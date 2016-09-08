"""
Module for recoding matrices.
date: 25/08/2016
"""
import sys
import os
import types
import subprocess
import re
import copy
import warnings
from Bio.Data import CodonTable
from p4 import Alignment
from p4 import func
from p4 import var
from p4 import P4Error
from p4 import read
from p4.code_utils import \
        codon_position, \
        codon_position_is_degenerate, \
        codon_slice_has_aas, \
        codon_slice_is_constant_aa, \
        codon_slice_is_degenerate, \
        codons_from_triplet_slice, \
        degenerate_codon_slice, \
        getBiopythonCode, nuc2val, \
        recode_sequence, \
        reduce_by_or, \
        val2nuc, \
        Code, R, Y

def formatwarning(message, category, filename, lineno, line):
    return "%s:%s: %s:%s" % (filename, lineno, category.__name__, message)
warnings.formatwarning = formatwarning

CAT = "".join

# This method uses a generalized codon size, but this is not the case of everything in this module.
def getCodonPositionMask(self, pos, codon_size=3):
    """This method returns a mask corresponding to sites at *pos*-th codon position.
    *codon_size* can be set to specify that the codons have a different size than 3."""
    def pos_mask(i):
        if (i + codon_size - pos) % codon_size == 0:
            return "1"
        else:
            return "0"
    return CAT([pos_mask(i+1) for i in range(self.length)])
Alignment.getCodonPositionMask = getCodonPositionMask


def getDegenerateSitesMask(self, transl_table=1, code=None, all_3rd_positions=False):
    """This method returns a mask corresponding to sites contributing to codon degeneracy.
    This is intended to be used for submatrix extraction using the noLRSall3 method,
    using :meth:`Alignment.getSubsetUsingMask` (with the option *inverse=True*, to get
    the degeneracy-free sites).

    If *all_3rd_positions* is set to True, then the mask includes all 3rd codon positions
    regardless of their effective contribution to codon degeneracy.

    The matrix is expected to start at a first codon position and stop at a third
    codon position.

    *transl_table* is an integer used to determine under which genetic code
    the codons are to be interpreted. The default value of 1 corresponds to the
    standard genetic code. Other values can be found in p4.GeneticCode.py

    Alternatively, the genetic code can be provided directly, using a
    dictionnary *code* whose keys are codons, and the values are the
    corresponding amino-acids. All triplets present in the matrix should also
    be present in the code dictionnary, except triplets of indels. Codons and
    the corresponding amino-acids are expected to be in lower case.
    If such a code is provided, the value of *transl_table* is ignored.

    The name of this method noLRSall3 comes from its effect in the case of the
    standard genetic code: it discards the sites participating in first
    position degeneracy for leucine (L) and arginine (R), first and second
    position degeneracy for serine (S), as well as all third codon positions
    where degeneracy is observed (or all of them if *all_3rd_positions* is True).
    Depending on the genetic code used, the type of amino-acid affected could
    be different.

    The goal of the submatrix extraction using the produced mask is to remove
    the sites that could have been affected by composition bias: mutations
    within a set of synonymous codons are more likely to favour the codons that
    conform to the general nucleotide composition.  However, one could argue
    that this bias is less likely to have played when the observed codons
    differ by more than one nucleotide and at least a non-synonymous mutation
    has to occur to bridge the gap. With the standard genetic code, this occurs
    for serine codons.  Indeed, the minimal mutation paths connecting the
    serine AGY and TCN codon categories are
    AGY (serine) <-> TGY (cysteine) <-> TCY (serine)
    and
    AGY (serine) <-> ACY (threonine) <-> TCY (serine)

    The current implementation (as of june 2012) does not check that a
    mutational path between synonymous codons exists, that consists only in
    synonymous point mutations. This may be considered as a bug, because you
    may not want AGY and TCN (or other similar cases that could occur with
    different genetic codes) to be considered as a single degeneracy continuum.
    """

    gm = ["Alignment.getDegenerateSitesMask()"]

    if code is None:
        #code = GeneticCode(transl_table).code
        # Use the generalized Code class defined in code_utils.py
        code = Code(transl_table).code

    n_codons = self.length / 3
    mask = ""
    # Loop over the successive triplets of sites.
    for c in range(n_codons):
        # 3 alignment slices. One for each codon position.
        slices = [self.sequenceSlice((3 * c) + pos-1) for pos in [1, 2, 3]]
        # The different codons found for the current triplet of sites.
        codons = set([codon.lower() for codon in ["%s%s%s" % nnn for nnn in zip(
            slices[0], slices[1], slices[2])]])
        # These are not Codon instances, this probably doesn't deal properly with ambiguity codes.
        # Record the amino-acids coded at the 3 nucleotides site, and the codons used for this aa.
        aas_codons = {}
        for codon in codons:
            # Determine the corresponding amino-acid.
            if codon == '---':
                aa = '-'
            elif code.has_key(codon):
                aa = code[codon]
            elif 'n' in codon:
                # This is a simplification. Some "degenerate" codons
                # can still code an unambiguous amino-acid.
                aa = 'x'
            else:
                gm.append("Codon %s is not defined in the chosen code "
                          "or translation table." % codon)
                gm.append("%s" % str(code))
                raise P4Error(gm)
            # Record the codon used for the aa.
            if aas_codons.has_key(aa):
                aas_codons[aa].append(codon)
            else:
                aas_codons[aa] = [codon]
        # Determine which positions in the triplet are degenerate.
        codon_mask = [False, False, False]
        # Loop over the recorded amino-acids.
        for aa in aas_codons.keys():
            if len(aas_codons[aa]) > 1:
                # Several codons have been found at this triplet for the amino-acid aa.
                # For each position, count the number of different nucleotides
                # present in the used codons.
                degeneracy = [len(set([cod[0] for cod in aas_codons[aa]])),
                              len(set([cod[1] for cod in aas_codons[aa]])),
                              len(set([cod[2] for cod in aas_codons[aa]]))]
                if all_3rd_positions:
                    # Put a position in the mask if it is already in the mask
                    # or if it is degenerate, or if it is a 3rd position.
                    codon_mask = [codon_mask[pos-1] or (degeneracy[pos-1] > 1)
                                  for pos in [1, 2]] + [True]
                else:
                    # Put a position in the mask if it is already in the mask
                    # or if it is degenerate.
                    codon_mask = [codon_mask[pos-1] or (degeneracy[pos-1] > 1)
                                  for pos in [1, 2, 3]]
            if all(codon_mask):
                # All positions of the triplet have been found to contribute to
                # some codon degeneracy somewhere in the alignment.
                # There is no need to search further.
                break
        # Append the codon mask to the mask.
        mask += CAT(map(lambda b: "1" if b else "0", codon_mask))
    return mask
Alignment.getDegenerateSitesMask = getDegenerateSitesMask


def pseudoTranslate(self, transl_table=1, out_type="standard", code=None):
    """Returns a pseudo protein alignment from *self*, a DNA alignment.
    The result is of datatype standard instead of protein, which allows
    the use of special recodings, like distinguishing between two types
    of serines, like in :meth:`Alignment.recode23aa()`.

    *self* is translated using :attribute:`Code(transl_table).code`.

    Alternatively, the genetic code can be provided through the parameter *code*.
    If such a code is provided, the value of *transl_table* is ignored.
    The parameter *code* can take to types of values:
    1) It can be a string naming the code to use, as defined in Biopython's
    `CodonTable.unambiguous_dna_by_name.keys()`
    2) It can be a dictionnary *code* whose keys are codons, and the values are
    the corresponding amino-acids. All triplets present in the matrix should
    also be present in the code dictionnary, except triplets of indels. Codons
    and the corresponding amino-acids are expected to be in lower case.
    It may be possible to use a code based on another codon length as 3,
    but this has not been tested as of June 2012.


    At the moment, we can only do translations where the sequences are phased
    with the coding frame, ie the first sequence position is the first position
    of the codon, and the last sequence position should be a last codon position.

    The default behaviour is to use translation table 1, that is the standard genetic code.
    Other available translation tables, this week::

        if transl_table == 1: # standard
        elif transl_table == 2: # vertebrate mito
        elif transl_table == 4: # Mold, Protozoan,
                                # and Coelenterate Mitochondrial Code
                                # and the Mycoplasma/Spiroplasma Code
        elif transl_table == 5: # invertebrate mito
        elif transl_table == 9: # echinoderm mito

        and now 6, 10, 11, 12, 13, 14, 21.

    (These are found in p4.GeneticCode.py or in :class:`Code`)

    *transl_table* may also be provided as text consisting in blank-separated elements.
    Each elements consists in n characters, where n is the number of defined codons.
    The first element lists the coded (pseudo-)amino-acids.
    The second elements describes whether a codon can be a start codon ('M') or not ('-').
    The other elements correspond to the (pseudo-)nucleotides at the successive codon positions.
    Example::
        FFJJZZZZYY**CC*WBBBBPPPPHHQQUUUUIIIMTTTTNNKKXXOOVVVVAAAADDEEGGGG
        ---M---------------M------------MMMM---------------M------------
        TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
        TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
        TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

    """

    gm = ['p4.alignment_recoding.pseudoTranslate()']
    if self.dataType != 'dna':
        gm.append("Self should be a DNA alignment")
        raise P4Error(gm)

    if code is None:
        #from GeneticCode import Code
        code = Code(transl_table, in_type="dna", out_type=out_type).code
        codelength = Code(transl_table).codelength
    else:
        if isinstance(code, types.StringType):
            code = getBiopythonCode(code) # defined in code_utils.py
        # We assume that the "codons" have all the same length,
        # and we look at the first codon in the dictionary to know this length.
        codelength = len(code.keys()[0])
        # We use standard type, because, depending on the code used to make the translation,
        # we may get something that contains symbols not corresponding to normal amino-acids.
        out_type = "standard"

    if self.length % codelength != 0:
        gm.append("The length of self should be a multiple of %i" % codelength)
        raise P4Error(gm)

    ali = self.dupe()
    ali.dataType = out_type
    ali.length = self.length / codelength
    ali.symbols = CAT(sorted(set(code.values())))
    ali.equates = {}
    ali.dim = len(ali.symbols)
    ali.nexusSets = None
    ali.parts = []
    ali.excludeDelete = None
    for seq in ali.sequences:
        # Initialize an all-gap sequence.
        seq.sequence = ['-'] * ali.length
        seq.dataType = out_type

    for i in range(len(self.sequences)):
        # the original sequence
        dnaSeq = self.sequences[i].sequence
        # the future pseudo-translation
        pseudoProtSeq = ali.sequences[i].sequence
        for j in range(ali.length):
            theCodon = dnaSeq[(j * codelength):((j+1) * codelength)]
            if code.has_key(theCodon):
                pseudoProtSeq[j] = code[theCodon]
            elif theCodon == '-' * codelength:
                # full indel
                pseudoProtSeq[j] = '-'
            elif theCodon.count('-'):
                # partial indel
                gm.append("    seq %i, position %4i, dnaSeq %4i, codon '%s' is incomplete" % (
                    i, j, (j*codelength), theCodon))
                raise P4Error(gm)
            else:
                # Should we use a CodonTranslationError (defined in code_utils.py) here ?
                gm.append("    seq %i position %4i, dnaSeq %4i, codon '%s' is not a known codon" % (
                    i, j, (j*codelength), theCodon))
                raise P4Error(gm)

    for seq in ali.sequences:
        # Convert from list to string.
        #s.sequence = string.join(s.sequence, '')
        seq.sequence = CAT(seq.sequence)
        #print s.sequence
    return ali
Alignment.pseudoTranslate = pseudoTranslate

def recode23aa(self):
    """
    This method gives a pseudo-translation of *self* where leucine, arginine and
    serine are coded differently depending on the codon category
    (CTN -> B, TTR -> J, AGR -> O, CGN -> U, AGY -> X, TCN -> Z)

    The original letters R, S and L are not used.

    Current implementation is based on the plastid/bacteria genetic code.
    The results should be valid also if the standard genetic code is assumed.
    """
    table = """FFJJZZZZYY**CC*WBBBBPPPPHHQQUUUUIIIMTTTTNNKKXXOOVVVVAAAADDEEGGGG
---M---------------M------------MMMM---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""
    return self.pseudoTranslate(transl_table=table)
Alignment.recode23aa = recode23aa


def triplet_slice(self, pos):
    """This method returns a list of 3 successive sequence slices from *self*,
    starting at position *pos*."""
    return [self.sequenceSlice(pos+i) for i in range(3)]
Alignment.triplet_slice = triplet_slice

# generalization of triplet_slice, not tested
def nuplet_slice(self, pos, n=3):
    """This method returns a list of n successive sequence slices from *self*,
    starting at position *pos*."""
    return [self.sequenceSlice(pos+i) for i in range(n)]
Alignment.nuplet_slice = nuplet_slice


def iter_codon_slices(self, code):
    """This method iterates over columns of codons of *self*."""
    for i in xrange(0, self.length, 3):
        triple = self.triplet_slice(i)
        yield codons_from_triplet_slice(triple, code)
Alignment.iter_codon_slices = iter_codon_slices

# generalization, not tested
#def iter_codon_slices(self, code, n=3):
#    """This method iterates over columns of codons of *self*."""
#    for i in xrange(0, self.length, n):
#        nuple = self.nuplet_slice(i, n)
#        yield codons_from_nuplet_slice(nuple, code)
#Alignment.iter_codon_slices = iter_codon_slices


def getHasAAsMask(self, aas, transl_table=1, code=None):
    """

    This method returns a mask corresponding to the triplets of sites where
    at least one amino-acid present in *aas* is coded. This is intended
    to be used for submatrix extraction using :meth:`Alignment.getSubsetUsingMask`
    (with the option *inverse=True*, to get the sites with none of the amino-acids
    present in *aas*).

    The matrix is expected to start at a first codon position and stop at a third
    codon position.

    *transl_table* is an integer used to determine under which genetic code
    the codons are to be interpreted. The default value of 1 corresponds to the
    standard genetic code. Other values can be found in p4.GeneticCode.py

    Alternatively, the genetic code can be provided through the parameter *code*.
    If such a code is provided, the value of *transl_table* is ignored.
    The parameter *code* can take to types of values:
    1) It can be a string naming the code to use, as defined in Biopython's
    `CodonTable.unambiguous_dna_by_name.keys()`
    2) It can be a dictionnary *code* whose keys are codons, and the values are
    the corresponding amino-acids. All triplets present in the matrix should
    also be present in the code dictionnary, except triplets of indels. Codons
    and the corresponding amino-acids are expected to be in lower case.
    Alternatively, the genetic code can be provided directly, using a
    dictionnary *code* whose keys are codons, and the values are the
    corresponding amino-acids. All triplets present in the matrix should also
    be present in the code dictionnary, except triplets of indels. Codons and
    the corresponding amino-acids are expected to be in lower case.

    """

    gm = ["Alignment.getHasAAsMask()"]

    if code is None:
        code = Code(transl_table).code
    elif isinstance(code, types.StringType):
        code = getBiopythonCode(code)
    else:
        msg = "code must be a dictionary, or a string naming the code in Biopython."
        assert isinstance(code, dict), msg

    return CAT(map(
        lambda c_slice: "111" if codon_slice_has_aas(c_slice, aas) else "000",
        self.iter_codon_slices(code)))
Alignment.getHasAAsMask = getHasAAsMask


def getConstantAAMask(self, transl_table=1, code=None, restrict_to=[]):
    """

    This method returns a mask corresponding to the triplets of sites where
    only one amino-acid is coded. This is intended to be used for submatrix
    extraction using :meth:`Alignment.getSubsetUsingMask` (with the option
    *inverse=True*, to get the non-constant sites).

    If *restrict_to* is not empty, only those sites that are constant and code
    for one of the amino-acids the one-letter code of which is in *restrict_to*
    will be considered constant.

    The matrix is expected to start at a first codon position and stop at a third
    codon position.

    *transl_table* is an integer used to determine under which genetic code
    the codons are to be interpreted. The default value of 1 corresponds to the
    standard genetic code. Other values can be found in p4.GeneticCode.py

    Alternatively, the genetic code can be provided through the parameter *code*.
    If such a code is provided, the value of *transl_table* is ignored.
    The parameter *code* can take to types of values:
    1) It can be a string naming the code to use, as defined in Biopython's
    `CodonTable.unambiguous_dna_by_name.keys()`
    2) It can be a dictionnary *code* whose keys are codons, and the values are
    the corresponding amino-acids. All triplets present in the matrix should
    also be present in the code dictionnary, except triplets of indels. Codons
    and the corresponding amino-acids are expected to be in lower case.
    Alternatively, the genetic code can be provided directly, using a
    dictionnary *code* whose keys are codons, and the values are the
    corresponding amino-acids. All triplets present in the matrix should also
    be present in the code dictionnary, except triplets of indels. Codons and
    the corresponding amino-acids are expected to be in lower case.

    """

    gm = ["Alignment.getConstantAAMask()"]

    if code is None:
        code = Code(transl_table).code
    elif isinstance(code, types.StringType):
        code = getBiopythonCode(code) # defined in code_utils.py
    else:
        msg = "code must be a dictionary, or a string naming the code in Biopython."
        assert isinstance(code, dict), msg

    return CAT(
        map(lambda c_slice: "111" if codon_slice_is_constant_aa(c_slice, restrict_to) else "000",
            self.iter_codon_slices(code)))
Alignment.getConstantAAMask = getConstantAAMask


def getDegenerateCodonsMask(self, transl_table=1, code=None, restrict_to=[], ignore=[]):
    """

    This method returns a mask corresponding to the triplets of sites where
    degeneracy has been observed. This is intended to be used for submatrix
    extraction using :meth:`Alignment.getSubsetUsingMask` (with the option
    *inverse=True*, to get the sites with no degenerate codons).

    If *restrict_to* is not empty, only those amino-acid the one-lettre code of
    which is in *restrict_to* are considered.

    If *ignore* is not empty, only those amino-acid the one-lettre code of
    which is not in *ignore* are considered.

    The matrix is expected to start at a first codon position and stop at a third
    codon position.

    *transl_table* is an integer used to determine under which genetic code
    the codons are to be interpreted. The default value of 1 corresponds to the
    standard genetic code. Other values can be found in p4.GeneticCode.py

    Alternatively, the genetic code can be provided through the parameter *code*.
    If such a code is provided, the value of *transl_table* is ignored.
    The parameter *code* can take to types of values:
    1) It can be a string naming the code to use, as defined in Biopython's
    `CodonTable.unambiguous_dna_by_name.keys()`
    2) It can be a dictionnary *code* whose keys are codons, and the values are
    the corresponding amino-acids. All triplets present in the matrix should
    also be present in the code dictionnary, except triplets of indels. Codons
    and the corresponding amino-acids are expected to be in lower case.
    Alternatively, the genetic code can be provided directly, using a
    dictionnary *code* whose keys are codons, and the values are the
    corresponding amino-acids. All triplets present in the matrix should also
    be present in the code dictionnary, except triplets of indels. Codons and
    the corresponding amino-acids are expected to be in lower case.

    """

    gm = ["Alignment.getDegenerateCodonsMask()"]

    if code is None:
        #code = GeneticCode(transl_table).code
        code = Code(transl_table).code
    elif isinstance(code, types.StringType):
        code = getBiopythonCode(code) # defined in code_utils.py
    else:
        msg = "code must be a dictionary, or a string naming the code in Biopython."
        assert isinstance(code, dict), msg

    # Experiments to test the speed of execution.
    #mask= "".join(("111" if codon_slice_is_degenerate(cod_slice, restrict_to) else "000") for cod_slice in self.iter_codon_slices(code))
    #return mask
    #return "".join(("111" if codon_slice_is_degenerate(cod_slice, restrict_to) else "000") for cod_slice in self.iter_codon_slices(code))
    #mask = ""
    # Loop over the successive triplets of sites.
    #for codon_slice in self.iter_codon_slices(code):
    #for codon_slice in [codons_from_triplet_slice(self.triplet_slice(i), code) for i in xrange(0, self.length, 3)]:
    #for i in xrange(0, self.length, 3):
    #    # i is positioned at the first codon position of the triplet.
    #    codon_slice = codons_from_triplet_slice(self.triplet_slice(i), code)
    #    if codon_slice_is_degenerate(codon_slice, restrict_to):
    #        mask += "111"
    #    else:
    #        mask += "000"
    #return "".join(("111" if codon_slice_is_degenerate(cod_slice, restrict_to) else "000") for cod_slice in [codons_from_triplet_slice(self.triplet_slice(i), code) for i in range(0, self.length, 3)])
    #mask= "".join(("111" if codon_slice_is_degenerate(cod_slice, restrict_to) else "000") for cod_slice in self.iter_codon_slices(code))
    #mask = "".join(map(lambda c_slice : "111" if codon_slice_is_degenerate(c_slice, restrict_to) else "000", self.iter_codon_slices(code)))
    #return "".join(map(lambda c_slice : "111" if codon_slice_is_degenerate(c_slice, restrict_to) else "000", self.iter_codon_slices(code)))
    #mask = "".join(("111" if codon_slice_is_degenerate(cod_slice, restrict_to) else "000") for cod_slice in [codons_from_triplet_slice(self.triplet_slice(i), code) for i in xrange(0, self.length, 3)])
    #mask = "".join(map(lambda c_slice : "111" if codon_slice_is_degenerate(c_slice, restrict_to) else "000", [codons_from_triplet_slice(self.triplet_slice(i), code) for i in xrange(0, self.length, 3)]))
    #return "".join(map(lambda c_slice : "111" if codon_slice_is_degenerate(c_slice, restrict_to) else "000", [codons_from_triplet_slice(self.triplet_slice(i), code) for i in range(0, self.length, 3)]))
    #return mask
    #return "".join(map(lambda c_slice : "111" if codon_slice_is_degenerate(c_slice, restrict_to) else "000", self.iter_codon_slices(code)))
    return CAT(map(
        lambda c_slice: "111" if codon_slice_is_degenerate(
            c_slice, restrict_to, ignore) else "000",
        self.iter_codon_slices(code)))
Alignment.getDegenerateCodonsMask = getDegenerateCodonsMask

def getDegenerateSiteMaskForPos(self, pos, transl_table=1, code=None, restrict_to=[], ignore=[]):
    """

    This method returns a mask corresponding to the sites where degeneracy has
    been observed if they correspond to a *pos*-th codon position.
    This is intended to be used for submatrix extraction using
    :meth:`Alignment.getSubsetUsingMask` (with the option *inverse=True*, to
    get the degeneracy-free sites).

    If *restrict_to* is not empty, only those amino-acids
    the one-lettre code of which is in *restrict_to* are considered.

    If *ignore* is not empty, only those amino-acid the one-lettre code of
    which is not in *ignore* are considered.

    The matrix is expected to start at a first codon position and stop at a third
    codon position.

    *transl_table* is an integer used to determine under which genetic code
    the codons are to be interpreted. The default value of 1 corresponds to the
    standard genetic code. Other values can be found in p4.GeneticCode.py

    Alternatively, the genetic code can be provided through the parameter *code*.
    If such a code is provided, the value of *transl_table* is ignored.
    The parameter *code* can take to types of values:
    1) It can be a string naming the code to use, as defined in Biopython's
    `CodonTable.unambiguous_dna_by_name.keys()`
    2) It can be a dictionnary *code* whose keys are codons, and the values are
    the corresponding amino-acids. All triplets present in the matrix should
    also be present in the code dictionnary, except triplets of indels. Codons
    and the corresponding amino-acids are expected to be in lower case.
    Alternatively, the genetic code can be provided directly, using a
    dictionnary *code* whose keys are codons, and the values are the
    corresponding amino-acids. All triplets present in the matrix should also
    be present in the code dictionnary, except triplets of indels. Codons and
    the corresponding amino-acids are expected to be in lower case.

    """

    gm = ["Alignment.getDegenerateSiteMaskForPos()"]

    if code is None:
        code = Code(transl_table).code
    elif isinstance(code, types.StringType):
        code = getBiopythonCode(code) # defined in code_utils.py
    else:
        msg = "code must be a dictionary, or a string naming the code in Biopython."
        assert isinstance(code, dict), msg

    def pos_mask(i):
        if i == pos:
            return "1"
        else:
            return "0"
    def triplet_mask(selected):
        if selected:
            return CAT(map(pos_mask, [1, 2, 3]))
        else:
            return "000"
    # Iterate over the slices, find the triplets that will be included in the mask
    # (those where degeneracy occurs), generate the corresponding mask portions,
    # and join the mask portions to make the matrix mask.
    return CAT(map(
        triplet_mask, [codon_position_is_degenerate(
            cod_slice, pos, restrict_to, ignore) for cod_slice in self.iter_codon_slices(code)]))
Alignment.getDegenerateSiteMaskForPos = getDegenerateSiteMaskForPos


def degenerate(self, code="Standard", positions=[1, 2, 3], restrict_to=[], ignore=[], sub_code=None):
    """
    This method returns a copy of *self* where the codons are replaced with degenerate versions.
    If *restrict_to* is not empty, only those codons that code amino-acids listed in *restrict_to*
    will be degenerated.
    If *ignore* is not empty, those codons that code amino-acids listed in *ignore*
    will not be degenerated.
    *positions* determines which codon positions are degenerated. By default, the whole codons are
    degenerated (if there is degeneracy of course).
    *code* is the Biopython name of the genetic code under which degeneracy has to be interpreted
    or a dictionary converting from codons to amino-acids (all in lower case).
    Default is to use the standard genetic code. Possible values for *code* are:
    %s
    *sub_code*, if provided, should be a dictionary associating amino-acids to codons
    (all in lower case). For the purpose of defining degeneracy groups, the codons present
    in *sub_code* will be considered as coding for the amino-acid defined there instead of
    the one defined by *code*. This can be used for instance to keep two distinct types of
    serine codons, with degeneracy only within each type. The codons still count as coding
    their original amino-acid with respect to the *restrict_to* and *ignore* options.
    """ % "\n".join(sorted(CodonTable.unambiguous_dna_by_name.keys()))
    if isinstance(code, types.StringType):
        code = getBiopythonCode(code) # defined in code_utils.py
    else:
        msg = "code must be a dictionary, or a string naming the code in Biopython."
        assert isinstance(code, dict), msg
    # codons belonging to different sub-groups of codons for one amino-acid
    # can be considered as coding different amino-acids
    # (sub-amino-acids of the "normal" amino-acid, for instance two types of serine)
    if sub_code is None:
        sub_code = {}
    else:
        assert isinstance(sub_code, dict), "sub_code must be a dictionary."
        sub_code = copy.copy(sub_code) # otherwise there are side effects:
                                       # the content of sub_code can be modified
                                       # in the calling context
        if any([sub_aa in code.values() for sub_aa in sub_code.values()]):
            msg = CAT(["Note that at least one sub-aminoacid provided in sub_code ",
                       "is identical to an amino-acid provided by the chosen genetic code.\n",
                       "The sub-amino-acids are:\n%s\n" % ", ".join(
                           [str(aa) for aa in sub_code.values()])])
            warnings.warn(msg)
    # Ensure the amino-acids are in lowercase.
    restrict_to = set([aa.lower() for aa in restrict_to])
    ignored_aas = set([aa.lower() for aa in ignore])
    # Find the groups of synonymous codons.
    # The keys are amino-acids, the values are lists of codons that code the amino-acid.
    aas_codons = {}
    for codon in code.keys():
        aa = code[codon]
        if not sub_code.has_key(codon):
            sub_code[codon] = aa # sub_aa will be the same as aa
        sub_aa = sub_code[codon]
        # Only consider codons that are compatible with the restriction rule, if there is one.
        if (len(restrict_to) == 0 or aa.lower() in restrict_to) and not (aa.lower() in ignored_aas) :
            #if aas_codons.has_key(aa):
            if aas_codons.has_key(sub_aa):
                #aas_codons[aa].append(codon)
                aas_codons[sub_aa].append(codon)
            else:
                #aas_codons[aa] = [codon]
                aas_codons[sub_aa] = [codon]
    # Build a conversion dictionary.
    # The keys are the codons, the values their degenerate replacements.
    cod2degen = {}
    for codons in aas_codons.values():
        # Compute degeneracy values at the 3 positions
        # The degenerate value at a position is the binary union
        # of the values of the nucleotides found at that position.
        # nuc2val and reduce_by_or are defined in code_utils.py
        degen1 = reduce_by_or([nuc2val[cod[0]] for cod in codons])
        degen2 = reduce_by_or([nuc2val[cod[1]] for cod in codons])
        degen3 = reduce_by_or([nuc2val[cod[2]] for cod in codons])
        # Compute the string representation of the resulting degenerate codon.
        # val2nuc is defined in code_utils.py
        degenerate_codon = val2nuc[degen1] + val2nuc[degen2] + val2nuc[degen3]
        # Associate this representation to all the synonymous codons it represents.
        for cod in codons:
            cod2degen[cod.lower()] = degenerate_codon.lower()
            # If restrict_to is not empty, it is likely that not all codons
            # are present in cod2degen, but the code_utils.recode_sequence function
            # will just keep those codons as is.
    # Make a copy of self.
    ali = self.dupe()
    for seq in ali.sequences:
        # Recode the sequence using the conversion dictionary built previously.
        # recode_sequence is defined in Codon_utils.py
        seq.sequence = recode_sequence(seq.sequence, cod2degen, positions, code=code)
    return ali
Alignment.degenerate = degenerate


def recodeRY(self, positions=[1, 2, 3]):
    """
    This method returns a copy of *self* where purines are replaced with IUPAC ambiguity code R
    and pyrimidines are replaced with IUPAC ambiguity code Y.
    *positions* determines which codon positions are degenerated. By default, the whole codons are
    degenerated (if there is degeneracy of course).
    """
    recode_table = {}
    # Make a copy of self.
    ali = self.dupe()
    for seq in ali.sequences:
        new_seq = []
        pos = 0
        while pos < len(seq.sequence):
            # codon_position is defined in code_utils.py
            if codon_position(pos + 1) in positions:
                letter = seq.sequence[pos]
                if not recode_table.has_key(letter):
                    # nuc2val, R and Y are defined in code_utils.py
                    val = nuc2val[letter]
                    if (val & R) and (val & Y):
                        # letter is an ambiguity code representing both purines and pyrimidines.
                        recode_table[letter] = "n"
                    elif val & R:
                        recode_table[letter] = "r"
                    elif val & Y:
                        recode_table[letter] = "y"
                    else:
                        msg = "Letter %s should be the code for a gap ('-')." % letter
                        assert letter == "-", msg
                        recode_table[letter] = "-"
                new_seq.append(recode_table[letter])
            else:
                # For consistency, all characters re-written in lower case.
                new_seq.append(seq.sequence[pos].lower())
            pos += 1
        seq.sequence = CAT(new_seq)
    return ali
Alignment.recodeRY = recodeRY


def degenerateByCodonColumn(self, code="Standard", restrict_to=[]):
    """
    This method returns a copy of *self* where codons coding for the same
    amino-acid in a given column of the matrix are replaced by their union
    (i.e. degenerated), contrary to `degenerate` which does this regardless
    of the codons present, just using the degeneracy observed in the genetic code.
    If *restrict_to* is not empty, only those codons that code amino-acids listed
    in *restrict_to* will be degenerated.
    *code* is the Biopython name of the genetic code under which degeneracy has
    to be interpreted or a dictionary converting from codons to amino-acids.
    Default is to use the standard genetic code. Possible values for *code* are:
    %s
    """ % "\n".join(sorted(CodonTable.unambiguous_dna_by_name.keys()))
    if isinstance(code, types.StringType):
        code = getBiopythonCode(code) # defined in code_utils.py
    else:
        msg = "code must be a dictionary, or a string naming the code in Biopython."
        assert isinstance(code, dict), msg
    # Ensure the amino-acids are in lowercase.
    restrict_to = set([aa.lower() for aa in restrict_to])
    # The matrix will be rebuilt column-wise,
    # and then the sequences will be rebuilt from these columns of codons.
    new_slices = []
    for cod_slice in iter_codon_slices(self, code):
        new_slices.append(degenerate_codon_slice(cod_slice, restrict_to))
    # Make a copy of self.
    ali = self.dupe()
    # Loop over the sequences.
    for i in range(self.nChar):
        ali.sequences[i].sequence = CAT([str(cod_slice[i]) for cod_slice in new_slices])
    return ali
Alignment.degenerateByCodonColumn = degenerateByCodonColumn


def indelizeCodons(self, aas, code="Standard"):
    """
    This method returns a copy of *self* where the codons corresponding to
    amino-acids listed in *aas* are replaced by indels.
    *code* is the Biopython name of the genetic code under which degeneracy has to be interpreted
    or a dictionary converting from codons to amino-acids.
    Default is to use the standard genetic code. Possible values for *code* are:
    %s
    """ % "\n".join(sorted(CodonTable.unambiguous_dna_by_name.keys()))
    if isinstance(code, types.StringType):
        code = getBiopythonCode(code) # defined in code_utils.py
    else:
        msg = "code must be a dictionary, or a string naming the code in Biopython."
        assert isinstance(code, dict), msg
    # Ensure the amino-acids are in lowercase.
    aas = set([aa.lower() for aa in aas])
    # Build a conversion dictionary.
    # The keys are the codons, the values their "indelized" replacements.
    cod2indels = {}
    for codon in code.keys():
        if code[codon] in aas:
            cod2indels[codon] = "---"
        else:
            cod2indels[codon] = codon
    # Make a copy of self.
    ali = self.dupe()
    for seq in ali.sequences:
        # Recode the sequence using the conversion dictionary built previously.
        # recode_sequence is defined in Codon_utils.py
        seq.sequence = recode_sequence(seq.sequence, cod2indels, code=code)
    return ali
Alignment.indelizeCodons = indelizeCodons

def blend_matrices(*matrices):
    """
    This function returns an alignment that is made by taking its column in the
    matrices provided as arguments, one column each and cycling between the
    matrices, starting with the first one. All matrices should have the same length,
    and this length should be a multiple of the number of matrices.
    They should also have the same taxa, and the taxa should be in the same order
    in the different matrices.
    """
    mat_len = matrices[0].nChar
    msg = "All matrices should have the same length."
    assert all([mat.nChar == mat_len for mat in matrices[1:]]), msg
    msg = "The length of the matrices should be a multiple of the number of matrices to blend."
    assert mat_len % len(matrices) == 0, msg
    # Start with a copy of the first matrix.
    ali = matrices[0].dupe()
    for i in range(ali.nTax):
        ali.sequences[i].sequence = blend_sequences([m.sequences[i].sequence for m in matrices])
    return ali

def blend_sequences(sequences):
    """This function returns a chain of characters made by taking characters in turn from the
    chains provided as arguments. These chains should have the same length and this length
    should be a multiple of the number of chains."""
    n_seq = len(sequences)
    #seq = ""
    #i = 0
    #while i < len(sequences[-1]):
    #    seq += sequences[i % n_seq][i]
    #    i += 1
    #return seq
    return CAT([sequences[i % n_seq][i] for i in range(len(sequences[-1]))])

def treeFinderMAPAnalysis(alignment, groups,
                          gamma=True, invariant=True, bootstrap=False,
                          nreplicates=100,
                          remove_files=False, run_analysis=True, verbose=False):
    """
    Uses TreeFinder to estimate a Maximum Likelihood tree using the MAP
    substitution model for grouped amino-acids.

    - *alignment*: p4 alignment object of original (un-recoded) protein data from
      which the "groups" are derived
    - *groups*: list of grouped amino-acids, possibly resuling from
      :meth:`Alignment.getKosiolAISGroups()` or :meth:`Alignment.getMinmaxChiSqGroups()`
    - *gamma*: include gamma distribution of among-site rate variation
    - *bootstrap*: run bootstrap analysis
    - *nreplicates*: number of bootstrap replicates
    - *invariant*: include a proportion of invariant sites
    - *run_analysis*: run the analysis if TreeFinder in $PATH, else just write the
      control file
    - *remove_files*: remove analysis files. Only available if run_analysis=True

    """

    gm = ["p4.alignment_recoding.treeFinderMAPAnalysis()"]

    if not isinstance(alignment, Alignment):
        msg = "alignment must be a Alignment object"
        gm.append(msg)
        raise P4Error(gm)

    if alignment.dataType != "protein":
        msg = "alignment should be the original protein data from" + \
              "which the groups were defined. Doing nothing."
        gm.append(msg)
        raise P4Error(gm)

    for param in [gamma, invariant, bootstrap,
                  remove_files, run_analysis, verbose]:
        if not isinstance(param, types.BooleanType):
            msg = "%s value must be either True or False" % param
            gm.append(msg)
            raise P4Error(gm)

    if not isinstance(nreplicates, types.IntType):
        msg = "nreplictes must be an integer"
        gm.append(msg)
        raise P4Error(gm)

    if run_analysis:
        if not func.which2("tf"):
            msg = "tf (treefinder) is not in your $PATH" + \
                  "Cannot run analysis"
            gm.append(msg)
            raise P4Error(gm)

    datafile_name = "tf_data.phy"

    #tf commands
    tls = """ReconstructPhylogeny[
             "%(datafile)s",
             SubstitutionModel->MAP[%(map)s][Optimum,Optimum]%(ifH)s,
             WithEdgeSupport->%(bootstrap)s%(nreplicates)s
             ],
             "%(outfile)s",SaveReport"""
    od = {}
    od["datafile"] = datafile_name
    if gamma:
        if invariant:
            od["ifH"] = ":GI[Optimum]"
        else:
            od["ifH"] = ":G[Optimum]"
    else:
        if invariant:
            od["ifH"] = ":I[Optimum]"
        else:
            od["ifH"] = ""
    if bootstrap:
        od["bootstrap"] = "True"
        od["nreplicates"] = ",NReplicates->%i" % nreplicates
    else:
        od["bootstrap"] = "False"
        od["nreplicates"] = ""
    od["outfile"] = "tf_reconstruction.output"
    od["map"] = ",".join(['"%s"' % i for i in [group.upper() for group in groups]])

    if run_analysis:

        #Write data file
        alignment.writePhylip(datafile_name)

        #Write control file
        tl_file = "tf_control.tl"
        fh = open(tl_file, "w")
        fh.write(tls % od)
        fh.close()

        if verbose:
            direct = subprocess.STDOUT
        else:
            direct = open("/dev/null", "w")

        child = subprocess.Popen("tf tf_control.tl", stderr=direct, shell=True)

        if verbose:
            print "Running TreeFinder, this could take some time...",
            sys.stdout.flush()

        child.communicate()

        if verbose:
            print "done."
            sys.stdout.flush()

        #This doesnt seem to work, why?
        #while child.poll() is None:
        #    time.sleep(60)
        #    if verbose:
        #        sys.stdout.write(".")
        #        sys.stdout.flush()

        if child.returncode != 0:
            msg = "TreeFinder returned error code %s"
            gm.append(msg % (child.returncode))
            raise P4Error(gm)

        fh = open(od["outfile"], "r")
        line = fh.readlines()[1]
        fh.close()

        rd = {}
        #Likelihood
        rd["Likelihood"] = float(line[line.index("Likelihood->")+12:line.index(",")])
        #Tree
        ts = line[line.index("Phylogeny->")+11:line.index("SubstitutionModel->")-1]
        rd["Phylogeny"] = ts
        #SubstitutionModel
        sm = line[line.index("SubstitutionModel->")+19:line.index("OSubstitutionModel->")-1]
        rd["SubstitutionModel"] = sm
        #OSubstitutionModel
        osm = line[line.index("OSubstitutionModel->")+20:line.index("OEdgeOptimizationOff->")-1]
        rd["OSubstitutionModel"] = osm
        #NSites
        ns = line[line.index("NSites->")+8:line.index("NParameters->")-1]
        rd["Nsites"] = int(ns)
        #NParameters
        np = line[line.index("NParameters->")+13:line.index("AIC->")-1]
        rd["NParameters"] = int(np)
        #AIC
        rd["AIC"] = float(line[line.index("AIC->")+5:line.index("AICc->")-1])
        #AICc->
        rd["AICc"] = float(line[line.index("AICc->")+6:line.index("HQ->")-1])
        #HQ
        rd["HQ"] = float(line[line.index("HQ->")+4:line.index("BIC->")-1])
        #BIC
        rd["BIC"] = float(line[line.index("BIC->")+5:line.index("Checksum->")-1])
        #LikelihoodTime
        lt = line[line.index("LikelihoodTime->")+16:line.index("LikelihoodMemory->")-1]
        rd["LikelihoodTime"] = float(lt)
        #LikelihoodMemory
        lm = line[line.index("LikelihoodMemory->")+18:-3]
        rd["LikelihoodMemory"] = int(lm)

        #Make a tree object
        tree = rd["Phylogeny"].replace("{", "(")
        tree = tree.replace("}", ")")
        tree = tree.replace("\"", "")
        tree = tree + ";"
        if bootstrap:
            #Tree viewer has the brlen before bootstrap value plus an extra colon
            # turn "xxx):0.00001:87.999,yyy" into "xxx)87.999:0.00001,yyy"
            patt = re.compile(r"\):([0-9]+\.[0-9e-]+):([0-9]+\.[0-9e-]*)")
            repl = r")\2:\1"
            tree = re.sub(patt, repl, tree)
        origw = var.warnReadNoFile
        var.warnReadNoFile = False
        read(tree)
        var.warnReadNoFile = origw
        result_tree = var.trees.pop()
        if bootstrap:
            #Round up floats to percentages
            for node in result_tree.iterInternalsNoRoot():
                node.name = "%2.f" % float(node.name)

        if remove_files:
            os.remove("tf_control.tl")
            os.remove("tf_data.phy")
            os.remove("tf_reconstruction.output")

        if verbose:
            print "\n"
            result_tree.draw()
            print "\nLikelihood: %.4f\n" % rd["Likelihood"]

        return result_tree, rd

    else:
        print tls % od
        return (None, None)

