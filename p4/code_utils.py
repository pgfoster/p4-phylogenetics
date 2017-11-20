"""
Utilities relating to genetic codes.
date: 25/08/2016
"""
import sys
import copy
import warnings
from itertools import combinations, product
import operator
from Bio.Data import CodonTable
from p4 import P4Error

def formatwarning(message, category, filename, lineno, line):
    return "%s:%s: %s:%s" % (filename, lineno, category.__name__, message)
warnings.formatwarning = formatwarning

CAT = "".join

def reduce_by_or(the_list):
    """This function returns the 'or' of the elements of the list *the_list*.
    This can be used to obtain a degenerate codon representing a list of codons,
    or to merge a list of sets of elements into a single set."""
    return reduce(operator.or_, the_list)

# numerical values for nucleotides to help with degeneracy and codon sets representation.
A = 1
C = 2
G = 4
T = 8
# ambiguity codes get the values corresponding to the bitwise 'or' of the nucleotides
B = C | G | T
D = A | G | T
H = A | C | T
K = G | T
M = A | C
N = A | C | G | T
R = A | G
S = C | G
V = A | C | G
W = A | T
Y = C | T
nucleotides = (A, C, G, T, 0) # 0 is for the gap
# gives the value of a nucleotide letter, upper case and lower case are in this dictionary.
nuc2val = {'A':A, 'C':C, 'G':G, 'T':T,
           'B':B, 'D':D, 'H':H, 'K':K,
           'M':M, 'N':N, 'R':R, 'S':S,
           'V':V, 'W':W, 'Y':Y, '-':0,
           'a':A, 'c':C, 'g':G, 't':T,
           'b':B, 'd':D, 'h':H, 'k':K,
           'm':M, 'n':N, 'r':R, 's':S,
           'v':V, 'w':W, 'y':Y}
# reverse table: gives the upper case letter corresponding to a 'nucleotidic value'
val2nuc = {}
for k in nuc2val.keys():
    val2nuc[nuc2val[k]] = k.upper()

class CodonTranslationError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def getBiopythonCode(code_name):
    """
    returns a dictionary translating codons to amino-acids using the translation
    table extracted from Biopython's CodonTable using the name *code_name*.
    """
    table = CodonTable.unambiguous_dna_by_name[code_name]
    # Complete the forward table with stop codons.
    for stop in table.stop_codons:
        table.forward_table[stop] = '*'
    # Also add a 'gap codon' (but partially gapped codons are not included in the code)
    table.forward_table["---"] = '-'
    # Convert to lower case (because Biopython and p4 don't use the same conventions).
    for codon in table.forward_table.keys():
        table.forward_table[codon.lower()] = table.forward_table[codon].lower()
        if codon != codon.lower():
            # Delete the old entry that has just been converted.
            del table.forward_table[codon]
        # else it is a gap: don't delete it
    return table.forward_table


# Adapted from GeneticCode.GeneticCode
class Code(object):
    """A generalized container for translation and recoding tables,
    adapted from :class:`GeneticCode`.
    This provides:
    
    - ``code`` A dictionary.  So you can ask for eg myCode.code['ggg']
    - ``codonsForAA`` Another dictionary, where you can ask for eg myCode.codonsForAA['v']
    - ``startList`` A list of start codons
    - ``in_type`` The data type that is used to make the codons
    - ``out_type`` The data type into which the codons are translated
    - ``codelength`` The number of in_type elements needed to make a codon
    """

    def __init__(self, transl_table=1, in_type="dna", out_type="standard"):
        """*transl-table* is an integer determining which predefined code to use.
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
        
        Limitation: This way of providing the translation table only permits a
        pseudo-translation into a space of one-letter states, like proteins.
        
        """

        self.transl_table = transl_table
        self.code = {}
        self.codonsForAA = {}
        self.startList = []
        self.in_type = in_type
        self.out_type = out_type
        self.codelength = None

        if transl_table == 1: # standard
            AAs    = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '---M---------------M---------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"

        elif transl_table == 2: # vertebrate mito
            AAs      = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG'
            Starts   = '--------------------------------MMMM---------------M------------'
            Base1    = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2    = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3    = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"

        elif transl_table == 4: # Mold, Protozoan,
                                # and Coelenterate Mitochondrial Code
                                # and the Mycoplasma/Spiroplasma Code
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '--MM---------------M------------MMMM---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"

        elif transl_table == 5: # invertebrate mito
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG'
            Starts = '---M----------------------------MMMM---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"

        elif transl_table == 6: # The Ciliate, Dasycladacean
                                # and Hexamita Nuclear Code (transl_table=6)
            AAs    = 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '-----------------------------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"

        # tables 7 and 8 have been deleted from NCBI.
        
        elif transl_table == 9: # echinoderm and flatworm mito
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG'
            Starts = '-----------------------------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"

        elif transl_table == 10: # The Euplotid Nuclear Code (transl_table=10)
            AAs    = 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '-----------------------------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"


        elif transl_table == 11: # The Bacterial and Plant Plastid Code (transl_table=11)
            AAs    = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '---M---------------M------------MMMM---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"

        elif transl_table == 12: # The Alternative Yeast Nuclear Code (transl_table=12)
            AAs    = 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '-------------------M---------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"


        elif transl_table == 13: # The Ascidian Mitochondrial Code (transl_table=13)
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG'
            Starts = '---M------------------------------MM---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"

        elif transl_table == 14: # The Alternative Flatworm Mitochondrial Code (transl_table=14)
            AAs    = 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG'
            Starts = '-----------------------------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"
            
        elif transl_table == 21: # Trematode Mitochondrial Code (transl_table=21)
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG'
            Starts = '-----------------------------------M---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            Bases = [Base1, Base2, Base3]
            symbols = "ACGT"


        else:
            lines = transl_table.split()
            AAs = lines[0]
            Starts = lines[1]
            Bases = lines[2:]
            symbols = CAT(sorted(reduce_by_or(
                [set([c for c in Bases[i]]) for i in range(len(Bases))])))

        # Find the length of a codon.
        self.codelength = len(Bases)
        # Define the 'gap codon'.
        self.code['-' * self.codelength] = "-"

        for i in range(len(symbols) ** self.codelength):
            # Generalized version (can deal with any code length):
            theCodon = reduce(operator.add, [Bases[j][i] for j in range(self.codelength)]).lower()
            theAA = AAs[i].lower()
            self.code[theCodon] = theAA
            if theAA in self.codonsForAA:
                self.codonsForAA[theAA].append(theCodon)
            else:
                self.codonsForAA[theAA] = [theCodon]
            if Starts[i] == 'M':
                self.startList.append(theCodon)

# Not generalized: uses 3 codon positions
class Codon(object):
    """This object represents a (possibly degenerate) codon."""
    def __init__(self, codon, code="Standard"):
        """A Codon object is initiated by giving its string representation to __init__.
        A genetic code can be specified according to which the codon should be translated.
        *code* is the Biopython name of the genetic code under which degeneracy has to be
        interpreted, or a dictionary giving the translation of each non-degenerate codon.
        By default, the standard genetic code will be used.
        """
        for letter in codon:
            assert letter.upper() in nuc2val.keys(), "Unrecognised symbol %s.\n" % letter
        self.codon = codon.upper()
        # The internal representation uses one number for each position.
        # A position is degenerate if it is not represented by a power of 2.
        self.v1 = nuc2val[self.codon[0]]
        self.v2 = nuc2val[self.codon[1]]
        self.v3 = nuc2val[self.codon[2]]
        # Is the codon a degenerate codon ?
        if self.v1 in nucleotides and self.v2 in nucleotides and self.v3 in nucleotides:
            self.degenerate = False
            self._decomposition = [self]
        else:
            self.degenerate = True
            self._decomposition = None # The decomposition into non-degenerate
                                       # codons will have to be computed.
        # The maximum value at a given position is A + C + G + T = 15
        # The following formula should therefore give a unique id to a given codon.
        self.id = self.v3 + 16 * (self.v2 + 16 * self.v1)
        # TODO: Since codons are hashable / uniquely defined, time might be saved by checking
        # wether a given codon has already been created when an operation requires the creation
        # of a codon, and use that already defined codon. Some clever programming pattern might
        # apply here. Maybe use the MutationGraph object to store the codons?
        
        # dictionary to record mutational distances
        # The keys are codon id values.
        # The values are mutational distances between self and the codon
        # The dictionary starts with the distance to self (0)
        self._distances = {self.id : 0}
        # genetic code according to which the codon should be translated
        if isinstance(code, str):
            self.code = getBiopythonCode(code)
        else:
            msg = "code must be a dictionary, or a string naming the code in Biopython."
            assert isinstance(code, dict), msg
            self.code = code
        # set of amino-acid coded by the codon
        self.aas = None
        self.set_aas(self.code)

    def __str__(self):
        # The lower case vs upper case thing could use some cleaning of the code.
        # We use lower case for p4 compatibility.
        return self.codon.lower()

    def __or__(self, other):
        """The union of two codons is the degenerate codon representing
        the union of the sets of codons that the two codons represent."""
        or1 = self.v1 | other.v1
        or2 = self.v2 | other.v2
        or3 = self.v3 | other.v3
        return Codon(val2nuc[or1] + val2nuc[or2] + val2nuc[or3], self.code)

    def __and__(self, other):
        """The and operation is somewhat artificial in that it returns indels
        at positions with no nucleotide in common."""
        and1 = self.v1 & other.v1
        and2 = self.v2 & other.v2
        and3 = self.v3 & other.v3
        return Codon(val2nuc[and1] + val2nuc[and2] + val2nuc[and3], self.code)

    def __contains__(self, other):
        """A codon 'contains' another codon if it is a degenerate
        or identical version of the other."""
        cont1 = self.v1 & other.v1 == other.v1
        cont2 = self.v2 & other.v2 == other.v2
        cont3 = self.v3 & other.v3 == other.v3
        # This defines a partial order relationship on the codons:
        # not(c1.__contains__(c2)) does not imply c2.__contains__(c1)
        return cont1 and cont2 and cont3

    def __eq__(self, other):
        return self.id == other.id

    def __ne__(self, other):
        return self.id != other.id

    def __hash__(self):
        return self.id

    def __getitem__(self, i):
        """returns the *i*-th nucleotide value of *self*."""
        if i == 1:
            return self.v1
        elif i == 2:
            return self.v2
        elif i == 3:
            return self.v3
        else:
            raise IndexError, "A codon has only 3 positions."

    def decomposition(self):
        """returns the list of non-degenerate codons that are contained in *self*."""
        if self._decomposition is None:
            # It is not defined yet: Compute it.
            # Find the non-degenerate nucleotides at each position.
            # The indel does not count as a nucleotide.
            # (I use list and set conversions, that are maybe a little ugly, to remove duplicates.)
            dec1 = sorted(set([self.v1 & nuc for nuc in [A, C, G, T]]) - set([0]))
            dec2 = sorted(set([self.v2 & nuc for nuc in [A, C, G, T]]) - set([0]))
            dec3 = sorted(set([self.v3 & nuc for nuc in [A, C, G, T]]) - set([0]))
            # Generate the list of implied codons.
            decomposition = [Codon(
                val2nuc[cod[0]] + val2nuc[cod[1]] + val2nuc[cod[2]],
                self.code) for cod in product(dec1, dec2, dec3)]
            if decomposition:
                self._decomposition = decomposition
            else:
                # We don't want to return an empty list, so we return a 'gap codon' instead.
                self._decomposition = [Codon("---", self.code)]
        #else:
        #    pass # no need to re-compute it
        return self._decomposition

    def set_aas(self, code=None):
        """This method determines the set of amino-acids coded by the non-degenerate codons
        implied by self.
        A genetic code can be specified according to which the codon should be translated.
        *code* is the Biopython name of the genetic code under which degeneracy has to be
        interpreted, or a Code object. By default, the genetic code of the codon will be used.
        Alternatively, a dictionary can be provided. It should take codon as keys and their
        translation as values.
        """
        if code is None:
            # self.code should be something valid.
            code = self.code
        elif isinstance(code, str):
            code = getBiopythonCode(code)
        else:
            msg = "code must be a dictionary, or a string naming the code in Biopython."
            assert isinstance(code, dict), msg
        # We assume that the "codons" have all the same length,
        # and we look at the first codon in the dictionary to know this length.
        codelength = len(code.keys()[0])
        assert codelength == 3, "Amino-acids should be coded by triplets of nucleotides."
        try:
            if self.degenerate:
                self.aas = set([code[str(codon)] for codon in self.decomposition()])
            else:
                self.aas = set([code[str(self)]])
        except KeyError:
            raise CodonTranslationError("The code does not say what the translation "
                                        "of %s should be.\n" % str(self))

    def distance(self, other):
        """This method returns the mutational distance between *self* and *other*.
        It is defined as the minimum number of nucleotide substitutions to transform
        *self* into *other*."""
        if other.id in self._distances:
            # The distance has already been calculated.
            # It should be known also on the other side.
            # other._distances[self.id] = self._distances[other.id]
            msg = "The distances should be symmetrical."
            assert other._distances[self.id] == self._distances[other.id], msg
            #pass
        elif self.id in other._distances:
            # The distance is not known, but strangely it is known on the other side.
            warnings.warn("It is unexpected that the distance is already known "
                          "on the other side but not on the side of self.\n")
            self._distances[other.id] = other._distances[self.id]
        # With degenerate codons, we use the shortest mutational path between the represented
        # sets of non-degenerate codons.
        elif self.degenerate:
            # Recursion may happen during the calculation of other.distance().
            self._distances[other.id] = min([
                other.distance(cod.id) for cod in self.decomposition()])
            # This is not yet know. Fix it now.
            other._distances[self.id] = self._distances[other.id]
        elif other.degenerate:
            # self is not degenerate.
            self._distances[other.id] = min([
                self.distance(cod.id) for cod in other.decomposition()])
            # This is not yet know. Fix it now.
            other._distances[self.id] = self._distances[other.id]
        else:
            # Neither self nor other are degenerate. Simply count the necessary mutations.
            d = 0
            if self.v1 != other.v1:
                d += 1
            if self.v2 != other.v2:
                d += 1
            if self.v3 != other.v3:
                d += 1
            self._distances[other.id] = d
            # This is not yet know. Fix it now.
            other._distances[self.id] = d
        return self._distances[other.id]

def codon_position(n, codon_length=3):
    """
    This function returns the codon position corresponding to position *n* in a sequence
    (starting at 1), assuming that the codons have a length of *codon_length*.
    """
    pos = n % codon_length
    if pos == 0:
        return codon_length
    else:
        return pos

#TODO: compute mutational paths between codons
class MutationGraph(object):
    """This object represents a graph of codons (and contains also other
    useful stuff that I should document / reorganize.).
    Non-degenerate codons are neighbours if they differ by only one nucleotide.
    Degenerate codons are sets of non-degenerate codons."""
    def __init__(self, code="Standard"):
        if isinstance(code, str):
            self.code = getBiopythonCode(code)
        else:
            msg = "code must be a dictionary, or a string naming the code in Biopython."
            assert isinstance(code, dict), msg
            self.code = code
        # dictionary of non-degenerate codons
        # The keys are codon id values.
        # The values are non-degenerate codons.
        self.non_degen = {}
        # dictionary of distances between non-degenerate codons
        # The keys are pairs (cardinality 2 frozensets) of codons id values.
        # The values are the number of different nucleotides between the two codons.
        self.distances = {}
        # dictionary of edges (links between non-degenerate codons differing by one nucleotide)
        # The keys are pairs (cardinality 2 frozensets) of codons id values.
        # The values are 1,
        # no entry exists where distance is more or less than one nucleotide substitution.
        self.edges = {}
        for n1 in ['A', 'C', 'G', 'T']:
            for n2 in ['A', 'C', 'G', 'T']:
                for n3 in ['A', 'C', 'G', 'T']:
                    codon = Codon(n1+n2+n3, self.code)
                    self.non_degen[codon.id] = codon
                    for cod in self.non_degen.values():
                        pair = frozenset([cod.id, codon.id])
                        d = codon.distance(cod)
                        self.distances[pair] = d
                        if d == 0:
                            self.edges[pair] = 1
        # To give colours to certain codons, for drawing purposes.
        self.codon_colours = {}
        # To group codons by synonymy classes.
        self.degen_groups = {}
        for codon in self.non_degen.values():
            aa = list(codon.aas)[0]
            n1 = codon[1]
            n2 = codon[2]
            n3 = codon[3]
            if aa == 'l':
                if n1 == C:
                    self.codon_colours[codon] = "blue"
                elif n1 == T:
                    self.codon_colours[codon] = "red"
                else:
                    pass
            if aa == 'r':
                if n1 == A:
                    self.codon_colours[codon] = "yellow"
                elif n1 == C:
                    self.codon_colours[codon] = "green"
                else:
                    pass
            if aa == 's':
                if n1 == A:
                    self.codon_colours[codon] = "orange"
                elif n1 == T:
                    self.codon_colours[codon] = "violet"
                else:
                    pass
            if aa not in self.degen_groups:
                # Create an entry for this amino-acid.
                # Each codon position has its own four sets of synonymous codons.
                # The codons are grouped on the basis of the nucleotide they have
                # at the considered position.
                self.degen_groups[aa] = {
                    1 : {A : set([]), C : set([]), G : set([]), T : set([])},
                    2 : {A : set([]), C : set([]), G : set([]), T : set([])},
                    3 : {A : set([]), C : set([]), G : set([]), T : set([])}}
            self.degen_groups[aa][1][n1].add(codon)
            self.degen_groups[aa][2][n2].add(codon)
            self.degen_groups[aa][3][n3].add(codon)
        # Examples with the standard code:
        # self.degen_groups["l"][1]["T"] should be:
        # set([Codon("TTA"), Codon("TTG")])
        #
        # self.degen_groups["s"][2]["G"] should be:
        # set([Codon("AGT"), Codon("AGC")])
        #
        # self.degen_groups["s"][3]["C"] should be:
        # set([Codon("TCC"), Codon("AGC")])
        #
        # Now create the degenerate codons for the synonymy classes.
        #self.degen_by_aa = {"-":{1:[Codon("---")], 2:[Codon("---")], 3:[Codon("---")]}}
        self.degen_by_aa = {}
        # Why did I exclude the stop codons?
        #for aa in set(self.code.values()) - set(['-', '*']):
        for aa in set(self.code.values()) - set(['-']):
        #for aa in set(self.code.values()):
            self.degen_by_aa[aa] = {}
            for pos in [1, 2, 3]:
                self.degen_by_aa[aa][pos] = []
                for nuc in [A, C, G, T]:
                    # The set of codons for amino-acid aa that have nucleotide nuc at position pos.
                    codons = self.degen_groups[aa][pos][nuc]
                    if len(codons) != 0:
                        codon = reduce_by_or(codons)
                        if pos == 1 and aa == 'l':
                            if nuc == C:
                                self.codon_colours[codon] = "blue"
                            elif nuc == T:
                                self.codon_colours[codon] = "red"
                            else:
                                pass
                        if pos == 1 and aa == 'r':
                            if nuc == A:
                                self.codon_colours[codon] = "yellow"
                            elif nuc == C:
                                self.codon_colours[codon] = "green"
                            else:
                                pass
                        if pos == 1 and aa == 's':
                            if nuc == A:
                                self.codon_colours[codon] = "orange"
                            elif nuc == T:
                                self.codon_colours[codon] = "violet"
                            else:
                                pass
                        # self.degen_by_aa[aa][pos] will be the list of degenerate codons
                        # representing the codons for amino-acid aa
                        # having the same nucleotide at position pos.
                        self.degen_by_aa[aa][pos].append(reduce_by_or(codons))
        # Examples with the standard code:
        # self.degen_by_aa["l"][1] should be:
        # [Codon("CTN"), Codon("TTR")]
        #
        # self.degen_by_aa["s"][2] should be:
        # set([Codon("TCN"), Codon("AGY")])
        #
        # self.degen_by_aa["s"][3] should be:
        # set([Codon("TCA"), Codon("WSC"), Codon("TCG"), Codon("WST")])
        #
        # dictionary giving, for each position, a dictionary giving the degenerate codon
        # representing a given codon and that have the same nucleotide at the position.
        # The keys are positions, the values start as (almost) empty dictionaries
        # and will be filled if necessary when self.give_degen is used.
        self.cod2degen = {1 : {Codon("---") : Codon("---")},
                          2 : {Codon("---") : Codon("---")},
                          3 : {Codon("---") : Codon("---")}}

    def give_degen(self, codon, pos=1):
        """This method returns the degenerate codon representing the group of codons
        synonymous to *codon* and having the same nucleotide at position *pos* as *codon*."""
        if codon not in self.cod2degen[pos]:
            if len(codon.aas) == 1:
                aa = list(codon.aas)[0]
            else:
                msg = CAT(["%s is already a degenerate codon." % codon,
                           "More thinking will be necessary in order to decide ",
                           "how to deal with such a case.\n"])
                raise NotImplementedError, msg
            for degenerate in list(self.degen_by_aa[aa][pos]):
                if codon in degenerate:
                    self.cod2degen[pos][codon] = degenerate
                    break
        try:
            return self.cod2degen[pos][codon]
        except KeyError:
            # The loop on the list of degenerate codons
            # finished without a suitable codon to be found.
            msg = CAT(["The degenerate codon for %s " % codon,
                       "having the same nucleotide at position %d " % pos,
                       "cannot be found.\nIt could be that %s " % codon,
                       "is already degenerate and spans several degeneracy classes.\n"])
            raise NotImplementedError, msg
    
    def colour(self, codon):
        """This method returns the colour to be associated to the codon *codon*.
        If the codon is not yet recorded in self.codon_colours, its colour is inferred from
        the colours of the degenerate codons already recorded."""
        if codon not in self.codon_colours:
            for cod in self.codon_colours.keys():
                if codon in cod:
                    self.codon_colours[codon] = self.codon_colours[cod]
                    break
        return self.codon_colours[codon]

    def calculate_distance(self, cod1, cod2):
        """This method calculates the distance between codons that can be degenerated."""
        # TODO
        raise NotImplementedError, "This method has not been implemented yet."

standard_mutation_graph = MutationGraph()


def codon_from_triplet_slice(trps, n, code="Standard"):
    """returns a Codon object corresponding to the triplet at position *n*
    in the triplet slice *trps*. A triplet slice is as returned by
    Alignment.triplet_slice, defined in alignment_recoding.py"""
    return Codon(trps[0][n] + trps[1][n] + trps[2][n], code)


def codons_from_triplet_slice(trps, code="Standard"):
    """returns the list of the codons corresponding to the triplet slice *trps*.
    A triplet slice is as returned by Alignment.triplet_slice, defined in
    alignment_recoding.py"""
    return [Codon(codon, code) for codon in ["%s%s%s" % nnn for nnn in zip(
        trps[0], trps[1], trps[2])]]

# generalization of codons_from_triplet_slice, will not work while Codon has not been generalized.
#def codons_from_nuplet_slice(nps, code="Standard", n=3):
#    """returns the list of the codons corresponding to the n-uplet slice *nps*.
#    A n-uplet slice is as returned by Alignment.nuplet_slice, defined in
#    alignment_recoding.py"""
#    return [Codon(codon, code) for codon in map(lambda c: "%s" * n % c, zip(*nps))]

def codon_slice_has_aas(cod_slice, aas):
    """returns True if the codon slice *cod_slice* contains at least one codon
    that codes an amino-acid listed in *aas*. *cod_slice* should be a list of Codon
    objects."""

    #aas_at_site = reduce(lambda s1, s2 : s1 | s2, [codon.aas for codon in set(cod_slice)])
    #aas_at_site = reduce(operator.or_, [codon.aas for codon in set(cod_slice)])
    aas_at_site = reduce_by_or([codon.aas for codon in set(cod_slice)])
    return len(aas_at_site & set([aa.lower() for aa in aas])) != 0

def codon_slice_is_constant_aa(cod_slice, restrict_to=[]):
    """returns True if the codon slice *cod_slice* contains only codons that
    code the same amino-acid. If *restrict_to* is not empty, only those
    sites that are constant and code for one of the amino-acids the
    one-letter code of which is in *restrict_to* will be considered constant.
    *cod_slice* should be a list of Codon objects."""
    # set of codons observed in the slice
    codons = set(cod_slice)
    if len(codons - set([Codon("---")])) == 0:
        # There are only indels here.
        warnings.warn("A slice of the matrix was found for which there were only indels. "
                      "It will be considered constant unless some amino-acids "
                      "were specified with the restrict_to option.\n")
        if restrict_to == set(["-"]) or len(restrict_to) == 0:
            return True
        else:
            return False
    #aas_at_site = reduce(lambda s1, s2 : s1 | s2, [codon.aas for codon in codons])
    #aas_at_site = reduce(operator.or_, [codon.aas for codon in codons])
    aas_at_site = reduce_by_or([codon.aas for codon in codons])
    if len(aas_at_site - set(["-"])) > 1: # ignore indels
        return False
    elif len(aas_at_site) == 1:
        #print aas_at_site
        if len(restrict_to):
            if aas_at_site <= set([aa.lower() for aa in restrict_to]):
                return True
            else:
                # does not code a correct amino-acid
                return False
        else:
            return True
    else:
        return False

#TODO: document ignore and implement in more functions
def codon_slice_is_degenerate(cod_slice, restrict_to=[], ignore=[]):
    """returns True if the codon slice *cod_slice* contains two different codons
    that code the same amino-acid. If *restrict_to* is not empty, only those
    amino-acid the one-lettre code of which is in *restrict_to* are considered.
    *cod_slice* should be a list of Codon objects."""
    ignored_aas = set([aa.lower() for aa in ignore])
    if len(restrict_to):
        # Convert to lower case and remove duplicates.
        restricted_aas = set([aa.lower() for aa in restrict_to])
        # Filter the codon slice by checking that the codons
        # code at least one of the selected amino-acids.
        codons = set([cod for cod in cod_slice if len(cod.aas & restricted_aas)])
    else:
        codons = set(cod_slice)
    # Is a set an iterable ?
    #codons = sorted(codons)
    # It seems so: no need to convert to a list.
    for (cod1, cod2) in combinations(codons, 2):
        if len((cod1.aas & cod2.aas) - ignored_aas):
            # There are common amino-acids coded by the two different codons.
            return True
    # If we still haven't returned at this point, this means that no degeneracy has been found.
    return False

#TODO: check that the behaviour is correct when there are degenerate codons in the matrix.
def codon_position_is_degenerate(cod_slice, pos, restrict_to=[], ignore=[]):
    """returns True if the codon slice *cod_slice* contains two codons that
    code the same amino-acid and differ at position *pos*.
    *pos* should be 1, 2 or 3.
    If *restrict_to* is not empty, only those amino-acids
    the one-lettre code of which is in *restrict_to* are considered.
    *cod_slice* should be a list of Codon objects."""
    ignored_aas = set([aa.lower() for aa in ignore])
    if len(restrict_to):
        # Convert to lower case and remove duplicates.
        restricted_aas = set([aa.lower() for aa in restrict_to])
        # Filter the codon slice by checking that the codons
        # code at least one of the selected amino-acids.
        codons = set([cod for cod in cod_slice if len(cod.aas & restricted_aas)])
    else:
        codons = set(cod_slice)
    # Is a set an iterable ?
    #codons = sorted(codons)
    # It seems so: no need to convert to a list.
    for (cod1, cod2) in combinations(codons, 2):
        #inter_aas = cod1.aas & cod2.aas
        if len((cod1.aas & cod2.aas) - ignored_aas):
            # There are common amino-acids coded by the two different codons.
            if cod1[pos] & cod2[pos] == 0:
                # The codons do not share a nucleotide at the considered position.
                return True
                # But what happens if the amino-acids present in restricted_aas
                # for cod1 and for cod2 have no intersection with cod1 & cod2?
    # If we still haven't returned at this point, this means that no degeneracy has been found.
    return False

# Note that this does not discriminate between types of Serine codons:
# Ser1 and Ser2 can be degenerated together.
def degenerate_codon_slice(cod_slice, restrict_to=[]):
    """
    returns a recoding of the codon slice *cod_slice* where codons coding
    for the same amino-acid are replaced by their union (i.e. are degenerated).
    If *restrict_to* is not empty, only those codons that code amino-acids
    listed in *restrict_to* are degenerated.
    """
    msg = CAT(["The programmer was too lazy to make sure the code makes sense ",
               "with already degenerate codons, so such codons are not allowed."])
    assert not any([cod.degenerate for cod in cod_slice]), msg
    # Ensure the amino-acids are in lowercase.
    restrict_to = set([aa.lower() for aa in restrict_to])
    n_seq = len(cod_slice)
    # Not sure the copy is necessary.
    new_slice = [copy.copy(cod) for cod in cod_slice]
    # Loop over the slice.
    for i in range(n_seq):
        cod1 = cod_slice[i]
        # Only consider codons that are compatible with the restriction rule, if there is one.
        if (len(restrict_to) == 0) or (len(restrict_to & cod1.aas) != 0):
            # Compare with the codons further in the slice.
            j = i+1
            while j < n_seq:
                cod2 = cod_slice[j]
                # Error discovered the 04/06/2012
                #common_aas = cod1.aas & cod1.aas
                # What does it affect? Nothing, because degenerate_codon_slice is used
                # in degenerateByCodonColumn in alignment_recoding.py,
                # and this is not used in recode_matrix.py
                common_aas = cod1.aas & cod2.aas
                # Only consider codons that are compatible
                # with the restriction rule, if there is one.
                if (len(restrict_to) == 0 and len(common_aas) != 0) \
                   or (len(restrict_to & common_aas) != 0):
                    # There is at least one amino-acid coded in common.
                    new_slice[i] = new_slice[i] | cod2
                    new_slice[j] = new_slice[j] | cod1
                j += 1
    return new_slice

def recode_sequence(sequence, converter, positions=None, code="Standard"):
    """uses the correspondence rules provided by the dictionary *converter*
    to produce a recoded version of *sequence*, and returns it.
    *positions* determines which codon positions are recoded.
    By default, all positions are recoded.
    """
    gm = ['p4.code_utils.recode_sequence()']
    if isinstance(code, str):
        code = getBiopythonCode(code)
    else:
        msg = "code must be a dictionary, or a string naming the code in Biopython."
        assert isinstance(code, dict), msg
    # To get the size of the motifs being substituted, we look at the first one in the dictionary.
    subst_size = len(converter.keys()[0])
    if len(sequence) % subst_size != 0:
        gm.append("The length of the sequence should be a multiple of %i" % subst_size)
        raise P4Error(gm)
    if positions is not None:
        # Filter the converter.
        for codon in converter.keys():
            convert = converter[codon]
            # Replace the positions to be recoded by the converted codon, but keep the others.
            converter[codon] = CAT(
                [(convert[i-1] if i in positions else codon[i-1]) for i in range(
                    1, subst_size+1)])
    # Build the recoded version of the sequence.
    new_seq = ""
    # Loop over the codons (triplets, if subst_size == 3).
    for i in range(len(sequence) / subst_size):
        try:
            # Make a Codon instance (to convert it afterwards).
            codon = Codon(sequence[(subst_size * i):(subst_size * (i+1))], code)
        except CodonTranslationError, e:
            sys.stderr.write(
                "%s\nProblem at sequence slice %i:%i\n" % (
                    e, subst_size * i, subst_size * (i+1)))
            warnings.warn("We will replace the codon by indels.\n")
            try:
                codon = Codon("-" * subst_size, code)
            except CodonTranslationError, e:
                sys.stderr.write("We still don't know how to translate the codon. "
                                 "Bad implementation?\n")
                sys.exit(1)
        # Convert the codon.
        # If the converter has no entry for the codon, we don't convert it,
        # hence the converter.get() syntax, using a default value.
        if codon.degenerate:
            # The codon is decomposed into non-degenerate codons.
            # These codons are converted, and the resulting conversions
            # are "recomposed" into a new codon.
            # Can it be done more efficiently ?
            new_seq += str(reduce_by_or(
                [Codon(converter.get(
                    motif, motif), code) for motif in [str(cod) for cod in codon.decomposition()]]))
        else:
            #motif = str(codon)
            #new_seq += str(Codon(converter.get(motif, motif)))
            new_seq += str(Codon(converter.get(str(codon), str(codon)), code))
    return new_seq

def codon_usage(sequence, code="Standard"):
    """returns a dictionary where the keys are amino-acids represented by
    one-letter codes, and the values are dictionaries where keys are codons and
    values are number of occurrences found in the sequence *sequence*. *code*
    should be the Biopython name of the genetic code under which the sequence is to
    be interpreted. Alternatively, it can be the dictionary giving the
    correspondence between codons and amino-acids."""
    if isinstance(code, str):
        code = getBiopythonCode(code)
    else:
        msg = "code must be a dictionary, or a string naming the code in Biopython."
        assert isinstance(code, dict), msg
    # We assume that the "codons" have all the same length,
    # and we look at the first codon in the dictionary to know this length.
    codelength = len(code.keys()[0])
    #assert codelength == 3, "Amino-acids should be coded by triplets of nucleotides."
    # Number of codons.
    n_codons = len(sequence) / codelength
    if n_codons * codelength !=  len(sequence):
        msg = CAT(["The sequence does not contain an integral number of codons.\n",
                   "It will be assumed that it starts at a first position, "
                   "and the extra nucleotides will not be taken into account.\n"])
        warnings.warn(msg)
    aa_stats = {}
    for i in range(n_codons):
        codon = Codon(sequence[i*codelength:(i+1)*codelength], code)
        # A degenerate codon can represent several amino-acids.
        # Does it really make sense to count the codon for several
        # amino-acids in the codon usage statistics ?
        for aa in codon.aas:
            if aa not in aa_stats:
                # Create the entry if it doesn't exist.
                aa_stats[aa] = {}
            if codon in aa_stats[aa]:
                aa_stats[aa][codon] += 1
            else:
                aa_stats[aa][codon] = 1
    return aa_stats

