from p4.sequence import Sequence
from p4.sequencelist import SequenceList
from p4.nexussets import NexusSets
from p4.p4exceptions import P4Error
import string
import copy
import os
import math
import p4.func
import re
import sys
import random
from p4.nexussets import CharSet
import subprocess
from p4.distancematrix import DistanceMatrix
from p4.var import var
from p4.part import Part
import numpy
import numpy.linalg
import statistics
from p4.geneticcode import GeneticCode
import p4.pf as pf

longMessage1 = """
 You may want to do the alignment method 

   checkForDuplicateSequences(removeDupes=True, makeDict=True)

 which keeps the first of the duplicate sequences, and removes the
 subsequent duplicates from the alignment.  If makeDict is set (it is
 set by default), it renames the first to p4Dupe1 (or p4Dupe2, ...),
 and saves a dictionary with all the names of the dupes in a file (by
 default 'p4DupeSeqRenameDict.py', but you can change that file name)
 to facilitate programmatic restoration of names to any output that
 you may get (generally a tree).

 To restore the taxa to your resulting tree you can use the Tree
 method restoreDupeTaxa().

 (You are probably seeing this message because the variable
 var.doCheckForDuplicateSequences is set.  If you do not want to
 automatically check for duplicate sequences when you read in an
 alignment, then turn that variable off.)

"""


class ExcludeDelete(object):

    """A class for Alignment exclude/include and delete/restore.
    
       - Exclude         for character sites.   Works.
       - Include         for character sites.   No Workee!
       - Delete/restore  for taxa.              No Workee!
    
    In Nexus-speak, when you exclude character sites, the remaining
    sites are not renumbered.  So if we have an alignment of one short
    sequence, as::
    
        123456789    1-based char numbers
        acgtacgta    the sequence
    
    and I exclude chars 4, 5, and 6, I get::
    
        123789    1-based char numbers
        acggta    the sequence
    
    Now what if I want to exclude sites 2 and 8 from the original as
    well?  It would not work if I re-numbered the characters; I need
    access to the original numbers.  This class provides that access.
    
    It involves making a copy of the original Alignment.  This is only
    done if needed-- ie it is not done routinely-- it is only done if
    a call to exclude or delete is made.  A mask (list) is used; that
    mask starts out as all 1's.  If a character is to be excluded, the
    mask is set to zero for that position.  So in the above example,
    we could do this::
    
        a=alignment()
        a.excludeCharSet('e4_5_6')
        a.excludeCharSet('e2_8')
        print(a.excludeDelete.mask)      # [1, 0, 1, 0, 0, 0, 1, 0, 1]
        print(a.sequences[0].sequence)   # agga
        print(a.excludeDelete.sequences[0].sequence)   # acgtacgta
    
    The Alignment interface stuff is only partly written.  We have::
    
        Alignment.excludeCharSet()   # works
    
    but we do not have an ``include()`` yet, and we do not have either
    ``deleteTaxaSet()`` or ``restore()`` yet.
    
    """

    def __init__(self, alignment):
        self.alignment = alignment
        self.length = alignment.length
        self.sequences = copy.deepcopy(alignment.sequences)
        #self.taxNames = copy.deepcopy(alignment.taxNames)
        self.mask = [1] * self.length
        self.excludedCharSets = []

    def dump(self):
        print("\n  ExcludeDelete dump.")
        print("        length = %i" % self.length)
        print("        mask = ", end=' ')
        upper = 30
        if self.length < upper:
            upper = self.length
        for i in range(upper):
            print("%1i" % self.mask[i], end=' ')
        if upper == self.length:
            print('')
        else:
            print(" ...")
        print("        %i excludedCharSets:" % len(self.excludedCharSets))
        for cs in self.excludedCharSets:
            print("            %s" % cs.name)

    def _resetMask(self):
        self.mask = [1] * self.length
        for cs in self.excludedCharSets:
            if not cs.mask:
                cs.setMask()
            for i in range(cs.aligNChar):
                if cs.mask[i] == '1':
                    self.mask[i] = 0

    def _resetSequences(self):
        for sNum in range(len(self.alignment.sequences)):
            s = []
            for i in range(len(self.mask)):
                if self.mask[i]:
                    s.append(self.sequences[sNum].sequence[i])
            self.alignment.sequences[sNum].sequence = ''.join(s)
        self.alignment.length = len(self.alignment.sequences[0].sequence)


class Alignment(SequenceList):

    """A SequenceList where all the sequences are the same length and dataType.

    From SequenceList, this inherits

    * ``sequences``, a list of :class:`~p4.sequence.Sequence` objects
    * ``fName``, the file name, by default None
    * ``sequenceForNameDict``, a dictionary which allows you to get a
      Sequence given its name.  This week it is not made by default
      --- you need to make it explicitly with the inherited method 
      :meth:`~p4.sequencelist.SequenceList.makeSequenceForNameDict`.
    * and the other methods from  :class:`~p4.sequencelist.SequenceList`

    Alignment objects additionally have the ``length``, which has a
    synonym ``nChar``, and is also accessible via ``len(self)``.

    Also:

    * ``nTax``, the number of sequences
    * ``taxNames``, a list of the names of the sequences
    * ``dataType``, a NEXUS-centric string --- 'dna', 'protein', or 'standard'
    * ``symbols``, the character symbols, not including equates or ambiguities
    * ``equates``, a dictionary of NEXUS-style ambiguities, eg ``'r':['a','g']`` for DNA
    * ``dim``, the dimension, ie the number of symbols
    * ``nexusSets``, a :class:`~p4.nexussets.NexusSets` object, usually
      made from a copy of ``var.nexusSets``, but made specific to self

    **Various checks are made when alignments are read in from files**


     .. currentmodule:: p4.alignment

     .. autosummary::
        :nosignatures:

        p4.sequencelist.SequenceList.checkNamesForDupes
        ~Alignment.checkForAllGapColumns
        ~Alignment.checkForBlankSequences
        ~Alignment.checkForDuplicateSequences
        ~Alignment.checkLengthsAndTypes

    The checks above are under the control of some variables in
    :py:class:`~p4.var.Var`, always available to you as the instance ``var``.

    * ``var.doCheckForDuplicateSequenceNames``
    * ``var.doRepairDupedTaxonNames``
    * ``var.doCheckForAllGapColumns``
    * ``var.doCheckForBlankSequences``
    * ``var.doCheckForDuplicateSequences``

    **Writing**
     .. autosummary::
        :nosignatures:

        ~Alignment.writeNexus
        ~Alignment.writePhylip
        p4.sequencelist.SequenceList.writeFasta
        ~Alignment.writeMolphy

    **Copy self**
     .. autosummary::
        :nosignatures:

        ~Alignment.dupe
        ~Alignment.noGapsOrAmbiguitiesCopy
        ~Alignment.bootstrap
        p4.data.Data.bootstrap
        ~Alignment.sequenceSlice

    **Masks are strings that are nChar long, usually just 0s and 1s**
     .. autosummary::
        :nosignatures:

        ~Alignment.constantMask
        ~Alignment.simpleConstantMask
        ~Alignment.gappedMask
        ~Alignment.getAllGapsMask
        ~Alignment.getEnoughCharsMask
        ~Alignment.getLikelihoodTopologyInformativeSitesMask
        ~Alignment.getMaskForAutapomorphies
        ~Alignment.getMaskForCharDiversity


    You can also make masks with :func:`p4.func.maskFromNexusCharacterList`.
    You can combine masks bitwise with the Alignment methods :meth:`~Alignment.andMasks` and :meth:`~Alignment.orMasks`

    **CharSets and CharPartitions**

    A NEXUS-style charSet or charPartition is usually defined using a
    NEXUS sets block, as explained in :class:`p4.nexussets.NexusSets`.  A
    NEXUS sets block might be something like::

        begin sets;
            charset pos1 = 1 - .\3;
            charset pos2 = 2 - .\3;
            charset pos3 = 3 - .\3;
        end;

    or ::

        begin sets;
            charset gene1 = 1 - 103;
            charset gene2 = 104 - 397;
            charpartition by_gene = gene1:gene1, gene2:gene2;
        end;

    However, you can also make new charSets from an
    Alignment.nexusSets with a mask using :meth:`p4.nexussets.NexusSets.newCharSet`

    CharPartitions must cover the entire length of self, without
    overlap of the character partition subsets.

    Whenever you do something like
    :meth:`Alignment.subsetUsingCharSet` that requires
    ``self.nexusSets``, if ``self.nexusSets`` does not exist yet then
    :meth:`Alignment.setNexusSets` is automatically called to make it
    --- so usually you do not need to do
    :meth:`Alignment.setNexusSets`.

     .. autosummary::
        :nosignatures:

        Alignment.excludeCharSet
        Alignment.setCharPartition
        Alignment.setGBlocksCharSet
        Alignment.setNexusSets
        p4.nexussets.NexusSets.newCharSet
        p4.nexussets.NexusSets.dupeCharSet

    **Extracting subsets**
     .. autosummary::
        :nosignatures:

        Alignment.subsetUsingCharSet
        Alignment.subsetUsingMask

    **Recoding the data into groups**
     .. autosummary::
        :nosignatures:

        Alignment.recodeDayhoff
        Alignment.recodeProteinIntoGroups
        Alignment.recodeRY

    **Translating DNA to protein**
     .. autosummary::
        :nosignatures:

        Alignment.translate
        Alignment.checkTranslation

    """

    def __init__(self):

        SequenceList.__init__(self)
        # Inherited from SequenceList:
        #self.sequences = []
        #self.fName = None
        #self.sequenceForNameDict = None

        self.length = 0

        #: Nexus-centric data type.  One of dna, rna, protein, or standard.
        self.dataType = None

        #self.sequences = []
        #self.fName = None
        #self.sequenceForNameDict = None

        #: Lowercase string, eg 'acgt'.  The order is all-important.
        self.symbols = None

        #: The number of symbols, eg 4 for DNA.
        self.dim = None

        #: A dictionary of NEXUS-style equates, eg r=[a,g] in DNA
        self.equates = {}   # A hash

        #: A :class:`~p4.nexussets.NexusSets` object, perhaps copied from
        #: var.nexusSets and made specific to self.  You can do a
        #: :meth:`p4.nexussets.NexusSets.dump` on it to see what is in
        #: there.
        self.nexusSets = None

        #: A list of Part objects, encapsulating data partitions in
        #: :class:`Data` objects.  There would be one or more parts in
        #: an Alignment.
        self.parts = []
        self.excludeDelete = None  # An ExcludeDelete object

    @property
    def nTax(self):
        """Return the number of sequences"""
        return len(self.sequences)

    @property
    def nChar(self):
        """Return the length of the alignment"""
        return self.length

    @property
    def nEquates(self):
        """Return the number of equates"""
        return len(self.equates)

    def _getTaxNames(self):
        theTaxNames = []
        for s in self.sequences:
            theTaxNames.append(s.name)
        return theTaxNames

    def _setTaxNames(self, theArg=None):
        gm = ["Alignment._setTaxNames()"]
        gm.append("Attempt to set Alignment taxNames.")
        gm.append("However, it is a property, so don't do that.")
        raise P4Error(gm)

    def _delTaxNames(self):
        gm = ["Alignment._delTaxNames()"]
        gm.append("Caught an attempt to delete self.taxNames, but")
        gm.append("self.taxNames is a property, so you shouldn't delete it.")
        raise P4Error(gm)

    taxNames = property(_getTaxNames, _setTaxNames, _delTaxNames)
    """A list of the names of the Sequences. """

    def _getRows(self):
        return self.sequences

    def _getColumns(self):
        columns = []
        for pos in range(self.nChar):
            column = []
            for seqNum in range(0, self.nTax):
                column.append(self.sequences[seqNum].sequence[pos])
            columns.append(column)
        return columns

    # Here I over-ride __bool__().  If the self.length len is zero, and I don't
    # have __bool__() redefined as below, then "assert self" will raise an
    # AssertionError, basing that response on the result of len(self).  Having
    # __bool__() redefined here makes "assert self" work, even with no sequence
    # length.  Previously, python 2 only, I had to over-ride __nonzero__() for
    # the same reason --- but that does not work with Python 3. 
    # Checked July 2020, this is still needed.  See similar in Sequence class. 

    def __bool__(self):
        return True

    def __len__(self):
        return self.length

    def checkLengthsAndTypes(self):
        """Last checks after reading an Alignment.

        Make sure the sequence lengths and dataType are all the same.
        Set self.length, self.dataType, self.symbols, self.dim, and self.equates
        """

        gm = ["Alignment.checkLengthsAndTypes()"]
        if self.fName:
            gm.append("fName = %s" % self.fName)
        if len(self.sequences) > 1:
            len0 = len(self.sequences[0].sequence)
            type0 = self.sequences[0].dataType
            for i in range(len(self.sequences))[1:]:
                if len0 != len(self.sequences[i].sequence):
                    gm.append("(Zero-based) sequence %i (%s) length (%i)," %
                              (i, self.sequences[i].name, len(self.sequences[i].sequence)))
                    gm.append(
                        "is not the same as the first sequence length (%i)." % len0)
                    raise P4Error(gm)
                if type0 != self.sequences[i].dataType:
                    gm.append("Type of (zero-based) sequence %i (%s)," %
                              (i, self.sequences[i].dataType))
                    gm.append(
                        "is not the same as the first sequence type (%s)." % type0)
                    raise P4Error(gm)
            self.length = len0
            self.dataType = type0
        elif len(self.sequences) == 1:
            self.length = len(self.sequences[0].sequence)
            self.dataType = self.sequences[0].dataType
        elif len(self.sequences) == 0:
            print(gm[0])
            print("    The alignment has no sequences!")
        if self.dataType == 'dna':
            self.symbols = 'acgt'
            self.dim = 4
            if not self.equates:
                self.equates = {'n': 'acgt', 'm': 'ac', 'k': 'gt',  # 'x': 'acgt',
                                'h': 'act', 'y': 'ct', 'v': 'acg',
                                'w': 'at', 'd': 'agt', 'b': 'cgt',
                                'r': 'ag', 's': 'cg'}
        elif self.dataType == 'protein':
            self.symbols = 'arndcqeghilkmfpstwyv'
            self.dim = 20
            if not self.equates:
                self.equates = {
                    'b': 'dn', 'x': 'arndcqeghilkmfpstwyv', 'z': 'eq'}

        elif self.dataType == 'standard':
            if not self.symbols:
                gm.append("symbols are missing.")
                raise P4Error(gm)
            if not self.dim:
                self.dim = len(self.symbols)
            if not self.equates:
                self.equates = {}
        elif self.dataType == 'rna':
            self.symbols = 'acgu'
            self.dim = 4
            if not self.equates:
                self.equates = {'n': 'acgu', 'm': 'ac', 'k': 'gu',  # 'x': 'acgu',
                                'h': 'acu', 'y': 'cu', 'v': 'acg',
                                'w': 'au', 'd': 'agu', 'b': 'cgu',
                                'r': 'ag', 's': 'cg'}
        else:
            gm.append("unknown dataType %s." % self.dataType)
            raise P4Error(gm)

    def composition(self, sequenceNumberList=None):
        """Returns a list of compositions.

        This returns a list of floats, the composition of the
        sequence(s) in the sequenceNumberList.  If the
        sequenceNumberList=None (the default), then the overall
        composition is given, which is the mean of the individual
        sequences (the sequence comps are not weighted by the sequence
        length, ie the non-gap sequence length).  For DNA, the order
        is acgt.  For protein, the order is arndcqeghilkmfpstwyv.  For
        standard, its the order in symbols.  Gaps and questionmarks
        are ignored.  Equates are handled properly, iterating to the
        final comp.
        """

        gm = ['Alignment.composition(sequenceNumberList=%s).' % sequenceNumberList]
        dbug = 0

        # symbolFreq and equateFreq are hashes for the raw counts
        symbolFreq = {}
        for symb in self.symbols:
            symbolFreq[symb] = 0.0
        if self.equates:
            equateFreq = {}
            for equate in self.equates:
                equateFreq[equate] = 0.0
        hasEquates = 0

        if dbug:
            print("Alignment.composition() sequenceNumberList = %s" % sequenceNumberList)
        if sequenceNumberList:
            if not isinstance(sequenceNumberList, list):
                gm.append("The sequenceNumberList should be a list, ok?")
                raise P4Error(gm)
            if len(sequenceNumberList) == 0:
                gm.append("The sequenceNumberList should have something in it, ok?")
                raise P4Error(gm)
        else:
            sequenceNumberList = range(len(self.sequences))

        result = [0.0] * self.dim
        grandNSites = 0
        import math
        epsilon = 1.0e-12
        maxIterations = 1000

        for i in sequenceNumberList:
            if not isinstance(i, int):
                gm.append("The sequenceNumberList should be integers, ok?")
                raise P4Error(gm)
            if i < 0 or i > len(self.sequences) - 1:
                gm.append("Item '%i' in sequenceNumberList is out of range" % i)
                raise P4Error(gm)

            seq = self.sequences[i]

            nGapsMissings = 0
            for j in seq.sequence:
                if j == '-' or j == '?':
                    nGapsMissings += 1
            nSites = self.length - nGapsMissings
            grandNSites = grandNSites + nSites

            for symb in self.symbols:
                #symbolFreq[symb] = symbolFreq[symb] + float(seq.sequence.count(symb))
                symbolFreq[symb] = float(seq.sequence.count(symb))
            if self.equates:
                for equate in self.equates:
                    # equateFreq[equate] = equateFreq[equate] + float(seq.sequence.count(equate))
                    equateFreq[equate] = float(seq.sequence.count(equate))
                    if equateFreq[equate]:
                        hasEquates = 1

            if dbug:
                print("symbolFreq = ", symbolFreq)
                print("equateFreq = ", equateFreq)

            initComp = 1.0 / self.dim
            comp = {}
            symbSum = {}
            for symb in self.symbols:
                comp[symb] = initComp
                symbSum[symb] = 0.0

            for dummy in range(maxIterations):
                for symb in self.symbols:
                    symbSum[symb] = symbolFreq[symb]
                if hasEquates:
                    for equate in self.equates:
                        if equateFreq[equate]:
                            factor = 0.0
                            for symb in self.equates[equate]:
                                factor = factor + comp[symb]
                            for symb in self.equates[equate]:
                                symbSum[symb] = symbSum[
                                    symb] + (equateFreq[equate] * (comp[symb] / factor))
                factor = 0.0
                for symb in self.symbols:
                    factor = factor + symbSum[symb]
                if not factor:
                    gm.append('(Zero-based) sequence %i. Empty?' % i)
                    gm.append(
                        'Perhaps exclude it by specifying arg sequenceNumberList.')
                    gm.append('(Or use the Data/Part composition.)')
                    raise P4Error(gm)
                diff = 0.0
                for symb in self.symbols:
                    oldComp = comp[symb]
                    comp[symb] = symbSum[symb] / factor
                    diff = diff + math.fabs(comp[symb] - oldComp)
                if 0:
                    print("diff=%8.5f %3i  " % (diff, i), end=' ')
                    for symb in self.symbols:
                        print("%s: %.6f  " % (symb, comp[symb]), end=' ')
                    print('')
                if diff < epsilon:
                    # print("did %i iterations" % dummy)
                    break

            for j in range(len(self.symbols)):
                result[j] = result[j] + (comp[self.symbols[j]] * nSites)

        for j in range(len(self.symbols)):
            result[j] = result[j] / grandNSites
        return result


# def subsetUsingCharPartition(self, charPartitionName, inverse=0):
# """Return a subset of self based on a character partition.

# A charpartition has one or more subsets, which together need
# not span the length of self.  This method would of course only
# be useful if the charpartition does not span the entire
# alignment-- otherwise you would get the whole alignment.  This
# method makes a new mask (using the CharPartition.mask()
# method) which has a 1 wherever any subset of the charPartition
# has a 1, and a zero otherwise.  Returns an alignment.  """

##        gm = ['Alignment.subsetUsingCharPartition()']
# if not self.nexusSets:
# self.setNexusSets()

# if not len(self.nexusSets.charPartitions):
##            gm.append("This alignment has no charPartitions")
##            raise P4Error(gm)
##        theCP = None
##        lowName = charPartitionName.lower()
# for cp in self.nexusSets.charPartitions:
# if cp.lowName == lowName:
##                theCP = cp
# break
# if theCP == None:
##            gm.append("This alignment has no charPartition named '%s'" % charPartitionName)
##            raise P4Error(gm)

# Get the mask
##        m = theCP.mask(self.nexusSets, self)

# print("The mask is: %s" % m)
# if not inverse:
##            a = self.subsetUsingMask(m, theMaskChar='1')
# else:
##            a = self.subsetUsingMask(m, theMaskChar='0')
# return a

# def subsetUsingCharPartitionSubset(self, charPartitionName, charPartitionSubsetName, inverse=0):
# """Return a subset of self based on a charPartition subset.

# A charPartition has one or more subsets, each of which can
# have a mask.  This method uses that mask to subset self.
# Returns an alignment.  """

##        gm = 'Alignment.subsetUsingCharPartitionSubset(charPartitionName=\'%s\'' % charPartitionName
##        gm += ', charPartitionSubsetName=\'%s\', inverse=%s)' % (charPartitionSubsetName, inverse)
##        gm = [gm]

# if not self.nexusSets:
# self.setNexusSets()

# if not len(self.nexusSets.charPartitions):
##            gm.append("This alignment has no charPartitions")
##            raise P4Error(gm)

# Find the charPartition
##        theCP = None
##        lowName = charPartitionName.lower()
# for cp in self.nexusSets.charPartitions:
# if cp.lowName == lowName:
##                theCP = cp
# break
# if theCP == None:
##            gm.append("This alignment has no charPartition named '%s'" % charPartitionName)
##            raise P4Error(gm)

# Find the charPartitionSubset
##        theCPsubset = None
##        lowName = charPartitionSubsetName.lower()
# for cps in theCP.subsets:
# if cps.lowName == lowName:
##                theCPsubset = cps
# break
# if theCPsubset == None:
# gm.append("The charPartition '%s' has no charPartitionSubset named '%s'" % \
# (charPartitionName, charPartitionSubsetName))
##            raise P4Error(gm)

# Prepare the mask
##        theCP.setSubsetMasks(self.nexusSets, self)

# if not inverse:
##            a = self.subsetUsingMask(theCPsubset.mask, theMaskChar='1')
# else:
##            a = self.subsetUsingMask(theCPsubset.mask, theMaskChar='0')
# return a

    def subsetUsingCharSet(self, charSetName, inverse=0):
        """Return a subset of self based on a charSet.

        A charset has a mask, composed of zeros and ones, which is
        used to subset self.  Returns an alignment.

        For example::

          read("myAlignment.nex")
          a = var.alignments[0]
          read('myNexusSets.nex') # with a charset named 'foo'
          b = a.subsetUsingCharSet('foo')
          b.writeNexus('myFooSubset.nex')

        """

        gm = ['Alignment.subsetUsingCharSet(charSetName=\'%s\', inverse=%s)' % (
            charSetName, inverse)]

        if not self.nexusSets:
            self.setNexusSets()

        theCS = None
        lowName = charSetName.lower()
        if lowName not in self.nexusSets.predefinedCharSetLowNames and lowName not in self.nexusSets.charSetLowNames:
            gm.append("This alignment has no charset named '%s'" % charSetName)
            raise P4Error(gm)
        if lowName in self.nexusSets.predefinedCharSetLowNames:
            if lowName == 'constant':
                theCS = self.nexusSets.constant
            elif lowName == 'gapped':
                theCS = self.nexusSets.gapped
        else:
            for cs in self.nexusSets.charSets:
                if cs.lowName == lowName:
                    theCS = cs
                    break
        if theCS == None:
            gm.append(
                "This should not happen -- alignment has no charset named '%s'" % charSetName)
            raise P4Error(gm)
        assert theCS.aligNChar
        assert theCS.mask

        # prepare the mask
        # if not theCS.mask:
        #    theCS.setMask()
        # print("The mask is: %s" % theCS.mask)
        if len(theCS.mask) != self.length:
            gm.append("The length of the mask is %i, the length of the alignment is %i" % (
                len(theCS.mask), self.length))
            raise P4Error(gm)

        if not inverse:
            a = self.subsetUsingMask(theCS.mask, theMaskChar='1')
        else:
            a = self.subsetUsingMask(theCS.mask, theMaskChar='0')
        return a

    def subsetUsingMask(self, theMask, theMaskChar='1', inverse=0):
        """Returns a subset of self based on a mask.

        This makes a copy of the alignment based on theMask.  Arg
        theMask is a string, the same length as the alignment.  The
        character theMaskChar determines which positions are included
        in the copy.  If the character in theMask is theMaskChar, that
        position is included in the copy.  Otherwise, no.  If inverse
        is set, then theMaskChar determines which positions are
        excluded, and all other positions are included.

        """

        gm = ['Alignment.subsetUsingMask()']

        # Make sure theMask is the right length, depending on whether
        # characters have been excluded or not.
        if self.excludeDelete:
            if len(theMask) != self.excludeDelete.length:
                gm.append("The mask length (%i) does not" % len(theMask))
                gm.append(
                    "equal the (pre-exclude charSets) alignment length (%i)" % self.excludeDelete.length)
                raise P4Error(gm)
        else:
            if len(theMask) != self.length:
                gm.append("The mask length (%i) does not" % len(theMask))
                gm.append("equal the alignment length (%i)" % self.length)
                raise P4Error(gm)

        if not isinstance(theMaskChar, str) or len(theMaskChar) != 1:
            gm.append("theMaskChar needs to be a single-character string")
            raise P4Error(gm)

        # first, make a copy
        a = copy.deepcopy(self)
        a.excludeDelete = None
        if len(a.parts):
            for i in a.parts:
                i.cPart = None
        a.parts = []
        a.nexusSets = None

        theMask2 = [0] * len(theMask)
        if not inverse:
            for i in range(len(theMask2)):
                if theMask[i] == theMaskChar:
                    theMask2[i] = 1
        else:
            for i in range(len(theMask2)):
                if theMask[i] != theMaskChar:
                    theMask2[i] = 1

        if self.excludeDelete:
            for i in range(self.excludeDelete.length):
                if self.excludeDelete.mask[i] == 0:
                    theMask2[i] = 0

        # From now on, inverse or not does not apply to theMask2.
        # Nor does theMatchChar apply-- its just 1's and zeros.
        # If it is 0, exclude it.  If it is 1, include it.

        a.length = theMask2.count(1)
        if a.length == 0:
            if not var.allowEmptyCharSetsAndTaxSets:
                gm.append("The mask has a length of zero.")
                gm.append("(Allow by turning var.allowEmptyCharSetsAndTaxSets on.)")
                raise P4Error(gm)

        # make a 2D array the same size as the sequences, filled.
        newList = []
        for i in range(len(self.sequences)):
            one = ['a'] * a.length
            newList.append(one)
        # print("newList = ", newList)

        # fill the array with slices from the sequences
        k = 0
        # print("self.length = %i, a.length = %i" % (self.length, a.length))

        if self.excludeDelete:
            for i in range(len(theMask2)):
                if theMask2[i]:
                    for j in range(len(self.excludeDelete.sequences)):
                        newList[j][k] = self.excludeDelete.sequences[
                            j].sequence[i]
                    k = k + 1
        else:
            for i in range(len(theMask2)):
                if theMask2[i]:
                    for j in range(len(self.sequences)):
                        newList[j][k] = self.sequences[j].sequence[i]
                    k = k + 1

        # replace the sequences
        for i in range(len(self.sequences)):
            a.sequences[i].sequence = ''.join(newList[i])

        a.checkLengthsAndTypes()
        if 0:
            from p4.nexussets import NexusSets
            a.nexusSets = NexusSets()
            a.nexusSets.aligNChar = a.length
            a.nexusSets.setPredefinedCharSets(a)
        return a

    def checkForDuplicateSequences(self, removeDupes=False, makeDict=True, dictFileName='p4DupeSeqRenameDict.py', dupeBaseName='p4Dupe'):
        """Like it says, with verbose output.

        If the p4 variable var.doCheckForDuplicateSequences is set,
        this method, with removeDupes=False, is called automatically
        every time an alignment is read in.  If you would rather it
        not do that, turn var.doCheckForDuplicateSequences off.

        When these automatic checks are done, if any dupe sequences
        are found then a verbose warning is issued, inviting you to
        run it again with removeDupes turned on.

        If removeDupes is set, the duplicate sequences after the first
        are removed.  This option is not set by default.

        If there are duplicate sequences pairs, then the list of
        sequence dupe pairs is returned, as a list of 2-tuples, of the
        two sequence objects.

        If both removeDupes and makeDict are set, then it will rename
        the first sequence to p4Dupe1 (or p4Dupe2, and so on --- the
        dupeBaseName is 'p4Dupe' by default, but it can be set as an
        arg) and make a dictionary to hold the other names, and write
        that dictionary to a file (by default p4DupeSeqRenameDict.py).
        The option makeDict is set by default, but it won't happen
        unless removeDupes is also set, and there are dupes to be
        removed.
        """

        #gm = ['Alignment.checkForDuplicateSequences()']
        dupeNumPairs = []
        dupes = []
        firsts = {}
        doneDupeSeqNums = set([])

        if removeDupes and makeDict:
            if os.path.isfile(dictFileName):
                gm = ['Alignment.checkForDuplicateSequences()']
                gm.append("file '%s' already exists" % dictFileName)
                raise P4Error(gm)

        theRange = range(self.length)
        for i in range(len(self.sequences))[:-1]:
            if i not in doneDupeSeqNums:
                si = self.sequences[i].sequence
                for j in range(len(self.sequences))[i + 1:]:
                    if j not in doneDupeSeqNums:
                        sj = self.sequences[j].sequence
                        # print("trying", i, j)
                        isSame = True
                        for k in theRange:
                            if si[k] != sj[k]:
                                isSame = False
                                break

                        # xfasta, with RNA structure line.
                        if hasattr(self.sequences[i], 'parens'):
                            sip = self.sequences[i].parens
                            sjp = self.sequences[j].parens
                            p_isSame = True
                            for k in theRange:
                                if sip[k] != sjp[k]:
                                    p_isSame = False
                                    break
                            if isSame and not p_isSame:
                                print("(One-based) Sequences %i and %i are the same," % (i + 1, j + 1), end=' ')
                                print("but the structures differ.")
                            # else:
                            #    print('ok')

                        if isSame:
                            dupeNumPairs.append([i, j])
                            doneDupeSeqNums.add(i)
                            doneDupeSeqNums.add(j)

        if dupeNumPairs:
            # print(dupeNumPairs)
            if not removeDupes:
                sequencePairs = []
                print()
                print("=" * 50)
                if self.fName:
                    print(" Alignment from file '%s'" % self.fName)
                print(" This alignment has duplicate sequences!")
                print(" Sequence numbers below are 1-based.")
                for dp in dupeNumPairs:
                    i = dp[0]
                    j = dp[1]
                    print("    sequence %i (%s) is the same as sequence %i (%s)." % (
                        i + 1, self.sequences[i].name, j + 1, self.sequences[j].name))
                    sequencePairs.append((self.sequences[i], self.sequences[j]))
                print(longMessage1)
                print("=" * 50)
                print()
                return sequencePairs

            else:     # ie do removeDupes 
                if makeDict:
                    myDict = {}
                    newNameCounter = 1
                    dpNum = 0  # dupe pair index
                    while 1:
                        dp = None
                        try:
                            dp = dupeNumPairs[dpNum]
                        except IndexError:
                            break
                        i = dp[0]
                        j = dp[1]
                        iName = self.sequences[i].name
                        jName = self.sequences[j].name
                        newIName = "%s%i" % (dupeBaseName, newNameCounter)
                        myDict[newIName] = [iName, jName]
                        self.sequences[i].name = newIName
                        newNameCounter += 1

                        # Get the other j's for the same i.
                        while 1:
                            dpNum += 1
                            dp = None
                            try:
                                dp = dupeNumPairs[dpNum]
                            except IndexError:
                                break
                            if dp[0] == i:
                                myDict[newIName].append(self.sequences[dp[1]].name)
                            else:
                                break
                    f = open(dictFileName, 'w')
                    f.write("p4DupeSeqRenameDict = %s\n" % myDict)
                    f.close()

                # Remove the dupe sequences.
                toRemove = []
                sequencePairs = []
                for dp in dupeNumPairs:
                    i = dp[0]
                    j = dp[1]   # remove the second, not the first
                    toRemove.append(self.sequences[j])
                    sequencePairs.append((self.sequences[i], self.sequences[j]))
                    #self.sequences[dp[0]].name += "_%s" % self.sequences[j].name
                for s in toRemove:
                    self.sequences.remove(s)

                if self.nexusSets and self.nexusSets.taxSets:
                    print()
                    print("-" * 50)
                    print("There are tax sets, possibly affected by dupe removal.")
                    print("So I am removing those taxSets.")
                    print("-" * 50)
                    self.nexusSets.taxSets = []
                return sequencePairs

    def checkForBlankSequences(self, removeBlanks=False, includeN=True, listSeqNumsOfBlanks=False):
        """Like it says, with verbose output.

        If the p4 variable var.doCheckForBlankSequences is set,
        this method, with removeBlanks=False, is called automatically
        every time an alignment is read in.  If you would rather it
        not do that, turn var.doCheckForBlankSequences off.

        When these automatic checks are done, if any blank sequences
        are found then a verbose P4Error is raised, inviting you to
        run it again with removeBlanks turned on.

        Sometimes you just want to know what sequences are blank; get
        the sequence numbers by turning arg listSeqNumsOfBlanks on.
        That returns a list of sequence numbers, without removing
        blanks or raising a P4Error.

        If removeBlanks is set, the blank sequences are removed.  This
        option is not set by default.

        If this method removes sequences, it returns the number of
        blank sequences removed.

        Blank sequences are defined as sequences wholly composed of
        '?' and '-'.  If includeN is turned on (which it is by
        default) then 'n' is included for DNA, and 'x' for protein.
        If that is done, that means, for DNA, that it is blank if it
        is wholly composed of '?', '-', and 'n'.
        """

        gm = ['Alignment.checkForBlankSequences()']
        if self.fName:
            gm.append("fName = %s" % self.fName)
        bChars = ['-', '?']
        if includeN:
            if self.dataType == 'dna':
                bChars.append('n')
            elif self.dataType == 'protein':
                bChars.append('x')
        blankSeqs = []
        seqNums = []
        for seqNum in range(len(self.sequences)):
            seqObj = self.sequences[seqNum]
            isBlank = True
            for c in seqObj.sequence:
                if c not in bChars:
                    isBlank = False
                    break
            if isBlank:
                blankSeqs.append(seqObj)
                seqNums.append(seqNum)
                # print("=" * 50)
                # seqObj.write()
        if listSeqNumsOfBlanks:
            return seqNums

        if blankSeqs and not removeBlanks:
            gm.append("This alignment has %i blank sequences," %
                      len(blankSeqs))
            gm.append("wholly composed of %s." % bChars)
            gm.append("To remove them, re-run this method, with arg")
            gm.append("removeBlanks turned on.")
            gm.append(
                "To prevent checking, turn var.doCheckForBlankSequences off.")
            raise P4Error(gm)

        if blankSeqs and removeBlanks:
            for s in blankSeqs:
                self.sequences.remove(s)

            if self.nexusSets and self.nexusSets.taxSets:
                print()
                print("-" * 50)
                print("There are tax sets, possibly affected by blank sequence removal.")
                print("So I am removing them.")
                print("-" * 50)
                self.nexusSets.taxSets = []
            return len(blankSeqs)
        return 0

    def checkForAllGapColumns(self, returnMask=False):
        """Check for alignment columns that are made of only gap or ? chars.

        By default, p4 does this on new alignments that are read.  It
        is under the control of var.doCheckForAllGapColumns.

        If there are all gap columns then a verbose output is printed
        and a P4Error is raised, unless ``returnMask`` is set, in which
        case no output is printed, no P4Error is raised, but the mask
        is returned.
        """
        gm = ['Alignment.checkForAllGapColumns()']
        if self.fName:
            gm.append("fName = %s" % self.fName)
        allGapPositions = []
        firstSeqSequence = self.sequences[0].sequence
        if '-' not in firstSeqSequence and '?' not in firstSeqSequence:
            return
        for pos in range(self.nChar):
            if firstSeqSequence[pos] == '-' or firstSeqSequence[pos] == '?':
                allGaps = True
                for seqNum in range(1, self.nTax):
                    c = self.sequences[seqNum].sequence[pos]
                    if c == '-' or c == '?':
                        pass
                    else:
                        allGaps = False
                        break
                if allGaps:
                    allGapPositions.append(pos)
        if allGapPositions:
            if returnMask:
                m = ['0'] * self.nChar
                for i in allGapPositions:
                    m[i] = '1'
                return ''.join(m)
            else:
                gm.append("The following %i positions were composed solely of '-' or '?'" % len(
                    allGapPositions))
                gm.append("Zero-based numbering - %s" % allGapPositions)
                for i in range(len(allGapPositions)):
                    allGapPositions[i] = allGapPositions[i] + 1
                gm.append("One-based numbering - %s" % allGapPositions)
                nxsString = ' '.join(["%i" % i for i in allGapPositions])
                gm.append('nexus charSet %s' % nxsString)
                gm.append("(To turn off auto-checking for all-gap columns,")
                gm.append("turn var.doCheckForAllGapColumns off.)")
                raise P4Error(gm)

    def dump(self):
        """Print rubbish about self."""

        print("\nAlignment dump:")
        if self.fName:
            print("  File name '%s'" % self.fName)
        if self.length:
            print("  Length is %s" % self.length)
        # if hasattr(self, 'nTax'):
        #    print("    nTax is %i" % self.nTax)
        # if hasattr(self, 'nChar'):
        #    print("    nChar is %i" % self.nChar)
        if hasattr(self, 'dataType'):
            print("  dataType is '%s'" % self.dataType)
        if hasattr(self, 'symbols'):
            print("  symbols are '%s'" % self.symbols)
        if self.equates:
            print("  equates")
            theKeys = list(self.equates.keys())
            theKeys.sort()
            for k in theKeys:
                print("%20s  %-30s" % (k, self.equates[k]))
        if self.nexusSets:
            if self.nexusSets.charSets:
                if len(self.nexusSets.charSets) == 1:
                    print("  There is 1 charSet")
                else:
                    print("  There are %i charSets" % len(self.nexusSets.charSets))
                for cp in self.nexusSets.charSets:
                    print("          %s" % cp.name)
            if self.nexusSets.charPartitions:
                if len(self.nexusSets.charPartitions) == 1:
                    print("  There is 1 charPartition")
                else:
                    print("  There are %i charPartitions" % len(self.nexusSets.charPartitions))
                for cp in self.nexusSets.charPartitions:
                    print("          %s" % cp.name)
            if self.nexusSets.charPartition:
                print("  The current charPartition is %s" % self.nexusSets.charPartition.name)
            else:
                print("  There is no current charPartition.")

        if self.excludeDelete:
            self.excludeDelete.dump()
        if len(self.parts):
            if len(self.parts) == 1:
                print("  There is %i part" % len(self.parts))
            else:
                print("  There are %i parts" % len(self.parts))
            for p in self.parts:
                print("          %s, length %i" % (p.name, p.nChar))
        # SequenceList.dump(self)
        print("  There are %i sequences" % len(self.sequences))
        upper = len(self.sequences)
        if upper > 5:
            upper = 5
        for i in range(upper):
            print("  %4i  %s" % (i, self.sequences[i].name))
        if len(self.sequences) > upper:
            print("       <... and more ...>")

    def setNexusSets(self):
        """Set self.nexusSets from var.nexusSets.

        A deepcopy is made of var.nexusSets, and then attached to
        self.  Sometimes other Nexus-set related methods trigger this.

        If var.nexusSets does not yet exist, a new blank one is made.
        """

        gm = ["Alignment.setNexusSets()"]
        if not var.nexusSets:
            var.nexusSets = NexusSets()
        self.nexusSets = copy.deepcopy(var.nexusSets)
        self.nexusSets.taxNames = self.taxNames
        self.nexusSets.nTax = self.nTax
        self.nexusSets.aligNChar = self.nChar

        self.nexusSets.constant.setAligNChar(self.nChar)
        self.nexusSets.gapped.setAligNChar(self.nChar)
        self.nexusSets.constant.mask = self.constantMask()
        self.nexusSets.gapped.mask = self.gappedMask()

        if self.nexusSets.charSets:
            for cs in self.nexusSets.charSets:
                cs.setAligNChar(self.nChar)
                cs.setMask()

        if self.nexusSets.taxSets:
            # print("%s. There are %i taxSets." % (gm[0], len(self.nexusSets.taxSets)))
            # Check that no taxSet name is a taxName
            lowSelfTaxNames = [txName.lower()
                               for txName in self.taxNames]
            for ts in self.nexusSets.taxSets:
                if ts.lowName in lowSelfTaxNames:
                    gm.append(
                        "Can't have taxSet names that are the same (case-insensitive) as a tax name")
                    gm.append(
                        "Lowercased taxSet name '%s' is the same as a lowcased taxName." % ts.name)
                    raise P4Error(gm)
            self.nexusSets.lowTaxNames = lowSelfTaxNames

            # If it is standard format,
            # convert triplets to numberTriplets, and then mask
            for ts in self.nexusSets.taxSets:
                if ts.format == 'standard':
                    ts.setNumberTriplets()
                    ts.setMask()
                    # print(ts.mask)
                elif ts.format == 'vector':
                    assert ts.mask
                    if len(ts.mask) != self.nTax:
                        gm.append("taxSet %s" % ts.name)
                        gm.append(
                            "It is vector format, but the length is wrong.")
                        gm.append(
                            "taxSet mask is length %i, but self nTax is %i" % (len(ts.mask), self.nTax))
                        raise P4Error(gm)
                else:
                    gm.append("taxSet %s" % ts.name)
                    gm.append("unknown format %s" % ts.format)
                    raise P4Error(gm)

                # Now set ts.taxNames from the mask.  Note self.taxNames is a property
                ts.taxNames = [self.taxNames[i] for i,c in enumerate(ts.mask) if c == '1']


        if self.nexusSets.charPartitions:
            pass

    def setCharPartition(self, charPartitionName):
        """Partition self into Parts based on charPartitionName.

        You need to do this before you ask for a Data object.

        You can also un-partition an already partitioned alignment by
        feeding this charPartitionName = None."""

        gm = ["Alignment.setCharPartition('%s')" % charPartitionName]
        if not self.nexusSets:
            if not charPartitionName:
                print(gm[0])
                print("setNexusSets() has not been done -- self has no nexusSets")
                print("yet we are doing setCharPartition, with no partition!")
                print("Its not an error, but are we a little confused?")
                return
            self.setNexusSets()
        self.nexusSets.charPartition = None
        self.parts = []

        if not charPartitionName:
            return

        for cp in self.nexusSets.charPartitions:
            if cp.name == charPartitionName:
                self.nexusSets.charPartition = cp
        if not self.nexusSets.charPartition:
            gm.append(
                "Could not find a CharPartition with the name '%s'" % charPartitionName)
            raise P4Error(gm)
        self.nexusSets.charPartition.setSubsetMasks()
        self.nexusSets.charPartition.checkForOverlaps()
        if 0:
            print(gm[0])
            self.nexusSets.charPartition.dump()
            self.dump()

    def changeDataTypeTo(self, newDataType, newSymbols=None, newEquates=None):
        """Coerce the alignment to be a new datatype.

        This would be good for pathological cases where eg DNA with
        lots of ambigs is mistaken for protein.

        It is not sufficient to simply change the dataType -- we must
        change the symbols and equates as well.  And the dataType for
        all the sequences in self.  If you are changing to 'standard'
        dataType, then you need to specify the symbols and the
        equates.  The symbols is a string, and the equates is a
        dictionary (see eg yourAlignment.equates for a DNA alignment
        to see the format).

        """

        gm = ['Alignment.changeDataTypeTo(%s, newSymbols=%s)' % (
            newDataType, newSymbols)]
        if newDataType not in ['dna', 'protein', 'standard']:
            gm.append(
                "newDataType must be one of 'dna', 'protein', 'standard'")
            raise P4Error(gm)

        if newDataType == self.dataType:
            gm.append("Self is already dataType %s" % self.dataType)
            raise P4Error(gm)

        if newDataType == 'standard':
            validChars = newSymbols + '-?' + ''.join(newEquates.keys())
            # print("standard datatype: got validChars '%s'" % validChars)

        for s in self.sequences:
            if newDataType == 'dna':
                for c in s.sequence:
                    if c not in var.validDnaChars:
                        gm.append(
                            "Sequence %s, char %s not a valid DNA character." % (s.name, c))
                        raise P4Error(gm)
                s.dataType = newDataType

            elif newDataType == 'protein':
                for c in s.sequence:
                    if c not in var.validProteinChars:
                        gm.append(
                            "Sequence %s, char %s not a valid protein character." % (s.name, c))
                        raise P4Error(gm)
                s.dataType = newDataType

            if newDataType == 'standard':
                for c in s.sequence:
                    if c not in validChars:
                        gm.append(
                            "Sequence %s, char '%s' not in valid chars '%s'." % (s.name, c, validChars))
                        raise P4Error(gm)
                s.dataType = newDataType

        self.dataType = newDataType
        if newDataType == 'dna':
            self.symbols = 'acgt'
            self.dim = 4
            self.equates = {'n': 'acgt', 'm': 'ac', 'k': 'gt',  # 'x': 'acgt',
                            'h': 'act', 'y': 'ct', 'v': 'acg',
                            'w': 'at', 'd': 'agt', 'b': 'cgt',
                            'r': 'ag', 's': 'cg'}
        elif newDataType == 'protein':
            self.symbols = 'arndcqeghilkmfpstwyv'
            self.dim = 20
            self.equates = {'b': 'dn', 'x': 'arndcqeghilkmfpstwyv', 'z': 'eq'}
        elif newDataType == 'standard':
            self.symbols = newSymbols
            self.dim = len(newSymbols)
            self.equates = newEquates

        # Is this needed?  Probably not.
        self.checkLengthsAndTypes()


    def entropyOfSites(self):
        """Simple entropy calc. Gaps and ambigs are ignored.

        Pure python.

        Returns:
            a list of entropy values, one for each site

        """

        siteEntropies = []
        counters = {}  # so that I can index it with a symbol rather than an int
        for c in self.symbols:
            counters[c] = 0
        for pos in range(len(self)):
            nCharsAtSite = 0
            keysHere = set()
            for seq in self.sequences:
                c = seq.sequence[pos]
                try:
                    counters[c] += 1  # ie ignore gaps and ambiguities
                    nCharsAtSite += 1
                    keysHere.add(c)
                except KeyError:
                    pass
            assert nCharsAtSite, f"a pos {pos} had no symbol chars, only gaps and ambigs."

            siteEntropy = 0.0
            totalCountsAtThisPos = 0
            for k in keysHere:
                v = counters[k]
                assert v
                totalCountsAtThisPos += v
                p = v / nCharsAtSite
                siteEntropy += p * math.log2(p)
                counters[k] = 0                # Initialize for next pos
            assert totalCountsAtThisPos, f"b pos {pos} had no symbol chars, only gaps and ambigs."
            assert totalCountsAtThisPos == nCharsAtSite
            if siteEntropy == 0.0:
                siteEntropies.append(0.0)     # to avoid -0.0
            else:
                siteEntropies.append(-siteEntropy)

        assert len(siteEntropies) == self.nChar
        return siteEntropies

    def entropyMeanAndVariance(self):
        """Calculate the mean and variance of the site entropy"""

        sll = self.entropyOfSites()
        m = statistics.mean(sll)
        v = statistics.variance(sll, m)

        return m,v


    #################################################### manip
    ####################################################

    def simpleConstantMask(self, ignoreGapQ=True, invert=False):
        """Returns a mask with 1 at simple constant sites, and 0 otherwise.

        See also :meth:`~p4.alignment.Alignment.constantMask`

        A simple constant site is a site that has only 1 kind of symbol, and no
        ambiguities.

        If ignoreGapQ (ignore '-' and '?') is turned on, then a site is a
        simple constant site if it has gaps and Q-marks, but otherwise
        only one kind of symbol.  If ignoreGapQ is turned off, then the
        presence of any '-' and '?' makes it a non-constant site.  If
        it is all gaps or Q-marks, then this method throws a P4Error.

        ::

          seq1  ac---cgt-
          seq2  acgtacgt-
          seq3  rcgtacgta
          mask  011111111   if ignoreGapQ=True
          mask  010001110   if ignoreGapQ=False


        All-gap columns are really undefined as to whether they are
        constant or not.  It gets worse if we invert to ask for
        potentially variable sites -- if an all-gap site is not constant,
        then it becomes potentially a variable site, which does not make
        sense to me.  Best to strip out all-gap sites first.

        """
        nTaxRange = range(self.nTax)
        mask = ['0'] * self.length
        for seqPos in range(self.length):
            theSlice = self.sequenceSlice(seqPos)
            # print("%2i: %s" % (seqPos, theSlice))

            # Is it all gaps and missing?  If so, its constant.
            nGapMiss = theSlice.count('-') + theSlice.count('?')
            if nGapMiss == len(theSlice):
                # print("    All miss-gap, ==> constant")
                gm = ["Alignment.simpleConstantMask()"]
                gm.append("All-gap site, position %i." % seqPos)
                gm.append("Get rid of it.")
                raise P4Error(gm)
            elif nGapMiss and not ignoreGapQ:
                mask[seqPos] = '0'
            else:
                firstSymb = None
                isSet = False
                for i in nTaxRange:
                    c = theSlice[i]
                    if c in self.symbols:
                        if not firstSymb:
                            firstSymb = c
                        else:
                            if c != firstSymb:
                                mask[seqPos] = '0'
                                isSet = True
                                break
                            else:
                                pass
                    elif c == '-' or c == '?':
                        pass
                    else:
                        # an ambig --> not simple constant
                        mask[seqPos] = '0'
                        isSet = True
                        break
                if not isSet:
                    mask[seqPos] = '1'

        if invert:
            for seqPos in range(self.length):
                if mask[seqPos] == '0':
                    mask[seqPos] = '1'
                elif mask[seqPos] == '1':
                    mask[seqPos] = '0'

        return ''.join(mask)

    def constantMask(self, invert=None):
        """Returns a mask string with 1 at constant sites, and 0 otherwise.

        See also :meth:`~p4.alignment.Alignment.simpleConstantMask`

        Constant sites are defined in a PAUP-like manner.  If, when all
        possibilities of ambiguities are tried, the site could possibly be
        constant, then it is constant.

        So if a site is all gaps or missing (or a combination), or if all
        the non-gap chars are the same, then it is a constant site.

        If a DNA site contains both r and y, it could not possibly be
        constant.  A site containing only a and r is constant, but a site
        containing only a and y cannot be constant.

        """

        #gm = ['Alignment.constantMask()']
        nSeq = len(self.sequences)

        # Make these lists once, so we do not spend all our time
        # malloc'ing and deleting the same lists over and over.
        symbolsSlice = ['a'] * nSeq
        equatesSlice = ['a'] * nSeq
        # Since these lists will not necessarily be filled, we need an
        # index for each.  These are nSymbolChars and nEquateChars, zeroed
        # below.

        # The general strategy here is to try to be fairly fast.  Options
        # for deciding whether a site is constant or not that are simple
        # and fast are tried first.  Only if there are ambiguities
        # (equates) does it get time-consuming.

        mask = ['0'] * self.length
        for seqPos in range(self.length):
            theSlice = self.sequenceSlice(seqPos)
            # print("%2i: %s" % (seqPos, theSlice))

            # Is it all gaps and missing?  If so, its constant.
            nGapMiss = theSlice.count('-') + theSlice.count('?')
            if nGapMiss == len(theSlice):
                # print("    All miss-gap, ==> constant")
                mask[seqPos] = '1'

            # if there is only 1 char that is not a gap or missing, then it is
            # constant
            elif nGapMiss == len(theSlice) - 1:
                # print("    Only 1 non-miss-gap, ==> constant")
                mask[seqPos] = '1'

            else:  # There are at least two non-gap chars

                # Divide up the non-miss-gap chars into symbols and equates
                nSymbolChars = 0
                nEquateChars = 0
                for i in range(nSeq):
                    aChar = theSlice[i]
                    # maybe I should include N? -- but I don't know if it is
                    # DNA.
                    if aChar not in '-?':
                        if aChar in self.symbols:
                            symbolsSlice[nSymbolChars] = aChar
                            nSymbolChars += 1
                        elif aChar in self.equates:
                            equatesSlice[nEquateChars] = aChar
                            nEquateChars += 1

                # If all the non-gap chars are symbols, then its easy.
                if not nEquateChars:
                    # print("    All (non-miss-gap) chars are symbols")
                    firstSymbol = symbolsSlice[0]
                    symbolsAreAllTheSame = 1
                    for i in range(1, nSymbolChars):
                        if symbolsSlice[i] != firstSymbol:
                            # print("        ...different symbols ==> not)
                            # constant"
                            symbolsAreAllTheSame = 0
                            break
                    if symbolsAreAllTheSame:
                        # print("        ... symbols all the same ==> constant")
                        mask[seqPos] = '1'

                else:  # We have equates
                    # print("    Some (non-miss-gap) chars are equates.")

                    if nSymbolChars:
                        # If there are different symbols, then it can't be
                        # constant
                        symbolsAreAllTheSame = 1
                        firstSymbol = symbolsSlice[0]
                        if nSymbolChars > 1:
                            for i in range(1, nSymbolChars):
                                if symbolsSlice[i] != firstSymbol:
                                    # print("        ...different symbols ==>)
                                    # not constant"
                                    symbolsAreAllTheSame = 0
                                    break
                        if symbolsAreAllTheSame:
                            # print("        ...symbols are all the same")
                            # But we cannot conclude that it is a constant
                            # site until we check the equates.  Which we
                            # now do.
                            symbolIndex = self.symbols.index(firstSymbol)
                            # print("firstSymbol=%s, symbolIndex = %s" % (firstSymbol, symbolIndex))
                            # print(self.equates)

                            # Here we make an array, eqArray, that
                            # contains coded info about what equates
                            # contain what symbols.  So for example, in
                            # DNA, n would be [1,1,1,1], and r would be
                            # [1,0,1,0], and so on.
                            eqArray = []
                            for eqNum in range(nEquateChars):
                                eq = equatesSlice[eqNum]
                                val = list(self.equates[eq])
                                # print("eq: %s  %s" % (eq, val))
                                oneLine = [0] * self.dim
                                for symbNum in range(self.dim):
                                    if self.symbols[symbNum] in val:
                                        oneLine[symbNum] = 1
                                eqArray.append(oneLine)
                            # print(eqArray)

                            allEquatesContainSymbol = 1  # to start
                            for i in range(nEquateChars):
                                # print(eqArray[i][symbolIndex])
                                if not eqArray[i][symbolIndex]:
                                    allEquatesContainSymbol = 0
                                    break
                            if allEquatesContainSymbol:
                                # print "        the equates all contain %s ==>
                                # constant" % firstSymbol
                                mask[seqPos] = '1'

                    else:
                        # print("        No symbols-- its all equates.")
                        firstEquate = equatesSlice[0]
                        equatesAreAllTheSame = 1
                        for i in range(1, nEquateChars):
                            if equatesSlice[i] != firstEquate:
                                # print("        ...different equates")
                                equatesAreAllTheSame = 0
                                break
                        if equatesAreAllTheSame:
                            # print "        ... equates all the same ==>
                            # constant"
                            mask[seqPos] = '1'
                        else:
                            # Need to test whether all equates "contain"
                            # the the same (any) symbol.  So, as above,
                            # make an eqArray.
                            eqArray = []
                            for eqNum in range(nEquateChars):
                                eq = equatesSlice[eqNum]
                                val = list(self.equates[eq])
                                # print("eq: %s  %s" % (eq, val))
                                oneLine = [0] * self.dim
                                for symbNum in range(self.dim):
                                    if self.symbols[symbNum] in val:
                                        oneLine[symbNum] = 1
                                eqArray.append(oneLine)
                            # print(eqArray)

                            # Now we go thru the eqArray column by column,
                            # and ask whether any column is all ones.
                            aColumnOfOnesExists = 0
                            for colNum in range(self.dim):
                                thisColIsAllOnes = 1
                                for eqNum in range(nEquateChars):
                                    if not eqArray[eqNum][colNum]:
                                        thisColIsAllOnes = 0
                                        break
                                if thisColIsAllOnes:
                                    aColumnOfOnesExists = 1
                                    break
                            if aColumnOfOnesExists:
                                # print "        the equates could be constant
                                # ==> constant"
                                mask[seqPos] = '1'

        if invert:
            for seqPos in range(self.length):
                if mask[seqPos] == '0':
                    mask[seqPos] = '1'
                elif mask[seqPos] == '1':
                    mask[seqPos] = '0'

        return ''.join(mask)

    def gappedMask(self, invert=None):
        """Returns a mask string with 1 at positions with any gaps, and 0 otherwise."""

        mask = ['0'] * self.length
        for i in range(self.length):
            theSlice = self.sequenceSlice(i)
            # if theSlice.count('-') == len(self.sequences):
            if theSlice.count('-'):
                mask[i] = '1'

        if invert:
            for i in range(self.length):
                if mask[i] == '0':
                    mask[i] = '1'
                elif mask[i] == '1':
                    mask[i] = '0'

        return ''.join(mask)

    def getLikelihoodTopologyInformativeSitesMask(self):
        """Make and return a mask for those sites that are likelihood informative about the topology.

        Mostly this means no constant and no autapomorphic sites.  
        `Autapomorphies <https://en.wikipedia.org/wiki/Autapomorphy>`_
        are sites that are all one character except for one taxon that
        has another character.

        The rules, this week--

        Its not informative if:

        - If there are no ambigs or gaps, then constants and autapomorphies are not informative.
        - If there are gaps but no ambigs,

            + if there are only 2 characters or fewer -- not informative
            + if constant + gaps -- not informative
            + if autapomorphy + gaps -- not informative

        - If there are ambigs but no gaps,

            + if constant except for a single ambig, then not informative
            + (an autapomorphy plus a single ambig can sometimes be informative)

        Otherwise, I'm saying that it is informative.  That includes

        - sites with gaps plus ambigs, and
        - sites with more than one ambig

        I am conservatively calling some sites as informative, eg sites
        with more than one ambiguity, even though sites like that are
        often not topologically informative, as they might be under some
        models-- I have not checked thoroughly.

        Because of the above, we cannot be sure that each site indicated
        by a 1 in the mask is informative, but we can be sure that each 0
        in the mask is uninformative. 

        Returns:
            a mask string
        """

        assert self.nTax > 2
        mask = ['0'] * len(self)
        if self.equates:
            equateKeys = list(self.equates)
        else:
            equateKeys = []
        counts = [0] * self.dim  # re-zero every loop
        hits = [0] * self.dim
        for sPos in range(len(self)):
            sl = self.sequenceSlice(sPos)
            for chNum in range(self.dim):
                counts[chNum] = 0
                # hits[chNum] = 0 done below, only if needed
            nGaps = 0
            nAmbigs = 0
            for tNum in range(self.nTax):
                c = sl[tNum]
                if c in self.symbols:
                    chNum = self.symbols.index(sl[tNum])
                    counts[chNum] += 1
                elif c in ['-', '?', 'n', 'x']:
                    nGaps += 1
                elif c in equateKeys:  # not n or x -- they have been done.
                    nAmbigs += 1
                else:
                    raise P4Error(
                        "can't deal with character '%s' at pos %i tax %i" % (c, sPos, tNum))

            if not nAmbigs and not nGaps:  # The simple, common case
                assert sum(counts) == self.nTax
                for chNum in range(self.dim):
                    hits[chNum] = 0
                for chNum in range(self.dim):
                    if counts[chNum]:
                        hits[chNum] = 1
                nDifferentChars = sum(hits)
                if nDifferentChars == 0:
                    raise P4Error(
                        "no nDifferentChars.  This should not happen.")
                elif nDifferentChars == 1:
                    # constant, mask stays as zero
                    pass
                elif nDifferentChars == 2:
                    if 1 in counts:  # a singleton
                        pass
                    else:
                        mask[sPos] = '1'
                else:
                    mask[sPos] = '1'
            elif nGaps and not nAmbigs:  # We got gaps.
                if 0:
                    mask[sPos] = '1'
                else:
                    totalNCharsInSlice = sum(counts)
                    # If 2 or less, mask stays as zero.
                    if totalNCharsInSlice > 2:
                        for chNum in range(self.dim):
                            hits[chNum] = 0
                        for chNum in range(self.dim):
                            if counts[chNum]:
                                hits[chNum] = 1
                        nDifferentChars = sum(hits)
                        if nDifferentChars == 0:
                            raise P4Error(
                                "no nDifferentChars.  This should not happen.")
                        elif nDifferentChars == 1:
                            # constant + gaps, mask stays as zero
                            pass
                        elif nDifferentChars == 2:
                            if 1 in counts:  # a singleton, plus gaps
                                pass
                            else:
                                mask[sPos] = '1'
                        else:
                            mask[sPos] = '1'
            elif nAmbigs and not nGaps:  # We got ambigs
                if 0:
                    mask[sPos] = '1'
                else:
                    totalNCharsInSlice = sum(counts)
                    if totalNCharsInSlice == self.nTax - 1:
                        assert nAmbigs == 1
                        for chNum in range(self.dim):
                            hits[chNum] = 0
                        for chNum in range(self.dim):
                            if counts[chNum]:
                                hits[chNum] = 1
                        nDifferentChars = sum(hits)
                        if nDifferentChars == 0:
                            raise P4Error(
                                "no nDifferentChars.  This should not happen.")
                        elif nDifferentChars == 1:
                            # constant + 1 ambig, mask stays as zero
                            # print "site %i, got constant + 1 ambig-- not
                            # topologically informative." % sPos
                            pass
                        elif nDifferentChars == 2:
                            # Sometimes a singleton + an ambig can be
                            # informative.
                            mask[sPos] = '1'
                            # if 1 in counts: # a singleton, plus 1 ambig, stays at zero
                            #    #print("site %i, got singleton + 1 ambig -- not topologically informative" % sPos)
                            #    mask[sPos] = '1'
                            # else:
                            #    mask[sPos] = '1'
                        else:
                            mask[sPos] = '1'
                    else:
                        mask[sPos] = '1'

            # We got ambigs and gaps.  Thats complicated.  So say that it is
            # informative.
            else:
                mask[sPos] = '1'

        mask = "".join(mask)
        return mask

    def getMaskForAutapomorphies(self):
        """Return a mask for autapomorphies 

        Returns a string of zeros and 1s, where the 1s are 
        `autapomorphic sites <https://en.wikipedia.org/wiki/Autapomorphy>`_.

        Gaps and ambiguities are ignored.

        In an alignment site, if there are two kinds of character states (ie the
        character diversity is 2) and one of the character states occurs once only,
        then it will be called an autapomorphy (or singleton).

        Returns:
            a string nChar long, of zeros and 1s.

        See also :meth:`~Alignment.getMaskForCharDiversity` and 
        :meth:`~Alignment.getLikelihoodTopologyInformativeSitesMask`

        """
        singletons = ['0'] * self.nChar
        counters = {}  # so that I can index it with a symbol rather than an int
        for c in self.symbols:
            counters[c] = 0
        for pos in range(len(self)):
            for seq in self.sequences:
                c = seq.sequence[pos]
                try:
                    counters[c] += 1  # ie not gaps and ambiguities
                except KeyError:
                    pass
            hitsAtThisPos = 0
            thereIsASingleton = False  
            for k,v in counters.items():
                if v:
                    hitsAtThisPos += 1
                    if v == 1:
                        thereIsASingleton = True
                counters[k] = 0                # Initialize for next pos
            assert hitsAtThisPos, "pos %i had no symbol chars, only gaps and ambigs."

            if hitsAtThisPos == 2 and thereIsASingleton:
                singletons[pos] = '1'
        return ''.join(singletons)



    def getMaskForCharDiversity(self, diversity=1):
        """Return a mask for a given character diversity

        Gaps and ambiguities are ignored.

        This makes a mask string with zeros and ones, with 1's showing the sites with
        the chosen diversity, and zero for other sites.

        So for example, if diversity is 1, that means constant sites, and this
        method returns a simple constant sites mask, with gaps and ambigs ignored.

        If the diversity is 2, it will return a mask showing sites with two kinds of
        characters.  And so on.

        Returns:
            a string nChar long, of zeros and 1s.

        See also :meth:`~Alignment.getMaskForAutapomorphies`
        """
        counters = {}  # so that I can index it with a symbol rather than an int
        myMask = []
        for c in self.symbols:
            counters[c] = 0
        for pos in range(len(self)):
            for seq in self.sequences:
                c = seq.sequence[pos]
                try:
                    counters[c] += 1  # ie ignore gaps and ambiguities
                except KeyError:
                    pass
            hitsAtThisPos = 0
            for k,v in counters.items():
                if v:
                    hitsAtThisPos += 1
                counters[k] = 0                # Initialize for next pos
            if hitsAtThisPos == diversity:
                myMask.append("1")
            else:
                myMask.append("0")
            assert hitsAtThisPos, "pos %i had no symbol chars, only gaps and ambigs."

        return ''.join(myMask)


    def getCharDiversityOverSites(self):
        """Return a list of character diversity values

        Gaps and ambiguities are ignored.

        Returns:
            a list of ints, nChar long

        See also :meth:`~Alignment.getMaskForCharDiversity`
        """
        counters = {}  # so that I can index it with a symbol rather than an int
        myCDList = []
        for c in self.symbols:
            counters[c] = 0
        for pos in range(len(self)):
            for seq in self.sequences:
                c = seq.sequence[pos]
                try:
                    counters[c] += 1  # ie ignore gaps and ambiguities
                except KeyError:
                    pass
            siteCD = 0
            for k,v in counters.items():
                if v:
                    siteCD += 1
                counters[k] = 0                # Initialize for next pos
            myCDList.append(siteCD)
            # assert siteCD, "pos %i had no symbol chars, only gaps and ambigs."

        return myCDList


    def getCharDiversityDistribution(self):
        """Distribution of different chars per site, as a distribution.

        Gaps and ambiguities are ignored.

        A tuple is returned, composed of

        +----------+---------------------------------------------------------+
        |    index |                                                        |
        +==========+=========================================================+
        |        0 | num of all-gap columns                                  |
        +----------+---------------------------------------------------------+
        |        1 | num of simple constant sites                            |
        +----------+---------------------------------------------------------+
        |        2 | num of sites with 2 kinds of char                       |
        +----------+---------------------------------------------------------+
        |        3 | num of sites with 3 kinds of char                       |
        +----------+---------------------------------------------------------+
        |      ... | ...                                                     |
        +----------+---------------------------------------------------------+
        | nSymbols | num of sites with some of each char                     |
        +----------+---------------------------------------------------------+


        So if there are 
        - no all-gap sites
        - 30 simple constant sites 
        - 10 sites with 2 kinds of char
        - and 20 sites with 3 kinds of character, 
        we would have a distro like this --- (0, 30, 10, 20, ...)

        This is in pure Python, and can be used for Alignments with one
        part.

        Returns:
            a tuple of ints, showing counts of sites with different diversities.

        See also :meth:`~Alignment.getCharDiversityOverSites`

        """
        #assert isinstance(self, Alignment)
        distro = [0] * (len(self.symbols) + 1)
        myCDList = self.getCharDiversityOverSites()

        for div in myCDList:
            distro[div] += 1

        # Check
        assert self.nChar == sum(distro), \
            "Something is wrong. The distribution does not add up to nChar"
        return tuple(distro)



    def orMasks(self, maskA, maskB):
        """Given two masks, this logically or's the string chars.

        | Only zero and '1' chars are allowed, and returned.
        | Eg. 0010 with 1000 will return 1010.
        | and 0010 with 1010 will return 1010.
        """

        gm = ["Alignment.orMasks()"]

        # check for silliness
        if not isinstance(maskA, str):
            gm.append("Alignment: orMasks(). Masks must be strings.")
            raise P4Error(gm)
        if not isinstance(maskB, str):
            gm.append("Alignment: orMasks(). Masks must be strings.")
            raise P4Error(gm)
        if len(maskA) != self.length or len(maskB) != self.length:
            gm.append("Masks must be the same length as the alignment.")
            raise P4Error(gm)
        l = self.length
        orMask = ['0'] * self.length
        for i in range(l):
            try:
                iA = int(maskA[i])
                iB = int(maskB[i])
            except ValueError:
                gm.append("All mask characters must be convertable to integers")
                raise P4Error(gm)
            if iA not in [0, 1] or iB not in [0, 1]:
                gm.append("All mask characters must be zero or 1")
                raise P4Error(gm)
            if iA or iB:
                orMask[i] = '1'
        return ''.join(orMask)

    def andMasks(self, maskA, maskB):
        """Given two masks, this logically and's the string chars.

        Only zero and '1' chars are allowed, and returned.

        Eg. 0010 with 1010 will return 0010.
        and 0010 with 1000 will return 0000.
        """
        gm = ['Alignment.andMasks']

        # check for silliness
        if not isinstance(maskA, str):
            gm.append("Masks must be strings.")
            raise P4Error(gm)
        if not isinstance(maskB, str):
            gm.append("Masks must be strings.")
            raise P4Error(gm)
        if len(maskA) != self.length or len(maskB) != self.length:
            gm.append("Masks must be the same length as the alignment.")
            raise P4Error(gm)

        andMask = ['0'] * self.length
        for i in range(self.length):
            try:
                iA = int(maskA[i])
                iB = int(maskB[i])
            except ValueError:
                gm.append("All mask characters must be convertable to integers")
                raise P4Error(gm)

            if iA not in [0, 1] or iB not in [0, 1]:
                gm.append("All mask characters must be zero or 1")
                raise P4Error(gm)
            if iA and iB:
                andMask[i] = '1'
        return ''.join(andMask)

    def sequenceSlice(self, pos):
        """Returns a list composed of the characters from the alignment at position pos.

        Args:
            pos (int): Zero-based position

        Returns:
            A list of the characters at that position

        """
        if self.length:
            if pos >= 0 and pos < self.length:
                sList = []
                for i in range(len(self.sequences)):
                    sList.append(self.sequences[i].sequence[pos])
                # return ''.join(sList)
                return sList
            else:
                raise P4Error("Alignment.sequenceSlice().  pos out of range")

    def bluntEndLigate(self, alig, allowDifferentDataTypes=False):
        """Attaches alig to the end of self.

        Various checks are made.  If the sequences are in a different
        order in the two alignments then it will not work.

        See the method :meth:`Alignment.concatenate`, which can handle missing
        and out-of-order sequences.

        """

        gm = ['Alignment.bluntEndLigate()']
        from p4.alignment import Alignment
        if not isinstance(alig, Alignment):
            gm.append("Arg must be an Alignment instance")
            raise P4Error(gm)
        if len(self.sequences) != len(alig.sequences):
            gm.append("Unequal number of sequences in the two alignments")
            raise P4Error(gm)
        elif self.length > 0 and alig.length == None:
            gm.append("self has sequence, but arg does not")
            raise P4Error(gm)
        elif self.length == None and alig.length > 0:
            gm.append("Arg has sequence, but self does not")
            raise P4Error(gm)
        if not allowDifferentDataTypes:
            if self.dataType != alig.dataType:
                gm.append("Arg and self dataTypes are different. (%s and %s)" % (
                    alig.dataType, self.dataType))
                gm.append(
                    "(If you really want this, you can set the arg allowDifferentDataTypes to True.)")
                raise P4Error(gm)
        for i in range(len(self.sequences)):
            if self.sequences[i].name != alig.sequences[i].name:
                gm.append("Name mismatch at zero-based sequence %i:" % i)
                gm.append("'%s' and '%s' don't match." %
                          (self.sequences[i].name, alig.sequences[i].name))
                raise P4Error(gm)
        if self.parts and len(self.parts):
            self.resetSequencesFromParts()
            self.parts = []
        if alig.parts and len(alig.parts):
            alig.resetSequencesFromParts()
        for i in range(len(self.sequences)):
            self.sequences[i].sequence = self.sequences[
                i].sequence + alig.sequences[i].sequence

        self.nexusSets = None
        self.length = len(self.sequences[0].sequence)
        self.checkLengthsAndTypes()

    def concatenate(self, alig, sNames):
        """Attaches alig to the end of self.

        You need to provide a list, argument sNames, of taxon names that are
        found in self and alig.  This will determine the order of the taxa in
        self, the result.

        It will still work if the sequences are in a different
        order in the two alignments.

        The order of the taxa in sNames need not be the same order as either
        alignment.  So this method can change the order of the sequences in self.

        It will still work if there are missing sequences in either self or
        alig.  It will add blank sequences as needed.

        Args:
            alig (Alignment):  The other Alignment object.
            sNames (list):  A list of all the taxon names found in self and alig

        """

        gm = ['Alignment.concatenate()']
        from p4.alignment import Alignment
        if not isinstance(alig, Alignment):
            gm.append("Arg must be an Alignment instance")
            raise P4Error(gm)
        elif self.length > 0 and not alig.length:
            gm.append("self has sequence, but arg does not")
            raise P4Error(gm)
        # elif alig.length > 0 and not self.length:
        #     gm.append("Arg has sequence, but self does not")
        #     raise P4Error(gm)
        if self.parts and len(self.parts):
            self.resetSequencesFromParts()
            self.parts = []
        if alig.parts and len(alig.parts):
            alig.resetSequencesFromParts()
        for tName in self.taxNames:
            assert tName in sNames, "self name %s is not in sNames." % tName
        for tName in alig.taxNames:
            assert tName in sNames, "other alig name %s is not in sNames." % tName
        if not self.sequenceForNameDict:
            self.makeSequenceForNameDict()
        if not alig.sequenceForNameDict:
            alig.makeSequenceForNameDict()
        newSequences = []
        for sName in sNames:
            # print sName,
            selfSeq = self.sequenceForNameDict.get(sName)
            # print selfSeq,
            if not selfSeq:
                selfSeq = self.sequences[0].dupe()
                selfSeq.name = sName
                selfSeq.sequence = '-' * self.nChar
                self.sequenceForNameDict[sName] = selfSeq
            # print selfSeq.sequence,
            aligSeq = alig.sequenceForNameDict.get(sName)
            # print aligSeq,
            if not aligSeq:
                aligSeq = alig.sequences[0].dupe()
                aligSeq.name = sName
                aligSeq.sequence = '-' * alig.nChar
            # print(aligSeq.sequence)
            selfSeq.sequence += aligSeq.sequence
            newSequences.append(selfSeq)
        self.sequences = newSequences
        self.nexusSets = None
        self.length = len(self.sequences[0].sequence)
        self.checkLengthsAndTypes()

    def constantSitesProportion(self):
        """Returns the proportion of the alignment that have possible constant sites.

        It uses constantSitesCount()
        """

        return float(self.constantSitesCount()) / float(self.length)

    def constantSitesCount(self):
        """Counts the sites that potentially have the same thing in each sequence.

        Constant sites are defined in a PAUP-like manner.  See the doc
        string for Alignment.constantMask().
        """

        if 0:
            # This works if there is is no var.nexusSets.  However, this
            # goes boom if there is an inappropriate var.nexusSets, maybe
            # leftover from something else.
            if self.nexusSets and self.nexusSets.charSets:
                pass
            else:
                self.setNexusSets()

            assert self.nexusSets
            assert self.nexusSets.constant
            constCS = self.nexusSets.constant
            if not constCS:
                gm = ['Alignment.constantSitesCount()']
                gm.append("Could not find the charSet 'constant'.  Fix me.")
                raise P4Error(gm)

            if not constCS.mask:
                constCS.setMask(self.nexusSets, self)
            return constCS.mask.count('1')

        return self.constantMask().count('1')

    def noGapsOrAmbiguitiesCopy(self):
        """Returns a new Alignment with sites with gaps or ambiguities removed.

        This makes another Alignment instance by copying over the
        sequences site by site as long as there are no gaps or
        ambiguities.  """

        from p4.alignment import Alignment
        dbug = 0
        seqCount = len(self.sequences)
        newAlig = Alignment()
        newAlig.dataType = self.dataType
        newAlig.symbols = self.symbols
        newAlig.equates = self.equates

        for s in self.sequences:
            newSeq = Sequence()
            newSeq.name = s.name
            newSeq.comment = s.comment
            #newSeq.transl_table = s.transl_table
            newSeq.dataType = s.dataType
            newSeq.sequence = []        # a list!
            newAlig.sequences.append(newSeq)

        for i in range(self.length):
            theSlice = self.sequenceSlice(i)
            useIt = 1
            if self.sequences[0].dataType == 'dna':
                for j in theSlice:
                    if j not in 'acgt':
                        if dbug:
                            print(j)
                        useIt = 0
                        if not dbug:
                            break
                if useIt:
                    for j in range(seqCount):
                        newAlig.sequences[j].sequence.append(theSlice[j])
            elif self.sequences[0].dataType == 'protein':
                for j in theSlice:
                    if j not in 'acdefghiklmnpqrstvwy':
                        if dbug:
                            print(j)
                        useIt = 0
                        if not dbug:
                            break
                if useIt:
                    for j in range(seqCount):
                        newAlig.sequences[j].sequence.append(theSlice[j])
            elif self.sequences[0].dataType == 'standard':
                for j in theSlice:
                    if j not in self.symbols:
                        if dbug:
                            print(j)
                        useIt = 0
                        if not dbug:
                            break
                if useIt:
                    for j in range(seqCount):
                        newAlig.sequences[j].sequence.append(theSlice[j])
        for j in range(seqCount):
            newAlig.sequences[j].sequence = ''.join(newAlig.sequences[j].sequence)
        newAlig.checkLengthsAndTypes()
        return newAlig

    def hasGapsOrAmbiguities(self):
        """Asks whether self has any gaps or ambiguities."""

        ambigs = list(self.equates.keys())
        ambigs.append('-')
        ambigs = ''.join(ambigs)
        # print("got ambigs = '%s'" % ambigs)

        for s in self.sequences:
            for c in s.sequence:
                if c in ambigs:
                    return True
        return False

    def bootstrap(self):
        """Returns a new Alignment, a bootstrap resampling of self.

        This is done only in Python (ie no cParts involved), and does
        not handle partitioned data.  If you want to bootstrap
        partitioned data, use the :meth:`p4.data.Data.bootstrap`
        method.

        To make this exactly repeatable, set the seed for the Python
        random module::

            import random
            random.seed(mySeed)
        """

        gm = ["Alignment.bootstrap()"]

        # Check that there is only one partition
        # its all one part
        if not self.nexusSets or not self.nexusSets.charPartition:
            pass
        # its partitioned.  Bad.
        elif self.nexusSets.charPartition and len(self.nexusSets.charPartition.subsets) > 1:
            gm.append(
                "This only works with Alignments having only one data partition.")
            raise P4Error(gm)

        # although we will be replacing the sequences...
        a = copy.deepcopy(self)
        n = len(self.sequences)

        # make a 2D array the same size as the sequences, filled.
        newList = []
        for i in range(n):
            one = ['a'] * self.length
            newList.append(one)
        # fill the array with random slices from the sequences
        for j in range(self.length):
            r = int(self.length * random.random())
            for i in range(n):
                newList[i][j] = self.sequences[i].sequence[r]
        # replace the sequences
        for i in range(n):
            a.sequences[i].sequence = ''.join(newList[i])
        return a

    def compositionEuclideanDistanceMatrix(self):
        """This returns a DistanceMatrix based on composition.

        The formula is as given in Lockhart et al 94, the logDet paper.
        Its equation 4 there, page 608.  One pairwise distance is the
        square root of the sum of the squares of the differences between
        the frequencies of character states of one pair of sequences.
        """

        d = DistanceMatrix()
        d.setDim(len(self.sequences))
        d.names = []
        for s in self.sequences:
            d.names.append(s.name)
        compList = []
        for i in range(len(self.sequences)):
            compList.append(self.composition([i]))
        # print(compList)
        for i in range(len(self.sequences)):
            for j in range(len(self.sequences))[i:]:
                s = 0.0
                for k in range(self.dim):
                    diff = compList[i][k] - compList[j][k]
                    s += diff * diff
                s = math.sqrt(s)
                d.matrix[i][j] = s
                d.matrix[j][i] = s
        return d

    def covarionStats(self, listA, listB, verbose=True):
        """Calculates some covarion statistics.

        After Pete Lockhart's covarion paper.  Reference?  Looks at the
        two groups of sequences defined by listA and listB, and determines
        the number of sites that are constant overall, constant in both
        groups but in a different character, constant in one but variable
        in the other, variable in one and constant in the other, and
        variable overall.  The lists can be either sequence numbers
        (zero-based) or names.  If verbose, (the default) it prints a nice
        summary, with a tiny explanation.  Returns a tuple."""

        gm = ['Alignment.covarionStats()']

        # listA and listB should be lists
        if not isinstance(listA, list) or not isinstance(listB, list):
            gm.append("The args should be lists of sequences numbers or names")
            raise P4Error(gm)
        lstA = []
        lstB = []
        for i in listA:
            if isinstance(i, str):
                it = None
                for s in range(len(self.sequences)):
                    if self.sequences[s].name == i:
                        it = s
                if it == None:
                    gm.append("Name '%s' is not in self." % i)
                    raise P4Error(gm)
                lstA.append(it)
            elif isinstance(i, int) and i >= 0 and i < len(self.sequences):
                lstA.append(i)
            else:
                gm.append(
                    "The args should be lists of sequences numbers or names")
                raise P4Error(gm)

        for i in listB:
            if isinstance(i, str):
                it = None
                for s in range(len(self.sequences)):
                    if self.sequences[s].name == i:
                        it = s
                if it == None:
                    gm.append("Name '%s' is not in self." % i)
                    raise P4Error(gm)
                lstB.append(it)
            elif isinstance(i, int) and i >= 0 and i < len(self.sequences):
                lstB.append(i)
            else:
                gm.append(
                    "The args should be lists of sequences numbers or names")
                raise P4Error(gm)

        for i in lstA:
            if i in lstB:
                gm.append(
                    "The arg lists overlap: sequence %i appears in both" % i)
                raise P4Error(gm)
        for i in lstB:
            if i in lstA:
                gm.append(
                    "The arg lists overlap: sequence %i appears in both" % i)
                raise P4Error(gm)

        # c1 = sames overall
        # c2 = sames in a, sames in b, but a is not the same as b
        # c3 = sames in a, differents in b
        # c4 = differents in a, sames in b
        # c5 = differents in both
        c1 = 0
        c2 = 0
        c3 = 0
        c4 = 0
        c5 = 0
        for seqPos in range(self.length):
            a = []
            b = []
            for i in lstA:
                a.append(self.sequences[i].sequence[seqPos])
            for i in lstB:
                b.append(self.sequences[i].sequence[seqPos])
            if a.count(a[0]) == len(a):
                if b.count(b[0]) == len(b):
                    if a[0] == b[0]:
                        c1 = c1 + 1
                    else:
                        c2 = c2 + 1
                else:
                    c3 = c3 + 1
            else:
                if b.count(b[0]) == len(b):
                    c4 = c4 + 1
                else:
                    c5 = c5 + 1

        if verbose:
            format = '%45s %i'
            print("\nCovarion stats")
            print(format % ('each group has sames, with the same char', c1))
            print(format % ('each group has sames, but a different char', c2))
            print(format % ('sames in A, differents in B', c3))
            print(format % ('differents in A, sames in B', c4))
            print(format % ('differents in both groups', c5))

        return (c1, c2, c3, c4, c5)

    def pDistances(self, ignoreGaps=True, divideByNPositionsCompared=True):
        """Returns a DistanceMatrix of mean character distances or pDistances.

        The default behaviour is that only pairwise positions are compared
        in which both sequences have a character.  The sum of differences
        is then divided by the number of positions compared.

        Its not the same as paup p-dists when there are gaps or
        ambiguities.  For example, the pDistances between the sequences
        'acgaa' and 'aaa--' is 0.4 in paup, but 0.667 here.  Paup would
        call these p4 distances 'mean character distances' (dset
        distance=mean).

        Ambiguities are not handled intelligently.  Ambigs are just
        another character.  So r-a is considered to be different, and so
        contributes to the distance.

        If ignoreGaps is set, then any position in a sequence pair that
        contains a gap in either sequence is ignored.  If ignoreGaps is
        set to False, then those positions are looked at, and gap-a would
        be considered a difference.

        If divideByNPositionsCompared is turned off, then the number of
        differences is divided by the length of the alignment.  """

        d = DistanceMatrix()
        d.setDim(len(self.sequences))
        d.names = []
        for s in self.sequences:
            d.names.append(s.name)

        for i in range(len(self.sequences) - 1):
            for j in range(i + 1, len(self.sequences)):
                nDiffs = 0
                nPositions = 0
                for k in range(self.length):
                    if ignoreGaps:
                        if self.sequences[i].sequence[k] != '-' and self.sequences[j].sequence[k] != '-':
                            nPositions += 1
                            if self.sequences[i].sequence[k] != self.sequences[j].sequence[k]:
                                nDiffs += 1
                    else:
                        nPositions += 1
                        if self.sequences[i].sequence[k] != self.sequences[j].sequence[k]:
                            nDiffs += 1
                if divideByNPositionsCompared:
                    if not nPositions:
                        print("No shared positions between (zero-based) seqs %i and %i.  Setting to 1.0" % (i, j))
                        fDiffs = 1.0
                    else:
                        fDiffs = float(nDiffs) / float(nPositions)
                else:
                    fDiffs = float(nDiffs) / float(self.length)
                # Uncomment the following to return the number, for just the first 2 sequences.
                # return fDiffs
                d.matrix[i][j] = fDiffs
                d.matrix[j][i] = fDiffs
        return d

    def recodeDayhoff(self, firstLetter=False, symbols="123456"):
        """Recode protein data into Dayhoff groups, in place.

            1. c

            2. stpag

            3. ndeq

            4. hrk

            5. milv

            6. fyw

        The ambiguity character 'x' is recoded as a gap.

        It does not make a new alignment-- it does the re-coding
        'in-place'.

        If arg *firstLetter* is set, then the character is recoded as
        the first letter of its group rather than as a number.  Eg k
        would be recoded as h rather than as 4.

        By default, the six symbols are "123456".  Setting
        *firstLetter* is equivalent to setting the symbols to be
        "csnhmf".  The symbols can be set to other values; for example
        for RAxML input you would want to set the symbols to "012345".
        """

        gm = ['Alignment.recodeDayhoff()']
        if self.dataType != 'protein':
            gm.append("This is only for protein alignments.")
            raise P4Error(gm)

        if firstLetter:
            symbols = "csnhmf"
        assert len(symbols) == 6
            

        for s in self.sequences:
            s.dataType = 'standard'
            s.sequence = list(s.sequence)
            for i in range(len(s.sequence)):
                c = s.sequence[i]
                if c == 'c':
                    s.sequence[i] = symbols[0]
                elif c in 'stpag':
                    s.sequence[i] = symbols[1]
                elif c in 'ndeq':
                    s.sequence[i] = symbols[2]
                elif c in 'hrk':
                    s.sequence[i] = symbols[3]
                elif c in 'milv':
                    s.sequence[i] = symbols[4]
                elif c in 'fyw':
                    s.sequence[i] = symbols[5]
                elif c in ['-', '?']:
                    pass  # They stay as they are.
                elif c == 'x':
                    s.sequence[i] = '-'
                else:
                    # Maybe this should raise a P4Error?
                    print("skipping character '%s'" % c)
            s.sequence = ''.join(s.sequence)
        self.dataType = 'standard'
        self.equates = {}
        self.dim = 6
        self.symbols = symbols

    def recodeProteinIntoGroups(self, groups, firstLetter=False):
        """Recode protein data into user-specified groups, in place.

        A generalization of :meth:`p4.alignment.Alignment.recodeDayhoff`

        The arg *groups* should be a list of strings indicating the
        groupings, with all AAs present.  Case does not matter.

        The ambiguity character 'x' is recoded as a gap.

        It does not make a new alignment-- it does the re-coding 'in-place'.

        If arg *firstLetter* is set, then the character is recoded as the
        first letter of its group rather than as a number. 
        """

        gm = ['Alignment.recodeProteinIntoGroups()']
        if self.dataType != 'protein':
            gm.append("This is only for protein alignments.")
            raise P4Error(gm)

        assert isinstance(groups, list)
        nGroups = len(groups)
        assert nGroups > 1
        assert nGroups < 20
        for gr in groups:
            assert isinstance(gr, str)
        myGroups = [gr.lower() for gr in groups]
        theseSymbols = ''.join(myGroups)
        assert len(theseSymbols) == 20
        for s in theseSymbols:
            assert s in self.symbols
            assert theseSymbols.count(s) == 1
        numeralSymbols = ['%i' % (i + 1) for i in range(nGroups)]
        firstLetters = [gr[0] for gr in myGroups]

        for s in self.sequences:
            s.dataType = 'standard'
            s.sequence = list(s.sequence)
            for i in range(len(s.sequence)):
                c = s.sequence[i]
                gotIt = False
                for grNum in range(nGroups):
                    gr = myGroups[grNum]
                    if c in gr:
                        if firstLetter:
                            s.sequence[i] = firstLetters[grNum]
                        else:
                            s.sequence[i] = numeralSymbols[grNum]
                        gotIt = True
                        break
                if not gotIt:
                    if c in ['-', '?']:
                        pass  # They stay as they are.
                    elif c == 'x':
                        s.sequence[i] = '-'
                    else:
                        # Maybe this should raise a P4Error?
                        print("skipping character '%s'" % c)
            s.sequence = ''.join(s.sequence)
        self.dataType = 'standard'
        self.equates = {}
        self.dim = nGroups
        if firstLetter:
            self.symbols = ''.join(firstLetters)
        if not firstLetter:
            self.symbols = ''.join(numeralSymbols)

    def recodeRY(self, ambigsBecomeGaps=True, use01=False):
        """Recode DNA data into purines and pyrimidines, in place.

        So A and G becomes R, and C and T becomes Y.  Gaps remain
        gaps, missing (?) remains missing, and Rs and Ys remain as
        they are.  Depending on the setting of ambigsBecomeGaps, Ns
        may also remain unmodified, but any other characters (DNA
        ambiguity characters) become either gaps or N, depending on
        the setting of ambigsBecomeGaps.  So if ambigsBecomeGaps is
        turned on, as it is by default, then N, S, M, and so on become
        '-' ie gap characters.  If ambigsBecomeGaps is turned off,
        then they all become Ns.  If ambigsBecomeGaps is turned off,
        then the alignment gets an equate, of N for R or Y.

        If use01 is turned on (it is off by default) then we use 0 and
        1 for characters, rather than R and Y.  This is a bit buggy
        because any existing R or Y characters are ignored, and other
        ambigs are encoded as N.

        The dataType becomes 'standard' and the dim becomes 2.

        It does not make a new alignment-- it does the re-coding
        'in-place'.
        """

        gm = ['Alignment.recodeRY()']
        if self.dataType != 'dna':
            gm.append("This is only for dna alignments.")
            raise P4Error(gm)

        if use01:
            theR = '0'
            theY = '1'
        else:
            theR = 'r'
            theY = 'y'

        for s in self.sequences:
            s.dataType = 'standard'
            s.sequence = list(s.sequence)
            for i in range(len(s.sequence)):
                c = s.sequence[i]
                if c in 'ag':
                    s.sequence[i] = theR
                elif c in 'ct':
                    s.sequence[i] = theY
                elif c in 'ry-?':
                    pass
                else:
                    if ambigsBecomeGaps:
                        s.sequence[i] = '-'
                    else:
                        s.sequence[i] = 'n'
            s.sequence = ''.join(s.sequence)
        self.dataType = 'standard'
        if use01:
            self.symbols = '01'
        else:
            self.symbols = 'ry'
        if not ambigsBecomeGaps:
            self.equates = {'n': 'ry'}
        else:
            self.equates = {}
        self.dim = 2

    def checkTranslation(self, theProteinAlignment, transl_table=1, checkStarts=False):
        """Check that self translates to theProteinAlignment.

        Self should be a DNA alignment.  It is translated using
        :meth:`p4.geneticcode.GeneticCode.translate` (so it should handle
        ambiguities) and compared against theProteinAlignment.  The
        theProteinAlignment sequence order, names, and gap pattern should
        be the same as in the DNA alignment.  The default transl_table is
        the standard (or so-called universal) genetic code.

        Other available translation tables, this week::

            if transl_table == 1:   # standard
            elif transl_table == 2: # vertebrate mito
            elif transl_table == 4: # Mold, Protozoan,
                                    # and Coelenterate Mitochondrial Code
                                    # and the Mycoplasma/Spiroplasma Code
            elif transl_table == 5: # invertebrate mito
            elif transl_table == 9: # echinoderm mito

            # and now 6, 10, 11, 12, 13, 14, 21.

        (These are found in :class:`~p4.geneticcode.GeneticCode`)

        See also :meth:`p4.alignment.Alignment.translate`

        If the arg *checkStarts* is turned on (by default it is not turned
        on) then this method checks whether the first codon is a start
        codon.

        """

        gm = ['Alignment.checkTranslation()']
        if self.dataType != 'dna':
            gm.append("Self should be a DNA alignment.")
            raise P4Error(gm)
        if not theProteinAlignment or \
                not isinstance(theProteinAlignment, p4.alignment.Alignment) or \
                theProteinAlignment.dataType != 'protein':
            gm.append("Something wrong with theProteinAlignment")
            raise P4Error(gm)

        if len(self.sequences) != len(theProteinAlignment.sequences):
            gm.append(
                "Self and theProteinAlignment have different numbers of sequences")
            raise P4Error(gm)

        for seqNum in range(len(self.sequences)):
            s1 = self.sequences[seqNum]
            s2 = theProteinAlignment.sequences[seqNum]
            if s1.name != s2.name:
                gm.append(
                    "The sequence names of self and theProteinAlignment are not the same")
                raise P4Error(gm)

        if self.length != (3 * theProteinAlignment.length):
            gm.append(
                "The length of the DNA alignment should be 3 times that of theProteinAlignment")
            gm.append("DNA alignment (self):  %i" % self.length)
            gm.append("Protein alignment:     %i  ( * 3 = %i)" %
                      (theProteinAlignment.length, (3 * theProteinAlignment.length)))
            raise P4Error(gm)

        gc = GeneticCode(transl_table)

        pLen = theProteinAlignment.length
        for i in range(len(self.sequences)):
            s1 = self.sequences[i]
            s2 = theProteinAlignment.sequences[i]
            print("Checking %s ..." % s1.name)
            crimes = 0
            for j in range(pLen):
                theCodon = s1.sequence[(3 * j) + 0] + \
                    s1.sequence[(3 * j) + 1] + \
                    s1.sequence[(3 * j) + 2]
                if theCodon == '---':
                    if s2.sequence[j] != '-':
                        print("    position %4i, codon '---' is '%s', should be '-'" % (j, s2.sequence[j]))
                        crimes += 1
                elif theCodon.count('-'):
                    print("    position %4i, codon '%s' is incomplete" % (j, theCodon))
                    crimes += 1
                # elif theCodon in gc.code:
                #     if gc.code[theCodon] != s2.sequence[j]:
                #         print "    position %4i, codon '%s' is '%s', should be '%s'" % (
                #             j, theCodon, s2.sequence[j], gc.code[theCodon])
                #         crimes += 1
                # else:
                #     print("    position %4i, codon '%s' is not a known codon" % (j, theCodon))
                #     crimes += 1
                else:
                    tr = gc.translate(theCodon)
                    if tr != s2.sequence[j]:
                        print("    position %4i, codon '%s' is '%s', should be '%s'" % (
                            j, theCodon, s2.sequence[j], gc.code[theCodon]))
                        crimes += 1

                    # If arg checkStarts is turned on -- Is the first
                    # codon a start?  -- if not, it is not a crime
                    if checkStarts and j == 0:
                        if theCodon in gc.startList:
                            print("    Seq %i (%s). The first codon, '%s', is a start codon" % (i, s1.name, theCodon))
                        else:
                            print("    Seq %i (%s). The first codon, '%s', is not a start codon" % (i, s1.name, theCodon))
                if crimes > 6:
                    break
            if crimes > 6:
                print("    ... and possibly others, skipped.")

    def translate(self, transl_table=1, checkStarts=False, nnn_is_gap=False):
        """Returns a protein alignment from self, a DNA alignment.

        Self is translated using
        :meth:`GeneticCode.GeneticCode.translate`, so it handles
        ambiguities.  At the moment, we can only do translations where the
        frame of the codon is 123, ie the first sequence position is the
        first position of the codon.  The default transl_table is the
        standard (or so-called universal) genetic code, but you can change
        it.

        Other available translation tables, this week::

            if transl_table == 1: # standard
            elif transl_table == 2: # vertebrate mito
            elif transl_table == 4: # Mold, Protozoan,
                                    # and Coelenterate Mitochondrial Code
                                    # and the Mycoplasma/Spiroplasma Code
            elif transl_table == 5: # invertebrate mito
            elif transl_table == 9: # echinoderm mito

            and now 6, 10, 11, 12, 13, 14, 21.

        (These are found in :class:`p4.geneticcode.GeneticCode`)

        See also :meth:`p4.alignment.Alignment.checkTranslation`.

        If the arg *checkStarts* is turned on (by default it is not turned
        on) then this method checks whether the first codon is a start
        codon.

        Arg *nnn_is_gap* is for odd alignments where there are long
        stretches of 'nnn' codons, which probably should be gaps.
        Probably best to correct those elsewise.

        """

        gm = ['Alignment.translate()']
        if self.dataType != 'dna':
            gm.append("Self should be a DNA alignment")
            raise P4Error(gm)

        if self.length % 3 != 0:
            gm.append("The length of self should be a multiple of 3")
            raise P4Error(gm)

        a = self.dupe()
        a.dataType = 'protein'
        a.length = self.length / 3
        a.symbols = 'arndcqeghilkmfpstwyv'
        a.equates = {'b': 'dn', 'x': 'arndcqeghilkmfpstwyv', 'z': 'eq'}
        a.dim = 20
        a.nexusSets = None
        a.parts = []
        a.excludeDelete = None
        for s in a.sequences:
            s.sequence = ['-'] * a.length
            s.dataType = 'protein'

        gc = GeneticCode(transl_table)

        # dnaEquates = self.equates.keys()
        # print dnaEquates  # ['b', 'd', 'h', 'k', 'm', 'n', 's', 'r', 'w', 'v', 'y']

        for i in range(len(self.sequences)):
            dnaSeq = self.sequences[i].sequence
            # self.sequences[i].writeFasta()
            protSeq = a.sequences[i].sequence
            for j in range(a.length):
                theCodon = dnaSeq[(j * 3):(j * 3) + 3]
                # print(theCodon)
                if theCodon == '---':
                    protSeq[j] = '-'
                elif theCodon.count('-'):
                    print("    seq %i, position %4i, dnaSeq %4i, codon '%s' is incomplete" % (i, j, (j * 3), theCodon))
                elif theCodon == 'nnn':
                    if nnn_is_gap:
                        print("    seq %i, position %4i, dnaSeq %4i, codon '%s' translating to a gap ('-')" % (i, j, (j * 3), theCodon))
                        protSeq[j] = '-'
                    else:
                        protSeq[j] = 'x'
                else:
                    protSeq[j] = gc.translate(theCodon)
                    # print("    seq %i position %4i, dnaSeq %4i, codon '%s' is not a known codon -- using x" % (i, j, (j*3), theCodon))
                    #protSeq[j] = 'x'
                    if checkStarts and j == 0:
                        if theCodon in gc.startList:
                            print("    Seq %i (%s). The first codon, '%s', is a start codon" % (
                                i, self.sequences[i].name, theCodon))
                        else:
                            print("    Seq %i (%s). The first codon, '%s', is not a start codon" % (
                                i, self.sequences[i].name, theCodon))

        for s in a.sequences:
            s.sequence = ''.join(s.sequence)
            # print(s.sequence)
        return a

    def excludeCharSet(self, theCharSetName):
        """Exclude a CharSet."""

        gm = ['Alignment.excludeCharSet()']
        if not self.nexusSets:
            self.setNexusSets()
        lowName = theCharSetName.lower()
        theCS = None

        # We have either a pre-defined or non pre-defined char set
        if lowName in ['constant', 'gapped']:
            if lowName == 'constant':
                theCS = self.nexusSets.constant
            else:
                theCS = self.nexusSets.gapped
        else:
            if not len(self.nexusSets.charSets):
                gm.append("This alignment has no non-pre-defined charSets")
                raise P4Error(gm)
            for cs in self.nexusSets.charSets:
                if cs.lowName == lowName:
                    theCS = cs
                    break

        if theCS == None:
            gm.append("This alignment has no charset named '%s'" %
                      theCharSetName)
            raise P4Error(gm)
        if theCS.aligNChar == None:
            if self.excludeDelete:
                theCS.setAligNChar(self.excludeDelete.length)
            else:
                theCS.setAligNChar(self.length)

        # prepare the mask
        if not theCS.mask:
            theCS.setMask(self.nexusSets, self)
        # print("The mask is: %s" % theCS.mask)

        if not self.excludeDelete:
            self.excludeDelete = ExcludeDelete(self)

        if theCS not in self.excludeDelete.excludedCharSets:
            self.excludeDelete.excludedCharSets.append(theCS)
            self.excludeDelete._resetMask()
            self.excludeDelete._resetSequences()
            # self.excludeDelete.dump()
        else:
            print(gm[0])
            print("    %s has already been excluded." % theCharSetName)
        self.parts = []

    def dupe(self):
        """Duplicates self, with no c-pointers.  And no parts"""
        theDupe = copy.deepcopy(self)
        for p in theDupe.parts:
           p.alignment = theDupe
           p.cPart = None
        for p in theDupe.parts:
            del(p)
        theDupe.parts = []
        return theDupe

    def putGaps(self, theDnaSequenceList):
        """Insert gaps in theDnaSequenceList based on gaps in self.

        Like James O. McInerney's 'putgaps' program.  Self should be a
        protein alignment.  The DNA is input as a SequenceList object,
        'theDnaSequenceList'.  It creates and returns a new DNA alignment.

        For example::

            read('myProteinAlignment.nex')
            pAlign = var.alignments[0]
            read('myUnalignedDna.fasta')
            sl = var.sequenceLists[0]
            newDnaAlign = pAlign.putGaps(sl)
            newDnaAlign.writeNexus('d.nex')

        Args:
            a DNA SequenceList object

        Returns:
            a new Alignment object

        """

        gm = ['Alignment.putGaps()']
        if self.dataType != 'protein':
            gm.append("self should be a protein alignment.")
            raise P4Error(gm)

        for s in theDnaSequenceList.sequences:
            if s.sequence.count('-'):
                gm.append("DNA sequence %s already has gaps." % s.name)
                raise P4Error(gm)
            if len(s.sequence) % 3 != 0:
                gm.append("DNA sequence %s" % s.name)
                gm.append("Length %i is not evenly divisible by 3." %
                          len(s.sequence))
                raise P4Error(gm)
            if (len(s.sequence) / 3) > self.length:
                gm.append("DNA sequence %s" % s.name)
                gm.append("Length %i is more than 3 times the protein length %i." % (
                    len(s.sequence), self.length))
                raise P4Error(gm)
            if s.dataType != 'dna':
                gm.append("DNA(?!?) sequence %s" % s.name)
                gm.append(
                    "Appears to not be a DNA sequence.  DataType %s." % s.dataType)
                raise P4Error(gm)

        for i in range(len(self.sequences)):
            dnaSeq = theDnaSequenceList.sequences[i]
            protSeq = self.sequences[i]
            if dnaSeq.name != protSeq.name:
                gm.append(
                    "Names for (zero-based) sequence %i don't match." % i)
                gm.append("Got protein %s, DNA %s." %
                          (protSeq.name, dnaSeq.name))
                raise P4Error(gm)

        from p4.alignment import Alignment
        a = Alignment()
        a.dataType = 'dna'
        a.symbols = 'acgt'
        a.dim = 4
        a.equates = {'n': 'acgt', 'm': 'ac', 'k': 'gt',  # 'x': 'acgt',
                     'h': 'act', 'y': 'ct', 'v': 'acg',
                     'w': 'at', 'd': 'agt', 'b': 'cgt',
                     'r': 'ag', 's': 'cg'}

        for i in range(len(self.sequences)):
            dnaSeq = theDnaSequenceList.sequences[i]
            protSeq = self.sequences[i]
            s = Sequence()
            s.dataType = 'dna'
            s.name = dnaSeq.name
            s.sequence = ['---'] * self.length
            dnaPos = 0
            for j in range(self.length):
                # print "protSeq[%i] = %s" % (j, protSeq.sequence[j]),
                if protSeq.sequence[j] != '-':
                    s.sequence[j] = dnaSeq.sequence[dnaPos:dnaPos + 3]
                    dnaPos = dnaPos + 3
                # print(", codon %s" %  s.sequence[j])
            s.sequence = ''.join(s.sequence)
            # print(s.sequence)
            a.sequences.append(s)

        a.checkLengthsAndTypes()
        # a.writePhylipFile(sys.stdout)
        return a

    def setGBlocksCharSet(self, b1=None, b2=None, b3=8, b4=10, b5='n', pathToGBlocks="Gblocks", verbose=False, deleteFiles=True):
        """Find conserved regions of an alignment using Gblocks

        Gblocks <http://molevol.ibmb.csic.es/Gblocks.html>

        This method sets a charset for the gblocks-included sites.  And it
        makes a duplicate charset.  Both are put in self.nexusSets, so
        they are written if you write self in nexus format.

            'b1' Minimum number of sequences for a conserved position.
                 Default- 50% of number of sequences + 1
            'b2' Minimum number of sequences for a flank position
                 Default- 85% of the number of sequences
            'b3' Maximum Number Of Contiguous Nonconserved Positions
                 Default 8
            'b4' Minimum Length Of A Block
                 Default 10
            'b5' Allowed Gap Positions (None, With Half, All) n,h,a
                 Default 'n'

        Gblocks removes all-gap sites from the alignment (before doing
        anything, I think), and the mask that it delivers is based on the
        (possibly) shortened alignment.  The present method modifies that
        mask to make it full-length again, so that it applies to self
        correctly.

        The default b5 of n is generally too conservative, and you may
        rather want b5='h' for half, or even b5='a'.

        A duplicate charset is also made and attached, so that you can
        start with a gblocks mask and modify it from there, while keeping
        the gblocks charset.

        So you might use it like this::

          a = func.readAndPop('myAlignment.fasta')
          a.setGBlocksCharSet(b5='h')
          a.writeNexus('align_withGblocksCharset.nex')  # Nexus format to get sets block

        Then read ``align_withGblocksCharset.nex`` with Seaview, pull down
        the Sites menu, and note that there are 2 charsets -- gblocks and
        myblocks.  They are identical.  You will often find that the
        gblocks alignment is not quite what you want, but I assume you
        want to keep it intact for reference, so set and modify the
        myblocks char set in Seaview.  Then save it in Seaview as
        ``align_myblocks_handEdited.nex``.  Then, back in p4::

          a = func.readAndPop('align_myblocks_handEdited.nex')
          b = a.subsetUsingCharSet('myblocks')
          b.writePhylip('align_goodSitesOnly.phy')

        """

        # This is right from the online documentation.  Block parameter 5
        #################################################################
        # b5 -- toggles among three different possibilities for treating
        # gap positions:

        # None: no gap positions are allowed in the final alignment.  All
        # positions with a single gap or more are treated as a gap
        # position for the block selection procedure, and they and the
        # adjacent nonconserved positions are eliminated.

        # With Half: only positions where 50% or more of the sequences
        # have a gap are treated as a gap position. Thus, positions with a
        # gap in less than 50% of the sequences can be selected in the
        # final alignment if they are within an appropriate block.

        # All: all gap positions can be selected. Positions with gaps are
        # not treated differently from other positions.
        ##################################################################

        # Some of this code is originally from Cymon, but its been changed
        # considerably.

        gm = ['Alignment.setGBlocksCharSet()']
        #assert self.dataType == 'protein'
        errors = []
        # if sequenceType not in ['p', 'd', 'c']:
        #    errors.append("\tsequenceType must be either p(rotein), d(na), or c(odons)")
        # if sequenceType not in ['p']:
        #    errors.append("\tsequenceType must be p(rotein)")
        # if sequenceType == 'p':
        #    if self.dataType != 'protein':
        #        errors.append("Gblocks sequenceType is p, but this alignment is dataType %s" % self.dataType)
        if b1 is None:
            if b2 is not None:
                errors.append(
                    "\tYou must set b1 and b2 together or not at all")
                errors.append(
                    "\tb1 = Minimum number of sequences for a conserved position")
                errors.append(
                    "\tb2 = Minimum number of sequences for a flank position")
            pass
        else:
            if not isinstance(b1, int):
                errors.append(
                    "\tb1 (Minimum number of sequences for a conserved position)")
                errors.append("\tmust be None or an integer")
            else:
                if b1 > self.nTax:
                    errors.append(
                        "\tb1 (Minimum number of sequences for a conserved position)")
                    errors.append(
                        "\tmust be <= to the number of taxa in matrix")
                elif b1 < self.nTax / 2:
                    errors.append(
                        "\tb1 (Minimum number of sequences for a conserved position)")
                    errors.append(
                        "\tmust be > than the number of taxa in matrix /2")
        if b2 is None:
            pass
        else:
            if not isinstance(b2, int):
                errors.append(
                    "\tb2 (Minimum number of sequences for a flank position)")
                errors.append("\tmust be None or an integer")
            elif b2 < b1:
                errors.append(
                    "\tb2 (Minimum number of sequences for a flank position) must be >= b1")
        if not isinstance(b3, int):
            errors.append(
                "\tb3 (Maximum Number Of Contiguous Nonconserved Positions) must be an integer")
        if not isinstance(b4, int):
            errors.append(
                "\tb4 (Minimum Length Of A Block) must be an integer")
        if b5 not in ['n', 'h', 'a']:
            errors.append(
                "\tb5 (Allowed Gap Positions (None, With Half, All) n,h,a")
            errors.append("\tmust be either n(one), h(alf), a(ll)")
        if pathToGBlocks != 'Gblocks':
            if not os.path.exists(pathToGBlocks):
                errors.append(
                    "\tUnable to locate Gblocks at %s" % pathToGBlocks)
        # if b5 in ['n', 'h']:
        #    if self.allGapsSites():
        #        errors.append("\tSome site have all gaps. Deleted these before using Gblocks")

        if errors != []:
            raise P4Error(gm) + errors

        self.renameForPhylip()

        for seqNum in range(len(self.sequences)):
            seq = self.sequences[seqNum]
            seq.comment = None

        # Write the fasta file
        #fastaFileName = sha.new(str(os.getpid())).hexdigest()[-10:] + ".fasta"
        fastaFileName = 'myFaStA_oUtPuT.fasta'
        assert not os.path.isfile(fastaFileName)
        self.writeFasta(fastaFileName)

        cmd = "%s %s -t=p %s%s -b3=%i -b4=%i -b5=%s -p=s -k=y"
        cmdLine = cmd % (pathToGBlocks, fastaFileName,
                         b1 and "-b1=%i " % b1 or "", b2 and "-b2=%i" % b2 or "", b3, b4, b5)
        if verbose:
            print(cmdLine)
        if not verbose:
            cmdLine = cmdLine + " >/dev/null"
        os.system(cmdLine)

        self.restoreNamesFromRenameForPhylip()

        outputTextFileName = fastaFileName + "-gb.txts"
        outputAlignFileName = fastaFileName + "-gb"
        outputMaskFileName = fastaFileName + "-gbMask"

        try:
            fh = open(outputTextFileName, 'r')
        except IOError:
            gm.append("Unable to read output from GBlocks")
            gm.append("Check that Gblocks is in your $PATH or set 'pathToGBlocks'")
            os.remove(fastaFileName)
            raise P4Error(gm)

        if verbose:
            print(fh.read())
        fh.close()

        theAllGapsMask = self.getAllGapsMask()

        # Get the gblocks mask
        f = open(outputMaskFileName)
        fLines = f.readlines()
        f.close()
        # Find the line number of the last line starting with ">"
        for pos in range(len(fLines)):
            if fLines[pos].startswith(">"):
                spot = pos
        # make sure its the Gblocks line
        aLine = fLines[spot].rstrip()
        if not aLine.endswith("Gblocks"):
            gm.append("Something wrong with reading the mask file.  No Gblocks line.")
            raise P4Error(gm)

        # collect lines until the end of the fLines
        mStrings = []
        while 1:
            spot += 1
            try:
                aLine = fLines[spot]
            except IndexError:
                break
            mStrings.append(aLine.strip())
        mString = ''.join(mStrings)
        # print(mString)

        # The mask from gblocks is based on an alignment that has had its
        # all-gap columns removed.  We can use theAllGapsMask from self to
        # restore the gblocks mask to full length.
        myMask = []
        mStringPos = 0
        for c1 in theAllGapsMask:
            if c1 == '1':  # all sequences are gaps
                myMask.append('0')
            else:
                c2 = None
                while c2 not in ['.', '#']:
                    c2 = mString[mStringPos]
                    # print("xxx got c2='%s'" % c2)
                    mStringPos += 1
                # print("    got c2='%s'" % c2)
                if c2 == '.':
                    myMask.append('0')
                elif c2 == '#':
                    myMask.append('1')
                else:
                    gm.append(
                        "Programming error getting the mask. c2='%s'" % c2)
                    raise P4Error(gm)

        myMask = ''.join(myMask)
        # print(myMask)
        # print "mask length is %i" % len(myMask)

        if len(myMask) != self.length:
            gm.append('mask length is %i.  self.length is %i.  Bad.' %
                      (len(myMask), self.length))
            raise P4Error(gm)

        if deleteFiles:
            os.remove(fastaFileName)
            os.remove(outputTextFileName)
            os.remove(outputAlignFileName)
            os.remove(outputMaskFileName)
            if os.path.isfile("p4_renameForPhylip_dict.py"):
                os.remove("p4_renameForPhylip_dict.py")

        # return myMask
        if not self.nexusSets:
            self.setNexusSets()

        # Choose a name -- I am allowing more than one gblocks charset --
        # gblocks, gblocks_1, gblocks_2, etc
        gNames = []
        for cs in self.nexusSets.charSets:
            if cs.lowName.startswith('gblocks'):
                gNames.append(cs.lowName)
        theGName = None
        if not gNames:
            theGName = 'gblocks'
        elif len(gNames) == 1:
            theGName = 'gblocks_1'
        else:
            endBits = [int(gName[8:]) for gName in gNames if len(gName) > 7]
            endBits.sort()
            lastOne = endBits.pop()
            theGName = 'gblocks_%i' % (lastOne + 1)

        cs = CharSet(self.nexusSets)
        cs.name = theGName
        cs.lowName = theGName
        cs.num = len(self.nexusSets.charSets)
        cs.format = 'vector'
        cs.mask = myMask
        cs.aligNChar = self.length
        cs.standardize()
        self.nexusSets.charSets.append(cs)
        self.nexusSets.charSetsDict[theGName] = cs
        self.nexusSets.charSetLowNames.append(cs.lowName)
        # self.nexusSets.dump()

        myNames = []
        for cs in self.nexusSets.charSets:
            if cs.lowName.startswith('myblocks'):
                myNames.append(cs.lowName)
        theMyName = None
        if not myNames:
            theMyName = 'myblocks'
        elif len(myNames) == 1:
            theMyName = 'myblocks_1'
        else:
            endBits = [int(gName[9:]) for gName in myNames if len(gName) > 8]
            endBits.sort()
            lastOne = endBits.pop()
            theMyName = 'myblocks_%i' % (lastOne + 1)

        self.nexusSets.dupeCharSet(theGName, theMyName)


    def meanNCharsPerSite(self, includeConstantSites=True, showDistribution=True):
        """Mean number of different chars per site.

        Gaps and ambiguities are ignored.

        Constant sites can optionally be ignored.

        This is in pure Python, and can be used for Alignments with one
        part.  It is also implemented in C in the Data class, which allows
        more than one part (but no distribution).

        """

        myCDList = self.getCharDiversityOverSites()
        distro = [0] * (len(self.symbols) + 1)
        for div in myCDList:
            distro[div] += 1


        if showDistribution:
            print("\nnChars count")
            print("------ -----")
            for dk,dVal in enumerate(distro):
                print("%2i  %3i" % (dk, dVal))
            print()

        if includeConstantSites:
            dSum = 0
            nSites = 0
            for dIndx in range(1, len(distro)):   # skip all-gap sites
                #print(dIndx, distro[dIndx])
                nSites += distro[dIndx]
                dSum += dIndx * distro[dIndx]
            #print(nSites, dSum)
            return float(dSum) / nSites
        else:
            dSum = 0
            nSites = 0
            for dIndx in range(2, len(distro)): # skip all-gap sites, and constant sites
                #print(dIndx, distro[dIndx])
                nSites += distro[dIndx]
                dSum += dIndx * distro[dIndx]
            #print(nSites, dSum)
            return float(dSum) / nSites


    def getAllGapsMask(self, andMissing=True):
        """If its all gaps in a site, its 1, else 0.

        If 'andMissing' is turned on, then if a site is all gaps or '?',
        then the mask is 1.  If 'andMissing' is turned off, then it is
        strictly gaps.
        """

        m = ['0'] * self.length
        gChars = ['-']
        if andMissing:
            gChars.append('?')
        for sPos in range(self.length):
            isAllGaps = True
            for s in self.sequences:
                if s.sequence[sPos] not in gChars:
                    isAllGaps = False
                    break
            if isAllGaps:
                m[sPos] = '1'
        return ''.join(m)

    def getEnoughCharsMask(self, perCent=85):
        """Mask out sites that have too many gaps.

        This returns a mask -- a string of 1's and 0's.

        If perCent or more fraction of chars in a site are good, its a 1.
        Else 0.

        Good means not gaps, missing, or ambigs.
        Good means in self.symbols.
        """
        m = ['0'] * self.length
        threshold = perCent * 0.01 * self.nTax
        for sPos in range(self.length):
            symbCount = 0
            for s in self.sequences:
                if s.sequence[sPos] in self.symbols:
                    symbCount += 1
            if symbCount >= threshold:
                m[sPos] = '1'
        return ''.join(m)

    # These methods following are from Ababneh, Jermiin, Ma, and Robinson,
    # 2006.  Matched-pairs tests of homogeneity with applications to
    # homologous nucleotide sequences.  Bioinformatics 22:1225--1231.

    # The methods are implemented in R code in Functions.txt, available at
    # http://www.maths.usyd.edu.au/u/johnr/testsym/

    # However, that code is only for nucleotide data, and the input format
    # is restrictive.

    # Additionally, in the Ababneh et al implementation, if any site
    # contains any character (or gap) that is not one of A, C, G, or T,
    # then the entire site is ignored; its not clear why this is so.

    def simpleCharCounts(self, seqNum=None):
        """Counts of chars that are symbols, only.

        Returns a Numpy array of floats.

        Gaps or ambigs are not counted.

        If seqNum is given, then the counts are for that sequence only.
        Otherwise the counts are for the whole alignment.

        """

        scc = numpy.zeros(self.dim, dtype=numpy.double)
        # Of course we cannot test "if seqNum:", as the seqNum might be zero.
        if seqNum == None:
            for seq in self.sequences:
                for posn in range(self.nChar):
                    c = seq.sequence[posn]
                    if c in self.symbols:
                        scc[self.symbols.index(c)] += 1.
            return scc
        seq = self.sequences[seqNum]
        for posn in range(self.nChar):
            c = seq.sequence[posn]
            if c in self.symbols:
                scc[self.symbols.index(c)] += 1.
        return scc

    def getSimpleBigF(self, iA, iB):
        """Make and return a simple bigF between two sequences.

        Eg::

          seqA acgtacgt
          seqB acgtcgta

          bigF =
          [[ 1.  1.  0.  0.]
           [ 0.  1.  1.  0.]
           [ 0.  0.  1.  1.]
           [ 1.  0.  0.  1.]]


        Sites with gaps and ambigs are ignored -- thats what makes it 'simple'.

        It returns a 2-D Numpy array of floats.

        """

        bigF = numpy.zeros([self.dim, self.dim], dtype=numpy.double)
        sA = self.sequences[iA].sequence
        sB = self.sequences[iB].sequence
        for seqPos in range(self.nChar):
            cA = sA[seqPos]
            cB = sB[seqPos]
            if cA in self.symbols and cB in self.symbols:
                bigF[self.symbols.index(cA), self.symbols.index(cB)] += 1.
        return bigF

    # A function, not an Alignment method.
    def _ababnehEtAlStatsAndProbs(bigF, dim, txNumA, txNumB):
        """Re-write of Ababneh et al 2006 Testpairs function for a single pair.

        Args txNumA and txNumB are not used in the calculations -- they
        are only used to help isolate problems.

        February 2020, changed degrees of freedom for PR to match symtest.
        """

        # Calculate 3 stats -- Bowker's stat QB, Stuart's stat QS, and
        # "Internal", QR.  And calculate the probs of those stats,
        # PB, PS, PR.

        gm = ["_ababnehEtAlStatsAndProbs()"]

        # First calculate QB, Bowker's stat.  If there are any double
        # blanks (ie both bigF[i,j] = 0 and bigF[j,i] = 0) in the bigF,
        # those do not contribute to the dof.
        QB = 0.0
        zcount = 0                       # The number of double zeros in bigF
        # dof to start.  Half of off-diags of bigF
        dof = int(((dim * dim) - dim) / 2)
        for i in range(dim - 1):
            for j in range(i + 1, dim):
                # Any double zeros?
                mySum = bigF[i, j] + bigF[j, i]
                if mySum == 0:
                    zcount += 1
                else:
                    diff = (bigF[i, j] - bigF[j, i])
                    QB += (diff * diff) / mySum
        dof2 = dof - zcount
        assert dof2

        # if dof2 == 0:
        #     print(f"Between (zero-based) seqs {txNumA} and {txNumB}, got zero degrees of freedom for Bowkers.  Identical sequences?")
        #     QB = None
        #     QS = None
        #     QR = None
        #     PB = 1.0
        #     PS = None
        #     PR = None
        #     return (QB, QS, QR, PB, PS, PR)

        # Now calculate QS, Stuart's stat.  For that we need the d vector
        # and the V matrix.  Oddly, d is dim-1 long, and V is (dim-1,dim-1).

        # Lesson in sum().  sum() is a Python function, not in numpy.
        # numpy.sum() is a numpy function.  aNumpyArray.sum() is a numpy
        # method of a numpy array.

        # Ababneh does not pre-compute these values.
        sumOfColumns = bigF.sum(axis=0)
        sumOfRows = bigF.sum(axis=1)

        d = numpy.zeros((dim - 1), dtype=numpy.double)
        V = numpy.zeros((dim - 1, dim - 1), dtype=numpy.double)

        # Ababneh iterated over all i and j (one-based, i from 1:3, and j
        # from 1:3), but the result is a symmetric matrix, so half those
        # off-diag calcs are not needed.
        for i in range(dim - 1):  # yes, dim-1
            d[i] = sumOfRows[i] - sumOfColumns[i]
            for j in range(i, dim - 1):  # yes, i to dim-1
                if i == j:
                    V[i, j] = (sumOfRows[i] + sumOfColumns[i]) - (2. * bigF[i, i])
                else:
                    theItem = -1 * (bigF[i, j] + bigF[j, i])
                    V[i, j] = theItem
                    V[j, i] = theItem
        #print("d = ", d)
        #print("V = ")
        #print(V)

        try:
            Vinv = numpy.linalg.inv(V)
            QS = numpy.dot(numpy.dot(d, Vinv), d)
            QR = QB - QS
        except numpy.linalg.linalg.LinAlgError:
            if 0:
                print("Got a numpy LinAlgError while trying to invert the V matrix between (zero-based)")
                print('sequences %i and %i' % (txNumA, txNumB), end=' ')
                print("as part of the test for marginal symmetry (Stuart's test).")
            if 0:
                print("V matrix =")
                print(V)
                print()
                print('bigF =')
                print(bigF)
                print('sumOfColumns = ', sumOfColumns)
                print('sumOfRows = ', sumOfRows)
            QS = None
            QR = None
            
        #print(f"Got QS:{QS} and QR:{QR}")

        # We are finished calculating the 3 stats.  Now do the P-values
        # The dof should be an int.
        dof_S = dim - 1
        dof_R = dof2 - dof_S
        # print(f"{txNumA} {txNumB} got dof2 {dof2}, dof_S {dof_S}, and dof_R {dof_R}, QB {QB}, QS {QS}, QR {QR}")

        PB = p4.func.chiSquaredProb(QB, dof2)
        if QS is not None:
            PS = p4.func.chiSquaredProb(QS, dof_S)
            PR = p4.func.chiSquaredProb(QR, dof_R)
        else:
            PS = None
            PR = None
        return (QB, QS, QR, PB, PS, PR)

    def matchedPairsTests(self, mostSignificantOnly=False, verbose=False):
        """Get all Ababneh et al 2006 matched pairs stats and probabilies.

        Args:

            mostSignificantOnly (bool):  False by default, which gives the
              full matrices.  Setting this to True returns the three most
              significant values only, as a tuple.

        Returns:

            By default it returns six :class:`~p4.distancematrix.DistanceMatrix` 
            objects.   QB, QS, and QR (QR=Internal) matrices contain the 
            statistics, and PB, PS, and PR contain the P-values.

        For example::

          a = var.alignments[0]
          QB, QS, QR, PB, PS, PR = a.matchedPairsTests()

        The tests are pairwise on all pairs of sequences.  Note that
        it may fail to do the calculations for a pair.  If so it will
        return None for that test for that pair, and that will end up
        in DistanceMatrix objects that are returned.

        """

        QBB = DistanceMatrix()
        QSS = DistanceMatrix()
        QRR = DistanceMatrix()
        for dm in [QBB, QSS, QRR]:
            dm.setDim(self.nTax)
            dm.names = self.taxNames
        PBB = DistanceMatrix()
        PSS = DistanceMatrix()
        PRR = DistanceMatrix()
        for dm in [PBB, PSS, PRR]:
            dm.setDim(self.nTax)
            dm.names = self.taxNames

        smallestPB = None
        smallestPS = None
        smallestPR = None
        badVInverts = 0
        identicals = 0
        

        for txNumA in range(self.nTax - 1):
            for txNumB in range(txNumA + 1, self.nTax):
                # print(txNumA, txNumB)
                bigF = self.getSimpleBigF(txNumA, txNumB)
                sumAll = bigF.sum()
                sumOffDiags = sumAll - bigF.diagonal().sum()

                if 0:
                    print("bigF is\n", bigF) # it is a numpy array
                    #print("sum of columns", bigF.sum(axis=0))
                    #print("sum of rows", bigF.sum(axis=1))
                    print(f"sumOffDiags = {sumOffDiags}")

                if sumOffDiags == 0.0:
                    identicals += 1
                    QB, QS, QR, PB, PS, PR = (None, None, None, None, None, None)
                else:
                    QB, QS, QR, PB, PS, PR = _ababnehEtAlStatsAndProbs(bigF, self.dim, txNumA, txNumB)
                    # print(f"QB:{QB} QS:{QS} QR:{QR} | PB:{PB} PS:{PS} PR:{PR}")

                if PB != None:
                    if smallestPB != None:
                        if PB < smallestPB:
                            smallestPB = PB
                    else:
                        smallestPB = PB

                if PS != None:
                    if smallestPS != None:
                        if PS < smallestPS:
                            smallestPS = PS
                    else:
                        smallestPS = PS

                if PR != None:
                    if smallestPR != None:
                        if PR < smallestPR:
                            smallestPR = PR
                    else:
                        smallestPR = PR

                if QS == None:
                    badVInverts += 1

                QBB.matrix[txNumA][txNumB] = QB
                QBB.matrix[txNumB][txNumA] = QB
                QSS.matrix[txNumA][txNumB] = QS
                QSS.matrix[txNumB][txNumA] = QS
                QRR.matrix[txNumA][txNumB] = QR
                QRR.matrix[txNumB][txNumA] = QR
                PBB.matrix[txNumA][txNumB] = PB
                PBB.matrix[txNumB][txNumA] = PB
                PSS.matrix[txNumA][txNumB] = PS
                PSS.matrix[txNumB][txNumA] = PS
                PRR.matrix[txNumA][txNumB] = PR
                PRR.matrix[txNumB][txNumA] = PR

        if verbose:
            nPairs = ((self.nTax * self.nTax) - self.nTax) / 2
            print("Matched Pairs Tests: ")
            print(f"{int(nPairs)} pairs tested, ")

            if identicals or badVInverts:
                if identicals:
                    print(f"{identicals} pairs identical, ")
                if badVInverts:
                    print(f"{badVInverts} failed matrix inversions ",)
            print(f"smallest PB (Bowker's) {smallestPB} \nsmallest PS (Stuart's, Marginal) {smallestPS} \nsmallest PR (Ababneh, Internal) {smallestPR}")
        
        if mostSignificantOnly:
            return (smallestPB, smallestPS, smallestPR)
        else:
            return QBB, QSS, QRR, PBB, PSS, PRR

    def symtestAsInIQTreeNaserKhdour(self, verbose=True):
        """Matched-pairs tests of one pair, as in IQTree

        This has appeared in IQTree betas from about 1.7beta onwards,
        and is part of ``iqtree2``.  (There was a small bug in
        ``--symtest`` that was fixed in v 2.0.6, June 2020.)

        See Naser-Khdour et al GBE 2019
        https://doi.org/10.1093/gbe/evz193

        This is my attempt to replicate it.  That implementation
        chooses the sequence pair with the biggest divergence, and
        only reports stats for that pair.  This does not necessarily
        show the biggest stats with the smallest P-values.

        If it is verbose, it speaks the results to stdout.

        It returns the 3 probabilities (not stats) -- PB, PS, PR,
        which are Bowker's, Stuart's and Ababneh's.
        """

        # Make a list to hold the bigFs with the biggest divergence.
        bigFs = []
        txNumPairs = []
        biggestDivergence = 0.0

        for txNumA in range(self.nTax - 1):
            for txNumB in range(txNumA + 1, self.nTax):
                bigF = self.getSimpleBigF(txNumA, txNumB)
                sumAll = bigF.sum()
                if not sumAll:        # eg if there are no chars that line up
                    continue
                sumOffDiags = sumAll - bigF.diagonal().sum()
                divergence = sumOffDiags / sumAll
                if 0:  # extreme debug
                    print("=" * 50)
                    print(txNumA, txNumB)
                    print(bigF)
                    print(f"sumAll {sumAll}")
                    print(f"bigF.diagonal().sum() {bigF.diagonal().sum()}")
                    print(f"sumOffDiags = {sumOffDiags}")
                    print(f"divergence {divergence}")
                    print("-" * 50, "\n\n")
                if divergence < biggestDivergence:
                    continue
                elif divergence == biggestDivergence:
                    bigFs.append(bigF)
                    txNumPairs.append((txNumA, txNumB))
                else:
                    # its bigger, so wipe previous results
                    bigFs = [bigF]
                    txNumPairs = [(txNumA, txNumB)]
                    biggestDivergence = divergence
                    
        if verbose:
            print("divergence matrices (F-matrix):")
            print(bigFs)
            print("biggest divergences")
            print(biggestDivergence)
            print("... between these (zero-based) sequence numbers:")
            print(txNumPairs)
        assert biggestDivergence > 0.0, "Got zero divergence between the sequences"

        myBigF = None
        myTxNumPair = None
        if len(bigFs) == 1:
            myBigF = bigFs[0]
            myTxNumPair = txNumPairs[0]
        else:
            print(f"Got {len(bigFs)} pairs with the same divergence.  Choosing randomly.", file=sys.stderr)
            myIndex = random.randrange(len(bigFs))
            myBigF = bigFs[myIndex]
            myTxNumPair = txNumPairs[myIndex]
        
        # We want all the pairwise p-vals.  These PB, PS, and PR are DistanceMatrix objects
        QB, QS, QR, PB, PS, PR = self.matchedPairsTests()
        for PM in [PB,PS,PR]:
            PM.ft = PM.flatTriangle()
            PM.nSig = len([it for it in PM.ft if it is not None  and it <= 0.05])
            PM.nNonSig = len([it for it in PM.ft if it is not None and it > 0.05])
            PM.forBiggestDivergence = PM.matrix[myTxNumPair[0]][myTxNumPair[1]]

        if verbose:
            print(f"Bowker's test of symmetry")
            print(f"  Significant:     {PB.nSig}")
            print(f"  Non-significant: {PB.nNonSig}")
            print(f"  P-value for most-diverged pair: {PB.forBiggestDivergence}")
            print(f"Stuart's test of marginal symmetry")
            print(f"  Significant:     {PS.nSig}")
            print(f"  Non-significant: {PS.nNonSig}")
            print(f"  P-value for most-diverged pair: {PS.forBiggestDivergence}")
            print(f"Ababneh's test of internal symmetry")
            print(f"  Significant:     {PR.nSig}")
            print(f"  Non-significant: {PR.nNonSig}")
            print(f"  P-value for most-diverged pair: {PR.forBiggestDivergence}")
            
        return PB.forBiggestDivergence, PS.forBiggestDivergence, PR.forBiggestDivergence
        
        

    def testOverallFromAbabnehEtAl2006(self):
        """Marginal symmetry (Stuarts's) test for more than two matched sequences

        As in Ababneh et al 2006 Page 1227, calculated as in the R function Testoverall.

        A tuple composed of the Ts stat, the degrees of freedom, and the probability, is returned.
        """

        # Make the J matrix.
        J = numpy.zeros(
            (self.dim * (self.nTax - 1), self.dim * (self.nTax - 1)), dtype=numpy.double)
        for seqNum in range(self.nTax - 1):
            start = seqNum * self.dim
            for i in range(self.dim):
                for j in range(self.dim):
                    J[start + i, start + j] = 1.
        # print J

        # Make the L matrix.
        L = numpy.zeros(
            (self.dim * (self.nTax - 1), self.dim * self.nTax), dtype=numpy.double)
        for seqNum in range(self.nTax - 1):
            start = self.dim * seqNum
            for i in range(self.dim):
                L[start + i, i] = 1.
        for seqNum in range(self.nTax - 1):
            start = self.dim * seqNum
            for i in range(self.dim):
                L[start + i, start + self.dim + i] = -1
        # print L

        # Make the v matrix, and the n vector
        v = numpy.zeros(
            (self.dim * self.nTax, self.dim * self.nTax), dtype=numpy.double)
        n = numpy.zeros(self.dim * self.nTax, dtype=numpy.double)

        # This is inefficient!
        for seqNumA in range(self.nTax):
            for seqNumB in range(self.nTax):
                if seqNumA == seqNumB:
                    scc = self.simpleCharCounts(seqNum=seqNumA)
                    fi = scc / scc.sum()
                    start = seqNumA * self.dim
                    for i in range(self.dim):
                        v[start + i, start + i] = fi[i]
                        n[start + i] = fi[i]
                else:
                    fij = self.getSimpleBigF(seqNumA, seqNumB)
                    fij /= fij.sum()
                    startA = seqNumA * self.dim
                    startB = seqNumB * self.dim
                    for i in range(self.dim):
                        for j in range(self.dim):
                            pass
                            v[startA + i, startB + j] = fij[i, j]
        # print v
        # print n
        firstBit = numpy.dot(L, n)
        myIdentity = numpy.identity(self.dim * (self.nTax - 1), dtype=numpy.double)
        secondBit = numpy.linalg.solve(
            (numpy.dot(numpy.dot(L, v), L.transpose()) + J), myIdentity)
        # print secondBit
        Ts = numpy.dot(numpy.dot(firstBit, secondBit), firstBit) * self.nChar
        # print Ts
        df = (self.dim - 1) * (self.nTax - 1)
        pval = p4.func.chiSquaredProb(Ts, df)
        return Ts, df, pval

    def getMinmaxChiSqGroups(self, percent_cutoff=0.05, min_bins=2, max_bins=20,
                             n_choices=1000, seed=42, verbose=False):
        """An interface for Susko and Roger's minmax-chisq program.

        Susko, E. and Roger, A.J. (2007). On reduced amino acid alphabets for
        phylogenetic inference.  Mol. Biol. Evol. 24:2139-2150.

        * *percent_cutoff* - P-value percentage above which chi-squared test is considered homogeneous. Default is 0.05.
        * *min_bins* - The minimum number of bins to consider. The default value is 2.
        * *max_bins* - The maximum number of bins to consider. The default value is 20.
        * *n_choices* - number. The number of random choices of bins to consider for each bin size. The default is 1000.
        * *seed* - An integer giving the starting seed for any random generation done by the program. By default this is 42.

        All in memory -- no files written.

        This is a very basic interface.  To be able to do multiple runs 
        with different random seeds or to specify a desired number of 
        bins requires re-jigging.  I don't use this method anymore.

        """
        gm = ["Alignment.getMinmaxChiSqGroups()"]

        # Unfortunately minmax-chsq only accepts 10 char taxon names and only accepts
        # data in CAPITALS!?!
        if self.dataType != "protein":
            gm.append("The data are not protein.")
            raise P4Error(gm)
        maxNameLen = max([len(s.name) for s in self.sequences])
        if maxNameLen > var.phylipDataMaxNameLength:
            gm.append(
                'The longest name length in this alignment is %i' % maxNameLen)
            arg = 'var.phylipDataMaxNameLength is now %i' % var.phylipDataMaxNameLength
            gm.append(arg)
            gm.append('Sequence names will be truncated.  Fix it.')
            gm.append("You may want to use the 'renameForPhylip()' method.")
            raise P4Error(gm)
        if not p4.func.which2("minmax-chisq"):
            gm.append("minmax-chisq is not in your path")
            raise P4Error(gm)
        the_data = []
        the_data.append("%i %i" % (self.nTax, self.nChar))
        theFormat = "%-10s"
        for s in self.sequences:
            the_data.append("%-10s %s" % (s.name, s.sequence.upper()))
        data_file = "\n".join(the_data)
        cmd = "minmax-chisq -l%i -u%i -n%i -s%i" % (
            min_bins, max_bins, n_choices, seed)
        child = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, shell=True)
        stdout, stderr = child.communicate(data_file)
        if stderr != "":
            gm.append("Minmax-chisq returned an error: %s" % stderr)
            raise P4Error(gm)
        results = zip(*[iter(stdout.split("\n"))] * 2)
        # Need the binning before pvalue falls below percent_cutoff
        # They could be all <= 0.05 or all >= 0.05 because of the binning range
        pvalues = [float(bin[0].split()[1]) for bin in results]
        if not pvalues[0] >= percent_cutoff:
            print("No p-value <= %s" % percent_cutoff)
            return None
        if not any(pv for pv in pvalues if pv >= percent_cutoff):
            print("No p-value >= %s" % percent_cutoff)
            return None
        for i, bin_result in enumerate(results):
            nbins, pvalue = bin_result[0].split()
            if float(pvalue) <= percent_cutoff:
                opt_bin = results[i - 1]
                break
        nbins, pvalue = opt_bin[0].split()
        scores = opt_bin[1].split()
        amino_acid_order = "A R N D C Q E G H I L K M F P S T W Y V".lower().split()
        c = zip(scores, amino_acid_order)
        groups = {}
        for bin, amino in c:
            sbin = str(bin)
            groups[sbin] = groups.get(sbin, "") + amino
        if verbose:
            print("\nminmax-chisq output:\n%s" % stdout)
            print("\nMaximum number of bins that maintains homogeneity: %s" % nbins)
            print("\nGroups: %s" % ", ".join(groups.values()))

        return groups.values()

    def getKosiolAISGroups(self, tree, n_bins, remove_files=False, verbose=True):
        """An interface for Kosiol's program AIS, for grouping amino acids.

        Writes three files: equi, q, and evec that are required as input to the
        software AIS (Almost Invariant Sets).

        - *tree* = tree object with optimised model compatible with the alignment
        - *nBins* = the number of groups required

        Kosiol, C., Goldman, N. and Buttimore, N. (2004). A new criterion and method
        for amino acid classification. J Theo Biol 228: 97-106.

        From share/examples/kosiol_ais.html
        This writes files.
        """
        gm = ["Alignment.writeKosiolAISFiles() "]
        if not tree.model:
            gm.append("Requires a tree with an optimimised model attached.")
            raise P4Error(gm)
        if not tree.data:
            d = Data()
            tree.data = d
        if self.dataType != "protein":
            gm.append("The data are not protein.")
            raise P4Error(gm)
        if not p4.func.which2("ais"):
            gm.append("ais is not in your path")
            raise P4Error(gm)
        if not n_bins > 1 or not n_bins < 20:
            gm.append("n_bins must be > 1 and < 20")
            raise P4Error(gm)
        # Init the model
        tree.calcLogLike(verbose=0)
        # Write the aa freqs
        f = open('equi', 'w')
        f.write('20\n')
        for i in range(20):
            f.write("%f\n" % tree.model.parts[0].comps[0].val[i])
        f.close()

        bigQ = tree.model.getBigQ()

        # Write the bigQ
        f = open('q', 'w')
        f.write('20\n')
        for i in range(20):
            for j in range(20):
                f.write("%5g  " % bigQ[i][j])
            f.write('\n')
        f.close()

        # Get the eigensystem
        evals, evecs = numpy.linalg.eig(bigQ)

        # According to the web page, "The right eigenvectors should be ordered
        # according to the absolute value of their eigenvalues."  Well, the
        # output from numpy, which uses lapack, is not so ordered.  So do it.
        sorter = numpy.argsort(evals)
        sorter = sorter[::-1]   # reverse
        # print sorter

        f = open('evec', 'w')
        f.write('20\n')
        f.write('20\n')
        for colNum in sorter:
            for rowNum in range(20):
                f.write("%5g\n" % evecs[rowNum][colNum])
        f.close()
        commands = b"equi\nq\n\nevec\n%i" % n_bins
        n = 0
        while n < 10:
            # Sometimes ais fails with the error "Warning: Eigenvectors are
            # perturbed etc" sometime it doesnt we're going to try a maximum of 10
            # times
            child = subprocess.Popen("ais", stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True)
            stdout, stderr = child.communicate(commands)
            #returncode = child.returncode
            if child.returncode != 0:
                n += 1
                continue
            else:
                if verbose:
                    print(stdout)
                fh = open("%isets.txt" % n_bins, "r")
                lines = fh.readlines()
                fh.close()
                sets = []
                for line in lines:
                    if line.startswith(" Set"):
                        aas = line.split("{")[1].split("}")[0]
                        sets.append("".join(aas.strip().lower().split()))
                if len(sets) != n_bins:
                    gm.append("Error: only recovered %i sets" % len(sets))
                    raise P4Error(gm)
                if remove_files:
                    for f in ["equi", "q", "evec", "%isets.txt" % n_bins,
                              "graph.txt", "coloring.txt"]:
                        if os.path.exists(f):
                            os.remove(f)
                return sets
        print("Unable to run ais after 10 attempts. Giving up.")

    def mrpSlice(self, pos, zeroBasedNumbering=True):
        """Pretty-print a mrp site, with no '?' positions.

        Zero-based numbering, unless arg zeroBasedNumbering is set to False.
        """
        if zeroBasedNumbering:
            ss = self.sequenceSlice(pos)
            print("mrp matrix position (zero-based) %i" % pos)
        else:
            pos2 = pos - 1
            ss = self.sequenceSlice(pos2)
            print("mrp matrix position (1-based) %i" % pos2)
        for txNum in range(self.nTax):
            if ss[txNum] == '?':
                pass
            else:
                print("  %25s %s" % (self.taxNames[txNum], ss[txNum]))

    #################################################### read-write
    ####################################################

    def readOpenPhylipFile(self, flob, nTax, nChar):
        """Read flob to get data in phylip format.

        The user would generally not need to call this method directly.
        It is called by read() etc.  """

        gm = ['Alignment.readOpenPhylipFile()']

        dbug = False

        if dbug:
            print("\nreadOpenPhylipFile here")
            if hasattr(flob, 'name'):
                print("    fileName is %s" % flob.name)
            print("    nTax is", nTax)
            print("    nChar is", nChar)
            print("    var.phylipDataMaxNameLength is", var.phylipDataMaxNameLength)

        # It is difficult to tell whether it is sequential or interleaved.
        # For example, is the following alignment sequential or interleaved?

        # 2 14
        # acctg     aaaa
        # gaattc    aaaa
        # gaattc    cccc
        # cccccc    cccc

        # It might be sequential:
        # acctg    aaaagaattcaaaa
        # gaattc    cccccccccccccc

        # or it might be interleaved:
        # acctg    aaaagaattccccc
        # gaattc    aaaacccccccccc

        # And it is impossible to tell!

        # Here is the top of a file that I had to read.

        # 36  209
        # Dvi        flnsfnakleqpvrqhlknvyaclamstmsaalgaaagflsaigalvfff
        # Gsp        finsfnskleqpvrqhlknvyacltmatmaaavgasagflsgigalvffg
        # Tca        flnsfsnsleapvrqhlknvyaclamstmaaaigasagflsgigaliffg
        # Ban        finsfqnrlespvrqhlknvygtlmmtcgaasagvyagilsaiagaalml
        # Bmo        fvnsfqnrleppvrqhlknvyatlmmtcvsasagvyagflsaivgaglml

        # After reading in Dvi, and then the next 3 lines, I had a
        # sequence name (Dvi) and exactly 209 characters, all of which
        # were aa characters.  So it was decided that it was Sequential!
        # Wrong!

        # Slurp the file into a list of lines, except for the header line.
        # Skip blank lines.
        # If a line startswith '//', stop (its a paml file)
        # If a line is composed wholly of numerals, don't collect it (paml again)
        # If a line starts with  '[', stop collecting until after a line that
        # starts with ']' (paml again)

        # Sometimes paml files have stuff at the end, without any marking.
        # So later we need to stop collecting sequences after nTax.

        # First, get all the lines in a list.  The headerLine is separate.
        # Blank lines and lines that have no potential characters are ignored.
        myPotentialChars = string.ascii_letters + '?-'

        flob.seek(0, 0)
        headerLine = None
        theLines = []
        inComment = False
        while 1:
            aLine = flob.readline()
            if not aLine:
                break
            aLine = aLine.strip()
            if not aLine:  # blank line
                pass
            elif aLine.startswith("//"):
                break
            else:
                if not headerLine:
                    headerLine = aLine
                else:
                    if aLine.startswith('['):
                        inComment = True
                    elif inComment:
                        if aLine.startswith(']'):
                            inComment = False
                    else:
                        # Check if it has letters.
                        hasLetters = False
                        for c in aLine:
                            if c in myPotentialChars:
                                hasLetters = True
                                break
                        if hasLetters:
                            theLines.append(aLine)
        assert theLines

        # Now we have the headerLine separate, and all the non-blank lines
        # in theLines list.  Check the headerLine.
        # If a paml file is interleaved, it seems to have an I in the header
        # line.
        splFirstLine = headerLine.split()
        assert len(splFirstLine) >= 2   # More than 2 sometimes
        try:
            firstInt = int(splFirstLine[0])
            secondInt = int(splFirstLine[1])
        # This should have been caught before, so this should never happen.
        except ValueError:
            gm.append('bad first line %s' % headerLine)
            raise P4Error(gm)
        assert firstInt == nTax
        assert secondInt == nChar

        # Look for an I as the third symbol in the header line.
        # Not used.
        #gotPamlInterleavedI = False
        # if len(splFirstLine) > 2:
        #    thirdBit = splFirstLine[2]
        #    if thirdBit == 'I':
        #        gotPamlInterleavedI = True

        # There are 2 conflicting strategies.  The strict phylip file has
        # a set maximum tax name size, and all names must fit in that size
        # -- usually 10 (but it is a compile-time option).  If the tax
        # name is shorter, it must make up for it with blanks.  And there
        # is no requirement for a space between the tax name and the
        # sequence.  The tax name may have spaces.

        # The other strategy sets off the tax name from the sequence with
        # a space, and so that space is needed.  That implies, I think,
        # that the taxname may not have spaces.  The taxname is often
        # longer than 10.  There is no requirement for the taxname to fill
        # a certain width, and so short names do not need to be filled
        # with multiple blanks -- a single blank will do for both long and
        # short tax names.

        # I am calling the former 'strict' and the latter
        # 'whitespaceSeparatesNames'

        # First I make a half-assed attempt to determine if it is
        # sequential or interleaved.  The decision could easily be wrong,
        # and thats ok.

        # The first line might be on a name by itself, which implies
        # sequential.  If it has no spaces, then that will be assumed.

        # The first line might be something like
        # SomeLongName    acgtacgtacgt
        # which probably means interleaved, with a space.  Not strict
        # phylip format -- the name is longer than 10.

        # The first line might (rarely) be
        # SomeNameABacgtacgtacgtacgt
        # which probably means interleaved, proper phylip format

        # Or the first line might be
        # Short Nameacgtacgtacgt
        # which again implies proper phylip format.

        isSequential = None  # undecided to start.
        moduloRemainderIsZero = None
        whitespaceSeparatesNames = True

        # Check whether the number of lines is some multiple of nTax -- if
        # so then it is probably interleaved.  If not, it cannot be
        # interleaved.
        remainder = len(theLines) % nTax
        if remainder == 0:
            moduloRemainderIsZero = True
        else:
            moduloRemainderIsZero = False
            isSequential = True

        if dbug:
            print("    moduloRemainderIsZero is %s" % moduloRemainderIsZero)

        if isSequential == None:
            # Look at the first line of stuff, after the numbers.
            # If the total length of the line is var.phylipDataMaxNameLength
            # (usually 10) or less, then I assume it is all name.
            firstLine = theLines[0]
            if len(firstLine) <= (var.phylipDataMaxNameLength):
                if dbug:
                    print("Got the name %s by virtue of it being a short name on the first data line." % firstLine)
                    print("Setting 'isSequential' to True, based on that.")
                isSequential = True
            elif not firstLine.count(' '):
                isSequential = True

        if isSequential == None:
            isSequential = True
        haveTriedSequential_strict = False
        haveTriedSequential = False
        haveTriedInterleaved_strict = False
        haveTriedInterleaved = False
        # (None, True, or False, depending on don't know, success, or failure)
        sequentialResult = None
        interleavedResult = None
        gotIt = False

        #######################################################################
        #######################################################################
        while not gotIt:
            if dbug:
                if hasattr(flob, 'name'):
                    print("  fileName is %s" % flob.name)
                print("  isSequential is %s" % isSequential)
                print("  whitespaceSeparatesNames is %s" % whitespaceSeparatesNames)
                print("  haveTriedSequential = %s" % haveTriedSequential)
                print("  haveTriedSequential_strict = %s" % haveTriedSequential_strict)
                print("  seqentialResult is %s" % sequentialResult)
                print("  haveTriedInterleaved is %s" % haveTriedInterleaved)
                print("  haveTriedInterleaved_strict is %s" % haveTriedInterleaved_strict)
                print("  interleavedResult is %s" % interleavedResult)

            if 0:
                if len(theLines) >= 5:
                    theRange = range(5)
                else:
                    theRange = range(len(theLines))
                for lineNum in theRange:
                    print(theLines[lineNum])

            if whitespaceSeparatesNames:
                if isSequential and not haveTriedSequential:
                    ret = self._readPhylipSequential(nTax, nChar, theLines)
                    # print "a ret = %s" % ret
                    if ret:
                        gotIt = True
                    else:
                        haveTriedSequential = True
                        self.sequences = []
                        if not haveTriedInterleaved:
                            isSequential = False
                        else:
                            if not haveTriedSequential_strict:
                                whitespaceSeparatesNames = False
                            elif not haveTriedInterleaved_strict:
                                whitespaceSeparatesNames = False
                                isSequential = False
                elif not isSequential and not haveTriedInterleaved:
                    ret = self._readPhylipInterleaved(nTax, nChar, theLines)
                    # print "b ret = %s" % ret
                    if ret:
                        gotIt = True
                    else:
                        haveTriedInterleaved = True
                        self.sequences = []
                        if not haveTriedSequential:
                            isSequential = True
                        else:
                            if not haveTriedSequential_strict:
                                whitespaceSeparatesNames = False
                            elif not haveTriedInterleaved_strict:
                                whitespaceSeparatesNames = False
                                isSequential = False
            else:
                if isSequential:
                    ret = self._readPhylipSequentialStrict(
                        nTax, nChar, theLines)
                    # print "c ret = %s" % ret
                    if ret:
                        gotIt = True
                    else:
                        haveTriedSequential_strict = True
                        self.sequences = []
                        if not haveTriedInterleaved_strict:
                            isSequential = True
                        else:
                            if not haveTriedSequential:
                                whitespaceSeparatesNames = True
                            elif not haveTriedInterleaved:
                                whitespaceSeparatesNames = True
                                isSequential = False
                else:
                    ret = self._readPhylipInterleavedStrict(
                        nTax, nChar, theLines)
                    # print "d ret = %s" % ret
                    if ret:
                        gotIt = True
                    else:
                        haveTriedInterleaved_strict = True
                        self.sequences = []
                        if not haveTriedSequential_strict:
                            isSequential = True
                        else:
                            if not haveTriedSequential:
                                whitespaceSeparatesNames = True
                            elif not haveTriedInterleaved:
                                whitespaceSeparatesNames = True
                                isSequential = False
            # print "x gotIt is now %s" % gotIt
            if gotIt:
                break
            if haveTriedSequential and haveTriedInterleaved and haveTriedSequential_strict and haveTriedInterleaved_strict:
                gm.append("Failed to read the phylip or phylip-like file.")
                if not var.verboseRead:
                    gm.append("(For more info, turn var.verboseRead on)")
                raise P4Error(gm)
        #######################################################################
        #######################################################################

        # Paml interleaved files sometimes have dots after the first sequence.
        firstSequence = self.sequences[0]
        for s in self.sequences[1:]:
            if s.sequence.count('.'):
                s.sequence = list(s.sequence)
                for seqPos in range(nChar):
                    if s.sequence[seqPos] == '.':
                        s.sequence[seqPos] = firstSequence.sequence[seqPos]
                s.sequence = ''.join(s.sequence)

        # Set the dataType
        for s in self.sequences:
            ret = None
            # returns 1,2 or 0, respectively
            ret = p4.func.isDnaRnaOrProtein(s.sequence)
            if ret == 1:
                s.dataType = 'dna'
                s.symbols = 'acgt'
            elif ret == 0:
                s.dataType = 'protein'
                s.symbols = 'arndcqeghilkmfpstwyv'
            else:
                raise P4Error("Got rna sequence.  Fix me.")

        # Having read in all the sequences, check for valid characters.
        if len(self.sequences) > 0:
            bads = 0
            if self.sequences[0].dataType == 'dna':
                for s in self.sequences:
                    j = 0
                    while j < nChar:
                        if s.sequence[j] not in var.validDnaChars:
                            print("Got bad character '%s' in (zero-based) dna sequence %s " % \
                                (s.sequence[j], self.sequences.index(s)))
                            print("          sequence name: %s" % s.name)
                            print("          at (zero-based) position %s" % j)
                            bads = bads + 1
                            if bads > 10:
                                print("...and possibly others")
                                break
                        j = j + 1
                    if bads > 10:
                        break
                if bads:
                    gm.append("See the list of bad chars above")
                    raise P4Error(gm)
            elif self.sequences[0].dataType == 'protein':
                for s in self.sequences:
                    j = 0
                    while j < nChar:
                        if s.sequence[j] not in var.validProteinChars:
                            print("Got bad character '%s' in (zero-based) protein sequence %s " % \
                                (s.sequence[j], self.sequences.index(s)))
                            print("          sequence name: %s" % s.name)
                            print("          at (zero-based) position %s" % j)
                            bads = bads + 1
                            if bads > 10:
                                print("...and possibly others")
                                break
                        j = j + 1
                    if bads > 10:
                        break
                if bads:
                    gm.append("See the list of bad chars above")
                    raise P4Error(gm)

        # for s in self.sequences:
        #    print s.name
        #    print s.dataType
        # sys.exit()

    def _readPhylipSequential(self, nTax, nChar, theLines):
        
        if var.verboseRead:
            print("Attempting to read the phylip file assuming")
            print("  sequential format, and that whitespaceSeparatesNames")

        lineNum = 0

        for taxNum in range(nTax):
            theSequenceBits = []
            seqLenSoFar = 0
            aLine = theLines[lineNum]
            lineNum += 1
            splitLine = aLine.split()
            theName = splitLine[0]
            if len(splitLine) > 1:
                for aBit in splitLine[1:]:
                    bBit = p4.func.stringZapWhitespaceAndDigits(aBit)
                    theSequenceBits.append(bBit)
                    seqLenSoFar += len(bBit)
            while seqLenSoFar < nChar:
                try:
                    aLine = theLines[lineNum]
                    lineNum += 1
                except IndexError:
                    if var.verboseRead:
                        print("    Early termination")
                    return False
                theSequenceBits.append(
                    p4.func.stringZapWhitespaceAndDigits(aLine))
                seqLenSoFar += len(theSequenceBits[-1])

            if seqLenSoFar != nChar:
                if var.verboseRead:
                    print("    Did not get exactly nChar.")
                    print("    Expected: ", seqLenSoFar)
                    print("    Got: ", nChar)
                return False

            s = Sequence()
            s.name = theName
            s.sequence = ''.join(theSequenceBits)
            s.sequence = s.sequence.lower()
            self.sequences.append(s)

        if len(self.sequences) != nTax:
            if var.verboseRead:
                print("    Did not get correct nTax.")
            return False

        isBad = False
        for s in self.sequences:
            if len(s.sequence) != nChar:
                isBad = True
                break
        if isBad:
            if var.verboseRead:
                print("    At least one sequence has the wrong nChar.")
            return False

        if var.verboseRead:
            print("    Got a candidate alignment.")
        return True

    def _readPhylipInterleaved(self, nTax, nChar, theLines):
        
        if var.verboseRead:
            print("Attempting to read the phylip file assuming")
            print("  interleaved format, and that whitespaceSeparatesNames")

        lineNum = 0
        for i in range(nTax):
            s = Sequence()
            aLine = theLines[lineNum]
            lineNum += 1
            splitLine = aLine.split()
            s.name = splitLine[0]
            s.theSequenceBits = []
            s.seqLenSoFar = 0
            if len(splitLine) > 1:
                for aBit in splitLine[1:]:
                    bBit = p4.func.stringZapWhitespaceAndDigits(aBit)
                    s.theSequenceBits.append(bBit)
                    s.seqLenSoFar += len(bBit)
            # print "%15s %s" % (s.name, bBit)
            self.sequences.append(s)
        if len(self.sequences) != nTax:
            if var.verboseRead:
                print("    Did not get exactly nTax sequences")
            return False

        # Subsequent cycles don't have names-- just sequence
        seqNum = 0
        isDone = False
        while not isDone:
            try:
                aLine = theLines[lineNum]
                lineNum += 1
            except IndexError:
                break
            # print "aLine a: %s" % aLine
            segment = p4.func.stringZapWhitespaceAndDigits(aLine)
            if len(segment):
                self.sequences[seqNum].theSequenceBits.append(segment)
                self.sequences[seqNum].seqLenSoFar += len(segment)
                seqNum += 1
                if seqNum == nTax:
                    seqNum = 0
                    # print "aLine b: %s" % aLine
                    # print "=" * 50
                    # Finished a cycle.  If we have enough sequence, break
                    if self.sequences[seqNum].seqLenSoFar >= nChar:
                        isDone = True
            else:
                print('_readPhylipInterleaved() bad line?: %s' % aLine)

        for s in self.sequences:
            s.sequence = ''.join(s.theSequenceBits)
            s.sequence = s.sequence.lower()
            del s.theSequenceBits
            del s.seqLenSoFar

        if len(self.sequences) != nTax:
            if var.verboseRead:
                print("    Did not get correct nTax.")
            return False

        isBad = False
        for s in self.sequences:
            # s.write()
            # print len(s.sequence), nChar
            if len(s.sequence) != nChar:
                isBad = True
                break
        if isBad:
            if var.verboseRead:
                print("    At least one sequence has the wrong nChar.")
                print("    Expected: ", nChar)
                print("    Got: ", len(s.sequence))
            return False

        if var.verboseRead:
            print("    Got a candidate alignment.")
        return True

    def _readPhylipSequentialStrict(self, nTax, nChar, theLines):
        
        if var.verboseRead:
            print("Attempting to read the phylip file assuming")
            print("  sequential format, with strict phylip format")

        lineNum = 0

        for taxNum in range(nTax):
            theSequenceBits = []
            seqLenSoFar = 0
            aLine = theLines[lineNum]
            lineNum += 1
            theName = aLine[:var.phylipDataMaxNameLength].strip()
            aBit = p4.func.stringZapWhitespaceAndDigits(
                aLine[var.phylipDataMaxNameLength:])
            theSequenceBits.append(aBit)
            seqLenSoFar = len(aBit)
            while seqLenSoFar < nChar:
                try:
                    aLine = theLines[lineNum]
                    lineNum += 1
                except IndexError:
                    if var.verboseRead:
                        print("    Early termination")
                    return False
                aBit = p4.func.stringZapWhitespaceAndDigits(aLine)
                theSequenceBits.append(aBit)
                seqLenSoFar += len(aBit)

            if seqLenSoFar != nChar:
                if var.verboseRead:
                    print("    Did not get exactly nChar.")
                return False

            s = Sequence()
            s.name = theName
            s.sequence = ''.join(theSequenceBits)
            s.sequence = s.sequence.lower()
            self.sequences.append(s)

        if len(self.sequences) != nTax:
            if var.verboseRead:
                print("    Did not get correct nTax.")
            return False

        isBad = False
        for s in self.sequences:
            if len(s.sequence) != nChar:
                isBad = True
                break
        if isBad:
            if var.verboseRead:
                print("    At least one sequence has the wrong nChar.")
            return False

        if var.verboseRead:
            print("    Got a candidate alignment.")
        return True

    def _readPhylipInterleavedStrict(self, nTax, nChar, theLines):
        if var.verboseRead:
            print("Attempting to read the phylip file assuming")
            print("  interleaved format, with strict phylip format")

        lineNum = 0
        for i in range(nTax):
            s = Sequence()
            aLine = theLines[lineNum]
            lineNum += 1

            s.name = aLine[:var.phylipDataMaxNameLength].strip()
            s.theSequenceBits = []
            aBit = p4.func.stringZapWhitespaceAndDigits(
                aLine[var.phylipDataMaxNameLength:])
            s.theSequenceBits.append(aBit)
            s.seqLenSoFar = len(aBit)
            self.sequences.append(s)
        if len(self.sequences) != nTax:
            if var.verboseRead:
                print("    Did not get exactly nTax sequences")
            return False

        # Subsequent cycles don't have names-- just sequence
        seqNum = 0
        isDone = False
        while not isDone:
            try:
                aLine = theLines[lineNum]
                lineNum += 1
            except IndexError:
                break
            # print "aLine a: %s" % aLine
            segment = p4.func.stringZapWhitespaceAndDigits(aLine)
            if len(segment):
                self.sequences[seqNum].theSequenceBits.append(segment)
                self.sequences[seqNum].seqLenSoFar += len(segment)
                seqNum += 1
                if seqNum == nTax:
                    seqNum = 0
                    # Finished a cycle.  If we have enough sequence, break
                    if self.sequences[seqNum].seqLenSoFar >= nChar:
                        isDone = True
            else:
                print('_readPhylipInterleavedStrict() bad line?: %s' % aLine)

        for s in self.sequences:
            s.sequence = ''.join(s.theSequenceBits)
            s.sequence = s.sequence.lower()
            del s.theSequenceBits
            del s.seqLenSoFar

        if len(self.sequences) != nTax:
            if var.verboseRead:
                print("    Did not get correct nTax.")
            return False

        isBad = False
        for s in self.sequences:
            if len(s.sequence) != nChar:
                isBad = True
                break
        if isBad:
            if var.verboseRead:
                print("    At least one sequence has the wrong nChar.")
            return False

        if var.verboseRead:
            print("    Got a candidate alignment.")
        return True

    def _readOpenClustalwFile(self, flob):
        
        gm = ['Alignment._readOpenClustalwFile()']

        dbug = 0
        if dbug:
            print("\n_readOpenClustalwFile() here")

        # readline up to the first line of sequence
        aLine = flob.readline()
        if not aLine:
            gm.append("No sequence?")
            raise P4Error(gm)
        while len(aLine) <= 1:
            aLine = flob.readline()
            if not aLine:
                gm.append("No sequence?")
                raise P4Error(gm)
        if dbug:
            print("a- got aLine: '%s'" % aLine)

        # Do the first cycle:
        while len(aLine) > 1 and aLine[0] not in string.whitespace:
            s = Sequence()
            splitLine = aLine.split()
            if len(splitLine) != 2:
                gm.append("Odd line:\n%s" % aLine)
                raise P4Error(gm)
            s.name = splitLine[0]
            s.temp = []
            s.temp.append(splitLine[1].strip().lower())
            if dbug:
                print("b got name: %s" % s.name)
                print("b got a line of sequence:\n          %s" % s.temp)
            self.sequences.append(s)
            aLine = flob.readline()
            if not aLine:
                break

        #nSeqs = len(self.sequences)

        nSeqs = len(self.sequences)
        if not nSeqs:
            gm.append("No sequences?")
            raise P4Error(gm)
        if dbug:
            print("Got %i sequences (after first cycle)" % nSeqs)

        # Do subsequent cycles
        while 1:
            while len(aLine) <= 1 or aLine[0] in string.whitespace:
                aLine = flob.readline()
                if not aLine:
                    break
            if aLine:
                if dbug:
                    print("begin subsequent cycle with aLine:\n    %s" % aLine)
                for s in self.sequences:
                    if aLine[0] in string.whitespace or len(aLine) <= 1:
                        gm.append("Bad line:\n    %s" % aLine)
                        #raise P4Error(gm)
                    splitLine = aLine.split()
                    if len(splitLine) != 2:
                        gm.append("Odd line:\n%s" % aLine)
                        raise P4Error(gm)
                    if s.name != splitLine[0]:
                        gm.append(
                            "Problem: existing name %s does not match new name %s" % (s.name, splitLine[0]))
                        raise P4Error(gm)
                    s.temp.append(splitLine[1].strip().lower())
                    if dbug:
                        print("got name: %s" % s.name)
                        print("got a line of sequence:\n          %s" % s.temp)
                    aLine = flob.readline()
                    if not aLine:
                        break
            else:
                break

        # sys.exit()

        nChar = len(''.join(self.sequences[0].temp))

        for s in self.sequences:
            if dbug:
                print("got name: %s" % s.name)
                print("got sequence:\n          %s" % s.temp)
            s.temp = ''.join(s.temp)
            if len(s.temp) != nChar:
                gm.append(
                    "Something is wrong with the sequence length of the clustalw file.")
                gm.append("(zero-based) sequence %i" % self.sequences.index(s))
                gm.append("expected %i, got %i" % (nChar, len(s.temp)))
                raise P4Error(gm)

            # returns 1,2 or 0, respectively
            if p4.func.isDnaRnaOrProtein(s.temp):
                s.sequence = s.temp
                s.dataType = 'dna'
                s.symbols = 'acgt'
            else:
                s.sequence = s.temp
                s.dataType = 'protein'
                s.symbols = 'arndcqeghilkmfpstwyv'
            del(s.temp)

        # Having read in all the sequences, check for valid characters
        if len(self.sequences) > 0:
            bads = 0
            if self.sequences[0].dataType == 'dna':
                for s in self.sequences:
                    j = 0
                    while j < nChar:
                        if s.sequence[j] not in var.validDnaChars:
                            print("bad character '%s' in (zero-based) dna sequence %s " % \
                                (s.sequence[j], self.sequences.index(s)))
                            print("          sequence name: %s" % s.name)
                            print("          at (zero-based) position %s" % j)
                            bads = bads + 1
                            if bads > 10:
                                print("...and possibly others")
                                break
                        j = j + 1
                    if bads > 10:
                        break
                if bads:
                    raise P4Error(gm)
            elif self.sequences[0].dataType == 'protein':
                for s in self.sequences:
                    j = 0
                    while j < nChar:
                        if s.sequence[j] not in var.validProteinChars:
                            print("bad character '%s' in (zero-based) protein sequence %s " % \
                                (s.sequence[j], self.sequences.index(s)))
                            print("          sequence name: %s" % s.name)
                            print("          at (zero-based) position %s" % j)
                            bads = bads + 1
                            if bads > 10:
                                print("...and possibly others")
                                break
                        j = j + 1
                    if bads > 10:
                        break
                if bads:
                    raise P4Error(gm)

        if dbug:
            for t in self.sequences:
                print('%20s  %-30s' % ('name', t.name))
                print('%20s  %-30s' % ('comment', t.comment))
                print('%20s  %-30s' % ('sequence', t.sequence))
                print('%20s  %-30s' % ('type', t.dataType))
                print('')

    def writeNexus(self, fName=None, writeDataBlock=0,  interleave=0, flat=0, append=0, userText=''):
        """Write self in Nexus format.

        If writeDataBlock=1, then a data block is written, rather than the
        default, which is to write a taxa and a characters block.

        Flat gives sequences all on one line.
        Append, if 0, writes #NEXUS first.  If 1, does not write #NEXUS.

        userText is anything, eg a comment or another Nexus block, that
        you might want to add.  It goes at the end. """

        gm = ["Alignment.writeNexus()"]
        if not fName:
            fName = sys.stdout
        if fName == sys.stdout:
            f = sys.stdout
            if append:
                pass
            else:
                f.write('#NEXUS\n\n')
        else:
            if append:
                import os
                if os.path.isfile(fName):
                    try:
                        f = open(fName, 'a')
                    except IOError:
                        gm.append("Can't open %s for appending." % fName)
                        raise P4Error(gm)
                else:
                    print("Alignment: writeNexusFile() 'append' is requested,")
                    print("    but '%s' is not a regular file (maybe it doesn't exist?)." % fName)
                    print("    Writing to a new file instead.")
                    try:
                        f = open(fName, 'w')
                        f.write('#NEXUS\n\n')
                    except IOError:
                        gm.append("Can't open %s for writing." % fName)
                        raise P4Error(gm)

            else:
                try:
                    f = open(fName, 'w')
                    f.write('#NEXUS\n\n')
                except IOError:
                    gm.append("Can't open %s for writing." % fName)
                    raise P4Error(gm)

        if not writeDataBlock:
            f.write('begin taxa;\n')
            f.write('  dimensions ntax=%s;\n' % len(self.sequences))
            f.write('  taxlabels')
            for i in range(len(self.sequences)):
                f.write(
                    ' %s' % p4.func.nexusFixNameIfQuotesAreNeeded(self.sequences[i].name))
            f.write(';\n')
            f.write('end;\n\n')
        else:  # ie writeDataBlock=1
            f.write('begin data;\n')
            if self.excludeDelete:
                if self.length < self.excludeDelete.length:
                    f.write('  [%i characters have been excluded]\n' %
                            (self.excludeDelete.length - self.length))
            f.write('  dimensions ntax=%s' % len(self.sequences))
        if not writeDataBlock:
            f.write('begin characters;\n')
            if self.excludeDelete:
                if self.length < self.excludeDelete.length:
                    f.write('  [%i characters have been excluded]\n' %
                            (self.excludeDelete.length - self.length))
            f.write('  dimensions')

        f.write(' nChar=%s;\n' % self.length)

        f.write('  format')
        if self.dataType == 'dna':
            f.write(' datatype=dna')
        elif self.dataType == 'protein':
            f.write(' datatype=protein')
        elif self.dataType == 'rna':
            f.write(' datatype=rna')
        elif self.dataType == 'standard':
            f.write(' datatype=standard')
            f.write(' symbols=\"%s\"' % self.symbols)
        if self.equates:
            # We generally do not want to write the usual set of equates.
            if self.dataType == 'dna' or self.dataType == 'rna':
                usualset = set(list('nrykmswbdhv'))
            elif self.dataType == 'protein':
                usualset = set(['x', 'b', 'z'])
            else:
                usualset = set([])  # no usual equates for standard dataType

            # kset is the currently defined set. Now we want to get ks, the
            # ones to write.
            kset = set(self.equates.keys())
            if usualset == kset:
                ks = []   # write none
            # Does the usualset have items not in the kset?  This would be odd,
            # unexpected, but possible.
            elif usualset.difference(kset):
                ks = list(kset)  # write all, as it is so unusual.
            elif kset.difference(usualset):
                # inefficient calculating it twice ...
                ks = list(kset.difference(usualset))
            # print " [ks=%s]" % ''.join(ks),

            if ks:
                f.write(' equate=\"')
                for k in ks[:-1]:
                    if len(self.equates[k]) == 1:
                        f.write('%s=%s ' % (k, self.equates[k]))
                    else:
                        f.write('%s={%s} ' % (k, self.equates[k]))
                if len(self.equates[ks[-1]]) == 1:
                    f.write('%s=%s\"' % (ks[-1], self.equates[ks[-1]]))
                else:
                    f.write('%s={%s}\"' % (ks[-1], self.equates[ks[-1]]))
        if interleave:
            f.write(' interleave')
        f.write(' gap=-')
        f.write(' missing=?')
        f.write(';\n')
        f.write('  matrix\n')
        if interleave:
            if flat:
                gm.append(
                    "'interleave' option does not make sense with 'flat' option.")
                if f != sys.stdout:
                    f.close()
                raise P4Error(gm)
            else:
                # first, get the length of the longest name
                longest = 0
                for i in range(len(self.sequences)):
                    s = self.sequences[i]
                    if len(s.name) > longest:
                        longest = len(
                            p4.func.nexusFixNameIfQuotesAreNeeded(s.name))
                # formatString = '    %' + `-longest` + 's '  # boring
                # left-justified
                # cool right-justified
                formatString = '    %' + '%i' % longest + 's '
                # print "format string is '%s'" % formatString
                if longest > 10:
                    wid = 50
                else:
                    wid = 60
                pos = 0
                left = len(self.sequences[0].sequence)
                while left > 0:
                    for i in range(len(self.sequences)):
                        s = self.sequences[i]
                        f.write(formatString %
                                p4.func.nexusFixNameIfQuotesAreNeeded(s.name))
                        if left >= wid:
                            f.write('%s\n' % s.sequence[pos: pos + wid])
                        elif left > 0:
                            f.write('%s\n' % s.sequence[pos:])
                    pos = pos + wid
                    left = left - wid
                    if left > 0:
                        f.write('\n')

        if not interleave:
            if flat:
                # first, get the length of the longest name
                longest = 0
                for i in range(len(self.sequences)):
                    s = self.sequences[i]
                    if len(s.name) > longest:
                        longest = len(
                            p4.func.nexusFixNameIfQuotesAreNeeded(s.name))
                # formatString = '    %' + `-longest` + 's '  # boring
                # left-justified
                # cool right-justified
                formatString = '    %' + "%i" % longest + 's '
                # print "format string is '%s'" % formatString
                for i in range(len(self.sequences)):
                    s = self.sequences[i]
                    f.write(formatString % s.name)
                    f.write('%s\n' % s.sequence)
            if not flat:
                wid = 60
                for i in range(len(self.sequences)):
                    s = self.sequences[i]
                    f.write('    %s\n' %
                            p4.func.nexusFixNameIfQuotesAreNeeded(s.name))
                    left = len(s.sequence)
                    pos = 0
                    while left >= wid:
                        f.write('      %s\n' % s.sequence[pos: pos + wid])
                        pos = pos + wid
                        left = left - wid
                    if left > 0:
                        f.write('      %s\n' % s.sequence[pos:])
        f.write('  ;\n')
        f.write('end;\n\n')

        if self.nexusSets:
            # print "self.nexusSets = %s" % self.nexusSets
            if self.excludeDelete:
                if self.length < self.excludeDelete.length:
                    f.write('[skipping out-of-sync nexus sets block]\n')
                else:
                    self.nexusSets.writeNexusToOpenFile(f)
            else:
                self.nexusSets.writeNexusToOpenFile(f)

        if userText:
            f.write(userText)
            f.write("\n")

        if f != sys.stdout:
            f.close()

    def writePhylip(self, fName=None, interleave=True, whitespaceSeparatesNames=True, flat=False, offset=1, uppercase=False):
        """Write the alignment in Phylip format.

        Arg fName can be None, the default, in which case the
        alignment is written to stdout, or the fName can be an open
        file-like object, or fName can be a string, the name of a file
        to write to.  If it is a file-like object it will be re-wound
        to the beginning, but it is is up to the user to open and
        close it.

        If interleave is turned off (it is on by default), then
        sequences are written sequentially.

        Phylip and phylip-like formats are too varied.  The strict
        Phylip format has a set number of spaces for the taxon name,
        and there may not necessarily be a space between the name and
        the sequence.  The name size is commonly 10 spaces, but it
        need not be -- it is a compile-time option in Phylip.

        Other programs, eg phyml and PAML, use a phylip-like format
        where the tax name is set off from the sequence by whitespace.
        There is no set number of spaces that the sequence needs to
        occupy, and there may not be spaces in the tax name.  'offset'
        is the number of spaces from the end of the name to the
        beginning of the sequence.

        This method used to write strict, real phylip format by
        default, where there is a set number of spaces for the taxon
        name, and where there may not necessarily be a space between
        the name and the sequence.  The name size is commonly 10
        spaces, but it need not be-- it is set by
        var.phylipDataMaxNameLength (default 10).

        However, it no longer does that by default, as people use the
        format too loosely.  So now 'whitespaceSeparatesNames' is
        turned on by default.  It accommodates names longer than 10
        chars.

        If you want to write strict Phylip format, turn
        'whitespaceSeparatesNames' off.  Note that in that format,
        described `here
        <http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles>`_,
        it is ok to have blanks in the sequence, and so p4 puts a
        blank at the beginning of the sequence when writing this way,
        even though the meaning of 'whitespaceSeparatesNames' would
        suggest not.  Here the main effect of
        whitespaceSeparatesNames=False is that it is assumed that
        whitespace will not be used to trigger the end of the name in
        the downstream program.  Rather it is assumed that in the
        downstream program the end of the name will triggered by
        reading in var.phylipDataMaxNameLength characters.  If I were
        to set whitespaceSeparatesNames=True when I write (because I
        expect that the downstream program that will read it expects
        that) then I check for internal blank spaces in the names
        (which would be disallowed), and also I then allow more than
        var.phylipDataMaxNameLength characters.  However, in both
        cases there is at least one blank after the name.

        Arg 'flat' puts it all on one line.  This is not compatible
        with 'interleave', of course.
        """

        gm = ['Alignment.writePhylip(fName=%s, interleave=%s, whitespaceSeparatesNames=%s, flat=%s, uppercase=%s)' % (
            fName, interleave, whitespaceSeparatesNames, flat, uppercase)]

        # Find the maxNameLen.
        maxNameLen = 0
        namesHaveSpaces = False
        for s in self.sequences:
            theNameLen = len(s.name)
            if theNameLen > maxNameLen:
                maxNameLen = theNameLen
            if s.name.count(' '):
                namesHaveSpaces = True
        # print 'The longest name length in this alignment is %i' % maxNameLen

        if whitespaceSeparatesNames and namesHaveSpaces:
            crimes = []
            for s in self.sequences:
                if s.name.count(' '):
                    crimes.append("has space in the name: %s" % s.name)
            gm.append("whitespaceSeparatesNames is set, but some tax names")
            gm.append("have spaces.  That won't work. -- Fix it.")
            for crime in crimes:
                gm.append(crime)
            raise P4Error(gm)

        doStrictOkAnyway = False
        if maxNameLen < var.phylipDataMaxNameLength:
            doStrictOkAnyway = True

        if interleave and flat:
            gm.append("Both 'interleave' and 'flat' are turned on -- does not work.")
            raise P4Error(gm)

        # Check and complain if any taxNames will be truncated.
        if not whitespaceSeparatesNames and maxNameLen > var.phylipDataMaxNameLength:
            gm.append('The longest name length in this alignment is %i' % maxNameLen)
            gm.append('var.phylipDataMaxNameLength is now %i' %
                      var.phylipDataMaxNameLength)
            gm.append('Sequence names will be truncated.  Fix it.')
            gm.append("You may want to use the 'renameForPhylip()' method.")
            raise P4Error(gm)

        nameWid = var.phylipDataMaxNameLength + 1
        spacer1 = var.phylipDataMaxNameLength + 1
        if whitespaceSeparatesNames and (maxNameLen >= nameWid):
            nameWid = maxNameLen + offset
            spacer1 = 11
        # print 'The nameWid is %i' % nameWid

        fileNeedsClosing = True    # for usual files on disk, fName is a string
        isUserFlob = False
        if fName == None or fName == sys.stdout:
            f = sys.stdout
            fileNeedsClosing = False
        elif hasattr(fName, 'read'):          # a file-like object
            f = fName
            fileNeedsClosing = False
            isUserFlob = True
        else:
            try:
                f = open(fName, 'w')
            except IOError:
                gm.append("Can't open '%s' for writing." % fName)
                raise P4Error(gm)
        f.write(' %i  %i\n' % (len(self.sequences), self.length))

        if interleave:
            #wid1 = 50
            if nameWid < 50:
                wid1 = 61 - nameWid
            else:
                wid1 = 0

            # do the first row
            if self.length >= wid1:
                upper = wid1
            else:
                upper = self.length
            for k in self.sequences:
                theFormat = "%-" + "%is" % (nameWid)
                f.write(theFormat % k.name[0:nameWid])
                if uppercase:
                    f.write('%s\n' % k.sequence[0:upper].upper())
                else:
                    f.write('%s\n' % k.sequence[0:upper])
            lower = upper
            f.write('\n')
            # do subsequent rows
            wid2 = 50
            while lower < self.length:
                if self.length >= lower + wid2:
                    upper = lower + wid2
                else:
                    upper = self.length
                for k in self.sequences:
                    f.write(' ' * spacer1)
                    if uppercase:
                        f.write('%s\n' % k.sequence[lower:upper].upper())
                    else:
                        f.write('%s\n' % k.sequence[lower:upper])
                lower = upper
                f.write('\n')
        if not interleave:
            if flat:
                # print "nameWid = ", nameWid
                theFormat = "%-" + "%is" % nameWid
                for s in self.sequences:
                    f.write(theFormat % s.name[:nameWid])
                    if uppercase:
                        f.write('%s\n' % s.sequence.upper())
                    else:
                        f.write('%s\n' % s.sequence)
            else:
                #wid = 50
                if nameWid < 50:
                    wid1 = 61 - nameWid
                else:
                    wid1 = 0
                wid2 = 50
                theFormat = "%-" + "%is" % nameWid
                for s in self.sequences:
                    f.write(theFormat % s.name[:nameWid])
                    left = len(s.sequence)
                    pos = 0
                    if left >= wid1:
                        if uppercase:
                            f.write('%s\n' % s.sequence[pos: pos + wid1].upper())
                        else:
                            f.write('%s\n' % s.sequence[pos: pos + wid1])
                        pos = pos + wid1
                        left = left - wid1
                    else:
                        if uppercase:
                            f.write('%s\n' % s.sequence[pos:].upper())
                        else:
                            f.write('%s\n' % s.sequence[pos:])
                    while left >= wid2:
                        f.write(' ' * spacer1)
                        if uppercase:
                            f.write('%s\n' % s.sequence[pos: pos + wid2].upper())
                        else:
                            f.write('%s\n' % s.sequence[pos: pos + wid2])
                        pos = pos + wid2
                        left = left - wid2
                    if left > 0:
                        f.write(' ' * spacer1)
                        if uppercase:
                            f.write('%s\n' % s.sequence[pos:].upper())
                        else:
                            f.write('%s\n' % s.sequence[pos:])
                    f.write('\n')

        if isUserFlob:
            f.seek(0)

        if fileNeedsClosing:
            f.close()

    def writeMolphy(self, fName=None):
        """Another phylogenetics program, another format.

        Does anybody use this anymore?
        """

        gm = ["Alignment.writeMolphy()"]
        if fName == None or fName == sys.stdout:
            f = sys.stdout
        else:
            try:
                f = open(fName, 'w')
            except IOError:
                gm.append("Can't open %s for writing." % fName)
                raise P4Error(gm)
        f.write(' %i  %i\n' % (len(self.sequences), self.length))

        wid = 50
        for i in range(len(self.sequences)):
            s = self.sequences[i]
            f.write('%s\n' % s.name)
            left = len(s.sequence)
            pos = 0
            if left >= wid:
                f.write('%s\n' % s.sequence[pos: pos + wid])
                pos = pos + wid
                left = left - wid
            else:
                f.write('%s\n' % s.sequence[pos:])
            while left >= wid:
                f.write('%s\n' % s.sequence[pos: pos + wid])
                pos = pos + wid
                left = left - wid
            if left > 0:
                f.write('%s\n' % s.sequence[pos:])
            f.write('\n')

        if f != sys.stdout:
            f.close()

    def _readOpenGdeFile(self, flob, inverseMasks=0):
        dbug = 0
        gm = ['Alignment._readOpenGdeFile()']
        sList = []
        maskList = []
        aLine = flob.readline()
        aLine = aLine.strip()
        if dbug:
            print(aLine)
        if aLine[0] != '{':
            gm.append("The first character is not '{' --- not a GDE file??")
            gm.append(
                "(Note that p4 only reads proper gde files, not gde flat files.)")
            raise P4Error(gm)
        while aLine:
            s = Sequence()
            if aLine[0] == '{':
                aLine = flob.readline().strip()
            else:
                gm.append("Problem with GDE file at line:\n\t%s" % aLine)
                raise P4Error(gm)
            if aLine[:4] == 'name':
                # the name line might be like: name    "D.ethanoge"
                # or it might be like: name    "C.acetobA ", ie with a spurious
                # space.
                if dbug:
                    print("got name line: \n%s" % aLine)
                splitLine = aLine.split('\"')
                s.name = splitLine[1].strip()
                # print s.name
                if not p4.func.nexusCheckName(s.name):
                    gm.append("Bad name '%s'" % s.name)
                    raise P4Error(gm)
            else:
                gm.append(
                    "I was expecting a name line: instead I got line:\n\t%s" % aLine)
                raise P4Error(gm)
            aLine = flob.readline().strip()
            if aLine[:4] == 'type':
                if dbug:
                    print("get type line: \n%s" % aLine)
                splitLine = aLine.split()
                type = splitLine[1][1:-1]
                # print type
                if type not in ['DNA', 'PROTEIN', 'MASK', 'TEXT', 'RNA']:
                    gm.append("Unknown type '%s'" % type)
                    raise P4Error(gm)
            else:
                gm.append(
                    "I was expecting a type line: instead I got line:\n\t%s" % aLine)
                raise P4Error(gm)

            if type == 'DNA' or type == 'RNA':
                s.dataType = 'dna'
            elif type == 'PROTEIN':
                s.dataType = 'protein'
            elif type == 'MASK':
                # a dodgy manoever, letting a Sequence object hold a mask or
                # text.
                s.dataType = 'mask'
            elif type == 'TEXT':
                s.dataType = 'text'

            if dbug:
                print("dataType set to: %s" % s.dataType)

            # Next, we are on the lookout for an 'offset' line, which
            # may or may not be there.  If we come across the first
            # line of the sequence, break
            offset = 0
            while 1:
                aLine = flob.readline()
                # print "a Reading line: %s" % aLine
                if aLine[:6] == 'offset':
                    aLine = aLine.strip()
                    if dbug:
                        print("got offset line: \n%s" % aLine)
                    splitLine = aLine.split()
                    try:
                        offset = int(splitLine[1])
                        if dbug:
                            print("get offset: %i" % offset)
                    except ValueError:
                        gm.append("Bad offset in line '%s'" % aLine)
                        raise P4Error(gm)
                elif aLine[:8] == 'sequence' and aLine[:11] != 'sequence-ID':
                    aLine = aLine.strip()
                    break

            # Now aLine has the first line of the sequence
            if dbug:
                print("\nGot first line of sequence: %s" % aLine)
            # Get the rest of the sequence.
            splitLine = aLine.split()
            thisSeqList = []
            thisSeqList.append(splitLine[1])
            aLine = flob.readline().strip()
            while aLine[0] != '}':
                thisSeqList.append(aLine)
                aLine = flob.readline().strip()
            # print "Got thisSeqList: %s" % thisSeqList
            s.sequence = (offset * '-') + ''.join(thisSeqList)[1:-1]
            # print "Finished: name:'%s', dataType %s, sequence = '%s'" % (s.name, s.dataType, s.sequence)
            # print "Got last line of group: '%s'" % aLine\
            if s.dataType == 'mask':
                maskList.append(s)
            else:
                sList.append(s)

            while 1:
                aLine = flob.readline()
                if not aLine:
                    break
                if aLine[0] == '{':
                    break

        # At this point we have read in to the end of the file, and
        # have a sList containing Sequences, and a maskList containing
        # masks.

        # Make all the sequences the same length.  GDE certainly does not do
        # this.
        if len(sList) == 0:
            gm.append(
                "Finished reading file, but got no sequences.  What gives?")
            raise P4Error(gm)
        maxSeqLen = 0
        for s in sList:
            if len(s.sequence) > maxSeqLen:
                maxSeqLen = len(s.sequence)
        for s in maskList:
            if len(s.sequence) > maxSeqLen:
                maxSeqLen = len(s.sequence)
        for s in sList:
            if len(s.sequence) != maxSeqLen:
                s.sequence = s.sequence + ((maxSeqLen - len(s.sequence)) * '-')
        for s in maskList:
            if len(s.sequence) != maxSeqLen:
                s.sequence = s.sequence + ((maxSeqLen - len(s.sequence)) * '-')

        # Now we have a bunch of sequences.  Fill in the alignment (ie self)
        self.sequences = sList
        self.checkNamesForDupes()
        if 0:
            self.checkLengthsAndTypes()
            self.nexusSets = NexusSets()
            self.nexusSets.nChar = self.length
            self.nexusSets.setPredefinedCharSets(self)
        else:
            # check lengths
            topLen = len(self.sequences[0].sequence)
            for s in self.sequences:
                # print "xxx dataType=%s" % s.dataType
                if len(s.sequence) != topLen:
                    gm.append("Got sequences of unequal lengths")
                    raise P4Error(gm)
            self.length = topLen

            # check dataType's
            topType = self.sequences[0].dataType
            for s in self.sequences:
                if s.dataType != topType:
                    gm.append("Mixed dataTypes.")
                    raise P4Error(gm)
            self.dataType = topType
            if self.dataType == 'dna':
                self.symbols = 'acgt'
                self.dim = 4
                if not self.equates:
                    self.equates = {'n': 'acgt', 'm': 'ac', 'k': 'gt',  # 'x': 'acgt',
                                    'h': 'act', 'y': 'ct', 'v': 'acg',
                                    'w': 'at', 'd': 'agt', 'b': 'cgt',
                                    'r': 'ag', 's': 'cg'}
            elif self.dataType == 'protein':
                self.symbols = 'arndcqeghilkmfpstwyv'
                self.dim = 20
                if not self.equates:
                    self.equates = {
                        'b': 'dn', 'x': 'arndcqeghilkmfpstwyv', 'z': 'eq'}
            else:
                gm.append("Can't deal with dataType '%s'." % self.dataType)
                raise P4Error(gm)

        # self.dump()
        # Is the following really needed?
        # for s in self.sequences:
        #    if s.name not in self.taxNames:
        #        self.taxNames.append(s.name)
        #    else:
        #        gm.append("Duplicated name '%s'" % s.name)
        #        raise P4Error(gm)

        if 1:
            if len(maskList):
                self.nexusSets = NexusSets()
            for s in maskList:
                c = CharSet()
                c.nexusSets = self.nexusSets
                c.nChar = self.length
                c.name = s.name
                c.lowName = s.name.lower()
                c.format = 'vector'
                if inverseMasks:
                    s.sequence = list(s.sequence)
                    for i in range(len(s.sequence)):
                        if s.sequence[i] == '0':
                            s.sequence[i] = '1'
                        else:
                            s.sequence[i] = '0'
                    s.sequence = ''.join(s.sequence)
                c.mask = s.sequence
                for c1 in self.nexusSets.charSets:
                    if c1.lowName == c.lowName:
                        gm.append(
                            "Duplicated charSet (ie mask) name %s" % c.name)
                        raise P4Error(gm)
                self.nexusSets.charSets.append(c)



    ####################################################### part
    #######################################################

    def _initParts(self):
        gm = ['Alignment._initParts()']

        if len(self.parts):
            for p in self.parts:
                del(p)
        self.parts = []
        if self.equates:
            eqSymb = list(self.equates.keys())
            eqSymb.sort()
            eqSymb = ''.join(eqSymb)
        else:
            eqSymb = ''

        if len(self.sequences) and self.length and self.symbols and self.dim:
            pass
        else:
            gm.append("Can't allocate part.")
            if not len(self.sequences):
                gm.append("-no sequences.")
            elif not self.length:
                gm.append("-the sequences have no length")
            elif not self.symbols:
                gm.append("-no symbols")
            elif not self.dim:
                gm.append("-dim not set")
            raise P4Error(gm)

        # its all one part
        if not self.nexusSets or not self.nexusSets.charPartition:
            aPart = Part()
            aPart.alignment = self
            aPart.name = 'all'
            aPart.lowName = 'all'
            aPart.dataType = self.dataType
            aPart.dim = self.dim
            aPart.symbols = self.symbols
            aPart.equates = self.equates
            aPart.nTax = len(self.sequences)
            aPart.nChar = self.length
            assert aPart.nChar

            if 0:
                print(gm[0])
                print("    symbols=%s" % self.symbols)

            aPart.cPart = pf.newPart(len(self.sequences), self.length,
                                     eqSymb, len(eqSymb), self.symbols, len(self.symbols))
            if not aPart or not aPart.cPart:
                gm.append("Failed to get memory for part.")
                raise P4Error(gm)

            # Make the equates table
            verbose = 0
            equatesTable = []
            if verbose:
                print("equates is %s" % self.equates)
                print("eqSymb is %s" % eqSymb)  # the keys
                print("symbols is %s" % self.symbols)
            for i in range(len(eqSymb)):
                if verbose:
                    print("%3s: " % eqSymb[i], end=' ')
                e = self.equates[eqSymb[i]]
                if verbose:
                    print("%8s : " % e, end=' ')
                for s in self.symbols:
                    if s in e:
                        if verbose:
                            print("%1i" % 1, end=' ')
                        equatesTable.append('1')
                    else:
                        if verbose:
                            print("%1i" % 0, end=' ')
                        equatesTable.append('0')
                if verbose:
                    print('')
            equatesTable = ''.join(equatesTable)
            if verbose:
                print("\n\nequatesTable:")
                print(equatesTable)
            pf.pokeEquatesTable(aPart.cPart, equatesTable)

            sList = []
            for s in self.sequences:
                sList.append(s.sequence)
            if 0:
                print(gm[0])
                print("sList = %s" % sList)
                print("joined = %s" % ''.join(sList))
            pf.pokeSequences(aPart.cPart, ''.join(sList))
            # print "about to makePatterns ..."
            pf.makePatterns(aPart.cPart)
            # print "about to setInvar"
            pf.setGlobalInvarSitesVec(aPart.cPart)

            # pf.dumpPart(aPart.cPart)
            self.parts.append(aPart)

        elif self.nexusSets.charPartition:
            for cpp in self.nexusSets.charPartition.subsets:
                # print "Doing subset '%s', mask: %s" % (cpp.name, cpp.mask)
                # print "About to subsetUsingMask (self length is %i)" %
                # self.length
                b = self.subsetUsingMask(cpp.mask)
                # This very method, but now there are no charPartitions in b.
                b._initParts()
                b.parts[0].name = cpp.name
                b.parts[0].lowName = cpp.name.lower()
                self.parts.append(b.parts[0])
                b.parts = []  # so we don't try free-ing it twice

    def initDataParts(self):
        gm = ['Alignment.initDataParts()']

        if len(self.parts):
            for p in self.parts:
                del(p)
        self.parts = []

        if len(self.sequences) and self.length and self.symbols and self.dim:
            pass
        else:
            gm.append("Can't allocate part.")
            if not len(self.sequences):
                gm.append("-no sequences.")
            elif not self.length:
                gm.append("-the sequences have no length")
            elif not self.symbols:
                gm.append("-no symbols")
            elif not self.dim:
                gm.append("-dim not set")
            raise P4Error(gm)

        # its all one part
        if not self.nexusSets or not self.nexusSets.charPartition:
            aPart = DataPart(self)
            self.parts.append(aPart)

        elif self.nexusSets.charPartition:
            for cpp in self.nexusSets.charPartition.subsets:
                # print "Doing subset '%s', mask: %s" % (cpp.name, cpp.mask)
                # print "About to subsetUsingMask (self length is %i)" %
                # self.length
                b = self.subsetUsingMask(cpp.mask)
                # This very method, but now there are no charPartitions in b.
                b.initDataParts()
                b.parts[0].alignment = self
                b.parts[0].name = cpp.name
                b.parts[0].lowName = cpp.name.lower()
                self.parts.append(b.parts[0])
                b.parts = []  # so we don't try free-ing the new part twice

    def resetSequencesFromParts(self):
        """Gets the sequences from Part.cPart, and installs them in self."""

        # print "Alignment.resetSequencesFromParts() here."
        if (not self.parts) or len(self.parts) == 0:
            gm = ["Alignment.resetSequencesFromParts()"]
            gm.append("No parts.")
            raise P4Error(gm)

        if not var.doDataPart:
            if len(self.parts) == 1 and self.parts[0].name == 'all':
                allSeq = pf.symbolSequences(self.parts[0].cPart)
                # print "allSeq[0:20] = %s" % allSeq[0:20]
                for i in range(len(self.sequences)):
                    self.sequences[i].sequence = allSeq[
                        (i * self.length): ((i + 1) * self.length)]
            else:
                for i in range(len(self.sequences)):
                    self.sequences[i].sequence = list(
                        self.sequences[i].sequence)
                for i in range(len(self.parts)):
                    partSeq = pf.symbolSequences(self.parts[i].cPart)
                    # print partSeq
                    spot = 0
                    m = self.nexusSets.charPartition.subsets[i].mask
                    for s in self.sequences:
                        for k in range(self.length):
                            if m[k] == '1':
                                s.sequence[k] = partSeq[spot]
                                spot += 1
                for i in range(len(self.sequences)):
                    self.sequences[i].sequence = ''.join(self.sequences[i].sequence)
        else:
            if len(self.parts) == 1:
                for i in range(len(self.sequences)):
                    self.sequences[i].sequence = self.parts[
                        0].sequenceString(i)
            else:
                for i in range(len(self.sequences)):
                    self.sequences[i].sequence = list(
                        self.sequences[i].sequence)
                for pNum in range(len(self.parts)):
                    for sNum in range(len(self.sequences)):
                        partSeq = self.parts[pNum].sequenceString(sNum)
                        print(partSeq)
                        spot = 0
                        m = self.nexusSets.charPartition.subsets[pNum].mask
                        s = self.sequences[sNum]
                        for k in range(self.length):
                            if m[k] == '1':
                                s.sequence[k] = partSeq[spot]
                                spot += 1
                for i in range(len(self.sequences)):
                    self.sequences[i].sequence = ''.join(self.sequences[i].sequence)

    def resetPartsContentFromSequences(self):
        """Reset Part.cPart sequences from self.sequences.

        It then makes patterns, and sets the global invariant sites
        array.  """

        gm = ['Alignment.resetPartsContentFromSequences()']
        if len(self.parts) == 1:  # its all one part
            aPart = self.parts[0]
            if not var.doDataPart:
                sList = []
                for s in self.sequences:
                    sList.append(s.sequence)
                pf.pokeSequences(aPart.cPart, ''.join(sList))
                # are the following necessary?
                pf.makePatterns(aPart.cPart)
                pf.setGlobalInvarSitesVec(aPart.cPart)
            else:
                for sNum in range(len(self.sequences)):
                    s = self.sequences[sNum]
                    for cNum in range(self.length):
                        theChar = s.sequence[cNum]
                        if theChar == '-':
                            aPart.seq[sNum, cNum] = var.GAP_CODE
                        if theChar == '?':
                            aPart.seq[sNum, cNum] = var.QMARK_CODE
                        elif aPart.equateSymbols and theChar in aPart.equateSymbols:
                            aPart.seq[sNum, cNum] = var.EQUATES_BASE + \
                                aPart.equateSymbols.index(theChar)
                        else:
                            aPart.seq[sNum, cNum] = aPart.symbols.index(
                                theChar)

        elif len(self.parts) > 1:
            # the number of parts is also the length of the subsets list
            if self.nexusSets and self.nexusSets.charPartition and \
                    self.nexusSets.charPartition.subsets and \
                    len(self.nexusSets.charPartition.subsets) == len(self.parts):
                pass
            else:
                gm.append(
                    'Something is wrong with the nexusSets or its charPartition')
                raise P4Error(gm)
            for i in range(len(self.parts)):
                cpSubset = self.nexusSets.charPartition.subsets[i]
                aPart = self.parts[i]
                b = self.subsetUsingMask(cpSubset.mask)
                if not var.doDataPart:
                    sList = []
                    for s in b.sequences:
                        sList.append(s.sequence)
                    pf.pokeSequences(aPart.cPart, ''.join(sList))
                    # are the following necessary?
                    pf.makePatterns(aPart.cPart)
                    pf.setGlobalInvarSitesVec(aPart.cPart)
                else:
                    for sNum in range(len(b.sequences)):
                        s = b.sequences[sNum]
                        for cNum in range(b.length):
                            theChar = s.sequence[cNum]
                            if theChar == '-':
                                aPart.seq[sNum, cNum] = var.GAP_CODE
                            if theChar == '?':
                                aPart.seq[sNum, cNum] = var.QMARK_CODE
                            elif aPart.equateSymbols and theChar in aPart.equateSymbols:
                                aPart.seq[
                                    sNum, cNum] = var.EQUATES_BASE + aPart.equateSymbols.index(theChar)
                            else:
                                aPart.seq[sNum, cNum] = aPart.symbols.index(
                                    theChar)

        else:
            gm.append("No parts.")
            raise P4Error(gm)

