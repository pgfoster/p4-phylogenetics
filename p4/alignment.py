from p4.sequence import Sequence
from p4.sequencelist import SequenceList
from p4.nexussets import NexusSets
from p4.p4exceptions import P4Error
import string
import copy
import os
import math
import string
import re
import sys
from p4.nexussets import CharSet
import subprocess
from p4.distancematrix import DistanceMatrix
from p4.var import var
from p4.part import Part
import numpy
import numpy.linalg

cListPat = re.compile('(\d+)-?(.+)?')
cList2Pat = re.compile('(.+)\\\\(\d+)')
cListAllPat = re.compile('all\\\\?(\d+)?')

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

    from p4.alignment_manip import simpleConstantMask, constantMask,gappedMask,getLikelihoodTopologyInformativeSitesMask,getMaskForAutapomorphies,getMaskForCharDiversity,getCharDiversityDistribution,orMasks,andMasks,sequenceSlice,bluntEndLigate, concatenate, constantSitesProportion, constantSitesCount, noGapsOrAmbiguitiesCopy, hasGapsOrAmbiguities, bootstrap, compositionEuclideanDistanceMatrix, covarionStats, pDistances, recodeDayhoff, recodeProteinIntoGroups, recodeRY, checkTranslation, translate, excludeCharSet, dupe, putGaps, setGBlocksCharSet, meanNCharsPerSite, getAllGapsMask, getEnoughCharsMask, simpleCharCounts, getSimpleBigF, matchedPairsTests, symtestAsInIQTreeNaserKhdour, testOverallFromAbabnehEtAl2006, getMinmaxChiSqGroups, getKosiolAISGroups, mrpSlice
    from p4.alignment_readwrite import readOpenPhylipFile, _readPhylipSequential, _readPhylipInterleaved, _readPhylipSequentialStrict, _readPhylipInterleavedStrict, _readOpenClustalwFile,writeNexus,writePhylip,writeMolphy,_readOpenGdeFile
    from p4.alignment_part import _initParts,initDataParts,resetSequencesFromParts,resetPartsContentFromSequences
    # from p4.alignment_logdet import logDet


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

