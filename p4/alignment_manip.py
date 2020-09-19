from p4.sequence import Sequence
from p4.sequencelist import SequenceList
from p4.nexussets import NexusSets
from p4.p4exceptions import P4Error
import string
import copy
import os
import math
import string
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
from p4.alignment import ExcludeDelete
from p4.geneticcode import GeneticCode

if True:
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


    def getCharDiversityDistribution(self):
        """Distribution of different chars per site, as a distribution.

        Gaps and ambiguities are ignored.

        A tuple is returned, composed of

        +----------+---------------------------------------------------------+
        |    index | \                                                       |
        +==========+=========================================================+
        |        0 | num of autapomorphies                                   |
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

        Autapomorphies are  sites with two kinds of characters, where there is
        only one of one of the characters.

        So if there are 10 sites with 2 kinds of char, and 4 of those are
        autapomorphies, 30 constant sites, and 20 sites with 3 kinds of character, we
        would have a distro like this --- (4, 30, 10, 20, ...)

        This is in pure Python, and can be used for Alignments with one
        part.

        Returns:
            a tuple of ints, showing counts of sites with different diversities.

        See also :meth:`~Alignment.getMaskForCharDiversity`

        """
        #assert isinstance(self, Alignment)
        distro = [0] * (len(self.symbols) + 1)
        counters = {}  # so that I can index it with a symbol rather than an int
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
            thereIsASingleton = False  
            for k,v in counters.items():
                if v:
                    hitsAtThisPos += 1
                    if v == 1:
                        thereIsASingleton = True
                counters[k] = 0                # Initialize for next pos
            assert hitsAtThisPos, "pos %i had no symbol chars, only gaps and ambigs."

            distro[hitsAtThisPos] += 1
            if hitsAtThisPos == 2:
                if thereIsASingleton:
                    distro[0] += 1
        # Check
        assert self.nChar == sum(distro[1:]), "Something is wrong. The distribution does not add up to nChar"
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

        This is done only in Python (ie no cParts involved), and does not
        handle partitioned data.  If you want to bootstrap partitioned
        data, use the :meth:`p4.data.Data.bootstrap` method.  """

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

    def recodeDayhoff(self, firstLetter=False):
        """Recode protein data into Dayhoff groups, in place.

        1.  c
        2.  stpag
        3.  ndeq
        4.  hrk
        5.  milv
        6.  fyw

        The ambiguity character 'x' is recoded as a gap.

        It does not make a new alignment-- it does the re-coding 'in-place'.

        If arg *firstLetter* is set, then the character is recoded as the
        first letter of its group rather than as a number.  Eg k would be
        recoded as h rather than as 4.
        """

        gm = ['Alignment.recodeDayhoff()']
        if self.dataType != 'protein':
            gm.append("This is only for protein alignments.")
            raise P4Error(gm)

        for s in self.sequences:
            s.dataType = 'standard'
            s.sequence = list(s.sequence)
            for i in range(len(s.sequence)):
                c = s.sequence[i]
                if c == 'c':
                    if firstLetter:
                        pass
                    else:
                        s.sequence[i] = '1'
                elif c in 'stpag':
                    if firstLetter:
                        s.sequence[i] = 's'
                    else:
                        s.sequence[i] = '2'
                elif c in 'ndeq':
                    if firstLetter:
                        s.sequence[i] = 'n'
                    else:
                        s.sequence[i] = '3'
                elif c in 'hrk':
                    if firstLetter:
                        s.sequence[i] = 'h'
                    else:
                        s.sequence[i] = '4'
                elif c in 'milv':
                    if firstLetter:
                        s.sequence[i] = 'm'
                    else:
                        s.sequence[i] = '5'
                elif c in 'fyw':
                    if firstLetter:
                        s.sequence[i] = 'f'
                    else:
                        s.sequence[i] = '6'
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
        if firstLetter:
            self.symbols = 'csnhmf'
        else:
            self.symbols = '123456'

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

        Constant sites can optionally be ignored.  Gaps and ambiguities are ignored.

        This is in pure Python, and can be used for Alignments with one
        part.  It is also implemented in c in the Data class, which allows
        more than one part (but no distribution).

        """
        #assert isinstance(self, Alignment)
        if showDistribution:
            distro = {}
        counters = {}
        for c in self.symbols:
            counters[c] = 0
        hits = 0
        nPos = 0  # the number of sites compared
        for pos in range(len(self)):
            for seq in self.sequences:
                c = seq.sequence[pos]
                if c in self.symbols:   # ie ignore gaps and ambiguities
                    counters[c] += 1
            #hits = 0
            hitsAtThisPos = 0
            for k, v in counters.items():
                if v:
                    hitsAtThisPos += 1
                    counters[k] = 0
            if includeConstantSites:
                hits += hitsAtThisPos
                nPos += 1
            else:
                if hitsAtThisPos > 1:   # ignore constant sites
                    hits += hitsAtThisPos
                    nPos += 1
            if showDistribution:
                if hitsAtThisPos in distro:
                    distro[hitsAtThisPos] += 1
                else:
                    distro[hitsAtThisPos] = 1
            # print pos, hits
        #print(hits, nPos)
        if showDistribution:
            kk = list(distro.keys())
            kk.sort()
            theMin = kk[0]
            theMax = kk[-1]
            print("\nnChars count")
            print("------ -----")
            # for k in kk:
            #    print "%2i  %3i" % (k, distro[k])
            # print
            for k in range(theMin, theMax + 1):
                # Tom Williams 7 Feb 2011 bug report and fix (Thanks!). Was 'if
                # distro[k]:' -- no workee if k is not a key.
                if k in distro:
                    print("%4i   %4i" % (k, distro[k]))
                else:
                    print("%4i   %4i" % (k, 0))
        return float(hits) / nPos

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

        scc = numpy.zeros(self.dim, numpy.float)
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

        bigF = numpy.zeros([self.dim, self.dim], numpy.float)
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

        d = numpy.zeros((dim - 1), numpy.float)
        V = numpy.zeros((dim - 1, dim - 1), numpy.float)

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

    def matchedPairsTests(self, mostSignificantOnly=False):
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

        The tests are pair-wize on all pairs of sequences.  Note that
        it may fail to do the calculations for a pair.  If so it will
        return None for that test for that pair, and that will end up
        in ``DistanceMatrix`` objects that are returned.

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
        and is part of ``iqtree2``.  (There was a small bug in ``--symtest``
        that was fixed in v 2.0.6, June 2020.)

        See Naser-Khdour et al GBE 2019 https://doi.org/10.1093/gbe/evz193

        This is my attempt to replicate it.  That implementation
        chooses the sequence pair with the biggest divergence, and
        only reports stats for that pair.  This does not necessarily
        show the smallest stats.

        If it is verbose, it speaks the results to stdout.

        It returns the 3 stats.

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
            print("... between these (zero-based)sequence numbers:")
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
        
        QB, QS, QR, PB, PS, PR = _ababnehEtAlStatsAndProbs(myBigF, self.dim, myTxNumPair[0], myTxNumPair[1])
        if verbose:
            print(f"PB:{PB} PS:{PS} PR:{PR}")
            print("(Turn off this verbose output by setting arg vebose=False)")
        return PB, PS, PR
        
        

    def testOverallFromAbabnehEtAl2006(self):
        """Marginal symmetry (Stuarts's) test for more than two matched sequences

        As in Ababneh et al 2006 Page 1227, calculated as in the R function Testoverall.

        A tuple composed of the Ts stat, the degrees of freedom, and the probability, is returned.
        """

        # Make the J matrix.
        J = numpy.zeros(
            (self.dim * (self.nTax - 1), self.dim * (self.nTax - 1)), numpy.float)
        for seqNum in range(self.nTax - 1):
            start = seqNum * self.dim
            for i in range(self.dim):
                for j in range(self.dim):
                    J[start + i, start + j] = 1.
        # print J

        # Make the L matrix.
        L = numpy.zeros(
            (self.dim * (self.nTax - 1), self.dim * self.nTax), numpy.float)
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
            (self.dim * self.nTax, self.dim * self.nTax), numpy.float)
        n = numpy.zeros(self.dim * self.nTax, numpy.float)

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
        myIdentity = numpy.identity(self.dim * (self.nTax - 1), numpy.float)
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

