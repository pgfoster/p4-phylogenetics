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
from p4.nexussets import CharSet
import subprocess
from p4.distancematrix import DistanceMatrix
from p4.var import var
from p4.part import Part
import numpy
import numpy.linalg
import p4.pf as pf

if True:
    def logDet(self, pInvar=None, missingCharacterStrategy='refuse', verbose=True):
        """logDet calculations.  Returns a DistanceMatrix object.

        Heavily influenced by paup and by LDDist.  Thanks to Swofford
        and Thollesson, the authors.

        The basis of the log det is the Fxy matrix, the matrix of
        observed changes between 2 sequences.  It is from the Fxy that
        the determinant is calculated.  The negative of the log of the
        determinant ('logdet') is related to the evolutionary distance
        separating the 2 sequences.

        **Corrections**

        The logdet can be corrected so that it expresses the distance
        between 2 sequences in the usual terms of mutations per site.
        There are 2 different corrections-- that from Lockhart94 Eqn3
        (pg 606), and Tamura and Kumar 2002.  Which one is given by
        the arg 'correction', which should be one of 'L94' or 'TK02'.
        The former is used in paup, and the latter is used in LDDist.
        The current implementation of this method only uses L94.

        **Ambiguities**

        In sequence pairs, when there are 2 N-like chars (n, gap, ?),
        they are completely ignored.  Other combinations are not
        ignored.  Consider this alignment::

            acgtr
            acgta

        Now you might think that r-a would mean either a-a or g-a, but
        in the present implementation it is slightly different-- the
        ambiguity is added to Fxy in proportion to the observed
        unambiguous counts.  So in this case there was an observed
        a-a, but no g-a, so r-a is considered to be wholly a-a.  Now
        consider the alignment::

            aagcgtr
            aaacgta

        Here r-a is considered to be both a-a and g-a, in proportion
        to their occurrence, ie 2/3 a-a and 1/3 g-a.

        **pInvar**

        The pInvar should not exceed the proportion of invariant
        sites.

        **Dealing with missing (or low-frequency) characters**

        If any character in any pairwise comparison in either sequence
        is missing then the calculation fails because of the
        correction.  (The correction involves the log of the product
        of all the compositions; if any comp is zero then the product
        will be zero, and the log will be undefined --- Boom!).  There
        are 2 strategies implemented for dealing with missing chars--

            1. refuse This is the way paup does it.  Refuse to
               calculate.  Make up a 'big distance' later.  By default
               in paup, the 'big distance' is twice the biggest
               defined distance (ie defined elsewhere in the distance
               matrix, from some other pairwise comparison).  Of
               course this will fail if all the pairs are refused,
               leaving no defined distances.

            2. fudge This is the way that LDDist does it.  LDDist
               replaces zeros with 0.5 in the un-normalized Fxy
               martrix.  The present implementation replaces zeros on
               the diagonal of Fxy with either 0.5 or half of the
               smallest positive Fii, whichever is smaller.

        **Non-positive determinants**

        Sometimes for very diverged sequences the determinant of the
        Fxy matrix is zero or less, which means that the logarithm is
        undefined.  (This can happen even if all Fii are positive).
        The current implementation does not deal with them except by
        putting a -1.0 in the resulting DistanceMatrix.  The user is
        expected to deal with them after, as they see fit.  Paup
        replaces the values with twice the biggest valid logDet, and
        something like that seems at least workable.
        """

        complaintHead = 'Alignment.logDet()'
        gm = [complaintHead]

        correction = 'L94'
        fastFillFxy = True  # in c

        # Check that there is only one partition
        # its all one part
        if not self.nexusSets or not self.nexusSets.charPartition:
            pass
        # its partitioned.  Bad.
        elif self.nexusSets.charPartition and len(self.nexusSets.charPartition.subsets) > 1:
            gm = [complaintHead]
            gm.append("This only works with Alignments having only one data partition.")
            raise P4Error(gm)

        # Check the corrections arg
        goodCorrections = ['L94']
        if correction not in goodCorrections:
            gm.append("The corrections arg should be one of: %s" %goodCorrections)
            gm.append("Got %s" % correction)
            raise P4Error(gm)

        # Check the pInvar or pInvarOfConstants args.  If zero, set to None.
        try:
            pInvar = float(pInvar)
            if math.fabs(pInvar) < 1.0e-10:
                pInvar = None
        except ValueError:
            pInvar = None
        except TypeError:
            pInvar = None
        if pInvar and (pInvar < 0.0 or pInvar > 1.0):
            gm.append("pInvar, if set, should be between zero and 1.0, inclusive.")
            raise P4Error(gm)

        alignmentComposition = self.composition()

        # Check the missingCharacterStrategy arg
        goodMissingCharacterStrategies = ['refuse', 'fudge']  # 'reduce' is turned off
        if missingCharacterStrategy in goodMissingCharacterStrategies:
            pass
        else:
            gm.append("Arg missingCharacterStrategy should be one of %s" %
                      goodMissingCharacterStrategies)
            gm.append("Got %s" % missingCharacterStrategy)
            raise P4Error(gm)


        # Equates
        if self.equates:
            self.equateSymbols = list(self.equates)
            self.equateSymbols.sort()

            # Make equatesArray
            equatesArray = numpy.zeros((self.nEquates, self.dim), numpy.int32)
            for i in range(self.nEquates):
                e = self.equates[self.equateSymbols[i]]
                for j in range(self.dim):
                    s = self.symbols[j]
                    if s in e:
                        equatesArray[i, j] = 1
            # print "equatesArray = "
            # print equatesArray
        else:
            equatesArray = numpy.array([-1], numpy.int32)

        # We need a charsLikeN list, showing characters that are fully
        # ambiguous.
        charsLikeN = []
        if self.equates:
            for i in range(self.nEquates):
                e = self.equateSymbols[i]
                isNLike = True
                for j in range(self.dim):
                    if not equatesArray[i, j]:
                        isNLike = False
                        break
                if isNLike:
                    charsLikeN.append(e)
        charsLikeN.append('-')
        charsLikeN.append('?')
        # print "charsLikeN = %s" % charsLikeN

        # Make a seq, a numpy array, in which to hold the recoded sequences.
        seq = numpy.zeros((self.nTax, self.nChar), numpy.int8)
        for sNum in range(self.nTax):
            s = self.sequences[sNum]
            for cNum in range(self.nChar):
                theChar = s.sequence[cNum]
                if theChar in charsLikeN:
                    seq[sNum, cNum] = var.N_LIKE
                elif self.equates and theChar in self.equates:
                    seq[sNum, cNum] = var.EQUATES_BASE + \
                        self.equateSymbols.index(theChar)
                else:
                    seq[sNum, cNum] = self.symbols.index(theChar)
        # print seq

        # The 'refUnambigCountMatrix' array is raw counts of changes between
        # the two sequences.
        refUnambigCountMatrix = numpy.zeros((self.dim, self.dim), numpy.int32)
        normUnambig = numpy.zeros((self.dim, self.dim), numpy.float)
        allSymbolNums = list(range(self.dim))
        allSymbolNums += list(range(var.EQUATES_BASE,
                               var.EQUATES_BASE + self.nEquates))
        allSymbolNums.append(var.N_LIKE)
        allSymbolNums = numpy.array(allSymbolNums, numpy.int32)
        # print "allSymbolNums = ", allSymbolNums
        bigDim = self.dim + 1 + len(self.equates)
        # print "bigDim = %i" % bigDim
        # Int, to keep track of ambiguous changes
        refAmbigCountMatrix = numpy.zeros((bigDim, bigDim), numpy.int32)
        # normUnambig = zeros((self.dim, self.dim), Float) # 'double' type
        bigFxy = numpy.zeros((self.dim, self.dim), numpy.float)
        fudgeCount = 0
        invertCount = 0

        nUnambig = numpy.zeros((1), numpy.int32)
        nAmbig = numpy.zeros((1), numpy.int32)
        nDoubleGap = numpy.zeros((1), numpy.int32)

        # #######################
        # The main pairwise loop.
        # #######################

        # Make a DistanceMatrix in which to put the results
        d = DistanceMatrix()
        d.setDim(self.nTax)
        d.names = self.taxNames
        nUndefinedLogDets = 0

        for sNum1 in range(self.nTax - 1):
            firstSeq = seq[sNum1]
            for sNum2 in range(sNum1 + 1, self.nTax):
                secondSeq = seq[sNum2]

                # There are 3 possibilities only.  Double gaps might
                # include "gap-gap" or "n-gap" or "n-n".  Ambigs only
                # involve non-n-like ambiguous chars eg "r-gap", "r-a",
                # "r-n".

                # Initialize refUnambigCountMatrix and refAmbigCountMatrix
                # to zeros
                if 0:  # slow, in python
                    for i in range(self.dim):
                        for j in range(self.dim):
                            refUnambigCountMatrix[i, j] = 0
                    for i in range(bigDim):
                        for j in range(bigDim):
                            refAmbigCountMatrix[i, j] = 0
                else:  # fast, in c
                    pf.zeroNumPyInts(
                        refUnambigCountMatrix, (self.dim * self.dim))
                    pf.zeroNumPyInts(
                        refAmbigCountMatrix, (bigDim * bigDim))

                # print(f"xzz FillFxy. sNum1 {sNum1} {self.sequences[sNum1].name}  sNum2 {sNum2} {self.sequences[sNum2].name}  ") 
                if fastFillFxy:
                    pf.logDetFillFxy(nUnambig,
                                     nAmbig,
                                     nDoubleGap,
                                     seq,
                                     sNum1,
                                     sNum2,
                                     refUnambigCountMatrix,
                                     refAmbigCountMatrix,
                                     allSymbolNums,
                                     self.nTax,
                                     self.nChar,
                                     self.dim,
                                     bigDim,
                                     normUnambig,
                                     bigFxy,
                                     equatesArray)
                else:
                    nUnambig[0] = 0
                    nAmbig[0] = 0
                    nDoubleGap[0] = 0

                    # Fill the refUnambigCountMatrix and
                    # refAmbigCountMatrix arrays, based on the sequences.
                    for i in range(self.nChar):
                        firstChar = firstSeq[i]
                        secondChar = secondSeq[i]
                        if firstChar >= 0 and secondChar >= 0:
                            refUnambigCountMatrix[
                                firstChar, secondChar] += 1
                            nUnambig[0] += 1
                        else:
                            if firstChar == var.N_LIKE and secondChar == var.N_LIKE:
                                nDoubleGap[0] += 1
                            else:
                                firstIndx = 0
                                secondIndx = 0
                                for j1 in range(bigDim):
                                    if allSymbolNums[j1] == firstChar:
                                        firstIndx = j1
                                        break
                                for j1 in range(bigDim):
                                    if allSymbolNums[j1] == secondChar:
                                        secondIndx = j1
                                        break
                                refAmbigCountMatrix[
                                    firstIndx, secondIndx] += 1
                                nAmbig[0] += 1

                    if 0:
                        print("nUnambig=%s, refUnambigCountMatrix=" % nUnambig)
                        print(refUnambigCountMatrix)
                        print("nAmbig=%s, refAmbigCountMatrix=" % nAmbig)
                        print(refAmbigCountMatrix)
                        print("nChar=%s, nAmbig=%s, nDoubleGap=%s, nUnambig=%s" % (
                            self.nChar, nAmbig, nDoubleGap, nUnambig))
                    assert nAmbig[0] + nDoubleGap[0] + \
                        nUnambig[0] == self.nChar

                    # normUnambig is just the refUnambigCountMatrix array, normalized to sum
                    # to 1.0.  It is only used for reference purposes-- it
                    # does not change.
                    if 1:
                        normUnambig = refUnambigCountMatrix / \
                            float(nUnambig[0])
                        # print "normUnambig="
                        # print normUnambig
                    if 0:
                        for i in range(self.dim):
                            for j in range(self.dim):
                                normUnambig[i][j] = refUnambigCountMatrix[
                                    i][j] / float(nUnambig[0])
                        # print "normUnambig="
                        # print normUnambig

                    # Initialize the "bigFxy" array.  The bigFxy array starts
                    # out as a Float version of the refUnambigCountMatrix array.  Then it is
                    # ajusted upwards with ambigs, and then constant sites are
                    # removed from it.  Finally it is normilized so that it
                    # sums to 1.0, and then it is used to make the logDet.
                    for i in range(self.dim):
                        for j in range(self.dim):
                            try:
                                bigFxy[i, j] = float(
                                    refUnambigCountMatrix[i, j])
                            except:
                                print("xxx i=%i, j=%i" % (i, j))
                                print("xxx", refUnambigCountMatrix[i, j])
                                print("xxx", float(refUnambigCountMatrix[i, j]))
                                print("xxx", bigFxy[i, j])
                                raise P4Error(gm)

                    if nAmbig[0]:
                        # Long section on resolving ambiguities.

                        # print "equatesArray = "
                        # print equatesArray

                        for i in range(bigDim):
                            for j in range(bigDim):
                                na = refAmbigCountMatrix[i, j]
                                if na:
                                    # print "firstChar=%i, secondChar=%i" %
                                    # (i,j)
                                    fsum = 0.0
                                    # the number of possible combinations
                                    # of resolved ambigs.
                                    nSlots = 0
                                    # firstChar is a symbol
                                    if i < self.dim:
                                        firstChar = self.symbols[i]
                                        # print "a firstChar = %s" % firstChar
                                        # secondChar must be N_LIKE or an
                                        # equateSymbol
                                        if j == bigDim - 1:  # N_LIKE
                                            # print "a secondChar is
                                            # N_LIKE"
                                            for j1 in range(self.dim):
                                                fsum += normUnambig[i, j1]
                                                nSlots += 1
                                        elif j >= self.dim:  # an equate
                                            #secondChar = self.equateSymbols[j - self.dim]
                                            # print "a secondChar = %s" %
                                            # secondChar
                                            for j1 in range(self.dim):
                                                if equatesArray[j - self.dim, j1]:
                                                    # print "    %s %i" %
                                                    # (secondChar, j1)
                                                    fsum += normUnambig[i,
                                                                        j1]
                                                    nSlots += 1
                                        else:
                                            raise P4Error(
                                                "This shouldn't happen. Ambig site, but both chars are non-ambig.")
                                    # firstChar is N_LIKE
                                    elif i == bigDim - 1:
                                        # print "b firstChar is N_LIKE"
                                        # secondChar must be either a
                                        # symbol or an equate
                                        # secondChar is a symbol
                                        if j < self.dim:
                                            secondChar = self.symbols[j]
                                            # print "b secondChar = %s" %
                                            # secondChar
                                            for i1 in range(self.dim):
                                                fsum += normUnambig[i1, j]
                                                nSlots += 1
                                        elif j == bigDim - 1:
                                            raise P4Error(
                                                "This shouldn't happen.  Ambig site with 2 N_LIKEs.  Should be a double gap.")
                                        else:  # secondChar is an equate
                                            secondChar = self.equateSymbols[
                                                j - self.dim]
                                            # print "b secondChar = %s" %
                                            # secondChar
                                            for i1 in range(self.dim):
                                                for j1 in range(self.dim):
                                                    if equatesArray[j - self.dim, j1]:
                                                        # print "    %s %i"
                                                        # % (secondChar,
                                                        # j1)
                                                        fsum += normUnambig[
                                                            i1, j1]
                                                        nSlots += 1

                                    else:  # firstChar an equateSymbol
                                        firstChar = self.equateSymbols[
                                            i - self.dim]
                                        # print "c firstChar = %s" %
                                        # firstChar

                                        # secondChar could be anything
                                        # secondChar is a symbol
                                        if j < self.dim:
                                            secondChar = self.symbols[j]
                                            # print "c secondChar = %s" %
                                            # secondChar
                                            for i1 in range(self.dim):
                                                if equatesArray[i - self.dim, i1]:
                                                    fsum += normUnambig[
                                                        i1, j]
                                                    nSlots += 1

                                        # secondChar is N_LIKE
                                        elif j == bigDim - 1:
                                            # print "c secondChar is
                                            # N_LIKE"
                                            for i1 in range(self.dim):
                                                if equatesArray[i - self.dim, i1]:
                                                    # print "%s %i" %
                                                    # (firstChar, i1)
                                                    for j1 in range(self.dim):
                                                        fsum += normUnambig[
                                                            i1, j1]
                                                        nSlots += 1

                                        else:  # secondChar is an equate
                                            secondChar = self.equateSymbols[
                                                j - self.dim]
                                            # print "c secondChar = %s" %
                                            # secondChar
                                            for i1 in range(self.dim):
                                                if equatesArray[i - self.dim, i1]:
                                                    # print "%s %i" %
                                                    # (firstChar, i1)
                                                    for j1 in range(self.dim):
                                                        if equatesArray[j - self.dim, j1]:
                                                            # print "    %s
                                                            # %i" %
                                                            # (secondChar,
                                                            # j1)
                                                            fsum += normUnambig[
                                                                i1, j1]
                                                            nSlots += 1

                                    # print "fsum=%f, nSlots=%i" % (fsum,
                                    # nSlots)

                                    if fsum == 0.0:
                                        # print "bigFxy= "
                                        # print bigFxy
                                        # print "nSlots = %i" % nSlots
                                        oneOverNSlots = 1.0 / nSlots
                                        # firstChar is a symbol
                                        if i < self.dim:
                                            # secondChar must be N_LIKE or
                                            # an equateSymbol
                                            if j == bigDim - 1:  # N_LIKE
                                                for j1 in range(self.dim):
                                                    bigFxy[
                                                        i, j1] += na * oneOverNSlots
                                            # an equate
                                            elif j >= self.dim:
                                                for j1 in range(self.dim):
                                                    if equatesArray[j - self.dim, j1]:
                                                        bigFxy[
                                                            i, j1] += na * oneOverNSlots
                                            else:
                                                raise P4Error(
                                                    "This shouldn't happen")
                                        # firstChar is N_LIKE
                                        elif i == bigDim - 1:
                                            # secondChar must be either a
                                            # symbol or an equate
                                            # secondChar is a symbol
                                            if j < self.dim:
                                                for i1 in range(self.dim):
                                                    bigFxy[
                                                        i1, j] += na * oneOverNSlots
                                            elif j == bigDim - 1:
                                                raise P4Error(
                                                    "This shouldn't happen.")
                                            # secondChar is an equate
                                            else:
                                                for i1 in range(self.dim):
                                                    for j1 in range(self.dim):
                                                        if equatesArray[j - self.dim, j1]:
                                                            bigFxy[
                                                                i1, j1] += na * oneOverNSlots

                                        else:  # firstChar an equateSymbol
                                            # secondChar could be anything
                                            # secondChar is a symbol
                                            if j < self.dim:
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        bigFxy[
                                                            i1, j] += na * oneOverNSlots
                                            # secondChar is N_LIKE
                                            elif j == bigDim - 1:
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        for j1 in range(self.dim):
                                                            bigFxy[
                                                                i1, j1] += na * oneOverNSlots

                                            # secondChar is an equate
                                            else:
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        for j1 in range(self.dim):
                                                            if equatesArray[j - self.dim, j1]:
                                                                bigFxy[
                                                                    i1, j1] += na * oneOverNSlots

                                    else:  # fsum is not zero
                                        # firstChar is a symbol
                                        if i < self.dim:
                                            # secondChar must be N_LIKE or
                                            # an equateSymbol
                                            if j == bigDim - 1:  # N_LIKE
                                                for j1 in range(self.dim):
                                                    bigFxy[
                                                        i, j1] += na * normUnambig[i, j1] / fsum
                                            # an equate
                                            elif j >= self.dim:
                                                for j1 in range(self.dim):
                                                    if equatesArray[j - self.dim, j1]:
                                                        bigFxy[
                                                            i, j1] += na * normUnambig[i, j1] / fsum
                                            else:
                                                raise P4Error(
                                                    "This shouldn't happen")
                                        # firstChar is N_LIKE
                                        elif i == bigDim - 1:
                                            # secondChar must be either a
                                            # symbol or an equate
                                            # secondChar is a symbol
                                            if j < self.dim:
                                                for i1 in range(self.dim):
                                                    bigFxy[
                                                        i1, j] += na * normUnambig[i1, j] / fsum
                                            elif j == bigDim - 1:
                                                raise P4Error(
                                                    "This shouldn't happen.")
                                            # secondChar is an equate
                                            else:
                                                for i1 in range(self.dim):
                                                    for j1 in range(self.dim):
                                                        if equatesArray[j - self.dim, j1]:
                                                            bigFxy[
                                                                i1, j1] += na * normUnambig[i1, j1] / fsum

                                        else:  # firstChar an equateSymbol
                                            # secondChar could be anything
                                            # secondChar is a symbol
                                            if j < self.dim:
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        bigFxy[
                                                            i1, j] += na * normUnambig[i1, j] / fsum
                                            # secondChar is N_LIKE
                                            elif j == bigDim - 1:
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        for j1 in range(self.dim):
                                                            bigFxy[
                                                                i1, j1] += na * normUnambig[i1, j1] / fsum

                                            # secondChar is an equate
                                            else:
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        for j1 in range(self.dim):
                                                            if equatesArray[j - self.dim, j1]:
                                                                bigFxy[
                                                                    i1, j1] += na * normUnambig[i1, j1] / fsum
                                    if 0:
                                        print("bigFxy= (after partial ambig resolution)")
                                        print(bigFxy)

                    # End of the long section on resolving ambiguities.
                # End of the long "else" clause to "if fastFillFxy:"

                if 0:
                    print(f"bigFxy=  (sNum1 {sNum1} sNum2 {sNum2} after ambig resolution)")
                    print(bigFxy)

                # pInvar stuff
                if pInvar != None:
                    nSitesToRemove = pInvar * self.nChar
                    for i in range(self.dim):
                        bigFxy[i][i] -= alignmentComposition[i] * nSitesToRemove

                if 0:
                    print(f"bigFxy=  (sNum1 {sNum1} sNum2 {sNum2} after pInvar removal)")
                    print(bigFxy)

                if missingCharacterStrategy == 'fudge':
                    # Replace zeros on the diagonal of Fxy with either 0.5
                    # or half of the smallest positive Fii, whichever is
                    # smaller.
                    minPositiveOnDiag = 1.0
                    for i in range(self.dim):
                        if bigFxy[i, i] > 0.0 and bigFxy[i, i] < minPositiveOnDiag:
                            minPositiveOnDiag = bigFxy[i, i]
                    minPositiveOnDiag /= 2.0
                    for i in range(self.dim):
                        if bigFxy[i, i] < minPositiveOnDiag:
                            if verbose:
                                print(f"zyx missingCharacterStrategy is 'fudge', sNum1 {sNum1} sNum2 {sNum2} ", end='')
                                print(f"char {i}({self.symbols[i]}) was {bigFxy[i, i]}, too small, set to {minPositiveOnDiag}")
                            bigFxy[i, i] = minPositiveOnDiag
                            fudgeCount += 1
                    theFxy = bigFxy

                else:
                    theFxy = bigFxy

                sumTheFxy = sum(sum(theFxy))

                # Normalize Fxy to 1.0
                theFxy /= float(sumTheFxy)
                if 0:
                    print("nSites=%f,  final normalized bigFxy (ie Fxy) = " % sumTheFxy)
                    print(theFxy)

                # Calculate the logDet, from the theFxy and from the bigPi
                theDet = numpy.linalg.det(theFxy)
                # print("theDet from theFxy is %f" % theDet)

                if theDet > 0.0:
                    bigPiX = theFxy.sum(axis=0)   # sum of columns
                    bigPiY = theFxy.sum(axis=1)   # sum of rows
                    if 0:
                        print("theFxy = ")
                        print(theFxy, type(theFxy))
                        print("bigPiX = ")
                        print(bigPiX, type(bigPiX))
                        print("bigPiY =")
                        print(bigPiY, type(bigPiY))
                    det_bigPiX = numpy.multiply.reduce(bigPiX)
                    det_bigPiY = numpy.multiply.reduce(bigPiY)
                    # print det_bigPiX, det_bigPiY
                    if det_bigPiX <= 0.0 or det_bigPiY <= 0.0:
                        if 1:
                            if det_bigPiX <= 0.0:
                                print(f"sNum1 {sNum1} sNum2 {sNum2} det_bigPiX {det_bigPiX} ")
                            if det_bigPiY <= 0.0:
                                print(f"sNum1 {sNum1} sNum2 {sNum2} det_bigPiY {det_bigPiY} ")
                        if missingCharacterStrategy == 'refuse':
                            if verbose:
                                print("zxx detPi is non-positive, missingCharacterStrategy is 'refuse', so log det is undefined.")
                            d.matrix[sNum1][sNum2] = -1.0
                            d.matrix[sNum2][sNum1] = -1.0
                            nUndefinedLogDets += 1
                        else:
                            if 1:
                                gm.append("sumTheFxy = %f" % sumTheFxy)
                                gm.append("sum(bigPiX)=%s, sum(bigPiY)=%s" % (numpy.sum(bigPiX), numpy.sum(bigPiY)))
                                gm.append("symbols = %s" % self.symbols)
                                gm.append("sNum1=%i, sNum2=%i" % (sNum1, sNum2))
                                gm.append("bigPiX = %s" % bigPiX)
                                gm.append("bigPiY=%s" % bigPiY)
                                gm.append("det_bigPiX = %s, det_bigPiY=%s" % (det_bigPiX, det_bigPiY))
                                gm.append("Got bad Pi det due to missing char(s).")
                                gm.append("This should not happen-- programming error.")
                                raise P4Error(gm)


                    else:  # det_bigPiX and det_bigPiY are over zero
                        # If we have been using a reduced Fxy, then
                        # self.dim is no longer appropriate.
                        reducedDim = len(bigPiX)

                        if 0:
                            # This section works, but is not very clear.
                            # Re-written below
                            theLogDet = numpy.log(theDet) - 0.5 * numpy.log(det_bigPiX) - 0.5 * numpy.log(det_bigPiY)

                            if correction == 'L94':
                                theCorrection = reducedDim
                            elif correction == 'TK02':  # from LDDist
                                squareSum = 0.0
                                for i in range(reducedDim):
                                    squareSum += (bigPiX[i] + bigPiY[i]
                                                  ) * (bigPiX[i] + bigPiY[i])
                                theCorrection = (
                                    reducedDim - 1) / (1.0 - squareSum / 4.0)
                            if 0:
                                print("theDet = %g" % theDet)
                                print("theLogDet = %f" % theLogDet)
                                print("theCorrection = %s" % theCorrection)

                            theLogDet /= theCorrection
                            theLogDet = -theLogDet

                        if 1:
                            # equation 3, pg 606, in L94
                            if correction == 'L94':
                                theLogDet = -numpy.log(theDet) + (numpy.log(det_bigPiX) + numpy.log(det_bigPiY)) / 2.0
                                theLogDet /= reducedDim
                                # print(f"sNum1 {sNum1} sNum2 {sNum2} theLogDet after L94 correction {theLogDet}") 
                            # equation 11, page 1729, in TK02
                            elif correction == 'TK02':
                                theLogDet = numpy.log(theDet) - (0.5 * (numpy.log(det_bigPiX) + numpy.log(det_bigPiY)))
                                squareSum = 0.0
                                for i in range(reducedDim):
                                    thePi = (bigPiX[i] + bigPiY[i]) / 2.0
                                    squareSum += thePi * thePi
                                theLogDet = - ((1.0 - squareSum) / (reducedDim - 1)) * theLogDet
                            elif correction == 'TK02_eqn10':
                                theLogDet = - (1. / reducedDim) * numpy.log(theDet) - numpy.log(reducedDim)

                        # dset allsitesmean=yes
                        #if doPInvarOfConstants == True and pInvarOfConstants != None:
                        #    theLogDet *= 1.0 - (pInvarOfConstants * nConstants) / self.nChar
                        if pInvar != None:
                            theLogDet *= 1.0 - pInvar

                        if theLogDet <= 0.0:
                            # gm = ["Alignment.logDet()"]
                            # gm.append(f"Between (zero-based) sNum1 {sNum1} {self.sequences[sNum1].name}")
                            # gm.append(f"and sNum2 {sNum2} {self.sequences[sNum2].name}") 
                            # gm.append(f"Got negative logDet ({theLogDet}). ")
                            # print(gm)
                            # #raise P4Error(gm)
                            if verbose:
                                print(f"xxy sNum1 {sNum1} sNum2 {sNum2} Got non-positive corrected logDet. {theLogDet}")
                            d.matrix[sNum1][sNum2] = -1.0
                            d.matrix[sNum2][sNum1] = -1.0
                            nUndefinedLogDets += 1
                        else:
                            d.matrix[sNum1][sNum2] = theLogDet
                            d.matrix[sNum2][sNum1] = theLogDet
                else: # theDet is less than or equal to zero

                    if verbose:
                        print(f"xxx sNum1 {sNum1} sNum2 {sNum2} Got non-positive logDet {theDet} on initial calculation")
                    d.matrix[sNum1][sNum2] = -1.0
                    d.matrix[sNum2][sNum1] = -1.0
                    nUndefinedLogDets += 1

        # End of the main pairwise loop

        dMessage = ["    "]
        dMessage.append("Log det distances from p4.")
        dMessage.append("Correction from %s" % correction)
        dMessage.append("The missingCharacterStrategy is set to '%s'." % missingCharacterStrategy)
        if missingCharacterStrategy == 'fudge':
            dMessage.append("    Did %i fudges." % fudgeCount)
        if nUndefinedLogDets:
                dMessage.append(f"There were {nUndefinedLogDets} undefined log dets")
                dMessage.append(f"given as -1.0 in the distance matrix.  You will want to deal with those.")
        d.message = '\n    '.join(dMessage)

        print(d.message)

        return d

