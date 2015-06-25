from Var import var
from DistanceMatrix import DistanceMatrix
from Glitch import Glitch
import numpy, numpy.linalg
import math,string
#import func  # temp
import pf

def logDet(self, correction='TK02', doPInvarOfConstants=True, pInvar=None, pInvarOfConstants=None, missingCharacterStrategy='fudge', minCompCount=1, nonPositiveDetStrategy='invert'):
    """logDet calculations.  Returns a DistanceMatrix object.

    Heavily influenced by paup and by LDDist.  Thanks to Swofford and
    Thollesson, the authors.

    The basis of the log det is the Fxy matrix, the matrix of observed
    changes between 2 sequences.  It is from the Fxy that the
    determinant is calculated.  The negative of the log of the
    determinant ('logdet') is related to the evolutionary distance
    separating the 2 sequences.

    **Corrections**
    
    The logdet can be corrected so that it expresses the distance
    between 2 sequences in the usual terms of mutations per site.
    There are 2 different corrections-- that from Lockhart94 Eqn3 (pg
    606), and Tamura and Kumar 2002.  Which one is given by the arg
    'correction', which should be one of 'L94' or 'TK02'.  The former
    is used in paup, and the latter is used in LDDist.

    **Ambiguities**
    

    In sequence pairs, when there are 2 N-like chars (n, gap, ?), they
    are completely ignored.  Other combinations are not ignored.
    Consider this alignment::

        acgtr
        acgta

    Now you might think that r-a would mean either a-a or g-a, but in
    the present implementation it is slightly different-- the
    ambiguity is added to Fxy in proportion to the observed
    unambiguous counts.  So in this case there was an observed a-a,
    but no g-a, so r-a is considered to be wholly a-a.  Now consider
    the alignment::

        aagcgtr
        aaacgta

    Here r-a is considered to be both a-a and g-a, in proportion to
    their occurrence, ie 2/3 a-a and 1/3 g-a.

    **pInvar**
    

    We can deal with pInvar in a couple of ways.  In both ways, we
    subtract some numbers from the diagonal of the un-normalized Fxy
    matrix, in proportion to the composition of the constant sites.
    The paup-like way is simply to remove a pInvar * nSites * comp[i]
    from each un-normalized Fii.  The LDDist-like way is to identify
    which sites are constant, and then to remove some
    pInvarOfConstants from the Fii.  Do the former by setting
    doPInvarOfConstants=False, and do the latter by setting
    doPInvarOfConstants=True.  To explain the difference,
    consider this alignment of 10 chars and 3 sequences::

        A acgtacgtta 
        B acgtacgttg
        C acgtacgtag
          ++++++++--  constant

    The un-normalized Fxy matrix for the A-B pair is::
    
       [[ 2.  0.  1.  0.]
        [ 0.  2.  0.  0.]
        [ 0.  0.  2.  0.]
        [ 0.  0.  0.  3.]]

    There are 8 constant sites, and those constant sites have equal
    base freqs.  Say I want pInvar=0.8, and I use a paup-like strategy
    by setting doPInvarOfConstantSites=False.  So I remove 0.8 * 10 *
    0.25 = 2. from each un-normalized Fii.  That leaves zeros on the
    diagonal and so it fails.  I can have pInvar=0.79 and that works,
    but if I try pInvar=0.8 or more than that then it fails.

    Now for more realistic sequences the Fii will be more than the
    pInvar*nSites*comp[i] because there will be some Fii contribution
    from non-constant sites; so often the Fii won't go to zero.  When
    that is the case, you will be able to successfully set pInvar to
    be equal to or even more than the proportion of constant sites.
    It seems strange to me to be able to set the pInvar higher than
    the proportion of constant sites, but paup allows you to do that,
    and p4 imitates it when doPInvarOfConstants=False.

    The LDDist-like strategy is to first identify the constant sites
    in the alignment, and only allow removal of that many counts from
    the un-normalized Fii.  I do that by setting
    doPInvarOfConstants=True, and I set pInvarOfConstants from 0-1.0.
    So in the alignment above, only if I set pInvarOfConstants=1.0 I
    will get some Fii = 0, but if pInvarOfConstants is anything less
    than 1.0 its ok. (For more realistic alignments there should be no
    trouble with pInvarOfConstants=1.0, as some Fii will come from
    non-constant sites).

    **Constant sites**
    

    Here constant sites are defined in a very simple way, different
    from the Alignment.constantSites*() methods.  Here, a site is
    constant if it is constant and unambiguous; if it has any ambigs
    or gaps then it is not constant.  In the
    Alignment.constantSites*() methods, a site is constant if it could
    possibly be constant in any resolution of ambigs or gaps.

    It appears that paup has a slightly different definition of a
    constant site, which comes in to play when there are ambigs.  It
    also appears that paup has a different way of calculating the
    composition of those constant sites.

    **Dealing with missing (or low-frequency) characters**
    

    If any character in any pairwise comparison in either sequence is
    missing then the calculation fails because of the correction. (The
    correction involves the log of the product of all the
    compositions; if any comp is zero then the product will be zero,
    and the log will be undefined --- Boom!).  There are 3 strategies
    implemented for dealing with missing chars--

    1.  refuse
                 This is the way paup does it.  Refuse to calculate.
                 Make up a 'big distance' later.  By default in paup,
                 the 'big distance' is twice the biggest defined
                 distance (ie defined elsewhere in the distance
                 matrix, from some other pairwise comparison).  Of
                 course this will fail if all the pairs are refused,
                 leaving no defined distances.
    2.  fudge
                 This is the way that LDDist does it.  LDDist replaces
                 zeros with 0.5 in the un-normalized Fxy martrix.  The
                 present implementation replaces zeros on the diagonal
                 of Fxy with either 0.5 or half of the smallest
                 positive Fii, whichever is smaller.
    3.  reduce
                 Robert Hirt's idea is that if the comp of any char in
                 either of the pair of sequences is zero or very
                 small, then probably all of the Fxy data that
                 involve that particular character are going to
                 unreliable, and so should not be used.  So in this
                 case the entire row and entire column for that
                 character is removed if the comp is less than the arg
                 'minCompCount' (default is 1).  If a reduced Fxy
                 matrix, made by removal of a certain row-column, is
                 needed for one pairwise comparison, then that
                 row-column is removed from all pairwise compares in
                 the alignment, even if the char in that row-column is
                 in high frequency in other sequence pairs.

    **Non-positive determinants**
    

    Sometimes for very diverged sequences the determinant of the Fxy
    matrix is zero or less, which means that the logarithm is
    undefined.  (This can happen even if all Fii are positive).  There
    are a couple of strategies for dealing with that.

    1.  refuse :: This is the way that paup does it.  Same as above.
    2.  invert ::  This is the way that LDDist does it.  If you get a
                negative determinant, the sign is simply inverted.
                This is a hack, without any reason that I can see
                except to make the distance calculable.  Fortunately
                when the det goes below zero, it only goes a little
                below zero, and so it does not matter much.  Now this
                strategy of course won't do anything for zero dets
                (for which LDDist throws an error), so for zero dets I
                simply make it a very small positive number.

    
    """

    complaintHead = 'Alignment.logDet()'
    gm = [complaintHead]

    fastFillFxy = True # in c

    # Check that there is only one partition
    if not self.nexusSets or not self.nexusSets.charPartition: # its all one part
        pass
    elif self.nexusSets.charPartition and len(self.nexusSets.charPartition.subsets) > 1: # its partitioned.  Bad.
        gm = [complaintHead]
        gm.append("This only works with Alignments having only one data partition.")
        raise Glitch, gm

    # Check the corrections arg
    goodCorrections = ['L94', 'TK02', 'TK02_eqn10']
    if correction not in goodCorrections:
        gm.append("The corrections arg should be one of: %s" % goodCorrections)
        gm.append("Got %s" % correction)
        raise Glitch, gm

    # Check the doPInvarOfConstants arg
    if doPInvarOfConstants not in [True, False]:
        gm.append("doPInvarOfConstants should be set to either True or False")
        raise Glitch, gm

    # Check the pInvar or pInvarOfConstants args.  If zero, set to None.
    if doPInvarOfConstants:
        if pInvar:
            gm.append("doPInvarOfConstants is set, which means that pInvar does not apply.")
            gm.append("To prove that you are not mixed up, set it to None.")
            raise Glitch, gm
        try:
            pInvarOfConstants = float(pInvarOfConstants)
            if math.fabs(pInvarOfConstants) < 1.0e-10:
                pInvarOfConstants = None
        except ValueError:
            pInvarOfConstants = None
        except TypeError:
            pInvarOfConstants = None
        if pInvarOfConstants and (pInvarOfConstants < 0.0 or pInvarOfConstants > 1.0):
            gm.append("pInvarOfConstants, if set, should be between zero and 1.0, inclusive.")
            raise Glitch, gm
    else:
        if pInvarOfConstants:
            gm.append("doPInvarOfConstants is off, which means that pInvarOfConstants does not apply.")
            gm.append("To prove that you are not mixed up, set it to None.")
            raise Glitch, gm
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
            raise Glitch, gm

    # Check the missingCharacterStrategy arg
    goodMissingCharacterStrategies = ['refuse', 'fudge', 'reduce']
    if missingCharacterStrategy in goodMissingCharacterStrategies:
        pass
    else:
        gm.append("Arg missingCharacterStrategy should be one of %s" % goodMissingCharacterStrategies)
        gm.append("Got %s" % missingCharacterStrategy)
        raise Glitch, gm
    
    # Check the nonPositiveDetStrategy arg
    goodNonPositiveDetStrategies = ['refuse', 'invert']
    if nonPositiveDetStrategy in goodNonPositiveDetStrategies:
        pass
    else:
        gm.append("Arg nonPositiveDetStrategy should be one of %s" % goodNonPositiveDetStrategies)
        gm.append("Got %s" % nonPositiveDetStrategy)
        raise Glitch, gm
    
    # Equates
    if self.equates:
        self.equateSymbols = self.equates.keys()
        self.equateSymbols.sort()

        # Make equatesArray
        equatesArray = numpy.zeros((self.nEquates, self.dim), numpy.int32)
        for i in range(self.nEquates):
            e = self.equates[self.equateSymbols[i]]
            for j in range(self.dim):
                s = self.symbols[j]
                if s in e:
                    equatesArray[i,j] = 1
        #print "equatesArray = "
        #print equatesArray
    else:
        equatesArray = numpy.array([-1], numpy.int32)

    # We need a charsLikeN list, showing characters that are fully ambiguous.
    charsLikeN = []
    if self.equates:
        for i in range(self.nEquates):
            e = self.equateSymbols[i]
            isNLike = True
            for j in range(self.dim):
                if not equatesArray[i,j]:
                    isNLike = False
                    break
            if isNLike:
                charsLikeN.append(e)
    charsLikeN.append('-')
    charsLikeN.append('?')
    #print "charsLikeN = %s" % charsLikeN

    # Make a seq, a numpy array, in which to hold the recoded sequences. 
    seq = numpy.zeros((self.nTax, self.nChar), numpy.int8)
    for sNum in range(self.nTax):
        s = self.sequences[sNum]
        for cNum in range(self.nChar):
            theChar = s.sequence[cNum]
            if theChar in charsLikeN:
                seq[sNum, cNum] = var.N_LIKE
            elif self.equates and theChar in self.equates:
                seq[sNum, cNum] = var.EQUATES_BASE + self.equateSymbols.index(theChar)
            else:
                seq[sNum, cNum] = self.symbols.index(theChar)
    #print seq

    constComps = None
    constCounts = None
    # Deal with pInvar or pInvarOfConstants, if needed.
    #print "pInvarOfConstants=%s, %s" % (pInvarOfConstants, pInvarOfConstants != None)
    if (doPInvarOfConstants==False and pInvar != None) or (doPInvarOfConstants==True and pInvarOfConstants != None):
        # If we are going to do something with constants, then we need
        # to know what sites are constant, if only to get the
        # composition.  If the site contains any gaps or ambigs then
        # it is not constant.  Only constant unambigs are constant.
        constants = numpy.ones(self.nChar, numpy.int32)
        constCounts = numpy.zeros(self.dim, numpy.int32)
        for j in range(self.nChar):
            theRefChar = seq[0,j]
            if theRefChar < 0:
                constants[j] = 0
            else:
                for i in range(1,self.nTax):
                    theChar = seq[i,j]
                    if theRefChar != theChar:
                        constants[j] = 0
                        break
        for j in range(self.nChar):
            if constants[j]:
                #print "j=%i, seq[0][j]=%s" % (j, seq[0][j])
                constCounts[seq[0][j]] +=1
        nConstants = numpy.sum(constants)
        sumOfCounts = float(numpy.sum(constCounts))
        if sumOfCounts:
            constComps = constCounts / float(numpy.sum(constCounts))
        else:
            constComps = None
        if 0:
            print "constants =", constants
            print "nConstants = %i" % nConstants
            print "constCounts = %s" % constCounts
            print "constComps = %s" % constComps
        # not used further
        del(constants)

    # The 'refUnambigCountMatrix' array is raw counts of changes between the two sequences.
    refUnambigCountMatrix = numpy.zeros((self.dim, self.dim), numpy.int32)
    normUnambig = numpy.zeros((self.dim, self.dim), numpy.float)
    allSymbolNums = range(self.dim)
    allSymbolNums += range(var.EQUATES_BASE, var.EQUATES_BASE + self.nEquates)
    allSymbolNums.append(var.N_LIKE)
    allSymbolNums = numpy.array(allSymbolNums, numpy.int32)
    #print "allSymbolNums = ", allSymbolNums
    bigDim = self.dim + 1 + len(self.equates)
    #print "bigDim = %i" % bigDim
    refAmbigCountMatrix = numpy.zeros((bigDim, bigDim), numpy.int32) # Int, to keep track of ambiguous changes
    #normUnambig = zeros((self.dim, self.dim), Float) # 'double' type
    bigFxy = numpy.zeros((self.dim, self.dim), numpy.float)
    fudgeCount = 0
    invertCount = 0


    if missingCharacterStrategy == 'reduce':
        # When this is set, I may use a reduced matrix, where low-freq
        # chars in a pair are ignored.  But if a pair has a low-freq
        # character, then that character should also be ignored in all the
        # other pairs (even if its at high freq, and so would otherwise be
        # ok).  So here we ask for what chars in all pairs need to be
        # ignored.  The answer is in the ignores array, so for DNA if the
        # ignores array ends up being [0,1,0,0], then we us a 3x3 Fxy and
        # ignore all sites containing a C.

        # The slow version can be very slow.
        hasIgnores = False
        if minCompCount:
            # slow
            
            #ignores = self._logDetSetReduceIgnores_slow(doPInvarOfConstants, pInvar, pInvarOfConstants, minCompCount, seq, constComps, constCounts, refUnambigCountMatrix, allSymbolNums, bigDim, refAmbigCountMatrix, bigFxy)
            ignores = self._logDetSetReduceIgnores(doPInvarOfConstants, pInvar, pInvarOfConstants, minCompCount, seq, constComps, constCounts)
            if numpy.sum(ignores):
                hasIgnores = True
            if 0:
                print "hasIgnores=%s" % hasIgnores
                print 'ignores= %s' % ignores
            totalNoIgnores = 0
            if hasIgnores:
                for i in range(self.dim):
                    if not ignores[i]:
                        totalNoIgnores += 1
                #print "totalNoIgnores = %i" % totalNoIgnores
                if totalNoIgnores < 2:
                    gm.append("The arg 'missingCharacterStrategy' is set to 'reduce'")
                    gm.append("The arg 'minCompCount' is turned on, and set to %i." % minCompCount)
                    gm.append("There is not enough variation in these sequences to make a valid distance.")
                    gm.append("There are too many sites that will be ignored because of low frequency characters.")
                    raise Glitch, gm
            if 0 and hasIgnores:
                print "missingCharacterStrategy is set to 'reduce'"
                print "The following char(s) will be ignored:"
                for i in range(self.dim):
                    if ignores[i]:
                        print " %s" % self.symbols[i],
                print "\n"


    nUnambig = numpy.zeros((1), numpy.int32)
    nAmbig = numpy.zeros((1), numpy.int32)
    nDoubleGap = numpy.zeros((1), numpy.int32)

    # The main pairwise loop.  If missingCharacterStrategy is
    # 'reduce', we may have to start it all over again if new ignores
    # are found to be needed.
    hasNewIgnores = True
    while hasNewIgnores:
        hasNewIgnores = False

        # Make a DistanceMatrix in which to put the results
        d = DistanceMatrix()
        d.setDim(self.nTax)
        d.names = self.taxNames
        nUndefinedLogDets = 0

        for sNum1 in range(self.nTax - 1):
            if hasNewIgnores:
                break
            firstSeq = seq[sNum1]
            for sNum2 in range(sNum1 + 1, self.nTax):
                if hasNewIgnores:
                    break
                secondSeq = seq[sNum2]

                # There are 3 possibilities only.  Double gaps might
                # include "gap-gap" or "n-gap" or "n-n".  Ambigs only
                # involve non-n-like ambiguous chars eg "r-gap", "r-a",
                # "r-n".

                # Initialize refUnambigCountMatrix and refAmbigCountMatrix to zeros
                if 0: # slow, in python
                    for i in range(self.dim):
                        for j in range(self.dim):
                            refUnambigCountMatrix[i,j] = 0
                    for i in range(bigDim):
                        for j in range(bigDim):
                            refAmbigCountMatrix[i,j] = 0
                else: # fast, in c
                    pf.zeroNumPyInts(refUnambigCountMatrix, (self.dim * self.dim))
                    pf.zeroNumPyInts(refAmbigCountMatrix, (bigDim * bigDim))



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

                    # Fill the refUnambigCountMatrix and refAmbigCountMatrix arrays, based on the sequences.
                    for i in range(self.nChar):
                        firstChar = firstSeq[i]
                        secondChar = secondSeq[i]
                        if firstChar >= 0 and secondChar >=0:
                            refUnambigCountMatrix[firstChar, secondChar] += 1
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
                                refAmbigCountMatrix[firstIndx, secondIndx] += 1
                                nAmbig[0] += 1

                    if 0:
                        print "nUnambig=%s, refUnambigCountMatrix=" % nUnambig
                        print refUnambigCountMatrix
                        print "nAmbig=%s, refAmbigCountMatrix=" % nAmbig
                        print refAmbigCountMatrix
                        print "nChar=%s, nAmbig=%s, nDoubleGap=%s, nUnambig=%s" % (
                            self.nChar, nAmbig, nDoubleGap, nUnambig)
                    assert nAmbig[0]+nDoubleGap[0]+nUnambig[0] == self.nChar

                    # normUnambig is just the refUnambigCountMatrix array, normalized to sum
                    # to 1.0.  It is only used for reference purposes-- it
                    # does not change.
                    if 1:
                        normUnambig = refUnambigCountMatrix / float(nUnambig[0])
                        #print "normUnambig="
                        #print normUnambig
                    if 0:
                        for i in range(self.dim):
                            for j in range(self.dim):
                                normUnambig[i][j] = refUnambigCountMatrix[i][j] / float(nUnambig[0])
                        #print "normUnambig="
                        #print normUnambig

                    # Initialize the "bigFxy" array.  The bigFxy array starts
                    # out as a Float version of the refUnambigCountMatrix array.  Then it is
                    # ajusted upwards with ambigs, and then constant sites are
                    # removed from it.  Finally it is normilized so that it
                    # sums to 1.0, and then it is used to make the logDet.
                    for i in range(self.dim):
                        for j in range(self.dim):
                            try:
                                bigFxy[i,j] = float(refUnambigCountMatrix[i,j])
                            except:
                                print "xxx i=%i, j=%i" % (i,j)
                                print "xxx", refUnambigCountMatrix[i,j]
                                print "xxx", float(refUnambigCountMatrix[i,j])
                                print "xxx", bigFxy[i,j]
                                raise Glitch, gm


                    if nAmbig[0]:
                        # Long section on resolving ambiguities.

                        #print "equatesArray = "
                        #print equatesArray

                        for i in range(bigDim):
                            for j in range(bigDim):
                                na = refAmbigCountMatrix[i,j]
                                if na:
                                    #print "firstChar=%i, secondChar=%i" % (i,j)
                                    fsum = 0.0
                                    nSlots = 0  # the number of possible combinations of resolved ambigs.
                                    if i < self.dim: # firstChar is a symbol
                                        firstChar = self.symbols[i]
                                        #print "a firstChar = %s" % firstChar
                                        # secondChar must be N_LIKE or an equateSymbol
                                        if j == bigDim - 1:  # N_LIKE
                                            #print "a secondChar is N_LIKE"
                                            for j1 in range(self.dim):
                                                fsum += normUnambig[i, j1]
                                                nSlots += 1
                                        elif j >= self.dim:  # an equate
                                            #secondChar = self.equateSymbols[j - self.dim]
                                            #print "a secondChar = %s" % secondChar
                                            for j1 in range(self.dim):
                                                if equatesArray[j - self.dim, j1]:
                                                    #print "    %s %i" % (secondChar, j1)
                                                    fsum += normUnambig[i, j1]
                                                    nSlots += 1
                                        else:
                                            raise Glitch, "This shouldn't happen. Ambig site, but both chars are non-ambig."
                                    elif i == bigDim - 1: # firstChar is N_LIKE
                                        #print "b firstChar is N_LIKE"
                                        # secondChar must be either a symbol or an equate
                                        if j < self.dim: # secondChar is a symbol
                                            secondChar = self.symbols[j]
                                            #print "b secondChar = %s" % secondChar
                                            for i1 in range(self.dim):
                                                fsum += normUnambig[i1, j]
                                                nSlots += 1
                                        elif j == bigDim - 1:
                                            raise Glitch, "This shouldn't happen.  Ambig site with 2 N_LIKEs.  Should be a double gap."
                                        else: # secondChar is an equate
                                            secondChar = self.equateSymbols[j - self.dim]
                                            #print "b secondChar = %s" % secondChar
                                            for i1 in range(self.dim):
                                                for j1 in range(self.dim):
                                                    if equatesArray[j - self.dim, j1]:
                                                        #print "    %s %i" % (secondChar, j1)
                                                        fsum += normUnambig[i1, j1]
                                                        nSlots += 1

                                    else:  # firstChar an equateSymbol
                                        firstChar = self.equateSymbols[i - self.dim]
                                        #print "c firstChar = %s" % firstChar

                                        # secondChar could be anything
                                        if j < self.dim: # secondChar is a symbol
                                            secondChar = self.symbols[j]
                                            #print "c secondChar = %s" % secondChar
                                            for i1 in range(self.dim):
                                                if equatesArray[i - self.dim, i1]:
                                                    fsum += normUnambig[i1, j]
                                                    nSlots += 1

                                        elif j == bigDim - 1:  # secondChar is N_LIKE
                                            #print "c secondChar is N_LIKE"
                                            for i1 in range(self.dim):
                                                if equatesArray[i - self.dim, i1]:
                                                    #print "%s %i" % (firstChar, i1)
                                                    for j1 in range(self.dim):
                                                        fsum += normUnambig[i1,j1]
                                                        nSlots += 1

                                        else:  # secondChar is an equate
                                            secondChar = self.equateSymbols[j - self.dim]
                                            #print "c secondChar = %s" % secondChar
                                            for i1 in range(self.dim):
                                                if equatesArray[i - self.dim, i1]:
                                                    #print "%s %i" % (firstChar, i1)
                                                    for j1 in range(self.dim):
                                                        if equatesArray[j - self.dim, j1]:
                                                            #print "    %s %i" % (secondChar, j1)
                                                            fsum += normUnambig[i1,j1]
                                                            nSlots += 1

                                    #print "fsum=%f, nSlots=%i" % (fsum, nSlots)

                                    if fsum == 0.0:
                                        #print "bigFxy= "
                                        #print bigFxy
                                        #print "nSlots = %i" % nSlots
                                        oneOverNSlots = 1.0 / nSlots
                                        if i < self.dim: # firstChar is a symbol
                                            # secondChar must be N_LIKE or an equateSymbol
                                            if j == bigDim - 1:  # N_LIKE
                                                for j1 in range(self.dim):
                                                    bigFxy[i,j1] += na * oneOverNSlots
                                            elif j >= self.dim:  # an equate
                                                for j1 in range(self.dim):
                                                    if equatesArray[j - self.dim, j1]:
                                                        bigFxy[i,j1] += na * oneOverNSlots
                                            else:
                                                raise Glitch, "This shouldn't happen"
                                        elif i == bigDim - 1: # firstChar is N_LIKE
                                            # secondChar must be either a symbol or an equate
                                            if j < self.dim: # secondChar is a symbol
                                                for i1 in range(self.dim):
                                                    bigFxy[i1,j] += na * oneOverNSlots
                                            elif j == bigDim - 1:
                                                raise Glitch, "This shouldn't happen."
                                            else: # secondChar is an equate
                                                for i1 in range(self.dim):
                                                    for j1 in range(self.dim):
                                                        if equatesArray[j - self.dim, j1]:
                                                            bigFxy[i1,j1] += na * oneOverNSlots

                                        else:  # firstChar an equateSymbol
                                            # secondChar could be anything
                                            if j < self.dim: # secondChar is a symbol
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        bigFxy[i1,j] += na * oneOverNSlots
                                            elif j == bigDim - 1:  # secondChar is N_LIKE
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        for j1 in range(self.dim):
                                                            bigFxy[i1,j1] += na * oneOverNSlots

                                            else:  # secondChar is an equate
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        for j1 in range(self.dim):
                                                            if equatesArray[j - self.dim, j1]:
                                                                bigFxy[i1,j1] += na * oneOverNSlots

                                    else: # fsum is not zero
                                        if i < self.dim: # firstChar is a symbol
                                            # secondChar must be N_LIKE or an equateSymbol
                                            if j == bigDim - 1:  # N_LIKE
                                                for j1 in range(self.dim):
                                                    bigFxy[i,j1] += na * normUnambig[i, j1] / fsum
                                            elif j >= self.dim:  # an equate
                                                for j1 in range(self.dim):
                                                    if equatesArray[j - self.dim, j1]:
                                                        bigFxy[i,j1] += na * normUnambig[i, j1] / fsum
                                            else:
                                                raise Glitch, "This shouldn't happen"
                                        elif i == bigDim - 1: # firstChar is N_LIKE
                                            # secondChar must be either a symbol or an equate
                                            if j < self.dim: # secondChar is a symbol
                                                for i1 in range(self.dim):
                                                    bigFxy[i1,j] += na * normUnambig[i1, j] / fsum
                                            elif j == bigDim - 1:
                                                raise Glitch, "This shouldn't happen."
                                            else: # secondChar is an equate
                                                for i1 in range(self.dim):
                                                    for j1 in range(self.dim):
                                                        if equatesArray[j - self.dim, j1]:
                                                            bigFxy[i1,j1] += na * normUnambig[i1, j1] / fsum

                                        else:  # firstChar an equateSymbol
                                            # secondChar could be anything
                                            if j < self.dim: # secondChar is a symbol
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        bigFxy[i1,j] += na * normUnambig[i1, j] / fsum
                                            elif j == bigDim - 1:  # secondChar is N_LIKE
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        for j1 in range(self.dim):
                                                            bigFxy[i1,j1] += na * normUnambig[i1, j1] / fsum

                                            else:  # secondChar is an equate
                                                for i1 in range(self.dim):
                                                    if equatesArray[i - self.dim, i1]:
                                                        for j1 in range(self.dim):
                                                            if equatesArray[j - self.dim, j1]:
                                                                bigFxy[i1,j1] += na * normUnambig[i1, j1] / fsum
                                    if 0:
                                        print "bigFxy= (after partial ambig resolution)"
                                        print bigFxy


                    # End of the long section on resolving ambiguities.
                # End of the long "else" clause to "if fastFillFxy:"

                if 0:
                    print "bigFxy=  (after ambig resolution)"
                    print bigFxy

                # pInvar stuff
                if doPInvarOfConstants==False and pInvar != None: # paup-like
                    nSitesToRemove = pInvar * bigFxy.sum() # sums over both axes
                    #print "pInvar=%s, nSitesToRemove=%s" % (pInvar, nSitesToRemove)
                    for i in range(self.dim):
                        bigFxy[i][i] -= constComps[i] * nSitesToRemove
                elif doPInvarOfConstants==True and pInvarOfConstants != None: # LDDist-like
                    for i in range(self.dim):
                        bigFxy[i][i] -= constCounts[i] * pInvarOfConstants

                if 0:
                    print "bigFxy=  (after pInvarOfConstants removal)"
                    print bigFxy

                if missingCharacterStrategy == 'fudge':
                    # Replace zeros on the diagonal of Fxy with either 0.5
                    # or half of the smallest positive Fii, whichever is
                    # smaller.
                    minPositiveOnDiag = 1.0
                    for i in range(self.dim):
                        if bigFxy[i,i] > 0.0 and bigFxy[i,i] < minPositiveOnDiag:
                            minPositiveOnDiag = bigFxy[i,i]
                    minPositiveOnDiag /= 2.0
                    for i in range(self.dim):
                        if bigFxy[i,i] < minPositiveOnDiag:
                            bigFxy[i,i] = minPositiveOnDiag
                            fudgeCount += 1
                    theFxy = bigFxy

                elif missingCharacterStrategy == 'reduce':
                    if minCompCount and hasIgnores:
                        # make a smaller matrix
                        smallerFxy = numpy.zeros((totalNoIgnores, totalNoIgnores), numpy.float)
                        i2 = 0
                        for i in range(self.dim):
                            if ignores[i]:
                                pass
                            else:
                                j2 = 0
                                for j in range(self.dim):
                                    if ignores[j]:
                                        pass
                                    else:
                                        #print "(%i, %i) -> (%i, %i)" % (i,j,i2,j2)
                                        smallerFxy[i2,j2] = bigFxy[i,j]
                                        j2 += 1
                                i2 += 1
                        #print "smallerFxy ="
                        #print smallerFxy
                        theFxy = smallerFxy
                    else:
                        theFxy = bigFxy
                else:
                    theFxy = bigFxy

                sumTheFxy = sum(sum(theFxy))

                # Normalize Fxy to 1.0
                theFxy /= float(sumTheFxy)
                if 0:
                    print "nSites=%f,  final normalized bigFxy (ie Fxy) = " % sumTheFxy
                    print theFxy

                # Calculate the logDet, from the theFxy and from the bigPi
                theDet = numpy.linalg.det(theFxy)
                #print "theDet from theFxy is %f" % theDet


                #if theDet <= 0.0:
                #    print "Got non-positive logDet."

                if nonPositiveDetStrategy == 'invert':
                    theDet = numpy.fabs(theDet)
                    if theDet < 1e-50: # eg zero
                        theDet = 1e-50
                        invertCount += 1

                if theDet > 0.0:
                    bigPiX = theFxy.sum(axis=0)   # sum of columns
                    bigPiY = theFxy.sum(axis=1)   # sum of rows
                    if 0:
                        print "theFxy = "
                        print theFxy, type(theFxy)
                        print "bigPiX = "
                        print bigPiX, type(bigPiX)
                        print "bigPiY ="
                        print bigPiY, type(bigPiY)
                    det_bigPiX = numpy.multiply.reduce(bigPiX)
                    det_bigPiY = numpy.multiply.reduce(bigPiY)
                    #print det_bigPiX, det_bigPiY
                    if det_bigPiX <= 0.0 or det_bigPiY <= 0.0:
                        if missingCharacterStrategy == 'refuse':
                            d.matrix[sNum1][sNum2] = -1.0
                            d.matrix[sNum2][sNum1] = -1.0
                            nUndefinedLogDets += 1
                        else:
                            if 0:
                                gm.append("sumTheFxy = %f" % sumTheFxy)
                                if missingCharacterStrategy == 'reduce':
                                    gm.append("reduce is on. hasIgnores=%s" % hasIgnores)
                                    gm.append("ignores = %s" % ignores)
                                gm.append("sum(bigPiX)=%s, sum(bigPiY)=%s" % (numpy.sum(bigPiX), numpy.sum(bigPiY)))
                                gm.append("symbols = %s" % self.symbols)
                                gm.append("sNum1=%i, sNum2=%i" % (sNum1, sNum2))
                                gm.append("bigPiX = %s" % bigPiX)
                                gm.append("bigPiY=%s" % bigPiY)
                                gm.append("det_bigPiX = %s, det_bigPiY=%s" % (det_bigPiX, det_bigPiY))
                                gm.append("Got bad Pi det due to missing char(s).")
                                gm.append("This should not happen-- programming error.")
                                raise Glitch, gm
                            if 1:
                                i2 = 0
                                for i1 in range(self.dim):
                                    if ignores[i1]:
                                        pass
                                    else:
                                        if bigPiX[i2] < 1.e-10:
                                            ignores[i1] = 1
                                            hasNewIgnores = True
                                        if bigPiY[i2] < 1.e-10:
                                            ignores[i1] = 1
                                            hasNewIgnores = True
                                        i2 += 1
                                assert hasNewIgnores, "If not hasNewIgnores, why are we here?"
                                #print "Added a new ignore, starting over ..."
                                totalNoIgnores = 0
                                for i in range(self.dim):
                                    if not ignores[i]:
                                        totalNoIgnores += 1
                                #print "totalNoIgnores = %i" % totalNoIgnores
                                if totalNoIgnores < 2:
                                    if 0:
                                        gm.append("The arg 'missingCharacterStrategy' is set to 'reduce'")
                                        gm.append("The arg 'minCompCount' is turned on, and set to %i." % minCompCount)
                                        gm.append("There is not enough variation in these sequences to make a valid distance.")
                                        gm.append("There are too many sites that will be ignored because of low frequency characters.")
                                        raise Glitch, gm
                                    else:
                                        return None
                                break



                                    
                    else:  # det_bigPiX and det_bigPiY are over zero
                        # If we have been using a reduced Fxy, then self.dim is no longer appropriate.
                        reducedDim = len(bigPiX)

                        if 0:
                            # This section works, but is not very clear.  Re-written below
                            theLogDet = numpy.log(theDet) - 0.5 * numpy.log(det_bigPiX) - 0.5 * numpy.log(det_bigPiY)

                            if correction == 'L94':
                                theCorrection = reducedDim
                            elif correction == 'TK02': # from LDDist
                                squareSum = 0.0
                                for i in range(reducedDim):
                                    squareSum += (bigPiX[i] + bigPiY[i]) * (bigPiX[i] + bigPiY[i])
                                theCorrection = (reducedDim - 1) / (1.0 - squareSum/4.0)
                            if 0:
                                print "theDet = %g" % theDet
                                print "theLogDet = %f" % theLogDet
                                print "theCorrection = %s" % theCorrection

                            theLogDet /= theCorrection
                            theLogDet = -theLogDet

                        if 1:
                            if correction == 'L94':  # equation 3, pg 606, in L94
                                theLogDet = -numpy.log(theDet) + (numpy.log(det_bigPiX) + numpy.log(det_bigPiY))/2.0
                                theLogDet /= reducedDim
                            elif correction == 'TK02': # equation 11, page 1729, in TK02
                                theLogDet = numpy.log(theDet) - (0.5 *(numpy.log(det_bigPiX) + numpy.log(det_bigPiY)))
                                squareSum = 0.0
                                for i in range(reducedDim):
                                    thePi = (bigPiX[i] + bigPiY[i])/2.0
                                    squareSum += thePi * thePi
                                theLogDet = -((1.0 - squareSum)/(reducedDim - 1)) * theLogDet
                            elif correction == 'TK02_eqn10':
                                theLogDet = -(1./reducedDim)*numpy.log(theDet) - numpy.log(reducedDim)
                                
                        # dset allsitesmean=yes
                        if doPInvarOfConstants==True and pInvarOfConstants != None:
                            theLogDet *= 1.0 - (pInvarOfConstants * nConstants)/self.nChar
                        if doPInvarOfConstants==False and pInvar != None:
                            theLogDet *= 1.0 - pInvar

                        if theLogDet < 0.0:
                            gm.append("Got negative logDet (%f).  This should not happen." % theLogDet)
                            raise Glitch, gm

                        #return theLogDet
                        d.matrix[sNum1][sNum2] = theLogDet
                        d.matrix[sNum2][sNum1] = theLogDet
                else:
                    if nonPositiveDetStrategy == 'refuse':
                        d.matrix[sNum1][sNum2] = -1.0
                        d.matrix[sNum2][sNum1] = -1.0
                        nUndefinedLogDets += 1
                    else:
                        gm.append("This should never happen.  Programming error.")
                        raise Glitch, gm
                    
    # End of the main pairwise loop

    
    if (missingCharacterStrategy == 'refuse' or nonPositiveDetStrategy == 'refuse') and nUndefinedLogDets:
        if nUndefinedLogDets == ((d.dim * d.dim) - d.dim) / 2:
            if 0:
                gm.append("All distances were undefined.")
                raise Glitch, gm
            else:
                return None
        #print "xyz There were %i undefined distances." % nUndefinedLogDets
        biggest = 0.0
        for sNum1 in range(self.nTax - 1):
            for sNum2 in range(sNum1 + 1, self.nTax):
                if d.matrix[sNum1][sNum2] > biggest:
                    biggest = d.matrix[sNum1][sNum2]
        biggest *= 2.0
        for sNum1 in range(self.nTax - 1):
            for sNum2 in range(sNum1 + 1, self.nTax):
                if numpy.fabs(d.matrix[sNum1][sNum2] -  -1.0) < 1e-10:
                    d.matrix[sNum1][sNum2] = biggest
                    d.matrix[sNum2][sNum1] = biggest

    dMessage = ["    "]
    dMessage.append("Log det distances from p4.")
    dMessage.append("Correction from %s" % correction)
    if doPInvarOfConstants:
        dMessage.append("doPInvarOfConstants is set, and pInvarOfConstants is %s" % (pInvarOfConstants))
    else:
        dMessage.append("doPInvarOfConstants is off, and pInvar is %s" % (pInvar))
    dMessage.append("The missingCharacterStrategy is set to '%s'." % missingCharacterStrategy)
    if missingCharacterStrategy == 'fudge':
        dMessage.append("    Did %i fudges." % fudgeCount)
    if missingCharacterStrategy == 'reduce':
        dMessage.append("    minCompCount = %i" % minCompCount)
        if hasIgnores:
            theIgnored = [self.symbols[i] for i in range(self.dim) if ignores[i]]
            dMessage.append("    These symbols were ignored: %s" % theIgnored)
        else:
            dMessage.append("    No chars were ignored.")
    dMessage.append("The nonPositiveDetStrategy is set to '%s'." % nonPositiveDetStrategy)
    if nonPositiveDetStrategy == 'invert':
        dMessage.append("    %i dets were inverted." % invertCount)
    if nonPositiveDetStrategy == 'refuse' or missingCharacterStrategy == 'refuse':
        if nUndefinedLogDets:
            dMessage.append("There were %i undefined log dets that were set to twice the largest defined distance (%f)"
                            % (nUndefinedLogDets, biggest))
        else:
            dMessage.append("There were no undefined log dets.")
    d.message = string.join(dMessage, '\n    ')

    return d

def  _logDetSetReduceIgnores(self, doPInvarOfConstants, pInvar, pInvarOfConstants, minCompCount, seq, constComps, constCounts):
    """Called by logDet when the missingCharactersStrategy is 'reduce'

    Make the ignores vector by going thru each sequence, not thru each
    pair.  Since it is a compromize, and does not resolve ambigs, it
    can miss ignores, but it is a lot faster than the slow version of
    this method.  """

    
    ignores = numpy.zeros(self.dim, numpy.int32)
    counts = numpy.zeros(self.dim, numpy.int32)
    for sNum in range(self.nTax):
        #print "sNum = %i" % sNum
        theSeq = seq[sNum]
        for i in range(self.dim):
            counts[i] = 0
        for cNum in range(self.nChar):
            if theSeq[cNum] >= 0:
                counts[theSeq[cNum]] += 1
        #print "counts = %s" % counts
        
        # pInvar stuff
        if doPInvarOfConstants==False and pInvar != None: # paup-like
            nSitesToRemove = pInvar * numpy.sum(counts)
            #print "pInvar=%s, nSitesToRemove=%s" % (pInvar, nSitesToRemove)
            for i in range(self.dim):
                counts[i] -= constComps[i] * nSitesToRemove
        elif doPInvarOfConstants==True and pInvarOfConstants != None: # LDDist-like
            for i in range(self.dim):
                counts[i] -= constCounts[i] * pInvarOfConstants

        #print "counts = %s" % counts

        for i in range(self.dim):
            if counts[i] < minCompCount:
                ignores[i] = 1
    return ignores



