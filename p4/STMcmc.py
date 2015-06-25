# This is STMcmc, for super tree mcmc.
# Started 18 March 2011, first commit 22 March 2011.

import pf,func
from Var import var
import math,random,string,sys,time,copy,os,cPickle,types,glob
import numpy as np
from Glitch import Glitch
from TreePartitions import TreePartitions
from Constraints import Constraints
from Tree import Tree
import datetime
import itertools

try:
    import bitarray
except ImportError:
    pass


def choose(n, k):
    """
    A fast way to calculate binomial coefficients 
    by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

# def nSplits(n):
#     mySum = 0
#     for k in range(2, n-1):
#         mySum += choose(n-1, k)
#     return mySum

def bForN(n):
    # This is the log version of this function.  The max diff (in
    # log(result)) between this and the non-log function seems to be
    # about 2.5e-10 for n up to 10000.

    prodLog = 0.0
    if n > 3:
        for k in range(4, n + 1):
            prodLog += math.log((2 * k) - 5)
    return prodLog




def BS2009_Eqn30_ZTApprox(n, beta, cT):
    # This log version of this function differs from from the non-log
    # version (in log(result)) by at most 6.82e-13 for n up to 150,
    # over a wide range of beta (0.001 -- 1000) and cT (2 -- n/2)

    myLambda = cT/(2.0*n)
    tester = 0.5 * math.log((n - 3.)/myLambda)

    epsilon = math.exp(-2. * beta)
    bigANEpsilon = 1 + (((2. * n) - 3.) * epsilon) + (2. * ((n * n) - (4. * n) - 6.) * epsilon * epsilon)
    termA = math.log(bigANEpsilon + 6 * cT * epsilon * epsilon)

    if beta < tester:
        termB = -(2. * beta) * (n - 3.) + (myLambda * (math.exp(2. * beta) - 1.))
        termB += bForN(n)
        if termA > termB:
            return termA
        else:
            return termB
    else:
        return termA


def popcountA(k, nBits):
    count = 0
    for i in range(nBits):
        tester = 1L << i
        if tester > k:
            return count
        if tester & k:
            count += 1
    return count

    


def bitReduce(bk, txBits, lLen, sLen, allOnes):
    #print "bitReduce: bk %i, txBits %i, lLen %i, sLen %i, allOnes %i" % (bk, txBits, lLen, sLen, allOnes)
    newBk = 0L
    counter = 0
    pops = 0
    for pos in range(lLen):
        tester = 1L << pos
        #print "pos %2i, tester: %3i" % (pos, tester)
        if tester & txBits:
            #print "    tester & txBits -- True"
            if tester & bk:
                adder = 1L << counter
                #print "        adding:", adder
                newBk += adder
                pops += 1
            else:
                #print "        not adding"
                pass
            counter += 1
    if (1 & newBk):
        #print "flipping"
        newBk = allOnes ^ newBk
        pops = sLen - pops
    #print "returning newBk %i, pops %i" % (newBk, pops)
    return newBk, pops

if 0: # test bitReduce
    sk = 6   # always at least 2 bits, even
    txBits = 30
    lLen = 5
    sLen = 4
    allOnes = 15
    print "     sk: %3i  %s" % (sk, func.getSplitStringFromKey(sk, lLen))
    print "taxBits: %3i  %s" % (txBits, func.getSplitStringFromKey(txBits, lLen))
                          
    rsk, popcount = bitReduce(sk, txBits, lLen, sLen, allOnes)
    print "    rsk: %3i  %s" % (rsk, func.getSplitStringFromKey(rsk, sLen))
    print "   popcount %i" % popcount
    #     sk:   6  .**..
    #     taxBits:  30  .****
    #     rsk:  12  ..**
    #     popcount 2


def maskedSymmetricDifference(skk, skSet, taxBits, longLen, shortLen, allOnes):
    if 0:
        print "-" * 50
        print "skk (skk_ppy1 from the current supertree)"
        for sk in skk:
            print func.getSplitStringFromKey(sk, longLen)
        print "skSet (from input tree)"
        for sk in skSet:
            print func.getSplitStringFromKey(sk, shortLen)
        print "taxBits:", taxBits, func.getSplitStringFromKey(taxBits, longLen)

    newSkk = []
    for sk in skk:
        reducedSk, popcount = bitReduce(sk, taxBits, longLen, shortLen, allOnes)
        if 0:
            print "taxBits: %s  " % func.getSplitStringFromKey(taxBits, longLen),
            print "%4i %s  " % (sk, func.getSplitStringFromKey(sk, longLen)),
            print "%4i %s  %i" % (reducedSk, func.getSplitStringFromKey(reducedSk, shortLen), popcount)
        if popcount <= 1 or popcount >= (shortLen - 1):
            pass
        else:
            newSkk.append(reducedSk)
    newSkkSet = set(newSkk)
    #print newSkkSet, skSet
    #print "reduced supertree splits =  newSkkSet = %s" % newSkkSet
    ret = len(newSkkSet.symmetric_difference(skSet))
    #print "symmetric difference %i" % ret
    nCherries = 0
    for sk in newSkkSet:
        popcount = popcountA(sk, shortLen)
        if popcount == 2:
            nCherries += 1
        if popcount == (shortLen - 2):  # not "elif", because they might both be True
            nCherries += 1
    #print "nCherries %i" % nCherries 
    return ret, nCherries

def slowQuartetDistance(st, inputTree):
    dst = st.dupe()
    toRemove = []
    for n in dst.iterLeavesNoRoot():
        if n.name not in inputTree.taxNames:
            toRemove.append(n)
    for n in toRemove:
        dst.removeNode(n)
    qd = dst.topologyDistance(inputTree, metric='scqdist')
    return qd
    

class STChain(object):

    def __init__(self, aSTMcmc):
        gm = ['STChain.__init__()']

        self.stMcmc = aSTMcmc
        self.tempNum = -1 # 'temp'erature, not 'temp'orary

        self.curTree = aSTMcmc.tree.dupe()
        self.propTree = aSTMcmc.tree.dupe()

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0

        self.frrf = None
        self.nInTreeSplits = 0

        if self.stMcmc.modelName.startswith('SR2008_rf'):
            self.curTree.beta = self.stMcmc.beta
            self.propTree.beta = self.stMcmc.beta

            if self.stMcmc.stRFCalc == 'purePython1':
                self.getTreeLogLike_ppy1()  
                
            elif self.stMcmc.stRFCalc == 'fastReducedRF':
                self.startFrrf()
                self.getTreeLogLike_fastReducedRF()

            elif self.stMcmc.stRFCalc == 'bitarray':
                self.setupBitarrayCalcs()
                self.getTreeLogLike_bitarray()
            
            self.curTree.logLike = self.propTree.logLike
        elif self.stMcmc.modelName.startswith('SPA'):
            self.curTree.spaQ= self.stMcmc.spaQ
            self.propTree.spaQ = self.stMcmc.spaQ

            for t in self.stMcmc.trees:
                self.nInTreeSplits += len(t.splSet)
            #print "Got nInTreeSplits %s" % self.nInTreeSplits
            self.setupBitarrayCalcs()
            self.getTreeLogLike_spa_bitarray()
            self.curTree.logLike = self.propTree.logLike
        elif self.stMcmc.modelName.startswith('QPA'):
            self.curTree.spaQ= self.stMcmc.spaQ
            self.propTree.spaQ = self.stMcmc.spaQ
            self.nPossibleQuartets = choose(self.stMcmc.tree.nTax, 4) * 3
            self.getTreeLogLike_qpa_slow()
            self.curTree.logLike = self.propTree.logLike

        else:
            gm.append('Unknown modelName %s' % self.stMcmc.modelName)
            raise Glitch, gm

        if 0:
            print "STChain init()"
            self.curTree.draw()
            print "logLike is %f" % self.curTree.logLike


    def getTreeLogLike_qpa_slow(self):
        gm = ["STChain.getTreeLogLike_qpa_slow()"]
        if self.propTree.spaQ > 1. or self.propTree.spaQ <= 0.0:
            gm.append("bad propTree.spaQ value %f" % self.propTree.spaQ)
            raise Glitch, gm
        
        for n in self.propTree.iterInternalsPostOrder():
            if n == self.propTree.root:
                break
            n.stSplitKey = n.leftChild.stSplitKey
            p = n.leftChild.sibling
            while p:
                n.stSplitKey |= p.stSplitKey    # "or", in-place
                p = p.sibling
        self.propTree.skk = [n.stSplitKey for n in self.propTree.iterInternalsNoRoot()]
        self.propTree.qSet = set()
        for sk in self.propTree.skk:
            ups = [txBit for txBit in self.propTree.taxBits if (sk & txBit)]
            downs = [txBit for txBit in self.propTree.taxBits if not (sk & txBit)]
            for down in itertools.combinations(downs, 2):
                if down[0] > down[1]:
                    down = (down[1], down[0])
                for up in itertools.combinations(ups, 2):
                    if up[0] > up[1]:
                        up = (up[1], up[0])
                    if down[0] < up[0]:
                        self.propTree.qSet.add(down+up)
                    else:
                        self.propTree.qSet.add(up+down)
        #print self.propTree.qSet
        self.propTree.nQuartets = len(self.propTree.qSet) 
        
        if self.propTree.nQuartets:
            q = self.propTree.spaQ / self.propTree.nQuartets
            R = 1. - self.propTree.spaQ
            r = R / (self.nPossibleQuartets - self.propTree.nQuartets)
            logq = math.log(q)
        else:
            R = 1.
            r = R / self.nPossibleQuartets
        logr = math.log(r)
        self.propTree.logLike = 0.0
        for it in self.stMcmc.trees:
            for qu in it.qSet:
                if qu in self.propTree.qSet:
                    self.propTree.logLike += logq
                else:
                    self.propTree.logLike += logr
        

    def getTreeLogLike_spa_bitarray(self):
        gm = ["STChain.getTreeLogLike_spa_bitarray"]
        if self.propTree.spaQ > 1. or self.propTree.spaQ <= 0.0:
            gm.append("bad propTree.spaQ value %f" % self.propTree.spaQ)
            raise Glitch, gm
        slowCheck = False
        if slowCheck:
            slowCheckLogLike = 0.0
            for it in self.stMcmc.trees:
                it.makeSplitKeys()
                it.skk =  [n.br.splitKey for n in it.iterInternalsNoRoot()]

        self.propTree.logLike = 0.0
        for it in self.stMcmc.trees:
            if 0:
                print "-" * 50
                it.draw()
                print "baTaxBits %s" % it.baTaxBits
                print "firstTax at %i" % it.firstTax

            if slowCheck:
                stDupe = self.propTree.dupe()
                toRemove = []
                for n in stDupe.iterLeavesNoRoot():
                    if n.name not in it.taxNames:
                        toRemove.append(n)
                for n in toRemove:
                    stDupe.removeNode(n)
                stDupe.taxNames = it.taxNames
                stDupe.makeSplitKeys(makeNodeForSplitKeyDict=True)

            # No need to consider (masked) splits with less than two
            # 1s or more than it.nTax - 2 1s.  
            upperGood = it.nTax - 2
            relevantStSplits = []
            for n in self.propTree.iterInternalsNoRoot():
                # Choose which spl (spl or spl2) based on it.firstTax)
                if n.ss.spl[it.firstTax]:
                    n.ss.theSpl = n.ss.spl
                else:
                    n.ss.theSpl = n.ss.spl2
                n.ss.maskedSplitWithTheFirstTaxOne = n.ss.theSpl & it.baTaxBits
                n.ss.onesCount = n.ss.maskedSplitWithTheFirstTaxOne.count()
                if 0:
                    print "bigT node %i" % n.nodeNum
                    print "  theSpl is %s" % n.ss.theSpl
                    print "  maskedSplitWithTheFirstTaxOne %s" %  n.ss.maskedSplitWithTheFirstTaxOne
                    print "  onesCount %i" % n.ss.onesCount
                    if n.ss.onesCount >= 2 and n.ss.onesCount <= upperGood:
                        print "    -> relevant"
                    else:
                        print "    -> not relevant"
                if n.ss.onesCount >= 2 and n.ss.onesCount <= upperGood:
                    relevantStSplits.append(n.ss)

            nonRedundantStSplits = []
            for ss in relevantStSplits:
                alreadyIn = False
                for ssB in nonRedundantStSplits:
                    if ss.maskedSplitWithTheFirstTaxOne == ssB.maskedSplitWithTheFirstTaxOne:
                        alreadyIn = True
                        break
                if alreadyIn == False:
                    nonRedundantStSplits.append(ss)
            if 0:
                for ss in relevantStSplits:
                    ss.dump()
                print "There are %i relevant splits in the st for this it." % len(relevantStSplits)
                for ss in nonRedundantStSplits:
                    ss.dump()
                print "There are %i non-redundant splits in the st for this it." % len(nonRedundantStSplits)


            S_st = len(nonRedundantStSplits)  # S_st is the number of splits in the reduced supertree
            if slowCheck:
                #stDupe.draw()
                #print "the drawing above is stDupe"
                slowCheckS_st = len([n for n in stDupe.iterInternalsNoRoot()])
                assert S_st == slowCheckS_st
                
            S = 2**(it.nTax - 1) - (it.nTax + 1)     # S is the number of possible splits in an it-sized tree
            #print "S=%i, S_st=%i" % (S, S_st)
            if S_st:
                q = self.propTree.spaQ / S_st
                R = 1. - self.propTree.spaQ
                r = R/(S - S_st)
                #print "q=%f" % q
                logq = math.log(q)
            else:
                R = 1.
                r = R/S
            #print "r=%f" % r
            logr = math.log(r)

            # for ss in nonRedundantStSplits:
            #     ss.bytes = ss.maskedSplitWithTheFirstTaxOne.tobytes()
            #     ret = ss.bytes in it.splSet
            #     if ret:
            #         print "  iT has reduced split %s" % ss.bytes
            #         self.propTree.logLike += logq
            #     else:
            #         print "  iT does not have reduced split %s" % ss.bytes
            #         self.propTree.logLike += logr
            mySSForBytesDict = {}
            for ss in nonRedundantStSplits:
                ss.bytes = ss.maskedSplitWithTheFirstTaxOne.tobytes()
                mySSForBytesDict[ss.bytes] = ss
            for spl in it.splSet:
                ret = mySSForBytesDict.get(spl)
                if ret:
                    #print "  st has reduced split %s" % spl
                    self.propTree.logLike += logq
                else:
                    #print "  st does not have reduced split %s" % spl
                    self.propTree.logLike += logr
                
            
            if slowCheck:
                for sk in it.skk:
                    ret = stDupe.nodeForSplitKeyDict.get(sk)
                    if ret:
                        slowCheckLogLike += logq
                    else:
                        slowCheckLogLike += logr

                myDiff = self.propTree.logLike - slowCheckLogLike
                if math.fabs(myDiff) > 1.e-12:
                    gm.append("Bad like calc. slowCheck %f, bitarray %f, diff %g" % (
                        slowCheckLogLike, self.propTree.logLike, myDiff))
                    raise Glitch, gm
            

        


    def setupBitarrayCalcs(self):
        # Prepare self.propTree (ie bigT).  First make n.stSplitKeys.  These are temporary.
        for n in self.propTree.iterPostOrder():
            if n == self.propTree.root:
                break
            if n.isLeaf:
                spot = self.stMcmc.taxNames.index(n.name)
                self.stMcmc.tBits[spot] = True
                n.stSplitKey = bitarray.bitarray(self.stMcmc.tBits)
                self.stMcmc.tBits[spot] = False
            else:
                n.stSplitKey = n.leftChild.stSplitKey.copy()
                p = n.leftChild.sibling
                while p:
                    n.stSplitKey |= p.stSplitKey    # "or", in-place
                    p = p.sibling

        # Next transfer the internal node split keys to BigTSplitStuff objects
        for n in self.propTree.iterInternalsNoRoot():
            n.ss = BigTSplitStuff()
            n.ss.spl = n.stSplitKey
            n.ss.spl2 = n.ss.spl.copy()
            n.ss.spl2.invert()
        # This next one will be empty, not used immediately, but will
        # be used after supertree rearrangements.
        self.propTree.root.ss = BigTSplitStuff()

    def refreshBitarrayPropTree(self):
        # Refresh self.propTree (ie bigT) after a topology change.  
        for n in self.propTree.iterPostOrder():
            if n == self.propTree.root:
                break
            if n.isLeaf:
                pass
            else:
                n.stSplitKey = n.leftChild.stSplitKey.copy()
                p = n.leftChild.sibling
                while p:
                    n.stSplitKey |= p.stSplitKey    # "or", in-place
                    p = p.sibling

        # Next transfer the internal node split keys to BigTSplitStuff objects
        for n in self.propTree.iterInternalsNoRoot():
            n.ss.spl = n.stSplitKey
            n.ss.spl2 = n.ss.spl.copy()
            n.ss.spl2.invert()


    def startFrrf(self):
        # if using self.stMcmc.stRFCalc= 'fastReducedRF'
        self.frrf = self.stMcmc.Frrf(len(self.stMcmc.taxNames))
        self.bigTr = self.frrf.setBigT(len(self.propTree.nodes), self.propTree.nTax, self.propTree.postOrder)

        for n in self.propTree.nodes:
            if n.parent:
                self.bigTr.setParent(n.nodeNum, n.parent.nodeNum)
            if n.leftChild:
                self.bigTr.setLeftChild(n.nodeNum, n.leftChild.nodeNum)
            else:
                self.bigTr.setNodeTaxNum(n.nodeNum, self.stMcmc.taxNames.index(n.name))
            if n.sibling:
                self.bigTr.setSibling(n.nodeNum, n.sibling.nodeNum)

        if 1:
            for t in self.stMcmc.trees:
                tr = self.frrf.appendInTree(len(t.nodes), t.nTax, t.postOrder)
                for n in t.nodes:
                    if n.parent:
                        tr.setParent(n.nodeNum, n.parent.nodeNum)
                    if n.leftChild:
                        tr.setLeftChild(n.nodeNum, n.leftChild.nodeNum)
                    else:
                        tr.setNodeTaxNum(n.nodeNum, self.stMcmc.taxNames.index(n.name))
                    if n.sibling:
                        tr.setSibling(n.nodeNum, n.sibling.nodeNum)
        self.frrf.setInTreeTaxBits()
        self.frrf.setInTreeInternalBits()
        self.frrf.maybeFlipInTreeBits()
        self.frrf.setBigTInternalBits()
        #self.frrf.dump()
        

    def getTreeLogLike_ppy1(self):
        gm = ['STChain.getTreeLogLike_pp1']
        self.propTree.makeSplitKeys()
        self.propTree.skk = [n.br.splitKey for n in self.propTree.iterInternalsNoRoot()]
        self.propTree.logLike = 0.0
        for t in self.stMcmc.trees:

            # Get the distance
            thisDist = None
            if self.stMcmc.modelName.startswith('SR2008_rf'):
                thisDist, nCherries = maskedSymmetricDifference(self.propTree.skk, t.skSet, 
                                              t.taxBits, self.stMcmc.nTax, t.nTax, t.allOnes)
            else:
                raise Glitch, "STChain.getTreeLogLike_ppy1() unknown model '%s'" % self.stMcmc.modelName

            # Now multiply by beta, and do approximate Z_T
            assert thisDist != None
            beta_distance = self.propTree.beta * thisDist
            if self.stMcmc.modelName == 'SR2008_rf_ia':
                self.propTree.logLike -= beta_distance
            elif self.stMcmc.modelName.startswith('SR2008_rf_aZ'):
                log_approxZT = BS2009_Eqn30_ZTApprox(t.nTax, self.propTree.beta, nCherries)
                if 0:
                    # Testing, testing ...
                    assert self.propTree.beta == 0.1
                    assert t.nTax == 6
                    if nCherries == 2:
                        log_approxZT = 4.13695897651  # exact
                    elif nCherries == 3:
                        log_approxZT = 4.14853562562
                self.propTree.logLike -= log_approxZT
                self.propTree.logLike -= beta_distance
            else:
                gm.append("Unknown modelName %s" % self.stMcmc.modelName)
                raise Glitch, gm


    def getTreeLogLike_fastReducedRF(self):
        slowCheck = False
        if slowCheck:
            self.getTreeLogLike_ppy1()
            savedLogLike = self.propTree.logLike

        self.frrf.wipeBigTPointers()
        for n in self.propTree.nodes:
            if n.parent:
                self.bigTr.setParent(n.nodeNum, n.parent.nodeNum)
            if n.leftChild:
                self.bigTr.setLeftChild(n.nodeNum, n.leftChild.nodeNum)
            #else:
            #    bigTr.setNodeTaxNum(n.nodeNum, tNames.index(n.name))
            if n.sibling:
                self.bigTr.setSibling(n.nodeNum, n.sibling.nodeNum)
        self.frrf.setBigTInternalBits()
        if self.stMcmc.modelName == 'SR2008_rf_ia':
            sd = self.frrf.getSymmDiff()
            self.propTree.logLike = -sd * self.propTree.beta
        elif self.stMcmc.modelName.startswith('SR2008_rf_aZ'):
            self.propTree.logLike = self.frrf.getLogLike(self.propTree.beta)
        if slowCheck:
            if self.propTree.logLike != savedLogLike:
                gm = ['STChain.getTreeLogLike_fastReducedRF()']
                gm.append("Slow likelihood %f" % savedLogLike)
                gm.append("Fast likelihood %f" % self.propTree.logLike)
                raise Glitch, gm

                      
    def getTreeLogLike_bitarray(self):
        self.propTree.logLike = 0.0
        slowCheck = False
        if slowCheck:
            self.propTree.makeSplitKeys()
            self.propTree.skk = [n.br.splitKey for n in self.propTree.iterInternalsNoRoot()]
        for t in self.stMcmc.trees:
            if 0:
                print "-" * 50
                t.draw()
                print "baTaxBits %s" % t.baTaxBits
                print "firstTax at %i" % t.firstTax
            usables = []     # splitStuff objects with onesCount >= 2 and <= t.nTax = 2
            # No need to consider (masked) splits with less than two
            # 1s or more than nTax - 2 1s.  The nTax depends on the
            # input tree.
            upperGood = t.nTax - 2
            for n in self.propTree.iterInternalsNoRoot():
                # Choose which spl (spl or spl2) based on t.firstTax)
                if n.ss.spl[t.firstTax]:
                    n.ss.theSpl = n.ss.spl
                else:
                    n.ss.theSpl = n.ss.spl2
                n.ss.maskedSplitWithTheFirstTaxOne = n.ss.theSpl & t.baTaxBits
                n.ss.onesCount = n.ss.maskedSplitWithTheFirstTaxOne.count()
                if 0:
                    print "bigT node %i" % n.nodeNum
                    print "  theSpl is %s" % n.ss.theSpl
                    print "  maskedSplitWithTheFirstTaxOne %s" %  n.ss.maskedSplitWithTheFirstTaxOne
                    print "  onesCount %i" % n.ss.onesCount
                    if n.ss.onesCount >= 2 and n.ss.onesCount <= upperGood:
                        print "    -> used"
                    else:
                        print "    -> not used"
                if n.ss.onesCount >= 2 and n.ss.onesCount <= upperGood:
                    usables.append(n.ss)
            usablesDict = {}
            for usable in usables:
                usable.bytes = usable.maskedSplitWithTheFirstTaxOne.tobytes()
                usablesDict[usable.bytes] = usable
            splSet = set()   # bytes, for RF calculation
            for usable in usables:
                # splSet.add(n.ss.maskedSplitWithTheFirstTaxOne.tobytes())
                splSet.add(usable.bytes)
            thisBaRF =  len(splSet.symmetric_difference(t.splSet))
            if slowCheck:  # with purePython1
                thisPPyRF, thisPPyNCherries = maskedSymmetricDifference(self.propTree.skk, t.skSet, 
                                                      t.taxBits, self.stMcmc.nTax, t.nTax, t.allOnes)
                if thisBaRF != thisPPyRF:
                    raise Glitch, "bitarray and purePython1 RF calcs differ."
            beta_distance = self.propTree.beta * thisBaRF
            if self.stMcmc.modelName == 'SR2008_rf_ia':
                self.propTree.logLike -= beta_distance
            elif self.stMcmc.modelName.startswith('SR2008_rf_aZ'):
                nCherries = 0
                for ba in splSet:
                    theSS = usablesDict[ba]
                    #theSS.dump()
                    if theSS.onesCount == 2:
                        nCherries += 1
                    if theSS.onesCount == upperGood:
                        nCherries += 1
                if slowCheck:
                    if nCherries != thisPPyNCherries:
                        raise Glitch, "bitarray and purePython1 nCherries calcs differ."
                log_approxZT = BS2009_Eqn30_ZTApprox(t.nTax, self.propTree.beta, nCherries)
                self.propTree.logLike -= log_approxZT
                self.propTree.logLike -= beta_distance
            else:
                gm.append("Unknown model %s" % self.stMcmc.modelName)
                raise Glitch, gm

    def proposePolytomy(self, theProposal):
        theProposal.doAbort = False
        dbug = False
        if dbug:
            #print "proposePolytomy() starting with this tree ..."
            #self.propTree.draw(width=80, addToBrLen=0.2)
            print "j There are %i internal nodes." % self.propTree.nInternalNodes
            if self.propTree.nInternalNodes == 1:
                print "-> so its a star tree -> proposeDeleteEdge is not possible."
            elif self.propTree.nInternalNodes == self.propTree.nTax - 2:
                print "-> so its a fully-resolved tree, so proposeAddEdge is not possible."

        if self.propTree.nInternalNodes == 1: # a star tree
            self.proposeAddEdge(theProposal)
        elif self.propTree.nInternalNodes == self.propTree.nTax - 2:
            candidateNodes = self._getCandidateNodesForDeleteEdge()
            if candidateNodes:
                self.proposeDeleteEdge(theProposal, candidateNodes)
            else:
                #gm = ["proposePolytomy()"]
                #gm.append("The tree is fully resolved, so I can't proposeAddEdge()")
                #gm.append("But there are no suitable nodes to remove.")
                #raise Glitch, gm
                theProposal.doAbort = True
                self.curTree._nInternalNodes = self.propTree._nInternalNodes
                return
        else:
            r = random.random()
            #r = 0.4
            if r < 0.5:
                self.proposeAddEdge(theProposal)
            else:
                candidateNodes = self._getCandidateNodesForDeleteEdge()
                if candidateNodes:
                    self.proposeDeleteEdge(theProposal, candidateNodes)
                else:
                    self.proposeAddEdge(theProposal)
        #if self.mcmc.constraints:
        #    print "checkSplitKeys() at the end of polytomy"
        #    self.propTree.checkSplitKeys()
    

    def proposeAddEdge(self, theProposal):
        gm = ["STChain.proposeAddEdge()"]
        #print "proposeAddEdge() here"
        dbug = False
        pTree = self.propTree
        if 0:
            print "proposeAddEdge(), starting with this tree ..."
            pTree.draw()
            print "k There are %i internal nodes." % pTree.nInternalNodes
            print "root is node %i" % pTree.root.nodeNum
        allPolytomies = []
        for n in pTree.iterInternalsNoRoot():
            if n.getNChildren() > 2:
                allPolytomies.append(n)
        if pTree.root.getNChildren() > 3:
            allPolytomies.append(pTree.root)

        theChosenPolytomy = random.choice(allPolytomies)

        # We want to choose one of the possible ways to add a node.  See
        # Lewis et al page 246, left top.  "The number of distinct ways of
        # dividing k edges into two groups, making sure that at least 3
        # edges are attached to each node afterwards, is 2^{k-1} - k - 1".
        # For non-root polytomies (with 3 or more children), it is
        # straightforward, but for root polytomies (ie with 4 or more
        # children) it is different.  I think in the case of root
        # polytomies that they will be equivalent to non-root polytomies
        # if I arbitrarily consider one randomly chosen child node to
        # take the role that the parent takes in the non-root-polytomies.
        # So a 4-child root will be considered to have a parent-like node
        # and 3 children.
        if theChosenPolytomy != pTree.root:
            nChildren = theChosenPolytomy.getNChildren()
            k = nChildren + 1
            childrenNodeNums = pTree.getChildrenNums(theChosenPolytomy) 
        else:
            # Its the root.  So we say that a random child takes the role
            # of the "parent", for purposes of these calculations.
            nChildren = theChosenPolytomy.getNChildren() - 1  # n - 1 children
            k = nChildren + 1
            childrenNodeNums = pTree.getChildrenNums(theChosenPolytomy) # Yes, all children.

        nPossibleWays = math.pow(2, k-1) - k - 1
        if dbug:
            print "These nodes are polytomies: %s" % [n.nodeNum for n in allPolytomies]
            print "We randomly choose to do node %i" % theChosenPolytomy.nodeNum
            print "It has %i children, so k=%i, so there are %i possible ways to add a node." % (
                nChildren, k, nPossibleWays)

        # We want to choose one of the possible ways to add a node, but we
        # want to choose it randomly.  I'll describe it for the case with
        # nChildren=5, so k is 6.  We know already that there are
        # nPossibleWays=25 different ways to add a node.  The complication
        # is that we could make a new group of 2, 3, or 4 nInNewGroup, and it will be
        # different numbers of possible ways in each.  The numbers of each are given by
        # func.nChoosek(), so there are 10 ways to make a group of 2 from 5
        # children, 10 ways to make a group of 3 from 5 children, and 5
        # ways to make a group of 4 from 5 children.  So thats [10, 10,
        # 5], which sums to 25 (nPossibleWays).  So we can make a
        # cumulative sum list ie [10, 20, 25], and use it to choose one
        # group randomly.
        nChooseKs = []
        for i in range(2, nChildren):
            nChooseKs.append(func.nChooseK(nChildren, i))
        cumSum = [nChooseKs[0]]
        for i in range(len(nChooseKs))[1:]:
            cumSum.append(nChooseKs[i] + cumSum[i-1])
        ran = random.randrange(nPossibleWays)
        for i in range(len(cumSum)):
            if ran < cumSum[i]:
                break
        nInNewGroup = i + 2
        # Ok, so we have decided that of the nChildren of
        # theChosenPolytomy, we will make a new node with a group of
        # nInNewGroup of them.  For that, we can use random.sample().
        newChildrenNodeNums = random.sample(childrenNodeNums, nInNewGroup)


        if dbug:
            print "The nChooseKs are %s" % nChooseKs
            print "The cumSum is %s" % cumSum
            print "Since there are nPossibleWays=%i, we choose a random number from 0-%i" % (
                nPossibleWays, nPossibleWays-1)
            print "->We chose a random number: %i" % ran
            print "So we choose the group at index %i, which means nInNewGroup=%i" % (i, nInNewGroup)
            print "So we make a new node with newChildrenNodeNums %s" % newChildrenNodeNums
            #sys.exit()

        # Choose to add a node between theChosenPolytomy and the first in
        # the list of newChildrenNodeNums.  The node that we add will be
        # chosen from pTree.nodes for the first node where both the parent
        # and the leftChild are None.
        firstNode = pTree.nodes[newChildrenNodeNums[0]]
        for newNode in pTree.nodes:
            if not newNode.parent and not newNode.leftChild:
                break
        #print "Got newNode = %i" % newNode.nodeNum

        # Add the newNode between theChosenPolytomy and firstNode
        newNode.parent = theChosenPolytomy
        newNode.leftChild = firstNode
        firstNode.parent = newNode
        if theChosenPolytomy.leftChild == firstNode:
            theChosenPolytomy.leftChild = newNode
        else:
            oldCh = theChosenPolytomy.leftChild
            while oldCh.sibling != firstNode:
                oldCh = oldCh.sibling
            oldCh.sibling = newNode
        if firstNode.sibling:
            newNode.sibling = firstNode.sibling
            firstNode.sibling = None
        pTree.setPreAndPostOrder()
        pTree._nInternalNodes += 1

        if 0:
            #pTree.setPreAndPostOrder()
            pTree.draw()

        for nodeNum in newChildrenNodeNums[1:]:
            n = pTree.pruneSubTreeWithoutParent(nodeNum)
            pTree.reconnectSubTreeWithoutParent(n, newNode)


        # Calculate the rawSplitKey and splitKey.
        # if self.mcmc.constraints:
        #     children = [n for n in newNode.iterChildren()]
        #     x = children[0].br.rawSplitKey
        #     for n in children[1:]:
        #         y = n.br.rawSplitKey
        #         x = x | y  # '|' is bitwise "OR".
        #     newNode.br.rawSplitKey = x
        #     if 1 & newNode.br.rawSplitKey: # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
        #         if self.mcmc.constraints:
        #             newNode.br.splitKey = self.mcmc.constraints.allOnes ^ newNode.br.rawSplitKey # "^" is xor, a bit-flipper.
        #         else:
        #             allOnes = 2L**(self.propTree.nTax) - 1
        #             newNode.br.splitKey = allOnes ^ newNode.br.rawSplitKey
        #     else:
        #         newNode.br.splitKey = newNode.br.rawSplitKey

        # Its a newly-added node, possibly in a new context.  We need to
        # deal with model stuff if it isHet.  The model.isHet if any part
        # isHet.
        if dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()

        # Now the Hastings ratio.  First calculate gamma_B.  If the
        # current tree is a star tree (nInternalNodes == 1) and the
        # proposed tree is not fully resolved (ie is less than
        # len(self.propTree.nodes) - 2), then gamma_B is 0.5.
        if (self.curTree.nInternalNodes == 1) and (pTree.nInternalNodes < (len(pTree.nodes) - 2)):
            gamma_B = 0.5
        # If the proposed tree is fully resolved and the current tree is not the star tree
        elif (pTree.nInternalNodes == (len(pTree.nodes) - 2)) and (self.curTree.nInternalNodes > 1):
            gamma_B = 2.0
        else:
            gamma_B = 1.0

        # n_e is number of internal edges present before the Add-edge move.  That would be self.curTree.nInternalNodes - 1
        n_e = float(self.curTree.nInternalNodes - 1)
        # n_p is the number of polytomies present before the move, len(allPolytomies)
        n_p = float(len(allPolytomies))
        hastingsRatio = (gamma_B * n_p * float(nPossibleWays)) / (1.0 + n_e)

        if dbug:
            print "The new node is given a random branch length of %f" % newNode.br.len
            print "For the Hastings ratio ..."
            print "gamma_B is %.1f" % gamma_B
            print "n_e is %.0f" % n_e
            print "k is (still) %i, and (2^{k-1} - k - 1) = nPossibleWays is still %i" % (k, nPossibleWays)
            print "n_p = %.0f is the number of polytomies present before the move." % n_p
            print "So the hastings ratio is %f" % hastingsRatio

        self.logProposalRatio = math.log(hastingsRatio)


        if 0:
            priorRatio = self.mcmc.tunings.brLenPriorLambda * math.exp(- self.mcmc.tunings.brLenPriorLambda * newNode.br.len)
            if dbug:
                print "The self.mcmc.tunings.brLenPriorLambda is %f" % self.mcmc.tunings.brLenPriorLambda
                print "So the prior ratio is %f" % priorRatio

            self.logPriorRatio = math.log(priorRatio)

            # The Jacobian
            jacobian = 1.0 / (self.mcmc.tunings.brLenPriorLambda * math.exp(- self.mcmc.tunings.brLenPriorLambda * newNode.br.len))
            self.logJacobian = math.log(jacobian)
            print "logPriorRatio = %f, logJacobian = %f" % (self.logPriorRatio, self.logJacobian)

        # Here I pull a fast one, as explained in Lewis et al.  The
        # priorRatio and the Jacobian terms cancel out.  So the logs might
        # as well be zeros.
        self.logPriorRatio = 0.0
        #self.logJacobian = 0.0
        # That was easy, wasn't it?
        if self.stMcmc.tunings.doPolytomyResolutionClassPrior:
            # We are gaining a node.  So the prior ratio is T_{n,m + 1} /
            # (T_{n,m} * C) .  We have the logs, and the result is the
            # log.
            if 0:
                print "-" * 30
                print 'curTree.nInternalNodes', self.curTree.nInternalNodes
                print 'pTree.nInternalNodes', pTree.nInternalNodes
                print 'logBigT[curTree.nInternalNodes]', theProposal.logBigT[self.curTree.nInternalNodes]
                #print math.exp(theProposal.logBigT[self.curTree.nInternalNodes])
                print 'C ', self.stMcmc.tunings.polytomyPriorLogBigC
                print 'logBigT[pTree.nInternalNodes]', theProposal.logBigT[pTree.nInternalNodes]
                #print math.exp(theProposal.logBigT[pTree.nInternalNodes])
                print "-" * 30
            self.logPriorRatio = (theProposal.logBigT[self.curTree.nInternalNodes] -
                                  (self.stMcmc.tunings.polytomyPriorLogBigC +
                                  theProposal.logBigT[pTree.nInternalNodes]))

        else:
            if self.stMcmc.tunings.polytomyPriorLogBigC:
                self.logPriorRatio =  -self.stMcmc.tunings.polytomyPriorLogBigC
            else:
                self.logPriorRatio = 0.0
        #print "gaining a node, m %2i->%2i. logPriorRatio is %f" % (self.curTree.nInternalNodes,
        #                                                              pTree.nInternalNodes, self.logPriorRatio)


    def _getCandidateNodesForDeleteEdge(self):
        pTree = self.propTree
        nodesWithInternalEdges = [n for n in pTree.iterInternalsNoRoot()]

        # Remove any that might violate constraints.
        # if self.mcmc.constraints:
        #     nodesToRemove = []
        #     for n in nodesWithInternalEdges:
        #         if n.br.splitKey in self.mcmc.constraints.constraints:
        #             nodesToRemove.append(n)
        #     for n in nodesToRemove:
        #         nodesWithInternalEdges.remove(n)

        return nodesWithInternalEdges

    def proposeDeleteEdge(self, theProposal, candidateNodes):

        dbug = False
        pTree = self.propTree
        #print "doing proposeDeleteEdge()"
        if 0:
            print "proposeDeleteEdge(), starting with this tree ..."
            pTree.draw()
            print "m There are %i internal nodes (before deleting the edge)." % pTree.nInternalNodes

        if not candidateNodes:
            raise Glitch, "proposeDeleteEdge() could not find a good node to attempt to delete."

        theChosenNode = random.choice(candidateNodes)
        if dbug:
            print "There are %i candidateNodes." % len(candidateNodes)
            print "node nums %s" % [n.nodeNum for n in candidateNodes]
            print "Randomly choose node %s" % theChosenNode.nodeNum



        theNewParent = theChosenNode.parent
        theRightmostChild = theChosenNode.rightmostChild()
        theLeftSib = theChosenNode.leftSibling()
        if theLeftSib:
            theLeftSib.sibling = theChosenNode.leftChild
        else:
            theNewParent.leftChild = theChosenNode.leftChild
        for n in theChosenNode.iterChildren():
            n.parent = theNewParent
        theRightmostChild.sibling = theChosenNode.sibling
        theChosenNode.wipe()
        pTree.setPreAndPostOrder()
        pTree._nInternalNodes -= 1
        #print pTree.preOrder
        #if dbug:
        #    pTree.draw()


        # Hastings ratio.  First calculate the gamma_D.  If the current
        # tree is fully resolved and the proposed tree is not the star
        # tree, then gamma_D is 0.5
        if (self.curTree.nInternalNodes == len(pTree.nodes) - 2) and pTree.nInternalNodes != 1:
            gamma_D = 0.5
        # If the proposed tree is the star tree and the current tree is not fully resolved
        elif (self.curTree.nInternalNodes < len(pTree.nodes) - 2) and pTree.nInternalNodes == 1:
            gamma_D = 2.
        else:
            gamma_D = 1.

        # n_e is the number of internal edges in existence before the move, which would be nInternalNodes - 1
        n_e = float(self.curTree.nInternalNodes - 1)
        # nStar_p is the number of polytomies in the tree after the move.
        nStar_p = 0
        for n in pTree.iterInternalsNoRoot():
            if n.getNChildren() > 2:
                nStar_p += 1
        if pTree.root.getNChildren() > 3:
            nStar_p += 1
        nStar_p = float(nStar_p)
        # kStar is the number of edges emanating from the polytomy created (or enlarged) by the move.
        kStar = theNewParent.getNChildren()
        if theNewParent.parent:
            kStar += 1

        hastingsRatio = (gamma_D * n_e) / (nStar_p * (2**(kStar - 1) - kStar - 1))
        self.logProposalRatio = math.log(hastingsRatio)

        if 0:
            # Now the prior ratio.  The prior probability density f(nu) for a
            # branch length is lambda * exp(-lambda * nu).  To a first
            # approximation, with equal priors on topologies, the prior ratio
            # is 1/f(nu)
            priorRatio = 1.0/(self.mcmc.tunings.brLenPriorLambda * math.exp(- self.mcmc.tunings.brLenPriorLambda * theChosenNode.br.len))
            if dbug:
                print "The self.mcmc.tunings.brLenPriorLambda is %f" % self.mcmc.tunings.brLenPriorLambda
                print "So the prior ratio is %f" % priorRatio

            self.logPriorRatio = math.log(priorRatio)    

            # The Jacobian
            jacobian = self.mcmc.tunings.brLenPriorLambda * math.exp(- self.mcmc.tunings.brLenPriorLambda * theChosenNode.br.len)
            self.logJacobian = math.log(jacobian)
            print "logPriorRatio = %f, logJacobian = %f" % (self.logPriorRatio, self.logJacobian)

        # Here I pull a fast one, as explained in Lewis et al.  The
        # priorRatio and the Jacobian terms cancel out.  So the logs might
        # as well be zeros.
        self.logPriorRatio = 0.0
        #self.logJacobian = 0.0
        # That was easy, wasn't it?

        if self.stMcmc.tunings.doPolytomyResolutionClassPrior:
            # We are losing a node.  So the prior ratio is (T_{n,m} * C) /
            # T_{n,m - 1}.  We have the logs, and the result is the log.
            if 0:
                print "-" * 30
                print 'curTree.nInternalNodes', self.curTree.nInternalNodes
                print 'pTree.nInternalNodes', pTree.nInternalNodes
                print 'logBigT[curTree.nInternalNodes]', theProposal.logBigT[self.curTree.nInternalNodes]
                #print math.exp(theProposal.logBigT[self.curTree.nInternalNodes])
                print 'C ', self.stMcmc.tunings.polytomyPriorLogBigC
                print 'logBigT[pTree.nInternalNodes]', theProposal.logBigT[pTree.nInternalNodes]
                #print math.exp(theProposal.logBigT[pTree.nInternalNodes])
                print "-" * 30
            self.logPriorRatio = ((theProposal.logBigT[self.curTree.nInternalNodes] +
                                   self.stMcmc.tunings.polytomyPriorLogBigC) -
                                  theProposal.logBigT[pTree.nInternalNodes])

        else:
            if self.stMcmc.tunings.polytomyPriorLogBigC:
                self.logPriorRatio =  self.stMcmc.tunings.polytomyPriorLogBigC
            else:
                self.logPriorRatio = 0.0

        #print " losing a node, m %2i->%2i. logPriorRatio is %f" % (self.curTree.nInternalNodes,
        #                                                           pTree.nInternalNodes, self.logPriorRatio)

        
        
    def propose(self, theProposal):
        gm = ['STChain.propose()']
        #print "propose() About to propose %s" % theProposal.name


        if theProposal.name == 'nni':
            #self.proposeNni(theProposal)
            self.propTree.nni()             # this does setPreAndPostOrder()
            if theProposal.doAbort:
                pass
            #else:
            #    if not self.propTree.preAndPostOrderAreValid:    # not needed
            #        self.propTree.setPreAndPostOrder()
        elif theProposal.name == 'spr':
            self.propTree.randomSpr()
            if theProposal.doAbort:
                pass
            else:
                if not self.propTree.preAndPostOrderAreValid:
                    self.propTree.setPreAndPostOrder()

        elif theProposal.name == 'SR2008beta_uniform':
            mt = self.propTree.beta

            # Slider proposal
            mt += (random.random() - 0.5) * theProposal.tuning

            # Linear reflect
            isGood = False
            myMIN = 1.e-10
            myMAX = 1.e+10
            while not isGood:
                if mt < myMIN:
                    mt = (myMIN - mt) + myMIN
                elif mt > myMAX:
                    mt = myMAX - (mt - myMAX)
                else:
                    isGood = True
            self.propTree.beta = mt
            self.logProposalRatio = 0.0
            self.logPriorRatio = 0.0
        elif theProposal.name == 'spaQ_uniform':
            mt = self.propTree.spaQ
            #originally = mt
            # Slider proposal
            mt += (random.random() - 0.5) * theProposal.tuning

            # Linear reflect
            isGood = False
            myMIN = 1.e-10
            myMAX = 1.
            while not isGood:
                if mt < myMIN:
                    mt = (myMIN - mt) + myMIN
                elif mt > myMAX:
                    mt = myMAX - (mt - myMAX)
                else:
                    isGood = True
            self.propTree.spaQ = mt
            self.logProposalRatio = 0.0
            self.logPriorRatio = 0.0
            #print "proposing mt from %.3f to %.3f, diff=%g" % (originally, mt, mt-originally)

        elif theProposal.name == 'polytomy':
            self.proposePolytomy(theProposal)
            if not self.propTree.preAndPostOrderAreValid:
                self.propTree.setPreAndPostOrder()
            #self.propTree.draw()
            
        else:
            gm.append('Unlisted proposal.name=%s  Fix me.' % theProposal.name)
            raise Glitch, gm

        #return 0.0
        if theProposal.doAbort:
            return 0.0
        else:
            #print "...about to calculate the likelihood of the propTree.  Model %s" % self.stMcmc.modelName
            if self.stMcmc.modelName.startswith('SR2008_rf'):
                if self.stMcmc.stRFCalc == 'fastReducedRF':
                    self.getTreeLogLike_fastReducedRF()
                elif self.stMcmc.stRFCalc == 'purePython1':
                    self.getTreeLogLike_ppy1()
                elif self.stMcmc.stRFCalc == 'bitarray':
                    self.refreshBitarrayPropTree()
                    self.getTreeLogLike_bitarray()
            elif self.stMcmc.modelName == 'SPA':
                self.refreshBitarrayPropTree()
                self.getTreeLogLike_spa_bitarray()
            elif self.stMcmc.modelName == 'QPA':
                self.getTreeLogLike_qpa_slow()
            else:
                gm.append('Unknown model %s' % self.stMcmc.modelName)
                raise Glitch, gm

            #if theProposal.name == 'polytomy':
            #print "propTree logLike is %f, curTree logLike is %f" % (
            #    self.propTree.logLike, self.curTree.logLike)
            #myDist = self.propTree.topologyDistance(self.curTree)
            #print "myDist %2i, propTree.logLike %.3f  curTree.logLike %.3f " % (myDist, self.propTree.logLike, self.curTree.logLike)

            logLikeRatio = self.propTree.logLike - self.curTree.logLike
            #print logLikeRatio
            #logLikeRatio = 0.0

            theSum = logLikeRatio + self.logProposalRatio + self.logPriorRatio
            #theSum = self.logProposalRatio + self.logPriorRatio
            #if theProposal.name == 'polytomy':
            #    print "%f  %f  %f  %f" % (theSum, logLikeRatio, self.logProposalRatio, self.logPriorRatio)
            return theSum

    
    def gen(self, aProposal):
        gm = ['STChain.gen()']

        # doAborts means that it was not a valid generation,
        # neither accepted or rejected.  Give up, by returning True.

        acceptMove = False

        #print "Doing %s" % aProposal.name
        pRet = self.propose(aProposal)

        #print "pRet = %.6f" % pRet,
        if not aProposal.doAbort:
            if pRet < -100.0:  # math.exp(-100.) is 3.7200759760208361e-44
                r = 0.0
            elif pRet >= 0.0:
                r = 1.0
            else:
                r = math.exp(pRet)

            if r == 1.0:
                acceptMove = True
            elif random.random() < r:
                acceptMove = True
        
        #if aProposal.name == 'polytomy':
        #print "acceptMove = %s" % acceptMove
        #print "------------"
        #print " %6.0f" % pRet
        if 0 and acceptMove:
            d1 = self.propTree.topologyDistance(self.curTree, metric='scqdist')
            d2 = self.stMcmc.tree.topologyDistance(self.propTree, metric='scqdist')
            print " %6.0f    %5i   %5i  %5s" % (pRet, d1, d2, acceptMove)

        aProposal.nProposals[self.tempNum] += 1
        if acceptMove:
            aProposal.accepted = True
            aProposal.nAcceptances[self.tempNum]  += 1
            
        #if not aProposal.doAbort:
        if acceptMove:
            a = self.propTree
            b = self.curTree
        else:
            a = self.curTree
            b = self.propTree


        if aProposal.name in ['nni', 'spr', 'polytomy']:
            b.logLike = a.logLike
            a.copyToTree(b)
        elif aProposal.name in ['SR2008beta_uniform']:
            b.logLike = a.logLike
            b.beta = a.beta
        elif aProposal.name in ['spaQ_uniform']:
            b.logLike = a.logLike
            b.spaQ = a.spaQ

        else:
            gm.append('Unlisted proposal.name = %s  Fix me.' % aProposal.name)
            raise Glitch, gm



# for proposal probs
fudgeFactor = {}
fudgeFactor['local'] = 1.5



class STMcmcTunings(object):
    def __init__(self):
        object.__setattr__(self, 'chainTemp', 0.15)  # was 0.2
        object.__setattr__(self, 'nni', None)  
        object.__setattr__(self, 'spr', None)  
        object.__setattr__(self, 'SR2008beta_uniform', 0.2)  
        object.__setattr__(self, 'spaQ_uniform', 0.1)  
        object.__setattr__(self, 'doPolytomyResolutionClassPrior', False)
        object.__setattr__(self, 'polytomyPriorLogBigC', 0.0)


    def __setattr__(self, item, val):
        #print "Got request to set %s to %s" % (item, val)
        if item in self.__dict__.keys():
            # Here is where I should do the sanity checking of the new vals.  Some day.
            #print "    Setting tuning '%s' to %s" % (item, val)
            object.__setattr__(self, item, val)
        else:
            print self.dump()
            gm = ["\nSTMcmcTunings.__setattr__()"]
            gm.append("Can't set tuning '%s'-- no such tuning." % item)
            raise Glitch, gm
    
    def reprString(self, advice=True):
        lst = ["\nSTMcmc.tunings:"]
        spacer = ' ' * 4
        lst.append("%s%20s: %s" % (spacer, 'chainTemp', self.chainTemp))
        lst.append("%s%20s: %s" % (spacer, 'nni', self.nni))
        lst.append("%s%20s: %s" % (spacer, 'spr', self.spr))
        lst.append("%s%20s: %s" % (spacer, 'SR2008beta_uniform', self.SR2008beta_uniform))
        lst.append("%s%20s: %s" % (spacer, 'spaQ_uniform', self.spaQ_uniform))
        return string.join(lst, '\n')

    def dump(self):
        print self.reprString()

    def __repr__(self):
        return  self.reprString()


class STMcmcProposalProbs(dict):
    """User-settable relative proposal probabilities.

    An instance of this class is made as STMcmc.prob, where you can
    do, for example,
        yourSTMcmc.prob.nni = 2.0

    These are relative proposal probs, that do not sum to 1.0, and
    affect the calculation of the final proposal probabilities (ie the
    kind that do sum to 1).  It is a relative setting, and the default
    is 1.0.  Setting it to 0 turns it off.  For small
    probabilities, setting it to 2.0 doubles it.  For bigger
    probabilities, setting it to 2.0 makes it somewhat bigger.

    Check the effect that it has by doing a
        yourSTMcmc.writeProposalIntendedProbs()
    which prints out the final calculated probabilities. 
    """
    
    def __init__(self):
        object.__setattr__(self, 'nni', 1.0)
        object.__setattr__(self, 'spr', 1.0)
        object.__setattr__(self, 'SR2008beta_uniform', 1.0)
        object.__setattr__(self, 'spaQ_uniform', 1.0)
        object.__setattr__(self, 'polytomy', 0.0)


    def __setattr__(self, item, val):
        # complaintHead = "\nSTMcmcProposalProbs.__setattr__()"
        gm = ["\nSTMcmcProposalProbs(). (set %s to %s)" % (item, val)]
        theKeys = self.__dict__.keys()
        if item in theKeys:
            try:
                val = float(val)
                if val < 1e-9:
                    val = 0
                object.__setattr__(self, item, val)
            except:
                gm.append("Should be a float.  Got '%s'" % val)
                raise Glitch, gm
                
        else:
            self.dump()
            gm.append("    Can't set '%s'-- no such proposal." % item)
            raise Glitch, gm

    def reprString(self):
        stuff = ["\nUser-settable relative proposal probabilities, from yourMcmc.prob"]
        stuff.append("  To change it, do eg ")
        stuff.append("    yourMcmc.prob.comp = 0.0 # turns comp proposals off")
        stuff.append("  Current settings:")
        theKeys = self.__dict__.keys()
        theKeys.sort()
        for k in theKeys:
            stuff.append("        %20s: %s" % (k, getattr(self, k)))
        return string.join(stuff, '\n')

    def dump(self):
        print self.reprString()

    def __repr__(self):
        return  self.reprString()


                           
class STProposal(object):
    def __init__(self, theSTMcmc=None):
        self.name = None
        self.stMcmc = theSTMcmc            # reference loop!
        self.nChains = theSTMcmc.nChains
        self.pNum = -1
        self.mtNum = -1
        self.weight = 1.0
        self.nProposals = [0] * self.nChains
        self.nAcceptances = [0] * self.nChains
        self.accepted = 0
        self.doAbort = False
        self.nAborts = [0] * self.nChains

    def dump(self):
        print "proposal name=%-10s pNum=%2i, mtNum=%2i, weight=%5.1f, tuning=%7.2f" % (
            '%s,' % self.name, self.pNum, self.mtNum, self.weight, self.tuning)
        print "    nProposals   by temperature:  %s" % self.nProposals
        print "    nAcceptances by temperature:  %s" % self.nAcceptances
        
    def _getTuning(self):
        if self.name in ['nni', 'spr', 'SR2008beta_uniform', 'spaQ_uniform']:
            #print "getting tuning for %s, returning %f" % (self.name, getattr(self.mcmc.tunings, self.name))
            #print self.stMcmc.tunings
            return getattr(self.stMcmc.tunings, self.name)
        else:
            return None
        
    def _setTuning(self, whatever):
        raise Glitch, "Can't set tuning this way."
    def _delTuning(self):
        raise Glitch, "Can't del tuning."
    
    tuning = property(_getTuning, _setTuning, _delTuning) 
            
class BigTSplitStuff(object):
    # An organizer for splits on STMcmc.tree (ie bigT) internal nodes, only for use with bitarray
    def __init__(self):
        self.spl = None
        self.spl2 = None
        self.theSpl = None
        self.maskedSplitWithFirstTaxOne = None
        self.onesCount = None
        self.bytes = None
    def dump(self):
        print "ss: spl=%s, spl2=%s, masked=%s, onesCount=%s" % (
            self.spl, self.spl2, self.maskedSplitWithFirstTaxOne, self.onesCount)



class STMcmc(object):
    """An MCMC for making supertrees from a set of input trees.

This week, it implements the Steel and Rodrigo 2008 model, with the
alpha calculation using the approximation in Bryant and Steel 2009.

**Arguments**

inTrees
    A list of p4 tree objects.  You could just use ``var.trees``.

modelName
    The SR2008 models implemented here are based on the Steel and
    Rodrigo 2008 description of a likelihood model, "Maximum
    likelihood supertrees" Syst. Biol. 57(2):243--250, 2008.  At
    the moment, they are all SR2008_rf, meaning that they use
    Robinson-Foulds distances.

    SR2008_rf_ia 

        Here 'ia' means 'ignore alpha'.  The alpha values are not
        calculated at all, as they are presumed (erroneously, but
        not too badly) to cancel out.

    SR2008_rf_aZ

        This uses the approximation for Z_T = alpha^{-1} as described
        in Equation 30 in the Bryant and Steel paper "Computing the
        distribution of a tree metric" in IEEE/ACM Transactions on
        computational biology and bioinformatics, VOL. 6, 2009.

    SR2008_rf_aZ_fb

        This is as SR2008_rf_aZ above, but additionally it allows
        beta to be a free parameter, and it is sampled.  Samples
        are written to mcmc_prams* files.

beta
    This only applies to SR2008.  The beta is the weight as
    given in Steel and Rodrigo 2008. By default it is 1.0.


stRFCalc 

    There are three ways to calculate the RF distances and
    likelihood, for these SR2008_rf models above --- all giving
    the same answer.

    1.  purePython1.  Slow.

    2.  bitarray, using the bitarray module.  About twice as fast
        as purePython1

    3.  fastReducedRF, written in C++ using boost and ublas.
        About 10 times faster than purePython1, but perhaps a bit
        of a bother to get going.  It needs the fastReducedRF
        module, included in the p4 source code.

    It is under control of the argument stRFCalc, which can be one
    of 'purePython1', 'bitarray', and 'fastReducedRF'.  By default
    it is purePython1, so you may want to at least install
    bitarray.

runNum

    You may want to do more than one 'run' in the same directory,
    to facilitate convergence testing.  The first runNum would be
    0, and samples, likelihoods, and checkPoints are written to
    files with that number.

sampleInterval

    Interval at which the chain is sampled, including writing a tree,
    and the logLike.  Plan to get perhaps 1000 samples; so if you are
    planning to make a run of 10000 generations then you might set
    sampleInterval=10.

checkPointInterval

    Interval at which checkpoints are made.  If set to None (the
    default) it means don't make checkpoints.  My taste is to aim to
    make perhaps 2 to 4 per run.  So if you are planning to start out
    with a run of 10000 generations, you could set
    checkPointInterval=5000, which will give you 2 checkpoints.  See
    more about checkpointing below.

To prepare for a run, instantiate an Mcmc object, for example::

    m = STMcmc(treeList, modelName='SR2008_rf_aZ_fb', stRFCalc='fastReducedRF', sampleInterval=10)

To start it running, do this::

    # Tell it the number of generations to do
    m.run(10000)

As it runs, it saves trees and likelihoods at sampleInterval
intervals (actually whenever the current generation number is
evenly divisible by the sampleInterval).

**CheckPoints**

Whenever the current generation number is evenly divisible by the
checkPointInterval it will write a checkPoint file.  A checkPoint
file is the whole MCMC, pickled.  Using a checkPoint, you can
re-start an STMcmc from the point you left off.  Or, in the event
of a crash, you can restart from the latest checkPoint.  But the
most useful thing about them is that you can query checkPoints to
get information about how the chain has been running, and about
convergence diagnostics.

In order to restart the MCMC from the end of a previous run:: 

    # read the last checkPoint file
    m = func.unPickleStMcmc(0)  # runNum 0
    m.run(20000)

Its that easy if your previous run finished properly.  However, if
your previous run has crashed and you want to restart it from a
checkPoint, then you will need to repair the sample output files
to remove samples that were taken after the last checkPoint, but
before the crash.  Fix the trees, likelihoods, prams, and sims.
(You probably do not need to beware of confusing gen (eg 9999) and
gen+1 (eg 10000) issues.)  When you remove trees from the tree
files be sure to leave the 'end;' at the end-- p4 needs it, and
will deal with it.

The checkPoints can help with convergence testing.  To help with
that, you can use the STMcmcCheckPointReader class.  It will print
out a table of average standard deviations of split supports
between 2 runs, or between 2 checkPoints from the same run.  It
will print out tables of proposal acceptances to show whether they
change over the course of the MCMC.

**Making a consensus tree**

See :class:`TreePartitions`.

    """

    
    def __init__(self, inTrees, bigT=None, modelName='SR2008_rf_aZ', beta=1.0, spaQ=0.5, stRFCalc='purePython1', runNum=0, sampleInterval=100, checkPointInterval=None):
        gm = ['STMcmc.__init__()']

        assert inTrees
        for t in inTrees:
            assert isinstance(t, Tree)
        if bigT:
            assert isinstance(bigT, Tree)
            assert bigT.taxNames
            bigT.stripBrLens()
            for n in bigT.iterInternalsNoRoot():
                n.name = None

        goodModelNames = ['SR2008_rf_ia', 'SR2008_rf_aZ', 'SR2008_rf_aZ_fb', 'SPA', 'QPA']
        if modelName not in goodModelNames:
            gm.append("Arg modelName '%s' is not recognized. " % modelName)
            gm.append("Good modelNames are %s" % goodModelNames)
            raise Glitch, gm
        self.modelName = modelName

        self.stRFCalc = None
        if modelName.startswith("SR2008"):
            try:
                fBeta = float(beta)
            except ValueError:
                gm.append("Arg beta (%s) should be a float" % beta)
                raise Glitch, gm
            self.beta = fBeta

            for t in inTrees:
                if t.isFullyBifurcating():
                    pass
                else:
                    gm.append("At the moment STMcmc wants trees that are fully bifurcating.")
                    raise Glitch, gm

            goodSTRFCalcNames = ['purePython1', 'bitarray', 'fastReducedRF']
            if stRFCalc not in goodSTRFCalcNames:
                gm.append("Arg stRFCalc '%s' is not recognized. " % modelName)
                gm.append("Good stRFCalc names are %s" % goodSTRFCalcNames)
                raise Glitch, gm
            self.stRFCalc = stRFCalc

        if modelName in ['SPA', 'QPA']:
            try:
                fspaQ = float(spaQ)
            except ValueError:
                gm.append("Arg spaQ (%s) should be a float" % spaQ)
                raise Glitch, gm
            self.spaQ = fspaQ

        nChains = 1 # mcmcmc is off, temporarily
        try:
            nChains = int(nChains)
        except (ValueError,TypeError):
            gm.append("nChains should be an int, 1 or more.  Got %s" % nChains)
            raise Glitch, gm
        if nChains < 1:
            gm.append("nChains should be an int, 1 or more.  Got %s" % nChains)
            raise Glitch, gm
        self.nChains = nChains
        self.chains = []
        self.gen = -1
        self.startMinusOne = -1
        self.constraints = None
        self.simulate = None

        try:
            runNum = int(runNum)
        except (ValueError, TypeError):
            gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
            raise Glitch, gm
        if runNum < 0:
            gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
            raise Glitch, gm
        self.runNum = runNum

        # Check that we are not going to over-write good stuff
        ff = os.listdir(os.getcwd())
        hasPickle = False
        for fName in ff:
            if fName.startswith("mcmc_checkPoint_%i." % self.runNum):
                hasPickle = True
                break
        if hasPickle:
            gm.append("runNum is set to %i" % self.runNum)
            gm.append("There is at least one mcmc_checkPoint_%i.xxx file in this directory." % self.runNum)
            gm.append("This is a new STMcmc, and I am refusing to over-write exisiting files.")
            gm.append("Maybe you want to re-start from the latest mcmc_checkPoint_%i file?" % self.runNum)
            gm.append("Otherwise, get rid of the existing mcmc_xxx_%i.xxx files and start again." % self.runNum)
            raise Glitch, gm

        if var.strictRunNumberChecking:
            # We want to start runs with number 0, so if runNum is more than that, check that there are other runs.
            if self.runNum > 0:
                for runNum2 in range(self.runNum):
                    hasTrees = False
                    for fName in ff:
                        if fName.startswith("mcmc_trees_%i" % runNum2):
                            hasTrees = True
                            break
                    if not hasTrees:
                        gm.append("runNum is set to %i" % self.runNum)
                        gm.append("runNums should go from zero up.")
                        gm.append("There are no mcmc_trees_%i.nex files to show that run %i has been done." % (runNum2, runNum2))
                        gm.append("Set the runNum to that, first.")
                        raise Glitch, gm
                              
        
        
        self.sampleInterval = sampleInterval
        self.checkPointInterval = checkPointInterval

        self.proposals = []
        self.proposalsHash = {}
        self.propWeights = []
        self.cumPropWeights = []
        self.totalPropWeights = 0.0

        self.treePartitions = None
        self.likesFileName = "mcmc_likes_%i" % runNum
        self.treeFileName = "mcmc_trees_%i.nex" % runNum
        self.pramsFileName = "mcmc_prams_%i" % runNum
        self.writePrams = False
        if self.modelName in ['SR2008_rf_aZ_fb', "SPA", "QPA"]:
            self.writePrams = True

        self.lastTimeCheck = None

        if self.nChains > 1:
            self.swapMatrix = []
            for i in range(self.nChains):
                self.swapMatrix.append([0] * self.nChains)
        else:
            self.swapMatrix = None

        self.tunings = STMcmcTunings() 
        self.prob = STMcmcProposalProbs()                
        if self.modelName in ['SPA', 'QPA']:   
            self.prob.polytomy = 1.0
            self.prob.spr = 0.0
                
        # Zap internal node names
        # for n in aTree.root.iterInternals():
        #     if n.name:
        #         n.name = None

        if not bigT:
            allNames = []
            for t in inTrees:
                t.unsorted_taxNames = [n.name for n in t.iterLeavesNoRoot()]
                allNames += t.unsorted_taxNames          # Efficient?  Probably does not matter.
            self.taxNames = list(set(allNames))
            self.taxNames.sort()                         # not needed, but nice for debugging
        else:
            for t in inTrees:
                t.unsorted_taxNames = [n.name for n in t.iterLeavesNoRoot()]
            self.taxNames = bigT.taxNames
        #print self.taxNames
        self.nTax = len(self.taxNames)


        if self.modelName in ['SPA'] or self.stRFCalc == 'bitarray':
            #print "self.taxNames = ", self.taxNames
            for t in inTrees:
                #print "-" * 50
                #t.draw()
                sorted_taxNames = []
                t.baTaxBits = []
                for tNum in range(self.nTax):
                    tN = self.taxNames[tNum]
                    if tN in t.unsorted_taxNames:
                        sorted_taxNames.append(tN)
                        t.baTaxBits.append(True)
                    else:
                        t.baTaxBits.append(False)
                t.taxNames = sorted_taxNames
                t.baTaxBits = bitarray.bitarray(t.baTaxBits)
                t.firstTax = t.baTaxBits.index(1)
                #print "intree baTaxBits is %s" % t.baTaxBits
                #print "intree firstTax is %i" % t.firstTax

                # Can't use Tree.makeSplitKeys(), unfortunately.  So
                # make split keys here.  STMcmc.tBits is only used for
                # the leaves, here and in
                # STChain.setupBitarrayCalcs(), and there only once,
                # during STChain.__init__().  So probably does not
                # need to be an instance attribute.  Maybe delete?
                self.tBits = [False] * self.nTax
                for n in t.iterPostOrder():
                    if n == t.root:
                        break
                    if n.isLeaf:
                        spot = self.taxNames.index(n.name)
                        self.tBits[spot] = True
                        n.stSplitKey = bitarray.bitarray(self.tBits)
                        self.tBits[spot] = False
                    else:
                        n.stSplitKey = n.leftChild.stSplitKey.copy()
                        p = n.leftChild.sibling
                        while p:
                            n.stSplitKey |= p.stSplitKey    # "or", in-place
                            p = p.sibling
                        #print "setting node %i stSplitKey to %s" % (n.nodeNum, n.stSplitKey)
                t.splSet = set()
                for n in t.iterInternalsNoRoot():
                    if not n.stSplitKey[t.firstTax]:   # make sure splitKey[firstTax] is a '1'
                        n.stSplitKey.invert()
                        n.stSplitKey &= t.baTaxBits     # 'and', in-place
                        #print "inverting and and-ing node %i stSplitKey to %s" % (n.nodeNum, n.stSplitKey)
                    t.splSet.add(n.stSplitKey.tobytes()) # bytes so that I can use it as a set element

        if self.modelName in ['QPA']:
            for t in inTrees:
                sorted_taxNames = []
                t.taxBits = []
                for tNum in range(self.nTax):
                    tN = self.taxNames[tNum]
                    if tN in t.unsorted_taxNames:
                        sorted_taxNames.append(tN)
                        t.taxBits.append(1L << tNum)
                    else:
                        t.taxBits.append(0)
                t.taxNames = sorted_taxNames
                #print "intree taxBits is %s" % t.taxBits

                # Can't use Tree.makeSplitKeys(), unfortunately.  So
                # make split keys here.  STMcmc.tBits is only used for
                # the leaves, here and in
                # STChain.setupBitarrayCalcs(), and there only once,
                # during STChain.__init__().  So probably does not
                # need to be an instance attribute.  Maybe delete?
                #self.tBits = [False] * self.nTax
                for n in t.iterPostOrder():
                    if n == t.root:
                        break
                    if n.isLeaf:
                        spot = self.taxNames.index(n.name)
                        #self.tBits[spot] = True
                        n.stSplitKey = 1L << spot
                        #self.tBits[spot] = False
                    else:
                        n.stSplitKey = n.leftChild.stSplitKey
                        p = n.leftChild.sibling
                        while p:
                            n.stSplitKey |= p.stSplitKey    # "or", in-place
                            p = p.sibling
                        #print "setting node %i stSplitKey to %s" % (n.nodeNum, n.stSplitKey)
                # t.splSet = set()
                # for n in t.iterInternalsNoRoot():
                #     if not n.stSplitKey[t.firstTax]:   # make sure splitKey[firstTax] is a '1'
                #         n.stSplitKey.invert()
                #         n.stSplitKey &= t.baTaxBits     # 'and', in-place
                #         #print "inverting and and-ing node %i stSplitKey to %s" % (n.nodeNum, n.stSplitKey)
                #     t.splSet.add(n.stSplitKey.tobytes()) # bytes so that I can use it as a set element
                t.skk = [n.stSplitKey for n in t.iterInternalsNoRoot()]
                t.qSet = set()
                for sk in t.skk:
                    ups = [txBit for txBit in t.taxBits if (sk & txBit)]
                    downs = [txBit for txBit in t.taxBits if not (sk & txBit)]
                    for down in itertools.combinations(downs, 2):
                        if down[0] > down[1]:
                            down = (down[1], down[0])
                        for up in itertools.combinations(ups, 2):
                            if up[0] > up[1]:
                                up = (up[1], up[0])
                            if down[0] < up[0]:
                                t.qSet.add(down+up)
                            else:
                                t.qSet.add(up+down)
                #print t.qSet
                t.nQuartets = len(t.qSet)
            
        self.trees = inTrees
        if bigT:
            self.tree = bigT
        else:
            self.tree = func.randomTree(taxNames=self.taxNames, name='stTree', randomBrLens=False)
                

        if self.stRFCalc in ['purePython1', 'fastReducedRF']:
            for t in inTrees:
                sorted_taxNames = []
                t.taxBits = 0L
                for tNum in range(self.nTax):
                    tN = self.taxNames[tNum]
                    if tN in t.unsorted_taxNames:
                        sorted_taxNames.append(tN)
                        adder = 1L << tNum
                        t.taxBits += adder
                t.taxNames = sorted_taxNames
                t.allOnes = 2L**(t.nTax) - 1
                t.makeSplitKeys()
                t.skSet = set([n.br.splitKey for n in t.iterInternalsNoRoot()])

        if self.stRFCalc in  ['purePython1', 'fastReducedRF']:
            self.tree.makeSplitKeys()

            self.Frrf = None
            if self.stRFCalc == 'fastReducedRF':
                try:
                    import fastReducedRF
                    self.Frrf = fastReducedRF.Frrf
                    import pyublas # not explicitly used--but makes converters available
                except ImportError:
                    gm.append("var.stRFCalc is set to 'fastReducedRF', but I could not import")
                    gm.append("at least one of fastReducedRF or pyublas.")
                    gm.append("Make sure they are installed.")
                    raise Glitch, gm

        if self.modelName in ['QPA']:
            self.tree.taxBits = [1L << i for i in range(self.tree.nTax)]
            t = self.tree
            for n in t.iterPostOrder():
                if n == t.root:
                    break
                if n.isLeaf:
                    spot = self.taxNames.index(n.name)
                    n.stSplitKey = 1L << spot
                else:
                    n.stSplitKey = n.leftChild.stSplitKey
                    p = n.leftChild.sibling
                    while p:
                        n.stSplitKey |= p.stSplitKey    # "or", in-place
                        p = p.sibling
            t.skk = [n.stSplitKey for n in t.iterInternalsNoRoot()]
            t.qSet = set()
            for sk in t.skk:
                ups = [txBit for txBit in t.taxBits if (sk & txBit)]
                downs = [txBit for txBit in t.taxBits if not (sk & txBit)]
                for down in itertools.combinations(downs, 2):
                    assert down[0] < down[1]   # probably not needed
                    for up in itertools.combinations(ups, 2):
                        assert up[0] < up[1]  # probably not needed
                        if down[0] < up[0]:
                            t.qSet.add(down+up)
                        else:
                            t.qSet.add(up+down)
            #print t.qSet
            t.nQuartets = len(t.qSet)            

        print "Initializing STMcmc"
        print "%-10s: %s" % ('modelName', modelName)
        if self.modelName.startswith("SR2008"):
            print "%-10s: %s" % ('stRFCalc', self.stRFCalc)
        print "%-10s: %s" % ('inTrees', len(self.trees))
        print "%-10s: %s" % ('nTax', self.nTax)



    def _makeProposals(self):
        """Make proposals for the STMcmc."""

        gm = ['STMcmc._makeProposals()']


        # nni
        if self.prob.nni:
            p = STProposal(self)
            p.name = 'nni'
            p.weight = self.prob.nni # * (len(self.tree.nodes) - 1) * fudgeFactor['nni']
            self.proposals.append(p)

        if self.prob.spr:
            p = STProposal(self)
            p.name = 'spr'
            p.weight = self.prob.spr # * (len(self.tree.nodes) - 1) * fudgeFactor['nni']
            self.proposals.append(p)

        if self.modelName in ['SR2008_rf_aZ_fb']:
            if self.prob.SR2008beta_uniform:
                p = STProposal(self)
                p.name = 'SR2008beta_uniform'
                p.weight = self.prob.SR2008beta_uniform # * (len(self.tree.nodes) - 1) * fudgeFactor['nni']
                self.proposals.append(p)
                #object.__setattr__(self.tuningsUsage, 'local', p)

        if self.modelName in ['SPA', 'QPA']:
            if self.prob.spaQ_uniform:
                p = STProposal(self)
                p.name = 'spaQ_uniform'
                p.weight = self.prob.spaQ_uniform # * (len(self.tree.nodes) - 1) * fudgeFactor['nni']
                self.proposals.append(p)
                #object.__setattr__(self.tuningsUsage, 'local', p)
            if self.prob.polytomy:
                p = STProposal(self)
                p.name = 'polytomy'
                p.weight = self.prob.polytomy
                self.proposals.append(p)
                


        if not self.proposals:
            gm.append("No proposals?")
            raise Glitch, gm
        self.propWeights = []
        for p in self.proposals:
            self.propWeights.append(p.weight)
        self.cumPropWeights = [self.propWeights[0]]
        for i in range(len(self.propWeights))[1:]:
            self.cumPropWeights.append(self.cumPropWeights[i - 1] + self.propWeights[i])
        self.totalPropWeights = sum(self.propWeights)
        if self.totalPropWeights < 1e-9:
            gm.append("No proposal weights?")
            raise Glitch, gm
        for p in self.proposals:
            self.proposalsHash[p.name] = p

    def _refreshProposalProbsAndTunings(self):
        """Adjust proposals after a restart."""

        gm = ['STMcmc._refreshProposalProbsAndTunings()']

        for p in self.proposals:
            # nni
            if p.name == 'nni':
                #p.weight = self.prob.local * (len(self.tree.nodes) - 1) * fudgeFactor['local']
                p.weight = self.prob.nni


        self.propWeights = []
        for p in self.proposals:
            self.propWeights.append(p.weight)
        self.cumPropWeights = [self.propWeights[0]]
        for i in range(len(self.propWeights))[1:]:
            self.cumPropWeights.append(self.cumPropWeights[i - 1] + self.propWeights[i])
        self.totalPropWeights = sum(self.propWeights)
        if self.totalPropWeights < 1e-9:
            gm.append("No proposal weights?")
            raise Glitch, gm




    def writeProposalAcceptances(self):
        """Pretty-print the proposal acceptances."""

        if (self.gen - self.startMinusOne) <= 0:
            print "\nSTMcmc.writeProposalAcceptances()  There is no info in memory. "
            print " Maybe it was just emptied after writing to a checkpoint?  "
            print "If so, read the checkPoint and get the proposalAcceptances from there."
        else:

            spacer = ' ' * 8
            print "\nProposal acceptances, run %i, for %i gens, from gens %i to %i, inclusive." % (
                self.runNum, (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen)
            print "%s %20s %10s %13s%8s" % (spacer, 'proposal', 'nProposals', 'acceptance(%)', 'tuning')
            for p in self.proposals:
                print "%s" % spacer,
                print "%20s" % p.name,
                print "%10i" % p.nProposals[0],

                if p.nProposals[0]: # Don't divide by zero
                    print "       %5.1f " % (100.0 * float(p.nAcceptances[0]) / float(p.nProposals[0])),
                else:
                    print "           - ",

                if p.tuning == None:
                    print "      -",
                elif p.tuning < 2.0:
                    print "  %5.3f" % p.tuning,
                else:
                    print "%7.1f" % p.tuning,
                print

            # # Tabulate topology changes, if any were attempted.
            # doTopol = 0
            # p = None
            # try:
            #     p = self.proposalsHash['local']
            # except KeyError:
            #     pass
            # if p:
            #     for tNum in range(self.nChains):
            #         if p.nTopologyChangeAttempts[tNum]:
            #             doTopol = 1
            #             break
            #     if doTopol:
            #         p = self.proposalsHash['local']
            #         print "'Local' proposal-- attempted topology changes"
            #         print "%s tempNum   nProps nAccepts percent nTopolChangeAttempts nTopolChanges percent" % spacer
            #         for tNum in range(self.nChains):
            #             print "%s" % spacer,
            #             print "%4i " % tNum,
            #             print "%9i" % p.nProposals[tNum],
            #             print "%8i" % p.nAcceptances[tNum],
            #             print "  %5.1f" % (100.0 * float(p.nAcceptances[tNum]) / float(p.nProposals[tNum])),
            #             print "%20i" % p.nTopologyChangeAttempts[tNum],
            #             print "%13i" % p.nTopologyChanges[tNum],
            #             print "  %5.1f" % (100.0 * float(p.nTopologyChanges[tNum])/float(p.nTopologyChangeAttempts[tNum]))
            #     else:
            #         print "%sFor the 'local' proposals, there were no attempted" % spacer
            #         print "%stopology changes in any of the chains." % spacer


            # Check for aborts.
            # p = None
            # try:
            #     p = self.proposalsHash['local']
            # except KeyError:
            #     pass
            # if p:
            #     if hasattr(p, 'nAborts'):
            #         if p.nAborts[0]:
            #             print "The 'local' proposal had %i aborts." % p.nAborts[0]
            #             print "(Aborts might be due to brLen proposals too big or too small)"
            #             if self.constraints:
            #                 print "(Or, more likely, due to violated constraints.)"
            #         else:
            #             print "The 'local' proposal had no aborts (either due to brLen proposals"
            #             print "too big or too small, or due to violated constraints)."
            # for pN in ['polytomy', 'compLocation', 'rMatrixLocation', 'gdasrvLocation']:
            #     p = None
            #     try:
            #         p = self.proposalsHash[pN]
            #     except KeyError:
            #         pass
            #     if p:
            #         if hasattr(p, 'nAborts'):
            #             print "The %15s proposal had %5i aborts." % (p.name, p.nAborts[0])
            
            
            

    def writeSwapMatrix(self):
        print "\nChain swapping, for %i gens, from gens %i to %i, inclusive." % (
            (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen)
        print "    Swaps are presented as a square matrix, nChains * nChains."
        print "    Upper triangle is the number of swaps proposed between two chains."
        print "    Lower triangle is the percent swaps accepted."
        print "    The current tunings.chainTemp is %5.3f\n" % self.tunings.chainTemp
        print " " * 10,
        for i in range(self.nChains):
            print "%7i" % i,
        print
        print " " * 10,
        for i in range(self.nChains):
            print "   ----",
        print
        for i in range(self.nChains):
            print " " * 7, "%2i" % i,
            for j in range(self.nChains):
                if i < j: # upper triangle
                    print "%7i" % self.swapMatrix[i][j],
                elif i == j:
                    print "      -",
                else:
                    if self.swapMatrix[j][i] == 0: # no proposals
                        print "      -",
                    else:
                        print "  %5.1f" % (100.0 * float(self.swapMatrix[i][j]) / float(self.swapMatrix[j][i])),
            print
        

    def _makeChainsAndProposals(self):
        """Make chains and proposals."""

        gm = ['STMcmc._makeChainsAndProposals()']

        #random.seed(0)

        # Make chains, if needed
        if not self.chains:
            self.chains = []
            for chNum in range(self.nChains):
                aChain = STChain(self)
                aChain.tempNum = chNum  # Temperature.  Set this way to start, but it changes.
                self.chains.append(aChain)
        if not self.proposals:
            self._makeProposals()

            # If we are going to be doing the resolution class prior
            # in the polytomy move, we want to pre-compute the logs of
            # T_{n,m}.  Its a vector with indices (ie m) from zero to
            # nTax-2 inclusive.
            # if self.proposalsHash.has_key('polytomy') and self.tunings.doPolytomyResolutionClassPrior:
            #     p = self.proposalsHash['polytomy']
            #     bigT = func.nUnrootedTreesWithMultifurcations(self.tree.nTax)
            #     p.logBigT = [0.0] * (self.tree.nTax - 1)
            #     for i in range(1, self.tree.nTax - 1):
            #         p.logBigT[i] = math.log(bigT[i])
            #     #print p.logBigT

    def _setOutputTreeFile(self):
        """Setup the (output) tree file for the STMcmc."""

        gm = ['STMcmc._setOutputTreeFile()']

        # Write the preamble for the trees outfile.
        self.treeFile = file(self.treeFileName, 'w')
        self.treeFile.write('#nexus\n\n')
        self.treeFile.write('begin taxa;\n')
        self.treeFile.write('  dimensions ntax=%s;\n' % self.tree.nTax)
        self.treeFile.write('  taxlabels')
        for tN in self.tree.taxNames:
            self.treeFile.write(' %s' % func.nexusFixNameIfQuotesAreNeeded(tN))
        self.treeFile.write(';\nend;\n\n')

        self.treeFile.write('begin trees;\n')
        self.translationHash = {}
        i = 1
        for tName in self.tree.taxNames:
            self.translationHash[tName] = i
            i += 1

        self.treeFile.write('  translate\n')
        for i in range(self.tree.nTax - 1):
            self.treeFile.write('    %3i %s,\n' % (
                i + 1, func.nexusFixNameIfQuotesAreNeeded(self.tree.taxNames[i])))
        self.treeFile.write('    %3i %s\n' % (
            self.tree.nTax, func.nexusFixNameIfQuotesAreNeeded(self.tree.taxNames[-1])))
        self.treeFile.write('  ;\n')
        self.treeFile.write('  [Tree numbers are gen+1]\n')
        self.treeFile.close()


    def run(self, nGensToDo, verbose=True):
        """Start the STMcmc running."""

        gm = ['STMcmc.run()']

        #Keep track of the first gen of this call to run(), maybe restart
        firstGen = self.gen + 1

        if self.checkPointInterval:
            # We want a couple of things:
            #  1.  The last gen should be on checkPointInterval.  For
            #      example, if the checkPointInterval is 200, then doing
            #      100 or 300 generations will not be allowed cuz the
            #      chain would continue past the checkPoint-- bad.  Or if
            #      you re-start after 500 gens and change to a
            #      checkPointInterval of 200, then you won't be allowed to
            #      do 500 gens.
            #if ((self.gen + 1) + nGensToDo) % self.checkPointInterval == 0:
            if nGensToDo % self.checkPointInterval == 0:
                pass
            else:
                gm.append("With the current settings, the last generation won't be on a checkPointInterval.")
                gm.append("self.gen+1=%i, nGensToDo=%i, checkPointInterval=%i" % ((self.gen + 1),
                                                                                  nGensToDo, self.checkPointInterval)) 
                raise Glitch, gm
            #  2.  We also want the checkPointInterval to be evenly
            #      divisible by the sampleInterval.
            if self.checkPointInterval % self.sampleInterval == 0:
                pass
            else:
                gm.append("The checkPointInterval (%i) should be evenly divisible" % self.checkPointInterval)
                gm.append("by the sampleInterval (%i)." % self.sampleInterval)
                raise Glitch, gm


        if self.proposals:
            # Its either a re-start, or it has been thru autoTune().
            # I can tell the difference by self.gen, which is -1 after
            # autoTune()
            if self.gen == -1:
                self._makeChainsAndProposals()
                self._setOutputTreeFile()
                #if self.simulate:
                #    self.writeSimFileHeader(self.tree)
            # The probs and tunings may have been changed by the user.
            self._refreshProposalProbsAndTunings()

            # This stuff below should be the same as is done after pickling, see below.
            self.startMinusOne = self.gen

            # Start the tree partitions over.
            self.treePartitions = None
            # Zero the proposal counts
            for p in self.proposals:
                p.nProposals = [0] * self.nChains
                p.nAcceptances = [0] * self.nChains
                p.nTopologyChangeAttempts = [0] * self.nChains
                p.nTopologyChanges = [0] * self.nChains
            # Zero the swap matrix
            if self.nChains > 1:
                self.swapMatrix = []
                for i in range(self.nChains):
                    self.swapMatrix.append([0] * self.nChains)
            
        else:
            self._makeChainsAndProposals()
            self._setOutputTreeFile()
            #if self.simulate:
            #    self.writeSimFileHeader(self.tree)
        if verbose:
            self.writeProposalIntendedProbs()
            sys.stdout.flush()

        coldChainNum = 0
        
        # If polytomy is turned on, then it is possible to get a star
        # tree, in which case local will not work.  So if we have both
        # polytomy and local proposals, we should also have brLen.
        # if self.proposalsHash.has_key("polytomy") and self.proposalsHash.has_key("local"):
        #     if not self.proposalsHash.has_key('brLen'):
        #         gm.append("If you have polytomy and local proposals, you should have a brLen proposal as well.")
        #         gm.append("It can have a low proposal probability, but it needs to be there.")
        #         gm.append("Turn it on by eg yourMcmc.prob.brLen = 0.001")
        #         raise Glitch, gm

        
        if self.gen > -1:
            # it is a re-start, so we need to back over the "end;" in the tree files.
            f2 = file(self.treeFileName, 'a+')
            pos = -1
            while 1:
                f2.seek(pos, 2)
                c = f2.read(1)
                if c == ';':
                    break
                pos -= 1
            #print "pos now %i" % pos
            pos -= 3 # end;
            f2.seek(pos, 2)
            c = f2.read(4)
            #print "got c = '%s'" % c
            if c != "end;":
                gm.append("Mcmc.run().  Failed to find and remove the 'end;' at the end of the tree file.")
                raise Glitch, gm
            else:
                f2.seek(pos, 2)
                f2.truncate()
            f2.close()
                
            if verbose:
                print
                print "Re-starting the MCMC run %i from gen=%i" % (self.runNum, self.gen)
                print "Set to do %i more generations." % nGensToDo
                #if self.writePrams:
                #    if self.chains[0].curTree.model.nFreePrams == 0:
                #        print "There are no free prams in the model, so I am turning writePrams off."
                #        self.writePrams = False
                sys.stdout.flush()

            self.startMinusOne = self.gen
        else:
            if verbose:
                print "Starting the MCMC %s run %i" % ((self.constraints and "(with constraints)" or ""), self.runNum)
                print "Set to do %i generations." % nGensToDo
                if self.writePrams:
                    # if self.chains[0].curTree.model.nFreePrams == 0:
                    #     print "There are no free prams in the model, so I am turning writePrams off."
                    #     self.writePrams = False
                    # else:
                    pramsFile = file(self.pramsFileName, 'a')
                    if self.modelName.startswith("SR2008"):
                        pramsFile.write("    genPlus1     beta\n")
                    elif self.modelName.startswith("SPA"):
                        pramsFile.write("    genPlus1     spaQ\n")
                    elif self.modelName.startswith("QPA"):
                        pramsFile.write("    genPlus1     spaQ\n")
                    pramsFile.close()
                sys.stdout.flush()
            
        if verbose:
            print "Sampling every %i." % self.sampleInterval
            if self.checkPointInterval:
                print "CheckPoints written every %i." % self.checkPointInterval
            if nGensToDo <= 20000:
                print "One dot is 100 generations."
            else:
                print "One dot is 1000 generations."
            sys.stdout.flush()

        self.treePartitions = None
        realTimeStart = time.time()
        self.lastTimeCheck = time.time()

        abortableProposals = ['nni', 'spr', 'polytomy']

        for gNum in range(nGensToDo):
            self.gen += 1
            #Do an initial time estimate based on 100 gens 
            if nGensToDo > 100 and self.gen-firstGen == 100:
                diff_secs = time.time() - realTimeStart
                total_secs = (float(nGensToDo)/float(100))*float(diff_secs)
                deltaTime = datetime.timedelta(seconds = int(round(total_secs)))
                print "Estimated completion time: %s days, %s" % (
                    deltaTime.days, time.strftime("%H:%M:%S",time.gmtime(deltaTime.seconds)))

            # Above is a list of proposals where it is possible to abort.
            # When a gen(aProposal) is made, below, aProposal.doAbort
            # might be set, in which case we want to skip it for this
            # gen.  But we want to start each 'for chNum' loop with
            # doAborts all turned off.
            
            for chNum in range(self.nChains):
                failure = True
                nAttempts = 0
                while failure:
                    # Get the next proposal
                    gotIt = False
                    safety = 0
                    while not gotIt:
                        theRan = random.uniform(0.0, self.totalPropWeights)
                        for i in range(len(self.cumPropWeights)):
                            if theRan < self.cumPropWeights[i]:
                                break
                        aProposal = self.proposals[i]
                        gotIt = True
                        if aProposal.name == 'nni':
                            if self.chains[chNum].curTree.nInternalNodes == 1:  # Can't do nni on a star tree.
                                aProposal = self.proposalsHash['polytomy']
                        #elif aProposal.name == 'root3':
                        #    if self.chains[chNum].curTree.nInternalNodes == 1:  # Can't do root3 on a star tree.
                        #        gotIt = False
                        if aProposal.doAbort:
                            gotIt = False
                        safety += 1
                        if safety > 1000:
                            gm.append("Could not find a proposal after %i attempts." % safety)
                            gm.append("Possibly a programming error.")
                            gm.append("Or possibly it is just a pathologically frustrating Mcmc.")
                            raise Glitch, gm

                    #if gNum % 2:
                    #    aProposal = self.proposalsHash['brLen']
                    #else:
                    #    aProposal = self.proposalsHash['comp']

                    if 0:
                        print "==== gNum=%i, chNum=%i, aProposal=%s (part %i)" % (
                            gNum, chNum, aProposal.name, aProposal.pNum),
                        sys.stdout.flush()
                        #print gNum,

                    failure = self.chains[chNum].gen(aProposal)  # success returns None

                    if 0:
                        if failure:
                            print "    failure"
                        else:
                            print
                        
                    nAttempts += 1
                    if nAttempts > 1000:
                        gm.append("Was not able to do a successful generation after %i attempts." % nAttempts)
                        raise Glitch, gm
                #print "   Mcmc.run(). finished a gen on chain %i" % (chNum)
                for pr in abortableProposals:
                    if self.proposalsHash.has_key(pr):
                        self.proposalsHash[pr].doAbort = False
                        
            # Do swap, if there is more than 1 chain.
            if self.nChains == 1:
                coldChain = 0
            else:
                # Chain swapping stuff was lifted from MrBayes.  Thanks again.
                chain1,chain2 = random.sample(self.chains, 2)

                # Use the upper triangle of swapMatrix for nProposed's
                if chain1.tempNum < chain2.tempNum:
                    self.swapMatrix[chain1.tempNum][chain2.tempNum] += 1
                else:
                    self.swapMatrix[chain2.tempNum][chain1.tempNum] += 1

                lnR =  (1.0 / (1.0 + (self.tunings.chainTemp * chain1.tempNum))) * chain2.curTree.logLike
                lnR += (1.0 / (1.0 + (self.tunings.chainTemp * chain2.tempNum))) * chain1.curTree.logLike
                lnR -= (1.0 / (1.0 + (self.tunings.chainTemp * chain1.tempNum))) * chain1.curTree.logLike
                lnR -= (1.0 / (1.0 + (self.tunings.chainTemp * chain2.tempNum))) * chain2.curTree.logLike


                if lnR < -100.0:
                    r = 0.0
                elif lnR >= 0.0:
                    r = 1.0
                else:
                    r = math.exp(lnR)

                acceptSwap = 0
                if random.random() < r:
                    acceptSwap = 1

                if acceptSwap:
                    # Use the lower triangle of swapMatrix to keep track of nAccepted's
                    if chain1.tempNum < chain2.tempNum:
                        self.swapMatrix[chain2.tempNum][chain1.tempNum] += 1
                    else:
                        self.swapMatrix[chain1.tempNum][chain2.tempNum] += 1

                    # Do the swap
                    chain1.tempNum, chain2.tempNum = chain2.tempNum, chain1.tempNum

                # Find the cold chain, the one where tempNum is 0
                coldChainNum = -1
                for i in range(len(self.chains)):
                    if self.chains[i].tempNum == 0:
                        coldChainNum = i
                        break
                if coldChainNum == -1:
                    gm.append("Unable to find which chain is the cold chain.  Bad.")
                    raise Glitch, gm

            # If it is a writeInterval, write stuff
            if (self.gen + 1) % self.sampleInterval == 0:
                if 1:
                    likesFile = file(self.likesFileName, 'a')
                    likesFile.write('%11i %f\n' % (self.gen + 1, self.chains[coldChainNum].curTree.logLike))
                    likesFile.close()
                    treeFile = file(self.treeFileName, 'a')
                    treeFile.write("  tree t_%i = [&U] " % (self.gen + 1))
                    self.chains[coldChainNum].curTree.writeNewick(treeFile,
                                                               withTranslation=1,
                                                               translationHash=self.translationHash,
                                                               doMcmcCommandComments=False)
                    treeFile.close()

                if self.writePrams:
                    pramsFile = file(self.pramsFileName, 'a')
                    #pramsFile.write("%12i " % (self.gen + 1))
                    pramsFile.write("%12i" % (self.gen + 1))
                    if self.modelName.startswith("SR2008"):
                        pramsFile.write("  %f\n" % self.chains[coldChainNum].curTree.beta)
                    elif self.modelName in ["SPA", "QPA"]:
                        pramsFile.write("  %f\n" % self.chains[coldChainNum].curTree.spaQ)
                    pramsFile.close()

                # Do a simulation
                if self.simulate:
                    #print "about to simulate..."
                    self.doSimulate(self.chains[coldChainNum].curTree)
                    #print "...finished simulate."
                
                # Do other stuff.
                if hasattr(self, 'hook'):
                    self.hook(self.chains[coldChainNum].curTree)

                if 0 and self.constraints:
                    print "Mcmc x1c"
                    print self.chains[0].verifyIdentityOfTwoTreesInChain()
                    print "b checking curTree .."
                    self.chains[0].curTree.checkSplitKeys()
                    print "b checking propTree ..."
                    self.chains[0].propTree.checkSplitKeys()
                    print "Mcmc xxx"

                # Add curTree to treePartitions
                if self.treePartitions:
                    self.treePartitions._getSplitsFromTree(self.chains[coldChainNum].curTree)
                else:
                    self.treePartitions = TreePartitions(self.chains[coldChainNum].curTree)
                # After _getSplitsFromTree, need to follow, at some point,
                # with _finishSplits().  Do that when it is pickled, or at the end of the run.

                # Checking and debugging constraints
                if 0 and self.constraints:
                    print "Mcmc x1d"
                    print self.chains[coldChainNum].verifyIdentityOfTwoTreesInChain()
                    print "c checking curTree ..."
                    self.chains[coldChainNum].curTree.checkSplitKeys()
                    print "c checking propTree ..."
                    self.chains[coldChainNum].propTree.checkSplitKeys()
                    #print "c checking that all constraints are present"
                    #theSplits = [n.br.splitKey for n in self.chains[0].curTree.iterNodesNoRoot()]
                    #for sk in self.constraints.constraints:
                    #    if sk not in theSplits:
                    #        gm.append("split %i is not present in the curTree." % sk)
                    #        raise Glitch, gm
                    print "Mcmc zzz"

                # Check that the curTree has all the constraints
                if self.constraints:
                    splitsInCurTree = [n.br.splitKey for n in self.chains[coldChainNum].curTree.iterInternalsNoRoot()]
                    for sk in self.constraints.constraints:
                        if sk not in splitsInCurTree:
                            gm.append("Programming error.")
                            gm.append("The current tree (the last tree sampled) does not contain constraint")
                            gm.append("%s" % func.getSplitStringFromKey(sk, self.tree.nTax))
                            raise Glitch, gm

                # If it is a checkPointInterval, pickle
                if self.checkPointInterval and (self.gen + 1) % self.checkPointInterval == 0:
                    self.checkPoint()

                    # The stuff below needs to be done in a re-start as well.  See above "if self.proposals:"
                    self.startMinusOne = self.gen

                    # Start the tree partitions over.
                    self.treePartitions = None
                    # Zero the proposal counts
                    for p in self.proposals:
                        p.nProposals = [0] * self.nChains
                        p.nAcceptances = [0] * self.nChains
                        p.nTopologyChangeAttempts = [0] * self.nChains
                        p.nTopologyChanges = [0] * self.nChains
                        p.nAborts = [0] * self.nChains
                    # Zero the swap matrix
                    if self.nChains > 1:
                        self.swapMatrix = []
                        for i in range(self.nChains):
                            self.swapMatrix.append([0] * self.nChains)
                

            # Reassuring pips ...
            if firstGen != self.gen: #We want to skip the first gen of every call to run()
                if nGensToDo <= 20000:
                    if (self.gen-firstGen) % 1000 == 0:
                        if verbose:
                            deltaTime = self._doTimeCheck(nGensToDo, firstGen, 1000)
                            if deltaTime.days:
                                timeString = "%s days, %s" % (
                                    deltaTime.days, time.strftime("%H:%M:%S",time.gmtime(deltaTime.seconds)))
                            else:
                                timeString = time.strftime("%H:%M:%S",time.gmtime(deltaTime.seconds))
                            print "%10i - %s" % (self.gen, timeString)

                        else:
                            sys.stdout.write(".")
                            sys.stdout.flush()
                    elif (self.gen-firstGen) % 100 == 0:
                        sys.stdout.write(".")
                        sys.stdout.flush()
                else:
                    if (self.gen-firstGen) % 50000 == 0:
                        if verbose:
                            deltaTime = self._doTimeCheck(nGensToDo, firstGen, 50000)
                            if deltaTime.days:
                                timeString = "%s days, %s" % (
                                    deltaTime.days, time.strftime("%H:%M:%S",time.gmtime(deltaTime.seconds)))
                            else:
                                timeString = time.strftime("%H:%M:%S",time.gmtime(deltaTime.seconds))
                            print "%10i - %s" % (self.gen, timeString)
                        else:
                            sys.stdout.write(".")
                            sys.stdout.flush()
                    elif (self.gen-firstGen) % 1000 == 0:
                        sys.stdout.write(".")
                        sys.stdout.flush()

        # Gens finished.  Clean up.
        print
        if verbose:
            print "Finished %s generations." % nGensToDo

        treeFile = file(self.treeFileName, 'a')
        treeFile.write('end;\n\n')
        treeFile.close()
        
    def _doTimeCheck(self, nGensToDo, firstGen, genInterval):
        """Time check 
        
        firstGen is the first generation of this call to Mcmc.run() else
        timing fails on restart"""
        nowTime = time.time()
        diff_secs = nowTime - self.lastTimeCheck
        total_secs = (float(nGensToDo-(self.gen-firstGen))/float(genInterval))*float(diff_secs)
        deltaTime = datetime.timedelta(seconds = int(round(total_secs)))
        self.lastTimeCheck = nowTime
        return deltaTime


    def checkPoint(self):
        # Maybe we should not save the inTrees? -- would make it more lightweight.
        if 0:
            for chNum in range(self.nChains):
                ch = self.chains[chNum]
                print "chain %i ==================" % chNum
                ch.curTree.summarizeModelThingsNNodes()

        # the Frrf object does not pickle
        savedFrrfs = []
        savedBigTrs = []
        if self.stRFCalc == 'fastReducedRF':
            for chNum in range(self.nChains):
                ch = self.chains[chNum]
                savedFrrfs.append(ch.frrf)
                ch.frrf = None
                savedBigTrs.append(ch.bigTr)
                ch.bigTr = None
            
        theCopy = copy.deepcopy(self)

        theCopy.treePartitions._finishSplits()
        theCopy.likesFile = None
        theCopy.treeFile = None
        #theCopy.treePartitions = None    # this can be the biggest part of the pickle.

        # Pickle it.
        fName = "mcmc_checkPoint_%i.%i" % (self.runNum, self.gen + 1)
        f = file(fName, 'w')
        cPickle.dump(theCopy, f, 1)
        f.close()

        if self.stRFCalc == 'fastReducedRF':
            for chNum in range(self.nChains):
                ch = self.chains[chNum]
                ch.frrf = savedFrrfs[chNum]
                ch.bigTr = savedBigTrs[chNum]
       

    

    def writeProposalProbs(self):
        """(Another) Pretty-print the proposal probabilities.

        See also STMcmc.writeProposalAcceptances().
        """

        nProposals = len(self.proposals)
        if not nProposals:
            print "STMcmc.writeProposalProbs().  No proposals (yet?)."
            return
        #intended = self.propWeights[:]
        #for i in range(len(intended)):
        #    intended[i] /= self.totalPropWeights
        #if math.fabs(sum(intended) - 1.0 > 1e-15):
        #    raise Glitch, 'bad sum of intended proposal probs. %s' % sum(intended)

        nAttained = [0] * nProposals
        nAccepted = [0] * nProposals
        for i in range(nProposals):
            nAttained[i] = self.proposals[i].nProposals[0]
            nAccepted[i] = self.proposals[i].nAcceptances[0]
        sumAttained = float(sum(nAttained)) # should be zero or nGen
        if not sumAttained:
            print "STMcmc.writeProposalProbs().  No proposals have been made."
            print "Possibly, due to it being a checkPoint interval, nProposals have all been set to zero."
            return
        #assert int(sumAttained) == self.gen + 1, "sumAttained is %i, should be gen+1, %i." % (
        #    int(sumAttained), self.gen + 1)
        probAttained = []
        for i in range(len(nAttained)):
            probAttained.append(100.0 * float(nAttained[i]) / sumAttained)
        if math.fabs(sum(probAttained) - 100.0 > 1e-13):
            raise Glitch, 'bad sum of attained proposal probs. %s' % sum(probAttained)

        spacer = ' ' * 4
        print "\nProposal probabilities (%)"
        #print "There are %i proposals" % len(self.proposals)
        print "For %i gens, from gens %i to %i, inclusive." % (
            (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen)
        print "%2s %11s %11s  %11s %10s %23s %5s %5s" % ('', 'nProposals', 'proposed(%)',
                                                    'accepted(%)', 'tuning', 'proposal', 'part', 'num')
        for i in range(len(self.proposals)):
            print "%2i" % i,
            p = self.proposals[i]
            print "   %7i " % self.proposals[i].nProposals[0],
            print "   %5.1f    " % probAttained[i],
            if nAttained[i]:
                print "   %5.1f   " % (100.0 * float(nAccepted[i]) / float(nAttained[i])),
            else:
                print "       -   ",

            if p.tuning == None:
                print "       -    ",
            elif p.tuning < 2.0:
                print "   %7.3f  " % p.tuning,
            else:
                print "   %7.1f  " % p.tuning,

            print " %20s" % p.name,
            if p.pNum != -1:
                print " %3i " % p.pNum,
            else:
                print "   - ",
            if p.mtNum != -1:
                print " %3i " % p.mtNum,
            else:
                print "   - ",
            print


    def writeProposalIntendedProbs(self):
        """Tabulate the intended proposal probabilities.
        """

        nProposals = len(self.proposals)
        if not nProposals:
            print "STMcmc.writeProposalIntendedProbs().  No proposals (yet?)."
            return
        intended = self.propWeights[:]
        for i in range(len(intended)):
            intended[i] /= self.totalPropWeights
        if math.fabs(sum(intended) - 1.0 > 1e-14):
            raise Glitch, 'bad sum of intended proposal probs. %s' % sum(intended)

        spacer = ' ' * 4
        print "\nIntended proposal probabilities (%)"
        #print "There are %i proposals" % len(self.proposals)
        print "%2s %11s %23s %5s %5s" % ('', 'intended(%)', 'proposal', 'part', 'num')
        for i in range(len(self.proposals)):
            print "%2i" % i,
            p = self.proposals[i]
            print "   %6.2f    " % (100. * intended[i]),

            print " %20s" % p.name,
            if p.pNum != -1:
                print " %3i " % p.pNum,
            else:
                print "   - ",
            if p.mtNum != -1:
                print " %3i " % p.mtNum,
            else:
                print "   - ",
            print

class STMcmcCheckPointReader(object):
    """Read in and display mcmc_checkPoint files.

    Three options--
    
    To read in a specific checkpoint file, specify the file name by
    fName=whatever

    To read in the most recent (by os.path.getmtime()) checkpoint
    file, say last=True
    
    If you specify neither of the above, it will read in all the
    checkPoint files that it finds.

    Where it looks is determined by theGlob, which by default is '*',
    ie everything in the current directory.  If you want to look
    somewhere else, you can specify eg

        theGlob='SomeWhereElse/*' 

    or, if it is unambiguous, just

        theGlob='S*/*' 

    So you might say

        cpr = STMcmcCheckPointReader(theGlob='*_0.*')

    to get all the checkpoints from the first run, run 0.  Then, you
    can tell the cpr object to do various things.  Eg

        cpr.writeProposalAcceptances()

    But perhaps the most powerful thing about it is that it allows
    easy access to the checkpointed Mcmc objects, in the list mm.  Eg
    to get the first one, ask for

        m = cpr.mm[0]

    and m is an STMcmc object, complete with all its records of
    proposals and acceptances and so on.  And the TreePartitions
    object.  

    (Sorry!  -- Lazy documentation.  See the source code for more that it can do.)
    """
    
    def __init__(self, fName=None, theGlob='*', last=False, verbose=True):
        self.mm = []
        if not fName:
            #fList = [fName for fName in os.listdir(os.getcwd()) if fName.startswith("mcmc_checkPoint")]
            #fList = glob.glob(theGlob)
            #print "Full glob = %s" % fList
            fList = [fName for fName in glob.glob(theGlob) if 
                     os.path.basename(fName).startswith("mcmc_checkPoint")]
            #print fList
            if not fList:
                raise Glitch, "No checkpoints found in this directory."
            if last:
                # Find the most recent
                mostRecent = os.path.getmtime(fList[0])
                mostRecentFileName = fList[0]
                if len(fList) > 1:
                    for fName in fList[1:]:
                        mtime = os.path.getmtime(fName)
                        if mtime > mostRecent:
                            mostRecent = mtime
                            mostRecentFileName = fName
                f = file(mostRecentFileName)
                m = cPickle.load(f)
                f.close()
                self.mm.append(m)
                
            else:
                # get all the files
                for fName in fList:
                    f = file(fName)
                    m = cPickle.load(f)
                    f.close()
                    self.mm.append(m)

                self.mm = func.sortListOfObjectsOn2Attributes(self.mm, "gen", 'runNum')
        else:
            # get the file by name
            f = file(fName)
            m = cPickle.load(f)
            f.close()
            self.mm.append(m)
        if verbose:
            self.dump()

    def dump(self):
        print "STMcmcCheckPoints (%i checkPoints read)" % len(self.mm)
        print "%12s %12s %12s %12s" % (" ", "index", "run", "gen+1")
        print "%12s %12s %12s %12s" % (" ", "-----", "---", "-----")
        for i in range(len(self.mm)):
            m = self.mm[i]
            #print "    %2i    run %2i,  gen+1 %11i" % (i, m.runNum, m.gen+1)
            print "%12s %12s %12s %12s" % (" ", i, m.runNum, m.gen+1)
            

    def compareSplits(self, mNum1, mNum2, verbose=True, minimumProportion=0.1):
        """Should we be only looking at splits within the 95% ci of the topologies?"""
        m1 = self.mm[mNum1]
        m2 = self.mm[mNum2]
        tp1 = m1.treePartitions
        tp2 = m2.treePartitions
        
        if verbose:
            print "\nSTMcmcCheckPointReader.compareSplits(%i,%i)" % (mNum1, mNum2)
            print "%12s %12s %12s %12s %12s" % ("mNum", "runNum", "start", "gen+1", "nTrees")
            for i in range(5):
                print "   ---------",
            print
            for mNum in [mNum1, mNum2]:
                print " %10i " % mNum,
                m = self.mm[mNum]            
                print " %10i " % m.runNum,
                print " %10i " % (m.startMinusOne + 1),
                print " %10i " % (m.gen + 1),
                #for i in m.splitCompares:
                #    print i
                print " %10i " % m.treePartitions.nTrees

        asdos = self.compareSplitsBetweenTwoTreePartitions(tp1, tp2, minimumProportion, verbose=verbose)
        if asdos == None and verbose:
                print "No splits > %s" % minimumProportion
        return asdos

        
    def compareSplitsBetweenTwoTreePartitions(tp1, tp2, minimumProportion, verbose=False):
        ret = tp1.compareSplits(tp2, minimumProportion=minimumProportion)
        if ret != []:
            sumOfStdDevs = 0.0
            diffs = []
            if ret and len(ret):
                nSplits = len(ret)
                for i in ret:
                    #print "            %.3f  %.3f    " % (i[2][0], i[2][1]),
                    stdDev = math.sqrt(func.variance(i[2]))
                    #print "%.5f" % stdDev
                    sumOfStdDevs += stdDev
                    diffs.append(math.fabs(i[2][0] - i[2][1]))
                if verbose:
                    #print "  %f " % sumOfStdDevs,
                    print "     nSplits=%i, average of std devs of splits %.4f " % (nSplits, sumOfStdDevs/nSplits)
                    print "     max difference %f, mean difference %f" % (max(diffs), sum(diffs)/nSplits)
            return sumOfStdDevs/nSplits
        else:
            return None

    compareSplitsBetweenTwoTreePartitions = staticmethod(compareSplitsBetweenTwoTreePartitions)

    def compareSplitsAll(self):
        nM = len(self.mm)
        nItems = ((nM * nM) - nM)/2
        results = np.zeros((nM, nM), np.float)
        vect = np.zeros(nItems, np.float)
        vCounter = 0
        for mNum1 in range(1, nM):
            for mNum2 in range(mNum1):
                ret = self.compareSplits(mNum1, mNum2, verbose=False)
                #print "+++ ret = %s" % ret
                if ret == None:
                    ret = 0.0
                results[mNum1][mNum2] = ret
                results[mNum2][mNum1] = ret
                vect[vCounter] = ret
                vCounter += 1
                if 0:
                    print " %10i " % mNum1,
                    print " %10i " % mNum2,
                    print "%.3f" % ret
        print results
        
        print "For the %i values in one triangle," % nItems
        print "max =  ", vect.max()
        print "min =  ", vect.min()
        print "mean = ", vect.mean()
        print "var =  ", vect.var()
            
        
    def writeProposalAcceptances(self):
        for m in self.mm:
            m.writeProposalAcceptances()

    def writeSwapMatrices(self):
        for m in self.mm:
            if m.nChains > 1:
                m.writeSwapMatrix()
            
    def writeProposalProbs(self):
        for m in self.mm:
            m.writeProposalProbs()
