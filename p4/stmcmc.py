# This is STMcmc, for super tree mcmc.
# Started 18 March 2011, first commit 22 March 2011.

import p4.pf as pf
import p4.func
from p4.var import var
import math
import random
import string
import sys
import time
import copy
import os
import pickle
import glob
import numpy as np
from p4.p4exceptions import P4Error
from p4.treepartitions import TreePartitions
from p4.constraints import Constraints
from p4.tree import Tree
import datetime
import itertools
from scipy.optimize import minimize
import logging
import bitarray


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

    myLambda = cT / (2.0 * n)
    tester = 0.5 * math.log((n - 3.) / myLambda)

    epsilon = math.exp(-2. * beta)
    bigANEpsilon = 1 + (((2. * n) - 3.) * epsilon) + \
        (2. * ((n * n) - (4. * n) - 6.) * epsilon * epsilon)
    termA = math.log(bigANEpsilon + 6 * cT * epsilon * epsilon)

    if beta < tester:
        termB = -(2. * beta) * (n - 3.) + \
            (myLambda * (math.exp(2. * beta) - 1.))
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
        tester = 1 << i
        if tester > k:
            return count
        if tester & k:
            count += 1
    return count


def bitReduce(bk, txBits, lLen, sLen, allOnes):
    # print "bitReduce: bk %i, txBits %i, lLen %i, sLen %i, allOnes %i" % (bk,
    # txBits, lLen, sLen, allOnes)
    newBk = 0
    counter = 0
    pops = 0
    for pos in range(lLen):
        tester = 1 << pos
        # print "pos %2i, tester: %3i" % (pos, tester)
        if tester & txBits:
            # print "    tester & txBits -- True"
            if tester & bk:
                adder = 1 << counter
                # print "        adding:", adder
                newBk += adder
                pops += 1
            else:
                # print "        not adding"
                pass
            counter += 1
    if (1 & newBk):
        # print "flipping"
        newBk = allOnes ^ newBk
        pops = sLen - pops
    # print "returning newBk %i, pops %i" % (newBk, pops)
    return newBk, pops

if 0:  # test bitReduce
    sk = 6   # always at least 2 bits, even
    txBits = 30
    lLen = 5
    sLen = 4
    allOnes = 15
    print("     sk: %3i  %s" % (sk, p4.func.getSplitStringFromKey(sk, lLen)))
    print("taxBits: %3i  %s" % (txBits, p4.func.getSplitStringFromKey(txBits, lLen)))

    rsk, popcount = bitReduce(sk, txBits, lLen, sLen, allOnes)
    print("    rsk: %3i  %s" % (rsk, p4.func.getSplitStringFromKey(rsk, sLen)))
    print("   popcount %i" % popcount)
    #     sk:   6  .**..
    #     taxBits:  30  .****
    #     rsk:  12  ..**
    #     popcount 2


def maskedSymmetricDifference(skk, skSet, taxBits, longLen, shortLen, allOnes):
    if 0:
        print("-" * 50)
        print("skk (skk_ppy1 from the current supertree)")
        for sk in skk:
            print(p4.func.getSplitStringFromKey(sk, longLen))
        print("skSet (from input tree)")
        for sk in skSet:
            print(p4.func.getSplitStringFromKey(sk, shortLen))
        print("taxBits:", taxBits, p4.func.getSplitStringFromKey(taxBits, longLen))

    newSkk = []
    for sk in skk:
        reducedSk, popcount = bitReduce(
            sk, taxBits, longLen, shortLen, allOnes)
        if 0:
            print("taxBits: %s  " % p4.func.getSplitStringFromKey(taxBits, longLen), end=' ')
            print("%4i %s  " % (sk, p4.func.getSplitStringFromKey(sk, longLen)), end=' ')
            print("%4i %s  %i" % (reducedSk, p4.func.getSplitStringFromKey(reducedSk, shortLen), popcount))
        if popcount <= 1 or popcount >= (shortLen - 1):
            pass
        else:
            newSkk.append(reducedSk)
    newSkkSet = set(newSkk)
    # print newSkkSet, skSet
    # print "reduced supertree splits =  newSkkSet = %s" % newSkkSet
    ret = len(newSkkSet.symmetric_difference(skSet))
    # print "symmetric difference %i" % ret
    nCherries = 0
    for sk in newSkkSet:
        popcount = popcountA(sk, shortLen)
        if popcount == 2:
            nCherries += 1
        # not "elif", because they might both be True
        if popcount == (shortLen - 2):
            nCherries += 1
    # print "nCherries %i" % nCherries
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

    def __init__(self, aSTMcmc, chNum):
        gm = ['STChain.__init__()']

        self.stMcmc = aSTMcmc
        self.chNum = chNum    # Does not change.  Used with fastspa only
        self.tempNum = chNum  # 'temp'erature, not 'temp'orary;  changes when swapped.

        if chNum == 0:
            self.curTree = aSTMcmc.tree.dupe()
            self.propTree = aSTMcmc.tree.dupe()
        else:      # heated chains are randomized, unless var.mcmc_sameBigTToStartOnAllChains is set.
            if var.mcmc_sameBigTToStartOnAllChains:  # False by default
                self.curTree = aSTMcmc.tree.dupe()
                self.propTree = aSTMcmc.tree.dupe()
            else:
                rTree = aSTMcmc.tree.dupe()
                rTree.randomizeTopology(randomBrLens = False)
                rTree.stripBrLens()
                rTree.setPreAndPostOrder()
                self.curTree = rTree
                self.propTree = rTree.dupe()

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0

        self.frrf = None
        #self.nInTreeSplits = 0

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
            self.curTree.spaQ = np.array([self.stMcmc.spaQ])
            self.propTree.spaQ = np.array([self.stMcmc.spaQ])

            # for t in self.stMcmc.trees:
            #    self.nInTreeSplits += len(t.splSet)
            # print "Got nInTreeSplits %s" % self.nInTreeSplits
            self.setupBitarrayCalcs()
            self.getTreeLogLike_spa_bitarray()
            if var.stmcmc_useFastSpa:
                #print("Here E.  bitarray propTree.logLike is %f" % self.propTree.logLike)
                fspaLike = self.stMcmc.fspa.calcLogLike(self.chNum)
                diff = math.fabs(self.propTree.logLike - fspaLike)
                #print("Got fspaLike %f, diff %g" % (fspaLike, diff))
                if diff > 1e-13:
                    gm.append("bad fastspa likelihood calc, %f vs %f, diff %f" % (self.propTree.logLike, fspaLike, diff))
                    raise P4Error(gm)

            self.curTree.logLike = self.propTree.logLike
        elif self.stMcmc.modelName.startswith('QPA'):
            self.curTree.spaQ = self.stMcmc.spaQ
            self.propTree.spaQ = self.stMcmc.spaQ
            self.nPossibleQuartets = choose(self.stMcmc.tree.nTax, 4) * 3
            self.getTreeLogLike_qpa_slow()
            self.curTree.logLike = self.propTree.logLike

        else:
            gm.append('Unknown modelName %s' % self.stMcmc.modelName)
            raise P4Error(gm)

        if 0:
            print("STChain init()")
            self.curTree.draw()
            print("logLike is %f" % self.curTree.logLike)

    def getTreeLogLike_qpa_slow(self):
        gm = ["STChain.getTreeLogLike_qpa_slow()"]
        if self.propTree.spaQ > 1. or self.propTree.spaQ <= 0.0:
            gm.append("bad propTree.spaQ value %f" % self.propTree.spaQ)
            raise P4Error(gm)

        
        for n in self.propTree.iterInternalsPostOrder():
            if n == self.propTree.root:
                break
            n.stSplitKey = n.leftChild.stSplitKey
            p = n.leftChild.sibling
            while p:
                n.stSplitKey |= p.stSplitKey    # "or", in-place
                p = p.sibling
        self.propTree.skk = [
            n.stSplitKey for n in self.propTree.iterInternalsNoRoot()]
        self.propTree.qSet = set()
        for sk in self.propTree.skk:
            ups = [txBit for txBit in self.propTree.taxBits if (sk & txBit)]
            downs = [
                txBit for txBit in self.propTree.taxBits if not (sk & txBit)]
            for down in itertools.combinations(downs, 2):
                if down[0] > down[1]:
                    down = (down[1], down[0])
                for up in itertools.combinations(ups, 2):
                    if up[0] > up[1]:
                        up = (up[1], up[0])
                    if down[0] < up[0]:
                        self.propTree.qSet.add(down + up)
                    else:
                        self.propTree.qSet.add(up + down)
        # print self.propTree.qSet
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
        if self.propTree.spaQ[0] > 1. or self.propTree.spaQ[0] <= 0.0:
            gm.append("bad propTree.spaQ value %f" % self.propTree.spaQ)
            raise P4Error(gm)
        slowCheck = False
        if slowCheck:
            print("\n", "-" * 30)
            print("Super tree: ", end=' ')
            self.propTree.write()
            slowCheckLogLike = 0.0
            for it in self.stMcmc.trees:
                it.makeSplitKeys()
                it.inbb = [n.br for n in it.iterInternalsNoRoot()]

        self.propTree.logLike = 0.0
        #sumOfLogqs = 0.0
        #sumOfLogrs = 0.0
        # self.propTree.draw()
        for it in self.stMcmc.trees:
            if 0:
                print("-" * 50)
                it.draw()
                print("baTaxBits %s" % it.baTaxBits)
                print("firstTax at %i" % it.firstTax)
            if 0:
                print("  input tree: ", end=' ')
                it.write()

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
                    print("bigT node %i" % n.nodeNum)
                    print("  theSpl is %s" % n.ss.theSpl)
                    print("  maskedSplitWithTheFirstTaxOne %s" % n.ss.maskedSplitWithTheFirstTaxOne)
                    print("  onesCount %i" % n.ss.onesCount)
                    if n.ss.onesCount >= 2 and n.ss.onesCount <= upperGood:
                        print("    -> relevant")
                    else:
                        print("    -> not relevant")
                if n.ss.onesCount >= 2 and n.ss.onesCount <= upperGood:
                    relevantStSplits.append(n.ss)

            nonRedundantStSplitDict = {}
            for ss in relevantStSplits:
                ss.bytes = ss.maskedSplitWithTheFirstTaxOne.tobytes()
                nonRedundantStSplitDict[ss.bytes] = ss

            if 0:
                for ss in relevantStSplits:
                    ss.dump()
                print("There are %i relevant splits in the st for this it." % len(relevantStSplits))
                for ss in nonRedundantStSplitDict:
                    ss.dump()
                print("There are %i non-redundant splits in the st for this it." % len(nonRedundantStSplitDict))
            if 0:
                if(len(relevantStSplits) != len(nonRedundantStSplitDict)):
                    print("Gen %12i: Got %i relevantStSplits; %i in nonRedundantStSplitDict" % (
                        self.stMcmc.gen, len(relevantStSplits),len(nonRedundantStSplitDict)))

            # S_st is the number of splits in the reduced supertree
            S_st = len(nonRedundantStSplitDict)
            if slowCheck:
                # stDupe.draw()
                # print "the drawing above is stDupe"
                slowCheckS_st = len([n for n in stDupe.iterInternalsNoRoot()])
                assert S_st == slowCheckS_st

            # S is the number of possible splits in an it-sized tree
            S = 2 ** (it.nTax - 1) - (it.nTax + 1)
            # print("    S=%i, S_st=%i, S_x=%i" % (S, S_st, S-S_st))
            if S_st:
                q = self.propTree.spaQ[0] / S_st
                R = 1. - self.propTree.spaQ[0]
                r = R / (S - S_st)
                #print("q=%g" % q)
                logq = math.log(q)
            else:
                R = 1.
                r = R / S
            #print ("r=%g" % r)
            logr = math.log(r)

            for n in it.internals:
                ret = nonRedundantStSplitDict.get(n.stSplitKeyBytes)
                if ret:
                    if self.stMcmc.useSplitSupport and n.br.support != None:
                        thisSplitLike = math.log(r + (n.br.support * (q - r)))
                        self.propTree.logLike += thisSplitLike
                        # sumOfLogqs += thisSplitLike
                    else:
                        self.propTree.logLike += logq
                        # sumOfLogqs += logq
                else:
                    # If we are here when S_st is zero, then q is undefined.
                    # So fall into the else clause
                    if 0:
                        # This has the dodgy assumption that 1-support is in the supertree.
                        # Might be zero support for a split in the supertree.
                        if self.stMcmc.useSplitSupport and S_st and n.br.support != None:
                            self.propTree.logLike += math.log(
                                r + ((1. - n.br.support) * (q - r)))
                        else:
                            self.propTree.logLike += logr
                    else:
                        # This does not make the assumption above.  Safer.
                        self.propTree.logLike += logr
                        # sumOfLogrs += logr

            if slowCheck:
                for inb in it.inbb:
                    splitString = p4.func.getSplitStringFromKey(
                        inb.splitKey, it.nTax)
                    print("    %s " % splitString, end=' ')
                    ret = stDupe.nodeForSplitKeyDict.get(inb.splitKey)
                    # Here we need to check that S_st is not zero.  If it is,
                    # then q is undefined.
                    if self.stMcmc.useSplitSupport and inb.support != None and S_st:
                        if ret:
                            thisLogQc = math.log(r + (inb.support * (q - r)))
                            print("qc %.3f" % thisLogQc)
                            slowCheckLogLike += thisLogQc
                        else:
                            thisLogRc = math.log(
                                r + ((1. - inb.support) * (q - r)))
                            print("rc %.3f" % thisLogRc)
                            slowCheckLogLike += thisLogRc
                    else:
                        if ret:
                            print("q %.3f" % q)
                            slowCheckLogLike += logq
                        else:
                            print("r %.3f" % r)
                            slowCheckLogLike += logr
        if 1:
            if slowCheck:
                # print self.propTree.logLike, slowCheckLogLike
                myDiff = self.propTree.logLike - slowCheckLogLike
                if math.fabs(myDiff) > 1.e-12:
                    gm.append("Bad like calc. slowCheck %f, bitarray %f, diff %g" % (
                        slowCheckLogLike, self.propTree.logLike, myDiff))
                    raise P4Error(gm)
        # print("sumOfLogqs = %g, sumOfLogrs = %g" % (sumOfLogqs, sumOfLogrs))

    def setupBitarrayCalcs(self):
        # Prepare self.propTree (ie bigT).  First make n.stSplitKeys.  These
        # are temporary; the info is held more permanently in n.ss, a BigTSplitStuff object
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

        if self.stMcmc.modelName.startswith('SPA') and var.stmcmc_useFastSpa:
            self.stMcmc.fspa.setBigT(len(self.propTree.nodes), self.propTree.nTax, 
                                     self.propTree.postOrder, self.propTree.spaQ)
            for nNum in self.propTree.postOrder:
                if nNum == -10000:
                    break
                n = self.propTree.nodes[nNum]
                if n == self.propTree.root or n.isLeaf:
                    theSpl = '0'
                    theSpl2 = '0'
                else:
                    theSpl = n.ss.spl.to01()
                    theSpl2 = n.ss.spl2.to01()
                self.stMcmc.fspa.setBigTNoSpl(self.chNum, nNum, theSpl, theSpl2)


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

        if self.stMcmc.modelName.startswith('SPA') and var.stmcmc_useFastSpa:
            for nNum in self.propTree.postOrder:
                if nNum == -10000:
                    break
                n = self.propTree.nodes[nNum]
                if n == self.propTree.root or n.isLeaf:
                    theSpl = '0'
                    theSpl2 = '0'
                else:
                    theSpl = n.ss.spl.to01()
                    theSpl2 = n.ss.spl2.to01()
                self.stMcmc.fspa.setBigTNoSpl(self.chNum, nNum, theSpl, theSpl2)

    def startFrrf(self):
        # if using self.stMcmc.stRFCalc= 'fastReducedRF'
        self.frrf = self.stMcmc.Frrf(len(self.stMcmc.taxNames))
        self.bigTr = self.frrf.setBigT(
            len(self.propTree.nodes), self.propTree.nTax, self.propTree.postOrder)

        for n in self.propTree.nodes:
            if n.parent:
                self.bigTr.setParent(n.nodeNum, n.parent.nodeNum)
            if n.leftChild:
                self.bigTr.setLeftChild(n.nodeNum, n.leftChild.nodeNum)
            else:
                self.bigTr.setNodeTaxNum(
                    n.nodeNum, self.stMcmc.taxNames.index(n.name))
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
                        tr.setNodeTaxNum(
                            n.nodeNum, self.stMcmc.taxNames.index(n.name))
                    if n.sibling:
                        tr.setSibling(n.nodeNum, n.sibling.nodeNum)
        self.frrf.setInTreeTaxBits()
        self.frrf.setInTreeInternalBits()
        self.frrf.maybeFlipInTreeBits()
        self.frrf.setBigTInternalBits()
        # self.frrf.dump()

    def getTreeLogLike_ppy1(self):
        gm = ['STChain.getTreeLogLike_pp1']
        self.propTree.makeSplitKeys()
        self.propTree.skk = [
            n.br.splitKey for n in self.propTree.iterInternalsNoRoot()]
        self.propTree.logLike = 0.0
        for t in self.stMcmc.trees:

            # Get the distance
            thisDist = None
            if self.stMcmc.modelName.startswith('SR2008_rf'):
                thisDist, nCherries = maskedSymmetricDifference(self.propTree.skk, t.skSet,
                                                                t.taxBits, self.stMcmc.nTax, t.nTax, t.allOnes)
            else:
                raise P4Error(
                    "STChain.getTreeLogLike_ppy1() unknown model '%s'" % self.stMcmc.modelName)

            # Now multiply by beta, and do approximate Z_T
            assert thisDist != None
            beta_distance = self.propTree.beta * thisDist
            if self.stMcmc.modelName == 'SR2008_rf_ia':
                self.propTree.logLike -= beta_distance
            elif self.stMcmc.modelName.startswith('SR2008_rf_aZ'):
                log_approxZT = BS2009_Eqn30_ZTApprox(
                    t.nTax, self.propTree.beta, nCherries)
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
                raise P4Error(gm)

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
            # else:
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
                raise P4Error(gm)

    def getTreeLogLike_bitarray(self):
        self.propTree.logLike = 0.0
        slowCheck = False
        if slowCheck:
            self.propTree.makeSplitKeys()
            self.propTree.skk = [
                n.br.splitKey for n in self.propTree.iterInternalsNoRoot()]
        for t in self.stMcmc.trees:
            if 0:
                print("-" * 50)
                t.draw()
                print("baTaxBits %s" % t.baTaxBits)
                print("firstTax at %i" % t.firstTax)
            # splitStuff objects with onesCount >= 2 and <= t.nTax = 2
            usables = []
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
                    print("bigT node %i" % n.nodeNum)
                    print("  theSpl is %s" % n.ss.theSpl)
                    print("  maskedSplitWithTheFirstTaxOne %s" % n.ss.maskedSplitWithTheFirstTaxOne)
                    print("  onesCount %i" % n.ss.onesCount)
                    if n.ss.onesCount >= 2 and n.ss.onesCount <= upperGood:
                        print("    -> used")
                    else:
                        print("    -> not used")
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
            thisBaRF = len(splSet.symmetric_difference(t.splSet))
            if slowCheck:  # with purePython1
                thisPPyRF, thisPPyNCherries = maskedSymmetricDifference(self.propTree.skk, t.skSet,
                                                                        t.taxBits, self.stMcmc.nTax, t.nTax, t.allOnes)
                if thisBaRF != thisPPyRF:
                    raise P4Error("bitarray and purePython1 RF calcs differ.")
            beta_distance = self.propTree.beta * thisBaRF
            if self.stMcmc.modelName == 'SR2008_rf_ia':
                self.propTree.logLike -= beta_distance
            elif self.stMcmc.modelName.startswith('SR2008_rf_aZ'):
                nCherries = 0
                for ba in splSet:
                    theSS = usablesDict[ba]
                    # theSS.dump()
                    if theSS.onesCount == 2:
                        nCherries += 1
                    if theSS.onesCount == upperGood:
                        nCherries += 1
                if slowCheck:
                    if nCherries != thisPPyNCherries:
                        raise P4Error(
                            "bitarray and purePython1 nCherries calcs differ.")
                log_approxZT = BS2009_Eqn30_ZTApprox(
                    t.nTax, self.propTree.beta, nCherries)
                self.propTree.logLike -= log_approxZT
                self.propTree.logLike -= beta_distance
            else:
                gm.append("Unknown model %s" % self.stMcmc.modelName)
                raise P4Error(gm)

    def proposePolytomy(self, theProposal):
        theProposal.doAbort = False
        dbug = False
        if dbug:
            # print "proposePolytomy() starting with this tree ..."
            #self.propTree.draw(width=80, addToBrLen=0.2)
            print("j There are %i internal nodes." % self.propTree.nInternalNodes)
            if self.propTree.nInternalNodes == 1:
                print("-> so its a star tree -> proposeDeleteEdge is not possible.")
            elif self.propTree.nInternalNodes == self.propTree.nTax - 2:
                print("-> so its a fully-resolved tree, so proposeAddEdge is not possible.")

        if self.propTree.nInternalNodes == 1:  # a star tree
            self.proposeAddEdge(theProposal)
        elif self.propTree.nInternalNodes == self.propTree.nTax - 2:
            candidateNodes = self._getCandidateNodesForDeleteEdge()
            if candidateNodes:
                self.proposeDeleteEdge(theProposal, candidateNodes)
            else:
                #gm = ["proposePolytomy()"]
                #gm.append("The tree is fully resolved, so I can't proposeAddEdge()")
                #gm.append("But there are no suitable nodes to remove.")
                #raise P4Error(gm)
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

    def proposeAddEdge(self, theProposal):
        gm = ["STChain.proposeAddEdge()"]
        # print "proposeAddEdge() here"
        dbug = False
        pTree = self.propTree
        if 0:
            print("proposeAddEdge(), starting with this tree ...")
            pTree.draw()
            print("k There are %i internal nodes." % pTree.nInternalNodes)
            print("root is node %i" % pTree.root.nodeNum)
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
            # Yes, all children.
            childrenNodeNums = pTree.getChildrenNums(theChosenPolytomy)

        nPossibleWays = math.pow(2, k - 1) - k - 1
        if dbug:
            print("These nodes are polytomies: %s" % [n.nodeNum for n in allPolytomies])
            print("We randomly choose to do node %i" % theChosenPolytomy.nodeNum)
            print("It has %i children, so k=%i, so there are %i possible ways to add a node." % (
                nChildren, k, nPossibleWays))

        # We want to choose one of the possible ways to add a node, but we
        # want to choose it randomly.  I'll describe it for the case with
        # nChildren=5, so k is 6.  We know already that there are
        # nPossibleWays=25 different ways to add a node.  The complication
        # is that we could make a new group of 2, 3, or 4 nInNewGroup, and it will be
        # different numbers of possible ways in each.  The numbers of each are given by
        # p4.func.nChoosek(), so there are 10 ways to make a group of 2 from 5
        # children, 10 ways to make a group of 3 from 5 children, and 5
        # ways to make a group of 4 from 5 children.  So thats [10, 10,
        # 5], which sums to 25 (nPossibleWays).  So we can make a
        # cumulative sum list ie [10, 20, 25], and use it to choose one
        # group randomly.
        nChooseKs = []
        for i in range(2, nChildren):
            nChooseKs.append(p4.func.nChooseK(nChildren, i))
        cumSum = [nChooseKs[0]]
        for i in range(len(nChooseKs))[1:]:
            cumSum.append(nChooseKs[i] + cumSum[i - 1])
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
            print("The nChooseKs are %s" % nChooseKs)
            print("The cumSum is %s" % cumSum)
            print("Since there are nPossibleWays=%i, we choose a random number from 0-%i" % (
                nPossibleWays, nPossibleWays - 1))
            print("->We chose a random number: %i" % ran)
            print("So we choose the group at index %i, which means nInNewGroup=%i" % (i, nInNewGroup))
            print("So we make a new node with newChildrenNodeNums %s" % newChildrenNodeNums)
            # sys.exit()

        # Choose to add a node between theChosenPolytomy and the first in
        # the list of newChildrenNodeNums.  The node that we add will be
        # chosen from pTree.nodes for the first node where both the parent
        # and the leftChild are None.
        firstNode = pTree.nodes[newChildrenNodeNums[0]]
        for newNode in pTree.nodes:
            if not newNode.parent and not newNode.leftChild:
                break
        # print "Got newNode = %i" % newNode.nodeNum

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
            # pTree.setPreAndPostOrder()
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
        #             allOnes = 2**(self.propTree.nTax) - 1
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
        # If the proposed tree is fully resolved and the current tree is not
        # the star tree
        elif (pTree.nInternalNodes == (len(pTree.nodes) - 2)) and (self.curTree.nInternalNodes > 1):
            gamma_B = 2.0
        else:
            gamma_B = 1.0

        # n_e is number of internal edges present before the Add-edge move.
        # That would be self.curTree.nInternalNodes - 1
        n_e = float(self.curTree.nInternalNodes - 1)
        # n_p is the number of polytomies present before the move,
        # len(allPolytomies)
        n_p = float(len(allPolytomies))
        hastingsRatio = (gamma_B * n_p * float(nPossibleWays)) / (1.0 + n_e)

        if dbug:
            print("The new node is given a random branch length of %f" % newNode.br.len)
            print("For the Hastings ratio ...")
            print("gamma_B is %.1f" % gamma_B)
            print("n_e is %.0f" % n_e)
            print("k is (still) %i, and (2^{k-1} - k - 1) = nPossibleWays is still %i" % (k, nPossibleWays))
            print("n_p = %.0f is the number of polytomies present before the move." % n_p)
            print("So the hastings ratio is %f" % hastingsRatio)

        self.logProposalRatio = math.log(hastingsRatio)

        if 0:
            priorRatio = theProposal.brLenPriorLambda * \
                math.exp(- theProposal.brLenPriorLambda * newNode.br.len)
            if dbug:
                print("The theProposal.brLenPriorLambda is %f" % theProposal.brLenPriorLambda)
                print("So the prior ratio is %f" % priorRatio)

            self.logPriorRatio = math.log(priorRatio)

            # The Jacobian
            jacobian = 1.0 / (theProposal.brLenPriorLambda *
                              math.exp(- theProposal.brLenPriorLambda * newNode.br.len))
            self.logJacobian = math.log(jacobian)
            print("logPriorRatio = %f, logJacobian = %f" % (self.logPriorRatio, self.logJacobian))

        # Here I pull a fast one, as explained in Lewis et al.  The
        # priorRatio and the Jacobian terms cancel out.  So the logs might
        # as well be zeros.
        self.logPriorRatio = 0.0
        #self.logJacobian = 0.0
        # That was easy, wasn't it?
        if theProposal.polytomyUseResolutionClassPrior:
            # We are gaining a node.  So the prior ratio is T_{n,m + 1} /
            # (T_{n,m} * C) .  We have the logs, and the result is the
            # log.
            if 0:
                print("-" * 30)
                print('curTree.nInternalNodes', self.curTree.nInternalNodes)
                print('pTree.nInternalNodes', pTree.nInternalNodes)
                print('logBigT[curTree.nInternalNodes]', theProposal.logBigT[self.curTree.nInternalNodes])
                # print
                # math.exp(theProposal.logBigT[self.curTree.nInternalNodes])
                print('C ', theProposal.polytomyPriorLogBigC)
                print('logBigT[pTree.nInternalNodes]', theProposal.logBigT[pTree.nInternalNodes])
                # print math.exp(theProposal.logBigT[pTree.nInternalNodes])
                print("-" * 30)
            self.logPriorRatio = (theProposal.logBigT[self.curTree.nInternalNodes] -
                                  (theProposal.polytomyPriorLogBigC +
                                   theProposal.logBigT[pTree.nInternalNodes]))

        else:
            if theProposal.polytomyPriorLogBigC:
                self.logPriorRatio = -theProposal.polytomyPriorLogBigC
            else:
                self.logPriorRatio = 0.0
        # print "gaining a node, m %2i->%2i. logPriorRatio is %f" % (self.curTree.nInternalNodes,
        # pTree.nInternalNodes, self.logPriorRatio)

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
        # print "doing proposeDeleteEdge()"
        if 0:
            print("proposeDeleteEdge(), starting with this tree ...")
            pTree.draw()
            print("m There are %i internal nodes (before deleting the edge)." % pTree.nInternalNodes)

        if not candidateNodes:
            raise P4Error(
                "proposeDeleteEdge() could not find a good node to attempt to delete.")

        theChosenNode = random.choice(candidateNodes)
        if dbug:
            print("There are %i candidateNodes." % len(candidateNodes))
            print("node nums %s" % [n.nodeNum for n in candidateNodes])
            print("Randomly choose node %s" % theChosenNode.nodeNum)

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
        # print pTree.preOrder
        # if dbug:
        #    pTree.draw()

        # Hastings ratio.  First calculate the gamma_D.  If the current
        # tree is fully resolved and the proposed tree is not the star
        # tree, then gamma_D is 0.5
        if (self.curTree.nInternalNodes == len(pTree.nodes) - 2) and pTree.nInternalNodes != 1:
            gamma_D = 0.5
        # If the proposed tree is the star tree and the current tree is not
        # fully resolved
        elif (self.curTree.nInternalNodes < len(pTree.nodes) - 2) and pTree.nInternalNodes == 1:
            gamma_D = 2.
        else:
            gamma_D = 1.

        # n_e is the number of internal edges in existence before the move,
        # which would be nInternalNodes - 1
        n_e = float(self.curTree.nInternalNodes - 1)
        # nStar_p is the number of polytomies in the tree after the move.
        nStar_p = 0
        for n in pTree.iterInternalsNoRoot():
            if n.getNChildren() > 2:
                nStar_p += 1
        if pTree.root.getNChildren() > 3:
            nStar_p += 1
        nStar_p = float(nStar_p)
        # kStar is the number of edges emanating from the polytomy created (or
        # enlarged) by the move.
        kStar = theNewParent.getNChildren()
        if theNewParent.parent:
            kStar += 1

        hastingsRatio = (gamma_D * n_e) / \
            (nStar_p * (2 ** (kStar - 1) - kStar - 1))
        self.logProposalRatio = math.log(hastingsRatio)

        if 0:
            # Now the prior ratio.  The prior probability density f(nu) for a
            # branch length is lambda * exp(-lambda * nu).  To a first
            # approximation, with equal priors on topologies, the prior ratio
            # is 1/f(nu)
            priorRatio = 1.0 / (theProposal.brLenPriorLambda *
                                math.exp(- theProposal.brLenPriorLambda * theChosenNode.br.len))
            if dbug:
                print("The theProposal.brLenPriorLambda is %f" % theProposal.brLenPriorLambda)
                print("So the prior ratio is %f" % priorRatio)

            self.logPriorRatio = math.log(priorRatio)

            # The Jacobian
            jacobian = theProposal.brLenPriorLambda * \
                math.exp(- theProposal.brLenPriorLambda *
                         theChosenNode.br.len)
            self.logJacobian = math.log(jacobian)
            print("logPriorRatio = %f, logJacobian = %f" % (self.logPriorRatio, self.logJacobian))

        # Here I pull a fast one, as explained in Lewis et al.  The
        # priorRatio and the Jacobian terms cancel out.  So the logs might
        # as well be zeros.
        self.logPriorRatio = 0.0
        #self.logJacobian = 0.0
        # That was easy, wasn't it?

        if theProposal.polytomyUseResolutionClassPrior:
            # We are losing a node.  So the prior ratio is (T_{n,m} * C) /
            # T_{n,m - 1}.  We have the logs, and the result is the log.
            if 0:
                print("-" * 30)
                print('curTree.nInternalNodes', self.curTree.nInternalNodes)
                print('pTree.nInternalNodes', pTree.nInternalNodes)
                print('logBigT[curTree.nInternalNodes]', theProposal.logBigT[self.curTree.nInternalNodes])
                # print
                # math.exp(theProposal.logBigT[self.curTree.nInternalNodes])
                print('C ', theProposal.polytomyPriorLogBigC)
                print('logBigT[pTree.nInternalNodes]', theProposal.logBigT[pTree.nInternalNodes])
                # print math.exp(theProposal.logBigT[pTree.nInternalNodes])
                print("-" * 30)
            self.logPriorRatio = ((theProposal.logBigT[self.curTree.nInternalNodes] +
                                   theProposal.polytomyPriorLogBigC) -
                                  theProposal.logBigT[pTree.nInternalNodes])

        else:
            if theProposal.polytomyPriorLogBigC:
                self.logPriorRatio = theProposal.polytomyPriorLogBigC
            else:
                self.logPriorRatio = 0.0

        # print " losing a node, m %2i->%2i. logPriorRatio is %f" % (self.curTree.nInternalNodes,
        # pTree.nInternalNodes, self.logPriorRatio)

    def propose(self, theProposal):
        gm = ['STChain.propose()']
        #print("propose() About to propose %s" % theProposal.name)

        if theProposal.name == 'nni':
            # self.proposeNni(theProposal)
            self.propTree.nni()             # this does setPreAndPostOrder()
            if theProposal.doAbort:
                pass
            # else:
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
            mt += (random.random() - 0.5) * theProposal.tuning[self.tempNum]

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
            mt = self.propTree.spaQ[0]
            originally = mt
            # Slider proposal
            mt += (random.random() - 0.5) * theProposal.tuning[self.tempNum]

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
            self.propTree.spaQ[0] = mt
            self.logProposalRatio = 0.0
            if theProposal.spaQPriorType == 'flat':
                self.logPriorRatio = 0.0
            elif theProposal.spaQPriorType == 'exponential':
                self.logPriorRatio = theProposal.spaQExpPriorLambda * (originally - mt)
            else:
                raise P4Error("this should not happen! wxyzz")
            # print "proposing mt from %.3f to %.3f, diff=%g" % (originally,
            # mt, mt-originally)

        elif theProposal.name == 'polytomy':
            self.proposePolytomy(theProposal)
            if not self.propTree.preAndPostOrderAreValid:
                self.propTree.setPreAndPostOrder()
            # self.propTree.draw()

        else:
            gm.append('Unlisted proposal.name=%s  Fix me.' % theProposal.name)
            raise P4Error(gm)

        if theProposal.doAbort:
            return 0.0

        # print "...about to calculate the likelihood of the propTree.
        # Model %s" % self.stMcmc.modelName
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
            if var.stmcmc_useFastSpa:
                if 0:  # check
                    self.getTreeLogLike_spa_bitarray()
                    #print("Here F.  bitarray propTree.logLike is %f" % self.propTree.logLike)
                    fspaLike = self.stMcmc.fspa.calcLogLike(self.chNum)
                    diff = math.fabs(self.propTree.logLike - fspaLike)
                    #print("Got fspaLike %f, diff %g" % (fspaLike, diff))
                    if diff > 1e-13:
                        gm.append("gen %i bad fastspa likelihood calc, %f vs %f, diff %f" % (
                            self.stMcmc.gen, self.propTree.logLike, fspaLike, diff))
                        raise P4Error(gm)
                else:
                    self.propTree.logLike = self.stMcmc.fspa.calcLogLike(self.chNum)
            else:
                self.getTreeLogLike_spa_bitarray()

        elif self.stMcmc.modelName == 'QPA':
            self.getTreeLogLike_qpa_slow()
        else:
            gm.append('Unknown model %s' % self.stMcmc.modelName)
            raise P4Error(gm)

        # if theProposal.name == 'polytomy':
        #print("propTree logLike is %f, curTree logLike is %f" % (
        #    self.propTree.logLike, self.curTree.logLike))
        #myDist = self.propTree.topologyDistance(self.curTree)
        # print "myDist %2i, propTree.logLike %.3f  curTree.logLike %.3f "
        # % (myDist, self.propTree.logLike, self.curTree.logLike)

        logLikeRatio = self.propTree.logLike - self.curTree.logLike
        # print logLikeRatio

        # To run "without the data", which shows the effect of priors.
        #logLikeRatio = 0.0

        # Mcmcmc
        if self.stMcmc.nChains > 1:
            if self.stMcmc.swapVector:
                heatBeta = 1.0 / (1.0 + self.stMcmc.chainTemps[self.tempNum])
            else:
                heatBeta = 1.0 / (1.0 + self.stMcmc.chainTemp * self.tempNum)
            logLikeRatio *= heatBeta
            self.logPriorRatio *= heatBeta
            #print("propose().  chainTemp=%s, heatBeta=%f" % (self.stMcmc.chainTemps, heatBeta))

        # Experimental Heating hack
        if self.stMcmc.doHeatingHack: # and theProposal.name in self.stMcmc.heatingHackProposalNames:
            heatFactor = 1.0 / (1.0 + self.stMcmc.heatingHackTemperature)
            logLikeRatio *= heatFactor
            self.logPriorRatio *= heatFactor

        theSum = logLikeRatio + self.logProposalRatio + self.logPriorRatio
        #theSum = self.logProposalRatio + self.logPriorRatio
        # if theProposal.name == 'polytomy':
        # print "%f  %f  %f  %f" % (theSum, logLikeRatio,
        # self.logProposalRatio, self.logPriorRatio)
        return theSum

    def gen(self, aProposal):
        gm = ['STChain.gen()']

        # doAborts means that it was not a valid generation,
        # neither accepted or rejected.  Give up, by returning True.

        acceptMove = False

        # print "Doing %s" % aProposal.name
        pRet = self.propose(aProposal)
        #if self.tempNum == 0:
        #    print(self.propTree.postOrder)

        # print "pRet = %.6f" % pRet,
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

        # if aProposal.name == 'polytomy':
        # print "acceptMove = %s" % acceptMove
        # print "------------"
        # print " %6.0f" % pRet
        if 0 and acceptMove:
            d1 = self.propTree.topologyDistance(self.curTree, metric='scqdist')
            d2 = self.stMcmc.tree.topologyDistance(
                self.propTree, metric='scqdist')
            print(" %6.0f    %5i   %5i  %5s" % (pRet, d1, d2, acceptMove))

        aProposal.nProposals[self.tempNum] += 1
        aProposal.tnNSamples[self.tempNum] += 1
        if acceptMove:
            aProposal.accepted = True
            aProposal.nAcceptances[self.tempNum] += 1
            aProposal.tnNAccepts[self.tempNum] += 1

        # if not aProposal.doAbort:
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
            b.spaQ[0] = a.spaQ[0]

        else:
            gm.append('Unlisted proposal.name = %s  Fix me.' % aProposal.name)
            raise P4Error(gm)


# for proposal probs
fudgeFactor = {}
fudgeFactor['nni'] = 1.0
fudgeFactor['spr'] = 1.0
fudgeFactor['SR2008beta_uniform'] = 0.1
fudgeFactor['spaQ_uniform'] = 0.1
fudgeFactor['polytomy'] = 0.5


class STMcmcTunings(object):
    def __init__(self):
        self.default = {}
        self.default['SR2008beta_uniform'] = 0.2
        self.default['spaQ_uniform'] = 0.1
        

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
                raise P4Error(gm)

        else:
            self.dump()
            gm.append("    Can't set '%s'-- no such proposal." % item)
            raise P4Error(gm)

    def reprString(self):
        stuff = ["\nUser-settable relative proposal probabilities, from yourStMcmc.prob"]
        stuff.append("  To change it, do eg ")
        stuff.append("    yourSTMcmc.prob.spaQ_uniform = 0.0 # turns spaQ_uniform proposals off")
        stuff.append("  Current settings:")
        theKeys = list(self.__dict__.keys())
        theKeys.sort()
        for k in theKeys:
            stuff.append("        %20s: %s" % (k, getattr(self, k)))
        return '\n'.join(stuff)

    def dump(self):
        print(self.reprString())

    def __repr__(self):
        return self.reprString()


class STProposal(object):

    def __init__(self, theSTMcmc=None):
        self.name = None
        self.stMcmc = theSTMcmc            # reference loop
        self.nChains = theSTMcmc.nChains
        self.pNum = -1
        self.weight = 1.0
        self.tuning = None
        self.tunings = {}
        self.nProposals = [0] * self.nChains
        self.nAcceptances = [0] * self.nChains
        self.accepted = 0
        self.doAbort = False
        self.nAborts = [0] * self.nChains

        self.tnSampleSize = 250
        self.tnNSamples = [0] * theSTMcmc.nChains
        self.tnNAccepts = [0] * theSTMcmc.nChains
        self.tnAccVeryHi = None
        self.tnAccHi = None
        self.tnAccLo = None
        self.tnAccVeryLo = None
        self.tnFactorVeryHi = None
        self.tnFactorHi = None
        self.tnFactorLo = None
        self.tnFactorVeryLo = None
        self.tnFactorZero = None

    def dump(self):
        print("proposal name=%-10s pNum=%2i, weight=%5.1f, tuning=%s" % (
            self.name, self.pNum, self.weight, self.tuning))
        #print("    nProposals   by temperature:  %s" % self.nProposals)
        #print("    nAcceptances by temperature:  %s" % self.nAcceptances)

    # def _getTuning(self):
    #     if self.name in ['nni', 'spr', 'SR2008beta_uniform', 'spaQ_uniform']:
    #         # print "getting tuning for %s, returning %f" % (self.name, getattr(self.mcmc.tunings, self.name))
    #         # print self.stMcmc.tunings
    #         return getattr(self.stMcmc.tunings, self.name)
    #     else:
    #         return None

    # def _setTuning(self, whatever):
    #     raise P4Error("Can't set tuning this way.")

    # def _delTuning(self):
    #     raise P4Error("Can't del tuning.")

    # tuning = property(_getTuning, _setTuning, _delTuning)

    def tune(self, tempNum):
        assert self.tnSampleSize >= 100.
        assert self.tnNSamples[tempNum] >= self.tnSampleSize
        acc = float(self.tnNAccepts[tempNum]) / self.tnNSamples[tempNum]   # float() for Py2
        doMessage = False
        if acc > self.tnAccHi:
            oldTn = self.tuning[tempNum]
            if acc > self.tnAccVeryHi:
                self.tuning[tempNum] *= self.tnFactorVeryHi
            else:
                self.tuning[tempNum] *= self.tnFactorHi
            doMessage = True
        elif acc < self.tnAccLo:
            oldTn = self.tuning[tempNum]
            if acc < self.tnAccVeryLo:
                self.tuning[tempNum] *= self.tnFactorVeryLo
            else:
                self.tuning[tempNum] *= self.tnFactorLo
            doMessage = True
        self.tnNSamples[tempNum] = 0
        self.tnNAccepts[tempNum] = 0
        if doMessage:
            message = "%s tune  gen=%i tempNum=%i acceptance=%.3f " % (self.name, self.stMcmc.gen, tempNum, acc)
            message += "(target %.3f -- %.3f) " % (self.tnAccLo, self.tnAccHi)
            message += "Adjusting tuning from %g to %g" % (oldTn, self.tuning[tempNum])
            #print(message)
            self.stMcmc.logger.info(message)




class Proposals(object):
    def __init__(self):
        self.proposals = []
        self.proposalsDict = {}
        self.propWeights = []
        self.cumPropWeights = []
        self.totalPropWeights = 0.0
        self.intended = None

    def summary(self):
        print("There are %i proposals" % len(self.proposals))
        for p in self.proposals:
            print("proposal name=%-10s pNum=%2s, weight=%s, tuning=%s" % (
                '%s,' % p.name, p.pNum, p.weight, p.tuning))
            
    def calculateWeights(self):
        gm = ["Proposals.calculateWeights()"]
        self.propWeights = []
        for p in self.proposals:
            #print("%s: %s" % (p.name, p.weight))
            self.propWeights.append(p.weight)
        #print(self.propWeights)
        self.cumPropWeights = [self.propWeights[0]]
        for i in range(len(self.propWeights))[1:]:
            self.cumPropWeights.append(
                self.cumPropWeights[i - 1] + self.propWeights[i])
        self.totalPropWeights = sum(self.propWeights)
        if self.totalPropWeights < 1e-9:
            gm.append("No proposal weights?")
            raise P4Error(gm)
        self.intended = self.propWeights[:]
        for i in range(len(self.intended)):
            self.intended[i] /= self.totalPropWeights
        if math.fabs(sum(self.intended) - 1.0 > 1e-14):
            raise P4Error("bad sum of intended proposal probs. %s" % sum(self.intended))
        #print(self.intended)

    def chooseProposal(self, equiProbableProposals):
        if equiProbableProposals:
            return random.choice(self.proposals)
        else:
            theRan = random.uniform(0.0, self.totalPropWeights)
            for i in range(len(self.cumPropWeights)):
                if theRan < self.cumPropWeights[i]:
                    break
            return self.proposals[i]
        

    def writeProposalIntendedProbs(self):
        """Tabulate the intended proposal probabilities"""

        spacer = ' ' * 4
        print("\nIntended proposal probabilities (%)")
        print("There are %i proposals" % len(self.proposals))
        print("%2s %11s %30s %5s %12s" % ('', 'intended(%)', 'proposal', 'part', 'tuning'))
        for i in range(len(self.proposals)):
            print("%2i" % i, end=' ')
            p = self.proposals[i]
            print("   %6.2f    " % (100. * self.intended[i]), end=' ')

            print(" %27s" % p.name, end=' ')

            if p.pNum != -1:
                print(" %3i " % p.pNum, end=' ')
            else:
                print("   - ", end=' ')

            if p.tuning == None:
                print(" %12s "% '    -   ', end=' ')
            else:
                if p.tuning[0] < 0.1:
                    print(" %12.4g" % p.tuning[0], end=' ')
                elif p.tuning[0] < 1.0:
                    print(" %12.4f" % p.tuning[0], end=' ')
                elif p.tuning[0] < 10.0:
                    print(" %12.3f" % p.tuning[0], end=' ')
                elif p.tuning[0] < 1000.0:
                    print(" %12.1f" % p.tuning[0], end=' ')
                else:
                    print(" %12.2g " % p.tuning[0], end=' ')
            print()

    def writeTunings(self):
        print("Proposal tunings:")
        print("%20s %12s" % ("proposal name", "tuning"))
        for p in self.proposals:
            print("%20s" % p.name, end=' ')
            if p.tuning:
                # if p.tuning < 10.0:
                #     print("%12.3f" % p.tuning, end=' ')
                # else:
                #     print("%12.1f" % p.tuning, end=' ')
                print(p.tuning)
            else:
                print("    %4s    " % '-', end=' ')
            print()



class SwapTuner(object):
    """Continuous tuning for swap temperature"""

    def __init__(self, sampleSize):
        assert sampleSize >= 100
        self.sampleSize = sampleSize
        self.swaps01_nAttempts = 0
        self.swaps01_nSwaps = 0

        self.tnAccVeryHi = 0.18
        self.tnAccHi = 0.12
        self.tnAccLo = 0.04
        self.tnAccVeryLo = 0.01
        self.tnFactorVeryHi = 1.4
        self.tnFactorHi = 1.2
        self.tnFactorLo = 0.9
        self.tnFactorVeryLo = 0.6
        self.tnFactorZero = 0.4


    def tune(self, theMcmc):
        assert self.swaps01_nAttempts >= self.sampleSize
        acc = float(self.swaps01_nSwaps) / self.swaps01_nAttempts    # float() for Py2
        #print("SwapTuner.tune() nSwaps %i, nAttemps %i, acc %s" % (
        #    self.swaps01_nSwaps, self.swaps01_nAttempts, acc))
        doMessage = False
        direction = None
        if acc > self.tnAccHi:
            oldTn = theMcmc.chainTemp
            if acc > self.tnAccVeryHi:
                theMcmc.chainTemp *= self.tnFactorVeryHi
            else:
                theMcmc.chainTemp *= self.tnFactorHi
            doMessage = True
            direction = 'Increase'
        elif acc < self.tnAccLo:
            oldTn = theMcmc.chainTemp
            if acc == 0.0:   # no swaps at all
                theMcmc.chainTemp *= self.tnFactorZero
            elif acc < self.tnAccVeryLo:
                theMcmc.chainTemp *= self.tnFactorVeryLo
            else:
                theMcmc.chainTemp *= self.tnFactorLo
            doMessage = True
            direction = 'Decrease'
        self.swaps01_nAttempts = 0
        self.swaps01_nSwaps = 0
        if doMessage:
            message = "%s tune  gen=%i acceptance=%.3f " % ('chainTemp', theMcmc.gen, acc)
            message += "(target %.3f -- %.3f) " % (self.tnAccLo, self.tnAccHi)
            message += "%s chainTemp from %g to %g" % (direction, oldTn, theMcmc.chainTemp)
            #print(message)
            theMcmc.logger.info(message)

class STSwapTunerV(object):
    """Continuous tuning for swap temperature"""

    def __init__(self, theMcmc):
        assert var.mcmc_swapTunerSampleSize >= 100
        self.mcmc = theMcmc
        self.nChains = self.mcmc.nChains

        # These are for adjacent pairs. Eg for attempts between chains 0 and 1,
        # we increment self.nAttempts[0], ie it is indexed with the lower number
        # in the pair.
        self.nAttempts = [0] * self.nChains
        self.nSwaps = [0] * self.nChains

        self.tnAccVeryHi = 0.30
        self.tnAccHi = 0.25
        self.tnAccLo = 0.10
        self.tnAccVeryLo = 0.05
        self.tnFactorVeryHi = 1.4
        self.tnFactorHi = 1.2
        self.tnFactorLo = 0.9
        self.tnFactorVeryLo = 0.6
        self.tnFactorZero = 0.4
        self.tnLimitHi = 10.0
        self.tnLimitLo = 0.2

    def tune(self, theTempNum):
        assert self.nAttempts[theTempNum] >= var.mcmc_swapTunerSampleSize
        acc = float(self.nSwaps[theTempNum]) / self.nAttempts[theTempNum]    # float() for Py2
        # print("STSwapTunerV.tune() theTempNum %i, nSwaps %i, nAttemps %i, acc %s" % (
        #     theTempNum, self.nSwaps[theTempNum], self.nAttempts[theTempNum], acc))
        # print("tempDiffs %s" % self.mcmc.chainTempDiffs)
        # print("temps     %s" % self.mcmc.chainTemps)

        doMessage = False
        direction = None
        oldTn = self.mcmc.chainTempDiffs[theTempNum]
        if acc > self.tnAccHi:
            if self.mcmc.chainTempDiffs[theTempNum] >= self.tnLimitHi:
                direction = "no change"
            else:
                if acc > self.tnAccVeryHi:
                    self.mcmc.chainTempDiffs[theTempNum] *= self.tnFactorVeryHi
                else:
                    self.mcmc.chainTempDiffs[theTempNum] *= self.tnFactorHi
                doMessage = True
                direction = 'Increase'
        elif acc < self.tnAccLo:
            if self.mcmc.chainTempDiffs[theTempNum] <= self.tnLimitLo:
                direction = "no change"
            else:
                if acc == 0.0:   # no swaps at all
                    self.mcmc.chainTempDiffs[theTempNum] *= self.tnFactorZero
                elif acc < self.tnAccVeryLo:
                    self.mcmc.chainTempDiffs[theTempNum] *= self.tnFactorVeryLo
                else:
                    self.mcmc.chainTempDiffs[theTempNum] *= self.tnFactorLo
                doMessage = True
                direction = 'Decrease'
        self.nAttempts[theTempNum] = 0
        self.nSwaps[theTempNum] = 0
        if direction != "no change":
            if doMessage:
                message = "%s tune  gen=%i tempNum=%i acceptance=%.3f " % ('chainTemp', self.mcmc.gen, theTempNum, acc)
                message += "(target %.3f -- %.3f) " % (self.tnAccLo, self.tnAccHi)
                message += "%s chainTempDiff from %g to %g" % (direction, oldTn, self.mcmc.chainTempDiffs[theTempNum])
                #print(message)
                self.mcmc.logger.info(message)
            # Make chainTemps from chainTempDiffs
            self.mcmc.chainTemps = [0.0]
            for dNum in range(self.mcmc.nChains - 1):
                self.mcmc.chainTemps.append(self.mcmc.chainTempDiffs[dNum] + self.mcmc.chainTemps[-1])
            if doMessage:
                message = "new chainTemps gen=%i " % (self.mcmc.gen)
                for cT in self.mcmc.chainTemps:
                    message += "%10.2f" % cT
                self.mcmc.logger.info(message)


class BigTSplitStuff(object):
    # An organizer for splits on STMcmc.tree (ie bigT) internal nodes, only
    # for use with bitarray

    def __init__(self):
        self.spl = None
        self.spl2 = None
        self.theSpl = None
        self.maskedSplitWithFirstTaxOne = None
        self.onesCount = None
        self.bytes = None

    def dump(self):
        print("ss: spl=%s, spl2=%s, masked=%s, onesCount=%s" % (
            self.spl, self.spl2, self.maskedSplitWithFirstTaxOne, self.onesCount))


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
        m = func.unPickleSTMcmc(0)  # runNum 0
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

    def __init__(self, inTrees, bigT=None, modelName='SR2008_rf_aZ',
                 beta=1.0, spaQ=0.5, stRFCalc='purePython1',
                 nChains=1, runNum=0, sampleInterval=100,
                 checkPointInterval=None, useSplitSupport=False, verbose=True,
                 checkForOutputFiles=True, swapTuner=250):

        import p4.func  # This should not be needed, but it is.  Why?
        #print(p4.func)
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

        goodModelNames = ['SR2008_rf_ia', 'SR2008_rf_aZ', 'SR2008_rf_aZ_fb',
                          'SPA', 'QPA']
        if modelName not in goodModelNames:
            gm.append("Arg modelName '%s' is not recognized. " % modelName)
            gm.append("Good modelNames are %s" % goodModelNames)
            raise P4Error(gm)
        self.modelName = modelName
        self.tree = None

        self.stRFCalc = None
        if modelName.startswith("SR2008"):
            try:
                fBeta = float(beta)
            except ValueError:
                gm.append("Arg beta (%s) should be a float" % beta)
                raise P4Error(gm)
            self.beta = fBeta

            for t in inTrees:
                if t.isFullyBifurcating():
                    pass
                else:
                    gm.append("The SR2008 model wants trees that are fully bifurcating.")
                    raise P4Error(gm)

            goodSTRFCalcNames = ['purePython1', 'bitarray', 'fastReducedRF']
            if stRFCalc not in goodSTRFCalcNames:
                gm.append("Arg stRFCalc '%s' is not recognized. " % modelName)
                gm.append("Good stRFCalc names are %s" % goodSTRFCalcNames)
                raise P4Error(gm)
            self.stRFCalc = stRFCalc


        try:
            nChains = int(nChains)
        except (ValueError, TypeError):
            gm.append("nChains should be an int, 1 or more.  Got %s" % nChains)
            raise P4Error(gm)
        if nChains < 1:
            gm.append("nChains should be an int, 1 or more.  Got %s" % nChains)
            raise P4Error(gm)
        self.nChains = nChains
        self.chains = []
        self.gen = -1
        self.startMinusOne = -1
        self.chainTemp = 1.0

        self.constraints = None
        self.simulate = None


        # spaQ is a property.  Whenever it is set, it is propagated to all the chains.
        self._spaQ = None
        if modelName in ['SPA', 'QPA']:
            try:
                self._spaQ = float(spaQ)
            except ValueError:
                gm.append("Arg spaQ (%s) should be a float" % spaQ)
                raise P4Error(gm)
            self.spaQ = self._spaQ

        try:
            runNum = int(runNum)
        except (ValueError, TypeError):
            gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
            raise P4Error(gm)
        if runNum < 0:
            gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
            raise P4Error(gm)
        self.runNum = runNum

        self._setLogger()

        if checkForOutputFiles:
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
                raise P4Error(gm)

            if var.strictRunNumberChecking:
                # We want to start runs with number 0, so if runNum is more than
                # that, check that there are other runs.
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
                            gm.append("There are no mcmc_trees_%i.nex files to show that run %i has been done." % (
                                runNum2, runNum2))
                            gm.append("Set the runNum to that, first.")
                            gm.append("Or else turn var.strictRunNumberChecking off to prevent checking.")
                            raise P4Error(gm)

        self.sampleInterval = sampleInterval
        self.checkPointInterval = checkPointInterval

        self.props = Proposals()
        self.tunableProps = """SR2008beta_uniform spaQ_uniform""".split()
        # maybeTunableButNotNow  polytomy


        self.treePartitions = None
        self.likesFileName = "mcmc_likes_%i" % runNum
        self.treeFileName = "mcmc_trees_%i.nex" % runNum
        #self.simFileName = "mcmc_sims_%i" % runNum
        self.pramsFileName = "mcmc_prams_%i" % runNum
        self.writePrams = False
        if self.modelName in ['SR2008_rf_aZ_fb', "SPA", "QPA"]:
            self.writePrams = True

        self.lastTimeCheck = None

        if self.nChains > 1:
            self.swapMatrix = []
            for i in range(self.nChains):
                self.swapMatrix.append([0] * self.nChains)
            # if self.swapVector:
            #     self.swapTuner = STSwapTunerV(self)
            # else:
            #     if swapTuner:             # a kwarg
            #         myST = int(swapTuner)
            #         if myST >= 100:
            #             self.swapTuner = SwapTuner(myST)
            #         else:
            #             gm.append("The swapTuner kwarg, the sample size, should be at least 100.  Got %i." % myST)
            #             raise P4Error(gm)
            #     else:
            #         self.swapTuner = None

        else:
            self.swapMatrix = None
        self.swapTuner = None

        self._tunings = STMcmcTunings()
        self.polytomyUseResolutionClassPrior = False
        self.polytomyPriorLogBigC = 0.0

        self.prob = STMcmcProposalProbs()
        if self.modelName in ['SPA', 'QPA']:
            self.prob.polytomy = 1.0
            self.prob.spr = 0.0

        # Zap internal node names
        # for n in aTree.root.iterInternals():
        #     if n.name:
        #         n.name = None

        if not bigT:
            allNames = set()
            for t in inTrees:
                t.unsorted_taxNames = [n.name for n in t.iterLeavesNoRoot()]
                # Get the union of a set and other stuff using set.update(stuff).
                allNames.update(t.unsorted_taxNames)
            self.taxNames = list(allNames)
            # not needed, but nice for debugging
            self.taxNames.sort()
        else:
            for t in inTrees:
                t.unsorted_taxNames = [n.name for n in t.iterLeavesNoRoot()]
            self.taxNames = bigT.taxNames
        # print self.taxNames
        self.nTax = len(self.taxNames)

        if self.modelName in ['SPA'] or self.stRFCalc == 'bitarray':
            # print "self.taxNames = ", self.taxNames
            for t in inTrees:
                # print "-" * 50
                # t.draw()
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
                # print "intree baTaxBits is %s" % t.baTaxBits
                # print "intree firstTax is %i" % t.firstTax

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
                        # print "setting node %i stSplitKey to %s" %
                        # (n.nodeNum, n.stSplitKey)
                if self.stRFCalc == 'bitarray':
                    t.splSet = set()
                    for n in t.iterInternalsNoRoot():
                        # make sure splitKey[firstTax] is a '1'
                        if not n.stSplitKey[t.firstTax]:
                            n.stSplitKey.invert()
                            n.stSplitKey &= t.baTaxBits     # 'and', in-place
                            # print "inverting and and-ing node %i stSplitKey
                            # to %s" % (n.nodeNum, n.stSplitKey)
                        # bytes so that I can use it as a set element
                        t.splSet.add(n.stSplitKey.tobytes())
                if self.modelName in ['SPA']:
                    t.internals = []
                    for n in t.iterInternalsNoRoot():
                        # make sure splitKey[firstTax] is a '1'
                        if not n.stSplitKey[t.firstTax]:
                            n.stSplitKey.invert()
                            n.stSplitKey &= t.baTaxBits     # 'and', in-place
                            # print "inverting and and-ing node %i stSplitKey
                            # to %s" % (n.nodeNum, n.stSplitKey)
                        # bytes so that I can use it as a set element
                        n.stSplitKeyBytes = n.stSplitKey.tobytes()
                        t.internals.append(n)

        self.fspa = None
        if self.modelName == 'SPA' and var.stmcmc_useFastSpa:
            import p4.fastspa as fastspa
            self.fspa = fastspa.FastSpa(useSplitSupport)
            for tNum, t in enumerate(inTrees):
                self.fspa.setInTr(tNum, t.nTax, self.nTax, t.baTaxBits.to01(), t.firstTax)
                for n in t.internals:
                    if n.br and hasattr(n.br, "support"):
                        support = n.br.support
                    else:
                        support = -1.0
                    self.fspa.setInTrNo(tNum, n.stSplitKey.to01(), support)
            #self.fspa.summarizeInTrs()

        if self.modelName in ['QPA']:
            for t in inTrees:
                sorted_taxNames = []
                t.taxBits = []
                for tNum in range(self.nTax):
                    tN = self.taxNames[tNum]
                    if tN in t.unsorted_taxNames:
                        sorted_taxNames.append(tN)
                        t.taxBits.append(1 << tNum)
                    else:
                        t.taxBits.append(0)
                t.taxNames = sorted_taxNames
                # print "intree taxBits is %s" % t.taxBits

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
                        n.stSplitKey = 1 << spot
                        #self.tBits[spot] = False
                    else:
                        n.stSplitKey = n.leftChild.stSplitKey
                        p = n.leftChild.sibling
                        while p:
                            n.stSplitKey |= p.stSplitKey    # "or", in-place
                            p = p.sibling
                        # print "setting node %i stSplitKey to %s" % (n.nodeNum, n.stSplitKey)
                # t.splSet = set()
                # for n in t.iterInternalsNoRoot():
                #     if not n.stSplitKey[t.firstTax]:   # make sure splitKey[firstTax] is a '1'
                #         n.stSplitKey.invert()
                #         n.stSplitKey &= t.baTaxBits     # 'and', in-place
                #         #print "inverting and and-ing node %i stSplitKey to %s" % (n.nodeNum, n.stSplitKey)
                # t.splSet.add(n.stSplitKey.tobytes()) # bytes so that I can
                # use it as a set element
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
                                t.qSet.add(down + up)
                            else:
                                t.qSet.add(up + down)
                # print t.qSet
                t.nQuartets = len(t.qSet)

        self.trees = inTrees
        if bigT:
            self.tree = bigT
        else:
            self.tree = p4.func.randomTree(taxNames=self.taxNames, name='stTree', randomBrLens=False)

        if self.stRFCalc in ['purePython1', 'fastReducedRF']:
            for t in inTrees:
                sorted_taxNames = []
                t.taxBits = 0
                for tNum in range(self.nTax):
                    tN = self.taxNames[tNum]
                    if tN in t.unsorted_taxNames:
                        sorted_taxNames.append(tN)
                        adder = 1 << tNum
                        t.taxBits += adder
                t.taxNames = sorted_taxNames
                t.allOnes = 2 ** (t.nTax) - 1
                t.makeSplitKeys()
                t.skSet = set([n.br.splitKey for n in t.iterInternalsNoRoot()])

        if self.stRFCalc in ['purePython1', 'fastReducedRF']:
            self.tree.makeSplitKeys()

            self.Frrf = None
            if self.stRFCalc == 'fastReducedRF':
                try:
                    import p4.fastReducedRF
                    self.Frrf = p4.fastReducedRF.Frrf
                    # not explicitly used--but makes converters available
                    import pyublas
                except ImportError:
                    gm.append("var.stRFCalc is set to 'fastReducedRF', but I could not import")
                    gm.append("at least one of fastReducedRF or pyublas.")
                    gm.append("Make sure they are installed.")
                    raise P4Error(gm)

        if self.modelName in ['QPA']:
            t = self.tree
            t.taxBits = [1 << i for i in range(t.nTax)]
            for n in t.iterPostOrder():
                if n == t.root:
                    break
                if n.isLeaf:
                    spot = self.taxNames.index(n.name)
                    n.stSplitKey = 1 << spot
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
                            t.qSet.add(down + up)
                        else:
                            t.qSet.add(up + down)
            # print t.qSet
            t.nQuartets = len(t.qSet)

        self.useSplitSupport = False
        if useSplitSupport:
            if self.modelName.startswith("SR2008"):
                gm.append(
                    "Arg useSplitSupport is turned on, but it is not implemented with SR2008")
                raise P4Error(gm)
            assert useSplitSupport in [True, 'percent']
            self.useSplitSupport = True
            hasSplitInfo = False
            for it in self.trees:
                for n in it.iterInternalsNoRoot():
                    if not hasattr(n.br, 'support'):
                        n.br.support = None
                    else:
                        assert n.br.support == None
                    if n.name:
                        flName = float(n.name)
                        hasSplitInfo = True
                        if useSplitSupport == 'percent':
                            flName *= 0.01
                        if flName < 0.0 or flName > 1.0:
                            gm.append("Input tree %s" %
                                      it.writeNewick(toString=True))
                            gm.append(
                                "Got support value %s, outside of range 0 to 1" % n.name)
                            if flName > 1.0 and useSplitSupport == True:
                                gm.append(
                                    "Maybe it is percent support?  If so, set useSplitSupport to 'percent' rather than True")
                            raise P4Error(gm)
                        n.br.support = flName
                        #n.br.logSupport = math.log(n.br.support)
                    else:
                        #n.br.logSupport = None
                        pass
            if not hasSplitInfo:
                gm.append(
                    "Arg useSplitSupport is turned on, but none of the trees seem to have split info.")
                raise P4Error(gm)


        splash = p4.func.splash2(verbose=False)
        for aLine in splash:
            self.logger.info(aLine)

        self.swapVector = True
        if self.nChains > 1:
            self.swapTuner = STSwapTunerV(self)


        # Hidden experimental hacking
        self.doHeatingHack = False
        self.heatingHackTemperature = 5.0
        #self.heatingHackProposalNames = ['nni', 'spr']

        if verbose:
            self.loggerPrinter.info("Initializing STMcmc")
            self.loggerPrinter.info("%-16s: %s" % ('modelName', modelName))
            if self.modelName.startswith("SR2008"):
                self.loggerPrinter.info("%-16s: %s" % ('stRFCalc', self.stRFCalc))
            if self.modelName in ["SPA", "QPA"]:
                self.loggerPrinter.info("%-16s: %s" % ('useSplitSupport', self.useSplitSupport))
            self.loggerPrinter.info("%-16s: %s" % ('inTrees', len(self.trees)))
            self.loggerPrinter.info("%-16s: %s" % ('nTax', self.nTax))
            if self.nChains == 1:
                self.loggerPrinter.info("%-16s: %s" % ('mcmcmc', "off: 1 chain"))
            elif self.nChains > 1:
                self.loggerPrinter.info("%-16s: %s" % ('mcmcmc', "on -- %i chains" % self.nChains))
                if self.swapVector:
                    self.loggerPrinter.info("%-16s: %s" % ('swapVector', "on"))
                    self.loggerPrinter.info("%-16s: %s" % ('swapTuner', "on"))

    def _del_nothing(self):
        gm = ["Don't/Can't delete this property."]
        raise P4Error(gm)

    def _get_spaQ(self):
        return self._spaQ

    def _set_spaQ(self, newVal):
        try:
            newVal = float(newVal)
        except:
            gm = ['This property should be set to a float.']
            raise P4Error(gm)
        self._spaQ = newVal
        if self.chains:
            for ch in self.chains:
                ch.propTree.spaQ = newVal
                ch.curTree.spaQ = newVal

    spaQ = property(_get_spaQ, _set_spaQ, _del_nothing)
    """(property) The current spaQ"""

    def _setLogger(self):
        """Make two loggers; one that writes to a file and to stderr, and one that writes only to a file."""

        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s %(message)s',
                            datefmt='[%Y-%m-%d %H:%M]',
                            filename="mcmc_log_%i" % self.runNum,
                            filemode='a')

        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(message)s')
        # tell the handler to use this format
        console.setFormatter(formatter)
        # add the handler to the root logger
        # Using named loggers allows me to keep them separate.
        self.loggerPrinter = logging.getLogger('withPrint')
        self.loggerPrinter.addHandler(console)
        # This logger only logs to the file, not to stderr.
        self.logger = logging.getLogger("logFileOnly")



    def _makeProposals(self):
        """Make proposals for the STMcmc."""

        gm = ['STMcmc._makeProposals()']

        # nni
        if self.prob.nni:
            p = STProposal(self)
            p.name = 'nni'
            # * (len(self.tree.nodes) - 1) * fudgeFactor['nni']
            #print(self.prob.nni)
            #print(fudgeFactor)
            
            p.weight = self.prob.nni * fudgeFactor['nni']
            self.props.proposals.append(p)

        if self.prob.spr:
            p = STProposal(self)
            p.name = 'spr'
            # * (len(self.tree.nodes) - 1) * fudgeFactor['spr']
            p.weight = self.prob.spr * fudgeFactor['spr']
            self.props.proposals.append(p)

        if self.modelName in ['SR2008_rf_aZ_fb']:
            if self.prob.SR2008beta_uniform:
                p = STProposal(self)
                p.name = 'SR2008beta_uniform'
                p.tuning = [self._tunings.default[p.name]] * self.nChains
                # * (len(self.tree.nodes) - 1) * fudgeFactor['SR2008beta_uniform']
                p.weight = self.prob.SR2008beta_uniform * fudgeFactor['SR2008beta_uniform']

                p.tnAccVeryHi = 0.7
                p.tnAccHi = 0.6
                p.tnAccLo = 0.05
                p.tnAccVeryLo = 0.03

                p.tnFactorVeryHi = 1.6
                p.tnFactorHi = 1.2
                p.tnFactorLo = 0.8
                p.tnFactorVeryLo = 0.7

                self.props.proposals.append(p)

        if self.modelName in ['SPA', 'QPA']:
            if self.prob.spaQ_uniform:
                p = STProposal(self)
                p.name = 'spaQ_uniform'
                p.tuning = [self._tunings.default[p.name]] * self.nChains
                p.spaQPriorType = 'flat'
                p.spaQExpPriorLambda = 100.
                # * (len(self.tree.nodes) - 1) * fudgeFactor['spaQ_uniform']
                p.weight = self.prob.spaQ_uniform * fudgeFactor['spaQ_uniform']
  
                p.tnAccVeryHi = 0.7
                p.tnAccHi = 0.6
                p.tnAccLo = 0.05
                p.tnAccVeryLo = 0.03

                p.tnFactorVeryHi = 1.6
                p.tnFactorHi = 1.2
                p.tnFactorLo = 0.8
                p.tnFactorVeryLo = 0.7

                self.props.proposals.append(p)

            if self.prob.polytomy:
                p = STProposal(self)
                p.name = 'polytomy'
                p.polytomyUseResolutionClassPrior = self.polytomyUseResolutionClassPrior
                p.polytomyPriorLogBigC = self.polytomyPriorLogBigC
                p.weight = self.prob.polytomy * fudgeFactor['polytomy']
                self.props.proposals.append(p)

        if not self.props.proposals:
            gm.append("No proposals?")
            raise P4Error(gm)
        for p in self.props.proposals:
            self.props.proposalsDict[p.name] = p
        self.props.calculateWeights()

    def _refreshProposalProbsAndTunings(self):
        """Adjust proposals after a restart."""

        gm = ['STMcmc._refreshProposalProbsAndTunings()']

        for p in self.props.proposals:
            # nni
            if p.name == 'nni':
                #p.weight = self.prob.local * (len(self.tree.nodes) - 1) * fudgeFactor['local']
                p.weight = self.prob.nni

        self.propWeights = []
        for p in self.props.proposals:
            self.propWeights.append(p.weight)
        self.cumPropWeights = [self.propWeights[0]]
        for i in range(len(self.propWeights))[1:]:
            self.cumPropWeights.append(
                self.cumPropWeights[i - 1] + self.propWeights[i])
        self.totalPropWeights = sum(self.propWeights)
        if self.totalPropWeights < 1e-9:
            gm.append("No proposal weights?")
            raise P4Error(gm)

    def writeProposalAcceptances(self):
        """Pretty-print the proposal acceptances."""

        if (self.gen - self.startMinusOne) <= 0:
            print("\nSTMcmc.writeProposalAcceptances()  There is no info in memory. ")
            print(" Maybe it was just emptied after writing to a checkpoint?  ")
            print("If so, read the checkPoint and get the proposalAcceptances from there.")
            return

        spacer = ' ' * 8
        print("\nProposal acceptances, run %i, for %i gens, from gens %i to %i, inclusive." % (
            self.runNum, (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen))
        print("%s %20s %10s %13s%8s" % (spacer, 'proposal', 'nProposals', 'acceptance(%)', 'tuning'))
        for p in self.props.proposals:
            print("%s" % spacer, end=' ')
            print("%20s" % p.name, end=' ')
            print("%10i" % p.nProposals[0], end=' ')

            if p.nProposals[0]:  # Don't divide by zero
                print("       %5.1f " % (100.0 * float(p.nAcceptances[0]) / float(p.nProposals[0])), end=' ')
            else:
                print("           - ", end=' ')

            if p.tuning == None:
                print("      -", end=' ')
            elif p.tuning[0] < 2.0:
                print("  %8.4f" % p.tuning[0], end=' ')
            elif p.tuning[0] < 20.0:
                print("  %8.3f" % p.tuning[0], end=' ')
            elif p.tuning[0] < 200.0:
                print("  %8.1f" % p.tuning[0], end=' ')
            else:
                print("  %8.3g" % p.tuning[0], end=' ')
            print()


        # # Tabulate topology changes by temperature
        if self.nChains > 1:
            for propName in ['nni', 'spr']:
                p = self.props.proposalsDict.get(propName)
                if p:
                    print("'%s' proposal-- topology changes by temperature" % (propName))
                    print("%s tempNum   nProps nAccepts percent" % spacer)
                    for tNum in range(self.nChains):
                        print("%s" % spacer, end=' ')
                        print("%4i " % tNum, end=' ')
                        print("%9i" % p.nProposals[tNum], end=' ')
                        print("%8i" % p.nAcceptances[tNum], end=' ')
                        print("  %5.1f" % (100.0 * float(p.nAcceptances[tNum]) / float(p.nProposals[tNum])))

        #     # Check for aborts.
        #     p = self.proposalsHash.get(propName)
        #     if p:
        #         if hasattr(p, 'nAborts'):
        #             if p.nAborts[0]:
        #                 print("The '%s' proposal had %i aborts in the cold chain." % (propName, p.nAborts[0]))
        #                 if self.constraints:
        #                     print("(Aborts might be due to violated constraints.)")
        #             else:
        #                 print("The '%s' proposal had no aborts in the cold chain" % propName)
        for pN in ['polytomy']:
            p = None
            try:
                p = self.props.proposalsDict[pN]
            except KeyError:
                pass
            if p:
                if hasattr(p, 'nAborts'):
                    print("The %s proposal had %5i aborts." % (p.name, p.nAborts[0]))

        if self.nChains > 1:
            print("\n\nAcceptances and tunings by temperature")
            print("%s %30s %5s %5s %10s %13s%10s" % (
                spacer, 'proposal', 'part', 'tempNum', 'nProposals', 'acceptance(%)', 'tuning'))
            for p in self.props.proposals:
                for tempNum in range(self.nChains):
                    print("%s" % spacer, end=' ')
                    print("%30s" % p.name, end=' ')
                    if p.pNum != -1:
                        print(" %3i " % p.pNum, end=' ')
                    else:
                        print("   - ", end=' ')
                    print(" %3i " % tempNum, end=' ')
                    print("%10i" % p.nProposals[tempNum], end=' ')

                    if p.nProposals[tempNum]:  # Don't divide by zero
                        print("       %5.1f " % (
                            100.0 * float(p.nAcceptances[tempNum]) / float(p.nProposals[tempNum])), end=' ')
                    else:
                        print("           - ", end=' ')

                    if p.tuning == None:
                        print("      -", end=' ')
                    elif p.tuning[tempNum] < 2.0:
                        print("  %8.4f" % p.tuning[tempNum], end=' ')
                    elif p.tuning[tempNum] < 20.0:
                        print("  %8.3f" % p.tuning[tempNum], end=' ')
                    elif p.tuning[tempNum] < 200.0:
                        print("  %8.1f" % p.tuning[tempNum], end=' ')
                    else:
                        print("  %8.3g" % p.tuning[tempNum], end=' ')
                    print()


    def writeSwapMatrix(self):
        print("\nChain swapping, for %i gens, from gens %i to %i, inclusive." % (
            (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen))
        #print("    Swaps are presented as a square matrix, nChains * nChains.")
        print("    Upper triangle is the number of swaps proposed between two chains.")
        print("    Lower triangle is the percent swaps accepted.")
        #print("    The current tunings.chainTemp is %5.3f\n" % self.chainTemp)
        #if var.mcmc_swapVector:
        #    print("    The chainTemp is continuously tuned for each chain\n")
        #else:
        #print("    The chainTemp is %f.\n" % self.chainTemp)

        print(" " * 10, end=' ')
        for i in range(self.nChains):
            print("%7i" % i, end=' ')
        print()
        print(" " * 10, end=' ')
        for i in range(self.nChains):
            print("   ----", end=' ')
        print()
        for i in range(self.nChains):
            print(" " * 7, "%2i" % i, end=' ')
            for j in range(self.nChains):
                if i < j:  # upper triangle
                    print("%7i" % self.swapMatrix[i][j], end=' ')
                elif i == j:
                    print("      -", end=' ')
                else:
                    if self.swapMatrix[j][i] == 0:  # no proposals
                        print("      -", end=' ')
                    else:
                        print("  %5.1f" % (100.0 * float(self.swapMatrix[i][j]) / float(self.swapMatrix[j][i])), end=' ')
            print()

    def _makeChainsAndProposals(self):
        """Make chains and proposals."""

        gm = ['STMcmc._makeChainsAndProposals()']

        # random.seed(0)

        # Make chains, if needed
        if not self.chains:
            self.chains = []
            # chNum is used by fastspa; it is also a starting point for tempNum (which changes in swaps)
            for chNum in range(self.nChains):
                aChain = STChain(self, chNum)
                self.chains.append(aChain)
        if not self.props.proposals:
            self._makeProposals()

            # If we are going to be doing the resolution class prior
            # in the polytomy move, we want to pre-compute the logs of
            # T_{n,m}.  Its a vector with indices (ie m) from zero to
            # nTax-2 inclusive.
            p = self.props.proposalsDict.get('polytomy')
            if p and self.polytomyUseResolutionClassPrior:
                bigT = p4.func.nUnrootedTreesWithMultifurcations(self.tree.nTax)
                p.logBigT = [0.0] * (self.tree.nTax - 1)
                for i in range(1, self.tree.nTax - 1):
                    p.logBigT[i] = math.log(bigT[i])
                #print p.logBigT

    def _setOutputTreeFile(self):
        """Setup the (output) tree file for the STMcmc."""

        gm = ['STMcmc._setOutputTreeFile()']

        # Write the preamble for the trees outfile.
        treeFile = open(self.treeFileName, 'w')
        treeFile.write('#nexus\n\n')
        treeFile.write('begin taxa;\n')
        treeFile.write('  dimensions ntax=%s;\n' % self.tree.nTax)
        treeFile.write('  taxlabels')
        for tN in self.tree.taxNames:
            treeFile.write(' %s' % p4.func.nexusFixNameIfQuotesAreNeeded(tN))
        treeFile.write(';\nend;\n\n')

        treeFile.write('begin trees;\n')
        self.translationHash = {}
        i = 1
        for tName in self.tree.taxNames:
            self.translationHash[tName] = i
            i += 1

        treeFile.write('  translate\n')
        for i in range(self.tree.nTax - 1):
            treeFile.write('    %3i %s,\n' % (
                i + 1, p4.func.nexusFixNameIfQuotesAreNeeded(self.tree.taxNames[i])))
        treeFile.write('    %3i %s\n' % (
            self.tree.nTax, p4.func.nexusFixNameIfQuotesAreNeeded(self.tree.taxNames[-1])))
        treeFile.write('  ;\n')
        treeFile.write('  [Tree numbers are gen+1]\n')
        treeFile.close()

    def run(self, nGensToDo, verbose=True, equiProbableProposals=False, writeSamples=True):
        """Start the STMcmc running."""

        gm = ['STMcmc.run()']

        # Hidden experimental hack
        if self.doHeatingHack:
            print("Heating hack is turned on.")
            assert self.nChains == 1, "MCMCMC does not work with the heating hack"
            print("Heating hack temperature is %.2f" % self.heatingHackTemperature)
            #print("Heating hack affects proposals %s" % self.heatingHackProposalNames)

        # Keep track of the first gen of this call to run(), maybe restart
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
            # if ((self.gen + 1) + nGensToDo) % self.checkPointInterval == 0:
            if nGensToDo % self.checkPointInterval == 0:
                pass
            else:
                gm.append(
                    "With the current settings, the last generation won't be on a checkPointInterval.")
                gm.append("self.gen+1=%i, nGensToDo=%i, checkPointInterval=%i" % ((self.gen + 1),
                                                                                  nGensToDo, self.checkPointInterval))
                raise P4Error(gm)
            #  2.  We also want the checkPointInterval to be evenly
            #      divisible by the sampleInterval.
            if self.checkPointInterval % self.sampleInterval == 0:
                pass
            else:
                gm.append(
                    "The checkPointInterval (%i) should be evenly divisible" % self.checkPointInterval)
                gm.append("by the sampleInterval (%i)." % self.sampleInterval)
                raise P4Error(gm)


        if self.props.proposals:
            # Its either a re-start, or it has been thru autoTune().
            # I can tell the difference by self.gen, which is -1 after
            # autoTune()
            if self.gen == -1:
                self._makeChainsAndProposals()
                self._setOutputTreeFile()
                # if self.simulate:
                #    self.writeSimFileHeader(self.tree)
            # The probs and tunings may have been changed by the user.
            self._refreshProposalProbsAndTunings()

            # This stuff below should be the same as is done after pickling,
            # see below.
            self.startMinusOne = self.gen

            # Start the tree partitions over.
            self.treePartitions = None
            # Zero the proposal counts
            for p in self.props.proposals:
                p.nProposals = [0] * self.nChains
                p.nAcceptances = [0] * self.nChains
                #p.nTopologyChangeAttempts = [0] * self.nChains
                #p.nTopologyChanges = [0] * self.nChains
            # Zero the swap matrix
            if self.nChains > 1:
                self.swapMatrix = []
                for i in range(self.nChains):
                    self.swapMatrix.append([0] * self.nChains)

        else:
            self._makeChainsAndProposals()
            self._setOutputTreeFile()
            # if self.simulate:
            #    self.writeSimFileHeader(self.tree)

            # The swap vector is just the diagonal of the swap matrix
            if self.swapVector and self.nChains > 1:
                # These are differences in temperatures between adjacent chains.  The last one is not used.
                self.chainTempDiffs = [self.chainTemp] * self.nChains 
                # These are cumulative, summed over the diffs.  This needs to be done whenever the diffs change
                self.chainTemps = [0.0]
                for dNum in range(self.nChains - 1):
                    self.chainTemps.append(self.chainTempDiffs[dNum] + self.chainTemps[-1])



        if verbose:
            self.props.writeProposalIntendedProbs()
            sys.stdout.flush()

        coldChainNum = 0            

        # If polytomy is turned on, then it is possible to get a star
        # tree, in which case local will not work.  So if we have both
        # polytomy and local proposals, we should also have brLen.
        # if "polytomy" in self.proposalsHash and 'local' in self.proposalsHash:
        #     if 'brLen' not in self.proposalsHash:
        #         gm.append("If you have polytomy and local proposals, you should have a brLen proposal as well.")
        #         gm.append("It can have a low proposal probability, but it needs to be there.")
        #         gm.append("Turn it on by eg yourMcmc.prob.brLen = 0.001")
        #         raise P4Error(gm)

        if self.gen > -1:
            # it is a re-start, so we need to back over the "end;" in the tree
            # files.
            f2 = open(self.treeFileName, 'r+b')
            pos = -1
            while 1:
                f2.seek(pos, 2)
                c = f2.read(1)
                if c == b';':
                    break
                pos -= 1
            # print "pos now %i" % pos
            pos -= 3  # end;
            f2.seek(pos, 2)
            c = f2.read(4)
            # print "got c = '%s'" % c
            if c != b"end;":
                gm.append("Stmcmc.run().  Failed to find and remove the 'end;' at the end of the tree file.")
                raise P4Error(gm)
            else:
                f2.seek(pos, 2)
                f2.truncate()
            f2.close()

            self.logger.info("Re-starting the ST MCMC run %i from gen=%i" % (self.runNum, self.gen))

            if verbose:
                print()
                print("Re-starting the ST MCMC run %i from gen=%i" % (self.runNum, self.gen))
                if not writeSamples:
                    print("Arg 'writeSamples' is off" )
                print("Set to do %i more generations." % nGensToDo)
                # if self.writePrams:
                #    if self.chains[0].curTree.model.nFreePrams == 0:
                #        print "There are no free prams in the model, so I am turning writePrams off."
                #        self.writePrams = False
                sys.stdout.flush()

            self.startMinusOne = self.gen
        else:
            self.logger.info("Starting the ST MCMC %s run %i" % ((self.constraints and "(with constraints)" or ""), self.runNum))
            self.logger.info("nChains %i" % (self.nChains))

            ret = self.props.proposalsDict.get('polytomy')
            if ret:
                message = "Doing polytomy proposal, with polytomyUseResolutionClassPrior=%s" % ret.polytomyUseResolutionClassPrior
                self.loggerPrinter.info(message)
                message = "polytomy: polytomyPriorLogBigC=%f" % ret.polytomyPriorLogBigC
                self.loggerPrinter.info(message)

            if verbose:
                if self.nChains > 1:
                    print("Using Metropolis-coupled MCMC, with %i chains.  Temperature %f." % (self.nChains, self.chainTemp))
                else:
                    print("Not using Metropolis-coupled MCMC.")
                self.loggerPrinter.info("Starting the ST MCMC %s run %i" % ((self.constraints and "(with constraints)" or ""), self.runNum))
                self.loggerPrinter.info("Set to do %i generations." % nGensToDo)
                if self.writePrams:
                    # if self.chains[0].curTree.model.nFreePrams == 0:
                    #     print "There are no free prams in the model, so I am turning writePrams off."
                    #     self.writePrams = False
                    # else:
                    pramsFile = open(self.pramsFileName, 'a')
                    if self.modelName.startswith("SR2008"):
                        pramsFile.write("    genPlus1     beta\n")
                    elif self.modelName.startswith("SPA"):
                        pramsFile.write("    genPlus1     spaQ\n")
                    elif self.modelName.startswith("QPA"):
                        pramsFile.write("    genPlus1     spaQ\n")
                    pramsFile.close()
                sys.stdout.flush()

        if not writeSamples:
            self.logger.info("STMcmc.run() arg 'writeSamples' is off, so samples are not being written")
        if equiProbableProposals:
            self.logger.info("STMcmc.run() arg 'equiProbableProposals' is turned on")



        if verbose:
            if writeSamples:
                print("Sampling every %i." % self.sampleInterval)
            else:
                print("Arg 'writeSamples' is off, so samples are not being written")
            if equiProbableProposals:
                print("Arg 'equiProbableProposals' is turned on")
            if self.checkPointInterval:
                print("CheckPoints written every %i." % self.checkPointInterval)
            if nGensToDo <= 20000:
                print("One dot is 100 generations.")
            else:
                print("One dot is 1000 generations.")
            sys.stdout.flush()

        self.treePartitions = None
        realTimeStart = time.time()
        self.lastTimeCheck = time.time()

        abortableProposals = ['nni', 'spr', 'polytomy']

        ##################################################
        ############### Main loop ########################
        ##################################################

        self.swapInterval = 2

        for gNum in range(nGensToDo):
            self.gen += 1

            # Do an initial time estimate based on 100 gens
            if nGensToDo > 100 and self.gen - firstGen == 100:
                diff_secs = time.time() - realTimeStart
                total_secs = (float(nGensToDo) / float(100)) * float(diff_secs)
                deltaTime = datetime.timedelta(seconds=int(round(total_secs)))
                print("Estimated completion time: %s days, %s" % (
                    deltaTime.days, time.strftime("%H:%M:%S", time.gmtime(deltaTime.seconds))))

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
                        # equiProbableProposals is True or False.  Usually False.
                        aProposal = self.props.chooseProposal(equiProbableProposals)
                        if aProposal:
                            gotIt = True

                        if aProposal.name == 'nni':
                            # Can't do nni on a star tree.
                            if self.chains[chNum].curTree.nInternalNodes == 1:
                                #aProposal = self.props.proposalsDict['polytomy']
                                gotIt = False

                        if aProposal.doAbort:
                            gotIt = False

                        safety += 1
                        if safety > 1000:
                            gm.append(
                                "Could not find a proposal after %i attempts." % safety)
                            gm.append("Possibly a programming error.")
                            gm.append(
                                "Or possibly it is just a pathologically frustrating Mcmc.")
                            raise P4Error(gm)

                    if 0:
                        print("==== gNum=%i, chNum=%i, aProposal=%s" % (
                            gNum, chNum, aProposal.name), end=' ')
                        sys.stdout.flush()
                        # print gNum,

                    # success returns None
                    failure = self.chains[chNum].gen(aProposal)

                    if failure:
                        myWarn = "STMcmc.run() main loop.  Proposal %s generated a 'failure'.  Why?" % aProposal.name
                        self.logger.warning(myWarn)

                    if 0:
                        if failure:
                            print("    failure")
                        else:
                            print()

                    nAttempts += 1
                    if nAttempts > 1000:
                        gm.append("Was not able to do a successful generation after %i attempts." % nAttempts)
                        raise P4Error(gm)


                    # Continuous tuning.  We have a tuning, and propose/accept
                    # tallies for each temperature, kept separately.  Note that
                    # since chNum does not equal tempNum, the most recently
                    # incremented values (and the ones we want to tune now) will
                    # be aProposal.tnNSamples[tempNum] and
                    # aProposal.tnNAccepts[tempNum], where tempNum will most
                    # likely not be chNum.  So we get the tempNum from this
                    # chNum, and tune it.

                    if aProposal.name in self.tunableProps:
                        tempNum = self.chains[chNum].tempNum
                        if aProposal.tnNSamples[tempNum] >= aProposal.tnSampleSize:
                            aProposal.tune(tempNum)



                # print "   Mcmc.run(). finished a gen on chain %i" % (chNum)
                for prNm in abortableProposals:
                    ret = self.props.proposalsDict.get(prNm)
                    if ret:
                        ret.doAbort = False

            # Do swap, if there is more than 1 chain.
            if (self.gen + 1) % self.swapInterval == 0:
                if self.nChains == 1:
                    coldChain = 0
                else:
                    if self.swapVector:
                        rTempNum1 = random.randrange(self.nChains - 1)
                        rTempNum2 = rTempNum1 + 1
                        chain1 = None
                        chain2 = None
                        for ch in self.chains:
                            if ch.tempNum == rTempNum1:
                                chain1 = ch
                            elif ch.tempNum == rTempNum2:
                                chain2 = ch
                        assert chain1 and chain2

                        # Use the upper triangle of swapMatrix for nAttempts
                        self.swapMatrix[chain1.tempNum][chain2.tempNum] += 1

                        lnR = (1.0 / (1.0 + (self.chainTemps[chain1.tempNum]))
                                ) * chain2.curTree.logLike
                        lnR += (1.0 / (1.0 + (self.chainTemps[chain2.tempNum]))
                                ) * chain1.curTree.logLike
                        lnR -= (1.0 / (1.0 + (self.chainTemps[chain1.tempNum]))
                                ) * chain1.curTree.logLike
                        lnR -= (1.0 / (1.0 + (self.chainTemps[chain2.tempNum]))
                                ) * chain2.curTree.logLike

                        # # An alternative calculation
                        # heatBeta1 = 1.0 / (1.0 + self.chainTemps[chain1.tempNum])
                        # heatBeta2 = 1.0 / (1.0 + self.chainTemps[chain2.tempNum])
                        # likeRatio12 = (chain2.curTree.logLike - chain1.curTree.logLike) * heatBeta1
                        # likeRatio21 = (chain1.curTree.logLike - chain2.curTree.logLike) * heatBeta2
                        # lnR2 = likeRatio12 + likeRatio21
                        # rDiff = math.fabs(lnR - lnR2)
                        # if rDiff > 1e-12:
                        #     print("bad swap rDiff %f (%g)   lnR=%f, lnR2=%f" % (rDiff, rDiff, lnR, lnR2))
                        # lnR = lnR2

                        if lnR < -100.0:
                            r = 0.0
                        elif lnR >= 0.0:
                            r = 1.0
                        else:
                            r = math.exp(lnR)

                        acceptSwap = 0
                        if random.random() < r:
                            acceptSwap = 1

                        # self.logger.info("swap proposed gen=%i between tempNum1=%i chNum1=%i temp1=%f and tempNum2=%i chNum2=%i temp2=%f acceptSwap=%s" % (
                        #     self.gen, rTempNum1, chain1.chNum, self.chainTemps[chain1.tempNum], 
                        #     rTempNum2, chain2.chNum, self.chainTemps[chain2.tempNum], acceptSwap))

                        # for continuous temperature tuning with self.swapTuner
                        if self.swapTuner:
                            # Index the nAttempts and nSwaps with the lower of the two tempNum's, which would be chain1.tempNum
                            self.swapTuner.nAttempts[chain1.tempNum] += 1
                            if acceptSwap:
                                self.swapTuner.nSwaps[chain1.tempNum] += 1
                            if self.swapTuner.nAttempts[chain1.tempNum] >= var.mcmc_swapTunerSampleSize:
                                self.swapTuner.tune(chain1.tempNum)
                                # tune() zeros nAttempts and nSwaps counters

                        if acceptSwap:
                            # Use the lower triangle of swapMatrix to keep track of
                            # nAccepted's
                            assert chain1.tempNum < chain2.tempNum
                            self.swapMatrix[chain2.tempNum][chain1.tempNum] += 1

                            # Do the swap
                            chain1.tempNum, chain2.tempNum = chain2.tempNum, chain1.tempNum


                    else:     # swap matrix
                        # Chain swapping stuff was lifted from MrBayes.  Thanks again.
                        chain1, chain2 = random.sample(self.chains, 2)

                        # Use the upper triangle of swapMatrix for nProposed's
                        if chain1.tempNum < chain2.tempNum:
                            self.swapMatrix[chain1.tempNum][chain2.tempNum] += 1
                            thisCh1Temp = chain1.tempNum
                            thisCh2Temp = chain2.tempNum
                        else:
                            self.swapMatrix[chain2.tempNum][chain1.tempNum] += 1
                            thisCh1Temp = chain2.tempNum
                            thisCh2Temp = chain1.tempNum

                        lnR = (1.0 / (1.0 + (self.chainTemp * chain1.tempNum))
                                ) * chain2.curTree.logLike
                        lnR += (1.0 / (1.0 + (self.chainTemp * chain2.tempNum))
                                ) * chain1.curTree.logLike
                        lnR -= (1.0 / (1.0 + (self.chainTemp * chain1.tempNum))
                                ) * chain1.curTree.logLike
                        lnR -= (1.0 / (1.0 + (self.chainTemp * chain2.tempNum))
                                ) * chain2.curTree.logLike

                        if lnR < -100.0:
                            r = 0.0
                        elif lnR >= 0.0:
                            r = 1.0
                        else:
                            r = math.exp(lnR)

                        acceptSwap = 0
                        if random.random() < r:
                            acceptSwap = 1

                        # for continuous temperature tuning with self.swapTuner
                        if self.swapTuner and thisCh1Temp == 0 and thisCh2Temp == 1:
                            self.swapTuner.swaps01_nAttempts += 1
                            if acceptSwap:
                                self.swapTuner.swaps01_nSwaps += 1
                            if self.swapTuner.swaps01_nAttempts >= self.swapTuner.sampleSize:
                                self.swapTuner.tune(self)
                                # tune() zeros nAttempts and nSwaps counters

                        if acceptSwap:
                            # Use the lower triangle of swapMatrix to keep track of
                            # nAccepted's
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
                        raise P4Error(gm)

            # If it is a writeInterval, write stuff
            if (self.gen + 1) % self.sampleInterval == 0:
                if writeSamples:
                    likesFile = open(self.likesFileName, 'a')
                    likesFile.write(
                        '%11i %f\n' % (self.gen + 1, self.chains[coldChainNum].curTree.logLike))
                    likesFile.close()
                    treeFile = open(self.treeFileName, 'a')
                    treeFile.write("  tree t_%i = [&U] " % (self.gen + 1))
                    self.chains[coldChainNum].curTree.writeNewick(treeFile,
                                                                  withTranslation=1,
                                                                  translationHash=self.translationHash,
                                                                  doMcmcCommandComments=False)
                    treeFile.close()

                if writeSamples and self.writePrams:
                    pramsFile = open(self.pramsFileName, 'a')
                    #pramsFile.write("%12i " % (self.gen + 1))
                    pramsFile.write("%12i" % (self.gen + 1))
                    if self.modelName.startswith("SR2008"):
                        pramsFile.write(
                            "  %f\n" % self.chains[coldChainNum].curTree.beta)
                    elif self.modelName in ["SPA", "QPA"]:
                        pramsFile.write(
                            "  %f\n" % self.chains[coldChainNum].curTree.spaQ)
                    pramsFile.close()

                # Do a simulation
                if self.simulate:
                    # print "about to simulate..."
                    self.doSimulate(self.chains[coldChainNum].curTree)
                    # print "...finished simulate."

                # Do other stuff.
                if hasattr(self, 'hook'):
                    self.hook(self.chains[coldChainNum].curTree)

                if 0 and self.constraints:
                    print("Mcmc x1c")
                    print(self.chains[0].verifyIdentityOfTwoTreesInChain())
                    print("b checking curTree ..")
                    self.chains[0].curTree.checkSplitKeys()
                    print("b checking propTree ...")
                    self.chains[0].propTree.checkSplitKeys()
                    print("Mcmc xxx")

                # Add curTree to treePartitions
                if self.treePartitions:
                    self.treePartitions._getSplitsFromTree(
                        self.chains[coldChainNum].curTree)
                else:
                    self.treePartitions = TreePartitions(
                        self.chains[coldChainNum].curTree)
                # After _getSplitsFromTree, need to follow, at some point,
                # with _finishSplits().  Do that when it is pickled, or at the
                # end of the run.

                # Checking and debugging constraints
                if 0 and self.constraints:
                    print("Mcmc x1d")
                    print(self.chains[coldChainNum].verifyIdentityOfTwoTreesInChain())
                    print("c checking curTree ...")
                    self.chains[coldChainNum].curTree.checkSplitKeys()
                    print("c checking propTree ...")
                    self.chains[coldChainNum].propTree.checkSplitKeys()
                    # print "c checking that all constraints are present"
                    #theSplits = [n.br.splitKey for n in self.chains[0].curTree.iterNodesNoRoot()]
                    # for sk in self.constraints.constraints:
                    #    if sk not in theSplits:
                    #        gm.append("split %i is not present in the curTree." % sk)
                    #        raise P4Error(gm)
                    print("Mcmc zzz")

                # Check that the curTree has all the constraints
                if self.constraints:
                    splitsInCurTree = [
                        n.br.splitKey for n in self.chains[coldChainNum].curTree.iterInternalsNoRoot()]
                    for sk in self.constraints.constraints:
                        if sk not in splitsInCurTree:
                            gm.append("Programming error.")
                            gm.append(
                                "The current tree (the last tree sampled) does not contain constraint")
                            gm.append(
                                "%s" % p4.func.getSplitStringFromKey(sk, self.tree.nTax))
                            raise P4Error(gm)

                # If it is a checkPointInterval, pickle
                if self.checkPointInterval and (self.gen + 1) % self.checkPointInterval == 0:
                    self.checkPoint()

                    # The stuff below needs to be done in a re-start as well.
                    # See above "if self.proposals:"
                    self.startMinusOne = self.gen

                    # Start the tree partitions over.
                    self.treePartitions = None
                    # Zero the proposal counts
                    for p in self.props.proposals:
                        p.nProposals = [0] * self.nChains
                        p.nAcceptances = [0] * self.nChains
                        #p.nTopologyChangeAttempts = [0] * self.nChains
                        #p.nTopologyChanges = [0] * self.nChains
                        p.nAborts = [0] * self.nChains
                    # Zero the swap matrix
                    if self.nChains > 1:
                        self.swapMatrix = []
                        for i in range(self.nChains):
                            self.swapMatrix.append([0] * self.nChains)

            # Reassuring pips ...
            # We want to skip the first gen of every call to run()
            if firstGen != self.gen:
                if nGensToDo <= 20000:
                    if (self.gen - firstGen) % 1000 == 0:
                        if verbose:
                            deltaTime = self._doTimeCheck(
                                nGensToDo, firstGen, 1000)
                            if deltaTime.days:
                                timeString = "%s days, %s" % (
                                    deltaTime.days, time.strftime("%H:%M:%S", time.gmtime(deltaTime.seconds)))
                            else:
                                timeString = time.strftime(
                                    "%H:%M:%S", time.gmtime(deltaTime.seconds))
                            print("%10i - %s" % (self.gen, timeString))

                        else:
                            sys.stdout.write(".")
                            sys.stdout.flush()
                    elif (self.gen - firstGen) % 100 == 0:
                        sys.stdout.write(".")
                        sys.stdout.flush()
                else:
                    if (self.gen - firstGen) % 50000 == 0:
                        if verbose:
                            deltaTime = self._doTimeCheck(
                                nGensToDo, firstGen, 50000)
                            if deltaTime.days:
                                timeString = "%s days, %s" % (
                                    deltaTime.days, time.strftime("%H:%M:%S", time.gmtime(deltaTime.seconds)))
                            else:
                                timeString = time.strftime(
                                    "%H:%M:%S", time.gmtime(deltaTime.seconds))
                            print("%10i - %s" % (self.gen, timeString))
                        else:
                            sys.stdout.write(".")
                            sys.stdout.flush()
                    elif (self.gen - firstGen) % 1000 == 0:
                        sys.stdout.write(".")
                        sys.stdout.flush()

        # Gens finished.  Clean up.
        print()
        if verbose:
            print("Finished %s generations." % nGensToDo)

        treeFile = open(self.treeFileName, 'a')
        treeFile.write('end;\n\n')
        treeFile.close()

    def _doTimeCheck(self, nGensToDo, firstGen, genInterval):
        """Time check 

        firstGen is the first generation of this call to Mcmc.run() else
        timing fails on restart"""
        nowTime = time.time()
        diff_secs = nowTime - self.lastTimeCheck
        total_secs = (float(nGensToDo - (self.gen - firstGen)) /
                      float(genInterval)) * float(diff_secs)
        deltaTime = datetime.timedelta(seconds=int(round(total_secs)))
        self.lastTimeCheck = nowTime
        return deltaTime

    def checkPoint(self):
        # Maybe we should not save the inTrees? -- would make it more
        # lightweight.
        if 0:
            for chNum in range(self.nChains):
                ch = self.chains[chNum]
                print("chain %i ==================" % chNum)
                ch.curTree.summarizeModelComponentsNNodes()

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

        # The logger does not pickle
        savedLogger = self.logger
        self.logger = None
        savedLoggerPrinter = self.loggerPrinter
        self.loggerPrinter = None

        # The FastSpa stuff does not pickle
        if self.modelName.startswith("SPA") and var.stmcmc_useFastSpa:
            savedFspa = self.fspa
            self.fspa = None


        # _io.TextIOWrapper objects (as returned by open(fileName)) cannot be
        # copied ("serialized"), even if they are closed.  So save them and
        # restore them.
        # Except that this is not needed, as it is not really used.  But I will
        # leave this comment here for when I do start to use it.
        # savedTreeFile = self.treeFile

        theCopy = copy.deepcopy(self)

        self.logger = savedLogger
        self.loggerPrinter = savedLoggerPrinter

        theCopy.treePartitions._finishSplits()
        # assert theCopy.treeFile == None
        # theCopy.treePartitions = None    # this can be the biggest part of
        # the pickle.

        # Pickle it.
        fName = "mcmc_checkPoint_%i.%i" % (self.runNum, self.gen + 1)
        f = open(fName, 'wb')
        pickle.dump(theCopy, f, pickle.HIGHEST_PROTOCOL)
        f.close()

        if self.stRFCalc == 'fastReducedRF':
            for chNum in range(self.nChains):
                ch = self.chains[chNum]
                ch.frrf = savedFrrfs[chNum]
                ch.bigTr = savedBigTrs[chNum]

        if self.modelName.startswith("SPA") and var.stmcmc_useFastSpa:
            self.fspa = savedFspa



    def writeProposalProbs(self):
        """(Another) Pretty-print the proposal probabilities.

        See also STMcmc.writeProposalAcceptances().
        """

        nProposals = len(self.props.proposals)
        if not nProposals:
            print("STMcmc.writeProposalProbs().  No proposals (yet?).")
            return

        nAttained = [0] * nProposals
        nAccepted = [0] * nProposals
        for i in range(nProposals):
            nAttained[i] = self.props.proposals[i].nProposals[0]
            nAccepted[i] = self.props.proposals[i].nAcceptances[0]
        sumAttained = float(sum(nAttained))  # should be zero or nGen
        if not sumAttained:
            print("STMcmc.writeProposalProbs().  No proposals have been made.")
            print("Possibly, due to it being a checkPoint interval, nProposals have all been set to zero.")
            return
        # assert int(sumAttained) == self.gen + 1, "sumAttained is %i, should be gen+1, %i." % (
        #    int(sumAttained), self.gen + 1)
        probAttained = []
        for i in range(len(nAttained)):
            probAttained.append(100.0 * float(nAttained[i]) / sumAttained)
        if math.fabs(sum(probAttained) - 100.0 > 1e-13):
            raise P4Error(
                "bad sum of attained proposal probs. %s" % sum(probAttained))

        spacer = ' ' * 4
        print("\nProposal probabilities (%)")
        # print "There are %i proposals" % len(self.proposals)
        print("For %i gens, from gens %i to %i, inclusive." % (
            (self.gen - self.startMinusOne), self.startMinusOne + 1, self.gen))
        print("%2s %11s %11s  %11s %10s %23s" % ('', 'nProposals', 'proposed(%)',
                                                         'accepted(%)', 'tuning', 'proposal'))
        for i in range(len(self.props.proposals)):
            print("%2i" % i, end=' ')
            p = self.props.proposals[i]
            print("   %7i " % self.props.proposals[i].nProposals[0], end=' ')
            print("   %5.1f    " % probAttained[i], end=' ')
            if nAttained[i]:
                print("   %5.1f   " % (100.0 * float(nAccepted[i]) / float(nAttained[i])), end=' ')
            else:
                print("       -   ", end=' ')

            if p.tuning == None:
                print("      -   ", end=' ')
            elif p.tuning[0] < 2.0:
                print("  %8.4f" % p.tuning[0], end=' ')
            elif p.tuning[0] < 20.0:
                print("  %8.3f" % p.tuning[0], end=' ')
            elif p.tuning[0] < 200.0:
                print("  %8.1f" % p.tuning[0], end=' ')
            else:
                print("  %8.3g" % p.tuning[0], end=' ')
            print("%23s " % p.name, end=' ')
            print()

    # def writeProposalIntendedProbs(self):
    #     """Tabulate the intended proposal probabilities.
    #     """

    #     nProposals = len(self.proposals)
    #     if not nProposals:
    #         print("STMcmc.writeProposalIntendedProbs().  No proposals (yet?).")
    #         return
    #     intended = self.propWeights[:]
    #     for i in range(len(intended)):
    #         intended[i] /= self.totalPropWeights
    #     if math.fabs(sum(intended) - 1.0 > 1e-14):
    #         raise P4Error(
    #             "bad sum of intended proposal probs. %s" % sum(intended))

    #     spacer = ' ' * 4
    #     print("\nIntended proposal probabilities (%)")
    #     # print "There are %i proposals" % len(self.proposals)
    #     print("%2s %11s %23s %5s %5s" % ('', 'intended(%)', 'proposal', 'part', 'num'))
    #     for i in range(len(self.proposals)):
    #         print("%2i" % i, end=' ')
    #         p = self.proposals[i]
    #         print("   %6.2f    " % (100. * intended[i]), end=' ')

    #         print(" %20s" % p.name, end=' ')
    #         if p.pNum != -1:
    #             print(" %3i " % p.pNum, end=' ')
    #         else:
    #             print("   - ", end=' ')
    #         if p.mtNum != -1:
    #             print(" %3i " % p.mtNum, end=' ')
    #         else:
    #             print("   - ", end=' ')
    #         print()



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
            # print "Full glob = %s" % fList
            fList = [fName for fName in glob.glob(theGlob) if
                     os.path.basename(fName).startswith("mcmc_checkPoint")]
            # print fList
            if not fList:
                raise P4Error("No checkpoints found in this directory.")
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
                f = open(mostRecentFileName, 'rb')
                m = pickle.load(f)
                f.close()
                self.mm.append(m)

            else:
                # get all the files
                for fName in fList:
                    f = open(fName, 'rb')
                    m = pickle.load(f)
                    f.close()
                    self.mm.append(m)

                self.mm = p4.func.sortListOfObjectsOn2Attributes(
                    self.mm, "gen", 'runNum')
        else:
            # get the file by name
            f = open(fName, 'rb')
            m = pickle.load(f)
            f.close()
            self.mm.append(m)
        if verbose:
            self.dump()

    def dump(self):
        print("STMcmcCheckPoints (%i checkPoints read)" % len(self.mm))
        print("%12s %12s %12s %12s" % (" ", "index", "run", "gen+1"))
        print("%12s %12s %12s %12s" % (" ", "-----", "---", "-----"))
        for i in range(len(self.mm)):
            m = self.mm[i]
            # print "    %2i    run %2i,  gen+1 %11i" % (i, m.runNum, m.gen+1)
            print("%12s %12s %12s %12s" % (" ", i, m.runNum, m.gen + 1))

    def compareSplits(self, mNum1, mNum2, verbose=True, minimumProportion=0.1):
        """Do the TreePartitions.compareSplits() method between two checkpoints 

        Args:
            mNum1, mNum2 (int): indices to STMcmc checkpoints in self

        Returns:
            a tuple of asdoss and the maximum difference in split supports

        """

        # Should we be only looking at splits within the 95% ci of the topologies?
        m1 = self.mm[mNum1]
        m2 = self.mm[mNum2]
        tp1 = m1.treePartitions
        tp2 = m2.treePartitions

        if verbose:
            print("\nSTMcmcCheckPointReader.compareSplits(%i,%i)" % (mNum1, mNum2))
            print("%12s %12s %12s %12s %12s" % ("mNum", "runNum", "start", "gen+1", "nTrees"))
            for i in range(5):
                print("   ---------", end=' ')
            print()
            for mNum in [mNum1, mNum2]:
                print(" %10i " % mNum, end=' ')
                m = self.mm[mNum]
                print(" %10i " % m.runNum, end=' ')
                print(" %10i " % (m.startMinusOne + 1), end=' ')
                print(" %10i " % (m.gen + 1), end=' ')
                # for i in m.splitCompares:
                #    print i
                print(" %10i " % m.treePartitions.nTrees)

        asdos, maxDiff, meanDiff = self.compareSplitsBetweenTwoTreePartitions(
            tp1, tp2, minimumProportion, verbose=verbose)
        asdos2, maxDiff2, meanDiff2 = self.compareSplitsBetweenTwoTreePartitions(
            tp2, tp1, minimumProportion, verbose=verbose)
        if math.fabs(asdos - asdos2) > 0.000001:
            print("Reciprocal assdos differs:  %s  %s" % (asdos, asdos2))

        if asdos == None and verbose:
            print("No splits > %s" % minimumProportion)
        return asdos, maxDiff, meanDiff

    def compareSplitsBetweenTwoTreePartitions(tp1, tp2, minimumProportion, verbose=False):
        """Returns a tuple of asdoss, maximum of the differences and mean of the differences

        This calls the method TreePartitions.compareSplits(), and digests the
        results returned from that.

        Args:
            tp1, tp2 (TreePartition): TreePartition objects
            minimumProportion (float): passed to TreePartitions.compareSplits()
        
        Returns:
            (asdoss, maxOfDiffs, meanOfDiffs)

        """

        ret = tp1.compareSplits(tp2, minimumProportion=minimumProportion)

        #print(ret)  # a list of 3-item lists
        #  1. The split key
        #  2. The split string
        #  3. A list of the 2 supports

        if not ret:
            return None

        sumOfStdDevs = 0.0
        nSplits = len(ret)
        diffs = []
        for i in ret:
            # print "            %.3f  %.3f    " % (i[2][0], i[2][1]),
            stdDev = math.sqrt(p4.func.variance(i[2]))
            # print "%.5f" % stdDev
            sumOfStdDevs += stdDev
            diffs.append(math.fabs(i[2][0] - i[2][1]))
        asdoss = sumOfStdDevs / nSplits
        maxOfDiffs = max(diffs)
        meanOfDiffs = sum(diffs) / nSplits
        if verbose:
            print("     nSplits=%i, average of std devs of split supports %.4f " % (nSplits, asdoss))
            print("     max of differences %f, mean of differences %f" % (maxOfDiffs, meanOfDiffs))
        return (asdoss, maxOfDiffs, meanOfDiffs)  

    compareSplitsBetweenTwoTreePartitions = staticmethod(
        compareSplitsBetweenTwoTreePartitions)

    def compareSplitsAll(self, precision=3, linewidth=120):
        """Do the compareSplits() method between all pairs

        Output is verbose.  Shows 
        - average standard deviation of split frequencies (or supports), like MrBayes
        - maximum difference between split supports from each pair of checkpoints, like PhyloBayes

        Returns:
            None

        """
        nM = len(self.mm)
        nItems = int(((nM * nM) - nM) / 2)
        asdosses = np.zeros((nM, nM), dtype=np.float64)
        vect = np.zeros(nItems, dtype=np.float64)
        maxDiffs = np.zeros((nM, nM), dtype=np.float64)

        vCounter = 0
        for mNum1 in range(1, nM):
            for mNum2 in range(mNum1):
                thisAsdoss, thisMaxDiff, thisMeanDiff = self.compareSplits(mNum1, mNum2, verbose=False)
                #print("+++ thisAsdoss = %s  thisMaxDiff=%f, mNum1=%i, mNum2=%i" % (
                #      thisAsdoss, thisMaxDiff, mNum1, mNum2))
                if thisAsdoss == None:
                    thisAsdoss = 0.0
                asdosses[mNum1][mNum2] = thisAsdoss
                asdosses[mNum2][mNum1] = thisAsdoss
                vect[vCounter] = thisAsdoss
                vCounter += 1
                maxDiffs[mNum1][mNum2] = thisMaxDiff
                maxDiffs[mNum2][mNum1] = thisMaxDiff

                if 0:
                    print(" %10i " % mNum1, end=' ')
                    print(" %10i " % mNum2, end=' ')
                    print("%.3f" % thisAsdoss)

        # Save current numpy printoptions, and restore, below.
        curr = np.get_printoptions()
        np.set_printoptions(precision=precision, linewidth=linewidth)
        print("Pairwise asdoss values ---")
        print(asdosses)
        print()
        print("For the %i values in one triangle," % nItems)
        print("max =  ", vect.max())
        print("min =  ", vect.min())
        print("mean = ", vect.mean())
        print("var =  ", vect.var())

        print()
        print("Pairwise maximum differences in split supports between the two runs ---")
        print(maxDiffs)

        # Reset printoptions back to what it was
        np.set_printoptions(
            precision=curr['precision'], linewidth=curr['linewidth'])


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


class QpaML(object):
    """Uses STMcmc to do likelihood calcs and Q optimization."""
    def __init__(self, inTrees, bigT):
        assert inTrees
        ttDupes = []
        for t in inTrees:
            ttDupes.append(t.dupe())
        

        assert bigT
        bigTDupe = bigT.dupe()
        if bigTDupe.taxNames:
            pass
        else:
            raise P4Error('The bigT needs taxNames')

        stm = STMcmc(ttDupes, bigT=bigTDupe, modelName='QPA',
                     beta=1.0, spaQ=0.5, stRFCalc='purePython1',
                     nChains=1, runNum=0, sampleInterval=100,
                     checkPointInterval=None, useSplitSupport=False, verbose=False,
                     checkForOutputFiles=False)
        self.ch = STChain(stm, 0)



    def setSuperTree(self, st):
        assert self.ch.propTree.taxNames
        st = st.dupe()
        if st.taxNames:
            assert st.taxNames == self.ch.propTree.taxNames
        else:
            st.taxNames = self.ch.propTree.taxNames
        st.setPreAndPostOrder()
        #st.draw()

        st.taxBits = [1 << i for i in range(st.nTax)]
        for n in st.iterPostOrder():
            if n == st.root:
                break
            if n.isLeaf:
                spot = st.taxNames.index(n.name)
                n.stSplitKey = 1 << spot
            else:
                n.stSplitKey = n.leftChild.stSplitKey
                p = n.leftChild.sibling
                while p:
                    n.stSplitKey |= p.stSplitKey    # "or", in-place
                    p = p.sibling
        st.skk = [n.stSplitKey for n in st.iterInternalsNoRoot()]
        st.qSet = set()
        for sk in st.skk:
            ups = [txBit for txBit in st.taxBits if (sk & txBit)]
            downs = [txBit for txBit in st.taxBits if not (sk & txBit)]
            for down in itertools.combinations(downs, 2):
                assert down[0] < down[1]   # probably not needed
                for up in itertools.combinations(ups, 2):
                    assert up[0] < up[1]  # probably not needed
                    if down[0] < up[0]:
                        st.qSet.add(down + up)
                    else:
                        st.qSet.add(up + down)
        # print st.qSet
        st.nQuartets = len(st.qSet)

        self.ch.propTree = st
        

                
    def calcP(self, Q):
        if Q >= 1.:
            return 10000000.
        if Q <= 0.:
            return 10000000.
        self.ch.propTree.spaQ = Q
        self.ch.getTreeLogLike_qpa_slow()
        return -self.ch.propTree.logLike


    def optimizeQ(self, x0=0.3):
        res = minimize(self.calcP, x0, method='Nelder-Mead')
        return (res.x, res.fun)



class SpaML(object):
    """Using STMcmc."""

    def __init__(self, inTrees, bigT):
        assert inTrees

        ttDupes = []
        for t in inTrees:
            ttDupes.append(t.dupe())
        

        assert bigT
        bigTDupe = bigT.dupe()
        if bigTDupe.taxNames:
            pass
        else:
            raise P4Error('The bigT needs taxNames')

        stm = STMcmc(ttDupes, bigT=bigTDupe, modelName='SPA',
                     beta=1.0, spaQ=0.5, stRFCalc='purePython1',
                     nChains=1, runNum=0, sampleInterval=100,
                     checkPointInterval=None, useSplitSupport=False, verbose=False, 
                     checkForOutputFiles=False)
        self.ch = STChain(stm, 0)


    def setSuperTree(self, st):
        assert self.ch.propTree.taxNames
        st = st.dupe()
        if st.taxNames:
            assert st.taxNames == self.ch.propTree.taxNames
        else:
            st.taxNames = self.ch.propTree.taxNames
        st.setPreAndPostOrder()

        self.ch.propTree = st
        self.ch.setupBitarrayCalcs()
        
                
    def calcP(self, Q):
        #if not isinstance(Q, np.ndarray):
        #    Q = np.array([Q])
        if Q >= 1.:
            return 10000000.
        if Q <= 0.:
            return 10000000.

        if not isinstance(Q, np.ndarray):
            Q = np.array([Q])
        self.ch.propTree.spaQ = Q
        self.ch.getTreeLogLike_spa_bitarray()
        return -self.ch.propTree.logLike


    def optimizeQ(self, x0=0.3):
        res = minimize(self.calcP, x0, method='Nelder-Mead')
        return (res.x, res.fun)


class SR2008ML(object):
    """Using STMcmc, with SR2008_rf_aZ_fb."""

    def __init__(self, inTrees, bigT):
        assert inTrees

        ttDupes = []
        for t in inTrees:
            ttDupes.append(t.dupe())
        

        assert bigT
        bigTDupe = bigT.dupe()
        if bigTDupe.taxNames:
            pass
        else:
            raise P4Error('SR2008ML: The bigT needs taxNames')

        stm = STMcmc(ttDupes, bigT=bigTDupe, modelName='SR2008_rf_aZ_fb',
                     beta=1.0, spaQ=0.5, 
                     stRFCalc='purePython1',
                     #stRFCalc='bitarray',
                     nChains=1, runNum=0, sampleInterval=100,
                     checkPointInterval=None, useSplitSupport=False, verbose=False, 
                     checkForOutputFiles=False)
        self.ch = STChain(stm, 0)


    def setSuperTree(self, st):
        assert self.ch.propTree.taxNames
        st = st.dupe()
        if st.taxNames:
            assert st.taxNames == self.ch.propTree.taxNames
        else:
            st.taxNames = self.ch.propTree.taxNames
        st.setPreAndPostOrder()

        self.ch.propTree = st
        #self.ch.setupBitarrayCalcs()
        
                
    def calcP(self, beta):
        myMIN = 1.e-10
        myMAX = 1.e+10
        if beta >= myMAX:
            return 10000000.
        if beta <= myMIN:
            return 10000000.
        self.ch.propTree.beta = beta
        self.ch.getTreeLogLike_ppy1()  # pure python
        return -self.ch.propTree.logLike


    def optimizeBeta(self, x0=1.0):
        res = minimize(self.calcP, x0, method='Nelder-Mead')
        return (res.x, res.fun)
    
