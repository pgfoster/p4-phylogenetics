from __future__ import print_function
import p4.func
import p4.pf as pf
from p4.var import var
import math
import random
import copy
import numpy
import scipy
import scipy.stats
from p4.p4exceptions import P4Error
import sys

#localCalls = 0


class Chain(object):

    from p4.chain_topol import proposeRoot3, proposeBrLen, proposeAllBrLens, proposeLocal, proposeETBR_Blaise, proposeESPR_Blaise, proposeETBR, proposeESPR, proposePolytomy, proposeAddEdge, _getCandidateNodesForDeleteEdge, proposeDeleteEdge


    def __init__(self, aMcmc):
        self.mcmc = aMcmc
        #self.num = -1
        self.tempNum = -1  # 'temp'erature, not 'temp'orary

        self.curTree = aMcmc.tree.dupe()
        self.curTree.data = aMcmc.tree.data
        self.curTree.calcLogLike(verbose=0)

        self.propTree = aMcmc.tree.dupe()
        self.propTree.data = aMcmc.tree.data
        self.propTree.calcLogLike(verbose=0)

        #print("Chain.init() curTree %f, propTree %f" % (
        #    self.curTree.logLike, self.propTree.logLike))

        # Oddly, the curTree and the propTree can have slightly
        # different condLikes and bigPDecks at this point.  The
        # difference being only a bit more than 1e-15.  However, that
        # is the epsilon that I use to test whether floating point
        # numbers are the same.  Rather than changing my epsilon, I
        # will simply copy the numbers over, so that we are starting
        # with identical numbers.
        # 1 means do all
        pf.p4_copyCondLikes(self.curTree.cTree, self.propTree.cTree, 1)
        # 1 means do all
        pf.p4_copyBigPDecks(self.curTree.cTree, self.propTree.cTree, 1)
        pf.p4_copyModelPrams(self.curTree.cTree, self.propTree.cTree)

        if 0:
            self.testTree = aMcmc.tree.dupe()
            self.testTree.data = aMcmc.tree.data
            self.testTree.calcLogLike(verbose=0)

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0
        #self.logJacobian = 0.0

        #self.lastProposal = None
        ret = self.verifyIdentityOfTwoTreesInChain()
        if ret == var.DIFFERENT:
            raise P4Error(
                "Chain.init().  Programming error. The prop tree should be identical to the cur tree, and it is not.")

    def propose(self, theProposal):
        gm = ['Chain.propose()']
        # print "propose().  gen %i, About to propose %s" % (self.mcmc.gen,
        # theProposal.name)

        if theProposal.name == 'comp':
            # print "theProposal.name = comp, pNum=%i" % theProposal.pNum
            self.proposeCompWithSlider(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'compDir':
            # print "theProposal.name = comDir, pNum=%i" % theProposal.pNum
            self.proposeCompWithDirichlet(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'rjComp':
            # print "theProposal.name = rjComp, pNum=%i" % theProposal.pNum
            self.proposeRjComp(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            self.propTree.setCStuff()  # for model usage info
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'rMatrix':
            self.proposeRMatrixWithSlider(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'rjRMatrix':
            self.proposeRjRMatrix(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            self.propTree.setCStuff()  # for model usage info
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'gdasrv':
            self.proposeGdasrv(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'pInvar':
            self.proposePInvar(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'compLocation':
            self.proposeCompLocation(theProposal)
            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'rMatrixLocation':
            self.proposeRMatrixLocation(theProposal)
            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'gdasrvLocation':
            self.proposeGdasrvLocation(theProposal)
            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'local':
            self.proposeLocal(theProposal)
            if theProposal.doAbort:
                pass
            else:
                if not self.propTree.preAndPostOrderAreValid:
                    self.propTree.setPreAndPostOrder()
                self.propTree.setCStuff()
                pf.p4_setPrams(self.propTree.cTree, -1)

        elif theProposal.name == 'brLen':
            self.proposeBrLen(theProposal)
            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, -1)

        elif theProposal.name == 'eTBR':
            self.proposeETBR_Blaise(theProposal)
            #self.proposeETBR(theProposal)
            if theProposal.doAbort:
                pass
            else:
                if not self.propTree.preAndPostOrderAreValid:
                    self.propTree.setPreAndPostOrder()
                self.propTree.setCStuff()
                pf.p4_setPrams(self.propTree.cTree, -1)

        elif theProposal.name == 'polytomy':
            self.proposePolytomy(theProposal)
            if not self.propTree.preAndPostOrderAreValid:
                self.propTree.setPreAndPostOrder()
            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, -1)

        elif theProposal.name == 'root3':
            self.proposeRoot3(theProposal)
            if not self.propTree.preAndPostOrderAreValid:
                self.propTree.setPreAndPostOrder()
            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, -1)

        elif theProposal.name == 'relRate':
            self.proposeRelRate(theProposal)
            self.propTree.model.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, -1)

        else:
            gm.append('Unlisted proposal.name=%s  Fix me.' % theProposal.name)
            raise P4Error(gm)

        if theProposal.doAbort:
            return 0.0
        else:

            # Topology moves may set n.br.lenChanged or n.flag, which
            # are not used in this slow method.  But they cause
            # problems later in checking and debugging.  So un-set
            # them.
            if 1:
                if theProposal.name in ['local', 'brLen', 'eTBR', 'polytomy']:
                    for n in self.propTree.iterNodes():
                        if n.br:
                            n.br.lenChanged = 0
                        n.flag = 0

            # print "...about to calculate the likelihood of the propTree."
            self.propTree.logLike = pf.p4_treeLogLike(self.propTree.cTree, 0)
            #self.propTree.logLike = self.getTreeLogLike()
            # print "propTree logLike is", self.propTree.logLike

            # slow check
            if 0:
                # _commonCStuff() has these ...
                #    self.model.setCStuff()
                #    self.setCStuff()
                #    #print "about to p4_setPrams()..."
                #    pf.p4_setPrams(self.cTree, -1) # "-1" means do all parts
                if 1:
                    firstCalc = self.propTree.logLike
                    # with _commonCStuff()
                    self.propTree.calcLogLike(verbose=0)
                    theDiff = math.fabs(firstCalc - self.propTree.logLike)
                if 0:
                    self.propTree.copyToTree(self.testTree)
                    self.propTree.model.copyValsTo(self.testTree.model)
                    self.testTree.calcLogLike(verbose=0)
                    theDiff = math.fabs(
                        self.testTree.logLike - self.propTree.logLike)
                # print "%g" % theDiff
                if theDiff > 1.e-9:
                    gm.append("Chain.propose().  '%s' Bad like calc.  theDiff = %g" % (
                        theProposal.name, theDiff))
                    raise P4Error(gm)

            logLikeRatio = self.propTree.logLike - self.curTree.logLike

            # To run "without the data", which shows the effect of priors.
            #logLikeRatio = 0.0

            if self.mcmc.nChains > 1:
                heatBeta = 1.0 / (1.0 + self.mcmc.tunings.chainTemp * self.tempNum)
                logLikeRatio *= heatBeta
                self.logPriorRatio *= heatBeta

            theSum = logLikeRatio + self.logProposalRatio + self.logPriorRatio
            if theProposal.name in ['rjComp', 'rjRMatrix']:
                theSum += self.logJacobian
            # if theProposal.name in ['rjComp']:
            #    print "%s: %10.2f %10.2f %10.2f %10.2f" % (theProposal.name, logLikeRatio,
            #    self.logPriorRatio, self.logProposalRatio, self.logJacobian)
            # if theProposal.name == 'polytomy':
            #    theSum += self.logJacobian
            #    self.logJacobian = 0.0
            # print "logLikeRatio = %f" % logLikeRatio
            return theSum

    def proposeSp(self, theProposal):
        gm = ['Chain.proposeSp()']
        # if self.mcmc.gen > 1300:
        #print "proposeSp().  gen %i, About to propose %s" % (self.mcmc.gen,
        #theProposal.name)

        if theProposal.name == 'comp':
            # print "theProposal.name = comp, pNum=%i" % theProposal.pNum
            self.proposeCompWithSlider(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'compDir':
            # print "theProposal.name = compDir, pNum=%i" % theProposal.pNum
            self.proposeCompWithDirichlet(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'allCompsDir':
            #print("theProposal.name = allCompsDir, pNum=%i" % theProposal.pNum)
            self.proposeAllCompsDir(theProposal)
            # This next line is not needed because comps are numpy arrays
            # self.propTree.model.setCStuff(partNum=theProposal.pNum)
            if 0:
                mpProp = self.propTree.model.parts[theProposal.pNum]
                for mtProp in mpProp.comps:
                    mySum = numpy.sum(mtProp.val)
                    myDiff = math.fabs(1.0 - mySum)
                    if myDiff > 1e-15:
                        print("Chain.proposeSp() gen %i, allCompsDir, x myDiff is %g" % (self.mcmc.gen, myDiff)) 
                for mtPropNum in range(len(mpProp.comps)):
                    mtProp = mpProp.comps[mtPropNum]
                    for chNum in range(mpProp.dim):
                        thisNp = mtProp.val[chNum]
                        thatNp = pf.test(self.propTree.cTree, theProposal.pNum, mtPropNum, chNum)
                        if math.fabs(thisNp - thatNp) > 1e-14:
                           print('++++++ gen %i comp %2i %2i' % (self.mcmc.gen, mtPropNum,chNum), "%17.15f %17.15f %g" % (thisNp, thatNp, (thisNp - thatNp))) 

            # This next line is needed, and it needs to go here.  At least
            # because of p4_calculateBigPDecksPart()
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)  # "-1" means do all parts
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)


        elif theProposal.name == 'ndch2_leafCompsDir':
            #print("theProposal.name = ndch2_leafCompsDir, pNum=%i" % theProposal.pNum)
            self.proposeNdch2_leafCompsDir(theProposal)
            # This next line is not needed because comps are numpy arrays
            # self.propTree.model.setCStuff(partNum=theProposal.pNum)
            if 0:
                mpProp = self.propTree.model.parts[theProposal.pNum]
                for mtProp in mpProp.comps:
                    mySum = numpy.sum(mtProp.val)
                    myDiff = math.fabs(1.0 - mySum)
                    if myDiff > 1e-15:
                        print("Chain.proposeSp() gen %i, ndch2_leafCompsDir, x myDiff is %g" % (self.mcmc.gen, myDiff)) 
                for mtPropNum in range(len(mpProp.comps)):
                    mtProp = mpProp.comps[mtPropNum]
                    for chNum in range(mpProp.dim):
                        thisNp = mtProp.val[chNum]
                        thatNp = pf.test(self.propTree.cTree, theProposal.pNum, mtPropNum, chNum)
                        if math.fabs(thisNp - thatNp) > 1e-14:
                           print('++++++ gen %i comp %2i %2i' % (self.mcmc.gen, mtPropNum,chNum), "%17.15f %17.15f %g" % (thisNp, thatNp, (thisNp - thatNp))) 
            
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)  # "-1" means do all parts
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)


        elif theProposal.name == 'ndch2_internalCompsDir':
            #print("theProposal.name = ndch2_internalCompsDir, pNum=%i" % theProposal.pNum)
            self.proposeNdch2_internalCompsDir(theProposal)
            # This next line is not needed because comps are numpy arrays
            # self.propTree.model.setCStuff(partNum=theProposal.pNum)
            if 0:
                mpProp = self.propTree.model.parts[theProposal.pNum]
                for mtProp in mpProp.comps:
                    mySum = numpy.sum(mtProp.val)
                    myDiff = math.fabs(1.0 - mySum)
                    if myDiff > 1e-15:
                        print("Chain.proposeSp() gen %i, ndch2_internalCompsDir, x myDiff is %g" % (self.mcmc.gen, myDiff)) 
                for mtPropNum in range(len(mpProp.comps)):
                    mtProp = mpProp.comps[mtPropNum]
                    for chNum in range(mpProp.dim):
                        thisNp = mtProp.val[chNum]
                        thatNp = pf.test(self.propTree.cTree, theProposal.pNum, mtPropNum, chNum)
                        if math.fabs(thisNp - thatNp) > 1e-14:
                           print('++++++ gen %i comp %2i %2i' % (self.mcmc.gen, mtPropNum,chNum), "%17.15f %17.15f %g" % (thisNp, thatNp, (thisNp - thatNp))) 
            
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)  # "-1" means do all parts
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)


        elif theProposal.name == 'ndch2_leafCompsDirAlpha':
            # print "theProposal.name = ndch2_leafCompsDirAlpha, pNum=%i" % theProposal.pNum
            self.proposeNdch2_leafCompsDirAlpha(theProposal)
            # No likelihood calcs!
            

        elif theProposal.name == 'ndch2_internalCompsDirAlpha':
            # print "theProposal.name = ndch2_internalCompsDirAlpha, pNum=%i" % theProposal.pNum
            self.proposeNdch2_internalCompsDirAlpha(theProposal)
            # No likelihood calcs!


        # elif theProposal.name == 'rjComp':
        #     # print "theProposal.name = rjComp, pNum=%i" % theProposal.pNum
        #     self.proposeRjComp(theProposal)
        #     if theProposal.doAbort:
        #         # print "abort rjComp"
        #         return 0.0
        #     # This next line transfers the newComp.val to C
        #     self.propTree.model.setCStuff(partNum=theProposal.pNum)
        #     self.propTree.setCStuff()  # for model usage info
        #     # print "about to p4_setPrams() ..."
        #     pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
        #     for n in self.propTree.iterPostOrder():
        #         if not n.isLeaf:
        #             pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
        #                 n.cNode, theProposal.pNum)
        #     pf.p4_partLogLike(self.propTree.cTree,
        #                       self.propTree.data.parts[theProposal.pNum].cPart,
        #                       theProposal.pNum, 0)

        # elif theProposal.name == 'cmd1_compDir':
        #     # print "theProposal.name = cmd1_compDir, pNum=%i" %
        #     # theProposal.pNum
        #     self.proposeCmd1CompDir(theProposal)
        #     if 1:
        #         self.propTree.model.setCStuff(partNum=theProposal.pNum)
        #         pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
        #         for n in self.propTree.iterPostOrder():
        #             if not n.isLeaf:
        #                 pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
        #                     n.cNode, theProposal.pNum)
        #         pf.p4_partLogLike(self.propTree.cTree,
        #                           self.propTree.data.parts[
        #                               theProposal.pNum].cPart,
        #                           theProposal.pNum, 0)

        # elif theProposal.name == 'cmd1_allCompDir':
        #     # print "theProposal.name = cmd1_allCompDir, pNum=%i" %
        #     # theProposal.pNum
        #     self.proposeCmd1AllCompDir(theProposal)
        #     if 1:
        #         self.propTree.model.setCStuff(partNum=theProposal.pNum)
        #         pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
        #         for n in self.propTree.iterPostOrder():
        #             if not n.isLeaf:
        #                 pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
        #                     n.cNode, theProposal.pNum)
        #         pf.p4_partLogLike(self.propTree.cTree,
        #                           self.propTree.data.parts[
        #                               theProposal.pNum].cPart,
        #                           theProposal.pNum, 0)

        # elif theProposal.name == 'cmd1_comp0Dir':
        #     # print "theProposal.name = cmd1_comp0Dir, pNum=%i" %
        #     # theProposal.pNum
        #     self.proposeCmd1Comp0Dir(theProposal)

        # elif theProposal.name == 'cmd1_alpha':
        #     # print "theProposal.name = cmd1_alpha, pNum=%i" % theProposal.pNum
        #     self.proposeCmd1Alpha(theProposal)

        elif theProposal.name in ['rMatrix', 'rMatrixDir']:
            if theProposal.name == 'rMatrix':
                self.proposeRMatrixWithSlider(theProposal)
            else:
                self.proposeRMatrixDirichlet(theProposal)

            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        # elif theProposal.name == 'rjRMatrix':
        #     # print "theProposal.name = rjRMatrix, pNum=%i" % theProposal.pNum
        #     self.proposeRjRMatrix(theProposal)
        #     if theProposal.doAbort:
        #         # print "abort rjRMatrix"
        #         return 0.0
        #     # This next line transfers the newRMatrix.val to C
        #     self.propTree.model.setCStuff(partNum=theProposal.pNum)
        #     self.propTree.setCStuff()  # for model usage info
        #     # print "about to p4_setPrams() ..."
        #     pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
        #     for n in self.propTree.iterPostOrder():
        #         if not n.isLeaf:
        #             pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
        #                 n.cNode, theProposal.pNum)
        #     pf.p4_partLogLike(self.propTree.cTree,
        #                       self.propTree.data.parts[theProposal.pNum].cPart,
        #                       theProposal.pNum, 0)

        elif theProposal.name == 'gdasrv':
            self.proposeGdasrv(theProposal)
            # THis next line is not needed because gdasrv is a numpy array
            #self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'pInvar':
            self.proposePInvar(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'compLocation':
            self.proposeCompLocation(theProposal)
            if theProposal.doAbort:
                return 0.0
            self.propTree.setCStuff()

            # p4_setPrams can affect bQETneedsReset (a shared numpy
            # array).  First, all cNum * rNum slots are set.  Then
            # only the combinations of cNum and rNum that are actually
            # used triggers the expensive reset.  That means that some
            # slots might remain in the state of saying that they need
            # resetting -- but they are not used.
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'rMatrixLocation':
            self.proposeRMatrixLocation(theProposal)
            if theProposal.doAbort:
                return 0.0
            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'gdasrvLocation':
            self.proposeGdasrvLocation(theProposal)
            if theProposal.doAbort:
                return 0.0
            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'local':
            self.proposeLocal(theProposal)
            if theProposal.doAbort:
                return 0.0
            else:
                if not self.propTree.preAndPostOrderAreValid:
                    self.propTree.setPreAndPostOrder()
                self.propTree.setCStuff()

                for n in self.propTree.iterNodesNoRoot():
                    if n.br.lenChanged:
                        pf.p4_calculateBigPDecks(n.cNode)
                        p = n
                        while p != self.propTree.root:
                            p = p.parent
                            p.flag = 1
                        n.br.lenChanged = False
                # for n in self.propTree.iterNodesNoRoot():
                #    if n.flag:
                #        print "    node %2i flag" % n.nodeNum
                for n in self.propTree.iterPostOrder():
                    if not n.isLeaf:
                        if n.flag:
                            for pNum in range(self.propTree.model.nParts):
                                pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                                    n.cNode, pNum)
                        n.flag = 0
                for pNum in range(self.propTree.model.nParts):
                    pf.p4_partLogLike(
                        self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)

                if 0 and self.mcmc.gen == 0:
                    self.propTree.calcLogLike()
                    self.curTree.calcLogLike()
                    self.curTree.draw()
                    self.propTree.draw()

        elif theProposal.name == 'brLen':
            self.proposeBrLen(theProposal)
            self.propTree.setCStuff()
            for n in self.propTree.iterNodesNoRoot():
                if n.br.lenChanged:
                    # print "node %i br.lenChanged" % n.nodeNum
                    pf.p4_calculateBigPDecks(n.cNode)
                    p = n
                    while p != self.propTree.root:
                        p = p.parent
                        p.flag = 1
                    n.br.lenChanged = False
                    break
            for pNum in range(self.propTree.model.nParts):
                for n in self.propTree.iterPostOrder():
                    if not n.isLeaf:
                        if n.flag:
                            pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                                n.cNode, pNum)
                pf.p4_partLogLike(
                    self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)
            for n in self.propTree.iterInternalsNoRoot():
                n.flag = 0
            self.propTree.root.flag = 0

            if 0:
                logLike1 = sum(self.propTree.partLikes)
                pf.p4_setPrams(self.propTree.cTree, -1)
                logLike2 = pf.p4_treeLogLike(self.propTree.cTree, 0)
                if math.fabs(logLike1 - logLike2) > 0.001:
                    print("propose brLen bad likes calc. %f %f" % (logLike1, logLike2))
                else:
                    print("propose brLen likes ok --  %f" % logLike1)
                sys.exit()

        elif theProposal.name == 'allBrLens':
            self.proposeAllBrLens(theProposal)
            self.propTree.setCStuff()
            for n in self.propTree.iterNodesNoRoot():  
                # all branch lengths have changed, 
                # so no need to check whether n.br.lenChanged
                pf.p4_calculateBigPDecks(n.cNode)
            for pNum in range(self.propTree.model.nParts):
                for n in self.propTree.iterPostOrder():
                    if not n.isLeaf:
                        # Normally we would check whether n.flag is set
                        # but here they all need to recalc cond likes
                        pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                            n.cNode, pNum)
                pf.p4_partLogLike(
                    self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)
            #for n in self.propTree.iterInternalsNoRoot():
            #    n.flag = 0
            #self.propTree.root.flag = 0

        # elif theProposal.name == 'treeScale':
        #     self.proposeTreeScale(theProposal)
        #     self.propTree.setCStuff()
        #     for n in self.propTree.iterNodesNoRoot():  
        #         # all branch lengths have changed, 
        #         # so no need to check whether n.br.lenChanged
        #         pf.p4_calculateBigPDecks(n.cNode)
        #     for pNum in range(self.propTree.model.nParts):
        #         for n in self.propTree.iterPostOrder():
        #             if not n.isLeaf:
        #                 # Normally we would check whether n.flag is set
        #                 # but here they all need to recalc cond likes
        #                 pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
        #                     n.cNode, pNum)
        #         pf.p4_partLogLike(
        #             self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)
        #     #for n in self.propTree.iterInternalsNoRoot():
        #     #    n.flag = 0
        #     #self.propTree.root.flag = 0

        #     if 0:
        #         logLike1 = sum(self.propTree.partLikes)
        #         pf.p4_setPrams(self.propTree.cTree, -1)
        #         logLike2 = pf.p4_treeLogLike(self.propTree.cTree, 0)
        #         if math.fabs(logLike1 - logLike2) > 0.001:
        #             print "propose treeScale bad likes calc. %f %f" % (logLike1, logLike2)
        #         else:
        #             print "propose treeScale likes ok --  %f" % logLike1
        #         sys.exit()

        elif theProposal.name == 'eTBR':
            self.proposeETBR_Blaise(theProposal)
            #self.proposeETBR(theProposal)
            #self.proposeESPR_Blaise(theProposal)
            if theProposal.doAbort:
                return 0.0
            else:
                if not self.propTree.preAndPostOrderAreValid:
                    self.propTree.setPreAndPostOrder()

                self.propTree.setCStuff()

                # Debugging litter ...
                if 0 and self.mcmc.gen == 270:
                    print()
                    self.propTree.draw()
                    #self.propTree.node(16).br.lenChanged = 1
                    #self.propTree.node(17).br.lenChanged = 1
                    for n in self.propTree.iterNodesNoRoot():
                        if n.br.lenChanged:
                            print("    node %2i br.lenChanged" % n.nodeNum)
                            #n.br.textDrawSymbol = 'C'
                    # self.propTree.draw()

                for n in self.propTree.iterNodesNoRoot():
                    if n.br.lenChanged:
                        # if 1:
                        # This next line can generate
                        # p4_calculateBigPDecksPart() pNum=0, compNum=1,
                        # rMatrixNum=0, needsReset. Fix me.
                        pf.p4_calculateBigPDecks(n.cNode)
                        p = n
                        while p != self.propTree.root:
                            p = p.parent
                            p.flag = 1
                        n.br.lenChanged = False

                # More debugging litter ...
                if 0 and self.mcmc.gen == 270:
                    for n in self.propTree.iterNodes():
                        if n.flag:
                            print("    node %2i flag" % n.nodeNum)
                            # if n.br:
                            #    n.br.textDrawSymbol = 'f'
                    # self.propTree.draw()
                    # for  n in self.propTree.iterNodesNoRoot():
                    #    n.br.textDrawSymbol = '-'

                # Recalculate condLikes for only the flagged nodes, in post
                # order down to the root.
                for n in self.propTree.iterPostOrder():
                    if not n.isLeaf:
                        if n.flag:
                            for pNum in range(self.propTree.model.nParts):
                                pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                                    n.cNode, pNum)
                        n.flag = 0
                for pNum in range(self.propTree.model.nParts):
                    pf.p4_partLogLike(
                        self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)

        elif theProposal.name == 'polytomy':
            self.proposePolytomy(theProposal)
            if theProposal.doAbort:
                return 0.0
            if not self.propTree.preAndPostOrderAreValid:
                self.propTree.setPreAndPostOrder()
            self.propTree.setCStuff()
            #pf.p4_setPrams(self.propTree.cTree, -1)
            # node.flag's were set by proposeDeleteEdge(), but no n.br.lenChanged.
            # But n.br.lenChanged set in proposeAddEdge(), but no n.flag's.
            if 0:
                #pf.p4_setPrams(self.propTree.cTree, -1)
                for n in self.propTree.iterNodesNoRoot():
                    pf.p4_calculateBigPDecks(n.cNode)
                for pNum in range(self.propTree.model.nParts):
                    for n in self.propTree.iterNodesNoRoot():
                        n.br.lenChanged = False
                        raise P4Error(
                            "this doesn't work for more than one part")
                        n.flag = 0
                    self.propTree.root.flag = 0
                    for n in self.propTree.iterPostOrder():
                        if not n.isLeaf:
                            pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                                n.cNode, pNum)
                    pf.p4_partLogLike(
                        self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)
            else:
                for n in self.propTree.iterNodesNoRoot():
                    if n.br.lenChanged:
                        pf.p4_calculateBigPDecks(n.cNode)
                        # Need to recalculate condLikes of the newly
                        # added node, as well as its ancestors.
                        # (usually just need to do the ancestors of
                        # the changed brLen).
                        n.flag = 1
                        p = n
                        while p != self.propTree.root:
                            p = p.parent
                            p.flag = 1
                        n.br.lenChanged = False
                        break
                for pNum in range(self.propTree.model.nParts):
                    for n in self.propTree.iterPostOrder():
                        if not n.isLeaf:
                            # print "node %i" % n.nodeNum
                            if n.flag:
                                pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                                    n.cNode, pNum)
                    pf.p4_partLogLike(
                        self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)
                for n in self.propTree.iterInternalsNoRoot():
                    n.flag = 0
                self.propTree.root.flag = 0

        elif theProposal.name == 'root3':
            self.proposeRoot3(theProposal)
            if not self.propTree.preAndPostOrderAreValid:
                self.propTree.setPreAndPostOrder()
            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, -1)
            for pNum in range(self.propTree.model.nParts):
                for n in self.propTree.iterPostOrder():
                    if not n.isLeaf:
                        pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                            n.cNode, pNum)
                pf.p4_partLogLike(
                    self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)

        elif theProposal.name == 'relRate':
            # print "theProposal.name = relRate, pNum=%i" % theProposal.pNum
            self.proposeRelRate(theProposal)

            # Model.setCStuff() transfers comp, rMatrix, pInvar, and
            # relRate info.  We don't need that much.
            # self.propTree.model.setCStuff()
            for pNum in range(self.propTree.model.nParts):
                mp = self.propTree.model.parts[pNum]
                pf.p4_setRelRateVal(
                    self.propTree.model.cModel, mp.num, mp.relRate)

            # pf.p4_setPrams(self.propTree.cTree, -1) recalculates all
            # Q-matrices, resets eigensystems, and recalculates
            # P-matrices for all the nodes.  We don't need the former
            # -- just need to recalculate P-matrices for all the
            # nodes.
            #pf.p4_setPrams(self.propTree.cTree, -1)
            # But it turns out the following is slower than p4_setPrams()!
            # for n in self.propTree.iterNodesNoRoot():
            #    pf.p4_calculateBigPDecks(n.cNode)
            # But the following is faster
            pf.p4_calculateAllBigPDecksAllParts(self.propTree.cTree)

            # This is the time-consuming part.
            for pNum in range(self.propTree.model.nParts):
                for n in self.propTree.iterPostOrder():
                    if not n.isLeaf:
                        pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                            n.cNode, pNum)
                pf.p4_partLogLike(
                    self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)

        else:
            gm.append('Unlisted proposal.name=%s  Fix me.' % theProposal.name)
            raise P4Error(gm)

        if theProposal.doAbort:
            raise P4Error(
                "programming error.  we should not be here.  proposal %s" % theProposal.name)

        proposalsWithNoLikeCalcs = ['ndch2_leafCompsDirAlpha', 'ndch2_internalCompsDirAlpha']

        if theProposal.name not in proposalsWithNoLikeCalcs:
            self.propTree.logLike = sum(self.propTree.partLikes)

        # Slow check.
        if 0:
            # Doing a pf.p4_treeLogLike() is expensive
            if 1:
                x = sum(self.curTree.partLikes)
                #pf.p4_partLogLike(self.curTree.cTree,self.curTree.data.parts[1].cPart, 1, 0)
                pf.p4_treeLogLike(self.curTree.cTree, 0)
                y = sum(self.curTree.partLikes)
                if math.fabs(x - y) > 0.00001:
                    print("***************************** gen %i, bad curTree here b" % self.mcmc.gen)
                    print(x, y)
                # else:
                # print "***************************** no difference to
                # curTree here b"
            if 1:
                x = sum(self.propTree.partLikes)
                #pf.p4_partLogLike(self.propTree.cTree,self.propTree.data.parts[1].cPart, 1, 0)
                pf.p4_treeLogLike(self.propTree.cTree, 0)
                y = sum(self.propTree.partLikes)
                if math.fabs(x - y) > 0.00001:
                    print("***************************** gen %i, bad propTree here b" % self.mcmc.gen)
                    print(x, y)
                # else:
                # print "***************************** no difference to
                # propTree here b"

        # slow check
        if 0:
            # _commonCStuff() has these ...
            #    self.model.setCStuff()
            #    self.setCStuff()
            #    #print "about to p4_setPrams()..."
            #    pf.p4_setPrams(self.cTree, -1) # "-1" means do all parts
            if 0:
                firstCalc = self.curTree.logLike
                # with _commonCStuff()
                self.curTree.calcLogLike(verbose=0)
                theDiff = math.fabs(firstCalc - self.curTree.logLike)
            if 1:
                firstCalc = self.propTree.logLike
                # with _commonCStuff()
                self.propTree.calcLogLike(verbose=0)
                theDiff = math.fabs(firstCalc - self.propTree.logLike)
            if 0:
                self.propTree.copyToTree(self.testTree)
                self.propTree.model.copyValsTo(self.testTree.model)
                self.testTree.calcLogLike(verbose=0)
                theDiff = math.fabs(
                    self.testTree.logLike - self.propTree.logLike)
            # print "%g" % theDiff
            if theDiff > 1.e-9:
                gm.append("gen %i, Bad like calc.  '%s', theDiff = %g" % (
                    self.mcmc.gen, theProposal.name, theDiff))
                raise P4Error(gm)

        if theProposal.name in proposalsWithNoLikeCalcs:
            logLikeRatio = 0.0
        else:
            logLikeRatio = self.propTree.logLike - self.curTree.logLike

        # To run "without the data", which shows the effect of priors.
        #logLikeRatio = 0.0

        if self.mcmc.nChains > 1:
            heatBeta = 1.0 / (1.0 + self.mcmc.tunings.chainTemp * self.tempNum)
            logLikeRatio *= heatBeta
            # print "logPriorRatio is %s, heatBeta is %s" %
            # (self.logPriorRatio, heatBeta)
            self.logPriorRatio *= heatBeta

        # Experimental Heating hack
        if self.mcmc.doHeatingHack: # and theProposal.name in self.mcmc.heatingHackProposalNames:
            heatFactor = 1.0 / (1.0 + self.mcmc.heatingHackTemperature)
            logLikeRatio *= heatFactor
            self.logPriorRatio *= heatFactor

        theSum = logLikeRatio + self.logProposalRatio + self.logPriorRatio
        if theProposal.name in ['rjComp', 'rjRMatrix']:
            theSum += self.logJacobian

        # if theProposal.name in ['ndch2_leafCompsDir']:
        #     print("%20s: %10.2f %10.2f %10.2f %10.2f" % (theProposal.name, logLikeRatio,
        #                                                  self.logPriorRatio, self.logProposalRatio, theSum), end=' ')

        # if theProposal.name in ['rjComp', 'rjRMatrix']:
        #    print "%12s: %10.2f %10.2f %10.2f %10.2f" % (theProposal.name, logLikeRatio,
        # self.logPriorRatio, self.logProposalRatio, self.logJacobian)

        # if theProposal.name == 'rMatrixLocation':
        #    print "logLikeRatio=%10.4f, logPriorRatio=%10.4f, logPosteriorRatio=%10.4f" % (
        #        logLikeRatio, self.logPriorRatio, theSum),
        # if theProposal.name == 'cmd1_allCompDir':
        #     print "logLikeRatio=%10.4f, logPriorRatio=%10.4f, logProposalRatio=%10.4f, logPosteriorRatio=%10.4f" % (
        # logLikeRatio, self.logPriorRatio, self.logProposalRatio, theSum),

        # if theProposal.name == 'polytomy':
        #    theSum += self.logJacobian
        #    self.logJacobian = 0.0
        # print "logLikeRatio = %f" % logLikeRatio
        # print "  %.2f  %.2f" % (self.logProposalRatio,
        # self.logPriorRatio)
        return theSum

    def gen(self, aProposal):
        gm = ['Chain.gen()']
        #print(gm[0], 'gen %s' % self.mcmc.gen, end=' ')


        # doAborts means that it was not a valid generation,
        # neither accepted or rejected.  Give up, by returning True.
        checkRj = False

        if checkRj:
            # Check rj stuff
            pNum = 0
            for mtNum in range(self.curTree.model.parts[pNum].nComps):
                c = self.curTree.model.parts[pNum].comps[mtNum]
                thisNNodes = 0
                for n in self.curTree.iterNodes():
                    if n.parts[pNum].compNum == c.num:
                        thisNNodes += 1
                if c.nNodes != thisNNodes:
                    gm.append(
                        "curTree  comp.nNodes=%i, but thisNNodes=%i" % (c.nNodes, thisNNodes))
                    raise P4Error(gm)
            for mtNum in range(self.propTree.model.parts[pNum].nComps):
                c = self.propTree.model.parts[pNum].comps[mtNum]
                thisNNodes = 0
                for n in self.propTree.iterNodes():
                    if n.parts[pNum].compNum == c.num:
                        thisNNodes += 1
                if c.nNodes != thisNNodes:
                    gm.append(
                        "propTree  comp.nNodes=%i, but thisNNodes=%i" % (c.nNodes, thisNNodes))
                    raise P4Error(gm)
            this_k = 0
            for mtNum in range(self.curTree.model.parts[pNum].nComps):
                c = self.curTree.model.parts[pNum].comps[mtNum]
                if c.rj_isInPool:
                    this_k += 1
            if self.curTree.model.parts[pNum].rjComp_k != this_k:
                gm.append("curTree. rjComp_k=%i, this_k=%i" %
                          (self.curTree.model.parts[pNum].rjComp_k, this_k))
                raise P4Error(gm)
            this_k = 0
            for mtNum in range(self.propTree.model.parts[pNum].nComps):
                c = self.propTree.model.parts[pNum].comps[mtNum]
                if c.rj_isInPool:
                    this_k += 1
            if self.propTree.model.parts[pNum].rjComp_k != this_k:
                gm.append("propTree rjComp_k=%i, this_k=%i" %
                          (self.propTree.model.parts[pNum].rjComp_k, this_k))
                raise P4Error(gm)

        # Same for rjRMatrix
        checkRjR = False

        if checkRjR:
            # Check rjRMatrix stuff
            for pNum in range(self.curTree.model.nParts):
                if self.curTree.model.parts[pNum].nRMatrices > 1:
                    thisK0_cur = 0
                    for mtNum in range(self.curTree.model.parts[pNum].nRMatrices):
                        c = self.curTree.model.parts[pNum].rMatrices[mtNum]
                        thisNNodes = 0
                        for n in self.curTree.iterNodesNoRoot():
                            if n.br.parts[pNum].rMatrixNum == c.num:
                                thisNNodes += 1
                        if c.nNodes != thisNNodes:
                            gm.append(
                                "curTree  rMatrix.nNodes=%i, but thisNNodes=%i" % (c.nNodes, thisNNodes))
                            raise P4Error(gm)
                        if c.nNodes:
                            thisK0_cur += 1
                    thisK0_prop = 0
                    for mtNum in range(self.propTree.model.parts[pNum].nRMatrices):
                        c = self.propTree.model.parts[pNum].rMatrices[mtNum]
                        thisNNodes = 0
                        for n in self.propTree.iterNodesNoRoot():
                            if n.br.parts[pNum].rMatrixNum == c.num:
                                thisNNodes += 1
                        if c.nNodes != thisNNodes:
                            gm.append(
                                "propTree  rMatrix.nNodes=%i, but thisNNodes=%i" % (c.nNodes, thisNNodes))
                            raise P4Error(gm)
                        if c.nNodes:
                            thisK0_prop += 1
                    if thisK0_cur != thisK0_prop:
                        gm.append("part %i, checkRjR: thisK0_cur %i, thisK0_prop %i" % (
                            pNum, thisK0_cur, thisK0_prop))
                        raise P4Error(gm)
                    this_k_cur = 0
                    for mtNum in range(self.curTree.model.parts[pNum].nRMatrices):
                        c = self.curTree.model.parts[pNum].rMatrices[mtNum]
                        if c.rj_isInPool:
                            this_k_cur += 1
                    if self.curTree.model.parts[pNum].rjRMatrix_k != this_k_cur:
                        gm.append("curTree. rjRMatrix_k=%i, this_k=%i" % (
                            self.curTree.model.parts[pNum].rjRMatrix_k, this_k_cur))
                        raise P4Error(gm)
                    this_k_prop = 0
                    for mtNum in range(self.propTree.model.parts[pNum].nRMatrices):
                        c = self.propTree.model.parts[pNum].rMatrices[mtNum]
                        if c.rj_isInPool:
                            this_k_prop += 1
                    if self.propTree.model.parts[pNum].rjRMatrix_k != this_k_prop:
                        gm.append("propTree rjRMatrix_k=%i, this_k=%i" % (
                            self.propTree.model.parts[pNum].rjRMatrix_k, this_k_prop))
                        raise P4Error(gm)

                    if this_k_cur != this_k_prop:
                        gm.append("part %i, checkRjR: this_k_cur %i, this_k_prop %i" % (
                            pNum, this_k_cur, this_k_prop))
                        raise P4Error(gm)
                    if thisK0_cur > this_k_cur:
                        gm.append("part %i, checkRjR: thisK0_cur %i, this_k_cur %i" % (
                            pNum, thisK0_cur, this_k_cur))
                        raise P4Error(gm)

        if 0:
            ret = self.verifyIdentityOfTwoTreesInChain(
                doSplitKeys=self.mcmc.constraints)
            if ret == var.DIFFERENT:
                gm.append("Trees differ at start of chain.")
                raise P4Error(gm)
            else:
                print("trees are the same -- ok")
                pass

        acceptMove = False

        if var.doMcmcSp:  # the speedy version
            if 0:
                # print "before proposal. curTree %f, %s   propTree %f, %s" % (
                # self.curTree.logLike, self.curTree.partLikes,
                # self.propTree.logLike, self.propTree.partLikes)
                if math.fabs(self.curTree.logLike - self.propTree.logLike) > 0.0001:
                    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Differs before proposal")
                if math.fabs(self.curTree.logLike - sum(self.curTree.partLikes)) > 0.0001:
                    print("7777777777777777777777777777777777777777777777777777777777 bad Cur Tree")
                if math.fabs(self.propTree.logLike - sum(self.propTree.partLikes)) > 0.0001:
                    print("8888888888888888888888888888888888888888888888888888888888 bad Prop Tree")
                # print(self.propTree.partLikes, type(self.propTree.partLikes))
                assert isinstance(self.propTree.partLikes, numpy.ndarray)
                assert isinstance(self.curTree.partLikes, numpy.ndarray)
            pRet = self.proposeSp(aProposal)
        else:
            pRet = self.propose(aProposal)

        # slow check
        if 0:
            # There should be no n.br.lenChanged or n.flag set at this point.
            if aProposal.name in ['local', 'brLen', 'eTBR', 'polytomy']:
                nnBrLenChanged = []
                nnFlags = []
                for n in self.propTree.iterNodes():
                    if n.br:
                        if n.br.lenChanged:
                            nnBrLenChanged.append(n.nodeNum)
                    if n.flag:
                        nnFlags.append(n.nodeNum)

                if nnBrLenChanged or nnFlags:
                    gm.append('gen %i, proposal %s' %
                              (self.mcmc.gen, aProposal.name))
                    gm.append('nnBrLenChanged %s, nnFlags %s' %
                              (nnBrLenChanged, nnFlags))
                    raise P4Error(gm)

                # The curTree should never have these set
                nnBrLenChanged = []
                nnFlags = []
                for n in self.curTree.iterNodes():
                    if n.br:
                        if n.br.lenChanged:
                            nnBrLenChanged.append(n.nodeNum)
                    if n.flag:
                        nnFlags.append(n.nodeNum)
                if nnBrLenChanged or nnFlags:
                    gm.append('curTree, gen %i' % self.mcmc.gen)
                    gm.append('nnBrLenChanged %s, nnFlags %s' %
                              (nnBrLenChanged, nnFlags))
                    raise P4Error(gm)

        if aProposal.name in ['local', 'eTBR'] and aProposal.doAbort:
            # local can abort because brLens were too short or too
            # long, or a constraint was violated.  There is every
            # reason to try again, up to a point.
            aProposal.nAborts[self.tempNum] += 1
            self.curTree.copyToTree(self.propTree)
            self.curTree.model.copyBQETneedsResetTo(self.propTree.model)

            safety = 0
            while 1:
                safety += 1
                if safety > 100:
                    print("Attempted %i '%s' proposals, and they all failed." % (safety, aProposal.name))
                    print("Giving up.")
                    return True  # ie failure

                if var.doMcmcSp:  # the speedy version
                    pRet = self.proposeSp(aProposal)
                else:
                    pRet = self.propose(aProposal)
                if not aProposal.doAbort:
                    break
                aProposal.nAborts[self.tempNum] += 1
                self.curTree.copyToTree(self.propTree)
                self.curTree.model.copyBQETneedsResetTo(self.propTree.model)
                # Slow check.
                if 1:
                    ret = self.verifyIdentityOfTwoTreesInChain(
                        doSplitKeys=self.mcmc.constraints)
                    if ret == var.DIFFERENT:
                        gm.append("Bad restore of propTree after doAbort.")
                        raise P4Error(gm)
                    else:
                        # print "ok"
                        pass

        if aProposal.doAbort and aProposal.name in ['compLocation', 'polytomy', 'rMatrixLocation',
                                                    'gdasrvLocation', 'rjComp', 'rjRMatrix']:
            # xxxLocation aborts if the move is impossible
            # because there are no xxx's with nNodes more than 1--
            # so its impossible, and there is no point in trying
            # again.  Best to do another move.

            # polytomy proposal aborts if it is fully resolved, and so
            # the only possible thing to do is to delete an edge, but
            # there are no suitable edges to delete.  Eg all comps
            # have nNodes < 1, or a constraint would be violated.
            # Again, there is no point in trying again, as it will
            # again be impossible.
            aProposal.nAborts[self.tempNum] += 1

            if aProposal.name in ['compLocation', 'rMatrixLocation', 'gdasrvLocation', 'rjComp', 'rjRMatrix']:
                a = self.curTree
                b = self.propTree
                a.copyToTree(b)
                a.model.parts[aProposal.pNum].copyValsTo(
                    b.model.parts[aProposal.pNum])
                b.model.setCStuff(partNum=aProposal.pNum)
                a.model.parts[aProposal.pNum].copyNNodesTo(
                    b.model.parts[aProposal.pNum])  # only one part
                a.model.parts[aProposal.pNum].copyBQETneedsResetTo(
                    b.model.parts[aProposal.pNum])
                b.setCStuff()
                # We did not do a likelihood calculation
                # pf.p4_copyCondLikes(a.cTree, b.cTree, 1) # 1 means do all
                # pf.p4_copyBigPDecks(a.cTree, b.cTree, 1) # 1 means do all
                pf.p4_copyModelPrams(a.cTree, b.cTree)

            # Slow check.
            if 1:
                ret = self.verifyIdentityOfTwoTreesInChain(
                    doSplitKeys=self.mcmc.constraints)
                if ret == var.DIFFERENT:
                    gm.append("Trees differ after doAbort.")
                    raise P4Error(gm)
                else:
                    # print "ok"
                    pass

            return True  # ie failure

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

        # print()
        # print(aProposal.name)
        # if aProposal.name == 'ndch2_internalCompsDirAlpha':
        #     print(" %f %f " %  (
        #         self.curTree.model.parts[aProposal.pNum].ndch2_internalAlpha,
        #         self.propTree.model.parts[aProposal.pNum].ndch2_internalAlpha),
        #         end=' ')
        #     print(" acceptMove = %s" % acceptMove)
        # if aProposal.name == 'ndch2_leafCompsDir':
        #     print(" acceptMove = %s" % acceptMove)

        # if aProposal.name in ['rMatrix', 'comp', 'gdasrv']:
        #    acceptMove = False

        #if self.mcmc.gen > 0 and self.mcmc.gen < 200:
        #    print("-------------- (gen %5i, %20s) acceptMove = %s" % (self.mcmc.gen, aProposal.name, acceptMove))

        aProposal.nProposals[self.tempNum] += 1
        if acceptMove:
            aProposal.accepted = True
            aProposal.nAcceptances[self.tempNum] += 1
            if aProposal.name in ['local', 'eTBR']:
                if aProposal.topologyChanged:
                    # print "zzz topologyChanged"
                    aProposal.nTopologyChangeAttempts[self.tempNum] += 1
                    aProposal.nTopologyChanges[self.tempNum] += 1
                    # aProposal.topologyChanged is (or should be) reset to zero
                    # by changeLocal() et al.
                else:
                    # print "zzz topology not changed"
                    pass
        else:
            if aProposal.name in ['local', 'eTBR']:
                if aProposal.topologyChanged:
                    aProposal.nTopologyChangeAttempts[self.tempNum] += 1
            aProposal.accepted = False

        if not aProposal.doAbort:
            if acceptMove:
                a = self.propTree
                b = self.curTree
            else:
                a = self.curTree
                b = self.propTree

            # Model values for one partition only.
            if aProposal.name in ['comp', 'compDir', 'allCompsDir',
                                  'ndch2_leafCompsDir', 'ndch2_internalCompsDir',
                                  'ndch2_leafCompsDirAlpha', 'ndch2_internalCompsDirAlpha', 'rMatrix', 
                                  'rMatrixDir', 'gdasrv', 'pInvar', 
                                  'cmd1_compDir', 'cmd1_allCompDir']:
                b.logLike = a.logLike
                pNum = aProposal.pNum
                b.partLikes[pNum] = a.partLikes[pNum]
                a.model.parts[pNum].copyValsTo(b.model.parts[pNum])
                if aProposal.name not in ['ndch2_leafCompsDir', 
                                          'ndch2_internalCompsDir', 
                                          'ndch2_leafCompsDirAlpha', 
                                          'ndch2_internalCompsDirAlpha', 
                                          'allCompsDir', 
                                          'gdasrv']:  # numpy arrays and hyperparameters
                    b.model.setCStuff(partNum=pNum)

                # Occasionally, pf.p4_setPrams() will change the bQETneedsReset
                if not (a.model.parts[pNum].bQETneedsReset == b.model.parts[pNum].bQETneedsReset).all():
                    a.model.parts[pNum].copyBQETneedsResetTo(
                        b.model.parts[pNum])  # only one part

                # These three could be faster, but they need to be re-written
                # to be part-specific.
                pf.p4_copyCondLikes(a.cTree, b.cTree, 1)  # 1 means do all
                pf.p4_copyBigPDecks(a.cTree, b.cTree, 1)
                pf.p4_copyModelPrams(a.cTree, b.cTree)

                if 0 and self.mcmc.gen == 35:  # slow check
                    previousA = a.logLike
                    a.calcLogLike(verbose=0)
                    diff = math.fabs(previousA - a.logLike)
                    if diff > 1.e-15:
                        gm.append("Chain.gen(%i).  LogLikes (a) do not match.  diff=%f (%g)" % (
                            self.mcmc.gen, diff, diff))
                        raise P4Error(gm)
                    previousB = b.logLike
                    b.calcLogLike(verbose=0)
                    diff = math.fabs(previousB - b.logLike)
                    if diff > 1.e-15:
                        gm.append("Chain.gen(%i).  LogLikes (b) do not match. diff=%f (%g)" % (
                            self.mcmc.gen, diff, diff))
                        raise P4Error(gm)
                if 0 and self.mcmc.gen == 34:
                    a.calcLogLike()
                    b.calcLogLike()

            elif aProposal.name in ['relRate']:
                b.logLike = a.logLike
                for pNum in range(self.propTree.model.nParts):  # do all parts
                    b.partLikes[pNum] = a.partLikes[pNum]
                a.model.copyValsTo(b.model)
                pf.p4_copyCondLikes(a.cTree, b.cTree, 1)  # 1 means do all
                pf.p4_copyBigPDecks(a.cTree, b.cTree, 1)  # 1 means do all
                pf.p4_copyModelPrams(a.cTree, b.cTree)

                # The propTree has, in the like calc, had
                # pf.p4_setPrams(self.propTree.cTree, -1) done to it.
                # That can put BQETneedsReset out of sync.  Easy, and
                # fast enough to simply do it to the other tree.  In
                # which case we do not need pf.p4_copyBigPDecks()
                # above, as they will all be recalculated.

                if 0 and self.mcmc.gen == 251:  # slow check
                    previousA = a.logLike
                    a.calcLogLike(verbose=0)
                    diff = math.fabs(previousA - a.logLike)
                    if diff > 1.e-15:
                        gm.append(
                            "Chain.gen().  LogLikes (a) do not match.  diff=%f (%g)" % (diff, diff))
                        raise P4Error(gm)
                    previousB = b.logLike
                    b.calcLogLike(verbose=0)
                    diff = math.fabs(previousB - b.logLike)
                    if diff > 1.e-15:
                        gm.append(
                            "Chain.gen().  LogLikes (b) do not match. diff=%f (%g)" % (diff, diff))
                        raise P4Error(gm)

            # This group is one part only.
            elif aProposal.name in ['compLocation', 'rMatrixLocation', 'gdasrvLocation']:
                b.logLike = a.logLike
                pNum = aProposal.pNum

                if 0 and self.mcmc.gen == 136:
                    print("k curTree:")
                    print(self.curTree.model.parts[pNum].bQETneedsReset)
                    print("propTree:")
                    print(self.propTree.model.parts[pNum].bQETneedsReset)

                b.partLikes[pNum] = a.partLikes[pNum]
                a.copyToTree(b)

                # Check for out-of-sync bigQET
                if acceptMove:
                    # a = self.propTree
                    # b = self.curTree
                    if b.model.parts[pNum].isHet:
                        # We are looking for combos of comp and rMatrix that
                        # were reset by the proposal and accepted
                        needsReset = b.model.parts[
                            pNum].bQETneedsReset - a.model.parts[pNum].bQETneedsReset
                        # print needsReset
                        if needsReset.any():
                            # print "Chain.gen()  fixing out-of-sync bQET after
                            # %s" % aProposal.name
                            for cNum in range(a.model.parts[pNum].nComps):
                                for rMatrixNum in range(a.model.parts[pNum].nRMatrices):
                                    if needsReset[cNum][rMatrixNum]:
                                        # print "reset cNum=%i, rMatrixNum=%i)"
                                        # % (cNum, rMatrixNum)
                                        pf.p4_resetBQET(
                                            b.model.cModel, pNum, cNum, rMatrixNum)

                if aProposal.name in ['compLocation', 'rMatrixLocation', 'gdasrvLocation']:
                    a.model.parts[pNum].copyNNodesTo(
                        b.model.parts[pNum])  # only one part
                    a.model.parts[pNum].copyBQETneedsResetTo(
                        b.model.parts[pNum])  # only one part
                b.setCStuff()

                # These next 3 could be made part-specific
                pf.p4_copyCondLikes(a.cTree, b.cTree, 1)  # 1 means do all
                pf.p4_copyBigPDecks(a.cTree, b.cTree, 1)  # 1 means do all
                pf.p4_copyModelPrams(a.cTree, b.cTree)

            # Tree topology, so all parts
            elif aProposal.name in ['local', 'eTBR', 'root3', 'brLen', 'polytomy', 'treeScale', 'allBrLens']:
                b.logLike = a.logLike
                for pNum in range(self.propTree.model.nParts):
                    b.partLikes[pNum] = a.partLikes[pNum]
                a.copyToTree(b)
                a.model.copyNNodesTo(b.model)  # all parts

                # Check for out-of-sync bigQET
                if acceptMove:
                    # a = self.propTree
                    # b = self.curTree
                    for pNum in range(a.model.nParts):
                        if b.model.parts[pNum].isHet:
                            # We are looking for combos of comp and rMatrix
                            # that were reset by the proposal and accepted
                            needsReset = b.model.parts[
                                pNum].bQETneedsReset - a.model.parts[pNum].bQETneedsReset
                            # print needsReset
                            if needsReset.any():
                                # print "Chain.gen()  fixing out-of-sync bQET
                                # after %s" % aProposal.name
                                for cNum in range(a.model.parts[pNum].nComps):
                                    for rMatrixNum in range(a.model.parts[pNum].nRMatrices):
                                        if needsReset[cNum][rMatrixNum]:
                                            # print "reset cNum=%i,
                                            # rMatrixNum=%i)" % (cNum,
                                            # rMatrixNum)
                                            pf.p4_resetBQET(
                                                b.model.cModel, pNum, cNum, rMatrixNum)

                a.model.copyBQETneedsResetTo(b.model)
                b.setCStuff()
                pf.p4_copyCondLikes(a.cTree, b.cTree, 1)  # 1 means do all
                pf.p4_copyBigPDecks(a.cTree, b.cTree, 1)  # 1 means do all
                pf.p4_copyModelPrams(a.cTree, b.cTree)

            elif aProposal.name in ['rjComp', 'rjRMatrix']:
                b.logLike = a.logLike
                pNum = aProposal.pNum
                b.partLikes[pNum] = a.partLikes[pNum]
                a.model.parts[pNum].copyValsTo(b.model.parts[pNum])
                a.copyToTree(b)
                a.model.parts[pNum].copyNNodesTo(
                    b.model.parts[pNum])  # only one part
                a.model.parts[pNum].copyBQETneedsResetTo(
                    b.model.parts[pNum])  # only one part
                b.model.setCStuff(partNum=pNum)
                b.setCStuff()

                # These three could be faster, but they need to be re-written
                # to be part-specific.
                pf.p4_copyCondLikes(a.cTree, b.cTree, 1)  # 1 means do all
                pf.p4_copyBigPDecks(a.cTree, b.cTree, 1)
                pf.p4_copyModelPrams(a.cTree, b.cTree)

            elif aProposal.name == 'cmd1_comp0Dir':
                b.model.parts[aProposal.pNum].cmd1_pi0 = a.model.parts[
                    aProposal.pNum].cmd1_pi0
            elif aProposal.name == 'cmd1_alpha':
                b.model.parts[aProposal.pNum].cmd1_alpha = a.model.parts[
                    aProposal.pNum].cmd1_alpha

            else:
                gm.append('Unlisted proposal.name = %s  Fix me.' %
                          aProposal.name)
                raise P4Error(gm)

        if 0:
            # This needs testTree.  copyToTree() is in Tree.py, and
            # copies node relations, branch lengths, and several other
            # things.
            self.curTree.copyToTree(self.testTree)
            self.curTree.model.copyValsTo(self.testTree.model)
            self.testTree.calcLogLike(verbose=False)
            myDiff = self.curTree.logLike - self.testTree.logLike
            print("diff = %f" % myDiff)
            if 1 and self.mcmc.gen == 13:
                # print "Too big!"
                print("Comparing topology stuff with Tree.verifyIdentityWith() ...")
                # python level only, false for 'doSplitKeys'
                ret = self.curTree.verifyIdentityWith(self.testTree, False)
                if ret == var.DIFFERENT:
                    print("verifyIdentityOfTwoTreesInChain() tree topology stuff differs")
                else:
                    print("topology stuff seems to be the same")
                print("Python-level: Verify model prams, with Model.verifyValsWith.")
                ret = self.curTree.model.verifyValsWith(
                    self.testTree.model)  # python level only
                if ret == var.DIFFERENT:
                    print("verifyIdentityOfTwoTreesInChain() model stuff differs")
                else:
                    print("model stuff appears to be the same")

                # cStuff.  This does model prams, tree and node stuff.
                print("about to pf.p4_verifyIdentityOfTwoTrees(self.curTree.cTree, self.testTree.cTree)")
                ret = pf.p4_verifyIdentityOfTwoTrees(
                    self.curTree.cTree, self.testTree.cTree)
                print("got ret %s" % ret)
            diffEpsi = 0.01
            if myDiff > diffEpsi or myDiff < -diffEpsi:
                raise P4Error("diff too big")

        if 0:
            if 1 and not var.doMcmcSp:
                for n in self.propTree.iterNodesNoRoot():
                    n.br.lenChanged = False
                    n.flag = False
                for n in self.curTree.iterNodesNoRoot():
                    n.br.lenChanged = False
                    n.flag = False
            isBad = False
            for n in self.propTree.iterNodesNoRoot():
                if n.br.lenChanged:
                    print("p node %2i, br.lenChanged" % n.nodeNum)
                    isBad = True
                if n.flag:
                    print("p node %2i, flag" % n.nodeNum)
                    isBad = True
            for n in self.curTree.iterNodesNoRoot():
                if n.br.lenChanged:
                    print("c node %2i, br.lenChanged" % n.nodeNum)
                    isBad = True
                if n.flag:
                    print("c node %2i, flag" % n.nodeNum)
                    isBad = True
            if isBad:
                gm.append(
                    "br.lenChanged or flag should not be set at this point.")
                raise P4Error(gm)

        # if 1:
        if (self.mcmc.gen + 1) % 100 == 0:  # every hundred gens
            ret = self.verifyIdentityOfTwoTreesInChain(
                doSplitKeys=self.mcmc.constraints)
            if ret == var.DIFFERENT:
                gm.append("gen %i" % self.mcmc.gen)
                gm.append(
                    "The two trees in the chain (ie the current tree and the proposed tree) differ.")
                gm.append("That is a programming error.")
                # if self.lastProposal: # a tuple, see a few lines below
                #    gm.append("Last proposal: %s, accepted=%s, topologyChanged=%s" % self.lastProposal)
                # else:
                #    gm.append("This appears to be the first proposal.")
                gm.append("This proposal: %s, accepted=%s, topologyChanged=%s" % (
                    aProposal.name, aProposal.accepted, aProposal.topologyChanged))
                raise P4Error(gm)
            # else:
            #    print "trees are the same at bottom of gen(), gen %i" % self.mcmc.gen
            #    print "x curTree ...."
            #    self.curTree.checkSplitKeys()
            #    print "x propTree ...."
            #    self.propTree.checkSplitKeys()

            #self.lastProposal = aProposal.name, aProposal.accepted, aProposal.topologyChanged
            # sys.exit()

        if 0:
            # no fix 397, 398
            # fix 399, curTree needed, propTree not needed.
            gNums = [136]  # random seed 3, 135 doesn't fix, 136 does
            if 1 and self.mcmc.gen in gNums:
                # self.curTree.calcLogLike()
                # self.curTree._commonCStuff(resetEmpiricalComps=False)

                # self.curTree.model.setCStuff()
                # self.curTree.setCStuff()
                # pf.p4_setPrams(self.curTree.cTree, -1) # "-1" means do all
                # parts   # sufficient to fix

                if 1:
                    self.curTree.draw(model=True)

                    print()
                    print(self.curTree.model.parts[pNum].bQETneedsReset)

                    for pNum in range(self.curTree.model.nParts):
                        for compNum in [0, 1]:
                            for rMatrixNum in [0, 1]:
                                pf.p4_resetBQET(
                                    self.curTree.model.cModel, pNum, compNum, rMatrixNum)
                    #pf.p4_resetBQET(self.curTree.model.cModel, 0, 0, 0)

                    # for n in self.curTree.iterPostOrder():
                    #    if n != self.curTree.root:
                    #        print "about to calculateBigPDecks for node %i" % n.nodeNum
                    #        pf.p4_calculateBigPDecks(n.cNode)
                        #p = n
                        # while p != self.curTree.root:
                        #    p = p.parent
                        #    p.flag = 1
                        #n.br.lenChanged = False
                    # for pNum in range(self.curTree.model.nParts):
                    #    for n in self.curTree.iterPostOrder():
                    #        if not n.isLeaf:
                    #            pf.p4_setConditionalLikelihoodsOfInteriorNodePart(n.cNode, pNum)
                    #    pf.p4_partLogLike(self.curTree.cTree, self.curTree.data.parts[pNum].cPart, pNum, 0)
                pass

            if 0 and self.mcmc.gen in gNums:
                self.propTree.calcLogLike()

        if checkRj:
            # Check rj stuff
            pNum = 0
            for mtNum in range(self.curTree.model.parts[pNum].nComps):
                c = self.curTree.model.parts[pNum].comps[mtNum]
                thisNNodes = 0
                for n in self.curTree.iterNodes():
                    if n.parts[pNum].compNum == c.num:
                        thisNNodes += 1
                if c.nNodes != thisNNodes:
                    gm.append(
                        "curTree  comp.nNodes=%i, but thisNNodes=%i" % (c.nNodes, thisNNodes))
                    raise P4Error(gm)
            for mtNum in range(self.propTree.model.parts[pNum].nComps):
                c = self.propTree.model.parts[pNum].comps[mtNum]
                thisNNodes = 0
                for n in self.propTree.iterNodes():
                    if n.parts[pNum].compNum == c.num:
                        thisNNodes += 1
                if c.nNodes != thisNNodes:
                    gm.append(
                        "propTree  comp.nNodes=%i, but thisNNodes=%i" % (c.nNodes, thisNNodes))
                    raise P4Error(gm)
            this_k = 0
            for mtNum in range(self.curTree.model.parts[pNum].nComps):
                c = self.curTree.model.parts[pNum].comps[mtNum]
                if c.rj_isInPool:
                    this_k += 1
            if self.curTree.model.parts[pNum].rjComp_k != this_k:
                gm.append("curTree. rjComp_k=%i, this_k=%i" %
                          (self.curTree.model.parts[pNum].rjComp_k, this_k))
                raise P4Error(gm)
            this_k = 0
            for mtNum in range(self.propTree.model.parts[pNum].nComps):
                c = self.propTree.model.parts[pNum].comps[mtNum]
                if c.rj_isInPool:
                    this_k += 1
            if self.propTree.model.parts[pNum].rjComp_k != this_k:
                gm.append("propTree rjComp_k=%i, this_k=%i" %
                          (self.propTree.model.parts[pNum].rjComp_k, this_k))
                raise P4Error(gm)

        if checkRjR:
            # Check rjRMatrix stuff
            for pNum in range(self.curTree.model.nParts):
                if self.curTree.model.parts[pNum].nRMatrices > 1:
                    thisK0_cur = 0
                    for mtNum in range(self.curTree.model.parts[pNum].nRMatrices):
                        c = self.curTree.model.parts[pNum].rMatrices[mtNum]
                        thisNNodes = 0
                        for n in self.curTree.iterNodesNoRoot():
                            if n.br.parts[pNum].rMatrixNum == c.num:
                                thisNNodes += 1
                        if c.nNodes != thisNNodes:
                            gm.append(
                                "curTree  rMatrix.nNodes=%i, but thisNNodes=%i" % (c.nNodes, thisNNodes))
                            raise P4Error(gm)
                        if c.nNodes:
                            thisK0_cur += 1
                    thisK0_prop = 0
                    for mtNum in range(self.propTree.model.parts[pNum].nRMatrices):
                        c = self.propTree.model.parts[pNum].rMatrices[mtNum]
                        thisNNodes = 0
                        for n in self.propTree.iterNodesNoRoot():
                            if n.br.parts[pNum].rMatrixNum == c.num:
                                thisNNodes += 1
                        if c.nNodes != thisNNodes:
                            gm.append(
                                "propTree  rMatrix.nNodes=%i, but thisNNodes=%i" % (c.nNodes, thisNNodes))
                            raise P4Error(gm)
                        if c.nNodes:
                            thisK0_prop += 1
                    if thisK0_cur != thisK0_prop:
                        gm.append("part %i, checkRjR: thisK0_cur %i, thisK0_prop %i" % (
                            pNum, thisK0_cur, thisK0_prop))
                        raise P4Error(gm)
                    this_k_cur = 0
                    for mtNum in range(self.curTree.model.parts[pNum].nRMatrices):
                        c = self.curTree.model.parts[pNum].rMatrices[mtNum]
                        if c.rj_isInPool:
                            this_k_cur += 1
                    if self.curTree.model.parts[pNum].rjRMatrix_k != this_k_cur:
                        gm.append("curTree. rjRMatrix_k=%i, this_k=%i" % (
                            self.curTree.model.parts[pNum].rjRMatrix_k, this_k_cur))
                        raise P4Error(gm)
                    this_k_prop = 0
                    for mtNum in range(self.propTree.model.parts[pNum].nRMatrices):
                        c = self.propTree.model.parts[pNum].rMatrices[mtNum]
                        if c.rj_isInPool:
                            this_k_prop += 1
                    if self.propTree.model.parts[pNum].rjRMatrix_k != this_k_prop:
                        gm.append("propTree rjRMatrix_k=%i, this_k=%i" % (
                            self.propTree.model.parts[pNum].rjRMatrix_k, this_k_prop))
                        raise P4Error(gm)

                    if this_k_cur != this_k_prop:
                        gm.append("part %i, checkRjR: this_k_cur %i, this_k_prop %i" % (
                            pNum, this_k_cur, this_k_prop))
                        raise P4Error(gm)
                    if thisK0_cur > this_k_cur:
                        gm.append("part %i, checkRjR: thisK0_cur %i, this_k_cur %i" % (
                            pNum, thisK0_cur, this_k_cur))
                        raise P4Error(gm)

    def verifyIdentityOfTwoTreesInChain(self, doSplitKeys=False):
        #gm = ['Chain.verifyIdentityOfTwoTreesInChain()']

        # print "Chain.verifyIdentityOfTwoTreesInChain().  Gen=%s" %
        # self.mcmc.gen

        # print "Python-level. Verify node relations, root, br.lens, model
        # usage, pre- and post-order."
        ret = self.curTree.verifyIdentityWith(
            self.propTree, doSplitKeys)  # python level only
        if ret == var.DIFFERENT:
            # print "verifyIdentityOfTwoTreesInChain() tree topology stuff
            # differs"
            return ret

        # print "Python-level: Verify model prams."
        ret = self.curTree.model.verifyValsWith(
            self.propTree.model)  # python level only
        if ret == var.DIFFERENT:
            # print "verifyIdentityOfTwoTreesInChain() model stuff differs"
            return ret

        # cStuff.  This does model prams, tree and node stuff.
        # print "about to pf.p4_verifyIdentityOfTwoTrees(self.curTree.cTree,
        # self.propTree.cTree)"
        ret = pf.p4_verifyIdentityOfTwoTrees(
            self.curTree.cTree, self.propTree.cTree)
        # print "ret = %s" % ret
        # if ret == var.DIFFERENT:
        #    print "verifyIdentityOfTwoTreesInChain() c stuff differs"
        return ret

    def proposeCompWithSlider(self, theProposal):
        gm = ['Chain.proposeCompWithSlider()']

        nComps = self.propTree.model.parts[theProposal.pNum].nComps
        mtNum = random.randrange(0, stop=nComps)
        mt = self.propTree.model.parts[theProposal.pNum].comps[mtNum]
        dim = self.propTree.model.parts[theProposal.pNum].dim

        # mt.val is a list, not a numpy array
        #assert isinstance(mt.val, numpy.ndarray)

        indxs = random.sample(range(dim), 2)
        currentAplusB = mt.val[indxs[0]] + mt.val[indxs[1]]
        thisMin = var.PIVEC_MIN / currentAplusB
        thisMax = 1. - thisMin

        minToMaxDiff = thisMax - thisMin
        thisTuning = theProposal.tuning

        # It is possible if A and B char states are missing that both A
        # and B values are very close to var.PIVEC_MIN, in which case
        # thisMin and thisMax will both be close to 0.5, and so the tuning
        # will be too much, requiring too many reflections.  In that case,
        # just change the tuning temporarily.
        if thisTuning > minToMaxDiff:
            thisTuning = minToMaxDiff
            # print "temporarily changing the tuning for comp proposal, to",
            # thisTuning

        x = mt.val[indxs[0]] / currentAplusB
        y = x + (thisTuning * (random.random() - 0.5))

        # reflection
        safety = -1
        while 1:
            safety += 1
            if safety > 100:
                gm.append(
                    "Did more than 100 reflections -- something is wrong.")
                raise P4Error(gm)
            if y < thisMin:
                y = thisMin + (thisMin - y)
            elif y > thisMax:
                y = thisMax - (y - thisMax)
            else:
                break
        # if safety > 1:
        #    print "comp reflections: ", safety
        mt.val[indxs[0]] = y * currentAplusB
        mt.val[indxs[1]] = currentAplusB - mt.val[indxs[0]]

        # The following normalization to a sum of 1 may not be needed.
        mySum = 0.0
        for stNum in range(dim):
            mySum += mt.val[stNum]
        for stNum in range(dim):
            mt.val[stNum] /= mySum

        self.logProposalRatio = 0.0
        # The prior here is a flat Dirichlet, ie Dirichlet(1, 1, 1, ...,
        # 1).  If it is informative, then the prior is affected.
        self.logPriorRatio = 0.0

    def proposeCompWithDirichlet(self, theProposal):
        gm = ['Chain.proposeCompWithDirichlet()']

        nComps = self.propTree.model.parts[theProposal.pNum].nComps
        mtNum = random.randrange(0, stop=nComps)
        mt = self.propTree.model.parts[theProposal.pNum].comps[mtNum]
        dim = self.propTree.model.parts[theProposal.pNum].dim

        # mt.val is a list of floats, not a numpy.ndarray
        # print(type(mt.val), type(mt.val[0]))

        # The tuning is the Dirichlet alpha.
        # print theProposal.tuning

        # This method previously used p4.func.dirichlet1, which is for lists not numpy
        # arrays.  A copy of inSeq is made, and the copy is modified and
        # returned.
        #dirichlet1(inSeq, alpha, theMin, theMax)
        #newVal = p4.func.dirichlet1(
        #    mt.val, theProposal.tuning, var.PIVEC_MIN, 1 - var.PIVEC_MIN)
        # Now it uses scipy.
        mtVal = numpy.array(mt.val)
        myProposer = scipy.stats.dirichlet(theProposal.tuning * mtVal)
        newVal = myProposer.rvs(size=1)[0]
        while  newVal.min() < var.PIVEC_MIN:
            for i in range(dim):
                if newVal[i] < var.PIVEC_MIN:
                    newVal[i] += (1.0 + random.random()) * var.PIVEC_MIN
            newVal = newVal / newVal.sum()

        # # Old way
        # newValList = newVal.tolist()
        # rangeDim = range(dim)
        # mySum = 0.0
        # for stNum in rangeDim:
        #     mySum += newValList[stNum] * theProposal.tuning
        # x = pf.gsl_sf_lngamma(mySum)
        # for stNum in rangeDim:
        #     x -= pf.gsl_sf_lngamma(newValList[stNum] * theProposal.tuning)
        # for stNum in rangeDim:
        #     x += ((newValList[stNum] * theProposal.tuning) - 1.) * \
        #         math.log(mt.val[stNum])

        # mySum = 0.0
        # for stNum in rangeDim:
        #     mySum += mt.val[stNum] * theProposal.tuning
        # y = pf.gsl_sf_lngamma(mySum)
        # for stNum in rangeDim:
        #     y -= pf.gsl_sf_lngamma(mt.val[stNum] * theProposal.tuning)
        # for stNum in rangeDim:
        #     y += ((mt.val[stNum] * theProposal.tuning) - 1.) * \
        #         math.log(newValList[stNum])
        # logProposalRatio = x - y

        # Calculate the proposal ratio
        # We can re-use myProposer to get the log pdf
        forwardLnPdf = myProposer.logpdf(newVal)
        # Another dirichlet distribution for the reverse
        spDist = scipy.stats.dirichlet(theProposal.tuning * newVal)
        reverseLnPdf = spDist.logpdf(mtVal)
        self.logProposalRatio = reverseLnPdf - forwardLnPdf

        # assert math.fabs(logProposalRatio - self.logProposalRatio) < 1.e-12 
        
        mt.val = newVal.tolist()

        # The prior here is a flat Dirichlet, ie Dirichlet(1, 1, 1, ..., 1).  If
        # it was not flat, then we would need to do some calculation here.
        self.logPriorRatio = 0.0


    def proposeRMatrixWithSlider(self, theProposal):

        # print "rMatrix proposal. the tuning is %s" % theProposal.tuning
        nRMatrices = self.propTree.model.parts[theProposal.pNum].nRMatrices
        mtNum = random.randrange(0, stop=nRMatrices)
        

        assert var.rMatrixNormalizeTo1
        mtCur = self.curTree.model.parts[theProposal.pNum].rMatrices[mtNum]
        mtProp = self.propTree.model.parts[theProposal.pNum].rMatrices[mtNum]
        if mtProp.spec == '2p':
            # For 2p, its actually a Dirichlet, not a slider.  All this is
            # stolen from MrBayes, where the default tuning is 50.  In
            # MrBayes, the "alphaDir" is a 2-item list of Dirichlet
            # parameters (not the multiplier) but they are both by default
            # 1, which makes the prior ratio 1.0 and the logPriorRatio
            # zero.

            # note that mtCur.val and mtProp.val are both ndarrays 
            # print(mtCur.val, type(mtCur.val))
            old = [0.0, 0.0]
            old[0] = mtCur.val[0] / (mtCur.val[0] + 1.0)
            old[1] = 1.0 - old[0]
            new = p4.func.dirichlet1(
                old, theProposal.tuning, var.KAPPA_MIN, var.KAPPA_MAX)
            mtProp.val[0] = new[0] / new[1]

            theSum = 0.0
            for i in range(2):
                theSum += new[i] * theProposal.tuning
            x = pf.gsl_sf_lngamma(theSum)
            for i in range(2):
                x -= pf.gsl_sf_lngamma(new[i] * theProposal.tuning)
            for i in range(2):
                x += ((new[i] * theProposal.tuning) - 1.0) * math.log(old[i])
            theSum = 0.0
            for i in range(2):
                theSum += old[i] * theProposal.tuning
            y = pf.gsl_sf_lngamma(theSum)
            for i in range(2):
                y -= pf.gsl_sf_lngamma(old[i] * theProposal.tuning)
            for i in range(2):
                y += ((old[i] * theProposal.tuning) - 1.0) * math.log(new[i])
            self.logProposalRatio = x - y

        else:  # specified, ones, eg gtr
            mt = self.propTree.model.parts[theProposal.pNum].rMatrices[mtNum]

            # mt.val is a numpy array
            assert isinstance(mt.val, numpy.ndarray)

            nRates = len(mt.val)  # eg 6 for dna gtr, not 5
            indxs = random.sample(range(nRates), 2)
            currentAplusB = mt.val[indxs[0]] + mt.val[indxs[1]]
            thisMin = var.RATE_MIN / currentAplusB
            thisMax = 1. - thisMin

            minToMaxDiff = thisMax - thisMin
            thisTuning = theProposal.tuning

            # It is possible that both A
            # and B values are very close to var.RATE_MIN, in which case
            # thisMin and thisMax will both be close to 0.5, and so the tuning
            # will be too much, requiring too many reflections.  In that case,
            # just change the tuning temporarily.
            if thisTuning > minToMaxDiff:
                thisTuning = minToMaxDiff
                # print "temporarily changing the tuning for rMatrix proposal,
                # to", thisTuning

            x = mt.val[indxs[0]] / currentAplusB
            y = x + (thisTuning * (random.random() - 0.5))

            # reflect
            safety = -1
            while 1:
                safety += 1
                if safety > 20:
                    gm.append(
                        "Did more than 20 reflections -- something is wrong.")
                    raise P4Error(gm)
                if y < thisMin:
                    y = thisMin + (thisMin - y)
                elif y > thisMax:
                    y = thisMax - (y - thisMax)
                else:
                    break
            # if safety > 1:
            #    print "rMatrix reflections: ", safety
            mt.val[indxs[0]] = y * currentAplusB
            mt.val[indxs[1]] = currentAplusB - mt.val[indxs[0]]

            mySum = 0.0
            for stNum in range(nRates):
                mySum += mt.val[stNum]
            for stNum in range(nRates):
                mt.val[stNum] /= mySum

            self.logProposalRatio = 0.0

        self.logPriorRatio = 0.0

    def proposeRMatrixDirichlet(self, theProposal):
        gm = ['Chain.proposeRMatrixDirichlet()']
        # print "rMatrix proposal. the tuning is %s" % theProposal.tuning

        assert var.rMatrixNormalizeTo1
        nRMatrices = self.propTree.model.parts[theProposal.pNum].nRMatrices
        mtNum = random.randrange(0, stop=nRMatrices)
        mtCur = self.curTree.model.parts[theProposal.pNum].rMatrices[mtNum]
        mtProp = self.propTree.model.parts[theProposal.pNum].rMatrices[mtNum]
        if mtProp.spec == '2p':

            # This is derived from MrBayes, where the default tuning is 50.  In
            # MrBayes, the "alphaDir" is a 2-item list of Dirichlet parameters
            # (not the multiplier) but they are both by default 1, which makes
            # the prior ratio 1.0 and the logPriorRatio zero.

            # note that mtCur.val and mtProp.val are both ndarrays, shape (1,) 
            # print(mtCur.val, type(mtCur.val), mtCur.val.shape)
            oldVal = numpy.array([0.0, 0.0])
            oldVal[0] = mtCur.val[0] / (mtCur.val[0] + 1.0)
            oldVal[1] = 1.0 - oldVal[0]
            myProposer = scipy.stats.dirichlet(theProposal.tuning * oldVal)
            newVal = myProposer.rvs(size=1)[0]

            safety = 0
            while 1:
                newVal = myProposer.rvs(size=1)[0]
                if newVal.min() > var.KAPPA_MIN and newVal.max() < var.KAPPA_MAX:
                    break
                safety += 1
                if safety > 100:
                    gm.append("Unable to draw a good proposal within var.KAPPA_MIN and var.KAPPA_MAX")
                    raise P4Error(gm)




            mtProp.val[0] = newVal[0] / newVal[1]
            # print(mtProp.val, type(mtProp.val), mtProp.val.shape)

            # Calculate the proposal ratio
            # We can re-use myProposer to get the log pdf
            forwardLnPdf = myProposer.logpdf(newVal)
            # Another dirichlet distribution for the reverse
            spDist = scipy.stats.dirichlet(theProposal.tuning * newVal)
            reverseLnPdf = spDist.logpdf(oldVal)
            self.logProposalRatio = reverseLnPdf - forwardLnPdf


        else:  # specified, ones, eg gtr
            mt = self.propTree.model.parts[theProposal.pNum].rMatrices[mtNum]

            # mt.val is a numpy array
            assert isinstance(mt.val, numpy.ndarray)
            myProposer = scipy.stats.dirichlet(theProposal.tuning * mt.val)

            safety = 0
            newVal = myProposer.rvs(size=1)[0]
            while newVal.min() < var.RATE_MIN or newVal.max() > var.RATE_MAX:
                for i in range(len(newVal)):
                    if newVal[i] < var.RATE_MIN:
                        newVal[i] += (1.0 + random.random()) * var.RATE_MIN
                    if newVal[i] > var.RATE_MAX:
                        newVal[i] = var.RATE_MAX - ((1.0 + random.random()) * var.RATE_MIN)
                    newVal = newVal / newVal.sum()
                    
            # Calculate the proposal ratio
            # We can re-use myProposer to get the log pdf
            forwardLnPdf = myProposer.logpdf(newVal)
            # Another dirichlet distribution for the reverse
            spDist = scipy.stats.dirichlet(theProposal.tuning * newVal)
            reverseLnPdf = spDist.logpdf(mt.val)
            self.logProposalRatio = reverseLnPdf - forwardLnPdf

            mtProp = self.propTree.model.parts[theProposal.pNum].rMatrices[mtNum]
            for i,val in enumerate(newVal):
                mtProp.val[i] = val

        self.logPriorRatio = 0.0


    def proposeGdasrv(self, theProposal):

        # This is a multiplier proposal.

        gm = ["Chain.proposeGdasrv()"]
        assert self.propTree.model.parts[theProposal.pNum].nGdasrvs == 1
        mt = self.propTree.model.parts[theProposal.pNum].gdasrvs[0]

        # We can't have alpha less than about 1.e-16, or DiscreteGamma hangs.
        # But that is moot, as var.GAMMA_SHAPE_MIN is much bigger

        # Dont' do something like the following, cuz mt.val is a property that invokes a function.
        #mt.val = newVal
        #mt.val /= theProposal.tuning

        # mt.val is a numpy.ndarray type, an array with 1 element.
        assert isinstance(mt.val, numpy.ndarray)
        oldVal = mt.val
        newVal = oldVal * math.exp(theProposal.tuning * (random.random() - 0.5))

        isGood = False
        while not isGood:
            if newVal < var.GAMMA_SHAPE_MIN:
                newVal = var.GAMMA_SHAPE_MIN * var.GAMMA_SHAPE_MIN / newVal
            elif newVal > var.GAMMA_SHAPE_MAX:
                newVal = var.GAMMA_SHAPE_MAX * var.GAMMA_SHAPE_MAX / newVal
            else:
                isGood = True

        # print(type(self.logProposalRatio), type(self.logPriorRatio), end=' ')
        self.logProposalRatio = math.log(newVal / oldVal)

        self.logPriorRatio = 0.0
        # as in proposeBrLen()
        #self.logPriorRatio = self.mcmc.tunings.parts[theProposal.pNum].gdasrvPriorLambda * float(oldVal - newVal)
        mt.val = newVal
        assert type(mt.val) == numpy.ndarray
        # print(type(self.logProposalRatio), type(self.logPriorRatio), end= ' ')
        # print(self.logProposalRatio, self.logPriorRatio)

    def proposePInvar(self, theProposal):
        mt = self.propTree.model.parts[theProposal.pNum].pInvar

        # Slider proposal
        mt.val += (random.random() - 0.5) * theProposal.tuning

        # Linear reflect
        isGood = False
        # while (mt.val < var.PINVAR_MIN) or (mt.val > var.PINVAR_MAX):
        while not isGood:
            if mt.val < var.PINVAR_MIN:
                mt.val = (var.PINVAR_MIN - mt.val) + var.PINVAR_MIN
            elif mt.val > var.PINVAR_MAX:
                mt.val = var.PINVAR_MAX - (mt.val - var.PINVAR_MAX)
            else:
                isGood = True

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0

    def proposeRelRate(self, theProposal):
        for pNum in range(self.propTree.model.nParts):
            mp = self.propTree.model.parts[pNum]
            ran = (random.random() - 0.5) * theProposal.tuning
            mp.relRate += ran
            isGood = False
            # while (mp.relRate < var.RELRATE_MIN) or (mp.relRate >
            # var.RELRATE_MAX):
            while not isGood:
                if mp.relRate < var.RELRATE_MIN:
                    mp.relRate = (
                        var.RELRATE_MIN - mp.relRate) + var.RELRATE_MIN
                elif mp.relRate > var.RELRATE_MAX:
                    mp.relRate = var.RELRATE_MAX - \
                        (mp.relRate - var.RELRATE_MAX)
                else:
                    isGood = True

        totDataLen = 0
        for p in self.propTree.data.parts:
            totDataLen += p.nChar
        fact = 0.0
        for pNum in range(self.propTree.model.nParts):
            fact += (self.propTree.model.parts[pNum].relRate *
                     self.propTree.data.parts[pNum].nChar)
        fact = float(totDataLen) / fact
        for p in self.propTree.model.parts:
            p.relRate *= fact

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0

    def proposeCompLocation(self, theProposal):
        gm = ["Chain.proposeCompLocation()"]
        mp = self.propTree.model.parts[theProposal.pNum]
        #nMT = self.propTree.model.parts[theProposal.pNum].nComps
        if mp.rjComp:
            pool = [c for c in mp.comps if c.rj_isInPool]
            # We need at least 2 comps, because one of them will be the
            # current comp for the node chosen below, and we need at least
            # one other to change to.
            if len(pool) < 2:
                theProposal.doAbort = True
                return True
        validNodeNums = [
            n for n in self.propTree.preOrder if n != var.NO_ORDER]
        validNodes = [self.propTree.nodes[n] for n in validNodeNums]
        if mp.rjComp:
            validNodes = [n for n in validNodes if (
                mp.comps[n.parts[theProposal.pNum].compNum].rj_isInPool)]
        else:
            validNodes = [n for n in validNodes if (
                mp.comps[n.parts[theProposal.pNum].compNum].nNodes > 1)]
        if not validNodes:
            theProposal.doAbort = True
            return True
        theNode = random.choice(validNodes)
        currentNum = theNode.parts[theProposal.pNum].compNum
        if mp.rjComp:
            validCompNumbers = [c.num for c in pool if c.num is not currentNum]
            if not validCompNumbers:
                theProposal.doAbort = True
                return True
        else:
            validCompNumbers = [
                c.num for c in mp.comps if c.num is not currentNum]
        #proposedNum = currentNum
        # while proposedNum == currentNum:
        #    proposedNum = random.randrange(nMT)
        proposedNum = random.choice(validCompNumbers)
        if 0 and self.mcmc.gen == 399:
            self.propTree.draw()
            print("proposeCompLocation().  node %i, before=%i, new=%s" % (theNode.nodeNum, currentNum, proposedNum))
        self.propTree.model.parts[theProposal.pNum].comps[
            currentNum].nNodes -= 1
        self.propTree.model.parts[theProposal.pNum].comps[
            proposedNum].nNodes += 1
        theNode.parts[theProposal.pNum].compNum = proposedNum
        self.logProposalRatio = 0.0
        #self.logPriorRatio = 0.0
        self.logPriorRatio = theProposal.tuning

    def proposeRMatrixLocation(self, theProposal):
        #gm = ["proposeRMatrixLocation()"]
        mp = self.propTree.model.parts[theProposal.pNum]
        #nMT = mp.nRMatrices
        if mp.rjRMatrix:
            pool = [c for c in mp.rMatrices if c.rj_isInPool]
            # We need at least 2 rMatrices, because one of them will be the
            # current rMatrix for the node chosen below, and we need at least
            # one other to change to.
            if len(pool) < 2:
                theProposal.doAbort = True
                return True
        validNodeNums = [
            n for n in self.propTree.preOrder if n != var.NO_ORDER]
        validNodeNums.remove(self.propTree.root.nodeNum)
        validNodes = [self.propTree.nodes[n] for n in validNodeNums]

        if mp.rjRMatrix:
            validNodes = [n for n in validNodes if (
                mp.rMatrices[n.br.parts[theProposal.pNum].rMatrixNum].rj_isInPool)]
        else:
            validNodes = [n for n in validNodes if (
                mp.rMatrices[n.br.parts[theProposal.pNum].rMatrixNum].nNodes > 1)]

        if not validNodes:
            theProposal.doAbort = True
            return True
        theNode = random.choice(validNodes)
        currentNum = theNode.br.parts[theProposal.pNum].rMatrixNum

        if mp.rjRMatrix:
            validRMatrixNumbers = [
                c.num for c in pool if c.num is not currentNum]
            if not validRMatrixNumbers:
                theProposal.doAbort = True
                return True
        else:
            validRMatrixNumbers = [
                c.num for c in mp.rMatrices if c.num is not currentNum]

        #proposedNum = currentNum
        # while proposedNum == currentNum:
        #    proposedNum = random.randrange(nMT)
        # print "proposeRMatrixLocation().  node %i, before=%i, new=%s" %
        # (theNode.nodeNum, currentNum, proposedNum)
        proposedNum = random.choice(validRMatrixNumbers)

        self.propTree.model.parts[theProposal.pNum].rMatrices[
            currentNum].nNodes -= 1
        self.propTree.model.parts[theProposal.pNum].rMatrices[
            proposedNum].nNodes += 1
        theNode.br.parts[theProposal.pNum].rMatrixNum = proposedNum
        self.logProposalRatio = 0.0
        #self.logPriorRatio = 0.0
        self.logPriorRatio = theProposal.tuning

    def proposeGdasrvLocation(self, theProposal):
        #gm = ["proposeGdasrvLocation()"]
        nMT = self.propTree.model.parts[theProposal.pNum].nGdasrvs
        validNodeNums = [
            n for n in self.propTree.preOrder if n != var.NO_ORDER]
        validNodeNums.remove(self.propTree.root.nodeNum)
        validNodes = [self.propTree.nodes[n] for n in validNodeNums]
        validNodes = [n for n in validNodes if (
            self.propTree.model.parts[theProposal.pNum].gdasrvs[n.br.parts[theProposal.pNum].gdasrvNum].nNodes > 1)]
        if not validNodes:
            theProposal.doAbort = True
            return True
        theNode = random.choice(validNodes)
        currentNum = theNode.br.parts[theProposal.pNum].gdasrvNum
        proposedNum = currentNum
        while proposedNum == currentNum:
            proposedNum = random.randrange(nMT)
        # print "proposeGdasrvLocation().  node %i, before=%i, new=%s" %
        # (theNode.nodeNum, currentNum, proposedNum)
        self.propTree.model.parts[theProposal.pNum].gdasrvs[
            currentNum].nNodes -= 1
        self.propTree.model.parts[theProposal.pNum].gdasrvs[
            proposedNum].nNodes += 1
        theNode.br.parts[theProposal.pNum].gdasrvNum = proposedNum
        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0



    def proposeAllCompsDir(self, theProposal):
        gm = ['Chain.proposeAllCompsDir()']
        # all the comps in one go.

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]

        assert not mpCur.ndch2, "allCompsDir proposal is not for ndch2"

        # Make proposals, accumulate log proposal ratios in the same loop
        self.logProposalRatio = 0.0
        for cNum in range(mpCur.nComps):
            mtCur = mpCur.comps[cNum]
            mtProp = mpProp.comps[cNum]
            # Result of the proposal goes into mtProp.val
            p4.func.gsl_ran_dirichlet(theProposal.tuning * mtCur.val, mtProp.val)
            while  mtProp.val.min() < var.PIVEC_MIN:
                for i in range(mpCur.dim):
                    if mtProp.val[i] < var.PIVEC_MIN:
                        mtProp.val[i] += (1.0 + random.random()) * var.PIVEC_MIN
                thisSum = mtProp.val.sum()
                mtProp.val /= thisSum

            forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, theProposal.tuning * mtCur.val, mtProp.val)
            reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, theProposal.tuning * mtProp.val, mtCur.val)
            self.logProposalRatio += reverseLnPdf - forwardLnPdf

        # prior ratio
        self.logPriorRatio = 0.0

    def proposeNdch2_leafCompsDir(self, theProposal):
        gm = ['Chain.proposeNdch2_leafCompsDir()']

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]
        assert mpCur.ndch2

        # Does this work for polytomies?  With that in mind, I iterate over
        # nodes rather than comps.

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0

        for nCur in self.curTree.iterLeavesNoRoot():
            mtNum = nCur.parts[theProposal.pNum].compNum
            mtCur = mpCur.comps[mtNum]
            mtProp = mpProp.comps[mtNum]

            # Make proposals. Result of the proposal goes into mtProp.val
            p4.func.gsl_ran_dirichlet(theProposal.tuning * mtCur.val, mtProp.val)
            while  mtProp.val.min() < var.PIVEC_MIN:
                for i in range(mpCur.dim):
                    if mtProp.val[i] < var.PIVEC_MIN:
                        mtProp.val[i] += (1.0 + random.random()) * var.PIVEC_MIN
                thisSum = mtProp.val.sum()
                mtProp.val /= thisSum

            # log proposal ratios
            forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, theProposal.tuning * mtCur.val, mtProp.val)
            reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, theProposal.tuning * mtProp.val, mtCur.val)
            self.logProposalRatio += reverseLnPdf - forwardLnPdf

            # prior ratio
            thisComp = mtCur.val
            dirPrams = mpCur.ndch2_leafAlpha * thisComp
            lnPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, dirPrams, mtCur.val)
            lnPdfProps = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, dirPrams, mtProp.val)
            self.logPriorRatio += lnPdfProps - lnPdfCurrs 
            # self.logPriorRatio = 0.0


    def proposeNdch2_internalCompsDir(self, theProposal):
        gm = ['Chain.proposeNdch2_internalCompsDir()']

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]
        assert mpCur.ndch2

        # Does this work for polytomies?  With that in mind, I iterate over
        # nodes rather than comps.

        # At the moment, the prior only looks at the global comp, and so I can
        # make the proposals, proposal ratios, and prior ratio all in one
        # iterInternals() loop.  If I were to look at local comps (ie neighbours)
        # then I would need to separate out the prior calc into its own loop.

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0

        for nCur in self.curTree.iterInternals():
            mtNum = nCur.parts[theProposal.pNum].compNum
            mtCur = mpCur.comps[mtNum]
            mtProp = mpProp.comps[mtNum]

            # Make proposals. Result of the proposal goes into mtProp.val
            p4.func.gsl_ran_dirichlet(theProposal.tuning * mtCur.val, mtProp.val)
            while  mtProp.val.min() < var.PIVEC_MIN:
                for i in range(mpCur.dim):
                    if mtProp.val[i] < var.PIVEC_MIN:
                        mtProp.val[i] += (1.0 + random.random()) * var.PIVEC_MIN
                thisSum = mtProp.val.sum()
                mtProp.val /= thisSum

            # log proposal ratios
            forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, theProposal.tuning * mtCur.val, mtProp.val)
            reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, theProposal.tuning * mtProp.val, mtCur.val)
            self.logProposalRatio += reverseLnPdf - forwardLnPdf

            # prior ratio
            dirPrams = mpCur.ndch2_internalAlpha * mpCur.ndch2_globalComp  # this is set in Mcmc.__init__()
            lnPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, dirPrams, mtCur.val)
            lnPdfProps = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, dirPrams, mtProp.val)
            self.logPriorRatio += lnPdfProps - lnPdfCurrs            




    def proposeNdch2_leafCompsDirAlpha(self, theProposal):
        # The Dirichlet hyperparameter alphaL, for leaves. 
        gm = ['Chain.proposNdch2_leafCompsDirAlpha()']

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]

        # if the max (below) gets too big, then the prior ratio can be so bad
        # that it will never be accepted.
        NDCH2_ALPHAL_MIN = 1.0 
        NDCH2_ALPHAL_MAX = 10000.   

        if 0:
            # Slider proposal
            myTuning = 50.
            oldVal = mpCur.ndch2_leafAlpha
            newVal = oldVal + (random.random() - 0.5) * myTuning

            # Linear reflect
            isGood = False
            while not isGood:
                if newVal < NDCH2_ALPHAL_MIN:
                    newVal = (NDCH2_ALPHAL_MIN - newVal) + NDCH2_ALPHAL_MIN
                elif newVal > NDCH2_ALPHAL_MAX:
                    newVal = NDCH2_ALPHAL_MAX - (newVal - NDCH2_ALPHAL_MAX)
                else:
                    isGood = True
            mpProp.ndch2_leafAlpha = newVal
            self.logProposalRatio = 0.0

        if 1: 
            # Multiplier proposal
            #myTuning = 2.0 * math.log(3.0)
            myTuning = 2.0 * math.log(1.2)
            oldVal = mpCur.ndch2_leafAlpha
            newVal = oldVal * math.exp((random.random() - 0.5) * myTuning)

            # Log reflect
            isGood = False
            while not isGood:
                if newVal < NDCH2_ALPHAL_MIN:
                    newVal = NDCH2_ALPHAL_MIN * NDCH2_ALPHAL_MIN  / newVal
                elif newVal > NDCH2_ALPHAL_MAX:
                    newVal = NDCH2_ALPHAL_MAX * NDCH2_ALPHAL_MAX / newVal
                else:
                    isGood = True
            mpProp.ndch2_leafAlpha = newVal
            self.logProposalRatio = math.log(newVal / oldVal)
        
        # Now the prior
        self.logPriorRatio = 0.0
        if 0:
            meanNeighbors = mpCur.ndch2_globalComp
            #print("proposeAllCompsDirAlphaL() meanNeighbors ", meanNeighbors)

        for nCur in self.curTree.iterNodes():
            mtNum = nCur.parts[theProposal.pNum].compNum
            mtCur = mpCur.comps[mtNum]

            # Leaf nodes use their own comp, only.
            if nCur.isLeaf:
                thisComp = mtCur.val
                #thisComp = mpCur.ndch2_globalComp
                lnPdfProp = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, newVal * thisComp, mtCur.val)
                lnPdfCur = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, oldVal * thisComp, mtCur.val)
                self.logPriorRatio += lnPdfProp - lnPdfCur

    def proposeNdch2_internalCompsDirAlpha(self, theProposal):
        # The Dirichlet hyperparameter alpha 
        gm = ['Chain.proposeNdch2_internalCompsDirAlpha()']
        # all the comps in one go.

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]

        NDCH2_ALPHAI_MIN = 1.0 
        NDCH2_ALPHAI_MAX = 2000.

        self.logProposalRatio = 0.0
        if 0:
            # Slider proposal
            myTuning = 50.
            oldVal = mpCur.ndch2_internalAlpha
            newVal = oldVal + (random.random() - 0.5) * myTuning

            # Linear reflect
            isGood = False
            while not isGood:
                if newVal < NDCH2_ALPHAI_MIN:
                    newVal = (NDCH2_ALPHAI_MIN - newVal) + NDCH2_ALPHAI_MIN
                elif newVal > NDCH2_ALPHAI_MAX:
                    newVal = NDCH2_ALPHAI_MAX - (newVal - NDCH2_ALPHAI_MAX)
                else:
                    isGood = True
            mpProp.ndch2_internalAlpha = newVal

        if 1: 
            # Multiplier proposal
            myTuning = 2.0 * math.log(3.0)
            oldVal = mpCur.ndch2_internalAlpha
            newVal = oldVal * math.exp((random.random() - 0.5) * myTuning)

            # Log reflect
            isGood = False
            while not isGood:
                if newVal < NDCH2_ALPHAI_MIN:
                    newVal = NDCH2_ALPHAI_MIN * NDCH2_ALPHAI_MIN  / newVal
                elif newVal > NDCH2_ALPHAI_MAX:
                    newVal = NDCH2_ALPHAI_MAX * NDCH2_ALPHAI_MAX / newVal
                else:
                    isGood = True
            mpProp.ndch2_internalAlpha = newVal
            self.logProposalRatio = math.log(newVal / oldVal)
        
        # Now the prior
        self.logPriorRatio = 0.0
        thisComp = mpCur.ndch2_globalComp
        for nCur in self.curTree.iterInternals():
            mtNum = nCur.parts[theProposal.pNum].compNum
            mtCur = mpCur.comps[mtNum]

            lnPdfProp = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, newVal * thisComp, mtCur.val)
            lnPdfCur = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, oldVal * thisComp, mtCur.val)
            self.logPriorRatio += lnPdfProp - lnPdfCur




