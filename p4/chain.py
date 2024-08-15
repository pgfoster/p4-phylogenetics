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

        # One more likelihood calculation, so that
        # self.verifyIdentityOfTwoTreesInChain() does not fail due to
        # part likes differing, below.

        # July 2020, comment out, see if it causes any problems.  It
        # did, with a similar problem.  This needs a better fix.  In
        # the meantime, un-comment-out.
        self.propTree.calcLogLike(verbose=0)
        
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
            eMessg = "Chain.init().  Programming error. The prop tree "
            eMessg += "should be identical to the cur tree, and it is not."
            raise P4Error(eMessg)

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

        elif theProposal.name == 'allCompsDir':
            # print "theProposal.name = allComsDir, pNum=%i" % theProposal.pNum
            self.proposeAllCompsDir(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'ndch2_leafCompsDir':
            self.proposeNdch2_leafCompsDir(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'ndch2_internalCompsDir':
            self.proposeNdch2_internalCompsDir(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'ndch2_leafCompsDirAlpha':
            # print "theProposal.name = ndch2_leafCompsDirAlpha, pNum=%i" % theProposal.pNum
            self.proposeNdch2_leafCompsDirAlpha(theProposal)
            # No likelihood calcs!

        elif theProposal.name == 'ndch2_internalCompsDirAlpha':
            # print "theProposal.name = ndch2_internalCompsDirAlpha, pNum=%i" % theProposal.pNum
            self.proposeNdch2_internalCompsDirAlpha(theProposal)
            # No likelihood calcs!

        ######################

        elif theProposal.name == 'ndrh2_leafRatesDir':
            self.proposeNdrh2_leafRatesDir(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'ndrh2_internalRatesDir':
            self.proposeNdrh2_internalRatesDir(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)

        elif theProposal.name == 'ndrh2_leafRatesDirAlpha':
            # print "theProposal.name = ndrh2_leafRatesDirAlpha, pNum=%i" % theProposal.pNum
            self.proposeNdrh2_leafRatesDirAlpha(theProposal)
            # No likelihood calcs!

        elif theProposal.name == 'ndrh2_internalRatesDirAlpha':
            # print "theProposal.name = ndrh2_internalRatesDirAlpha, pNum=%i" % theProposal.pNum
            self.proposeNdrh2_internalRatesDirAlpha(theProposal)
            # No likelihood calcs!

        ##########################

        elif theProposal.name in ['rMatrix', 'rMatrixDir', 'allRMatricesDir']:
            if theProposal.name == 'rMatrix':
                self.proposeRMatrixWithSlider(theProposal)
            elif theProposal.name == 'rMatrixDir':
                self.proposeRMatrixDirichlet(theProposal)
            else:
                self.proposeAllRMatricesDir(theProposal)

            self.propTree.model.setCStuff(partNum=theProposal.pNum)
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

        elif theProposal.name == 'allBrLens':
            self.proposeAllBrLens(theProposal)
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

        elif theProposal.name == 'root3n':
            self.proposeRoot3n(theProposal)
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
            raise P4Error("programming error.  we should not be here.  proposal %s" % theProposal.name)



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

        proposalsWithNoLikeCalcs = ['ndch2_leafCompsDirAlpha', 
                                    'ndch2_internalCompsDirAlpha',
                                    'ndch2_priorRefCompDir',
                                    'ndrh2_leafRatesDirAlpha', 
                                    'ndrh2_internalRatesDirAlpha']

        #print("...about to calculate the likelihood of the propTree.")
        self.propTree.logLike = pf.p4_treeLogLike(self.propTree.cTree, 0)    # second arg is getSiteLikes
        #self.propTree.calcLogLike()                                         # sets self.propTree.logLike
        #print("propTree logLike is", self.propTree.logLike)

        # slow check
        if 1:
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

        if theProposal.name in proposalsWithNoLikeCalcs:
            logLikeRatio = 0.0
        else:
            logLikeRatio = self.propTree.logLike - self.curTree.logLike


        # To run "without the data", which shows the effect of priors.
        #logLikeRatio = 0.0

        if self.mcmc.nChains > 1:
            if self.mcmc.swapVector:
                heatBeta = 1.0 / (1.0 + self.mcmc.chainTemps[self.tempNum])
            else:
                heatBeta = 1.0 / (1.0 + self.mcmc.chainTemp * self.tempNum)
            logLikeRatio *= heatBeta
            try:
                self.logPriorRatio *= heatBeta
            except TypeError:
                gm.append("logPriorRatio is %s, heatBeta is %s" % (self.logPriorRatio, heatBeta))
                gm.append("proposal name %s" % theProposal.name)
                raise P4Error(gm)

        theSum = logLikeRatio + self.logProposalRatio + self.logPriorRatio
        return theSum

    def proposeSp(self, theProposal):
        gm = ['Chain.proposeSp()']

        if theProposal.name == 'comp':
            # print "theProposal.name = comp, pNum=%i" % theProposal.pNum
            self.proposeCompWithSlider(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'compDir':
            # print "theProposal.name = compDir, pNum=%i" % theProposal.pNum
            self.proposeCompWithDirichlet(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

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
                           print('++++++ gen %i comp %2i %2i' % (
                               self.mcmc.gen, mtPropNum,chNum), 
                                 "%17.15f %17.15f %g" % (thisNp, thatNp, (thisNp - thatNp))) 

            # This next line is needed, and it needs to go here.  At least
            # because of p4_calculateBigPDecksPart()
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)  # "-1" means do all parts
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)


        elif theProposal.name == 'ndch2_leafCompsDir':
            # print("theProposal.name = ndch2_leafCompsDir, pNum=%i" % theProposal.pNum)
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
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

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
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

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

        elif theProposal.name == 'ndch2_priorRefCompDir':
            self.proposeNdch2_priorRefCompDir(theProposal)
            # No likelihood calcs!



        elif theProposal.name in ['rMatrix', 'rMatrixDir', 'allRMatricesDir']:
            if theProposal.name == 'rMatrix':
                self.proposeRMatrixWithSlider(theProposal)
            elif theProposal.name == 'rMatrixDir':
                self.proposeRMatrixDirichlet(theProposal)
            else:
                self.proposeAllRMatricesDir(theProposal)

            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        ######################################
        ######################################

        elif theProposal.name == 'ndrh2_leafRatesDir':
            #print("theProposal.name = ndrh2_leafRatesDir, pNum=%i" % theProposal.pNum)
            self.proposeNdrh2_leafRatesDir(theProposal)
            # This next line is needed because rMatrix is not a numpy array
            self.propTree.model.setCStuff(partNum=theProposal.pNum) 
            if 0:
                mpProp = self.propTree.model.parts[theProposal.pNum]
                for mtProp in mpProp.rMatrices:
                    mySum = sum(mtProp.val)                             #### rMatrix objects are not numpy arrays, this week
                    myDiff = math.fabs(1.0 - mySum)
                    if myDiff > 1e-15:
                        print("Chain.proposeSp() gen %i, ndrh2_leafRatesDir, x myDiff is %g" % (self.mcmc.gen, myDiff)) 
            
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)  # "-1" means do all parts
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)


        elif theProposal.name == 'ndrh2_internalRatesDir':
            #print("theProposal.name = ndrh2_internalRatesDir, pNum=%i" % theProposal.pNum)
            self.proposeNdrh2_internalRatesDir(theProposal)
            # This next line is needed because rMatrix is not a numpy array 
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            if 0:
                mpProp = self.propTree.model.parts[theProposal.pNum]
                for mtProp in mpProp.rMatrices:
                    mySum = sum(mtProp.val)
                    myDiff = math.fabs(1.0 - mySum)
                    if myDiff > 1e-15:
                        print("Chain.proposeSp() gen %i, ndrh2_internalRatesDir, x myDiff is %g" % (self.mcmc.gen, myDiff)) 
            
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)  # "-1" means do all parts
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'ndrh2_leafRatesDirAlpha':
            # print "theProposal.name = ndrh2_leafRatesDirAlpha, pNum=%i" % theProposal.pNum
            self.proposeNdrh2_leafRatesDirAlpha(theProposal)
            # No likelihood calcs!

        elif theProposal.name == 'ndrh2_internalRatesDirAlpha':
            # print "theProposal.name = ndrh2_internalRatesDirAlpha, pNum=%i" % theProposal.pNum
            self.proposeNdrh2_internalRatesDirAlpha(theProposal)
            # No likelihood calcs!

        elif theProposal.name == 'ndrh2_priorRefRMatrixDir':
            self.proposeNdrh2_priorRefRMatrixDir(theProposal)
            # No likelihood calcs!


        ######################################
        ######################################

        elif theProposal.name == 'gdasrv':
            self.proposeGdasrv(theProposal)
            # This next line is not needed because gdasrv is a numpy array
            #self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf  or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'pInvar':
            self.proposePInvar(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

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
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

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
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

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
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, theProposal.pNum)

            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'local':
            if self.propTree.root.getNChildren() == 2:
                rooter = self.propTree.attachRooter()
                self.propTree.makeSplitKeys()
                #print("rooter br is %s, rooter.nodeNum=%i" % (rooter.br, rooter.nodeNum))
                #print("Starting local proposal ...")
                self.proposeLocal(theProposal)
                #print("  ... finished local proposal")

                if rooter == self.propTree.root:
                    #self.propTree.draw()
                    assert rooter.leftChild.leftChild.sibling
                    assert not rooter.leftChild.leftChild.sibling.sibling
                    self.propTree.reRoot(rooter.leftChild)

                if rooter.br and rooter.br.parts:
                    for n in self.propTree.root.iterChildren():
                        if not n.br.parts:
                            n.br.parts = rooter.br.parts
                            rooter.br.parts = []
                            break

                if rooter.br.parts:
                    print("xx rooter is node %i, br=%s, parts=%s" % (rooter.nodeNum, rooter.br, rooter.br.parts))
                    for n in self.propTree.iterNodesNoRoot():
                        if n.br.parts:
                            pass
                        else:
                            print("xx Node %i br (%s) has no parts." % (n.nodeNum, n.br))
                    raise P4Error(gm)


                self.propTree.reRoot(rooter.parent)
                self.propTree.detachRooter()

            else:
                self.proposeLocal(theProposal)
            if theProposal.doAbort:
                return 0.0
            
            if not self.propTree.preAndPostOrderAreValid:
                self.propTree.setPreAndPostOrder()

            # Are we violating constraints?
            if self.mcmc.constraints and theProposal.topologyChanged:
                self.propTree.makeSplitKeys(makeNodeForSplitKeyDict=True)  # makeNodeForSplitKeyDict=True is default
                ret = self.mcmc.constraints.areConsistentWithTree(self.propTree)
                if not ret:
                    theProposal.doAbort = True
                    return 0.0
                if self.mcmc.constraints.rooting:
                    ret = self.mcmc.constraints.areConsistentWithTreeRoot(self.propTree)
                    if not ret:
                        theProposal.doAbort = True
                        return 0.0


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
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    if n.flag:
                        for pNum in range(self.propTree.model.nParts):
                            pf.p4_setConditionalLikelihoodsOfInternalNodePart(
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
                    if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                        if n.flag:
                            pf.p4_setConditionalLikelihoodsOfInternalNodePart(
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
                    if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                        # Normally we would check whether n.flag is set
                        # but here they all need to recalc cond likes
                        pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, pNum)

                pf.p4_partLogLike(
                    self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)
            #for n in self.propTree.iterInternalsNoRoot():
            #    n.flag = 0
            #self.propTree.root.flag = 0

        elif theProposal.name == 'eTBR':
            # doAbort is set, below, if constraints are violated

            if self.propTree.root.getNChildren() == 2:
                    rooter = self.propTree.attachRooter()
                    # self.propTree.makeSplitKeys()
                    #print("rooter br is %s, rooter.nodeNum=%i" % (rooter.br, rooter.nodeNum))
                    #print("Starting eTBR proposal ...")
                    self.proposeETBR_Blaise(theProposal)
                    #print("  ... finished eTBR proposal")

                    if rooter == self.propTree.root:
                        #self.propTree.draw()
                        assert rooter.leftChild.leftChild.sibling
                        assert not rooter.leftChild.leftChild.sibling.sibling
                        self.propTree.reRoot(rooter.leftChild)

                    if rooter.br and rooter.br.parts:
                        for n in self.propTree.root.iterChildren():
                            if not n.br.parts:
                                n.br.parts = rooter.br.parts
                                rooter.br.parts = []
                                break

                    if rooter.br.parts:
                        print("xx rooter is node %i, br=%s, parts=%s" % (rooter.nodeNum, rooter.br, rooter.br.parts))
                        for n in self.propTree.iterNodesNoRoot():
                            if n.br.parts:
                                pass
                            else:
                                print("xx Node %i br (%s) has no parts." % (n.nodeNum, n.br))
                        raise P4Error(gm)
                    

                    self.propTree.reRoot(rooter.parent)
                    self.propTree.detachRooter()
                    
            
            else:
                self.proposeETBR_Blaise(theProposal)

            
            if not self.propTree.preAndPostOrderAreValid:
                self.propTree.setPreAndPostOrder()

            # Are we violating constraints?
            if self.mcmc.constraints and theProposal.topologyChanged:
                self.propTree.makeSplitKeys(makeNodeForSplitKeyDict=True)  # makeNodeForSplitKeyDict=True is default
                ret = self.mcmc.constraints.areConsistentWithTree(self.propTree)
                if not ret:
                    #print("proposeSp() after eTBR, doAbort set due to violated constraints")
                    theProposal.doAbort = True
                    return 0.0
                if self.mcmc.constraints.rooting:
                    ret = self.mcmc.constraints.areConsistentWithTreeRoot(self.propTree)
                    if not ret:
                        #print("proposeSp() after eTBR, doAbort set due to violated constraint for root")
                        theProposal.doAbort = True
                        return 0.0

            # Bugfix, moved from proposeETBR_Blaise() to here, *after* checking for violated constraints.
            if theProposal.topologyChanged:
                for pNum in range(self.propTree.model.nParts):
                    # print("\n")
                    # print(self.propTree.model.parts[pNum].bQETneedsReset)
                    for cNum in range(self.propTree.model.parts[pNum].nComps):
                        for rNum in range(self.propTree.model.parts[pNum].nRMatrices):
                            if self.propTree.model.parts[pNum].bQETneedsReset[cNum][rNum]:
                                pf.p4_resetBQET(self.propTree.model.cModel, pNum, cNum, rNum)
                    # print("after p4_resetBQET, ...")
                    # print(self.propTree.model.parts[pNum].bQETneedsReset)


            self.propTree.setCStuff()

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


            # Recalculate condLikes for only the flagged nodes, in post
            # order down to the root.
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                    if n.flag:
                        for pNum in range(self.propTree.model.nParts):
                            pf.p4_setConditionalLikelihoodsOfInternalNodePart(
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
                        if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                            pf.p4_setConditionalLikelihoodsOfInternalNodePart(
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
                        if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                            # print "node %i" % n.nodeNum
                            if n.flag:
                                pf.p4_setConditionalLikelihoodsOfInternalNodePart(
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
                    if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                        pf.p4_setConditionalLikelihoodsOfInternalNodePart(
                            n.cNode, pNum)
                pf.p4_partLogLike(
                    self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)

        elif theProposal.name == 'root3n':
            self.proposeRoot3n(theProposal)
            if not self.propTree.preAndPostOrderAreValid:
                self.propTree.setPreAndPostOrder()
            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, -1)
            for pNum in range(self.propTree.model.nParts):
                for n in self.propTree.iterPostOrder():
                    if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                        pf.p4_setConditionalLikelihoodsOfInternalNodePart(
                            n.cNode, pNum)
                pf.p4_partLogLike(
                    self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)

        elif theProposal.name == 'root2':
            self.proposeRoot2(theProposal)

            # The root may have changed
            if not self.propTree.preAndPostOrderAreValid:
                self.propTree.setPreAndPostOrder()
            self.propTree.makeSplitKeys()

            if self.mcmc.constraints:   # root constraint?
                # It does not make a lot of sense to check constraints for this move.
                pass


            self.propTree.setCStuff()
            pf.p4_setPrams(self.propTree.cTree, -1)
            for pNum in range(self.propTree.model.nParts):
                for n in self.propTree.iterPostOrder():
                    if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                        pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, pNum)
                pf.p4_partLogLike(self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)

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
                    if not n.isLeaf or (n == self.propTree.root and n.isLeaf):  # possibly leaf root
                        pf.p4_setConditionalLikelihoodsOfInternalNodePart(
                            n.cNode, pNum)
                pf.p4_partLogLike(
                    self.propTree.cTree, self.propTree.data.parts[pNum].cPart, pNum, 0)

        else:
            gm.append('Unlisted proposal.name=%s  Fix me.' % theProposal.name)
            raise P4Error(gm)

        if theProposal.doAbort:
            raise P4Error("programming error.  we should not be here.  proposal %s" % theProposal.name)

        proposalsWithNoLikeCalcs = ['ndch2_leafCompsDirAlpha', 'ndch2_internalCompsDirAlpha',
                                    'ndrh2_leafRatesDirAlpha', 'ndrh2_internalRatesDirAlpha', 
                                    'ndch2_priorRefCompDir', 'ndrh2_priorRefRMatrixDir']


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

        if self.mcmc.ssBeta is not None:
            logLikeRatio *= self.mcmc.ssBeta


        # Heating via Heating hack, 
        if self.mcmc.doHeatingHack:
            assert self.mcmc.swapVector
            if self.mcmc.heatingHackRamper:
                hHTemp = self.mcmc.originalHeatingHackTemperature
                myG = self.mcmc.gen
                if myG <= self.mcmc.heatingHackRamper[0]:
                    pass
                elif myG >= self.mcmc.heatingHackRamper[1]:
                    hHTemp = 0.0
                else:
                    rdiff = self.mcmc.heatingHackRamper[1] - self.mcmc.heatingHackRamper[0]
                    factor = (self.mcmc.heatingHackRamper[1] - myG) / rdiff 
                    hHTemp = self.mcmc.originalHeatingHackTemperature * factor
                    #print("hHTemp. factor %.3f, new hHTemp %.3f" % (factor, hHTemp))
                self.mcmc.heatingHackTemperature = hHTemp
            
            if self.mcmc.nChains > 1:
                heatFactor = 1.0 / (1.0 + self.mcmc.heatingHackTemperature + self.mcmc.chainTemps[self.tempNum])
            else:
                heatFactor = 1.0 / (1.0 + self.mcmc.heatingHackTemperature)
            logLikeRatio *= heatFactor
            self.logPriorRatio *= heatFactor

        # or simulated tempering, .. 
        elif self.mcmc.simTemp:
            myTmpO = self.mcmc.simTemp.temps[self.mcmc.simTemp.curTemp]  # SimTempTemp object
            myTmp = myTmpO.temp
            heatFactor = 1.0 / (1.0 + myTmp)
            logLikeRatio *= heatFactor
            self.logPriorRatio *= heatFactor

        else:
            #  ... or MCMCMC
            # See Altekar et al 2004.  Parallel Metropolis coupled Markov chain MonteCarlo for 
            # Bayesian phylogenetic inference Bioinformatics 20:407 https://doi.org/10.1093/bioinformatics/btg427
            # Both the likelihood ratio and the prior ratio are powered.  Not the proposal ratio.
            if self.mcmc.nChains > 1:
                if self.mcmc.swapVector:
                    heatBeta = 1.0 / (1.0 + self.mcmc.chainTemps[self.tempNum])
                else:
                    heatBeta = 1.0 / (1.0 + self.mcmc.chainTemp * self.tempNum)
                logLikeRatio *= heatBeta
                try:
                    self.logPriorRatio *= heatBeta
                except TypeError:
                    gm.append("logPriorRatio is %s, heatBeta is %s" % (self.logPriorRatio, heatBeta))
                    gm.append("proposal name %s" % theProposal.name)
                    raise P4Error(gm)

        theSum = logLikeRatio + self.logProposalRatio + self.logPriorRatio
        # if theProposal.name in ['ndrh2_leafRatesDirAlpha']:
        #     print("%20s: %10.2f %10.2f %10.2f %10.2f" % (theProposal.name, logLikeRatio,
        #                                                  self.logPriorRatio, self.logProposalRatio, theSum), end=' ')

        # if theProposal.name == 'rMatrixLocation':
        #    print "logLikeRatio=%10.4f, logPriorRatio=%10.4f, logPosteriorRatio=%10.4f" % (
        #        logLikeRatio, self.logPriorRatio, theSum),
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
        #sys.stdout.flush()


        # doAborts means that it was not a valid generation,
        # neither accepted or rejected.  Give up, by returning True.

        if 0:
            ret = self.verifyIdentityOfTwoTreesInChain(
                doSplitKeys=self.mcmc.constraints)
            if ret == var.DIFFERENT:
                gm.append("Trees differ at start of chain.gen(), before making a proposal")
                raise P4Error(gm)
            else:
                # print("trees are the same -- ok")
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
        # print("finished proposal")

        # slow check
        if 0:
            print("doing slow check ...")
            # There should be no n.br.lenChanged or n.flag set at this point.
            if aProposal.name in ['local', 'brLen', 'eTBR', 'polytomy', 'allBrLens']:
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

        if aProposal.doAbort:
            if aProposal.name in ['local', 'eTBR']:
                # local can abort because brLens were too short or too
                # long, or a constraint was violated.  There is every
                # reason to try again, up to a point.
                safety = 0
                while 1:
                    safety += 1

                    aProposal.nAborts[self.tempNum] += 1
                    self.curTree.copyToTree(self.propTree)
                    self.curTree.model.copyBQETneedsResetTo(self.propTree.model)

                    if safety >= 100:
                        print("Attempted %i '%s' proposals, and they all failed." % (safety, aProposal.name))
                        print("Giving up.")
                        return True  # ie failure

                    # Slow check after restore
                    if 0:
                        ret = self.verifyIdentityOfTwoTreesInChain(
                            doSplitKeys=self.mcmc.constraints)
                        if ret == var.DIFFERENT:
                            gm.append(f"{aProposal.name} - Bad restore of propTree after doAbort.")
                            raise P4Error(gm)
                        else:
                            # print "ok"
                            pass

                    # Try again
                    if var.doMcmcSp:  # the speedy version
                        # print(gm[0], f"proposing {aProposal.name} again, after doAbort True")
                        pRet = self.proposeSp(aProposal)
                    else:
                        pRet = self.propose(aProposal)
                    if not aProposal.doAbort:
                        break
                    


            elif aProposal.name in ['compLocation', 'polytomy', 'rMatrixLocation',
                                                        'gdasrvLocation']:
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

                if aProposal.name in ['compLocation', 'rMatrixLocation', 'gdasrvLocation']:
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

                    # slow check
                    if 0:
                        ret = self.verifyIdentityOfTwoTreesInChain(doSplitKeys=self.mcmc.constraints)
                        if ret == var.DIFFERENT:
                            gm.append(f"{aProposal.name} Trees differ after doAbort.")
                            raise P4Error(gm)
                        else:
                            # print "ok"
                            pass

                    # Giving up on this proposal, returning True to
                    # the mcmc.run(), the calling method where it is assigned to the variable  "failure".
                    return True

                else:
                    gm.append(f"doAbort was set for proposal {aProposal.name}")
                    gm.append("but it was not handled.  This is a programming error.  Fix me.")
                    raise P4Error(gm)

        assert not aProposal.doAbort, f"proposal {aProposal.name} doAbort set.  Fix me." 
        # if aProposal.name == 'ndrh2_leafRatesDir':
        #     print("pRet = %10.6f" % pRet, end=' ')


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


        if 0: 
            # if aProposal.name == 'root2': #and self.mcmc.gen >= 0 and self.mcmc.gen < 1000:
            print("-------------- (gen %5i, %30s) acceptMove = %6s" % (self.mcmc.gen, aProposal.name, acceptMove), end='\n')

        aProposal.nProposals[self.tempNum] += 1
        aProposal.tnNSamples[self.tempNum] += 1
        if acceptMove:
            aProposal.accepted = True
            aProposal.nAcceptances[self.tempNum] += 1
            aProposal.tnNAccepts[self.tempNum] += 1
            if aProposal.name in ['local', 'eTBR']:
                if aProposal.topologyChanged:
                    # print(gm[0], "zzz topologyChanged")
                    aProposal.nTopologyChangeAttempts[self.tempNum] += 1
                    aProposal.nTopologyChanges[self.tempNum] += 1
                    # aProposal.topologyChanged is (or should be) reset to zero
                    # by changeLocal() et al.
                else:
                    # print(gm[0], "zzz topology not changed")
                    pass
        else:
            if aProposal.name in ['local', 'eTBR']:
                if aProposal.topologyChanged:
                    aProposal.nTopologyChangeAttempts[self.tempNum] += 1
                    # but no nTopologyChanges
            aProposal.accepted = False


        if acceptMove:
            a = self.propTree
            b = self.curTree
        else:
            a = self.curTree
            b = self.propTree

        # Model values for one partition only.
        if aProposal.name in ['comp', 'compDir', 'allCompsDir',
                              'ndch2_leafCompsDir', 'ndch2_internalCompsDir',
                              'ndch2_leafCompsDirAlpha', 'ndch2_internalCompsDirAlpha', 
                              'ndch2_priorRefCompDir',
                              'rMatrix', 'rMatrixDir', 'allRMatricesDir', 
                              'ndrh2_leafRatesDir', 'ndrh2_internalRatesDir',
                              'ndrh2_leafRatesDirAlpha', 'ndrh2_internalRatesDirAlpha', 
                              'ndrh2_priorRefRMatrixDir',
                              'gdasrv', 'pInvar']:
            b.logLike = a.logLike
            pNum = aProposal.pNum
            b.partLikes[pNum] = a.partLikes[pNum]
            a.model.parts[pNum].copyValsTo(b.model.parts[pNum])
            if aProposal.name not in ['ndch2_leafCompsDir', 
                                      'ndch2_internalCompsDir', 
                                      'ndch2_leafCompsDirAlpha', 
                                      'ndch2_internalCompsDirAlpha', 
                                      'ndch2_priorRefCompDir',
                                      'ndrh2_leafRatesDirAlpha', 
                                      'ndrh2_internalRatesDirAlpha',
                                      'ndrh2_priorRefRMatrixDir',
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
        elif aProposal.name in ['local', 'eTBR', 'root3', 'root3n', 'root2',
                                'brLen', 'polytomy', 'treeScale', 
                                'allBrLens']:
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
                ret = self.curTree.verifyIdentityWith(self.testTree, doSplitKeys=False)
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
        #print("curTree.logLike %12.6f" % self.curTree.logLike)

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
            # Check ndch2 and ndrh2; make sure all comps or rMatrices
            # are there, and that there is no more than one node for
            # each comp or rMatrix.

            # Note that checking the length of the rNums or cNums
            # against self.curTree.nodes will not work in the case of
            # bi-rooted trees, which have an extra node (not seen in
            # tree traversal methods).

            if self.curTree.model.parts[0].ndrh2:
                rNums = []
                for n in self.curTree.iterNodesNoRoot():
                    rNum = n.br.parts[0].rMatrixNum
                    rNums.append(rNum)
                for rNum in rNums:
                    assert rNums.count(rNum) == 1

            if self.curTree.model.parts[0].ndch2:
                cNums = []
                for n in self.curTree.iterNodes():
                    cNum = n.parts[0].compNum
                    cNums.append(cNum)
                for cNum in cNums:
                    assert cNums.count(cNum) == 1
               
            
        #if self.mcmc.gen % 1000 == 0:
        #    print(self.curTree.model.parts[0].ndch2_priorRefComp)



        if 0:
            # Leftover debugging stuff.
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
                    #        if not n.isLeaf or (n == propTree.root and n.isLeaf):  # possibly leaf root
                    #            pf.p4_setConditionalLikelihoodsOfInternalNodePart(n.cNode, pNum)
                    #    pf.p4_partLogLike(self.curTree.cTree, self.curTree.data.parts[pNum].cPart, pNum, 0)
                pass

            if 0 and self.mcmc.gen in gNums:
                self.propTree.calcLogLike()



    def verifyIdentityOfTwoTreesInChain(self, doSplitKeys=False):
        #gm = ['Chain.verifyIdentityOfTwoTreesInChain()']

        # print "Chain.verifyIdentityOfTwoTreesInChain().  Gen=%s" %
        # self.mcmc.gen

        # print "Python-level. Verify node relations, root, br.lens, model
        # usage, pre- and post-order."
        ret = self.curTree.verifyIdentityWith(self.propTree, doSplitKeys=doSplitKeys)  # python level only
        if ret == var.DIFFERENT:
            # print("verifyIdentityOfTwoTreesInChain() tree topology stuff differs")
            return ret

        # print "Python-level: Verify model prams."
        ret = self.curTree.model.verifyValsWith(
            self.propTree.model)  # python level only
        if ret == var.DIFFERENT:
            # print "verifyIdentityOfTwoTreesInChain() model stuff differs"
            return ret

        # cStuff.  This does model prams, tree and node stuff.
        # print("about to pf.p4_verifyIdentityOfTwoTrees(self.curTree.cTree,self.propTree.cTree)")
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
        # print theProposal.tuning[self.tempNum]

        # This method previously used p4.func.dirichlet1, which is for lists not numpy
        # arrays.  A copy of inSeq is made, and the copy is modified and
        # returned.
        #dirichlet1(inSeq, alpha, theMin, theMax)
        #newVal = p4.func.dirichlet1(
        #    mt.val, theProposal.tuning, var.PIVEC_MIN, 1 - var.PIVEC_MIN)
        # Now it uses scipy.
        mtVal = numpy.array(mt.val)
        myProposer = scipy.stats.dirichlet(theProposal.tuning[self.tempNum] * mtVal)
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
        spDist = scipy.stats.dirichlet(theProposal.tuning[self.tempNum] * newVal)
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
            myProposer = scipy.stats.dirichlet(theProposal.tuning[self.tempNum] * oldVal)
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
            spDist = scipy.stats.dirichlet(theProposal.tuning[self.tempNum] * newVal)
            reverseLnPdf = spDist.logpdf(oldVal)
            self.logProposalRatio = reverseLnPdf - forwardLnPdf


        else:  # specified, ones, eg gtr
            mt = self.propTree.model.parts[theProposal.pNum].rMatrices[mtNum]

            # mt.val is a numpy array
            assert isinstance(mt.val, numpy.ndarray)
            myProposer = scipy.stats.dirichlet(theProposal.tuning[self.tempNum] * mt.val)

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
            spDist = scipy.stats.dirichlet(theProposal.tuning[self.tempNum] * newVal)
            reverseLnPdf = spDist.logpdf(mt.val)
            self.logProposalRatio = reverseLnPdf - forwardLnPdf

            mtProp = self.propTree.model.parts[theProposal.pNum].rMatrices[mtNum]
            for i,val in enumerate(newVal):
                mtProp.val[i] = val

        self.logPriorRatio = 0.0


    def proposeGdasrv(self, theProposal):
        # This is a "multiplier" proposal.

        gm = ["Chain.proposeGdasrv()"]
        assert self.propTree.model.parts[theProposal.pNum].nGdasrvs == 1
        mt = self.propTree.model.parts[theProposal.pNum].gdasrvs[0]

        # We can't have alpha less than about 1.e-16, or DiscreteGamma hangs.
        # But that is moot, as var.GAMMA_SHAPE_MIN is much bigger

        # Don't do something like the following, cuz mt.val is a property 
        # that invokes a function, and we don't want to do that until later 
        # (below).  See the Gdasrv class.
        #mt.val = newVal
        #mt.val /= theProposal.tuning

        # mt.val is a numpy.ndarray type, an array with 1 element.
        assert isinstance(mt.val, numpy.ndarray)
        oldVal = mt.val
        newVal = oldVal * math.exp(theProposal.tuning[self.tempNum] * (random.random() - 0.5))

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

        # prior
        if theProposal.prior == None:
            # Default exponential, lambda=1.0
            self.logPriorRatio = oldVal - newVal
        elif isinstance(theProposal.prior, float):
            # Exponential with specified lambda
            self.logPriorRatio = theProposal.prior * (oldVal - newVal)
        elif theProposal.prior == 'uniform':
            self.logPriorRatio = 0.0
        elif hasattr(theProposal.prior, "logpdf"):  # a scipy.stats distribution
            oldPr = theProposal.prior.logpdf(oldVal)
            newPr = theProposal.prior.logpdf(newVal)
            self.logPriorRatio = newPr - oldPr
        else:
            gm.append(f"prior '{theProposal.prior}' is not understood.")
            raise P4Error(gm)

        # Setting mt.val triggers recalculation of discrete gamma values
        mt.val = newVal
        assert type(mt.val) == numpy.ndarray
        # print(type(self.logProposalRatio), type(self.logPriorRatio), end= ' ')
        # print(self.logProposalRatio, self.logPriorRatio)

    def proposePInvar(self, theProposal):
        mt = self.propTree.model.parts[theProposal.pNum].pInvar

        # Slider proposal
        mt.val += (random.random() - 0.5) * theProposal.tuning[self.tempNum]

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
            ran = (random.random() - 0.5) * theProposal.tuning[self.tempNum]
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
        validNodeNums = [n for n in self.propTree.preOrder if n != var.NO_ORDER]
        validNodes = [self.propTree.nodes[n] for n in validNodeNums]
        validNodes = [n for n in validNodes if (mp.comps[n.parts[theProposal.pNum].compNum].nNodes > 1)]
        if not validNodes:
            theProposal.doAbort = True
            return True
        theNode = random.choice(validNodes)
        currentNum = theNode.parts[theProposal.pNum].compNum
        validCompNumbers = [c.num for c in mp.comps if c.num is not currentNum]
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
        if theProposal.tuning:
            self.logPriorRatio = theProposal.tuning[self.tempNum]  # hack alert!
        else:
            self.logPriorRatio = 0.0

    def proposeRMatrixLocation(self, theProposal):
        #gm = ["proposeRMatrixLocation()"]
        mp = self.propTree.model.parts[theProposal.pNum]
        #nMT = mp.nRMatrices
        validNodeNums = [
            n for n in self.propTree.preOrder if n != var.NO_ORDER]
        validNodeNums.remove(self.propTree.root.nodeNum)
        validNodes = [self.propTree.nodes[n] for n in validNodeNums]

        validNodes = [n for n in validNodes if (
            mp.rMatrices[n.br.parts[theProposal.pNum].rMatrixNum].nNodes > 1)]

        if not validNodes:
            theProposal.doAbort = True
            return True
        theNode = random.choice(validNodes)
        currentNum = theNode.br.parts[theProposal.pNum].rMatrixNum

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
        if theProposal.tuning:
            self.logPriorRatio = theProposal.tuning[self.tempNum]  # hack alert!
        else:
            self.logPriorRatio = 0.0

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
            p4.func.gsl_ran_dirichlet(theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            while  mtProp.val.min() < var.PIVEC_MIN:
                for i in range(mpCur.dim):
                    if mtProp.val[i] < var.PIVEC_MIN:
                        mtProp.val[i] += (1.0 + random.random()) * var.PIVEC_MIN
                thisSum = mtProp.val.sum()
                mtProp.val /= thisSum

            forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(
                mpCur.dim, theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(
                mpCur.dim, theProposal.tuning[self.tempNum] * mtProp.val, mtCur.val)
            self.logProposalRatio += reverseLnPdf - forwardLnPdf

        # prior ratio
        self.logPriorRatio = 0.0

    def proposeAllRMatricesDir(self, theProposal):
        gm = ['Chain.proposeAllRMatricesDir()']
        # all the rMatrices in one go.

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]

        #assert not mpCur.ndch2, "allRMatricesDir proposal is not for ndch2"  # Why not?

        # Make proposals, accumulate log proposal ratios in the same loop
        self.logProposalRatio = 0.0
        for mtNum in range(mpCur.nRMatrices):
            mtCur = mpCur.rMatrices[mtNum]
            mtProp = mpProp.rMatrices[mtNum]
            assert mtProp.spec != '2p', "proposeAllRMatricesDir is not set up for 2p models yet.  To do."
            # Result of the proposal goes into mtProp.val
            assert isinstance(mtProp.val, numpy.ndarray)
            p4.func.gsl_ran_dirichlet(theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            #assert isinstance(mtProp.val, numpy.ndarray)
            while  mtProp.val.min() < var.RATE_MIN or mtProp.val.max() > var.RATE_MAX:
                for i in range(len(mtProp.val)):
                    if mtProp.val[i] < var.RATE_MIN:
                        mtProp.val[i] += (1.0 + random.random()) * var.RATE_MIN
                    if mtProp.val[i] > var.RATE_MAX:
                        mtProp.val[i] = var.RATE_MAX - ((1.0 + random.random()) * var.RATE_MIN)
                thisSum = mtProp.val.sum()
                mtProp.val /= thisSum

            forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(len(mtCur.val), theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(len(mtCur.val), theProposal.tuning[self.tempNum] * mtProp.val, mtCur.val)
            self.logProposalRatio += reverseLnPdf - forwardLnPdf

        # prior ratio
        self.logPriorRatio = 0.0



    ####################################  NDRH2
    ####################################

    def proposeNdrh2_leafRatesDir(self, theProposal):
        gm = ['Chain.proposeNdrh2_leafRatesDir()']

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]
        assert mpCur.ndrh2

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0

        ratesLen = int(((mpCur.dim * mpCur.dim) - mpCur.dim) / 2)

        for nCur in self.curTree.iterLeavesNoRoot():
            mtNum = nCur.br.parts[theProposal.pNum].rMatrixNum
            mtCur = mpCur.rMatrices[mtNum]
            mtProp = mpProp.rMatrices[mtNum]

            # Make proposals. Result of the proposal goes into mtProp.val
            p4.func.gsl_ran_dirichlet(theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            while  mtProp.val.min() < var.RATE_MIN:
                for i in range(ratesLen):
                    if mtProp.val[i] < var.RATE_MIN:
                        mtProp.val[i] += (1.0 + (1.1 * random.random())) * var.RATE_MIN
                mtProp.val /= mtProp.val.sum()

            # log proposal ratios
            forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(
                ratesLen, theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(
                ratesLen, theProposal.tuning[self.tempNum] * mtProp.val, mtCur.val)
            self.logProposalRatio += reverseLnPdf - forwardLnPdf

            # prior ratio
            # dirPrams = mpCur.ndch2_leafAlpha * mtCur.empiricalComp
            dirPrams = mpCur.ndrh2_leafAlpha * mpCur.ndrh2_priorRefRMatrix
            lnPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(ratesLen, dirPrams, mtCur.val)
            lnPdfProps = pf.gsl_ran_dirichlet_lnpdf(ratesLen, dirPrams, mtProp.val)
            self.logPriorRatio += lnPdfProps - lnPdfCurrs 
        # self.logPriorRatio = 0.0   # to turn off prior 


    def proposeNdrh2_internalRatesDir(self, theProposal):
        gm = ['Chain.proposeNdrh2_internalRatesDir()']

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]
        assert mpCur.ndrh2

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0

        ratesLen = int(((mpCur.dim * mpCur.dim) - mpCur.dim) / 2)

        for nCur in self.curTree.iterInternalsNoRoot():
            mtNum = nCur.br.parts[theProposal.pNum].rMatrixNum
            mtCur = mpCur.rMatrices[mtNum]
            mtProp = mpProp.rMatrices[mtNum]

            # Make proposals. Result of the proposal goes into mtProp.val
            p4.func.gsl_ran_dirichlet(theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            while  mtProp.val.min() < var.RATE_MIN:
                for i in range(ratesLen):
                    if mtProp.val[i] < var.RATE_MIN:
                        mtProp.val[i] += (1.0 + random.random()) * var.RATE_MIN
                thisSum = mtProp.val.sum()
                mtProp.val /= thisSum

            # log proposal ratios
            forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(
                ratesLen, theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(
                ratesLen, theProposal.tuning[self.tempNum] * mtProp.val, mtCur.val)
            self.logProposalRatio += reverseLnPdf - forwardLnPdf

            # prior ratio
            dirPrams = mpCur.ndrh2_internalAlpha * mpCur.ndrh2_priorRefRMatrix
            lnPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(ratesLen, dirPrams, mtCur.val)
            lnPdfProps = pf.gsl_ran_dirichlet_lnpdf(ratesLen, dirPrams, mtProp.val)
            self.logPriorRatio += lnPdfProps - lnPdfCurrs            




    def proposeNdrh2_leafRatesDirAlpha(self, theProposal):
        # The Dirichlet hyperparameter alphaL, for leaves. 
        gm = ['Chain.proposNdrh2_leafRatesDirAlpha()']

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]

        # if the max (below) gets too big, then the prior ratio can be so bad
        # that it will never be accepted.
        NDRH2_ALPHAL_MIN = 1.0 
        NDRH2_ALPHAL_MAX = 10000. 

        ratesLen = int(((mpCur.dim * mpCur.dim) - mpCur.dim) / 2)

        if 0:
            # Slider proposal
            myTuning = 50.
            oldVal = mpCur.ndrh2_leafAlpha
            newVal = oldVal + (random.random() - 0.5) * myTuning

            # Linear reflect
            isGood = False
            while not isGood:
                if newVal < NDRH2_ALPHAL_MIN:
                    newVal = (NDRH2_ALPHAL_MIN - newVal) + NDRH2_ALPHAL_MIN
                elif newVal > NDRH2_ALPHAL_MAX:
                    newVal = NDRH2_ALPHAL_MAX - (newVal - NDRH2_ALPHAL_MAX)
                else:
                    isGood = True
            mpProp.ndrh2_leafAlpha = newVal
            self.logProposalRatio = 0.0

        if 1: 
            # Multiplier proposal
            #myTuning = 2.0 * math.log(1.2)
            myTuning = theProposal.tuning[self.tempNum]
            oldVal = mpCur.ndrh2_leafAlpha
            newVal = oldVal * math.exp((random.random() - 0.5) * myTuning)

            # Log reflect
            isGood = False
            while not isGood:
                if newVal < NDRH2_ALPHAL_MIN:
                    newVal = NDRH2_ALPHAL_MIN * NDRH2_ALPHAL_MIN  / newVal
                elif newVal > NDRH2_ALPHAL_MAX:
                    newVal = NDRH2_ALPHAL_MAX * NDRH2_ALPHAL_MAX / newVal
                else:
                    isGood = True
            mpProp.ndrh2_leafAlpha = newVal
            self.logProposalRatio = math.log(newVal / oldVal)
        
        # prior
        self.logPriorRatio = 0.0

        for nCur in self.curTree.iterLeavesNoRoot():
            mtNum = nCur.br.parts[theProposal.pNum].rMatrixNum
            mtCur = mpCur.rMatrices[mtNum]

            lnPdfProp = pf.gsl_ran_dirichlet_lnpdf(ratesLen, newVal * mpCur.ndrh2_priorRefRMatrix, mtCur.val)
            lnPdfCur  = pf.gsl_ran_dirichlet_lnpdf(ratesLen, oldVal * mpCur.ndrh2_priorRefRMatrix, mtCur.val)
            self.logPriorRatio += lnPdfProp - lnPdfCur
        # self.logPriorRatio = 0.0   # to turn off prior


    def proposeNdrh2_internalRatesDirAlpha(self, theProposal):
        # The Dirichlet hyperparameter alpha 
        gm = ['Chain.proposeNdrh2_internalRatesDirAlpha()']
        # all the rMatrices in one go.

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]

        NDRH2_ALPHAI_MIN = 1.0 
        NDRH2_ALPHAI_MAX = 10000.

        ratesLen = int(((mpCur.dim * mpCur.dim) - mpCur.dim) / 2)

        self.logProposalRatio = 0.0
        if 0:
            # Slider proposal
            myTuning = 50.
            oldVal = mpCur.ndrh2_internalAlpha
            newVal = oldVal + (random.random() - 0.5) * myTuning

            # Linear reflect
            isGood = False
            while not isGood:
                if newVal < NDRH2_ALPHAI_MIN:
                    newVal = (NDRH2_ALPHAI_MIN - newVal) + NDRH2_ALPHAI_MIN
                elif newVal > NDRH2_ALPHAI_MAX:
                    newVal = NDRH2_ALPHAI_MAX - (newVal - NDRH2_ALPHAI_MAX)
                else:
                    isGood = True
            mpProp.ndrh2_internalAlpha = newVal

        if 1: 
            # Multiplier proposal
            # myTuning = 2.0 * math.log(3.0)
            myTuning = theProposal.tuning[self.tempNum]
            oldVal = mpCur.ndrh2_internalAlpha
            newVal = oldVal * math.exp((random.random() - 0.5) * myTuning)

            # Log reflect
            isGood = False
            while not isGood:
                if newVal < NDRH2_ALPHAI_MIN:
                    newVal = NDRH2_ALPHAI_MIN * NDRH2_ALPHAI_MIN  / newVal
                elif newVal > NDRH2_ALPHAI_MAX:
                    newVal = NDRH2_ALPHAI_MAX * NDRH2_ALPHAI_MAX / newVal
                else:
                    isGood = True
            mpProp.ndrh2_internalAlpha = newVal
            self.logProposalRatio = math.log(newVal / oldVal)
        
        # prior
        self.logPriorRatio = 0.0
        for nCur in self.curTree.iterInternalsNoRoot():
            mtNum = nCur.br.parts[theProposal.pNum].rMatrixNum
            mtCur = mpCur.rMatrices[mtNum]

            lnPdfProp = pf.gsl_ran_dirichlet_lnpdf(ratesLen, newVal * mpCur.ndrh2_priorRefRMatrix, mtCur.val)
            lnPdfCur  = pf.gsl_ran_dirichlet_lnpdf(ratesLen, oldVal * mpCur.ndrh2_priorRefRMatrix, mtCur.val)
            self.logPriorRatio += lnPdfProp - lnPdfCur


    def proposeNdrh2_priorRefRMatrixDir(self, theProposal):
        gm = ['Chain.proposeNdrh2_priorRefRMatrixDir()']

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]
        assert mpCur.ndrh2

        ratesLen = int(((mpCur.dim * mpCur.dim) - mpCur.dim) / 2)

        mtCur = mpCur.ndrh2_priorRefRMatrix
        mtProp = mpProp.ndrh2_priorRefRMatrix
        assert len(mtProp) == ratesLen

        # Make proposals. Result of the proposal goes into mtProp
        p4.func.gsl_ran_dirichlet(theProposal.tuning[self.tempNum] * mtCur, mtProp)

        while  mtProp.min() < var.RATE_MIN:
            for i in range(ratesLen):
                if mtProp[i] < var.RATE_MIN:
                    mtProp[i] += (1.0 + (1.1 * random.random())) * var.RATE_MIN
            mtProp /= mtProp.sum()

        # log proposal ratios
        self.logProposalRatio = 0.0
        forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(
            ratesLen, theProposal.tuning[self.tempNum] * mtCur, mtProp)
        reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(
            ratesLen, theProposal.tuning[self.tempNum] * mtProp, mtCur)
        self.logProposalRatio = reverseLnPdf - forwardLnPdf

        # prior ratio. All the RMatrices.  RMatrices and alphas will be the same (cur and prop).
        self.logPriorRatio = 0.0
        for nCur in self.curTree.iterNodesNoRoot():
            mtNum = nCur.br.parts[theProposal.pNum].rMatrixNum
            mtRMatrix = mpCur.rMatrices[mtNum]

            if nCur.isLeaf:
                thisAlpha = mpCur.ndrh2_leafAlpha
            else:
                thisAlpha = mpCur.ndrh2_internalAlpha

            lnPdfCur = pf.gsl_ran_dirichlet_lnpdf(ratesLen, thisAlpha * mtCur, mtRMatrix.val)
            lnPdfProp = pf.gsl_ran_dirichlet_lnpdf(ratesLen, thisAlpha * mtProp, mtRMatrix.val)
            self.logPriorRatio += lnPdfProp - lnPdfCur
            
        #self.logPriorRatio = 0.0   # to turn off prior 


    ####################################
    ####################################


    def proposeNdch2_leafCompsDir(self, theProposal):
        gm = ['Chain.proposeNdch2_leafCompsDir()']

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]
        assert mpCur.ndch2

        # Does this work for polytomies?  With that in mind, I iterate over
        # nodes rather than comps.

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0

        #print(f"proposeNdch2_leafCompsDir() mpCur.ndch2_priorRefComp is {mpCur.ndch2_priorRefComp}")
        #sys.stdout.flush()

        for nCur in self.curTree.iterLeavesNoRoot():
            mtNum = nCur.parts[theProposal.pNum].compNum
            mtCur = mpCur.comps[mtNum]
            mtProp = mpProp.comps[mtNum]

            # Make proposals. Result of the proposal goes into mtProp.val
            p4.func.gsl_ran_dirichlet(theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            while  mtProp.val.min() < var.PIVEC_MIN:
                for i in range(mpCur.dim):
                    if mtProp.val[i] < var.PIVEC_MIN:
                        mtProp.val[i] += (1.0 + (1.1 * random.random())) * var.PIVEC_MIN
                mtProp.val /= mtProp.val.sum()

            # log proposal ratios
            forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(
                mpCur.dim, theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(
                mpCur.dim, theProposal.tuning[self.tempNum] * mtProp.val, mtCur.val)
            self.logProposalRatio += reverseLnPdf - forwardLnPdf

            # prior ratio
            dirPrams = mpCur.ndch2_leafAlpha * mpCur.ndch2_priorRefComp
            lnPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, dirPrams, mtCur.val)
            lnPdfProps = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, dirPrams, mtProp.val)
            self.logPriorRatio += lnPdfProps - lnPdfCurrs 
            # self.logPriorRatio = 0.0   # to turn off prior 



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
            p4.func.gsl_ran_dirichlet(theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            while  mtProp.val.min() < var.PIVEC_MIN:
                for i in range(mpCur.dim):
                    if mtProp.val[i] < var.PIVEC_MIN:
                        mtProp.val[i] += (1.0 + random.random()) * var.PIVEC_MIN
                thisSum = mtProp.val.sum()
                mtProp.val /= thisSum

            # log proposal ratios
            forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(
                mpCur.dim, theProposal.tuning[self.tempNum] * mtCur.val, mtProp.val)
            reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(
                mpCur.dim, theProposal.tuning[self.tempNum] * mtProp.val, mtCur.val)
            self.logProposalRatio += reverseLnPdf - forwardLnPdf

            # prior ratio
            dirPrams = mpCur.ndch2_internalAlpha * mpCur.ndch2_priorRefComp  
            lnPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, dirPrams, mtCur.val)
            lnPdfProps = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, dirPrams, mtProp.val)
            self.logPriorRatio += lnPdfProps - lnPdfCurrs
            # self.logPriorRatio = 0.0   # to turn off prior 


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
            #myTuning = 2.0 * math.log(1.2)
            myTuning = theProposal.tuning[self.tempNum]
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
        # thisComp = mtCur.empiricalComp
        thisComp = mpCur.ndch2_priorRefComp
        for nCur in self.curTree.iterLeavesNoRoot():
            mtNum = nCur.parts[theProposal.pNum].compNum
            mtCur = mpCur.comps[mtNum]

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
        NDCH2_ALPHAI_MAX = 10000.

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
            # myTuning = 2.0 * math.log(3.0)
            myTuning = theProposal.tuning[self.tempNum]
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
        thisComp = mpCur.ndch2_priorRefComp
        for nCur in self.curTree.iterInternals():
            mtNum = nCur.parts[theProposal.pNum].compNum
            mtCur = mpCur.comps[mtNum]

            lnPdfProp = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, newVal * thisComp, mtCur.val)
            lnPdfCur = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, oldVal * thisComp, mtCur.val)
            self.logPriorRatio += lnPdfProp - lnPdfCur




    def proposeNdch2_priorRefCompDir(self, theProposal):
        gm = ['Chain.proposeNdch2_priorRefCompDir()']

        mpCur = self.curTree.model.parts[theProposal.pNum]
        mpProp = self.propTree.model.parts[theProposal.pNum]
        assert mpCur.ndch2

        mtCur = mpCur.ndch2_priorRefComp
        mtProp = mpProp.ndch2_priorRefComp

        # Make proposals. Result of the proposal goes into mtProp
        p4.func.gsl_ran_dirichlet(theProposal.tuning[self.tempNum] * mtCur, mtProp)

        while  mtProp.min() < var.PIVEC_MIN:
            for i in range(mpCur.dim):
                if mtProp[i] < var.PIVEC_MIN:
                    mtProp[i] += (1.0 + (1.1 * random.random())) * var.PIVEC_MIN
            mtProp /= mtProp.sum()

        # log proposal ratios
        self.logProposalRatio = 0.0
        forwardLnPdf = pf.gsl_ran_dirichlet_lnpdf(
            mpCur.dim, theProposal.tuning[self.tempNum] * mtCur, mtProp)
        reverseLnPdf = pf.gsl_ran_dirichlet_lnpdf(
            mpCur.dim, theProposal.tuning[self.tempNum] * mtProp, mtCur)
        self.logProposalRatio = reverseLnPdf - forwardLnPdf

        # prior ratio. All the comps.  Comps and alphas will be the same (cur and prop).
        self.logPriorRatio = 0.0
        for nCur in self.curTree.iterNodes():
            mtNum = nCur.parts[theProposal.pNum].compNum
            mtComp = mpCur.comps[mtNum]

            if nCur.isLeaf:
                thisAlpha = mpCur.ndch2_leafAlpha
            else:
                thisAlpha = mpCur.ndch2_internalAlpha

            lnPdfCur = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, thisAlpha * mtCur, mtComp.val)
            lnPdfProp = pf.gsl_ran_dirichlet_lnpdf(mpCur.dim, thisAlpha * mtProp, mtComp.val)
            self.logPriorRatio += lnPdfProp - lnPdfCur
            
        #self.logPriorRatio = 0.0   # to turn off prior 


    ############################################ topology
    ############################################

    def proposeRoot3(self, theProposal):
        """For non-biRooted trees.  Root on another internal node."""

        candidates = [n for n in self.propTree.iterInternalsNoRoot()]
        if len(candidates):
            oldRoot = self.propTree.root
            newRoot = random.choice(candidates)
            self.propTree.reRoot(
                newRoot, moveInternalName=False, 
                fixRawSplitKeys=self.mcmc.constraints)
            if self.mcmc.stickyRootComp:
                # move the compNum's too
                nParts = len(oldRoot.parts)
                for mpNum in range(nParts):
                    mpOld = oldRoot.parts[mpNum]
                    mpNew = newRoot.parts[mpNum]
                    # ... the switch
                    savedOldCompNum = mpOld.compNum
                    mpOld.compNum = mpNew.compNum
                    mpNew.compNum = savedOldCompNum
        else:
            print("Chain.proposeRoot3().  No other internal nodes.  Fix me.")
        
        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0



    def proposeRoot3n(self, theProposal):
        """For non-biRooted trees.  Root on a neighbouring internal node."""

        candidates = [n for n in self.propTree.root.iterChildren() if not n.isLeaf]
        if len(candidates):
            oldRoot = self.propTree.root
            newRoot = random.choice(candidates)
            self.propTree.reRoot(
                newRoot, moveInternalName=False, 
                fixRawSplitKeys=self.mcmc.constraints)
            if self.mcmc.stickyRootComp:
                # move the compNum's too
                nParts = len(oldRoot.parts)
                for mpNum in range(nParts):
                    mpOld = oldRoot.parts[mpNum]
                    mpNew = newRoot.parts[mpNum]
                    # ... the switch
                    savedOldCompNum = mpOld.compNum
                    mpOld.compNum = mpNew.compNum
                    mpNew.compNum = savedOldCompNum
        else:
            print("Chain.proposeRoot3().  No other internal nodes.  Fix me.")

        # There is a proposal imbalance, so we need to calculate the
        # proposal ratio.  It has to do with the number of non-leaves, as
        # leaves are not root candidates. 
        oldRootNNonLeaves = len(candidates)
        newRootNNonLeaves = len([n for n in self.propTree.root.iterChildren() if not n.isLeaf])
        proposalRatio = oldRootNNonLeaves/newRootNNonLeaves
        self.logProposalRatio = math.log(proposalRatio)
        self.logPriorRatio = 0.0



    def proposeRoot2(self, theProposal):
        """For bi-rooted trees; slides the node, perhaps past other nodes.

        Branch lengths are not affected.
        """

        gm = ["Chain.proposeRoot2()"]

        # Check that we have a bifurcating root.
        assert self.propTree.root.leftChild and self.propTree.root.leftChild.sibling
        assert not self.propTree.root.leftChild.sibling.sibling

        mvDist = theProposal.tuning[self.tempNum] * random.random()
        # print(gm[0], "Setting mvDist ")
        # mvDist = 0.15
        oldRoot = self.propTree.root
        curNode = self.propTree.root
        while True:
            dests = [n for n in curNode.iterChildren()]
            dest = random.choice(dests)
            #print("Got dests %s, choose dest %i" % ([n.nodeNum for n in dests], dest.nodeNum))

            if mvDist > dest.br.len:
                #print("go to node %i" % dest.nodeNum)
                mvDist -= dest.br.len
                curNode = dest
                if curNode.isLeaf:
                    #print("reRoot to node %i" % curNode.nodeNum)
                    self.propTree.reRoot(curNode, checkBiRoot=False)
            else:
                break
        # print(gm[0], "finished loop with curNode %i, dest=%i, and mvDist %4f" % (curNode.nodeNum, dest.nodeNum, mvDist))

        if dest.parent == curNode:     # curNode ----- dest
            pass
        else:
            print("curNode %i, dest %i" % (curNode.nodeNum, dest.nodeNum))
            self.propTree.draw()
            print("xx fix me")
            raise P4Error()

        if curNode == oldRoot:   # then it is easy
            # print(gm[0], "curNode is oldRoot")
            # as there may have been some rearrangements in the loop above, reRoot
            self.propTree.reRoot(oldRoot, checkBiRoot=False) 
            if dest == self.propTree.root.leftChild:
                dest.br.len -= mvDist
                if dest.br.len <= var.BRLEN_MIN:
                    dest.br.len += var.BRLEN_MIN
                dest.sibling.br.len += mvDist
            elif dest == self.propTree.root.leftChild.sibling:
                dest.br.len -= mvDist
                if dest.br.len <= var.BRLEN_MIN:
                    dest.br.len += var.BRLEN_MIN
                self.propTree.root.leftChild.br.len += mvDist
            else:
                raise P4Error("This should not happen.")

        elif dest == oldRoot:   # also easy
            # print(gm[0], "dest is oldRoot")
            pNode = oldRoot.leftChild
            # curNode---------oldRoot------pNode
            #                  dest
            pNode.br.len += (oldRoot.br.len - mvDist)
            oldRoot.br.len = mvDist
            if oldRoot.br.len <= var.BRLEN_MIN:
                oldRoot.br.len += var.BRLEN_MIN
            self.propTree.reRoot(oldRoot, checkBiRoot=False)
            self.propTree.preAndPostOrderAreValid = False   # check if this is necessary?

        else:
            #print(gm[0], "Neither curNode nor dest is the root. Need to rearrange.")

            # Remove the current root.  
            if oldRoot.parent:
                qNode = oldRoot.parent
            else:
                # reRoot
                self.propTree.reRoot(self.propTree.root.leftChild, checkBiRoot=False)
                qNode = oldRoot.parent

            #self.propTree.draw()
            #print("in the tree above, the root should have a parent.)
            #print("qNode is %i, dest = %i" % (qNode.nodeNum, dest.nodeNum))

            pNode = oldRoot.leftChild      #  qNode --- oldRoot ----pNode
            sib = oldRoot.sibling
            leftSib = oldRoot.leftSibling()

            # remove the oldRoot from the tree.

            #             +----------leftSib
            #             |
            #    ---------qNode              +----------
            #             +---------0--------pNode
            #             |     oldRoot      +----------
            #             |
            #             +---------sib

            pNode.parent = qNode
            pNode.sibling = sib  # might be None
            if leftSib:
                leftSib.sibling = pNode
            else:
                qNode.leftChild = pNode
            oldRoot.sibling = None
            oldRoot.parent = None
            oldRoot.leftChild = None
            pNode.br.len += oldRoot.br.len

            if dest.parent == curNode:     # curNode ----- dest
                pass
            else:
                print("curNode %i, dest %i" % (curNode.nodeNum, dest.nodeNum))
                self.setPreAndPostOrder()
                self.draw()
                print("fix me")
                raise P4Error

            #oldRoot.nodeNum = var.NO_ORDER  # out of the tree
            #self.propTree.setPreAndPostOrder()
            #self.propTree.dump(node=True)
            #self.propTree.draw()

            #             +----------leftSib (may not exist)
            #             |
            #    ---------qNode       +----------
            #             +-----------pNode
            #             |           +----------
            #             | 
            #             +---------sib (may be None)

            # now re-attach the old root
            qNode = curNode
            pNode = dest                  #  to be:  qNode --- root ----pNode

            assert pNode.parent == qNode

            sib = pNode.sibling
            leftSib = pNode.leftSibling()

            pNode.parent = oldRoot
            oldRoot.parent = qNode
            oldRoot.leftChild = pNode
            pNode.sibling = None
            oldRoot.sibling = sib   # may be None
            if leftSib:
                leftSib.sibling = oldRoot
            else:
                qNode.leftChild = oldRoot
            oldRoot.br.len = mvDist
            pNode.br.len -= mvDist
            for n in [oldRoot, pNode]:
                if n.br.len < var.BRLEN_MIN:
                    n.br.len += var.BRLEN_MIN
            self.propTree.reRoot(oldRoot)
            self.propTree.preAndPostOrderAreValid = False

            #oldRoot.nodeNum = 0
            #self.propTree.dump(node=True)
            #self.propTree.setPreAndPostOrder()
            #self.propTree.draw()


        assert self.propTree.root.leftChild and self.propTree.root.leftChild.sibling
        assert not self.propTree.root.leftChild.sibling.sibling

        # checking for constraint violation is done in Chain.proposeSp() (which calls this method).

        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0



    def proposeBrLen(self, theProposal):
        #gm = ['Chain.proposeBrLen']

        # Choose a node.
        nodesNoRoot = [n for n in self.propTree.iterNodesNoRoot()]
        theNode = random.choice(nodesNoRoot)
        #theNode = self.propTree.nodes[1]
        oldBrLen = theNode.br.len

        if 1:  # "Multiplier" proposal
            newBrLen = oldBrLen * \
                math.exp(theProposal.tuning[self.tempNum] * (random.random() - 0.5))

            # Logarithmic reflect if needed
            while (newBrLen < var.BRLEN_MIN) or (newBrLen > var.BRLEN_MAX):
                if newBrLen < var.BRLEN_MIN:
                    newBrLen = var.BRLEN_MIN * var.BRLEN_MIN / newBrLen
                elif newBrLen > var.BRLEN_MAX:
                    newBrLen = var.BRLEN_MAX * var.BRLEN_MAX / newBrLen
            theNode.br.len = newBrLen
            self.logProposalRatio = math.log(newBrLen / oldBrLen)

        else:  # Sliding window.
            newBrLen = oldBrLen + (theProposal.tuning[self.tempNum] * (random.random() - 0.5))
            #newBrLen = oldBrLen + (2.0 * (random.random() - 0.5))

            # Linear reflect if needed
            while (newBrLen < var.BRLEN_MIN) or (newBrLen > var.BRLEN_MAX):
                if newBrLen < var.BRLEN_MIN:
                    newBrLen = (var.BRLEN_MIN - newBrLen) + var.BRLEN_MIN
                elif newBrLen > var.BRLEN_MAX:
                    newBrLen = var.BRLEN_MAX - (newBrLen - var.BRLEN_MAX)
            theNode.br.len = newBrLen
            self.logProposalRatio = 0.0

        if 0 and hasattr(self.mcmc.tunings, 'doInternalBrLenPrior') and self.mcmc.tunings.doInternalBrLenPrior:
            if theNode.isLeaf:
                self.logPriorRatio = self.mcmc.tunings.brLenPriorLambda * \
                    (oldBrLen - newBrLen)
            else:
                self.logPriorRatio = self.mcmc.tunings.brLenPriorLambdaForInternals * \
                    (oldBrLen - newBrLen)
        else:
            if theProposal.brLenPriorType == 'exponential':
                self.logPriorRatio = theProposal.brLenPriorLambda * \
                    (oldBrLen - newBrLen)
            else:
                self.logPriorRatio = 0.

        if var.doMcmcSp:
            theNode.br.lenChanged = True

    def proposeAllBrLens(self, theProposal):
        gm = ['Chain.proposeAllBrLens']
        pTree = self.propTree
        self.logPriorRatio = 0.0
        self.logProposalRatio = 0.0

        if 0:    # multiplier proposal, does not work
            nBranches = 0
            oldBrLens = []
            newBrLens = []
            for n in pTree.iterNodesNoRoot():
                oldBrLen = n.br.len
                oldBrLens.append(oldBrLen)
                newBrLen = oldBrLen *  math.exp(theProposal.tuning * (random.random() - 0.5))

                # Logarithmic reflect if needed
                while (newBrLen < var.BRLEN_MIN) or (newBrLen > var.BRLEN_MAX):
                    if newBrLen < var.BRLEN_MIN:
                        newBrLen = var.BRLEN_MIN * var.BRLEN_MIN / newBrLen
                    elif newBrLen > var.BRLEN_MAX:
                        newBrLen = var.BRLEN_MAX * var.BRLEN_MAX / newBrLen

                newBrLens.append(newBrLen)
                nBranches += 1
                n.br.len = newBrLen
                if theProposal.brLenPriorType == 'exponential':
                    self.logPriorRatio += theProposal.brLenPriorLambda * \
                                          (oldBrLen - newBrLen)
                else:
                    pass # log prior remains zero for uniform prior
            # This hasting ratio appears to be wrong.
            self.logProposalRatio = nBranches * math.log(sum(newBrLens)/sum(oldBrLens))

        else:  # sliding window
            for n in pTree.iterNodesNoRoot():
                oldBrLen = n.br.len

                # sliding window
                newBrLen = oldBrLen + (theProposal.tuning[self.tempNum] * (random.random() - 0.5))

                # Linear reflect if needed
                while (newBrLen < var.BRLEN_MIN) or (newBrLen > var.BRLEN_MAX):
                    if newBrLen < var.BRLEN_MIN:
                        newBrLen = (var.BRLEN_MIN - newBrLen) + var.BRLEN_MIN
                    elif newBrLen > var.BRLEN_MAX:
                        newBrLen = var.BRLEN_MAX - (newBrLen - var.BRLEN_MAX)

                n.br.len = newBrLen
                if theProposal.brLenPriorType == 'exponential':
                    self.logPriorRatio += theProposal.brLenPriorLambda * \
                                          (oldBrLen - newBrLen)
                else:
                    pass # log prior remains zero for uniform prior

            

    def proposeLocal(self, theProposal):
        """from BAMBE and MrBayes."""
        
        # doAbort is set if brLens are too long or too short, or if a
        # constraint is violated.

        #global localCalls
        #localCalls += 1
        # print 'localCalls %i' % localCalls

        gm = ['Chain.proposeLocal()']
        theProposal.topologyChanged = 0
        theProposal.doAbort = False
        pTree = self.propTree
        dbug = False
        if 0 and self.mcmc.gen == 0:
            dbug = True

        if dbug:
            for n in pTree.iterNodesNoRoot():
                n.br.oldLen = n.br.len
                n.br.oldNode = n

        if dbug:
            print("=" * 80)
            # print "proposeLocal() starting with this tree ..."
            #pTree.draw(width=80, addToBrLen=0.0)
            # for n in pTree.iterInternalsNoRoot():
            #    print n.nodeNum, n.br.splitKey
            for n in pTree.iterNodes():
                n.oldName = n.name
            for n in pTree.iterLeavesNoRoot():
                n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
            for n in pTree.iterInternalsNoRoot():
                n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
            print("proposeLocal() starting with this tree ...")
            pTree.draw(width=100, addToBrLen=0.2, model=True)
            if self.mcmc.constraints:
                pTree.checkSplitKeys(useOldName=True)

        if pTree.root.getNChildren() == 2:
            gm.append("This method is not working for biRoot'd trees yet.")
            raise P4Error(gm)

        assert pTree.nInternalNodes > 1, "For local, we need trees with more than 1 internal node."

        # We want to find a node with a great-grandparent, but that might be
        # impossible!, eg if it is a 4-taxon tree rooted on an internal
        # node.  So at least find a node with a grandparent.
        usedNodes = [n for n in pTree.iterNodes()]
        candidateC = random.choice(usedNodes)
        safety = 0
        while 1:
            if candidateC.parent and candidateC.parent.parent:
                break
            else:
                safety += 1
                if safety > 100:
                    pTree.draw()
                    gm.append(
                        "Unable to find a node with a grandparent after 100 tries.")
                    gm.append("The propTree has %i internal nodes." %
                              pTree.nInternalNodes)
                    del(pTree.nInternalNodes)
                    gm.append(
                        "Recalculated: the propTree has %i internal nodes." % pTree.nInternalNodes)
                    raise P4Error(gm)
                candidateC = random.choice(usedNodes)

        # Check whether candidateC has a great-grandparent.  (We know it
        # has a grandparent).  If so, then fine.  If not, we need to
        # re-root.
        oldRoot = None
        if candidateC.parent.parent.parent:
            pass
        else:
            oldRoot = pTree.root
            possibleRoots = []
            for n in oldRoot.iterChildren():
                if n != candidateC.parent:
                    possibleRoots.append(n)
            if not possibleRoots:
                print("=" * 50)
                pTree.draw()
                gm.append(
                    "Programming error. Could not find any possibleRoots")
                raise P4Error(gm)

            newRoot = random.choice(possibleRoots)
            pTree.reRoot(
                newRoot, moveInternalName=False, fixRawSplitKeys=self.mcmc.constraints)

        if 0 and dbug:
            print("candidateC is node %i" % candidateC.nodeNum)
            if oldRoot:
                print("I had to reRoot to node %i." % newRoot.nodeNum)
            pTree.draw(width=80, showInternalNodeNames=1, addToBrLen=0.0)

        # We want a tree like this:
        #                +------c
        #        +-------|(v)
        # +------|(u)    +------d
        # |      |
        # |(a)   +-------b
        # |
        # +------X (which might be None)

        # set up the nodes as in Larget and Simon MBE, pg 754, fig 4
        c = candidateC
        v = c.parent
        safety = 0
        while c != v.leftChild:
            pTree.rotateAround(v)
            safety += 1
            if safety > 100:
                gm.append("Unable to make c as v's leftChild.")
                raise P4Error(gm)
        assert c.sibling
        d = c.sibling
        u = v.parent
        safety = 0
        while v != u.leftChild:
            pTree.rotateAround(u)
            safety += 1
            if safety > 100:
                gm.append("Unable to make v as u's leftChild.")
                raise P4Error(gm)
        assert v.sibling
        b = v.sibling
        a = u.parent
        safety = 0
        while a.leftChild != u:
            pTree.rotateAround(a)
            safety += 1
            if safety > 100:
                gm.append("Unable to make u as a's leftChild.")
                raise P4Error(gm)

        if dbug:
            #v.oldName = v.name
            #u.oldName = u.name
            #c.oldName = c.name
            #d.oldName = d.name
            #a.oldName = a.name
            #b.oldName = b.name
            for n in pTree.iterNodes():
                n.name = n.oldName
            for n in pTree.iterLeavesNoRoot():
                n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
            for n in pTree.iterInternalsNoRoot():
                n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)

            if v.name:
                v.name += '(v)'
            else:
                v.name = 'v'
            if u.name:
                u.name += '(u)'
            else:
                u.name = 'u'
            if c.name:
                c.name += '(c)'
            else:
                c.name = 'c'
            if d.name:
                d.name += '(d)'
            else:
                d.name = 'd'
            if a.name:
                a.name += '(a)'
            else:
                a.name = 'a'
            if b.name:
                b.name += '(b)'
            else:
                b.name = 'b'

            print("Label the nodes a,b,c,d,u,v, and arrange into a 'standard form' ...")
            pTree.draw(
                width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)

        # At this point, the tree should look like this:
        #                +------c
        #        +-------|(v)
        # +------|(u)    +------d
        # |      |
        # |(a)   +-------b
        # |
        # +------X (which might be None)

        m = c.br.len + v.br.len + u.br.len
        x = u.br.len
        y = x + v.br.len
        # by default, 0.909 to 1.1
        newMRatio = math.exp(theProposal.tuning[self.tempNum] * (random.random() - 0.5))
        newM = m * newMRatio

        # Hopefully these checks will not be needed forever.
        ##tooShort = False
        # if c.br.len < var.BRLEN_MIN:
        ##    gm.append("c.br.len (%g) is too short" % c.br.len)
        ##    tooShort = True
        # elif v.br.len < var.BRLEN_MIN:
        ##    gm.append("v.br.len (%g) is too short" % v.br.len)
        ##    tooShort = True
        # elif u.br.len < var.BRLEN_MIN:
        ##    gm.append("u.br.len (%g) is too short" % u.br.len)
        ##    tooShort = True
        # elif newM < (3.0 * var.BRLEN_MIN):
        ##    gm.append("newM (%g) is too short." % newM)
        ##    tooShort = True
        # if tooShort:
        # if self.lastProposal: # a tuple, see a few lines below
        ##        gm.append("Last proposal: %s, accepted=%s, topologyChanged=%s" % self.lastProposal)
        # else:
        ##        gm.append("This appears to be the first proposal.")
        ##    raise P4Error(gm)

        if 0 and dbug:
            print()
            print("m, the sum of brLens from a up to c, is %f" % m)
            print("x, from a to u, is %f" % x)
            print("y, from a to v, is %f" % y)
            print("newMRatio is %f" % newMRatio)
            print("newM is %f" % newM)

        #################################################################
        # Detach either u or v, then re-attach somewhere between a and c.
        #################################################################

        if random.random() < 0.5:
            # detach u
            #                +------c
            #        +-------|(v)
            # +------|(u)    +------d
            # |      |
            # |(a)   +-------b
            # |
            # +------X
            ##
            #            +----------c
            # +----------|(v)
            # |(a)       +----------d
            # |
            # +----------X

            newY = y * newMRatio
            newX = random.random() * newM

            # newX should be at least var.BRLEN_MIN away from newY
            safety = 0
            while math.fabs(newY - newX) < var.BRLEN_MIN:
                newX = random.random() * newM
                safety += 1
                if safety > 100:
                    if dbug:
                        print("Unable to place newX sufficiently far away from newY")
                    theProposal.doAbort = True
                    return

            if 0 and dbug:
                print("Choose to detach node u (not v)")
                print("newY is (%f * %f =) %f" % (y, newMRatio, newY))
                print("newX, a random spot along newM, is %f" % newX)
                if newX < newY:
                    print("-> Since newX is still less than newY, there will be no topology change.")
                else:
                    print("-> Since newX is now more than newY, there will be a topology change.")

            a.leftChild = v
            v.parent = a
            v.sibling = u.sibling  # which might be None

            # now re-attach at newX
            if newX < newY:
                # no topology change, set up the same as above
                #            +----------c
                # +----------|(v)
                # |(a)       +----------d
                # |
                # +----------X
                ##
                #        newX    newY   newM
                #        +       +      +
                ##
                #                +------c
                #        +-------|(v)
                # +------|(u)    +------d
                # |      |
                # |(a)   +-------b
                # |
                # +------X
                ##
                a.leftChild = u
                u.parent = a
                u.leftChild = v
                u.sibling = v.sibling
                v.parent = u
                v.sibling = b
                u.br.len = newX
                v.br.len = newY - newX
                c.br.len = newM - newY
                if 0 and dbug:
                    print("-> detach u, reattach between a and v, so no topology change")
                    pTree.draw(
                        width=80, showInternalNodeNames=1, addToBrLen=0.0)

            else:
                # a topology change
                ##
                #            +----------c
                # +----------|(v)
                # |(a)       +----------d
                # |
                # +----------X
                ##
                ##
                #        newY    newX   newM
                #        +       +      +
                ##
                #                +------c
                #        +-------|(u)
                # +------|(v)    +------b
                # |      |
                # |(a)   +-------d
                # |
                # +------X
                v.leftChild = u
                u.parent = v
                u.leftChild = c
                u.sibling = d
                c.parent = u
                c.sibling = b
                v.br.len = newY
                u.br.len = newX - newY
                c.br.len = newM - newX
                if self.mcmc.constraints:
                    tempRawSplitKey = v.br.rawSplitKey
                    tempSplitKey = v.br.splitKey
                    v.br.rawSplitKey = u.br.rawSplitKey
                    v.br.splitKey = u.br.splitKey
                    u.br.rawSplitKey = tempRawSplitKey
                    u.br.splitKey = tempSplitKey
                pTree.preAndPostOrderAreValid = 0
                theProposal.topologyChanged = 1
                if dbug:
                    print("-> detach u, re-attach between v and c, so there is a topology change")
                    for n in pTree.iterNodes():
                        n.name = n.oldName
                    for n in pTree.iterLeavesNoRoot():
                        n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
                    for n in pTree.iterInternalsNoRoot():
                        n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)

                    if v.name:
                        v.name += '(v)'
                    else:
                        v.name = 'v'
                    if u.name:
                        u.name += '(u)'
                    else:
                        u.name = 'u'
                    if c.name:
                        c.name += '(c)'
                    else:
                        c.name = 'c'
                    if d.name:
                        d.name += '(d)'
                    else:
                        d.name = 'd'
                    if a.name:
                        a.name += '(a)'
                    else:
                        a.name = 'a'
                    if b.name:
                        b.name += '(b)'
                    else:
                        b.name = 'b'
                    pTree.draw(
                        width=100, showInternalNodeNames=1, addToBrLen=0.2)

        else:
            # detach v
            #                +------c
            #        +-------|(v)
            # +------|(u)    +------d
            # |      |
            # |(a)   +-------b
            # |
            # +------X
            ##
            ##
            #            +----------c
            # +----------|(u)
            # |(a)       +----------b
            # |
            # +----------X
            ##

            newX = x * newMRatio
            newY = random.random() * newM

            # newY should be at least var.BRLEN_MIN away from newX
            safety = 0
            while math.fabs(newY - newX) < var.BRLEN_MIN:
                newY = random.random() * newM
                safety += 1
                if safety > 100:
                    if dbug:
                        print("Unable to place newY sufficiently far away from newX")
                    theProposal.doAbort = True
                    return

            if 0 and dbug:
                print("Choose to detach node v (not u)")
                print("newX is (%f * %f =) %f" % (x, newMRatio, newX))
                print("newY, a random spot along newM, is %f" % newY)
                if newY < newX:
                    print("-> Since newY is now less than newX, there will be a topology change.")
                else:
                    print("-> Since newY is still more than newX, there will not be a topology change.")

            u.leftChild = c
            c.parent = u
            c.sibling = b

            # now reattach at newY
            if newX < newY:
                # no topology change
                ##
                #            +----------c
                # +----------|(u)
                # |(a)       +----------b
                # |
                # +----------X
                ##
                #        newX    newY   newM
                #        +       +      +
                ##
                #                +------c
                #        +-------|(v)
                # +------|(u)    +------d
                # |      |
                # |(a)   +-------b
                # |
                # +------X
                ##
                u.leftChild = v
                v.parent = u
                v.leftChild = c
                v.sibling = b
                c.parent = v
                c.sibling = d
                u.br.len = newX
                v.br.len = newY - newX
                c.br.len = newM - newY
                if 0 and dbug:
                    print("-> detach v, reattach between u and c, so no topology change")
                    pTree.draw(
                        width=80, showInternalNodeNames=1, addToBrLen=0.0)
            else:
                # with a topology change
                ##
                #            +----------c
                # +----------|(u)
                # |(a)       +----------b
                # |
                # +----------X
                ##
                ##
                #        newY    newX   newM
                #        +       +      +
                ##
                #                +------c
                #        +-------|(u)
                # +------|(v)    +------b
                # |      |
                # |(a)   +-------d
                # |
                # +------X
                a.leftChild = v
                v.parent = a
                v.leftChild = u
                v.sibling = u.sibling
                u.parent = v
                u.sibling = d
                v.br.len = newY
                u.br.len = newX - newY
                c.br.len = newM - newX
                if self.mcmc.constraints:
                    tempRawSplitKey = v.br.rawSplitKey
                    tempSplitKey = v.br.splitKey
                    v.br.rawSplitKey = u.br.rawSplitKey
                    v.br.splitKey = u.br.splitKey
                    u.br.rawSplitKey = tempRawSplitKey
                    u.br.splitKey = tempSplitKey
                pTree.preAndPostOrderAreValid = 0
                theProposal.topologyChanged = 1
                if dbug:
                    print("-> detach v, re-attach between a and u, so there is a topology change")
                    print("   (splitKeys are wrong on nodes v and u, in the figure below)")
                    for n in pTree.iterNodes():
                        n.name = n.oldName
                    for n in pTree.iterLeavesNoRoot():
                        n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
                    for n in pTree.iterInternalsNoRoot():
                        n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)

                    if v.name:
                        v.name += '(v)'
                    else:
                        v.name = 'v'
                    if u.name:
                        u.name += '(u)'
                    else:
                        u.name = 'u'
                    if c.name:
                        c.name += '(c)'
                    else:
                        c.name = 'c'
                    if d.name:
                        d.name += '(d)'
                    else:
                        d.name = 'd'
                    if a.name:
                        a.name += '(a)'
                    else:
                        a.name = 'a'
                    if b.name:
                        b.name += '(b)'
                    else:
                        b.name = 'b'
                    pTree.draw(
                        width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)
                    if self.mcmc.constraints:
                        pTree.checkSplitKeys(useOldName=True, glitch=False)

        # Check that new brLens are not too short or too long.  Ronquist
        # suggests that if that is the case, then just abort, rather than
        # fussing with reflections.
        if 0:
            if c.br.len < var.BRLEN_MIN or v.br.len < var.BRLEN_MIN or u.br.len < var.BRLEN_MIN:
                # if dbug:
                print("At least 1 brLen is too short.")
                theProposal.doAbort = True
                return
            elif c.br.len > var.BRLEN_MAX or v.br.len > var.BRLEN_MAX or u.br.len > var.BRLEN_MAX:
                # if dbug:
                print("At least 1 brLen is too long.  Aborting. (No big deal ...)")
                theProposal.doAbort = True
                return

        if 1:
            complain = False
            if c.br.len < var.BRLEN_MIN:
                if complain:
                    self.mcmc.logger.info("proposeLocal() tempNum=%i gen=%i  br c too short" % (self.tempNum, self.mcmc.gen))
                theProposal.doAbort = True
                return
            if v.br.len < var.BRLEN_MIN:
                if complain:
                    self.mcmc.logger.info("proposeLocal() tempNum=%i gen=%i  br v too short" % (self.tempNum, self.mcmc.gen))
                theProposal.doAbort = True
                return
            if u.br.len < var.BRLEN_MIN:
                if complain:
                    self.mcmc.logger.info("proposeLocal() tempNum=%i gen=%i  br u too short" % (self.tempNum, self.mcmc.gen))
                theProposal.doAbort = True
                return
            if c.br.len > var.BRLEN_MAX:
                if complain:
                    self.mcmc.logger.info("proposeLocal() tempNum=%i gen=%i  br c too long" % (self.tempNum, self.mcmc.gen))
                theProposal.doAbort = True
                return
            if v.br.len > var.BRLEN_MAX:
                if complain:
                    self.mcmc.logger.info("proposeLocal() tempNum=%i gen=%i  br v too long" % (self.tempNum, self.mcmc.gen))
                theProposal.doAbort = True
                return
            if u.br.len > var.BRLEN_MAX:
                if complain:
                    self.mcmc.logger.info("proposeLocal() tempNum=%i gen=%i  br u too long" % (self.tempNum, self.mcmc.gen))
                theProposal.doAbort = True
                return

        if self.mcmc.constraints and theProposal.topologyChanged:
            # Check whether any constraints have been involved, and if so abort.
            ##
            # It was like this:
            #                +------c
            #        +-------|(v)
            # +------|(u)    +------d
            # |      |
            # |(a)   +-------b
            # |
            # +------X
            ##
            # And now its like this:
            #                +------c
            #        +-------|(u)
            # +------|(v)    +------b
            # |      |
            # |(a)   +-------d
            # |
            # +------X
            #
            # So splits on branches on all nodes except u are unaffected.
            # The rawSplitKey on node v might be affected, if there was a
            # re-rooting, but the splitKey is still ok.  So the
            # rawSplitKey and the splitKey for node u in its new position
            # needs to be re-calculated.

            oldUSplitKey = u.br.splitKey
            pTree.recalculateSplitKeysOfNodeFromChildren(
                u, self.mcmc.constraints.allOnes)
            # print "u, node %i recalculated br.rawSplitKey=%s, br.splitKey = %s" % (
            #    u.nodeNum, u.br.rawSplitKey, u.br.splitKey)

            if dbug:
                for n in pTree.iterNodes():
                    n.name = n.oldName
                pTree.checkTaxNames()
                for n in pTree.iterLeavesNoRoot():
                    n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
                for n in pTree.iterInternalsNoRoot():
                    n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
                pTree.draw(width=100, showInternalNodeNames=1, addToBrLen=0.2)
                if self.mcmc.constraints:
                    # If v is wrong, it gets righted later.
                    pTree.checkSplitKeys(useOldName=True, glitch=False)

            for sk in self.mcmc.constraints.constraints:
                if oldUSplitKey == sk:
                    # print "oldUSplitKey was constrained to %i, and now it is
                    # gone. --> abort!" % sk
                    theProposal.doAbort = True
                    if dbug:
                        for n in pTree.iterNodes():
                            n.name = n.oldName
                    return

        if var.doMcmcSp:
            c.br.lenChanged = True
            u.br.lenChanged = True
            v.br.lenChanged = True

        # Now evaluate the prior ratio.  If doInternalBrLenPrior is set,
        # then we do it like Yang and Rannala Syst Biol 54(3):455-470,
        # 2005.  Do it before reRooting (does it make any difference?)
        # print self.mcmc.tunings.doInternalBrLenPrior
        if 0 and hasattr(self.mcmc.tunings, 'doInternalBrLenPrior') and self.mcmc.tunings.doInternalBrLenPrior:
            # originally this
            #                +------c
            #        +-------|(v)
            # +------|(u)    +------d
            # |      |
            # |(a)   +-------b
            # |
            # +------X (which might be None)
            # Either u or v was detached, and reattached somewhere between a
            # and c.

            # See the "long way" calculation below, which shows how it can
            # be done with a single prior.  The complication here is that
            # there are two priors, depending on whether it is a leaf or
            # not.
            theSum = 0.0
            if newX < newY:  # no topology change
                # Nodes a or c might be leaf nodes.  Nodes u and v cannot be.
                # Do the 3 edges in turn.  First the a edge, then the internal
                # edge, then the c edge.
                if a.isLeaf:
                    theSum += (theProposal.brLenPriorLambda * x) - \
                        (theProposal.brLenPriorLambda * newX)
                else:
                    theSum += (theProposal.brLenPriorLambdaForInternals * x) - (
                        theProposal.brLenPriorLambdaForInternals * newX)
                theSum += (theProposal.brLenPriorLambdaForInternals * (y - x)) - (
                    theProposal.brLenPriorLambdaForInternals * (newY - newX))
                if c.isLeaf:
                    theSum += (theProposal.brLenPriorLambda * (m - y)) - \
                        (theProposal.brLenPriorLambda * (newM - newY))
                else:
                    theSum += (theProposal.brLenPriorLambdaForInternals * (m - y)) - (
                        theProposal.brLenPriorLambdaForInternals * (newM - newY))
            else:  # with topology change
                if a.isLeaf:
                    theSum += (theProposal.brLenPriorLambda * x) - \
                        (theProposal.brLenPriorLambda * newY)
                else:
                    theSum += (theProposal.brLenPriorLambdaForInternals * x) - (
                        theProposal.brLenPriorLambdaForInternals * newY)
                theSum += (theProposal.brLenPriorLambdaForInternals * (y - x)) - (
                    theProposal.brLenPriorLambdaForInternals * (newX - newY))
                if c.isLeaf:
                    theSum += (theProposal.brLenPriorLambda * (m - y)) - \
                        (theProposal.brLenPriorLambda * (newM - newX))
                else:
                    theSum += (theProposal.brLenPriorLambdaForInternals * (m - y)) - (
                        theProposal.brLenPriorLambdaForInternals * (newM - newX))

            if 0:
                # Slow check, via priorDensities.
                theta1 = theProposal.brLenPriorLambda
                theta2 = theProposal.brLenPriorLambdaForInternals

                if newX < newY:
                    if a.isLeaf:
                        theta = theProposal.brLenPriorLambda
                    else:
                        theta = theProposal.brLenPriorLambdaForInternals
                    prDensNu1 = theta * math.exp(-theta * x)
                    prDensNuStar1 = theta * math.exp(-theta * newX)

                    theta = theProposal.brLenPriorLambdaForInternals
                    prDensNu2 = theta * math.exp(-theta * (y - x))
                    prDensNuStar2 = theta * math.exp(-theta * (newY - newX))
                    if c.isLeaf:
                        theta = theProposal.brLenPriorLambda
                    else:
                        theta = theProposal.brLenPriorLambdaForInternals
                    prDensNu3 = theta * math.exp(-theta * (m - y))
                    prDensNuStar3 = theta * math.exp(-theta * (newM - newY))
                else:
                    if a.isLeaf:
                        theta = theProposal.brLenPriorLambda
                    else:
                        theta = theProposal.brLenPriorLambdaForInternals
                    prDensNu1 = theta * math.exp(-theta * x)
                    prDensNuStar1 = theta * math.exp(-theta * newY)

                    theta = theProposal.brLenPriorLambdaForInternals
                    prDensNu2 = theta * math.exp(-theta * (y - x))
                    prDensNuStar2 = theta * math.exp(-theta * (newX - newY))
                    if c.isLeaf:
                        theta = theProposal.brLenPriorLambda
                    else:
                        theta = theProposal.brLenPriorLambdaForInternals
                    prDensNu3 = theta * math.exp(-theta * (m - y))
                    prDensNuStar3 = theta * math.exp(-theta * (newM - newX))
                prRat = (prDensNuStar1 * prDensNuStar2 * prDensNuStar3) / \
                    (prDensNu1 * prDensNu2 * prDensNu3)
                logPrRat = math.log(prRat)

                if math.fabs(logPrRat - theSum) > 1.e-10:
                    print("xxzz differs.  logPrRat=%g, theSum=%g" % (logPrRat, theSum))
                # else:
                #    print "s",

            self.logPriorRatio = theSum

        else:  # Do not doInternalBrLenPrior
            if theProposal.brLenPriorType == 'uniform':
                self.logPriorRatio = 0.0
            elif theProposal.brLenPriorType == 'exponential':
                self.logPriorRatio = theProposal.brLenPriorLambda *  (m - newM)
                if 0:  # Do the same calculation the long way, edge by edge.
                    # print "logPriorRatio = %+.4f" % self.logPriorRatio,
                    foo0 = (theProposal.brLenPriorLambda * m) - \
                        (theProposal.brLenPriorLambda * newM)
                    # print "%+.4f" % foo0,

                    foo = 0.0
                    if newX < newY:  # no topology change
                        #        newX    newY   newM
                        #        +       +      +
                        ##
                        #                +------c
                        #        +-------|(v)
                        # +------|(u)    +------d
                        # |      |
                        # |(a)   +-------b
                        # |
                        # +------X
                        foo += (theProposal.brLenPriorLambda * x) - \
                            (theProposal.brLenPriorLambda * newX)
                        foo += (theProposal.brLenPriorLambda * (y - x)) - \
                            (theProposal.brLenPriorLambda *
                             (newY - newX))
                        foo += (theProposal.brLenPriorLambda * (m - y)) - \
                            (theProposal.brLenPriorLambda *
                             (newM - newY))
                    else:  # with topology change
                        #        newY    newX   newM
                        #        +       +      +
                        ##
                        #                +------c
                        #        +-------|(u)
                        # +------|(v)    +------b
                        # |      |
                        # |(a)   +-------d
                        # |
                        # +------X
                        foo += (theProposal.brLenPriorLambda * x) - \
                            (theProposal.brLenPriorLambda * newY)
                        foo += (theProposal.brLenPriorLambda * (y - x)) - \
                            (theProposal.brLenPriorLambda *
                             (newX - newY))
                        foo += (theProposal.brLenPriorLambda * (m - y)) - \
                            (theProposal.brLenPriorLambda *
                             (newM - newX))
                    # print "%+.4f" % foo
                    if (math.fabs(self.logPriorRatio - foo0) > 1.e-10):
                        print("differs-- foo0, %g %g" % (self.logPriorRatio, foo0))
                    if (math.fabs(self.logPriorRatio - foo) > 1.e-10):
                        print("differs-- foo, %g %g" % (self.logPriorRatio, foo))
            else:
                raise P4Error("This should not happen.")

        if oldRoot:
            if dbug:
                print('-------------------- about to reRoot -----------')
                pTree.draw(
                    width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)
                if self.mcmc.constraints:
                    pTree.checkSplitKeys(useOldName=True, glitch=False)
                for n in pTree.iterNodes():
                    n.name = n.oldName

            pTree.reRoot(
                oldRoot, moveInternalName=False, fixRawSplitKeys=self.mcmc.constraints)

            if dbug:
                print('--------------after reRoot --------------')
                for n in pTree.iterLeavesNoRoot():
                    n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
                for n in pTree.iterInternalsNoRoot():
                    n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
                pTree.draw(
                    width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)
                if self.mcmc.constraints:
                    pTree.checkSplitKeys(useOldName=True)

            if self.mcmc.constraints and theProposal.topologyChanged:
                # Node v in its new position might now need to have its
                # rawSplitKey updated.  The splitKey should still be ok.
                children = [n for n in v.iterChildren()]
                # print "leaves above v, node %i: %s" % (v.nodeNum, [n.nodeNum for n in children])
                # for n in children:
                # print "    %2i  %4i  %4i" % (n.nodeNum, n.br.rawSplitKey,
                # n.br.splitKey)
                x = children[0].br.rawSplitKey
                # print "rawSplitKeys: ", x,
                for n in children[1:]:
                    y = n.br.rawSplitKey
                    # print y,
                    x = x | y  # '|' is bitwise "OR".
                v.br.rawSplitKey = x
                # print '=>', x
                # print "v, node %i recalculated br.rawSplitKey=%s" %
                # (v.nodeNum, v.br.rawSplitKey)

        if dbug:
            if theProposal.topologyChanged:
                print("The topology CHANGED")
            else:
                print("Topology -- no change.")
            # pTree.draw(width=80)

            for n in pTree.iterNodes():
                n.name = n.oldName
            pTree.checkTaxNames()
            for n in pTree.iterLeavesNoRoot():
                n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
            for n in pTree.iterInternalsNoRoot():
                n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
            pTree.draw(
                width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)
            if self.mcmc.constraints:
                pTree.checkSplitKeys(useOldName=True)
            #c.name = c.oldName
            #u.name = u.oldName
            #v.name = v.oldName
            #d.name = d.oldName
            #a.name = a.oldName
            #b.name = b.oldName
            if 0:
                del(c.oldName)
                del(u.oldName)
                del(v.oldName)
                del(d.oldName)
                del(a.oldName)
                del(b.oldName)
            for n in pTree.iterNodes():
                n.name = n.oldName

        self.logProposalRatio = 3.0 * math.log(newMRatio)

        if 0:
            for n in pTree.iterNodesNoRoot():
                if n.br.lenChanged:
                    print("l node %2i br.lenChanged" % n.nodeNum)
        if dbug:
            if self.mcmc.constraints:
                pTree.checkSplitKeys()
        # if self.mcmc.constraints:
        #    print "checkSplitKeys() at the end of local"
        #    pTree.checkSplitKeys()

        # Check if we have a new combo of comp and rMatrix
        if theProposal.topologyChanged:
            for pNum in range(pTree.model.nParts):
                # print("\n")
                # print(pTree.model.parts[pNum].bQETneedsReset)
                for cNum in range(pTree.model.parts[pNum].nComps):
                    for rNum in range(pTree.model.parts[pNum].nRMatrices):
                        if pTree.model.parts[pNum].bQETneedsReset[cNum][rNum]:
                            pf.p4_resetBQET(pTree.model.cModel, pNum, cNum, rNum)
                # print("after p4_resetBQET, ...")
                # print(pTree.model.parts[pNum].bQETneedsReset)
                
        if dbug:
            for n in pTree.iterNodesNoRoot():
                if math.fabs(n.br.len - n.br.oldLen) > 0.0000001:
                    # print "Node %2i br len changed" % n.nodeNum
                    assert n.br.lenChanged
                else:
                    if n.br.lenChanged:
                        print("Node %2i lenChanged set, but its the same length." % n.nodeNum)
                        raise P4Error

                # If the branch has been inverted, we will want to recalculate
                # the bigPDecks, even if the length has not really changed.
                # Trigger that intent by setting lenChanged.
                if n.br.oldNode != n:
                    # print "Node %2i branch: oldNode %2i" % (n.nodeNum,
                    # n.br.oldNode.nodeNum)
                    n.br.lenChanged = 1

    def proposeETBR_Blaise(self, theProposal):
        """Adapted and modified from the Evans and Sullivan version

        See the docstring for proposeETBR.  I have shamelessly and extensively
        borrowed code from Evans' Crux.

        Many thanks to Blaise Li who convincingly pointed out that LOCAL
        was not enough, and then pointed out that even the Crux version of
        eTBR did not appear to be reversible when there were polytomies.
        This version has a modification from Blaise as described
        in his poster 'An eTBR proposal for non-binary trees in MCMC
        Bayesian phylogeny', B Li and P Foster, presented (by BL) at the
        15th Evolutionary Biology Meeting, September 27-30, 2011
        Marseilles.

        """

        # doAbort is set if brLens are too long or too short, or if a
        # constraint is violated.  Constraints are not checked here, in this method.
        gm = ['Chain.proposeETBR_Blaise()']

        theProposal.topologyChanged = 0
        theProposal.doAbort = False
        pTree = self.propTree
        dbug = False

        for n in pTree.iterNodesNoRoot():
            n.br.oldLen = n.br.len
            n.br.oldNode = n

        if 0 and self.mcmc.gen == 217:
            dbug = True
            print("------------- eTBR gen %i -----------------" % self.mcmc.gen)
            if 0:
                currentLogLike = self.propTree.logLike
                self.propTree.calcLogLike(verbose=0)  # with _commonCStuff()
                theDiff = math.fabs(currentLogLike - self.propTree.logLike)
                if theDiff > 1.e-9:
                    gm.append("propTree like diff %f (%g)" %
                              (theDiff, theDiff))
                    raise P4Error(gm)

        oldRoot = pTree.root
        eA = None
        eX = None
        eY = None

        etbrLambda = theProposal.etbrLambda
        etbrPExt = theProposal.etbrPExt

        if 0 and dbug:
            print("=" * 50)
            pTree.draw()
            print("starting with the tree above.")

        # Choose a node, not the root.  In crux the choice for the original
        # branch can be any branch at all, including leaf branches.  Here, since
        # I choose a random node that is not the root, then the parent, the X0
        # end, will never be a leaf, but the Y0 end may be.  It will have edge
        # eA in Jason's diagram.  y0 will be the asterisk node in Jason's
        # diagram.  In eTBR, it (y0) may be extended, below, as the second of
        # the two eSPR moves.

        y0 = None
        while not y0:
            nNum = random.choice(pTree.preOrder)
            if nNum != var.NO_ORDER and nNum != pTree.root.nodeNum:
                y0 = pTree.node(nNum)
        x0 = y0.parent

        # The identity of x0 and y0 should be random.  That is, the first
        # direction of extension, that is the first eSPR, should be random.  The
        # way it is coded it below is that the x-direction is extended first,
        # then the y-direction.  However, the way I have it above, it is x0 =
        # y0.parent.  To make it properly random, switch ends with probability
        # 0.5.
        myRan = random.random()
        if myRan >= 0.5:
            x0 = y0
            y0 = x0.parent
            
        if dbug:
            print("y0 is node %i" % y0.nodeNum)
            y0.br.textDrawSymbol = '='
            if y0.name:
                y0.name += '_y0'
            else:
                y0.name = 'y0'
            if x0.name:
                x0.name += '_x0'
            else:
                x0.name = 'x0'

        # Name the edge here, used below when we modify the br.len
        if x0 == y0.parent:
            eA = y0.br
        else:
            eA = x0.br

        # #########
        # Extend x
        # #########
        # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
        r1 = x0

        #  If x0 is not a leaf, it is 'unconstrained'
        x0Degree = pTree.getDegree(x0)
        if x0Degree == 1:
            x0Uncon = False   # a leaf
        else:
            x0Uncon = True
        if x0Uncon:
            myRan = random.randrange(x0Degree - 1)
            if y0.parent == x0:  # x-side is down
                r0 = pTree.nextNode(y0, x0)
            elif x0.parent == y0: # x-side is up
                r0 = pTree.nextNode(x0, x0)
            else:
                gm.append("This shouldn't happen.")
                raise P4Error(gm)

            for i in range(myRan):
                r0 = pTree.nextNode(r0, x0)
            if r0 == x0:  # points down 
                r0 = r0.parent
            x1 = r0
            if dbug:
                if x1.name:
                    x1.name += '_x1'
                else:
                    x1.name = 'x1'

            # So we are set up like this ...
            # 
            #        eA        eX
            #   y0--------x0--------x1
            #             r1        r0
            #
            # ... and below, if r0 is unconstrained, 
            # it may turn into this (n? is just some other node) ---
            #
            #        eA        eX        eR
            #   y0--------x0--------x1--------n?
            #                       r1        r0
            #
            # and then maybe
            #        eA        eX                  eR
            #   y0--------x0--------x1--------n?--------n?
            #                                 r1        r0
            
            # We name the edge here, and use the name below when we modify the br.len
            if x0.parent == x1:
                eX = x0.br
            elif x1.parent == x0:
                eX = x1.br
            else:
                raise P4Error("This should not happen")

            #  if r0 is not a leaf, it is 'unconstrained'
            r0Degree = pTree.getDegree(r0)
            r0Uncon = (r0Degree > 1)
            while r0Uncon:
                # randomly determine whether to extend
                myRan = random.random()
                if etbrPExt < myRan:
                    break
                # choose an edge (direction) to go, and step to the next node
                myRan = random.randrange(r0Degree - 1)
                # We now call nextNode(spoke, hub) where hub is r0.  If r0
                # is the parent of r1, then the spoke of nextNode(spoke,
                # hub) is r1.  But if its the other way around, if r1 is
                # the parent of r0, then spoke is r0.
                if r1.parent == r0:
                    r0new = pTree.nextNode(r1, r0)
                elif r0.parent == r1:
                    r0new = pTree.nextNode(r0, r0)
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)

                for i in range(myRan):
                    r0new = pTree.nextNode(r0new, r0)
                if r0new == r0:
                    r0new = r0.parent
                r1 = r0
                r0 = r0new
                r0Degree = pTree.getDegree(r0)
                r0Uncon = (r0Degree > 1)
            if dbug:
                if r1.name:
                    r1.name += '_r1'
                else:
                    r1.name = 'r1'
                if r0.name:
                    r0.name += '_r0'
                else:
                    r0.name = 'r0'
                # pTree.draw()

            # Perform rearrangement unless it would be a no-op.  It would be a
            # no-op if r0 was still x1.
            if r0 == x1:
                # We did not extend, at all.
                if dbug:
                    print("No extension from x1 was done (because x1=r0), so no rearrangement on the x side.")
                pass
            else:
                # Do the rearrangement.  
                pTree.reRoot(r1, moveInternalName=False)
                assert x0.parent == x1
                assert y0.parent == x0
                assert r0.parent == r1

                if dbug:
                    x0.br.textDrawSymbol = 'X'
                    pTree.setPreAndPostOrder()
                    pTree.draw()
                    print("The drawing above is just before the rearrangement.")

                # Get children of x0 that are not y0
                ch = [n for n in x0.iterChildren() if n != y0]
                assert ch
                if len(ch) == 1:
                    myChoice = ch[0]
                else:
                    myChoice = random.choice(ch)

                # Three spr moves.
                pTree.pruneSubTreeWithoutParent(
                    myChoice, allowSingleChildNode=True)
                pTree.reconnectSubTreeWithoutParent(myChoice, x1)

                pTree.pruneSubTreeWithoutParent(x0)
                pTree.reconnectSubTreeWithoutParent(x0, r1)

                pTree.pruneSubTreeWithoutParent(r0)
                pTree.reconnectSubTreeWithoutParent(r0, x0)

        if 1 and dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()
            if x0Uncon and r0 != x1:
                print("The drawing above shows that X extended")
            else:
                print("The drawing above shows that X did not extend.")

        # #########
        # Extend y
        # #########
        # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
        s1 = y0

        #  If y0 is not a leaf, it is 'unconstrained'
        y0Degree = pTree.getDegree(y0)
        if y0Degree == 1:
            y0Uncon = False    # If y0 is a leaf
        else:
            y0Uncon = True

        if y0Uncon:
            myRan = random.randrange(y0Degree - 1)
            if x0.parent == y0:  # y-side is down
                s0 = pTree.nextNode(x0, y0)
            elif y0.parent == x0:  # y-side is up
                s0 = pTree.nextNode(y0, y0)
            else:
                gm.append("This shouldn't happen.")
                raise P4Error(gm)

            for i in range(myRan):
                s0 = pTree.nextNode(s0, y0)
            if s0 == y0:   # points down
                s0 = s0.parent
            y1 = s0
            if dbug:
                if y1.name:
                    y1.name += '_y1'
                else:
                    y1.name = 'y1'

            # name the edge here, to be used below when we modify the br.len
            if y0.parent == y1:
                eY = y0.br
            elif y1.parent == y0:
                eY = y1.br
            else:
                raise P4Error("This should not happen (name eY)")

            # So we are set up like this ...
            #   nS1      nS0
            #   nY0==eY==nY1

            #  if s0 is not a leaf, it is 'unconstrained'
            s0Degree = pTree.getDegree(s0)
            s0Uncon = (s0Degree > 1)
            while s0Uncon:
                # randomly determine whether to extend
                myRan = random.random()
                if etbrPExt < myRan:
                    break
                # choose an edge (direction) to go, and step to the next node
                myRan = random.randrange(s0Degree - 1)
                # We now call nextNode(spoke, hub) where hub is s0.  If s0
                # is the parent of s1, then the spoke of nextNode(spoke,
                # hub) is s1.  But if its the other way around, if s1 is
                # the parent of s0, then spoke is s0.
                if s1.parent == s0:
                    s0new = pTree.nextNode(s1, s0)
                elif s0.parent == s1:
                    s0new = pTree.nextNode(s0, s0)
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)

                for i in range(myRan):
                    s0new = pTree.nextNode(s0new, s0)
                if s0new == s0:
                   s0new = s0.parent
                s1 = s0
                s0 = s0new
                s0Degree = pTree.getDegree(s0)
                s0Uncon = (s0Degree > 1)

            if dbug:
                if s1.name:
                    s1.name += '_s1'
                else:
                    s1.name = 's1'
                if s0.name:
                    s0.name += '_s0'
                else:
                    s0.name = 's0'
                y1.br.textDrawSymbol = 'Y'
                # pTree.draw()

            # Perform rearrangement unless it would be a no-op.  It would be a
            # no-op if s0 was still y1.
            if s0 == y1:
                # We did not extend, at all.
                if dbug:
                    print("No Y extension was made (because y1=s0), so nothing to do.")
            else:
                # Do the rearrangement.  
                pTree.reRoot(s1, moveInternalName=False)
                assert y0.parent == y1
                assert x0.parent == y0
                assert s0.parent == s1

                if dbug:
                    y0.br.textDrawSymbol = 'Y'
                    pTree.setPreAndPostOrder()
                    pTree.draw()
                    print("The drawing above is the tree just before rearrangement on the Y side.")

                # Get children of y0 that are not x0
                ch = [n for n in y0.iterChildren() if n != x0]
                assert ch
                if len(ch) == 1:
                    myChoice = ch[0]
                else:
                    myChoice = random.choice(ch)

                # Three spr moves.
                # myChoice is returned, also
                pTree.pruneSubTreeWithoutParent(
                    myChoice, allowSingleChildNode=True)
                pTree.reconnectSubTreeWithoutParent(myChoice, y1)

                pTree.pruneSubTreeWithoutParent(y0)
                pTree.reconnectSubTreeWithoutParent(y0, s1)

                pTree.pruneSubTreeWithoutParent(s0)
                pTree.reconnectSubTreeWithoutParent(s0, y0)

        if dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()
            if y0Uncon and s0 != y1:
                print("The drawing above shows that Y extended")
            else:
                print("The drawing above shows that Y did not extend.")

        if oldRoot != pTree.root:
            pTree.reRoot(oldRoot, moveInternalName=False)
        if dbug:
            pTree.draw()
            print("The above is back to the original root.")

        if dbug:
            for n in pTree.nodes:
                if n.isLeaf:
                    while '_' in n.name:
                        n.name = n.name[:-3]
                elif not n.isLeaf:
                    n.name = None
                if n.br:
                    n.br.textDrawSymbol = '-'

        xRearranged = False
        if x0Uncon:
            if r0 is not x1:
                xRearranged = True
        yRearranged = False
        if y0Uncon:
            if s0 is not y1:
                yRearranged = True
        if xRearranged or yRearranged:
            theProposal.topologyChanged = True

        if dbug:
            print("-" * 20)
            print("xRearranged = %s" % xRearranged)
            print("yRearranged = %s" % yRearranged)
            if eA:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eA:
                        print("eA is from node", n.nodeNum)
            else:
                print("eA is None")
            if eX:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eX:
                        print("eX is from node", n.nodeNum)
            else:
                print("eX is None")
            if eY:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eY:
                        print("eY is from node", n.nodeNum)
            else:
                print("eY is None")
            print("x0Uncon is", x0Uncon)
            print("y0Uncon is", y0Uncon)


        # n.flag is set if the condLikes need recalculating.  Edges eA,
        # eX, and eY will have their bigPDecks recalculated, and all the
        # nodes below those will have their flags set (in
        # Chain.proposeSp()).  However, there are additional nodes that
        # need the flag set.
        if xRearranged:
            p = x1
            while 1:
                p.flag = 1
                p = p.parent
                if not p:
                    break
        if yRearranged:
            p = y1
            while 1:
                p.flag = 1
                p = p.parent
                if not p:
                    break

        # Do the branch length changes.  Taken nearly verbatim from Crux.
        # Thanks Jason!  The 3 edges to modify are eA, eX, and eY, which
        # were named above, as NodeBranch objects.

        lnProp = 0.0

        # Generate branch length multipliers and set new branch lengths.
        # eA.
        u = random.random()
        lnMA = etbrLambda * (u - 0.5)
        lnProp += lnMA
        mA = math.exp(lnMA)
        vA0 = eA.len
        vA1 = vA0 * mA
        eA.len = vA1
        eA.lenChanged = True
        # eX.
        if x0Uncon:
            u = random.random()
            lnMX = etbrLambda * (u - 0.5)
            lnProp += lnMX
            mX = math.exp(lnMX)
            vX0 = eX.len
            vX1 = vX0 * mX
            eX.len = vX1
            eX.lenChanged = True
        # eY.
        if y0Uncon:
            u = random.random()
            lnMY = etbrLambda * (u - 0.5)
            lnProp += lnMY
            mY = math.exp(lnMY)
            vY0 = eY.len
            vY1 = vY0 * mY
            eY.len = vY1
            eY.lenChanged = True

        # More from Crux ...
        # The prior ratio is the product of the prior ratios for each modified
        # branch length.  The number of internal branches does not change
        # (though the number of polytomies may change), so the topology prior
        # ratio is always 1.
        if theProposal.brLenPriorType == 'exponential':
            lnPrior = -theProposal.brLenPriorLambda * (vA1 - vA0)
            if x0Uncon:
                lnPrior += -theProposal.brLenPriorLambda * (vX1 - vX0)
            if y0Uncon:
                lnPrior += -theProposal.brLenPriorLambda * (vY1 - vY0)
        elif theProposal.brLenPriorType == 'uniform':
            lnPrior = 0.0

        # The proposal ratio is the product of the proposal ratios for
        # extension of each end of eA, as well as the branch multipliers.  The
        # ratio is 1 for the constrained/constrained and
        # unconstrained/unconstrained extension cases.
        #
        # lnProp for branch len changes were done above.

        # extension on the X-side
        if x0Uncon and r0 is not x1:  # has there been a rearrangement on the X-side?
            if y0Uncon:
                if not r0Uncon:
                    # x-side constrained, y-side unconstrained
                    lnProp += math.log(1.0 - etbrPExt)
            else:           # so y-side is constrained
                if r0Uncon: 
                    # x-side unconstrained, y-side constrained
                    lnProp += math.log(1.0 / (1.0 - etbrPExt))

        # extension on the Y-side
        if y0Uncon and s0 is not y1:  # has there been a rearrangement on the Y-side?
            if x0Uncon:
                if not s0Uncon:
                    # y-side constrained, x-side unconstrained
                    lnProp += math.log(1.0 - etbrPExt)
            else:           # so y-side is constrained
                if s0Uncon: 
                    # y-side unconstrained, x-side constrained
                    lnProp += math.log(1.0 / (1.0 - etbrPExt))

        self.propTree.preAndPostOrderAreValid = False

        self.logPriorRatio = lnPrior
        self.logProposalRatio = lnProp

        # Check if we have a new combo of comp and rMatrix.  This might be
        # more efficient if I cleverly only look at the affected nodes.
        # if 0:
        #     if theProposal.topologyChanged:
        #         for n in pTree.iterNodesNoRoot():
        #             for pNum in range(pTree.model.nParts):
        #                 theCompNum = n.parts[pNum].compNum
        #                 theRMatrixNum = n.br.parts[pNum].rMatrixNum
        #                 if pTree.model.parts[pNum].bQETneedsReset[theCompNum][theRMatrixNum]:
        #                     pf.p4_resetBQET(pTree.model.cModel, pNum, theCompNum, theRMatrixNum)

        # if 1:
        #     if theProposal.topologyChanged:
        #         for pNum in range(pTree.model.nParts):
        #             # print("\n")
        #             # print(pTree.model.parts[pNum].bQETneedsReset)
        #             for cNum in range(pTree.model.parts[pNum].nComps):
        #                 for rNum in range(pTree.model.parts[pNum].nRMatrices):
        #                     if pTree.model.parts[pNum].bQETneedsReset[cNum][rNum]:
        #                         pf.p4_resetBQET(pTree.model.cModel, pNum, cNum, rNum)
        #             # print("after p4_resetBQET, ...")
        #             # print(pTree.model.parts[pNum].bQETneedsReset)


        # This stuff below could probably be done more cleverly, but this
        # works.
        for n in pTree.iterNodesNoRoot():
            if math.fabs(n.br.len - n.br.oldLen) > 0.0000001:
                # print "Node %2i br len changed" % n.nodeNum
                assert n.br.lenChanged
            # else:
            #    if n.br.lenChanged:
            #        print "Node %2i lenChanged set, but its the same length." % n.nodeNum
            #        raise P4Error

            # If the branch has been inverted, we will want to recalculate
            # the bigPDecks, even if the length has not really changed.
            # Trigger that intent by setting lenChanged.
            if n.br.oldNode != n:
                # print "Node %2i branch: oldNode %2i" % (n.nodeNum,
                # n.br.oldNode.nodeNum)
                n.br.lenChanged = 1

    def proposeESPR_Blaise(self, theProposal):
        """SPR version of proposeETBR_Blaise

        See the docstring for proposeETBR_Blaise.  This is just the first half
        of the eTBR move.

        """
        
        # doAbort is set if brLens are too long or too short, or if a
        # constraint is violated.
        gm = ['Chain.proposeESPR_Blaise()']

        # Now works with constraints.
        theProposal.topologyChanged = 0
        theProposal.doAbort = False
        pTree = self.propTree
        dbug = False

        for n in pTree.iterNodesNoRoot():
            n.br.oldLen = n.br.len
            n.br.oldNode = n

        if 0 and self.mcmc.gen == 217:
            dbug = True
            print("------------- eSPR gen %i -----------------" % self.mcmc.gen)
            if 0:
                currentLogLike = self.propTree.logLike
                self.propTree.calcLogLike(verbose=0)  # with _commonCStuff()
                theDiff = math.fabs(currentLogLike - self.propTree.logLike)
                if theDiff > 1.e-9:
                    gm.append("propTree like diff %f (%g)" %
                              (theDiff, theDiff))
                    raise P4Error(gm)

        oldRoot = pTree.root
        eA = None
        eX = None
        eY = None

        etbrLambda = theProposal.tuning
        etbrPExt = theProposal.etbrPExt

        if 0 and dbug:
            print("=" * 50)
            pTree.draw()
            print("starting with the tree above.")

        # Choose a node, not the root.  In crux the choice for the original
        # branch can be any branch at all, including leaf branches.  Here, since
        # I choose a random node that is not the root, then the parent, the X0
        # end, will never be a leaf, but the Y0 end may be.  It will have edge
        # eA in Jason's diagram.  y0 will be the asterisk node in Jason's
        # diagram.  In eTBR, it (y0) may be extended, below, as the second of
        # the two eSPR moves.

        y0 = None
        while not y0:
            nNum = random.choice(pTree.preOrder)
            if nNum != var.NO_ORDER and nNum != pTree.root.nodeNum:
                y0 = pTree.node(nNum)
        x0 = y0.parent

        # The identity of x0 and y0 should be random.  That is, the first
        # direction of extension, that is the first eSPR, should be random.  The
        # way it is coded it below is that the x-direction is extended first,
        # then the y-direction.  However, the way I have it above, it is x0 =
        # y0.parent.  To make it properly random, switch ends with probability
        # 0.5.
        myRan = random.random()
        if myRan >= 0.5:
            x0 = y0
            y0 = x0.parent
            
        if 0:                                               # for debugging and demos
            y0 = pTree.node(2)
            x0 = pTree.node(1)
            #y0 = pTree.node(3)
            #x0 = pTree.node(2)

        if dbug:
            print("y0 is node %i" % y0.nodeNum)
            y0.br.textDrawSymbol = '='
            if y0.name:
                y0.name += '_y0'
            else:
                y0.name = 'y0'
            if x0.name:
                x0.name += '_x0'
            else:
                x0.name = 'x0'

        # Name the edge here, used below when we modify the br.len
        if x0 == y0.parent:
            eA = y0.br
        else:
            eA = x0.br

        # #########
        # Extend x
        # #########
        # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
        r1 = x0

        #  If x0 is not a leaf, it is 'unconstrained'
        x0Degree = pTree.getDegree(x0)
        if x0Degree == 1:
            x0Uncon = False   # a leaf
        else:
            x0Uncon = True
        if x0Uncon:
            myRan = random.randrange(x0Degree - 1)
            if y0.parent == x0:  # x-side is down
                r0 = pTree.nextNode(y0, x0)
            elif x0.parent == y0: # x-side is up
                r0 = pTree.nextNode(x0, x0)
            else:
                gm.append("This shouldn't happen.")
                raise P4Error(gm)

            for i in range(myRan):
                r0 = pTree.nextNode(r0, x0)
            if r0 == x0:  # points down 
                r0 = r0.parent

            if 0:
                r0 = pTree.node(0)                            # More for debugging and demos

            x1 = r0
            if dbug:
                if x1.name:
                    x1.name += '_x1'
                else:
                    x1.name = 'x1'

            # So we are set up like this ...
            # 
            #        eA        eX
            #   y0--------x0--------x1
            #             r1        r0
            #
            # ... and below, if r0 is unconstrained, 
            # it may turn into this (n? is just some other node) ---
            #
            #        eA        eX        eR
            #   y0--------x0--------x1--------n?
            #                       r1        r0
            #
            # and then maybe
            #        eA        eX                  eR
            #   y0--------x0--------x1--------n?--------n?
            #                                 r1        r0
            
            # We name the edge here, and use the name below when we modify the br.len
            if x0.parent == x1:
                eX = x0.br
            elif x1.parent == x0:
                eX = x1.br
            else:
                raise P4Error("This should not happen")

            #  if r0 is not a leaf, it is 'unconstrained'
            r0Degree = pTree.getDegree(r0)
            r0Uncon = (r0Degree > 1)
            while r0Uncon:
                # randomly determine whether to extend
                myRan = random.random()
                if etbrPExt < myRan:
                    break
                # choose an edge (direction) to go, and step to the next node
                myRan = random.randrange(r0Degree - 1)
                # We now call nextNode(spoke, hub) where hub is r0.  If r0
                # is the parent of r1, then the spoke of nextNode(spoke,
                # hub) is r1.  But if its the other way around, if r1 is
                # the parent of r0, then spoke is r0.
                if r1.parent == r0:
                    r0new = pTree.nextNode(r1, r0)
                elif r0.parent == r1:
                    r0new = pTree.nextNode(r0, r0)
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)

                for i in range(myRan):
                    r0new = pTree.nextNode(r0new, r0)
                if r0new == r0:
                    r0new = r0.parent
                r1 = r0
                r0 = r0new
                r0Degree = pTree.getDegree(r0)
                r0Uncon = (r0Degree > 1)
            if dbug:
                if r1.name:
                    r1.name += '_r1'
                else:
                    r1.name = 'r1'
                if r0.name:
                    r0.name += '_r0'
                else:
                    r0.name = 'r0'
                # pTree.draw()

            # Perform rearrangement unless it would be a no-op.  It would be a
            # no-op if r0 was still x1.
            if r0 == x1:
                # We did not extend, at all.
                if dbug:
                    print("No extension from x1 was done (because x1=r0), so no rearrangement on the x side.")
                pass
            else:
                # Do the rearrangement.  
                pTree.reRoot(r1, moveInternalName=False)
                assert x0.parent == x1
                assert y0.parent == x0
                assert r0.parent == r1

                if dbug:
                    x0.br.textDrawSymbol = 'X'
                    pTree.setPreAndPostOrder()
                    pTree.draw()
                    print("The drawing above is just before the rearrangement.")

                # Get children of x0 that are not y0
                ch = [n for n in x0.iterChildren() if n != y0]
                assert ch
                if len(ch) == 1:
                    myChoice = ch[0]
                else:
                    myChoice = random.choice(ch)

                # Three spr moves.
                pTree.pruneSubTreeWithoutParent(
                    myChoice, allowSingleChildNode=True)
                pTree.reconnectSubTreeWithoutParent(myChoice, x1)

                pTree.pruneSubTreeWithoutParent(x0)
                pTree.reconnectSubTreeWithoutParent(x0, r1)

                pTree.pruneSubTreeWithoutParent(r0)
                pTree.reconnectSubTreeWithoutParent(r0, x0)

        if 1 and dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()
            if x0Uncon and r0 != x1:
                print("The drawing above shows that X extended")
            else:
                print("The drawing above shows that X did not extend.")


        if oldRoot != pTree.root:
            pTree.reRoot(oldRoot, moveInternalName=False)
        if dbug:
            pTree.draw()
            print("The above is back to the original root.")

            for n in pTree.nodes:
                if n.isLeaf:
                    while '_' in n.name:
                        n.name = n.name[:-3]
                elif not n.isLeaf:
                    n.name = None
                if n.br:
                    n.br.textDrawSymbol = '-'

        xRearranged = False
        if x0Uncon:
            if r0 != x1:
                xRearranged = True
        if xRearranged:
            theProposal.topologyChanged = True

        if dbug:
            print("-" * 20)
            print("xRearranged = %s" % xRearranged)
            if eA:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eA:
                        print("eA is from node", n.nodeNum)
            else:
                print("eA is None")
            if eX:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eX:
                        print("eX is from node", n.nodeNum)
            else:
                print("eX is None")
            print("x0Uncon is", x0Uncon)


        # n.flag is set if the condLikes need recalculating.  Edges eA,
        # eX, and eY will have their bigPDecks recalculated, and all the
        # nodes below those will have their flags set (in
        # Chain.proposeSp()).  However, there are additional nodes that
        # need the flag set.
        if xRearranged:
            p = x1
            while 1:
                p.flag = 1
                p = p.parent
                if not p:
                    break

        # Do the branch length changes.  Taken nearly verbatim from Crux.
        # Thanks Jason!  The 2 edges to modify are eA, and eX, which
        # were named above, as NodeBranch objects.

        lnProp = 0.0
        lnPrior = 0.0

        if 1:
            # Generate branch length multipliers and set new branch lengths.
            # eA.
            u = random.random()
            lnMA = etbrLambda * (u - 0.5)
            lnProp += lnMA
            mA = math.exp(lnMA)
            vA0 = eA.len
            vA1 = vA0 * mA
            eA.len = vA1
            eA.lenChanged = True
            # eX.
            if x0Uncon:
                u = random.random()
                lnMX = etbrLambda * (u - 0.5)
                lnProp += lnMX
                mX = math.exp(lnMX)
                vX0 = eX.len
                vX1 = vX0 * mX
                eX.len = vX1
                eX.lenChanged = True

            # More from Crux ...  The prior ratio is the product of the prior ratios
            # for each modified branch length.  The number of internal branches does
            # not change, so the topology prior ratio is always 1.
            if theProposal.brLenPriorType == 'exponential':
                lnPrior = -theProposal.brLenPriorLambda * (vA1 - vA0)
                if x0Uncon:
                    lnPrior += -theProposal.brLenPriorLambda * (vX1 - vX0)
            elif theProposal.brLenPriorType == 'uniform':
                lnPrior = 0.0

        # The proposal ratio is the product of the proposal ratios for
        # extension of the x- end of eA, as well as the branch multipliers.  The
        # ratio is 1 for the constrained/constrained and
        # unconstrained/unconstrained extension cases.
        #
        # nY0/nR0.

        # proposal ratio
        # For this, we need to know if y0 is a leaf (constrained) or not.
        y0Degree = pTree.getDegree(y0)
        if y0Degree == 1:
            y0Uncon = False   # a leaf
        else:
            y0Uncon = True

        # extension on the x-side
        if x0Uncon and r0 is not x1:  # has there been a rearrangement on the x-side?
            if y0Uncon:
                if not r0Uncon:
                    # x-side constrained, y-side unconstrained
                    lnProp += math.log(1.0 - etbrPExt)
            else:           # so y-side is constrained
                if r0Uncon: 
                    # x-side unconstrained, y-side constrained
                    lnProp += math.log(1.0 / (1.0 - etbrPExt))

        self.propTree.preAndPostOrderAreValid = False

        self.logPriorRatio = lnPrior
        self.logProposalRatio = lnProp

        # Check if we have a new combo of comp and rMatrix.  This might be
        # more efficient if I cleverly only look at the affected nodes.
        if 1:
            if theProposal.topologyChanged:
                for n in pTree.iterNodesNoRoot():
                    for pNum in range(pTree.model.nParts):
                        theCompNum = n.parts[pNum].compNum
                        theRMatrixNum = n.br.parts[pNum].rMatrixNum
                        if pTree.model.parts[pNum].bQETneedsReset[theCompNum][theRMatrixNum]:
                            pf.p4_resetBQET(
                                pTree.model.cModel, pNum, theCompNum, theRMatrixNum)

        # This stuff below could probably be done more cleverly, but this
        # works.
        for n in pTree.iterNodesNoRoot():
            if math.fabs(n.br.len - n.br.oldLen) > 0.0000001:
                # print "Node %2i br len changed" % n.nodeNum
                assert n.br.lenChanged
            # else:
            #    if n.br.lenChanged:
            #        print "Node %2i lenChanged set, but its the same length." % n.nodeNum
            #        raise P4Error

            # If the branch has been inverted, we will want to recalculate
            # the bigPDecks, even if the length has not really changed.
            # Trigger that intent by setting lenChanged.
            if n.br.oldNode != n:
                # print "Node %2i branch: oldNode %2i" % (n.nodeNum,
                # n.br.oldNode.nodeNum)
                n.br.lenChanged = 1

    def proposeETBR(self, theProposal):
        """Evans and Sullivan eTBR move

        Described in Evans, J and J Sullivan 2012.  Generalized Mixture Models
        for Molecular Phylogenetic Estimation.  Syst Biol 61:12-21
        http://sysbio.oxfordjournals.org/content/61/1/12.short

        Adapted from Jason Evans' excellent Crux v 1.2.0.

        Many thanks to JE, who wrote such clear code.  The Crux version
        came from Lakner et al, but was modified by JE so that it works on
        polytomies, and so that it works on leaf nodes.  Also many thanks
        to Blaise Li who convincingly pointed out that LOCAL was not enough.

        It does not work with constraints.
        This proposal is not reversible if there are polytomies.
        """

        # doAbort is set if brLens are too long or too short, or if a
        # constraint is violated.
        gm = ['Chain.proposeETBR()']

        if self.mcmc.constraints:
            gm.append(
                "Sorry, due to lazy programming, proposeETBR() does not work with constraints yet.")
            raise P4Error(gm)
        theProposal.topologyChanged = 0
        theProposal.doAbort = False
        pTree = self.propTree
        dbug = False

        if 0 and self.mcmc.gen == 434:
            dbug = True
            if 0:
                currentLogLike = self.propTree.logLike
                self.propTree.calcLogLike(verbose=0)  # with _commonCStuff()
                theDiff = math.fabs(currentLogLike - self.propTree.logLike)
                if theDiff > 1.e-9:
                    gm.append("propTree like diff %f (%g)" %
                              (theDiff, theDiff))
                    raise P4Error(gm)

        # The figure below is from Evans' Crux v 1.2.0, 
        # file crux-1.2.0/pkg/Crux/Mc3/Chain.pyx
        # in the method defined by cdef etbrPropose

        # For each end of eA, choose a random starting direction away from eA.
        # With probability etbrPExt, "extend" (move the end point past the next
        # node along that path.  Randomly choose a direction away from the
        # previous edge.  Continue iterative extension until extension failure
        # due to randomness (the "unconstrained" case), or a leaf edge is
        # reached (the "constrained" case).  Re-arrange the tree by extracting
        # *==eA==nX0==eX== and moving it to the extension point.  The following
        # figure depicts the transformation due to extending the left end of eA
        # to eR:
        #
        #         cH               cI                            cI            #
        #          \              /    ===>                        \           #
        #           \            /                                  \          #
        #           nR1--------nR0                                  nR0--cJ    #
        #    cG     /     eR     \                cH                /          #
        #      \   /              \    ===>         \            eR/           #
        #       \ /                cJ                \            /            #
        #   cF--nX1                                  nR1========nX0            #
        #        \\               cA          cG     /     eX     \\           #
        #         \\eX            /    ===>     \   /            eA\\          #
        #          \\     eA     /               \ /                \\         #
        #       cE--nX0=========*            cF--nX1--cC             *--cA     #
        #           / \          \               / \                /          #
        #          /   \          \    ===>     /   \              /           #
        #         cD    cC         cB          cE    cD          cB            #

        #=======================================================================


        oldRoot = pTree.root
        eA = None
        eX = None
        eY = None

        etbrLambda = theProposal.tuning
        etbrPExt = theProposal.etbrPExt

        if 1 and dbug:
            print("=" * 80)
            print("=" * 80)
            pTree.draw()
            print("starting with the tree above.")

        # Choose a node, not the root.  In crux the choice for the original
        # branch can be any branch at all, including leaf branches.  Here, since
        # I choose a random node that is not the root, then the parent, the X0
        # end, will never be a leaf, but the Y0 end may be.  It will have edge
        # eA in Jason's diagram.  y0 will be the asterisk node in Jason's
        # diagram.  In eTBR, it (y0) may be extended, below, as the second of
        # the two eSPR moves.

        y0 = None        # The "*" node in Jason's diagram
        while not y0:
            nNum = random.choice(pTree.preOrder)
            if nNum != var.NO_ORDER and nNum != pTree.root.nodeNum:
                y0 = pTree.node(nNum)
        x0 = y0.parent
        assert x0

        # The identity of x0 and y0 should be random.  That is, the first
        # direction of extension, that is the first eSPR, should be random.  The
        # way it is coded it below is that the x-direction is extended first,
        # then the y-direction.  However, the way I have it above, it is x0 =
        # y0.parent.  To make it properly random, switch ends with probability
        # 0.5.
        myRan = random.random()
        if myRan >= 0.5:
            x0 = y0
            y0 = x0.parent
            

        if dbug:
            print("y0 is node %i" % y0.nodeNum)
            if x0 == y0.parent:
                y0.br.textDrawSymbol = '='
            else:
                x0.br.textDrawSymbol = '='
            if y0.name:
                y0.name += '_y0'
            else:
                y0.name = 'y0'
            if x0.name:
                x0.name += '_x0'
            else:
                x0.name = 'x0'
            # pTree.draw()
        # Name the edge here, used below when we modify the br.len
        if x0 == y0.parent:
            eA = y0.br
        else:
            eA = x0.br

        # #########
        # Extend x
        # #########
        # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
        r1 = x0

        # If x0 is not a leaf, it is 'unconstrained'
        x0Degree = pTree.getDegree(x0)
        if x0Degree == 1:
            x0Uncon = False   # a leaf
        else:
            x0Uncon = True
        if x0Uncon:
            myRan = random.randrange(x0Degree - 1)
            if y0.parent == x0:  # x-side is down
                r0 = pTree.nextNode(y0, x0)
            elif x0.parent == y0: # x-side is up
                r0 = pTree.nextNode(x0, x0)
            else:
                gm.append("This shouldn't happen.")
                raise P4Error(gm)

            for i in range(myRan):
                r0 = pTree.nextNode(r0, x0)
            if r0 == x0:  # points down 
                r0 = r0.parent
            x1 = r0
            if dbug:
                if x1.name:
                    x1.name += '_x1'
                else:
                    x1.name = 'x1'

            # So we are set up like this ...
            # 
            #        eA        eX
            #   y0--------x0--------x1
            #             r1        r0
            #
            # ... and below, if r0 is unconstrained, 
            # it may turn into this (n? is just some other node) ---
            #
            #        eA        eX        eR
            #   y0--------x0--------x1--------n?
            #                       r1        r0
            #
            # and then maybe
            #        eA        eX                  eR
            #   y0--------x0--------x1--------n?--------n?
            #                                 r1        r0
            
            # We name the edge here, and use the name below when we modify the
            # br.len
            if x0.parent == x1:
                eX = x0.br
            elif x1.parent == x0:
                eX = x1.br
            else:
                raise P4Error("This should not happen")

            #  if r0 is not a leaf, it is 'unconstrained'
            r0Degree = pTree.getDegree(r0)
            r0Uncon = (r0Degree > 1)
            while r0Uncon:
                # randomly determine whether to extend
                myRan = random.random()
                if etbrPExt < myRan:
                    break
                # choose an edge (direction) to go, and step to the next node
                myRan = random.randrange(r0Degree - 1)
                # We now call nextNode(spoke, hub) where hub is r0.  If r0
                # is the parent of r1, then the spoke of nextNode(spoke,
                # hub) is r1.  But if its the other way around, if r1 is
                # the parent of r0, then spoke is r0.
                if r1.parent == r0:
                    r0new = pTree.nextNode(r1, r0)
                elif r0.parent == r1:
                    r0new = pTree.nextNode(r0, r0)
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)

                for i in range(myRan):
                    r0new = pTree.nextNode(r0new, r0)
                if r0new == r0:
                    r0new = r0.parent
                r1 = r0
                r0 = r0new
                r0Degree = pTree.getDegree(r0)
                r0Uncon = (r0Degree > 1)
            if dbug:
                if r1.name:
                    r1.name += '_r1'
                else:
                    r1.name = 'r1'
                if r0.name:
                    r0.name += '_r0'
                else:
                    r0.name = 'r0'
                #pTree.draw()

            # Perform rearrangement unless it would be a no-op.  It would be a
            # no-op if r0 was still x1.
            if r0 == x1:
                # We did not extend, at all.
                pass
            else:
                # Do the rearrangement.  
                pTree.reRoot(r1, moveInternalName=False)
                assert x0.parent == x1
                assert y0.parent == x0
                assert r0.parent == r1

                if dbug:
                    x0.br.textDrawSymbol = 'X'
                    pTree.setPreAndPostOrder()
                    pTree.draw()

                # Collapse the node between x0 and x1, ...
                theRightmostChild = x0.rightmostChild()
                theLeftSib = x0.leftSibling()
                if theLeftSib:
                    theLeftSib.sibling = x0.leftChild
                else:
                    x1.leftChild = x0.leftChild
                for n in x0.iterChildren():
                    n.parent = x1
                theRightmostChild.sibling = x0.sibling
                x0.wipe()  # needed?

                # ... and insert that node (x0) between r0 and r1.
                x0.parent = r1
                x0.leftChild = r0
                r0.parent = x0
                if r1.leftChild == r0:
                    r1.leftChild = x0
                else:
                    oldCh = r1.leftChild
                    while oldCh.sibling != r0:
                        oldCh = oldCh.sibling
                    oldCh.sibling = x0
                if r0.sibling:
                    x0.sibling = r0.sibling
                    r0.sibling = None

                if 1 and dbug:
                    pTree.setPreAndPostOrder()
                    pTree.draw()

                # this method returns y0 as well.
                pTree.pruneSubTreeWithoutParent(y0)
                pTree.reconnectSubTreeWithoutParent(y0, x0)

                # if oldRoot != pTree.root:
                #    pTree.reRoot(oldRoot, moveInternalName=False)

        if 1 and dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()
            if x0Uncon and r0 != x1:
                print("The drawing above shows that X extended")
            else:
                print("The drawing above shows that X did not extend.")

        # #########
        # Extend y
        # #########
        # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
        s1 = y0

        #  If y0 is not a leaf, it is 'unconstrained'
        y0Degree = pTree.getDegree(y0)
        if y0Degree == 1:
            y0Uncon = False    # If y0 is a leaf
        else:
            y0Uncon = True
        if y0Uncon:
            myRan = random.randrange(y0Degree - 1)
            if x0.parent == y0:  # y-side is down
                s0 = pTree.nextNode(x0, y0)
            elif y0.parent == x0:  # y-side is up
                s0 = pTree.nextNode(y0, y0)
            else:
                gm.append("This shouldn't happen.")
                raise P4Error(gm)

            for i in range(myRan):
                s0 = pTree.nextNode(s0, y0)
            if s0 == y0:   # points down
                s0 = s0.parent
            y1 = s0
            if dbug:
                if y1.name:
                    y1.name += '_y1'
                else:
                    y1.name = 'y1'

            # name the edge here, to be used below when we modify the br.len
            if y0.parent == y1:
                eY = y0.br
            elif y1.parent == y0:
                eY = y1.br
            else:
                raise P4Error("This should not happen (name eY)")

            # So we are set up like this ...
            #   nS1      nS0
            #   nY0==eY==nY1

            #  if s0 is not a leaf, it is 'unconstrained'
            s0Degree = pTree.getDegree(s0)
            s0Uncon = (s0Degree > 1)
            while s0Uncon:
                # randomly determine whether to extend
                myRan = random.random()
                if etbrPExt < myRan:
                    break
                # choose an edge (direction) to go, and step to the next node
                myRan = random.randrange(s0Degree - 1)
                # We now call nextNode(spoke, hub) where hub is s0.  If s0
                # is the parent of s1, then the spoke of nextNode(spoke,
                # hub) is s1.  But if its the other way around, if s1 is
                # the parent of s0, then spoke is s0.
                if s1.parent == s0:
                    s0new = pTree.nextNode(s1, s0)
                elif s0.parent == s1:
                    s0new = pTree.nextNode(s0, s0)
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)

                for i in range(myRan):
                    s0new = pTree.nextNode(s0new, s0)
                if s0new == s0:
                   s0new = s0.parent
                s1 = s0
                s0 = s0new
                s0Degree = pTree.getDegree(s0)
                s0Uncon = (s0Degree > 1)
            if dbug:
                if s1.name:
                    s1.name += '_s1'
                else:
                    s1.name = 's1'
                if s0.name:
                    s0.name += '_s0'
                else:
                    s0.name = 's0'
                if y1.parent == y0:
                    y1.br.textDrawSymbol = 'Y'
                elif y0.parent == y1:
                    y0.br.textDrawSymbol = 'Y'
                #pTree.draw()

            # Perform rearrangement unless it would be a no-op.  It would be a
            # no-op if s0 was still y1.
            if s0 == y1:
                # We did not extend, at all.
                pass
            else:
                # Do the rearrangement.  
                pTree.reRoot(s1, moveInternalName=False)
                assert y0.parent == y1
                assert x0.parent == y0
                assert s0.parent == s1

                if dbug:
                    y0.br.textDrawSymbol = 'Y'
                    pTree.setPreAndPostOrder()
                    pTree.draw()

                # Collapse the node between y0 and y1, ...
                theRightmostChild = y0.rightmostChild()
                theLeftSib = y0.leftSibling()
                if theLeftSib:
                    theLeftSib.sibling = y0.leftChild
                else:
                    y1.leftChild = y0.leftChild
                for n in y0.iterChildren():
                    n.parent = y1
                theRightmostChild.sibling = y0.sibling
                y0.wipe()  # needed?

                # ... and insert that node (y0) between s0 and s1.
                y0.parent = s1
                y0.leftChild = s0
                s0.parent = y0
                if s1.leftChild == s0:
                    s1.leftChild = y0
                else:
                    oldCh = s1.leftChild
                    while oldCh.sibling != s0:
                        oldCh = oldCh.sibling
                    oldCh.sibling = y0
                if s0.sibling:
                    y0.sibling = s0.sibling
                    s0.sibling = None

                if dbug:
                    pTree.setPreAndPostOrder()
                    pTree.draw()

                # this method returns x0 as well.
                pTree.pruneSubTreeWithoutParent(x0)
                pTree.reconnectSubTreeWithoutParent(x0, y0)

                if 0 and dbug:
                    pTree.setPreAndPostOrder()
                    pTree.draw()

        if oldRoot != pTree.root:
            pTree.reRoot(oldRoot, moveInternalName=False)

        if dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()
            if y0Uncon and s0 != y1:
                print("The drawing above shows that Y extended")
            else:
                print("The drawing above shows that Y did not extend.")

        if dbug:
            for n in pTree.nodes:
                if n.isLeaf:
                    while '_' in n.name:
                        n.name = n.name[:-3]
                elif not n.isLeaf:
                    n.name = None
                if n.br:
                    n.br.textDrawSymbol = '-'

        xRearranged = False
        if x0Uncon:
            if r0 != x1:
                xRearranged = True
        yRearranged = False
        if y0Uncon:
            if s0 != y1:
                yRearranged = True
        if xRearranged or yRearranged:
            theProposal.topologyChanged = True

        if dbug:
            print("-" * 20)
            print("xRearranged = %s" % xRearranged)
            print("yRearranged = %s" % yRearranged)
            if eA:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eA:
                        print("eA is from node", n.nodeNum)
            else:
                print("eA is None")
            if eX:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eX:
                        print("eX is from node", n.nodeNum)
            else:
                print("eX is None")
            if eY:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eY:
                        print("eY is from node", n.nodeNum)
            else:
                print("eY is None")
            print("x0Uncon is", x0Uncon)
            print("y0Uncon is", y0Uncon)

        # n.flag is set if the condLikes need recalculating.  Edges eA,
        # eX, and eY will have their bigPDecks recalculated, and all the
        # nodes below those will have their flags set (in
        # Chain.proposeSp()).  However, there are additional nodes that
        # need the flag set.

        # Overkill
        if 0:
            for n in pTree.iterNodesNoRoot():
                if not n.isLeaf:
                    n.flag = 1
                n.br.lenChanged = True
            pTree.root.flag = 1
        

        if 1:
            if x0.br:
                x0.br.lenChanged = True
                n = x0
                while 1:
                    if not n.isLeaf:
                        n.flag = 1
                    n = n.parent
                    if not n:
                        break
                
            if y0.br:
                y0.br.lenChanged = True
                n = y0
                while 1:
                    if not n.isLeaf:
                        n.flag = 1
                    n = n.parent
                    if not n:
                        break
            if x0Uncon:
                #for n in [x1, r0, r1]:
                for n in [x1]:
                    if n and n.br:
                        n.br.lenChanged = True
                    while 1:
                        if not n.isLeaf:
                            n.flag = 1
                        n = n.parent
                        if not n:
                            break
            if y0Uncon:
                #for n in [y1, s0, s1]:
                for n in [y1]:
                    if n and n.br:
                        n.br.lenChanged = True
                    while 1:
                        if not n.isLeaf:
                            n.flag = 1
                        n = n.parent
                        if not n:
                            break

            if 0 and self.mcmc.gen == 434:
                for n in pTree.iterNodesNoRoot():
                    print("%2i  %5s  %i" % (n.nodeNum, n.br.lenChanged, n.flag))
                n = pTree.root
                print("%2i  %5s  %i" % (n.nodeNum, "-", n.flag))
                #for nNum in [11]:
                #    n = pTree.node(nNum)
                #    n.br.lenChanged = True


        # odd things with rearrangement due to the re-rooting, so we need to set
        # node.br.lenChanged from x1 down to r1, and from y1 down to s1.

        if xRearranged:
            if x0.isAncestorOf(r0):
                if r1.isAncestorOf(x1):
                    p = x1
                    while 1:
                        p.br.lenChanged = True
                        p = p.parent
                        if p == r1:
                            break
        if yRearranged:
            if y0.isAncestorOf(s0):
                if s1.isAncestorOf(y1):
                    p = y1
                    while 1:
                        p.br.lenChanged = True
                        p = p.parent
                        if p == s1:
                            break

        # Do the branch length changes.  Taken nearly verbatim from Crux.
        # Thanks Jason!  The 3 edges to modify are eA, eX, and eY, which
        # were named above, as NodeBranch objects.

        lnPrior = 0.0
        lnProp = 0.0

        if 1:
            # Generate branch length multipliers and set new branch lengths.
            # eA.
            u = random.random()
            lnMA = etbrLambda * (u - 0.5)
            lnProp += lnMA
            mA = math.exp(lnMA)
            vA0 = eA.len
            vA1 = vA0 * mA
            eA.len = vA1
            eA.lenChanged = True
            # eX.
            if x0Uncon:
                u = random.random()
                lnMX = etbrLambda * (u - 0.5)
                lnProp += lnMX
                mX = math.exp(lnMX)
                vX0 = eX.len
                vX1 = vX0 * mX
                eX.len = vX1
                eX.lenChanged = True
            # eY.
            if y0Uncon:
                u = random.random()
                lnMY = etbrLambda * (u - 0.5)
                lnProp += lnMY
                mY = math.exp(lnMY)
                vY0 = eY.len
                vY1 = vY0 * mY
                eY.len = vY1
                eY.lenChanged = True

            # More from Crux ...
            # The prior ratio is the product of the prior ratios for each modified
            # branch length.  The number of internal branches does not change
            # (though the number of polytomies may change), so the topology prior
            # ratio is always 1.
            lnPrior = -theProposal.brLenPriorLambda * (vA1 - vA0)
            if x0Uncon:
                lnPrior += -theProposal.brLenPriorLambda * (vX1 - vX0)
            if y0Uncon:
                lnPrior += -theProposal.brLenPriorLambda * (vY1 - vY0)

        # The proposal ratio is the product of the proposal ratios for
        # extension of each end of eA, as well as the branch multipliers.  The
        # ratio is 1 for the constrained/constrained and
        # unconstrained/unconstrained extension cases.
        #

        # In crux, the code was as shown here (with similar names)
        # y0/r0.
        # if x0Uncon and r0 is not x1:
        #     if y0Uncon:
        #         if not r0Uncon:
        #             lnProp += math.log(1.0 - etbrPExt)
        #     elif r0Uncon:
        #         lnProp += math.log(1.0 / (1.0 - etbrPExt))
        # # x0/s0.
        # if y0Uncon and s0 is not y1:
        #     if x0Uncon:
        #         if not s0Uncon:
        #             lnProp += math.log(1.0 - etbrPExt)
        #     elif s0Uncon:
        #         lnProp += math.log(1.0 / (1.0 - etbrPExt))
        #
        # The code above is saying ---
        # if x0 is an internal node and it extended on the x-side:
        #     if y0 is an internal node:
        #         if r0 is a leaf:  # so x-side is constrained
        #             # x-side is constrained, y is not
        #             proposalRatio *= 1 - p_ext
        #     elif r0 is an internal node:   # so x-side is unconstrained
        #         # we know that y-side is a leaf, so it is constrained
        #         # x is unconstrained, y is constrained
        #         proposalRatio *= 1/(1-p_ext)
        # if y0 is an internal node and it extended on the y-side:
        #     if x0 is an internal node:
        #         if s0 is a leaf:  # so y-side is constrained
        #             # y-side constrained, x-side unconstrained
        #             proposalRatio *= 1 - p_ext
        #     elif s0 is an internal node:  # so y-side is unconstrained
        #         # We know that the x-side is a leaf, so it is constrained
        #         # y-side is unconstrained, x-side constrained
        #         proposalRatio *= 1/(1-p_ext)
       
        # The code above works, but I do not like the obfuscated logic, with no
        # obvious "else" clause.  So I recode it a bit more simply, as

        # y0/r0.
        if x0Uncon and r0 is not x1:
            if y0Uncon:
                if not r0Uncon:
                    # x-side constrained, y-side unconstrained
                    lnProp += math.log(1.0 - etbrPExt)
            else:           # so y-side is constrained
                if r0Uncon: # x-side is not
                    # x-side unconstrained, y-side constrained
                    lnProp += math.log(1.0 / (1.0 - etbrPExt))
        # x0/s0.
        if y0Uncon and s0 is not y1:
            if x0Uncon:
                if not s0Uncon:
                    # y-side constrained, x-side unconstrained
                    lnProp += math.log(1.0 - etbrPExt)
            else:            # x-side is constrained
                if s0Uncon:  # y-side is not
                    # y-side unconstrained, x-side is constrained
                    lnProp += math.log(1.0 / (1.0 - etbrPExt))

        self.propTree.preAndPostOrderAreValid = False

        self.logPriorRatio = lnPrior
        self.logProposalRatio = lnProp

        # Check if we have a new combo of comp and rMatrix.  This might be
        # more efficient if I cleverly only look at the affected nodes.
        if 1:
            if theProposal.topologyChanged:
                for n in pTree.iterNodesNoRoot():
                    for pNum in range(pTree.model.nParts):
                        theCompNum = n.parts[pNum].compNum
                        theRMatrixNum = n.br.parts[pNum].rMatrixNum
                        if pTree.model.parts[pNum].bQETneedsReset[theCompNum][theRMatrixNum]:
                            pf.p4_resetBQET(
                                pTree.model.cModel, pNum, theCompNum, theRMatrixNum)

    def proposeESPR(self, theProposal):
        """Evans and Sullivan eSPR move

        See the description for proposeETBR.  This is just the first half of an
        eTBR.

        This seems to be working.  But there is no user interface for it.
        """
        # doAbort is set if brLens are too long or too short, or if a
        # constraint is violated.
        gm = ['Chain.proposeESPR()']

        if self.mcmc.constraints:
            gm.append(
                "Sorry, due to lazy programming, proposeSPR() does not work with constraints yet.")
            raise P4Error(gm)
        theProposal.topologyChanged = 0
        theProposal.doAbort = False
        pTree = self.propTree
        dbug = False

        if 0 and self.mcmc.gen == 434:
            dbug = True
            if 0:
                currentLogLike = self.propTree.logLike
                self.propTree.calcLogLike(verbose=0)  # with _commonCStuff()
                theDiff = math.fabs(currentLogLike - self.propTree.logLike)
                if theDiff > 1.e-9:
                    gm.append("propTree like diff %f (%g)" %
                              (theDiff, theDiff))
                    raise P4Error(gm)

        oldRoot = pTree.root
        eA = None
        eX = None
        eY = None

        etbrLambda = theProposal.tuning
        etbrPExt = theProposal.etbrPExt

        if 1 and dbug:
            print("=" * 80)
            print("=" * 80)
            pTree.draw()
            print("starting with the tree above.")

        # Choose a node, not the root.  In crux the choice for the original
        # branch can be any branch at all, including leaf branches.  Here, since
        # I choose a random node that is not the root, then the parent, the X0
        # end, will never be a leaf, but the Y0 end may be.  It will have edge
        # eA in Jason's diagram.  y0 will be the asterisk node in Jason's
        # diagram.  

        y0 = None        # The "*" node in Jason's diagram
        while not y0:
            nNum = random.choice(pTree.preOrder)
            if nNum != var.NO_ORDER and nNum != pTree.root.nodeNum:
                y0 = pTree.node(nNum)
        x0 = y0.parent
        assert x0

        # The identity of x0 and y0 should be random.  That is, the first
        # direction of extension, and the first eSPR, should be random.  The way
        # I have it, it is always x0 = y0.parent, and x0 is extended (down)
        # first, and then y0 is extended up.  So randomly switch ends.

        myRan = random.random()
        if myRan >= 0.5:
            x0 = y0
            y0 = x0.parent

        if 0:                                              # for debugging and demos
            y0 = pTree.node(2)
            x0 = pTree.node(1)

        if dbug:
            print("y0 is node %i" % y0.nodeNum)
            if x0 == y0.parent:
                y0.br.textDrawSymbol = '='
            else:
                x0.br.textDrawSymbol = '='
            if y0.name:
                y0.name += '_y0'
            else:
                y0.name = 'y0'
            if x0.name:
                x0.name += '_x0'
            else:
                x0.name = 'x0'
            # pTree.draw()
        # Name the edge here, used below when we modify the br.len
        if x0 == y0.parent:
            eA = y0.br
        else:
            eA = x0.br

        # #########
        # Extend x
        # #########
        # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
        r1 = x0

        # If x0 is not a leaf, it is 'unconstrained'
        x0Degree = pTree.getDegree(x0)
        if x0Degree == 1:
            x0Uncon = False   # a leaf
        else:
            x0Uncon = True
        if x0Uncon:
            myRan = random.randrange(x0Degree - 1)
            if y0.parent == x0:  # x-side is down
                r0 = pTree.nextNode(y0, x0)
            elif x0.parent == y0: # x-side is up
                r0 = pTree.nextNode(x0, x0)
            else:
                gm.append("This shouldn't happen.")
                raise P4Error(gm)

            for i in range(myRan):
                r0 = pTree.nextNode(r0, x0)
            if r0 == x0:  # points down 
                r0 = r0.parent
            if 0:
                r0 = pTree.node(0)                            # more for debugging and demo
            x1 = r0
            if dbug:
                if x1.name:
                    x1.name += '_x1'
                else:
                    x1.name = 'x1'

            # So we are set up like this ...
            # 
            #        eA        eX
            #   y0--------x0--------x1
            #             r1        r0
            #
            # ... and below, if r0 is unconstrained, 
            # it may turn into this (n? is just some other node) ---
            #
            #        eA        eX        eR
            #   y0--------x0--------x1--------n?
            #                       r1        r0
            #
            # and then maybe
            #        eA        eX                  eR
            #   y0--------x0--------x1--------n?--------n?
            #                                 r1        r0
            
            # We name the edge here, and use the name below when we modify the
            # br.len
            if x0.parent == x1:
                eX = x0.br
            elif x1.parent == x0:
                eX = x1.br
            else:
                raise P4Error("This should not happen")

            #  if r0 is not a leaf, it is 'unconstrained'
            r0Degree = pTree.getDegree(r0)
            r0Uncon = (r0Degree > 1)
            while r0Uncon:
                # randomly determine whether to extend
                myRan = random.random()
                if etbrPExt < myRan:
                    break
                # choose an edge (direction) to go, and step to the next node
                myRan = random.randrange(r0Degree - 1)
                # We now call nextNode(spoke, hub) where hub is r0.  If r0
                # is the parent of r1, then the spoke of nextNode(spoke,
                # hub) is r1.  But if its the other way around, if r1 is
                # the parent of r0, then spoke is r0.
                if r1.parent == r0:
                    r0new = pTree.nextNode(r1, r0)
                elif r0.parent == r1:
                    r0new = pTree.nextNode(r0, r0)
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)

                for i in range(myRan):
                    r0new = pTree.nextNode(r0new, r0)
                if r0new == r0:
                    r0new = r0.parent
                r1 = r0
                r0 = r0new

                r0Degree = pTree.getDegree(r0)
                r0Uncon = (r0Degree > 1)
            if dbug:
                if r1.name:
                    r1.name += '_r1'
                else:
                    r1.name = 'r1'
                if r0.name:
                    r0.name += '_r0'
                else:
                    r0.name = 'r0'
                #pTree.draw()

            # Perform rearrangement unless it would be a no-op.  It would be a
            # no-op if r0 was still x1.
            if r0 == x1:
                # We did not extend, at all.
                pass
            else:
                # Do the rearrangement.  
                pTree.reRoot(r1, moveInternalName=False)
                assert x0.parent == x1
                assert y0.parent == x0
                assert r0.parent == r1

                if dbug:
                    x0.br.textDrawSymbol = 'X'
                    pTree.setPreAndPostOrder()
                    pTree.draw()

                # Collapse the node between x0 and x1, ...
                theRightmostChild = x0.rightmostChild()
                theLeftSib = x0.leftSibling()
                if theLeftSib:
                    theLeftSib.sibling = x0.leftChild
                else:
                    x1.leftChild = x0.leftChild
                for n in x0.iterChildren():
                    n.parent = x1
                theRightmostChild.sibling = x0.sibling
                x0.wipe()  # needed?

                # ... and insert that node (x0) between r0 and r1.
                x0.parent = r1
                x0.leftChild = r0
                r0.parent = x0
                if r1.leftChild == r0:
                    r1.leftChild = x0
                else:
                    oldCh = r1.leftChild
                    while oldCh.sibling != r0:
                        oldCh = oldCh.sibling
                    oldCh.sibling = x0
                if r0.sibling:
                    x0.sibling = r0.sibling
                    r0.sibling = None

                if 1 and dbug:
                    pTree.setPreAndPostOrder()
                    pTree.draw()

                # this method returns y0 as well.
                pTree.pruneSubTreeWithoutParent(y0)
                pTree.reconnectSubTreeWithoutParent(y0, x0)

                # if oldRoot != pTree.root:
                #    pTree.reRoot(oldRoot, moveInternalName=False)

        if 1 and dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()
            if x0Uncon and r0 != x1:
                print("The drawing above shows that X extended")
            else:
                print("The drawing above shows that X did not extend.")


        if oldRoot != pTree.root:
            pTree.reRoot(oldRoot, moveInternalName=False)

        if dbug:
            pTree.draw()
            for n in pTree.nodes:
                if n.isLeaf:
                    while '_' in n.name:
                        n.name = n.name[:-3]
                elif not n.isLeaf:
                    n.name = None
                if n.br:
                    n.br.textDrawSymbol = '-'

        xRearranged = False
        if x0Uncon:
            if r0 != x1:
                xRearranged = True
        if xRearranged:
            theProposal.topologyChanged = True

        if dbug:
            print("-" * 20)
            print("xRearranged = %s" % xRearranged)
            if eA:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eA:
                        print("eA is from node", n.nodeNum)
            else:
                print("eA is None")
            if eX:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eX:
                        print("eX is from node", n.nodeNum)
            else:
                print("eX is None")
            print("x0Uncon is", x0Uncon)

        # n.flag is set if the condLikes need recalculating.  Edges eA,
        # eX, will have their bigPDecks recalculated, and all the
        # nodes below those will have their flags set (in
        # Chain.proposeSp()).  However, there are additional nodes that
        # need the flag set.

        # Overkill
        if 0:
            for n in pTree.iterNodesNoRoot():
                if not n.isLeaf:
                    n.flag = 1
                n.br.lenChanged = True
            pTree.root.flag = 1
        

        if 1:
            if x0.br:
                x0.br.lenChanged = True
                n = x0
                while 1:
                    if not n.isLeaf:
                        n.flag = 1
                    n = n.parent
                    if not n:
                        break
                
            if x0Uncon:
                #for n in [x1, r0, r1]:
                for n in [x1]:
                    if n and n.br:
                        n.br.lenChanged = True
                    while 1:
                        if not n.isLeaf:
                            n.flag = 1
                        n = n.parent
                        if not n:
                            break

            if 0 and self.mcmc.gen == 253:
                for n in pTree.iterNodesNoRoot():
                    print("%2i  %5s  %i" % (n.nodeNum, n.br.lenChanged, n.flag))
                n = pTree.root
                print("%2i  %5s  %i" % (n.nodeNum, "-", n.flag))
                #for nNum in [11]:
                #    n = pTree.node(nNum)
                #    n.br.lenChanged = True


        # odd things with rearrangement due to the re-rooting, so we need to set
        # node.br.lenChanged from x1 down to r1, and from y1 down to s1.

        if xRearranged:
            if x0.isAncestorOf(r0):
                if r1.isAncestorOf(x1):
                    p = x1
                    while 1:
                        p.br.lenChanged = True
                        p = p.parent
                        if p == r1:
                            break

        # Do the branch length changes.  Taken nearly verbatim from Crux.
        # Thanks Jason!  The edges to modify are eA, and eX, which
        # were named above, as NodeBranch objects.

        lnPrior = 0.0
        lnProp = 0.0

        if 1:
            # Generate branch length multipliers and set new branch lengths.
            # eA.
            u = random.random()
            lnMA = etbrLambda * (u - 0.5)
            lnProp += lnMA
            mA = math.exp(lnMA)
            vA0 = eA.len
            vA1 = vA0 * mA
            eA.len = vA1
            eA.lenChanged = True
            # eX.
            if x0Uncon:
                u = random.random()
                lnMX = etbrLambda * (u - 0.5)
                lnProp += lnMX
                mX = math.exp(lnMX)
                vX0 = eX.len
                vX1 = vX0 * mX
                eX.len = vX1
                eX.lenChanged = True

            # More from Crux ...
            # The prior ratio is the product of the prior ratios for each modified
            # branch length.  The number of internal branches does not change
            # (though the number of polytomies may change), so the topology prior
            # ratio is always 1.

            # This needs fixing to accommodate uniform/exponential brLenPriorType. 
            lnPrior = -theProposal.brLenPriorLambda * (vA1 - vA0)
            if x0Uncon:
                lnPrior += -theProposal.brLenPriorLambda * (vX1 - vX0)

        # proposal ratio
        # For this, we need to know if y0 is a leaf (constrained) or not.
        y0Degree = pTree.getDegree(y0)
        if y0Degree == 1:
            y0Uncon = False   # a leaf
        else:
            y0Uncon = True


        # y0/r0.
        if x0Uncon and r0 is not x1:
            if y0Uncon:
                if not r0Uncon:
                    # x-side constrained, y-side unconstrained
                    lnProp += math.log(1.0 - etbrPExt)
            else:           # so y-side is constrained
                if r0Uncon: # x-side is not
                    # x-side unconstrained, y-side constrained
                    lnProp += math.log(1.0 / (1.0 - etbrPExt))

        self.propTree.preAndPostOrderAreValid = False

        self.logPriorRatio = lnPrior
        self.logProposalRatio = lnProp

        # Check if we have a new combo of comp and rMatrix.  This might be
        # more efficient if I cleverly only look at the affected nodes.
        if 1:
            if theProposal.topologyChanged:
                for n in pTree.iterNodesNoRoot():
                    for pNum in range(pTree.model.nParts):
                        theCompNum = n.parts[pNum].compNum
                        theRMatrixNum = n.br.parts[pNum].rMatrixNum
                        if pTree.model.parts[pNum].bQETneedsReset[theCompNum][theRMatrixNum]:
                            pf.p4_resetBQET(
                                pTree.model.cModel, pNum, theCompNum, theRMatrixNum)



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
        # if self.mcmc.constraints:
        #    print "checkSplitKeys() at the end of polytomy"
        #    self.propTree.checkSplitKeys()

    def proposeAddEdge(self, theProposal):
        gm = ["Chain.proposeAddEdge()"]
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

        # Choose a branch length for newNode.  See LewisHolderHolsinger eqn 6.
        newNode.br.len = - \
            (1.0 / theProposal.brLenPriorLambda) * \
            math.log(1. - random.random())
        if newNode.br.len < var.BRLEN_MIN or newNode.br.len > var.BRLEN_MAX:
            safety = 0
            while newNode.br.len < var.BRLEN_MIN or newNode.br.len > var.BRLEN_MAX:
                newNode.br.len = - \
                    (1.0 / theProposal.brLenPriorLambda) * \
                    math.log(1. - random.random())
                safety += 1
                if safety > 20:
                    gm.append(
                        "Unable to find a good branch length for the new edge.")
                    gm.append("Probably a programming error.")
                    raise P4Error(gm)
        if var.doMcmcSp:
            newNode.br.lenChanged = True

        # Calculate the rawSplitKey and splitKey.
        if self.mcmc.constraints:
            children = [n for n in newNode.iterChildren()]
            x = children[0].br.rawSplitKey
            for n in children[1:]:
                y = n.br.rawSplitKey
                x = x | y  # '|' is bitwise "OR".
            newNode.br.rawSplitKey = x
            # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
            if 1 & newNode.br.rawSplitKey:
                if self.mcmc.constraints:
                    # "^" is xor, a bit-flipper.
                    newNode.br.splitKey = self.mcmc.constraints.allOnes ^ newNode.br.rawSplitKey
                else:
                    allOnes = 2 ** (self.propTree.nTax) - 1
                    newNode.br.splitKey = allOnes ^ newNode.br.rawSplitKey
            else:
                newNode.br.splitKey = newNode.br.rawSplitKey

        # Its a newly-added node, possibly in a new context.  We need to
        # deal with model stuff if it isHet.  The model.isHet if any part
        # isHet.
        if pTree.model.isHet:
            for pNum in range(pTree.model.nParts):
                mp = pTree.model.parts[pNum]
                if mp.isHet:
                    # We want to make the model prams appropriate, and not
                    # just use what was left over from the last time the
                    # node was in the tree.
                    if mp.nComps > 1:
                        # Pick a child, and use the same comp.
                        oneChildNum = random.choice(newChildrenNodeNums)
                        newNode.parts[pNum].compNum = pTree.nodes[
                            oneChildNum].parts[pNum].compNum
                        mp.comps[newNode.parts[pNum].compNum].nNodes += 1
                    if mp.nRMatrices > 1:
                        oneChildNum = random.choice(newChildrenNodeNums)
                        newNode.br.parts[pNum].rMatrixNum = pTree.nodes[
                            oneChildNum].br.parts[pNum].rMatrixNum
                        mp.rMatrices[
                            newNode.br.parts[pNum].rMatrixNum].nNodes += 1
                    if mp.nGdasrvs > 1:
                        oneChildNum = random.choice(newChildrenNodeNums)
                        newNode.br.parts[pNum].gdasrvNum = pTree.nodes[
                            oneChildNum].br.parts[pNum].gdasrvNum
                        mp.gdasrvs[
                            newNode.br.parts[pNum].gdasrvNum].nNodes += 1

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
        # print "\ngaining a node, m %2i->%2i. logPriorRatio is %f" % (self.curTree.nInternalNodes,
        # pTree.nInternalNodes, self.logPriorRatio)

    def _getCandidateNodesForDeleteEdge(self):
        pTree = self.propTree
        nodesWithInternalEdges = [n for n in pTree.iterInternalsNoRoot()]

        # Remove any that might violate constraints.
        if self.mcmc.constraints:
            nodesToRemove = []
            for n in nodesWithInternalEdges:
                if n.br.splitKey in self.mcmc.constraints.constraints:
                    nodesToRemove.append(n)
            for n in nodesToRemove:
                nodesWithInternalEdges.remove(n)

        # We need to check that we will not be deleting model components (eg comps) that are
        # only on one node.
        if pTree.model.isHet:
            nodesToRemove = []
            for pNum in range(pTree.model.nParts):
                if pTree.model.parts[pNum].isHet:
                    mp = pTree.model.parts[pNum]
                    if mp.nComps > 1:
                        for mtNum in range(mp.nComps):
                            mt = mp.comps[mtNum]
                            # These model components are on only one node
                            if mt.nNodes <= 1:
                                for nNum in range(len(nodesWithInternalEdges)):
                                    n = nodesWithInternalEdges[nNum]
                                    if n.parts[pNum].compNum == mtNum:
                                        if n not in nodesToRemove:
                                            nodesToRemove.append(n)
                    if mp.nRMatrices > 1:
                        for mtNum in range(mp.nRMatrices):
                            mt = mp.rMatrices[mtNum]
                            # These model components are on only one node
                            if mt.nNodes <= 1:
                                for nNum in range(len(nodesWithInternalEdges)):
                                    n = nodesWithInternalEdges[nNum]
                                    if n.br.parts[pNum].rMatrixNum == mtNum:
                                        if n not in nodesToRemove:
                                            nodesToRemove.append(n)
                    if mp.nGdasrvs > 1:
                        for mtNum in range(mp.nGdasrvs):
                            mt = mp.gdasrvs[mtNum]
                            # These model components are on only one node
                            if mt.nNodes <= 1:
                                for nNum in range(len(nodesWithInternalEdges)):
                                    n = nodesWithInternalEdges[nNum]
                                    if n.br.parts[pNum].gdasrvNum == mtNum:
                                        if n not in nodesToRemove:
                                            nodesToRemove.append(n)
            # print("There are %i nodesWithInternalEdges, and I need to remove %i nodes" % (
            #     len(nodesWithInternalEdges) ,len(nodesToRemove)))
            for n in nodesToRemove:
                nodesWithInternalEdges.remove(n)
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

        if pTree.model.isHet:
            for pNum in range(pTree.model.nParts):
                if pTree.model.parts[pNum].isHet:
                    mp = pTree.model.parts[pNum]
                    if mp.nComps > 1:
                        mp.comps[theChosenNode.parts[pNum].compNum].nNodes -= 1
                    if mp.nRMatrices > 1:
                        mp.rMatrices[
                            theChosenNode.br.parts[pNum].rMatrixNum].nNodes -= 1
                    if mp.nGdasrvs > 1:
                        mp.gdasrvs[
                            theChosenNode.br.parts[pNum].gdasrvNum].nNodes -= 1

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

        if var.doMcmcSp:
            p = theNewParent.leftChild
            while p != pTree.root:
                p = p.parent
                p.flag = 1
                # print "setting flag of node %i" % p.nodeNum

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
