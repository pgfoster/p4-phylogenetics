import func
import pf
from var import var
import math
import random
import copy
import numpy
from p4exceptions import P4Error
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

        # print "Chain.init() curTree %f, propTree %f" % (self.curTree.logLike,
        # self.propTree.logLike)

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
                heatBeta = 1.0 / \
                    (1.0 + self.mcmc.tunings.chainTemp * self.tempNum)
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
        # print "proposeSp().  gen %i, About to propose %s" % (self.mcmc.gen,
        # theProposal.name)

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

        elif theProposal.name == 'rjComp':
            # print "theProposal.name = rjComp, pNum=%i" % theProposal.pNum
            self.proposeRjComp(theProposal)
            if theProposal.doAbort:
                # print "abort rjComp"
                return 0.0
            # This next line transfers the newComp.val to C
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            self.propTree.setCStuff()  # for model usage info
            # print "about to p4_setPrams() ..."
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'cmd1_compDir':
            # print "theProposal.name = cmd1_compDir, pNum=%i" %
            # theProposal.pNum
            self.proposeCmd1CompDir(theProposal)
            if 1:
                self.propTree.model.setCStuff(partNum=theProposal.pNum)
                pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
                for n in self.propTree.iterPostOrder():
                    if not n.isLeaf:
                        pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                            n.cNode, theProposal.pNum)
                pf.p4_partLogLike(self.propTree.cTree,
                                  self.propTree.data.parts[
                                      theProposal.pNum].cPart,
                                  theProposal.pNum, 0)

        elif theProposal.name == 'cmd1_allCompDir':
            # print "theProposal.name = cmd1_allCompDir, pNum=%i" %
            # theProposal.pNum
            self.proposeCmd1AllCompDir(theProposal)
            if 1:
                self.propTree.model.setCStuff(partNum=theProposal.pNum)
                pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
                for n in self.propTree.iterPostOrder():
                    if not n.isLeaf:
                        pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                            n.cNode, theProposal.pNum)
                pf.p4_partLogLike(self.propTree.cTree,
                                  self.propTree.data.parts[
                                      theProposal.pNum].cPart,
                                  theProposal.pNum, 0)

        elif theProposal.name == 'cmd1_comp0Dir':
            # print "theProposal.name = cmd1_comp0Dir, pNum=%i" %
            # theProposal.pNum
            self.proposeCmd1Comp0Dir(theProposal)

        elif theProposal.name == 'cmd1_alpha':
            # print "theProposal.name = cmd1_alpha, pNum=%i" % theProposal.pNum
            self.proposeCmd1Alpha(theProposal)

        elif theProposal.name == 'rMatrix':
            self.proposeRMatrixWithSlider(theProposal)

            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'rjRMatrix':
            # print "theProposal.name = rjRMatrix, pNum=%i" % theProposal.pNum
            self.proposeRjRMatrix(theProposal)
            if theProposal.doAbort:
                # print "abort rjRMatrix"
                return 0.0
            # This next line transfers the newRMatrix.val to C
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
            self.propTree.setCStuff()  # for model usage info
            # print "about to p4_setPrams() ..."
            pf.p4_setPrams(self.propTree.cTree, theProposal.pNum)
            for n in self.propTree.iterPostOrder():
                if not n.isLeaf:
                    pf.p4_setConditionalLikelihoodsOfInteriorNodePart(
                        n.cNode, theProposal.pNum)
            pf.p4_partLogLike(self.propTree.cTree,
                              self.propTree.data.parts[theProposal.pNum].cPart,
                              theProposal.pNum, 0)

        elif theProposal.name == 'gdasrv':
            self.proposeGdasrv(theProposal)
            self.propTree.model.setCStuff(partNum=theProposal.pNum)
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
                    print "propose brLen bad likes calc. %f %f" % (logLike1, logLike2)
                else:
                    print "propose brLen likes ok --  %f" % logLike1
                sys.exit()

        elif theProposal.name == 'eTBR':
            self.proposeETBR_Blaise(theProposal)
            if theProposal.doAbort:
                return 0.0
            else:
                if not self.propTree.preAndPostOrderAreValid:
                    self.propTree.setPreAndPostOrder()

                self.propTree.setCStuff()

                # Debugging litter ...
                if 0 and self.mcmc.gen == 270:
                    print
                    self.propTree.draw()
                    #self.propTree.node(16).br.lenChanged = 1
                    #self.propTree.node(17).br.lenChanged = 1
                    for n in self.propTree.iterNodesNoRoot():
                        if n.br.lenChanged:
                            print "    node %2i br.lenChanged" % n.nodeNum
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
                            print "    node %2i flag" % n.nodeNum
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

        else:
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
                        print "***************************** gen %i, bad curTree here b" % self.mcmc.gen
                        print x, y
                    # else:
                    # print "***************************** no difference to
                    # curTree here b"
                if 1:
                    x = sum(self.propTree.partLikes)
                    #pf.p4_partLogLike(self.propTree.cTree,self.propTree.data.parts[1].cPart, 1, 0)
                    pf.p4_treeLogLike(self.propTree.cTree, 0)
                    y = sum(self.propTree.partLikes)
                    if math.fabs(x - y) > 0.00001:
                        print "***************************** gen %i, bad propTree here b" % self.mcmc.gen
                        print x, y
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

            logLikeRatio = self.propTree.logLike - self.curTree.logLike

            # To run "without the data", which shows the effect of priors.
            logLikeRatio = 0.0

            if self.mcmc.nChains > 1:
                heatBeta = 1.0 / \
                    (1.0 + self.mcmc.tunings.chainTemp * self.tempNum)
                logLikeRatio *= heatBeta
                # print "logPriorRatio is %s, heatBeta is %s" %
                # (self.logPriorRatio, heatBeta)
                self.logPriorRatio *= heatBeta

            theSum = logLikeRatio + self.logProposalRatio + self.logPriorRatio
            if theProposal.name in ['rjComp', 'rjRMatrix']:
                theSum += self.logJacobian

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
                print "trees are the same -- ok"
                pass

        acceptMove = False

        if var.doMcmcSp:  # the speedy version
            if 0:
                # print "before proposal. curTree %f, %s   propTree %f, %s" % (
                # self.curTree.logLike, self.curTree.partLikes,
                # self.propTree.logLike, self.propTree.partLikes)
                if math.fabs(self.curTree.logLike - self.propTree.logLike) > 0.0001:
                    print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Differs before proposal"
                if math.fabs(self.curTree.logLike - sum(self.curTree.partLikes)) > 0.0001:
                    print "7777777777777777777777777777777777777777777777777777777777 bad Cur Tree"
                if math.fabs(self.propTree.logLike - sum(self.propTree.partLikes)) > 0.0001:
                    print "8888888888888888888888888888888888888888888888888888888888 bad Prop Tree"
                # print self.propTree.partLikes, type(self.propTree.partLikes)
                assert type(self.propTree.partLikes) == type(
                    self.propTree.preOrder)
                assert type(self.curTree.partLikes) == type(
                    self.curTree.preOrder)
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
                    print "Attempted %i '%s' proposals, and they all failed." % (safety, aProposal.name)
                    print "Giving up."
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

        # if aProposal.name == 'rMatrixLocation':
        #    print " acceptMove = %s" % acceptMove
        # if aProposal.name == 'cmd1_allCompDir':
        #     print " acceptMove = %s" % acceptMove

        # if aProposal.name in ['rMatrix', 'comp', 'gdasrv']:
        #    acceptMove = False

        # if self.mcmc.gen > 130 and self.mcmc.gen < 140:
        # print "-------------- (gen %5i, %20s) acceptMove = %s" %
        # (self.mcmc.gen, aProposal.name, acceptMove)

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
            if aProposal.name in ['comp', 'compDir', 'rMatrix', 'gdasrv', 'pInvar', 'cmd1_compDir', 'cmd1_allCompDir']:
                b.logLike = a.logLike
                pNum = aProposal.pNum
                b.partLikes[pNum] = a.partLikes[pNum]
                a.model.parts[pNum].copyValsTo(b.model.parts[pNum])
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
                    print "k curTree:"
                    print self.curTree.model.parts[pNum].bQETneedsReset
                    print "propTree:"
                    print self.propTree.model.parts[pNum].bQETneedsReset

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
            elif aProposal.name in ['local', 'eTBR', 'root3', 'brLen', 'polytomy']:
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
            print "diff = %f" % myDiff
            if 1 and self.mcmc.gen == 13:
                # print "Too big!"
                print "Comparing topology stuff with Tree.verifyIdentityWith() ..."
                # python level only, false for 'doSplitKeys'
                ret = self.curTree.verifyIdentityWith(self.testTree, False)
                if ret == var.DIFFERENT:
                    print "verifyIdentityOfTwoTreesInChain() tree topology stuff differs"
                else:
                    print "topology stuff seems to be the same"
                print "Python-level: Verify model prams, with Model.verifyValsWith."
                ret = self.curTree.model.verifyValsWith(
                    self.testTree.model)  # python level only
                if ret == var.DIFFERENT:
                    print "verifyIdentityOfTwoTreesInChain() model stuff differs"
                else:
                    print "model stuff appears to be the same"

                # cStuff.  This does model prams, tree and node stuff.
                print "about to pf.p4_verifyIdentityOfTwoTrees(self.curTree.cTree, self.testTree.cTree)"
                ret = pf.p4_verifyIdentityOfTwoTrees(
                    self.curTree.cTree, self.testTree.cTree)
                print "got ret %s" % ret
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
                    print "p node %2i, br.lenChanged" % n.nodeNum
                    isBad = True
                if n.flag:
                    print "p node %2i, flag" % n.nodeNum
                    isBad = True
            for n in self.curTree.iterNodesNoRoot():
                if n.br.lenChanged:
                    print "c node %2i, br.lenChanged" % n.nodeNum
                    isBad = True
                if n.flag:
                    print "c node %2i, flag" % n.nodeNum
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

                    print
                    print self.curTree.model.parts[pNum].bQETneedsReset

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

    def proposeRoot3(self, theProposal):
        """For non-biRooted trees.  Root on another internal node."""
        internalsNoRoot = [n for n in self.propTree.iterInternalsNoRoot()]
        if len(internalsNoRoot):
            newRoot = random.choice(internalsNoRoot)
            self.propTree.reRoot(
                newRoot, moveInternalName=False, fixRawSplitKeys=self.mcmc.constraints)
        else:
            print "Chain.proposeRoot3().  No other internal nodes.  Fix me."
        self.logProposalRatio = 0.0
        self.logPriorRatio = 0.0
        # if self.mcmc.constraints:
        #    print "checkSplitKeys() at the end of root3"
        #    self.propTree.checkSplitKeys()

    def proposeBrLen(self, theProposal):
        #gm = ['Chain.proposeBrLen']

        # Choose a node.
        nodesNoRoot = [n for n in self.propTree.iterNodesNoRoot()]
        theNode = random.choice(nodesNoRoot)
        #theNode = self.propTree.nodes[1]
        oldBrLen = theNode.br.len

        if 1:  # "Multiplier" proposal
            newBrLen = oldBrLen * \
                math.exp(theProposal.tuning * (random.random() - 0.5))

            # Logarithmic reflect if needed
            while (newBrLen < var.BRLEN_MIN) or (newBrLen > var.BRLEN_MAX):
                if newBrLen < var.BRLEN_MIN:
                    newBrLen = var.BRLEN_MIN * var.BRLEN_MIN / newBrLen
                elif newBrLen > var.BRLEN_MAX:
                    newBrLen = var.BRLEN_MAX * var.BRLEN_MAX / newBrLen
            theNode.br.len = newBrLen
            self.logProposalRatio = math.log(newBrLen / oldBrLen)
        else:  # Sliding window.
            newBrLen = oldBrLen + \
                (theProposal.tuning * (random.random() - 0.5))
            #newBrLen = oldBrLen + (2.0 * (random.random() - 0.5))

            # Linear reflect if needed
            while (newBrLen < var.BRLEN_MIN) or (newBrLen > var.BRLEN_MAX):
                if newBrLen < var.BRLEN_MIN:
                    newBrLen = (var.BRLEN_MIN - newBrLen) + var.BRLEN_MIN
                elif newBrLen > var.BRLEN_MAX:
                    newBrLen = var.BRLEN_MAX - (newBrLen - var.BRLEN_MAX)
            theNode.br.len = newBrLen
            self.logProposalRatio = 0.0

        if hasattr(self.mcmc.tunings, 'doInternalBrLenPrior') and self.mcmc.tunings.doInternalBrLenPrior:
            if theNode.isLeaf:
                self.logPriorRatio = self.mcmc.tunings.brLenPriorLambda * \
                    (oldBrLen - newBrLen)
            else:
                self.logPriorRatio = self.mcmc.tunings.brLenPriorLambdaForInternals * \
                    (oldBrLen - newBrLen)
        else:
            if self.mcmc.tunings.brLenPriorType == 'exponential':
                self.logPriorRatio = self.mcmc.tunings.brLenPriorLambda * \
                    (oldBrLen - newBrLen)
            else:
                self.logPriorRatio = 0.

        if var.doMcmcSp:
            theNode.br.lenChanged = True

    def proposeLocal(self, theProposal):  # from BAMBE and MrBayes.

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
            print "=" * 80
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
            print "proposeLocal() starting with this tree ..."
            pTree.draw(width=100, addToBrLen=0.2, model=True)
            if self.mcmc.constraints:
                pTree.checkSplitKeys(useOldName=True)

        if pTree.root.getNChildren() == 2:
            isBiRoot = True
            gm.append("This method is not working for biRoot'd trees yet.")
            raise P4Error(gm)
        else:
            isBiRoot = False

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
                print "=" * 50
                pTree.draw()
                gm.append(
                    "Programming error. Could not find any possibleRoots")
                raise P4Error(gm)

            newRoot = random.choice(possibleRoots)
            pTree.reRoot(
                newRoot, moveInternalName=False, fixRawSplitKeys=self.mcmc.constraints)

        if 0 and dbug:
            print "candidateC is node %i" % candidateC.nodeNum
            if oldRoot:
                print "I had to reRoot to node %i." % newRoot.nodeNum
            pTree.draw(width=80, showInternalNodeNames=1, addToBrLen=0.0)

        # We want a tree like this:
        # +------c
        # +-------|(v)
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

            print "Label the nodes a,b,c,d,u,v, and arrange into a 'standard form' ..."
            pTree.draw(
                width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)

        # At this point, the tree should look like this:
        # +------c
        # +-------|(v)
        # +------|(u)    +------d
        # |      |
        # |(a)   +-------b
        # |
        # +------X (which might be None)

        m = c.br.len + v.br.len + u.br.len
        x = u.br.len
        y = x + v.br.len
        # by default, 0.909 to 1.1
        newMRatio = math.exp(theProposal.tuning * (random.random() - 0.5))
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
            print
            print "m, the sum of brLens from a up to c, is %f" % m
            print "x, from a to u, is %f" % x
            print "y, from a to v, is %f" % y
            print "newMRatio is %f" % newMRatio
            print "newM is %f" % newM

        #################################################################
        # Detach either u or v, then re-attach somewhere between a and c.
        #################################################################

        if random.random() < 0.5:
            # detach u
            # +------c
            # +-------|(v)
            # +------|(u)    +------d
            # |      |
            # |(a)   +-------b
            # |
            # +------X
            ##
            # +----------c
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
                        print "Unable to place newX sufficiently far away from newY"
                    theProposal.doAbort = True
                    return

            if 0 and dbug:
                print "Choose to detach node u (not v)"
                print "newY is (%f * %f =) %f" % (y, newMRatio, newY)
                print "newX, a random spot along newM, is %f" % newX
                if newX < newY:
                    print "-> Since newX is still less than newY, there will be no topology change."
                else:
                    print "-> Since newX is now more than newY, there will be a topology change."

            a.leftChild = v
            v.parent = a
            v.sibling = u.sibling  # which might be None

            # now re-attach at newX
            if newX < newY:
                # no topology change, set up the same as above
                # +----------c
                # +----------|(v)
                # |(a)       +----------d
                # |
                # +----------X
                ##
                # newX    newY   newM
                # +       +      +
                ##
                # +------c
                # +-------|(v)
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
                    print "-> detach u, reattach between a and v, so no topology change"
                    pTree.draw(
                        width=80, showInternalNodeNames=1, addToBrLen=0.0)

            else:
                # a topology change
                ##
                # +----------c
                # +----------|(v)
                # |(a)       +----------d
                # |
                # +----------X
                ##
                ##
                # newY    newX   newM
                # +       +      +
                ##
                # +------c
                # +-------|(u)
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
                    print "-> detach u, re-attach between v and c, so there is a topology change"
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
            # +------c
            # +-------|(v)
            # +------|(u)    +------d
            # |      |
            # |(a)   +-------b
            # |
            # +------X
            ##
            ##
            # +----------c
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
                        print "Unable to place newY sufficiently far away from newX"
                    theProposal.doAbort = True
                    return

            if 0 and dbug:
                print "Choose to detach node v (not u)"
                print "newX is (%f * %f =) %f" % (x, newMRatio, newX)
                print "newY, a random spot along newM, is %f" % newY
                if newY < newX:
                    print "-> Since newY is now less than newX, there will be a topology change."
                else:
                    print "-> Since newY is still more than newX, there will not be a topology change."

            u.leftChild = c
            c.parent = u
            c.sibling = b

            # now reattach at newY
            if newX < newY:
                # no topology change
                ##
                # +----------c
                # +----------|(u)
                # |(a)       +----------b
                # |
                # +----------X
                ##
                # newX    newY   newM
                # +       +      +
                ##
                # +------c
                # +-------|(v)
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
                    print "-> detach v, reattach between u and c, so no topology change"
                    pTree.draw(
                        width=80, showInternalNodeNames=1, addToBrLen=0.0)
            else:
                # with a topology change
                ##
                # +----------c
                # +----------|(u)
                # |(a)       +----------b
                # |
                # +----------X
                ##
                ##
                # newY    newX   newM
                # +       +      +
                ##
                # +------c
                # +-------|(u)
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
                    print "-> detach v, re-attach between a and u, so there is a topology change"
                    print "   (splitKeys are wrong on nodes v and u, in the figure below)"
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
                print "At least 1 brLen is too short."
                theProposal.doAbort = True
                return
            elif c.br.len > var.BRLEN_MAX or v.br.len > var.BRLEN_MAX or u.br.len > var.BRLEN_MAX:
                # if dbug:
                print "At least 1 brLen is too long.  Aborting. (No big deal ...)"
                theProposal.doAbort = True
                return

        if 1:
            complain = True
            if c.br.len < var.BRLEN_MIN:
                if complain:
                    print "c  %i  too short" % self.mcmc.gen
                theProposal.doAbort = True
                return
            if v.br.len < var.BRLEN_MIN:
                if complain:
                    print "v  %i  too short" % self.mcmc.gen
                theProposal.doAbort = True
                return
            if u.br.len < var.BRLEN_MIN:
                if complain:
                    print "u  %i  too short" % self.mcmc.gen
                theProposal.doAbort = True
                return
            if c.br.len > var.BRLEN_MAX:
                if complain:
                    print "c  %i  too long" % self.mcmc.gen
                theProposal.doAbort = True
                return
            if v.br.len > var.BRLEN_MAX:
                if complain:
                    print "v  %i  too long" % self.mcmc.gen
                theProposal.doAbort = True
                return
            if u.br.len > var.BRLEN_MAX:
                if complain:
                    print "u  %i  too long" % self.mcmc.gen
                theProposal.doAbort = True
                return

        if self.mcmc.constraints and theProposal.topologyChanged:
            # Check whether any constraints have been involved, and if so abort.
            ##
            # It was like this:
            # +------c
            # +-------|(v)
            # +------|(u)    +------d
            # |      |
            # |(a)   +-------b
            # |
            # +------X
            ##
            # And now its like this:
            # +------c
            # +-------|(u)
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
        if hasattr(self.mcmc.tunings, 'doInternalBrLenPrior') and self.mcmc.tunings.doInternalBrLenPrior:
            # originally this
            # +------c
            # +-------|(v)
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
                    theSum += (self.mcmc.tunings.brLenPriorLambda * x) - \
                        (self.mcmc.tunings.brLenPriorLambda * newX)
                else:
                    theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * x) - (
                        self.mcmc.tunings.brLenPriorLambdaForInternals * newX)
                theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * (y - x)) - (
                    self.mcmc.tunings.brLenPriorLambdaForInternals * (newY - newX))
                if c.isLeaf:
                    theSum += (self.mcmc.tunings.brLenPriorLambda * (m - y)) - \
                        (self.mcmc.tunings.brLenPriorLambda * (newM - newY))
                else:
                    theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * (m - y)) - (
                        self.mcmc.tunings.brLenPriorLambdaForInternals * (newM - newY))
            else:  # with topology change
                if a.isLeaf:
                    theSum += (self.mcmc.tunings.brLenPriorLambda * x) - \
                        (self.mcmc.tunings.brLenPriorLambda * newY)
                else:
                    theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * x) - (
                        self.mcmc.tunings.brLenPriorLambdaForInternals * newY)
                theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * (y - x)) - (
                    self.mcmc.tunings.brLenPriorLambdaForInternals * (newX - newY))
                if c.isLeaf:
                    theSum += (self.mcmc.tunings.brLenPriorLambda * (m - y)) - \
                        (self.mcmc.tunings.brLenPriorLambda * (newM - newX))
                else:
                    theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * (m - y)) - (
                        self.mcmc.tunings.brLenPriorLambdaForInternals * (newM - newX))

            if 0:
                # Slow check, via priorDensities.
                theta1 = self.mcmc.tunings.brLenPriorLambda
                theta2 = self.mcmc.tunings.brLenPriorLambdaForInternals

                if newX < newY:
                    if a.isLeaf:
                        theta = self.mcmc.tunings.brLenPriorLambda
                    else:
                        theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                    prDensNu1 = theta * math.exp(-theta * x)
                    prDensNuStar1 = theta * math.exp(-theta * newX)

                    theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                    prDensNu2 = theta * math.exp(-theta * (y - x))
                    prDensNuStar2 = theta * math.exp(-theta * (newY - newX))
                    if c.isLeaf:
                        theta = self.mcmc.tunings.brLenPriorLambda
                    else:
                        theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                    prDensNu3 = theta * math.exp(-theta * (m - y))
                    prDensNuStar3 = theta * math.exp(-theta * (newM - newY))
                else:
                    if a.isLeaf:
                        theta = self.mcmc.tunings.brLenPriorLambda
                    else:
                        theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                    prDensNu1 = theta * math.exp(-theta * x)
                    prDensNuStar1 = theta * math.exp(-theta * newY)

                    theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                    prDensNu2 = theta * math.exp(-theta * (y - x))
                    prDensNuStar2 = theta * math.exp(-theta * (newX - newY))
                    if c.isLeaf:
                        theta = self.mcmc.tunings.brLenPriorLambda
                    else:
                        theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                    prDensNu3 = theta * math.exp(-theta * (m - y))
                    prDensNuStar3 = theta * math.exp(-theta * (newM - newX))
                prRat = (prDensNuStar1 * prDensNuStar2 * prDensNuStar3) / \
                    (prDensNu1 * prDensNu2 * prDensNu3)
                logPrRat = math.log(prRat)

                if math.fabs(logPrRat - theSum) > 1.e-10:
                    print "xxzz differs.  logPrRat=%g, theSum=%g" % (logPrRat, theSum)
                # else:
                #    print "s",

            self.logPriorRatio = theSum

        else:  # Do not doInternalBrLenPrior
            if self.mcmc.tunings.brLenPriorType == 'uniform':
                self.logPriorRatio = 0.0
            elif self.mcmc.tunings.brLenPriorType == 'exponential':
                self.logPriorRatio = self.mcmc.tunings.brLenPriorLambda * \
                    (m - newM)
                if 0:  # Do the same calculation the long way, edge by edge.
                    # print "logPriorRatio = %+.4f" % self.logPriorRatio,
                    foo0 = (self.mcmc.tunings.brLenPriorLambda * m) - \
                        (self.mcmc.tunings.brLenPriorLambda * newM)
                    # print "%+.4f" % foo0,

                    foo = 0.0
                    if newX < newY:  # no topology change
                        # newX    newY   newM
                        # +       +      +
                        ##
                        # +------c
                        # +-------|(v)
                        # +------|(u)    +------d
                        # |      |
                        # |(a)   +-------b
                        # |
                        # +------X
                        foo += (self.mcmc.tunings.brLenPriorLambda * x) - \
                            (self.mcmc.tunings.brLenPriorLambda * newX)
                        foo += (self.mcmc.tunings.brLenPriorLambda * (y - x)) - \
                            (self.mcmc.tunings.brLenPriorLambda *
                             (newY - newX))
                        foo += (self.mcmc.tunings.brLenPriorLambda * (m - y)) - \
                            (self.mcmc.tunings.brLenPriorLambda *
                             (newM - newY))
                    else:  # with topology change
                        # newY    newX   newM
                        # +       +      +
                        ##
                        # +------c
                        # +-------|(u)
                        # +------|(v)    +------b
                        # |      |
                        # |(a)   +-------d
                        # |
                        # +------X
                        foo += (self.mcmc.tunings.brLenPriorLambda * x) - \
                            (self.mcmc.tunings.brLenPriorLambda * newY)
                        foo += (self.mcmc.tunings.brLenPriorLambda * (y - x)) - \
                            (self.mcmc.tunings.brLenPriorLambda *
                             (newX - newY))
                        foo += (self.mcmc.tunings.brLenPriorLambda * (m - y)) - \
                            (self.mcmc.tunings.brLenPriorLambda *
                             (newM - newX))
                    # print "%+.4f" % foo
                    if (math.fabs(self.logPriorRatio - foo0) > 1.e-10):
                        print "differs-- foo0, %g %g" % (self.logPriorRatio, foo0)
                    if (math.fabs(self.logPriorRatio - foo) > 1.e-10):
                        print "differs-- foo, %g %g" % (self.logPriorRatio, foo)
            else:
                raise P4Error("This should not happen.")

        if oldRoot:
            if dbug:
                print '-------------------- about to reRoot -----------'
                pTree.draw(
                    width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)
                if self.mcmc.constraints:
                    pTree.checkSplitKeys(useOldName=True, glitch=False)
                for n in pTree.iterNodes():
                    n.name = n.oldName

            pTree.reRoot(
                oldRoot, moveInternalName=False, fixRawSplitKeys=self.mcmc.constraints)

            if dbug:
                print '--------------after reRoot --------------'
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
                print "The topology CHANGED"
            else:
                print "Topology -- no change."
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
                    print "l node %2i br.lenChanged" % n.nodeNum
        if dbug:
            if self.mcmc.constraints:
                pTree.checkSplitKeys()
        # if self.mcmc.constraints:
        #    print "checkSplitKeys() at the end of local"
        #    pTree.checkSplitKeys()

        # Check if we have a new combo of comp and rMatrix
        if theProposal.topologyChanged:
            for pNum in range(pTree.model.nParts):
                if 0 and self.mcmc.gen == 14:
                    print
                    print pTree.model.parts[pNum].bQETneedsReset
                    print "a is node %i" % a.nodeNum
                    print "u is node %i" % u.nodeNum
                    print "v is node %i" % v.nodeNum
                for n in [a, u, v]:
                    if n.br:
                        if pTree.model.parts[pNum].bQETneedsReset[n.parts[pNum].compNum][n.br.parts[pNum].rMatrixNum]:
                            # print "bQETneedsReset is set for %i, %i." % (
                            # n.parts[pNum].compNum,
                            # n.br.parts[pNum].rMatrixNum)
                            pf.p4_resetBQET(
                                pTree.model.cModel, pNum, n.parts[pNum].compNum, n.br.parts[pNum].rMatrixNum)
                if 0 and self.mcmc.gen == 14:
                    print
                    print pTree.model.parts[pNum].bQETneedsReset

        if dbug:
            for n in pTree.iterNodesNoRoot():
                if math.fabs(n.br.len - n.br.oldLen) > 0.0000001:
                    # print "Node %2i br len changed" % n.nodeNum
                    assert n.br.lenChanged
                else:
                    if n.br.lenChanged:
                        print "Node %2i lenChanged set, but its the same length." % n.nodeNum
                        raise P4Error

                # If the branch has been inverted, we will want to recalculate
                # the bigPDecks, even if the length has not really changed.
                # Trigger that intent by setting lenChanged.
                if n.br.oldNode != n:
                    # print "Node %2i branch: oldNode %2i" % (n.nodeNum,
                    # n.br.oldNode.nodeNum)
                    n.br.lenChanged = 1

    def proposeETBR_Blaise(self, theProposal):
        """Adapted and modified from Jason Evans' excellent Crux v 1.2.0

        Many thanks to JE, who wrote such clear code.  The Lakner et al
        version was modified by JE for Crux so that it works on
        polytomies, and so that it works on leaf nodes.

        Many thanks to Blaise Li who convincingly pointed out that LOCAL
        was not enough, and then pointed out that even the Crux version of
        eTBR did not appear to be reversible when there were polytomies.
        This version has a suggested modification from Blaise as described
        in his poster 'An eTBR proposal for non-binary trees in MCMC
        Bayesian phylogeny', B Li and P Foster, presented (by BL) at the
        15th Evolutionary Biology Meeting, September 27-30, 2011
        Marseilles.

        """

        # doAbort is set if brLens are too long or too short, or if a
        # constraint is violated.
        gm = ['Chain.proposeETBR_Blaise()']

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
            print "------------- eTBR gen %i -----------------" % self.mcmc.gen
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
        etbrPExt = self.mcmc.tunings.etbrPExt

        if 0 and dbug:
            print "=" * 50
            pTree.draw()
            print "starting with the tree above."

        # Choose a node, not the root.  It will have edge eA in Jason's diagram.
        # y0 will be the asterisk node in Jason's diagram.  It may be extended,
        # below.
        y0 = None
        while not y0:
            nNum = random.choice(pTree.preOrder)
            if nNum != var.NO_ORDER and nNum != pTree.root.nodeNum:
                y0 = pTree.node(nNum)
        x0 = y0.parent
        if dbug:
            print "y0 is node %i" % y0.nodeNum
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
        eA = y0.br

        # Extend x
        # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
        r1 = x0

        #  If x0 is not a leaf, it is 'unconstrained', and since x0 is the
        #  parent of y0, it will always be so in p4.
        x0Degree = pTree.getDegree(x0)
        x0Uncon = (x0Degree > 1)  # I think this will always be so.
        if x0Uncon:
            myRan = random.randrange(x0Degree - 1)
            r0 = pTree.nextNode(y0, x0)
            for i in range(myRan):
                r0 = pTree.nextNode(r0, x0)
            if r0 == x0:
                r0 = r0.parent
            x1 = r0
            if dbug:
                if x1.name:
                    x1.name += '_x1'
                else:
                    x1.name = 'x1'

            # So we are set up like this ...
            #   nR1      nR0
            #   nX0==eX==nX1

            # We name the edge here, and use it below when we modify the br.len
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
                    print "No extension from x1 was done (because x1=r0), so no rearrangement on the x side."
                pass
            else:
                # Do the rearrangement.  We need to have x0 above x1, and
                # r0 above r1.  We could guarrantee that by re-rooting on
                # r1, but that may not be needed.
                if x1.parent == x0:
                    pTree.reRoot(x1, moveInternalName=False)
                elif x0.parent == x1:
                    pass
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)
                if dbug:
                    x0.br.textDrawSymbol = 'X'

                # The r1-r0 branch might be pointing up or down, but we need
                # it to be pointing up, such that r0 is the child of r1 and so
                # r0.parent = r1.  If it is not, then we need to re-root.

                if r0.parent == r1:
                    # no need to re-root, the r branch points up
                    pass
                elif r1.parent == r0:
                    pTree.reRoot(r1, moveInternalName=False)
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)

                assert r0.parent == r1
                assert x0.parent == x1

                if dbug:
                    pTree.draw()
                    print "The drawing above is just before the rearrangement."

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

                if x0 != r1:
                    pTree.pruneSubTreeWithoutParent(x0)
                    pTree.reconnectSubTreeWithoutParent(x0, r1)

                pTree.pruneSubTreeWithoutParent(r0)
                pTree.reconnectSubTreeWithoutParent(r0, x0)

        if 1 and dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()
            if x0Uncon and r0 != x1:
                print "The drawing above shows that X extended"
            else:
                print "The drawing above shows that X did not extend."

        # Extend y
        # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
        s1 = y0

        #  If y0 is not a leaf, it is 'unconstrained'
        y0Degree = pTree.getDegree(y0)
        y0Uncon = (y0Degree > 1)  # y0 will sometimes be a leaf
        if y0Uncon:
            myRan = random.randrange(y0Degree - 1)
            s0 = pTree.nextNode(y0, y0)
            for i in range(myRan):
                s0 = pTree.nextNode(s0, y0)
            assert s0.parent == y0
            y1 = s0
            if dbug:
                if y1.name:
                    y1.name += '_y1'
                else:
                    y1.name = 'y1'

            # name the edge here, to be used below when we modify the br.len
            eY = y1.br

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

                # Since we are extending Y here, and this entire subtree goes
                # up, it should always be that s0.parent is s1.
                if s1.parent == s0:
                    gm.append("s1.parent is s0.  This should not happen.")
                    raise P4Error(gm)
                elif s0.parent == s1:
                    s0new = pTree.nextNode(s0, s0)
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)

                for i in range(myRan):
                    s0new = pTree.nextNode(s0new, s0)

                assert s0new != s0
                # if s0new == s0:
                #    s0new = s0.parent
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
                    print "No Y extension was made (because y1=s0), so nothing to do."
            else:
                # Do the rearrangement.  Since it is Y that is being
                # extended, it should always be pointing such that
                # y1.parent is y0.
                assert y1.parent == y0

                # Because we are extending Y, the s1-s0 branch should always be
                # be pointing such that s0 is the child of s1 and so
                # s0.parent = s1.
                assert s0.parent == s1

                # We want y0.parent to be y1, and we want s0.parent to be
                # s1; and only the latter is currently true.  So reRoot to
                # y1.

                pTree.reRoot(y1, moveInternalName=False)
                assert y0.parent == y1
                assert s0.parent == s1

                if dbug:
                    pTree.draw()
                    print "The drawing above is the tree just before rearrangement on the Y side."

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

                if y0 != s1:
                    pTree.pruneSubTreeWithoutParent(y0)
                    pTree.reconnectSubTreeWithoutParent(y0, s1)

                pTree.pruneSubTreeWithoutParent(s0)
                pTree.reconnectSubTreeWithoutParent(s0, y0)

        if dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()
            if y0Uncon and s0 != y1:
                print "The drawing above shows that Y extended"
            else:
                print "The drawing above shows that Y did not extend."

        if oldRoot != pTree.root:
            pTree.reRoot(oldRoot, moveInternalName=False)
        if dbug:
            pTree.draw()
            print "The above is back to the original root."

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
            print "-" * 20
            print "xRearranged = %s" % xRearranged
            print "yRearranged = %s" % yRearranged
            if eA:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eA:
                        print "eA is from node", n.nodeNum
            else:
                print "eA is None"
            if eX:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eX:
                        print "eX is from node", n.nodeNum
            else:
                print "eX is None"
            if eY:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eY:
                        print "eY is from node", n.nodeNum
            else:
                print "eY is None"
            print "x0Uncon is", x0Uncon
            print "y0Uncon is", y0Uncon

        # Are we violating constraints?
        if 1:
            if self.mcmc.constraints and theProposal.topologyChanged:
                pTree.makeSplitKeys()
                pTreeSKSet = set(
                    [n.br.splitKey for n in pTree.iterInternalsNoRoot()])
                isViolating = False
                for sk in self.mcmc.constraints.constraints:
                    if sk not in pTreeSKSet:
                        isViolating = True
                        break
                if isViolating:
                    theProposal.doAbort = True
                    return

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
        if self.mcmc.tunings.brLenPriorType == 'exponential':
            lnPrior = -self.mcmc.tunings.brLenPriorLambda * (vA1 - vA0)
            if x0Uncon:
                lnPrior += -self.mcmc.tunings.brLenPriorLambda * (vX1 - vX0)
            if y0Uncon:
                lnPrior += -self.mcmc.tunings.brLenPriorLambda * (vY1 - vY0)
        elif self.mcmc.tunings.brLenPriorType == 'uniform':
            lnPrior = 0.0

        # The proposal ratio is the product of the proposal ratios for
        # extension of each end of eA, as well as the branch multipliers.  The
        # ratio is 1 for the constrained/constrained and
        # unconstrained/unconstrained extension cases.
        #
        # nY0/nR0.
        if x0Uncon and r0 is not x1:
            if y0Uncon:
                if not r0Uncon:
                    lnProp += math.log(1.0 - etbrPExt)
            elif r0Uncon:
                lnProp += math.log(1.0 / (1.0 - etbrPExt))
        # nX0/nS0.
        if y0Uncon and s0 is not y1:
            if x0Uncon:
                if not s0Uncon:
                    lnProp += math.log(1.0 - etbrPExt)
            elif s0Uncon:
                lnProp += log(1.0 / (1.0 - etbrPExt))

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
        """Adapted from Jason Evans' excellent Crux v 1.2.0

        Many thanks to JE, who wrote such clear code.  The Crux version
        came from Lakner et al, but was modified by JE so that it works on
        polytomies, and so that it works on leaf nodes.  Also many thanks
        to Blaise Li who convincingly pointed out that LOCAL was not enough.

        It does not work with constraints yet.
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

        if 0 and self.mcmc.gen == 404:
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
        etbrPExt = self.mcmc.tunings.etbrPExt

        if 1 and dbug:
            print "=" * 50
            # pTree.draw()
            # print "starting with the tree above."

        # Choose a node, not the root.  It will have edge eA in Jason's diagram.
        # y0 will be the asterisk node in Jason's diagram.  It may be extended,
        # below.
        y0 = None
        while not y0:
            nNum = random.choice(pTree.preOrder)
            if nNum != var.NO_ORDER and nNum != pTree.root.nodeNum:
                y0 = pTree.node(nNum)
        x0 = y0.parent
        if dbug:
            print "y0 is node %i" % y0.nodeNum
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
        eA = y0.br

        # Extend x
        if dbug:
            xUp = None
        # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
        r1 = x0

        #  If x0 is not a leaf, it is 'unconstrained', and since x0 is the
        #  parent of y0, it will always be so in p4.
        x0Degree = pTree.getDegree(x0)
        x0Uncon = (x0Degree > 1)  # I think this will always be so.
        if x0Uncon:
            myRan = random.randrange(x0Degree - 1)
            r0 = pTree.nextNode(y0, x0)
            for i in range(myRan):
                r0 = pTree.nextNode(r0, x0)
            if r0 == x0:
                r0 = r0.parent
            x1 = r0
            if dbug:
                if x1.name:
                    x1.name += '_x1'
                else:
                    x1.name = 'x1'

            # So we are set up like this ...
            #   nR1      nR0
            #   nX0==eX==nX1

            # We name the edge here, and used it below when we modify the
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
            if dbug:
                pTree.draw()

            # Perform rearrangement unless it would be a no-op.  It would be a
            # no-op if r0 was still x1.
            if r0 == x1:
                # We did not extend, at all.
                # print "No extension from x1 was done, so no rearrangement on
                # the x side."
                pass
            else:
                # Do the rearrangement.  We want to collapse edge x1-x0, and so we
                # need to know if it is pointing up or down.
                if x1.parent == x0:
                    xUp = x1
                    xDn = x0
                elif x0.parent == x1:
                    xUp = x0
                    xDn = x1
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)
                assert xUp.parent == xDn
                if dbug:
                    #xUp.name += '_xUp'
                    #xDn.name += '_xDn'
                    xUp.br.textDrawSymbol = 'X'

                # The r1-r0 branch might be pointing up or down, but we need
                # it to be pointing up, such that r0 is the child of r1 and so
                # r0.parent = r1.  If it is not, then we need to re-root.

                if r0.parent == r1:
                    # no need to re-root, the r branch points up
                    pass
                elif r1.parent == r0:
                    pTree.reRoot(r1, moveInternalName=False)
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)
                assert r0.parent == r1
                assert xUp.parent == xDn

                # Collapse the node between xUp and xDn, ...
                if r1 == xUp:
                    r1 = xDn
                    if dbug:
                        # fix the internal node names
                        x1.name = x1.name[:-3]
                        x0.name += "_r1"
                theRightmostChild = xUp.rightmostChild()
                theLeftSib = xUp.leftSibling()
                if theLeftSib:
                    theLeftSib.sibling = xUp.leftChild
                else:
                    xDn.leftChild = xUp.leftChild
                for n in xUp.iterChildren():
                    n.parent = xDn
                theRightmostChild.sibling = xUp.sibling
                xUp.wipe()  # needed?
                assert r0.parent == r1

                # ... and insert that node (xUp) between r0 and r1.
                xUp.parent = r1
                xUp.leftChild = r0
                r0.parent = xUp
                if r1.leftChild == r0:
                    r1.leftChild = xUp
                else:
                    oldCh = r1.leftChild
                    while oldCh.sibling != r0:
                        oldCh = oldCh.sibling
                    oldCh.sibling = xUp
                if r0.sibling:
                    xUp.sibling = r0.sibling
                    r0.sibling = None

                if dbug:
                    pTree.setPreAndPostOrder()
                    pTree.draw()

                # this method returns y0 as well.
                pTree.pruneSubTreeWithoutParent(y0)
                pTree.reconnectSubTreeWithoutParent(y0, xUp)

                # if oldRoot != pTree.root:
                #    pTree.reRoot(oldRoot, moveInternalName=False)

        if 1 and dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()
            if x0Uncon and r0 != x1:
                print "The drawing above shows that X extended"
            else:
                print "The drawing above shows that X did not extend."

        # Extend y
        # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
        s1 = y0

        #  If y0 is not a leaf, it is 'unconstrained'
        y0Degree = pTree.getDegree(y0)
        y0Uncon = (y0Degree > 1)  # y0 will sometimes be a leaf
        if y0Uncon:
            myRan = random.randrange(y0Degree - 1)
            s0 = pTree.nextNode(y0, y0)
            for i in range(myRan):
                s0 = pTree.nextNode(s0, y0)
            assert s0.parent == y0
            y1 = s0
            if dbug:
                if y1.name:
                    y1.name += '_y1'
                else:
                    y1.name = 'y1'

            # name the edge here, to be used below when we modify the br.len
            eY = y1.br

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

                # Since we are extending Y here, and this entire subtree goes
                # up, it should always be that s0.parent is s1.
                if s1.parent == s0:
                    gm.append("s1.parent is s0.  This should not happen.")
                    raise P4Error(gm)
                elif s0.parent == s1:
                    s0new = pTree.nextNode(s0, s0)
                else:
                    gm.append("This shouldn't happen.")
                    raise P4Error(gm)

                for i in range(myRan):
                    s0new = pTree.nextNode(s0new, s0)

                assert s0new != s0
                # if s0new == s0:
                #    s0new = s0.parent
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
                pTree.draw()

            # Perform rearrangement unless it would be a no-op.  It would be a
            # no-op if s0 was still y1.
            if s0 == y1:
                # We did not extend, at all.
                pass
            else:
                # Do the rearrangement.  We want to collapse edge y1-y0.
                # Since it is Y that is being extended, it should always
                # be pointing such that y1.parent is y0.
                assert y1.parent == y0

                # Because we are extending Y, the s1-s0 branch should always be
                # be pointing such that s0 is the child of s1 and so
                # s0.parent = s1.
                assert s0.parent == s1

                # We can't easily prune off the subtree downwards, so reRoot to
                # y0
                pTree.reRoot(y0, moveInternalName=False)
                if dbug:
                    pTree.draw()

                # Collapse the node between y1 and y0, ...
                if s1 == y1:
                    s1 = y0
                    if dbug:
                        # fix the internal node names
                        y1.name = y1.name[:-3]
                        y0.name += "_s1"
                theRightmostChild = y1.rightmostChild()
                theLeftSib = y1.leftSibling()
                if theLeftSib:
                    theLeftSib.sibling = y1.leftChild
                else:
                    y0.leftChild = y1.leftChild
                for n in y1.iterChildren():
                    n.parent = y0
                theRightmostChild.sibling = y1.sibling
                y1.wipe()  # needed?
                assert s0.parent == s1

                # ... and insert that node (y1) between s0 and s1.
                y1.parent = s1
                y1.leftChild = s0
                s0.parent = y1
                if s1.leftChild == s0:
                    s1.leftChild = y1
                else:
                    oldCh = s1.leftChild
                    while oldCh.sibling != s0:
                        oldCh = oldCh.sibling
                    oldCh.sibling = y1
                if s0.sibling:
                    y1.sibling = s0.sibling
                    s0.sibling = None

                if dbug:
                    pTree.setPreAndPostOrder()
                    pTree.draw()

                # Since the tree is rooted on y0 at the moment, we prune
                # off the x-subtree.  Usually this will be x0, but it can
                # happen that due to a rearrangement above that x1 is now
                # between y0 and x0.

                theX = None
                if x0.parent == y0:
                    theX = x0
                elif x1.parent == y0:
                    theX = x1
                    # assert x0.parent == x1  Nope, not always.
                else:
                    gm.append("Fix me.")
                    raise P4Error(gm)
                # this method returns x0 as well.
                pTree.pruneSubTreeWithoutParent(theX)
                pTree.reconnectSubTreeWithoutParent(theX, y1)

                if dbug:
                    pTree.setPreAndPostOrder()
                    pTree.draw()

        if oldRoot != pTree.root:
            pTree.reRoot(oldRoot, moveInternalName=False)

        if dbug:
            pTree.setPreAndPostOrder()
            pTree.draw()
            if y0Uncon and s0 != y1:
                print "The drawing above shows that Y extended"
            else:
                print "The drawing above shows that Y did not extend."

        if dbug:
            for n in pTree.nodes:
                if n.isLeaf:
                    while '_' in n.name:
                        n.name = n.name[:-3]
                    # if n.name.endswith('_xUp'):
                    #    n.name = n.name[:-4]
                    # elif n.name.endswith('_xDn'):
                    #    n.name = n.name[:-4]
                    # else:
                    #    n.name = n.name[:-3]
                elif not n.isLeaf:
                    n.name = None
                if n.br:
                    n.br.textDrawSymbol = '-'
            # if xUp:
            #    xUp.br.textDrawSymbol = '-'
            # pTree.draw()

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
            print "-" * 20
            print "xRearranged = %s" % xRearranged
            print "yRearranged = %s" % yRearranged
            if eA:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eA:
                        print "eA is from node", n.nodeNum
            else:
                print "eA is None"
            if eX:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eX:
                        print "eX is from node", n.nodeNum
            else:
                print "eX is None"
            if eY:
                for n in pTree.iterNodesNoRoot():
                    if n.br == eY:
                        print "eY is from node", n.nodeNum
            else:
                print "eY is None"
            print "x0Uncon is", x0Uncon
            print "y0Uncon is", y0Uncon

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
            p = y0
            while 1:
                p.flag = 1
                p = p.parent
                if not p:
                    break

        # if r1 is above r0, unusually, then eX points 'up' towards r1.
        # In that case, the branch on node x0 needs to have its lenChanged
        # set.
        if xRearranged:
            if r1.parent == x0:
                if x0.parent == r0:
                    # by virtue of it being upside down.
                    x0.br.lenChanged = True

        # if s0 is above y1, then there was a rearrangement due to the
        # re-rooting, so we need to set node.br.lenChanged from y0 down to s1.

        if yRearranged:
            if y1.isAncestorOf(s0):
                if s1.isAncestorOf(y0):
                    p = y0
                    while 1:
                        p.br.lenChanged = True
                        p = p.parent
                        if p == s1:
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
        lnPrior = -self.mcmc.tunings.brLenPriorLambda * (vA1 - vA0)
        if x0Uncon:
            lnPrior += -self.mcmc.tunings.brLenPriorLambda * (vX1 - vX0)
        if y0Uncon:
            lnPrior += -self.mcmc.tunings.brLenPriorLambda * (vY1 - vY0)

        # The proposal ratio is the product of the proposal ratios for
        # extension of each end of eA, as well as the branch multipliers.  The
        # ratio is 1 for the constrained/constrained and
        # unconstrained/unconstrained extension cases.
        #
        # nY0/nR0.
        if x0Uncon and r0 is not x1:
            if y0Uncon:
                if not r0Uncon:
                    lnProp += math.log(1.0 - etbrPExt)
            elif r0Uncon:
                lnProp += math.log(1.0 / (1.0 - etbrPExt))
        # nX0/nS0.
        if y0Uncon and s0 is not y1:
            if x0Uncon:
                if not s0Uncon:
                    lnProp += math.log(1.0 - etbrPExt)
            elif s0Uncon:
                lnProp += log(1.0 / (1.0 - etbrPExt))

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
            print "j There are %i internal nodes." % self.propTree.nInternalNodes
            if self.propTree.nInternalNodes == 1:
                print "-> so its a star tree -> proposeDeleteEdge is not possible."
            elif self.propTree.nInternalNodes == self.propTree.nTax - 2:
                print "-> so its a fully-resolved tree, so proposeAddEdge is not possible."

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
            # Yes, all children.
            childrenNodeNums = pTree.getChildrenNums(theChosenPolytomy)

        nPossibleWays = math.pow(2, k - 1) - k - 1
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
            print "The nChooseKs are %s" % nChooseKs
            print "The cumSum is %s" % cumSum
            print "Since there are nPossibleWays=%i, we choose a random number from 0-%i" % (
                nPossibleWays, nPossibleWays - 1)
            print "->We chose a random number: %i" % ran
            print "So we choose the group at index %i, which means nInNewGroup=%i" % (i, nInNewGroup)
            print "So we make a new node with newChildrenNodeNums %s" % newChildrenNodeNums
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
            (1.0 / self.mcmc.tunings.brLenPriorLambda) * \
            math.log(1. - random.random())
        if newNode.br.len < var.BRLEN_MIN or newNode.br.len > var.BRLEN_MAX:
            safety = 0
            while newNode.br.len < var.BRLEN_MIN or newNode.br.len > var.BRLEN_MAX:
                newNode.br.len = - \
                    (1.0 / self.mcmc.tunings.brLenPriorLambda) * \
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
                    allOnes = 2L ** (self.propTree.nTax) - 1
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
            print "The new node is given a random branch length of %f" % newNode.br.len
            print "For the Hastings ratio ..."
            print "gamma_B is %.1f" % gamma_B
            print "n_e is %.0f" % n_e
            print "k is (still) %i, and (2^{k-1} - k - 1) = nPossibleWays is still %i" % (k, nPossibleWays)
            print "n_p = %.0f is the number of polytomies present before the move." % n_p
            print "So the hastings ratio is %f" % hastingsRatio

        self.logProposalRatio = math.log(hastingsRatio)

        if 0:
            priorRatio = self.mcmc.tunings.brLenPriorLambda * \
                math.exp(- self.mcmc.tunings.brLenPriorLambda * newNode.br.len)
            if dbug:
                print "The self.mcmc.tunings.brLenPriorLambda is %f" % self.mcmc.tunings.brLenPriorLambda
                print "So the prior ratio is %f" % priorRatio

            self.logPriorRatio = math.log(priorRatio)

            # The Jacobian
            jacobian = 1.0 / (self.mcmc.tunings.brLenPriorLambda *
                              math.exp(- self.mcmc.tunings.brLenPriorLambda * newNode.br.len))
            self.logJacobian = math.log(jacobian)
            print "logPriorRatio = %f, logJacobian = %f" % (self.logPriorRatio, self.logJacobian)

        # Here I pull a fast one, as explained in Lewis et al.  The
        # priorRatio and the Jacobian terms cancel out.  So the logs might
        # as well be zeros.
        self.logPriorRatio = 0.0
        #self.logJacobian = 0.0
        # That was easy, wasn't it?
        if self.mcmc.tunings.doPolytomyResolutionClassPrior:
            # We are gaining a node.  So the prior ratio is T_{n,m + 1} /
            # (T_{n,m} * C) .  We have the logs, and the result is the
            # log.
            if 0:
                print "-" * 30
                print 'curTree.nInternalNodes', self.curTree.nInternalNodes
                print 'pTree.nInternalNodes', pTree.nInternalNodes
                print 'logBigT[curTree.nInternalNodes]', theProposal.logBigT[self.curTree.nInternalNodes]
                # print
                # math.exp(theProposal.logBigT[self.curTree.nInternalNodes])
                print 'C ', self.mcmc.tunings.polytomyPriorLogBigC
                print 'logBigT[pTree.nInternalNodes]', theProposal.logBigT[pTree.nInternalNodes]
                # print math.exp(theProposal.logBigT[pTree.nInternalNodes])
                print "-" * 30
            self.logPriorRatio = (theProposal.logBigT[self.curTree.nInternalNodes] -
                                  (self.mcmc.tunings.polytomyPriorLogBigC +
                                   theProposal.logBigT[pTree.nInternalNodes]))

        else:
            if self.mcmc.tunings.polytomyPriorLogBigC:
                self.logPriorRatio = -self.mcmc.tunings.polytomyPriorLogBigC
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

        # We need to check that we will not be deleting modelThings that are
        # only on one node.
        if pTree.model.isHet:
            nodesToRemove = []
            for pNum in range(pTree.model.nParts):
                if pTree.model.parts[pNum].isHet:
                    mp = pTree.model.parts[pNum]
                    if mp.nComps > 1:
                        for mtNum in range(mp.nComps):
                            mt = mp.comps[mtNum]
                            # These modelThings are on only one node
                            if mt.nNodes <= 1:
                                for nNum in range(len(nodesWithInternalEdges)):
                                    n = nodesWithInternalEdges[nNum]
                                    if n.parts[pNum].compNum == mtNum:
                                        if n not in nodesToRemove:
                                            nodesToRemove.append(n)
                    if mp.nRMatrices > 1:
                        for mtNum in range(mp.nRMatrices):
                            mt = mp.rMatrices[mtNum]
                            # These modelThings are on only one node
                            if mt.nNodes <= 1:
                                for nNum in range(len(nodesWithInternalEdges)):
                                    n = nodesWithInternalEdges[nNum]
                                    if n.br.parts[pNum].rMatrixNum == mtNum:
                                        if n not in nodesToRemove:
                                            nodesToRemove.append(n)
                    if mp.nGdasrvs > 1:
                        for mtNum in range(mp.nGdasrvs):
                            mt = mp.gdasrvs[mtNum]
                            # These modelThings are on only one node
                            if mt.nNodes <= 1:
                                for nNum in range(len(nodesWithInternalEdges)):
                                    n = nodesWithInternalEdges[nNum]
                                    if n.br.parts[pNum].gdasrvNum == mtNum:
                                        if n not in nodesToRemove:
                                            nodesToRemove.append(n)
            # print "There are %i nodesWithInternalEdges, and I need to remove %i nodes" % (
            #    len(nodesWithInternalEdges) ,len(nodesToRemove))
            for n in nodesToRemove:
                nodesWithInternalEdges.remove(n)
        return nodesWithInternalEdges

    def proposeDeleteEdge(self, theProposal, candidateNodes):

        dbug = False
        pTree = self.propTree
        # print "doing proposeDeleteEdge()"
        if 0:
            print "proposeDeleteEdge(), starting with this tree ..."
            pTree.draw()
            print "m There are %i internal nodes (before deleting the edge)." % pTree.nInternalNodes

        if not candidateNodes:
            raise P4Error(
                "proposeDeleteEdge() could not find a good node to attempt to delete.")

        theChosenNode = random.choice(candidateNodes)
        if dbug:
            print "There are %i candidateNodes." % len(candidateNodes)
            print "node nums %s" % [n.nodeNum for n in candidateNodes]
            print "Randomly choose node %s" % theChosenNode.nodeNum

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
            priorRatio = 1.0 / (self.mcmc.tunings.brLenPriorLambda *
                                math.exp(- self.mcmc.tunings.brLenPriorLambda * theChosenNode.br.len))
            if dbug:
                print "The self.mcmc.tunings.brLenPriorLambda is %f" % self.mcmc.tunings.brLenPriorLambda
                print "So the prior ratio is %f" % priorRatio

            self.logPriorRatio = math.log(priorRatio)

            # The Jacobian
            jacobian = self.mcmc.tunings.brLenPriorLambda * \
                math.exp(- self.mcmc.tunings.brLenPriorLambda *
                         theChosenNode.br.len)
            self.logJacobian = math.log(jacobian)
            print "logPriorRatio = %f, logJacobian = %f" % (self.logPriorRatio, self.logJacobian)

        # Here I pull a fast one, as explained in Lewis et al.  The
        # priorRatio and the Jacobian terms cancel out.  So the logs might
        # as well be zeros.
        self.logPriorRatio = 0.0
        #self.logJacobian = 0.0
        # That was easy, wasn't it?

        if self.mcmc.tunings.doPolytomyResolutionClassPrior:
            # We are losing a node.  So the prior ratio is (T_{n,m} * C) /
            # T_{n,m - 1}.  We have the logs, and the result is the log.
            if 0:
                print "-" * 30
                print 'curTree.nInternalNodes', self.curTree.nInternalNodes
                print 'pTree.nInternalNodes', pTree.nInternalNodes
                print 'logBigT[curTree.nInternalNodes]', theProposal.logBigT[self.curTree.nInternalNodes]
                # print
                # math.exp(theProposal.logBigT[self.curTree.nInternalNodes])
                print 'C ', self.mcmc.tunings.polytomyPriorLogBigC
                print 'logBigT[pTree.nInternalNodes]', theProposal.logBigT[pTree.nInternalNodes]
                # print math.exp(theProposal.logBigT[pTree.nInternalNodes])
                print "-" * 30
            self.logPriorRatio = ((theProposal.logBigT[self.curTree.nInternalNodes] +
                                   self.mcmc.tunings.polytomyPriorLogBigC) -
                                  theProposal.logBigT[pTree.nInternalNodes])

        else:
            if self.mcmc.tunings.polytomyPriorLogBigC:
                self.logPriorRatio = self.mcmc.tunings.polytomyPriorLogBigC
            else:
                self.logPriorRatio = 0.0

        # print " losing a node, m %2i->%2i. logPriorRatio is %f" % (self.curTree.nInternalNodes,
        # pTree.nInternalNodes, self.logPriorRatio)

    def proposeCompWithSlider(self, theProposal):
        gm = ['Chain.proposeCompWithSlider()']

        mt = self.propTree.model.parts[
            theProposal.pNum].comps[theProposal.mtNum]
        dim = self.propTree.model.parts[theProposal.pNum].dim

        # mt.val is a list, not a numpy array
        #assert type(mt.val) == numpy.ndarray

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

        mt = self.propTree.model.parts[
            theProposal.pNum].comps[theProposal.mtNum]
        dim = self.propTree.model.parts[theProposal.pNum].dim

        # mt.val is a list of floats, not a numpy.ndarray
        # print type(mt.val), type(mt.val[0])

        # The tuning is the Dirichlet alpha.
        # print theProposal.tuning

        # This method uses func.dirichlet1, which is for lists not numpy
        # arrays.  A copy of inSeq is made, and the copy is modified and
        # returned.
        #dirichlet1(inSeq, alpha, theMin, theMax)
        newVal = func.dirichlet1(
            mt.val, theProposal.tuning, var.PIVEC_MIN, 1 - var.PIVEC_MIN)

        self.logProposalRatio = 0.0

        rangeDim = range(dim)
        mySum = 0.0
        for stNum in rangeDim:
            mySum += newVal[stNum] * theProposal.tuning
        x = pf.gsl_sf_lngamma(mySum)
        for stNum in rangeDim:
            x -= pf.gsl_sf_lngamma(newVal[stNum] * theProposal.tuning)
        for stNum in rangeDim:
            x += ((newVal[stNum] * theProposal.tuning) - 1.) * \
                math.log(mt.val[stNum])

        mySum = 0.0
        for stNum in rangeDim:
            mySum += mt.val[stNum] * theProposal.tuning
        y = pf.gsl_sf_lngamma(mySum)
        for stNum in rangeDim:
            y -= pf.gsl_sf_lngamma(mt.val[stNum] * theProposal.tuning)
        for stNum in rangeDim:
            y += ((mt.val[stNum] * theProposal.tuning) - 1.) * \
                math.log(newVal[stNum])
        self.logProposalRatio = x - y
        mt.val = newVal

        # The prior here is a flat Dirichlet, ie Dirichlet(1, 1, 1, ...,
        # 1).  If it is informative, then the prior is affected.
        self.logPriorRatio = 0.0

    def proposeRjComp(self, theProposal):
        gm = ['Chain.proposeRjComp()']
        theProposal.doAbort = False
        mp = self.propTree.model.parts[theProposal.pNum]
        assert mp.rjComp
        assert mp.rjComp_k >= 1
        assert mp.rjComp_k <= mp.nComps

        # If the pool size k, mp.rjComp_k, is only 1, then we can only split
        if mp.rjComp_k == 1:
            self.proposeSplitComp(theProposal)

        # If k, mp.rjComp_k is mp.nComps, then we can only merge
        if mp.rjComp_k == mp.nComps:
            self.proposeMergeComp(theProposal)

        # Otherwise, choose randomly.
        if random.random() < 0.5:
            self.proposeSplitComp(theProposal)
        else:
            self.proposeMergeComp(theProposal)

    def proposeSplitComp(self, theProposal):
        gm = ['Chain.proposeSplitComp()']

        # var.rjCompUniformAllocationPrior  True by default
        # theProposal.tuning 200.  becomes p0 below

        mp = self.propTree.model.parts[theProposal.pNum]
        assert mp.rjComp
        dim = mp.dim

        # Check that k is less than k_max, which is nComps.  This should have been checked before, but check again.
        # print gm[0], "rjComp_k is currently %i, with %i comps" %
        # (mp.rjComp_k, mp.nComps)
        assert mp.rjComp_k < mp.nComps

        # Select an existing comp vector from the pool
        pool = [c for c in mp.comps if c.rj_isInPool]
        notInPool = [c for c in mp.comps if not c.rj_isInPool]
        assert notInPool  # or else we can't split
        pi0 = random.choice(pool)

        # The nodes currently associated with pi0
        beta0 = [n for n in self.propTree.iterNodes(
        ) if n.parts[theProposal.pNum].compNum == pi0.num]
        b0 = float(len(beta0))
        # print gm[0], "comp %i is chosen, f=%f, currently on nodes" %
        # (pi0.num, pi0.rj_f), [n.nodeNum for n in beta0]

        # Divvy up the contents of beta0 into (new) beta1 and beta2, based on
        # probability u
        if var.rjCompUniformAllocationPrior:
            u = 0.5
        else:
            u = random.random()

        beta1 = []
        beta2 = []
        for it in beta0:
            r = random.random()
            if r < u:
                beta1.append(it)
            else:
                beta2.append(it)
        b1 = float(len(beta1))
        b2 = float(len(beta2))
        bPrime0 = b0 + 2.
        bPrime1 = b1 + 1.
        bPrime2 = b2 + 1.

        # Calculation of f1 and f2 depends on u
        f0 = pi0.rj_f
        f1 = u * f0
        f2 = (1.0 - u) * f0

        uu = [random.normalvariate(0., 1.) for i in range(dim)]
        p0 = theProposal.tuning
        s0 = random.gammavariate(p0, 1.)
        m0 = [s0 * it for it in pi0.val]
        # print m0

        # I get a math range error here -- needs debugging.
        # m1 = [m0[j] * math.exp((bPrime0 * uu[j])/(bPrime1 * math.sqrt(m0[j])))
        #      for j in range(dim)]
        # m2 = [m0[j] * math.exp((-bPrime0 * uu[j])/(bPrime2 * math.sqrt(m0[j])))
        #      for j in range(dim)]

        safety = 0
        while 1:
            try:
                m1 = [m0[j] * math.exp((bPrime0 * uu[j]) / (bPrime1 * math.sqrt(m0[j])))
                      for j in range(dim)]
                m2 = [m0[j] * math.exp((-bPrime0 * uu[j]) / (bPrime2 * math.sqrt(m0[j])))
                      for j in range(dim)]
                break
            except OverflowError:
                print "Overflow error in splitComp() (%2i)" % safety
                safety += 1
                if safety >= 100:
                    theProposal.doAbort = True
                    # print "Too many overflows in splitComp.  Aborting!"
                    return
                uu = [random.normalvariate(0., 1.) for i in range(dim)]

        if 0:
            # testing ...
            m1 = [m0[j] * math.exp((bPrime0 * uu[j]) / bPrime1)
                  for j in range(dim)]
            m2 = [m0[j] * math.exp((-bPrime0 * uu[j]) / bPrime2)
                  for j in range(dim)]

        # Long form of the above for debugging --
        if 0:
            m1 = []
            for j in range(dim):
                top = (bPrime0 * uu[j])
                bottom = (bPrime1 * math.sqrt(m0[j]))
                quot = top / bottom
                try:
                    myexp = math.exp(quot)
                except OverflowError:
                    gm.append(
                        "Got overflow error for m1 exp(%f) at j=%i" % (quot, j))
                    gm.append("s0 is %f" % s0)
                    gm.append("bPrime0 = %f" % bPrime0)
                    gm.append("uu[j] = %f" % uu[j])
                    gm.append("bPrime1 = %f" % bPrime1)
                    gm.append("m0[j] = %f, sqrt=%f" %
                              (m0[j], math.sqrt(m0[j])))
                    gm.append("m0 is %s" % m0)
                    gm.append("top = %f" % top)
                    gm.append("bottom = %f" % bottom)
                    raise P4Error(gm)
                m1.append(m0[j] * myexp)

            m2 = []
            for j in range(dim):
                top = (-bPrime0 * uu[j])
                bottom = (bPrime2 * math.sqrt(m0[j]))
                quot = top / bottom
                try:
                    myexp = math.exp(quot)
                except OverflowError:
                    gm.append(
                        "Got overflow error for m2 exp(%f) at j=%i" % (quot, j))
                    gm.append("s0 is %f" % s0)
                    gm.append("-bPrime0 = %f" % -bPrime0)
                    gm.append("uu[j] = %f" % uu[j])
                    gm.append("bPrime2 = %f" % bPrime2)
                    gm.append("m0[j] = %f, sqrt=%f" %
                              (m0[j], math.sqrt(m0[j])))
                    gm.append("m0 is %s" % m0)
                    gm.append("top = %f" % top)
                    gm.append("bottom = %f" % bottom)
                    raise P4Error(gm)
                m2.append(m0[j] * myexp)

        if 0:
            # Loggified version, as in Gowri-Shankar and Rattray, eqn 7.
            log_m1 = [math.log(m0[j]) + ((bPrime0 * uu[j]) / (bPrime1 * math.sqrt(m0[j])))
                      for j in range(dim)]
            log_m2 = [math.log(m0[j]) - ((bPrime0 * uu[j]) / (bPrime2 * math.sqrt(m0[j])))
                      for j in range(dim)]

            try:
                m1 = [math.exp(it) for it in log_m1]
            except OverflowError:
                gm.append("m0 = %s" % m0)
                gm.append("log_m1 = %s" % log_m1)
                gm.append('overflow m1')
                raise P4Error(gm)
            try:
                m2 = [math.exp(it) for it in log_m2]
            except OverflowError:
                gm.append("m0 = %s" % m0)
                gm.append("log_m2 = %s" % log_m2)
                gm.append("overflow m2")
                raise P4Error(gm)

        if 0:
            print m0
            print m1
            print m2
            print

        s1 = sum(m1)
        s2 = sum(m2)
        newVal1 = [it / s1 for it in m1]
        newVal2 = [it / s2 for it in m2]

        # print newVal1
        # print newVal2

        if 1:
            # Peter adds, the following few lines to make sure the vals are
            # more than var.PIVEC_MIN
            isChanged = False
            for vNum in range(len(newVal1)):
                isGood = False
                while not isGood:
                    # print "gen %i" % self.mcmc.gen
                    if newVal1[vNum] < var.PIVEC_MIN:
                        newVal1[vNum] = (
                            var.PIVEC_MIN - newVal1[vNum]) + var.PIVEC_MIN
                        isChanged = True
                    else:
                        isGood = True
            if isChanged:
                s1 = sum(newVal1)
                newVal1 = [it / s1 for it in newVal1]

            isChanged = False
            for vNum in range(len(newVal2)):
                isGood = False
                while not isGood:
                    # print "y gen %i" % self.mcmc.gen
                    if newVal2[vNum] < var.PIVEC_MIN:
                        newVal2[vNum] = (
                            var.PIVEC_MIN - newVal2[vNum]) + var.PIVEC_MIN
                        isChanged = True
                    else:
                        isGood = True
            if isChanged:
                s2 = sum(newVal2)
                newVal2 = [it / s2 for it in newVal2]

        # print newVal1
        # print newVal2

        # Log prior ratio
        # We could have a prior on the pool size, reflected in t1.  If all pool
        # sizes are equally probable, then t1 = 0
        t1 = 0.

        if var.rjCompUniformAllocationPrior:
            b = len([n for n in self.propTree.iterNodes()])
            t2 = b * (math.log(mp.rjComp_k) - math.log(mp.rjComp_k + 1))
        else:
            t2 = (b1 * math.log(f1)) + \
                (b2 * math.log(f2)) - (b0 * math.log(f0))

        # t3 is for the prior on comp vectors.  With the Dirichlet prior alpha
        # values all 1, t3 is log Gamma dim
        t3 = pf.gsl_sf_lngamma(dim)

        # t4 is for the f values.
        if var.rjCompUniformAllocationPrior:
            t4 = 0.0
        else:
            # If its a uniform Dirichlet, then t4 = log k, where k is from
            # before the split
            t4 = math.log(mp.rjComp_k)

        self.logPriorRatio = t1 + t2 + t3 + t4

        # Log proposal ratio
        if mp.rjComp_k == 1:
            t1 = math.log(0.5)
        else:
            t1 = 0.
        if var.rjCompUniformAllocationPrior:
            t2 = b0 * math.log(2.)
        else:
            t2 = (b0 * math.log(f0)) - (b1 * math.log(f1)) - \
                (b2 * math.log(f2))          # this was changed 26 sept

        # for t3, below, do some pre-calculations
        sum_uu2 = sum([u * u for u in uu])
        lastTerm = -pf.gsl_sf_lngamma(p0) + (0.5 * sum_uu2) + \
            ((dim / 2.) * math.log(2 * math.pi))
        t3 = (s0 - s1 - s2) + ((p0 - 1.) *
                               (math.log(s1) + math.log(s2) - math.log(s0))) + lastTerm
        self.logProposalRatio = t1 + t2 + t3
        # print t1,t2,t3,s0,s1,s2

        #self.logProposalRatio = 0.

        # The Jacobian
        lastTerm = 0.5 * sum([math.log(v) for v in pi0.val])    # added 26 sept
        t1 = ((((3. * dim) - 2.) / 2.) * math.log(s0)) - \
            ((dim - 1.) * (math.log(s1) + math.log(s2))) + lastTerm
        t2 = (2. * dim * math.log(bPrime0)) - \
            (dim * (math.log(bPrime1) + math.log(bPrime2)))
        t3 = sum([uu[j] / (math.sqrt(s0 * pi0.val[j])) for j in range(dim)])
        t3 = ((bPrime0 * (bPrime2 - bPrime1)) / (bPrime1 * bPrime2)) * t3

        if var.rjCompUniformAllocationPrior:
            self.logJacobian = t1 + t2 + t3
        else:
            self.logJacobian = t1 + t2 + t3 + math.log(f0)

        # We will now make pi1 and pi2.  The pi1 will be made from pi0,
        # and pi2 will be popped from the notInPool list.  We have newVal1
        # and newVal2 which will be their vals, and we assign them to
        # nodes in beta1 and beta2, and give them rj_f values of f1 and
        # f2.
        pi1 = pi0
        pi1.val = newVal1
        pi1.rj_f = f1
        for n in beta1:
            # not needed, its already that.
            n.parts[theProposal.pNum].compNum = pi1.num
            pf.p4_setCompNum(n.cNode, theProposal.pNum, pi1.num)
        pi1.nNodes = b1

        pi2 = notInPool.pop()
        pi2.val = newVal2
        pi2.rj_f = f2
        for n in beta2:
            n.parts[theProposal.pNum].compNum = pi2.num
            pf.p4_setCompNum(n.cNode, theProposal.pNum, pi2.num)
        pi2.nNodes = b2
        pi2.rj_isInPool = True

        self.propTree.model.parts[theProposal.pNum].rjComp_k += 1
        # print "...finished proposeSplitComp()"

    def proposeMergeComp(self, theProposal):
        gm = ['Chain.proposeMergeComp()']

        mp = self.propTree.model.parts[theProposal.pNum]
        assert mp.rjComp
        dim = mp.dim
        p0 = theProposal.tuning

        # Check that k is more than 1.  This should have been checked before, but check again.
        # print "rjComp_k is currently %i, with %i comps" % (mp.rjComp_k,
        # mp.nComps)
        if mp.rjComp_k <= 1:
            gm.append("part %i, rjComp_k = %i" %
                      (theProposal.pNum, mp.rjComp_k))
            pool = [c for c in mp.comps if c.rj_isInPool]
            gm.append(
                'len of pool = %i (should be the same as rjComp_k)' % len(pool))
            gm.append(
                "rjComp_k, the pool size, should be more than 1 for a merge.  This isn't.")
            raise P4Error(gm)

        # Choose two comps (to make into one).  They must be in the pool.
        pool = [c for c in mp.comps if c.rj_isInPool]
        assert len(pool) == mp.rjComp_k
        pi1, pi2 = random.sample(pool, 2)
        # print "proposing to merge comps %i and %i" % (pi1.num, pi2.num)

        beta1 = []
        beta2 = []
        for n in self.propTree.iterNodes():
            theCompNum = n.parts[theProposal.pNum].compNum
            if theCompNum == pi1.num:
                beta1.append(n)
            elif theCompNum == pi2.num:
                beta2.append(n)
        beta0 = beta1 + beta2
        b1 = float(len(beta1))
        b2 = float(len(beta2))
        b0 = float(b1 + b2)
        assert len(beta0) == b0
        bPrime0 = b0 + 2.
        bPrime1 = b1 + 1.
        bPrime2 = b2 + 1.
        f1 = pi1.rj_f
        f2 = pi2.rj_f
        f0 = f1 + f2

        # Obtain composition vector proposal
        s1 = random.gammavariate(p0, 1.)
        s2 = random.gammavariate(p0, 1.)
        m1 = [v * s1 for v in pi1.val]
        m2 = [v * s2 for v in pi2.val]
        # print "m1 = ", m1
        # print "m2 = ", m2
        # print b0, b1, b2, bPrime0, bPrime1, bPrime2
        m0 = [math.exp(((bPrime1 / bPrime0) * math.log(m1[j])) + ((bPrime2 / bPrime0) * math.log(m2[j])))
              for j in range(dim)]

        # print "m0 = ", m0
        s0 = sum(m0)

        newVal0 = [m0k / s0 for m0k in m0]

        # Log prior ratio
        # We could have a prior on the pool size, reflected in t1.  If all pool
        # sizes are equally probable, then t1 = 0
        t1 = 0.

        if var.rjCompUniformAllocationPrior:
            b = len([n for n in self.propTree.iterNodes()])
            t2 = b * (math.log(mp.rjComp_k) - math.log(mp.rjComp_k - 1))
        else:
            t2 = (b0 * math.log(f0)) - \
                (b1 * math.log(f1)) - (b2 * math.log(f2))

        # t3 is for the prior on comp vectors.  With the Dirichlet prior alpha
        # values all 1, t3 is log Gamma dim
        t3 = -pf.gsl_sf_lngamma(dim)

        # t4 is for the f values.
        if var.rjCompUniformAllocationPrior:
            t4 = 0.0
        else:
            # If its a uniform Dirichlet, then t4 = - log (k - 1), where k is
            # from before the merge
            t4 = -math.log(mp.rjComp_k - 1)

        self.logPriorRatio = t1 + t2 + t3 + t4

        # Log proposal ratio
        if mp.rjComp_k == mp.nComps:  # nComps is k_max
            t1 = math.log(0.5)
        else:
            t1 = 0.

        if var.rjCompUniformAllocationPrior:
            t2 = - (b0 * math.log(2.))
        else:
            t2 = (b1 * math.log(f1)) + \
                (b2 * math.log(f2)) - (b0 * math.log(f0))

        # for t3, below, do some pre-calculations
        uu = [(bPrime1 / bPrime0) * math.sqrt(m0[j])
              * (math.log(m1[j]) - math.log(m0[j])) for j in range(dim)]
        sum_uu2 = sum([u * u for u in uu])
        lastTerm = pf.gsl_sf_lngamma(p0) - (0.5 * sum_uu2) - \
            ((dim / 2.) * math.log(2 * math.pi))
        t3 = (s1 + s2 - s0) - ((p0 - 1.) *
                               (math.log(s1) + math.log(s2) - math.log(s0))) + lastTerm

        #logSterm = ((1. - p0) * (math.log(s1) + math.log(s2) - math.log(s0)))
        # print "s1=%.1f s2=%.1f s0=%.1f    sum_uu2=%.1f  logGamma(p0)=%.1f, logSterm=%.1f" % (
        #    s1, s2, s0, sum_uu2, pf.gsl_sf_lngamma(p0), logSterm)
        self.logProposalRatio = t1 + t2 + t3

        #self.logProposalRatio = 20.

        # The Jacobian
        lastTerm = 0.5 * sum([math.log(v) for v in newVal0])  # new 26 sept
        t1 = ((((3. * dim) - 2.) / 2.) * math.log(s0)) - \
            ((dim - 1.) * (math.log(s1) + math.log(s2))) + lastTerm
        t2 = (2. * dim * math.log(bPrime0)) - \
            (dim * (math.log(bPrime1) + math.log(bPrime2)))
        t3 = sum([uu[j] / (math.sqrt(s0 * newVal0[j])) for j in range(dim)])
        t3 = ((bPrime0 * (bPrime2 - bPrime1)) / (bPrime1 * bPrime2)) * t3
        if var.rjCompUniformAllocationPrior:
            self.logJacobian = -(t1 + t2 + t3)
        else:
            self.logJacobian = -(t1 + t2 + t3 + math.log(f0))

        # Merge pi1 and pi2 => pi0, where pi0 is actually pi1, re-used, by
        # giving "0" values to pi1 = pi0
        pi1.rj_f = f0
        pi1.val = newVal0
        for n in beta0:
            n.parts[theProposal.pNum].compNum = pi1.num
            pf.p4_setCompNum(n.cNode, theProposal.pNum, pi1.num)
        pi1.nNodes = b0
        mp.rjComp_k -= 1
        pi2.rj_isInPool = False
        pi2.nNodes = 0

    def proposeRMatrixWithSlider(self, theProposal):

        # print "rMatrix proposal. the tuning is %s" % theProposal.tuning

        assert var.rMatrixNormalizeTo1
        mtCur = self.curTree.model.parts[
            theProposal.pNum].rMatrices[theProposal.mtNum]
        mtProp = self.propTree.model.parts[
            theProposal.pNum].rMatrices[theProposal.mtNum]
        if mtProp.spec == '2p':
            # For 2p, its actually a Dirichlet, not a slider.  All this is
            # stolen from MrBayes, where the default tuning is 50.  In
            # MrBayes, the "alphaDir" is a 2-item list of Dirichlet
            # parameters (not the multiplier) but they are both by default
            # 1, which makes the prior ratio 1.0 and the logPriorRatio
            # zero.

            old = [0.0, 0.0]
            old[0] = mtCur.val / (mtCur.val + 1.0)
            old[1] = 1.0 - old[0]
            new = func.dirichlet1(
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
            mt = self.propTree.model.parts[
                theProposal.pNum].rMatrices[theProposal.mtNum]

            # mt.val is a numpy array
            assert type(mt.val) == numpy.ndarray

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

    def proposeRjRMatrix(self, theProposal):
        gm = ['Chain.proposeRjRMatrix()']
        mp = self.propTree.model.parts[theProposal.pNum]
        assert mp.rjRMatrix
        assert mp.rjRMatrix_k >= 1
        assert mp.rjRMatrix_k <= mp.nRMatrices

        # If the pool size k, mp.rjRMatrix_k, is only 1, then we can only split
        if mp.rjRMatrix_k == 1:
            self.proposeSplitRMatrix(theProposal)

        # If k, mp.rjRMatrix_k is mp.nRMatrices, then we can only merge
        if mp.rjRMatrix_k == mp.nRMatrices:
            self.proposeMergeRMatrix(theProposal)

        # Otherwise, choose randomly.
        if random.random() < 0.5:
            self.proposeSplitRMatrix(theProposal)
        else:
            self.proposeMergeRMatrix(theProposal)

    def proposeSplitRMatrix(self, theProposal):
        gm = ['Chain.proposeSplitRMatrix()']
        # var.rjRMatrixUniformAllocationPrior  True by default
        # theProposal.tuning 300.  becomes p0 below

        mp = self.propTree.model.parts[theProposal.pNum]
        assert mp.rjRMatrix
        rDim = ((mp.dim * mp.dim) - mp.dim) / 2

        # Check that k is less than k_max, which is nRMatrices.  This should have been checked before, but check again.
        # print gm[0], "rjRMatrix_k is currently %i, with %i rMatrices" %
        # (mp.rjRMatrix_k, mp.nRMatrices)
        assert mp.rjRMatrix_k < mp.nRMatrices

        # Select an existing rMatrix from the pool
        pool = [c for c in mp.rMatrices if c.rj_isInPool]
        assert mp.rjRMatrix_k == len(pool)
        notInPool = [c for c in mp.rMatrices if not c.rj_isInPool]
        assert notInPool  # or else we can't split
        assert mp.nRMatrices == len(pool) + len(notInPool)
        rm0 = random.choice(pool)

        # The nodes currently associated with rm0
        beta0 = [n for n in self.propTree.iterNodesNoRoot(
        ) if n.br.parts[theProposal.pNum].rMatrixNum == rm0.num]
        b0 = float(len(beta0))
        # print gm[0], "rMatrix %i is chosen, f=%f, currently on nodes" %
        # (rm0.num, rm0.rj_f), [n.nodeNum for n in beta0]

        # Divvy up the contents of beta0 into (new) beta1 and beta2, based on
        # probability u
        if var.rjRMatrixUniformAllocationPrior:
            u = 0.5
        else:
            u = random.random()

        beta1 = []
        beta2 = []
        for it in beta0:
            r = random.random()
            if r < u:
                beta1.append(it)
            else:
                beta2.append(it)
        b1 = float(len(beta1))
        b2 = float(len(beta2))
        bPrime0 = b0 + 2.
        bPrime1 = b1 + 1.
        bPrime2 = b2 + 1.

        # Calculation of f1 and f2 depends on u
        f0 = rm0.rj_f
        f1 = u * f0
        f2 = (1.0 - u) * f0

        uu = [random.normalvariate(0., 1.) for i in range(rDim)]
        p0 = theProposal.tuning
        s0 = random.gammavariate(p0, 1.)
        m0 = [s0 * it for it in rm0.val]
        # print m0

        # I get a math range error here -- needs debugging.
        # m1 = [m0[j] * math.exp((bPrime0 * uu[j])/(bPrime1 * math.sqrt(m0[j])))
        #      for j in range(rDim)]
        # m2 = [m0[j] * math.exp((-bPrime0 * uu[j])/(bPrime2 * math.sqrt(m0[j])))
        #      for j in range(rDim)]
        safety = 0
        while 1:
            try:
                m1 = [m0[j] * math.exp((bPrime0 * uu[j]) / (bPrime1 * math.sqrt(m0[j])))
                      for j in range(rDim)]
                m2 = [m0[j] * math.exp((-bPrime0 * uu[j]) / (bPrime2 * math.sqrt(m0[j])))
                      for j in range(rDim)]
                break
            except OverflowError:
                print "Overflow error in splitRMatrix() (%2i)" % safety
                safety += 1
                if safety >= 100:
                    theProposal.doAbort = True
                    print "Too many overflows in splitComp.  Aborting!"
                    return
                uu = [random.normalvariate(0., 1.) for i in range(rDim)]

        if 0:
            # Long form of the above for debugging --
            m1 = []
            for j in range(rDim):
                top = (bPrime0 * uu[j])
                bottom = (bPrime1 * math.sqrt(m0[j]))
                quot = top / bottom
                try:
                    myexp = math.exp(quot)
                except OverflowError:
                    gm.append(
                        "Got overflow error for m1 exp(%f) at j=%i" % (quot, j))
                    gm.append("bPrime0 = %f" % bPrime0)
                    gm.append("uu[j] = %f" % uu[j])
                    gm.append("bPrime1 = %f" % bPrime1)
                    gm.append("m0[j] = %f, sqrt=%f" %
                              (m0[j], math.sqrt(m0[j])))
                    gm.append("m0 is %s" % m0)
                    gm.append("top = %f" % top)
                    gm.append("bottom = %f" % bottom)
                    raise P4Error(gm)
                m1.append(m0[j] * myexp)

            m2 = []
            for j in range(rDim):
                top = (-bPrime0 * uu[j])
                bottom = (bPrime2 * math.sqrt(m0[j]))
                quot = top / bottom
                try:
                    myexp = math.exp(quot)
                except OverflowError:
                    gm.append(
                        "Got overflow error for m2 exp(%f) at j=%i" % (quot, j))
                    gm.append("-bPrime0 = %f" % -bPrime0)
                    gm.append("uu[j] = %f" % uu[j])
                    gm.append("bPrime2 = %f" % bPrime2)
                    gm.append("m0[j] = %f, sqrt=%f" %
                              (m0[j], math.sqrt(m0[j])))
                    gm.append("m0 is %s" % m0)
                    gm.append("top = %f" % top)
                    gm.append("bottom = %f" % bottom)
                    raise P4Error(gm)
                m2.append(m0[j] * myexp)

        s1 = sum(m1)
        s2 = sum(m2)
        newVal1 = [it / s1 for it in m1]
        newVal2 = [it / s2 for it in m2]

        # print newVal1
        # print newVal2

        if 1:
            # Peter adds, the following few lines to get the vals more than
            # var.RATE_MIN
            isChanged = False
            for vNum in range(len(newVal1)):
                isGood = False
                while not isGood:
                    # print "gen %i" % self.mcmc.gen
                    if newVal1[vNum] < var.RATE_MIN:
                        newVal1[vNum] = (
                            var.RATE_MIN - newVal1[vNum]) + var.RATE_MIN
                        isChanged = True
                    else:
                        isGood = True
            if isChanged:
                s1 = sum(newVal1)
                newVal1 = [it / s1 for it in newVal1]

            isChanged = False
            for vNum in range(len(newVal2)):
                isGood = False
                while not isGood:
                    # print "y gen %i" % self.mcmc.gen
                    if newVal2[vNum] < var.RATE_MIN:
                        newVal2[vNum] = (
                            var.RATE_MIN - newVal2[vNum]) + var.RATE_MIN
                        isChanged = True
                    else:
                        isGood = True
            if isChanged:
                s2 = sum(newVal2)
                newVal2 = [it / s2 for it in newVal2]

        # print newVal1
        # print newVal2

        # Log prior ratio
        # We could have a prior on the pool size, reflected in t1.  If all pool
        # sizes are equally probable, then t1 = 0
        t1 = 0.

        if var.rjRMatrixUniformAllocationPrior:
            b = len([n for n in self.propTree.iterNodesNoRoot()])
            t2 = b * (math.log(mp.rjRMatrix_k) - math.log(mp.rjRMatrix_k + 1))
        else:
            t2 = (b1 * math.log(f1)) + \
                (b2 * math.log(f2)) - (b0 * math.log(f0))

        # t3 is for the prior on rMatrices.  With the Dirichlet prior alpha
        # values all 1, t3 is log Gamma rDim
        t3 = pf.gsl_sf_lngamma(rDim)

        # t4 is for the f values.
        if var.rjRMatrixUniformAllocationPrior:
            t4 = 0.0
        else:
            # If its a uniform Dirichlet, then t4 = log k, where k is from
            # before the split
            t4 = math.log(mp.rjRMatrix_k)

        self.logPriorRatio = t1 + t2 + t3 + t4

        # Log proposal ratio
        if mp.rjRMatrix_k == 1:
            t1 = math.log(0.5)
        else:
            t1 = 0.
        if var.rjRMatrixUniformAllocationPrior:
            t2 = b0 * math.log(2.)
        else:
            t2 = (b0 * math.log(f0)) - (b1 * math.log(f1)) - \
                (b2 * math.log(f2))          # this was changed 26 sept

        # for t3, below, do some pre-calculations
        sum_uu2 = sum([u * u for u in uu])
        lastTerm = -pf.gsl_sf_lngamma(p0) + (0.5 * sum_uu2) + \
            ((rDim / 2.) * math.log(2 * math.pi))
        t3 = (s0 - s1 - s2) + ((p0 - 1.) *
                               (math.log(s1) + math.log(s2) - math.log(s0))) + lastTerm
        self.logProposalRatio = t1 + t2 + t3
        # print t1,t2,t3,s0,s1,s2

        #self.logProposalRatio = 0.

        # The Jacobian
        lastTerm = 0.5 * sum([math.log(v) for v in rm0.val])    # added 26 sept
        t1 = ((((3. * rDim) - 2.) / 2.) * math.log(s0)) - \
            ((rDim - 1.) * (math.log(s1) + math.log(s2))) + lastTerm
        t2 = (2. * rDim * math.log(bPrime0)) - \
            (rDim * (math.log(bPrime1) + math.log(bPrime2)))
        t3 = sum([uu[j] / (math.sqrt(s0 * rm0.val[j])) for j in range(rDim)])
        t3 = ((bPrime0 * (bPrime2 - bPrime1)) / (bPrime1 * bPrime2)) * t3

        if var.rjRMatrixUniformAllocationPrior:
            self.logJacobian = t1 + t2 + t3
        else:
            self.logJacobian = t1 + t2 + t3 + math.log(f0)

        # We will now make rm1 and rm2.  The rm1 will be made from rm0,
        # and rm2 will be popped from the notInPool list.  We have newVal1
        # and newVal2 which will be their vals, and we assign them to
        # nodes in beta1 and beta2, and give them rj_f values of f1 and
        # f2.
        rm1 = rm0
        for rNum in range(rDim):
            rm1.val[rNum] = newVal1[rNum]
        rm1.rj_f = f1
        for n in beta1:
            # not needed, its already that.
            n.br.parts[theProposal.pNum].rMatrixNum = rm1.num
            pf.p4_setRMatrixNum(n.cNode, theProposal.pNum, rm1.num)
        rm1.nNodes = b1

        rm2 = notInPool.pop()
        for rNum in range(rDim):
            rm2.val[rNum] = newVal2[rNum]
        rm2.rj_f = f2
        for n in beta2:
            n.br.parts[theProposal.pNum].rMatrixNum = rm2.num
            pf.p4_setRMatrixNum(n.cNode, theProposal.pNum, rm2.num)
        rm2.nNodes = b2
        rm2.rj_isInPool = True

        self.propTree.model.parts[theProposal.pNum].rjRMatrix_k += 1
        # print "...finished proposeSplitRMatrix()"

    def proposeMergeRMatrix(self, theProposal):
        gm = ['Chain.proposeMergeRMatrix()']

        mp = self.propTree.model.parts[theProposal.pNum]
        assert mp.rjRMatrix
        rDim = ((mp.dim * mp.dim) - mp.dim) / 2
        p0 = theProposal.tuning

        # Check that k is more than 1.  This should have been checked before, but check again.
        # print "rjRMatrix_k is currently %i, with %i rMatrices" %
        # (mp.rjRMatrix_k, mp.nRMatrices)
        if mp.rjRMatrix_k <= 1:
            gm.append("part %i, rjRMatrix_k = %i" %
                      (theProposal.pNum, mp.rjRMatrix_k))
            pool = [c for c in mp.rMatrices if c.rj_isInPool]
            gm.append(
                'len of pool = %i (should be the same as rjRMatrix_k)' % len(pool))
            gm.append(
                "rjRMatrix_k, the pool size, should be more than 1 for a merge.  This isn't.")
            raise P4Error(gm)

        # Choose two rMatrices (to make into one).  They must be in the pool.
        pool = [c for c in mp.rMatrices if c.rj_isInPool]
        assert len(pool) == mp.rjRMatrix_k
        rm1, rm2 = random.sample(pool, 2)
        # print "proposing to merge rMatrices %i and %i" % (rm1.num, rm2.num)

        beta1 = []
        beta2 = []
        for n in self.propTree.iterNodesNoRoot():
            theRMatrixNum = n.br.parts[theProposal.pNum].rMatrixNum
            if theRMatrixNum == rm1.num:
                beta1.append(n)
            elif theRMatrixNum == rm2.num:
                beta2.append(n)
        beta0 = beta1 + beta2
        b1 = float(len(beta1))
        b2 = float(len(beta2))
        b0 = float(b1 + b2)
        assert len(beta0) == b0
        bPrime0 = b0 + 2.
        bPrime1 = b1 + 1.
        bPrime2 = b2 + 1.
        f1 = rm1.rj_f
        f2 = rm2.rj_f
        f0 = f1 + f2

        # Obtain rMatrix proposal
        s1 = random.gammavariate(p0, 1.)
        s2 = random.gammavariate(p0, 1.)
        m1 = [v * s1 for v in rm1.val]
        m2 = [v * s2 for v in rm2.val]
        # print "m1 = ", m1
        # print "m2 = ", m2
        # print b0, b1, b2, bPrime0, bPrime1, bPrime2
        m0 = [math.exp(((bPrime1 / bPrime0) * math.log(m1[j])) + ((bPrime2 / bPrime0) * math.log(m2[j])))
              for j in range(rDim)]

        # print "m0 = ", m0
        s0 = sum(m0)

        newVal0 = [m0k / s0 for m0k in m0]

        # Log prior ratio
        # We could have a prior on the pool size, reflected in t1.  If all pool
        # sizes are equally probable, then t1 = 0
        t1 = 0.

        if var.rjRMatrixUniformAllocationPrior:
            b = len([n for n in self.propTree.iterNodes()])
            t2 = b * (math.log(mp.rjRMatrix_k) - math.log(mp.rjRMatrix_k - 1))
        else:
            t2 = (b0 * math.log(f0)) - \
                (b1 * math.log(f1)) - (b2 * math.log(f2))

        # t3 is for the prior on rMatrices.  With the Dirichlet prior alpha
        # values all 1, t3 is log Gamma rDim
        t3 = -pf.gsl_sf_lngamma(rDim)

        # t4 is for the f values.
        if var.rjRMatrixUniformAllocationPrior:
            t4 = 0.0
        else:
            # If its a uniform Dirichlet, then t4 = - log (k - 1), where k is
            # from before the merge
            t4 = -math.log(mp.rjRMatrix_k - 1)

        self.logPriorRatio = t1 + t2 + t3 + t4

        # Log proposal ratio
        if mp.rjRMatrix_k == mp.nRMatrices:  # nRMatrices is k_max
            t1 = math.log(0.5)
        else:
            t1 = 0.

        if var.rjRMatrixUniformAllocationPrior:
            t2 = - (b0 * math.log(2.))
        else:
            t2 = (b1 * math.log(f1)) + \
                (b2 * math.log(f2)) - (b0 * math.log(f0))

        # for t3, below, do some pre-calculations
        uu = [(bPrime1 / bPrime0) * math.sqrt(m0[j])
              * (math.log(m1[j]) - math.log(m0[j])) for j in range(rDim)]
        sum_uu2 = sum([u * u for u in uu])
        lastTerm = pf.gsl_sf_lngamma(p0) - (0.5 * sum_uu2) - \
            ((rDim / 2.) * math.log(2 * math.pi))
        t3 = (s1 + s2 - s0) - ((p0 - 1.) *
                               (math.log(s1) + math.log(s2) - math.log(s0))) + lastTerm

        #logSterm = ((1. - p0) * (math.log(s1) + math.log(s2) - math.log(s0)))
        # print "s1=%.1f s2=%.1f s0=%.1f    sum_uu2=%.1f  logGamma(p0)=%.1f, logSterm=%.1f" % (
        #    s1, s2, s0, sum_uu2, pf.gsl_sf_lngamma(p0), logSterm)
        self.logProposalRatio = t1 + t2 + t3

        #self.logProposalRatio = 20.

        # The Jacobian
        lastTerm = 0.5 * sum([math.log(v) for v in newVal0])  # new 26 sept
        t1 = ((((3. * rDim) - 2.) / 2.) * math.log(s0)) - \
            ((rDim - 1.) * (math.log(s1) + math.log(s2))) + lastTerm
        t2 = (2. * rDim * math.log(bPrime0)) - \
            (rDim * (math.log(bPrime1) + math.log(bPrime2)))
        t3 = sum([uu[j] / (math.sqrt(s0 * newVal0[j])) for j in range(rDim)])
        t3 = ((bPrime0 * (bPrime2 - bPrime1)) / (bPrime1 * bPrime2)) * t3
        if var.rjRMatrixUniformAllocationPrior:
            self.logJacobian = -(t1 + t2 + t3)
        else:
            self.logJacobian = -(t1 + t2 + t3 + math.log(f0))

        # Merge rm1 and rm2 => rm0, where rm0 is actually rm1, re-used, by
        # giving "0" values to rm1 = rm0
        rm1.rj_f = f0
        for rNum in range(rDim):
            rm1.val[rNum] = newVal0[rNum]
        for n in beta0:
            n.br.parts[theProposal.pNum].rMatrixNum = rm1.num
            pf.p4_setRMatrixNum(n.cNode, theProposal.pNum, rm1.num)
        rm1.nNodes = b0
        mp.rjRMatrix_k -= 1
        rm2.rj_isInPool = False
        rm2.nNodes = 0

    def proposeGdasrv(self, theProposal):

        # This is a multiplier proposal.

        gm = ["Chain.proposeGdasrv()"]
        mt = self.propTree.model.parts[
            theProposal.pNum].gdasrvs[theProposal.mtNum]

        # We can't have alpha less than about 1.e-16, or DiscreteGamma hangs.
        # But that is moot, as var.GAMMA_SHAPE_MIN is much bigger

        # Dont' do something like the following, cuz mt.val is a property that invokes a function.
        #mt.val = newVal
        #mt.val /= theProposal.tuning

        # mt.val is a numpy.ndarray type, an array with 1 element.
        assert type(mt.val) == numpy.ndarray
        oldVal = mt.val
        newVal = oldVal * \
            math.exp(theProposal.tuning * (random.random() - 0.5))

        isGood = False
        while not isGood:
            if newVal < var.GAMMA_SHAPE_MIN:
                newVal = var.GAMMA_SHAPE_MIN * var.GAMMA_SHAPE_MIN / newVal
            elif newVal > var.GAMMA_SHAPE_MAX:
                newVal = var.GAMMA_SHAPE_MAX * var.GAMMA_SHAPE_MAX / newVal
            else:
                isGood = True

        # print type(self.logProposalRatio), type(self.logPriorRatio),
        self.logProposalRatio = math.log(newVal / oldVal)

        self.logPriorRatio = 0.0
        # as in proposeBrLen()
        #self.logPriorRatio = self.mcmc.tunings.parts[theProposal.pNum].gdasrvPriorLambda * float(oldVal - newVal)
        mt.val = newVal
        assert type(mt.val) == numpy.ndarray
        # print type(self.logProposalRatio), type(self.logPriorRatio),
        # print self.logProposalRatio, self.logPriorRatio

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
            print "proposeCompLocation().  node %i, before=%i, new=%s" % (theNode.nodeNum, currentNum, proposedNum)
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

    def proposeCmd1CompDir(self, theProposal):
        gm = ['Chain.proposeCmd1CompDir()']

        # print gm[0], theProposal.pNum, theProposal.mtNum

        mp = self.propTree.model.parts[theProposal.pNum]

        # The proposal mtNum is -1, meaning do all, or any
        assert theProposal.mtNum == -1
        mt = random.choice(mp.comps)

        # mt.val is a list of floats, not a numpy.ndarray
        # print type(mt.val), type(mt.val[0]), mt.val, mt.num

        # This method uses func.dirichlet1, which is for lists not numpy
        # arrays.  A copy of inSeq is made, and the copy is modified and
        # returned.
        #dirichlet1(inSeq, alpha, theMin, theMax)
        newVal = func.dirichlet1(
            mt.val, mp.cmd1_p, var.PIVEC_MIN, 1 - var.PIVEC_MIN)

        # proposal ratio
        dirPrams = [mp.cmd1_p * v for v in newVal]
        logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(mt.val))
        dirPrams = [mp.cmd1_p * v for v in mt.val]
        logPdfProps = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(newVal))
        self.logProposalRatio = logPdfProps - logPdfCurrs

        # prior ratio
        dirPrams = [mp.cmd1_alpha * v for v in mp.cmd1_pi0]
        logPdfProps = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(newVal))
        logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(mt.val))
        self.logPriorRatio = logPdfProps - logPdfCurrs

        mt.val = newVal

    def proposeCmd1Comp0Dir(self, theProposal):
        gm = ['Chain.proposeCmd1Comp0Dir()']

        # print gm[0], theProposal.pNum, theProposal.mtNum

        mp = self.propTree.model.parts[theProposal.pNum]

        # The proposal mtNum is -1, meaning do all, or any
        assert theProposal.mtNum == -1
        curVal = mp.cmd1_pi0

        # mt.val is a list of floats, not a numpy.ndarray
        # print type(mt.val), type(mt.val[0]), mt.val, mt.num

        # This method uses func.dirichlet1, which is for lists not numpy
        # arrays.  A copy of inSeq is made, and the copy is modified and
        # returned.
        #dirichlet1(inSeq, alpha, theMin, theMax)
        newVal = func.dirichlet1(
            curVal, mp.cmd1_q, var.PIVEC_MIN, 1 - var.PIVEC_MIN)

        # proposal ratio
        dirPrams = [mp.cmd1_q * v for v in newVal]
        logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(curVal))
        dirPrams = [mp.cmd1_q * v for v in curVal]
        logPdfProps = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(newVal))
        self.logProposalRatio = logPdfProps - logPdfCurrs

        # prior ratio
        dirPrams = [mp.cmd1_s * v for v in [1.0] * mp.dim]
        logPdfProps = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(newVal))
        logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(curVal))
        self.logPriorRatio = logPdfProps - logPdfCurrs

        for pi_i in mp.comps:
            dirPrams = [mp.cmd1_alpha * v for v in newVal]
            logPdfProps = pf.gsl_ran_dirichlet_lnpdf(
                mp.dim, numpy.array(dirPrams), numpy.array(pi_i.val))
            dirPrams = [mp.cmd1_alpha * v for v in curVal]
            logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(
                mp.dim, numpy.array(dirPrams), numpy.array(pi_i.val))
            diff = logPdfProps - logPdfCurrs
            self.logPriorRatio += diff

        mp.cmd1_pi0 = newVal

    def proposeCmd1AllCompDir(self, theProposal):
        gm = ['Chain.proposeCmd1AllCompDir()']
        # all the comps, including pi0, in one go.

        mp = self.propTree.model.parts[theProposal.pNum]

        # The proposal mtNum is -1, meaning do all, or any
        assert theProposal.mtNum == -1

        # First do proposal for pi0

        # mt.val is a list of floats, not a numpy.ndarray
        # print type(mt.val), type(mt.val[0]), mt.val, mt.num
        # This method uses func.dirichlet1, which is for lists not numpy
        # arrays.  A copy of inSeq is made, and the copy is modified and
        # returned.
        #dirichlet1(inSeq, alpha, theMin, theMax)
        myU = 0.0
        pi0_newVal = func.dirichlet1(
            mp.cmd1_pi0, mp.cmd1_q, var.PIVEC_MIN, 1 - var.PIVEC_MIN, u=myU)

        # Now do proposals for all the comps in mp.comps, now using u
        # added to the dirichlet prams within func.dirichlet1().  This of
        # course needs to be taken into account when calculating the
        # proposal ratio below.
        piNewVals = [func.dirichlet1(
            mt.val, mp.cmd1_p, var.PIVEC_MIN, 1 - var.PIVEC_MIN, u=myU) for mt in mp.comps]

        # proposal ratio for pi0
        #myU = 0.0
        dirPrams = [(mp.cmd1_q * v) + myU for v in pi0_newVal]
        logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(mp.cmd1_pi0))
        dirPrams = [(mp.cmd1_q * v) + myU for v in mp.cmd1_pi0]
        logPdfProps = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(pi0_newVal))
        self.logProposalRatio = logPdfProps - logPdfCurrs

        # proposal ratios for all the comps in mp.comps
        for mtNum in range(len(mp.comps)):
            mt = mp.comps[mtNum]
            newVal = piNewVals[mtNum]
            dirPrams = [(mp.cmd1_p * v) + myU for v in newVal]
            logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(
                mp.dim, numpy.array(dirPrams), numpy.array(mt.val))
            dirPrams = [(mp.cmd1_p * v) + myU for v in mt.val]
            logPdfProps = pf.gsl_ran_dirichlet_lnpdf(
                mp.dim, numpy.array(dirPrams), numpy.array(newVal))
            self.logProposalRatio += (logPdfProps - logPdfCurrs)

        # prior ratio for pi0 component
        dirPrams = [mp.cmd1_s * v for v in [1.0] * mp.dim]
        logPdfProps = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(pi0_newVal))
        logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(
            mp.dim, numpy.array(dirPrams), numpy.array(mp.cmd1_pi0))
        self.logPriorRatio = logPdfProps - logPdfCurrs

        # prior ratios for all the comps in mp.comps
        for mtNum in range(len(mp.comps)):
            mt = mp.comps[mtNum]
            newVal = piNewVals[mtNum]
            dirPrams = [mp.cmd1_alpha * v for v in pi0_newVal]
            logPdfProps = pf.gsl_ran_dirichlet_lnpdf(
                mp.dim, numpy.array(dirPrams), numpy.array(newVal))
            dirPrams = [mp.cmd1_alpha * v for v in mp.cmd1_pi0]
            logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(
                mp.dim, numpy.array(dirPrams), numpy.array(mt.val))
            diff = logPdfProps - logPdfCurrs
            self.logPriorRatio += diff

        # assign the proposals to the prop tree
        mp.cmd1_pi0 = pi0_newVal
        for mtNum in range(len(mp.comps)):
            mt = mp.comps[mtNum]
            newVal = piNewVals[mtNum]
            mt.val = newVal

    def proposeCmd1Alpha(self, theProposal):
        gm = ['Chain.proposeCmd1Alpha()']
        MIN = 1.
        MAX = 1000.

        # print gm[0], theProposal.pNum, theProposal.mtNum

        mp = self.propTree.model.parts[theProposal.pNum]

        # The proposal mtNum is -1, meaning do all, or any
        assert theProposal.mtNum == -1
        curVal = mp.cmd1_alpha

        assert curVal <= MAX
        assert curVal >= MIN

        # mt.val is a list of floats, not a numpy.ndarray
        # print type(mt.val), type(mt.val[0]), mt.val, mt.num

        if 0:
            # Do a log scale proposal
            if 0:
                # Make proposals, and if it is outside MIN, MAX, then try
                # again.
                while 1:
                    newVal = curVal * \
                        math.exp(
                            mp.cmd1_alphaLogScaleProposalTuning * (random.random() - 0.5))
                    if (newVal < MIN) or (newVal > MAX):
                        continue
                    else:
                        break
            else:
                # If it is outside MIN, MAX, then do logarithmic reflect
                newVal = curVal * \
                    math.exp(
                        mp.cmd1_alphaLogScaleProposalTuning * (random.random() - 0.5))
                if 1:
                    # Logarithmic reflect if needed
                    while (newVal < MIN) or (newVal > MAX):
                        if newVal < MIN:
                            newVal = MIN * MIN / newVal
                        elif newVal > MAX:
                            newVal = MAX * MAX / newVal
        else:
            # Do linear proposal
            newVal = curVal + \
                (mp.cmd1_alphaLinearScaleProposalTuning *
                 (random.random() - 0.5))

            # Linear reflect if needed
            while (newVal < MIN) or (newVal > MAX):
                if newVal < MIN:
                    newVal = (MIN - newVal) + MIN
                elif newVal > MAX:
                    newVal = MAX - (newVal - MAX)

        # print curVal, newVal

        # proposal ratio
        #self.logProposalRatio = math.log(newVal/curVal)
        self.logProposalRatio = 0.0

        # prior ratio
        if 0:
            # As in Tom's blurb
            # logNormal for self
            zeta = mp.cmd1_LN_a
            sigma = mp.cmd1_LN_t
            pdfProp = pf.gsl_ran_lognormal_pdf(newVal, zeta, sigma)
            pdfCurr = pf.gsl_ran_lognormal_pdf(curVal, zeta, sigma)
            priorRatio = pdfProp / pdfCurr
            self.logPriorRatio = math.log(priorRatio)
        else:
            self.logPriorRatio = 0.0  # flat!

        for pi_i in mp.comps:
            dirPrams = [newVal * v for v in mp.cmd1_pi0]
            logPdfProps = pf.gsl_ran_dirichlet_lnpdf(
                mp.dim, numpy.array(dirPrams), numpy.array(pi_i.val))
            dirPrams = [curVal * v for v in mp.cmd1_pi0]
            logPdfCurrs = pf.gsl_ran_dirichlet_lnpdf(
                mp.dim, numpy.array(dirPrams), numpy.array(pi_i.val))
            diff = logPdfProps - logPdfCurrs
            self.logPriorRatio += diff

        mp.cmd1_alpha = newVal
