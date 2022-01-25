import sys
import string
import math
import copy
import os
import p4.func
import time
import glob
from p4.var import var
from p4.p4exceptions import P4Error
from p4.node import Node, NodeBranch, NodePart, NodeBranchPart
from p4.distancematrix import DistanceMatrix

import numpy
import p4.pf as pf
from p4.model import Model
from p4.data import Data
from p4.alignment import Part
import random

if True:
    def __del__(self, freeTree=pf.p4_freeTree, freeNode=pf.p4_freeNode, mysys=sys):
        #mysys.stdout.write('Tree.__del__() here.\n')
        # mysys.stdout.flush()
        # Refers to nodes, which causes grief.
        if hasattr(self, "splitKeyHash"):
            del(self.splitKeyHash)
        self._data = None
        # self._model = None  # model is needed for freeNode()
        # If this is not here, then nodes tend to hang around forever ...
        if 1:
            for n in self.nodes:
                n.wipe()
            for n in self.nodes:
                if n.cNode:
                    #mysys.stdout.write('  Tree.__del__(), freeing node %i\n' % n.nodeNum)
                    # mysys.stdout.flush()
                    freeNode(n.cNode)
                    n.cNode = None
            for n in self.nodes:
                del(n)
            self.root = None
        self.nodes = None
        if self.cTree:
            if self.doDataPart:
                dp_freeTree(self.cTree)
            else:
                freeTree(self.cTree)
            self.cTree = None
        #mysys.stdout.write('Tree.__del__() finished.\n')
        # mysys.stdout.flush()

    def deleteCStuff(self):
        """Deletes c-pointers from nodes, self, and model, but not the data."""

        # print 'Tree.deleteCStuff() here.'
        for n in self.nodes:
            if n.cNode:
                # print '  about to free node %i, cNode %s' % (n.nodeNum,
                # n.cNode)
                pf.p4_freeNode(n.cNode)
                n.cNode = 0
        if self.cTree:
            # print '  about to free cTree'
            pf.p4_freeTree(self.cTree)
            self.cTree = 0
        # I need to delay deleting the cModel until after deleting the
        # self.cStuff, because free-ing self.cStuff (eg nodes)
        # requires the cModel.
        if self.model and self.model.cModel:
            # print '  about to free cModel'
            pf.p4_freeModel(self.model.cModel)
            self.model.cModel = 0

    def _allocCStuff(self, resetEmpiricalComps=True):
        """Allocate c-memory for self and its nodes."""

        gm = ['Tree._allocCStuff()']

        # Make sure the nodeNums go from zero to N-1
        for i,n in enumerate(self.iterNodes()):           # allow NO_ORDER nodes at the end
            if self.nodes[i].nodeNum != i:
                gm.append(
                    "Programming error: Problem with node number %i." % i)
                gm.append("Nodes should be numbered consecutively from zero.")
                raise P4Error(gm)

        self.modelSanityCheck(resetEmpiricalComps=resetEmpiricalComps)
        if not self.data.cData:
            self.data._setCStuff()
        if not self.model.cModel:
            self.model.allocCStuff()

        if var.doDataPart:
            # print 'about to dp_newTree'
            self.cTree = pf.dp_newTree(len(self.nodes), self.preOrder,
                                       self.postOrder, self.data.cData, self.model.cModel)
            self.doDataPart = 1
            if not self.cTree:
                gm.append("Unable to allocate a cTree")
                raise P4Error(gm)
            for n in self.iterNodes():
                n.doDataPart = 1
                # print 'about to dp_newNode (%i)' % n.nodeNum
                cNode = pf.dp_newNode(
                    n.nodeNum, self.cTree, n.seqNum, n.isLeaf)
                if not cNode:
                    gm.append("Unable to allocate a cNode.")
                    raise P4Error(gm)
                n.cNode = cNode
        else:
            nLeaves = 0
            for n in self.iterNodes():
                if n.isLeaf:
                    nLeaves += 1
            self.partLikes = numpy.zeros(self.model.nParts, numpy.float)
            self.cTree = pf.p4_newTree(len(list(self.iterNodes())), nLeaves, self.preOrder,
                                       self.postOrder, var._newtAndBrentPowellOptPassLimit, self.partLikes, 
                                       self.data.cData, self.model.cModel)
            if not self.cTree:
                gm.append("Unable to allocate a cTree")
                raise P4Error(gm)
            for i in range(len(self.nodes)):
                n = self.nodes[i]
                if n.nodeNum == var.NO_ORDER:
                    continue                 # seems a little mixed up.
                if i in self.preOrder:
                    inTree = 1
                else:
                    inTree = 0
                # We include the inTree as a flag for whether the node
                # is in the tree or not.  If the inTree flag is 0,
                # then the node is not actually part of the tree, and so
                # clNeedsUpdating is turned off.
                n.cNode = pf.p4_newNode(
                    n.nodeNum, self.cTree, n.seqNum, n.isLeaf, inTree)
                if not n.cNode:
                    gm.append("Unable to allocate a cNode")
                    raise P4Error(gm)


    def setCStuff(self):
        """Transfer info about self to c-language stuff.

        Transfer relationships among nodes, the root position, branch
        lengths, model usage info (ie what model attributes apply to what
        nodes), and pre- and post-order."""

        #gm = ['Tree.setCStuff()']

        # Set node relations, br.len, root, node modelNums, preOrder?,
        # postOrder

        # Set relations- parent, leftChild, sibling.  Here's the code for
        # pf.p4_setRelative(int theCNode, int relation, int relNum)
        # parent- relation = 0, leftChild- relation = 1, sibling- relation
        # = 2
        for n in self.iterNodes():
            if n.parent:
                pf.p4_setNodeRelation(n.cNode, 0, n.parent.nodeNum)
            else:
                pf.p4_setNodeRelation(n.cNode, 0, -1)  # "-1" gives NULL

            if n.leftChild:
                pf.p4_setNodeRelation(n.cNode, 1, n.leftChild.nodeNum)
            else:
                pf.p4_setNodeRelation(n.cNode, 1, -1)

            if n.sibling:
                pf.p4_setNodeRelation(n.cNode, 2, n.sibling.nodeNum)
            else:
                pf.p4_setNodeRelation(n.cNode, 2, -1)

        # Root
        pf.p4_setTreeRoot(self.cTree, self.root.cNode)

        # br.lens
        for n in self.iterNodesNoRoot():
            #pf.p4_setBrLen(n.cNode, n.br.len, n.br.lenChanged)
            pf.p4_setBrLen(n.cNode, n.br.len)

        # Model usage info
        if self.model.isHet:
            for pNum in range(self.model.nParts):
                if self.model.parts[pNum].isHet:
                    # print "setCStuff().  about to setCompNum"
                    for n in self.iterNodes():
                        pf.p4_setCompNum(n.cNode, pNum, n.parts[pNum].compNum)
                        if n != self.root:
                            pf.p4_setRMatrixNum(
                                n.cNode, pNum, n.br.parts[pNum].rMatrixNum)
                            pf.p4_setGdasrvNum(
                                n.cNode, pNum, n.br.parts[pNum].gdasrvNum)

        # pre- and postOrder
        if not self.preAndPostOrderAreValid:
            self.setPreAndPostOrder()

    def _commonCStuff(self, resetEmpiricalComps=True):
        """Allocate and set c-stuff, and setPrams."""
        if not self.data:
            if self.name:
                gm = ["Tree %s  (_commonCStuff)" % self.name]
            else:
                gm = ["Tree (_commonCStuff)"]
            gm.append(
                "This tree has no data attached.  Before doing an optimization, likelihood")
            gm.append(
                "calculation, or simulation, you need to do something like this:")
            gm.append("    theTree.data = theData")
            raise P4Error(gm)

        #print("self.cTree = %s" % self.cTree)
        if not self.cTree:
            # This calls self.modelSanityCheck(), which calls
            # self.setEmpiricalComps()
            self._allocCStuff(resetEmpiricalComps=resetEmpiricalComps)
        #print("About to self.model.setCStuff()")
        self.model.setCStuff()
        #print("About to self.setCStuff()")
        self.setCStuff()
        #print("about to p4_setPrams()...")
        pf.p4_setPrams(self.cTree, -1)  # "-1" means do all parts
        #print("finished _commonCStuff()")

    def calcLogLike(self, verbose=1, resetEmpiricalComps=True):
        """Calculate the likelihood of the tree, without optimization."""

        self._commonCStuff(resetEmpiricalComps=resetEmpiricalComps)
        # print("about to p4_treeLogLike()...")
        # second arg is getSiteLikes
        self.logLike = pf.p4_treeLogLike(self.cTree, 0)
        if verbose:
            print("Tree.calcLogLike(). %f" % self.logLike)


    def optLogLike(self, verbose=1, method="newtAndBrentPowell"):
        """Calculate the likelihood of the tree, with optimization.

        There are different optimization methods-- choose one.  I've made
        'newtAndBrentPowell' the default, as it is fast and seems to be
        working.  The 'allBrentPowell' optimizer used to be the default,
        as it seems to be the most robust, although it is slow.  It would
        be good for checking important calculations. 

        There are two new optimization methods from the nlopt library --- 'newtAndBOBYQA' and 'BOBYQA'.  They are fast and seem to work well.

        For difficult optimizations it may help to repeat the call to optLogLike(), perhaps with a different method.
        """

        gm = ["Tree.optLogLike()"]
        if verbose:
            theStartTime = time.time()

        for n in self.iterNodesNoRoot():
            if n.br.len < var.BRLEN_MIN:
                gm.append("All branch lengths should be greater than or equal to var.BRLEN_MIN,") 
                gm.append(f"    which at the moment is {var.BRLEN_MIN}")
                gm.append(f"Got a branch length of {n.br.len:.8f} {n.br.len:g}")
                gm.append("Either make the branch length bigger, or lower var.BRLEN_MIN.")
                gm.append("You could, for example, t.stripBrLens() which makes all br lens default 0.1")
                raise P4Error(gm)

        self._commonCStuff()

        if method == "newtAndBrentPowell":
            pf.p4_newtSetup(self.cTree)
            pf.p4_newtAndBrentPowellOpt(self.cTree)
        elif method == "allBrentPowell":
            pf.p4_allBrentPowellOptimize(self.cTree)
        elif method == "newtAndBOBYQA":
            pf.p4_newtSetup(self.cTree)
            pf.p4_newtAndBOBYQAOpt(self.cTree)
        elif method == "BOBYQA":
            pf.p4_allBOBYQAOptimize(self.cTree)
        else:
            gm.append('method should be one of "newtAndBrentPowell", "allBrentPowell", "newtAndBOBYQA", or "BOBYQA"')
            raise P4Error(gm)

        # Do a final like calc.  (second arg is getSiteLikes)
        self.logLike = pf.p4_treeLogLike(self.cTree, 0)

        # get the brLens
        brLens = pf.p4_getBrLens(self.cTree)
        for n in self.iterNodesNoRoot():
            n.br.len = brLens[n.nodeNum]

        # get the other free prams
        prams = pf.p4_getFreePrams(self.cTree)
        self.model.restoreFreePrams(prams)

        if verbose:
            print("optLogLike = %f" % self.logLike)
            theEndTime = time.time()
            print("cpu time %s seconds." % (theEndTime - theStartTime))


    def optTest(self):
        self._commonCStuff()
        theStartTime = time.time()
        doXfer = 0
        for i in range(1):
            if doXfer:
                self.model.setCStuff()
                self.setCStuff()
            pf.p4_setPrams(self.cTree, -1)
            self.logLike = pf.p4_treeLogLike(self.cTree, 0)

            if doXfer:
                # get the brLens
                brLens = pf.p4_getBrLens(self.cTree)
                for i in range(len(self.nodes)):
                    n = self.nodes[i]
                    if n != self.root:
                        n.br.len = brLens[i]

                # get the other free prams
                prams = pf.p4_getFreePrams(self.cTree)
                self.model.restoreFreePrams(prams)

        print("time %s seconds." % (time.time() - theStartTime))

    def simulate(self, calculatePatterns=True, resetSequences=True, resetNexusSetsConstantMask=True, refTree=None):
        """Simulate into the attached data.

        The tree self needs to have a data and model attached.

        Generation of random numbers uses the GSL random number
        generator.  The state is held in var.gsl_rng, which is None by
        default.  If you do a simulation using this method, it will
        use ``var.gsl_rng`` if it exists, or make it if it does not exist
        yet.  When it makes it, it seeds the state based on the
        current time.  That should give you lots of variation in the
        simulations.

        If on the other hand you want to make simulations that are the
        same you can reseed the randomizer with the same seed whenever
        you do it, like this::

            if not var.gsl_rng:
                var.gsl_rng = pf.gsl_rng_get()
            # unusually, set the seed with each simulation
            mySeed = 23    # your chosen int seed
            pf.gsl_rng_set(var.gsl_rng, mySeed)

        The usual way to simulate does not use reference data.  An unusual way to
        simulate comes from (inspired by?) PhyloBayes, where the simulation is
        conditional on the original data.  It uses conditional likelihoods of
        that reference data at the root.  To turn that on, set refTree to the
        tree+model+data that you would like to use.  Calculate a likelihood with
        that refTree before using it, so that conditional likelihoods are set.
        The tree and model for refTree should be identical to the tree and model
        for self.

        Args: 

            calculatePatterns (bool): True by default. Whether to "compress" the
                newly simulated data to facilitate a faster likelihood
                calculation.

            resetSequences (bool): True by default. whether to bring the
                simulated sequences in C back into Python

            resetNexusSetsConstantMask (bool): True by default.  When
                simulations are made, the constant mask in any associated nexus
                sets will get out of sync.  Setting this to True makes a new
                mask and sets it.

            refTree (Tree): None by default.  If supplied, a tree+model+data
                which has had its likelihood calculated, where the tree+model is
                identical to self.

        """

        if refTree:
            from p4.tree import Tree
            assert isinstance(refTree, Tree)
            assert refTree.model
            assert refTree.data
            if not refTree.cTree:
                refTree.calcLogLike(verbose=False)
            assert refTree.model.cModel
            assert refTree.data.cData
            
        if not var.gsl_rng:
            var.gsl_rng = pf.gsl_rng_get()
            pf.gsl_rng_set(var.gsl_rng, int(time.time()))

        self._commonCStuff()
        if refTree:
            assert refTree.data.cData != self.data.cData
            assert refTree.data.nParts == self.data.nParts
            assert refTree.data.nTax == self.data.nTax
            for i in range(self.data.nTax):
                assert refTree.data.taxNames[i] == self.data.taxNames[i]
            assert len(refTree.data.alignments) == len(self.data.alignments)
            assert refTree.logLike, "Do a likelihood calculation with the refTree before using it here."
            # could have some more checks ...
            

        # If there is a NexusSets object attached to any of the alignments
        # in the Data, the constant sites mask at least will become out of sync, but we can't just
        # delete the whole nexusSets object, as they define what the parts are.
        # for a in self.data.alignments:
        #
        #    if a.nexusSets:
        #        a.nexusSets = None

        # Probably better to do something like this
        # a.nexusSets.constant.mask = self.constantMask()
        # at the end.

        # print "About to pf.p4_simulate(self.cTree)"
        if refTree:
            pf.p4_simulate(self.cTree, refTree.cTree, var.gsl_rng)
        else:
            pf.p4_simulate(self.cTree, 0, var.gsl_rng)
        if calculatePatterns:
            for p in self.data.parts:
                pf.makePatterns(p.cPart)
                pf.setGlobalInvarSitesVec(p.cPart)
        if resetSequences:
            self.data.resetSequencesFromParts()
            if resetNexusSetsConstantMask:
                for a in self.data.alignments:
                    if a.nexusSets:
                        a.nexusSets.constant.mask = a.constantMask()
        else:
            if resetNexusSetsConstantMask:
                gm = ['Tree.simulate().']
                gm.append(
                    "resetSequences is not set, but resetNexusSetsConstantMask is set,")
                gm.append("which is probably not going to work as you want.")
                raise P4Error(gm)

    def ancestralStateDraw(self):
        """Make a draw from the inferred root character state distribution

        This method works on a tree with an attached model and data.

        Conditional on the tree, branch lengths, model, and data, this method
        infers the ancestral character states of the root node.  However, that
        inference is probabilistic, a distribution, and this method takes a
        single draw.  It returns a string.

        """

        gm = ['Tree.ancestralStateDraw().']
        self._commonCStuff()
        self.logLike = pf.p4_treeLogLike(self.cTree, 0)
        draw = numpy.empty(4, dtype=numpy.int32)
        ancSts = []
        for pNum in range(self.data.nParts):
            dp = self.data.parts[pNum]
            ancStsPart = []
            for seqPos in range(dp.nChar):
                pf.p4_drawAncState(self.cTree, pNum, seqPos, draw)
                if draw[1] >= 0:        # gamma cat if it is a variable site, else -1  
                    assert draw[2] == 0 # not invar
                    assert draw[0] >= 0 # char num
                    ancStsPart.append(dp.symbols[draw[0]])
                elif draw[2]:           # isInvar, zero if not
                    assert draw[0] == -1
                    assert draw[1] == -1
                    assert draw[3] >= 0    # invar char num
                    ancStsPart.append(dp.symbols[draw[3]])
                else:
                    gm.append("Problem with returned draw.  Got %s" % draw)
                    raise P4Error(gm)
            assert len(ancStsPart) == dp.nChar
            ancSts.append(''.join(ancStsPart))
        return ''.join(ancSts)
                    

    def getSiteLikes(self):
        """Likelihoods, not log likes. Placed in self.siteLikes, a list."""
        self._commonCStuff()
        # second arg is getSiteLikes
        self.logLike = pf.p4_treeLogLike(self.cTree, 1)
        self.siteLikes = []
        for p in self.data.parts:
            self.siteLikes += pf.getSiteLikes(p.cPart)

    # def getSiteRates(self):
    #     """Get posterior mean site rate, and gamma category.

    #     This says two things --
    #     1. The posterior mean site rate, calculated like PAML
    #     2. Which GDASRV category contributes most to the likelihood.

    #     The posterior mean site rate calculation requires that there be
    #     only one gdasrv over the tree, which will usually be the case.

    #     For placement in categories, if its a tie score, then it is placed
    #     in the first one.

    #     The list of site rates, and the list of categories, both with one
    #     value for each site, are put into separate numpy arrays, returned
    #     as a list, ie [siteRatesArray, categoriesArray]

    #     There is one of these lists for each data partition, and the results as a
    #     whole are returned as a list.  So if you only have one data
    #     partition, then you get a 1-item list, and that single item is a list with 2
    #     numpy arrays.  Ie [[siteRatesArray, categoriesArray]]

    #     If nGammaCat for a partition is 1, it will give that partition an
    #     array of ones for the site rates and zeros for the categories.

    #     """

    #     self._commonCStuff()
    #     # second arg is getSiteLikes
    #     self.logLike = pf.p4_treeLogLike(self.cTree, 0)
    #     #self.winningGammaCats = []
    #     # for p in self.data.parts:
    #     #    self.winningGammaCats += pf.getWinningGammaCats(p.cPart)
    #     results = []

    #     for partNum in range(len(self.data.parts)):
    #         if len(self.model.parts[partNum].gdasrvs) > 1:
    #             gm = ['Tree.getSiteRates()']
    #             gm.append("Part %i has %i gdasrvs.  Maximum 1 allowed." % (
    #                 partNum, len(self.model.parts[partNum].gdasrvs)))
    #             raise P4Error(gm)

    #     for partNum in range(len(self.data.parts)):
    #         p = self.data.parts[partNum]
    #         if self.model.parts[partNum].nGammaCat == 1:
    #             siteRates = numpy.ones(p.nChar, numpy.float)
    #             gammaCats = numpy.zeros(p.nChar, numpy.int32)
    #         elif self.model.parts[partNum].nGammaCat > 1:
    #             siteRates = numpy.zeros(p.nChar, numpy.float)
    #             gammaCats = numpy.zeros(p.nChar, numpy.int32)
    #             work = numpy.zeros(
    #                 self.model.parts[partNum].nGammaCat, numpy.float)
    #             for charNum in range(p.nChar):
    #                 gammaCats[charNum] = -1
    #             #pf.getWinningGammaCats(self.cTree, p.cPart, i, gammaCats, work)
    #             pf.getSiteRates(
    #                 self.cTree, p.cPart, partNum, siteRates, gammaCats, work)
    #             # print siteRates
    #             # print gammaCats
    #             # print work
    #             if 0:
    #                 counts = numpy.zeros(
    #                     self.model.parts[partNum].nGammaCat, numpy.int32)
    #                 for charNum in range(p.nChar):
    #                     counts[winningGammaCats[charNum]] += 1
    #                 print(counts)

    #         else:
    #             raise P4Error("This should not happen.")
    #         results.append([siteRates, gammaCats])
    #     return results


