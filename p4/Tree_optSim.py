import time,os,string,sys
import pf,func
from Var import var
from Glitch import Glitch
import numpy



#def __del__(self, freeTree=pf.p4_freeTree, freeNode=pf.p4_freeNode):
#def __del__(self, freeTree=pf.p4_freeTree, dp_freeTree = pf.dp_freeTree, mysys=sys):
#def __del__(self, freeTree=pf.p4_freeTree, dp_freeTree = pf.dp_freeTree):
def __del__(self, freeTree=pf.p4_freeTree, freeNode=pf.p4_freeNode, mysys=sys):
    #mysys.stdout.write('Tree.__del__() here.\n')
    #mysys.stdout.flush()
    if hasattr(self, "splitKeyHash"): # Refers to nodes, which causes grief.
        del(self.splitKeyHash)
    self._data = None
    #self._model = None  # model is needed for freeNode()
    if 1:       # If this is not here, then nodes tend to hang around forever ...
        for n in self.nodes:
            n.wipe()
        for n in self.nodes:
            if n.cNode:
                #mysys.stdout.write('  Tree.__del__(), freeing node %i\n' % n.nodeNum)
                #mysys.stdout.flush()
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
    #mysys.stdout.flush()


def deleteCStuff(self):
    """Deletes c-pointers from nodes, self, and model, but not the data."""

    #print 'Tree.deleteCStuff() here.'
    for n in self.nodes:
        if n.cNode:
            #print '  about to free node %i, cNode %s' % (n.nodeNum, n.cNode)
            pf.p4_freeNode(n.cNode)
            n.cNode = 0
    if self.cTree:
        #print '  about to free cTree'
        pf.p4_freeTree(self.cTree)
        self.cTree = 0
    # I need to delay deleting the cModel until after deleting the
    # self.cStuff, because free-ing self.cStuff (eg nodes)
    # requires the cModel.
    if self.model and self.model.cModel:
        #print '  about to free cModel'
        pf.p4_freeModel(self.model.cModel)
        self.model.cModel = 0





##Ignore
def allocCStuff(self, resetEmpiricalComps=True):
    """Allocate c-memory for self and its nodes."""

    gm = ['Tree.allocCStuff()']

    # Make sure the nodeNums go from zero to N-1
    for i in range(len(self.nodes)):
        if self.nodes[i].nodeNum != i:
            gm.append("Programming error: Problem with node number %i." % i)
            gm.append("Nodes should be numbered consecutively from zero.")
            raise Glitch, gm

    self.modelSanityCheck(resetEmpiricalComps=resetEmpiricalComps)
    if not self.data.cData:
        self.data._setCStuff()
    if not self.model.cModel:
        self.model.allocCStuff()

    if var.doDataPart:
        #print 'about to dp_newTree'
        self.cTree = pf.dp_newTree(len(self.nodes), self.preOrder,
                                   self.postOrder, self.data.cData, self.model.cModel)
        self.doDataPart = 1
        if not self.cTree:
            gm.append("Unable to allocate a cTree")
            raise Glitch, gm
        for n in self.nodes:
            n.doDataPart = 1
            #print 'about to dp_newNode (%i)' % n.nodeNum
            cNode = pf.dp_newNode(n.nodeNum, self.cTree, n.seqNum, n.isLeaf)
            if not cNode:
                gm.append("Unable to allocate a cNode.")
                raise Glitch, gm
            n.cNode = cNode
    else:
        nLeaves = 0
        for n in self.nodes:
            if n.isLeaf:
                nLeaves += 1
        self.partLikes = numpy.zeros(self.model.nParts, numpy.float)
        self.cTree = pf.p4_newTree(len(self.nodes), nLeaves, self.preOrder,
                                   self.postOrder, self.partLikes, self.data.cData, self.model.cModel)
        if not self.cTree:
            gm.append("Unable to allocate a cTree")
            raise Glitch, gm
        for i in range(len(self.nodes)):
            n = self.nodes[i]
            if i in self.preOrder:
                inTree = 1
            else:
                inTree = 0
            # We include the inTree as a flag for whether the node
            # is in the tree or not.  If the inTree flag is 0,
            # then the node is not actually part of the tree, and so
            # clNeedsUpdating is turned off.
            n.cNode = pf.p4_newNode(n.nodeNum, self.cTree, n.seqNum, n.isLeaf, inTree)
            if not n.cNode:
                gm.append("Unable to allocate a cNode")
                raise Glitch, gm
            
    #print "finished Tree.allocCStuff()"


def setCStuff(self):
    """Transfer info about self to c-language stuff.

    Transfer relationships among nodes, the root position, branch
    lengths, model usage info (ie what model attributes apply to what
    nodes), and pre- and post-order."""

    #gm = ['Tree.setCStuff()']

    # Set node relations, br.len, root, node modelNums, preOrder?, postOrder

    # Set relations- parent, leftChild, sibling.  Here's the code for
    # pf.p4_setRelative(int theCNode, int relation, int relNum)
    # parent- relation = 0, leftChild- relation = 1, sibling- relation
    # = 2
    for n in self.nodes:
        if n.parent:
            pf.p4_setNodeRelation(n.cNode, 0, n.parent.nodeNum)
        else:
            pf.p4_setNodeRelation(n.cNode, 0, -1) # "-1" gives NULL

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
                #print "setCStuff().  about to setCompNum"
                for n in self.nodes:
                    pf.p4_setCompNum(n.cNode, pNum, n.parts[pNum].compNum)
                    if n != self.root:
                        pf.p4_setRMatrixNum(n.cNode, pNum, n.br.parts[pNum].rMatrixNum)
                        pf.p4_setGdasrvNum(n.cNode, pNum, n.br.parts[pNum].gdasrvNum)

    # pre- and postOrder
    if not self.preAndPostOrderAreValid:
        self.setPreAndPostOrder()
    
    #for i in range(len(self.nodes)):
    #    pf.p4_setPreAndPostOrder(self.cTree, i, self.preOrder[i], self.postOrder[i]) # no longer needed

    #print "finished Tree.setCStuff()"


def _commonCStuff(self, resetEmpiricalComps=True):
    """Allocate and set c-stuff, and setPrams."""
    if not self.data:
        if self.name:
            gm = ["Tree %s  (_commonCStuff)" % self.name]
        else:
            gm = ["Tree (_commonCStuff)"]
        gm.append("This tree has no data attached.  Before doing an optimization, likelihood")
        gm.append("calculation, or simulation, you need to do something like this:")
        gm.append("    theTree.data = theData")
        raise Glitch, gm

    #print "self.cTree = %s" % self.cTree
    if not self.cTree:
        # This calls self.modelSanityCheck(), which calls self.setEmpiricalComps()
        self.allocCStuff(resetEmpiricalComps=resetEmpiricalComps)
    #print "About to self.model.setCStuff()"
    self.model.setCStuff()
    #print "About to self.setCStuff()"
    self.setCStuff()
    #print "about to p4_setPrams()..."
    pf.p4_setPrams(self.cTree, -1) # "-1" means do all parts


def calcLogLike(self, verbose=1, resetEmpiricalComps=True):
    """Calculate the likelihood of the tree, without optimization."""

    self._commonCStuff(resetEmpiricalComps=resetEmpiricalComps)
    #print "about to p4_treeLogLike()..."
    self.logLike = pf.p4_treeLogLike(self.cTree, 0) # second arg is getSiteLikes
    if verbose:
        print "Tree.calcLogLike(). %f" % self.logLike



def optLogLike(self, verbose=1, newtAndBrentPowell=1, allBrentPowell=0, simplex=0):
    """Calculate the likelihood of the tree, with optimization.

    There are 3 optimization methods-- choose one.  I've made
    'newtAndBrentPowell' the default, as it is fast and seems to be
    working.  The 'allBrentPowell' optimizer used to be the default,
    as it seems to be the most robust, although it is slow.  It would
    be good for checking important calculations.  The simplex
    optimizer is the slowest, and will sometimes find better optima
    for difficult data, but often fails to optimize (with no
    warning)."""


    if verbose:
        theStartTime = time.clock()
    self._commonCStuff()

    # We want only one opt method.
    if newtAndBrentPowell:
        newtAndBrentPowell = 1
    if allBrentPowell:
        allBrentPowell = 1
    if simplex:
        simplex = 1
    if (newtAndBrentPowell + allBrentPowell + simplex) != 1:
        gm = ['Tree.optLogLike()']
        gm.append("Choose 1 opt method.")
        raise Glitch, gm
    
    # Do the opt.
    if allBrentPowell:
        pf.p4_allBrentPowellOptimize(self.cTree)
    elif simplex:
        from Tree import Tree
        pf.p4_simplexOptimize(self.cTree, self, Tree.simplexDump)
    else:
        pf.p4_newtSetup(self.cTree)
        pf.p4_newtAndBrentPowellOpt(self.cTree)

    self.logLike = pf.p4_treeLogLike(self.cTree, 0) # second arg is getSiteLikes

    # get the brLens
    brLens = pf.p4_getBrLens(self.cTree)
    for n in self.iterNodesNoRoot():
        n.br.len = brLens[n.nodeNum]

    # get the other free prams
    prams = pf.p4_getFreePrams(self.cTree)
    self.model.restoreFreePrams(prams)

    if verbose:
        print "optLogLike = %f" % self.logLike
        theEndTime = time.clock()
        print "cpu time %s seconds." % (theEndTime - theStartTime)

##Ignore
def optTest(self):
    self._commonCStuff()
    theStartTime = time.clock()
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

    print "time %s seconds." % (time.clock() - theStartTime)

##def simplexDump(self):
##    """Take a timepoint in the simplex optimization.

##    This is called by the C-language function simplexOptimize.
##    It is for taking samples of the optimization process and
##    writing them to a file (called simplex_timepoint.nex).
##    """

##    print "Tree.simplexDump here.  Doesn't work, tho."

##Ignore
def simplexDump(self):
    pass



def simulate(self, calculatePatterns=True, resetSequences=True, resetNexusSetsConstantMask=True):
    """Simulate into the attached data.

    The tree self needs to have a data and model attached.

    This week, generation of random numbers uses the C language random
    function, which is in stdlib on Linux.  It will use the same
    series of random numbers over and over, unless you tell it
    otherwise.  That means that (unless you tell it otherwise) it will
    generate the same simulated data if you run it twice.  To reset
    the randomizer, you can use func.reseedCRandomizer(), eg

    func.reseedCRandomizer(os.getpid())

    
    """
    
    self._commonCStuff()

    # If there is a NexusSets object attached to any of the alignments
    # in the Data, the constant sites mask at least will become out of sync, but we can't just
    # delete the whole nexusSets object, as they define what the parts are.
    #for a in self.data.alignments:
    # 
    #    if a.nexusSets:
    #        a.nexusSets = None

    # Probably better to do something like this
    # a.nexusSets.constant.mask = self.constantMask()
    # at the end.
    
    #print "About to pf.p4_simulate(self.cTree)"
    pf.p4_simulate(self.cTree)
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
            gm.append("resetSequences is not set, but resetNexusSetsConstantMask is set,")
            gm.append("which is probably not going to work as you want.")
            raise Glitch, gm
        

def getSiteLikes(self):
    """Likelihoods, not log likes. Placed in self.siteLikes, a list."""
    self._commonCStuff()
    self.logLike = pf.p4_treeLogLike(self.cTree, 1) # second arg is getSiteLikes
    self.siteLikes = []
    for p in self.data.parts:
        self.siteLikes += pf.getSiteLikes(p.cPart)

#def getWinningGammaCats(self):
def getSiteRates(self):
    """Get posterior mean site rate, and gamma category.

    This says two things --
    1. The posterior mean site rate, calculated like PAML
    2. Which GDASRV category contributes most to the likelihood.

    The posterior mean site rate calculation requires that there be
    only one gdasrv over the tree, which will usually be the case.

    For placement in categories, if its a tie score, then it is placed
    in the first one.

    The list of site rates, and the list of categories, both with one
    value for each site, are put into separate numpy arrays, returned
    as a list, ie [siteRatesArray, categoriesArray]

    There is one of these lists for each data partition, and the results as a
    whole are returned as a list.  So if you only have one data
    partition, then you get a 1-item list, and that single item is a list with 2
    numpy arrays.  Ie [[siteRatesArray, categoriesArray]]

    If nGammaCat for a partition is 1, it will give that partition an
    array of ones for the site rates and zeros for the categories.

    """

    self._commonCStuff()
    self.logLike = pf.p4_treeLogLike(self.cTree, 0) # second arg is getSiteLikes
    #self.winningGammaCats = []
    #for p in self.data.parts:
    #    self.winningGammaCats += pf.getWinningGammaCats(p.cPart)
    results = []

    for partNum in range(len(self.data.parts)):
        if len(self.model.parts[partNum].gdasrvs) > 1:
            gm = ['Tree.getSiteRates()']
            gm.append("Part %i has %i gdasrvs.  Maximum 1 allowed." % (
                partNum, len(self.model.parts[partNum].gdasrvs)))
            raise Glitch, gm
                      
    for partNum in range(len(self.data.parts)):
        p = self.data.parts[partNum]
        if self.model.parts[partNum].nGammaCat == 1:
            siteRates = numpy.ones(p.nChar, numpy.float)
            gammaCats = numpy.zeros(p.nChar, numpy.int32)
        elif self.model.parts[partNum].nGammaCat >  1:
            siteRates = numpy.zeros(p.nChar, numpy.float)
            gammaCats = numpy.zeros(p.nChar, numpy.int32)
            work = numpy.zeros(self.model.parts[partNum].nGammaCat, numpy.float)
            for charNum in range(p.nChar):
                gammaCats[charNum] = -1
            #pf.getWinningGammaCats(self.cTree, p.cPart, i, gammaCats, work)
            pf.getSiteRates(self.cTree, p.cPart, partNum, siteRates, gammaCats, work)
            #print siteRates
            #print gammaCats
            #print work
            if 0:
                counts = numpy.zeros(self.model.parts[partNum].nGammaCat, numpy.int32)
                for charNum in range(p.nChar):
                    counts[winningGammaCats[charNum]] += 1
                print counts
            
        else:
            raise Glitch, "This should not happen."
        results.append([siteRates, gammaCats])
    return results
