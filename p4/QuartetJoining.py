from Glitch import Glitch
import func
from func import read
from Var import var
from Tree import Tree
from Node import Node,NodeBranch
import sys
import random,  copy
import types
from Numbers import Numbers
import math

# I suppose if there are input trees whose tax sets are equal to or subsets
# of the current bigT, then they need not ever be considered again and
# could be removed from the list of input trees also. Or rather put on a list of ignores.

# On Mon, Oct 19, 2009 at 2:17 PM, mark wilkinson <mw@bmnh.org> wrote:
# > Do you have any prioritisation of the order in which leaves are added.  For
# > example you might choose to add leaves that are present in the most input
# > trees first?   --Done.


BIG_T_COLOUR = 'RED'
INPUT_T_COLOUR = 'blue'

class QJTax(object):
    def __init__(self, taxName):
        self.num = None
        self.name = taxName
        self.bitKey = None
        self.node = None         # used by getQuartetFromFourTaxa()
        self.kBrs = None
        self.nTrees = 0          # ie how many input trees have this taxon

    def dump(self):
        print "Tax.dump().       name = %s" % self.name
        print "                bitKey = %s" % self.bitKey
        if not self.kBrs:
            print "                  kBrs = %s" % self.kBrs
        else:
            print "                  kBrs = %s" % [n.nodeNum for n in self.kBrs]

######################################################################
# What is a quartet?
#
# In a theoretical description of a quartet ij|kx, all of i, j, k, and
# x are taxa.  But that is not the case here!  Here x is indeed a
# taxon, but k does not need to be defined as such, and only needs to
# be defined as a subtree on a node on bigT.  There are 3 such
# subtrees, so k is one of 0, 1, or 2.  Confusingly, I treat k
# like a taxon in some descriptions below, and so I refer to "kx" for
# example.  The k there would mean one of the taxa on the
# k-subtree in bigT, and that subtree would be identified by the
# k-number.  
#
# Since we know the node in bigT, and we know x, then a quartet is
# adequately described by a number, one of 0, 1, or 2 -- that is the
# k-number, or just k, in the code below.  So a quartet is k.
#
# It gets a little more complicated when we turn on doAddXSubTree.
# Then we need to know more -- we need to know about the input tree
# that the quartet came from, and other bits of info, as well as the
# k-number.  That is packaged up in a QJQuartet object.  It still has
# a k-number, called k, as an attribute.
#
# Below, the method that finds quartets is getQuartetsForXForNode().
# It returns a list of quartets.  Those quartets are just numbers, or
# if doAddXSubTree is turned on, then the quartets are QJQuartet
# objects.
####################################################################



class QJQuartet(object):
    # This is used only when adding x-subtrees.

    def __init__(self, theTree, theNode, k, xAndKAreAbove):
        self.tree = theTree  # the input tree, not bigT
        self.node = theNode  # the node on self.tree, with quartet ij|kx, with ij on one side and kx on the other
        self.k = k           # The k-number.-- 0, 1, or 2, indicating a subtree in bigT
        # I am not quite sure that this next bit is useful.  Might be dropped
        self.xAndKAreAbove = xAndKAreAbove         # True or False.  If kx is above then it would be True.
        # These next 2 have to do with making a smaller x-subtree,
        # making sure that it contains no polytomies nor any other
        # taxa from bigT.
        self.xNode = None    # The node on self.tree at the origin of the x-subtree
        self.up = None       # Whether the x-subtree goes up or down from the xNode.
    def dump(self):
        print "QJQuartet.dump()."
        print "node is node %i," % self.node.nodeNum,
        print " k = %i," % self.k,
        print " xAndKAreAbove = %s" % self.xAndKAreAbove
        self.tree.draw()

class QJDiagnostics(object):
    # For debugging.  If the QuartetJoining object has diagn set, eg
    # as qj.diagn = 1, then each time the quartetJoining() method is
    # called it instantiates a new QJDiagnostics object.  Completely
    # flexible.  All this will be deleted when everything works.
    def __init__(self):
        pass
    def dump(self):
        print self.__dict__


class QuartetJoining(object):
    """Mark Wilkinson's fast supertree method.

    You feed this a list of p4 Tree objects.

    For example.
    read('myTreeFile.nex')
    qj = QuartetJoining(var.trees)
    qj.useAnInputTreeAsTheStartingTree = True
    qj.doBreakTies = True
    qj.doPips=True
    qj.maxQuartetsPerN = 0
    qj.doAddXSubTree = True
    qjt = qj.quartetJoining()
    qjt.writeNexus('myQJTree.nex')

    """
    
    def __init__(self, theTrees):

        gm = ['QuartetJoining.__init__()']
        # Check the input trees.
        assert type(theTrees) == type([]), "Arg theTrees should be a list of p4 tree objects."
        for t in theTrees:
            assert isinstance(t, Tree), "Input trees should be p4 Tree objects.  Got %s" % t
            rootNChildren = t.root.getNChildren()
            if rootNChildren == 2:
                t.write()
                gm.append('This input tree has a bifurcating root.  Not allowed.')
                raise Glitch, gm
            if not t.taxNames:
                t.taxNames = [n.name for n in t.iterLeavesNoRoot()]
        self.trees = theTrees


        # Make a list of taxa that are found in all the trees.  We
        # have a list of names, and a corresponding list of Tax
        # objects, that hold more info.
        self.taxa = []  # Tax objects
        self.taxNames = []
        self.taxaDict = {}  # key is the taxName
        if theTrees:
            t = theTrees[0]
            for tName in t.taxNames:
                self.taxNames.append(tName)
                theQJTax = QJTax(tName)
                theQJTax.nTrees = 1
                self.taxa.append(theQJTax)
                self.taxaDict[tName] = theQJTax
            for t in theTrees[1:]:
                for tName in t.taxNames:
                    ret = self.taxaDict.get(tName)
                    if ret:
                        ret.nTrees += 1
                    else:
                        self.taxNames.append(tName)
                        theQJTax = QJTax(tName)
                        theQJTax.nTrees = 1
                        self.taxa.append(theQJTax)
                        self.taxaDict[tName] = theQJTax
        
        # Set the bitKey for each Tax object.
        for i in range(len(self.taxa)):
            tx = self.taxa[i]
            tx.bitKey = 1L << i
            tx.num = i

        # Calculate bitKeys and taxBits for the input trees.  Each
        # internal branch has a bitKey, that says what leaves are
        # above it.  Each tree has a taxBits that says which taxa are
        # in it.
        for t in self.trees:
            self.decorateTreeWithBitKeys(t)
            self.setTreeTaxBits(t)
            #t.loserDict = {}
            #t.nLosers = 0

        # A useful binary number.
        self.allOnes = 2L**(len(self.taxNames)) - 1

        # The supertree is made by the quartetJoining() method.
        self.bigT = None
        # This list, bigTInternalNodes, needs to be kept up to date.
        self.bigTInternalNodes = []

        ##########
        # OPTIONS
        ##########

        # Option 'verbose' says how much info to spit to the screen as
        # it goes.  Zero is off, 1 is some, 2 is lots, and 3 is also
        # lots but it stops occasionally, allowing you to read the
        # output.  Zero by default.
        self.verbose = 0

        # This 'central node strategy' is the strategy of picking an
        # internal node that maximizes the product of the number of
        # leaves in the three subtrees.  It may be a speedup, but it
        # involves extra calculations.
        self.doCentralNodeStrategy = True

        # The starting tree can be a quartet or a full input tree,
        # controlled by the following.  True means use a randomly
        # chosen input tree as the starting tree, False means start
        # with a quartet from an input tree (or trees).  If an input
        # tree is used, it must be fully resolved, and if an input
        # tree that is fully resolved cannot be found then it dies
        # with an error.
        self.useAnInputTreeAsTheStartingTree = False

        # If we do opt to use a quartet as the starting tree, by
        # setting useAnInputTreeAsTheStartingTree=False, we can either
        # use a single quartet or use all the quartets in all the
        # input trees that have those four taxa (and make a consensus
        # quartet) to start.  False means use the quartet from the
        # first randomly chosen tree only.
        self.useAllTreesInStartingQuartet = True

        # If we have tried to add all the unused taxa and failed, we
        # can choose a random taxon with few choices in where it can
        # go in bigT and simply add it.  Only leaves can be added this
        # way, not subtrees.  When a leaf has been added, it goes back
        # to normal and cycles thru the taxa one by one trying to add
        # them.  Only if none of them can be added will it attempt to
        # break a tie, and then with only 1 leaf.
        self.doBreakTies = True

        self.ifUnfinishedAddUnusedTaxaIfTheyHaveThreeKBrs = False
        #self.ifUnfinishedCollapseNodesAndAddUnusedTaxa = False

        self.addHighFrequencyTaxaFirst = False

        # Write a dot every time a taxon gets added.
        self.doPips=True

        # We only get one quartet per tree-- that is a built-in
        # limitation.  However, we can get more than one quartet if
        # they exist in different input trees.  Getting quartets is
        # expensive, and so one way to make things go faster is to
        # limit the number of quartets that we get from any N.  If
        # zero, then there is no limit, and it will get all the
        # quartets that it can, up to a maximum of 1 per input tree.
        self.maxQuartetsPerN = 4

        # We can add an x leaf one at a time, or save time and add an
        # x-subtree.  But there is a computational cost to this, which
        # might be overcome by its benefits in speed and accuracy.  Or
        # not.
        self.doAddXSubTree = False
        
        self.kk = []
        self.votes = None
        self.dbug = False
        self.slowCheck = False
        self.diagn = False
        self.didBreakTie = False

        self.taxFrequenciesNTrees = []
        self.taxFrequenciesDict = {}


#        # If doCountTries is turned on, then every time a leaf or
#        # subtree is added it prints a line saying how many
#        # x-candidates it tried, how many quartets it attempted to
#        # find from input trees, and how many taxa it added.  It will
#        # also announce if there has been a tie-breaking step.  (You
#        # probably do not want to have both this and doPips turned on,
#        # as the pips will interfere with this output.)
#        self.doCountTries = False
#        self.count_xTries = 0
#        self.count_getQuartetOnTreeTries = 0
#        self.count_txAdded = 0

        # When you ask for a quartet from a tree that has polytomies,
        # the tree might not have it, even though it has the 4 taxa.
        # That becomes a problem if we repeatedly ask the same tree
        # for the same quartet that it does not have.  The solution is
        # to check whether a tree has been asked for a quartet before,
        # and if so, then don't bother trying to get one.  That is the
        # point of doAvoidLoserQuartets.  There is a small computation
        # cost for asking the question, and so it is possible to turn
        # it off.  It does not offer any advantage if the input trees
        # are all fully resolved, because there will never be any
        # loser quartets, but on the other hand it costs less than 1%
        # of the speed to have it on when it is not needed.  One very
        # nice thing about it is that when you do multiple calls to
        # quartetJoining() from the same QuartetJoining object, the
        # compilations of loserTrees persist between reps, so it they
        # do not have to be found and added to the compilation again.

        # I don't know if this is needed or works. .
        #self.doAvoidLoserQuartets = True
        #self.count_foundLosers = 0

        #print "Got %i input trees." % len(self.trees)
        #print "There are %i taxa in all the input trees." % len(self.taxa)

    nTrees = property(lambda self: len(self.trees))
    nTax = property(lambda self: len(self.taxNames))

    def dump(self, fName=None):
        """Dump info and options about self."""

        #if not fName:
        #    f = sys.stdout
        f = sys.stdout
        f.write("QuartetJoining.dump()\n")
        f.write("---------------------\n")
        f.write("nTrees: %i\n" % self.nTrees)
        f.write("nTax: %i\n" % self.nTax)
        f.write("verbose: %s\n" % self.verbose)
        f.write("doCentralNodeStrategy: %s\n" % self.doCentralNodeStrategy)
        f.write("useAnInputTreeAsTheStartingTree: %s\n" % self.useAnInputTreeAsTheStartingTree)
        if self.useAnInputTreeAsTheStartingTree:
            f.write("inapplicable -- useAllTreesInStartingQuartet: %s\n" % self.useAllTreesInStartingQuartet)
        else:
            f.write("useAllTreesInStartingQuartet: %s\n" % self.useAllTreesInStartingQuartet)
        f.write("doBreakTies: %s\n" % self.doBreakTies)
        if self.doBreakTies:
            f.write("inapplicable -- ifUnfinishedAddUnusedTaxaIfTheyHaveThreeKBrs: %s\n" %
                    self.ifUnfinishedAddUnusedTaxaIfTheyHaveThreeKBrs)
        else:
            f.write("ifUnfinishedAddUnusedTaxaIfTheyHaveThreeKBrs: %s\n" %
                    self.ifUnfinishedAddUnusedTaxaIfTheyHaveThreeKBrs)
        f.write("addHighFrequencyTaxaFirst: %s\n" % self.addHighFrequencyTaxaFirst)
        f.write("doPips: %s\n" % self.doPips)
        f.write("maxQuartetsPerN: %i\n" % self.maxQuartetsPerN)
        f.write("doAddXSubTree: %s\n" % self.doAddXSubTree)
        f.write('---------------------\n')

        if self.trees:
            #f.write('N Trees per taxon histogram --\n')
            taxaNTrees = [tx.nTrees for tx in self.taxa]
            n = Numbers(taxaNTrees)
            if math.fabs(n.max - n.min) < 0.001:
                f.write("All %i QJ-wide taxa are each found in %i input trees.\n" % (self.nTax, n.max))
            else:
                n.histo(verbose=False, binSize=1)
                f.write("The following table shows how many input trees have each QJ-wide taxon.\n")
                f.write('  taxa    nTrees\n')
                f.write('  ----    ------\n')
                for bin in n.bins:
                    if bin[1]:
                        f.write("%5i    %6.0f\n" % (bin[1], bin[0]))

            treeSizes = [t.nTax for t in self.trees]
            n = Numbers(treeSizes)
            if math.fabs(n.max - n.min) < 0.001:
                f.write("All %i input trees have %i taxa\n" % (self.nTrees, n.max))
            else:
                n.histo(verbose=False, binSize=1)
                f.write('\n')
                f.write("In the following table, the tree size is the number of taxa in that tree.\n")
                f.write('  tree size   nTrees\n')
                f.write('  ---------   ------\n')
                for bin in n.bins:
                    if bin[1]:
                        f.write("%9i    %5.0f\n" % (bin[0], bin[1]))


            # Characterize the input trees in terms of resolution.
            # Here I am making up an index to describe how well the
            # trees are resolved.  It is the number of nodes that
            # would be needed to make it fully resolved.  So a fully
            # resolved tree has an index of zero.  A fully resolved
            # tree has nTax - 2 internal nodes, and a tree with
            # polytomies has less.

            polytomyIndices = [(t.nTax - 2) - t.nInternalNodes for t in self.trees]
            #print polytomyIndices
            n = Numbers(polytomyIndices)
            f.write('\n')
            f.write("Here the 'polytomy index' is the number of nodes that must be added to make the tree fully resolved.\n")
            if math.fabs(n.max - n.min) < 0.001:
                f.write("All %i input trees have a polytomy index of %i\n" % (self.nTrees, n.max))
            else:
                n.histo(verbose=False, binSize=1)
                f.write('  polytomy index   nTrees\n')
                f.write('  -------------   ------\n')
                for bin in n.bins:
                    if bin[1]:
                        f.write("%10i        %5.0f\n" % (bin[0], bin[1]))
            
        if f != sys.stdout:
            f.close()
        
    def quartetJoining(self, verbose=True):
        """Do The quartet joining method.

        This returns a Tree object (bigT), which either has all the
        taxa, or not.  It has a list of unusedTaxa attached, and if
        that list is empty then it has all the taxa.  On the other
        hand if the supertree is incomplete then there are Tax objects
        in that list, each of which contains a list of nodes on bigT
        on which the leftover taxon could go.

        If you set verbose=True for this method (which is different
        from setting verbose for the whole QuartetJoining instance)
        then it will tell you if it gets a complete tree or not, and
        if not it will make a table of taxa that did not make it into
        the supertree, and the branches on which they could go.

        You can test if you have a complete tree by asking the tree
        that is returned whether it has any unusedTaxa, eg if you call
        your tree 'qt' you can ask
        
            if qt.unusedTaxa:
                <it is incomplete>
            else:
                <its a complete supertree>

        """
        gm = ['QuartetJoining.quartetJoining()']

        # It doesn't make sense if ...
        if self.ifUnfinishedAddUnusedTaxaIfTheyHaveThreeKBrs and self.doBreakTies:
            gm.append("If ifUnfinishedAddUnusedTaxaIfTheyHaveThreeKBrs is set, it does not make sense to doBreakTies as well.")
            raise Glitch, gm
        #if self.ifUnfinishedCollapseNodesAndAddUnusedTaxa and self.doBreakTies:
        #    gm.append("If ifUnfinishedCollapseNodesAndAddUnusedTaxa is set, it does not make sense to doBreakTies as well.")
        #    raise Glitch, gm

        if self.useAnInputTreeAsTheStartingTree and self.addHighFrequencyTaxaFirst:
            gm.append("Sorry - useAnIputTreeAsTheStartingTree and addHighFrequencyTaxaFirst do not work together at the moment.")
            raise Glitch, gm

        if self.addHighFrequencyTaxaFirst and not self.taxFrequenciesDict:
            taxaNTrees = [tx.nTrees for tx in self.taxa]
            n = Numbers(taxaNTrees)
            if math.fabs(n.max - n.min) < 0.001:
                #print "All %i QJ-wide taxa are each found in %i input trees." % (self.nTax, n.max))
                self.taxFrequenciesNTrees = [n.max]
                self.taxFrequenciesDict = {n.max: self.taxa[:]}
            else:
                n.histo(verbose=False, binSize=1)
                #print "The following table shows how many input trees have each QJ-wide taxon."
                #print '  taxa    nTrees'
                #print '  ----    ------'
                #for bin in n.bins:
                #    if bin[1]:
                #        print "%5i    %6.0f" % (bin[1], bin[0])
                self.taxFrequenciesNTrees = [int(b[0]) for b in n.bins if b[1]]
                self.taxFrequenciesDict = {}
                for tS in self.taxFrequenciesNTrees:
                    self.taxFrequenciesDict[tS] = []
                for tx in self.taxa:
                    self.taxFrequenciesDict[tx.nTrees].append(tx)

            #print
            #for nT in self.taxFrequenciesNTrees:
            #    print "%5i    %6.0f" % (nT, len(self.taxFrequenciesDict[nT]))
        
        # kBrs is an attribute of Tax objects.  If it exists, it is a
        # list of bigT nodes where that tax could be added to bigT as
        # an x.  As such, it is specific to a particular state of
        # bigT, and so we need to wipe it if bigT has been modified.
        # That is done below whenever an x tax is added to bigT, but
        # it also needs to be done here, as this might not be the
        # first time that quartetJoining() is called from self.  So
        # initialize ...
        
        for tx in self.taxa:
            tx.kBrs = None

        # A flexible object for diagnostics.
        if self.diagn:
            self.diagn = QJDiagnostics()

        # This next method sets self.bigT
        if self.useAnInputTreeAsTheStartingTree:
            self.getAnInputTreeForTheStartingTree()
        else:
            self.getStartingTreeFromQuartet()

        self.unusedTaxa = []
        self.usedTaxa = []
        for tx in self.taxa:
            if tx.name in self.bigT.taxNames:
                self.usedTaxa.append(tx)
            else:
                self.unusedTaxa.append(tx)

        # If we have self.addHighFrequencyTaxaFirst, we want to sort
        # the unusedTaxa in order of their coverage (tx.nTrees).  But
        # we want to be somewhat random, so that repeated calls to
        # this method will not give the same order of unusedTaxa.
        # Since we only need to do this once, a hack would be to
        # shuffle the unusedTaxa, and then sort them by nTrees.  Is
        # that too slow?
        if self.addHighFrequencyTaxaFirst:
            random.shuffle(self.unusedTaxa)
            self.unusedTaxa = func.sortListOfObjectsOnAttribute(self.unusedTaxa, 'nTrees')


        self.bigTInternalNodes = [n for n in self.bigT.iterInternalsNoRoot()]
        self.bigTInternalNodes.append(self.bigT.root)

        if self.doPips:
            for i in range(self.bigT.nTax): 
                sys.stdout.write(".")
                if i > 0 and i % 50 == 0:
                    print " %4i" % i
            sys.stdout.flush()

        # So now we have a starting bigT, either an input tree or a
        # quartet taken from one or more of the input trees.
        if self.verbose:
            if self.doPips:
                print
            func.setTerminalColour(BIG_T_COLOUR)
            self.bigT.draw()
            func.unsetTerminalColour()
            print "The tree above is the starting tree, becoming self.bigT"
            #for n in self.bigT.iterNodesNoRoot():
            #    print "  node %2i %10s  %s" % (n.nodeNum, n.name, self.getTaxBitsString(n.br.bitKey))
            #print self.taxNames
            self.pause()

        # Keep adding taxa to it, until the tree is finished, or taxa
        # cannot be added.  The following is the main tree-building
        # loop.
        while len(self.bigT.taxNames) < self.nTax:
            self.didBreakTie = False
            if self.verbose >= 2:
                func.setTerminalColour('BLUE')
                print "About to tryToAddAnyX() ..."
                func.unsetTerminalColour()
                print "self.usedTaxa = %s" % [tx.name for tx in self.usedTaxa]

            #if self.doCountTries:
            #    self.count_xTries = 0
            #    self.count_getQuartetOnTreeTries = 0
            #    self.count_txAdded = 0

            if self.slowCheck: 
                self.checkBitKeys(self.bigT)
                for tx in self.unusedTaxa:
                    if tx.kBrs:
                        gm.append("tx %s has a kBrs before going into tryToAddAnyX()" % tx.name)
                        raise Glitch, gm

            # This next method returns self if it manages to add an x,
            # and None if not.  If it returns None, it has tried every
            # unused taxon, and failed.
            ret = self.tryToAddAnyX()

            #if self.doCountTries:
            #    print "xTries=%3i  getQuartetTries=%4i  txAdded=%2i  nTx=%3i" % (
            #        self.count_xTries, self.count_getQuartetOnTreeTries,
            #        self.count_txAdded, len(self.usedTaxa)),
            #    if not ret and self.doBreakTies:
            #        print "breakTie"
            #    else:
            #        print

            if ret:  # An x was added.
                #self.bigT.draw()
                #print "the tree above is the latest bigT"
                for tx in self.unusedTaxa:
                    tx.kBrs = None
            else:
                # We were unable to add any x.  Our options now
                # are limited.  If doBreakTies is turned on, then
                # we can do that, and continue on.  If that is not
                # turned on, then we are done.
                if self.slowCheck:
                    for tx in self.unusedTaxa:
                        if tx.kBrs:
                            for n in tx.kBrs:
                                if n not in self.bigT.nodes:
                                    gm.append("unused tx %s kBr %i (%s) is not in self.bigT" % (
                                        tx.name, n.nodeNum, n))
                                    raise Glitch, gm

                if self.doBreakTies:
                    self.breakTie()
                    for tx in self.unusedTaxa:
                        tx.kBrs = None
                else:
                    if self.ifUnfinishedAddUnusedTaxaIfTheyHaveThreeKBrs:
                        addedLeaves = self.addUnusedTaxaIfTheyHaveThreeKBrs()
                        if 0 and addedLeaves:
                            tDupe = self.bigT.dupe()
                            self.decorateTreeWithBitKeys(tDupe)
                            self.setTreeTaxBits(tDupe)
                            self.trees.append(tDupe)
                    #if self.unusedTaxa and self.ifUnfinishedCollapseNodesAndAddUnusedTaxa:
                    #    self.collapseNodesAndAddUnusedTaxa()

                    leftOverTaxNames = [tx.name for tx in self.unusedTaxa]
                    if verbose and self.unusedTaxa:
                        self.summarizeUnusedTaxa(leftOverTaxNames=leftOverTaxNames)
                    tnn = [tName for tName in self.taxNames if tName not in leftOverTaxNames]
                    self.bigT.taxNames = tnn # triggers bigT.checkTaxNames()
                    self.bigT.unusedTaxa = self.unusedTaxa
                    self.unusedTaxa = None # needed?
                    return self.bigT  # Incomplete
            if self.slowCheck: 
                self.checkBitKeys(self.bigT)

        if self.doPips:
            print
        if verbose:
            print "quartetJoining() Got a completed tree."
        self.bigT.taxNames = self.taxNames[:]
        self.bigT.unusedTaxa = []  # So I do not need to ask hasattr(self.bigT, 'unusedTaxa')
        #self.tNums[self.diagn.treeNum] += 1

        #if self.doAvoidLoserQuartets:
            # Show how many loser quartets each tree had.
            #for t in self.trees:
            #    print t.nLosers,

            # Say how many times we avoided trying to get the same quartet more than once.
            #print "\ncount_foundLosers = %i" % self.count_foundLosers
        #    self.count_foundLosers = 0
        return self.bigT

    def summarizeUnusedTaxa(self, leftOverTaxNames=None):
        if not leftOverTaxNames:
            leftOverTaxNames = [tx.name for tx in self.unusedTaxa]
        print "\nquartetJoining() Did not make a complete tree."
        print "Still have %i taxa left: %s" % (len(leftOverTaxNames), leftOverTaxNames)
        print "LeftoverTaxon            could be on the branch on node numbers"
        print "-------------            --------------------------------------"
        for tx in self.unusedTaxa:
            #tx.dump()
            if tx.kBrs:
                print "%-25s%s" % (tx.name, [n.nodeNum for n in tx.kBrs])
            else:
                print "%-25s Any" % tx.name        

    def addUnusedTaxaIfTheyHaveThreeKBrs(self):
        gm = ['QuartetJoining.addUnusedTaxaIfTheyHaveThreeKBrs()']

        #          +---------2:A
        # +--------1
        # |        |         +--------4:B
        # |        +---------3
        # |                  +--------5:C
        # 0
        # |--------6:D
        # |
        # +--------7:E

        # So lets say that QJ has ended with an unfinished tree
        # (obviously without breaking ties), and we have
        # myQJ.ifUnfinishedAddUnusedTaxaIfTheyHaveThreeKBrs set,
        # meaning that this method will be called.  The leftover taxa
        # if they have 3 kBrs, will allways be in a relationship like
        # branches [1,2,3] or [3,4,5] or sometimes [1,6,7] -- either
        # going up with a parent and 2 children, or going down where
        # all 3 are children of the root.

        # Or it might be like [2,3,4].  Maybe like [6,1,3].  -- ie two
        # levels.  But in those cases we do not want to do anything,
        # as it would mean collapsing existing nodes, destroying
        # information.

        # Might be like [1,3,5], which would be the same sort of thing
        # as [6,1,3] -- also don't do anything, as it is two-level.

        toDo = []
        for uT in self.unusedTaxa:
            # kBrs will be either None, implying any location, or a list of nodes.
            if uT.kBrs:
                if len(uT.kBrs) == 3:
                    #print uT.name, [n.nodeNum for n in uT.kBrs]
                    parents = set()
                    parentsOutsideKbrs = set()
                    for n in uT.kBrs:
                        if n.parent in uT.kBrs:
                            parents.add(n.parent)
                        else:
                            parentsOutsideKbrs.add(n.parent)
                    #print "got parents", [n.nodeNum for n in parents]
                    #print "got parentsOutsideKbrs", [n.nodeNum for n in parentsOutsideKbrs]

                    if len(parents) == 1:
                        # The usual case for non-roots.  Like [1,2,3]
                        # or [3,4,5] above.  One parent in kBrs, and
                        # the others are children of that parent.
                        theParent = parents.pop()
                        othersAreChildren = True
                        for n in uT.kBrs:
                            if n != theParent:
                                if n.parent != theParent:
                                    othersAreChildren = False
                        if othersAreChildren:
                            #print "parent node %i, others are children" % theParent.nodeNum
                            # This would be the case where we have a stem and 2 children, a two-pronged fork.
                            toDo.append([uT.name, theParent, uT])
                            #pass
                        else:
                            # This would be the case like [2,3,4] or [6,1,3] above.
                            pass
                    elif len(parents) == 2:
                        # This would happen in cases like [1,3,5] above.  Two level, so don't do anything.
                        pass
                    elif len(parentsOutsideKbrs) == 1:
                        # Eg [1,6,7]
                        theParentOutsideKbrs = parentsOutsideKbrs.pop()
                        if theParentOutsideKbrs == self.bigT.root:
                            othersAreChildren = True
                            for n in uT.kBrs:
                                if n.parent != theParentOutsideKbrs:
                                    othersAreChildren = False
                            if othersAreChildren:
                                #print "parentOutsideKbrs is the root, all Kbrs are children."
                                toDo.append([uT.name, self.bigT.root, uT])
                        else:
                            # How would this happen?
                            gm.append("b something wrong! -- fix me!")
                            gm.append("attempting to add taxon %s" % uT.name)
                            gm.append("kBrs %s" % [n.nodeNum for n in uT.kBrs])
                            gm.append("got parents %s" % [n.nodeNum for n in parents])
                            gm.append("got parentsOutsideKbrs %s" % [n.nodeNum for n in parentsOutsideKbrs])
                            self.bigT.draw(width=140)
                            raise Glitch, gm
                    else:
                        self.bigT.draw(width=140)
                        gm.append("attempting to add taxon %s" % uT.name)
                        gm.append("kBrs %s" % [n.nodeNum for n in uT.kBrs])
                        gm.append("got parents %s" % [n.nodeNum for n in parents])
                        gm.append("got parentsOutsideKbrs %s" % [n.nodeNum for n in parentsOutsideKbrs])
                        gm.append("c something is wrong here ...")
                        raise Glitch, gm
        for l in toDo:
            self.bigT.addSibLeaf(l[1], l[0])
            self.unusedTaxa.remove(l[2])
        return len(toDo)
        

    def collapseNodesAndAddUnusedTaxa(self):


        # This is unfinished, obviously.

        gm = ['QuartetJoining.collapseNodesAndAddUnusedTaxa()']
        print "collapseNodesAndAddUnusedTaxa() here ..."
        self.summarizeUnusedTaxa()
        self.bigT.setPreAndPostOrder()
        self.bigT.draw()
        if 0:
            for tx in self.unusedTaxa:
                # tx may not have kBrs, indicating no preference for placement.
                if tx.kBrs:
                    nodeNums = [n.nodeNum for n in tx.kBrs]
                    # Here we ask whether all the nodes in nodeNums are
                    # connected.  If they are connected, they will all have
                    # the same parent closest to the root.
                    parents = []
                    #penultimate = None
                    for n in tx.kBrs:
                        p = n
                        while p.nodeNum in nodeNums:
                            #penultimate = p
                            p = p.parent
                        if p not in parents:
                            parents.append(p)
                    if len(parents) != 1:
                        gm.append("Appears to be discontiguous branches for tx %s" % tx.name)
                        raise Glitch
                    #assert penultimate and penultimate.nodeNum in nodeNums
                    #print "kBrs %s have parent %i, penultimate %i" % (nodeNums, parents[0].nodeNum, penultimate.nodeNum)
                    print "kBrs %s have parent %i" % (nodeNums, parents[0].nodeNum)
                    #toCollapse = []
                    #for n in self.bigT.iterInternalsNoRootPostOrder():
                    #    if n.nodeNum in nodeNums:
                    #        toCollapse.append(n)
                    #for n in toCollapse[:-1]:
                    #    self.bigT.collapseNode(n)
                    #for n in tx.kBrs:
                    #    if not n.isLeaf:
                    #        self.bigT.collapseNode(n)
        if 1:
            comb = False
            toCollapseSet = set()
            for tx in self.unusedTaxa:
                if not tx.kBrs:
                    print "Add tx %s anywhere -- comb" % tx.name
                    comb = True
                    break
                else:
                    #for kNode in tx.kBrs:
                    #    if not kNode.isLeaf:
                    #        print "internal node %i" % kNode.nodeNum
                    # The rules:

                    # First do the nodes where the node's parent is not the root.  If the node's parent is not the root, and if any of the children of the node are in tx.kBrs, then do collapse those children.  If none of the children are in tx.kBrs, then don't collapse.
                    
                    # If the node's parent is the root, and the other children of the root are not in tx.kBrs, then don't collapse it.  If none of the node's children are in tx.kBrs, then don't collapse.  Otherwise do.

                    theAttachNode = None
                    for kNode in tx.kBrs:
                        if kNode.isLeaf:
                            pass
                        elif kNode in toCollapseSet:
                            if not theAttachNode:
                                theAttachNode = kNode.parent
                            pass
                        else:
                            print "considering kNode %i" % kNode.nodeNum
                            if kNode.parent == self.bigT.root:
                                doCollapse = False
                                doCollapse2 = False
                                for ch in self.bigT.root.iterChildren():
                                    if ch != kNode and ch in tx.kBrs:
                                        doCollapse = True
                                        print "kNode %i: root has more than one child in kBrs" % kNode.nodeNum
                                        break
                                if doCollapse:
                                    doCollapse2 = False
                                    for ch in kNode.iterChildren():
                                        if ch in tx.kBrs:
                                            doCollapse2 = True
                                            print "kNode %i: has a child in kBrs" % kNode.nodeNum
                                            break
                                if doCollapse and doCollapse2:
                                    print "a kNode %i, adding to collapse set" % kNode.nodeNum
                                    toCollapseSet.add(kNode)
                                    theAttachNode = kNode.parent
                            else:
                                doCollapse = False
                                if kNode.parent in tx.kBrs:
                                    for ch in kNode.iterChildren():
                                        if ch in tx.kBrs:
                                            doCollapse = True
                                            break
                                if doCollapse:
                                    print "b adding kNode %i to collapse set" % kNode.nodeNum 
                                    toCollapseSet.add(kNode)
                                    theAttachNode = kNode.parent
                    if theAttachNode:
                        print "attach node for tx %s is %s" % (tx.name, theAttachNode.parent)
                    else:
                        print "attach node for tx %s not found" % tx.name
            if comb:
                print "Comb"
            elif toCollapseSet:
                for cNode in toCollapseSet:
                    print "collapse node %i" % cNode.nodeNum
                    
            
            
    def getAnInputTreeForTheStartingTree(self):

        gm = ['QuartetJoining.getAnInputTreeForTheStartingTree()']
        if self.verbose:
            print "\nSince useAnInputTreeAsTheStartingTree is set, we start from an input tree..."

        # Old code, from when polytomous trees were not allowed ...

        # # We want to randomly choose a tree.  We only want to get
        # # a fully bifurcating tree.  The slow way to do that would
        # # be to make a list of all the input trees that are fully
        # # bifurcating and then choose one randomly.  But it is
        # # probably faster to choose a random tree and then ask
        # # whether it is fully bifurcating, until I get one.  So
        # # make a list of indices, shuffle the list, and pop from
        # # it.
        # indices = range(len(self.trees))
        # random.shuffle(indices)
        # gotIt = False
        # while not gotIt:
        #     try:
        #         randomIndex = indices.pop()
        #     except IndexError:
        #         gm.append("useAnInputTreeAsTheStartingTree is set")
        #         gm.append("However, I was unable to find any fully bifurcating input tree.")
        #         raise Glitch, gm
        #     randTree = self.trees[randomIndex]
        #     if randTree.isFullyBifurcating():
        #         gotIt = True
        #         #self.diagn.treeNum = randomIndex
        # if self.dbug:
        #     print "dbug.  Using input tree 1"
        #     randTree = self.trees[1]

        # self.bigT = randTree.dupe()
        # if not self.bigT.taxNames:
        #     self.bigT.taxNames = [n.name for n in self.bigT.iterLeavesNoRoot()]
        # #self.countSubTreeTaxa(self.bigT)

        # There are various possible strategies.  One could simply
        # take a random tree.  Or one could take the biggest tree.  Or
        # somehow the 'best' tree, defined by some criterion.  That
        # criterion would be to use the highest frequency taxa first.
        # And especially avoid low frequency taxa.

        # This is tricky, and especially tricky for polytomous trees,
        # as we don't know the content of those trees at the moment.
        # I don't think there is a good easy solution, so at the
        # moment it is not allowed if addHighFrequencyTaxaFirst is
        # set.

        # As a first go, just a random tree.
        randTree = random.choice(self.trees)
        if randTree.isFullyBifurcating():
            self.bigT = randTree.dupe()
        else:
            self.bigT = self._qj_extractFullyBifurcatingTree(randTree)  # this method dupes the randTree first
        if not self.bigT.taxNames:
            self.bigT.taxNames = [n.name for n in self.bigT.iterLeavesNoRoot()]
        self.decorateTreeWithBitKeys(self.bigT)
        self.setTreeTaxBits(self.bigT)

    def getStartingTreeFromQuartet(self):
        if self.addHighFrequencyTaxaFirst:
            self.getStartingTreeFromQuartet_highFrequencyFirst()
        else:
            self.getStartingTreeFromQuartet_randomTaxa()

    def getStartingTreeFromQuartet_highFrequencyFirst(self):
        gm = ['QuartetJoining.getStartingTreeFromQuartet_highFrequencyFirst()']

        # self.taxFrequenciesNTrees, a series of ints
        # self.taxFrequeciesDict, key is an int from self.taxFrequenciesNTrees, vals are lists of Tax objects.

        nCandidatesWanted = 10
        
        candidates = []
        indx = -1
        while len(candidates) < nCandidatesWanted:
            # potential IndexError here
            thisNTrees = self.taxFrequenciesNTrees[indx]
            candidates += self.taxFrequenciesDict[thisNTrees]
            indx -= 1
        #print "Got %i candidates" % len(candidates)
        if len(candidates) > nCandidatesWanted:
            candidates = random.sample(candidates, nCandidatesWanted)
        allSubsetsOfFour = self.allChooseFour(candidates)
        #print "got %i subsetsOfFour" % len(allSubsetsOfFour)
        bigT = None
        while allSubsetsOfFour:
            oneSubsetOfFour = random.choice(allSubsetsOfFour)
            allSubsetsOfFour.remove(oneSubsetOfFour)
            allBits = 0L
            for tx in oneSubsetOfFour:
                allBits += tx.bitKey
                #print tx.name, tx.bitKey, tx.nTrees
            #print allBits
            tt = [t for t in self.trees if (t.taxBits & allBits) == allBits]
            #print len(tt)
            if tt:
                bigT = self.getQuartetForFourTaxaInTrees(oneSubsetOfFour, tt)
            if bigT:
                break
        if bigT:
            self.bigT = bigT
            if not self.bigT.taxNames:
                self.bigT.taxNames = [n.name for n in self.bigT.iterLeavesNoRoot()]
            self.decorateTreeWithBitKeys(self.bigT)
            self.setTreeTaxBits(self.bigT)
        else:
            raise Glitch, "no bigT!"

        
    def getStartingTreeFromQuartet_randomTaxa(self):
        gm = ['QuartetJoining.getStartingTreeFromQuartet_randomTaxa()']

        if self.verbose:
            print "\nSince useAnInputTreeAsTheStartingTree is not set, we need a quartet to start with..."
        safety = 0
        self.bigT = None

        # Choose a random tree and make sure that it has at least one split.
        indices = range(len(self.trees))
        random.shuffle(indices)
        gotIt = False
        while not gotIt:
            try:
                randomIndex = indices.pop()
                #if self.dbug:
                #    randomIndex = 0
            except IndexError:
                gm.append("useAnInputTreeAsTheStartingTree is not set")
                gm.append("However, I was unable to find any input trees with internal nodes.")
                raise Glitch, gm
            randTree = self.trees[randomIndex]
            if randTree.nInternalNodes > 1: # root is an internal, but we want more
                gotIt = True
        #print "randomIndex=%i" % randomIndex
        if self.diagn:
            self.diagn.randomIndex = randomIndex
        noRootInternals = [n for n in randTree.iterInternalsNoRoot()]
        assert noRootInternals
        aNode = random.choice(noRootInternals)
        nodesAbove = [n for n in aNode.iterLeaves()]
        nodesBelow = [n for n in aNode.iterDown() if n.isLeaf]
        if len(nodesAbove) < 2:
            print
            randTree.draw()
            gm.append("Something is wrong, maybe with the randomly chosen input tree.")
            gm.append("There are only %i leaf nodes above node %i" % (len(nodesAbove), aNode.nodeNum))
            raise Glitch, gm
        if len(nodesBelow) < 2:
            print
            randTree.draw()
            gm.append("Something is wrong, maybe with the randomly chosen input tree.")
            gm.append("There are only %i leaf nodes below node %i" % (len(nodesBelow), aNode.nodeNum))
            raise Glitch, gm

        twoNodesAbove = random.sample(nodesAbove, 2)
        twoNodesBelow = random.sample(nodesBelow, 2)

        fourTaxNames = [n.name for n in twoNodesAbove]
        fourTaxNames += [n.name for n in twoNodesBelow]
        #print "fourTaxNames = %s" % fourTaxNames
        #if self.dbug:
        #    fourTaxNames = ['A', 'B', 'C', 'X']

        fourTaxa = [tx for tx in self.taxa if tx.name in fourTaxNames]
        assert len(fourTaxa) == 4

        #if self.diagn:
        #    self.diagn.firstFourTaxa = ''.join([tx.name for tx in fourTaxa])

        if self.verbose >= 2:
            func.setTerminalColour(INPUT_T_COLOUR)
            randTree.draw()
            func.unsetTerminalColour()
            print "To get a quartet, I start with a randomly chosen input tree (above),"
            print "and chose a random internal node -- in this case node %i" % aNode.nodeNum
            print "Then I randomly choose 2 taxa on one side and 2 taxa on the other side of that split."
            print "So the starting quartet will involve %s" % fourTaxNames

        # If self.useAllTreesInStartingQuartet is not set, then we
        # do not need to call makeStartingQuartetFromAllTrees(),
        # as we have the quartet already.
        if not self.useAllTreesInStartingQuartet:
            savedWarnReadNoFile = var.warnReadNoFile
            var.warnReadNoFile = 0
            read("(A, B, (C, D));")
            var.warnReadNoFile = savedWarnReadNoFile
            self.bigT = var.trees.pop()
            self.bigT.nodes[1].name = twoNodesAbove[0].name
            self.bigT.nodes[2].name = twoNodesAbove[1].name
            self.bigT.nodes[4].name = twoNodesBelow[0].name
            self.bigT.nodes[5].name = twoNodesBelow[1].name

        else:
            self.bigT = self.makeStartingQuartetFromAllTrees(randTree, twoNodesAbove, twoNodesBelow, fourTaxa)

        assert self.bigT
        self.bigT.taxNames = [n.name for n in self.bigT.iterLeavesNoRoot()]
        self.decorateTreeWithBitKeys(self.bigT)
        self.setTreeTaxBits(self.bigT)
        #self.bigT.draw()
        

        
    def makeStartingQuartetFromAllTrees(self, theOriginalTree, twoNodesAbove, twoNodesBelow, fourTaxa):
        """Having chosen fourTaxa, get quartets from all input trees, and make a staring quartet."""
        
        assert  self.useAllTreesInStartingQuartet
        bk0 = fourTaxa[0].bitKey
        bk1 = fourTaxa[1].bitKey
        bk2 = fourTaxa[2].bitKey
        bk3 = fourTaxa[3].bitKey

        allBits = bk0 | bk1 | bk2 | bk3

        treesWithTheFourTaxa = []
        for t in self.trees:
            if (t.taxBits & allBits) == allBits:
                treesWithTheFourTaxa.append(t)
        assert treesWithTheFourTaxa
        #if self.diagn:
        #    self.diagn.nTreesWithTheFourTaxa = len(treesWithTheFourTaxa)

        if self.verbose:
            print "\nuseAllTreesInStartingQuartet is set"
            print "Search all the %i input trees for quartets with those 4 taxa." % len(self.trees)
            print "Got %i input trees that had the first four taxa." % len(treesWithTheFourTaxa)
            print "(first fourTaxa are %s)" % [tx.name for tx in fourTaxa]

        # If we only got one tree in treesWithTheFourTaxa, then it
        # must be theOriginalTree, and that makes it easy, as we have
        # the quartet already in the twoNodesAbove and twoNodesBelow.
        if len(treesWithTheFourTaxa) == 1:
            assert treesWithTheFourTaxa[0] == theOriginalTree
            savedWarnReadNoFile = var.warnReadNoFile
            var.warnReadNoFile = 0
            read("(A, B, (C, D));")
            var.warnReadNoFile = savedWarnReadNoFile
            bigT = var.trees.pop()
            bigT.nodes[1].name = twoNodesAbove[0].name
            bigT.nodes[2].name = twoNodesAbove[1].name
            bigT.nodes[4].name = twoNodesBelow[0].name
            bigT.nodes[5].name = twoNodesBelow[1].name
            return bigT
        else:
            # More than one tree has the chosen quartet.
            votes = [0] * 3
            for gT in treesWithTheFourTaxa:
                if self.verbose >= 2:
                    func.setTerminalColour(INPUT_T_COLOUR)
                    gT.draw()
                    func.unsetTerminalColour()
                gotIt = False
                tNode = gT.root
                while not gotIt:
                    biggestNBits = 0
                    n = tNode.leftChild
                    badBitsBelow = False
                    while n:
                        if n.isLeaf:
                            pass
                        elif allBits & n.br.bitKey: # at least one bit
                            bitHits = []
                            nBits = 0
                            for bk in [bk0, bk1, bk2, bk3]:
                                if bk & n.br.bitKey:
                                    nBits += 1
                                    bitHits.append(bk)
                                    if nBits > 2:
                                        #print "nBits is more than 2, so its further up this subtree."
                                        break
                            if nBits > biggestNBits:
                                biggestNBits = nBits
                            if nBits == 2:

                                # Now we have n, and n.br is the split, and there are
                                # 2 bits above it.  Presumably there are 2 bits below
                                # it as well, but it might be ambiguous, in which case
                                # it is not 2 bits.  So it is worth checking.  Maybe
                                # this should be optional, with a switch.
                                badBitsBelow = False
                                bitHitsBelow = []
                                bitsBelow = (self.allOnes ^ n.br.bitKey) & gT.taxBits
                                if 0:
                                    func.setTerminalColour('violet')
                                    gT.draw()
                                    print "taxa below node %i are %s" % (
                                        n.nodeNum, self.getTaxBitsString(bitsBelow))
                                    func.unsetTerminalColour()

                                for bk in [bk0, bk1, bk2, bk3]:
                                    if bk & bitsBelow:
                                        bitHitsBelow.append(bk)
                                if len(bitHitsBelow) == 2:
                                    gotIt = True
                                    break
                                else:
                                    # ambiguous!
                                    #func.writeInColour("Got %i bitHitsBelow\n" % len(bitHitsBelow), "CYAN")
                                    badBitsBelow = True

                            elif nBits > 2: # Its further up that subtree
                                tNode = n
                                nBits = 0
                                break
                        n = n.sibling  # Try the next subtree
                    if biggestNBits < 2: # useless to go further
                        break
                    if gotIt:
                        break
                    if badBitsBelow:
                        break
                if not gotIt:
                    continue
                else:
                    # Now we have n, and n.br is the split.
                    #print "makeStartingQuartetFromAllTrees() gotIt.  node num = %i" % n.nodeNum
                    assert len(bitHits) == 2   # 2 of [bk0, bk1, bk2, bk3]
                    assert len(bitHitsBelow) == 2   # 2 of [bk0, bk1, bk2, bk3]
                    # bk0 and bk1 => vote 0   (ie bk2 and bk3)
                    # bk0 and bk2 => vote 1   (ie bk1 and bk3)
                    # bk0 and bk3 => vote 2   (ie bk1 and bk2)
                    theVoteIndex = None
                    if 0:
                        print "bk0=%s, bk1=%s, bk2=%s, bk3=%s" % (bk0, bk1, bk2, bk3)
                        print "bitHits are [%s, %s]" % (bitHits[0], bitHits[1])
                    if (bk0 in bitHits and bk1 in bitHits) or (bk2 in bitHits and bk3 in bitHits):
                        theVoteIndex = 0
                    elif (bk0 in bitHits and bk2 in bitHits) or (bk1 in bitHits and bk3 in bitHits):
                        theVoteIndex = 1
                    elif (bk0 in bitHits and bk3 in bitHits) or (bk1 in bitHits and bk2 in bitHits):
                        theVoteIndex = 2
                    if self.verbose >= 2:
                        print "the vote for the tree above is %i" % theVoteIndex
                    assert theVoteIndex != None
                    votes[theVoteIndex] += 1
            theMaxVotes = max(votes)
            maxIndexes = []
            for i in range(3):
                if votes[i] == theMaxVotes:
                    maxIndexes.append(i)
            myMaxIndex = random.choice(maxIndexes)
            if self.verbose >= 2:
                print "The votes for all %i trees are" % len(treesWithTheFourTaxa), votes
                print "Choosing vote at index %s" % myMaxIndex

            #if self.dbug:
            #    myMaxIndex = 1
            #if self.diagn:
            #    self.diagn.myMaxIndex = myMaxIndex

            savedWarnReadNoFile = var.warnReadNoFile
            var.warnReadNoFile = 0
            read('(A, B, (C, D));')
            var.warnReadNoFile = savedWarnReadNoFile
            t = var.trees.pop()
            if myMaxIndex == 0: # 0 and 1 go together, so 2 and 3 go together
                t.nodes[1].name = fourTaxa[2].name
                t.nodes[2].name = fourTaxa[3].name
                t.nodes[4].name = fourTaxa[0].name
                t.nodes[5].name = fourTaxa[1].name
            elif myMaxIndex == 1: # 0 and 2 go together, so 1 and 3 go together
                t.nodes[1].name = fourTaxa[1].name
                t.nodes[2].name = fourTaxa[3].name
                t.nodes[4].name = fourTaxa[0].name
                t.nodes[5].name = fourTaxa[2].name
            elif myMaxIndex == 2: # 0 and 3 go together, so 1 and 2 go together
                t.nodes[1].name = fourTaxa[0].name
                t.nodes[2].name = fourTaxa[3].name
                t.nodes[4].name = fourTaxa[1].name
                t.nodes[5].name = fourTaxa[2].name
            
            if self.diagn:
                self.diagn.st = t.dupe()
            return t


    def breakTie(self):

        # This method is called by quartetJoining(), if tryToAddAnyX() fails.

        # Redundantly prevent adding a subTree, by zeroing the votes
        self.votes = [0] * 3
        
        smallestKBrs = None
        for tx in self.unusedTaxa:
            if tx.kBrs:
                if smallestKBrs == None:
                    smallestKBrs = len(tx.kBrs)
                else:
                    if (len(tx.kBrs) < smallestKBrs):
                        smallestKBrs = len(tx.kBrs)
        if smallestKBrs == None:
            txWithSmallestKBrs = self.unusedTaxa
        else:
            txWithSmallestKBrs = [tx for tx in self.unusedTaxa if tx.kBrs and (len(tx.kBrs) == smallestKBrs)]
        x = random.choice(txWithSmallestKBrs)

        # It is possible that x.kBrs is None, which will give an IndexError
        try:
            myBr = random.choice(x.kBrs)
        except TypeError:
            # any branch in bigT
            if self.verbose:
                print "Warning: no quartet info, so adding leaf %s randomly!" % x.name
            myBr = random.choice([n for n in self.bigT.iterNodesNoRoot()])
        x.kBrs = [myBr]

        # This is straight from tryToAddAnyX()
        if self.verbose == 1:
            print "  +Add it (x=%s), after breaking ties" % x.name
        if self.verbose >= 2:
            if self.doAddXSubTree:
                print "By breaking a tie, we have only 1 kBr.  We do not add a subtree-- just add the leaf."
            else:
                print "By breaking a tie, we have only 1 kBr, so we can go on to add the leaf"
        #if self.doAddXSubTree:
            # addXSubTreeToBigT() uses self.votes, but self.votes
            # probably does not apply to the chosen x.  It is too much
            # bother to recalculate.  So just add the leaf.
            
            #subTreeTaxNames = self.addXSubTreeToBigT(x.kBrs[0], x)
            #self.bigT.checkTaxNames()
            #taxaInSubTree = [tx for tx in self.unusedTaxa if tx.name in subTreeTaxNames]
            ##print "taxaInSubTree = %s" % [tx.name for tx in taxaInSubTree]
            #self.usedTaxa += taxaInSubTree
            #self.unusedTaxa = [tx for tx in self.unusedTaxa if tx not in taxaInSubTree]
            ##print "self.unusedTaxa is now %s" % [tx.name for tx in self.unusedTaxa]
        #else:
        #if not (x.kBrs[0] in self.bigT.nodes):
        #    raise Glitch, "breakTie().  node %i is not in bigT." % x.kBrs[0].nodeNum
        self.didBreakTie = True
        self.addLeafToBigT(x.kBrs[0], x)
        self.usedTaxa.append(x)
        self.unusedTaxa.remove(x)

        if self.slowCheck:
            self.checkBitKeys(self.bigT)
        return self # to the quartetJoining() method.

    

    def tryToAddAnyX(self):
        """Add any one of the unused taxa.

        This is called only by the main loop in the quartetJoining()
        method.

        Its an infinite loop, and the only way to get out of it is to
        add a taxon to bigT or to run out of choices.

        First choose a random x. Having chosen an x, we invoke the
        method tryToAddAParticularX().  This will return either None
        or a list, kBrs, with at least one item in it.  KBrs are
        branches in bigT where x might be added.  Ideally, the list of
        kBrs will have only one item, the branch on the k-subtree on
        bigT where the x should be added.  But often it returns a list
        with 2 or more, meaning that it is ambiguous.  If so, we move
        on to the next x.

        If tryToAddAParticularX() returns a single kBr, then we add
        that leaf (x) to bigT at kBr.  If doAddXSubTree is turned on,
        then an x-subtree is added.

        """

        gm = ["QuartetJoining.tryToAddAnyX()"]

        if self.verbose >= 2:
            #if self.dbug:
            #    self.unusedTaxa = func.sortListOfObjectsOnAttribute(self.unusedTaxa, 'name')
            print "tryToAddAnyX(), choices = %s" % [tx.name for tx in self.unusedTaxa]
        #if self.dbug:
        #    currTxNames = [tx.name for tx in self.unusedTaxa]

        if self.doCentralNodeStrategy:
            #self.checkBitKeys(self.bigT)
            for n in self.bigTInternalNodes:
                n.pc0 = self.popcount(n.leftChild.br.bitKey)
                n.pc1 = self.popcount(n.leftChild.sibling.br.bitKey)
                if n == self.bigT.root:
                    n.pc2 = self.popcount(n.leftChild.sibling.sibling.br.bitKey) 
                else:
                    n.pc2 = self.popcount(self.bigT.taxBits & (self.allOnes ^ n.br.bitKey))
                n.pc012 = n.pc0 * n.pc1 * n.pc2
            self.bigTInternalNodes = func.sortListOfObjectsOnAttribute(self.bigTInternalNodes, 'pc012')

            if self.verbose >= 2:
                print "self.taxa are", [tx.name for tx in self.taxa]
                func.setTerminalColour(BIG_T_COLOUR)
                self.bigT.draw()
                func.unsetTerminalColour()
                print "(The tree above is bigT)"
                for n in self.bigTInternalNodes:
                    print "  node %2i " % n.nodeNum,
                    print "pc0=%s, " % n.pc0,
                    print "pc1=%s, " % n.pc1,
                    print "pc2=%s, " % n.pc2,
                    print "pc012=%i"  % n.pc012

            
            

        indices = range(len(self.unusedTaxa))
        # don't shuffle-- its too expensive!
        while 1:
                
            try:
                if self.addHighFrequencyTaxaFirst:
                    # self.unusedTaxa are sorted -- use the ones at the end first, if possible
                    myIndx = indices.pop()
                else:
                    myIndx = random.choice(indices)
            except IndexError:
                #for tx in self.unusedTaxa:
                #    if tx.kBrs:
                #        for n in tx.kBrs:
                #            if n not in self.bigT.nodes:
                #                gm.append("unused tx %s kBr %i (%s) is not in self.bigT" % (
                #                            tx.name, n.nodeNum, n))
                #                raise Glitch, gm
                return None
            
            if 0: #self.dbug:
                if 'W' in currTxNames:
                    self.bigT.write()
                    if currTxNames.index('W') in indices:
                        myIndx = currTxNames.index('W')
                elif 'Z' in currTxNames:
                    if len(currTxNames) == 2:
                        if currTxNames.index('Z') in indices:
                            self.bigT.write()
                            myIndx = currTxNames.index('Z')
                    #print "dbug: considering adding tax 'X'"

            #print "indices =%s, myIndx = %s" % (indices, myIndx)

            if self.addHighFrequencyTaxaFirst:
                pass # its already pop'ed, above
            else:
                indices.remove(myIndx)
            x = self.unusedTaxa[myIndx]
            if self.verbose >= 2:
                print "Considering adding leaf x=%s" % x.name
                if self.verbose >= 3:
                    self.pause()

            # The return value of tryToAddAParticularX() is kBrs, a
            # list of nodes in bigT.  Ideally with only one item, but
            # of course often with more.
            kBrs = self.tryToAddAParticularX(x)

            if 0: # slow check for debug
                if kBrs:
                    for n in kBrs:
                        if n not in self.bigT.nodes:
                            raise Glitch, "tryToAddAnyX().  node %i (%s) is not in self.bigT" % (n.nodeNum, n)
                for tx in self.unusedTaxa:
                    if tx.kBrs:
                        for n in tx.kBrs:
                            if n not in self.bigT.nodes:
                                gm.append("current x=%s" % x.name)
                                gm.append("unusedTaxa=%s" % [tx2.name for tx2 in self.unusedTaxa]) 
                                gm.append("x unused tx %s kBr %i (%s) is not in self.bigT" % (
                                        tx.name, n.nodeNum, n))
                                raise Glitch, gm
            
            #if self.diagn:
            #    if x.name == 'W':
            #        if kBrs:
            #            self.diagn.nKbrsW = len(kBrs)
            #        else:
            #            self.diagn.nKbrsW = 0
            #    if x.name == 'Z':
            #        if kBrs:
            #            self.diagn.nKbrsZ = len(kBrs)
            #        else:
            #            self.diagn.nKbrsZ = 0
                
            if kBrs == None:
                if self.verbose >= 1:
                    if indices:
                        print "  +kBrs is None.  Try another x"
                    else:
                        print "  +kBrs is None.  There are no more x candidates. "
                continue
            if len(kBrs) == 0:
                raise Glitch, "len(kBrs) is zero.  This should not happen.  Programming error."
            if len(kBrs) == 1:
                if self.verbose == 1:
                    print "  +Add it."
                if self.verbose >= 2:
                    if self.doAddXSubTree:
                        print "Back in tryToAddAnyX().  We have only 1 kBr, so we can go on to add the subtree"
                    else:
                        print "Back in tryToAddAnyX().  We have only 1 kBr, so we can go on to add the leaf"
                if self.doAddXSubTree:
                    subTreeTaxNames = self.addXSubTreeToBigT(kBrs[0], x)
                    self.bigT.checkTaxNames()
                    taxaInSubTree = [tx for tx in self.unusedTaxa if tx.name in subTreeTaxNames]
                    if self.diagn:
                        if hasattr(self.diagn, 'afters'):
                            self.diagn.afters.append(''.join([tx.name for tx in taxaInSubTree]))
                        else:
                            self.diagn.afters = [''.join([tx.name for tx in taxaInSubTree])]
                    #print "taxaInSubTree = %s" % [tx.name for tx in taxaInSubTree]
                    self.usedTaxa += taxaInSubTree
                    self.unusedTaxa = [tx for tx in self.unusedTaxa if tx not in taxaInSubTree]
                    #print "self.unusedTaxa is now %s" % [tx.name for tx in self.unusedTaxa]
                else:
                    self.addLeafToBigT(kBrs[0], x)
                    self.usedTaxa.append(x)
                    self.unusedTaxa.remove(x)
                if self.slowCheck: 
                    self.checkBitKeys(self.bigT)
                return self # to the quartetJoining() method.
            else:
                if self.verbose == 1:
                    print "  +Can't add it (x=%s) now." % x.name
                if self.verbose >= 2:
                    print "We have more than 1 kBr, so attach that list to x, as x.kBrs"
                x.kBrs = kBrs

        
        
    def tryToAddAParticularX(self, x):
        """Given a particular x, try to add it.

        This method tries to add the particular x 'centered' on all
        the internal nodes in bigT, until success or until it runs out
        of internal nodes in bigT.

        It returns a list of branches on bigT where x could be added.
        The list contains at least 1 item, but may contain more.
        """


        #if self.doCountTries:
        #    self.count_xTries += 1
        #self.diagn.calls += 1

        goodInternalNodeNums = [n.nodeNum for n in self.bigTInternalNodes]
        if not self.doCentralNodeStrategy:
            random.shuffle(goodInternalNodeNums)

        #if self.dbug:
        #    if x.name == 'W':
        #        print "goodInternalNodeNums", goodInternalNodeNums
        if self.verbose >= 2:
            print "goodInternals (in bigT):"
            func.setTerminalColour(BIG_T_COLOUR)
            self.bigT.draw()
            func.unsetTerminalColour()
            print "(The tree above is bigT)"
            if self.doCentralNodeStrategy:
                for n in self.bigTInternalNodes:
                    print "  node %2i  pc0=%2i, pc1=%2i, pc2=%2i, pc012=%2i" % (
                        n.nodeNum, n.pc0, n.pc1, n.pc2, n.pc012)
            else:
                for i in goodInternalNodeNums:
                    print "  node %2i" % i
            if self.verbose >= 3:
                self.pause()
        #print "used=%s, choices=%s" % ([tx.name for tx in self.usedTaxa], [tx.name for tx in choices]) 

        # In case we need to refineKSubTree() below, we keep track of
        # the internalNodesWithNoQuartets, that is nodes from
        # goodInternals that we try to find quartets for but there are
        # none in the input trees.  refineKSubTree() also needs to
        # look for quartets centered on various nodes, and there is no
        # point in looking for them again if they do not exist.
        internalNodesWithNoQuartets = []

        # kk is a list of quartets, either as numbers (0, 1, or 2) or QJQuartet objects.
        self.kk = []
        while not self.kk:
            if self.verbose >= 2:
                print "  tryToAddAParticularX(x=%s).  " % x.name,
                print "goodInternalNodeNums are now %s" % goodInternalNodeNums
            if 0 and self.dbug and x.name == 'Z':
                theNodeNum = 3
                func.setTerminalColour('VIOLET')
                self.bigT.draw()
                func.unsetTerminalColour()
                print "dbug: goodInternalNodeNums", goodInternalNodeNums
                print "dbug: choosing node %i in bigT" % theNodeNum
                goodInternalNodeNums.remove(theNodeNum)
                theNode = self.bigT.nodes[theNodeNum]
                if self.diagn:
                    self.diagn.tZ = self.bigT.dupe()
            else:
                try:
                    theNodeNum = goodInternalNodeNums.pop()  # from the end, ie highest pc012.
                    theNode = self.bigT.nodes[theNodeNum]
                except IndexError:
                    #print "No more goodInternals.  Not an error. It just means choose another x"
                    return None  # to tryToAddAnyX(), above.
            if self.verbose >= 1:
                print "Try to add x=%s 'centered on' node %i" % (x.name, theNode.nodeNum)
                if self.verbose >= 3:
                    self.pause()
            #if self.dbug:
            #    if x.name == 'W':
            #        print "theNodeNum is", theNodeNum
            self.kk = self.getQuartetsForXForNode(x, theNode)
            #print "xxxxx got %i kk" % len(self.kk)
            #if self.dbug:
            #    if x.name == 'W':
            #        if self.kk:
            #            print "theNodeNum %i, len(self.kk) = %i, k=%i" % (theNodeNum, len(self.kk), self.kk[0].k)
            #        else:
            #            print "theNodeNum %i, self.kk = %s" % (theNodeNum, self.kk)
            if not self.kk: # unable to find an appropriate tree, so try another theNode
                internalNodesWithNoQuartets.append(theNode)
                if self.verbose >= 1:
                    print "  Can't get any quartets centered on that node.  Try another node."
                if self.verbose >= 2:
                    if not goodInternalNodeNums:
                        #print "  Failed to get a quartet, at all, for x=%s.  Bad! -- Glitch follows!" % x.name
                        print "  ...except that there are no other nodes.  ",
                        print "So we have failed to add x=%s on any node." % x.name
                        return None
                    else:
                        # This happens a lot, both for qj's that will be successful or will fail.
                        print "  Failed to get a quartet, so choose another theNode from goodInternals"
                        if self.verbose >= 3:
                            func.setTerminalColour(BIG_T_COLOUR)
                            self.bigT.draw()
                            func.unsetTerminalColour()
                            print "  (The tree above is the current bigT)"


        assert self.kk  # A list, one for each tree that had a k, unless limited by maxQuartetsPerN
        #print "got %i quartets" % len(self.kk)
        #print self.kk
        #for k in self.kk:
        #    k.dump()
        #sys.exit()
        
        # We have self.kk, so find kNodes and kBrs.
        if self.verbose >= 1:
            print "  Maybe able to add x=%s.  Got quartets from %i trees.  Now see how big the k-subtree is." % (
                x.name, len(self.kk))
            #print " ",
        kIntNodes, kBrs = self.getKNodes(theNode)
        if self.verbose >= 2:
            print "    After assessing the first quartet, we have %i kBrs, %s" % (
                len(kBrs), [n.nodeNum for n in kBrs])
            #print "    After assessing the first quartet, we have %i kBrs" % len(kBrs)
            print "    and %i kIntNodes, %s" % (len(kIntNodes), [n.nodeNum for n in kIntNodes])
            if self.verbose >= 3:
                func.setTerminalColour(BIG_T_COLOUR)
                self.bigT.draw()
                func.unsetTerminalColour()
                print "(The tree above is bigT, in which we find the kIntNodes %s)" % (
                    [n.nodeNum for n in kIntNodes])
            if len(kBrs) > 1:
                print "    So refineKSubTree() it."
        if self.verbose == 1:
            if len(kBrs) > 1:
                print "    k-subtree is %i kBrs.  Need to refine the k-subtree." % len(kBrs)
            elif len(kBrs) == 1:
                print "    k-subtree is %i kBrs.  Lucky-- can add it right now." % len(kBrs)

        #if self.dbug:
        #    if x.name == 'Z':
        #        pass
        #        #print "** len(kBrs)", len(kBrs), "len(kIntNodes)", len(kIntNodes)

        #print "len(kBrs) == %i." % len(kBrs)
        #sys.exit()
        if len(kBrs) == 1:
            return kBrs  # to tryToAddAnyX()
        elif len(kBrs) > 1 and not kIntNodes: # This can happen with a tie score in getKNodes()
            return kBrs  # to tryToAddAnyX()
        elif not kBrs:
            raise Glitch, "No kBrs.  %i kIntNodes.  Programming error." % len(kIntNodes)
        else:
            #if not kIntNodes or not kBrs:
            #    self.bigT.draw()
            #    print "kIntNodes %s" % [n.nodeNum for n in kIntNodes]
            #    print "kBrs %s" % [n.nodeNum for n in kBrs]
            #    raise Glitch, "There are %i kIntNodes and %i kBrs" % (len(kIntNodes), len(kBrs))
            kBrs = self.refineKSubTree(kIntNodes, kBrs, x, internalNodesWithNoQuartets)
            if self.verbose >= 2:
                print "    refineKSubTree() returned kBrs = %s" % [n.nodeNum for n in kBrs]
                self.pause()
            assert kBrs, "no kBrs.  Programming error. How did *that* happen?!?"
            return kBrs  # to tryToAddAnyX()
        



    def refineKSubTree(self, kIntNodes, kBrs, x, internalNodesWithNoQuartets):
        #print "refineKSubTree() start. got %i kIntNodes, %i kBrs" % (len(kIntNodes), len(kBrs))
        if self.doCentralNodeStrategy:
            kIntNodes = func.sortListOfObjectsOnAttribute(kIntNodes, 'pc012')
        else:
            random.shuffle(kIntNodes)
        #print "Starting refineKSubTree() with %i kIntNodes" % len(kIntNodes)
        while 1:

            if self.verbose == 1:
                print "    Refining. %i kBrs" % len(kBrs)
                
            if self.verbose >= 2:
                print "      refineKSubTree()  kIntNodes are now %s" % [n.nodeNum for n in kIntNodes]
            try:
                #print "x%i" % len(kIntNodes),
                #sys.stdout.flush()
                theIntNode = kIntNodes.pop()
            except IndexError:
                #print "No more kIntNodes.  Don't bother re-looping?"
                #print "Z"
                if self.verbose >= 2:
                    print "      refineKSubTree()  no more kIntNodes.  Returning kBrs."
                return kBrs
            if theIntNode in internalNodesWithNoQuartets:
                if self.verbose >= 2:
                    print "      Ignoring node %i -- it is in internalNodesWithNoQuartets" % theIntNode.nodeNum 
            else:
                if self.verbose >= 2:
                    #self.bigT.draw()
                    print "      Continue to try to add x=%s, this time via a quartet centered on %i" % (
                        x.name, theIntNode.nodeNum)
                self.kk = self.getQuartetsForXForNode(x, theIntNode)
                if not self.kk:
                    #print "y",
                    if self.verbose >= 2:
                        print "      Failed to get a quartet, so choose another theIntNode from kIntNodes."
                else:
                    oldKBrs = kBrs
                    oldKIntNodes = kIntNodes

                    # getKNodes gets all of the kBrs and kIntNodes on the k-subtree.
                    newKIntNodes, newKBrs = self.getKNodes(theIntNode)

                    # For both kBrs and kIntNodes, get the intersection of the old and new sets.
                    kBrs = []
                    for kN in oldKBrs:
                        if kN in newKBrs:
                            kBrs.append(kN)
                    kIntNodes = []
                    for kN in oldKIntNodes:
                        if kN in newKIntNodes:
                            kIntNodes.append(kN)

                    # We should be refining here-- making the list smaller. Is it smaller?
                    if len(kBrs) > len(oldKBrs):  # should it be just >?  No, often it stays the same.
                        print "refineKSubTree().  Added a quartet, but the number of kBrs did not decrease."

                    # This next 'else' happens a lot
                    #if len(kIntNodes) < len(oldKIntNodes):
                    #    pass
                    #else:
                    #    print "refineKSubTree().  Added a quartet, but the number of kIntNodes did not decrease."

                    if len(kBrs) == 1:
                        return kBrs
                    elif not kBrs:
                        print "refineKSubTree().  No kBrs.  How did *that* happen?!?"
                        sys.exit()
                    else:
                        if self.verbose >= 2:
                            print "      after getting set intersections, got %i kIntNodes and %i kBrs" % (
                                len(kIntNodes), len(kBrs))
                            print "      kBrs %s" % [n.nodeNum for n in kBrs]

    def addLeafToBigT(self, theNode, x):
        # Add a leaf on the branch on the chosenNode.  That causes
        # 2 nodes to be added -- one internal, and one leaf.  For
        # both, the bitKeys need to be added.
        #print "adding taxon %s ..." % x.name

        if self.verbose >= 2:
            print "Before adding the leaf, bigT is:"
            func.setTerminalColour(BIG_T_COLOUR)
            self.bigT.draw()
            func.unsetTerminalColour()

        # The Tree.addLeaf() method adds 2 nodes, one for the leaf,
        # and one on the branch leading from theNode.  The leaf node
        # is returned (here addedNode).  The other new node is
        # unambiguously its parent (here theParent).
        addedNode = self.bigT.addLeaf(theNode, x.name)
        addedNode.br.bitKey = x.bitKey

        theParent = addedNode.parent
        self.bigTInternalNodes.append(theParent)
        #self.checkBigTInternalNodes() # slow check
        p = theParent
        while p != self.bigT.root:
            self.calculateBitKeyForInternalNodeFromChildren(p)
            p = p.parent

        self.bigT.taxBits = self.bigT.taxBits | x.bitKey


        if self.verbose:
            print "    addLeafToBigT: Add a leaf x=%s on the branch for node %i." % (
                x.name, theNode.nodeNum)
            if self.verbose >= 2:
                func.setTerminalColour(BIG_T_COLOUR)
                self.bigT.draw()
                func.unsetTerminalColour()
                print "... resulting in the tree above, with %i leaves." % self.bigT.nTax
                self.pause()

        #if self.doCountTries:
        #    self.count_txAdded = 1

        if self.doPips:
            if self.didBreakTie:
                sys.stdout.write('*')
            else:
                sys.stdout.write(".")
            if self.bigT.nTax % 50 == 0:
                print " %4i" % self.bigT.nTax
            sys.stdout.flush()

    def addXSubTreeToBigT(self, theNode, x):

        # theNode is a node in bigT, where an x-subtree may now be added.

        gm = ["addXSubTreeToBigT(theNodeOnBigT=%i, x=%s)" % (theNode.nodeNum, x.name)]
        if self.verbose >= 2:
            print "Now doing %s" % gm[0]
        # Use only those trees that had the high vote in the choice of the kBr.
        maxVotes = max(self.votes)
        assert maxVotes
        nMaxs = self.votes.count(maxVotes)
        #maxs = []
            
        #for i in range(3):
        #    if self.votes[i] == maxVotes:
        #        maxs.append(i)
        if nMaxs > 1:
            # We only call addXSubTreeToBigT() if there is only 1 kBr,
            # and that would only happen if there was only 1 maxVotes
            # in self.votes.  So if we are here, with nMaxs more
            # than 1, we must have made a tie-breaking choice.  And
            # that should have been prevented.
            gm.append("nMaxs is %i.  Odd.  Programming error. Due to tie-breaking?" % nMaxs)
            gm.append("self.votes = %s" % self.votes)
            raise Glitch, gm
        theMaxNum = self.votes.index(maxVotes)
        
        # Find all the quartets that had the maxVote, where theQJQuartet.k == theMaxNum
        maxVoteQuartets = [q for q in self.kk if q.k == theMaxNum]


        if not maxVoteQuartets:
            for i in range(len(self.kk)):
                q = self.kk[i]
                gm.append("addXSubTreeToBigT(), q.x.name=%s, q.vote=%s" % (q.x.name, q.vote))
            gm.append("addXSubTreeToBigT()  self.votes=%s, maxVotes=%s" % (self.votes, maxVotes))
            gm.append("addXSubTreeToBigT() maxs = %s" % maxs)
            gm.append("no maxVoteQuartets -- how can that be?")
            raise Glitch, gm
                      
        if 0:
            for i in range(len(maxVoteQuartets)):
                q = maxVoteQuartets[i]
                q.tree.draw()
                print "Above is the tree for maxVoteQuartets number %i" % i
            print "got %i maxVoteQuartets" % len(maxVoteQuartets)
            print "adding taxon x=%s, and its subtree." % x.name
            #sys.exit()

        # Now we want to find the QJQuartet with biggest x-subtree.
        # There might be more than one QJQuartet with a tree with the
        # biggest x-subtree -- the way this loop is set up is it finds
        # and keeps the first one, and ignores subsequent ones the
        # same size.
        qWithBiggestSubTree = None
        nLeavesInBiggestSubTree = 0
        for q in maxVoteQuartets:
            if not qWithBiggestSubTree:  # the first one, so it must be the biggest so far
                qWithBiggestSubTree = q
                nLeavesInBiggestSubTree = self.locateAndCountLeavesInXSubTreeFromQuartet(q, x)
            else:
                thisNLeaves = self.locateAndCountLeavesInXSubTreeFromQuartet(q, x)
                if thisNLeaves > nLeavesInBiggestSubTree:
                    nLeavesInBiggestSubTree = thisNLeaves
                    qWithBiggestSubTree = q
        if self.verbose >= 2:
            print "The (first) tree with the biggest sub-tree (with %i leaves) tree follows:" % (
                nLeavesInBiggestSubTree)
            func.setTerminalColour(INPUT_T_COLOUR)
            qWithBiggestSubTree.tree.draw()
            func.unsetTerminalColour()

        theSubTree = qWithBiggestSubTree.tree.dupeSubTree(qWithBiggestSubTree.xNode, up=qWithBiggestSubTree.up)
        assert theSubTree.isFullyBifurcating()
        self.decorateTreeWithBitKeys(theSubTree)
        #theSubTree.draw()
        #print "The above is the subtree."
        #sys.exit()

        # Calculate the taxBits for theSubTree
        subTreeTaxBits = qWithBiggestSubTree.xNode.br.bitKey
        if not qWithBiggestSubTree.up:
            subTreeTaxBits = (self.allOnes ^ subTreeTaxBits) & qWithBiggestSubTree.tree.taxBits # flip bits

        # Slow check of theSubTree.taxBits
        if 1:
            checkTaxBits = 0L
            for n in theSubTree.iterLeavesNoRoot():
                checkTaxBits = checkTaxBits | n.br.bitKey
            if subTreeTaxBits != checkTaxBits:
                raise Glitch, "subTreeTaxBits are wrong.  Existing=%s, check=%s" % (
                    self.getTaxBitsString(subTreeTaxBits), self.getTaxBitsString(checkTaxBits))
        
        # theSubTree does not last long -- it is destroyed when it is
        # added to bigT.  But I need the taxNames and the taxBits.  So
        # save them.  And the internalNodes
        subTreeTaxNames = [n.name for n in theSubTree.iterLeavesNoRoot()]
        assert len(subTreeTaxNames) == nLeavesInBiggestSubTree
        #subTreeTaxBits = theSubTree.taxBits
        subTreeInternalNodes = [n for n in theSubTree.iterInternalsNoRoot()]
        subTreeInternalNodes.append(theSubTree.root)

        if self.verbose >= 2:
            #if not qWithBiggestSubTree.up:
            if 1:
                func.setTerminalColour(BIG_T_COLOUR)
                self.bigT.draw()
                func.unsetTerminalColour()
                print "addXSubTreeToBigT().  The tree above is the current bigT"
                func.setTerminalColour('violet')
                theSubTree.draw()
                func.unsetTerminalColour()
                print "addXSubTreeToBigT().  The tree above is the sub-tree"
                print "Add it to node %i in bigT" % theNode.nodeNum
                print "addXSubTreeToBigT().  nLeaves = %2i, up=%5s" % (
                    nLeavesInBiggestSubTree, qWithBiggestSubTree.up)
                print "subTreeTaxBits = %s" % self.getTaxBitsString(subTreeTaxBits)
                #sys.exit()
        
        if self.diagn:
            self.xSubTree = nLeavesInBiggestSubTree
        # Doing the following zaps theSubTree.
        oldBigTNTax = self.bigT.nTax
        self.bigT.addSubTree(theNode, theSubTree, subTreeTaxNames)
        self.bigTInternalNodes += subTreeInternalNodes # includes root of subtree
        #self.checkBigTInternalNodes() # slow check
        theNode.parent.name = None
        theNode.sibling.br.bitKey = subTreeTaxBits
        

        #self.calculateBitKeyForInternalNodeFromChildren(theNode.parent)
        p = theNode.parent
        while p != self.bigT.root:
            self.calculateBitKeyForInternalNodeFromChildren(p)
            p = p.parent


        if 0:
            if qWithBiggestSubTree.up == False:
                for n in self.bigT.nodes:
                    if hasattr(n, 'br') and hasattr(n.br, 'bitKey'):
                        print "node %2i  %s" % (n.nodeNum, self.getTaxBitsString(n.br.bitKey))
                    else:
                        print "node %2i  None" % n.nodeNum
                func.setTerminalColour(BIG_T_COLOUR)
                self.bigT.draw()
                func.unsetTerminalColour()
                print "The above is bigT, after adding %i nodes" % nLeavesInBiggestSubTree
            
        #if qWithBiggestSubTree.up == False:
        #    self.checkBitKeys(self.bigT)

        if self.verbose >= 2:
            func.setTerminalColour(BIG_T_COLOUR)
            self.bigT.draw()
            func.unsetTerminalColour()
            print "The tree above is the current bigT, after adding the sub-tree."

        # Need the taxBits from the subTree to add to self.bigT
        self.bigT.taxBits = self.bigT.taxBits | subTreeTaxBits

        #if self.doCountTries:
        #    self.count_txAdded = nLeavesInBiggestSubTree
        
        if self.doPips:
            pipsWrittenNow = 0
            for pip in range(nLeavesInBiggestSubTree):
                sys.stdout.write(".")
                pipsWrittenNow +=1
                if (oldBigTNTax + pipsWrittenNow) % 50 == 0:
                    print " %4i" % (oldBigTNTax + pipsWrittenNow)
            sys.stdout.flush()
        #print "adding %i leaves" % nLeavesInBiggestSubTree
        return subTreeTaxNames  # to tryToAddAnyX
    


    
    def getQuartetsForXForNode(self, x, theNode):

        # Our goal here is to identify k as the subtree on theNode on
        # bigT were x goes.  Since bigT is fully bifurcating, we can
        # identify k by one of the numbers 0, 1, or 2, for one of the
        # 3 subtrees.  Number 0 is the subtree on theNode.leftChild.
        # Number 1 is the subtree on theNode.leftChild.sibling.  For
        # usual internal nodes, number 2 would be the subtree below
        # theNode, but if theNode is the root of bigT then number 2
        # would be the subtree on theNode.sibling.sibling.


        if self.slowCheck:
            self.checkBitKeys(self.bigT)
 
        #if self.verbose >= 1:
        #    print "        getQuartetsForXForNode()"

        # First find all the 'goodTrees' in the input trees.  It uses
        # 'theNode', which is an internal node in bigT.  A good tree
        # is one that contains at least one taxon from each of the 3
        # sub-trees attached to theNode, plus x.

        if theNode == self.bigT.root:
            bk0 = theNode.leftChild.br.bitKey
            bk1 = theNode.leftChild.sibling.br.bitKey
            bk2 = theNode.leftChild.sibling.sibling.br.bitKey
            goodTrees = [t for t in self.trees if (
                t.taxBits & x.bitKey and \
                t.taxBits & bk0) and \
                (t.taxBits & bk1) and \
                (t.taxBits & bk2)]
        else:
            bk0 = theNode.leftChild.br.bitKey
            bk1 = theNode.leftChild.sibling.br.bitKey
            bk2 = (self.allOnes ^ theNode.br.bitKey) & self.bigT.taxBits
            goodTrees = [t for t in self.trees if (
                t.taxBits & x.bitKey and \
                t.taxBits & bk0) and \
                (t.taxBits & bk1) and \
                (t.taxBits & bk2)]

        if not goodTrees:
            return
        if goodTrees:
            if self.verbose >= 2:
                print "There are %i good input trees, so will attempt to find a quartet with x=%s" % (
                    len(goodTrees), x.name)

                if 0:
                    print "  The following are the %i input trees from which to get a quartet:" % len(goodTrees)
                    for t in goodTrees:
                        func.setTerminalColour(INPUT_T_COLOUR)
                        sL = t.textDrawList(showInternalNodeNames=False,
                                           addToBrLen=0.0, width=50,
                                           autoIncreaseWidth=True,
                                           showNodeNums=True, partNum=None, model=None)
                        for s in sL:
                            print "    ",
                            print s
                        func.unsetTerminalColour()
                #if self.verbose >= 3:
                #    self.pause()

                func.setTerminalColour(BIG_T_COLOUR)
                self.bigT.draw()
                func.unsetTerminalColour()
                print "(The tree above is bigT)"
                print "getQuartetsForXForNode()  theNode=%i, x=%s" % (theNode.nodeNum, x.name)
                print "theNode.pc0 = %s" % theNode.pc0
                print "theNode.pc1 = %s" % theNode.pc1
                print "theNode.pc2 = %s" % theNode.pc2
                print "bk0 = %s" % self.getTaxBitsString(bk0)
                print "bk1 = %s" % self.getTaxBitsString(bk1)
                print "bk2 = %s" % self.getTaxBitsString(bk2)
                
        #sys.exit()
        kk = []                # k's (quartets) from all the trees.  At most one per tree
        allBits = x.bitKey | bk0 | bk1 | bk2
        if 0: #self.dbug:
            if x.name == 'W':
                func.setTerminalColour(BIG_T_COLOUR)
                self.bigT.draw()
                func.unsetTerminalColour()
                print "(The tree above is bigT)"
                print "getQuartetsForXForNode()  theNode=%i, x=%s" % (theNode.nodeNum, x.name)
                print "bk0 = %s" % self.getTaxBitsString(bk0)
                print "bk1 = %s" % self.getTaxBitsString(bk1)
                print "bk2 = %s" % self.getTaxBitsString(bk2)
                print 
                print "allBits = %s" % self.getTaxBitsString(allBits)
                print "goodTree[0] is ", 
                goodTrees[0].write()

        for gT in goodTrees:
            if self.verbose >= 2:
                func.setTerminalColour(INPUT_T_COLOUR)
                sL = gT.textDrawList(showInternalNodeNames=False,
                                   addToBrLen=0.0, width=50,
                                   autoIncreaseWidth=True,
                                   showNodeNums=True, partNum=None, model=None)
                for s in sL:
                    print "    ",
                    print s
                func.unsetTerminalColour()
                print "Input tree containing the right bitKeys"
                
            gotIt = False
            badBitsBelow = False
            #mBits = 0
            tNode = gT.root
            while not gotIt:
                biggestNBits = 0
                n = tNode.leftChild
                while n:
                    #print "while n.  Doing node %i" % n.nodeNum 
                    if n.isLeaf:
                        pass
                    elif allBits & n.br.bitKey: # at least one bit
                        bitHits = []
                        nBits = 0
                        for bk in [bk0, bk1, bk2, x.bitKey]:
                            if bk & n.br.bitKey:
                                nBits += 1
                                bitHits.append(bk)
                                if nBits > 2:
                                    #print "nBits is more than 2, so its further up this subtree."
                                    break
                        if nBits > biggestNBits:
                            biggestNBits = nBits
                        if nBits == 2:

                            # Now we have n, and n.br is the split, and there are
                            # 2 bits above it.  Presumably there are 2 bits below
                            # it as well, but it might be ambiguous, in which case
                            # it is not 2 bits.  So it is worth checking.  Maybe
                            # this should be optional, with a switch.
                            badBitsBelow = False
                            bitHitsBelow = []
                            bitsBelow = (self.allOnes ^ n.br.bitKey) & gT.taxBits
                            if 0:
                                func.setTerminalColour('violet')
                                gT.draw()
                                print "taxa below node %i are %s" % (n.nodeNum, self.getTaxBitsString(bitsBelow))
                                func.unsetTerminalColour()

                            for bk in [bk0, bk1, bk2, x.bitKey]:
                                if bk & bitsBelow:
                                    bitHitsBelow.append(bk)
                            if len(bitHitsBelow) == 2:
                                gotIt = True
                                break
                            else:
                                # ambiguous!
                                #func.writeInColour("Got %i bitHitsBelow\n" % len(bitHitsBelow), "CYAN")
                                badBitsBelow = True
                                

                        elif nBits > 2: # Its further up that subtree
                            tNode = n
                            nBits = 0
                            break
                    n = n.sibling  # Try the next subtree

                if biggestNBits < 2: # useless to go further
                    break
                if gotIt:
                    break
                if badBitsBelow:
                    break
            if not gotIt:
                if self.verbose >= 2:
                    print "The tree above does not have a relevant quartet."
                continue

            else:
                assert len(bitHits) == 2   # 2 of [bk0, bk1, bk2, x.bitKey]
                assert len(bitHitsBelow) == 2   # 2 of [bk0, bk1, bk2, x.bitKey]
                    



                if self.verbose >= 2:
                    print "The tree above has a relevant quartet, at node num %i" % n.nodeNum
                #if self.dbug and x.name == 'W':
                #    print "The tree above has a relevant quartet, at node num %i" % n.nodeNum
                k = None
                xAndKAreAbove = None

                # If x.bitKey is in bitHits, then it must be the 2nd of 2.
                if bitHits[1] == x.bitKey:
                    # Then bitHits[0] is from the k-subtree
                    if bitHits[0] == bk0:
                        k = 0
                    elif bitHits[0] == bk1:
                        k = 1
                    elif bitHits[0] == bk2:
                        k = 2

                    # x is in bitHits, so must be above n, and so must k
                    xAndKAreAbove = True
                else:
                    # x.bitKey goes with the 3rd one, the one that is
                    # not in bitHits, so it will be the k-subtree
                    allThree = [bk0, bk1, bk2]
                    allThree.remove(bitHits[0])
                    allThree.remove(bitHits[1])
                    thirdOne = allThree[0]
                    if thirdOne == bk0:
                        k = 0
                    elif thirdOne == bk1:
                        k = 1
                    elif thirdOne == bk2:
                        k = 2 

                    xAndKAreAbove = False
                assert k != None
                assert xAndKAreAbove != None
                if self.verbose >= 2:
                    print "Got k = %i, and xAndKAreAbove=%s, which means that x=%s and one of" % (
                        k, xAndKAreAbove, x.name)
                    if k == 0:
                        print "%s" % self.getTaxBitsString(bk0)
                    elif k == 1:
                        print "%s" % self.getTaxBitsString(bk1)
                    elif k == 2:
                        print "%s" % self.getTaxBitsString(bk2)
                    if xAndKAreAbove:
                        print "are above the split at node %i" % n.nodeNum
                    else:
                        print "are below the split at node %i" % n.nodeNum
                    
                    
                if self.doAddXSubTree:
                    theQuart = QJQuartet(gT, n, k, xAndKAreAbove)
                else:
                    theQuart = k
                kk.append(theQuart)

                if self.maxQuartetsPerN and len(kk) >= self.maxQuartetsPerN:
                    return kk
        return kk
                    
            
       
    def getKNodes(self, theNode):
        """Find and return the kIntNodes and kBrs in the k-subtree."""

        # We have a list kk, which may be a list of Tax objects, or a
        # list of QJQuartet objects, or ints.  So vote.
        self.votes = [0] * 3  # assume only 3 subtrees on theNode.
        for item in self.kk:
            if self.doAddXSubTree: # item is a QJQuartet object
                self.votes[item.k] += 1

            else:                # item is an int
                self.votes[item] += 1

        if self.verbose >= 2:
            print "  getKNodes(). self.votes = %s" % self.votes
        
        # In the following, we take into account the possibility of there being a tie score.
        maxVotes = max(self.votes)   
        assert maxVotes

        # It may be a tie score.  It happens a fair amount.
        if 0:
            tieScore = False
            if self.votes.count(maxVotes) > 1:
                tieScore = True
            if tieScore:
                print "tieScore is %s" % tieScore

        # Now that we have voted, get the kIntNodes and kBrs.
        kIntNodes = []
        kBrs = []

        for theVoteIndx in [0, 1, 2]:
            if self.votes[theVoteIndx] == maxVotes:
                if theNode == self.bigT.root:
                    if theVoteIndx == 0:
                        ch = theNode.leftChild
                    elif theVoteIndx == 1:
                        ch = theNode.leftChild.sibling
                    elif theVoteIndx == 2:
                        ch = theNode.leftChild.sibling.sibling
                    for n in ch.iterPreOrder():
                        if not n.isLeaf:
                            kIntNodes.append(n)
                        kBrs.append(n)
                else: # theNode is not root.
                    if theVoteIndx < 2:
                        if theVoteIndx == 0:
                            ch = theNode.leftChild
                        elif theVoteIndx == 1:
                            ch = theNode.leftChild.sibling
                        for n in ch.iterPreOrder():
                            if not n.isLeaf:
                                kIntNodes.append(n)
                            kBrs.append(n)
                    else:                                   # theVoteIndx is 2, so its below
                        for n in theNode.iterDown():
                            if not n.isLeaf:
                                kIntNodes.append(n)
                            if n.br:
                                kBrs.append(n)
                        kIntNodes = kIntNodes[1:] 

        #print "getKNodes() got kIntNodes=%s, kBrs=%s" % (
        #                [n.nodeNum for n in kIntNodes], [n.nodeNum for n in kBrs])
        #assert kIntNodes
        assert kBrs
        return (kIntNodes, kBrs)


    def getQuartetForFourTaxaInTrees(self, taxa, theTrees):

        assert len(taxa) == 4
        votes = [0] * 3
        # **..  -- 0
        # *.*.  -- 1
        # *..*  -- 2

        bk0 = taxa[0].bitKey
        bk1 = taxa[1].bitKey
        bk2 = taxa[2].bitKey
        bk3 = taxa[3].bitKey
        allBits = bk0 | bk1 | bk2 | bk3
        
        for gT in theTrees:

            if self.verbose >= 2:
                func.setTerminalColour(INPUT_T_COLOUR)
                sL = gT.textDrawList(showInternalNodeNames=False,
                                   addToBrLen=0.0, width=50,
                                   autoIncreaseWidth=True,
                                   showNodeNums=True, partNum=None, model=None)
                for s in sL:
                    print "    ",
                    print s
                func.unsetTerminalColour()
                #print "Input tree containing the right bitKeys"
                
            gotIt = False
            badBitsBelow = False
            #mBits = 0
            tNode = gT.root
            while not gotIt:
                biggestNBits = 0
                n = tNode.leftChild
                while n:
                    #print "while n.  Doing node %i" % n.nodeNum 
                    if n.isLeaf:
                        pass
                    elif allBits & n.br.bitKey: # at least one bit
                        bitHits = []
                        nBits = 0
                        for bk in [bk0, bk1, bk2, bk3]:
                            if bk & n.br.bitKey:
                                nBits += 1
                                bitHits.append(bk)
                                if nBits > 2:
                                    #print "nBits is more than 2, so its further up this subtree."
                                    break
                        if nBits > biggestNBits:
                            biggestNBits = nBits
                        if nBits == 2:

                            # Now we have n, and n.br is the split, and there are
                            # 2 bits above it.  Presumably there are 2 bits below
                            # it as well, but it might be ambiguous, in which case
                            # it is not 2 bits.  So it is worth checking.  Maybe
                            # this should be optional, with a switch.
                            badBitsBelow = False
                            bitHitsBelow = []
                            bitsBelow = (self.allOnes ^ n.br.bitKey) & gT.taxBits
                            if 0:
                                func.setTerminalColour('violet')
                                gT.draw()
                                print "taxa below node %i are %s" % (n.nodeNum, self.getTaxBitsString(bitsBelow))
                                func.unsetTerminalColour()

                            for bk in [bk0, bk1, bk2, bk3]:
                                if bk & bitsBelow:
                                    bitHitsBelow.append(bk)
                            if len(bitHitsBelow) == 2:
                                gotIt = True
                                break
                            else:
                                # ambiguous!
                                #func.writeInColour("Got %i bitHitsBelow\n" % len(bitHitsBelow), "CYAN")
                                badBitsBelow = True
                                

                        elif nBits > 2: # Its further up that subtree
                            tNode = n
                            nBits = 0
                            break
                    n = n.sibling  # Try the next subtree

                if biggestNBits < 2: # useless to go further
                    break
                if gotIt:
                    break
                if badBitsBelow:
                    break
            if not gotIt:
                if self.verbose >= 2:
                    print "The tree above does not have a quartet with those four taxa."
                continue

            else:
                assert len(bitHits) == 2   # 2 of [bk0, bk1, bk2, bk3]
                assert len(bitHitsBelow) == 2   # 2 of [bk0, bk1, bk2, bk3]
                    



                if self.verbose >= 2:
                    print "The tree above has a quartet, at node num %i" % n.nodeNum

                # **..  -- 0
                # *.*.  -- 1
                # *..*  -- 2

                # bitHits are lists of bitKeys
                if bk0 in bitHits:
                    if bk1 in bitHits:
                        assert bk0 in bitHits and bk1 in bitHits
                        assert bk2 in bitHitsBelow and bk3 in bitHitsBelow
                        votes[0] += 1
                    elif bk2 in bitHits:
                        assert bk0 in bitHits and bk2 in bitHits
                        assert bk1 in bitHitsBelow and bk3 in bitHitsBelow
                        votes[1] += 1
                    else:
                        assert bk0 in bitHits and bk3 in bitHits
                        assert bk1 in bitHitsBelow and bk2 in bitHitsBelow
                        votes[2] += 1
                            
                            
                else:
                    if bk1 in bitHits:
                        if bk2 in bitHits:
                            assert bk1 in bitHits and bk2 in bitHits
                            assert bk0 in bitHitsBelow and bk3 in bitHitsBelow
                            votes[2] += 1
                        else:
                            assert bk1 in bitHits and bk3 in bitHits
                            assert bk0 in bitHitsBelow and bk2 in bitHitsBelow
                            votes[1] += 1
                    else:
                        assert bk2 in bitHits and bk3 in bitHits
                        assert bk0 in bitHitsBelow and bk1 in bitHitsBelow
                        votes[0] += 1
                        

        theMax = max(votes)
        if not theMax:
            return None
        nMaxs = votes.count(theMax)
        if nMaxs > 1:  # we don't want ambiguity
            return None
        myMaxIndex = votes.index(theMax)
        #print "taxa names are %s, " % [tx.name for tx in taxa],
        #print "myMaxIndex is", myMaxIndex
        
        savedWarnReadNoFile = var.warnReadNoFile
        var.warnReadNoFile = 0
        read('(A, B, (C, D));')
        var.warnReadNoFile = savedWarnReadNoFile
        t = var.trees.pop()
        if myMaxIndex == 0: # 0 and 1 go together, so 2 and 3 go together
            t.nodes[1].name = taxa[2].name
            t.nodes[2].name = taxa[3].name
            t.nodes[4].name = taxa[0].name
            t.nodes[5].name = taxa[1].name
        elif myMaxIndex == 1: # 0 and 2 go together, so 1 and 3 go together
            t.nodes[1].name = taxa[1].name
            t.nodes[2].name = taxa[3].name
            t.nodes[4].name = taxa[0].name
            t.nodes[5].name = taxa[2].name
        elif myMaxIndex == 2: # 0 and 3 go together, so 1 and 2 go together
            t.nodes[1].name = taxa[0].name
            t.nodes[2].name = taxa[3].name
            t.nodes[4].name = taxa[1].name
            t.nodes[5].name = taxa[2].name

        return t
        

    # *************************************************
    # *************************************************
    # *************************************************
    # *************************************************


    def pause(self):
        if self.verbose >= 3:
            func.setTerminalColour('blue')
            raw_input("Hit return to continue ...")
            print "=" * 90
            func.unsetTerminalColour()
        else:
            pass

    def popcount(self, n):
        count = 0
        for tx in self.taxa:
            if n < tx.bitKey:
                return count
            if tx.bitKey & n:
                count += 1
        return count



    def setTreeTaxBits(self, theTree):
        """This needs to have bitKey attributes for nodes already set.

        taxBits says which of self.taxa are in theTree.
        """

        theTree.taxBits = 0L
        for ch in theTree.root.iterChildren():
            #print ch.br.bitKey
            theTree.taxBits = theTree.taxBits | ch.br.bitKey # bitwise 'or'


    def decorateTreeWithBitKeys(self, theTree):
        """Put bitKey attributes on all branches of the tree.

        bitKey's are attributes of node branches, and are binary
        numbers that say what taxa are above that branch in the tree.
        """
        
        for n in theTree.iterPostOrder():
            if n != theTree.root:
                if n.isLeaf:
                    # order comes from self.taxNames, not from the tree
                    n.br.bitKey = 1L << self.taxNames.index(n.name)
                else:
                    childrenNums = theTree.getChildrenNums(n)
                    x = theTree.nodes[childrenNums[0]].br.bitKey
                    for i in childrenNums[1:]:
                        y = theTree.nodes[i].br.bitKey
                        x = x | y
                    n.br.bitKey = x

        if 0:
            theTree.draw()
            print "bitKeys for the tree above:"
            for n in theTree.iterInternalsNoRoot():
                print " %2i  %s" % (n.nodeNum, self.getTaxBitsString(n.br.bitKey))
            self.setTreeTaxBits(theTree)
            print "taxBits is %s" % self.getTaxBitsString(theTree.taxBits)
            print "order:",
            tNames = [tx.name for tx in self.taxa]
            for tName in tNames:
                print " %s" % tName,
            print
    

    def checkBitKeys(self, theTree):
        for n in theTree.iterPostOrder():
            if n != theTree.root:
                if n.isLeaf:
                    # order comes from self.taxNames, not from the tree
                    theBitKey = 1L << self.taxNames.index(n.name)
                    # Hard to imagine how this would get messed up.
                    if n.br.bitKey != theBitKey: #, "checkBitKeys().  node %i is wrong" % n.nodeNum
                        print
                        theTree.draw()
                        gm = ["checkBitKeys()"]
                        gm.append("leaf node %i bitKey is %s, should be %s" % (n.nodeNum, n.br.bitKey, theBitKey))
                        raise Glitch, gm
                else:
                    childrenNums = theTree.getChildrenNums(n)
                    x = theTree.nodes[childrenNums[0]].br.bitKey
                    for i in childrenNums[1:]:
                        y = theTree.nodes[i].br.bitKey
                        x = x | y
                    if n.br.bitKey != x:
                        print
                        theTree.draw()
                        gm = ['checkBitKeys()']
                        for ch in n.iterChildren():
                            gm.append("node %i bitKey is %s" % (ch.nodeNum, ch.br.bitKey))
                        gm.append("but node %i bitKey is %s -- wrong!" % (n.nodeNum, n.br.bitKey))
                        gm.append("should be %s" % x)
                        raise Glitch, gm
        # Check theTree.taxBits
        checkTaxBits = 0L
        for ch in theTree.root.iterChildren():
            #print ch.br.bitKey
            checkTaxBits = checkTaxBits | ch.br.bitKey # bitwise 'or'
        if checkTaxBits != theTree.taxBits:
            gm = ['checkBitKeys()  tcheck tree taxBits']
            gm.append('theTree.taxBits=%s, check=%s' % (theTree.taxBits, checkTaxBits))
            raise Glitch, gm

                        
    def checkBigTInternalNodes(self):
        recalculatedInternalNodesSet = set([n for n in self.bigT.iterInternalsNoRoot()])
        recalculatedInternalNodesSet.add(self.bigT.root)
        existingSet = set(self.bigTInternalNodes)
        if recalculatedInternalNodesSet != existingSet:
            func.setTerminalColour(BIG_T_COLOUR)
            self.bigT.draw()
            func.unsetTerminalColour()
            gm = ["checkBigTInternalNodes()"]
            gm.append("existingSet = %s" % [n.nodeNum for n in existingSet])
            gm.append("recalculated = %s" % [n.nodeNum for n in recalculatedInternalNodesSet])
            raise Glitch, gm
        
                            
        

    def calculateBitKeyForInternalNodeFromChildren(self, theNode):
        children = [ch for ch in theNode.iterChildren()]
        x = children[0].br.bitKey
        for ch in children[1:]:
            y = ch.br.bitKey
            x = x | y
        theNode.br.bitKey = x
        
        
        
    def getTaxBitsString(self, taxBits):
        sTax = []
        for tx in self.taxa:
            if (taxBits & tx.bitKey) == tx.bitKey:
                sTax.append(tx.name)
            else:
                #sTax.append('.')
                pass
        #return ' '.join(sTax) + " (small bits on the left)"
        return ' '.join(sTax)

    def locateAndCountLeavesInXSubTreeFromQuartet(self, q, x):
        """Find the node in the q.tree with the x-subtree, and count leaves.

        We also find the xNode, the node in q.tree at which the
        x-subtree can be taken, and 'up', which says whether that
        x-subtree goes up or down from that xNode.  Both of these
        become attributes of q, the QJQuartet object.  The leaf count
        is returned.
        """

        gm = ['locateAndCountLeavesInXSubTreeFromQuartet()']
        #func.printColour("Now in %s" % gm[0], 'BLUE')
        #q.dump()                # includes a tree drawing

        # Arg q is a QJQuartet object, and q.tree is an input tree
        # that has a relevant quartet ij|kx.  Here x is the leaf that
        # we want to add, and q.k is the number (one of 0, 1, or 2)
        # indicating which subtree in bigT that contains one or more
        # taxa in bigT that are on the x-side of the relevant quartet.
        # Sorry! - k is being used in two senses here -- in the
        # quartet ij|kx k is a taxon, but q.k is number indicating a
        # bigT-subtree that contains a k-taxon (or perhaps more than
        # one k-taxa) shared by bigT and q.tree.  But at this point
        # the k-taxon (or k-taxa) is not pinpointed, not identified,
        # and we don't really ever need to actually identify it as a
        # taxon; we just know that an identified subtree (the
        # k-subtree) in bigT and the x-side of the split at q.node in
        # q.tree share taxa (exclusively, ie there are no bigT-k-taxa
        # on the non-x side of the split in q.tree, and no
        # bigT-non-k-taxa on the x-side of the split in q.tree).

        # So we know that q.tree at q.node, either up or down, has a
        # k-taxon and x.  The k-taxon (or k-taxa) is in bigT, and we
        # do not want it in the x-subtree that we add, so we start
        # with the split at q.node and work toward x in q.tree until
        # no bigT taxa are in that smaller x-subtree.  The smaller
        # x-subtree might be just x by itself, but if we are lucky it
        # might have other taxa that are not yet in bigT.  That
        # smaller x-subtree may be added to bigT if the branch on the
        # k-subtree in bigT can be identified, but that is another
        # story.


        if q.xAndKAreAbove:
            theNode = q.node
            #if theNode == q.tree.root:
            #    for n in theNode.iterChildren():
            #        if n.br.bitKey & x.bitKey:
            #            theNode = n
            #            break
            assert theNode != q.tree.root
            #print "theNode %i, and x (%s) is above." % (theNode.nodeNum, x.name)
            if theNode.br.bitKey & self.bigT.taxBits:
                bigTBitsAreOnTheSameBranch = True
            else:
                bigTBitsAreOnTheSameBranch = False
            #print "bigTBitsAreOnTheSameBranch is %s" % bigTBitsAreOnTheSameBranch
            while bigTBitsAreOnTheSameBranch:
                if theNode.isLeaf:
                    bigTBitsAreOnTheSameBranch = False
                    break
                if not bigTBitsAreOnTheSameBranch:
                    break
                for n in theNode.iterChildren():
                    #print "   ** Looking at node %i" % n.nodeNum
                    if n.br.bitKey & x.bitKey:
                        #print "     ** x is above"
                        #print "     ** bigT bits above is %s" % (n.br.bitKey & self.bigT.taxBits)
                        if not (n.br.bitKey & self.bigT.taxBits):
                            bigTBitsAreOnTheSameBranch = False
                        theNode = n
                        break
            assert not bigTBitsAreOnTheSameBranch
            assert theNode.br.bitKey & x.bitKey
            #print "theNode %i, and x is above, and there are no bigT bits above" % theNode.nodeNum
            #sys.exit()

            # Now we want to make sure that the x-subtree that we have
            # found is fully bifurcating.  If it is not, then follow
            # nodes toward x until it is fully bifurcating.  If we are
            # very unlucky, we will need to follow all the way to x,
            # and the subtree will be a single leaf.

            while 1:
                if theNode.isLeaf:
                    break
                if q.tree.subTreeIsFullyBifurcating(theNode, up=True):
                    break
                else:
                    for n in theNode.iterChildren():
                        if n.br.bitKey & x.bitKey:
                            theNode = n
                            break
            assert q.tree.subTreeIsFullyBifurcating(theNode, up=True)

            #print "the x subtree (fully bifurcating, and no bigT taxa) is everything from node %i up." % (
            #    theNode.nodeNum)
            #sys.exit()
            q.xNode = theNode
            q.up = True
            if theNode.isLeaf:
                return 1
            else:
                return self.popcount(theNode.br.bitKey)
 
        else: # xAndK are below
            theNode = q.node
            xIsBelow = True
            bigTBitsAreBelow = True
            #print 'bigT.taxBits are %s' % self.getTaxBitsString(self.bigT.taxBits) 
            #print 'theNode %i, theNode.br.bitKey is %s' % (
            #    theNode.nodeNum, self.getTaxBitsString(theNode.br.bitKey))
            bitsBelow = (self.allOnes ^ theNode.br.bitKey) & q.tree.taxBits
            #print '    taxa below are %s' % (self.getTaxBitsString(bitsBelow))
            assert x.bitKey & bitsBelow
            assert self.bigT.taxBits & bitsBelow
            if 0:
                for qn in q.tree.iterLeavesNoRoot():
                    if qn.name in self.bigT.taxNames:
                        qn.oldName = qn.name
                        qn.name = "\033[1;36m%s ***\033[m" % qn.name
                    elif qn.name == x.name:
                        qn.oldName = qn.name
                        qn.name = "\033[1;31m%s\033[m" % qn.name
                q.tree.draw()
                for qn in q.tree.iterLeavesNoRoot():
                    if hasattr(qn, 'oldName'):
                        qn.name = qn.oldName
                        del(qn.oldName)
                print "x (%s) and at least some bigT taxa (%s)" % (x.name, self.getTaxBitsString(self.bigT.taxBits))
                print "   are both somewhere down from node %i on the QJQuartet tree" %  q.node.nodeNum

            while 1:
                theNode = theNode.parent
                #print "while 1: theNode=%i" % theNode.nodeNum
                if theNode == q.tree.root:
                    xIsBelow = False
                    #print "x (%s) is no longer below" % x.name
                    break
                bitsBelow = (self.allOnes ^ theNode.br.bitKey) & q.tree.taxBits
                #print "theNode %i, bitsBelow = %s" % (theNode.nodeNum, self.getTaxBitsString(bitsBelow))
                if not (x.bitKey & bitsBelow):
                    xIsBelow = False
                    #print "x (%s) is no longer below" % x.name
                if not (self.bigT.taxBits & bitsBelow):
                    bigTBitsAreBelow = False
                    #print "there are no longer any bigT bits below."
                if xIsBelow and bigTBitsAreBelow:
                    pass
                else:
                    break
            #print "*** theNode=%i, xIsBelow is %s, bigTBitsAreBelow = %s" % (
            #    theNode.nodeNum, xIsBelow, bigTBitsAreBelow)
            #sys.exit()
            
            # So we have now found the x-subtree.  It is off of theNode, either above it or below it.  

            # Now we want to make sure that the x-subtree that we have
            # found is fully bifurcating.  If it is not, then follow
            # nodes toward x until it is fully bifurcating.  If we are
            # very unlucky, we will need to follow all the way to x,
            # and the subtree will be a single leaf.

            if xIsBelow:
                while 1:
                    if not xIsBelow:
                        break
                    if q.tree.subTreeIsFullyBifurcating(theNode, up=False):
                        #print "theSubTree going down from theNode %i is fully bifurcating" % theNode.nodeNum
                        break
                    else:
                        #print "theSubTree going down from theNode %i is not fully bifurcating" % theNode.nodeNum
                        previousNode = theNode
                        theNode = theNode.parent
                        for n in theNode.iterChildren():
                            if n == previousNode:
                                pass
                            elif n.br.bitKey & x.bitKey:
                                xIsBelow = False
                                break
                if xIsBelow:
                    assert q.tree.subTreeIsFullyBifurcating(theNode, up=False)
                    


            #print "locateAndCountLeavesInXSubTreeFromQuartet()  x=%s, k=%s.  node with split is %i" % (
            #                                                   x.name, q.k.name, q.node.nodeNum)
            #print "theNode is now %i" % theNode.nodeNum
            #print "xIsBelow=%s, bigTBitsAreBelow=%s" % (xIsBelow, bigTBitsAreBelow)
            #if not xIsBelow:
            #    print "(However, the bigTBitsAreBelow value does not matter since x is not below)."
            #q.tree.draw()
            #print "Drawing above is q.tree.  theNode is %i" % theNode.nodeNum
            #print "bigT taxNames is %s" % self.bigT.taxNames
            #sys.exit()
            if xIsBelow and not bigTBitsAreBelow:
                #print "x is everything below node %i (now isolated from bigT taxa)." % theNode.nodeNum
                flippedBits = q.tree.taxBits  ^ theNode.br.bitKey
                #print "flippedBits = %s" % self.getTaxBitsString(flippedBits)
                assert not flippedBits & self.bigT.taxBits
                q.xNode = theNode
                q.up = False
                return self.popcount(flippedBits)
            elif not xIsBelow: # At this point, whether bigTBitsAreBelow or not does not matter.
                if theNode == q.tree.root:
                    for n in theNode.iterChildren():
                        if n.br.bitKey & x.bitKey:
                            theNode = n
                            break
                #print "theNode %i, and x is above." % theNode.nodeNum

                if theNode.br.bitKey & self.bigT.taxBits:
                    bigTBitsAreOnTheSameBranch = True
                else:
                    bigTBitsAreOnTheSameBranch = False
                #print "bigTBitsAreOnTheSameBranch is %s" % bigTBitsAreOnTheSameBranch
                while bigTBitsAreOnTheSameBranch:
                    if theNode.isLeaf:
                        bigTBitsAreOnTheSameBranch = False
                        break
                    if not bigTBitsAreOnTheSameBranch:
                        break
                    for n in theNode.iterChildren():
                        #print "   ** Looking at node %i" % n.nodeNum
                        if n.br.bitKey & x.bitKey:
                            #print "     ** x is above"
                            #print "     ** bigT bits above is %s" % (n.br.bitKey & self.bigT.taxBits)
                            if not (n.br.bitKey & self.bigT.taxBits):
                                bigTBitsAreOnTheSameBranch = False
                            theNode = n
                            break
                assert not bigTBitsAreOnTheSameBranch
                assert theNode.br.bitKey & x.bitKey
                #print "Found theNode %i, and x is above, and there are no bigT bits above" % theNode.nodeNum
                while 1:
                    if theNode.isLeaf:
                        break
                    if q.tree.subTreeIsFullyBifurcating(theNode, up=True):
                        break
                    else:
                        for n in theNode.iterChildren():
                            if n.br.bitKey & x.bitKey:
                                theNode = n
                                break
                assert q.tree.subTreeIsFullyBifurcating(theNode, up=True)
                
                #print "the x subtree (fully bifurcating, and no bigT taxa) is everything from node %i up." % (
                #    theNode.nodeNum)
                q.xNode = theNode
                q.up = True
                if theNode.isLeaf:
                    return 1
                else:
                    return self.popcount(theNode.br.bitKey)
 

    def _qj_extractFullyBifurcatingTree(self, aTree):
        """Dupe the tree, and remove leaves until it is fully bifurcating.

        Tries to make the biggest possible fully bifurcating tree.
        """

        theTree = aTree.dupe()

        PP = []
        for n in theTree.iterNodes():
            if n.isLeaf:
                pass
            elif n == theTree.root:
                nCh = n.getNChildren()
                if nCh > 3:
                    PP.append(n)
            else:
                nCh = n.getNChildren()
                if nCh > 2:
                    PP.append(n)
        #if PP:
        #    print "There are %i polytomies" % len(PP)
        if not PP:
            return theTree

        while PP:
            pN = PP.pop()
            if pN in theTree.nodes:
                #print "==================== %i ===========" % pN.nodeNum
                chNN = []
                sumOfAboves = 0
                for n in pN.iterChildren():
                    if n.isLeaf:
                        above = 1
                        sumOfAboves += 1
                    else:
                        above = self.popcount(n.br.bitKey)
                        sumOfAboves += above
                    chNN.append([n, above, 1])
                if pN != theTree.root:
                    chNN.append([pN.parent, theTree.nTax - sumOfAboves, 0])
                chNN = func.sortListOfListsOnListElementNumber(chNN, 1)
                chNN.reverse()

                #for chN in chNN:
                #    print "%2i (%i %i) " % (chN[0].nodeNum, chN[1], chN[2]) 


                if pN == theTree.root:
                    theNChildren = 3
                else:
                    theNChildren = 2

                while 1:
                    if pN.getNChildren() <= theNChildren:
                        break
                    theAction = chNN.pop()
                    if theAction[2]: # going up -- easy
                        pass
                    else:
                        # Going down -- not as easy.
                        # Does this ever happen? -- yes indeed.
                        #print "_qj_extractFullyBifurcatingTree() ============== ReRooting."
                        theTree.reRoot(pN)
                        theNChildren = 3
                    rNode = theAction[0]
                    # removeNode() renumbers the nodes and resets pre and post order.
                    theTree.removeNode(rNode)
                    #print "The tree following is after removal of the node, and re-numbering."
                    #theTree.draw()
                    self.decorateTreeWithBitKeys(theTree)
                    self.setTreeTaxBits(theTree)

        return theTree


    def allChooseFour(self, things):
        ttNN = []
        nTx = len(things)
        nTxMinus1 = nTx - 1
        nTxMinus2 = nTx - 2
        nTxMinus3 = nTx - 3
        w = 0
        x = 1
        y = 2
        z = 3
        counter = 0
        #print "%4i | %4i %4i %4i %4i" % (counter, w,x,y,z)
        tN = [things[w],things[x],things[y],things[z]]
        ttNN.append(tN)
        counter += 1
        while 1:
            z += 1
            if z >= nTx:
                y += 1
                if y >= nTxMinus1:
                    x += 1
                    if x >= nTxMinus2:
                        w += 1
                        if w >= nTxMinus3:
                            #print "w is now %i" % w
                            break
                        x = w + 1
                    y = x + 1
                z = y + 1
            #print "%4i | %4i %4i %4i %4i" % (counter, w,x,y,z)
            tN = [things[w],things[x],things[y],things[z]]
            ttNN.append(tN)
            #thisQNum = self.qNum([w,x,y,z])
            #assert counter == thisQNum
            counter += 1
        #print "There are %i possible quartets." % counter
        return ttNN
