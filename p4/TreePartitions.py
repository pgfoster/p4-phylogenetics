import types,string,cStringIO,sys,os
from Tree import Tree
from Node import Node,NodePart,NodeBranchPart
from Trees import Trees
from Nexus import Nexus,NexusData
from Glitch import Glitch
#from NexusToken import nextTok,safeNextTok,nexusSkipPastNextSemiColon
import NexusToken   # needed for cStrings
#import NexusToken2  # Faster, in c
import func
from Var import var

longMessage1 = """
This table shows, for splits that were used in the cons tree (ie not
necessarily all splits), which of those splits had the root on it in
the original bifurcating input trees.  The tree as returned from this
method is rooted on a node on the branch with the highest biRootCount
(or if there were several highest, it is rooted on a node arbitrarily
on the first of those).
"""

longMessage2 = """
This table shows, for nodes derived from splits that were used in the
cons tree (ie not necessarily all splits), which of those nodes were
roots in the original input trees.  The tree as returned from this
method is rooted on the node with the highest rootCount (or if there
are several highest, it is rooted arbitrarily on the first of
those).
"""


class SplitModelUsage(object):
    def __init__(self, TPMI):
        if not isinstance(TPMI, TreePartitionsModelInfo):
            gm = ["SplitModelUsage init."]
            gm.append("Expecting a TreePartitionsModelInfo instance.")
            raise Glitch, gm
        self.nParts = TPMI.nParts
        self.parts = []
        for i in range(self.nParts):
            SMTP = SplitModelUsagePart()
            SMTP.compCounts = [0] * TPMI.parts[i].nComps
            SMTP.rMatrixCounts = [0] * TPMI.parts[i].nRMatrices
            SMTP.gdasrvCounts = [0] * TPMI.parts[i].nGdasrvs
            self.parts.append(SMTP)

    def dump(self):
        print "                modelUsage nParts=%s" % self.nParts
        for pNum in range(self.nParts):
            print "                    part %s" % pNum
            print "%35s = %s" % ('compCounts', self.parts[pNum].compCounts)
            print "%35s = %s" % ('rMatrixCounts', self.parts[pNum].rMatrixCounts)
            print "%35s = %s" % ('gdasrvCounts', self.parts[pNum].gdasrvCounts)

    def zero(self):
        for p in self.parts:
            for i in range(len(p.compCounts)):
                p.compCounts[i] = 0
            for i in range(len(p.rMatrixCounts)):
                p.rMatrixCounts[i] = 0
            for i in range(len(p.gdasrvCounts)):
                p.gdasrvCounts[i] = 0
                


class SplitModelUsagePart(object):
    def __init__(self):
        self.compCounts = None
        self.rMatrixCounts = None
        self.gdasrvCounts = None


class Split(object):
    #def __init__(self, nParts, nModels):
    def __init__(self):
        # self.set is a 1-based set of integers representing taxa on
        # one side of the split.  So naturally the taxa have to be
        # ordered.  It is 1-based because when I make a consensus, the
        # first thing I do is make a comb tree, and the terminal node
        # numbers in that tree and the resulting cons tree will
        # correspond to these 1-based set numbers.
        self.set = None
        self.key = None
        self.string = None
        self.count = 0.0
        self.rootCount = 0.0
        # rootCount2 is needed only for non-biRoot trees that are
        # rooted on the parent of the first taxon.
        self.rootCount2 = 0.0  
        self.proportion = 0.0
        self.cumBrLen = 0.0
        self.biRootCumBrLen = 0.0
        self.modelUsage = None
        self.biRootModelUsage = None
        self.rootModelUsage = None

    def dump(self):
        print "    Split dump:"
        print "%25s = %s" % ('set', self.set)
        print "%25s = %s" % ('key', self.key)
        print "%25s = %s" % ('string', self.string)
        print "%25s = %s" % ('count', self.count)
        print "%25s = %s" % ('rootCount', self.rootCount)
        print "%25s = %s" % ('proportion', self.proportion)
        print "%25s = %s" % ('cumBrLen', self.cumBrLen)
        if self.modelUsage:
            self.modelUsage.dump()
        if self.rootModelUsage:
            print "                rootModelUsage.dump()"
            #self.rootModelUsage.dump()
            for pNum in range(self.rootModelUsage.nParts):
                print "                    part %s" % pNum
                print "%35s = %s" % ('compCounts', self.rootModelUsage.parts[pNum].compCounts)

    def combineWith(self, otherSplitObject):
        # no need to do proportion -- that is handled by _finishSplits()
        # model stuff is ignored -- actually, wiped! -- due to lazy programming.
        oso = otherSplitObject
        assert self.key == oso.key
        self.count += oso.count
        self.rootCount += oso.rootCount
        self.rootCount2 += oso.rootCount2
        self.cumBrLen += oso.cumBrLen
        self.biRootCumBrLen += oso.biRootCumBrLen
        self.modelUsage = None
        self.biRootModelUsage = None
        self.rootModelUsage = None
        



class TreePartitionsModelInfo(object):
    def __init__(self):
        self.nParts = 0
        self.parts = []  # TreePartitionsModelInfoPart objects

    def dump(self):
        print "TreePartitionsModelInfo.dump()  nParts=%s" % self.nParts
        for i in range(self.nParts):
            self.parts[i].dump()

    def check(self):
        gm = ['TreePartitionsModelInfo.check()']
        if self.nParts < 1:
            gm.append("No parts.")
            raise Glitch, gm
        isHetero = 0
        for pNum in range(self.nParts):
            p = self.parts[pNum]
            if p.partNum != pNum:
                self.dump()
                gm.append("partNum is wrong.")
                raise Glitch, gm
            if p.nComps < 1:
                gm.append("Incomplete model info.  Part %s, no nComps." % pNum)
                raise Glitch, gm
            elif p.nComps > 1:
                isHetero = 1

            if p.nRMatrices < 1:
                gm.append("Incomplete model info.  Part %s, no nRMatrices." % pNum)
                raise Glitch, gm
            elif p.nRMatrices > 1:
                isHetero = 1

            if p.nGdasrvs < 0:
                gm.append("Incomplete model info.  Part %s, nGdasrvs not specified." % pNum)
                raise Glitch, gm
            elif p.nGdasrvs > 1:
                isHetero = 1
        if isHetero:
            return 2
        else:
            return 1



class TreePartitionsModelInfoPart(object):
    def __init__(self):
        self.partNum = -1
        self.nComps = 0
        self.nRMatrices = 0
        self.nGdasrvs = -1

    def dump(self):
        #print "    TreePartitionsModelInfoPart.dump  partNum=%s" % self.partNum
        print "    partNum=%s" % self.partNum
        print "        nComps     = %s" % self.nComps
        print "        nRMatrices = %s" % self.nRMatrices
        print "        nGdasrvs   = %s" % self.nGdasrvs



def _getModelInfo(theComment):

    # Eg [&&p4 models p1 c0.2 r0.2 g0.1]
    # or [&&p4 models p2 c0.2 r0.2 g0.1 c1.1 r1.1 g1.1]
    gm = ['TreePartitions._getModelInfo()']
    #print "_getModelInfo(), from theComment = %s" % theComment
    flob = cStringIO.StringIO(theComment)
    flob.seek(1, 0) # past the [

    tok = NexusToken.nextTok(flob)
    if not tok or tok != '&&p4':
        print "a got tok=%s" % tok
        flob.close()
        return  # not an error, just the wrong kind of comment
    tok = NexusToken.nextTok(flob)
    if not tok or tok != 'models':
        print "b got tok=%s" % tok
        flob.close()
        gm.append("Expecting 'models'.  Got %s" % tok)
        raise Glitch, gm
    tok = NexusToken.nextTok(flob)
    if tok[0] != 'p':
        gm.append("Expecting 'pN'.  Got %s" % tok)
        raise Glitch, gm
    try:
        nParts = int(tok[1:])

        if nParts > 0:
            theModelInfo = TreePartitionsModelInfo()
            theModelInfo.nParts = nParts
            for i in range(nParts):
                mp = TreePartitionsModelInfoPart()
                mp.partNum = i
                theModelInfo.parts.append(mp)
        else:
            theModelInfo = None
    except ValueError:
        gm.append("Failed to parse model comment.")
        raise Glitch, gm

    tok = NexusToken.safeNextTok(flob)
    while 1:
        #print "top of loop, tok = %s" % tok
        if tok == ']':
            break
        elif tok[0] in ['c', 'r', 'g']:
            ending = tok[1:]
            splitEnding = string.split(ending, '.')
            #print "got splitEnding = %s" % splitEnding
            try:
                firstNum = int(splitEnding[0])
                secondNum = int(splitEnding[1])
            except ValueError:
                gm.append("Failed to parse numbers in model comment.")
                raise Glitch, gm
            if tok[0] == 'c':
                theModelInfo.parts[firstNum].nComps = secondNum
            elif tok[0] == 'r':
                theModelInfo.parts[firstNum].nRMatrices = secondNum
            elif tok[0] == 'g':
                theModelInfo.parts[firstNum].nGdasrvs = secondNum
        else:
            gm.append("Bad token %s")
            raise Glitch, gm
        tok = NexusToken.safeNextTok(flob)

    #theModelInfo.dump()
    return theModelInfo


class TP_TinyTax(object):
    def __init__(self, taxName):
        self.name = taxName
        self.rawSplitKey = None
        self.splitKey = None
    



class TreePartitions(object):
    """A container for tree bipartitions or splits.

    Start it up like this::

        tp = TreePartitions(inThing)

    where inThing can be a file name, a Trees object, or a Tree
    object.

    If you are reading from a file (generally a bootstrap or mcmc
    output), you can skip some trees at the beginning, and optionally
    after that read only a maximum number of trees, like this::

        tp = TreePartitions('myFile', skip=1000, max=500)

    Then you can do::

        tp.writeSplits()  # in dot-star notation, like .***..***
        t = tp.consensus()

    All the input trees need to have the same taxa.

    If your input trees are in more than one place, you can read in
    more trees with the method:: 
    
        read(whatever, skip, max)  
    
    The default setting in p4 is to use any tree weights supplied.
    This will cause consensus trees created by p4 to differ from ones
    created using Paup as the default in Paup is to not consider 
    tree weights. P4s behavior can be modified to mimic Paup 
    using the following statement before creating a TreePartitions object. :: 
    
        var.nexus_getWeightCommandComments = 0
        tp = TreePartitions(inThing)
        
    If you want to restore the default behavior issue this command::
    
        var.nexus_getWeightCommandComments = 1

**Making consensus trees**

P4 makes majority rule consensus trees with extra compatible splits.  It
is like PAUP does when you do the ``contree`` command like the following::

     contree all/strict=no majrule=yes percent=50 le50=yes useTreeWts=yes;

or what MrBayes does when you do a ``sumt contype=allcompat``.

Making cons trees uses this class, TreePartitions, which takes trees apart
into their 'splits', aka components or tree bipartitions.  So if you have
a tree file with trees from an MCMC output or from a bootstrap, you can
make a cons tree by the following::

     tp = TreePartitions('yourFile')
     t = tp.consensus()

If you want to skip some trees at the beginning of a file (often the
burn-in for an MCMC), or if you want to read in a maximum number of
trees (which might be useful for convergence testing when using an
MCMC), you do something like::

     tp = TreePartitions('yourFile', skip=100, max=200)

When you get a cons tree with the ``consensus()`` method, the support for
splits is placed in Node.br.support attributes.  This allows some
flexibility in how the support is displayed.  To see those support
values when you draw the tree to the screen, you will need to transfer
the support information to the node names, like this::

     for n in t.iterInternalsNoRoot():
         n.name = "%.2f" % n.br.support # you can specify the precision

For nice eps drawings, you might want to put the support on the tree as
the branch.name rather than the node.name, and you can do that with
something like this::

     for n in t.iterInternalsNoRoot():
         n.br.name = "%.0f" % (n.br.support * 100.) # convert to percent

    """

    def __init__(self, inThing=None, skip=0, max=None, taxNames=None):
        self.nTrees = 0
        self.taxNames = None
        self.nTax = 0        # Merely len(taxNames)
        self.splits = []
        self.splitsHash = {}
        # biSplits and biSplitsHash for the "even" side of the biRoot bifurcation.
        self.biSplits = []
        self.biSplitsHash = {}
        self.isBiRoot = None      # Is it a bifurcating root?
        self.doModelComments = 0
        self.modelInfo = None
        if inThing:
            self.read(inThing, skip, max, taxNames)
                      

    def read(self, inThing, skip=0, max=None, taxNames=None):
        """Read in a tree, or some trees, or a file of trees.

        Arg inThing can be a file name, a Tree instance, or a
        Trees instance.  Args skip and max only come into play if
        inThing is a file."""

        gm = ['TreePartitions.read()']

        if not inThing:
            gm.append("No input?")
            raise Glitch, gm
        #print "inThing = %s, type %s" % (inThing, type(inThing))
        if self.taxNames:
            if taxNames:
                for txNum in range(self.nTax):
                    assert self.taxNames[txNum] == taxNames[txNum], "Mismatched tax names."
        else:
            if taxNames:
                self.taxNames = taxNames
                self.nTax = len(taxNames)
        if type(inThing) == type('string'):
            if not os.path.isfile(inThing):
                gm.append("The inThing is a string, but does not appear to be a file name.")
                raise Glitch, gm
            assumeIsPhylip = None
            f = file(inThing)
            c = ' '
            while c in string.whitespace:
                c = f.read(1)
            f.close()
            if c == '(':
                assumeIsPhylip=True
            if assumeIsPhylip:
                if not taxNames:
                    gm.append("File inThing %s" % inThing)
                    gm.append("assumed to be a phylip tree file.")
                    gm.append("Needs arg taxNames-- an ordered taxNames list.")
                    raise Glitch, gm
                self._readPhylipTreeFile(inThing,skip,max,taxNames)
            else:
                self._readNexusFile(inThing, skip, max)
        elif isinstance(inThing, Tree):
            if skip or max:
                gm.append("Args skip and and max only come into play when reading from files.")
                raise Glitch, gm
            if not inThing.taxNames:
                gm.append("If inThing is a Tree object, it needs taxNames.")
                raise Glitch, gm
            self.taxNames = inThing.taxNames
            self.nTax = len(self.taxNames)
            self._getSplitsFromTree(inThing)
        elif isinstance(inThing, Trees):
            if skip or max:
                gm.append("Args skip and and max only come into play when reading from files.")
                raise Glitch, gm
            if not inThing.taxNames:
                gm.append("If inThing is a Trees object, it needs a taxNames attached.")
                gm.append("(Try 'theTrees.setTaxNames()')")
                raise Glitch, gm
            self.taxNames = inThing.taxNames
            self.nTax = len(self.taxNames)
            if hasattr(inThing.trees[0], "modelInfo"):
                self.doModelComments = 1
                self.modelInfo = inThing.trees[0].modelInfo
            for t in inThing.trees:
                self._getSplitsFromTree(t)
        else:
            gm.append("Sorry-- I can't grok your input.  What is '%s'?" % inThing)
            raise Glitch, gm

        self._finishSplits()

    def _readPhylipTreeFile(self, fName, skip, max, taxNames):
        gm = ['TreePartitions._readPhylipTreeFile()']

        f = file(fName)

        #print 'TreePartitions._readPhylipTreeFile().  var.nexus_doFastNextTok=%s' % var.nexus_doFastNextTok
        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok,safeNextTok,nexusSkipPastNextSemiColon
            from NexusToken2 import checkLineLengths
            checkLineLengths(f)
        else:
            from NexusToken import nextTok,safeNextTok,nexusSkipPastNextSemiColon

        # Skip over trees
        if skip:
            for i in range(skip):
                nexusSkipPastNextSemiColon(f)
            #tok = safeNextTok(f)
            #lowTok = string.lower(tok)
            

        # Read in the trees
        while 1:
            t = Tree()
            t.parseNewick(f, translationHash=None)
            if not t.nodes:
                break
            t.initFinish()
            t.taxNames = taxNames
            if 0:
                print t
                print "got tree:",
                t.write()
            self._getSplitsFromTree(t)
            # The line above increments self.nTrees by the weight.
            if max and self.nTrees >= max:
                break
        f.close()
        self.taxNames = taxNames
        self.nTax = len(self.taxNames)


    def _readNexusFile(self, fName, skip, max):
        gm = ['TreePartitions._readNexusFile()']

        f = file(fName)

        #print 'TreePartitions._readNexusFile().  var.nexus_doFastNextTok=%s' % var.nexus_doFastNextTok
        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok,safeNextTok,nexusSkipPastNextSemiColon
            from NexusToken2 import checkLineLengths
            checkLineLengths(f)
        else:
            from NexusToken import nextTok,safeNextTok,nexusSkipPastNextSemiColon

        
        nf = Nexus()  # provides nf.readBlock(), nf.readTranslateCommand()

        # Read the #nexus
        tok = nextTok(f)
        if tok:
            lowTok = string.lower(tok)
        else:
            gm.append("Empty file?!?")
            raise Glitch, gm
        if lowTok != '#nexus':
            gm.append("Hmmm..., this doesn't appear to be a nexus file. ")
            gm.append("The first token is not '#NEXUS'.")
            raise Glitch, gm


        # Get everything before the trees block.
        inBlock = 0
        while 1:
            tok = safeNextTok(f)
            lowTok = string.lower(tok)
            if not inBlock:
                if lowTok != 'begin':
                    gm.append("Expecting the 'begin' of a nexus block.")
                    raise Glitch, gm
                inBlock = 1
            else:
                if lowTok == 'trees':
                    break
                elif lowTok == 'taxa':
                    #print "TreePartitions._readNexusFile()  Entering taxa block ..."
                    nf.nexusData = NexusData()
                    nf.nexusData.readBlock(f, 'taxa')
                    nexusSkipPastNextSemiColon(f)
                    #print "TreePartitions._readNexusFile()  Finished taxa block."

                    # If we have a taxa block, we can use it to get taxNames
                    self.taxNames = nf.nexusData.taxNames
                    self.nTax = len(self.taxNames)
                    inBlock = 0
                else:
                    gm.append("Unexpected '%s' block." % tok)    # Can't I just skip un-usable blocks?
                    raise Glitch, gm

        # At this point we have read the 'trees' token in the 'begin trees;' line.
        translationHash = None
        nexusSkipPastNextSemiColon(f)
        tok = safeNextTok(f)
        #print "t got tok = %s" % tok
        lowTok = string.lower(tok)
        if lowTok == 'translate':
            translationHash = nf.readTranslateCommand(f)
            if not self.taxNames:
                # If we did not get taxNames from a taxa block, we can
                # get taxNames from the translate command.  News -- this fails for odd keys, anything other than 1-nTax.
                try:
                    self.taxNames = []
                    for i in range(len(translationHash)):
                        self.taxNames.append(translationHash['%s' % (i + 1)])
                except KeyError:
                    #self.taxNames = translationHash.values()
                    gm.append("No taxNames were supplied, so I am getting the taxNames from the translation.")
                    gm.append("However, the keys of the translation do not appear to be numbers from 1 to nTax.")
                    gm.append("I am easily confused, and I don't want to make taxNames in an arbitrary order.")
                    gm.append("So you should supply a taxNames list when you invoke a TreePartitions object.")
                    raise Glitch, gm
                self.nTax = len(self.taxNames)

            # We might expect a p4 command comment if there has been a
            # translate command.
            savedState = var.nexus_getP4CommandComments
            var.nexus_getP4CommandComments = 1
            while lowTok != 'tree':
                tok = safeNextTok(f)
                if tok[0] == '[':
                    #print "got comment: %s" % tok
                    self.modelInfo = _getModelInfo(tok)
                    if self.modelInfo:
                        # self.modelInfo.check() returns
                        #   1 for complete info, but no heterogeneity over the tree
                        #   2 for complete info, and heterogeneity over the tree
                        ret = self.modelInfo.check()
                        if ret == 1:
                            pass # self.doModelComments remains 0
                        elif ret == 2:
                            self.doModelComments = 1
                else:
                    lowTok = string.lower(tok)
            var.nexus_getP4CommandComments = savedState

        if not self.taxNames or not self.nTax:
            gm.append("No taxa block or translation command.")
            gm.append("This file is not suitable for file input.")
            gm.append("(Maybe input it as a Trees object.)")
            raise Glitch, gm

        # At this point, lowTok should be 'tree'
        if lowTok != 'tree':
            if lowTok in ['end', 'endblock']:
                gm.append("Got %s, expecting 'tree' command.  No trees in %s?" % (tok, fName))
                raise Glitch, gm
            else:
                gm.append("I can't grok %s.  I was expecting a 'trees' command." % tok)
                raise Glitch, gm

        # Skip over trees
        if skip:
            for i in range(skip):
                nexusSkipPastNextSemiColon(f)
            tok = safeNextTok(f)
            lowTok = string.lower(tok)

        # Read in the trees
        while 1:
            if lowTok == 'tree':
                t = Tree()
                #print 'initialized t'
                #if 0:
                #    print 'self.taxNames = %s' % self.taxNames
                #    print 'translationHash = %s' % translationHash
                #if self.taxNames and not translationHash:
                #    t._taxNames = self.taxNames 
                if self.doModelComments:
                    t.parseNexus(f, translationHash=translationHash, doModelComments=self.modelInfo.nParts)
                else:
                    t.parseNexus(f, translationHash=translationHash)
                #if t.taxNames:
                #    pass
                #elif self.taxNames:
                #    t.taxNames = self.taxNames
                if self.taxNames:
                    t.taxNames = self.taxNames
                if not t.taxNames:
                    gm.append('Input tree does not have a taxNames list.  Fix me.')
                    raise Glitch, gm
                #t.dump(tree=True)
                self._getSplitsFromTree(t)

                # The line above increments self.nTrees by the weight.
                
                if max and self.nTrees >= max:
                    break
            elif lowTok == 'end' or lowTok == 'endblock':
                break
            else:
                gm.append("I don't understand (lowercased) token '%s'" % lowTok)
                gm.append("I was expecting 'tree' or 'end'.")
                raise Glitch, gm
            tok = safeNextTok(f)
            lowTok = string.lower(tok)
            #print "xx got tok %s" % tok
        f.close()




    def _getSplitsFromTree(self, theTree):
        """Make split objects from theTree, append to self.splits"""


        gm = ['TreePartitions._getSplitsFromTree()']

        if self.doModelComments and theTree.recipWeight and theTree.recipWeight != 1:
            gm.append("doModelComments is set, but the tree has a weight-- should not have both.")
            raise Glitch, gm

        theWeight = 1.0
        if theTree.recipWeight:
            if theTree.recipWeight == 1:
                pass
            else:
                theWeight = 1.0 / theTree.recipWeight

        # In a virgin TreePartitions object, isBiRoot is set to
        # None.  On reading the first tree, it is set to 0 or 1,
        # depending.  Subsequent trees must agree with the first tree.
        nRootChildren = theTree.root.getNChildren()
        if not nRootChildren:
            gm.append("Root has no children.")
            raise Glitch, gm
        elif nRootChildren == 1:
            gm.append("Tree is rooted on a terminal node.  No workee.")
            raise Glitch, gm
        elif nRootChildren == 2:
            if self.isBiRoot == None:
                if 0:
                    gm.append("Got a tree with a bifurcating root.")
                    gm.append("TreePartitions for bi-rooted trees is not working yet.")
                    gm.append("Maybe you could use the Tree.removeRoot() method.")
                    raise Glitch, gm
                self.isBiRoot = 1
            elif self.isBiRoot == 1:
                pass
            else:
                gm.append("Self.isBiRoot has been previously turned off, but now we have a biRoot tree.")
                raise Glitch, gm
        else:
            if self.isBiRoot == None:
                self.isBiRoot = 0
            elif self.isBiRoot == 0:
                pass
            else:
                gm.append("Self.isBiRoot has been previously turned on, but now we have a non-biRoot tree.")
                raise Glitch, gm

        # This next method call, makeSplitKeys(), checks whether any
        # internal nodes have exactly 1 child.  That would be bad cuz
        # it would mean duped splits.
        theTree.makeSplitKeys()
        
        for n in theTree.iterNodesNoRoot():
            theSKey = n.br.splitKey
            #print "node %s, splitKey %s  parent=%s" % (n.nodeNum, theSKey, n.parent.nodeNum)
            
            # Either we have seen this splitKey before, in which
            # case we get it from self.splitsHash, or we have to
            # make a new one.

            if self.splitsHash.has_key(theSKey):
                theSplit = self.splitsHash[theSKey]
            else:
                theSplit = Split()
                self.splits.append(theSplit)
                self.splitsHash[theSKey] = theSplit
                theSplit.key = theSKey

            if not self.isBiRoot:
                theSplit.count += theWeight
                theSplit.cumBrLen += n.br.len * theWeight
            else:
                if n.parent != theTree.root:
                    theSplit.count += theWeight
                    theSplit.cumBrLen += n.br.len * theWeight
                    
                else:
                    # We will be here twice, cuz both children of the
                    # root have the same splitKey.  We want to only
                    # count it once.  So we count it only if the
                    # rawSplitKey is odd.
                    if 1 & n.br.rawSplitKey:
                        theSplit.count += theWeight

                    if self.biSplitsHash.has_key(theSKey):
                        theBiSplit = self.biSplitsHash[theSKey]
                    else:
                        theBiSplit = Split()
                        self.biSplits.append(theBiSplit)
                        self.biSplitsHash[theSKey] = theBiSplit
                        theBiSplit.key = theSKey
                    if 1 & n.br.rawSplitKey:
                        theBiSplit.count += theWeight

                    # We need to save the brLens.  If the rawSplitKey
                    # is odd, we save to theBiSplit.cumBrLen, and if
                    # it is even, we save to biRootCumBrLen.
                    if 1 & n.br.rawSplitKey:
                        theBiSplit.cumBrLen += n.br.len * theWeight
                    else:
                        theBiSplit.biRootCumBrLen += n.br.len * theWeight

                    if self.doModelComments:
                        # We save the model stuff on the "odd" branch
                        # to theBiSplit.modelUsage.  We save the model
                        # stuff on the "even" side to
                        # theBiSplit.biRootModelUsage.  We save the
                        # rootModelUsage to theBiSplit.rootModelUsage.
                        if 1 & n.br.rawSplitKey:
                            if theBiSplit.modelUsage:
                                pass
                            else:
                                theBiSplit.modelUsage = SplitModelUsage(self.modelInfo)
                            smt = theBiSplit.modelUsage
                            for pNum in range(self.modelInfo.nParts):
                                if n.parts[pNum].compNum >= 0:
                                    smt.parts[pNum].compCounts[n.parts[pNum].compNum] += 1
                                if n.br.parts[pNum].rMatrixNum >= 0:
                                    smt.parts[pNum].rMatrixCounts[n.br.parts[pNum].rMatrixNum] += 1
                                if n.br.parts[pNum].gdasrvNum >= 0:
                                    smt.parts[pNum].gdasrvCounts[n.br.parts[pNum].gdasrvNum] += 1

                            # If the rawSplitKey is odd, which it is, get rootModelUsage
                            if theBiSplit.rootModelUsage:
                                pass
                            else:
                                theBiSplit.rootModelUsage = SplitModelUsage(self.modelInfo)
                            smt = theBiSplit.rootModelUsage
                            for pNum in range(self.modelInfo.nParts):
                                if theTree.root.parts[pNum].compNum >= 0:
                                    smt.parts[pNum].compCounts[theTree.root.parts[pNum].compNum] += 1

                        else:
                            if theBiSplit.biRootModelUsage:
                                pass
                            else:
                                theBiSplit.biRootModelUsage = SplitModelUsage(self.modelInfo)
                            smt = theBiSplit.biRootModelUsage
                            for pNum in range(self.modelInfo.nParts):
                                if n.parts[pNum].compNum >= 0:
                                    smt.parts[pNum].compCounts[n.parts[pNum].compNum] += 1
                                if n.br.parts[pNum].rMatrixNum >= 0:
                                    smt.parts[pNum].rMatrixCounts[n.br.parts[pNum].rMatrixNum] += 1
                                if n.br.parts[pNum].gdasrvNum >= 0:
                                    smt.parts[pNum].gdasrvCounts[n.br.parts[pNum].gdasrvNum] += 1
                        

            if self.doModelComments:
                #print "node %i, self.isBiRoot=%s" % (n.nodeNum, self.isBiRoot)
                #print "node %i, self.modelInfo=%s" % (n.nodeNum, self.modelInfo)
                # We can be sure that theWeight is 1.0
                if theSplit.modelUsage:
                    pass
                else:
                    theSplit.modelUsage = SplitModelUsage(self.modelInfo)

                smt = theSplit.modelUsage
                for pNum in range(self.modelInfo.nParts):
                    if n.parts[pNum].compNum >= 0:
                        smt.parts[pNum].compCounts[n.parts[pNum].compNum] += 1
                    if n.br.parts[pNum].rMatrixNum >= 0:
                        smt.parts[pNum].rMatrixCounts[n.br.parts[pNum].rMatrixNum] += 1
                    if n.br.parts[pNum].gdasrvNum >= 0:
                        smt.parts[pNum].gdasrvCounts[n.br.parts[pNum].gdasrvNum] += 1

                    

        #self.dump()
        #sys.exit()

        if not self.isBiRoot:
            for n in theTree.iterNodesNoRoot():
                theSplit = self.splitsHash[n.br.splitKey]
                # If the split is off the root and its rawSplitKey contains a
                # 1 (ie has the first taxon), then the corresponding split will have
                # its rootCount incremented.

                # That means that when a consensus tree is made and
                # rooted on the first taxon, the root should be placed
                # on the node with the split with the highest rootCount.
                if n.parent == theTree.root and (1 & n.br.rawSplitKey):
                    if n.br.rawSplitKey == 1:
                        # We do *not* want to do this next line.
                        #theSplit.rootCount += theWeight

                        # Hack alert! This is an inelegant solution.  This
                        # tree is rooted on the parent of the first taxon.
                        # We want to save the rootCount, but saving it on
                        # this split is wrong, because this split is
                        # guarranteed to be in every cons tree.  If this
                        # particular rooting does not make it to the cons
                        # tree then we do not want its rootCount (and we
                        # *will* get its rootCount if we put its root
                        # count on this split).

                        # We (generally) increment the rootCount of the
                        # split of the node with the clade containing the
                        # first taxon, but the reason is arbitrary,
                        # ultimately.  We could increment the rootCount of
                        # another split, as long as we only do it for one
                        # node.  So that is what we do here-- we store the
                        # info in another split-- another split that is
                        # off the root.  But it is not that simple--
                        # instead of incrementing the rootCount of
                        # theSplit, we increment the rootCount2 of some
                        # other split off the root.  That at least saves
                        # the information that the parent of that split
                        # was a root.

                        # Later, after we make a cons tree, if that split
                        # does not make it into the cons tree, then its
                        # rootCount2 is lost, as it should be.  If it does
                        # make it to the cons tree, then we increment the
                        # rootCount of its parent.
                        aNode = None
                        for bNode in theTree.root.iterChildren():
                            # We want a child of the root, that is an internal node.
                            if bNode != n and not bNode.isLeaf:
                                aNode = bNode
                                break
                        if aNode:
                            self.splitsHash[aNode.br.splitKey].rootCount2 += theWeight
                        else:
                            pass
                            #gm.append("This tree appears to have no splits!")
                            #raise Glitch, gm



                    else:
                        theSplit.rootCount += theWeight

                    if self.doModelComments:
                        # We can be sure that theWeight is 1.0

                        # Now we want to get the root model info.
                        # That info is currently residing in
                        # theTree.root.modelUsage, but I want to store
                        # it in theSplit.rootModelUsage.

                        if theSplit.rootModelUsage:
                            pass
                        else:
                            theSplit.rootModelUsage = SplitModelUsage(self.modelInfo)
                        smt = theSplit.rootModelUsage
                        for pNum in range(self.modelInfo.nParts):
                            if theTree.root.parts[pNum].compNum >= 0:
                                smt.parts[pNum].compCounts[theTree.root.parts[pNum].compNum] += 1

        self.nTrees += theWeight





    def _finishSplits(self):
        """Having made all the splits, now do some finalizing adjustments.

        Make split.string
        Make split.set
        Make split.proportion
        Order the splits based on proportion and the string.
        Put biRoot counts on appropriate branches.
        """

        #gm = ['TreePartitions._finishSplits()']
        #print "self.nTrees = %s" % self.nTrees
        for spl in self.splits:
            spl.string = func.getSplitStringFromKey(spl.key, self.nTax)
            listForSet = []
            for i in range(len(spl.string)):
                if spl.string[i] == '*':
                    # It is 1-based because the first thing that I do
                    # when I make a consensus tree is to make a comb
                    # tree.  Numbers on that comb tree then correspond
                    # to the numbers in the 1-based sets.  Otherwise
                    # it does not matter--- it just makes the making
                    # of cons trees easier to debug.  If it is changed
                    # to 0-based, it still works.
                    listForSet.append(i + 1)  # 1-based

            spl.set = set(listForSet)
            spl.proportion = spl.count / float(self.nTrees)

        # We want the splits to be ordered on both the reverse
        # proportion, and the split string, with stars before dots
        # (arbitrarily, so it looks nice, and is consistent)
        for spl in self.splits:
            spl.proportion = -1.0 * spl.proportion
        self.splits = func.sortListOfObjectsOn2Attributes(self.splits, 'proportion', 'string')
        #self.splits = func.sortListOfObjectsOn2Attributes(self.splits, 'proportion', 'key')
        for spl in self.splits:
            spl.proportion = -1.0 * spl.proportion

        # Finish self.biSplits
        for spl in self.biSplits:
            spl.string = func.getSplitStringFromKey(spl.key, self.nTax)
            spl.proportion = spl.count / float(self.nTrees)


    def dump(self):
        #print "\nTreePartitions dump: \n\t(to get the full translationHash and splits, do 'writeSplits')"
        #tH = '%s' % self.translationHash
        #if len(tH) > 30:
        #    tH = tH[:30] + ' ...'
        #print "%25s = %s' % ('translationHash", tH)
        if 1:
            print "%25s = %s" % ('nTax', self.nTax)
            print "%25s = %s" % ('nTrees', self.nTrees)
            print "%25s = %s" % ('number of splits', len(self.splits))
            print "%25s = %s" % ('isBiRoot', self.isBiRoot)
        #print "%25s = %s' % ('nParts", self.nParts)
        #print "%25s = %s' % ('modelNames", self.modelNames)
        if 0:
            for s in self.splits:
                s.dump()
        if 0:
            print
            print "%12s %12s %12s %12s %12s %12s %12s" % (
                'string', 'key', 'count', 'rootCount', 'rootCount2', 'cumBrLen', 'biRtCumBrLen')
            for s in self.splits:
                print "%12s %12s %12s %12s %12s %12s %12s" % (
                    s.string, s.key, s.count, s.rootCount, s.rootCount2, s.cumBrLen, '-')
            if len(self.biSplits):
                print "biSplits"
                for s in self.biSplits:
                    print "%12s %12s %12s %12s %12s %12s %12s" % (
                        s.string, s.key, s.count, s.rootCount, s.rootCount2, s.cumBrLen, s.biRootCumBrLen)

        if 1:
            print
            print "%12s %12s %12s %12s %12s %12s" % (
                'string', 'key', 'count', 'modelUsage', 'biRtMdlUsg', 'rtModelUsage')
            for s in self.splits:
                print "%12s %12s %12s" % (
                    s.string, s.key, s.count),
                if s.modelUsage:
                    print "%12s" % s.modelUsage.nParts,
                else:
                    print "%12s" % "None",
                if s.biRootModelUsage:
                    print "%12s" % s.biRootModelUsage.nParts,
                else:
                    print "%12s" % "None",
                if s.rootModelUsage:
                    print "%12s" % s.rootModelUsage.nParts
                else:
                    print "%12s" % "None"
                    
            if len(self.biSplits):
                print "biSplits"
                for s in self.biSplits:
                    print "%12s %12s %12s" % (
                        s.string, s.key, s.count),
                    if s.modelUsage:
                        print "%12s" % s.modelUsage.nParts,
                    else:
                        print "%12s" % "None",
                    if s.biRootModelUsage:
                        print "%12s" % s.biRootModelUsage.nParts,
                    else:
                        print "%12s" % "None",
                    if s.rootModelUsage:
                        print "%12s" % s.rootModelUsage.nParts
                    else:
                        print "%12s" % "None"

        

    def writeSplits(self, fName=None, minimumProportion=0.05, doLeaves=True):
        """Write a table of splits, in dot-star notation.

        Writes all the splits with a proportion >= the minimumProportion."""
        
        #print "\nTree bipartitions."

        if not fName:
            f = sys.stdout
        else:
            f = file(fName, 'w')
            
        for i in range(len(self.taxNames)):
            f.write("  %5i   %s\n" % ((i + 1), self.taxNames[i]))

        # Print headings for the table
        f.write("\nSplits\n======\n")
        firstColWidth = self.nTax
        if firstColWidth < 12:
            firstColWidth = 12
        firstColSig = '%%-%is' % firstColWidth

        if not self.isBiRoot:
            f.write(firstColSig % ' ')
            f.write("    count proportion brLen\n")
            nSplits = len(self.splits)
            if 0:
                for i in range(nSplits):
                    p = self.splits[i]
                    f.write(firstColSig % p.string)
                    f.write("%9.3f    %5.3f   %5.3f\n" % (p.count, p.proportion, (p.cumBrLen/p.count)))
                print "-" * 50
            for i in range(nSplits):
                p = self.splits[i]
                if not doLeaves:
                    theStarCount = p.string.count('*')
                    if theStarCount == 1 or theStarCount == self.nTax - 1:
                        continue
                if p.proportion >= minimumProportion:
                    f.write(firstColSig % p.string)
                    f.write("%9.3f    %5.3f   %5.3f\n" % (p.count, p.proportion, (p.cumBrLen/p.count)))
                else:
                    break
            if i < nSplits - 1:
                f.write("(plus %i splits with proportion less than %.3f)\n" % ((nSplits - i), minimumProportion))

        else:
            # This needs to be fixed to take into account minimumProportion
            # The brLen depends on where its rooted, so don't report it.
            f.write(firstColSig % ' ')
            f.write("    count proportion\n")
            for p in self.splits:
                f.write(firstColSig % p.string)
                f.write("%9.3f    %5.3f\n" % (p.count, p.proportion))


            if len(self.biSplits):
                f.write("\nBiSplits (ie splits involving a bifurcating root)\n========\n")
                f.write(firstColSig % ' ')
                f.write("    count proportion\n")
                for p in self.biSplits:
                    f.write(firstColSig % p.string)
                    f.write("%9.3f    %5.3f\n" % (p.count, p.proportion))

        if 0 and self.modelInfo:  # This needs to be fixed up to show modelUsage for biSplits.
            f.write('\n')
            self.modelInfo.dump()
            f.write('\n')

            for pNum in range(self.modelInfo.nParts):
                f.write("partNum %i\n" % pNum)
                mip = self.modelInfo.parts[pNum]
                if mip.nComps > 1 or mip.nRMatrices > 1 or mip.nGdasrvs > 1:
                    for p in self.splits:
                        f.write("%s  " % p.string)
                        pmt = p.modelUsage
                        if pmt:
                            if mip.nComps > 1:
                                f.write("compCounts=%s " % pmt.parts[pNum].compCounts)
                            if mip.nRMatrices > 1:
                                f.write("rMatrixCounts=%s " % pmt.parts[pNum].rMatrixCounts)
                            if mip.nGdasrvs > 1:
                                f.write("gdasrvCounts=%s " % pmt.parts[pNum].gdasrvCounts)
                        else:
                            #f.write("no modelUsage for this split.")
                            pass
                        print
                else:
                    f.write("part %i is tree-homogeneous\n" % pNum)
                f.write('\n')
        if f != sys.stdout:
            f.close()
                            




    def compareSplits(self, otherTP, noLeaves=True, minimumProportion=0.1, bothMustMeetMinimum=False):
        """Compare split support with another TreePartitions.

        It returns a list of comparisons, each comparison being a list of 3 things:
          1. The split key
          2. The split string
          3. A list of the 2 supports
        
        Lets say you have 2 MCMC runs, and you want to compare the
        split supports for them.  You make a TreePartitions object for
        each, and then use this method to get a list of comparisons.
        A single comparison is the support from each of the two MCMC
        runs for a particular split.  If one MCMC run has a split but
        the other does not, it will still make a pair, with one of the
        elements being zero.

        noLeaves
                            means that trivial splits, ie for leaves,
                            are excluded.
                            
        minimumProportion
                            - If bothMustMeetMinimum is turned off,
                              which is now the default (and is more
                              informative), then a comparison is not
                              included unless one support value is
                              greater than or equal to the minimum
                              proportion
                            
                            - If bothMustMeetMinimum is turned on,
                              then a comparison is not included unless
                              both supports are at or above the
                              minimumProportion.

        """

        gm = ['TreePartitions.compareSplits()']
        
        if self.nTax != otherTP.nTax:
            gm.append("Mismatched taxa.")
            raise Glitch, gm
            
        for i in range(len(self.taxNames)):
            if self.taxNames[i] != otherTP.taxNames[i]:
                gm.append("Mismatched taxa.")
                raise Glitch, gm

        if 0:
            self.writeSplits(minimumProportion=minimumProportion, doLeaves=False)
            otherTP.writeSplits(minimumProportion=minimumProportion, doLeaves=False)

        pairs = []

        # Compare all the splits in self with otherTP
        for s in self.splits:
            theString = s.string
            if noLeaves:
                theStarCount = theString.count('*')
                if theStarCount == 1 or theStarCount == self.nTax - 1:
                    continue
            theKey = s.key
            prop1 = s.proportion
            
            if otherTP.splitsHash.has_key(theKey):
                prop2 = otherTP.splitsHash[theKey].proportion
            else:
                prop2 = 0.0

            if bothMustMeetMinimum:
                # The following is my old way, where both must meet the minimum.  Changed 18 May 2011.
                if (prop1 < minimumProportion) or (prop2 < minimumProportion):
                    continue
            else:
                # The following is more informative -- at least one split greater than minpartfreq
                if (prop1 < minimumProportion) and (prop2 < minimumProportion):
                    continue

            #print "%s  %6.3f %6.3f" % (theKey, prop1, prop2)
            pairs.append([theKey, theString, [prop1, prop2]])

        # otherTP might have splits that self does not have.  Compare those.
        for s in otherTP.splits:
            theString = s.string
            if noLeaves:
                theStarCount = theString.count('*')
                if theStarCount == 1 or theStarCount == self.nTax - 1:
                    continue
            theKey = s.key
            if not self.splitsHash.has_key(theKey):
                prop1 = 0.0
                prop2 = s.proportion
                if bothMustMeetMinimum:
                    if (prop1 < minimumProportion) or (prop2 < minimumProportion):
                        continue
                else:
                    if (prop1 < minimumProportion) and (prop2 < minimumProportion):
                        continue
                #print "%s  %6.3f %6.3f" % (theKey, prop1, prop2)
                pairs.append([theKey, theString, [prop1, prop2]])

        return pairs


    def consensus(self, conTreeName='consensus', showRootInfo=0, minimumProportion=None):
        """Make a consensus tree.

        This method assembles a Tree object from the tree partitions
        contained in self.splits.  It makes a 'majority-rule consensus
        with extra compatible splits'.  It will use splits that have
        more than minimumProportion support (note that
        minimumProportion is given as a number between 0 and 1.0,
        often it will be 0.5, meaning 50%).

        It returns the tree.  Do it like this:
            tp = TreePartitions('myTreeFile.nex')
            tp.writeSplits() # If you want ...
            t = tp.consensus()

        The support for the splits (often the Bayesian posterior
        support or bootstrap support, as given in the attribute
        split.proportion) is put on the tree as node.br.support as a
        float, not a string.  If you want to draw it later showing the
        support, you can format the float as you like (eg as percent)
        and put it on node.br.name (for example) for making a picture,
        or as node.name if you want to save it in Nexus format with
        that info included.  You can do that like this::

            for n in t.iterInternalsNoRoot():
                if n.br.support != None:
                    # The signature '%.0f' gives zero decimal places;
                    # change that if you want more precision.
                    n.name = '%.0f' % (100. * n.br.support)
            t.writeNexus('consTreeWithSupport.nex')

        A feature of this consensus method is that the resulting tree
        is rooted on the majority root position of the input trees.
        (If there is more than one majority root position, it is
        rooted on the first majority root position, arbitrarily.)  If
        showRootInfo is set, a table of root positions is printed out.

        The resulting con tree is decorated with the rootCounts.  With
        trifurcating roots, they get put in the node.rootCount, and
        that shows directly which nodes were roots.  If it is a
        bifurcating root the root counts get put in
        node.br.biRootCount, and shows how many times, for those
        splits that made it into the cons tree, that the root was
        found on that branch.

        Since the resulting tree often contains more info than can be
        accommodated in Newick format, you may want to save it by
        pickling it.  Use::
        
            t.tPickle()
        """


        gm = ['TreePartitions.consensus()']

        if self.nTax <= 3:
            gm.append("Only %i taxa? -- too few to consider for a consensus" % self.nTax)
            raise Glitch, gm

        #############################
        # Make a star tree
        #############################

        conTree = Tree()
        conTree.name = conTreeName

        # Make a comb tree.  Right after I finish making it, the tree
        # will then be rooted on the first taxon, but for construction
        # purposes here it will be rooted on the single internal node.
        conTree.root = Node()
        conTree.root.nodeNum = 0
        conTree.root.isLeaf = 0
        conTree.nodes.append(conTree.root)

        n = Node()
        n.nodeNum = 1
        conTree.nodes.append(n)
        conTree.root.leftChild = n
        n.parent = conTree.root
        n.isLeaf = 1
        n.name = self.taxNames[0]
        n.seqNum = 0
        n.br.splitKey = 2 ** self.nTax - 2
        previous = n

        for i in range(self.nTax - 1):
            n = Node()
            n.nodeNum = i + 2
            conTree.nodes.append(n)
            previous.sibling = n
            n.parent = conTree.root
            n.isLeaf = 1
            n.name = self.taxNames[i + 1]
            n.seqNum = i + 1
            n.br.splitKey = 2 ** (i + 1)
            previous = n

        conTree.reRoot(1, moveInternalName=False)
        # Now its a comb on a stick.

        if 0:
            conTree.preAndPostOrderAreValid = 0
            conTree.draw()
            for n in conTree.nodes:
                if n != conTree.root:
                    print "%3i  %s" % (n.nodeNum, n.br.splitKey)
                else:
                    print "%3i  root" % n.nodeNum
            conTree.preOrder = None
            conTreePostOrder = None
            conTree.preAndPostOrderAreValid = 0
            #sys.exit()


        # #####################################
        # Add splits
        # #####################################

        #nTaxSet = set(range(self.nTax+ 1)[1:])
        nInternalNodes = 1

        for spl in self.splits:
            # Skip the trivial splits-- they are already in the star tree.
            if len(spl.set) <= 1 or len(spl.set) == (self.nTax - 1):
                #print "...skipping split %s" % spl.set
                continue

            # If the tree is fully resolved, stop.
            if nInternalNodes >= self.nTax -2:
                #print "Tree is fully resolved, so stop."
                break

            # If the split proportion is less than minimumProportion ...
            if minimumProportion:
                if spl.proportion < minimumProportion:
                    #print " ... is less than minimumProportion, so break"
                    break

            #print "\nLooking at split %s, %s, %s" % (spl.string, spl.key, spl.set)

            # Check whether spl is compatible with the tree splits so
            # far.  (This was not in the Page code, that I could find,
            # and its absence leads to grief and awkwardness.)

            if nInternalNodes > 1: # The first one will be compatible for sure.
                #print "    Checking for compatibility."
                isCompatible = 1
                for n in conTree.nodes[self.nTax + 1:]: # No need to check the first trivial nodes.
                    nSet = self.splitsHash[n.br.splitKey].set
                    #print "      ...checking   %s       %s" % (n.br.splitKey, nSet)

                    # A split is incompatible (for our purposes here) if:
                    #
                    #    1. The sets overlap, ie the intersection is non-empty
                    #    2. nSet.difference(spl.set) is non-empty
                    #    3. spl.set.difference(nSet) is non-empty

                    # (For completeness, the other criterion for
                    # incompatibility is that both splits (nSet and
                    # spl.set) have at least one shared taxon missing.
                    # Well, both are missing taxon 1, so that is
                    # always fulfilled, and does not need to be
                    # checked.  Thanks to Mark Wilkinson for defining
                    # split compatibility like this.)

                    if len(nSet.intersection(spl.set)) \
                           and len(nSet.difference(spl.set)) \
                           and len(spl.set.difference(nSet)):
                        #print "        %s is incompatible with" % n.br.splitKey
                        #print "        %s" % spl.string
                        isCompatible = 0
                        break
                if not isCompatible:
                    continue        # next spl

            # If we made it this far, spl is compatible.
            # Find a place in the tree in which to insert the split.
            # This is the crux of Page's code.
            # Thanks, Rod!

            p = conTree.root
            q = p.leftChild
            isDone = 0
            addThisNode = 0

            while q and not isDone:
                relationship = None
                #print "q.br.splitKey = %s" % q.br.splitKey
                qSet = self.splitsHash[q.br.splitKey].set
                #print "qSet = %s" % qSet
                #print "spl.set = %s" % spl.set

                # 1. spl is identical to q, but that should never happen
                if q.br.splitKey == spl.key:
                    gm.append("Candidate split is identical to an existing split.")
                    gm.append("This should never happen.  Programming error.")
                    raise Glitch, gm

                # 2. spl is a subset of q, then choose a new p and q (as children of previous p and q) and try again
                elif spl.set.issubset(qSet):
                    relationship = 'SUBSET'
                    p = q
                    q = q.leftChild

                # 3. spl is disjoint to q, then choose a new p (child of p) and q (sib of q) and try again
                elif not len(spl.set.intersection(qSet)):
                    relationship = 'DISJOINT'
                    p = q
                    q = q.sibling

                # 4. spl is a superset of q, in which case make a new node between p and q
                elif spl.set.issuperset(qSet):
                    #print "spl.set (%s) is a superset of qSet (%s)" % (spl.set.__getstate__()[0].keys(),
                    #                                                   qSet.__getstate__()[0].keys())
                    relationship = 'SUPERSET'
                    isDone = 1
                    addThisNode = 1

                # 5. spl overlaps q, so spl cant be part of the tree
                elif len(spl.set.intersection(qSet)):
                    #relationship = 'OVERLAPPING'
                    #print "...Got overlapping..."
                    isDone = 1

                else:
                    gm.append("xxx Why would this happen?")
                    raise Glitch, gm

            #print "a relationship is %s" % relationship
            if addThisNode:
                #print "    Adding a new node. spl.string= %s" % spl.string
                if not q:
                    gm.append("Programming error. q does not exist.")
                    raise Glitch, gm
                n = Node()
                n.nodeNum = len(conTree.nodes)
                conTree.nodes.append(n)
                nInternalNodes = nInternalNodes + 1
                n.br.splitKey = spl.key
                # I don't think that either SUBSET nor DISJOINT will
                # happen.  It will never be SUBSET or DISJOINT if I
                # start with a star tree.  However, I'll leave this
                # code in, in case I ever decide to not start with a
                # star tree.
                if relationship == 'SUBSET':
                    #print "b relationship is SUBSET"
                    p.leftChild = n
                    n.parent = p
                elif relationship == 'DISJOINT':
                    #print "b relationship is DISJOINT"
                    p.sibling = n
                    n.parent = p.parent
                elif relationship == 'SUPERSET':
                    #print "b relationship is SUPERSET"
                    if q == p.leftChild:
                        #print "here 1"
                        p.leftChild = n
                        n.leftChild = q
                        n.parent = p
                        q.parent = n
                    else:
                        #print "here 2. p = %s, n = %s, q = %s" % (p.nodeNum, n.nodeNum, q.nodeNum)
                        p.sibling = n
                        n.leftChild = q
                        n.parent = p.parent
                        q.parent = n

                    tmp = q
                    s = q.sibling
                    t = q.parent

                    while s:
                        sSet = self.splitsHash[s.br.splitKey].set
                        if s.br.splitKey == spl.key:
                            print "why would this ever happen?"
                            raise Glitch, gm
                        elif sSet.issubset(spl.set):
                            #print "here 3"
                            s.parent = n
                            tmp = s
                            s = s.sibling
                        else:
                            #print "here 4"
                            t.sibling = s
                            tmp.sibling = s.sibling
                            t = s
                            t.sibling = None
                            s = tmp.sibling
                    if tmp == q:
                        gm.append("Programming error.  This shouldn't happen.")
                        gm.append("No further subsets have been found.  A useless node has been added")
                        raise Glitch, gm

        ###################################################
        ###################################################

        # At this point, we have a (perhaps fully resolved?) tree,
        # rooted on the first taxon.
        if 0:
            conTree.draw()
            print "%12s %12s %12s %12s %12s %12s" % ('nodeNum', '', 'count', 'rootCount', 'rootCount2', 'cumBrLen')
            for n in conTree.iterNodesNoRoot():
                s = self.splitsHash[n.br.splitKey]
                print "%12s %12s %12s %12s %12s %12s" % (
                    n.nodeNum, s.string, s.count, s.rootCount, s.rootCount2, s.cumBrLen)
            sys.exit()

        # For convenience, temporarily attach the split objects to
        # node.br's; then I don't need to keep looking up the split
        # in the self.splitsHash.  The n.br.split's are deleted at the
        # end of this method, below.
        conTree.setPreAndPostOrder()
        for n in conTree.iterNodesNoRoot():
            n.br.split = self.splitsHash[n.br.splitKey]

        # We need to deal with the rootCount2 counts.  When the tree
        # is rooted on the first taxon, which it is, then the parent
        # of any nodes with a rootCount2 should have its rootCount
        # incremented.  Its a hack originating in
        # _getSplitsFromTree(), above.
        if not self.isBiRoot:
            for n in conTree.iterNodesNoRoot():
                if n.br.split.rootCount2:  # The single child of the root will never have a rootCount2.
                    n.parent.br.split.rootCount += n.br.split.rootCount2

        # Get br.lens
        if self.isBiRoot:
            for n in conTree.iterNodesNoRoot():
                # This is a little tricky.  We are only dealing with
                # non-root splits here.  The cumBrLen's only come from
                # non-root splits, so we can use it.  However, the
                # counts are composed of non-root splits plus root
                # splits; we only want the non-root splits.  We can
                # get the contribution of the root splits from the
                # count from the biRootHash.
                biRootCountContribution = 0.0
                if self.biSplitsHash.has_key(n.br.split.key):
                    biRootCountContribution = self.biSplitsHash[n.br.split.key].count
                theCount = n.br.split.count - biRootCountContribution
                if theCount:
                    n.br.len = n.br.split.cumBrLen / theCount
        else:
            for n in conTree.iterNodesNoRoot():
                n.br.len = n.br.split.cumBrLen / n.br.split.count

        # Get splitSupport
        for n in conTree.nodes[2:]: # skip the (terminal) root, and the node 0
            if not n.isLeaf:
                #print "setting n %i support to %s" % (n.nodeNum, n.br.split.proportion)
                n.br.support = n.br.split.proportion
            conTree.nodes[0].br.support = -1.0

        # Get rootCount and biRootCount
        if self.isBiRoot:
            for n in conTree.iterNodesNoRoot():
                if self.biSplitsHash.has_key(n.br.split.key):
                    n.br.biRootCount = self.biSplitsHash[n.br.split.key].count
                else:
                    n.br.biRootCount = 0.0
                n.rootCount = -1.0
                n.br.support = -1.0
                    
        else:
            for n in conTree.iterNodesNoRoot():
                if n.br.split.rootCount:
                    n.rootCount = n.br.split.rootCount

        if 0:
            conTree.preAndPostOrderAreValid = 0
            conTree.draw()
            print "%12s %12s %12s %12s %12s %12s" % (
                'nodeNum', '', 'n.br.len', 'n.br.support', 'n.rootCount', 'n.br.biRtCnt')
            for n in conTree.iterNodesNoRoot():
                print "%12s %12s %12s %12s %12s %12s" % (
                    n.nodeNum, n.br.split.string, n.br.len, n.br.support, n.rootCount, n.br.biRootCount)
            #sys.exit()

        ######################
        # Re-root
        ######################

        # Put together a maxRootNodes list, in 2 steps.
        # First, without identifying which node or nodes have it, find
        # the the maxRootCount.  Then find which nodes in the conTree
        # have that maxRootCount.  There will be at least one, and
        # maybe more.  Then re-root.
        if self.isBiRoot:
            maxRootCount = 0
            for n in conTree.iterNodesNoRoot():
                if 0:
                    print "%3i   %i" % (n.nodeNum, n.br.biRootCount)
                if n.br.biRootCount > maxRootCount:
                    maxRootCount = n.br.biRootCount
            maxRootNodes = []
            for n in conTree.iterNodesNoRoot():
                if n.br.biRootCount == maxRootCount:
                    maxRootNodes.append(n)
            if 0:
                print "maxRootNodes = ",
                for n in maxRootNodes:
                    print n.nodeNum,
                print

            if 0 and showRootInfo:
                print "There are %s maxRootNodes.  Choosing the first to root on." % len(maxRootNodes)
            biRootChild = maxRootNodes[0]
            theBiSplit = self.biSplitsHash[biRootChild.br.splitKey]
            #print "    Max root node, ie biRootChild, is node %i" % biRootChild.nodeNum
            # Add a root node.
            theBiRoot = conTree.addNodeBetweenNodes(biRootChild, biRootChild.parent)
            theBiRoot.br.split = Split()
            
            
            theBiRoot.br.len = theBiSplit.cumBrLen / theBiSplit.count
            theBiRoot.br.split.modelUsage = theBiSplit.modelUsage
            
            biRootChild.br.len = theBiSplit.biRootCumBrLen / theBiSplit.count
            biRootChild.br.split.modelUsage = theBiSplit.biRootModelUsage

            if 0:
                print "biRootChild is node %i" % biRootChild.nodeNum
                print "theBiRoot.rootModelUsage = %s" % theBiRoot.rootModelUsage
                sys.exit()
            conTree.reRoot(theBiRoot, moveInternalName=False)
            theBiRoot.rootModelUsage = theBiSplit.rootModelUsage # Might be None, if there was no model info

        else:
            maxRootCount = 0
            for n in conTree.iterNodesNoRoot():
                if 0:
                    print "%3i   %i" % (n.nodeNum, n.br.split.rootCount)
                if n.br.split.rootCount > maxRootCount:
                    maxRootCount = n.br.split.rootCount
            maxRootNodes = []
            for n in conTree.iterNodesNoRoot():
                if n.br.split.rootCount == maxRootCount:
                    maxRootNodes.append(n)

            if 0 and showRootInfo:
                print gm[0]
                if len(maxRootNodes) > 1:
                    print "    There are %s maxRootNodes.  Choosing the first to root on." % len(maxRootNodes)
            maxRootNodes[0].rootModelUsage = maxRootNodes[0].br.split.rootModelUsage
            conTree.reRoot(maxRootNodes[0], moveInternalName=False)

        # Put the nodes in pre-order, only for cosmetic purposes.
        conTree.setPreAndPostOrder()
        newNodes = []
        for i in conTree.preOrder:
            newNodes.append(conTree.nodes[int(i)])
        for i in range(len(newNodes)):
            newNodes[int(i)].nodeNum = int(i)
        conTree.nodes = newNodes

        if 0:
            conTree.preAndPostOrderAreValid = 0
            conTree.draw()

        # Now tabulate model usage.
        if self.modelInfo:
            for pNum in range(self.modelInfo.nParts):
                print "\nPartition %i" % pNum

                if self.modelInfo.parts[pNum].nComps > 1:
                    print "\nc%i.%i" % (pNum, self.modelInfo.parts[pNum].nComps)
                    print "           Composition Numbers"
                    print "%8s" % 'NodeNum',
                    for cNum in range(self.modelInfo.parts[pNum].nComps):
                        print "%4s " % cNum,
                    print
                    print "%8s" % '=======',
                    for cNum in range(self.modelInfo.parts[pNum].nComps):
                        print " %4s" % '===',
                    print

                    for n in conTree.nodes:
                        if n != conTree.root:
                            print "%6i  " % n.nodeNum,
                            theCompCounts = n.br.split.modelUsage.parts[pNum].compCounts
                            for cNum in range(self.modelInfo.parts[pNum].nComps):
                                print "%4i " % theCompCounts[cNum],
                            print
                        else:
                            print "%6i  " % n.nodeNum,
                            #theCompCounts = maxRootNodes[0].br.split.rootModelUsage.parts[pNum].compCounts
                            theCompCounts = conTree.root.rootModelUsage.parts[pNum].compCounts
                            for cNum in range(self.modelInfo.parts[pNum].nComps):
                                print "%4i " % theCompCounts[cNum],
                            print
                else:
                    print "Comps in this partition are homogeneous."

                if self.modelInfo.parts[pNum].nRMatrices > 1:
                    print "\nr%i.%i" % (pNum, self.modelInfo.parts[pNum].nRMatrices)
                    print "           rMatrix Numbers"
                    print "%8s" % 'NodeNum',
                    for rNum in range(self.modelInfo.parts[pNum].nRMatrices):
                        print "%4s " % rNum,
                    print
                    print "%8s" % '=======',
                    for rNum in range(self.modelInfo.parts[pNum].nRMatrices):
                        print " %4s" % '===',
                    print

                    for n in conTree.nodes:
                        if n != conTree.root:
                            print "%6i  " % n.nodeNum,
                            theRMatrixCounts = n.br.split.modelUsage.parts[pNum].rMatrixCounts
                            for rNum in range(self.modelInfo.parts[pNum].nRMatrices):
                                print "%4i " % theRMatrixCounts[rNum],
                            print
                else:
                    print "rMatrices in this partition are homogeneous."


        # Find the majority compNum's and rMatrixNum's.
        if self.modelInfo:
            # First make a place to put them.  Initialize n.compNum
            # etc with -1's.
            for n in conTree.nodes:
                if n.parts:
                    gm.append('node already has parts.')
                    raise Glitch, gm
                for pNum in range(self.modelInfo.nParts):
                    n.parts.append(NodePart())
                if n != conTree.root:
                    for pNum in range(self.modelInfo.nParts):
                        n.br.parts.append(NodeBranchPart())

            # Now find the majority
            for pNum in range(self.modelInfo.nParts):
                if self.modelInfo.parts[pNum].nComps > 1:
                    for n in conTree.nodes:
                        if n != conTree.root:
                            theCompCounts = n.br.split.modelUsage.parts[pNum].compCounts
                            n.parts[pNum].compNum = theCompCounts.index(max(theCompCounts))
                        else:
                            theCompCounts = conTree.root.rootModelUsage.parts[pNum].compCounts
                            n.parts[pNum].compNum = theCompCounts.index(max(theCompCounts))
                if self.modelInfo.parts[pNum].nRMatrices > 1:
                    for n in conTree.iterNodesNoRoot():
                        theRMatrixCounts = n.br.split.modelUsage.parts[pNum].rMatrixCounts
                        n.br.parts[pNum].rMatrixNum = theRMatrixCounts.index(max(theRMatrixCounts))

            # Draw them
            conTree.setPreAndPostOrder() # Needed.
            for pNum in range(self.modelInfo.nParts):
                if self.modelInfo.parts[pNum].nComps > 1:
                    modelKeyHash = {}
                    for n in conTree.iterNodesNoRoot():
                        n.br.textDrawSymbol = var.modelSymbols[n.parts[pNum].compNum]
                        if not modelKeyHash.has_key(n.parts[pNum].compNum):
                            modelKeyHash[n.parts[pNum].compNum] = n.br.textDrawSymbol
                    print
                    conTree.draw()
                    print "\nThe drawing above shows majority comp numbers in partition %i" % pNum
                    kk = modelKeyHash.keys()
                    kk.sort()
                    for k in kk:
                        print "    comp %2i - %s" % (k, modelKeyHash[k])
                    print "    Root has comp %i" % conTree.root.parts[pNum].compNum
                if self.modelInfo.parts[pNum].nRMatrices > 1:
                    modelKeyHash = {}
                    for n in conTree.iterNodesNoRoot():
                        n.br.textDrawSymbol = var.modelSymbols[n.br.parts[pNum].rMatrixNum]
                        if not modelKeyHash.has_key(n.br.parts[pNum].rMatrixNum):
                            modelKeyHash[n.br.parts[pNum].rMatrixNum] = n.br.textDrawSymbol
                    print
                    conTree.draw()
                    print "\nThe drawing above shows majority rMatrix numbers in partition %i" % pNum
                    kk = modelKeyHash.keys()
                    kk.sort()
                    for k in kk:
                        print "    rMatrix %2i - %s" % (k, modelKeyHash[k])
            print

        if showRootInfo:
            # Print out a table showing what nodes or branches had the root.
            if self.isBiRoot:
                if 0:
                    conTree.draw()
                print gm[0]
                print longMessage1 # see top of file
                print "node  br.biRootCount"
                for n in conTree.iterNodesNoRoot():
                    if n.br.biRootCount == None:
                        pass
                    else:
                        if n == biRootChild:
                            # Point out that it has been rooted on that branch.
                            print "%4i  %6.1f <==" % (n.nodeNum, n.br.biRootCount)
                        else:
                            print "%4i  %6.1f" % (n.nodeNum, n.br.biRootCount)
            else:
                print gm[0]
                print longMessage2 # see top of file.
                print "node  rootCount"
                for n in conTree.nodes:
                    if hasattr(n, 'rootCount'):
                        if n.rootCount == None:
                            pass
                        else:
                            if n == conTree.root:
                                print "%4i  %6.1f  <==" % (n.nodeNum, n.rootCount)
                            else:
                                print "%4i  %6.1f" % (n.nodeNum, n.rootCount)


        # Attaching splits to node.br's was only temporary.
        for n in conTree.iterNodesNoRoot():
            if hasattr(n.br, 'split'): # theBiRoot doesnt have one
                del(n.br.split)
            n.br.splitKey = None
        # So was the root rootModelUsage
        if hasattr(conTree.root, 'rootModelUsage'):
            del(conTree.root.rootModelUsage)

        # At the moment, the textDrawSymbols are whatever was left
        # after drawing the model.  Put them back to default.
        for n in conTree.iterNodesNoRoot():
            n.br.textDrawSymbol = '-'

        conTree.preAndPostOrderAreValid = 0
        conTree.taxNames = self.taxNames
        return conTree

    def makeTreeFromPartitions(self, partitions, taxNames=None, zeroBasedNumbering=True):
        """Make a tree from a list of partitions.

        Each partition in the list of partitions is a list of taxon
        numbers (zero-based, or 1-based if the arg zeroBasedNumbering
        is set to False).

        It needs taxNames set, which can be done as an arg in this
        method, or set before.

        """
        
        gm = ['TreePartitions.makeTreeFromPartitions()']
        if type(partitions) != types.ListType:
            gm.append('The partitions should be a list of lists of ints.')
            raise Glitch, gm

        for it in partitions:
            if type(it) != types.ListType:
                gm.append('The partitions should be a list of lists of ints.')
                raise Glitch, gm
            for it2 in it:
                if type(it2) != types.IntType:
                    gm.append('The partitions should be a list of lists of ints.')
                    raise Glitch, gm

        if zeroBasedNumbering:
            zbParts = partitions
        else:
            zbParts = []
            for p in partitions:
                cp = p[:]
                for i in range(len(cp)):
                    cp[i] -= 1
                zbParts.append(cp)
        #print zbParts
        if taxNames:
            self.taxNames = taxNames
        self.nTax = len(self.taxNames)
        self.nTrees = 1
        allOnes = 2L**(self.nTax) - 1
        txx = []
        for tNum in range(self.nTax):
            tName = self.taxNames[tNum]
            tx = TP_TinyTax(tName)
            tx.rawSplitKey = 1L << tNum
            txx.append(tx)
        for tx in txx:
            spl = Split()
            self.splits.append(spl)
            spl.rawSplitKey = tx.rawSplitKey
            spl.count = 1.0
    
        for zbPart in zbParts:
            rsk = 0L
            for it in zbPart:
                rsk = rsk | txx[it].rawSplitKey
            spl = Split()
            self.splits.append(spl)
            spl.rawSplitKey = rsk
            spl.count = 1.0

        for spl in self.splits:
            if (1 & spl.rawSplitKey):
                spl.key = allOnes ^ spl.rawSplitKey
            else:
                spl.key = spl.rawSplitKey

        toRemoveFromSplits = []
        for spl in self.splits:
            if not self.splitsHash.get(spl.key):
                self.splitsHash[spl.key] = spl
            else:
                toRemoveFromSplits.append(spl)
                self.biSplits.append(spl)
                self.biSplitsHash[spl.key] = spl
        
        if toRemoveFromSplits:
            self.isBiRoot = True
        for spl in toRemoveFromSplits:
            self.splits.remove(spl)

        # print '------'
        # for spl in self.splits:
        #     print "%3i %3i" % (spl.rawSplitKey, spl.key)
        # print "biSplits ====="
        # for spl in self.biSplits:
        #     print "%3i %3i" % (spl.rawSplitKey, spl.key)
    
        if self.isBiRoot and len(self.biSplits) != 1:
            gm.append('It appears to be a bifurcationg root.')
            gm.append("But I got %i 'biSplits' -- should be 1." % len(self.biSplits))
            gm.append("Possibly a badly-formed matrix?")
            gm.append("Possibly a duplicated site?")
            raise Glitch, gm

        self._finishSplits()
        #self.dump()

        t = self.consensus()
        t.stripBrLens()
        t.name = 'treeFromPartitions'
        #t.write()
        return t


    def combineWith(self, otherTreePartitionsObject):
        """Combine the info from another TreePartitions object into self.

        Hack alert!

        Only for topology info, not for model info.  Attribute
        doModelComments from self needs to be set to zero, and
        attribute modelInfo from self needs to be set to None (by you,
        the user) or else it won't work.

        """

        otp = otherTreePartitionsObject
        gm = ["TreePartitions.combineWith()"]
        
        # Checks.
        if self.nTax != otp.nTax:
            gm.append("nTax differs.  self %i, other %i" % (self.nTax, otp.nTax))
            raise Glitch, gm
        for tNum in range(len(self.taxNames)):
            if self.taxNames[tNum] != otp.taxNames[tNum]:
                gm.append("taxNames differ.")
                raise Glitch, gm
        if self.isBiRoot or otp.isBiRoot:
            gm.append("self.isBiRoot = %s, other.isBiRoot = %s" % (self.isBiRoot, otp.isBiRoot))
            gm.append("not tested yet for biRoots")
            raise Glitch, gm
        if self.doModelComments or otp.doModelComments:
            gm.append("self.doModelComments = %s, other.doModelComments = %s" % (self.doModelComments, otp.doModelComments))
            raise Glitch, gm
        if self.modelInfo:
            gm.append("self.modelInfo exists -- set it to None.")
            raise Glitch, gm

        # Combine.
        self.nTrees += otp.nTrees

        for ospl in otp.splits:
            sspl = self.splitsHash.get(ospl.key)
            if sspl:
                sspl.combineWith(ospl)
            else:
                self.splits.append(ospl)
                self.splitsHash[ospl.key] = ospl
        for ospl in otp.biSplits:
            sspl = self.biSplitsHash.get(ospl.key)
            if sspl:
                sspl.combineWith(ospl)
            else:
                self.biSplits.append(ospl)
                self.biSplitsHash[ospl.key] = ospl
        self._finishSplits()

        
    def getSplitForTaxNames(self, txNames):
        k = func.getSplitKeyFromTaxNames(self.taxNames, txNames)
        return self.splitsHash.get(k)
