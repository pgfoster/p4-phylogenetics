import sys,string,types,cStringIO,math,copy
import func
from Var import var
from Glitch import Glitch
if var.usePfAndNumpy:
    import numpy
    from Model import Model
from Node import Node,NodePart,NodeBranchPart
import NexusToken


class Tree(object):
    """A phylogenetic tree.


**Some instance variables**

    * ``fName``, if the Tree was read from a file
    * ``name``, the name of the tree.
    * ``root``, the root node
    * ``nodes``, a list of nodes.
    * ``preOrder`` and ``postOrder``, lists of node numbers
    * ``recipWeight``, the weight, if it exists, is usually 1/something, so the reciprocal looks nicer ...
    * ``nexusSets``, if it exists, a NexusSets object.

**Properties**
    
    * ``taxNames``, a list of names.  Usually the order is important!
    * ``data``, a :class:`Data.Data` object
    * ``model``, a :class:`Model.Model` object
    * ``nTax``, the number of taxa
    * ``nInternalNodes``, the number of non-leaf nodes

**The node method**

    You often will want to refer to a node in a tree.  This can be
    done via its name, or its nodeNum, or as an object, via the method
    :meth:`Tree.Tree.node`.  For example, from a Tree ``t`` you can
    get node number 3 via::

        n = t.node(3)

    and you can get the node that is the parent of the Mastodon via::

        n = t.node('Mastodon').parent

    For many methods that require specifying a node, the method argument is *nodeSpecifier*, eg::

        t.reRoot(23)

    ``reRoots``'s the tree to node number 23.

**Describe, draw, and get information about the tree**

     .. autosummary::
     
        Tree.dump
        Tree.draw
        Tree.textDrawList
        Tree.tv
        Tree.btv
        Tree.isFullyBifurcating
        Tree.taxSetIsASplit
        Tree.getAllLeafNames
        Tree.getChildrenNums
        Tree.getDegree
        Tree.getLen
        Tree.getNodeNumsAbove
        Tree.getPreAndPostOrderAbove
        Tree.getPreAndPostOrderAboveRoot
        Tree.getSeqNumsAbove
        Tree.subTreeIsFullyBifurcating
        Tree.summarizeModelThingsNNodes
        Tree.verifyIdentityWith

**Write**

    .. autosummary::

       Tree.write
       Tree.writeNewick
       Tree.writeNexus
       Tree.writePhylip
       Tree.tPickle
        
    See also Trees methods :meth:`Trees.Trees.writeNexus` and
    :meth:`Trees.Trees.writeNewick` for doing trees by the bunch.

**Iteration over the nodes**

    Sometimes you don't want to just iterate over the self.nodes list,
    because after some manipulations a node might be in self.nodes but
    not actually in the tree; using these 'iter' methods takes care of
    that, skipping such nodes.

    .. autosummary::

       Tree.iterInternals
       Tree.iterInternalsNoRoot
       Tree.iterInternalsNoRootPostOrder
       Tree.iterInternalsNoRootPreOrder
       Tree.iterInternalsPostOrder
       Tree.iterLeavesNoRoot
       Tree.iterLeavesPostOrder
       Tree.iterLeavesPreOrder
       Tree.iterNodes
       Tree.iterNodesNoRoot
       Tree.iterPostOrder
       Tree.iterPreOrder
       Tree.nextNode

    See also Node methods that do similar things starting from a given node.
    
**Copy**

     .. autosummary::
     
        Tree.dupe
        Tree.copyToTree
        Tree.dupeSubTree

**In combination with Data and Model**

     .. autosummary::

        Tree.calcLogLike
        Tree.optLogLike
        Tree.simulate
        Tree.getSiteLikes
        Tree.getSiteRates
        Tree.bigXSquaredSubM
        Tree.compStatFromCharFreqs
        Tree.compoTestUsingSimulations
        Tree.modelFitTests
        Tree.modelSanityCheck
        Tree.simsForModelFitTests

**Setting a model**

     .. autosummary::

        Tree.newComp
        Tree.newRMatrix
        Tree.newGdasrv
        Tree.setPInvar
        Tree.setRelRate
        Tree.setModelThing
        Tree.setModelThingsRandomly
        Tree.setModelThingsNNodes
        Tree.summarizeModelThingsNNodes
        Tree.setNGammaCat
        Tree.setTextDrawSymbol


**Tree manipulation**

     .. autosummary::

        Tree.addLeaf
        Tree.addNodeBetweenNodes
        Tree.addSibLeaf
        Tree.addSubTree
        Tree.allBiRootedTrees
        Tree.collapseNode
        Tree.ladderize
        Tree.lineUpLeaves
        Tree.nni
        Tree.pruneSubTreeWithoutParent
        Tree.randomSpr
        Tree.randomizeTopology
        Tree.reRoot
        Tree.reconnectSubTreeWithoutParent
        Tree.removeEverythingExceptCladeAtNode
        Tree.removeNode
        Tree.removeAboveNode
        Tree.removeRoot
        Tree.renameForPhylip
        Tree.restoreDupeTaxa
        Tree.restoreNamesFromRenameForPhylip
        Tree.rotateAround
        Tree.spr
        Tree.stripBrLens

**Misc**

     .. autosummary::
     
        Tree.checkDupedTaxonNames
        Tree.checkSplitKeys
        Tree.checkTaxNames
        Tree.checkThatAllSelfNodesAreInTheTree
        Tree.inputTreesToSuperTreeDistances
        Tree.makeSplitKeys
        Tree.readBipartitionsFromPaupLogFile
        Tree.recalculateSplitKeysOfNodeFromChildren
        Tree.setNexusSets
        Tree.topologyDistance
        Tree.tvTopologyCompare
        Tree.patristicDistanceMatrix
        Tree.inputTreesToSuperTreeDistances


    
    """
    # Tree methods in other files.
    from Tree_muck import node, rotateAround, reRoot, removeRoot, removeNode, removeAboveNode, collapseNode, pruneSubTreeWithoutParent, reconnectSubTreeWithoutParent, addNodeBetweenNodes, allBiRootedTrees, ladderize, randomizeTopology, readBipartitionsFromPaupLogFile, renameForPhylip, restoreNamesFromRenameForPhylip, restoreDupeTaxa, lineUpLeaves, removeEverythingExceptCladeAtNode, dupeSubTree, addSubTree, addLeaf, addSibLeaf, subTreeIsFullyBifurcating, nni, checkThatAllSelfNodesAreInTheTree, spr, randomSpr, inputTreesToSuperTreeDistances
    from Tree_write import patristicDistanceMatrix, tPickle, writePhylip, writeNexus, write, writeNewick, _getMcmcCommandComment, draw, textDrawList, eps, svg
    if var.usePfAndNumpy:
        from Tree_model import data, model, _setData, _checkModelThing, newComp, newRMatrix, newGdasrv, setPInvar, setRelRate, setRjComp, setRjRMatrix, setModelThing, setModelThingsRandomly, setModelThingsNNodes, summarizeModelThingsNNodes, setTextDrawSymbol, setNGammaCat, modelSanityCheck, setEmpiricalComps
        from Tree_optSim import __del__, deleteCStuff, allocCStuff, setCStuff, _commonCStuff, calcLogLike, optLogLike, optTest, simplexDump, simulate, getSiteLikes, getSiteRates
        from Tree_fit import simsForModelFitTests, modelFitTests, compoTestUsingSimulations, bigXSquaredSubM, compStatFromCharFreqs
    #from Tree_pyLike import allocBigP, allocPartPatternsAndCondLikes, calcPartPatterns, calcBigP, calcBigPPart, calcCondLikes, calcCondLikesForPart, calcCondLikesForPartAtNode,partLogLike, calcLogLike2

    def __init__(self):
        self.fName = None   # The name of the file it came from
        self.name = None
        self.root = None
        self.nodes = []
        self.preOrder = None                # nodeNums of nodes root -> tips   A numpy array
        self.postOrder = None               # nodeNums of nodes tips -> root   A numpy array
        self.preAndPostOrderAreValid = 0
        self.recipWeight = None           # Usually weight is 1/N, so the reciprocal looks nicer
        #self.weight = None               # Only for floating point weights, so not usually ...
        self._taxNames = []               # An ordered list.  self.taxNames is a property
        self._data = None                 # A Data object.  self.data is a property
        self.cTree = None                 # A pointer to a c-struct
        self.logLike = None
        self.partLikes = None
        self._model = None                # A Model object.  self.model is a property
        # self.nTax, a property
        self._nTax = 0
        self._nInternalNodes = -1
        # self.nInternalNodes, a property
        self.doDataPart = 0
        self.nexusSets = None
        self.nodeForSplitKeyDict = None


    #########################################################
    # Properties: data, model in Tree_model
    #########################################################
    

    #########################################################
    # Properties: taxNames, nTax, nInternalNodes
    #########################################################
    

    def _setTaxNames(self, theTaxNames):
        gm = ['Tree._setTaxNames()']
        if type(theTaxNames) != type([]):
            gm.append("You can only set property 'taxNames' to a list.")
            gm.append("Got attempt to set to '%s'" % theTaxNames)
            raise Glitch, gm
        self._taxNames = theTaxNames
        if theTaxNames: # and not var.allowTreesWithDifferingTaxonSets:  # Peter commented out until it is sorted. Why is it here?
            self.checkTaxNames()
    
    def _setTaxNamesFromLeaves(self):
        tax = []
        for n in self.iterNodes():
            if n.isLeaf and n.name:
                tax.append(n.name)
            # This next line should not be needed, as root leaves should be leaves.
            elif n == self.root and n.name and n.getNChildren() < 2:  # terminal root that has a taxName
                tax.append(n.name)
            tax.sort()
        self._taxNames = tax
        if self._taxNames:
            self.checkTaxNames()
    
    def _delTaxNames(self):
        gm = ['Tree._delTaxNames()']
        gm.append("    Caught an attempt to delete self.taxNames, but")
        gm.append("self.taxNames is a property, so you can't delete it.")
        gm.append("But you can set it to an empty list if you like.")
        raise Glitch, gm

    taxNames = property(lambda self: self._taxNames, _setTaxNames, _delTaxNames)

    def _getNTax(self):
        # We can't rely on len(self.taxNames), cuz it might not exist.
        #if hasattr(self, '_nTax') and self._nTax:
        if self._nTax:
            return self._nTax
        else:
            nTax = 0
            if self.nodes:
                for n in self.iterNodes():
                    if n.isLeaf:
                        nTax += 1
            self._nTax = nTax
            return nTax
    
    nTax = property(_getNTax) 

    def _getNInternalNodes(self):
        if self._nInternalNodes >= 0:
            return self._nInternalNodes
        elif not self.nodes:
            return 0
        else:
            self._nInternalNodes = len([n for n in self.iterInternalsNoRoot()])
            if not self.root.isLeaf:
                self._nInternalNodes += 1
            return self._nInternalNodes

    def _setNInternalNodes(self, theNInternalNodes):
        gm = ['Tree._setNInternalNodes()']
        gm.append("Caught an attempt to set self.nInternalNodes, but")
        gm.append("self.nInternalNodes is a property, so you shouldn't do that.")
        raise Glitch, gm

    def _delNInternalNodes(self):
        self._nInternalNodes = -1

    nInternalNodes = property(_getNInternalNodes, _setNInternalNodes, _delNInternalNodes)



    ##################################################
    ##################################################

    def dupe(self):
        """Duplicates self, but with no c-pointers.  And no data object.

        If there is a model, it is duped.
        """

        if var.usePfAndNumpy:
            storedData = None
            if self.data:
                storedData = self.data
                self.data = None  # We don't want to have to copy a big data object, now do we?
        import copy
        dupe = copy.deepcopy(self)
        #print 'Tree.dupe()   self.root=%s, dupe.root=%s' % (self.root, dupe.root) 

        # Delete cPointers
        for n in dupe.nodes:
            if n.cNode:
                n.cNode = None
        if dupe.cTree:
            dupe.cTree = None
        if var.usePfAndNumpy:
            if dupe.model and dupe.model.cModel:
                dupe.model.cModel = None

        if var.usePfAndNumpy:
            if storedData:
                self.data = storedData

        return dupe


        


    def parseNexus(self, flob, translationHash=None, doModelComments=0):
        """Start parsing nexus format newick tree description.

        From just after the command word 'tree', to the first paren of
        the Newick part of the tree."""

        gm  = ['Tree.parseNexus()'] # re-defined below
        if 0:
            print 'Tree.parseNexus() translationHash = %s' % translationHash
            print '    doModelComments = %s (nParts)' % doModelComments
            print '    var.nexus_doFastNextTok = %s' % var.nexus_doFastNextTok

        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok,safeNextTok
        else:
            from NexusToken import nextTok,safeNextTok



        tok = safeNextTok(flob, 'Tree.parseNexus()')
        #print 'parseNexus() tok = %s' % tok
        tok = func.nexusUnquoteName(tok)
        if tok == '*':
            print gm[0]
            print "    Ignoring '*' in tree description"
            tok = safeNextTok(flob, 'Tree.parseNexus()')
        if not func.nexusCheckName(tok):
            gm.append("Bad tree name: '%s'" % tok)
            raise Glitch, gm
        self.name = tok
        #print "got name: '%s'" % tok
        #print "%s" % tok
        gm = ["Tree.parseNexus() '%s'" % self.name] # re-defining
        tok = safeNextTok(flob, gm[0])
        if tok != '=':
            gm.append("Tree name must be followed by '='")
            raise Glitch, gm

        # Generally this is the beginning of the newick tree
        # description.  But we have to look ahead to see if there is a
        # weight comment.        
        savedPos = flob.tell()
        while 1:
            tok = safeNextTok(flob, gm[0])
            #print "parseNexus: tok after '=' is '%s'" % tok

            # This next bit will only happen if either var.nexus_getWeightCommandComments
            # or var nexus_getAllCommandComments is set.
            if tok[0] == '[':
                self.getWeightCommandComment(tok)
            elif tok == '(':
                flob.seek(-1, 1)
                self.parseNewick(flob, translationHash, doModelComments)
                #self.initFinish()
                break
            elif tok == ';':
                gm.append("Got ';' before any tree description.")
                raise Glitch, gm
            elif tok[0] in string.letters + string.digits + '_' + '\'':
                flob.seek(savedPos, 0)
                self.parseNewick(flob, translationHash, doModelComments)
                #self.initFinish()
                break
            else:
                gm.append('Expecting a newick tree description.')
                raise Glitch, gm
        self.initFinish()
        #print 'finished Tree.parseNexus()'

    def getWeightCommandComment(self, tok):
        if 0:
            print 'var.nexus_getWeightCommandComments = %s' % var.nexus_getWeightCommandComments
            print 'var.nexus_getAllCommandComments = %s' % var.nexus_getAllCommandComments
            print "Got comment '%s', checking if it is a 'weight' comment." % tok
        gm = ["Tree.getWeightCommandComment()"]
        from NexusToken import nextTok,safeNextTok  # python, not c, so I can use StringIO
        cFlob = cStringIO.StringIO(tok)
        cFlob.seek(1) # The [
        cTok = nextTok(cFlob)
        if not cTok:
            #print "no cTok -- returning nothing"
            return
        lowCTok = string.lower(cTok)
        if lowCTok in ['&r', '&u']:
            #print "got %s -- returning nothing" % cTok
            return
        if lowCTok != '&w':
            gm.append('Expecting a weight comment.  Got %s' % tok)
            raise Glitch, gm
        cTok = nextTok(cFlob)
        # It might be a float, or the more usual 1/something
        if ("." in cTok):
            # A float?
            try:
                self.weight = float(cTok)
            except:
                gm.append("I can't grok '%s' in weight comment %s" % (cTok, tok))
                raise Glitch, gm

        # Should check for scientific notation?

        else:
            try:
                theNumerator = int(cTok)
                if theNumerator != 1:
                    gm.append('Expecting a numerator 1 in weight comment %s' % tok)
                    raise Glitch, gm
                #print 'got theNumerator %i' % theNumerator
            except ValueError:
                gm.append('Expecting a numerator 1 in weight comment %s' % tok)
                raise Glitch, gm
            cTok = nextTok(cFlob)
            if cTok == '/':
                cTok = safeNextTok(cFlob, 'Getting weight comment %s' % tok)
                try:
                    self.recipWeight = int(cTok)
                except ValueError:
                    gm.append('Bad denominator in weight comment %s' % tok)
                    raise Glitch, gm
            elif cTok == ']':
                #self.recipWeight = theNumerator # ie 1, might as well leave it as None
                pass
            else:
                gm.append("I can't grok '%s' in weight comment %s" % (cTok, tok))
                raise Glitch, gm
        cFlob.close()
        #print 'got recipWeight = %s' % self.recipWeight
        




##    ##Ignore
##    def printStack(self, theStack):  # only used for debugging parseNewick()
##        print 'stack = ',
##        for n in theStack:
##            print "%i['%s'] " % (n.nodeNum, n.name),
##        print ''

    def parseNewick(self, flob, translationHash, doModelComments=0):
        """Parse Newick tree descriptions.

        This is stack-based, and does not use recursion.
        """

        #print 'parseNewick here. var.nexus_doFastNextTok=%s' % var.nexus_doFastNextTok
        #print 'parseNewick here. doModelComments=%s' % doModelComments
        #print "parseNewick()  translationHash=%s, self.taxNames=%s" % (translationHash, self.taxNames)

        if self.name:
            gm = ["Tree.parseNewick(), tree '%s'" % self.name]
        else:
            gm = ['Tree.parseNewick()']

        if hasattr(flob, 'name') and flob.name:
            self.fName = flob.name
            gm[0] += ", File %s" % self.fName

        if doModelComments:
            savedP4Nexus_getAllCommandComments = var.nexus_getAllCommandComments # restore at end
            var.nexus_getAllCommandComments = 1

        stack = []
        isAfterParen = 1 # to start, even tho its not true
        isAfterComma = 0
        parenNestLevel = 0
        lastPopped = None

        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok,safeNextTok
        else:
            from NexusToken import nextTok,safeNextTok

        tok = nextTok(flob)
        if not tok:
            return
        
        tok = func.nexusUnquoteName(tok)  # Should generally be the opening paren, except if its a single-node tree.
        while tok != ';':
            #print "top of loop tok '%s', tok[0] is '%s'" % (tok, tok[0])
            
            if tok == '(':
                #print "Got '(': new node (%i)." % len(self.nodes)
                if not (isAfterParen or isAfterComma):
                    gm.append('Got badly-placed paren, not after a paren or comma.')
                    raise Glitch, gm
                newNode = Node()
                if doModelComments:
                    for pNum in range(doModelComments):
                        newNode.parts.append(NodePart())
                        newNode.br.parts.append(NodeBranchPart())
                        
                #self.printStack(stack)
                if len(stack):
                    newNode.parent = stack[-1]
                    if newNode.parent.leftChild == None:
                        newNode.parent.leftChild = newNode
                    else:
                        newNode.parent.rightmostChild().sibling = newNode
                else:
                    if len(self.nodes) == 0:
                        self.root = newNode
                        newNode.isLeaf = 1  # Sometimes. Generally not true-- corrected at the end.
                    else:
                        gm.append('Something is wrong. Stack is empty.')
                        gm.append('Extra paren?')
                        raise Glitch, gm
                newNode.nodeNum = len(self.nodes)
                self.nodes.append(newNode)
                stack.append(newNode)
                isAfterParen = 1
                parenNestLevel += 1

            elif tok == ',':
                if isAfterParen:
                    gm.append('Got comma after paren.')
                    raise Glitch, gm
                elif isAfterComma:
                    gm.append('Got comma after comma.')
                    raise Glitch, gm
                #self.printStack(stack)
                try:
                    lastPopped = stack.pop()
                except IndexError:
                    gm.append('Empty stack.  Out of place comma?')
                    raise Glitch, gm
                isAfterComma = 1
                if len(stack) == 0:
                    gm.append('Empty stack.  Out of place comma?')
                    raise Glitch, gm

            elif tok == ')':
                try:
                    lastPopped = stack.pop()
                except IndexError:
                    gm.append('Empty stack.  Out of place unparen?')
                    raise Glitch, gm

                isAfterParen = 0
                isAfterComma = 0
                parenNestLevel = parenNestLevel - 1
                if parenNestLevel < 0:
                    gm.append('Unmatched unparen.')
                    raise Glitch, gm
                if len(stack) == 0 and len(self.nodes) > 1:
                    gm.append('Empty stack.  Out of place unparen?')
                    raise Glitch, gm

            elif tok[0] in string.letters or tok[0] in string.digits or tok[0] == "'" or tok[0] in [
                '_', '#', '\\', '/', '"', '(', ')']:
                if len(self.nodes) == 0: # A single-node tree, not ()aName, rather just aName.
                    isAfterParen = 1
                if not (isAfterParen or isAfterComma):
                    # Probably a name of an internal node.
                    if len(stack):
                        #if stack[-1].isLeaf and stack[-1].name != '(':
                        if stack[-1].name:
                            if not var.newick_allowSpacesInNames:
                                # a second name after a node name, eg (A foo, B)   =>foo is bad
                                # or eg (A, B)foo bar    => bar is bad
                                gm.append("Badly placed token '%s'." % tok)
                                gm.append("Appears to be a second node name, after '%s'" % stack[-1].name)
                                gm.append('Missing comma maybe?  Or punctuation or spaces in an unquoted name?')
                                gm.append("To allow reading Newick (or Nexus) with spaces, ")
                                gm.append("turn var.newick_allowSpacesInNames on")
                                raise Glitch, gm
                            else:
                                stack[-1].name += ' '
                                stack[-1].name += tok
                        else:
                            # Usually this...
                            #print "naming node %i as '%s'" % (stack[-1].nodeNum, tok)
                            # We allow bad names on internal nodes, ie we do not nexusCheckName(tok)
                            stack[-1].name = tok

                    else:    # len(stack) == 0
                        if lastPopped and lastPopped.name == None: # ()A
                            #print "naming lastPopped node %i with '%s'" % (lastPopped.nodeNum, tok)
                            lastPopped.isLeaf = 1
                            #lastPopped.label = tok
                            lastPopped.name = tok
                        else:
                            gm.append("Badly placed token '%s' in tree description." % tok)
                            raise Glitch, gm


                else:
                    # A new terminal node.
                    if tok[0] in string.letters or tok[0] in ['_']:
                        if translationHash and translationHash.has_key(tok):
                            #print 'got key %s, val is %s' % (tok, translationHash[tok])
                            tok = translationHash[tok]

                    elif tok[0] in string.digits:
                        if var.nexus_allowAllDigitNames:
                            if translationHash and translationHash.has_key(tok):
                                #print 'got key %s, val is %s' % (tok, translationHash[tok])
                                tok = translationHash[tok]
                        else:
                            try:
                                tok = int(tok)
                                if translationHash and translationHash.has_key(`tok`):
                                    tok = translationHash[`tok`]
                                elif translationHash and not translationHash.has_key(`tok`):
                                    gm.append("There is a 'translation' for this tree, but the")
                                    gm.append("number '%i' in the tree description" % tok)
                                    gm.append('is not included in that translate command.')
                                    raise Glitch, gm
                                elif self.taxNames:
                                    try:
                                        tok = self.taxNames[tok - 1]
                                    except IndexError:
                                        gm.append("Can't make sense out of token '%s' for a new terminal node." % tok)
                                        gm.append('There is no translate command, and the taxNames does not')
                                        gm.append('have a value for that number.')
                                        raise Glitch, gm
                                else:
                                    gm.append("We have a taxon name '%s', composed only of numerals." % tok)
                                    gm.append(" ")
                                    gm.append('The Nexus format allows tree specifications with no')
                                    gm.append('translate command to use integers to refer to taxa.')
                                    gm.append('That is possible because in a proper Nexus file the')
                                    gm.append('taxa are defined before the trees.  P4, however, does')
                                    gm.append('not require definition of taxa before the trees, and in')
                                    gm.append('the present case no definition was made.  Deal with it.')
                                    raise Glitch, gm
                            except ValueError:
                                if translationHash and translationHash.has_key(`tok`):
                                    tok = translationHash[`tok`]
                                #else:  # starts with a digit, but it is not an int.
                                #    gm.append('Problem token %s' % tok)
                                #    raise Glitch, gm

                    #print "Got terminal node '%s'" % tok
                    newNode = Node()
                    if doModelComments:
                        for pNum in range(doModelComments):
                            newNode.parts.append(NodePart())
                            newNode.br.parts.append(NodeBranchPart())

                    newNode.isLeaf = 1
                    if func.nexusCheckName(tok):
                        newNode.name = tok
                        #print 'got newNode.name = %s' % tok
                    else:
                        gm.append("Bad name '%s'" % tok)
                        raise Glitch, gm
                    if len(stack):
                        newNode.parent = stack[-1]
                        if newNode.parent.leftChild == None:
                            newNode.parent.leftChild = newNode
                        else:
                            newNode.parent.rightmostChild().sibling = newNode
                    newNode.nodeNum = len(self.nodes)
                    if len(self.nodes) == 0:
                        self.root = newNode
                    self.nodes.append(newNode)
                    stack.append(newNode)
                    isAfterParen = 0
                    isAfterComma = 0


            elif tok == ':':
                #  Looking for a br.len number, which might be eg 0.234 or -1.23e-05
                #  It might be a multi-token operation.  Accumulate tok's in theNum
                theNum = nextTok(flob)
                if not theNum:
                    gm.append('Tree description ended with a colon.  Bad!')
                    raise Glitch, gm
                #print "  Got token after colon:  '%s'" % theNum
                if theNum == '-' or theNum == '+':
                    tok = nextTok(flob)
                    #print "  Got tok: '%s' after '%s'" % (tok, theNum)
                    if not tok:
                        gm.append('Trying to deal with a branch length.')
                        gm.append("It didn't work, tho.")
                        gm.append("Got this after colon:  '%s'" % theNum)
                        gm.append('followed by nothing.')
                        raise Glitch, gm
                    theNum += tok
                try:
                    # If it is a simple number like 0.123 or -23, then we are finished.
                    stack[-1].br.len = float(theNum) # Won't work if it ends in 'e'
                    #print '  Successfully got br.len %f' % stack[-1].br.len
                except ValueError:
                    # The first bit after the colon is hopefully something like +1.23e
                    if theNum[-1] not in ['e', 'E']:
                        gm.append('Trying to deal with a branch length after a colon, but I am totally confused.')
                        gm.append("Can't make sense out of '%s'" % theNum)
                        raise Glitch(gm, 'newick_badBranchLength')
                    try:
                        float(theNum[:-1])
                    except ValueError:
                        gm.append('Trying to deal with a branch length after a colon, but I am totally confused.')
                        gm.append("Can't make sense out of '%s'" % theNum)
                        raise Glitch(gm, 'newick_badBranchLength')

                    # Now we are sure that the first bit *is* something like +1.23e
                    # We do not allow spaces after the 'e', so we do not use nextTok().
                    # That introduces a bug, where comments inserted in the number don't get ignored.   <<== unfixed bug!
                    # The first thing must be a '+' or a '-'.
                    c = flob.read(1)
                    if not c:
                        gm.append('Trying to deal with a branch length, possibly in scientific notation.')
                        gm.append("Got '%s' after the colon, but then nothing." % theNum)
                        raise Glitch, gm
                    if c not in ['+', '-']:
                        gm.append('Trying to deal with a branch length, possibly in scientific notation.')
                        gm.append("Got '%s' after the colon." % theNum)
                        gm.append("Expecting a '+' or '-' after that (no spaces allowed).")
                        gm.append("Got '%s'." % c)
                        raise Glitch, gm
                    # Accumulate characters in 'theExp'.  We need at least one digit.
                    theExp = c
                    c = flob.read(1)
                    if not c:
                        gm.append('Trying to deal with a branch length, possibly in scientific notation.')
                        gm.append("Got '%s%s' after the colon, but then nothing." % (theNum, theExp))
                        raise Glitch, gm
                    if c not in string.digits:
                        gm.append("Trying to deal with a branch length, possibly in scientific notation.")
                        gm.append("Got '%s%s' after the colon." % (theNum, theExp))
                        gm.append('Expecting one or more digits.')
                        gm.append("Got '%s'" % c)
                        raise Glitch, gm
                    theExp += c
                    # So we got one good digit.  Are there any more?
                    while 1:
                        c = flob.read(1)
                        if not c:
                            gm.append('Trying to deal with a branch length, possibly in scientific notation.')
                            gm.append("Got '%s%s' after the colon, but then nothing." % (theNum, theExp))
                            raise Glitch, gm
                        # We got something.  If its a digit, add it to
                        # theExp.  If its anything else, back up one
                        # space and then break
                        if c in string.digits:
                            theExp += c
                        else:
                            flob.seek(-1, 1)
                            break
                    #print "  At this point, theNum='%s' and theExp='%s'" % (theNum, theExp)
                    try:
                        #print "  Trying to see if theExp '%s' can be converted to an int." % theExp
                        int(theExp)
                        try:
                            theBrLen = float(theNum + theExp)
                            #print '  Successfully got br.len %g (from %s%s)' % (theBrLen, theNum, theExp)
                            stack[-1].br.len = theBrLen
                        except ValueError:
                            gm.append('Trying to deal with a branch length, possibly in scientific notation.')
                            gm.append("It didn't work, tho.")
                            gm.append("Got these after colon: '%s' and '%s'" % (theNum, theExp))
                            gm.append('And they could not be converted to an exponential float.')
                            raise Glitch, gm
                    except ValueError:
                        gm.append('Trying to deal with a branch length, possibly in scientific notation.')
                        gm.append("It didn't work, tho.")
                        gm.append("Got these after colon: '%s' and '%s'." % (theNum, theExp))
                        gm.append('And the latter does not appear to be an int.')
                        raise Glitch, gm

            elif tok[0] == '[':
                #print "openSquareBracket.  Got tok '%s'" % tok
                # if doModelComments is set, it should be set to nParts.
                if doModelComments:
                    # eg [& c0.1 r0.0]
                    n = stack[-1]
                    #print 'got comment %s, node %i' % (tok, n.nodeNum)
                    cFlob = cStringIO.StringIO(tok)
                    cFlob.seek(2)
                    tok2 = NexusToken.safeNextTok(cFlob)
                    while 1:
                        if tok2 == ']':
                            break
                        elif tok2[0] in ['c', 'r', 'g']:
                            ending = tok2[1:]
                            splitEnding = string.split(ending, '.')
                            try:
                                firstNum = int(splitEnding[0])
                                secondNum = int(splitEnding[1])
                            except ValueError:
                                gm.append('Bad command comment %s' % tok)
                                raise Glitch, gm 
                            if tok2[0] == 'c':
                                n.parts[firstNum].compNum = secondNum
                            if tok2[0] == 'r':
                                n.br.parts[firstNum].rMatrixNum = secondNum
                            if tok2[0] == 'g':
                                n.br.parts[firstNum].gdasrvNum = secondNum
                        else:
                            gm.append('Bad command comment %s' % tok)
                            raise Glitch, gm
                        tok2 = NexusToken.safeNextTok(cFlob)
                elif 0:
                    # Ugly hack for RAxML trees with bootstrap
                    # supports in square brackets after the br len, on
                    # internal nodes.  First modify the eg [100] to be
                    # [&100], set var.nexus_getAllCommandComments =
                    # True, and turn on this elif section.
                    myNode = stack[-1]
                    assert not myNode.isLeaf
                    assert not myNode.name
                    mySupportString = tok[2:-1]
                    #print mySupportString
                    myNode.name = mySupportString
                elif var.nexus_readBeastTreeCommandComments:
                    n = stack[-1]
                    i = 2
                    while i < (len(tok) - 1):
                        j = i
                        inBraces = False
                        while 1:
                            j += 1
                            #print tok[j]
                            if tok[j] == ']':
                                break
                            if tok[j] == '{':
                                inBraces = True
                            if inBraces:
                                if tok[j] == '}':
                                    inBraces = False
                            else:
                                if tok[j] == ',':
                                    break
                        substring = tok[i:j]
                        #print substring
                        splSS = substring.split('=')
                        theNameString = splSS[0].strip()
                        if '%' in theNameString:
                            theNameString = theNameString.replace('%', '')
                        theValString = splSS[1]
                        if '{' in theValString:
                            theValString = theValString.replace('{', '(')
                            theValString = theValString.replace('}', ')')
                        if theValString == 'true':
                            theVal = True
                        elif theValString == 'false':
                            theVal = False
                        else:
                            theVal = eval(theValString)
                        assert type(theVal) in [types.FloatType, types.TupleType, types.BooleanType]
                        n.__setattr__(theNameString, theVal) 
                        i = j + 1
                    

            else:
                gm.append("I can't make sense of the token '%s'" % tok)
                if len(tok) == 1:
                    if tok[0] in var.punctuation:
                        gm.append("The token is in var.punctuation. If you don't think it should")
                        gm.append("be, you can modify what p4 thinks that punctuation is.")
                        if var.nexus_doFastNextTok:
                            gm.append("But to do that you can't use nexus_doFastNextToken, so you ")
                            gm.append("need to turn that off, temporarily.  ")
                            gm.append("(var.nexus_doFastNextTok is currently on.)")
                            gm.append("So you might do this:")
                            gm.append("var.nexus_doFastNextTok = False ")
                            gm.append("var.punctuation = var.phylip_punctuation")
                            gm.append("(or use your own definition -- see Var.py)")
                            gm.append("read('yourWackyTreeFile.phy')")
                            gm.append("That might work.")
                        else:
                            gm.append("(Doing that does not work if nexus_doFastNextToken is turned on,")
                            gm.append("but var.nexus_doFastNextTok is currently off.)")
                            gm.append("So you might do this:")
                            gm.append("var.punctuation = var.phylip_punctuation")
                            gm.append("(or use your own definition -- see Var.py)")
                            gm.append("read('yourWackyTreeFile.phy')")
                            gm.append("That might work.")
                    if tok[0] not in var.punctuation:
                        gm.append("The token is not in your current var.punctuation.")
                        if var.nexus_doFastNextTok:
                            gm.append("var.nexus_doFastNextTok is currently on.")
                            gm.append("It uses a hard-wired list of punctuation, and so you may need to ")
                            gm.append("need to turn var.nexus_doFastNextTok off, temporarily.  ")
                            gm.append("So you might do this:")
                            gm.append("var.nexus_doFastNextTok = False ")
                            gm.append("read('yourWackyTreeFile.phy')")
                            gm.append("That might work.")
                        else:
                            gm.append("var.nexus_doFastNextTok is currently off.")
                            
                            
                #gm.append("tok[0] is '%s'" % tok[0])
                raise Glitch, gm

            tok = func.nexusUnquoteName(safeNextTok(flob, 'Tree init, reading tree string'))
            #print 'got tok for next round = %s' % tok
            # This is the end of the "while tok != ';':" loop

        #print '\n*** Stack len = %i ***' % len(stack)

        if parenNestLevel > 0:
            gm.append('Unmatched paren.')
            raise Glitch, gm
        elif parenNestLevel < 0:
            gm.append('Unmatched unparen.')
            raise Glitch, gm

        if len(stack) == 0:
            if len(self.nodes) == 1:
                pass
            else:
                gm.append("Got an oddly-placed ';' in the tree %s description." % self.name)
                self.dump(treeInfo=0, nodeInfo=1)
                raise Glitch, gm
        elif len(stack) > 1:
            gm.append("Got an oddly-placed ';' in the tree %s description." % self.name)
            #gm.append('(stack len = %i)' % len(stack)
            #self.dump(tree=0, node=1)
            raise Glitch, gm

        if self.root.leftChild and self.root.leftChild.sibling:  # usually this
            self.root.isLeaf = 0

        # Should a root on a stick be a leaf?  If it is just for
        # display purposes, then it should be ok to not be a leaf.
        # But if you are going to re-Root, then it will cause trouble.
        # So by default, a root on a stick should be a leaf.  I think.

        # Hopefully if you are dealing with this, then you know what
        # you are doing and what you want, and how to modify things to
        # get it.

        # Uncomment this next line to make it always non-leaf, even if it is a leaf.
        
        #self.root.isLeaf = 0 # do this to always have a non-leaf root-on-a-stick   <-- Potential Trouble!!!
        self.root.br = None
        #self.draw()
        #self.dump(tree=0, node=1, treeModel=0)

        if doModelComments:
            # restore the value of var.nexus_getAllCommandComments, which was saved above.
            var.nexus_getAllCommandComments = savedP4Nexus_getAllCommandComments



    ##Ignore
    def initFinish(self):

        if self.name:
            gm = ["Tree.initFinish()   tree '%s'" % self.name]
        else:
            gm = ['Tree.initFinish()']

        # Checking for duped taxon names used to be here, but it was
        # pushed further ahead, so that self.taxNames can be corrected
        # also, if need be.  At this point, self does not have
        # taxNames set.

        # Check that all terminal nodes have names
        for item in self.nodes:
            if item.isLeaf:
                #print 'leaf name %s' % item.name
                if not item.name:
                    if item == self.root:
                        if var.warnAboutTerminalRootWithNoName:
                            print 'Tree.initFinish()'
                            print '    Non-fatal warning: the root is terminal, but has no name.'
                            print '    This may be what you wanted.  Or not?'
                            print '    (To get rid of this warning, turn off var.warnAboutTerminalRootWithNoName)'
                    else:
                        gm.append('Got a terminal node with no name.')
                        raise Glitch, gm

        if var.usePfAndNumpy:
            self.preOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
            self.postOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
        else:
            self.preOrder = [var.NO_ORDER] * len(self.nodes)
            self.postOrder = [var.NO_ORDER] * len(self.nodes)
            
        if len(self.nodes) > 1:
            self.setPreAndPostOrder()


    def checkDupedTaxonNames(self):

        # Called by func._tryToReadNexusFile() and func._tryToReadPhylipFile()
        # Check for duped names
        if self.name:
            gm = ["Tree.checkDupedTaxonNames()   tree '%s'" % self.name]
        else:
            gm = ['Tree.checkDupedTaxonNames()']
        if self.fName:
            gm[0] += ' file=%s' % self.fName
        
        hasDupedName = 0
        loNames = []
        for n in self.nodes:
            if n.isLeaf and n.name:
                loNames.append(n.name.lower())
        for loName in loNames:
            if loNames.count(loName) > 1:
                if var.allowDupedTaxonNames:
                    pass
                elif not var.doRepairDupedTaxonNames:
                    gm.append("Got duplicated taxon (lowercased) name '%s'." % loName)
                    gm.append('Since var.doRepairDupedTaxonNames is not turned on, p4 will not fix duplications.')
                    gm.append('To repair duplications verbosely, set ')
                    gm.append('var.doRepairDupedTaxonNames = 1')
                    gm.append('To repair duplications silently, set')
                    gm.append('var.doRepairDupedTaxonNames = 2')
                    raise Glitch, gm
                hasDupedName = 1
                break

        if hasDupedName:
            #print self.name
            if var.allowDupedTaxonNames:
                # more hacking ...
                if var.allowDupedTaxonNames == 2:  # ie silently.
                    pass
                else:
                    complainedAlready = []
                    for loName in loNames:
                        if loNames.count(loName) > 1 and loName not in complainedAlready:
                            if self.name:
                                print "        Tree %s. Duped tax name (lowercased) '%s'" % (
                                    self.name, loName)
                            else:
                                print "        Duped tax name (lowercased) '%s'" % loName
                            complainedAlready.append(loName)

            elif var.doRepairDupedTaxonNames:
                repairedNames = []
                for loName in loNames:
                    if loNames.count(loName) > 1 and n not in repairedNames:
                        repairCounter = 1
                        repairCounter2 = 1
                        for n in self.nodes:
                            if n.isLeaf:
                                if n.name and n.name.lower() == loName:
                                    newName = '%s_%i' % (n.name, repairCounter)
                                    if var.doRepairDupedTaxonNames == 1:
                                        if self.name:
                                            print "        Tree %s. Changing '%s' to '%s'" % (
                                                self.name, n.name, newName)
                                        else:
                                            print "        Changing '%s' to '%s'" % (n.name, newName)
                                    n.name = newName
                                    repairedNames.append(loName)
                                    repairCounter += 1
                        if self.taxNames:
                            for tNameNum in range(len(self.taxNames)):
                                tName = self.taxNames[tNameNum]
                                if tName.lower() == loName:
                                    newName = '%s_%i' % (tName, repairCounter2)
                                    self.taxNames[tNameNum] = newName
                                    repairCounter2 += 1
                            assert repairCounter == repairCounter2, "Got a problem with re-naming duped taxa."
                        

        

    ##############
    ##############
    #
    # dump()
    #
    ##############
    ##############

    def dump(self, tree=0, node=0, model=0, all=0):
        """Print rubbish about self.

        tree
            is the default, showing basic info about the tree.
        node
            shows info about all the nodes.
        model
            shows which modelThing number goes on which node.
            (which you can also get by drawing the tree)

        (If you want the info about the model itself, do a
        aTree.model.dump() instead.)

        """

        if all:
            self._doTreeInfo()
            self._doNodeInfo()
            self._doNodeModelInfo()
        elif not tree and not node and not model:
            self._doTreeInfo()
        else:
            if tree:
                self._doTreeInfo()
            if node:
                self._doNodeInfo()
            if model:
                self._doNodeModelInfo()

    def _doTreeInfo(self):
        if self.name:
            print "Tree '%s' dump" % self.name
        else:
            print 'Tree dump.  No name.'
        if self.fName:
            print "    From file '%s'" % self.fName
        else:
            print "    From an unknown file, or no file."
        if self.root:
            print '    Node %i is root' % self.root.nodeNum
        else:
            print '    There is no root'
        if self.recipWeight:
            print '    The tree recipWeight is %s' % self.recipWeight
        else:
            print '    There is no recipWeight'
        print '    There are %i nodes' % len(self.nodes)
        terminals = 0
        for i in self.nodes:
            if i.isLeaf:
                terminals += 1
        print '        of which %i are terminal nodes' % terminals
        if self.data:
            print '    There is a data object, with %i parts.' % self.data.nParts
        else:
            print '    There is no data object.'

        if self.data:
            print '    The data came from file(s):'
            for a in self.data.alignments:
                if a.fName:
                    print '        %s' % a.fName

        if self.model:
            print '    There is a model object, with %i parts.' % self.model.nParts
            if self.model.cModel:
                print '    model.cModel is %i' % self.model.cModel
            else:
                print '    There is no cModel.'
        else:
            print '    There is no model object.'

        if self.taxNames:
            print '    There is a taxNames list.'
        else:
            print '    There is no taxNames list.'
        if self.cTree:
            print '    cTree is %i' % self.cTree
        else:
            print '    There is no cTree.'


    def _doNodeInfo(self):
        """Basic rubbish about nodes."""

        print '\n-------- nodes -----------------------------------------'
        print '%7s %6s %6s %6s %6s %7s %6s  %4s' % ('nodeNum', 'isLeaf', 'parent', 'leftCh',
                                                    'siblng', 'br.len', 'seqNum', 'name')
        for n in self.nodes:
            print '%7s %6s' % (n.nodeNum, n.isLeaf),
            if n.parent:
                print '%6s' % n.parent.nodeNum,
            else:
                print '%6s' % 'None',

            if n.leftChild:
                print '%6s' % n.leftChild.nodeNum,
            else:
                print '%6s' % 'None',

            if n.sibling:
                print '%6s' % n.sibling.nodeNum,
            else:
                print '%6s' % 'None',

            if n.br and (n.br.len or n.br.len == 0.0):
                print '%7.3f' % n.br.len,
            else:
                print '%7s' % 'None',

            if n.seqNum or n.seqNum == 0:
                print '%6s' % n.seqNum,
            else:
                print '%6s' % 'None',

            if n.name:
                print '  %s' % n.name
            else:
                print '  %s' % 'None'
        print '--------------------------------------------------------\n'


        doMore = 0
        for n in self.iterNodesNoRoot():
            if hasattr(n.br, 'name') and n.br.name:
                doMore = 1
                break
            elif hasattr(n.br, 'uName') and n.br.uName:
                doMore = 1
                break
            elif hasattr(n.br, 'support') and n.br.support:
                doMore = 1
                break
            
        if doMore:
            print '\n-------- more node stuff -------------------------------'
            print '%7s   %10s %10s %10s   %4s' % ('nodeNum', 'br.name', 'br.uName', 'br.support', 'name')
            for n in self.nodes:
                print '%7s' % n.nodeNum,

                if n.br and hasattr(n.br, 'name') and n.br.name:
                    print '%10s' % n.br.name,
                else:
                    print '%10s' % '-',

                if n.br and hasattr(n.br, 'uName') and n.br.uName:
                    print '%10s' % n.br.uName,
                else:
                    print '%10s' % '-',

                if n.br and hasattr(n.br, 'support') and n.br.support:
                    print '%10.4f' % n.br.support,
                else:
                    print '%10s' % '-',

                if n.name:
                    print '    %s' % n.name
                else:
                    print '    %s' % 'None'
            print '--------------------------------------------------------\n'


        doMore = 0
        for n in self.nodes:
            if hasattr(n, 'rootCount') and n.rootCount:
                doMore = 1
                break
            if n.br:
                if hasattr(n.br, 'color') and n.br.color:
                    doMore = 1
                    break
                elif hasattr(n.br, 'biRootCount') and n.br.biRootCount:
                    doMore = 1
                    break

        if doMore:
            print '\n-------- even more node stuff --------------------------'
            print '%7s   %10s %14s %10s   %4s' % ('nodeNum', 'br.color', 'br.biRootCount', 'rootCount', 'name')
            for n in self.nodes:
                print '%7s' % n.nodeNum,

                if n.br and hasattr(n.br, 'color') and n.br.color:
                    print '%10s' % n.br.color,
                else:
                    print '%10s' % '-',

                if n.br and hasattr(n.br, 'biRootCount') and n.br.biRootCount:
                    print '%14s' % n.br.biRootCount,
                else:
                    print '%14s' % '-',

                if hasattr(n, 'rootCount') and n.rootCount:
                    print '%10s' % n.rootCount,
                else:
                    print '%10s' % '-',

                if n.name:
                    print '    %s' % n.name
                else:
                    print '    %s' % 'None'
            print '--------------------------------------------------------\n'









    def _doNodeModelInfo(self):
        if not self.model:
            print '\n****** Node Model Info.  No model.'
            if not self.data:
                print '(no data attached, either)'
        else:
            print '\n****** Node Model Info.  nParts=%s' % self.model.nParts
            if not self.data:
                print 'no data'
            if not self.model.nParts:
                return

            print '\nComps in the nodes:'
            print ' %13s' % 'nodeNum',
            for i in range(self.model.nParts):
                print ' %8s' % 'part%i' % i,
            print ''
            for n in self.nodes:
                print ' %13i' % n.nodeNum,
                for i in range(self.model.nParts):
                    print '%8i' % n.parts[i].compNum,
                print ''

            print '\nrMatrices in the nodes:'
            print ' %13s' % 'nodeNum',
            for i in range(self.model.nParts):
                print ' %8s' % 'part%i' % i,
            print ''
            for n in self.iterNodesNoRoot():
                print ' %13i' % n.nodeNum,
                for i in range(self.model.nParts):
                    print '%8i' % n.br.parts[i].rMatrixNum,
                print ''

            print '\ngdasrvs in the nodes:'
            print ' %13s' % '',
            for i in range(self.model.nParts):
                print ' %8s' % 'part%i' % i,
            print ''
            print ' %13s' % 'nGammaCats ->',
            for i in range(self.model.nParts):
                print '%8i' % self.model.parts[i].nGammaCat,
            print '\n'

            print ' %13s' % 'nodeNum',
            for i in range(self.model.nParts):
                print ' %8s' % 'part%i' % i,
            print ''
            for n in self.iterNodesNoRoot():
                print ' %13i' % n.nodeNum,
                for i in range(self.model.nParts):
                    print '%8i' % n.br.parts[i].gdasrvNum,
                print ''



    ###########################
    #
    # Set a NexusSets object, for taxSets ...
    #
    ###########################


    def setNexusSets(self):
        """Set self.nexusSets from var.nexusSets.

        A deepcopy is made of var.nexusSets, only if it exists.  If
        var.nexusSets does not yet exist, a new blank one is not made
        (cf this method in Alignment class, where it would be
        made).

        Important!  This method depends on a correct taxNames.
        
        """

        assert self.taxNames, "This method requires correct taxNames, in the correct order."
        gm = ["Tree.setNexusSets()"]
        if not var.nexusSets:
            return
        self.nexusSets = copy.deepcopy(var.nexusSets)
        self.nexusSets.taxNames = self.taxNames
        self.nexusSets.nTax = self.nTax
        self.nexusSets.constant = None
        self.nexusSets.gapped = None
        self.nexusSets.charSets = []
        self.nexusSets.charPartitions = []

        if self.nexusSets.taxSets:
            #print "%s. There are %i taxSets." % (gm[0], len(self.nexusSets.taxSets))
            # Check that no taxSet name is a taxName
            lowSelfTaxNames = [string.lower(txName) for txName in self.taxNames]
            for ts in self.nexusSets.taxSets:
                if ts.lowName in lowSelfTaxNames:
                    gm.append("Can't have taxSet names that are the same (case-insensitive) as a tax name")
                    gm.append("Lowercased taxSet name '%s' is the same as a lowcased taxName." % ts.name)
                    raise Glitch, gm
            self.nexusSets.lowTaxNames = lowSelfTaxNames
            
            # If it is standard format, 
            # convert triplets to numberTriplets, and then mask
            for ts in self.nexusSets.taxSets:
                if ts.format == 'standard':
                    ts.setNumberTriplets()
                    ts.setMask()
                    #print ts.mask
                elif ts.format == 'vector':
                    assert ts.mask
                    if len(ts.mask) != self.nTax:
                        gm.append("taxSet %s" % ts.name)
                        gm.append("It is vector format, but the length is wrong.")
                        gm.append("taxSet mask is length %i, but self nTax is %i" % (len(ts.mask), self.nTax))
                        raise Glitch, gm
                else:
                    gm.append("taxSet %s" % ts.name)
                    gm.append("unknown format %s" % ts.format)
                    raise Glitch, gm

                    



    ###########################
    #
    # Get lists of nodeNums ...
    #
    ###########################


    def getNodeNumsAbove(self, nodeSpecifier, leavesOnly=0):
        """Gets a list of nodeNums, in postOrder, above nodeSpecifier.

        The node specified is not included.
        """
        x, y = self.getPreAndPostOrderAbove(nodeSpecifier)
        if leavesOnly:
            tOnly = []
            for i in y[:-1]:
                if self.nodes[i].isLeaf:
                    tOnly.append(i)
            return tOnly
        return y[:-1]


    def getSeqNumsAbove(self, nodeSpecifier):
        """Gets a list of seqNums above nodeSpecifier."""
        x, y = self.getPreAndPostOrderAbove(nodeSpecifier)
        seqNums = []
        for nNum in y[:-1]:
            n = self.nodes[nNum]
            if n.isLeaf:
                seqNums.append(n.seqNum)
        return seqNums


    # def getAllChildrenNums(self, specifier):
    #     """Returns a list of the nodeNums of all children of the specified node
    #     Ambiguous, unused.
    #     """
    #     theNode = self.node(specifier)
    #     if not theNode:
    #         gm = ['Tree.getChildrenNums()']
    #         gm.append('Bad node specifier')
    #         raise Glitch, gm
    #     ret = []
    #     x, y = self.getPreAndPostOrderAbove(specifier)
    #     for i in y[:-1]:
    #         ret.append(self.nodes[i].nodeNum)
    #     return ret
    
    def getAllLeafNames(self, specifier):
        """Returns a list of the leaf names of all children 
        """
        theNode = self.node(specifier)
        if not theNode:
            gm = ['Tree.getChildrenNums()']
            gm.append('Bad node specifier')
            raise Glitch, gm
        ret = []
        x, y = self.getPreAndPostOrderAbove(specifier)
        for i in y[:-1]:
            if self.nodes[i].isLeaf:
                ret.append(self.nodes[i].name)
        return ret

    def getChildrenNums(self, specifier):
        """Returns a list of nodeNums of children of the specified node.

        See also Node.getNChildren()"""

        theNode = self.node(specifier)
        if not theNode:
            gm = ['Tree.getChildrenNums()']
            gm.append('Bad node specifier')
            raise Glitch, gm
        ret = []
        c = theNode.leftChild
        while 1:
            if c:
                ret.append(c.nodeNum)
                c = c.sibling
            else:
                return ret



    def setPreAndPostOrder(self):
        """Sets or re-sets self.preOrder and self.postOrder lists of node numbers.

        PreOrder starts from the root and goes to the tips; postOrder
        starts from the tips and goes to the root."""
        
        self.getPreAndPostOrderAboveRoot()
        self.preAndPostOrderAreValid = 1


    def getPreAndPostOrderAbove(self, nodeSpecifier):
        """Returns 2 lists of node numbers, preOrder and postOrder.

        This uses a stack, not recursion, so it should work for large
        trees without bumping into the recursion limit.  The 2 lists
        are relative to the node specified, and include the node
        specified.  PreOrder starts from theNode and goes to the tips;
        postOrder starts from the tips and goes to theNode."""

        gm = ['Tree.getPreAndPostOrderAbove()']
        theNode = self.node(nodeSpecifier)
        preOrder = [] # nodeNum's
        postOrder = []
        if not theNode.leftChild:
            preOrder.append(theNode.nodeNum)
            postOrder.append(theNode.nodeNum)
            return preOrder, postOrder
        stack = []    # nodes
        stack.append(theNode)
        preOrder.append(theNode.nodeNum)
        while len(stack):

            if stack[-1].leftChild:
                #print 'leftChild: %i' % stack[-1].leftChild.nodeNum
                theNodeNum = stack[-1].leftChild.nodeNum
                stack.append(stack[-1].leftChild)
                preOrder.append(theNodeNum)
            elif stack[-1].sibling:
                #print 'sibling: %i' % stack[-1].sibling.nodeNum
                theNodeNum = stack[-1].sibling.nodeNum
                theSib = stack[-1].sibling
                #print '                 postOrder appending u %i' % stack[-1].nodeNum
                postOrder.append(stack[-1].nodeNum)
                stack.pop()
                stack.append(theSib)
                preOrder.append(theNodeNum)
            else:
                #print '                 postOrder appending  v %i' % stack[-1].nodeNum
                postOrder.append(stack[-1].nodeNum)
                stack.pop()
                if len(stack) == 0:
                    break
                #print '                 postOrder appending  w %i' % stack[-1].nodeNum
                postOrder.append(stack[-1].nodeNum)
                theNode = stack.pop()
                while not theNode.sibling:
                    if len(stack) == 0:
                        break
                    #print '                 postOrder appending x %i' % stack[-1].nodeNum
                    postOrder.append(stack[-1].nodeNum)
                    theNode = stack.pop()
                if len(stack) == 0:
                    break
                if theNode.sibling:
                    stack.append(theNode.sibling)
                    preOrder.append(theNode.sibling.nodeNum)
                else:
                    gm.append('Problemo.')
                    gm.append('xxx got theNode %s' % theNode.nodeNum)
                    raise Glitch, gm

        return preOrder, postOrder

    def getPreAndPostOrderAboveRoot(self):
        """Sets self.preOrder and self.postOrder.

        This uses a stack, not recursion, so it should work for large
        trees without bumping into the recursion limit.  PreOrder
        starts from the root and goes to the tips; postOrder starts
        from the tips and goes to the root."""

        gm = ['Tree.getPreAndPostOrderAboveRoot()']
        theNode = self.root
        preOrdIndx = 0
        postOrdIndx = 0
        if type(self.preOrder) == types.NoneType or type(self.postOrder) == types.NoneType or len(self.preOrder) != len(self.nodes) or len(self.postOrder) != len(self.nodes):
            if var.usePfAndNumpy:
                self.preOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
                self.postOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
            else:
                self.preOrder = [var.NO_ORDER] * len(self.nodes)
                self.postOrder = [var.NO_ORDER] * len(self.nodes)

        #print "xx self.preOrder=%s" % self.preOrder
        if not theNode.leftChild:
            self.preOrder[preOrdIndx] = theNode.nodeNum
            preOrdIndx += 1
            self.postOrder[postOrdIndx] = theNode.nodeNum
            postOrdIndx += 1
        else:
            stack = []    # nodes
            stack.append(theNode)
            self.preOrder[preOrdIndx] = theNode.nodeNum
            preOrdIndx += 1
            while len(stack):

                if stack[-1].leftChild:
                    #print 'leftChild: %i  (%s)' % (stack[-1].leftChild.nodeNum, [n.nodeNum for n in stack])
                    theNodeNum = stack[-1].leftChild.nodeNum
                    stack.append(stack[-1].leftChild)
                    self.preOrder[preOrdIndx] = theNodeNum
                    preOrdIndx += 1
                elif stack[-1].sibling:
                    #print 'sibling: %i  (%s)' % (stack[-1].sibling.nodeNum, [n.nodeNum for n in stack])
                    theNodeNum = stack[-1].sibling.nodeNum
                    theSib = stack[-1].sibling
                    #print '                 postOrder appending u %i' % stack[-1].nodeNum
                    self.postOrder[postOrdIndx] = stack[-1].nodeNum
                    postOrdIndx += 1
                    stack.pop()
                    stack.append(theSib)
                    try:
                        self.preOrder[preOrdIndx] = theNodeNum
                    except IndexError:
                        gm.append("preOrdIndx=%s, theNodeNum=%i" % (preOrdIndx, theNodeNum))
                        gm.append("preOrder = %s" % self.preOrder)
                        raise Glitch, gm
                    preOrdIndx += 1
                else:
                    #print '                 postOrder appending  v %i' % stack[-1].nodeNum
                    self.postOrder[postOrdIndx] = stack[-1].nodeNum
                    postOrdIndx += 1
                    stack.pop()
                    if len(stack) == 0:
                        break
                    #print '                 postOrder appending  w %i' % stack[-1].nodeNum
                    self.postOrder[postOrdIndx] = stack[-1].nodeNum
                    postOrdIndx += 1
                    theNode = stack.pop()
                    while not theNode.sibling:
                        if len(stack) == 0:
                            break
                        #print '                 postOrder appending x %i' % stack[-1].nodeNum
                        self.postOrder[postOrdIndx] = stack[-1].nodeNum
                        postOrdIndx += 1
                        theNode = stack.pop()
                    if len(stack) == 0:
                        break
                    if theNode.sibling:
                        stack.append(theNode.sibling)
                        #print "self.preOrder = %s, preOrdIndx=%i" % (self.preOrder, preOrdIndx)
                        self.preOrder[preOrdIndx] = theNode.sibling.nodeNum
                        preOrdIndx += 1
                    else:
                        gm.append('Problemo.')
                        gm.append('xxx got theNode %s' % theNode.nodeNum)
                        raise Glitch, gm
        if 1:
            assert preOrdIndx == postOrdIndx
            #print "a preOrdIndx = %i, len(self.nodes) = %i" % (preOrdIndx, len(self.nodes))
            if preOrdIndx != len(self.nodes):
                pOI = preOrdIndx
                for i in range(preOrdIndx, len(self.nodes)):
                    self.preOrder[i] = var.NO_ORDER
                    self.postOrder[i] = var.NO_ORDER
                    preOrdIndx += 1
                    postOrdIndx += 1
            #print "b preOrdIndx = %i, len(self.nodes) = %i" % (preOrdIndx, len(self.nodes))
        assert preOrdIndx == len(self.nodes) and postOrdIndx == len(self.nodes)


    def iterPreOrder(self):
        """Node generator.  Assumes preAndPostOrderAreValid."""

        for i in self.preOrder:
            j = int(i)           # this speeds things up a lot!  Two-fold!
            if j != var.NO_ORDER:
                #yield self.nodes[int(i)]
                yield self.nodes[j]
    
    def iterPostOrder(self):
        """Node generator.  Assumes preAndPostOrderAreValid."""

        for i in self.postOrder:
            j = int(i)
            if j != var.NO_ORDER:
                yield self.nodes[j]
    
    def iterNodes(self):
        """Node generator, in preOrder.  Assumes preAndPostOrderAreValid."""

        for i in self.preOrder:
            j = int(i)
            if j != var.NO_ORDER:
                yield self.nodes[j]

    def iterNodesNoRoot(self):
        """Node generator, skipping the root. PreOrder."""

        for i in self.preOrder:
            j = int(i)
            if j != var.NO_ORDER:
                n = self.nodes[j]
                if n != self.root:
                    yield n
                    
    def iterLeavesNoRoot(self):
        """Leaf node generator, skipping the root.  PreOrder."""

        for i in self.preOrder:
            j = int(i)
            if j != var.NO_ORDER:
                n = self.nodes[j]
                if n != self.root and n.isLeaf:
                    yield n
                
    def iterLeavesPreOrder(self):

        for i in self.preOrder:
            j = int(i)
            if j != var.NO_ORDER:
                n = self.nodes[j]
                if n.isLeaf:         
                    yield n
                
    def iterLeavesPostOrder(self):

        for i in self.postOrder:
            j = int(i)
            if j != var.NO_ORDER:
                n = self.nodes[j]
                if n.isLeaf:
                    yield n
   

    def iterInternalsNoRootPreOrder(self):
        """Internal post order node generator, skipping the root.  Assumes preAndPostOrderAreValid."""

        for i in self.preOrder:
            j = int(i)
            if j != var.NO_ORDER:
                n = self.nodes[j]
                if n != self.root and not n.isLeaf:
                    yield n

    def iterInternalsNoRootPostOrder(self):
        """Internal post order node generator, skipping the root.  Assumes preAndPostOrderAreValid."""

        for i in self.postOrder:
            j = int(i)
            if j != var.NO_ORDER:
                n = self.nodes[j]
                if n != self.root and not n.isLeaf:
                    yield n
                
    def iterInternalsPostOrder(self):
        """Internal post order node generator.  Assumes preAndPostOrderAreValid."""

        for i in self.postOrder:
            j = int(i)
            if j != var.NO_ORDER:
                n = self.nodes[j]
                if not n.isLeaf:
                    yield n

    def iterInternalsNoRoot(self):
        """Internal node generator, skipping the root.  PreOrder"""

        for i in self.preOrder:
            j = int(i)
            if j != var.NO_ORDER:
                n = self.nodes[j]
                if n != self.root and not n.isLeaf:
                    yield n
            
    def iterInternals(self):
        """Internal node generator. PreOrder.  Including the root, if it is internal."""

        for i in self.preOrder:
            j = int(i)
            if j != var.NO_ORDER:
                n = self.nodes[j]
                if not n.isLeaf:
                    yield n

    ################################################
    #

    def stripBrLens(self):
        """Sets all node.br.len's to 0.1, the default in p4.

        Then, if you were to write it out in Nexus or Newick format,
        no branch lengths would be printed.
        """
        for n in self.iterNodesNoRoot():
            n.br.len = 0.1

    def getLen(self):
        """Return the sum of all br.len's."""
        if self.cTree:
            # About 0.01 msec for a tree with 1000 leaves
            return pf.p4_getTreeLen(self.cTree)
        else:
            # About 2 msec for a tree with 1000 leaves.
            total = 0.0
            for n in self.iterNodesNoRoot():
                total += n.br.len
            return total

    # def lenInternals(self):
    #     """Return the sum of all internal br.len's."""
    #     total = 0.0
    #     for n in self.iterInternalsNoRoot():
    #         total += n.br.len
    #     return total

    # def stemminess(self):
    #     """Ratio of internal branches to overall tree length.
        
    #     Also called 'treeness'.  Via Phillips and Penny, MPE 2003, but
    #     they cite Lanyon 1988."""

    #     total = self.len()
    #     internals = self.lenInternals()
    #     return internals/total
    
    def _makeRCSplitKeys(self, splitList=None):
        """Make long integer-valued split keys.
        This is dependent on that the leaf bitkeys are already set, ie this method should only 
        be used by the reduced consensus supertree class
        """
        if not self.preAndPostOrderAreValid:
            self.setPreAndPostOrder()

        #self.draw()

#        allOnes = 2L**(self.nTax) - 1
        #print 'nTax = %i, allOnes = %i' % (self.nTax, allOnes)

#        for n in self.iterInternalsPostOrder():
        for n in self.iterInternalsNoRootPostOrder():
#            if n != self.root:
                #print 'doing node %s' % n.nodeNum
#            if not n.isLeaf():
            if n.leftChild:

                childrenNums = self.getChildrenNums(n)
#                    x = childrenNums[0]
                x = self.nodes[childrenNums[0]].br.rc
                for i in childrenNums[1:]:
                    x = x | self.nodes[i].br.rc
#                        x = x | i
                n.br.rc = x
                if splitList != None:
                    splitList.append([x,0])
                n.br.rcList = [n.br.rc]

    def makeSplitKeys(self, makeNodeForSplitKeyDict=False):
        """Make long integer-valued split keys.

        This needs to have self.taxNames set.

        We make 2 kinds of split keys-- rawSplitKeys and splitKeys.
        Both are attributes of node.br, so we have eg node.br.splitKey.

        Raw split keys for terminal nodes are 2**n, where n is the index
        of the taxon name.  Eg for the first taxon, the rawSplitKey will
        be 1, for the 3rd taxon the rawSplitKey will be 4.

        RawSplitKeys for internal nodes are the rawSplitKey's for the
        children, bitwise OR'ed together.

        SplitKeys, cf rawSplitKeys, are in 'standard form', where the
        numbers are even, ie do not contain the 1-bit.  Having it in
        standard form means that you can compare splits among trees.
        If the rawSplitKey is even, then the splitKey is simply that,
        unchanged.  If, however, the rawSplitKey is odd, then the
        splitKey is the rawSplitKey bit-flipped.  For example, if
        there are 5 taxa, and one of the rawSplitKeys is 9 (odd), we
        can calculate the splitKey by bit-flipping, as::

            01001 =  9   rawSplitKey
            10110 = 22   splitKey

        (Bit-flipping is done by exclusive-or'ing (xor) with 11111.)

        The splitKey is readily converted to a splitString for
        display, as 22 becomes '.**.*' (note the '1' bit is now on the
        left).  It is conventional that the first taxon, on the left,
        is always a dot.  (I don't know where the convention comes
        from.)

        The root has no rawSplitKey or splitKey.

        For example, the tree::

                +-------2:B (rawSplitKey = 2)
            +---1
            |   +---------3:C  (rawSplitKey = 4)
            |
            0-------------4:E  (rawSplitKey = 16)
            |
            |    +-----6:A  (rawSplitKey = 1)
            +----5
                 +-----------7:D  (rawSplitKey = 8)

        has 2 internal splits, on nodes 1 and 5.

        ::
        
            Node      n.br.rawSplitKey     n.br.splitKey
            1             6                    6
            5             9                   22

        There should be no duplicated rawSplitKeys, but if the tree
        has a bifurcating root then there will be a duped splitKey.

        This method will fail for trees with internal nodes that have
        only one child, because that will make duplicated splits.

        If arg *makeNodeForSplitKeyDict* is set, then it will make a
        dictionary ``nodeForSplitKeyDict`` where the keys are the
        splitKeys and the values are the corresponding nodes.
        
        """

        gm = ['Tree.makeSplitKeys()']
        #raise Glitch, gm
        if not self.taxNames:
            gm.append('No taxNames.')
            raise Glitch, gm

        if not self.preAndPostOrderAreValid:
            self.setPreAndPostOrder()

        #self.draw()
        if makeNodeForSplitKeyDict:
            self.nodeForSplitKeyDict = {}

        allOnes = 2L**(self.nTax) - 1
        #print 'nTax = %i, allOnes = %i' % (self.nTax, allOnes)

        for n in self.iterPostOrder():
            if n != self.root:
                #print 'doing node %s' % n.nodeNum

                if not n.leftChild:
                    # A long int, eg 1L, has no upper limit on its value
                    try:
                        n.br.rawSplitKey = 1L << self.taxNames.index(n.name)  # "<<" is left-shift
                    except ValueError:
                        gm.append('node.name %s' % n.name)
                        gm.append('is not in taxNames %s' % self.taxNames)
                        raise Glitch, gm
                    #n.br.rawSplitKey = 1L << self.taxNames.index(n.name)  # "<<" is left-shift
                    
                    if n.br.rawSplitKey == 1:
                        n.br.splitKey = allOnes - 1
                    else:
                        n.br.splitKey = n.br.rawSplitKey

                    #print 'upward leaf   node %s, rawSplitKey %s, splitKey %s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
                else:
                    childrenNums = self.getChildrenNums(n)
                    if len(childrenNums) == 1:
                        gm.append('Internal node has only one child.  That will make a duped split.')
                        raise Glitch, gm
                    x = self.nodes[childrenNums[0]].br.rawSplitKey
                    for i in childrenNums[1:]:
                        y = self.nodes[i].br.rawSplitKey
                        x = x | y  # '|' is bitwise "OR".
                    n.br.rawSplitKey = x

                    # Make node.br.splitKey's in "standard form", ie
                    # without the first taxon, ie without a 1.  To do that
                    # we use the '&' operator, to bitwise "and" with 1.
                    if 1 & n.br.rawSplitKey: # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
                        n.br.splitKey = allOnes ^ n.br.rawSplitKey # "^" is xor, a bit-flipper.
                    else:
                        n.br.splitKey = n.br.rawSplitKey
                    #print 'intern node %s, rawSplitKey %s, splitKey %s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
                if makeNodeForSplitKeyDict:
                    self.nodeForSplitKeyDict[n.br.splitKey] = n


        if 0:
            # There should be no duped rawSplitKeys
            theKeys = []
            for n in self.iterNodesNoRoot():
                theKeys.append(n.br.rawSplitKey)
            for k in theKeys:
                if theKeys.count(k) > 1:
                    gm.append('Duped rawSplitKey %i.' % k)
                    for n in self.nodes:
                        if n != self.root:
                            print '%7s     %4s       %4s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
                    raise Glitch, gm

            # Any duped splitKeys?  There will be if the tree is bi-Rooted.
            if 0:
                theKeys = []
                for n in self.iterNodesNoRoot():
                    theKeys.append(n.br.splitKey)
                for k in theKeys:
                    if theKeys.count(k) > 1:
                        gm.append('Duped splitKey %i.' % k)
                        for n in self.iterNodesNoRoot():
                            print '%7s     %4s       %4s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
                        raise Glitch, gm

        if 0:
            print gm[0]
            print self.taxNames
            print 'nodeNum  rawSplitKey  splitKey'
            for n in self.iterNodesNoRoot():
                print '%7s     %4s       %4s        %s' % (
                    n.nodeNum, n.br.rawSplitKey, n.br.splitKey, func.getSplitStringFromKey(n.br.splitKey, self.nTax))


    def recalculateSplitKeysOfNodeFromChildren(self, aNode, allOnes):
        children = [n for n in aNode.iterChildren()]
        x = children[0].br.rawSplitKey
        for n in children[1:]:
            x = x | n.br.rawSplitKey   # '|' is bitwise "OR".
        aNode.br.rawSplitKey = x
        if 1 & aNode.br.rawSplitKey: # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
            aNode.br.splitKey = allOnes ^ aNode.br.rawSplitKey # "^" is xor, a bit-flipper.
        else:
            aNode.br.splitKey = aNode.br.rawSplitKey
        

    def checkSplitKeys(self, useOldName=False, glitch=True, verbose=True):
        gm = ['Tree.checkSplitKeys()']
        
        allOnes = 2L**(self.nTax) - 1
        #print 'nTax = %i, allOnes = %i' % (self.nTax, allOnes)

        isBad = False
        for n in self.iterPostOrder():
            if n != self.root:
                #print 'doing node %s' % n.nodeNum

                if not n.leftChild:
                    # A long int, eg 1L, has no upper limit on its value
                    try:
                        if useOldName:
                            rawSplitKey = 1L << self.taxNames.index(n.oldName)  # "<<" is left-shift
                        else:
                            rawSplitKey = 1L << self.taxNames.index(n.name)  # "<<" is left-shift
                    except ValueError:
                        if useOldName:
                            gm.append('node.name %s' % n.oldName)
                        else:
                            gm.append('node.name %s' % n.name)
                        gm.append('is not in taxNames %s' % self.taxNames)
                        raise Glitch, gm
                    #n.br.rawSplitKey = 1L << self.taxNames.index(n.name)  # "<<" is left-shift

                    if rawSplitKey == 1:
                        splitKey = allOnes - 1
                    else:
                        splitKey = rawSplitKey

                    #print 'upward leaf   node %s, rawSplitKey %s, splitKey %s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
                else:
                    childrenNums = self.getChildrenNums(n)
                    if len(childrenNums) == 1:
                        gm.append('Internal node has only one child.  That will make a duped split.')
                        raise Glitch, gm
                    x = self.nodes[childrenNums[0]].br.rawSplitKey
                    for i in childrenNums[1:]:
                        y = self.nodes[i].br.rawSplitKey
                        x = x | y  # '|' is bitwise "OR".
                    rawSplitKey = x

                    # Make node.br.splitKey's in "standard form", ie
                    # without the first taxon, ie without a 1.  To do that
                    # we use the '&' operator, to bitwise "and" with 1.
                    if 1 & rawSplitKey: # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
                        splitKey = allOnes ^ rawSplitKey # "^" is xor, a bit-flipper.
                    else:
                        splitKey = rawSplitKey
                    #print 'intern node %s, rawSplitKey %s, splitKey %s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
                if n.br.rawSplitKey != rawSplitKey:
                    print "checkSplitKeys node %2i rawSplitKey: existing %s, calculated %s" % (n.nodeNum, n.br.rawSplitKey, rawSplitKey)
                    isBad = True
                if n.br.splitKey != splitKey:
                    print "checkSplitKeys node %2i splitKey: existing %s, calculated %s" % (n.nodeNum, n.br.splitKey, splitKey)
                    isBad = True
        if glitch and isBad:
            raise Glitch, gm
        if verbose and not isBad:
            print "checkSplitKeys().  ok"


    def taxSetIsASplit(self, taxSetName):

        assert self.nexusSets
        assert self.taxNames
        assert self.nexusSets.taxSets
        lowArgTaxSetName = string.lower(taxSetName)
        theTS = None
        for ts in self.nexusSets.taxSets:
            if ts.lowName == lowArgTaxSetName:
                theTS = ts
                break
        assert theTS, "Could not find the taxSet named %s" % taxSetName
        #theTS.dump()
        assert theTS.mask
        rawSplitKey = 0L
        for i in range(len(theTS.mask)):
            #print i, theTS.mask[i]
            if theTS.mask[i] == '1':
                rawSplitKey += (1L << i)
        if 1 & rawSplitKey: # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
            allOnes = 2L**(self.nTax) - 1
            splitKey = allOnes ^ rawSplitKey # "^" is xor, a bit-flipper.
        else:
            splitKey = rawSplitKey
        #print "got splitKey %s" % splitKey

        for n in self.nodes:
            if n.br and not n.isLeaf:
                #print "  %2i  %s  %s" % (n.nodeNum, n.br.splitKey, splitKey)
                if n.br.splitKey == splitKey:
                    #self.draw()
                    #return n.nodeNum
                    return n
        return None    # Was -1 before, when n.nodeNum was returned for hits.  Now a node is returned.
        
        
                


    def checkTaxNames(self):
        """Check that all taxNames are in the tree, and vice versa."""

        #If var.allowTreesWithDifferingTaxonSets is set to True we will not check
        #the taxnames. This is primarily used to read trees for supertree and 
        #leafstability calcuations.

        # Peter comments that Tobias added this, but it is not needed,
        # and messes other things up -- so comment out until it is
        # sorted.
        #if var.allowTreesWithDifferingTaxonSets:
        #     return
         

        # self.taxNames is a property, so when it is set, it calls this method
        
        if self.name:
            gm = ['Tree.checkTaxNames()   tree %s' % self.name]
        else:
            gm = ['Tree.checkTaxNames()']

        if not self.taxNames:
            gm.append('No taxNames.')
            raise Glitch, gm

        tax = []
        for n in self.iterNodes():
            if n.isLeaf and n.name:
                tax.append(n.name)
            # This next line should not be needed, as root leaves should be leaves.
            elif n == self.root and n.name and n.getNChildren() < 2:  # terminal root that has a taxName
                tax.append(n.name)

        #print 'tax from tree = %s' % tax
        #print 'self.taxNames = %s' % self.taxNames

        # Check for same number of taxa
        if len(tax) != len(self.taxNames):
            #self.draw()
            #self.dump(node=1)
            gm.append('Mismatch.  Length of self.taxNames is wrong.')
            gm.append('The tree has %i leaves with names, but len(self.taxNames) = %i' % (
                len(tax), len(self.taxNames)))
            gm.append('leaves on the tree = %s' % tax)
            gm.append('self.taxNames = %s' % self.taxNames)
            gm.append('symmetric_difference is %s' % set(tax).symmetric_difference(set(self.taxNames)))
            #print '    Setting invalid taxNames to an empty list.'
            #self.taxNames = []
            raise Glitch(gm, 'tree_badTaxNamesLength')

        # Check for mis-matched taxNames
        isBad = 0
        taxSet = set(tax)
        selfTaxNamesSet = set(self.taxNames)

        s = taxSet.difference(selfTaxNamesSet)
        if len(s):
            isBad = 1
            print gm[0]
            print 'TaxName mismatch between the tree and self.taxNames.'
            print 'These taxa are found in the tree but not in self.taxNames:'
            print s
        s = selfTaxNamesSet.difference(taxSet)
        if len(s):
            isBad = 1
            print gm[0]
            print 'TaxName mismatch between the tree and self.taxNames.'
            print 'These taxa are found in the self.taxNames but not in the tree:'
            print s

        if isBad:
            raise Glitch(gm, 'tree_taxNamesMisMatch')


    ################################################
    # Copy and Verify

    def copyToTree(self, otherTree):

        gm = ['Tree.copyToTree()']

        if len(self.nodes) != len(otherTree.nodes):
            gm.append('Different number of nodes.')
            raise Glitch, gm

        # node relations (parent, child, sib)
        for nNum in range(len(self.nodes)):
            selfNode = self.nodes[nNum]
            otherNode = otherTree.nodes[nNum]

            # parent
            if selfNode.parent:
                otherNode.parent = otherTree.nodes[selfNode.parent.nodeNum]
            else:
                #print "otherNode.parent", otherNode.parent
                otherNode.parent = None

            # leftChild
            if selfNode.leftChild:
                otherNode.leftChild = otherTree.nodes[selfNode.leftChild.nodeNum]
            else:
                #print "otherNode.leftChild", otherNode.leftChild
                otherNode.leftChild = None

            # sibling
            if selfNode.sibling:
                otherNode.sibling = otherTree.nodes[selfNode.sibling.nodeNum]
            else:
                #print "otherNode.sibling", otherNode.sibling
                otherNode.sibling = None

        # root
        otherTree.root.br = otherTree.nodes[self.root.nodeNum].br
        otherTree.nodes[self.root.nodeNum].br = None
        otherTree.root = otherTree.nodes[self.root.nodeNum]

        # brLens and splitKeys
        for nNum in range(len(self.nodes)):
            if self.nodes[nNum] != self.root:
                otherTree.nodes[nNum].br.len = self.nodes[nNum].br.len
                otherTree.nodes[nNum].br.lenChanged = self.nodes[nNum].br.lenChanged
                #otherTree.nodes[nNum].br.flag = self.nodes[nNum].flag
                otherTree.nodes[nNum].br.splitKey = self.nodes[nNum].br.splitKey
                otherTree.nodes[nNum].br.rawSplitKey = self.nodes[nNum].br.rawSplitKey

        # model usage numbers
        if self.model:
            for nNum in range(len(self.nodes)):
                selfNode = self.nodes[nNum]
                otherNode = otherTree.nodes[nNum]
                for pNum in range(self.model.nParts):
                    otherNode.parts[pNum].compNum = selfNode.parts[pNum].compNum
                    if selfNode != self.root:
                        otherNode.br.parts[pNum].rMatrixNum = selfNode.br.parts[pNum].rMatrixNum
                        otherNode.br.parts[pNum].gdasrvNum = selfNode.br.parts[pNum].gdasrvNum

        # pre- and postOrder
        for i in range(len(self.preOrder)):
            otherTree.preOrder[i] = self.preOrder[i]
            otherTree.postOrder[i] = self.postOrder[i]

        # partLikes
        if self.model:
            for pNum in range(self.model.nParts):
                otherTree.partLikes[pNum] = self.partLikes[pNum]

        otherTree._nInternalNodes = self._nInternalNodes

    def verifyIdentityWith(self, otherTree, doSplitKeys):
        """For MCMC debugging.  Verifies that two trees are identical."""

        complaintHead = '\nTree.verifyIdentityWith()' # keep

        if len(self.nodes) != len(otherTree.nodes):
            print complaintHead
            print '    Different number of nodes.'
            return var.DIFFERENT

        # check node relations (parent, child, sib)
        isBad = 0
        for nNum in range(len(self.nodes)):
            selfNode = self.nodes[nNum]
            otherNode = otherTree.nodes[nNum]

            # parent
            if selfNode.parent:
                if otherNode.parent:
                    if otherNode.parent.nodeNum != selfNode.parent.nodeNum:
                        isBad = 1
                else:
                    isBad = 1
            else:
                if otherNode.parent:
                    isBad = 1

            # leftChild
            if selfNode.leftChild:
                if otherNode.leftChild:
                    if otherNode.leftChild.nodeNum != selfNode.leftChild.nodeNum:
                        isBad = 1
                else:
                    isBad = 1
            else:
                if otherNode.leftChild:
                    isBad = 1

            # sibling
            if selfNode.sibling:
                if otherNode.sibling:
                    if otherNode.sibling.nodeNum != selfNode.sibling.nodeNum:
                        isBad = 1
                else:
                    isBad = 1
            else:
                if otherNode.sibling:
                    isBad = 1

            if isBad:
                print complaintHead
                print '    Node %i, relations differ.' % nNum
                self.write()
                otherTree.write()
                return var.DIFFERENT

        if self.root.nodeNum != otherTree.root.nodeNum:
            print complaintHead
            print '    Roots differ.'
            return var.DIFFERENT

        # brLens, lenChanged, and node.flag. and splitKeys
        for nNum in range(len(self.nodes)):
            if self.nodes[nNum] != self.root:
                #if self.nodes[nNum].br.len != otherTree.nodes[nNum].br.len:
                if math.fabs(self.nodes[nNum].br.len - otherTree.nodes[nNum].br.len) > 1.e-8:
                    print complaintHead
                    print '    BrLens differ.'
                    return var.DIFFERENT
                if self.nodes[nNum].br.lenChanged != otherTree.nodes[nNum].br.lenChanged:
                    print complaintHead
                    print '    br.lenChanged differs.'
                    return var.DIFFERENT
                if self.nodes[nNum].flag != otherTree.nodes[nNum].flag:
                    print complaintHead
                    print '    flag differs, nodeNum %i.  %s vs %s' % (nNum, self.nodes[nNum].flag, otherTree.nodes[nNum].flag)
                    return var.DIFFERENT
                if doSplitKeys:
                    if self.nodes[nNum].br.splitKey != otherTree.nodes[nNum].br.splitKey:
                        print complaintHead
                        print '    SplitKeys differ.'
                        return var.DIFFERENT
                    if self.nodes[nNum].br.rawSplitKey != otherTree.nodes[nNum].br.rawSplitKey:
                        print complaintHead
                        print '    rawSplitKeys differ.'
                        return var.DIFFERENT

        # model usage numbers
        isBad = 0
        for pNum in range(self.model.nParts):
            for nNum in range(len(self.nodes)):
                selfNode = self.nodes[nNum]
                otherNode = otherTree.nodes[nNum]
                if selfNode.parts[pNum].compNum != otherNode.parts[pNum].compNum:
                    isBad = 1
                if self.nodes[nNum] != self.root:
                    if selfNode.br.parts[pNum].rMatrixNum != otherNode.br.parts[pNum].rMatrixNum:
                        isBad = 1
                    elif selfNode.br.parts[pNum].gdasrvNum != otherNode.br.parts[pNum].gdasrvNum:
                        isBad = 1

                if isBad:
                    print complaintHead
                    print '    Node %i, model usage info does not match.' % nNum
                    return var.DIFFERENT

        # pre- and postOrder
        isBad = 0
        for i in range(len(self.preOrder)):
            if self.preOrder[i] != otherTree.preOrder[i]:
                isBad = 1
                break
            elif self.postOrder[i] != otherTree.postOrder[i]:
                isBad = 1
                break
        if isBad:
            print complaintHead
            print '    Pre- or postOrder do not match.'
            return var.DIFFERENT

        if self._nInternalNodes != otherTree._nInternalNodes:
            print complaintHead
            print '    _nInternalNodes differ.'
            return var.DIFFERENT

        # partLikes
        for pNum in range(self.model.nParts):
            #if otherTree.partLikes[pNum] != self.partLikes[pNum]:
            if math.fabs(otherTree.partLikes[pNum] - self.partLikes[pNum]) > 1.e-8:
                print complaintHead
                print "    partLikes differ.  (%.5f, (%g) %.5f (%g)" % (
                    otherTree.partLikes[pNum], otherTree.partLikes[pNum], self.partLikes[pNum], self.partLikes[pNum])
                return var.DIFFERENT

        if 0: # some more
            for nNum in range(len(self.nodes)):
                selfNode = self.nodes[nNum]
                otherNode = otherTree.nodes[nNum]
                if selfNode.nodeNum != otherNode.nodeNum:
                    print complaintHead
                    print '    nodeNum differs'
                    return var.DIFFERENT
                if selfNode.seqNum != otherNode.seqNum:
                    print complaintHead
                    print '    seqNum differs'
                    return var.DIFFERENT
                if selfNode.name != otherNode.name:
                    print complaintHead
                    print '    name differs'
                    return var.DIFFERENT
                if selfNode.isLeaf != otherNode.isLeaf:
                    print complaintHead
                    print '    isLeaf differs'
                    return var.DIFFERENT

        
        return var.SAME


    ############################################


    def isFullyBifurcating(self, verbose=False):
        """Returns True if the tree is fully bifurcating.  Else False. """
        
        if self.root and self.root.leftChild and self.root.leftChild.sibling and self.root.leftChild.sibling.sibling:
            if self.root.leftChild.sibling.sibling.sibling:
                if verbose:
                    print "isFullyBifurcating() returning False, due to root with 4 or more children." 
                return False
        elif self.root and self.root.isLeaf:
            pass
        else:
            if verbose:
                print "isFullyBifurcating() returning False, due to (non-leaf) root not having 3 children."
            return False
        for n in self.iterInternalsNoRoot():
            if n.leftChild and n.leftChild.sibling:
                if n.leftChild.sibling.sibling:
                    if verbose:
                        print "isFullyBifurcating() returning False, due to node %i having 3 or more children." % n.nodeNum
                    return False
            else:
                if verbose:
                    print "isFullyBifurcating() returning False, due to non-leaf node %i having 1 or fewer children." % n.nodeNum
                return False
        return True

    # These next two are for the eTBR implementation that I got from Jason Evans' Crux.  Thanks Jason!
    def getDegree(self, nodeSpecifier):
        n = self.node(nodeSpecifier)
        if n.isLeaf:
            if n.parent:
                return 1
            else:
                return 0
        else:
            #assert n.leftChild
            deg = 1 # the leftChild
            if n.parent:
                deg += 1
            s = n.leftChild.sibling
            while s:
                deg += 1
                s = s.sibling
            return deg

    def nextNode(self, spokeSpecifier, hubSpecifier):
        """Get next node cycling around a hub node.

        A bit of a hack to make a p4 Node behave sorta like a
        Felsenstein node.  Imagine cycling around the branches
        emanating from a node like spokes on a hub, starting from
        anywhere, with no end.

        The hub node would usually be the parent of the spoke, or the
        spoke would be the hub itself.  Usually
        self.nextNode(spoke,hub) delivers spoke.sibling.  What happens
        when the siblings run out is that
        self.nextNode(rightmostSibling, hub) delivers hub itself, and
        of course its branch (spoke) points toward the hub.parent.
        (Unless hub is the root, of course, in which case
        self.nextNode(rightmostSibling, hub) delivers hub.leftChild.)
        In the usual case of the hub not being the root, the next node
        to be delivered by nextNode(spokeIsHub, hub) is usually the
        leftChild of the hub.  Round and round, clockwise.

        """

        spoke = self.node(spokeSpecifier)
        hub = self.node(hubSpecifier)
        
        if spoke.parent == hub or hub == spoke:
            if spoke == hub:
                if spoke.leftChild:
                    return spoke.leftChild
                else:
                    return hub
            else:
                if spoke.sibling:
                    return spoke.sibling
                else:
                    if hub.parent:
                        return hub
                    else:
                        return hub.leftChild
        else:
            print "*=" * 25
            self.draw()
            gm = ["Tree.nextNode() spoke=%i, hub=%i" % (spoke.nodeNum, hub.nodeNum)]
            gm.append("Need to have either spoke.parent == hub or hub == spoke.")
            raise Glitch, gm
            
        

    def topologyDistance(self, tree2, metric='sd', resetSplitKeySet=False):
        """Compares the topology of self with tree2.

        The complete list of metrics is given in var.topologyDistanceMetrics

        For most metrics using this method, taxNames needs to be set,
        to the same in the two trees.  If the taxa differ, this method
        simply returns -1

        The 'metric' can be one of 'sd' (symmetric difference), 'wrf'
        (weighted Robinson-Foulds), 'bld' (Felsenstein's branch-
        length distance), or 'diffs'.  The unweighted Robinson-Foulds
        metric would be the same as the symmetric difference.

        There is also an experimental scqdist, but that needs the
        scqdist.so module, in the QDist directory.

        See Felsenstein 2004 Inferring Phylogenies, Pg 529.

        The default metric is the very simple 'sd', symmetric
        difference.  Using this metric, if the 2 trees share the same
        set of splits, they are deemed to be the same topology; branch
        lengths are not compared.  This method returns the number of
        splits that are in self that are not in tree2 plus the number
        of splits that are in tree2 that are not in self.  So it would
        return 0 for trees that are the same.

        The 'wrf' and 'bld' metrics take branch lengths into account.
        Bifurcating roots complicate things, so they are not allowed
        for weighted distance calculations.

        In the unweighted case (ie metric='sd'), whether the trees
        compared have bifurcating roots or not is ignored.  So the
        trees (A,B,(C,D)) and ((A,B),(C,D)) will be deemed to have the
        same topology, since they have the same splits.

        The measurement 'diffs', which returns a tuple of 2 numbers --
        both are set differences.  The first is the number of splits
        in self that are not in tree2, and the second is the number of
        splits in tree2 that are not in self.  (Consider it as the the
        symmetric difference split into its 2 parts.)

        If you calculate a distance and then make a topology change, a
        subsequent sd topologyDistance calculation will be wrong, as it
        uses previous splits.  So then you need to 'resetSplitKeySet'.
 
        The 'scqdist' metric also gives quartet distances.  It was
        written by Anders Kabell Kristensen for his Masters degree at
        Aarhus University, 2010.  http://www.cs.au.dk/~dalko/thesis/
        It has two versions -- a pure Python version (that needs
        scipy) that I do not include here, and a fast C++ version,
        that I wrapped in python.  Its speedy -- the 'sc' in 'scqdist'
        is for 'sub-cubic', ie better than O(n^3).

        """

        gm = ['Tree.topologyDistance()']
        
        if metric not in var.topologyDistanceMetrics:
            gm.append("Got a request for unknown metric '%s'" % metric)
            gm.append("The 'metric' arg should be one of %s" % var.topologyDistanceMetrics)
            raise Glitch, gm
        if metric == 'scqdist': # no need for taxNames
            try:
                import scqdist
            except ImportError:
                gm.append("Could not find the 'scqdist' module needed for this metric.")
                gm.append("See the instructions for making it in the p4 source, in the Qdist directory.")
                raise Glitch, gm
            tsSelf = self.writeNewick(toString=True)
            tsTree2 = tree2.writeNewick(toString=True)
            return scqdist.qdist(tsSelf, tsTree2)
        if not self.taxNames or not tree2.taxNames:
            gm.append("This method requires taxNames to be set.")
            raise Glitch, gm
        if self.taxNames != tree2.taxNames:
            gm.append("The taxNames are different for the two trees.")
            gm.append("Self:  %s" % self.taxNames)
            gm.append("tree2: %s" % tree2.taxNames)
            raise Glitch, gm
        if (self.root.getNChildren() == 2 or tree2.root.getNChildren() == 2) and ( metric in ['wrf', 'bld']):
            gm.append('One of the input trees has a bifurcating root.')
            gm.append('Weighted tree distance calculations do not work on bi-rooted trees.')
            raise Glitch, gm


        # I might be doing a lot of these, and I don't want to make
        # splitKeys and make the splitKeys set over and over again.
        # So store it as a Tree.attribute.

        if resetSplitKeySet:
            if hasattr(self, 'splitKeySet'):
                del(self.splitKeySet)
            if hasattr(tree2, 'splitKeySet'):
                del(tree2.splitKeySet)

        if not hasattr(self, 'splitKeySet'):
            self.makeSplitKeys()
            self.splitKeySet = set([n.br.splitKey for n in self.iterNodesNoRoot()])

        if not hasattr(tree2, 'splitKeySet'):
            tree2.makeSplitKeys()
            tree2.splitKeySet = set([n.br.splitKey for n in tree2.iterNodesNoRoot()])

        if metric == 'sd':
            # Symmetric difference.  The symmetric_difference method
            # returns all elements that are in exactly one of the sets.
            theSD = len(self.splitKeySet.symmetric_difference(tree2.splitKeySet))
            return theSD

        # The difference method returns the difference of two sets as
        # a new Set.  I.e. all elements that are in self and not in
        # the other.
        selfHasButTree2DoesNot = self.splitKeySet.difference(tree2.splitKeySet)
        tree2HasButSelfDoesNot = tree2.splitKeySet.difference(self.splitKeySet)

        if metric == 'diffs':
            return len(selfHasButTree2DoesNot),len(tree2HasButSelfDoesNot)

        if metric in ['wrf', 'bld']:
            self.splitKeyHash = {}
            for n in self.iterNodesNoRoot():
                self.splitKeyHash[n.br.splitKey] = n
            tree2.splitKeyHash = {}
            for n in tree2.iterNodesNoRoot():
                tree2.splitKeyHash[n.br.splitKey] = n
            
        if metric == 'wrf':
            theSum = 0.0
            for k in self.splitKeySet.intersection(tree2.splitKeySet):
                #print '%s - %s' % (self.splitKeyHash[k].br.len, tree2.splitKeyHash[k].br.len)
                theSum += abs(self.splitKeyHash[k].br.len - tree2.splitKeyHash[k].br.len)
            for k in selfHasButTree2DoesNot:
                #print 'x %s' % self.splitKeyHash[k].br.len
                theSum += self.splitKeyHash[k].br.len
            for k in tree2HasButSelfDoesNot:
                #print 'y %s' % tree2.splitKeyHash[k].br.len
                theSum += tree2.splitKeyHash[k].br.len
            return theSum
        elif metric == 'bld':
            theSum = 0.0
            for k in self.splitKeySet.intersection(tree2.splitKeySet):
                theDiff = self.splitKeyHash[k].br.len - tree2.splitKeyHash[k].br.len
                theSum += theDiff * theDiff
            for k in selfHasButTree2DoesNot:
                theSum += self.splitKeyHash[k].br.len * self.splitKeyHash[k].br.len
            for k in tree2HasButSelfDoesNot:
                theSum += tree2.splitKeyHash[k].br.len * tree2.splitKeyHash[k].br.len
            #print 'branch score =', theSum
            return math.sqrt(theSum)


    ##################################################
    

    def tv(self):
        """Tree Viewer. Show the tree in a gui window.

        Needs Tkinter.

        If you have nexus taxsets defined, you can show them.
        """
        from BTV import TV
        #import os
        #os.environ['PYTHONINSPECT'] = '1'
        TV(self)
    
    def btv(self):
        """Big Tree Viewer. Show the tree in a gui window.

        This is for looking at big trees.  The viewer has 2 panels --
        one for an overall view of the whole tree, and one for a
        zoomed view, controlled by a selection rectangle on the whole
        tree view.

        Needs Tkinter.

        If you have nexus taxsets defined, you can show them.
        """
        from BTV import BTV
        #import os
        #os.environ['PYTHONINSPECT'] = '1'
        BTV(self)
    

    def tvTopologyCompare(self, treeB):
        """Graphically show topology differences.

        The taxNames need to be set, and need to be the same for both
        trees.

        (If the red lines don't show up right away, try adjusting the
        size of the windows slightly.)
        
        """
        
        sd = self.topologyDistance(treeB)
        if sd == 0:
            print "The trees are the same. No tv."
            return
        #for sk in self.splitKeyHash.iterkeys():
        #    if not treeB.splitKeyHash.has_key(sk):
        #        print "self has sk

        self.splitKeyHash = {}
        for n in self.iterInternalsNoRoot():
            self.splitKeyHash[n.br.splitKey] = n
        treeB.splitKeyHash = {}
        for n in treeB.iterNodesNoRoot():
            treeB.splitKeyHash[n.br.splitKey] = n
        assert self.splitKeySet
        assert treeB.splitKeySet
        selfHasButTreeBDoesnt = self.splitKeySet.difference(treeB.splitKeySet)
        treeBHasButSelfDoesnt = treeB.splitKeySet.difference(self.splitKeySet)

        from BTV import TV
        #import os
        #os.environ['PYTHONINSPECT'] = '1'
        TV(self, title='TV self')
        TV(treeB, title='TV treeB')
        for spl in selfHasButTreeBDoesnt:
            n = self.splitKeyHash[spl]
            n.br.color = 'red'
        for spl in treeBHasButSelfDoesnt:
            n = treeB.splitKeyHash[spl]
            n.br.color = 'orange'




