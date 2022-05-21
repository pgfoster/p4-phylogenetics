import sys
import string
import io
import math
import copy
import os
import p4.func
import time
import glob
from p4.var import var
from p4.p4exceptions import P4Error
from p4.node import Node, NodePart, NodeBranch, NodeBranchPart
from p4.nexustoken import nextTok, safeNextTok
from p4.distancematrix import DistanceMatrix

import numpy
import p4.pf as pf
from p4.model import Model
from p4.data import Data
from p4.alignment import Part
import random

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
    :meth:`~p4.tree.Tree.node`.  For example, from a Tree ``t`` you can
    get node number 3 via::

        n = t.node(3)

    and you can get the node that is the parent of the Mastodon via::

        n = t.node('Mastodon').parent

    For many methods that require specifying a node, the method argument is *nodeSpecifier*, eg::

        t.reRoot(23)

    ``reRoots``'s the tree to node number 23.

**Describe, draw, and get information about the tree**

    .. autosummary::
       :nosignatures:

       ~Tree.dump
       ~Tree.draw
       ~Tree.textDrawList
       ~Tree.tv
       ~Tree.btv
       ~Tree.isFullyBifurcating
       ~Tree.taxSetIsASplit
       ~Tree.getAllLeafNames
       ~Tree.getChildrenNums
       ~Tree.getDegree
       ~Tree.getLen
       ~Tree.getNodeNumsAbove
       ~Tree.getPreAndPostOrderAbove
       ~Tree.getPreAndPostOrderAboveRoot
       ~Tree.getSeqNumsAbove
       ~Tree.subTreeIsFullyBifurcating
       ~Tree.summarizeModelComponentsNNodes
       ~Tree.verifyIdentityWith

**Write**

    .. autosummary::
       :nosignatures:

       ~Tree.write
       ~Tree.writeNewick
       ~Tree.writeNexus
       ~Tree.writePhylip
       ~Tree.tPickle

    See also Trees methods :meth:`p4.trees.Trees.writeNexus` and
    :meth:`p4.trees.Trees.writeNewick` for doing trees by the bunch.

**Iteration over the nodes**

    Sometimes you don't want to just iterate over the self.nodes list,
    because after some manipulations a node might be in self.nodes but
    not actually in the tree; using these 'iter' methods takes care of
    that, skipping such nodes.

    .. autosummary::
       :nosignatures:

       ~Tree.iterInternals
       ~Tree.iterInternalsNoRoot
       ~Tree.iterInternalsNoRootPostOrder
       ~Tree.iterInternalsNoRootPreOrder
       ~Tree.iterInternalsPostOrder
       ~Tree.iterLeavesNoRoot
       ~Tree.iterLeavesPostOrder
       ~Tree.iterLeavesPreOrder
       ~Tree.iterNodes
       ~Tree.iterNodesNoRoot
       ~Tree.iterPostOrder
       ~Tree.iterPreOrder
       ~Tree.nextNode

    See also :class:`~p4.node.Node` methods that do similar things starting from a given node.

**Copy**

    .. autosummary::
       :nosignatures:

       ~Tree.dupe
       ~Tree.copyToTree
       ~Tree.dupeSubTree

**In combination with Data and Model**

     .. autosummary::
        :nosignatures:

        ~Tree.calcLogLike
        ~Tree.optLogLike
        ~Tree.simulate
        ~Tree.getSiteLikes
        ~Tree.ancestralStateDraw
        ~Tree.bigXSquaredSubM
        ~Tree.compStatFromCharFreqs
        ~Tree.compoTestUsingSimulations
        ~Tree.modelFitTests
        ~Tree.modelSanityCheck
        ~Tree.simsForModelFitTests
        ~Tree.getEuclideanDistanceFromSelfDataToExpectedComposition

**Setting a model**

     .. autosummary::
        :nosignatures:

        ~Tree.newComp
        ~Tree.newRMatrix
        ~Tree.newGdasrv
        ~Tree.setPInvar
        ~Tree.setRelRate
        ~Tree.setModelComponentOnNode
        ~Tree.setModelComponentsOnNodesRandomly
        ~Tree.setModelComponentsNNodes
        ~Tree.summarizeModelComponentsNNodes
        ~Tree.setNGammaCat
        ~Tree.setTextDrawSymbol


**Tree manipulation**

     .. autosummary::
        :nosignatures:

        ~Tree.addLeaf
        ~Tree.addNodeBetweenNodes
        ~Tree.addSibLeaf
        ~Tree.addSubTree
        ~Tree.allBiRootedTrees
        ~Tree.collapseNode
        ~Tree.collapseClade
        ~Tree.ladderize
        ~Tree.lineUpLeaves
        ~Tree.nni
        ~Tree.nni2
        ~Tree.pruneSubTreeWithoutParent
        ~Tree.pruneSubTreeWithParent
        ~Tree.randomSpr
        ~Tree.randomizeTopology
        ~Tree.reRoot
        ~Tree.reconnectSubTreeWithoutParent
        ~Tree.reconnectSubTreeWithParent
        ~Tree.removeEverythingExceptCladeAtNode
        ~Tree.removeNode
        ~Tree.removeAboveNode
        ~Tree.removeRoot
        ~Tree.renameForPhylip
        ~Tree.resolve
        ~Tree.resolvePolytomyAtNode
        ~Tree.restoreDupeTaxa
        ~Tree.restoreNamesFromRenameForPhylip
        ~Tree.rotateAround
        ~Tree.spr
        ~Tree.stripBrLens

**Misc**

     .. autosummary::
        :nosignatures:

        ~Tree.checkDupedTaxonNames
        ~Tree.checkSplitKeys
        ~Tree.checkTaxNames
        ~Tree.checkThatAllSelfNodesAreInTheTree
        ~Tree.inputTreesToSuperTreeDistances
        ~Tree.makeSplitKeys
        ~Tree.readBipartitionsFromPaupLogFile
        ~Tree.recalculateSplitKeysOfNodeFromChildren
        ~Tree.setNexusSets
        ~Tree.topologyDistance
        ~Tree.tvTopologyCompare
        ~Tree.patristicDistanceMatrix



    """

    from p4.tree_manip import node, rotateAround, reRoot, removeRoot, removeNode, removeAboveNode, collapseNode, collapseClade, pruneSubTreeWithoutParent, reconnectSubTreeWithoutParent, pruneSubTreeWithParent, reconnectSubTreeWithParent, addNodeBetweenNodes, allBiRootedTrees, ladderize, randomizeTopology, readBipartitionsFromPaupLogFile, renameForPhylip, restoreNamesFromRenameForPhylip, restoreDupeTaxa, lineUpLeaves, removeEverythingExceptCladeAtNode, dupeSubTree, addSubTree, addLeaf, addSibLeaf, subTreeIsFullyBifurcating, nni, nni2, checkThatAllSelfNodesAreInTheTree, spr, randomSpr, inputTreesToSuperTreeDistances, resolvePolytomyAtNode, resolve
    from p4.tree_optsim import __del__, deleteCStuff, _allocCStuff, setCStuff, _commonCStuff, calcLogLike, optLogLike, optTest, simulate, ancestralStateDraw, getSiteLikes
    from p4.tree_model import data, model, _checkModelThing, newComp, newRMatrix, newGdasrv, setPInvar, setRelRate, setModelComponentOnNode, setModelThingsRandomly, setModelComponentsOnNodesRandomly, setModelComponentsNNodes, summarizeModelComponentsNNodes, setTextDrawSymbol, setNGammaCat, modelSanityCheck, setEmpiricalComps
    from p4.tree_write import patristicDistanceMatrix, tPickle, writeNexus, write, writePhylip, writeNewick, _getMcmcCommandComment, draw, textDrawList, eps
    from p4.tree_fit import simsForModelFitTests, modelFitTests, compoTestUsingSimulations, bigXSquaredSubM, compStatFromCharFreqs, getEuclideanDistanceFromSelfDataToExpectedComposition

    def __init__(self):
        self.fName = None   # The name of the file it came from
        self.name = None
        self.root = None
        self.nodes = []
        # nodeNums of nodes root -> tips   A numpy array
        self.preOrder = None
        # nodeNums of nodes tips -> root   A numpy array
        self.postOrder = None
        self.preAndPostOrderAreValid = 0
        # Usually weight is 1/N, so the reciprocal looks nicer
        self.recipWeight = None
        # self.weight = None               # Only for floating point weights,
        # so not usually ...
        # An ordered list.  self.taxNames is a property
        self._taxNames = []
        # A Data object.  self.data is a property
        self._data = None
        self.cTree = None                 # A pointer to a c-struct
        self.logLike = None
        self.partLikes = None
        # A Model object.  self.model is a property
        self._model = None
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

    @property
    def taxNames(self):
        """(property) taxNames"""
        return self._taxNames

    @taxNames.setter
    def taxNames(self, theTaxNames):
        gm = ['Tree.taxNames']
        if not isinstance(theTaxNames, list):
            gm.append("You can only set property 'taxNames' to a list.")
            gm.append("Got attempt to set to '%s'" % theTaxNames)
            raise P4Error(gm)
        self._taxNames = theTaxNames
        if theTaxNames:
            self.checkTaxNames()
    
    @taxNames.deleter
    def taxNames(self):
        gm = ['Tree.taxNames']
        gm.append("    Caught an attempt to delete self.taxNames, but")
        gm.append("self.taxNames is a property, so you can't delete it.")
        gm.append("But you can set it to an empty list if you like.")
        raise P4Error(gm)


    def _setTaxNamesFromLeaves(self):
        tax = []
        for n in self.iterNodes():
            if n.isLeaf and n.name:
                tax.append(n.name)
            # This next line should not be needed, as root leaves should be
            # leaves.
            # terminal root that has a taxName
            elif n == self.root and n.name and n.getNChildren() < 2:
                tax.append(n.name)
            tax.sort()
        self._taxNames = tax
        if self._taxNames:
            self.checkTaxNames()

    def _getNTax(self):
        # We can't rely on len(self.taxNames), cuz it might not exist.
        # if hasattr(self, '_nTax') and self._nTax:
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
    """(property) nTax"""

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
        gm.append(
            "self.nInternalNodes is a property, so you shouldn't do that.")
        raise P4Error(gm)

    def _delNInternalNodes(self):
        self._nInternalNodes = -1

    nInternalNodes = property(_getNInternalNodes, _setNInternalNodes, _delNInternalNodes)
    """(property) nInternalNodes"""

    ##################################################
    ##################################################

    def dupe(self):
        """Duplicates self, but with no c-pointers.  And no data object.

        If there is a model, it is duped.

        Returns a copy of self.
        """

        storedData = None
        if self.data:
            storedData = self.data
            # We don't want to have to copy a big data object, now do we?
            self.data = None
        dupe = copy.deepcopy(self)
        # print 'Tree.dupe()   self.root=%s, dupe.root=%s' % (self.root,
        # dupe.root)

        # Delete cPointers
        for n in dupe.nodes:
            if n.cNode:
                n.cNode = None
        if dupe.cTree:
            dupe.cTree = None
        if dupe.model and dupe.model.cModel:
            dupe.model.cModel = None
        if storedData:
            self.data = storedData
        return dupe

    def parseNexus(self, flob, translationHash=None, doModelComments=0):
        """Start parsing nexus format newick tree description.

        From just after the command word 'tree', to the first paren of
        the Newick part of the tree.

        Args:
            flob: an open file or file-like object

            translationHash (dict): associates short names or numbers with long
            proper taxon names

            doModelComments (bool): whether to parse p4-specific model command
            comments in the tree description

        Returns:
            None

        """

        gm = ['Tree.parseNexus()']  # re-defined below
        if 0:
            print('Tree.parseNexus() translationHash = %s' % translationHash)
            print('    doModelComments = %s (nParts)' % doModelComments)
            print('    flob type is %s' % type(flob))

        tok = safeNextTok(flob, 'Tree.parseNexus()')
        # print 'parseNexus() tok = %s' % tok
        tok = p4.func.nexusUnquoteName(tok)
        if tok == '*':
            print(gm[0])
            print("    Ignoring '*' in tree description")
            tok = safeNextTok(flob, 'Tree.parseNexus()')
        if not p4.func.nexusCheckName(tok):
            gm.append("Bad tree name: '%s'" % tok)
            raise P4Error(gm)
        self.name = tok
        # print "got name: '%s'" % tok
        # print "%s" % tok
        gm = ["Tree.parseNexus() '%s'" % self.name]  # re-defining
        tok = safeNextTok(flob, gm[0])
        if tok != '=':
            gm.append("Tree name must be followed by '='")
            raise P4Error(gm)

        # Generally this is the beginning of the newick tree
        # description.  But we have to look ahead to see if there is a
        # weight comment.
        savedPos = flob.tell()
        while 1:
            beforeSafeNextTokPosn = flob.tell()
            #print('mnop var.nexus_getAllCommandComments = %s' % var.nexus_getAllCommandComments)

            tok = safeNextTok(flob, gm[0])          # skips [&U] if there is one
            #print("parseNexus: tok after '=' is '%s'.  Before pos %i, after pos %i" % (
            #    tok, beforeSafeNextTokPosn, flob.tell()))

            # This next bit will only happen if either var.nexus_getWeightCommandComments
            # or var nexus_getAllCommandComments is set.
            if tok[0] == '[':
                self.getWeightCommandComment(tok)
            elif tok == '(':
                flob.seek(beforeSafeNextTokPosn)
                #print("wxy var.nexus_getAllCommandComments is %s" % var.nexus_getAllCommandComments)
                self.parseNewick(flob, translationHash, doModelComments)
                # self._initFinish()
                break
            elif tok == ';':
                gm.append("Got ';' before any tree description.")
                raise P4Error(gm)
            elif tok[0] in string.ascii_letters + string.digits + '_' + "'":
                flob.seek(savedPos, 0)
                self.parseNewick(flob, translationHash, doModelComments)
                # self._initFinish()
                break
            else:
                gm.append('Expecting a newick tree description.')
                raise P4Error(gm)
        self._initFinish()
        # print 'finished Tree.parseNexus()'

    def getWeightCommandComment(self, tok):
        if 0:
            print('var.nexus_getWeightCommandComments = %s' % var.nexus_getWeightCommandComments)
            print('var.nexus_getAllCommandComments = %s' % var.nexus_getAllCommandComments)
            print("Got comment '%s' (type %s), checking if it is a 'weight' comment." % (tok, type(tok)))
        gm = ["Tree.getWeightCommandComment()"]
        cFlob = io.StringIO(tok)
        cFlob.seek(1)  # The [
        cTok = nextTok(cFlob)
        if not cTok:
            # print "no cTok -- returning nothing"
            return
        lowCTok = cTok.lower()
        if lowCTok in ['&r', '&u']:
            # print "got %s -- returning nothing" % cTok
            return
        if lowCTok != '&w':
            gm.append('Expecting a weight comment.  Got %s' % tok)
            raise P4Error(gm)
        cTok = nextTok(cFlob)
        # It might be a float, or the more usual 1/something
        if ("." in cTok):
            # A float?
            try:
                self.weight = float(cTok)
            except:
                gm.append("I can't grok '%s' in weight comment %s" %
                          (cTok, tok))
                raise P4Error(gm)

        # Should check for scientific notation?

        else:
            try:
                theNumerator = int(cTok)
                if theNumerator != 1:
                    gm.append(
                        'Expecting a numerator 1 in weight comment %s' % tok)
                    raise P4Error(gm)
                # print 'got theNumerator %i' % theNumerator
            except ValueError:
                gm.append('Expecting a numerator 1 in weight comment %s' % tok)
                raise P4Error(gm)
            cTok = nextTok(cFlob)
            if cTok == '/':
                cTok = safeNextTok(cFlob, 'Getting weight comment %s' % tok)
                try:
                    self.recipWeight = int(cTok)
                except ValueError:
                    gm.append('Bad denominator in weight comment %s' % tok)
                    raise P4Error(gm)
            elif cTok == ']':
                # self.recipWeight = theNumerator # ie 1, might as well leave
                # it as None
                pass
            else:
                gm.append("I can't grok '%s' in weight comment %s" %
                          (cTok, tok))
                raise P4Error(gm)
        cFlob.close()
        # print 'got recipWeight = %s' % self.recipWeight


# def printStack(self, theStack):  # only used for debugging parseNewick()
# print 'stack = ',
# for n in theStack:
# print "%i['%s'] " % (n.nodeNum, n.name),
# print ''

    def parseNewick(self, flob, translationHash, doModelComments=0):
        """Parse Newick tree descriptions.

        This is stack-based, and does not use recursion.
        """

        if 0:
            print('parseNewick here. doModelComments=%s' % doModelComments)
            print("    translationHash=%s, self.taxNames=%s" % (translationHash, self.taxNames))
            print("    flob type is %s, pos is %i" % (type(flob), flob.tell()))
            print("wyy var.nexus_getAllCommandComments is %s" % var.nexus_getAllCommandComments)
            print(f"    var.punctuation is {var.punctuation}")

        if self.name:
            gm = ["Tree.parseNewick(), tree '%s'" % self.name]
        else:
            gm = ['Tree.parseNewick()']

        if hasattr(flob, 'name') and flob.name:
            self.fName = flob.name
            gm[0] += ", File %s" % self.fName


        stack = []
        isAfterParen = 1  # to start, even tho its not true
        isAfterComma = 0
        parenNestLevel = 0
        lastPopped = None
        
        tok = nextTok(flob)          # skip over [&U] if there is one.
        # print("xx Got tok %s" % tok)
        if not tok:
            return
        isQuotedTok = False
        if tok.startswith("'"):
            isQuotedTok = True
        # Should generally be the opening paren, except if its a single-node
        # tree.
        tok = p4.func.nexusUnquoteName(tok)

        if doModelComments:
            # Turn on var.nexus_getAllCommandComments in order 
            # to be able to read model info on nodes, eg [& c0.1]
            # restore at end
            # We need this after the safeNextTok() above, in case there is a [&U] to be ignored.
            savedP4Nexus_getAllCommandComments = var.nexus_getAllCommandComments
            var.nexus_getAllCommandComments = 1

        while tok != ';':
            # print("top of loop tok '%s', isQuotedTok=%s, tok[0] is '%s'" % (tok, isQuotedTok, tok[0]))
            if tok == '(':
                # print "Got '(': new node (%i)." % len(self.nodes)
                if not (isAfterParen or isAfterComma):
                    gm.append(
                        'Got badly-placed paren, not after a paren or comma.')
                    raise P4Error(gm)
                newNode = Node()
                if doModelComments:
                    for pNum in range(doModelComments):
                        newNode.parts.append(NodePart())
                        newNode.br.parts.append(NodeBranchPart())

                # self.printStack(stack)
                if len(stack):
                    newNode.parent = stack[-1]
                    if newNode.parent.leftChild == None:
                        newNode.parent.leftChild = newNode
                    else:
                        newNode.parent.rightmostChild().sibling = newNode
                else:
                    if len(self.nodes) == 0:
                        self.root = newNode
                        # Sometimes. Generally not true-- corrected at the end.
                        newNode.isLeaf = 1
                    else:
                        gm.append('Something is wrong. Stack is empty.')
                        gm.append('Extra paren?')
                        raise P4Error(gm)
                newNode.nodeNum = len(self.nodes)
                self.nodes.append(newNode)
                stack.append(newNode)
                isAfterParen = 1
                parenNestLevel += 1

            elif tok == ',':
                if isAfterParen:
                    gm.append('Got comma after paren.')
                    raise P4Error(gm)
                elif isAfterComma:
                    gm.append('Got comma after comma.')
                    raise P4Error(gm)
                # self.printStack(stack)
                try:
                    lastPopped = stack.pop()
                except IndexError:
                    gm.append('Empty stack.  Out of place comma?')
                    raise P4Error(gm)
                isAfterComma = 1
                if len(stack) == 0:
                    gm.append('Empty stack.  Out of place comma?')
                    raise P4Error(gm)

            elif tok == ')':
                try:
                    lastPopped = stack.pop()
                except IndexError:
                    gm.append('Empty stack.  Out of place unparen?')
                    raise P4Error(gm)

                isAfterParen = 0
                isAfterComma = 0
                parenNestLevel = parenNestLevel - 1
                if parenNestLevel < 0:
                    gm.append('Unmatched unparen.')
                    raise P4Error(gm)
                if len(stack) == 0 and len(self.nodes) > 1:
                    gm.append('Empty stack.  Out of place unparen?')
                    raise P4Error(gm)

            elif tok[0] in string.ascii_letters or tok[0] in string.digits or tok[0] in var.nexus_safeChars \
                 or isQuotedTok or tok[0] in ['_', '#', '\\', '/', '"', '(', ')']:
                # A single-node tree, not ()aName, rather just aName.
                # print(f"parseNewick() here ACE tok={tok}")
                if len(self.nodes) == 0:
                    isAfterParen = 1
                if not (isAfterParen or isAfterComma):
                    # Probably a name of an internal node.
                    # print(f"parseNewick() here ACG tok={tok} isAfterParen={isAfterParen} len(stack)={len(stack)}")
                    if len(stack):
                        # print(f"parseNewick() here ACI stack[-1].name={stack[-1].name}")
                        # if stack[-1].isLeaf and stack[-1].name != '(':
                        if stack[-1].name:
                            if not var.newick_allowSpacesInNames:
                                # a second name after a node name, eg (A foo, B)   =>foo is bad
                                # or eg (A, B)foo bar    => bar is bad
                                gm.append("Badly placed token '%s'." % tok)
                                gm.append("Appears to be a second node name, after '%s'" % stack[-1].name)
                                gm.append('Missing comma maybe?  Or punctuation or spaces in an unquoted name?')
                                gm.append("To allow reading Newick (or Nexus) with spaces, ")
                                gm.append("turn var.newick_allowSpacesInNames on")
                                raise P4Error(gm)
                            else:
                                stack[-1].name += ' '
                                stack[-1].name += tok
                        else:
                            # Usually this...
                            # print("parseNewick() here ACK naming node %i as '%s'" % (stack[-1].nodeNum, tok))
                            # We allow bad names on internal nodes, ie we do
                            # not nexusCheckName(tok)
                            stack[-1].name = tok
                            # We may have just named the root node, and it may be a leaf, and perhaps should have a translation.
                            # But we do not know that yet.  Checked at the end.

                    else:    # len(stack) == 0
                        if lastPopped and lastPopped.name == None:  # ()A
                            # print "naming lastPopped node %i with '%s'" %
                            # (lastPopped.nodeNum, tok)
                            lastPopped.isLeaf = 1
                            #lastPopped.label = tok
                            lastPopped.name = tok
                        else:
                            gm.append("Badly placed token '%s' in tree description." % tok)
                            raise P4Error(gm)

                else:
                    # A new terminal node.
                    if tok[0] in string.ascii_letters or tok[0] in ['_']:
                        if translationHash and tok in translationHash:
                            # print 'got key %s, val is %s' % (tok,
                            # translationHash[tok])
                            tok = translationHash[tok]

                    elif tok[0] in string.digits:
                        if var.nexus_allowAllDigitNames:
                            if translationHash and tok in translationHash:
                                # print 'got key %s, val is %s' % (tok,
                                # translationHash[tok])
                                tok = translationHash[tok]
                        else:
                            try:
                                tok = int(tok)
                                if translationHash and repr(tok) in translationHash:
                                    tok = translationHash[repr(tok)]
                                elif translationHash and repr(tok) not in translationHash:
                                    gm.append("There is a 'translation' for this tree, but the")
                                    gm.append("number '%i' in the tree description" % tok)
                                    gm.append('is not included in that translate command.')
                                    raise P4Error(gm)
                                elif self.taxNames:
                                    try:
                                        tok = self.taxNames[tok - 1]
                                    except IndexError:
                                        gm.append("Can't make sense out of token '%s' for a new terminal node." % tok)
                                        gm.append('There is no translate command, and the taxNames does not')
                                        gm.append('have a value for that number.')
                                        raise P4Error(gm)
                                else:
                                    gm.append("We have a taxon name '%s', composed only of numerals." % tok)
                                    gm.append(" ")
                                    gm.append('The Nexus format allows tree specifications with no')
                                    gm.append('translate command to use integers to refer to taxa.')
                                    gm.append('That is possible because in a proper Nexus file the')
                                    gm.append('taxa are defined before the trees.  P4, however, does')
                                    gm.append('not require definition of taxa before the trees, and in')
                                    gm.append('the present case no definition was made.  Deal with it.')
                                    raise P4Error(gm)
                            except ValueError:
                                if translationHash and repr(tok) in translationHash:
                                    tok = translationHash[repr(tok)]
                                # else:  # starts with a digit, but it is not an int.
                                #    gm.append('Problem token %s' % tok)
                                #    raise P4Error(gm)

                    #print("Got terminal node '%s'" % tok)
                    newNode = Node()
                    if doModelComments:
                        for pNum in range(doModelComments):
                            newNode.parts.append(NodePart())
                            newNode.br.parts.append(NodeBranchPart())

                    newNode.isLeaf = 1
                    if p4.func.nexusCheckName(tok):
                        newNode.name = tok
                        # print 'got newNode.name = %s' % tok
                    else:
                        gm.append("Bad name '%s'" % tok)
                        raise P4Error(gm)
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
                # It might be a multi-token operation.  Accumulate tok's in
                # theNum
                theNum = nextTok(flob)
                if not theNum:
                    gm.append('Tree description ended with a colon.  Bad!')
                    raise P4Error(gm)
                # print "  Got token after colon:  '%s'" % theNum
                if theNum == '-' or theNum == '+':
                    tok = nextTok(flob)
                    # print "  Got tok: '%s' after '%s'" % (tok, theNum)
                    if not tok:
                        gm.append('Trying to deal with a branch length.')
                        gm.append("It didn't work, tho.")
                        gm.append("Got this after colon:  '%s'" % theNum)
                        gm.append('followed by nothing.')
                        raise P4Error(gm)
                    theNum += tok
                try:
                    # If it is a simple number like 0.123 or -23, then we are
                    # finished.
                    # Won't work if it ends in 'e'
                    stack[-1].br.len = float(theNum)
                    #print('  Successfully got br.len %f' % stack[-1].br.len)
                except ValueError:
                    # The first bit after the colon is hopefully something like
                    # +1.23e
                    if theNum[-1] not in ['e', 'E']:
                        gm.append(
                            'Trying to deal with a branch length after a colon, but I am totally confused.')
                        gm.append("Can't make sense out of '%s'" % theNum)
                        raise P4Error(gm, 'newick_badBranchLength')
                    try:
                        float(theNum[:-1])
                    except ValueError:
                        gm.append(
                            'Trying to deal with a branch length after a colon, but I am totally confused.')
                        gm.append("Can't make sense out of '%s'" % theNum)
                        raise P4Error(gm, 'newick_badBranchLength')

                    # Now we are sure that the first bit *is* something like +1.23e
                    # We do not allow spaces after the 'e', so we do not use nextTok().
                    # That introduces a bug, where comments inserted in the number don't get ignored.   <<== unfixed bug!
                    # The first thing must be a '+' or a '-'.
                    c = flob.read(1)
                    if not c:
                        gm.append(
                            'Trying to deal with a branch length, possibly in scientific notation.')
                        gm.append(
                            "Got '%s' after the colon, but then nothing." % theNum)
                        raise P4Error(gm)
                    if c not in ['+', '-']:
                        gm.append(
                            'Trying to deal with a branch length, possibly in scientific notation.')
                        gm.append("Got '%s' after the colon." % theNum)
                        gm.append(
                            "Expecting a '+' or '-' after that (no spaces allowed).")
                        gm.append("Got '%s'." % c)
                        raise P4Error(gm)
                    # Accumulate characters in 'theExp'.  We need at least one
                    # digit.
                    theExp = c
                    c = flob.read(1)
                    if not c:
                        gm.append('Trying to deal with a branch length, possibly in scientific notation.')
                        gm.append("Got '%s%s' after the colon, but then nothing." % (theNum, theExp))
                        raise P4Error(gm)
                    if c not in string.digits:
                        gm.append("Trying to deal with a branch length, possibly in scientific notation.")
                        gm.append("Got '%s%s' after the colon." % (theNum, theExp))
                        gm.append('Expecting one or more digits.')
                        gm.append("Got '%s'" % c)
                        raise P4Error(gm)
                    theExp += c
                    # So we got one good digit.  Are there any more?
                    while 1:
                        positionBeforeRead = flob.tell()
                        c = flob.read(1)
                        if not c:
                            gm.append('Trying to deal with a branch length, possibly in scientific notation.')
                            gm.append("Got '%s%s' after the colon, but then nothing." % (theNum, theExp))
                            raise P4Error(gm)
                        # We got something.  If its a digit, add it to
                        # theExp.  If its anything else, back up one
                        # space and then break
                        if c in string.digits:
                            theExp += c
                        else:
                            flob.seek(positionBeforeRead)
                            break
                    #print("  At this point, theNum='%s' and theExp='%s'" % (theNum, theExp))
                    try:
                        # print "  Trying to see if theExp '%s' can be
                        # converted to an int." % theExp
                        int(theExp)
                        try:
                            theBrLen = float(theNum + theExp)
                            # print '  Successfully got br.len %g (from %s%s)'
                            # % (theBrLen, theNum, theExp)
                            stack[-1].br.len = theBrLen
                        except ValueError:
                            gm.append('Trying to deal with a branch length, possibly in scientific notation.')
                            gm.append("It didn't work, tho.")
                            gm.append("Got these after colon: '%s' and '%s'" % (theNum, theExp))
                            gm.append('And they could not be converted to an exponential float.')
                            raise P4Error(gm)
                    except ValueError:
                        gm.append('Trying to deal with a branch length, possibly in scientific notation.')
                        gm.append("It didn't work, tho.")
                        gm.append("Got these after colon: '%s' and '%s'." % (theNum, theExp))
                        gm.append('And the latter does not appear to be an int.')
                        raise P4Error(gm)

            elif tok[0] == '[':
                # print("openSquareBracket.  Got tok '%s'" % tok)
                # if doModelComments is set, it should be set to nParts.
                if doModelComments:
                    # eg [& c0.1 r0.0]
                    n = stack[-1]
                    # print('got comment %s, node %i' % (tok, n.nodeNum))
                    cFlob = io.StringIO(tok)
                    cFlob.seek(2)
                    tok2 = safeNextTok(cFlob)
                    while 1:
                        if tok2 == ']':
                            break
                        elif tok2[0] in ['c', 'r', 'g']:
                            ending = tok2[1:]
                            splitEnding = ending.split('.')
                            try:
                                firstNum = int(splitEnding[0])
                                secondNum = int(splitEnding[1])
                            except ValueError:
                                gm.append('Bad command comment %s' % tok)
                                raise P4Error(gm)
                            if tok2[0] == 'c':
                                n.parts[firstNum].compNum = secondNum
                            if tok2[0] == 'r':
                                n.br.parts[firstNum].rMatrixNum = secondNum
                            if tok2[0] == 'g':
                                n.br.parts[firstNum].gdasrvNum = secondNum
                        else:
                            gm.append('Bad command comment %s' % tok)
                            raise P4Error(gm)
                        tok2 = safeNextTok(cFlob)
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
                    # print mySupportString
                    myNode.name = mySupportString
                elif var.nexus_readBeastTreeCommandComments:
                    n = stack[-1]
                    i = 2
                    while i < (len(tok) - 1):
                        j = i
                        inBraces = False
                        while 1:
                            j += 1
                            # print tok[j]
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
                        # print(substring)
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
                        assert isinstance(theVal, (float, tuple, bool))
                        n.__setattr__(theNameString, theVal)
                        i = j + 1
                elif 0:
                    # branch colours, after the branch length, eg '[&!color=#8caaf4]'
                    if tok.startswith("[&!color="):
                        myNode = stack[-1]
                        theColour = tok[9:-1]
                        #print(theColour)
                        myNode.br.colour = theColour

                    
                    
            else:
                gm.append("I can't make sense of the token '%s'" % tok)
                if len(tok) == 1:
                    if tok[0] in var.punctuation:
                        gm.append("The token is in var.punctuation. If you don't think it should")
                        gm.append("be, you can modify what p4 thinks that punctuation is.")
                        gm.append("So you might do this:")
                        gm.append("var.punctuation = var.phylip_punctuation")
                        gm.append("(or use your own definition -- see var.py)")
                        gm.append("read('yourWackyTreeFile.phy')")
                        gm.append("That might work.")
                    if tok[0] not in var.punctuation:
                        gm.append("The token is not in your current var.punctuation.")

                #gm.append("tok[0] is '%s'" % tok[0])
                raise P4Error(gm)

            sTok = safeNextTok(flob, 'Tree init, reading tree string')
            if sTok.startswith("'"):
                isQuotedTok = True
            else:
                isQuotedTok = False
            tok = p4.func.nexusUnquoteName(sTok)
            #print('got tok for next round = %s' % tok)
            # This is the end of the "while tok != ';':" loop

        # print '\n*** Stack len = %i ***' % len(stack)

        if parenNestLevel > 0:
            gm.append('Unmatched paren.')
            raise P4Error(gm)
        elif parenNestLevel < 0:
            gm.append('Unmatched unparen.')
            raise P4Error(gm)

        if len(stack) == 0:
            if len(self.nodes) == 1:
                pass
            else:
                gm.append(
                    "Got an oddly-placed ';' in the tree %s description." % self.name)
                self.dump(treeInfo=0, nodeInfo=1)
                raise P4Error(gm)
        elif len(stack) > 1:
            gm.append(
                "Got an oddly-placed ';' in the tree %s description." % self.name)
            # gm.append('(stack len = %i)' % len(stack)
            #self.dump(tree=0, node=1)
            raise P4Error(gm)

        if self.root.leftChild and self.root.leftChild.sibling:  # usually this
            self.root.isLeaf = 0
        else:
            # The root is a leaf
            assert self.root.isLeaf
            if self.root.name:
                if translationHash and self.root.name in translationHash:
                    self.root.name = translationHash[self.root.name]

        # Should a root on a stick be a leaf?  If it is just for
        # display purposes, then it should be ok to not be a leaf.
        # But if you are going to re-Root, then it will cause trouble.
        # So by default, a root on a stick should be a leaf.  I think.

        # Hopefully if you are dealing with this, then you know what
        # you are doing and what you want, and how to modify things to
        # get it.

        # Uncomment this next line to make it always non-leaf, even if it is a
        # leaf.

        # self.root.isLeaf = 0 # do this to always have a non-leaf
        # root-on-a-stick   <-- Potential Trouble!!!
        self.root.br = None
        # self.draw()
        #self.dump(tree=0, node=1, treeModel=0)

        if doModelComments:
            # restore the value of var.nexus_getAllCommandComments, which was
            # saved above.
            var.nexus_getAllCommandComments = savedP4Nexus_getAllCommandComments

    def _initFinish(self):

        if self.name:
            gm = ["Tree._initFinish()   tree '%s'" % self.name]
        else:
            gm = ['Tree._initFinish()']

        # Checking for duped taxon names used to be here, but it was
        # pushed further ahead, so that self.taxNames can be corrected
        # also, if need be.  At this point, self does not have
        # taxNames set.

        # Check that all terminal nodes have names
        for item in self.nodes:
            if item.isLeaf:
                # print 'leaf name %s' % item.name
                if not item.name:
                    if item == self.root:
                        if var.warnAboutTerminalRootWithNoName:
                            print('Tree._initFinish()')
                            print('    Non-fatal warning: the root is terminal, but has no name.')
                            print('    This may be what you wanted.  Or not?')
                            print('    (To get rid of this warning, turn off var.warnAboutTerminalRootWithNoName)')
                    else:
                        gm.append('Got a terminal node with no name.')
                        raise P4Error(gm)

        self.preOrder = numpy.array(
            [var.NO_ORDER] * len(self.nodes), numpy.int32)
        self.postOrder = numpy.array(
            [var.NO_ORDER] * len(self.nodes), numpy.int32)

        if len(self.nodes) > 1:
            self.setPreAndPostOrder()

    def checkDupedTaxonNames(self):

        # Called by p4.func._tryToReadNexusFile() and p4.func._tryToReadPhylipFile()
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
                    gm.append(
                        "Got duplicated taxon (lowercased) name '%s'." % loName)
                    gm.append(
                        'Since var.doRepairDupedTaxonNames is not turned on, p4 will not fix duplications.')
                    gm.append('To repair duplications verbosely, set ')
                    gm.append('var.doRepairDupedTaxonNames = 1')
                    gm.append('To repair duplications silently, set')
                    gm.append('var.doRepairDupedTaxonNames = 2')
                    raise P4Error(gm)
                hasDupedName = 1
                break

        if hasDupedName:
            # print self.name
            if var.allowDupedTaxonNames:
                # more hacking ...
                if var.allowDupedTaxonNames == 2:  # ie silently.
                    pass
                else:
                    complainedAlready = []
                    for loName in loNames:
                        if loNames.count(loName) > 1 and loName not in complainedAlready:
                            if self.name:
                                print("        Tree %s. Duped tax name (lowercased) '%s'" % (
                                    self.name, loName))
                            else:
                                print("        Duped tax name (lowercased) '%s'" % loName)
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
                                            print("        Tree %s. Changing '%s' to '%s'" % (
                                                self.name, n.name, newName))
                                        else:
                                            print("        Changing '%s' to '%s'" % (n.name, newName))
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
            shows which model component number goes on which node.
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
            print("Tree '%s' dump" % self.name)
        else:
            print('Tree dump.  No name.')
        if self.fName:
            print("    From file '%s'" % self.fName)
        else:
            print("    From an unknown file, or no file.")
        if self.root:
            print('    Node %i is root' % self.root.nodeNum)
        else:
            print('    There is no root')
        if self.recipWeight:
            print('    The tree recipWeight is %s' % self.recipWeight)
        else:
            print('    There is no recipWeight')
        print('    There are %i nodes' % len(self.nodes))
        terminals = 0
        for i in self.nodes:
            if i.isLeaf:
                terminals += 1
        print('        of which %i are terminal nodes' % terminals)
        if self.data:
            print('    There is a data object, with %i parts.' % self.data.nParts)
        else:
            print('    There is no data object.')

        if self.data:
            print('    The data came from file(s):')
            for a in self.data.alignments:
                if a.fName:
                    print('        %s' % a.fName)

        if self.model:
            print('    There is a model object, with %i parts.' % self.model.nParts)
            if self.model.cModel:
                print('    model.cModel is %i' % self.model.cModel)
            else:
                print('    There is no cModel.')
        else:
            print('    There is no model object.')

        if self.taxNames:
            print('    There is a taxNames list.')
        else:
            print('    There is no taxNames list.')
        if self.cTree:
            print('    cTree is %i' % self.cTree)
        else:
            print('    There is no cTree.')

    def _doNodeInfo(self):
        """Basic rubbish about nodes."""

        print('\n-------- nodes -----------------------------------------')
        print('%7s %6s %6s %6s %6s %7s %6s  %4s' % ('nodeNum', 'isLeaf', 'parent', 'leftCh',
                                                    'siblng', 'br.len', 'seqNum', 'name'))
        for n in self.nodes:
            print('%7s %6s' % (n.nodeNum, n.isLeaf), end=' ')
            if n.parent:
                print('%6s' % n.parent.nodeNum, end=' ')
            else:
                print('%6s' % 'None', end=' ')

            if n.leftChild:
                print('%6s' % n.leftChild.nodeNum, end=' ')
            else:
                print('%6s' % 'None', end=' ')

            if n.sibling:
                print('%6s' % n.sibling.nodeNum, end=' ')
            else:
                print('%6s' % 'None', end=' ')

            if n.br and (n.br.len or n.br.len == 0.0):
                print('%7.3f' % n.br.len, end=' ')
            else:
                print('%7s' % 'None', end=' ')

            if n.seqNum or n.seqNum == 0:
                print('%6s' % n.seqNum, end=' ')
            else:
                print('%6s' % 'None', end=' ')

            if n.name:
                print('  %s' % n.name)
            else:
                print('  %s' % 'None')
        print('--------------------------------------------------------\n')

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
            print('\n-------- more node stuff -------------------------------')
            print('%7s   %10s %10s %10s   %4s' % ('nodeNum', 'br.name', 'br.uName', 'br.support', 'name'))
            for n in self.nodes:
                print('%7s' % n.nodeNum, end=' ')

                if n.br and hasattr(n.br, 'name') and n.br.name:
                    print('%10s' % n.br.name, end=' ')
                else:
                    print('%10s' % '-', end=' ')

                if n.br and hasattr(n.br, 'uName') and n.br.uName:
                    print('%10s' % n.br.uName, end=' ')
                else:
                    print('%10s' % '-', end=' ')

                if n.br and hasattr(n.br, 'support') and n.br.support:
                    print('%10.4f' % n.br.support, end=' ')
                else:
                    print('%10s' % '-', end=' ')

                if n.name:
                    print('    %s' % n.name)
                else:
                    print('    %s' % 'None')
            print('--------------------------------------------------------\n')

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
            print('\n-------- even more node stuff --------------------------')
            print('%7s   %10s %14s %10s   %4s' % ('nodeNum', 'br.color', 'br.biRootCount', 'rootCount', 'name'))
            for n in self.nodes:
                print('%7s' % n.nodeNum, end=' ')

                if n.br and hasattr(n.br, 'color') and n.br.color:
                    print('%10s' % n.br.color, end=' ')
                else:
                    print('%10s' % '-', end=' ')

                if n.br and hasattr(n.br, 'biRootCount') and n.br.biRootCount:
                    print('%14s' % n.br.biRootCount, end=' ')
                else:
                    print('%14s' % '-', end=' ')

                if hasattr(n, 'rootCount') and n.rootCount:
                    print('%10s' % n.rootCount, end=' ')
                else:
                    print('%10s' % '-', end=' ')

                if n.name:
                    print('    %s' % n.name)
                else:
                    print('    %s' % 'None')
            print('--------------------------------------------------------\n')

    def _doNodeModelInfo(self):
        if not self.model:
            print('\n****** Node Model Info.  No model.')
            if not self.data:
                print('(no data attached, either)')
        else:
            print('\n****** Node Model Info.  nParts=%s' % self.model.nParts)
            if not self.data:
                print('no data')
            if not self.model.nParts:
                return

            print('\nComps in the nodes:')
            print(' %13s' % 'nodeNum', end=' ')
            for i in range(self.model.nParts):
                print(' %8s' % 'part%i' % i, end=' ')
            print('')
            for n in self.nodes:
                print(' %13i' % n.nodeNum, end=' ')
                for i in range(self.model.nParts):
                    print('%8i' % n.parts[i].compNum, end=' ')
                print('')

            print('\nrMatrices in the nodes:')
            print(' %13s' % 'nodeNum', end=' ')
            for i in range(self.model.nParts):
                print(' %8s' % 'part%i' % i, end=' ')
            print('')
            for n in self.iterNodesNoRoot():
                print(' %13i' % n.nodeNum, end=' ')
                for i in range(self.model.nParts):
                    print('%8i' % n.br.parts[i].rMatrixNum, end=' ')
                print('')

            print('\ngdasrvs in the nodes:')
            print(' %13s' % '', end=' ')
            for i in range(self.model.nParts):
                print(' %8s' % 'part%i' % i, end=' ')
            print('')
            print(' %13s' % 'nGammaCats ->', end=' ')
            for i in range(self.model.nParts):
                print('%8i' % self.model.parts[i].nGammaCat, end=' ')
            print('\n')

            print(' %13s' % 'nodeNum', end=' ')
            for i in range(self.model.nParts):
                print(' %8s' % 'part%i' % i, end=' ')
            print('')
            for n in self.iterNodesNoRoot():
                print(' %13i' % n.nodeNum, end=' ')
                for i in range(self.model.nParts):
                    print('%8i' % n.br.parts[i].gdasrvNum, end=' ')
                print('')

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
            # print "%s. There are %i taxSets." % (gm[0], len(self.nexusSets.taxSets))
            # Check that no taxSet name is a taxName
            lowSelfTaxNames = [txName.lower()
                               for txName in self.taxNames]
            for ts in self.nexusSets.taxSets:
                if ts.lowName in lowSelfTaxNames:
                    gm.append(
                        "Can't have taxSet names that are the same (case-insensitive) as a tax name")
                    gm.append(
                        "Lowercased taxSet name '%s' is the same as a lowcased taxName." % ts.name)
                    raise P4Error(gm)
            self.nexusSets.lowTaxNames = lowSelfTaxNames

            # If it is standard format,
            # convert triplets to numberTriplets, and then mask
            for ts in self.nexusSets.taxSets:
                if ts.format == 'standard':
                    ts.setNumberTriplets()
                    ts.setMask()
                elif ts.format == 'vector':
                    assert ts.mask
                    if len(ts.mask) != self.nTax:
                        gm.append("taxSet %s" % ts.name)
                        gm.append(
                            "It is vector format, but the length is wrong.")
                        gm.append(
                            "taxSet mask is length %i, but self nTax is %i" % (len(ts.mask), self.nTax))
                        raise P4Error(gm)
                else:
                    gm.append("taxSet %s" % ts.name)
                    gm.append("unknown format %s" % ts.format)
                    raise P4Error(gm)

                # Now set ts.taxNames from the mask
                ts.taxNames = [self.taxNames[i] for i,c in enumerate(ts.mask) if c == '1']


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
    #         raise P4Error(gm)
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
            raise P4Error(gm)
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
            raise P4Error(gm)
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
        preOrder = []  # nodeNum's
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
                # print 'leftChild: %i' % stack[-1].leftChild.nodeNum
                theNodeNum = stack[-1].leftChild.nodeNum
                stack.append(stack[-1].leftChild)
                preOrder.append(theNodeNum)
            elif stack[-1].sibling:
                # print 'sibling: %i' % stack[-1].sibling.nodeNum
                theNodeNum = stack[-1].sibling.nodeNum
                theSib = stack[-1].sibling
                # print '                 postOrder appending u %i' %
                # stack[-1].nodeNum
                postOrder.append(stack[-1].nodeNum)
                stack.pop()
                stack.append(theSib)
                preOrder.append(theNodeNum)
            else:
                # print '                 postOrder appending  v %i' %
                # stack[-1].nodeNum
                postOrder.append(stack[-1].nodeNum)
                stack.pop()
                if len(stack) == 0:
                    break
                # print '                 postOrder appending  w %i' %
                # stack[-1].nodeNum
                postOrder.append(stack[-1].nodeNum)
                theNode = stack.pop()
                while not theNode.sibling:
                    if len(stack) == 0:
                        break
                    # print '                 postOrder appending x %i' %
                    # stack[-1].nodeNum
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
                    raise P4Error(gm)

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
        if self.preOrder is None or self.postOrder is None or len(self.preOrder) != len(self.nodes) or len(self.postOrder) != len(self.nodes):
            self.preOrder = numpy.array(
                [var.NO_ORDER] * len(self.nodes), numpy.int32)
            self.postOrder = numpy.array(
                [var.NO_ORDER] * len(self.nodes), numpy.int32)

        # print "xx self.preOrder=%s" % self.preOrder
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
                    # print 'leftChild: %i  (%s)' %
                    # (stack[-1].leftChild.nodeNum, [n.nodeNum for n in stack])
                    theNodeNum = stack[-1].leftChild.nodeNum
                    stack.append(stack[-1].leftChild)
                    self.preOrder[preOrdIndx] = theNodeNum
                    preOrdIndx += 1
                elif stack[-1].sibling:
                    # print 'sibling: %i  (%s)' % (stack[-1].sibling.nodeNum,
                    # [n.nodeNum for n in stack])
                    theNodeNum = stack[-1].sibling.nodeNum
                    theSib = stack[-1].sibling
                    # print '                 postOrder appending u %i' %
                    # stack[-1].nodeNum
                    self.postOrder[postOrdIndx] = stack[-1].nodeNum
                    postOrdIndx += 1
                    stack.pop()
                    stack.append(theSib)
                    try:
                        self.preOrder[preOrdIndx] = theNodeNum
                    except IndexError:
                        gm.append("preOrdIndx=%s, theNodeNum=%i" %
                                  (preOrdIndx, theNodeNum))
                        gm.append("preOrder = %s" % self.preOrder)
                        raise P4Error(gm)
                    preOrdIndx += 1
                else:
                    # print '                 postOrder appending  v %i' %
                    # stack[-1].nodeNum
                    self.postOrder[postOrdIndx] = stack[-1].nodeNum
                    postOrdIndx += 1
                    stack.pop()
                    if len(stack) == 0:
                        break
                    # print '                 postOrder appending  w %i' %
                    # stack[-1].nodeNum
                    self.postOrder[postOrdIndx] = stack[-1].nodeNum
                    postOrdIndx += 1
                    theNode = stack.pop()
                    while not theNode.sibling:
                        if len(stack) == 0:
                            break
                        # print '                 postOrder appending x %i' %
                        # stack[-1].nodeNum
                        self.postOrder[postOrdIndx] = stack[-1].nodeNum
                        postOrdIndx += 1
                        theNode = stack.pop()
                    if len(stack) == 0:
                        break
                    if theNode.sibling:
                        stack.append(theNode.sibling)
                        # print "self.preOrder = %s, preOrdIndx=%i" %
                        # (self.preOrder, preOrdIndx)
                        self.preOrder[preOrdIndx] = theNode.sibling.nodeNum
                        preOrdIndx += 1
                    else:
                        gm.append('Problemo.')
                        gm.append('xxx got theNode %s' % theNode.nodeNum)
                        raise P4Error(gm)
        if 1:
            assert preOrdIndx == postOrdIndx
            # print "a preOrdIndx = %i, len(self.nodes) = %i" % (preOrdIndx,
            # len(self.nodes))
            if preOrdIndx != len(self.nodes):
                pOI = preOrdIndx
                for i in range(preOrdIndx, len(self.nodes)):
                    self.preOrder[i] = var.NO_ORDER
                    self.postOrder[i] = var.NO_ORDER
                    preOrdIndx += 1
                    postOrdIndx += 1
            # print "b preOrdIndx = %i, len(self.nodes) = %i" % (preOrdIndx,
            # len(self.nodes))
        assert preOrdIndx == len(self.nodes) and postOrdIndx == len(self.nodes)

    def iterPreOrder(self):
        """Node generator.  Assumes preAndPostOrderAreValid."""

        for i in self.preOrder:
            j = int(i)           # this speeds things up a lot!  Two-fold!
            if j != var.NO_ORDER:
                # yield self.nodes[int(i)]
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

        # self.draw()

#        allOnes = 2**(self.nTax) - 1
        # print 'nTax = %i, allOnes = %i' % (self.nTax, allOnes)

#        for n in self.iterInternalsPostOrder():
        for n in self.iterInternalsNoRootPostOrder():
            #            if n != self.root:
                # print 'doing node %s' % n.nodeNum
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
                    splitList.append([x, 0])
                n.br.rcList = [n.br.rc]

    def makeSplitKeys(self, makeNodeForSplitKeyDict=True):
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
        #raise P4Error(gm)
        if not self.taxNames:
            gm.append('No taxNames.')
            raise P4Error(gm)

        if not self.preAndPostOrderAreValid:
            self.setPreAndPostOrder()

        # self.draw()
        if makeNodeForSplitKeyDict:
            self.nodeForSplitKeyDict = {}

        allOnes = 2 ** (self.nTax) - 1
        # print 'nTax = %i, allOnes = %i' % (self.nTax, allOnes)

        for n in self.iterPostOrder():
            if n != self.root:
                # print 'doing node %s' % n.nodeNum

                if not n.leftChild:
                    try:
                        n.br.rawSplitKey = 1 << self.taxNames.index(
                            n.name)  # "<<" is left-shift
                    except ValueError:
                        gm.append('node.name %s' % n.name)
                        gm.append('is not in taxNames %s' % self.taxNames)
                        raise P4Error(gm)
                    # n.br.rawSplitKey = 1 << self.taxNames.index(n.name)  #
                    # "<<" is left-shift

                    if n.br.rawSplitKey == 1:
                        n.br.splitKey = allOnes - 1
                    else:
                        n.br.splitKey = n.br.rawSplitKey

                    # print 'upward leaf   node %s, rawSplitKey %s, splitKey
                    # %s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
                else:
                    childrenNums = self.getChildrenNums(n)
                    if len(childrenNums) == 1:
                        gm.append(
                            'Internal node has only one child.  That will make a duped split.')
                        raise P4Error(gm)
                    x = self.nodes[childrenNums[0]].br.rawSplitKey
                    for i in childrenNums[1:]:
                        y = self.nodes[i].br.rawSplitKey
                        x = x | y  # '|' is bitwise "OR".
                    n.br.rawSplitKey = x

                    # Make node.br.splitKey's in "standard form", ie
                    # without the first taxon, ie without a 1.  To do that
                    # we use the '&' operator, to bitwise "and" with 1.
                    # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey
                    # odd?"
                    if 1 & n.br.rawSplitKey:
                        # "^" is xor, a bit-flipper.
                        n.br.splitKey = allOnes ^ n.br.rawSplitKey
                    else:
                        n.br.splitKey = n.br.rawSplitKey
                    # print 'intern node %s, rawSplitKey %s, splitKey %s' %
                    # (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
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
                            print('%7s     %4s       %4s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey))
                    raise P4Error(gm)

            # Any duped splitKeys?  There will be if the tree is bi-Rooted.
            if 0:
                theKeys = []
                for n in self.iterNodesNoRoot():
                    theKeys.append(n.br.splitKey)
                for k in theKeys:
                    if theKeys.count(k) > 1:
                        gm.append('Duped splitKey %i.' % k)
                        for n in self.iterNodesNoRoot():
                            print('%7s     %4s       %4s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey))
                        raise P4Error(gm)

        if 0:
            print(gm[0])
            print(self.taxNames)
            print(f"root is node {self.root.nodeNum}")
            print('nodeNum  rawSplitKey  splitKey')
            # getSplitsFromKey would usually use self.nTax.  
            # Changed below to accommodate tempRoot for bi-rooted trees.
            for n in self.iterNodesNoRoot():
                print('%7s     %4s       %4s        %s' % (
                    n.nodeNum, n.br.rawSplitKey, n.br.splitKey, 
                    p4.func.getSplitStringFromKey(n.br.splitKey, len(self.taxNames))))

    def recalculateSplitKeysOfNodeFromChildren(self, aNode, allOnes):
        children = [n for n in aNode.iterChildren()]
        x = children[0].br.rawSplitKey
        for n in children[1:]:
            x = x | n.br.rawSplitKey   # '|' is bitwise "OR".
        aNode.br.rawSplitKey = x
        # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
        if 1 & aNode.br.rawSplitKey:
            # "^" is xor, a bit-flipper.
            aNode.br.splitKey = allOnes ^ aNode.br.rawSplitKey
        else:
            aNode.br.splitKey = aNode.br.rawSplitKey

    def checkSplitKeys(self, useOldName=False, glitch=True, verbose=True):
        gm = ['Tree.checkSplitKeys()']

        allOnes = 2 ** (self.nTax) - 1
        # print 'nTax = %i, allOnes = %i' % (self.nTax, allOnes)

        isBad = False
        for n in self.iterPostOrder():
            if n != self.root:
                # print 'doing node %s' % n.nodeNum

                if not n.leftChild:
                    # A long int, eg 1, has no upper limit on its value
                    try:
                        if useOldName:
                            rawSplitKey = 1 << self.taxNames.index(
                                n.oldName)  # "<<" is left-shift
                        else:
                            rawSplitKey = 1 << self.taxNames.index(
                                n.name)  # "<<" is left-shift
                    except ValueError:
                        if useOldName:
                            gm.append('node.name %s' % n.oldName)
                        else:
                            gm.append('node.name %s' % n.name)
                        gm.append('is not in taxNames %s' % self.taxNames)
                        raise P4Error(gm)
                    # n.br.rawSplitKey = 1 << self.taxNames.index(n.name)  #
                    # "<<" is left-shift

                    if rawSplitKey == 1:
                        splitKey = allOnes - 1
                    else:
                        splitKey = rawSplitKey

                    # print 'upward leaf   node %s, rawSplitKey %s, splitKey
                    # %s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
                else:
                    childrenNums = self.getChildrenNums(n)
                    if len(childrenNums) == 1:
                        gm.append(
                            'Internal node has only one child.  That will make a duped split.')
                        raise P4Error(gm)
                    x = self.nodes[childrenNums[0]].br.rawSplitKey
                    for i in childrenNums[1:]:
                        y = self.nodes[i].br.rawSplitKey
                        x = x | y  # '|' is bitwise "OR".
                    rawSplitKey = x

                    # Make node.br.splitKey's in "standard form", ie
                    # without the first taxon, ie without a 1.  To do that
                    # we use the '&' operator, to bitwise "and" with 1.
                    # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey
                    # odd?"
                    if 1 & rawSplitKey:
                        # "^" is xor, a bit-flipper.
                        splitKey = allOnes ^ rawSplitKey
                    else:
                        splitKey = rawSplitKey
                    # print 'intern node %s, rawSplitKey %s, splitKey %s' %
                    # (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
                if n.br.rawSplitKey != rawSplitKey:
                    print("checkSplitKeys node %2i rawSplitKey: existing %s, calculated %s" % (n.nodeNum, n.br.rawSplitKey, rawSplitKey))
                    isBad = True
                if n.br.splitKey != splitKey:
                    print("checkSplitKeys node %2i splitKey: existing %s, calculated %s" % (n.nodeNum, n.br.splitKey, splitKey))
                    isBad = True
        if glitch and isBad:
            raise P4Error(gm)
        if verbose and not isBad:
            print("checkSplitKeys().  ok")

    def taxSetIsASplit(self, taxSetName):
        """Asks whether a nexus taxset is a split in the tree.

        Args:
            taxSetName (str): The name of the taxset. Case does not matter.

        Returns:
            Node: the node in self if the taxset is a split, or else None

        """

        gm = ["Tree.taxSetIsASplit(%s)" % taxSetName]
        assert self.nexusSets
        assert self.taxNames
        assert self.nexusSets.taxSets
        lowArgTaxSetName = taxSetName.lower()
        theTS = None
        for ts in self.nexusSets.taxSets:
            if ts.lowName == lowArgTaxSetName:
                theTS = ts
                break
        assert theTS, "Could not find the taxSet named %s" % taxSetName
        # theTS.dump()
        assert theTS.mask

        # Make sure that self splitKey has been set. Check the first internal node.
        needsDoing = False
        for n in self.iterInternalsNoRoot():
            if not n.br.splitKey:
                needsDoing = True
            break
        if needsDoing:
            self.makeSplitKeys()

        rawSplitKey = 0
        for i in range(len(theTS.mask)):
            # print i, theTS.mask[i]
            if theTS.mask[i] == '1':
                rawSplitKey += (1 << i)
        # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
        if 1 & rawSplitKey:
            allOnes = 2 ** (self.nTax) - 1
            splitKey = allOnes ^ rawSplitKey  # "^" is xor, a bit-flipper.
        else:
            splitKey = rawSplitKey
        # print "got splitKey %s" % splitKey

        for n in self.nodes:
            if n.br and not n.isLeaf:
                assert n.br.splitKey, "Needs node.br.splitKey on all internal nodes."
                # print "  %2i  %s  %s" % (n.nodeNum, n.br.splitKey, splitKey)
                if n.br.splitKey == splitKey:
                    # self.draw()
                    # return n.nodeNum
                    return n
        # Was -1 before, when n.nodeNum was returned for hits.  Now a node is
        # returned.
        return None

    def checkTaxNames(self):
        """Check that all taxNames are in the tree, and vice versa."""

        # If var.allowTreesWithDifferingTaxonSets is set to True we will not check
        # the taxnames. This is primarily used to read trees for supertree and
        # leafstability calcuations.

        # Peter comments that Tobias added this, but it is not needed,
        # and messes other things up -- so comment out until it is
        # sorted.
        # if var.allowTreesWithDifferingTaxonSets:
        #     return

        # self.taxNames is a property, so when it is set, it calls this method

        if self.name:
            gm = ['Tree.checkTaxNames()   tree %s' % self.name]
        else:
            gm = ['Tree.checkTaxNames()']
        if self.fName:
            gm.append("From file '%s' " % self.fName)

        if not self.taxNames:
            gm.append('No taxNames.')
            raise P4Error(gm)

        tax = []
        for n in self.iterNodes():
            if n.isLeaf and n.name:
                tax.append(n.name)
            # This next line should not be needed, as root leaves should be
            # leaves.
            # terminal root that has a taxName
            elif n == self.root and n.name and n.getNChildren() < 2:
                tax.append(n.name)

        # print 'tax from tree = %s' % tax
        # print 'self.taxNames = %s' % self.taxNames

        # Check for same number of taxa
        if len(tax) != len(self.taxNames):
            # self.draw()
            # self.dump(node=1)
            gm.append('Mismatch.  Length of self.taxNames is wrong.')
            gm.append('The tree has %i leaves with names, but len(self.taxNames) = %i' % (
                len(tax), len(self.taxNames)))
            gm.append('leaves on the tree = %s' % tax)
            gm.append('self.taxNames = %s' % self.taxNames)
            gm.append('symmetric_difference is %s' %
                      set(tax).symmetric_difference(set(self.taxNames)))
            # print '    Setting invalid taxNames to an empty list.'
            #self.taxNames = []
            raise P4Error(gm, 'tree_badTaxNamesLength')

        # Check for mis-matched taxNames
        isBad = 0
        taxSet = set(tax)
        selfTaxNamesSet = set(self.taxNames)

        s = taxSet.difference(selfTaxNamesSet)
        if len(s):
            isBad = 1
            print(gm[0])
            print('TaxName mismatch between the tree and self.taxNames.')
            print('These taxa are found in the tree but not in self.taxNames:')
            print(s)
        s = selfTaxNamesSet.difference(taxSet)
        if len(s):
            isBad = 1
            print(gm[0])
            print('TaxName mismatch between the tree and self.taxNames.')
            print('These taxa are found in the self.taxNames but not in the tree:')
            print(s)

        if isBad:
            raise P4Error(gm, 'tree_taxNamesMisMatch')

    ################################################
    # Copy and Verify

    def copyToTree(self, otherTree):

        gm = ['Tree.copyToTree()']

        if len(self.nodes) != len(otherTree.nodes):
            gm.append('Different number of nodes.')
            raise P4Error(gm)

        # node relations (parent, child, sib)
        for nNum in range(len(self.nodes)):
            selfNode = self.nodes[nNum]
            otherNode = otherTree.nodes[nNum]

            # parent
            if selfNode.parent:
                otherNode.parent = otherTree.nodes[selfNode.parent.nodeNum]
            else:
                # print "otherNode.parent", otherNode.parent
                otherNode.parent = None

            # leftChild
            if selfNode.leftChild:
                otherNode.leftChild = otherTree.nodes[
                    selfNode.leftChild.nodeNum]
            else:
                # print "otherNode.leftChild", otherNode.leftChild
                otherNode.leftChild = None

            # sibling
            if selfNode.sibling:
                otherNode.sibling = otherTree.nodes[selfNode.sibling.nodeNum]
            else:
                # print "otherNode.sibling", otherNode.sibling
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
                otherTree.nodes[nNum].flag = self.nodes[nNum].flag
                otherTree.nodes[nNum].br.splitKey = self.nodes[nNum].br.splitKey
                otherTree.nodes[nNum].br.rawSplitKey = self.nodes[nNum].br.rawSplitKey

        if 0:
            print("Tree.copyToTree()")
            for nNum in range(len(self.nodes)):
                if self.nodes[nNum] != self.root:
                    print(f"{nNum:3} {self.nodes[nNum].br.splitKey} {otherTree.nodes[nNum].br.splitKey}")
                else:
                    print(f"{nNum:3} is root")

        # model usage numbers
        if self.model:
            for nNum in range(len(self.nodes)):
                selfNode = self.nodes[nNum]
                otherNode = otherTree.nodes[nNum]
                if selfNode.parts and otherNode.parts:
                    for pNum in range(self.model.nParts):
                        otherNode.parts[pNum].compNum = selfNode.parts[pNum].compNum
                        if selfNode != self.root:
                            if selfNode.br.parts and otherNode.br.parts:
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

    def verifyIdentityWith(self, otherTree, doSplitKeys=False, doMore=False):
        """For MCMC debugging.  Verifies that two trees are identical."""

        complaintHead = '\nTree.verifyIdentityWith()'  # keep

        if len(self.nodes) != len(otherTree.nodes):
            print(complaintHead)
            print('    Different number of nodes.')
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
                print(complaintHead)
                print('    Node %i, relations differ.' % nNum)
                self.write()
                otherTree.write()
                return var.DIFFERENT

        if self.root.nodeNum != otherTree.root.nodeNum:
            print(complaintHead)
            print('    Roots differ.')
            return var.DIFFERENT

        # brLens, lenChanged, and node.flag. and splitKeys
        if 0:
            print("Tree.verifyIdentityWith() flag")
            for nNum in range(len(self.nodes)):
                if self.nodes[nNum] != self.root:
                    print(f"{nNum:3} {self.nodes[nNum].flag} {otherTree.nodes[nNum].flag}")
                else:
                    print(f"{nNum:3} is root")


        for nNum in range(len(self.nodes)):
            if self.nodes[nNum] != self.root:
                # if self.nodes[nNum].br.len != otherTree.nodes[nNum].br.len:
                if math.fabs(self.nodes[nNum].br.len - otherTree.nodes[nNum].br.len) > 1.e-8:
                    print(complaintHead)
                    print('    BrLens differ.')
                    return var.DIFFERENT
                if self.nodes[nNum].br.lenChanged != otherTree.nodes[nNum].br.lenChanged:
                    print(complaintHead)
                    print('    br.lenChanged differs.')
                    return var.DIFFERENT
                if self.nodes[nNum].flag != otherTree.nodes[nNum].flag:
                    print(complaintHead)
                    print('    flag differs, nodeNum %i.  %s vs %s' % (nNum, self.nodes[nNum].flag, otherTree.nodes[nNum].flag))
                    return var.DIFFERENT
                if doSplitKeys:
                    if self.nodes[nNum].br.splitKey != otherTree.nodes[nNum].br.splitKey:
                        print(complaintHead)
                        print(f'    SplitKeys differ. nodeNum {nNum}')
                        return var.DIFFERENT
                    if self.nodes[nNum].br.rawSplitKey != otherTree.nodes[nNum].br.rawSplitKey:
                        print(complaintHead)
                        print(f'    rawSplitKeys differ. nodeNum {nNum}')
                        return var.DIFFERENT

        if self.model:
            # model usage numbers
            isBad = 0
            for pNum in range(self.model.nParts):
                for nNum in range(len(self.nodes)):
                    selfNode = self.nodes[nNum]
                    otherNode = otherTree.nodes[nNum]
                    if selfNode.parts and otherNode.parts:
                        if selfNode.parts[pNum].compNum != otherNode.parts[pNum].compNum:
                            isBad = 1
                        if self.nodes[nNum] != self.root:
                            if selfNode.br.parts and otherNode.br.parts:
                                if selfNode.br.parts[pNum].rMatrixNum != otherNode.br.parts[pNum].rMatrixNum:
                                    isBad = 1
                                elif selfNode.br.parts[pNum].gdasrvNum != otherNode.br.parts[pNum].gdasrvNum:
                                    isBad = 1

                    if isBad:
                        print(complaintHead)
                        print('    Node %i, model usage info does not match.' % nNum)
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
                print(complaintHead)
                print('    Pre- or postOrder do not match.')
                return var.DIFFERENT

            if self.nInternalNodes != otherTree.nInternalNodes:
                print(complaintHead)
                print('    nInternalNodes differ.')
                return var.DIFFERENT

            # partLikes
            for pNum in range(self.model.nParts):
                # if otherTree.partLikes[pNum] != self.partLikes[pNum]:
                myDiff = math.fabs(otherTree.partLikes[pNum] - self.partLikes[pNum])
                if myDiff > 1.e-5:
                    print(complaintHead)
                    print("    partLikes differ by %.8f (%g).  (%.8f, (%g) %.8f (%g)" % (
                        myDiff, myDiff,
                        otherTree.partLikes[pNum], otherTree.partLikes[pNum], self.partLikes[pNum], self.partLikes[pNum]))
                    return var.DIFFERENT

        if doMore:  # some more
            for nNum in range(len(self.nodes)):
                selfNode = self.nodes[nNum]
                otherNode = otherTree.nodes[nNum]
                if selfNode.nodeNum != otherNode.nodeNum:
                    print(complaintHead)
                    print('    nodeNum differs')
                    return var.DIFFERENT
                if selfNode.seqNum != otherNode.seqNum:
                    print(complaintHead)
                    print('    seqNum differs')
                    return var.DIFFERENT
                if selfNode.name != otherNode.name:
                    print(complaintHead)
                    print('    name differs')
                    return var.DIFFERENT
                if selfNode.isLeaf != otherNode.isLeaf:
                    print(complaintHead)
                    print('    isLeaf differs')
                    return var.DIFFERENT

        return var.SAME

    ############################################


    def isFullyBifurcating(self, verbose=False, biRoot=False):
        """Returns True if the tree is fully bifurcating.  Else False. 

        If arg biRoot is True, then it is required that the tree be
        biRoot'ed (ie have a bifurcating root)

        If arg biRoot is False, then biRoot'ed trees are not allowed,
        but trees rooted on a leaf are allowed.

        """

        rootNChildren = self.root.getNChildren()
        if rootNChildren > 3:
            if verbose:
                print("isFullyBifurcating() returning False, due to root with %i children." % rootNChildren)
            return False
        if not biRoot:
            if rootNChildren == 2:
                if verbose:
                    print("isFullyBifurcating() returning False, due to root with %i children." % rootNChildren)
                return False
            
            elif rootNChildren == 1 or rootNChildren == 3:   # rooting on a leaf is OK
                if rootNChildren == 1:
                    if not self.root.isLeaf:
                        if verbose:
                            print("isFullyBifurcating()", end=" ")
                            print("returning False, because the single-child root is not a leaf.")
                        return False

                for n in self.iterInternalsNoRoot():
                    if n.leftChild and n.leftChild.sibling:
                        if n.leftChild.sibling.sibling:
                            if verbose:
                                print("isFullyBifurcating()", end=" ")
                                print("returning False, due to node %i having 3 or more children." % n.nodeNum)
                            return False
                    else:
                        if verbose:
                            print("isFullyBifurcating()", end=" ")
                            print("returning False, due to non-leaf node %i having 1 or fewer children." % n.nodeNum)
                        return False
                return True
            else:
                raise P4Error("Tree.isFullyBifurcating() --- this should not happen; fix me.")

        if biRoot:
            # root on a leaf not allowed on a biRoot tree, as the root is explicitly not a leaf.
            if rootNChildren != 2:    
                if verbose:
                    print("isFullyBifurcating(biRoot=True)", end=" ")
                    print("returning False, due to root with %i children." % rootNChildren)
                return False

            for n in self.iterInternalsNoRoot():
                if n.leftChild and n.leftChild.sibling:
                    if n.leftChild.sibling.sibling:
                        if verbose:
                            print("isFullyBifurcating(biRoot=True)", end=" ")
                            print("returning False, due to node %i having 3 or more children." % n.nodeNum)
                        return False
                else:
                    if verbose:
                        print("isFullyBifurcating(biRoot=True)", end=" ")
                        print("returning False, due to non-leaf node %i having 1 or fewer children." % n.nodeNum)
                    return False
            return True

    # These next two are for the eTBR implementation that I got from Jason
    # Evans' Crux.  Thanks Jason!
    def getDegree(self, nodeSpecifier):
        n = self.node(nodeSpecifier)
        if n.isLeaf:
            if n.parent:
                return 1
            else:
                assert n == self.root
                if n.leftChild:   # a tree-on-a-stick
                    return 1
                else:             # a single isolated node
                    return 0
        else:
            #assert n.leftChild
            deg = 1  # the leftChild
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
            print("*=" * 25)
            self.draw()
            gm = ["Tree.nextNode() spoke=%i, hub=%i" %
                  (spoke.nodeNum, hub.nodeNum)]
            gm.append(
                "Need to have either spoke.parent == hub or hub == spoke.")
            raise P4Error(gm)

    def topologyDistance(self, tree2, metric='sd', resetSplitKeySet=True):
        """Compares the topology of self with tree2.

        The complete list of metrics is given in
        var.topologyDistanceMetrics

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
        subsequent sd topologyDistance calculation will be wrong, as
        it uses previous splits.  So then you need to
        'resetSplitKeySet'.

        The 'scqdist' metric calculates quartet distances.  The code
        was written by Anders Kabell Kristensen for his Masters degree
        at Aarhus University, 2010.
        http://www.cs.au.dk/~dalko/thesis/ It has two versions -- a
        pure Python version (that needs scipy) that I do not include
        here, and a fast C++ version, that I wrapped in python.  Its
        speedy -- the 'sc' in 'scqdist' is for 'sub-cubic', ie better
        than O(n^3).

        I have also incorporated the tqDist v1.0.2 code, from 2014,
        also for quartet distance calculations, from Christian N.  S.
        Pedersen and his group at BiRC in Aarhus.  See
        https://users-cs.au.dk/cstorm/software/tqdist/ 2014.
        It is available here via the 'tqdist' metric.

        """

        gm = ['Tree.topologyDistance()']

        if metric not in var.topologyDistanceMetrics:
            gm.append("Got a request for unknown metric '%s'" % metric)
            gm.append("The 'metric' arg should be one of %s" %
                      var.topologyDistanceMetrics)
            raise P4Error(gm)
        if metric == 'scqdist':  # no need for taxNames
            try:
                import p4.scqdist as scqdist
            except ImportError:
                gm.append("Could not find the 'scqdist' module needed for this metric.")
                gm.append("See the instructions for making it in the p4 source, in the Qdist directory.")
                raise P4Error(gm)
            tsSelf = self.writeNewick(toString=True)
            tsTree2 = tree2.writeNewick(toString=True)
            return scqdist.qdist(tsSelf, tsTree2)

        elif metric == 'tqdist':  # no need for taxNames
            try:
                import p4.pytqdist as pytqdist
            except ImportError:
                gm.append("Could not find the 'tqdist' module needed for this metric.")
                gm.append("See the instructions for making it in the p4 source, in the tqDist directory.")
                raise P4Error(gm)
            tsSelf = self.writeNewick(toString=True, spaceAfterComma=False)
            tsTree2 = tree2.writeNewick(toString=True, spaceAfterComma=False)
            return pytqdist.qdistFromStrings(tsSelf, tsTree2)

        if not self.taxNames or not tree2.taxNames:
            gm.append("This method requires taxNames to be set.")
            raise P4Error(gm)
        if self.taxNames != tree2.taxNames:
            gm.append("The taxNames are different for the two trees.")
            gm.append("Self:  %s" % self.taxNames)
            gm.append("tree2: %s" % tree2.taxNames)
            raise P4Error(gm)
        if (self.root.getNChildren() == 2 or tree2.root.getNChildren() == 2) and (metric in ['wrf', 'bld']):
            gm.append('One of the input trees has a bifurcating root.')
            gm.append(
                'Weighted tree distance calculations do not work on bi-rooted trees.')
            raise P4Error(gm)

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
            self.splitKeySet = set(
                [n.br.splitKey for n in self.iterNodesNoRoot()])

        if not hasattr(tree2, 'splitKeySet'):
            tree2.makeSplitKeys()
            tree2.splitKeySet = set(
                [n.br.splitKey for n in tree2.iterNodesNoRoot()])

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
            return len(selfHasButTree2DoesNot), len(tree2HasButSelfDoesNot)

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
                # print '%s - %s' % (self.splitKeyHash[k].br.len,
                # tree2.splitKeyHash[k].br.len)
                theSum += abs(self.splitKeyHash[k].br.len -
                              tree2.splitKeyHash[k].br.len)
            for k in selfHasButTree2DoesNot:
                # print 'x %s' % self.splitKeyHash[k].br.len
                theSum += self.splitKeyHash[k].br.len
            for k in tree2HasButSelfDoesNot:
                # print 'y %s' % tree2.splitKeyHash[k].br.len
                theSum += tree2.splitKeyHash[k].br.len
            return theSum
        elif metric == 'bld':
            theSum = 0.0
            for k in self.splitKeySet.intersection(tree2.splitKeySet):
                theDiff = self.splitKeyHash[
                    k].br.len - tree2.splitKeyHash[k].br.len
                theSum += theDiff * theDiff
            for k in selfHasButTree2DoesNot:
                theSum += self.splitKeyHash[k].br.len * \
                    self.splitKeyHash[k].br.len
            for k in tree2HasButSelfDoesNot:
                theSum += tree2.splitKeyHash[k].br.len * \
                    tree2.splitKeyHash[k].br.len
            # print 'branch score =', theSum
            return math.sqrt(theSum)

    ##################################################

    def tv(self):
        """Tree Viewer. Show the tree in a gui window.

        Needs Tkinter.

        If you have nexus taxsets defined, you can show them.
        """
        from p4.btv import TV
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
        from p4.btv import BTV
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
            print("The trees are the same. No tv.")
            return
        # for sk in self.splitKeyHash.keys():
        #    if sk not in treeB.splitKeyHash:
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

        from p4.btv import TV
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

    ##################################################

    def drawTopologyCompare(self, treeB, showNodeNums=False):
        """Graphically show topology differences between two trees.

        The two trees (self and treeB) are drawn as text, with differences
        highlighted with a thick branch.

        The taxNames need to be set, and need to be the same for both
        trees.

        """

        sd = self.topologyDistance(treeB)
        if sd == 0:
            print("The trees are the same. ")
            return

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

        for spl in selfHasButTreeBDoesnt:
            n = self.splitKeyHash[spl]
            n.br.textDrawSymbol = '='
        for spl in treeBHasButSelfDoesnt:
            n = treeB.splitKeyHash[spl]
            n.br.textDrawSymbol = '='

        self.draw(showNodeNums=showNodeNums)
        treeB.draw(showNodeNums=showNodeNums)


    def getNodeOnReferenceTreeCorrespondingToSelfRoot(self, refTree, verbose=True):
        """Find, on a ref tree, a node corresponding to where the self root is.

        This works for both bifurcating and non-bifurcating roots (of
        self).

        The refTree should not be bi-rooted.

        To facilitate doing a lot of (self) trees, a counter in the
        corresponding node in the refTree is incremented.  The refTree
        is rooted on the first taxon.  For bi-rooted trees, the
        node.br.biRootCount is incremented, and for non-bi-rooted
        trees the node.rootCount is incremented.  The root counts
        should therefore be immune to reRoot()'ing.

        """

        gm = ["Tree.getNodeOnReferenceTreeCorrespondingToSelfRoot()"]
        assert self is not refTree
        assert self.taxNames
        assert refTree.taxNames
        assert self.taxNames == refTree.taxNames

        isBiRoot = False
        # sr = self.root; srnChildren = nChildren of self.root
        srnChildren = self.root.getNChildren()
        if srnChildren == 1:
            gm.append("Self root has only one child; does not work")
            raise P4Error(gm)
        elif srnChildren == 2:
            isBiRoot = True
            for n in refTree.iterNodesNoRoot():
                if not hasattr(n.br, "biRootCount"):
                    n.br.biRootCount = 0
        else:
            for n in refTree.iterNodes():
                if not hasattr(n, "rootCount"):
                    n.rootCount = 0

        #self.draw()
        refTree.reRoot(self.taxNames[0])
        #refTree.draw()

        self.makeSplitKeys()  # default makeNodeForSplitKeyDict=True
        refTree.makeSplitKeys()

        for ch in self.root.iterChildren():
            if ch.br.rawSplitKey & 1:
                break
        # ch clade has the first taxon, ch.parent is self.root
        #print("self root child %i clade has the first taxon" % ch.nodeNum)
        refNode = refTree.nodeForSplitKeyDict.get(ch.br.splitKey)
        if refNode:
            if isBiRoot:
                refNode.br.biRootCount += 1
                if verbose:
                    #print("When the refTree is rooted on the first taxon,")
                    #print("the bi-root is on the branch on refTree node %i" % refNode.nodeNum)
                    print("Incrementing refTree node %i br.biRootCount by 1" % refNode.nodeNum)

            else:
                refNode.rootCount += 1
                if verbose:
                    #print("When the refTree is rooted on the first taxon,")
                    #print("the root is at refTree node %i" % refNode.nodeNum)
                    print("Incrementing refTree node %i rootCount by 1" % refNode.nodeNum)

            return refNode
        else:
            if verbose:
                print("There is no node in the refTree corresponding to the self.root")
            return None
                
                

    def readPhyloXmlFile(self, fName, verbose=False):
        """Start with an empty Tree, read in a phyloxml file"""

        assert not self.nodes, "The tree should be empty, with no nodes"
        import xml.etree.ElementTree as ET
        xAll = ET.parse(fName)
        for el in xAll.iter():
            if '}' in el.tag:
                el.tag = el.tag.split('}', 1)[1]  # strip all namespaces

        xroot = xAll.getroot()
        xphylogeny = xroot.find('phylogeny')
        assert xphylogeny
        rootclade = xphylogeny.find('clade')
        assert rootclade

        def recurseClade(theClade):
            global currentNode
            for el in theClade:
                #print(el.tag)
                if el.tag == 'clade':
                    n = Node()
                    n.nodeNum = len(self.nodes)
                    el.attrib['p4Node'] = n
                    n.parent = theClade.attrib['p4Node']
                    if n.parent:
                        if n.parent.leftChild == None:
                            n.parent.leftChild = n
                        else:
                            rmCh = n.parent.rightmostChild()
                            assert rmCh, "Can't get rightmostChild"
                            rmCh.sibling = n

                    self.nodes.append(n)
                    currentNode = n
                    recurseClade(el)
                elif el.tag == "branch_length":
                    currentNode.br.len = float(el.text)
                elif el.tag == 'name':
                    #print(el.text)
                    currentNode.name = el.text
                    if verbose:
                        print("Got node name", el.text)
                elif el.tag == 'width':
                    currentNode.br.width = float(el.text)
                    if verbose:
                        print("Got node width, as node.br.width", el.text)
                elif el.tag == 'color':
                    thiscol = {}
                    for el2 in el:
                        thiscol[el2.tag] = int(el2.text)
                    cols = ['red', 'green', 'blue']
                    for col in cols:
                        assert thiscol.get(col)
                    myHexString = ""
                    for col in cols:
                        bit = thiscol.get(col)
                        hx = hex(bit)[2:].upper()
                        myHexString += hx
                    #print(myHexString)
                    currentNode.br.hexColor = myHexString
                    if verbose:
                        print("Got node color (as attribute node.br.hexColor)", currentNode.br.hexColor)
                else:
                    print("Fixme. Unhandled/unused element with tag %s" % el.tag)


        n = Node()
        n.nodeNum = 0
        n.br = None
        self.root = n
        self.nodes.append(n)

        rootclade.attrib['p4Node'] = n

        recurseClade(rootclade)

        for n in self.nodes:
            if not n.leftChild:
                n.isLeaf = 1
        self.setPreAndPostOrder()

    def attachRooter(self):
        rNode = self.nodes[-1]
        assert rNode.name == 'tempRooter'
        assert rNode.nodeNum == var.NO_ORDER
        rtMostCh = self.root.rightmostChild()
        rtMostCh.sibling = rNode
        rNode.sibling = None
        rNode.parent = self.root
        rNode.nodeNum = len(self.nodes) - 1
        if self.taxNames:
            # taxnames can be shared among trees, so check whether it has already been added
            if rNode.name not in self.taxNames:
                self.taxNames.append(rNode.name)
        self.preAndPostOrderAreValid = False
        self.setPreAndPostOrder()
        #if self._nTax:
        #    self._nTax += 1
        return rNode


    def detachRooter(self):
        rNode = self.nodes[-1]
        assert rNode.parent == self.root
        safety = 0
        while self.root.leftChild.sibling.sibling != rNode:
            self.rotateAround(self.root)
            safety += 1
            if safety >= 5:
                self.draw()
                raise P4Error("Tree.detachRooter() Too many rotations.")

        assert self.root.leftChild.sibling.sibling == rNode
        self.root.leftChild.sibling.sibling = None
        rNode.parent = None
        rNode.nodeNum = var.NO_ORDER
        self.preAndPostOrderAreValid = False
        self.setPreAndPostOrder()
        if self.taxNames:
            # taxnames can be shared among trees, so check whether it is in
            if rNode.name == self.taxNames[-1]:
                self.taxNames.pop()
        #if self._nTax:
        #    self._nTax -= 1

    def isBiRoot(self):
        """Answers whether self has a bifurcating root"""

        rootNChildren = self.root.getNChildren()
        if rootNChildren == 2:
            return True
        return False

    def isTriRoot(self):
        """Answers whether self has a trifurcating root"""

        rootNChildren = self.root.getNChildren()
        if rootNChildren == 3:
            return True
        return False

    def isCompatibleWith(self, otherTree):
        """Determines whether self is compatible with otherTree

        Returns True or False.
        """

        assert otherTree.nTax == self.nTax
        assert otherTree.taxNames == self.taxNames

        # It is easier with sets, so make sets from splitKeys
        self.makeSplitKeys()
        for n in self.iterInternalsNoRoot():
            ss = set()
            for i in range(self.nTax):
                tester = 2 ** i
                if tester & n.br.splitKey:
                    ss.add(tester)
            n.br.splSet = ss

        otherTree.makeSplitKeys()
        for n in otherTree.iterInternalsNoRoot():
            ss = set()
            for i in range(self.nTax):
                tester = 2 ** i
                if tester & n.br.splitKey:
                    ss.add(tester)
            n.br.splSet = ss

        isCompatible = True
        for nA in self.iterInternalsNoRoot():
            for nB in otherTree.iterInternalsNoRoot():
                if len(nA.br.splSet.intersection(nB.br.splSet)) \
                   and len(nA.br.splSet.difference(nB.br.splSet)) \
                   and len(nB.br.splSet.difference(nA.br.splSet)):
                    isCompatible = False
                    break
            if not isCompatible:
                break
        # print(f"isCompatible is {isCompatible}")
        return isCompatible
