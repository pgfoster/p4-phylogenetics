import sys
import string
import types
import cStringIO
import math
import copy
import os
import func
import time
import glob
from var import var
from p4exceptions import P4Error
from node import Node, NodePart, NodeBranchPart
import nexustoken
from distancematrix import DistanceMatrix

import numpy
import pf
from model import Model
from data import Data
from alignment import Part
from node import NodeBranch, NodePart, NodeBranchPart
import random


def _fixFileName(fName):
    if fName.count('.') or fName.count(' '):
        fName = list(fName)
        for i in range(len(fName)):
            theChar = fName[i]
            if theChar == '.' or theChar == ' ':
                fName[i] = '_'
        fName = string.join(fName, '')
    return fName

longMessage1 = """

(Boring old) Chi-square test for compositional homogeneity
==========================================================

    The statistic is Sum[(Obs - Exp)^2 / Exp],
        where Exp comes from the data.
    The statistic is X^2 (X squared), in the sense used by
        Sokal & Rohlf in 'Biometry', to distinguish it from
        'chi square', which is a distribution, not a statistic.
    This is (mostly) the same as the test in PAUP, using the
        basefreq command.  There are small differences having
        to do with calculation of degrees of freedom.
    Significance is assessed by Chi-square.

"""

longMessage2 = """

Tree- and model-based composition fit test
==========================================

    The statistic is Sum[(Obs - Exp)^2 / Exp],
        like the statistic used in the Chi-squared test,
        except that Exp comes from the model, not from the data.
    I call the statistic X^2_m (X squared sub m) since it is like
        the X^2 statistic used in the chi-square test above, but uses
        the m subscript to say that the expected values come from the model.
    Significance is assessed by simulations on the tree and model.
    A critical point of 95% is used to decide if the data fits or not.
"""


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
        :nosignatures:

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

    def _setTaxNames(self, theTaxNames):
        gm = ['Tree._setTaxNames()']
        if type(theTaxNames) != type([]):
            gm.append("You can only set property 'taxNames' to a list.")
            gm.append("Got attempt to set to '%s'" % theTaxNames)
            raise P4Error(gm)
        self._taxNames = theTaxNames
        # and not var.allowTreesWithDifferingTaxonSets:  # Peter commented out
        # until it is sorted. Why is it here?
        if theTaxNames:
            self.checkTaxNames()

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

    def _delTaxNames(self):
        gm = ['Tree._delTaxNames()']
        gm.append("    Caught an attempt to delete self.taxNames, but")
        gm.append("self.taxNames is a property, so you can't delete it.")
        gm.append("But you can set it to an empty list if you like.")
        raise P4Error(gm)

    taxNames = property(
        lambda self: self._taxNames, _setTaxNames, _delTaxNames)

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

    nInternalNodes = property(
        _getNInternalNodes, _setNInternalNodes, _delNInternalNodes)

    ##################################################
    ##################################################

    def dupe(self):
        """Duplicates self, but with no c-pointers.  And no data object.

        If there is a model, it is duped.
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
        the Newick part of the tree."""

        gm = ['Tree.parseNexus()']  # re-defined below
        if 0:
            print 'Tree.parseNexus() translationHash = %s' % translationHash
            print '    doModelComments = %s (nParts)' % doModelComments
            print '    var.nexus_doFastNextTok = %s' % var.nexus_doFastNextTok

        if var.nexus_doFastNextTok:
            from nexustoken2 import nextTok, safeNextTok
        else:
            from nexustoken import nextTok, safeNextTok

        tok = safeNextTok(flob, 'Tree.parseNexus()')
        # print 'parseNexus() tok = %s' % tok
        tok = func.nexusUnquoteName(tok)
        if tok == '*':
            print gm[0]
            print "    Ignoring '*' in tree description"
            tok = safeNextTok(flob, 'Tree.parseNexus()')
        if not func.nexusCheckName(tok):
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
            tok = safeNextTok(flob, gm[0])
            # print "parseNexus: tok after '=' is '%s'" % tok

            # This next bit will only happen if either var.nexus_getWeightCommandComments
            # or var nexus_getAllCommandComments is set.
            if tok[0] == '[':
                self.getWeightCommandComment(tok)
            elif tok == '(':
                flob.seek(-1, 1)
                self.parseNewick(flob, translationHash, doModelComments)
                # self._initFinish()
                break
            elif tok == ';':
                gm.append("Got ';' before any tree description.")
                raise P4Error(gm)
            elif tok[0] in string.letters + string.digits + '_' + '\'':
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
            print 'var.nexus_getWeightCommandComments = %s' % var.nexus_getWeightCommandComments
            print 'var.nexus_getAllCommandComments = %s' % var.nexus_getAllCommandComments
            print "Got comment '%s', checking if it is a 'weight' comment." % tok
        gm = ["Tree.getWeightCommandComment()"]
        # python, not c, so I can use StringIO
        from nexustoken import nextTok, safeNextTok
        cFlob = cStringIO.StringIO(tok)
        cFlob.seek(1)  # The [
        cTok = nextTok(cFlob)
        if not cTok:
            # print "no cTok -- returning nothing"
            return
        lowCTok = string.lower(cTok)
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

        # print 'parseNewick here. var.nexus_doFastNextTok=%s' % var.nexus_doFastNextTok
        # print 'parseNewick here. doModelComments=%s' % doModelComments
        # print "parseNewick()  translationHash=%s, self.taxNames=%s" %
        # (translationHash, self.taxNames)

        if self.name:
            gm = ["Tree.parseNewick(), tree '%s'" % self.name]
        else:
            gm = ['Tree.parseNewick()']

        if hasattr(flob, 'name') and flob.name:
            self.fName = flob.name
            gm[0] += ", File %s" % self.fName

        if doModelComments:
            # restore at end
            savedP4Nexus_getAllCommandComments = var.nexus_getAllCommandComments
            var.nexus_getAllCommandComments = 1

        stack = []
        isAfterParen = 1  # to start, even tho its not true
        isAfterComma = 0
        parenNestLevel = 0
        lastPopped = None

        if var.nexus_doFastNextTok:
            from nexustoken2 import nextTok, safeNextTok
        else:
            from nexustoken import nextTok, safeNextTok

        tok = nextTok(flob)
        if not tok:
            return

        # Should generally be the opening paren, except if its a single-node
        # tree.
        tok = func.nexusUnquoteName(tok)
        while tok != ';':
            # print "top of loop tok '%s', tok[0] is '%s'" % (tok, tok[0])

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

            elif tok[0] in string.letters or tok[0] in string.digits or tok[0] == "'" or tok[0] in [
                    '_', '#', '\\', '/', '"', '(', ')']:
                # A single-node tree, not ()aName, rather just aName.
                if len(self.nodes) == 0:
                    isAfterParen = 1
                if not (isAfterParen or isAfterComma):
                    # Probably a name of an internal node.
                    if len(stack):
                        # if stack[-1].isLeaf and stack[-1].name != '(':
                        if stack[-1].name:
                            if not var.newick_allowSpacesInNames:
                                # a second name after a node name, eg (A foo, B)   =>foo is bad
                                # or eg (A, B)foo bar    => bar is bad
                                gm.append("Badly placed token '%s'." % tok)
                                gm.append(
                                    "Appears to be a second node name, after '%s'" % stack[-1].name)
                                gm.append(
                                    'Missing comma maybe?  Or punctuation or spaces in an unquoted name?')
                                gm.append(
                                    "To allow reading Newick (or Nexus) with spaces, ")
                                gm.append(
                                    "turn var.newick_allowSpacesInNames on")
                                raise P4Error(gm)
                            else:
                                stack[-1].name += ' '
                                stack[-1].name += tok
                        else:
                            # Usually this...
                            # print "naming node %i as '%s'" % (stack[-1].nodeNum, tok)
                            # We allow bad names on internal nodes, ie we do
                            # not nexusCheckName(tok)
                            stack[-1].name = tok

                    else:    # len(stack) == 0
                        if lastPopped and lastPopped.name == None:  # ()A
                            # print "naming lastPopped node %i with '%s'" %
                            # (lastPopped.nodeNum, tok)
                            lastPopped.isLeaf = 1
                            #lastPopped.label = tok
                            lastPopped.name = tok
                        else:
                            gm.append(
                                "Badly placed token '%s' in tree description." % tok)
                            raise P4Error(gm)

                else:
                    # A new terminal node.
                    if tok[0] in string.letters or tok[0] in ['_']:
                        if translationHash and translationHash.has_key(tok):
                            # print 'got key %s, val is %s' % (tok,
                            # translationHash[tok])
                            tok = translationHash[tok]

                    elif tok[0] in string.digits:
                        if var.nexus_allowAllDigitNames:
                            if translationHash and translationHash.has_key(tok):
                                # print 'got key %s, val is %s' % (tok,
                                # translationHash[tok])
                                tok = translationHash[tok]
                        else:
                            try:
                                tok = int(tok)
                                if translationHash and translationHash.has_key(`tok`):
                                    tok = translationHash[`tok`]
                                elif translationHash and not translationHash.has_key(`tok`):
                                    gm.append(
                                        "There is a 'translation' for this tree, but the")
                                    gm.append(
                                        "number '%i' in the tree description" % tok)
                                    gm.append(
                                        'is not included in that translate command.')
                                    raise P4Error(gm)
                                elif self.taxNames:
                                    try:
                                        tok = self.taxNames[tok - 1]
                                    except IndexError:
                                        gm.append(
                                            "Can't make sense out of token '%s' for a new terminal node." % tok)
                                        gm.append(
                                            'There is no translate command, and the taxNames does not')
                                        gm.append(
                                            'have a value for that number.')
                                        raise P4Error(gm)
                                else:
                                    gm.append(
                                        "We have a taxon name '%s', composed only of numerals." % tok)
                                    gm.append(" ")
                                    gm.append(
                                        'The Nexus format allows tree specifications with no')
                                    gm.append(
                                        'translate command to use integers to refer to taxa.')
                                    gm.append(
                                        'That is possible because in a proper Nexus file the')
                                    gm.append(
                                        'taxa are defined before the trees.  P4, however, does')
                                    gm.append(
                                        'not require definition of taxa before the trees, and in')
                                    gm.append(
                                        'the present case no definition was made.  Deal with it.')
                                    raise P4Error(gm)
                            except ValueError:
                                if translationHash and translationHash.has_key(`tok`):
                                    tok = translationHash[`tok`]
                                # else:  # starts with a digit, but it is not an int.
                                #    gm.append('Problem token %s' % tok)
                                #    raise P4Error(gm)

                    # print "Got terminal node '%s'" % tok
                    newNode = Node()
                    if doModelComments:
                        for pNum in range(doModelComments):
                            newNode.parts.append(NodePart())
                            newNode.br.parts.append(NodeBranchPart())

                    newNode.isLeaf = 1
                    if func.nexusCheckName(tok):
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
                    # print '  Successfully got br.len %f' % stack[-1].br.len
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
                        gm.append(
                            'Trying to deal with a branch length, possibly in scientific notation.')
                        gm.append(
                            "Got '%s%s' after the colon, but then nothing." % (theNum, theExp))
                        raise P4Error(gm)
                    if c not in string.digits:
                        gm.append(
                            "Trying to deal with a branch length, possibly in scientific notation.")
                        gm.append("Got '%s%s' after the colon." %
                                  (theNum, theExp))
                        gm.append('Expecting one or more digits.')
                        gm.append("Got '%s'" % c)
                        raise P4Error(gm)
                    theExp += c
                    # So we got one good digit.  Are there any more?
                    while 1:
                        c = flob.read(1)
                        if not c:
                            gm.append(
                                'Trying to deal with a branch length, possibly in scientific notation.')
                            gm.append(
                                "Got '%s%s' after the colon, but then nothing." % (theNum, theExp))
                            raise P4Error(gm)
                        # We got something.  If its a digit, add it to
                        # theExp.  If its anything else, back up one
                        # space and then break
                        if c in string.digits:
                            theExp += c
                        else:
                            flob.seek(-1, 1)
                            break
                    # print "  At this point, theNum='%s' and theExp='%s'" %
                    # (theNum, theExp)
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
                            gm.append(
                                'Trying to deal with a branch length, possibly in scientific notation.')
                            gm.append("It didn't work, tho.")
                            gm.append(
                                "Got these after colon: '%s' and '%s'" % (theNum, theExp))
                            gm.append(
                                'And they could not be converted to an exponential float.')
                            raise P4Error(gm)
                    except ValueError:
                        gm.append(
                            'Trying to deal with a branch length, possibly in scientific notation.')
                        gm.append("It didn't work, tho.")
                        gm.append(
                            "Got these after colon: '%s' and '%s'." % (theNum, theExp))
                        gm.append(
                            'And the latter does not appear to be an int.')
                        raise P4Error(gm)

            elif tok[0] == '[':
                # print "openSquareBracket.  Got tok '%s'" % tok
                # if doModelComments is set, it should be set to nParts.
                if doModelComments:
                    # eg [& c0.1 r0.0]
                    n = stack[-1]
                    # print 'got comment %s, node %i' % (tok, n.nodeNum)
                    cFlob = cStringIO.StringIO(tok)
                    cFlob.seek(2)
                    tok2 = nexustoken.safeNextTok(cFlob)
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
                        tok2 = nexustoken.safeNextTok(cFlob)
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
                        # print substring
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
                        assert type(theVal) in [
                            types.FloatType, types.TupleType, types.BooleanType]
                        n.__setattr__(theNameString, theVal)
                        i = j + 1

            else:
                gm.append("I can't make sense of the token '%s'" % tok)
                if len(tok) == 1:
                    if tok[0] in var.punctuation:
                        gm.append(
                            "The token is in var.punctuation. If you don't think it should")
                        gm.append(
                            "be, you can modify what p4 thinks that punctuation is.")
                        if var.nexus_doFastNextTok:
                            gm.append(
                                "But to do that you can't use nexus_doFastNextToken, so you ")
                            gm.append("need to turn that off, temporarily.  ")
                            gm.append(
                                "(var.nexus_doFastNextTok is currently on.)")
                            gm.append("So you might do this:")
                            gm.append("var.nexus_doFastNextTok = False ")
                            gm.append(
                                "var.punctuation = var.phylip_punctuation")
                            gm.append(
                                "(or use your own definition -- see Var.py)")
                            gm.append("read('yourWackyTreeFile.phy')")
                            gm.append("That might work.")
                        else:
                            gm.append(
                                "(Doing that does not work if nexus_doFastNextToken is turned on,")
                            gm.append(
                                "but var.nexus_doFastNextTok is currently off.)")
                            gm.append("So you might do this:")
                            gm.append(
                                "var.punctuation = var.phylip_punctuation")
                            gm.append(
                                "(or use your own definition -- see Var.py)")
                            gm.append("read('yourWackyTreeFile.phy')")
                            gm.append("That might work.")
                    if tok[0] not in var.punctuation:
                        gm.append(
                            "The token is not in your current var.punctuation.")
                        if var.nexus_doFastNextTok:
                            gm.append(
                                "var.nexus_doFastNextTok is currently on.")
                            gm.append(
                                "It uses a hard-wired list of punctuation, and so you may need to ")
                            gm.append(
                                "need to turn var.nexus_doFastNextTok off, temporarily.  ")
                            gm.append("So you might do this:")
                            gm.append("var.nexus_doFastNextTok = False ")
                            gm.append("read('yourWackyTreeFile.phy')")
                            gm.append("That might work.")
                        else:
                            gm.append(
                                "var.nexus_doFastNextTok is currently off.")

                #gm.append("tok[0] is '%s'" % tok[0])
                raise P4Error(gm)

            tok = func.nexusUnquoteName(
                safeNextTok(flob, 'Tree init, reading tree string'))
            # print 'got tok for next round = %s' % tok
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
                            print 'Tree._initFinish()'
                            print '    Non-fatal warning: the root is terminal, but has no name.'
                            print '    This may be what you wanted.  Or not?'
                            print '    (To get rid of this warning, turn off var.warnAboutTerminalRootWithNoName)'
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
            # print "%s. There are %i taxSets." % (gm[0], len(self.nexusSets.taxSets))
            # Check that no taxSet name is a taxName
            lowSelfTaxNames = [string.lower(txName)
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
                    # print ts.mask
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
        if type(self.preOrder) == types.NoneType or type(self.postOrder) == types.NoneType or len(self.preOrder) != len(self.nodes) or len(self.postOrder) != len(self.nodes):
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

#        allOnes = 2L**(self.nTax) - 1
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
        #raise P4Error(gm)
        if not self.taxNames:
            gm.append('No taxNames.')
            raise P4Error(gm)

        if not self.preAndPostOrderAreValid:
            self.setPreAndPostOrder()

        # self.draw()
        if makeNodeForSplitKeyDict:
            self.nodeForSplitKeyDict = {}

        allOnes = 2L ** (self.nTax) - 1
        # print 'nTax = %i, allOnes = %i' % (self.nTax, allOnes)

        for n in self.iterPostOrder():
            if n != self.root:
                # print 'doing node %s' % n.nodeNum

                if not n.leftChild:
                    # A long int, eg 1L, has no upper limit on its value
                    try:
                        n.br.rawSplitKey = 1L << self.taxNames.index(
                            n.name)  # "<<" is left-shift
                    except ValueError:
                        gm.append('node.name %s' % n.name)
                        gm.append('is not in taxNames %s' % self.taxNames)
                        raise P4Error(gm)
                    # n.br.rawSplitKey = 1L << self.taxNames.index(n.name)  #
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
                            print '%7s     %4s       %4s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
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
                            print '%7s     %4s       %4s' % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
                        raise P4Error(gm)

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
        # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
        if 1 & aNode.br.rawSplitKey:
            # "^" is xor, a bit-flipper.
            aNode.br.splitKey = allOnes ^ aNode.br.rawSplitKey
        else:
            aNode.br.splitKey = aNode.br.rawSplitKey

    def checkSplitKeys(self, useOldName=False, glitch=True, verbose=True):
        gm = ['Tree.checkSplitKeys()']

        allOnes = 2L ** (self.nTax) - 1
        # print 'nTax = %i, allOnes = %i' % (self.nTax, allOnes)

        isBad = False
        for n in self.iterPostOrder():
            if n != self.root:
                # print 'doing node %s' % n.nodeNum

                if not n.leftChild:
                    # A long int, eg 1L, has no upper limit on its value
                    try:
                        if useOldName:
                            rawSplitKey = 1L << self.taxNames.index(
                                n.oldName)  # "<<" is left-shift
                        else:
                            rawSplitKey = 1L << self.taxNames.index(
                                n.name)  # "<<" is left-shift
                    except ValueError:
                        if useOldName:
                            gm.append('node.name %s' % n.oldName)
                        else:
                            gm.append('node.name %s' % n.name)
                        gm.append('is not in taxNames %s' % self.taxNames)
                        raise P4Error(gm)
                    # n.br.rawSplitKey = 1L << self.taxNames.index(n.name)  #
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
                    print "checkSplitKeys node %2i rawSplitKey: existing %s, calculated %s" % (n.nodeNum, n.br.rawSplitKey, rawSplitKey)
                    isBad = True
                if n.br.splitKey != splitKey:
                    print "checkSplitKeys node %2i splitKey: existing %s, calculated %s" % (n.nodeNum, n.br.splitKey, splitKey)
                    isBad = True
        if glitch and isBad:
            raise P4Error(gm)
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
        # theTS.dump()
        assert theTS.mask
        rawSplitKey = 0L
        for i in range(len(theTS.mask)):
            # print i, theTS.mask[i]
            if theTS.mask[i] == '1':
                rawSplitKey += (1L << i)
        # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
        if 1 & rawSplitKey:
            allOnes = 2L ** (self.nTax) - 1
            splitKey = allOnes ^ rawSplitKey  # "^" is xor, a bit-flipper.
        else:
            splitKey = rawSplitKey
        # print "got splitKey %s" % splitKey

        for n in self.nodes:
            if n.br and not n.isLeaf:
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
                otherTree.nodes[nNum].br.lenChanged = self.nodes[
                    nNum].br.lenChanged
                #otherTree.nodes[nNum].br.flag = self.nodes[nNum].flag
                otherTree.nodes[nNum].br.splitKey = self.nodes[
                    nNum].br.splitKey
                otherTree.nodes[nNum].br.rawSplitKey = self.nodes[
                    nNum].br.rawSplitKey

        # model usage numbers
        if self.model:
            for nNum in range(len(self.nodes)):
                selfNode = self.nodes[nNum]
                otherNode = otherTree.nodes[nNum]
                for pNum in range(self.model.nParts):
                    otherNode.parts[pNum].compNum = selfNode.parts[
                        pNum].compNum
                    if selfNode != self.root:
                        otherNode.br.parts[pNum].rMatrixNum = selfNode.br.parts[
                            pNum].rMatrixNum
                        otherNode.br.parts[
                            pNum].gdasrvNum = selfNode.br.parts[pNum].gdasrvNum

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

        complaintHead = '\nTree.verifyIdentityWith()'  # keep

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
                # if self.nodes[nNum].br.len != otherTree.nodes[nNum].br.len:
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
            # if otherTree.partLikes[pNum] != self.partLikes[pNum]:
            if math.fabs(otherTree.partLikes[pNum] - self.partLikes[pNum]) > 1.e-8:
                print complaintHead
                print "    partLikes differ.  (%.5f, (%g) %.5f (%g)" % (
                    otherTree.partLikes[pNum], otherTree.partLikes[pNum], self.partLikes[pNum], self.partLikes[pNum])
                return var.DIFFERENT

        if 0:  # some more
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
            print "*=" * 25
            self.draw()
            gm = ["Tree.nextNode() spoke=%i, hub=%i" %
                  (spoke.nodeNum, hub.nodeNum)]
            gm.append(
                "Need to have either spoke.parent == hub or hub == spoke.")
            raise P4Error(gm)

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
            gm.append("The 'metric' arg should be one of %s" %
                      var.topologyDistanceMetrics)
            raise P4Error(gm)
        if metric == 'scqdist':  # no need for taxNames
            try:
                import scqdist
            except ImportError:
                gm.append(
                    "Could not find the 'scqdist' module needed for this metric.")
                gm.append(
                    "See the instructions for making it in the p4 source, in the Qdist directory.")
                raise P4Error(gm)
            tsSelf = self.writeNewick(toString=True)
            tsTree2 = tree2.writeNewick(toString=True)
            return scqdist.qdist(tsSelf, tsTree2)
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
            theSD = len(
                self.splitKeySet.symmetric_difference(tree2.splitKeySet))
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
        from btv import TV
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
        from btv import BTV
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
        # for sk in self.splitKeyHash.iterkeys():
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

        from btv import TV
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

    def node(self, specifier):
        """Get a node based on a specifier.

        The *specifier* can be a nodeNum, name, or node object.
        """

        nodeNum = None

        gm = ['Tree.node()']

        if type(specifier) == types.IntType:
            nodeNum = specifier
        elif type(specifier) == numpy.int32:
            nodeNum = specifier
        elif isinstance(specifier, Node):   # if its a node object
            if specifier in self.nodes:
                return specifier
            else:
                gm.append(
                    "The specifier is a node object, but is not part of self.")
                raise P4Error(gm)
        elif type(specifier) == types.StringType:   # if its a string
            for n in self.iterNodes():
                if n.name == specifier:
                    return n
            # if we haven't found a node matching the specier...
            if nodeNum == None:
                gm.append(
                    "Specifier string '%s' is not a node name.  What gives?" % specifier)
                raise P4Error(gm)

        else:
            gm.append("I don't understand the specifier '%s', type '%s'." % (
                specifier, type(specifier)))
            raise P4Error(gm)

        if nodeNum < 0 or nodeNum >= len(self.nodes):
            gm.append("The node number is out of range.")
            raise P4Error(gm)

        return self.nodes[nodeNum]

    def rotateAround(self, specifier):
        """Rotate a clade around a node.

        The *specifier* can be a nodeNum, name, or node object.
        """

        gm = ['Tree.rotateAround()']
        rotateNode = self.node(specifier)
        if rotateNode.isLeaf:
            print gm[0]
            print "    The rotateNode is a terminal node.  Not doing anything ..."
            return
        if rotateNode.getNChildren() == 1:
            print gm[0]
            print "    The rotateNode has only one child.  Not doing anything..."
            return

        # set up to unattach the rightmost child, and reattach at the left
        oldLeftChild = rotateNode.leftChild
        oldRightmostChild = rotateNode.rightmostChild()
        # we need to know the node second from the right, which will become the
        # newRightmostChild
        newRightmostChild = rotateNode.leftChild  # as a first guess...
        while newRightmostChild.sibling != oldRightmostChild:
            newRightmostChild = newRightmostChild.sibling
        # now the newRightmostChild is indeed the second from the right
        # now unattach the rightmost child, and reattach at the left
        newRightmostChild.sibling = None
        oldRightmostChild.sibling = oldLeftChild
        oldLeftChild.parent.leftChild = oldRightmostChild
        self.preAndPostOrderAreValid = 0

    def reRoot(self, specifier, moveInternalName=True, stompRootName=True, checkBiRoot=True, fixRawSplitKeys=False):
        """Re-root the tree to the node described by the specifier.

        The *specifier* can be a node.nodeNum, node.name, or node object.

        Here is a potential problem.  Lets say you start with this tree,
        with split support as shown::

            (((A, B)99, C)70, D, E);

                               +--------3:A
                     +---------2:99
            +--------1:70      +--------4:B
            |        |
            |        +---------5:C
            0
            |--------6:D
            |
            +--------7:E


        Now we want to reRoot it to node 2.  So we do that. Often,
        node.name's are really there to label the branch, not the node;
        you may be using node.name's to label your branches (splits), eg
        with support values.  If that is the case, then you want to keep
        the node name with the branch (not the node) as you reRoot().  Ie
        you want node labels to behave like branch labels.  If that is the
        case, we set moveInternalName=True; that is the default.  When
        that is done in the example above, we get::


            (A, B, (C, (D, E)70)99);

            +--------3:A
            |
            |--------4:B
            2
            |        +---------5:C
            +--------1:99
                     |         +--------6:D
                     +---------0:70
                               +--------7:E


        Another possibility is that the node names actually are there to
        name the node, not the branch, and you want to keep the node name
        with the node during the re-rooting process.  That can be done by
        setting moveInternalName=False, and the tree below is what happens
        when you do that.  Each internal node.name stays with its node in
        the re-rooting process::


            (A, B, (C, (D, E))70)99;

               +--------3:A
               |
               |--------4:B
            99:2
               |        +---------5:C
               +--------1:70
                        |         +--------6:D
                        +---------0
                                  +--------7:E

        Now if you had the default moveInternalName=True and you had a
        node name on the root, that would not work (draw it out to
        convince yourself ...).  So in that case you probably want
        stompRootName=True as well --- Both are default.  If you set
        stompRootName=True, it gives you a little warning as it does it.
        If you set stompRootName=2, it will do it silently.  If
        moveInternalName is not set, then stompRootName is not used.

        You probably do not want to ``reRoot()`` if the current root is
        bifurcating.  If you do that, you will get a node in the tree with
        a single child, which is strictly ok, but useless and confusing.
        So by default I checkBiRoot=True, and throw a P4Error if there is
        one.  If you want to draw such a pathological tree with a node
        with a single child, set checkBiRoot=False, and it will allow it.

        """

        gm = ['Tree.reRoot()']

        # The user probably does not want to reRoot if the current root is
        # bifurcating.  So check that.
        if checkBiRoot:
            if self.root.getNChildren() == 2:
                gm.append(
                    "The tree has a bifurcating root, so you probably do not")
                gm.append(
                    "want to reRoot() it.  You can remove the bifurcating root")
                gm.append(
                    "with yourTree.removeRoot().  If you really want to reRoot()")
                gm.append(
                    "with a bifurcating root, set checkBiRoot=False in the reRoot() args.")
                raise P4Error(gm)
        if self.root.isLeaf and not self.root.name:
            gm.append("The root is a leaf, but has no name.")
            gm.append(
                "So when you reRoot() it, some other leaf will have no name.")
            gm.append("That is a recipe for trouble, and is not allowed.")
            raise P4Error(gm)

        newRoot = self.node(specifier)
        oldRoot = self.root
        if newRoot == oldRoot:
            return

        if moveInternalName:
            if self.root.name and not self.root.isLeaf:
                if stompRootName:
                    if stompRootName != 2:
                        print "Notice.  Tree.reRoot(stompRootName) is set, so the root name '%s' is being zapped..." % self.root.name
                        print "(Set stompRootName=2 to do this silently ...)"
                    self.root.name = None
                else:
                    gm.append(
                        "Setting 'moveInternalName' implies keeping node names with their branches.")
                    gm.append(
                        "The root in this tree has a name, but has no branch.")
                    gm.append("So that does not work.")
                    gm.append("Set arg stompRootName to work around this.")
                    raise P4Error(gm)
            if self.root.isLeaf and self.root.leftChild.name:
                assert self.root.name
                gm.append("The current root is a leaf, with a name.")
                gm.append(
                    "Its sole child has a name also, which is an 'internal' node name.")
                gm.append("Arg moveInternalName is turned on.")
                gm.append(
                    "So when the tree gets re-rooted, that internal node name  should stay with its branch, not its node.")
                gm.append(
                    "But the rerooted branch will already have a name-- the current root name.")
                gm.append("So that does not work.")
                raise P4Error(gm)

        path = [newRoot]
        theParent = newRoot.parent
        while theParent != oldRoot:
            path.append(theParent)
            theParent = theParent.parent
        if 0:
            print "path from the newRoot to the oldRoot:"
            for i in path:
                print "   %i" % i.nodeNum

        while len(path):
            # We reverse the path above.  Its last entry is a child of the old root.
            # Start there-- its to be called the 'swivelNode'
            # print "Rooting on node %i..." % path[-1].nodeNum

            # +----------2:A
            # +----------1
            # |          +----------3:B
            # 0
            # +----------4:C
            # |
            # +----------5:D

            # +----------2:A
            # |
            # 1----------3:B
            # |
            # |          +----------4:C
            # +----------0
            # +----------5:D

            # in case this is not the first time around this loop...
            oldRoot = self.root
            swivelNode = path.pop()

            # Move the old root around the swivel node to the right.
            oldRoot.parent = swivelNode
            oldRoot.br = swivelNode.br
            swivelNode.br = None

            # If we are keeping the node.name with the branch, then do so now.
            if moveInternalName:
                if swivelNode.isLeaf:
                    pass
                elif oldRoot.isLeaf:
                    pass
                else:
                    oldRoot.name = swivelNode.name
                    swivelNode.name = None

            swivelNode.parent = None  # as befits a root
            if swivelNode.leftChild:
                swivelNode.rightmostChild().sibling = oldRoot
            else:
                swivelNode.leftChild = oldRoot

            # What we do next depends on whether the swivelNode is the
            # left, middle, or right child of the oldRoot
            if oldRoot.leftChild == swivelNode:
                oldRoot.leftChild = swivelNode.sibling
            else:
                leftSib = oldRoot.leftChild
                while leftSib.sibling != swivelNode:
                    leftSib = leftSib.sibling
                # At this point leftSib is the node that was left sib of the swivelNode
                # But the swivelNode is no longer there, so skip to the other side of it
                # (which may be None, if the swivelNode was rightmost child of the oldRoot).
                leftSib.sibling = swivelNode.sibling
            swivelNode.sibling = None  # as befits a root
            self.root = swivelNode

            # The splitKey is still good, but the rawSplitKey needs updating.
            if fixRawSplitKeys and oldRoot.br.rawSplitKey:
                oldRoot.br.rawSplitKey = 0L
                if oldRoot.leftChild:
                    for n in oldRoot.iterChildren():
                        oldRoot.br.rawSplitKey += n.br.rawSplitKey
                else:
                    oldRoot.br.rawSplitKey = 1L << self.taxNames.index(
                        oldRoot.name)  # "<<" is left-shift

        self.preAndPostOrderAreValid = 0

    def removeRoot(self):
        """Removes the root if self.root is mono- or bifurcating.

        This removes the root node if the tree is rooted on a terminal
        node, or if the tree is rooted on a bifurcating node.  Otherwise,
        it refuses to do anything.

        In the usual case of removing a bifurcating root, the branch
        length of one fork of the bifurcation is added to the other fork,
        so the tree length is preserved.

        In the unusual case of removing a monofurcating root (a root that
        is a terminal node, a tree-on-a-stick) then its branch length
        disappears.

        """

        gm = ['\nTree.removeRoot()']

        oldRoot = self.root
        newRoot = None
        if not oldRoot.leftChild:
            gm.append("The root has no children.")
            raise P4Error(gm)

        self.deleteCStuff()

        nRootChildren = self.root.getNChildren()
        if nRootChildren == 1:
            # The root has only one child.  Its rooted on a terminal node.
            # Its like this:
            #      +---A
            #  0---1
            #      +---B
            newRoot = oldRoot.leftChild
            newRoot.parent = None
            newRoot.br = None
        elif nRootChildren == 2:
            # The root has exactly two children.  This would be the usual,
            # expected, case.  We want to find the first child (of the
            # oldRoot) with more than one child (of the child).
            newRoot = oldRoot.leftChild  # Try this one first:
            # Ideally, and usually, the new root should have two or more
            # children.  Imagine what would happen if this tree
            #
            #  +-----A
            #  0   +----B
            #  +---2
            #      +----C
            #
            # had its root removed, and then
            # had A for its new root.  We get
            #      +---B
            #  A---2
            #      +---C

            # It would be better in this case if the new root was 2.  So we
            # want to search for candidate newRoots that have at least 2
            # children.  So is the current candidate, newRoot, good
            # enough?
            while newRoot:
                if newRoot.getNChildren() >= 2:
                    break  # its good enough, use it
                newRoot = newRoot.sibling  # its not good, try its sib

            # If we got this far an failed to find a good root, then
            # newRoot is None.  We are dealing with a tree where all the
            # root children are either leaves or have only one child.  Ok,
            # its unusual, but we do need a new root, even if it is a
            # terminal node.  So arbitrarily choose the leftChild to root
            # on.
            if not newRoot:
                newRoot = oldRoot.leftChild
                # print "x Re-rooting on node number %i" % newRoot.nodeNum
                if newRoot.leftChild:
                    newCh = newRoot.leftChild
                    newCh.sibling = newRoot.sibling
                else:
                    newCh = newRoot.sibling
                newRoot.leftChild = newCh
                newRoot.parent = None
                newRoot.sibling = None
                newCh.parent = newRoot
                newCh.br.len = newCh.br.len + newRoot.br.len
                newRoot.br = None
            else:   # We found a good new root, with two or more children
                # print "y reRooting on nodeNum %i" % newRoot.nodeNum
                if newRoot == oldRoot.leftChild:
                    newRoot.sibling.br.len = newRoot.sibling.br.len + \
                        newRoot.br.len
                    newRoot.br = None
                    newRoot.rightmostChild().sibling = newRoot.sibling
                    newRoot.sibling.parent = newRoot
                    newRoot.parent = None
                    newRoot.sibling = None
                else:
                    # the new root is the right child of the old
                    oldLeftCh = newRoot.leftChild
                    oldRoot.leftChild.br.len += newRoot.br.len
                    newRoot.br = None
                    newRoot.leftChild = oldRoot.leftChild
                    newRoot.leftChild.parent = newRoot
                    newRoot.leftChild.sibling = oldLeftCh
                    newRoot.parent = None
                    newRoot.sibling = None

        else:
            gm.append("The root has more than two children.")
            gm.append(
                "Removing the root with more than two children is not implemented.")
            gm.append("Are you even sure you want to do that?")
            raise P4Error(gm)

        oldRoot.wipe()
        self.nodes.remove(oldRoot)
        del oldRoot
        if newRoot:
            self.root = newRoot
        else:
            gm.append("No newRoot?   Programming error.")
            raise P4Error(gm)
        for i in range(len(self.nodes)):
            self.nodes[i].nodeNum = i
        self.preOrder = None
        self.postOrder = None
        self.preAndPostOrderAreValid = 0
        self.setPreAndPostOrder()
        self._nTax = 0

    def removeNode(self, specifier, alsoRemoveSingleChildParentNode=True, alsoRemoveBiRoot=True, alsoRemoveSingleChildRoot=True):
        """Remove a node, together with everything above it.

        Arg *specifier* can be a nodeNum, name, or node object.

        So lets say that we have a tree like this::

            +-------1:A
            0
            |       +--------3:B
            +-------2
                    +--------4:C

        and we remove node 4.  When it is removed, node 2 ends up having
        only one child.  Generally you would want to remove it as well (so
        that the parent of node 3 is node 0), so the option
        *alsoRemoveSingleChildParentNode* is turned on by default.  If
        *alsoRemoveSingleChildParentNode* is turned off, nodes like node 2
        get left in the tree.

        Removal of a node might cause the creation of a bifurcating root.
        I assume that is not desired, so alsoRemoveBiRoot is turned on by
        default.

        In the example above, if I were to remove node 4, by default node
        2 would also disappear, but by default node 0 would also disappear
        because it would then be a tree with a bifurcating root node.  So
        starting with a 5-node tree, by removing 1 node you would end up
        with a 2-node tree, with 1 branch.

        In the example here, if I were to simply remove node 1::

            +--------1:A
            |
            0        +---------3:B
            +--------2
                     |         +--------5:C
                     +---------4
                               +--------6:D

        Then node 0 would remain, as::

                     +---------2:B
            0--------1
                     |         +--------4:C
                     +---------3
                               +--------5:D

        Presumably that is not wanted, so the arg
        *alsoRemoveSingleChildRoot* is turned on by default.  When the root
        is removed, we are left with a bi-root, which (if the arg
        alsoRemoveBiRoot is set) would also be removed.  The resulting
        tree would be::

            +-------0:B
            |
            1-------2:C
            |
            +-------3:D

        The deleted nodes are really deleted, and do not remain in self.nodes.
        """

        gm = ['Tree.removeNode()']

        rNode = self.node(specifier)
        if rNode == self.root:
            print gm[0]
            print "    The specified node appears to be the root."
            print "    Removing everything above the root would leave nothing."
            print "    I assume that you do not want to do that."
            print "    So I'm not doing that."
            # self.draw()
            # print "the specifier was %s" % specifier
            # print "the specified node was node number %i" % rNode.nodeNum
            #raise P4Error
            return
        rNodeParnt = rNode.parent

        self.deleteCStuff()

        # For cases where the tree is originally with a single child root
        # -- we don't want to then delete that root below.
        assert self.root.leftChild
        isOriginallySingleChildRoot = False
        if not self.root.leftChild.sibling:
            isOriginallySingleChildRoot = True

        # For cases where we originally have a bi-Root, that we would like to
        # keep.
        isOriginallyBiRoot = False
        if self.root.getNChildren() == 2:
            isOriginallyBiRoot = True

        # print "Removing node number", rNode.nodeNum
        #hitList = []
        #self.recursivelyListNodeIndicesDownTo(hitList, rNode)

        # getNodeNumsAbove does not require setting
        # self.preAndPostOrderAreValid = 0.  Its irrelevant.
        # does not include rNode
        hitList = self.getNodeNumsAbove(rNode.nodeNum)
        hitList.append(rNode.nodeNum)
        # hitList.sort()
        # hitList.reverse()
        # print "Hit list is", hitList
        hitNodes = []
        for i in hitList:
            hitNodes.append(self.nodes[i])

        # Disconnect it from the tree.
        # the rNode is the left, middle, or right child of the parent
        if rNodeParnt.leftChild == rNode:
            # its the left child
            rNodeParnt.leftChild = rNode.sibling
            rNode.sibling = None
        else:
            # its a middle or a right child
            # If its a right child, then rNode.sibling = None
            leftSib = rNode.leftSibling()
            if 1:
                if leftSib:
                    # print "leftSib is node %i" % leftSib.nodeNum
                    leftSib.sibling = rNode.sibling
                else:
                    gm.append("leftSib is None.  This shouldn't happen")
                    gm.append("Programming error?")
                    raise P4Error(gm)

        haveRemovedSingleChildRoot = False
        if not isOriginallySingleChildRoot:
            # The tree may have been bifurcating, and left a single-child root.
            if alsoRemoveSingleChildRoot and self.root.leftChild and not self.root.leftChild.sibling:
                hitNodes.append(self.root)
                self.root = self.root.leftChild
                self.root.parent = None
                haveRemovedSingleChildRoot = True

        for n in hitNodes:
            n.wipe()
            self.nodes.remove(n)
            if n.isLeaf and self.taxNames and n.name and n.name in self.taxNames:
                self.taxNames.remove(n.name)
            del n

        self._nTax = 0
        if 1:
            for i in range(len(self.nodes)):
                self.nodes[i].nodeNum = i
            # self.dump(node=1)
            self.preOrder = None
            self.postOrder = None
            self.preAndPostOrderAreValid = 0
            # self.draw()  # This won't work unless preAndPostOrderAreValid set
            # to 0

        # the parent of the removed node may now only have one child,
        # in which case it (the parent) should be removed

        if alsoRemoveSingleChildParentNode or alsoRemoveBiRoot or haveRemovedSingleChildRoot:

            ignoreBrLens = True
            self.preOrder = numpy.array(
                [var.NO_ORDER] * len(self.nodes), numpy.int32)
            self.postOrder = numpy.array(
                [var.NO_ORDER] * len(self.nodes), numpy.int32)

            if len(self.nodes) > 1:
                self.setPreAndPostOrder()
            for n in self.iterNodesNoRoot():
                if n.br.len != 0.1:
                    ignoreBrLens = False
                    break

        if alsoRemoveSingleChildParentNode:

            if rNodeParnt and rNodeParnt.leftChild and not rNodeParnt.leftChild.sibling:
                if rNodeParnt.parent:
                    # rNodeParnt is a left, middle, or right child
                    if rNodeParnt.parent.leftChild == rNodeParnt:
                        # print 'rNodeParnt is a left child'
                        rNodeParnt.parent.leftChild = rNodeParnt.leftChild
                    else:
                        # print 'rNodeParnt is a middle or right child'
                        rNodeParnt.leftSibling().sibling = rNodeParnt.leftChild
                    rNodeParnt.leftChild.sibling = rNodeParnt.sibling
                    if not ignoreBrLens:
                        rNodeParnt.leftChild.br.len += rNodeParnt.br.len
                    rNodeParnt.leftChild.parent = rNodeParnt.parent
                    rNodeParnt.wipe()
                    self.nodes.remove(rNodeParnt)
                    del rNodeParnt

        if not isOriginallyBiRoot:
            if self.root.getNChildren() == 2 and alsoRemoveBiRoot:
                self.removeRoot()
                if ignoreBrLens:
                    for ch in self.root.iterChildren():
                        ch.br.len = 0.1

        for i in range(len(self.nodes)):
            self.nodes[i].nodeNum = i
        self.preOrder = numpy.array(
            [var.NO_ORDER] * len(self.nodes), numpy.int32)
        self.postOrder = numpy.array(
            [var.NO_ORDER] * len(self.nodes), numpy.int32)

        if len(self.nodes) > 1:
            self.setPreAndPostOrder()

    def removeAboveNode(self, specifier, newName):
        """Remove everything above an internal node, making it a leaf, and so needing a new name.
        """
        rNode = self.node(specifier)
        assert rNode != self.root
        assert not rNode.isLeaf
        toDelete = [n for n in rNode.iterPostOrder() if n != rNode]
        # print [n.nodeNum for n in toDelete]
        for n in toDelete:
            # print 'deleting node %i' % n.nodeNum
            self.removeNode(n, alsoRemoveSingleChildParentNode=False)
        rNode.name = newName
        rNode.isLeaf = 1

    def collapseNode(self, specifier):
        """Collapse the specified node to make a polytomy, and remove it from the tree.

        Arg *specifier*, as usual, can be a node, node number, or node name.

        The specified node remains in self.nodes, and is returned.
        """

        theNode = self.node(specifier)
        assert theNode in self.nodes, "The specified Node is not in the tree."
        #assert theNode.leftChild and theNode.leftChild.sibling, "The specified node must have at least 2 children."
        assert not theNode.isLeaf, "The specified node must not be a leaf."
        assert theNode is not self.root, "The specified node must not be the tree root."
        # print "Collapsing node %i" % theNode.nodeNum

        theNewParent = theNode.parent
        theRightmostChild = theNode.rightmostChild()
        theLeftSib = theNode.leftSibling()
        if theLeftSib:
            theLeftSib.sibling = theNode.leftChild
        else:
            theNewParent.leftChild = theNode.leftChild
        for n in theNode.iterChildren():
            n.parent = theNewParent
        theRightmostChild.sibling = theNode.sibling
        theNode.wipe()

        # The following does not work well.
        # self.nodes.remove(theNode)
        # del(theNode)

        self.setPreAndPostOrder()
        self._nInternalNodes -= 1

    def pruneSubTreeWithoutParent(self, specifier, allowSingleChildNode=False):
        """Remove and return a node, together with everything above it.

        Arg *specifier* can be a nodeNum, name, or node object.

        By default, the arg allowSingleChildNode is turned off, and is for
        those cases where the parent of the node has more than 2 children.
        So when the subTree is removed, the parent node that is left
        behind has more than one child.

        The stuff that is removed is returned.  The nodes are left in
        self; the idea being that the subTree will be added back to the
        tree again (via reconnectSubTreeWithoutParent()).
        """

        gm = ['Tree.pruneSubTreeWithoutParent()']

        rNode = self.node(specifier)
        if rNode == self.root:
            gm.append("The specified node is the root.")
            raise P4Error(gm)
        rNodeParnt = rNode.parent
        if not allowSingleChildNode:
            if rNodeParnt.getNChildren() < 3:
                # self.draw()
                gm.append("The arg allowSingleChildNode is turned off.")
                gm.append(
                    "This would be for those cases where the parent of the subTree has more than 2 children.")
                raise P4Error(gm)

        # self.deleteCStuff()

        # print "rNode is node %i" % rNode.nodeNum
        # print "rNodeParnt is node %i" % rNodeParnt.nodeNum

        # Disconnect it from the tree.
        # the rNode is the left, middle, or right child of the parent
        if rNodeParnt.leftChild == rNode:
            # print "its the left child"
            rNodeParnt.leftChild = rNode.sibling
            rNode.sibling = None
            rNode.parent = None
        else:
            # print "its a middle or a right child"
            leftSib = rNode.leftSibling()
            assert leftSib
            leftSib.sibling = rNode.sibling
            rNode.sibling = None
            rNode.parent = None

        self.preAndPostOrderAreValid = 0
        self._nTax = 0
        return rNode

    def reconnectSubTreeWithoutParent(self, stNode, newParent, beforeNode=None):
        """Attach subtree stNode to the rest of the tree at newParent.


        The beforeNode is by default None, and then the subtree is
        reconnected as the rightmost child of the new parent.  However, if
        you want it somewhere else, for example as the leftmost child, or
        between two existing child nodes, specify a beforeNode (specified
        as usual as a node, nodeNumber, or node name) and the subtree will
        be inserted there.
        """

        gm = ["Tree.reconnectSubTreeWithoutParent()"]

        newParent = self.node(newParent)
        if not newParent.leftChild:
            gm.append("Can't attach to a leaf.")
            raise P4Error(gm)
        stNode.parent = newParent
        if beforeNode == None:  # easy -- just add it to the rightmost child.
            theRMChildOfNewParent = newParent.rightmostChild()
            theRMChildOfNewParent.sibling = stNode
        else:
            bNode = self.node(beforeNode)
            if bNode.parent != newParent:
                gm.append(
                    "The parent of the 'beforeNode' should be the newParent.")
                raise P4Error(gm)
            if newParent.leftChild == bNode:
                newParent.leftChild = stNode
            else:
                lSib = newParent.leftChild
                while lSib.sibling != bNode:
                    lSib = lSib.sibling
                lSib.sibling = stNode
            stNode.sibling = bNode

        self._nTax = 0

    def addNodeBetweenNodes(self, specifier1, specifier2):
        """Add a node between 2 exisiting nodes, which should be parent-child.

        The *specifier* can be a nodeNum, name, or node object.

        Returns the new node object.
        """

        gm = ['Tree.addNodeBetweenNodes()']

        aNode1 = self.node(specifier1)
        aNode2 = self.node(specifier2)

        # aNode1 should be the parent of aNode2
        if aNode1 == aNode2.parent:
            pass
        elif aNode1 == aNode2:
            gm.append("The two specified nodes are the same.")
            raise P4Error(gm)
        elif aNode2 == aNode1.parent:
            temp = aNode1
            aNode1 = aNode2
            aNode2 = temp
        else:
            gm.append(
                "The 2 specified nodes should have a parent-child relationship")
            raise P4Error(gm)

        self.deleteCStuff()

        hasBrLens = False
        for n in self.iterNodes():
            if n.br and math.fabs(n.br.len - 0.1) > 1.e-15:
                hasBrLens = True
                break

        newNode = copy.deepcopy(aNode2)
        newNode.br.len = 0.1
        newNode.name = None
        newNode.isLeaf = False
        newNode.nodeNum = len(self.nodes)
        self.nodes.append(newNode)

        newNode.parent = aNode1
        newNode.leftChild = aNode2
        aNode2.parent = newNode
        if aNode1.leftChild == aNode2:
            aNode1.leftChild = newNode
        else:
            oldCh = aNode1.leftChild
            while oldCh.sibling != aNode2:
                oldCh = oldCh.sibling
            oldCh.sibling = newNode
        if aNode2.sibling:
            newNode.sibling = aNode2.sibling
            aNode2.sibling = None
        if hasBrLens:
            halfBrLen = aNode2.br.len / 2.0
            aNode2.br.len = halfBrLen
            newNode.br.len = halfBrLen

        if 1:
            self.preOrder = numpy.array(
                [var.NO_ORDER] * len(self.nodes), numpy.int32)
            self.postOrder = numpy.array(
                [var.NO_ORDER] * len(self.nodes), numpy.int32)
            if len(self.nodes) > 1:
                self.setPreAndPostOrder()
        return newNode

    def allBiRootedTrees(self):
        """Returns a Trees object containing all possible bi-rootings of self.

        Self should have a root node of degree > 2, but need not be fully
        resolved.

        Self needs a taxNames.
        """

        gm = ['Tree.allBiRootedTrees()']

        if self.root.getNChildren() < 3:
            gm.append("Self root should be of degree > 2.")
            raise P4Error(gm)
        if not self.taxNames:
            gm.append("Self (ie the tree) needs to have taxNames set.")
            raise P4Error(gm)

        tList = []
        for i in range(len(self.nodes)):
            if self.nodes[i] == self.root:
                pass
            else:
                t = self.dupe()
                n = t.nodes[i]
                x = t.addNodeBetweenNodes(n.parent, n)
                t.reRoot(x, moveInternalName=False)
                t.name = 'r%i' % n.nodeNum
                tList.append(t)

        from trees import Trees
        tt = Trees(trees=tList, taxNames=self.taxNames)
        return tt

    def ladderize(self, biggerGroupsOnBottom=True):
        """Rotate nodes for a staircase effect.

        This method, in its default biggerGroupsOnBottom way, will take a
        tree like this::

                                        +---------4:A
                               +--------3
                     +---------2        +---------5:B
                     |         |
                     |         +--------6:C
            +--------1
            |        |         +--------8:D
            |        +---------7
            |                  +--------9:E
            0
            |--------10:F
            |
            |        +---------12:G
            +--------11
                     +---------13:H


        and rearranges it so that it is like ... ::


            +--------10:F
            |
            |        +---------12:G
            |--------11
            |        +---------13:H
            0
            |                  +--------8:D
            |        +---------7
            |        |         +--------9:E
            +--------1
                     |         +--------6:C
                     +---------2
                               |        +---------4:A
                               +--------3
                                        +---------5:B

        Note that for each node, the more populated group is on the bottom,
        the secondmost populated second, and so on.

        To get it with the bigger groups on top, set
        biggerGroupsOnBottom=False.  I made the default with the bigger
        groups on the bottom so that it often makes room for a scale bar.

        The setting biggerGroupsOnBottom, the default here, would
        equivalent to set torder=right in paup; torder=left puts the
        bigger groups on the top.

        """

        # self.draw()
        self.root._ladderize(biggerGroupsOnBottom)
        self.preAndPostOrderAreValid = 0
        # self.draw()

    def randomizeTopology(self, randomBrLens=True):

        gm = ["Tree.randomizeTopology()"]
        if self.root.getNChildren() != 3 or not self.isFullyBifurcating():
            gm.append(
                "Should be a fully bifurcating tree, this week.  Fix me?")
            raise P4Error(gm)
        if self.cTree:
            self.deleteCStuff()
        nTax = self.nTax
        oldRoot = self.root
        leaves = []
        internals = []
        for n in self.nodes:
            if n.isLeaf:
                leaves.append(n)
            else:
                internals.append(n)
        random.shuffle(leaves)
        random.shuffle(internals)
        lIndx = 0
        iIndx = 0
        self.nodes = []
        nodeNum = 0

        # new root
        self.root = internals[iIndx]
        iIndx += 1
        self.nodes.append(self.root)
        self.root.parent = None
        self.root.leftChild = None
        self.root.sibling = None
        if self.root == oldRoot:
            pass
        else:
            oldRoot.br = self.root.br
            self.root.br = None
        self.root.nodeNum = nodeNum
        nodeNum += 1

        # add a leaf as left child to the root
        n = leaves[lIndx]
        lIndx += 1
        self.nodes.append(n)
        n.sibling = None
        n.nodeNum = nodeNum
        nodeNum += 1
        self.root.leftChild = n
        n.parent = self.root
        previousNode = n

        # add the rest of the leaves
        while lIndx < nTax:
            n = leaves[lIndx]
            lIndx += 1
            self.nodes.append(n)
            n.sibling = None
            n.nodeNum = nodeNum
            nodeNum += 1
            previousNode.sibling = n
            previousNode = n
            n.parent = self.root

        # Now we have a star tree.  Now add internal nodes until it is all
        # resolved, which needs nTax - 3 nodes
        for i in range(nTax - 3):
            ssNodes = []
            for n in self.nodes:
                if n.sibling and n.sibling.sibling:
                    ssNodes.append(n)
            lChild = random.choice(ssNodes)

            # print "lChild = node %i" % lChild.nodeNum

            # +----------1:oldLeftSib
            # |
            # +----------2:lChild
            # 0
            # +----------3:lChildSib
            # |
            # +----------4:oldLChildSibSib

            # +----------1:oldLeftSib
            # |
            # |          +----------3:lChild
            # 0----------2(n)
            # |          +----------4:lChildSib
            # |
            # +----------5:oldLChildSibSib

            n = internals[iIndx]
            iIndx += 1
            n.parent = None
            n.sibling = None
            n.leftChild = None
            n.nodeNum = nodeNum
            nodeNum += 1
            lChildSib = lChild.sibling  # guarranteed to have one
            oldLChildSibSib = lChildSib.sibling  # ditto
            oldLeftSib = lChild.parent.leftChild  # first guess ...
            if oldLeftSib != lChild:
                while oldLeftSib.sibling != lChild:
                    oldLeftSib = oldLeftSib.sibling
            else:
                oldLeftSib = None
            if 0:
                if oldLeftSib:
                    print "oldLeftSib = %i" % oldLeftSib.nodeNum
                else:
                    print "oldLeftSib = None"
                print "lChildSib = %i" % lChildSib.nodeNum
                if oldLChildSibSib:
                    print "oldLChildSibSib = %i" % oldLChildSibSib.nodeNum
                else:
                    print "oldLChildSibSib = None"

            if oldLeftSib:
                oldLeftSib.sibling = n
            else:
                lChild.parent.leftChild = n

            n.parent = lChild.parent
            lChild.parent = n
            n.leftChild = lChild
            lChildSib.parent = n
            n.sibling = oldLChildSibSib
            lChildSib.sibling = None
            self.nodes.append(n)

        # self.dump(all=True)
        # self.draw()

        # The way it is now, the root rightmost child is always a
        # leaf.  Not really random, then, right?  So choose a random
        # internal node, and re-root it there.
        # print "nTax=%i, len(t.nodes)=%i" % (nTax, len(t.nodes))
        if nTax > 3:
            n = self.nodes[random.randrange(nTax + 1, len(self.nodes))]
            self.reRoot(n, moveInternalName=False)

        # The default is to have randomBrLens, where internal nodes get
        # brLens of 0.02 - 0.05, and terminal nodes get brLens of 0.2 -
        # 0.5.  Branch lengths are all 0.1 if randomBrLens is turned
        # off.
        if randomBrLens:
            for n in self.nodes:
                if n != self.root:
                    if n.isLeaf:
                        n.br.len = 0.02 + (random.random() * 0.48)
                    else:
                        n.br.len = 0.02 + (random.random() * 0.03)

        self.preAndPostOrderAreValid = 0
        # self.dump(all=True)

    def readBipartitionsFromPaupLogFile(self, thePaupLogFileName):
        """Assigns support to the tree, from the PAUP bipartitions table.

        This needs to have self.taxNames set.

        This is useful if you want to make a consensus tree using PAUP,
        and get the support values.  When you make a cons tree with PAUP,
        the support values, usually bootstrap values, are unfortunately
        not saved with the tree.  That information is in the Bipartitions
        table, which can be saved to a PAUP log file so that p4 can get
        it.  This method will read thePaupLogFileName and extract the
        split (tree bipartition) supports, and assign those supports to self
        as node.br.support's (as a float, not a string).

        It also returns a hash with the split strings as keys and the
        split support as values, if you need it.
        """

        gm = ['Tree.readBipartitionsFromPaupLogFile()']

        if not self.taxNames:
            gm.append("This method needs self.taxNames.")
            raise P4Error(gm)
        self.checkTaxNames()

        f = open(thePaupLogFileName, 'r')
        while 1:
            aLine = f.readline()
            if not aLine:
                gm.append("No Bipartitions line in %s?" % thePaupLogFileName)
                raise P4Error(gm)
            if aLine.startswith('Bipartitions found'):
                # print aLine
                break

        # The splits table might be all in one piece, or it might be split
        # into sections (CYCLEs).  Only the last section has the Freq or
        # percent.  If the first section is the only section, then it is
        # the LAST_CYCLE, so that we get the Freq or percent.

        # Surprise-- the PAUP output contains both Freq and % support,
        # unless the sum of the weights (or if there are no weights, the
        # number of trees) is 100, in which case only the Freq is given.

        FIRST_CYCLE = -1
        MIDDLE_CYCLE = 0
        LAST_CYCLE = 1
        cycleType = None
        thisKeyLength = None
        accumulatedKeyLength = 0
        nKeys = 0
        keyCounter = 0
        theHash = {}
        theKeys = []
        hasFreqOnly = None

        while 1:
            prevLine = aLine
            aLine = f.readline()
            # print "a Got Line ==>%s<==" % aLine
            if not aLine:
                gm.append('Unexpected end of file.')
                raise P4Error(gm)
            aLine = aLine.rstrip()
            # print "b Got Line ==>%s<==" % aLine
            if len(aLine):
                if aLine[0] == '-':
                    # print prevLine
                    # print aLine
                    if prevLine.endswith('Freq') or prevLine.endswith('%'):
                        if prevLine.endswith('Freq'):
                            hasFreqOnly = 1
                        elif prevLine.endswith('%'):
                            hasFreqOnly = 0

                        cycleType = LAST_CYCLE
                        # print "We are now in the last cycle.  hasFreqOnly=%s"
                        # % hasFreqOnly
                        splitLine = string.split(prevLine)
                        thisKeyLength = len(splitLine[0])
                        keyCounter = 0
                    elif cycleType == MIDDLE_CYCLE:
                        # print "We are now in a middle cycle"
                        thisKeyLength = len(prevLine)
                        keyCounter = 0
                    elif cycleType == FIRST_CYCLE:
                        gm.append('This should never happen.')
                        raise P4Error(gm)
                    else:
                        cycleType = FIRST_CYCLE
                        # print "We are now in the first cycle"
                        thisKeyLength = len(prevLine)

                elif aLine[0] in ['.', '*']:
                    if cycleType in [FIRST_CYCLE, MIDDLE_CYCLE] and len(aLine) != thisKeyLength:
                        gm.append("Unequal key lengths?!?")
                        gm.append("thisKeyLength = %i" % thisKeyLength)
                        gm.append("%s" % aLine)
                        raise P4Error(gm)

                    if cycleType == FIRST_CYCLE:
                        theKeys.append(aLine)
                    elif cycleType == MIDDLE_CYCLE:
                        theKeys[keyCounter] = theKeys[keyCounter] + aLine
                        keyCounter = keyCounter + 1
                        if keyCounter > nKeys:
                            gm.append("Too many keys.")
                            raise P4Error(gm)
                    elif cycleType == LAST_CYCLE:
                        # It will usually be a line like: ...*....*.     67.83  67.9%
                        # But it will sometimes be a line like: ...*....*.
                        # 68
                        splitLine = string.split(aLine)
                        theKey = splitLine[0]
                        if hasFreqOnly:
                            theSupport = float(splitLine[-1])
                            theSupport /= 100.0
                        else:
                            # don't read the % sign
                            theSupport = float(splitLine[-1][:-1])
                            theSupport /= 100.0

                        if len(theKey) != thisKeyLength:
                            gm.append("Unequal key lengths?!?")
                            gm.append("thisKeyLength = %i" % thisKeyLength)
                            gm.append("%s" % theKey)
                            raise P4Error(gm)
                        if accumulatedKeyLength:
                            theKeys[keyCounter] = theKeys[keyCounter] + theKey
                            theHash[theKeys[keyCounter]] = theSupport
                            keyCounter = keyCounter + 1
                            if keyCounter > nKeys:
                                gm.append("Too many keys.")
                                raise P4Error(gm)
                        # LAST_CYCLE is also the FIRST_CYCLE, ie the table is
                        # in one section.
                        else:
                            theKeys.append(theKey)
                            theHash[theKey] = theSupport
                    else:
                        gm.append("This should never happen.")
                        raise P4Error(gm)
                else:
                    pass  # Skip lines with numbers
            else:  # a blank line
                if cycleType == FIRST_CYCLE:
                    cycleType = MIDDLE_CYCLE
                    nKeys = len(theKeys)
                    # print "At the end of the first cycle, got %i keys" %
                    # nKeys
                    accumulatedKeyLength = thisKeyLength
                elif cycleType == MIDDLE_CYCLE:
                    accumulatedKeyLength += thisKeyLength
                elif cycleType == LAST_CYCLE:
                    accumulatedKeyLength += thisKeyLength
                    nKeys = len(theKeys)
                    # print "finished"
                    # print theKeys[0]
                    break
                else:
                    pass
        f.close()

        if 0:
            print "Finished getting split strings, and supports"
            print "Got %i items in the hash" % len(theHash)
            print "accumulatedKeyLength = %i" % accumulatedKeyLength
            print "nKeys = %i" % nKeys
            # self.draw()

        # I will need to make splitStrings (in dot-star notation) from
        # splitKeys.  To do that, I can use the
        # func.getSplitStringFromKey(theKey, nTax) function.

        self.makeSplitKeys()
        for n in self.nodes:
            if n != self.root:
                if not n.isLeaf:
                    theNodeSplitString = func.getSplitStringFromKey(
                        n.br.splitKey, self.nTax)
                    if theHash.has_key(theNodeSplitString):
                        if hasattr(n.br, 'support') and n.br.support is not None:
                            gm.append(
                                "Node %i already has a br.support." % n.nodeNum)
                            gm.append(
                                "I am refusing to clobber it with the split support.")
                            gm.append(
                                "Either fix the tree or fix this method.")
                            raise P4Error(gm)
                        n.br.support = float(theHash[theNodeSplitString])
        return theHash

    def renameForPhylip(self, dictFName='p4_renameForPhylip_dict.py'):
        """Rename with phylip-friendly short boring names.

        It saves the old names (together with the new) in a python
        dictionary, in a file by default named p4_renameForPhylip_dict.py

        If self does not have taxNames set, it does not write
        originalNames to that file-- which may cause problems
        restoring names.  If you want to avoid that, be sure to set
        self.taxNames before you do this method.

        This method does not deal with internal node names, at all.  They
        are silently ignored.  If they are too long for phylip, they are
        still silently ignored, which might cause problems.
        """

        gm = ['Tree.renameForPhylip()']
        if os.path.exists(dictFName):
            gm.append("The dictionary file '%s' already exists." % dictFName)
            raise P4Error(gm)
        d = {}
        if self.taxNames:
            d2 = {}
            originalNames = self.taxNames[:]
            for i in range(len(self.taxNames)):
                oldName = self._taxNames[i]
                newName = 's%i' % i
                d[newName] = oldName
                d2[oldName] = newName
                self._taxNames[i] = newName
            for n in self.iterLeavesNoRoot():
                n.name = d2[n.name]
            if self.root.isLeaf and self.root.name:
                self.root.name = d2[self.root.name]

        else:
            originalNames = None
            i = 0
            for n in self.iterLeavesNoRoot():
                oldName = n.name
                newName = 's%i' % i
                d[newName] = oldName
                n.name = newName
                i += 1
            if self.root.isLeaf and self.root.name:
                oldName = self.root.name
                newName = 's%i' % i
                d[newName] = oldName
                self.root.name = newName

        f = file(dictFName, 'w')
        f.write("p4_renameForPhylip_originalNames = %s\np4_renameForPhylip_dict = %s\n" % (
            originalNames, d))
        f.close()

    def restoreNamesFromRenameForPhylip(self, dictFName='p4_renameForPhylip_dict.py'):
        """Given the dictionary file, restore proper names.

        The renaming is done by the Alignment method renameForPhylip(),
        which makes the dictionary file.  The dictionary file is by
        default named p4_renameForPhylip_dict.py
        """

        gm = ["Tree.restoreNamesFromRenameForPhylip()"]
        if os.path.exists(dictFName):
            import __main__
            execfile(dictFName, __main__.__dict__,  __main__.__dict__)
            from __main__ import p4_renameForPhylip_dict, p4_renameForPhylip_originalNames
        else:
            gm.append("The dictionary file '%s' can't be found." % dictFName)
            raise P4Error(gm)
        for n in self.iterNodes():
            if n.isLeaf:
                if p4_renameForPhylip_dict.has_key(n.name):
                    n.name = p4_renameForPhylip_dict[n.name]
                else:
                    gm.append(
                        "The dictionary does not contain a key for '%s'." % n.name)
                    raise P4Error(gm)
        if p4_renameForPhylip_originalNames:
            self.taxNames = p4_renameForPhylip_originalNames
        else:
            if self.taxNames:
                gm.append("self.taxNames is set, and should be replaced, but")
                gm.append("p4_renameForPhylip_originalNames is None. ?!?")
                raise P4Error(gm)
        del(__main__.p4_renameForPhylip_dict)
        del(__main__.p4_renameForPhylip_originalNames)

    def restoreDupeTaxa(self, dictFileName='p4DupeSeqRenameDict.py', asMultiNames=True):
        """Restore previously removed duplicate taxa from a dict file.

        The usual story would be like this: You read in your alignment and
        p4 tells you that you have duplicate sequences.  So you use the
        Alignment method checkForDuplicateSequences() to remove them,
        which makes a dictionary file, by default
        'p4DupeSeqRenameDict.py', to facilitate restoration of the names.
        You do your analysis on the reduced alignment, and get a tree.
        Then you use the dictionary file with this method to restore all
        the taxon names.

        If asMultiNames is turned on, the default, then the leaf nodes are
        not replicated, and the name is changed to be a long multi-name.

        If asMultiNames is turned off, then the restored taxa are made to
        be siblings, and the branch lengths are set to zero.

        """

        gm = ['Tree.restoreDupeTaxa()']
        if not os.path.isfile(dictFileName):
            gm.append("Can't find dict file '%s'" % dictFileName)
            raise P4Error(gm)
        loc = {}
        execfile(dictFileName, {}, loc)
        try:
            p4DupeSeqRenameDict = loc['p4DupeSeqRenameDict']
        except KeyError:
            gm.append(
                "Can't get the dictionary named 'p4DupeSeqRenameDict' from the dict file.")
        # print p4DupeSeqRenameDict

        kk = p4DupeSeqRenameDict.keys()
        # print kk

        if asMultiNames:
            for oldName in kk:
                rNode = self.node(oldName)
                newNames = p4DupeSeqRenameDict[oldName]
                multiName = ', '.join(newNames)
                rNode.name = multiName
        else:
            for oldName in kk:
                rNode = self.node(oldName)
                newNames = p4DupeSeqRenameDict[oldName]
                lCh = Node()
                lCh.isLeaf = 1
                lCh.name = newNames[0]
                lCh.br.len = 0.0
                lCh.parent = rNode
                lCh.nodeNum = len(self.nodes)
                self.nodes.append(lCh)
                rNode.isLeaf = 0
                rNode.name = None
                rNode.leftChild = lCh
                theLeftSib = rNode.leftChild
                for sibName in newNames[1:]:
                    n = Node()
                    n.parent = rNode
                    theLeftSib.sibling = n
                    n.isLeaf = 1
                    n.name = sibName
                    n.br.len = 0.0
                    n.nodeNum = len(self.nodes)
                    self.nodes.append(n)
                    theLeftSib = n
        self.taxNames = []
        self.preOrder = None
        self.postOrder = None
        self.setPreAndPostOrder()
        self._nTax = 0

    def lineUpLeaves(self, rootToLeaf=1.0, overWriteBrLens=True):
        """Make the leaves line up, as in a cladogram.

        This makes the rootToLeaf distance the same for all leaves.

        If overWriteBrLens is set, then the newly calculated br.lens
        replace the original br.lens.  If it is not set, then the new
        br.lens are placed in br.lenL, and does not over-write the
        original br.lens.  """

        levels = 0
        for n1 in self.iterLeavesNoRoot():
            nodeCount = 1
            p = n1.parent
            while p != self.root:
                p = p.parent
                nodeCount += 1
            if nodeCount > levels:
                levels = nodeCount
        # print "levels = %i" % levels
        for n1 in self.iterLeavesNoRoot():
            n1.lul_pos = levels
            ch = n1
            p = n1.parent
            while p:
                newPos = ch.lul_pos - 1
                if not hasattr(p, 'lul_pos'):
                    p.lul_pos = newPos
                else:
                    if p.lul_pos <= newPos:
                        pass
                    else:
                        p.lul_pos = newPos
                ch = p
                p = p.parent
        for n in self.iterNodesNoRoot():
            p = n.parent
            if p == self.root:
                n.br.lenL = float(n.lul_pos) / float(levels) * rootToLeaf
            else:
                n.br.lenL = float(n.lul_pos - p.lul_pos) / \
                    float(levels) * rootToLeaf

        # Clean up
        for n in self.iterNodes():
            if hasattr(n, 'lul_pos'):
                del(n.lul_pos)
            if overWriteBrLens:
                if n.br and hasattr(n.br, 'lenL'):
                    n.br.len = n.br.lenL
                    del(n.br.lenL)

    def removeEverythingExceptCladeAtNode(self, specifier):
        """Like it says.  Leaves a tree with a root-on-a-stick."""

        theNode = self.node(specifier)
        if theNode == self.root:
            return

        self.reRoot(theNode.parent, checkBiRoot=False)
        toRemoves = []
        for ch in self.root.iterChildren():
            if ch != theNode:
                toRemoves.append(ch)
        for ch in toRemoves:
            self.removeNode(
                ch, alsoRemoveSingleChildParentNode=False, alsoRemoveBiRoot=False)

    def dupeSubTree(self, dupeNodeSpecifier, up, doBrLens=True, doSupport=True):
        """Makes and returns a new Tree object, duping part of self.

        The dupeNodeSpecifier can be a node name, node number, or node
        object.

        Arg 'up' should be True or False.

        The returned subtree has a root-on-a-stick.

        BrLens are not duped -- brLens are default in the new subtree.

        So if the tree is like this::

                     +---------2:A
            +--------1
            |        +---------3:B
            |
            0--------4:C
            |
            |        +---------6:D
            +--------5:dupeNode
                     |         +--------8:E
                     +---------7
                               +--------9:F

        Then the subtree from node 5 up is::

                                 +---------2:D
            subTreeRoot:0--------1:dupeNode
                                 |         +--------4:E
                                 +---------3
                                           +--------5:F

        and the subtree from node 5 down is::

                                        +--------2:A
                              +---------1
            dupeNode:0--------5         +--------3:B
                              |
                              +---------4:C

        """

        dupeNode = self.node(dupeNodeSpecifier)
        if dupeNode == self.root and not up:
            print "The dupeNode is self.root, and you want a subtree below that?!?"
            sys.exit()
        from tree import Tree
        st = Tree()
        if up:
            n = Node()
            n.br = None
            n.nodeNum = 0
            n.name = 'subTreeRoot'
            n.isLeaf = 1
            st.root = n
            st.nodes.append(n)

            if dupeNode.isLeaf:
                # its easy -- just dupe the leaf
                n = Node()
                if doBrLens:
                    n.br.len = dupeNode.br.len
                if doSupport:
                    n.br.support = dupeNode.br.support
                n.parent = st.root
                st.root.leftChild = n
                n.name = dupeNode.name
                n.isLeaf = 1
                st.nodes.append(n)
            else:
                i = 1  # nodeNum counter
                nodeNumDict = {}
                for selfNode in dupeNode.iterPreOrder():
                    n = Node()
                    if doBrLens:
                        n.br.len = selfNode.br.len
                    if doSupport and hasattr(selfNode.br, 'support'):
                        n.br.support = selfNode.br.support
                    n.nodeNum = i
                    nodeNumDict[selfNode.nodeNum] = i
                    n.isLeaf = selfNode.isLeaf
                    n.name = selfNode.name
                    i += 1
                    st.nodes.append(n)
                for selfNode in dupeNode.iterPreOrder():
                    if selfNode.leftChild:
                        st.nodes[nodeNumDict[selfNode.nodeNum]].leftChild = \
                            st.nodes[nodeNumDict[selfNode.leftChild.nodeNum]]
                    if selfNode == dupeNode:
                        pass  # skip parent and sibling
                    else:
                        # parents and siblings.  There will always be a parent.
                        st.nodes[nodeNumDict[selfNode.nodeNum]].parent = \
                            st.nodes[nodeNumDict[selfNode.parent.nodeNum]]
                        if selfNode.sibling:
                            st.nodes[nodeNumDict[selfNode.nodeNum]].sibling = \
                                st.nodes[nodeNumDict[selfNode.sibling.nodeNum]]
                st.root.leftChild = st.nodes[nodeNumDict[dupeNode.nodeNum]]
                st.root.leftChild.parent = st.root
        else:  # down
            i = 0  # nodeNum counter
            nodeNumDict = {}
            mostRecentDown = None
            for selfNode in dupeNode.iterDown():
                n = Node()
                if selfNode != self.root:
                    if doBrLens:
                        n.br.len = selfNode.br.len
                    if not selfNode.isLeaf and doSupport:
                        if hasattr(selfNode.br, 'support'):
                            n.br.support = selfNode.br.support
                        else:
                            print "Warning: dupeSubTree() doSupport is turned on, but node %i has no support." % selfNode.nodeNum
                n.nodeNum = i
                nodeNumDict[selfNode.nodeNum] = i
                n.isLeaf = selfNode.isLeaf
                n.name = selfNode.name
                i += 1
                st.nodes.append(n)
            for selfNode in dupeNode.iterDown():
                # parents and siblings.
                if selfNode.parent:
                    st.nodes[nodeNumDict[selfNode.nodeNum]].parent = \
                        st.nodes[nodeNumDict[selfNode.parent.nodeNum]]
                else:
                    # Its the root
                    st.root = st.nodes[nodeNumDict[selfNode.nodeNum]]
                    st.root.br = None
                if selfNode.sibling:
                    st.nodes[nodeNumDict[selfNode.nodeNum]].sibling = \
                        st.nodes[nodeNumDict[selfNode.sibling.nodeNum]]
                if selfNode == dupeNode:
                    st.nodes[nodeNumDict[dupeNode.nodeNum]].leftChild = None
                    st.nodes[nodeNumDict[dupeNode.nodeNum]].isLeaf = 1
                else:
                    if selfNode.leftChild:
                        st.nodes[nodeNumDict[selfNode.nodeNum]].leftChild = \
                            st.nodes[nodeNumDict[selfNode.leftChild.nodeNum]]
            st.reRoot(nodeNumDict[dupeNode.nodeNum], moveInternalName=False)
            st.root.isLeaf = 1
            if not st.root.name:
                st.root.name = 'subTreeRoot'
        st.setPreAndPostOrder()
        return st

    def addSubTree(self, selfNode, theSubTree, subTreeTaxNames=None):
        """Add a subtree to a tree.

        The nodes from theSubTree are added to self.nodes, and theSubTree
        is deleted.

        If subTreeTaxNames is provided, fine, but if not this method can
        find them.  Providing them saves a bit of time, I assume.
        """

        #    subTreeRootNode:0-------1:oldSubTreeRootNodeLeftChild

        #                      +-------1:A
        #    oldSelfNodeParent:0
        #                      +-------2:selfNode

        #                    +--------1:A
        #  oldSelfNodeParent:0
        #                    |        +--------2:selfNode
        #                    +--------3:subTreeRootNode
        #                             +--------4:oldSubTreeRootNodeLeftChild

        # =========================================

        #                      +-------1:selfNode
        #    oldSelfNodeParent:0
        #                      +-------2:A

        #                             +--------1:selfNode
        #                    +--------3:subTreeRootNode
        #  oldSelfNodeParent:0        +--------4:oldSubTreeRootNodeLeftChild
        #                    |
        #                    +--------2:A

        assert selfNode in self.nodes
        assert selfNode.parent
        # its a root on a stick
        assert theSubTree.root.leftChild and not theSubTree.root.leftChild.sibling
        if not subTreeTaxNames:
            subTreeTaxNames = [n.name for n in theSubTree.iterLeavesNoRoot()]

        oldSelfNodeParent = selfNode.parent

        subTreeRootNode = theSubTree.root
        theSubTree.root = None
        oldSubTreeRootNodeLeftChild = subTreeRootNode.leftChild

        subTreeRootNode.parent = oldSelfNodeParent
        subTreeRootNode.leftChild = selfNode
        selfNode.parent = subTreeRootNode
        if oldSelfNodeParent.leftChild == selfNode:
            oldSelfNodeParent.leftChild = subTreeRootNode
        else:
            oldCh = oldSelfNodeParent.leftChild
            while oldCh.sibling != selfNode:
                oldCh = oldCh.sibling
            oldCh.sibling = subTreeRootNode
        if selfNode.sibling:
            subTreeRootNode.sibling = selfNode.sibling
        selfNode.sibling = oldSubTreeRootNodeLeftChild

        subTreeRootNode.isLeaf = 0

        nSelfNodes = len(self.nodes)
        self.nodes += theSubTree.nodes
        theSubTree.nodes = []

        selfNode.br.len = 0.1
        subTreeRootNode.br = NodeBranch()
        subTreeRootNode.br.len = 0.1
        self.taxNames += theSubTree.taxNames
        for i in range(nSelfNodes, len(self.nodes)):
            n = self.nodes[i]
            n.nodeNum = i
        # print
        # for n in self.nodes:
        #    print n.nodeNum
        # print "self.taxNames is %s" % self.taxNames
        # print "subTreeTaxNames %s" % subTreeTaxNames
        if self.taxNames:
            self._taxNames += subTreeTaxNames

        if self._nTax:
            self._nTax += len(subTreeTaxNames)
        self.preAndPostOrderAreValid = False
        self.preOrder = None
        self.postOrder = None
        self.setPreAndPostOrder()
        del(theSubTree)

    def addLeaf(self, attachmentNode, taxName):
        """Add a leaf to a tree.

        The leaf is added to the branch leading from the specified node.
        A new node is made on that branch, so actually 2 nodes are added
        to the tree.  The new leaf node is returned.
        """

        assert attachmentNode in self.nodes
        aNode = attachmentNode
        # self.draw()
        # print "attachmentNode is %i" % aNode.nodeNum
        #bNode = self.addNodeBetweenNodes(aNode, aNode.parent)
        aNodeP = attachmentNode.parent
        bNode = Node()
        bNode.nodeNum = len(self.nodes)
        self.nodes.append(bNode)
        bNode.parent = aNodeP
        bNode.leftChild = aNode
        aNode.parent = bNode
        if aNodeP.leftChild == aNode:
            aNodeP.leftChild = bNode
        else:
            oldCh = aNodeP.leftChild
            while oldCh.sibling != aNode:
                oldCh = oldCh.sibling
            oldCh.sibling = bNode
        if aNode.sibling:
            bNode.sibling = aNode.sibling
            aNode.sibling = None

        aNode.br.len = 0.1
        bNode.br.len = 0.1
        # self.draw()
        n = Node()
        n.nodeNum = len(self.nodes)
        n.name = taxName
        n.isLeaf = 1
        oldSibling = attachmentNode.sibling
        n.parent = bNode
        aNode.sibling = n
        n.sibling = oldSibling
        self.nodes.append(n)
        if self.taxNames:
            self.taxNames.append(n.name)
        self.preAndPostOrderAreValid = False
        self.preOrder = None
        self.postOrder = None
        self.getPreAndPostOrderAboveRoot()
        # self.dump(node=True)
        # self.draw()
        if self._nTax:
            self._nTax += 1
        return n

    def addSibLeaf(self, attachmentNode, taxName):
        """Add a leaf to a tree as a sibling, by specifying its parent.

        The leaf is added so that its parent is the specified node (ie
        attachmentNode), adding the node as a rightmost child to that
        parent.  The attachmentNode should not be a leaf -- it must have
        children nodes, to which the new leaf can be added as a sibling.

        The new node is returned.
        """

        assert attachmentNode in self.nodes
        assert not attachmentNode.isLeaf
        aNode = attachmentNode
        # self.draw()
        # print "attachmentNode is %i" % aNode.nodeNum
        n = Node()
        n.nodeNum = len(self.nodes)
        n.name = taxName
        n.isLeaf = 1
        n.br.len = 0.1
        rtMostCh = aNode.rightmostChild()
        rtMostCh.sibling = n
        n.sibling = None
        n.parent = aNode
        self.nodes.append(n)
        if self.taxNames:
            self.taxNames.append(n.name)
        self.preAndPostOrderAreValid = False
        self.preOrder = None
        self.postOrder = None
        self.getPreAndPostOrderAboveRoot()
        # self.dump(node=True)
        # self.draw()
        if self._nTax:
            self._nTax += 1
        return n

    def subTreeIsFullyBifurcating(self, theNode, up=True):
        """Is theNode and everything above it (or below it) bifurcating?

        Arg *up* says whether its above or below.
        """

        # don't use getNChildren() -- too slow!
        assert theNode != self.root
        if up:
            # Includes theNode, which we want
            for n in theNode.iterInternals():
                # if n.getNChildren() != 2:
                #    return False
                if n.leftChild and n.leftChild.sibling:
                    if n.leftChild.sibling.sibling:
                        return False
                else:
                    return False
            return True
        else:
            # Includes theNode, which we do not want
            for n in theNode.iterDown():
                if n != theNode:
                    if not n.isLeaf:
                        if n == self.root:
                            # if n.getNChildren() != 3:
                            #    return False
                            if n.leftChild and n.leftChild.sibling and n.leftChild.sibling.sibling:
                                if n.leftChild.sibling.sibling.sibling:
                                    return False
                            else:
                                return False
                        else:
                            # if n.getNChildren() != 2:
                            #    return False
                            if n.leftChild and n.leftChild.sibling:
                                if n.leftChild.sibling.sibling:
                                    return False
                            else:
                                return False
            return True

    def nni(self, upperNodeSpec=None):
        """Simple nearest-neighbor interchange.

        You specify an 'upper' node, via an upperNodeSpec, which as usual
        can be a node name, node number, or node object.  If you don't
        specify something, a random node will be chosen for you.  (This
        latter option might be a little slow if you are doing many of
        them, as it uses iterInternalsNoRoot(), but mostly it should be
        fast enough).

        The upper node has a parent -- the 'lower' node.  One subtree from
        the upper node and one subtree from the lower node are exchanged.
        Both subtrees are chosen randomly.

        This works on biRooted trees also, preserving the biRoot.
        """

        gm = ["Tree.nni()"]
        if upperNodeSpec:
            # This makes sure that upperNode is part of self.
            upperNode = self.node(upperNodeSpec)
        else:
            candidates = [n for n in self.iterInternalsNoRoot()]
            upperNode = random.choice(candidates)

        # Want the upperNode to have at least 2 children
        upperChildren = [n for n in upperNode.iterChildren()]
        if len(upperChildren) < 2:
            gm.append("upperNode needs to have at least 2 children.")
            raise P4Error(gm)
        if upperNode.parent:
            lowerNode = upperNode.parent
        else:
            gm.append("upperNode needs to have a parent node.")
            raise P4Error(gm)
        lowerChildren = [n for n in lowerNode.iterChildren() if n != upperNode]
        if 0:
            if lowerNode.parent:
                if len(lowerChildren) < 1:
                    gm.append("The lower node has a parent.")
                    gm.append(
                        "It needs at least one more child besides the upperNode.")
                    raise P4Error(gm)
            else:
                if len(lowerChildren) < 2:
                    gm.append("The lower node does not have a parent.")
                    gm.append(
                        "It needs at least 2 children besides the upperNode.")
                    raise P4Error(gm)
        if len(lowerChildren) < 1:
            gm.append(
                "The lower node needs at least one more child besides the upperNode.")
            raise P4Error(gm)

        upperSubTreeNode = random.choice(upperChildren)
        lowerSubTreeNode = random.choice(lowerChildren)

        upperSubTree = self.pruneSubTreeWithoutParent(
            upperSubTreeNode, allowSingleChildNode=True)
        lowerSubTree = self.pruneSubTreeWithoutParent(
            lowerSubTreeNode, allowSingleChildNode=True)

        self.reconnectSubTreeWithoutParent(upperSubTree, lowerNode)
        self.reconnectSubTreeWithoutParent(lowerSubTree, upperNode)

        self.setPreAndPostOrder()

    def checkThatAllSelfNodesAreInTheTree(self, verbose=False, andRemoveThem=False):
        """Check that all self.nodes are actually part of the tree.

        Arg *andRemoveThem* will remove those nodes, renumber the nodes,
        and reset pre- and postOrder, and return None

        If *andRemoveThem* is not set (the default is not set) then this
        method returns the list of nodes that are in self.nodes but not in
        the tree.
        """
        selfNodesSet = set(self.nodes)
        treeSet = set([n for n in self.iterNodes()])
        inSelfNodesButNotInTree = list(selfNodesSet.difference(treeSet))
        if verbose:
            if inSelfNodesButNotInTree:
                print "These nodes are in self.nodes, but not part of the tree."
                for n in inSelfNodesButNotInTree:
                    print n.nodeNum
            else:
                print "All nodes in self.nodes are also in the tree."
        if inSelfNodesButNotInTree and andRemoveThem:
            for n in inSelfNodesButNotInTree:
                self.nodes.remove(n)
            for i in range(len(self.nodes)):
                self.nodes[i].nodeNum = i
            self.preOrder = numpy.array(
                [var.NO_ORDER] * len(self.nodes), numpy.int32)
            self.postOrder = numpy.array(
                [var.NO_ORDER] * len(self.nodes), numpy.int32)

            if len(self.nodes) > 1:
                self.setPreAndPostOrder()

            return
        return list(inSelfNodesButNotInTree)

    def spr(self, pruneNode=None, above=True, graftNode=None):
        """Subtree pruning and reconnection.

        See also the Tree.randomSpr() method.  It uses this method to do a
        random spr move.

        This only works on fully bifurcating trees.  Doing spr moves would
        tend to break up polytomies anyway; pruning subtrees from a
        polytomy would require creation of new nodes.

        The subtree to be pruned might be pointing up or pointing down
        from a specified node.  If the subtree is pointing up, the subtree
        to be pruned is specified by the appropriate child of the root of
        the subtree; the subtree would have a root-on-a-stick (Is
        monofurcating a proper word?) with the subtree root's single child
        being the specified node.  If the subtree is pointing down, then
        the tree is re-rooted to the specified node to allow pruning of
        the subtree, now above the specified node, with the specified node
        as the root, including the subtree with the pre-re-rooting parent
        of the specified node.

        I'll draw that out.  Lets say we want to prune the subtree below
        node 2 in this tree.  That would include nodes 0, 1, 2, and 7. ::

            +--------1:A
            |
            |        +---------3:B
            |--------2
            0        |         +--------5:C
            |        +---------4
            |                  +--------6:D
            |
            +--------7:E

        The way it is done in this method is to re-root at node 2, which
        is the specified node.  Then the subtree including the
        pre-re-rooting parent of the specified node, ie node 0, is pruned.
        ::

            +--------3:B
            |
            |        +--------5:C
            2--------4
            |        +--------6:D
            |
            |        +--------1:A
            +--------0
                     +--------7:E

        """

        gm = ["Tree.spr()"]

        # This is only for fully bifurcating trees.
        if not self.isFullyBifurcating():
            gm.append("This method is only for fully bifurcating trees.")
            raise P4Error(gm)

        pnNode = self.node(pruneNode)
        grNode = self.node(graftNode)

        # self.draw()

        self.deleteCStuff()

        if not above:
            preReRootingParent = pnNode.parent
            self.reRoot(pnNode)
            pnNode = preReRootingParent

        # Check for silliness
        assert pnNode != self.root
        assert grNode != self.root
        assert pnNode != grNode
        subTreeNodes = [n for n in pnNode.iterPreOrder()]
        pnNodeParnt = pnNode.parent
        subTreeNodes.append(pnNodeParnt)
        # print [n.nodeNum for n in subTreeNodes]
        if grNode in subTreeNodes:
            if above:
                gm.append("grNode %i is part of the pruned subtree from %i-%i up.  No workee!" % (
                    grNode.nodeNum, pnNodeParnt.nodeNum, pnNode.nodeNum))
            else:
                gm.append("grNode %i is part of the pruned subtree below %i-%i.  No workee!" % (
                    grNode.nodeNum, pnNodeParnt.nodeNum, pnNode.nodeNum))

            raise P4Error(gm)

        # Prune it from the tree.
        if pnNodeParnt == self.root:
            if grNode.parent == self.root:  # as well,
                gm.append(
                    "prune node and graft node both have root as parent -- ie same origin and destination.")
                raise P4Error(gm)
            # Check if removal of the subtree will result in only 2 taxa
            singles = 0
            for ch in self.root.iterChildren():
                if ch != pnNode:
                    if not ch.leftChild:
                        singles += 1
            if singles == 2:
                gm.append(
                    "Removing subtree at %i will leave only 2 taxa." % pnNode.nodeNum)
                raise P4Error(gm)

            newRoot = None
            for ch in self.root.iterChildren():
                if ch != pnNode:
                    if ch.leftChild:
                        newRoot = ch
                        break
            # print "rerooting to node %i" % newRoot.nodeNum
            self.reRoot(newRoot)
            # self.draw()

        # pnNodeParnt is not the root, and so we can be sure pnNodeParnt has a parent.
        # the pnNode is the left, middle, or right child of the parent
        if pnNodeParnt.leftChild == pnNode:
            # its the left child
            # remove the pnNodeParnt along with the subtree
            # (pnNodeParnt.parent.leftChild, (pnNode, pnNode.sibling)pnNodeParnt, pnNodeParnt.sibling)pnNodeParnt.parent;
            #                   +----------1:pnNodeParnt.parent.leftChild
            #                   |
            #                   |          +-----------3:pnNode
            # pnNodeParnt.parent:0----------2:pnNodeParnt
            #                   |          +-----------4:pnNode.sibling
            #                   |
            #                   +----------5:pnNodeParnt.sibling

            pnNode.sibling.parent = pnNodeParnt.parent

            if pnNodeParnt.parent.leftChild == pnNodeParnt:
                pnNodeParnt.parent.leftChild = pnNode.sibling
            elif pnNodeParnt.parent.leftChild.sibling == pnNodeParnt:
                pnNodeParnt.parent.leftChild.sibling = pnNode.sibling
            # if pnNodeParnt.parent is the root, it can have 3 children, and
            # maybe this ...
            else:
                pnNodeParnt.parent.leftChild.sibling.sibling = pnNode.sibling
            pnNode.sibling.sibling = pnNodeParnt.sibling
            #pnNode.sibling = None
            #pnNodeParnt.parnt = None
        else:
            # its a right child
            # (pnNodeParnt.parent.leftChild, (pnNodeLeftSib, pnNode)pnNodeParnt, pnNodeParnt.sibling)pnNodeParnt.parent;
            #                    +-----------1:pnNodeParnt.parent.leftChild
            #                    |
            #                    |           +-----------3:pnNodeLeftSib
            # pnNodeParnt.parent:0-----------2:pnNodeParnt
            #                    |           +-----------4:pnNode
            #                    |
            #                    +-----------5:pnNodeParnt.sibling

            pnNodeLeftSib = pnNodeParnt.leftChild
            pnNodeLeftSib.parent = pnNodeParnt.parent
            pnNodeLeftSib.sibling = pnNodeParnt.sibling
            if pnNodeParnt.parent.leftChild == pnNodeParnt:
                pnNodeParnt.parent.leftChild = pnNodeLeftSib
            elif pnNodeParnt.parent.leftChild.sibling == pnNodeParnt:
                pnNodeParnt.parent.leftChild.sibling = pnNodeLeftSib
            # if pnNodeParnt.parent is the root, it can have 3 children, and
            # maybe this ...
            else:
                pnNodeParnt.parent.leftChild.sibling.sibling = pnNodeLeftSib

        # To look at the pruned tree, before grafting ...
        if 0:
            print "removal of subtree at %i-%i gives .." % (pnNodeParnt.nodeNum, pnNode.nodeNum)
            self._nTax = 0
            self.preAndPostOrderAreValid = 0
            # This won't work unless preAndPostOrderAreValid set to 0
            self.draw()

        # Now graft it back on ...

        #  (grNodeParnt.leftChild, grNode, grNode.sibling)grNodeParnt;
        #             +-------1:grNodeParnt.leftChild
        #             |
        # grNodeParnt:0-------2:grNode
        #             |
        #             +-------3:grNode.sibling

        # (grNode, grNode.sibling)grNodeParnt;

        #             +-------1:grNode
        # grNodeParnt:0
        #             +-------2:grNode.sibling

        # (grNodeParnt.leftChild, grNodeParnt.leftChild.sibling, grNode)grNodeParnt;
        #             +-------1:grNodeParnt.leftChild
        #             |
        # grNodeParnt:0-------2:grNodeParnt.leftChild.sibling
        #             |
        #             +-------3:grNode

        grNodeParnt = grNode.parent
        pnNodeParnt.parent = grNodeParnt
        grNode.parent = pnNodeParnt
        pnNodeParnt.sibling = grNode.sibling
        grNode.sibling = pnNode
        pnNode.sibling = None
        grNode.sibling = pnNode
        pnNodeParnt.leftChild = grNode
        if grNodeParnt.leftChild == grNode:
            grNodeParnt.leftChild = pnNodeParnt
        elif grNodeParnt.leftChild.sibling == grNode:
            grNodeParnt.leftChild.sibling = pnNodeParnt
        else:
            grNodeParnt.leftChild.sibling.sibling = pnNodeParnt

        # To look at the tree after grafting ...
        if 0:
            print "grafting the subtree at grNode %i gives ..." % grNode.nodeNum
            self._nTax = 0
            self.preAndPostOrderAreValid = 0
            # This won't work unless preAndPostOrderAreValid set to 0
            self.draw()

        self.preAndPostOrderAreValid = 0
        self.setPreAndPostOrder()

    def randomSpr(self):
        """Do a random spr move.

        """

        myAbove = random.choice([True, False])
        tNodes = [n for n in self.iterNodesNoRoot()]
        while 1:
            try:
                pNode = random.choice(tNodes)
            except IndexError:
                # we have run out of choices, all were unsuitable.
                return
            tNodes.remove(pNode)
            if pNode == self.root:
                continue
            if myAbove:
                # iterDown() will return pNode, the root, and whatever else.
                # If whatever else is only 2 nodes, it won't work.
                nodesDown = [n2 for n2 in pNode.iterDown()]
                if len(nodesDown) <= 4:
                    continue
            else:
                # iterPostOrder() ends with the pNode.  We need at least the
                # pNode and 3 others, so total 3 is too few.
                nodesUp = [n2 for n2 in pNode.iterPostOrder()]
                # self.draw()
                # print "--", pNode.nodeNum, [n2.nodeNum for n2 in nodesUp]
                if len(nodesUp) <= 3:
                    continue
            break

        if myAbove:
            #pNode = self.node(pNNum)
            possibles = [n2.nodeNum for n2 in pNode.iterDown()
                         if n2 is not pNode
                         and n2 is not pNode.parent
                         and n2.parent is not pNode.parent]
        else:
            #pNode = self.node(pNNum)
            possibles = [n2.nodeNum for n2 in pNode.iterPreOrder()
                         if n2.parent is not pNode and
                         pNode is not n2]
        if 0:
            self.draw()
            print "above=%s, pNode %i, " % (myAbove, pNode.nodeNum), subtreeNodeNums
            print possibles
            sys.exit()

        safety = 0
        giveUp = 20
        while 1:
            safety += 1
            if safety >= giveUp:
                break
            gNNum = random.choice(possibles)
            if gNNum == pNode.nodeNum:
                continue
            if gNNum == self.root.nodeNum:
                continue
            # print "===", "pNode=%i" % pNode.nodeNum, "gNNum=%i" % gNNum,
            # subtreeNodeNums
            if pNode.parent == self.root and self.node(gNNum).parent == self.root:
                continue
            if not myAbove and self.node(gNNum).parent == pNode:
                continue
            break
        if safety >= giveUp:
            print "pNode is %i, above=%s" % (pNode.nodeNum, myAbove)
            self.draw()
            raise P4Error
        # print "spr()  pruneNum %i, above=%s, graftNum %i" % (pNode.nodeNum,
        # myAbove, gNNum)

        self.spr(pruneNode=pNode, above=myAbove, graftNode=gNNum)

    def inputTreesToSuperTreeDistances(self, inputTrees, doSd=True, doScqdist=True):
        """Return the topology distance between input trees and pruned self.

        Either or both of two distances, RF ('sd') and quartet ('scqdist')
        distances, are returned.

        See also the Trees method of the same name that does a bunch of them.
        """

        for inTree in inputTrees:
            if not inTree.taxNames:
                inTree.taxNames = [n.name for n in inTree.iterLeavesNoRoot()]
            # inTree.makeSplitKeys()

        totalSd = 0
        totalScqdist = 0
        for inTree in inputTrees:
            sDupe = self.dupe()
            toRemove = []
            for n in sDupe.iterLeavesNoRoot():
                if n.name not in inTree.taxNames:
                    toRemove.append(n)
            for n in toRemove:
                sDupe.removeNode(n)
            sDupe.taxNames = inTree.taxNames
            if doSd:
                rfDist = sDupe.topologyDistance(
                    inTree, metric='sd', resetSplitKeySet=True)
                totalSd += rfDist
            if doScqdist:
                qd = sDupe.topologyDistance(inTree, metric='scqdist')
                totalScqdist += qd
        if doSd and doScqdist:
            return totalSd, totalScqdist
        elif doSd:
            return totalSd
        elif doScqdist:
            return totalScqdist

    def patristicDistanceMatrix(self):
        """Matrix of distances along tree path.

        This method sums the branch lengths between each pair of taxa, and
        puts the result in a DistanceMatrix object, which is returned.

        Self.taxNames is required.
        """

        gm = ['Tree.patristicDistanceMatrix()']

        if not self.taxNames:
            gm.append("No taxNames.")
            raise P4Error(gm)

        # The tree will be rearranged, so make a copy to play with, so
        # self is undisturbed.
        t = self.dupe()

        d = DistanceMatrix()
        d.names = self.taxNames
        d.dim = len(d.names)
        # print d.names
        d.matrix = []
        for i in range(d.dim):
            d.matrix.append([0.0] * d.dim)

        for i in range(d.dim):
            n1 = t.node(d.names[i])
            t.reRoot(n1, moveInternalName=False)
            for j in range(i + 1, d.dim):
                n2 = t.node(d.names[j])
                # sum up distances between n1 and n2
                sum = n2.br.len
                n2 = n2.parent
                while n2 != n1:
                    sum = sum + n2.br.len
                    n2 = n2.parent
                # print "Dist from %s to %s is %f" % (d.names[i], d.names[j],
                # sum)
                d.matrix[i][j] = sum
                d.matrix[j][i] = sum
        return d

    def tPickle(self, fName=None):
        """Pickle self to a file with a 'p4_tPickle' suffix.

        If there is an attached Data object, it is not pickled.  If there
        is an attached model object, it is pickled.  Pointers to c-structs
        are not pickled.

        If fName is supplied, the file name becomes fName.p4_tPickle,
        unless fName already ends with .p4_tPickle.  If fName is not
        supplied, self.name is used in fName's place.  If neither is
        supplied, the pid is used as fName.

        If a file with the chosen name already exists, it is silently
        over-written!

        p4 can read a p4_tPickle file from the command line or using the
        read() function, as usual.

        (This would not be a good option for long-term storage, because if
        you upgrade p4 and the p4 Classes change a lot then it may become
        impossible to unpickle it.  If that happens, you can use the old
        version of p4 to unpickle.)

        """

        import cPickle
        #suffix = '.p4_%s_tPickle' % var.versionString
        suffix = '.p4_tPickle'

        if not self.name and not fName:
            fName = '%s' % os.getpid()

        if fName:
            if fName.endswith(suffix):
                fN = fName
            else:
                fN = fName + suffix
        elif self.name:
            if self.name.endswith(suffix):
                fN = self.name
            else:
                fN = self.name + suffix
        f = file(fN, 'w')
        cPickle.dump(self.dupe(), f, 1)  # 1 for binary
        # cPickle.dump(self, f, 1) # Don't do this -- has pointers that would
        # not have been malloc'ed!  And data!
        f.close()

    def writeNexus(self, fName=None, append=0, writeTaxaBlockIfTaxNamesIsSet=1, message=None):
        """Write the tree out in Nexus format, in a trees block.

        If fName is None, the default, it is written to sys.stdout.

        #NEXUS is written unless we are appending-- set append=1.

        If you want to write with a translation, use a Trees object.
        """

        gm = ['Tree.writeNexus()']

        if fName == None or fName == sys.stdout:
            f = sys.stdout
            if not append:
                f.write('#NEXUS\n\n')
        else:
            if append:
                import os
                if os.path.isfile(fName):
                    try:
                        f = open(fName, 'a')
                    except IOError:
                        gm.append("Can't open %s for appending." % fName)
                        raise P4Error(gm)
                else:
                    if 0:
                        print "Tree.writeNexus()"
                        print "    'append' is requested,"
                        print "    but '%s' is not a regular file (doesn't exist?)." \
                              % fName
                        print "    Writing to a new file instead."
                    try:
                        f = open(fName, 'w')
                        f.write('#NEXUS\n\n')
                    except IOError:
                        gm.append("Can't open %s for writing." % fName)
                        raise P4Error(gm)

            else:
                try:
                    f = open(fName, 'w')
                    f.write('#NEXUS\n\n')
                except IOError:
                    gm.append("Can't open %s for writing." % fName)
                    raise P4Error(gm)

        if writeTaxaBlockIfTaxNamesIsSet and self.taxNames:
            f.write('begin taxa;\n')
            f.write('  dimensions ntax=%s;\n' % self.nTax)
            f.write('  taxlabels')
            for i in self.taxNames:
                f.write(' %s' % func.nexusFixNameIfQuotesAreNeeded(i))
            f.write(';\nend;\n\n')

        f.write('begin trees;\n')
        if message:
            f.write('  [%s\n  ]\n' % message)
        if self.logLike:
            f.write('  [logLike for tree %s is %f]\n' %
                    (self.name, self.logLike))

        f.write('  tree %s = [&U] ' %
                func.nexusFixNameIfQuotesAreNeeded(self.name))
        if self.recipWeight:
            if self.recipWeight == 1:
                f.write('[&W 1] ')
            else:
                f.write('[&W 1/%i] ' % self.recipWeight)
        if hasattr(self, 'weight'):
            f.write('[&W %f] ' % self.weight)

        self.writeNewick(f)
        f.write('end;\n\n')
        if f != sys.stdout:
            f.close()

    def write(self):
        """This writes out the Newick tree description to sys.stdout."""
        self.writeNewick(sys.stdout)

    def writePhylip(self, fName=None, withTranslation=0, translationHash=None, doMcmcCommandComments=0):
        """Write the tree in Phylip or Newick format.

        (This method is just a dupe of writeNewick().  Without the
        'toString' or 'append' args.)

        fName may also be an open file object.
        """
        self.writeNewick(
            fName, withTranslation, translationHash, doMcmcCommandComments)

    def writeNewick(self, fName=None, withTranslation=0, translationHash=None, doMcmcCommandComments=0, toString=False, append=False):
        """Write the tree in Newick, aka Phylip, format.

        This is done in a Nexus-oriented way.  If taxNames have spaces or
        odd characters, they are single-quoted.  There is no restriction
        on the length of the taxon names.  A translationHash can be used.

        fName may also be an open file object.

        If 'toString' is turned on, then 'fName' should be None, and a Newick
        representation of the tree is returned as a string.

        The method 'writePhylip()' is the same as this, with fewer arguments.
        """

        gm = ['Tree.writeNewick()']

        sList = []
        if withTranslation and not translationHash:
            gm.append("No translationHash.")
            raise P4Error(gm)

        if fName and toString:
            gm.append(
                "You cannot write to a file and string at the same time.")
            raise P4Error(gm)

        if doMcmcCommandComments:
            if not self.model:
                gm.append("No model attached to tree.")
                gm.append("Set doMcmcCommandComments=0")
                raise P4Error(gm)

        # print 'self.preAndPostOrderAreValid = %s' %
        # self.preAndPostOrderAreValid
        if not self.preAndPostOrderAreValid:
            self.setPreAndPostOrder()
        # print "self.preOrder = %s" % self.preOrder

        # Don't count un-used nodes.
        nNodes = len([n for n in self.iterNodes()])
        # print "nNodes = %i" % nNodes

        if nNodes == 1:
            # print "Single node.  isLeaf=%s, name=%s" % (self.root.isLeaf,
            # self.root.name)
            if self.root.isLeaf:
                if withTranslation:
                    sList.append('%s' % translationHash[self.root.name])
                elif self.root.name:
                    sList.append(
                        '%s' % func.nexusFixNameIfQuotesAreNeeded(self.root.name))
                else:
                    sList.append('()')
            else:
                # Will this ever happen?
                gm.append(
                    "Something is wrong.  There is only one node, and it is not terminal.")
                raise P4Error(gm)

        elif nNodes > 1:
            writeBrLens = 0
            for n in self.iterNodesNoRoot():
                if n.br.len != 0.1:
                    writeBrLens = 1
                    break
            stack = []
            for n in self.iterPreOrder():
                stack.append(n)

                if n.leftChild:
                    sList.append('(')
                    continue

                while len(stack):
                    n1 = stack.pop()
                    # print "stacklen=%i, n1 name=%s" % (len(stack), n1.name)
                    if n1.isLeaf:
                        if n1 == self.root:
                            sList.append(')')
                        if withTranslation:
                            sList.append('%s' % (translationHash[n1.name]))
                        else:
                            if n1.name:
                                sList.append(
                                    '%s' % func.nexusFixNameIfQuotesAreNeeded(n1.name))
                            else:
                                if n1 != self.root:
                                    gm.append("Terminal node with no name?")
                                    raise P4Error(gm)
                    else:
                        sList.append(')')
                        if n1.name:
                            sList.append(
                                '%s' % func.nexusFixNameIfQuotesAreNeeded(n1.name))
                    if writeBrLens:
                        if n1 != self.root:
                            sList.append(':%g' % n1.br.len)
                    if doMcmcCommandComments:
                        sList.append(self._getMcmcCommandComment(n1))
                    if n1.sibling:
                        sList.append(', ')
                        break
        sList.append(';\n')
        if toString:
            return "".join(sList)
        elif fName == None:
            print "".join(sList)
        elif type(fName) == type('string'):
            if append:
                fName2 = file(fName, 'a')
            else:
                fName2 = file(fName, 'w')
            fName2.write(string.join(sList, ''))
    #        fName2.write('\n')
            fName2.close()
        elif hasattr(fName, 'write'):
            fName.write(string.join(sList, ''))
            # fName.write('\n')
            # Somebody else opened the fName, so somebody else can close it.
        else:
            gm.append(
                "I don't understand (%s) passed to me to write to." % fName)
            raise P4Error(gm)

    def _getMcmcCommandComment(self, theNode):
        sList = [' [&']
        for pNum in range(self.model.nParts):
            if self.model.parts[pNum].nComps > 1:
                sList.append(' c%i.%i' % (pNum, theNode.parts[pNum].compNum))
            if theNode != self.root:
                if self.model.parts[pNum].nRMatrices > 1:
                    sList.append(' r%i.%i' %
                                 (pNum, theNode.br.parts[pNum].rMatrixNum))
                if self.model.parts[pNum].nGdasrvs > 1:
                    sList.append(' g%i.%i' %
                                 (pNum, theNode.br.parts[pNum].gdasrvNum))
        sList.append(']')
        return string.join(sList, '')

    def draw(self, showInternalNodeNames=1, addToBrLen=0.2, width=None, showNodeNums=1, partNum=0, model=None):
        """Draw the tree to the screen.

        This method makes a text drawing of the tree and writes it to sys.stdout.

        Arg addToBrLen adds, by default 0.2, to each branch length, to
        make the tree more legible.  If you want the branch lengths more
        realistic, you can set it to zero, or better, use vector graphics
        for drawing the trees.

        Setting arg model aids in drawing trees with tree-hetero models.
        If the model characteristic (usually composition or rMatrix)
        differs over the tree, this method can draw it for you.

        See also :meth:`Tree.Tree.textDrawList`, which returns the drawing as a list of strings.

        See the method :meth:`Tree.Tree.setTextDrawSymbol`, which facilitates
        drawing different branches with different symbols.
        """

        s = self.textDrawList(showInternalNodeNames=showInternalNodeNames,
                              addToBrLen=addToBrLen, width=width,
                              autoIncreaseWidth=True,
                              showNodeNums=showNodeNums, partNum=partNum, model=model)
        print string.join(s, '\n')

    def textDrawList(self, showInternalNodeNames=1, addToBrLen=0.2, width=None, autoIncreaseWidth=True, showNodeNums=1, partNum=0, model=False):

        if len(self.nodes) == 0:
            return ['']
        elif len(self.nodes) == 1:
            if showNodeNums:
                return ['%i:%s' % (self.nodes[0].nodeNum, self.nodes[0].name)]
            else:
                return ['%s' % self.nodes[0].name]
        gm = ['Tree.textDrawList()']
        if not self.preAndPostOrderAreValid:
            self.setPreAndPostOrder()

        from treepicture import TreePicture
        p = TreePicture(self)
        p.fName = None
        p.width = width
        p.xScale = None
        p.yScale = 1
        p.pointsPerLetter = 1
        p.addToBrLen = addToBrLen
        p.textShowNodeNums = showNodeNums
        p.showInternalNodeNames = showInternalNodeNames
        if showNodeNums:
            p.nameOffset = 0
        else:
            p.nameOffset = 1
        p.xOrigin = 0.0
        p.yOrigin = 0.0
        p.textSize = 1
        p.labelTextSize = 1

        # If the width is not specified, make a guess
        if width == None:
            tLen = 0  # Longest number of horizontal sections to draw
            longestNameLen = 0
            for n in self.nodes:
                if n.isLeaf and n != self.root:
                    if n.name:  # This assumes short internal node name lengths
                        if len(n.name) > longestNameLen:
                            longestNameLen = len(n.name)
                    thisLen = 0
                    n1 = n
                    # print "x n1 is node %i" % n1.nodeNum
                    while n1 != self.root:
                        n1 = n1.parent
                        # print "y n1 is node %i" % n1.nodeNum
                        thisLen += 1
                    if thisLen > tLen:
                        tLen = thisLen
            # print "tLen =", tLen
            # print "longestNameLen =", longestNameLen
            rootNameLen = 0
            if self.root.name:
                rootNameLen = len(self.root.name) + 1
            p.width = (10 * tLen) + rootNameLen + longestNameLen
            if p.width > 100:
                p.width = 100

        # print "p.width =", p.width
        #import sys; sys.exit()

        # Make sure the names fit.
        if not autoIncreaseWidth:
            for n in self.nodes:
                if n.isLeaf and n != self.root:
                    if n.name and len(n.name) > p.width:
                        gm.append(
                            "There are long names, and the given width is not enough.")
                        raise P4Error(gm)

        if model:
            try:
                partNum = int(partNum)
            except:
                gm.append("partNum arg should be an integer.")
                raise P4Error(gm)
            if not self.model:
                gm.append("If model arg is set, then self.model must exist.")
                raise P4Error(gm)
            if partNum < 0 or partNum >= self.model.nParts:
                gm.append(
                    "Zero-based partNum %i is out of range of %s parts." % (partNum, self.model.nParts))
                raise P4Error(gm)
            p.partNum = partNum
            p.doModel = 1

            doComps = 1
            doRMatrices = 1
            if self.model.parts[partNum].nComps < 2:
                doComps = 0
            if self.model.parts[partNum].nRMatrices < 2:
                doRMatrices = 0

            if not doComps and not doRMatrices:
                p._setPos(autoIncreaseWidth)
                s = p.textString(returnAsList=True)
                s.append(
                    "Both the composition of the model and the rate matrix are homogeneous in part %i.\n" % partNum)
                return s

            # First do compositions
            s = []
            if doComps:
                if 0:
                    if self.model.nParts > 1:
                        print "Compositions for part %i" % partNum
                    else:
                        print "\nCompositions\n------------"
                p.textDrawModelThing = var.TEXTDRAW_COMP
                p._setPos(autoIncreaseWidth)
                s = p.textString(returnAsList=True)
                # print s

            # Then do RMatrices
            if doRMatrices:
                if 0:
                    if self.model.nParts > 1:
                        print "RMatrices for part %i" % partNum
                    else:
                        print "\nRMatrices\n---------"
                p.textDrawModelThing = var.TEXTDRAW_RMATRIX
                p._setPos(autoIncreaseWidth)
                if s:
                    s += p.textString(returnAsList=True)
                else:
                    s = p.textString(returnAsList=True)
                if self.model.parts[partNum].nComps == 1:
                    s.append(
                        "The composition of the model is homogeneous in part %i\n" % partNum)
            elif self.model.parts[partNum].nRMatrices == 1:
                s.append(
                    "The rate matrix is homogeneous in part %i\n" % partNum)
            # Don't bother with GDASRV, yet
            return s

        else:
            p._setPos(autoIncreaseWidth)
            s = p.textString(returnAsList=True)
            return s

    # outFileName=None, width=500, heightFactor=0.85, pointsPerLetter=6.0,
    # textSize=11, labelSize=9, putInternalNodeNamesOnBranches=0)

    def eps(self, fName=None, width=500, putInternalNodeNamesOnBranches=0):
        """Make a basic eps drawing of self.

        The 'width' is in postscript points.

        By default, internal node names label the node, where the node
        name goes on the right of the node.  You can make the node name
        label the branch by setting 'putInternalNodeNamesOnBranches'.

        """

        gm = ['Tree.eps()']

        if not self.preAndPostOrderAreValid:
            self.setPreAndPostOrder()

        from treepicture import TreePicture
        p = TreePicture(self)
        p.addToBrLen = 0.0
        p.width = width
        p.yScale = 17.0
        p.xScale = None
        p.pointsPerLetter = 6.0
        p.textSize = 11
        p.labelTextSize = 8
        p.putInternalNodeNamesOnBranches = putInternalNodeNamesOnBranches
        p.xOrigin = 5.0
        p._setPos()
        s = p.vectorString()

        if not fName:
            if self.name:
                fName = '%s.eps' % self.name
            else:
                fName = '%i.eps' % os.getpid()
        # if not fName.endswith('.eps'):
        #    fName = '%s.eps' % fName
        f = file(fName, 'w')
        f.write(s)
        f.close()

    #############################################
    # Various Tree methods for defining models
    #############################################

    def _setData(self, theData):
        """Sets self.data, and self.nParts"""

        # Two cases.  Either
        #   1.  Self already has a data.
        #                - deleteCStuff(), and replace the old data with the new data.
        #   2.  Self has no data.
        #     a.  Has a model.  So just check for compatibility.
        #     b.  Has no model, so has never seen a data before.

        complaintHead = '\nTree._setData()'
        if self.name:
            complaintHead += ", tree '%s'" % self.name
        if self.fName:
            complaintHead += ", file %s" % self.fName
        gm = [complaintHead]

        if isinstance(theData, Data) or theData == None:
            pass
        else:
            gm.append("Set data only to a Data object, or None, ok?")
            raise P4Error(gm)

        if self.data or self.model:
            # Normally there won't be anything to delete, but you never know...
            self.deleteCStuff()

        if not theData:
            self._data = None
            return

        if self.model:

            # We have seen a data object before (otherwise we would not
            # have been able to set the model).  Check for compatibility.

            # print "_setData() here.   self.model exits."

            if not self.taxNames:
                gm.append(
                    "Self has model, but no taxNames.  Programming error.")
                raise P4Error(gm)

            # Check for same number of taxa
            treeNTax = len(self.taxNames)
            dataNTax = len(theData.taxNames)
            if self.nTax != dataNTax:
                gm.append("The number of taxa in the tree (%s)" % treeNTax)
                gm.append("is not the same as in the data (%s)" % dataNTax)
                raise P4Error(gm)

            # Check for mis-matched taxNames
            isBad = 0
            for tn in self.taxNames:
                if tn not in theData.taxNames:
                    isBad = 1
                    break
            for tn in theData.taxNames:
                if tn not in self.taxNames:
                    isBad = 1
                    break
            if isBad:
                gm.append("TaxName mismatch between the tree and the data.")
                gm.append("Here they are, sorted to show mis-matches.")
                gm.append("    %25s    %25s" % ('data', 'tree'))
                self.taxNames.sort()
                theData.taxNames.sort()
                for i in range(len(self.taxNames)):
                    if theData.taxNames[i] == self.taxNames[i]:
                        gm.append("    %25s    %25s" %
                                  (theData.taxNames[i], self.taxNames[i]))
                    else:
                        gm.append("    %25s    %25s  ***" %
                                  (theData.taxNames[i], self.taxNames[i]))
                raise P4Error(gm)

            # Same number of parts
            if len(theData.parts) != self.model.nParts:
                gm.append("nParts mis-match.  len(theData.parts)=%i, model.nParts=%i" % (
                    len(theData.parts), self.model.nParts))
                raise P4Error(gm)

            # Check dims and symbols in the parts
            for pNum in range(self.model.nParts):
                if theData.parts[pNum].dim != self.model.parts[pNum].dim:
                    gm.append("Parts dim mis-match.")
                    raise P4Error(gm)
                if theData.parts[pNum].symbols != self.model.parts[pNum].symbols:
                    gm.append("Parts symbols mis-match.")
                    raise P4Error(gm)

            # Set seqNum
            for n in self.iterLeavesNoRoot():
                if n.seqNum != self.taxNames.index(n.name):
                    gm.append("seqNums do not match up with taxNames.")
                    raise P4Error(gm)

            self._data = theData

        else:
            # When you do this method, _setData(), self gets a suitable
            # model.  We have here no model, so we may have never seen a
            # data before.  Or we might have just lost the model.

            # Check for same number of taxa
            treeNTax = 0
            treeTaxNames = []
            for n in self.nodes:
                if n.isLeaf:
                    treeNTax += 1
                    treeTaxNames.append(n.name)
            dataNTax = len(theData.taxNames)
            if treeNTax != dataNTax:
                gm.append("The number of taxa in the tree (%s)" % treeNTax)
                gm.append("is not the same as in the data (%s)" % dataNTax)
                raise P4Error(gm)

            # Check for mis-matched taxNames
            isBad = 0
            for tn in treeTaxNames:
                if tn not in theData.taxNames:
                    isBad = 1
                    break
            for tn in theData.taxNames:
                if tn not in treeTaxNames:
                    isBad = 1
                    break
            if isBad:
                gm.append("TaxName mismatch between the tree and the data.")
                gm.append("Here they are, sorted to show mis-matches.")
                gm.append("    %25s  %25s" % ('data', 'tree'))
                treeTaxNames.sort()
                theData.taxNames.sort()
                for i in range(len(treeTaxNames)):
                    if theData.taxNames[i] == treeTaxNames[i]:
                        gm.append("    %25s  %25s" %
                                  (theData.taxNames[i], treeTaxNames[i]))
                    else:
                        gm.append("*** %25s  %25s" %
                                  (theData.taxNames[i], treeTaxNames[i]))
                raise P4Error(gm)

            # attach
            self.taxNames = theData.taxNames
            self._data = theData
            #self.nParts = len(theData.parts)

            # Now that nParts is known ...
            # print "_setData.  len(theData.parts) = %s" % len(theData.parts)
            # calls self.deleteCStuff()
            self.model = Model(len(theData.parts))
            for n in self.nodes:
                if n.parts:
                    n.parts = []
                for i in range(self.model.nParts):
                    n.parts.append(NodePart())
            for n in self.nodes:
                if n != self.root:
                    n.br.parts = []
                    for i in range(self.model.nParts):
                        n.br.parts.append(NodeBranchPart())

            # Set modelPart dims and symbols
            for pNum in range(self.model.nParts):
                self.model.parts[pNum].dim = theData.parts[pNum].dim
                self.model.parts[pNum].symbols = theData.parts[pNum].symbols

            # There is intentionally no default pInvar, forcing the user to be
            # explicit.
            for p in self.model.parts:
                p.pInvar = None

            # Set seqNum
            for n in self.nodes:
                if n.isLeaf:
                    n.seqNum = self.taxNames.index(n.name)

    data = property(lambda self: self._data, _setData)

    def _setModel(self, theModel):
        gm = ['Tree._setModel()']
        # print gm[0]
        # print "    Got '%s'" % theModel
        if isinstance(theModel, Model) or theModel == None:
            pass
        else:
            gm.append("Attempt to set Tree.model to '%s'.  " % theModel)
            gm.append(
                "Don't set the model to anything other than 'None' or a Model, ok?  ")
            gm.append("(And generally the user only sets it to None.)  ")
            raise P4Error(gm)
        # if theModel and self._model:  # Why do I do this?
        #    gm.append("The tree already has a model object; I am refusing to clobber it.")
        #    gm.append("Perhaps use a (perhaps duplicate) tree with no model.")
        #    raise P4Error(gm)
        if self.model or self.data:
            self.deleteCStuff()
            # print 'Tree._setModel()  finished deleteCStuff()'
        self._model = theModel

    def _delModel(self):
        gm = ['Tree._delModel()']
        gm.append("Caught an attempt to delete self.model, but")
        gm.append("self.model is a property, so you shouldn't delete it.")
        gm.append("But you can set it to None if you like.")
        raise P4Error(gm)

    model = property(lambda self: self._model, _setModel, _delModel)

    def _checkModelThing(self, partNum, symbol, complaintHead):
        gm = [complaintHead]
        if not self.data:
            gm.append("No data.  Set the data first.")
            raise P4Error(gm)

        if not self.model:
            # When you _setData(), a model object of suitable dimensions
            # is made and attached to self.  If we have got here, it is
            # because the model has subsequently been lost.  So just
            # re-instate it.
            self._setData(self.data)

        if partNum < 0 or partNum >= self.model.nParts:
            gm.append("Zero-based partNum (%s) is out of range (of %s parts)" %
                      (partNum, self.model.nParts))
            raise P4Error(gm)

        if symbol:
            if type(symbol) != type('s') or len(symbol) != 1:
                gm.append("Symbols must be 1-length strings.")
                raise P4Error(gm)
            if symbol == '?':
                gm.append("Got assigned text drawing symbol '?'.")
                gm.append(
                    "Don't use it-- it is reserved for missing modelThings")
                raise P4Error(gm)

    def newComp(self, partNum=0, free=0, spec='empirical', val=None, symbol=None):
        """Make, attach, and return a new Comp object.

        The arg *spec* should be a string, one of::

          'equal'          no val
          'empirical'      no val
          'specified'      val=[aList]
          'wag', etc       no val
             (ie one of the empirical protein models, including
             cpREV, d78, jtt, mtREV24, mtmam, wag, etc)

        If spec='specified', then you specify dim or dim-1 values in a
        list as the 'val' arg.

        This method returns a Comp object, which you can ignore if it is a
        tree-homogeneous model.  However, if it is a tree-hetero model
        then you may want to get that Comp object so that you can place
        it on the tree explicitly with setModelThing(), like this::

            c0 = newComp(partNum=0, free=1, spec='empirical')
            c1 = newComp(partNum=0, free=1, spec='empirical')
            myTree.setModelThing(c0, node=myTree.root, clade=1)
            myTree.setModelThing(c1, node=5, clade=1)
            myTree.setModelThing(c1, node=18, clade=0)

        Alternatively, you can simply let p4 place them randomly::

            newComp(partNum=0, free=1, spec='empirical')
            newComp(partNum=0, free=1, spec='empirical')
            myTree.setModelThingsRandomly()

        Calculation of probability matrices for likelihood calcs etc are
        wrong when there are any comp values that are zero, so that is not
        allowed.  Any zeros are converted to var.PIVEC_MIN, which is 1e-18
        this week.  Hopefully close enough to zero for you.
        """

        gm = ['Tree.newComp()']

        self._checkModelThing(partNum, symbol, gm[0])
        if self.model.cModel:
            self.deleteCStuff()
        mt = Comp()
        mt.partNum = partNum
        #mt.dim = self.data.parts[partNum].dim
        mt.free = free

        # spec
        if spec not in var.compSpecs:
            gm.append("The spec should be one of %s" % var.compSpecs)
            raise P4Error(gm)
        mt.spec = spec

        mt.num = len(self.model.parts[partNum].comps)
        if symbol:
            mt.symbol = symbol
        else:
            mt.symbol = var.modelSymbols[mt.num]

        self.model.parts[partNum].comps.append(mt)

        # assign val
        dim = self.model.parts[partNum].dim
        if spec == 'equal':
            oneComp = 1.0 / dim
            mt.val = []
            for i in range(dim):
                mt.val.append(oneComp)
        elif spec == 'empirical':
            mt.val = None
        elif spec == 'specified':
            if not val:
                gm.append("Specified comp, but no val.")
                raise P4Error(gm)
            try:
                val = list(val)
            except TypeError:
                gm.append("The 'val' arg should be a list or tuple.")
                raise P4Error(gm)
            if len(val) == dim or len(val) == dim - 1:
                pass
            else:
                gm.append(
                    "Bad length for val arg.  Should be dim or dim-1 long.")
                gm.append("(Dim for this part is %i)" % dim)
                raise P4Error(gm)

            # I allow val's of dim or dim-1.
            if len(val) == dim - 1:
                lastVal = 1.0 - sum(val)
                if lastVal > 0.0:
                    val = val + [1.0 - sum(val)]
                else:
                    gm.append("Bad comp vals %s" % val)
                    gm.append("sum to 1.0 or more.")
                    raise P4Error(gm)
            else:  # len = dim
                theSum = sum(val)
                theDiff = math.fabs(theSum - 1.0)
                # How big to make the delta?  With reasonably good,
                # normalized protein comps (where all the values had just
                # been divided by the total, so it should have summed to
                # 1.0 at that point) I kept getting 1.1e-16.  So make it
                # 5.e-16
                if theDiff > 5.e-16:  # 1e-17 was too small for protein
                    gm.append("Bad comp vals %s" % val)
                    gm.append("does not sum to 1.0")
                    gm.append("The sum = %f" % theSum)
                    gm.append("abs(1.0 - theSum) = %g" % theDiff)
                    raise P4Error(gm)

            # Are any specified values less than PIVEC_MIN?
            needsNormalizing = 0
            for i in range(len(val)):
                thisVal = val[i]
                if thisVal < var.PIVEC_MIN:
                    print gm[0]
                    print "    Specifying a comp of zero for a character is not allowed."
                    print "    Re-setting to %g" % var.PIVEC_MIN
                    val[i] = var.PIVEC_MIN
                    needsNormalizing = 1

            if needsNormalizing:
                theSum = sum(val)
                for i in range(len(val)):
                    val[i] /= theSum
                if math.fabs(sum(val) - 1.0) > 5.e-16:
                    gm.append("Bad comp vals %s" % val)
                    gm.append("does not sum to 1.0")
                    raise P4Error(gm)
            # print "sum(val) - 1.0 = %f (%g)" % (sum(val) - 1.0, sum(val) -
            # 1.0)
            mt.val = val

        # Empirical protein comps are from the dat files in PAML.  Thanks,
        # Ziheng!
        elif spec == 'cpREV':
            mt.val = [0.0755, 0.0621, 0.0410, 0.0371, 0.0091,
                      0.0382, 0.0495, 0.0838, 0.0246, 0.0806,
                      0.1011, 0.0504, 0.0220, 0.0506, 0.0431,
                      0.0622, 0.0543, 0.0181, 0.0307, 0.0660]
            #theSum = sum(mt.val)
            # for i in range(len(mt.val)):
            #    mt.val[i] /= theSum
        elif spec == 'd78':
            # These first values have a couple more decimal places.  I
            # think I got these from some obscure code in a back corner of
            # the NCBI ftp site.  I believe it is obtainable by raising
            # the d78 matrix to a high power.  It is a more precise comp,
            # but is not used here because it is not standard.
            # mt.val = [0.08713, 0.04090, 0.04043, 0.04687, 0.03347,
            #          0.03826, 0.04953, 0.08861, 0.03362, 0.03689,
            #          0.08536, 0.08048, 0.01475, 0.03977, 0.05068,
            #          0.06958, 0.05854, 0.01049, 0.02992, 0.06472]

            # This next set of values is from her paper, and is the set
            # that everybody uses.
            # mt.val = [0.087, 0.041, 0.040, 0.047, 0.033,
            #          0.038, 0.05, 0.089, 0.034, 0.037,
            #          0.085, 0.08, 0.015, 0.04, 0.051,
            #          0.07, 0.058, 0.01, 0.03, 0.065]

            # These values are from Goldman's recommendations (Kosiol &
            # Goldman)
            mt.val = [0.087127, 0.040904, 0.040432, 0.046872, 0.033474,
                      0.038255, 0.049530, 0.088612, 0.033619, 0.036886,
                      0.085357, 0.080481, 0.014753, 0.039772, 0.050680,
                      0.069577, 0.058542, 0.010494, 0.029916, 0.064718]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'jtt':
            # mt.val = [0.077,0.051, 0.043, 0.052, 0.02,
            #          0.041, 0.062, 0.074, 0.023, 0.052,
            #          0.091, 0.059, 0.024, 0.04, 0.051,
            #          0.069, 0.059, 0.014, 0.032, 0.066]

            # again, a Goldman recommendation
            mt.val = [0.076862, 0.051057, 0.042546, 0.051269, 0.020279,
                      0.041061, 0.061820, 0.074714, 0.022983, 0.052569,
                      0.091111, 0.059498, 0.023414, 0.040530, 0.050532,
                      0.068225, 0.058518, 0.014336, 0.032303, 0.066374]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'mtREV24':
            mt.val = [0.072, 0.019, 0.039, 0.019, 0.006,
                      0.025, 0.024, 0.056, 0.028, 0.088,
                      0.168, 0.023, 0.054, 0.061, 0.054,
                      0.072, 0.086, 0.029, 0.033, 0.043]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'mtmam':
            mt.val = [0.0692, 0.0184, 0.0400, 0.0186, 0.0065,
                      0.0238, 0.0236, 0.0557, 0.0277, 0.0905,
                      0.1675, 0.0221, 0.0561, 0.0611, 0.0536,
                      0.0725, 0.0870, 0.0293, 0.0340, 0.0428]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'wag':
            mt.val = [0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078,
                      0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466,
                      0.086209,  0.0620286, 0.0195027, 0.0384319, 0.0457631,
                      0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'rtRev':
            mt.val = [0.0646, 0.0453, 0.0376, 0.0422, 0.0114, 0.0606,
                      0.0607, 0.0639, 0.0273, 0.0679, 0.1018, 0.0751,
                      0.0150, 0.0287, 0.0681, 0.0488, 0.0622, 0.0251,
                      0.0318, 0.0619]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum

        elif spec == 'tmjtt94':
            mt.val = [0.105068479, 0.015695291, 0.018494452, 0.008897331,
                      0.021893432, 0.014095771, 0.009697091, 0.075777267,
                      0.016794962, 0.118764371, 0.163450965, 0.011196641,
                      0.033290013, 0.077676697, 0.025992202, 0.056782965,
                      0.052284315, 0.022293312, 0.032390283, 0.119464161]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'tmlg99':
            mt.val = [0.100632, 0.014017, 0.014706, 0.010371, 0.030668,
                      0.015152, 0.011343, 0.069235, 0.017501, 0.107722,
                      0.155161, 0.009723, 0.038730, 0.086453, 0.031761,
                      0.064333, 0.044847, 0.028277, 0.036988, 0.112380]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'lg':
            mt.val = [0.079066, 0.055941, 0.041977, 0.053052, 0.012937,
                      0.040767, 0.071586, 0.057337, 0.022355, 0.062157,
                      0.099081, 0.064600, 0.022951, 0.042302, 0.044040,
                      0.061197, 0.053287, 0.012066, 0.034155, 0.069147]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'blosum62':
            mt.val = [0.074, 0.052, 0.045, 0.054, 0.025,
                      0.034, 0.054, 0.074, 0.026, 0.068,
                      0.099, 0.058, 0.025, 0.047, 0.039,
                      0.057, 0.051, 0.013, 0.032, 0.073]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'hivb':
            mt.val = [0.060490222, 0.066039665, 0.044127815, 0.042109048, 0.020075899,
                      0.053606488, 0.071567447, 0.072308239, 0.022293943, 0.069730629,
                      0.098851122, 0.056968211, 0.019768318, 0.028809447, 0.046025282,
                      0.05060433, 0.053636813, 0.033011601, 0.028350243, 0.061625237]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'mtart':
            mt.val = [0.054116, 0.018227, 0.039903, 0.020160, 0.009709,
                      0.018781, 0.024289, 0.068183, 0.024518, 0.092638,
                      0.148658, 0.021718, 0.061453, 0.088668, 0.041826,
                      0.091030, 0.049194, 0.029786, 0.039443, 0.057700]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        elif spec == 'mtzoa':
            mt.val = [0.068880,    0.021037,    0.030390,    0.020696,    0.009966,
                      0.018623,    0.024989,    0.071968,    0.026814,    0.085072,
                      0.156717,    0.019276,    0.050652,    0.081712,    0.044803,
                      0.080535,    0.056386,    0.027998,    0.037404,    0.066083]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum
        # Cymon was here!  gcpREV (green chloroplast plant REV seeing how
        # is estimated from green plant chloroplasts alone rather than all
        # chloroplasts).
        elif spec == 'gcpREV':
            mt.val = [0.079510, 0.056001, 0.040459, 0.033220, 0.009051,
                      0.037505, 0.049675, 0.080233, 0.021880, 0.080496,
                      0.107512, 0.049324, 0.020776, 0.047731, 0.039916,
                      0.073820, 0.053615, 0.016705, 0.030790, 0.071781]
            theSum = sum(mt.val)
            for i in range(len(mt.val)):
                mt.val[i] /= theSum

        return mt

    def newRMatrix(self, partNum=0, free=0, spec='ones', val=None, symbol=None):
        """Make, attach, and return a new RMatrix instance.

        spec should be one of:

        -   'ones'         - for JC, poisson, F81
        -   '2p'           - for k2p and hky
        -   'specified'
        -   'cpREV'
        -   'd78'
        -   'jtt'
        -   'mtREV24'
        -   'mtmam'
        -   'wag'
        -   'rtRev'
        -   'tmjtt94'
        -   'tmlg99'
        -   'lg'
        -   'blosum62'
        -   'hivb'
        -   'mtart'
        -   'mtzoa'

        You do not set the 'val' arg unless the spec is 'specified' or
        '2p'.  If spec='2p', then you set val to kappa.

        If the spec is 'specified', you specify all the numerical values
        in a list given as the 'val' arg.  The length of that list will be
        (((dim * dim) - dim) / 2) - 1, so for DNA, where dim=4, you would
        specify a list containing 5 numbers.  """

        # not implemented:
        # 'blosum62a'
        # 'blosum62b'
        # 'phat70'

        complaintHead = '\nTree.newRMatrix()'
        gm = [complaintHead]
        self._checkModelThing(partNum, symbol, complaintHead)
        if self.model.cModel:
            self.deleteCStuff()
        mt = RMatrix()
        mt.partNum = partNum
        #mt.dim = self.data.parts[partNum].dim
        mt.free = free
        if spec not in var.rMatrixSpecs:
            gm.append("Got unknown rMatrix spec '%s'." % spec)
            gm.append("Should be one of: %s" % var.rMatrixSpecs)
            raise P4Error(gm)
        mt.spec = spec
        mt.num = len(self.model.parts[partNum].rMatrices)
        if symbol:
            mt.symbol = symbol
        else:
            mt.symbol = var.modelSymbols[mt.num]

        self.model.parts[partNum].rMatrices.append(mt)

        # assign val
        dim = self.model.parts[partNum].dim
        if var.rMatrixNormalizeTo1:
            goodLen = (((dim * dim) - dim) / 2)
        else:
            goodLen = (((dim * dim) - dim) / 2) - 1

        v = None
        if spec == 'specified':
            if val:
                # should check that values are all floats
                if len(val) == goodLen:
                    v = numpy.array(val, numpy.float)
                    if var.rMatrixNormalizeTo1:
                        v /= v.sum()
                elif var.rMatrixNormalizeTo1 and len(val) == goodLen - 1:
                    gm.append("var.rMatrixNormalizeTo1 is set, val length should be %i, got %i" % (
                        goodLen, len(val)))
                    raise P4Error(gm)
                else:
                    gm.append(
                        "Bad length for arg val.  Length %i, should be %i" % (len(val), goodLen))
                    raise P4Error(gm)
            else:
                gm.append(
                    "spec is 'specified', but there are no specified rMatrix values.")
                gm.append(
                    "Specify rMatrix values by eg val=[2.0, 3.0, 4.0, 5.0,6.0]")
                raise P4Error(gm)
        elif spec == 'ones':
            v = numpy.array([1.0] * goodLen)
            if var.rMatrixNormalizeTo1:
                v /= v.sum()
        elif spec == '2p':
            try:
                v = float(val)
            except (ValueError, TypeError):
                gm.append(
                    "Kappa ('val' arg) should be a float.  Setting to 2.0")
                v = 2.0
            if v < var.KAPPA_MIN:
                gm.append("Kappa is too small.   Setting to %f" %
                          var.KAPPA_MIN)
                v = var.KAPPA_MIN
            elif v > var.KAPPA_MAX:
                gm.append("Kappa is too big.  Setting to %f" % var.KAPPA_MAX)
                v = var.KAPPA_MAX
            v = numpy.array([v], numpy.float)
        elif spec in var.rMatrixProteinSpecs:
            if self.data.parts[partNum].dataType != 'protein':
                gm.append("A protein matrix has been specified, but the dataType for part %i is %s." % (
                    partNum, self.data.parts[partNum].dataType))
                raise P4Error(gm)
            if free:
                gm.append(
                    'The rMatrix should not be free if it is an empirical protein matrix.')
                raise P4Error(gm)

        mt.val = v  # type numpy.ndarray, or None for protein
        return mt

    def newGdasrv(self, partNum=0, free=0, val=None, symbol=None):
        gm = ['Tree.newGdasrv()']

        if not self.model:
            gm.append("Set the data first.  Eg myTree.data = Data()")
            raise P4Error(gm)

        if self.model.cModel:
            self.deleteCStuff()

        # check if there is an nGammaCat > 1:
        if self.model.parts[partNum].nGammaCat == 1:
            gm.append(
                "For this part (%s), the number of nGammaCat has been set to 1." % partNum)
            gm.append("So gdasrv won't work.")
            gm.append(
                "You can set the nGammaCat with yourTree.setNGammaCat(partNum=x, nGammaCat=y)")
            raise P4Error(gm)

        # check val
        if val == None:
            gm.append("Please specify a val, a positive float.")
            raise P4Error(gm)
        try:
            v = float(val)
        except:
            gm.append("Arg val must be a float.  Got '%s'" % val)
            raise P4Error(gm)

        # This week, we have in defines.h
        # define GAMMA_SHAPE_MIN 0.000001
        # define GAMMA_SHAPE_MAX 300.0
        if v <= 0.000001 or v >= 300.0:
            gm.append("Arg val must be between 0.000001 and 300.  Got %f" % v)
            raise P4Error(gm)

        self._checkModelThing(partNum, symbol, gm[0])

        if self.model.parts[partNum].isMixture:
            gm.append("Don't do this if it is a mixture.")
            raise P4Error(gm)

        mt = Gdasrv()
        mt.nGammaCat = self.model.parts[partNum].nGammaCat
        mt.partNum = partNum
        mt.free = free
        # no spec or dim
        mt.num = len(self.model.parts[partNum].gdasrvs)
        if symbol:
            mt.symbol = symbol
        else:
            mt.symbol = var.modelSymbols[mt.num]

        self.model.parts[partNum].gdasrvs.append(mt)
        mt.freqs = numpy.zeros(mt.nGammaCat, numpy.float)
        mt.rates = numpy.zeros(mt.nGammaCat, numpy.float)
        mt._val[0] = v
        mt.calcRates()
        return mt

    def setPInvar(self, partNum=0, free=0, val=0.0):
        complaintHead = '\nTree.setPInvar()'
        gm = [complaintHead]

        # check val
        try:
            v = float(val)
        except:
            gm.append("Arg val must be a float.  Got '%s'" % val)
            raise P4Error(gm)

        if v < 0.0 or v >= 1.0:
            gm.append(
                "Arg val must be zero or more, and less than 1.  Got %f" % v)
            raise P4Error(gm)

        self._checkModelThing(partNum, None, complaintHead)
        if self.model.cModel:
            self.deleteCStuff()
        mt = PInvar()
        mt.partNum = partNum
        mt.free = free
        mt.val = v
        self.model.parts[partNum].pInvar = mt

    def setRelRate(self, partNum=0, val=0.0):
        complaintHead = '\nTree.setRelRate()'
        gm = [complaintHead]

        # check val
        try:
            v = float(val)
        except:
            gm.append("Arg val must be a float.  Got '%s'" % val)
            raise P4Error(gm)

        if v <= 0.0 or v >= 1000.0:
            gm.append(
                "Arg val must be more than zero, and less than 1000 (arbitrarily).  Got %f" % v)
            raise P4Error(gm)

        self._checkModelThing(partNum, None, complaintHead)
        if self.model.cModel:
            self.deleteCStuff()
        self.model.parts[partNum].relRate = v

    def setRjComp(self, partNum=0, val=True):
        if self.model.cModel:
            self.deleteCStuff()
        self.model.parts[partNum].rjComp = val

    def setRjRMatrix(self, partNum=0, val=True):
        if self.model.cModel:
            self.deleteCStuff()
        self.model.parts[partNum].rjRMatrix = val

    def setModelThing(self, theModelThing, node=None, clade=1):
        complaintHead = '\nTree.setModelThing()'
        gm = [complaintHead]

        if self.model.parts[theModelThing.partNum].isMixture:
            gm.append("Don't do this if the part uses a mixture model.")
            raise P4Error(gm)

        if theModelThing and \
            (isinstance(theModelThing, Comp) or
                isinstance(theModelThing, RMatrix) or
                isinstance(theModelThing, Gdasrv)):
            pass
        else:
            gm.append("Expecting a model thing instance of some sort.")
            gm.append("Ie a comp, rMatrix, or gdasrv, instance.")
            gm.append("Got theModelThing = %s" % theModelThing)
            raise P4Error(gm)

        if self.model.cModel:
            self.deleteCStuff()

        partNum = theModelThing.partNum

        if node == None:
            theNode = self.root
        else:
            theNode = self.node(node)

        isBad = 0
        if isinstance(theModelThing, Comp):
            if theModelThing != self.model.parts[partNum].comps[theModelThing.num]:
                isBad = 1
        elif isinstance(theModelThing, RMatrix):
            if theModelThing != self.model.parts[partNum].rMatrices[theModelThing.num]:
                isBad = 1
        elif isinstance(theModelThing, Gdasrv):
            if theModelThing != self.model.parts[partNum].gdasrvs[theModelThing.num]:
                isBad = 1
        else:  # This will never happen-- we checked above.  Overkill.
            gm.append("I don't recognise theModelThing.")
            raise P4Error(gm)
        if isBad:
            gm.append(
                "The modelThing can only be set on the tree that made it.")
            raise P4Error(gm)

        # For the root, we set comps and nothing else.  For other nodes we
        # set anything.
        if theNode == self.root:
            if isinstance(theModelThing, Comp):
                theNode.parts[partNum].compNum = theModelThing.num
        else:
            if isinstance(theModelThing, Comp):
                theNode.parts[partNum].compNum = theModelThing.num
            elif isinstance(theModelThing, RMatrix):
                theNode.br.parts[partNum].rMatrixNum = theModelThing.num
            elif isinstance(theModelThing, Gdasrv):
                theNode.br.parts[partNum].gdasrvNum = theModelThing.num

        if clade:
            aboves = self.getNodeNumsAbove(theNode, leavesOnly=0)
            for i in aboves:
                if isinstance(theModelThing, Comp):
                    self.nodes[i].parts[partNum].compNum = theModelThing.num
                elif isinstance(theModelThing, RMatrix):
                    self.nodes[i].br.parts[
                        partNum].rMatrixNum = theModelThing.num
                elif isinstance(theModelThing, Gdasrv):
                    self.nodes[i].br.parts[
                        partNum].gdasrvNum = theModelThing.num

    def setModelThingsRandomly(self, forceRepresentation=2):
        """Place model things (semi-)randomly on the tree.

        For example, if there are 2 compositions in model part partNum,
        this method will decorate each node of the tree with zeros and
        ones, randomly. The actual thing set is
        node.parts[partNum].compNum.  If the model thing is homogeneous,
        it will just put zeros in all the nodes.

        We want to have each model thing on the tree somewhere, and so it
        is not really randomly set.  If the model thing numbers were
        assigned randomly on the tree, it may occur that some model thing
        numbers by chance would not be represented.  This is not allowed,
        and you can set forceRepresentation to some positive integer, 1 or
        more.  That number will be the lower limit allowed on the number
        of nodes that get assigned the model thing number.  For example,
        if forceRepresentation is set to 2, then each model thing must get
        assigned to at least 2 nodes."""

        gm = ['Tree.setModelThingsRandomly()']

        if not self.model or not self.model.nParts:
            gm.append("No model parts?")
            raise P4Error(gm)

        if self.model.cModel:
            self.deleteCStuff()
        # self.model.dump()

        if type(forceRepresentation) != type(1) or forceRepresentation < 1:
            gm.append("Arg 'forceRepresentation' should be 1 or more.")
            gm.append("Got forceRepresentation = %s" % forceRepresentation)
            raise P4Error(gm)

        for i in self.preOrder:
            if i == var.NO_ORDER:
                gm.append(
                    "This method does not work if any nodes are not used in the tree.")
                raise P4Error(gm)

        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]

            # First do comps
            if mp.nComps == 1:
                for n in self.nodes:
                    n.parts[pNum].compNum = 0
            elif mp.nComps > 1:
                nNodes = len(self.nodes)
                if (mp.nComps * forceRepresentation) > nNodes:
                    gm.append("Part %i" % pNum)
                    gm.append(
                        "There are not enough nodes (%i) to put %i" % (nNodes, mp.nComps))
                    gm.append(
                        "comps on at least forceRepresentation (%i) nodes." % forceRepresentation)
                    raise P4Error(gm)
                nList = self.nodes[:]
                random.shuffle(nList)
                # get the forceRepresentation out of the way first
                for mtNum in range(mp.nComps):
                    for fr in range(forceRepresentation):
                        n = nList.pop()
                        n.parts[pNum].compNum = mtNum
                # Now do the rest
                for n in nList:
                    n.parts[pNum].compNum = random.randrange(mp.nComps)
            else:
                gm.append("No comps in part %i" % pNum)
                raise P4Error(gm)

            # Second do rMatrices
            if mp.nRMatrices == 1:
                for n in self.nodes:
                    if n != self.root:
                        n.br.parts[pNum].rMatrixNum = 0
            elif mp.nRMatrices > 1:
                nNodes = len(self.nodes) - 1
                if (mp.nRMatrices * forceRepresentation) > nNodes:
                    gm.append("Part %i" % pNum)
                    gm.append(
                        "There are not enough nodes (%i) to put %i" % (nNodes, mp.nRMatrices))
                    gm.append(
                        "rMatrices on at least forceRepresentation (%i) nodes." % forceRepresentation)
                    raise P4Error(gm)
                nList = self.nodes[:]
                nList.remove(self.root)
                random.shuffle(nList)
                # get the forceRepresentation out of the way first
                for mtNum in range(mp.nRMatrices):
                    for fr in range(forceRepresentation):
                        n = nList.pop()
                        n.br.parts[pNum].rMatrixNum = mtNum
                # Now do the rest
                for n in nList:
                    n.br.parts[pNum].rMatrixNum = random.randrange(
                        mp.nRMatrices)

            else:
                gm.append("No rMatrices in part %i" % pNum)
                raise P4Error(gm)

            # Third do gdasrvs
            if mp.nGammaCat > 1:
                if mp.nGdasrvs == 1:
                    for n in self.nodes:
                        if n != self.root:
                            n.br.parts[pNum].gdasrvNum = 0
                elif mp.nGdasrvs > 1:
                    nNodes = len(self.nodes) - 1
                    if (mp.nGdasrvs * forceRepresentation) > nNodes:
                        gm.append("Part %i" % pNum)
                        gm.append(
                            "There are not enough nodes (%i) to put %i" % (nNodes, mp.nGdasrvs))
                        gm.append(
                            "gdasrvs on at least forceRepresentation (%i) nodes." % forceRepresentation)
                        raise P4Error(gm)
                    nList = self.nodes[:]
                    nList.remove(self.root)
                    random.shuffle(nList)
                    # get the forceRepresentation out of the way first
                    for mtNum in range(mp.nGdasrvs):
                        for fr in range(forceRepresentation):
                            n = nList.pop()
                            n.br.parts[pNum].gdasrvNum = mtNum
                    # Now do the rest
                    for n in nList:
                        n.br.parts[pNum].gdasrvNum = random.randrange(
                            mp.nGdasrvs)
                else:
                    gm.append(
                        "No gdasrvs in part %i and yet nGammaCat > 1" % pNum)
                    raise P4Error(gm)

        # self.dump(model=True)

    def setModelThingsNNodes(self):
        """Set nNodes for all modelThings"""

        gm = ['Tree.setModelThingsNNodes()']

        if not self.model or not self.model.nParts:
            gm.append("No model parts?")
            raise P4Error(gm)

        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]

            if not mp.nComps:
                gm.append("No comps in model part %i." % pNum)
                raise P4Error(gm)
            elif not mp.nRMatrices:
                gm.append("No rMatrices in model part %i." % pNum)
                raise P4Error(gm)

        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]

            # First do comps
            if mp.nComps == 1:
                pass
            elif mp.nComps > 1:
                for mtNum in range(mp.nComps):
                    mp.comps[mtNum].nNodes = 0
                for n in self.iterNodes():
                    mp.comps[n.parts[pNum].compNum].nNodes += 1

            # Second do rMatrices
            if mp.nRMatrices == 1:
                pass
            elif mp.nRMatrices > 1:
                for mtNum in range(mp.nRMatrices):
                    mp.rMatrices[mtNum].nNodes = 0
                for n in self.iterNodesNoRoot():
                    mp.rMatrices[n.br.parts[pNum].rMatrixNum].nNodes += 1

            # Third do gdasrvs
            if mp.nGammaCat > 1:
                if mp.nGdasrvs == 1:
                    pass
                elif mp.nGdasrvs > 1:
                    for mtNum in range(mp.nGdasrvs):
                        mp.gdasrvs[mtNum].nNodes = 0
                    for n in self.iterNodesNoRoot():
                        mp.gdasrvs[n.br.parts[pNum].gdasrvNum].nNodes += 1
                else:
                    gm.append("No gdasrvs in part %i" % pNum)
                    raise P4Error(gm)

    def summarizeModelThingsNNodes(self):
        """Summarize nNodes for all modelThings if isHet"""

        gm = ['Tree.summarizeModelThingsNNodes()']

        if not self.model or not self.model.nParts:
            gm.append("No model parts?")
            raise P4Error(gm)
        if not self.model.isHet:
            gm.append("This method is for hetero models")
            raise P4Error(gm)

        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]

            if not mp.nComps:
                gm.append("No comps in model part %i." % pNum)
                raise P4Error(gm)
            elif not mp.nRMatrices:
                gm.append("No rMatrices in model part %i." % pNum)
                raise P4Error(gm)

        for pNum in range(self.model.nParts):
            print "\n%6s %s:" % ("Part", pNum)
            mp = self.model.parts[pNum]

            # First do comps
            if mp.nComps == 1:
                pass
            elif mp.nComps > 1:
                for mtNum in range(mp.nComps):
                    # print "  comp %i nNodes=%i" % (mtNum,
                    # mp.comps[mtNum].nNodes)
                    print "%16s %i %s = %i" % ("composition", mtNum, "nNodes",
                                               mp.comps[mtNum].nNodes)

            # Second do rMatrices
            if mp.nRMatrices == 1:
                pass
            elif mp.nRMatrices > 1:
                for mtNum in range(mp.nRMatrices):
                    print "%16s %i %s = %i" % ("rate matrix", mtNum,
                                               "nNodes", mp.rMatrices[mtNum].nNodes)

            # Third do gdasrvs
            if mp.nGammaCat > 1:
                if mp.nGdasrvs == 1:
                    pass
                elif mp.nGdasrvs > 1:
                    for mtNum in range(mp.nGdasrvs):
                        print "  gdasrv %i nNodes =%i" % (mtNum, mp.gdasrvs[mtNum].nNodes)
                else:
                    gm.append("No gdasrvs in part %i" % pNum)
                    raise P4Error(gm)

    def setTextDrawSymbol(self, theSymbol='-', node=None, clade=1):
        gm = ['\nTree.setTextDrawString()']

        if not theSymbol or type(theSymbol) != type('c') or len(theSymbol) != 1:
            gm.append("theSymbol should be a single character string.")
            raise P4Error(gm)

        if not node:
            theNode = self.root
        else:
            theNode = self.node(node)

        if theNode == self.root:
            pass
        else:
            theNode.br.textDrawSymbol = theSymbol

        if clade:
            aboves = self.getNodeNumsAbove(theNode, leavesOnly=0)
            for i in aboves:
                self.nodes[i].br.textDrawSymbol = theSymbol

    def setNGammaCat(self, partNum=0, nGammaCat=1):
        gm = ['\nTree.setNGammaCat()']
        if not self.data or not self.model:
            gm.append("No data?")
            raise P4Error(gm)
        if self.model.cModel:
            self.deleteCStuff()
        if partNum < 0 or partNum >= self.model.nParts:
            gm.append("PartNum %s is out of range of %s parts." %
                      (partNum, self.model.nParts))
            raise P4Error(gm)
        if self.model.parts[partNum].isMixture:
            gm.append("Don't do this if the part uses a mixture model.")
            raise P4Error(gm)

        try:
            x = int(nGammaCat)
        except ValueError:
            gm.append("'%s' does not appear to be an integer." % i)
            raise P4Error(gm)
        if x < 1:
            gm.append("nGammaCat should not be less than 1.")
            raise P4Error(gm)
        elif x > 16:
            gm.append("nGammaCat '%s' exceeds the arbitrary limit of 16." % x)
            raise P4Error(gm)
        self.model.parts[partNum].nGammaCat = nGammaCat

    # def setMixture(self, partNum=0, free=0, freqs=None, rates=None):
    #     complaintHead = '\nTree.setMixture()'
    #     gm = [complaintHead]

    #     if 1:
    #         gm.append("Sorry -- not turned on yet.")
    #         raise P4Error(gm)

    #     assert self.model and self.data
    #     if self.model.cModel:
    #         self.deleteCStuff()
    #     if partNum < 0 or partNum >= self.model.nParts:
    #         gm.append("PartNum %s is out of range of %s parts." % (partNum, self.model.nParts))
    #         raise P4Error(gm)
    #     mp = self.model.parts[partNum]

    #     isHet = 0
    #     for n in self.nodes:
    #         if n.parts[partNum].compNum >= 0:
    #             isHet = 1
    #             break
    #         if n != self.root:
    #             if n.br.parts[partNum].rMatrixNum >= 0 or n.br.parts[partNum].gdasrvNum >= 0:
    #                 isHet = 1
    #                 break
    #     if isHet:
    #         gm.append("Can't be both hetero on the tree and a mixture model.")
    #         gm.append("There seems to have been a previous setModelThing()")
    #         raise P4Error(gm)
    #     mp.isMixture = 1
    #     mp.mixture = Mixture()
    #     mp.mixture.free = free

    #     #if freqs and rates:
    #     nMix = mp.nComps * mp.nRMatrices
    #     mp.mixture.freqs = numpy.zeros(nMix, numpy.float)
    #     mp.mixture.rates = numpy.zeros(nMix, numpy.float)

    #     try:
    #         for i in range(nMix):
    #             mp.mixture.freqs[i] = float(freqs[i])
    #             mp.mixture.rates[i] = float(rates[i])
    #     except:
    #         gm.append('Args freqs and rates should be sequences of floats, each nComps * nMatrices long.')
    #         raise P4Error(gm)
    #     if len(freqs) != nMix or len(rates) != nMix:
    #         gm.append('Args freqs and rates should be sequences of floats, each nComps * nMatrices long.')
    #         raise P4Error(gm)

    def modelSanityCheck(self, resetEmpiricalComps=True):
        """Check that the tree, data, and model specs are good to go.

        Complain and exit if there is anything wrong that might prevent a
        likelihood evaluation from being done.  We are assuming that a
        data object exists and is attached, and that model stuff has been
        set.

        Check that each part has at least 1 each from comps, rMatrices,
        and gdasrvs (if nGammaCat is > 1).

        If it is not a mixture model for a particular part, check that
        each node has a comp, rMatrix, and gdasr.  Check that all comps,
        rMatrices, gdasrvs are used on a node somewhere.

        Here relRate, ie the relative rate of each data partition, is
        adjusted based on the size of the data partitions.

        newRelRate_p = oldRelRate_p * (Sum_p[oldRelRate_i * partLen_i] / Sum[partLen_i])

        That ensures that Sum(relRate_i * partLen_i) = totalDataLength, ie
        that the weighted mean of the rates is 1.0.

        This method also tallies up the number of free prams in the whole
        model, and sets self.model.nFreePrams.

        """

        complaintHead = '\nTree.modelSanityCheck()'
        gm = [complaintHead]
        # print "\nTree.modelSanityCheck() here. self.model.nParts=%s" % self.model.nParts
        # print "\nTree.modelSanityCheck() here.  resetEmpiricalComps=%s" %
        # resetEmpiricalComps
        isBad = 0
        complaints = []
        if not self.data:
            complaints.append('    No data.')
            isBad = 1
        if not self.model:
            complaints.append('    No model.')
            isBad = 1

        # Set isHet.
        for pNum in range(self.model.nParts):
            mp = self.model.parts[pNum]
            mp.isHet = 0
            if mp.nComps > 1 or mp.nRMatrices > 1:
                mp.isHet = 1
            if mp.nGammaCat > 1 and mp.nGdasrvs > 1:
                mp.isHet = 1

        # Check that all parts have all the required stuff.  Make a list
        # of errors.  If there is something missing or wrong, don't die
        # right away, but add the problem to the list, and write it all
        # out at the end.  It gives the user a chance to fix more than one
        # error at a time.
        for pNum in range(self.model.nParts):
            complaints.append('  Part %i' % pNum)
            partIsBad = 0
            mp = self.model.parts[pNum]

            # Check if essential things have been set
            if not mp.nComps:
                complaints.append('    No comps in part %s' % pNum)
                partIsBad = 1
            if not mp.nRMatrices:
                complaints.append('    No rMatrices in part %s' % pNum)
                partIsBad = 1
            if mp.nGammaCat > 1:
                if not mp.nGdasrvs:
                    complaints.append('    No gdasrvs in part %s' % pNum)
                    partIsBad = 1
            if mp.nGammaCat == 1:
                if mp.nGdasrvs:
                    complaints.append(
                        '    There should be no gdasrvs in part %s, with nGammaCat=1' % pNum)
                    partIsBad = 1
            if not mp.pInvar:
                complaints.append('    No pInvar in part %s' % pNum)
                partIsBad = 1

            if partIsBad:
                gm.append("  (Indices are zero-based.)")
                gm += complaints
                raise P4Error(gm)

            # Check if comp values have been set.
            for mt in mp.comps:
                if mt.spec != 'empirical' or not resetEmpiricalComps:
                    if not mt.val:
                        complaints.append(
                            '    No composition val in part %s' % pNum)
                        partIsBad = 1
                    if len(mt.val) != mp.dim:
                        complaints.append('    Composition val is wrong length (%i), but dim is %i' % (
                            len(mt.val), mp.dim))
                        partIsBad = 1

            # We don't want multiple rMatrices or free rMatrices if mp.dim is 2
            if mp.dim == 2:
                if mp.nRMatrices > 1:
                    complaints.append(
                        '    Part %s is dim 2, but we have more than one rMatrix' % pNum)
                    partIsBad = 1
                mt = mp.rMatrices[0]  # hopefully only one
                if mt.free:
                    complaints.append(
                        '    Part %s is dim 2, but rMatrix 0 is free' % pNum)
                    partIsBad = 1

            # If isMixture, then it may not have nGdasrvs
            if mp.isMixture:
                if mp.nGdasrvs:
                    complaints.append(
                        '    If it isMixture, then gdasrv may not be on.')
                    partIsBad = 1

            # If isMixture, then it cannot be isHet
            if mp.isMixture:
                for n in self.nodes:
                    n.parts[pNum].compNum = -1
                    if n != self.root:
                        n.br.parts[pNum].rMatrixNum = -1
                        if mp.nGammaCat > 1:
                            n.br.parts[pNum].gdasrvNum = -1

            if mp.isMixture:
                mt = mp.mixture
                if not mt.freqs or not mt.rates:
                    complaints.append(
                        '    This week, you must specify mixture freqs and rates.')
                    partIsBad = 1
                print 'mt.freqs = %s' % mt.freqs
                print 'mt.rates = %s' % mt.rates
                if len(mt.freqs) != len(mt.rates):
                    complaints.append(
                        '    Lengths of mixture freqs and rates differ.')
                    partIsBad = 1
                mp.nCat = mp.nComps * mp.nRMatrices
                # if nCompsTimesNRMatrices == 1:
                #    complaints.append('    nComps * nRMatrices = 1, no point in having a mixture.')
                #    partIsBad = 1
                if len(mt.freqs) != mp.nCat:
                    complaints.append(
                        '    Wrong length of mixture freqs and rates.  Should be %i' % mp.nCat)
                    partIsBad = 1

                # print "Freqs = %s" % mt.freqs
                # print "Rates = %s" % mt.rates
                theSum = sum(mt.freqs)
                if theSum != 1.0:
                    for i in range(len(mt.freqs)):
                        mt.freqs[i] /= theSum
                theSum = 0.0
                for i in range(len(mt.freqs)):
                    theSum += mt.freqs[i] * mt.rates[i]
                # print "Mixture mean = %f (un-normalized)" % theSum
                if theSum != 1.0:
                    for i in range(len(mt.freqs)):
                        mt.rates[i] /= theSum
                if 1:
                    theSum = 0
                    for i in range(len(mt.freqs)):
                        theSum += mt.freqs[i] * mt.rates[i]
                    if theSum < 1.0 - 1.0e-9 or theSum > 1.0 + 1.0e-9:
                        gm.append(
                            "Failed to normalize mixture rates.  Sum = %19.17f" % theSum)
                        raise P4Error(gm)
                    # else:
                    #    print "...successfully normalized mixture rates."
                    # print "Freqs = %s" % mt.freqs
                    # print "Rates = %s" % mt.rates

            else:  # not isMixture
                mp.nCat = mp.nGammaCat
                # If the model part isHet, we need to check that all nodes
                # have something assigned, and that all model things are
                # used.  If the model part is not het, we can skip that,
                # but we need to check that all the
                # node.parts[pNum].compNum are 0, and all the
                # node.br.parts[pNum].rMatrixNum and
                # node.br.parts[pNum].gdasrvNum are set to 0.
                if not mp.isHet:
                    # print "model part %i is not het" % pNum
                    for n in self.nodes:
                        # print "pNum = %i, n.nodeNum=%i, len n.parts = %i" %
                        # (pNum, n.nodeNum, len(n.parts))
                        n.parts[pNum].compNum = 0
                        if n != self.root:
                            n.br.parts[pNum].rMatrixNum = 0
                            if mp.nGammaCat > 1:
                                n.br.parts[pNum].gdasrvNum = 0
                else:  # isHet

                    # If there is only one comp, rMatrix, or gdasrv, then
                    # simply set it.
                    if mp.nComps == 1:
                        for n in self.nodes:
                            n.parts[pNum].compNum = 0
                    if mp.nRMatrices == 1:
                        for n in self.nodes:
                            if n != self.root:
                                n.br.parts[pNum].rMatrixNum = 0
                    if mp.nGammaCat > 1 and mp.nGdasrvs == 1:
                        for n in self.nodes:
                            if n != self.root:
                                n.br.parts[pNum].gdasrvNum = 0

                    # print "model part %i is het" % pNum
                    # New ad hoc attribute 'isUsed', to keep track of whether
                    # any node uses it.
                    for mt in mp.comps:
                        mt.isUsed = 0
                    for mt in mp.rMatrices:
                        mt.isUsed = 0
                    for mt in mp.gdasrvs:
                        mt.isUsed = 0

                    # Does every node have all required things?
                    for n in self.nodes:
                        mtNum = n.parts[pNum].compNum
                        if mtNum >= 0 and mtNum < mp.nComps:
                            mt = mp.comps[mtNum]
                            mt.isUsed = 1
                        else:
                            complaints.append(
                                '    Part %s, node %s has no comp.' % (pNum, n.nodeNum))
                            partIsBad = 1

                        if n != self.root:
                            mtNum = n.br.parts[pNum].rMatrixNum
                            if mtNum >= 0 and mtNum < mp.nRMatrices:
                                mt = mp.rMatrices[n.br.parts[pNum].rMatrixNum]
                                mt.isUsed = 1
                            else:
                                complaints.append(
                                    '    Part %s, node %s has no rMatrix.' % (pNum, n.nodeNum))
                                partIsBad = 1
                            if mp.nGammaCat > 1:
                                mtNum = n.br.parts[pNum].gdasrvNum
                                if mtNum >= 0 and mtNum < mp.nGdasrvs:
                                    mt = mp.gdasrvs[n.br.parts[pNum].gdasrvNum]
                                    mt.isUsed = 1
                                else:
                                    complaints.append('    Part %s, node %s has no gdasrv. nGammaCat=%s' % (
                                        pNum, n.nodeNum, mp.nGammaCat))
                                    partIsBad = 1
                            if mp.nGammaCat == 1:
                                if n.br.parts[pNum].gdasrvNum != -1:
                                    complaints.append('    Part %s, node %s has a gdasrv, but nGammaCat is 1.' % (
                                        pNum, n.nodeNum))
                                    partIsBad = 1

                    # Is every model thing used?
                    if not mp.rjComp:
                        for mt in mp.comps:
                            if not mt.isUsed:
                                complaints.append(
                                    '    Part %s, comp %s is not used.' % (pNum, mt.num))
                                partIsBad = 1
                    if not mp.rjRMatrix:
                        for mt in mp.rMatrices:
                            if not mt.isUsed:
                                complaints.append(
                                    '    Part %s, rMatrix %s is not used.' % (pNum, mt.num))
                                partIsBad = 1
                    for mt in mp.gdasrvs:
                        if not mt.isUsed:
                            complaints.append(
                                '    Part %s, gdasrv %s is not used.' % (pNum, mt.num))
                            partIsBad = 1

                    # Clean up ad hoc attr 'isUsed'
                    for mt in mp.comps:
                        del(mt.isUsed)
                    for mt in mp.rMatrices:
                        del(mt.isUsed)
                    for mt in mp.gdasrvs:
                        del(mt.isUsed)

            if partIsBad:
                isBad = 1
            else:
                complaints.append('    ok')

        # ##################################
        if resetEmpiricalComps:
            self.setEmpiricalComps()

        # self.model.isHet if any part isHet
        self.model.isHet = 0
        for pNum in range(self.model.nParts):
            if self.model.parts[pNum].isHet:
                self.model.isHet = 1
                break

        # relativeRates
        self.model.doRelRates = 0
        if self.model.nParts > 1:
            for p in self.model.parts:
                if p.relRate != 1.0:  # This week, the default relRate is 1.0
                    self.model.doRelRates = 1
                    break
        if self.model.relRatesAreFree:
            self.model.doRelRates = 1

        if self.model.doRelRates:
            totDataLen = 0
            for p in self.data.parts:
                totDataLen += p.nChar
            fact = 0.0
            for i in range(self.model.nParts):
                fact += (self.model.parts[i].relRate *
                         self.data.parts[i].nChar)
            fact = float(totDataLen) / fact
            for p in self.model.parts:
                p.relRate *= fact
            if 0:
                print "RelativeRates (adjusted for length)"
                for i in range(self.model.nParts):
                    p = self.model.parts[i]
                    print "  part %s,  nChar %5s, relRate %s" % (p.num, self.data.parts[i].nChar, p.relRate)
            if 1:
                total = 0.0
                for i in range(self.model.nParts):
                    total += (self.model.parts[i].relRate *
                              (float(self.data.parts[i].nChar) / float(totDataLen)))
                if abs(total - 1.0) > 1.0e-12:
                    gm.append(
                        'Error in relativeRate calculation (total=%s).' % total)
                    raise P4Error(gm)

        # print "modelSanityCheck. relRatesAreFree=%s, doRelRates=%s" %
        # (self.model.relRatesAreFree, self.model.doRelRates)

        # tSCovarion
        for p in self.model.parts:
            if p.tSCovarion:
                if p.nComps > 1 or p.nRMatrices > 1 or p.nGammaCat > 1:
                    gm.append(
                        "When tSCovarion is on, there should be no heterogeneity in comps, rMatrices, and gdasrv.")
                    raise P4Error(gm)
                if p.pInvar.val > 0.0 or p.pInvar.free:
                    gm.append(
                        "When tSCovarion is on, you can't use pInvar.  Turn it off.")
                    raise P4Error(gm)

        # model.nFreePrams
        self.model.nFreePrams = 0
        for mp in self.model.parts:
            for mt in mp.comps:
                if mt.free:
                    self.model.nFreePrams += mp.dim - 1
            for mt in mp.rMatrices:
                if mt.free:
                    if mt.spec == '2p':
                        self.model.nFreePrams += 1
                    else:
                        self.model.nFreePrams += (
                            ((mp.dim * mp.dim) - mp.dim) / 2) - 1
            for mt in mp.gdasrvs:
                if mt.free:
                    self.model.nFreePrams += 1
            if mp.pInvar.free:
                self.model.nFreePrams += 1
            if mp.doTSCovarion and mp.tSCovarion.free:
                self.model.nFreePrams += 2
            if mp.isMixture and mp.mixture.free:
                self.model.nFreePrams += 2 * (len(mp.mixture.freqs) - 1)
        # print "Tree.modelSanityCheck().  Counted %i free params." %
        # self.model.nFreePrams

        if self.model.doRelRates and self.model.relRatesAreFree:
            self.model.nFreePrams += self.model.nParts - 1

        if isBad:
            gm.append("(Indices are zero-based.)")
            gm += complaints
            raise P4Error(gm)

    def setEmpiricalComps(self):
        """Set any empirical model comps to the comp of the data.

        This is done by self.modelSanityCheck(), but sometimes you may
        want to do it at other times.  For example, do this after
        exchanging Data objects, or after simulating.  In those cases
        there does not seem to be a reasonable way to do it
        automatically."""
        complaintHead = '\nTree.setEmpiricalComps()'
        gm = [complaintHead]
        if not self.model:
            gm.append("This tree has no model.")
            raise P4Error(gm)
        if not self.data:
            gm.append("This tree has no data.")
            raise P4Error(gm)

        for mp in self.model.parts:
            for c in mp.comps:
                if c.spec == 'empirical':
                    # print "got empirical comp, comp %s in part %s. (nComps=%i, isHet=%s)" % (
                    #    c.num, mp.num, mp.nComps, mp.isHet)
                    if not mp.isHet:
                        seqNums = None
                    elif mp.nComps == 1:
                        seqNums = None
                    else:
                        seqNums = []
                        # for n in self.nodes:
                        #    print "node %2i seqNum=%3i n.parts[%i].compNum=%3i" % (
                        # n.nodeNum, n.seqNum, mp.num, n.parts[mp.num].compNum)

                        for n in self.nodes:
                            # Is the comp used by the node?
                            if n.parts[mp.num].compNum == c.num:
                                # print "comp %s is used by node %s" % (c.num,
                                # n.nodeNum)
                                if n.isLeaf:
                                    nodeNums = [n.nodeNum]
                                else:
                                    nodeNums = self.getNodeNumsAbove(
                                        n, leavesOnly=1)
                                # gm.append("nodeNums for %s = %s" %
                                # (n.nodeNum, nodeNums)
                                for i in nodeNums:
                                    seqNum = self.nodes[i].seqNum
                                    if seqNum not in seqNums:
                                        seqNums.append(seqNum)

                        # print "setEmpiricalComps() got seqNums = %s" %
                        # seqNums
                        if not seqNums:
                            gm.append(
                                "Something is wrong here.  part %i, comp %i." % (mp.num, c.num))
                            gm.append(
                                "This comp object has no sequences from which to get the empirical comp.")
                            gm.append(
                                "Maybe you need to yourTree.setModelThing() or ")
                            gm.append("yourTree.setModelThingsRandomly()")
                            gm.append(
                                "Or maybe its an extra comp in an RJ MCMC? -- If so, fix")
                            gm.append("the comp val to eg 'equal'.")
                            raise P4Error(gm)

                    c.val = self.data.parts[mp.num].composition(
                        seqNums)  # dim long, not dim - 1
                    # print "  seqNums=%s, c.val=%s" % (seqNums, c.val)

                    needsNormalizing = 0
                    theSum = 0.0
                    for i in range(len(c.val)):
                        if c.val[i] < var.PIVEC_MIN:
                            c.val[
                                i] = var.PIVEC_MIN + (0.2 * var.PIVEC_MIN) + (var.PIVEC_MIN * random.random())
                            needsNormalizing = 1
                        theSum += c.val[i]
                    # print "setEmpiricalComps().  Got theSum = %i" % theSum

                    # We may have asked for the comp of an empty sequence,
                    # in which case val is all zeros.  Check for that.
                    if abs(1.0 - theSum) > 0.1:
                        gm.append(
                            "Something is very wrong here.  Empirical comp vals should sum to 1.0")
                        gm.append("The sum of the comp vals for part %s, comp %s, is %s" % (
                            mp.num, c.num, theSum))
                        gm.append(
                            "Probably the sequences from which the composition was taken were blank.")
                        raise P4Error(gm)

                    if needsNormalizing or abs(theSum - 1.0) > 1e-16:
                        for i in range(len(c.val)):
                            c.val[i] /= theSum

    # def __del__(self, freeTree=pf.p4_freeTree, freeNode=pf.p4_freeNode):
    # def __del__(self, freeTree=pf.p4_freeTree, dp_freeTree = pf.dp_freeTree, mysys=sys):
    # def __del__(self, freeTree=pf.p4_freeTree, dp_freeTree = pf.dp_freeTree):
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
        for i in range(len(self.nodes)):
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
            for n in self.nodes:
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
            for n in self.nodes:
                if n.isLeaf:
                    nLeaves += 1
            self.partLikes = numpy.zeros(self.model.nParts, numpy.float)
            self.cTree = pf.p4_newTree(len(self.nodes), nLeaves, self.preOrder,
                                       self.postOrder, self.partLikes, self.data.cData, self.model.cModel)
            if not self.cTree:
                gm.append("Unable to allocate a cTree")
                raise P4Error(gm)
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
                n.cNode = pf.p4_newNode(
                    n.nodeNum, self.cTree, n.seqNum, n.isLeaf, inTree)
                if not n.cNode:
                    gm.append("Unable to allocate a cNode")
                    raise P4Error(gm)

        # print "finished Tree._allocCStuff()"

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
        for n in self.nodes:
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
                    for n in self.nodes:
                        pf.p4_setCompNum(n.cNode, pNum, n.parts[pNum].compNum)
                        if n != self.root:
                            pf.p4_setRMatrixNum(
                                n.cNode, pNum, n.br.parts[pNum].rMatrixNum)
                            pf.p4_setGdasrvNum(
                                n.cNode, pNum, n.br.parts[pNum].gdasrvNum)

        # pre- and postOrder
        if not self.preAndPostOrderAreValid:
            self.setPreAndPostOrder()

        # for i in range(len(self.nodes)):
        # pf.p4_setPreAndPostOrder(self.cTree, i, self.preOrder[i],
        # self.postOrder[i]) # no longer needed

        # print "finished Tree.setCStuff()"

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

        # print "self.cTree = %s" % self.cTree
        if not self.cTree:
            # This calls self.modelSanityCheck(), which calls
            # self.setEmpiricalComps()
            self._allocCStuff(resetEmpiricalComps=resetEmpiricalComps)
        # print "About to self.model.setCStuff()"
        self.model.setCStuff()
        # print "About to self.setCStuff()"
        self.setCStuff()
        # print "about to p4_setPrams()..."
        pf.p4_setPrams(self.cTree, -1)  # "-1" means do all parts

    def calcLogLike(self, verbose=1, resetEmpiricalComps=True):
        """Calculate the likelihood of the tree, without optimization."""

        self._commonCStuff(resetEmpiricalComps=resetEmpiricalComps)
        # print "about to p4_treeLogLike()..."
        # second arg is getSiteLikes
        self.logLike = pf.p4_treeLogLike(self.cTree, 0)
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
            raise P4Error(gm)

        # Do the opt.
        if allBrentPowell:
            pf.p4_allBrentPowellOptimize(self.cTree)
        elif simplex:
            pf.p4_simplexOptimize(self.cTree, self, Tree.simplexDump)
        else:
            pf.p4_newtSetup(self.cTree)
            pf.p4_newtAndBrentPowellOpt(self.cTree)

        # second arg is getSiteLikes
        self.logLike = pf.p4_treeLogLike(self.cTree, 0)

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
        # for a in self.data.alignments:
        #
        #    if a.nexusSets:
        #        a.nexusSets = None

        # Probably better to do something like this
        # a.nexusSets.constant.mask = self.constantMask()
        # at the end.

        # print "About to pf.p4_simulate(self.cTree)"
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
                gm.append(
                    "resetSequences is not set, but resetNexusSetsConstantMask is set,")
                gm.append("which is probably not going to work as you want.")
                raise P4Error(gm)

    def getSiteLikes(self):
        """Likelihoods, not log likes. Placed in self.siteLikes, a list."""
        self._commonCStuff()
        # second arg is getSiteLikes
        self.logLike = pf.p4_treeLogLike(self.cTree, 1)
        self.siteLikes = []
        for p in self.data.parts:
            self.siteLikes += pf.getSiteLikes(p.cPart)

    # def getWinningGammaCats(self):
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
        # second arg is getSiteLikes
        self.logLike = pf.p4_treeLogLike(self.cTree, 0)
        #self.winningGammaCats = []
        # for p in self.data.parts:
        #    self.winningGammaCats += pf.getWinningGammaCats(p.cPart)
        results = []

        for partNum in range(len(self.data.parts)):
            if len(self.model.parts[partNum].gdasrvs) > 1:
                gm = ['Tree.getSiteRates()']
                gm.append("Part %i has %i gdasrvs.  Maximum 1 allowed." % (
                    partNum, len(self.model.parts[partNum].gdasrvs)))
                raise P4Error(gm)

        for partNum in range(len(self.data.parts)):
            p = self.data.parts[partNum]
            if self.model.parts[partNum].nGammaCat == 1:
                siteRates = numpy.ones(p.nChar, numpy.float)
                gammaCats = numpy.zeros(p.nChar, numpy.int32)
            elif self.model.parts[partNum].nGammaCat > 1:
                siteRates = numpy.zeros(p.nChar, numpy.float)
                gammaCats = numpy.zeros(p.nChar, numpy.int32)
                work = numpy.zeros(
                    self.model.parts[partNum].nGammaCat, numpy.float)
                for charNum in range(p.nChar):
                    gammaCats[charNum] = -1
                #pf.getWinningGammaCats(self.cTree, p.cPart, i, gammaCats, work)
                pf.getSiteRates(
                    self.cTree, p.cPart, partNum, siteRates, gammaCats, work)
                # print siteRates
                # print gammaCats
                # print work
                if 0:
                    counts = numpy.zeros(
                        self.model.parts[partNum].nGammaCat, numpy.int32)
                    for charNum in range(p.nChar):
                        counts[winningGammaCats[charNum]] += 1
                    print counts

            else:
                raise P4Error("This should not happen.")
            results.append([siteRates, gammaCats])
        return results


    def simsForModelFitTests(self, reps=10, seed=None):
        """Do simulations for model fit tests.

        The model fit tests are the Goldman-Cox test, and the tree- and
        model-based composition fit test.  Both of those tests require
        simulations, optimization of the tree and model parameters on the
        simulated data, and extraction of statistics for use in the null
        distribution.  So might as well do them together.  The Goldman-Cox
        test is not possible if there are any gaps or ambiguities, and in
        that case Goldman-Cox simulation stats are not collected.

        Doing the simulations is therefore the time-consuming part, and so
        this method facilitates doing that job in sections.  If you do
        that, set the random number seed to different numbers.  If the
        seed is not set, the process id is used.  (So obviously you should
        explicitly set the seed if you are doing several runs in the same
        process.)  Perhaps you may want to do the simulations on different
        machines in a cluster.  The stats are saved to files.  The output
        files have the seed number attached to the end, so that different
        runs of this method will have different output file names.
        Hopefully.

        When your model uses empirical comps, simulation uses the
        empirical comp of the original data for simulation (good), then
        the optimization part uses the empirical comp of the
        newly-simulated data (also good, I think).  In that case, if it is
        tree-homogeneous, the X^2_m statistic would be identical to the
        X^2 statistic.

        You would follow this method with the modelFitTests() method,
        which uses all the stats files to make null distributions to
        assess significance of the same stats from self."""

        #gm = ['Tree.simsForModelFitTests()']

        # Make a new data object in which to do the sims, so we do not over-write self
        # print "a self.data = %s" % self.data
        # self.data.dump()
        savedData = self.data
        self.data = None  # This triggers self.deleteCStuff()
        self.data = savedData.dupe()

        # We need 2 trees, one for sims, and one for evaluations.  We can
        # use self for sims.  Make a copy for evaluations.
        evalTree = self.dupe()
        evalTree.data = self.data

        # make sure all the memory works ...
        self.calcLogLike(verbose=0)
        evalTree.calcLogLike(verbose=0)

        # We can't do the Goldman-Cox test if there are any gaps or
        # ambiguities.
        doGoldmanCox = True
        for a in self.data.alignments:
            if a.hasGapsOrAmbiguities():
                doGoldmanCox = False
                break
        # print "sims doGoldmanCox = %s" % doGoldmanCox

        # Collect info about the observed data
        statsHashList = []  # one for each data part
        for pNum in range(self.data.nParts):
            h = {}
            statsHashList.append(h)
            h['individualNSites'] = []
            h['observedIndividualCounts'] = []
            for j in range(self.data.nTax):
                h['individualNSites'].append(
                    pf.partSequenceSitesCount(self.data.parts[pNum].cPart, j))  # no gaps or qmarks
                h['observedIndividualCounts'].append(
                    self.data.parts[pNum].composition([j]))
                # (In the line above, its temporarily composition, not counts)
                # print "got seq %i comp = %s' % (j, h['observedIndividualCounts"][-1])
            # At the moment, h['observedIndividualCounts'] has composition,
            # not counts.  So multiply by h['individualNSites']
            for i in range(self.data.nTax):
                for j in range(self.data.parts[pNum].dim):
                    h['observedIndividualCounts'][i][
                        j] *= h['individualNSites'][i]

        # We will want to skip any sequences composed of all gaps
        skipTaxNums = []
        for pNum in range(self.data.nParts):
            stn = []
            for tNum in range(self.data.nTax):
                if not statsHashList[pNum]['individualNSites'][tNum]:
                    stn.append(tNum)
            skipTaxNums.append(stn)
        # print "skipTaxNums = %s" % skipTaxNums

        if seed == None:
            seed = os.getpid()
        pf.reseedCRandomizer(int(seed))

        # Open up some output files in which to put the sim data
        outfileBaseName = 'sims'  # Could be an argument, user-assignable.
        if doGoldmanCox:
            f2Name = outfileBaseName + '_GoldmanStats_%s' % seed
            f2 = open(f2Name, 'w')
            f2.write('# part\tunconstr L\t log like \tGoldman-Cox stat\n')

        f3Name = outfileBaseName + '_CompStats_%s' % seed
        f3 = open(f3Name, 'w')

        # When sims are done when the comp is empirical (whether or not
        # free) we need to re-set the comps based on the newly-simulated
        # data.  So first find out if any comps are empirical.
        hasEmpiricalComps = 0
        for mp in self.model.parts:
            for c in mp.comps:
                if c.spec == 'empirical':
                    hasEmpiricalComps = 1
                    break
        # print "hasEmpiricalComps=%s" % hasEmpiricalComps

        # Do the sims
        for i in range(reps):
            self.simulate()
            if hasEmpiricalComps:
                # Set empirical comps based on newly-simulated data
                evalTree.setEmpiricalComps()
            evalTree.optLogLike(verbose=0)     # The time-consuming part

            if doGoldmanCox:
                if self.data.nParts > 1:
                    self.data.calcUnconstrainedLogLikelihood2()
                    diff = self.data.unconstrainedLogLikelihood - \
                        evalTree.logLike
                    f2.write(
                        '-1\t%f\t%f\t%f\n' % (self.data.unconstrainedLogLikelihood, evalTree.logLike, diff))
                for pNum in range(self.data.nParts):
                    unc = pf.getUnconstrainedLogLike(
                        self.data.parts[pNum].cPart)
                    # 0 for getSiteLikes
                    like = pf.p4_partLogLike(
                        evalTree.cTree, self.data.parts[pNum].cPart, pNum, 0)
                    diff = unc - like
                    f2.write('%i\t%f\t%f\t%f\n' % (pNum, unc, like, diff))

            for pNum in range(self.data.nParts):
                h = statsHashList[pNum]
                # pf.p4_expectedCompositionCounts returns a tuple of tuples
                # representing the counts of the nodes in proper alignment
                # order.
                ret = pf.p4_expectedCompositionCounts(evalTree.cTree, pNum)
                h['expectedIndividualCounts'] = list(ret)  # alignment order
                # print h['expectedIndividualCounts']
                h['overallSimStat'] = 0.0
                h['individualSimStats'] = [0.0] * self.data.nTax
                for seqNum in range(self.data.nTax):
                    if seqNum in skipTaxNums[pNum]:
                        pass
                    else:
                        obsCounts = list(
                            pf.singleSequenceBaseCounts(self.data.parts[pNum].cPart, seqNum))
                        # obsCounts is the counts observed in the simulation.
                        # It assumes that there are no gaps.  If there are
                        # gaps, adjust the counts.
                        if h['individualNSites'][seqNum] != self.data.parts[pNum].nChar:
                            factor = float(
                                h['individualNSites'][seqNum]) / self.data.parts[pNum].nChar
                            # print "factor = %s" % factor
                            for j in range(self.data.parts[pNum].dim):
                                obsCounts[j] = float(obsCounts[j]) * factor
                        # print obsCounts
                        for j in range(self.data.parts[pNum].dim):
                            # Avoid dividing by Zero.
                            if h['expectedIndividualCounts'][seqNum][j]:
                                dif = obsCounts[
                                    j] - h['expectedIndividualCounts'][seqNum][j]
                                h['individualSimStats'][
                                    seqNum] += ((dif * dif) / h['expectedIndividualCounts'][seqNum][j])
                        h['overallSimStat'] += h['individualSimStats'][seqNum]

                f3.write('%i\t' % pNum)
                for seqNum in range(self.data.nTax):
                    f3.write('%f\t' % h['individualSimStats'][seqNum])
                f3.write('%f\n' % h['overallSimStat'])
                # print h['overallSimStat']

        if doGoldmanCox:
            f2.close()
        f3.close()

        # Replace the saved data
        # Since we are replacing an exisiting data, this triggers
        # self.deleteCStuff()
        self.data = savedData

    def modelFitTests(self, fName='model_fit_tests_out', writeRawStats=0):
        """Do model fit tests on the data.

        The two tests are the Goldman-Cox test, and the tree- and model-
        based composition fit test.  Both require simulations with
        optimizations in order to get a null distribution, and those
        simulations need to be done before this method.  The simulations
        should be done with the simsForModelFitTests() method.

        Self should have a data and a model attached, and be optimized.

        The Goldman-Cox test (Goldman 1993.  Statistical tests of models
        of DNA substitution.  J Mol Evol 36: 182-198.) is a test for
        overall fit of the model to the data.  It does not work if the
        data have gaps or ambiguities.

        The tree- and model-based composition test asks the question:
        'Does the composition implied by the model fit the data?'  If the
        model is homogeneous and empirical comp is used, then this is the
        same as the chi-square test except that the null distribution
        comes from simulations, not from the chi-square distribution.  In
        that case only the question is, additionally, 'Are the data
        homogeneous in composition?', ie the same question asked by the
        chi-square test.  However, the data might be heterogeneous, and
        the model might be heterogeneous over the tree; the tree- and
        model-based composition fit test can ask whether the heterogeneous
        model fits the heterogeneous data.  The composition is tested in
        each data partition, separately.  The test is done both overall,
        ie for all the sequences together, and for individual sequences.

        If you just want a compo homogeneity test with empirical
        homogeneous comp, try the compoTestUsingSimulations() method-- its
        way faster, because there are not optimizations in the sims part.

        Output is verbose, to a file."""

        gm = ['Tree.modelFitTests()']
        self.calcLogLike(verbose=0)
        # Usually True.  Set to False for debugging, experimentation, getting
        # individual stats, etc
        doOut = True

        # We can't do the Goldman-Cox test if there are any gaps or
        # ambiguities.
        doGoldmanCox = True
        for a in self.data.alignments:
            if a.hasGapsOrAmbiguities():
                doGoldmanCox = False
                break
        # print "test doGoldmanCox = %s" % doGoldmanCox

        rawFName = '%s_raw.py' % fName
        #flob = sys.stderr
        #fRaw = sys.stderr
        if doOut:
            flob = file(fName, 'w')
        else:
            flob = None
        if writeRawStats:
            fRaw = file(rawFName, 'w')
        else:
            fRaw = None

        #######################
        # Goldman-Cox stats
        #######################

        # For a two-part data analysis, the first few lines of the
        # sims_GoldmanStats_* file will be like the following.  Its in
        # groups of 3-- the first one for all parts together (part number
        # -1), and the next lines for separate parts.

        # part  unconstr L       log like       Goldman-Cox stat
        # -1      -921.888705     -1085.696919    163.808215
        # 0       -357.089057     -430.941958     73.852901
        # 1       -564.799648     -654.754962     89.955314
        # -1      -952.063037     -1130.195799    178.132761
        # 0       -362.164119     -439.709824     77.545705
        # ... and so on.

        # For a one-part analysis, it will be the same except that one sim
        # gets only one line, starting with zero.

        if doGoldmanCox:
            goldmanOverallSimStats = []
            if self.data.nParts > 1:
                goldmanIndividualSimStats = []
                for partNum in range(self.data.nParts):
                    goldmanIndividualSimStats.append([])

            import glob
            goldmanFNames = glob.glob('sims_GoldmanStats_*')
            # print "nParts=%s" % self.data.nParts
            # print goldmanFNames
            for fName1 in goldmanFNames:
                f2 = open(fName1)
                aLine = f2.readline()
                if not aLine:
                    gm.append("Empty file %s" % fName1)
                    raise P4Error(gm)
                if aLine[0] != '#':
                    gm.append(
                        "Expecting a '#' as the first character in file %s" % fName1)
                    raise P4Error(gm)
                aLine = f2.readline()
                # print "a got line %s" % aLine,
                while aLine:
                    if self.data.nParts > 1:
                        splitLine = aLine.split()
                        if len(splitLine) != 4:
                            gm.append(
                                "Bad line in Goldman stats file %s" % fName1)
                            gm.append("'%s'" % aLine)
                            raise P4Error(gm)
                        if int(splitLine[0]) != -1:
                            gm.append(
                                "Bad line in Goldman stats file %s" % fName1)
                            gm.append("First item should be -1")
                            gm.append("'%s'" % aLine)
                            raise P4Error(gm)
                        # print splitLine[-1]
                        goldmanOverallSimStats.append(float(splitLine[-1]))

                        aLine = f2.readline()
                        # print "b got line %s" % aLine,
                        if not aLine:
                            gm.append("Premature end to file %s" % fName1)
                            raise P4Error(gm)

                    for partNum in range(self.data.nParts):
                        splitLine = aLine.split()
                        # print "partNum %i, splitLine=%s" % (partNum,
                        # splitLine)
                        if len(splitLine) != 4:
                            gm.append(
                                "Bad line in Goldman stats file %s" % fName1)
                            gm.append("'%s'" % aLine)
                            raise P4Error(gm)
                        try:
                            splitLine[0] = int(splitLine[0])
                        except ValueError:
                            gm.append(
                                "Bad line in Goldman stats file %s" % fName1)
                            gm.append(
                                "First item should be the partNum %i" % partNum)
                            gm.append("'%s'" % aLine)
                            raise P4Error(gm)
                        if splitLine[0] != partNum:
                            gm.append(
                                "Bad line in Goldman stats file %s" % fName1)
                            gm.append(
                                "First item should be the partNum %i" % partNum)
                            gm.append("'%s'" % aLine)
                            raise P4Error(gm)
                        # for taxNum in range(self.data.nTax):
                        #    print splitLine[taxNum + 1]
                        # print splitLine[-1]
                        if self.data.nParts == 1:
                            goldmanOverallSimStats.append(float(splitLine[-1]))
                        else:
                            goldmanIndividualSimStats[
                                partNum].append(float(splitLine[-1]))

                        aLine = f2.readline()
                        # print "c got line %s" % aLine,
                f2.close()

            # print "goldmanOverallSimStats =", goldmanOverallSimStats
            # print "goldmanIndividualSimStats =", goldmanIndividualSimStats
            # sys.exit()

            if doOut:
                flob.write('Model fit tests\n===============\n\n')
                flob.write(
                    'The data that we are testing have %i taxa,\n' % self.data.nTax)

                if len(self.data.alignments) == 1:
                    flob.write('1 alignment, ')
                else:
                    flob.write('%i alignments, ' % len(self.data.alignments))
                if self.data.nParts == 1:
                    flob.write('and 1 data partition.\n')
                else:
                    flob.write('and %i data partitions.\n' % self.data.nParts)

                flob.write('The lengths of those partitions are as follows:\n')
                flob.write('                  partNum    nChar \n')
                for i in range(self.data.nParts):
                    flob.write('                      %3i   %5i\n' %
                               (i, self.data.parts[i].nChar))
            self.data.calcUnconstrainedLogLikelihood2()
            if doOut:
                flob.write("\nThe unconstrained likelihood is %f\n" %
                           self.data.unconstrainedLogLikelihood)
                flob.write(
                    '(This is the partition-by-partition unconstrained log likelihood, \n')
                flob.write(
                    'ie the sum of the unconstrained log likes from each partition separately, \n')
                flob.write(
                    'and so will not be the same as that given by PAUP, if the data are partitioned.)\n')

                flob.write('\n\nGoldman-Cox test for overall model fit\n')
                flob.write('======================================\n')
                flob.write(
                    'The log likelihood for these data for this tree is %f\n' % self.logLike)
                flob.write('The unconstrained log likelihood for these data is %f\n' %
                           self.data.unconstrainedLogLikelihood)
            originalGoldmanCoxStat = self.data.unconstrainedLogLikelihood - \
                self.logLike
            if doOut:
                flob.write(
                    'The Goldman-Cox statistic for the original data is the difference, %f\n' % originalGoldmanCoxStat)
                if self.data.nParts > 1:
                    flob.write(
                        '(The unconstrained log likelihood for these data is calculated partition by partition.)\n')
                flob.write('\n')

            if self.data.nParts > 1:
                originalGoldmanCoxStatsByPart = []
                if doOut:
                    flob.write('Stats by partition.\n')
                    flob.write(
                        'part\t unconstrLogL\t log like \tGoldman-Cox stat\n')
                    flob.write(
                        '----\t ----------\t -------- \t----------------\n')
                for partNum in range(self.data.nParts):
                    unc = pf.getUnconstrainedLogLike(
                        self.data.parts[partNum].cPart)
                    like = pf.p4_partLogLike(
                        self.cTree, self.data.parts[partNum].cPart, partNum, 0)
                    diff = unc - like
                    if doOut:
                        flob.write('  %i\t%f\t%f\t   %f\n' %
                                   (partNum, unc, like, diff))
                    originalGoldmanCoxStatsByPart.append(diff)

            # Do the overall stat
            nSims = len(goldmanOverallSimStats)
            if doOut:
                flob.write('\nThere were %i simulations.\n\n' % nSims)

            if writeRawStats:
                fRaw.write('# Goldman-Cox null distributions.\n')
                if self.data.nParts > 1:
                    fRaw.write(
                        '# Simulation stats for overall data, ie for all data partitions combined.\n')
                else:
                    fRaw.write('# Simulation stats.\n')
                fRaw.write('goldman_cox_overall = %s\n' %
                           goldmanOverallSimStats)
                if self.data.nParts > 1:
                    for partNum in range(self.data.nParts):
                        fRaw.write(
                            '# Simulation stats for data partition %i\n' % partNum)
                        fRaw.write('goldman_cox_part%i = %s\n' %
                                   (partNum, goldmanIndividualSimStats[partNum]))

            prob = func.tailAreaProbability(
                originalGoldmanCoxStat, goldmanOverallSimStats, verbose=0)
            if doOut:
                flob.write('\n              Overall Goldman-Cox test: ')
                if prob <= 0.05:
                    flob.write('%13s' % "Doesn't fit.")
                else:
                    flob.write('%13s' % 'Fits.')
                flob.write('    P = %5.3f\n' % prob)

            if self.data.nParts > 1:
                if doOut:
                    flob.write('  Tests for individual data partitions:\n')
                for partNum in range(self.data.nParts):
                    prob = func.tailAreaProbability(originalGoldmanCoxStatsByPart[partNum],
                                                    goldmanIndividualSimStats[partNum], verbose=0)
                    if doOut:
                        flob.write(
                            '                               Part %-2i: ' % partNum)
                        if prob <= 0.05:
                            flob.write('%13s' % 'Doesn\'t fit.')
                        else:
                            flob.write('%13s' % 'Fits.')
                        flob.write('    P = %5.3f\n' % prob)

        #########################
        # COMPOSITION
        #########################

        statsHashList = []
        for pNum in range(self.data.nParts):
            h = {}
            statsHashList.append(h)
            h['individualNSites'] = []
            h['observedIndividualCounts'] = []
            for j in range(self.data.nTax):
                # print pf.partSequenceSitesCount(self.data.parts[pNum].cPart,
                # j)
                h['individualNSites'].append(
                    pf.partSequenceSitesCount(self.data.parts[pNum].cPart, j))  # no gaps or qmarks
                # print self.data.parts[pNum].composition([j])
                h['observedIndividualCounts'].append(
                    self.data.parts[pNum].composition([j]))
                # The line above is temporarily composition, not counts
            # pf.expectedCompositionCounts returns a tuple of tuples
            # representing the counts of the nodes in proper alignment order.
            h['expectedIndividualCounts'] = list(
                pf.p4_expectedCompositionCounts(self.cTree, pNum))  # alignment order

            # At the moment, h['observedIndividualCounts'] has composition,
            # not counts.  So multiply by h['individualNSites']
            for i in range(self.data.nTax):
                for j in range(self.data.parts[pNum].dim):
                    h['observedIndividualCounts'][i][
                        j] *= h['individualNSites'][i]

        # We will want to skip any sequences composed of all gaps
        skipTaxNums = []
        for pNum in range(self.data.nParts):
            stn = []
            for tNum in range(self.data.nTax):
                if not statsHashList[pNum]['individualNSites'][tNum]:
                    stn.append(tNum)
            skipTaxNums.append(stn)
        # print "skipTaxNums = %s" % skipTaxNums

        # Do the boring old compo chi square test.
        if doOut:
            flob.write(longMessage1)  # explanation ...
        for pNum in range(self.data.nParts):
            h = statsHashList[pNum]
            # Can't use func.xSquared(), because there might be column
            # zeros.
            # print "observedIndividualCounts = %s' %
            # h['observedIndividualCounts"]
            nRows = len(h['observedIndividualCounts'])
            nCols = len(h['observedIndividualCounts'][0])
            # I could have just used nSites, above
            theSumOfRows = func._sumOfRows(h['observedIndividualCounts'])
            theSumOfCols = func._sumOfColumns(h['observedIndividualCounts'])
            # print theSumOfCols
            isOk = 1
            columnZeros = []
            # for j in range(len(theSumOfRows)):
            #    if theSumOfRows[j] == 0.0:
            #        gm.append("Zero in a row sum.  Programming error.")
            #        raise P4Error(gm)
            for j in range(len(theSumOfCols)):
                if theSumOfCols[j] <= 0.0:
                    columnZeros.append(j)
            theExpected = func._expected(theSumOfRows, theSumOfCols)
            # print "theExpected = %s" % theExpected
            # print "columnZeros = %s" % columnZeros
            xSq = 0.0
            for rowNum in range(nRows):
                if rowNum in skipTaxNums[pNum]:
                    pass
                else:
                    xSq_row = 0.0
                    for colNum in range(nCols):
                        if colNum in columnZeros:
                            pass
                        else:
                            theDiff = h['observedIndividualCounts'][rowNum][
                                colNum] - theExpected[rowNum][colNum]
                            xSq_row += (theDiff * theDiff) / \
                                theExpected[rowNum][colNum]
                    xSq += xSq_row
            dof = (nCols - len(columnZeros) - 1) * \
                (nRows - len(skipTaxNums[pNum]) - 1)
            prob = func.chiSquaredProb(xSq, dof)
            if doOut:
                flob.write(
                    '        Part %i: Chi-square = %f, (dof=%i) P = %f\n' % (pNum, xSq, dof, prob))

        for pNum in range(self.data.nParts):
            h = statsHashList[pNum]
            h['overallStat'] = 0.0
            h['individualStats'] = [0.0] * self.data.nTax
            for i in range(self.data.nTax):
                if i in skipTaxNums[pNum]:
                    pass  # h['individualStats'] stays at zeros
                else:
                    for j in range(self.data.parts[pNum].dim):
                        # Avoid dividing by Zero.
                        if h['expectedIndividualCounts'][i][j]:
                            dif = h['observedIndividualCounts'][i][
                                j] - h['expectedIndividualCounts'][i][j]
                            h['individualStats'][
                                i] += ((dif * dif) / h['expectedIndividualCounts'][i][j])
                    h['overallStat'] += h['individualStats'][i]

            h['overallSimStats'] = []
            h['individualSimStats'] = []
            for i in range(self.data.nTax):
                h['individualSimStats'].append([])

            if 0:
                print "h['individualNSites'] = %s" % h['individualNSites']
                print "h['observedIndividualCounts'] = %s" % h['observedIndividualCounts']
                print "h['expectedIndividualCounts'] = %s" % h['expectedIndividualCounts']
                print "h['overallStat'] = %s" % h['overallStat']
                print "h['individualStats'] = %s" % h['individualStats']
                raise P4Error(gm)

        import glob
        compoFNames = glob.glob('sims_CompStats_*')
        # print compoFNames
        for fName1 in compoFNames:
            f2 = open(fName1)
            aLine = f2.readline()
            if not aLine:
                gm.append("Empty file %s" % fName1)
                raise P4Error(gm)
            # print "a got line %s" % aLine,
            while aLine:
                for partNum in range(self.data.nParts):
                    h = statsHashList[partNum]
                    splitLine = aLine.split()
                    if len(splitLine) != (self.data.nTax + 2):
                        gm.append(
                            "Bad line in composition stats file %s" % fName1)
                        gm.append("'%s'" % aLine)
                        raise P4Error(gm)
                    if int(splitLine[0]) != partNum:
                        gm.append(
                            "Bad line in composition stats file %s" % fName1)
                        gm.append(
                            "First item should be the partNum %i" % partNum)
                        gm.append("'%s'" % aLine)
                        raise P4Error(gm)
                    # for taxNum in range(self.data.nTax):
                    #    print splitLine[taxNum + 1]
                    # print splitLine[-1]
                    h['overallSimStats'].append(float(splitLine[-1]))
                    for i in range(self.data.nTax):
                        h['individualSimStats'][i].append(
                            float(splitLine[i + 1]))
                    #raise P4Error(gm)

                    aLine = f2.readline()
                    if not aLine:
                        break
                    # print "b got line %s" % aLine,
            f2.close()

        nSims = len(statsHashList[0]['overallSimStats'])
        if doOut:
            # Explain tree- and model-based compo fit stat, X^2_m
            flob.write(longMessage2)
            flob.write('    %i simulation reps were used.\n\n' % nSims)

        spacer1 = ' ' * 10
        for partNum in range(self.data.nParts):
            h = statsHashList[partNum]
            if doOut:
                flob.write('Part %-2i:\n-------\n\n' % partNum)
                flob.write('Statistics from the original data\n')
                flob.write('%s%30s: %f\n' %
                           (spacer1, 'Overall observed stat', h['overallStat']))
                flob.write('%s%30s:\n' %
                           (spacer1, 'Stats for individual taxa'))
                for taxNum in range(self.data.nTax):
                    if taxNum not in skipTaxNums[partNum]:
                        flob.write('%s%30s: %f\n' % (
                            spacer1, self.data.taxNames[taxNum], h['individualStats'][taxNum]))
                    else:
                        flob.write('%s%30s: skipped\n' %
                                   (spacer1, self.data.taxNames[taxNum]))

                flob.write(
                    '\nAssessment of fit from null distribution from %i simulations\n' % nSims)
                flob.write('%s%30s:  ' % (spacer1, 'Overall'))
            prob = func.tailAreaProbability(
                h['overallStat'], h['overallSimStats'], verbose=0)
            if doOut:
                if prob <= 0.05:
                    flob.write('%13s' % 'Doesn\'t fit.')
                else:
                    flob.write('%13s' % 'Fits.')
                flob.write('    P = %5.3f\n' % prob)
            #############
            theRet = prob
            #############
            for taxNum in range(self.data.nTax):
                if doOut:
                    flob.write('%s%30s:  ' %
                               (spacer1, self.data.taxNames[taxNum]))
                if taxNum in skipTaxNums[partNum]:
                    if doOut:
                        flob.write('%13s\n' % 'skipped.')
                else:
                    prob = func.tailAreaProbability(h['individualStats'][taxNum],
                                                    h['individualSimStats'][taxNum], verbose=0)
                    if doOut:
                        if prob <= 0.05:
                            flob.write('%13s' % "Doesn't fit.")
                        else:
                            flob.write('%13s' % 'Fits.')
                        flob.write('    P = %5.3f\n' % prob)

            if writeRawStats:
                fRaw.write('#\n# Tree and model based composition fit test\n')
                fRaw.write('# =========================================\n')
                fRaw.write(
                    '# Simulation statistics, ie the null distributions\n\n')
                fRaw.write('# Part %i:\n' % partNum)
                fRaw.write('part%i_overall_compo_null = %s\n' %
                           (partNum, h['overallSimStats']))
                for taxNum in range(self.data.nTax):
                    fRaw.write('part%i_%s_compo_null = %s\n' % (partNum,
                                                                _fixFileName(
                                                                    self.data.taxNames[taxNum]),
                                                                h['individualSimStats'][taxNum]))

        # Yes, it is possible to close sys.stdout
        if flob and flob != sys.stdout:
            flob.close()
        if fRaw and fRaw != sys.stdout:
            fRaw.close()
        return theRet

    def compoTestUsingSimulations(self, nSims=100, doIndividualSequences=0, doChiSquare=0, verbose=1):
        """Compositional homogeneity test using a null distribution from simulations.

        This does a compositional homogeneity test on each data partition.
        The statistic used here is X^2, obtained via
        Data.compoChiSquaredTest().

        The null distribution of the stat is made using simulations, so of
        course you need to provide a tree with a model, with optimized
        branch lengths and model parameters.  This is a comp homogeneity
        test, so the model should be tree-homogeneous.

        The analysis usually tests all sequences in the data partition
        together (like paup), but you can also 'doIndividualSequences'
        (like puzzle).  Beware that the latter is a multiple simultaneous
        stats test, and so the power may be compromized.

        For purposes of comparison, this test can also do compo tests in
        the style of PAUP and puzzle, using chi-square to assess
        significance.  Do this by turning 'doChiSquare' on.  The compo test
        in PAUP tests all sequences together, while the compo test in
        puzzle tests all sequences separately.  There are advantages and
        disadvantages to the latter-- doing all sequences separately
        allows you to identify the worst offenders, but suffers due to the
        problems of multiple simultaneous stats tests.  There are slight
        differences between the computation of the Chi-square in PAUP and
        puzzle and the p4 version.  The compo test in PAUP (basefreq) does
        the chi-squared test, but if sequences are blank it still counts
        them in the degrees of freedom; p4 does not count blank sequences
        in the degrees of freedom.  Puzzle simply uses the row sums, ie
        the contributions of each sequence to the total X-squared, and
        assesses significance with chi-squared using the number of symbols
        minus 1 as the degrees of freedom.  Ie for DNA dof=3, for protein
        dof=19.  Puzzle correctly gets the composition from sequences with
        gaps, but does not do the right thing for sequences with
        ambiguities like r, y, and so on.  P4 does calculate the
        composition correctly when there are such ambiguities.  So p4 will
        give you the same numbers as paup and puzzle for the chi-squared
        part as long as you don't have blank sequences or ambiguities like
        r and y.

        This uses the Data.compoChiSquaredTest() method to get the
        stats. See the doc string for that method, where it describes how
        zero column sums (ie some character is absent) can be dealt with.
        Here, when that method is invoked, 'skipColumnZeros' is turned on,
        so that the analysis is robust against data with zero or low
        values for some characters.

        For example::

            # First, do a homog opt, and pickle the optimized tree.
            # Here I use a bionj tree, but you could use whatever.
            read('d.nex')
            a = var.alignments[0]
            dm = a.pDistances()
            t = dm.bionj()
            d = Data()
            t.data = d
            t.newComp(free=1, spec='empirical')
            t.newRMatrix(free=1, spec='ones')
            t.setNGammaCat(nGammaCat=4)
            t.newGdasrv(free=1, val=0.5)
            t.setPInvar(free=0, val=0.0)
            t.optLogLike()
            t.name = 'homogOpt'
            t.tPickle()

            # Then, do the test ...
            read('homogOpt.p4_tPickle')
            t = var.trees[0]
            read('d.nex')
            d = Data()
            t.data = d
            t.compoTestUsingSimulations()

            # Output would be something like ...
            # Composition homogeneity test using simulations.
            # P-values are shown.

            #             Part Num       0     
            #            Part Name      all    
            # --------------------    -------- 
            #        All Sequences     0.0000  

            # Or using more sims for more precision, and also doing the
            # Chi-square test for contrast ...

            t.compoTestUsingSimulations(nSims=1000, doChiSquare=True)

            # Output might be something like ...
            # Composition homogeneity test using simulations.
            # P-values are shown.
            # (P-values from Chi-Square are shown in parens.)

            #             Part Num       0     
            #            Part Name      all    
            # --------------------    -------- 
            #        All Sequences     0.0140  
            #   (Chi-Squared Prob)    (0.9933) 

        It is often the case, as above, that this test will show
        significance while the Chi-square test does not.

        """

        gm = ['Tree.compoTestUsingSimulations()']

        # print "inComp = %s" % self.model.parts[0].comps[0].val

        if not self.data:
            gm.append("No data.  Set the data first.")
            raise P4Error(gm)
        if not self.model:
            gm.append("No model.  You need to set the model first.")
            raise P4Error(gm)
        self.modelSanityCheck()
        if self.model.isHet:
            gm.append("The model for this tree is tree-heterogeneous.")
            gm.append("This test is not valid for tree-hetero models.")
            raise P4Error(gm)

        # Make a new data object in which to do the sims, so we do not over-write self
        # print "a self.data = %s" % self.data
        # self.data.dump()
        savedData = self.data
        self.data = None  # This triggers self.deleteCStuff()
        self.data = savedData.dupe()

        # print "b self.data = %s" % self.data
        # self.data.dump()
        #raise P4Error(gm)

        # Check for missing sequences in any of the parts.  Missing seq
        # nums go in skips, a list of lists.
        skips = []
        for pNum in range(self.data.nParts):
            skips.append([])
        for pNum in range(self.data.nParts):
            for tNum in range(self.data.nTax):
                nSites = pf.partSequenceSitesCount(
                    self.data.parts[pNum].cPart, tNum)  # no gaps, no missings
                if not nSites:
                    skips[pNum].append(tNum)

        # Get the original stats from self.data.
        # compoChiSquaredTest(self, verbose=1, skipColumnZeros=0, useConstantSites=1, skipTaxNums=None, getRows=0)
        original = self.data.compoChiSquaredTest(verbose=0,
                                                 skipColumnZeros=1,
                                                 skipTaxNums=skips,
                                                 getRows=doIndividualSequences)
        # print "original =", original

        # Make some empty lists in which to put our stats
        full = []
        if doIndividualSequences:
            rows = []
        for pNum in range(self.data.nParts):
            full.append([])
            if doIndividualSequences:
                onePartRows = []
                for i in range(self.data.nTax):
                    onePartRows.append([])
                rows.append(onePartRows)

        # Do the sims
        for i in range(nSims):
            # if i < 5:
            # print "%i simComp = %s" % (i, self.model.parts[0].comps[0].val)
            self.simulate()
            ret = self.data.compoChiSquaredTest(skipColumnZeros=1,
                                                skipTaxNums=skips, getRows=doIndividualSequences, verbose=0)
            # print "%i ret=%s" % (i, ret)
            for pNum in range(self.data.nParts):
                full[pNum].append(ret[pNum][0])
                if doIndividualSequences:
                    for tNum in range(self.data.nTax):
                        if tNum not in skips[pNum]:
                            rows[pNum][tNum].append(ret[pNum][3][tNum])

        # Find the longest part name length, and heading width, so the output
        # looks nice.
        partWid = 8
        for p in self.data.parts:
            if len(p.name) > partWid:
                partWid = len(p.name)
        partWid += 2

        headWid = 20
        for tN in self.data.taxNames:
            if len(tN) > headWid:
                headWid = len(tN)
        headWid += 2
        #headSig = '%-' + `headWid` + 's'
        headSig = '%' + `headWid - 2` + 's  '

        # Get the all-sequences tail area probs
        partTaps = []
        for pNum in range(self.data.nParts):
            partTaps.append(
                func.tailAreaProbability(original[pNum][0], full[pNum], verbose=0))

        # Intro
        if verbose:
            print "Composition homogeneity test using simulations."
            print "P-values are shown."
            if doChiSquare:
                print "(P-values from Chi-Square are shown in parens.)"
            print

        # Print the Part Nums and Part Names
        if verbose:
            print headSig % 'Part Num',
            for pNum in range(self.data.nParts):
                print string.center('%i' % pNum, partWid),
            print
            print headSig % 'Part Name',
            for pNum in range(self.data.nParts):
                print string.center('%s' % self.data.parts[pNum].name, partWid),
            print
            print headSig % ('-' * (headWid - 2)),
            for pNum in range(self.data.nParts):
                print string.center('%s' % ('-' * (partWid - 2)), partWid),
            print

        # Print the all-sequences results
        if verbose:
            print headSig % 'All Sequences',
            for pNum in range(self.data.nParts):
                print string.center('%6.4f' % partTaps[pNum], partWid),
            print
            if doChiSquare:
                print headSig % '(Chi-Squared Prob)',
                for pNum in range(self.data.nParts):
                    print string.center('(%6.4f)' % original[pNum][2], partWid),
                print

        if doIndividualSequences and verbose:
            print
            # print "Individual sequences"
            # print "--------------------"

            for tNum in range(self.data.nTax):
                print headSig % self.data.taxNames[tNum],
                for pNum in range(self.data.nParts):
                    if tNum not in skips[pNum]:
                        ret = func.tailAreaProbability(
                            original[pNum][3][tNum], rows[pNum][tNum], verbose=0)
                        print string.center('%6.4f' % ret, partWid),
                    else:
                        print string.center('%s' % ('-' * 4), partWid),
                print
                if doChiSquare:
                    print headSig % ' ',
                    for pNum in range(self.data.nParts):
                        # degrees of freedom
                        dof = self.data.parts[pNum].dim - 1
                        if tNum not in skips[pNum]:
                            ret = func.chiSquaredProb(
                                original[pNum][3][tNum], dof)
                            print string.center('(%6.4f)' % ret, partWid),
                        else:
                            print string.center('%s' % ('-' * 4), partWid),
                    print

        # Replace the saved data
        # Since we are replacing an exisiting data, this triggers
        # self.deleteCStuff()
        self.data = savedData

        return partTaps[0]

    def bigXSquaredSubM(self, verbose=False):
        """Calculate the X^2_m stat

        This can handle gaps and ambiguities.

        Column zeros in the observed is not a problem with this stat,
        as we are dividing by the expected composition, and that comes
        from the model, which does not allow compositions with values
        of zero.  """

        if not self.cTree:
            self._commonCStuff(resetEmpiricalComps=True)
        l = []
        for pNum in range(self.data.nParts):
            if verbose:
                print "Part %i" % pNum
                print "======"
            obs = []
            nSites = []  # no gaps or ?
            for taxNum in range(self.nTax):
                thisNSites = pf.partSequenceSitesCount(
                    self.data.parts[pNum].cPart, taxNum)
                comp = self.data.parts[pNum].composition([taxNum])
                for symbNum in range(self.data.parts[pNum].dim):
                    comp[symbNum] *= thisNSites
                nSites.append(thisNSites)
                obs.append(comp)
            if verbose:
                print "\n  Observed"
                print " " * 10,
                for symbNum in range(self.data.parts[pNum].dim):
                    print "%8s" % self.data.parts[pNum].symbols[symbNum],
                print
                for taxNum in range(self.nTax):
                    print "%10s" % self.taxNames[taxNum],
                    for symbNum in range(self.data.parts[pNum].dim):
                        print "%8.4f" % obs[taxNum][symbNum],
                    print "   n=%i" % nSites[taxNum]

            # pf.p4_expectedCompositionCounts returns a tuple of tuples
            # representing the counts of the nodes in proper alignment order.
            exp = list(pf.p4_expectedCompositionCounts(self.cTree, pNum))
            if verbose:
                print "\n  Expected"
                print " " * 10,
                for symbNum in range(self.data.parts[pNum].dim):
                    print "%8s" % self.data.parts[pNum].symbols[symbNum],
                print
                for taxNum in range(self.nTax):
                    print "%10s" % self.taxNames[taxNum],
                    for symbNum in range(self.data.parts[pNum].dim):
                        print "%8.4f" % exp[taxNum][symbNum],
                    print "   n=%i" % nSites[taxNum]

            # do the summation
            theSum = 0.0
            for taxNum in range(self.nTax):
                for symbNum in range(self.data.parts[pNum].dim):
                    x = obs[taxNum][symbNum] - exp[taxNum][symbNum]
                    theSum += (x * x) / exp[taxNum][symbNum]

            l.append(theSum)
            if verbose:
                print "The bigXSquaredSubM stat for this part is %.5f" % theSum
        return l

    def compStatFromCharFreqs(self, verbose=False):
        """Calculate a statistic from observed and model character frequencies.

        Call it c_m, little c sub m.

        It is calculated from observed character frequencies and character
        frequencies expected from the (possibly tree-heterogeneous) model.

        It would be the sum of abs(obs-exp)/exp
        """

        if not self.cTree:
            self._commonCStuff(resetEmpiricalComps=True)
        # pf.p4_expectedComposition returns a tuple of tuples of tuples
        # representing the counts of the nodes in proper alignment order.
        exp = pf.p4_expectedComposition(self.cTree)
        l = []
        for pNum in range(self.data.nParts):
            if verbose:
                print "Part %i" % pNum
                print "======"
            obs = []
            for taxNum in range(self.nTax):
                comp = self.data.parts[pNum].composition([taxNum])
                obs.append(comp)
            if verbose:
                print "\n  Observed"
                print " " * 10,
                for symbNum in range(self.data.parts[pNum].dim):
                    print "%8s" % self.data.parts[pNum].symbols[symbNum],
                print
                for taxNum in range(self.nTax):
                    print "%10s" % self.taxNames[taxNum],
                    for symbNum in range(self.data.parts[pNum].dim):
                        print "%8.4f" % obs[taxNum][symbNum],
                    print

            if verbose:
                print "\n  Expected"
                print " " * 10,
                for symbNum in range(self.data.parts[pNum].dim):
                    print "%8s" % self.data.parts[pNum].symbols[symbNum],
                print
                for taxNum in range(self.nTax):
                    print "%10s" % self.taxNames[taxNum],
                    for symbNum in range(self.data.parts[pNum].dim):
                        print "%8.4f" % exp[pNum][taxNum][symbNum],
                    print

            # do the summation
            theSum = 0.0
            for taxNum in range(self.nTax):
                for symbNum in range(self.data.parts[pNum].dim):
                    theSum += math.fabs(obs[taxNum][symbNum] - exp[pNum]
                                        [taxNum][symbNum]) / exp[pNum][taxNum][symbNum]

            l.append(theSum)
            if verbose:
                print "The c_m stat for this part is %.5f" % theSum
        return l


class RMatrix(object):

    def __init__(self):
        self.num = -1
        self.partNum = None
        self.free = None
        self.spec = None
        self.symbol = None
        self.val = None
        self.nNodes = 0
        self.rj_isInPool = False
        self.rj_f = 0.0


class Gdasrv(object):

    def __init__(self):
        self.num = -1
        self.partNum = None
        self.free = None
        self.symbol = None
        # self.val=None
        self._val = numpy.zeros(1, numpy.float)
        self.freqs = None
        self.rates = None
        self.nGammaCat = None
        self.c = None    # a p4_gdasrvStruct, if it exists.
        self.nNodes = 0

    def _setVal(self, theVal):
        if theVal < 1.e-16:
            gm = ["Gdasrv._setVal()"]
            gm.append("Attempt to set Gdasrv.val (ie alpha) to %g" % theVal)
            gm.append(
                "However, we cannot calculate the discrete categories with a value so low.")
            raise P4Error(gm)
        self._val[0] = theVal
        self.calcRates()

    val = property(lambda self: self._val, _setVal)

    def calcRates(self):
        # Use either the p4_gdasrvStruct, or just use the NumPy
        # array vals (np = NumPy).
        # print "self.c = %s" % self.c
        if self.c:
            pf.gdasrvCalcRates(self.c)
        else:
            pf.gdasrvCalcRates_np(
                self.nGammaCat, self._val[0], self.freqs, self.rates)
        # print 'xxx self.rates = %s, val=%s' % (self.rates, self._val[0])


class Comp(object):

    def __init__(self):
        self.num = -1
        self.partNum = None
        self.free = None
        self.spec = None
        self.symbol = None
        self.val = None
        self.nNodes = 0
        self.rj_isInPool = False
        self.rj_f = 0.0


class PInvar(object):

    def __init__(self):
        self.num = -1
        self.partNum = None
        self.free = None
        self.val = None


class TSCovarion(object):

    def __init__(self):
        self.partNum = None
        self.free = None
        self.s1 = None
        self.s2 = None


class Mixture(object):

    def __init__(self):
        self.partNum = None
        self.free = None
        self.freqs = None
        self.rates = None
        self.nMix = None
