from __future__ import print_function
import sys
import string
import types
import math
import copy
import os
import p4.func
import time
import glob
import pickle
from p4.var import var
from p4.p4exceptions import P4Error
from p4.node import Node, NodePart, NodeBranch, NodeBranchPart
from p4.distancematrix import DistanceMatrix

import numpy
import p4.pf as pf
from p4.model import Model
from p4.data import Data
from p4.alignment import Part
import random

if True:
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
        f = open(fN, 'wb')
        pickle.dump(self.dupe(), f, pickle.HIGHEST_PROTOCOL)  
        # pickle.dump(self, f, 1) # Don't do this -- has pointers that would
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
                        print("Tree.writeNexus()")
                        print("    'append' is requested,")
                        print("    but '%s' is not a regular file (doesn't exist?)." \
                              % fName)
                        print("    Writing to a new file instead.")
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
                f.write(' %s' % p4.func.nexusFixNameIfQuotesAreNeeded(i))
            f.write(';\nend;\n\n')

        f.write('begin trees;\n')
        if message:
            f.write('  [%s\n  ]\n' % message)
        if self.logLike:
            f.write('  [logLike for tree %s is %f]\n' %
                    (self.name, self.logLike))

        f.write('  tree %s = [&U] ' %
                p4.func.nexusFixNameIfQuotesAreNeeded(self.name))
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

    def writeNewick(self, fName=None, withTranslation=0, translationHash=None, doMcmcCommandComments=0, toString=False, append=False, spaceAfterComma=True):
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
                        '%s' % p4.func.nexusFixNameIfQuotesAreNeeded(self.root.name))
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
                                    '%s' % p4.func.nexusFixNameIfQuotesAreNeeded(n1.name))
                            else:
                                if n1 != self.root:
                                    gm.append("Terminal node with no name?")
                                    raise P4Error(gm)
                    else:
                        sList.append(')')
                        if n1.name:
                            sList.append(
                                '%s' % p4.func.nexusFixNameIfQuotesAreNeeded(n1.name))
                    if writeBrLens:
                        if n1 != self.root:
                            sList.append(':%g' % n1.br.len)
                    if doMcmcCommandComments:
                        sList.append(self._getMcmcCommandComment(n1))
                    if n1.sibling:
                        if spaceAfterComma:
                            sList.append(', ')
                        else:
                            sList.append(',')
                        break
        sList.append(';\n')
        if toString:
            return "".join(sList)
        elif fName == None:
            print("".join(sList))
        elif type(fName) == type('string'):
            if append:
                fName2 = open(fName, 'a')
            else:
                fName2 = open(fName, 'w')
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
        print(string.join(s, '\n'))

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
                        print("Compositions for part %i" % partNum)
                    else:
                        print("\nCompositions\n------------")
                p.textDrawModelThing = var.TEXTDRAW_COMP
                p._setPos(autoIncreaseWidth)
                s = p.textString(returnAsList=True)
                # print s

            # Then do RMatrices
            if doRMatrices:
                if 0:
                    if self.model.nParts > 1:
                        print("RMatrices for part %i" % partNum)
                    else:
                        print("\nRMatrices\n---------")
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
        f = open(fName, 'w')
        f.write(s)
        f.close()

    #############################################
    # Various Tree methods for defining models
    #############################################

