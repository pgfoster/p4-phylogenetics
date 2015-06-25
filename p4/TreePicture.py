import string,os
from Var import var
from Glitch import Glitch

"""This class is used by Tree.draw(), Tree.eps(), and Tree.svg().
This week, there is no 'user-interface' for it, other than those two
methods."""


class TreePicture(object):
    def __init__(self, theTree = None):
        gm = ['TreePicture.__init__()']
        self.tree = None
        if theTree:
            try:
                if len(theTree.nodes) < 2:
                    if len(theTree.nodes) == 1:
                        gm.append("Only 1 node in this tree, so it can't be drawn.")
                    else:
                        gm.append("No nodes.")
                    raise Glitch, gm
            except AttributeError:
                gm.append("Expecting a Tree instance.")
                raise Glitch, gm

            self.tree = theTree
            if self.tree:
                for n in self.tree.nodes:
                    n.xPos = None
                    n.yPos = None

        self.width = None
        self.fName = None
        self.heightFactor = None
        self.pointsPerLetter = None
        self.textSize = 11
        self.labelTextSize = 9
        self.nameOffset = self.textSize / 3.0
        self.nameDrop = self.textSize / 4.0

        self.showInternalNodeNames = 1
        self.textShowNodeNums = 1
        self.putInternalNamesOnBranches = 0
        self.addToBrLen = 0.0

        self.xOrigin = 0.0
        self.yOrigin = 0.0
        self.xScale = None
        self.yScale = None

        self.partNum = -1
        self.textDrawModelThing = None

        self.svg = False



    ##Ignore
    def setPos(self, autoIncreaseWidth=True):
        gm = ['TreePicture.setPos()']

        # set xPos
        self.tree.root.xPos = 0.0
        for n in self.tree.iterPreOrder():
            if n != self.tree.root:
                #print "n=%s, n.parent=%s" % (n.nodeNum, n.parent)
                #print "n.parent.xPos=%s, n.br.len=%s, self.addToBrLen=%s" % (n.parent.xPos, n.br.len, self.addToBrLen)
                if n.parent.xPos == None:
                    raise Glitch, "Programming error. node %i n.parent.xPos is None" % n.nodeNum
                if n.br == None:
                    raise Glitch, "Programming error. node %i n.br is None" % n.nodeNum
                n.xPos = n.parent.xPos + n.br.len + self.addToBrLen
                    

        # set yPos
        counter = len([n for n in self.tree.iterLeavesNoRoot()])
        for n in self.tree.iterPostOrder():
            if n.isLeaf:
                if n == self.tree.root:
                    pass
                else:
                    n.yPos = float(counter)
                    counter -= 1
            else:
                try:
                    n.yPos = (n.leftChild.yPos + n.rightmostChild().yPos) / 2.0
                except:
                    print "-" * 50
                    print "TreePicture.setPos()  got an exception."
                    print "self.postOrder is %s" % self.tree.postOrder
                    print "n.nodeNum=%i, n.name=%s " % (n.nodeNum, n.name)
                    
                    if not n.leftChild:
                        print "There is no n.leftChild"
                    else:
                        print "n.leftChild.nodeNum ", n.leftChild.nodeNum
                        print "n.leftChild.isLeaf=%s" % n.leftChild.isLeaf
                        print "n.leftChild.yPos ", n.leftChild.yPos
                    print "n.rightmostChild().nodeNum ", n.rightmostChild().nodeNum
                    print "n.rightmostChild().isLeaf=%s" % n.rightmostChild().isLeaf
                    print "n.rightmostChild().yPos ", n.rightmostChild().yPos
                    raise Glitch, gm

        if self.tree.root.isLeaf:
            self.tree.root.yPos = self.tree.root.leftChild.yPos

        if 0:
            print "index   isLeaf        xPos         yPos"
            for n in self.tree.iterNodes():
                print "%4i" % n.nodeNum,
                print "%8i" % n.isLeaf,
                print "%12s" % n.xPos,
                print "%12s" % n.yPos

        # Find maxX
        self.maxX = 0.0
        for n in self.tree.iterNodes():
            if n.xPos > self.maxX:
                self.maxX = n.xPos

        if self.maxX < 1.e-100:
            gm.append("maxX (%f, %g) is too small.  Can't draw it." % (self.maxX, self.maxX))
            raise Glitch, gm
        # Find maxY
        self.maxY = 0.0
        for n in self.tree.iterNodes():
            if n.yPos > self.maxY:
                self.maxY = n.yPos
        #print "maxX is %f, maxY is %f" % (self.maxX, self.maxY)
        if self.maxY < 1.e-100:
            gm.append("maxY (%f, %g) is too small.  How did that happen?" % (self.maxY, self.maxY))
            raise Glitch, gm
        

        


        # The variable 'nameLengths' is a list of 2-item tuples.
        # First item is the xPos, and second item is the length of the
        # node name (plus any extras -- nodeNum and colon).  Eg,
        # without any extras, for a tree ((xyz, C:0.2), wzyz:0.2),
        # nameLengths = [(0.2, 3), (0.4, 1), (0.2, 4)]
        nameLengths = []
        for n in self.tree.iterNodes():
            if n != self.tree.root:
                if n.name:
                    if not self.showInternalNodeNames and not n.isLeaf:
                        pass
                    else:
                        # The extras would include the nodeNum length
                        # plus a colon, but only if
                        # self.textShowNodeNums is turned on.
                        if n.isLeaf:
                            extras = 0
                            if self.textShowNodeNums:
                                extras += len(`n.nodeNum`) + 1 # 1 for the colon
                        else:
                            if self.textShowNodeNums:
                                extras = 0 # So we do over-write a vertical line in text draw
                                extras += len(`n.nodeNum`) + 1 # 1 for the colon
                            else:
                                extras = 1 # So we do not over-write a vertical line in text draw

                        #print "%i: '%s'" % (n.nodeNum, n.name)
                        nameLengths.append((n.xPos, len('%s' % n.name) + extras))
        #print "nameLengths = %s" % nameLengths

        # Same idea for the root name.  It will extend to the left.
        if self.tree.root.xPos != 0:
            gm.append("Root xPos is not zero.  Fix me.")
            raise Glitch, gm

        if not self.tree.root.name:
            rootNameLength = 0
        else:
            if not self.showInternalNodeNames and not self.tree.root.isLeaf:
                rootNameLength = 0
            else:
                # The node num (if it is placed at all), will extend
                # to the right, so I don't need to account for it.
                rootNameLength = len(self.tree.root.name)
                if self.textShowNodeNums:
                    rootNameLength += 1
        #print "rootNameLength is %s" % rootNameLength


        # In 12 point Helvetica-Oblique, the M is 10 points wide, and the
        # i is 2.75 points wide.

        # Generally, xScale is going to be None at this point.
        if self.xScale == None:
            # I don't know how to adjust the xScale analytically, so I
            # will do it incrementally.  Start by calculating an xScale
            # assuming that taxon names on the right have no width.  We
            # want a bit (xOrigin) of space on both the left
            # and the right.
            self.xScale = (self.width - (2 * self.xOrigin)) / self.maxX

            # Find the contribution of the root.name to the width.
            if self.tree.root.name:
                # We have rootNameLength from above
                thisPointsPerLetter = (float(self.labelTextSize) / float(self.textSize)) * self.pointsPerLetter
                self.rootNameContribution =  (rootNameLength * thisPointsPerLetter) + self.nameOffset
            else:
                self.rootNameContribution = 0.0
            #print "self.rootNameContribution = %s" % self.rootNameContribution

            # Check if the self.width is wide enough, even if
            # self.xScale is zero.  We need enough room for the
            # rootName and the longest name.
            for i in nameLengths:
                theSum = self.rootNameContribution + (i[1] * self.pointsPerLetter) + self.nameOffset
                if theSum >= self.width - (2.0 * self.xOrigin):
                    if autoIncreaseWidth:
                        self.width =  theSum + (2.0 * self.xOrigin) + 1
                    else:
                        gm.append("The picture is not wide enough to fit the node names.")
                        gm.append("Increase the width.")
                        raise Glitch, gm

            #print "before: xScale is %f" % self.xScale
            makeAdjustments = 1
            while makeAdjustments:
                adjusted = 0
                for i in nameLengths:
                    theSum = self.rootNameContribution + (i[0] * self.xScale) + (i[1] * self.pointsPerLetter) + self.nameOffset
                    if theSum > (self.width - (2.0 * self.xOrigin)):
                        self.xScale = self.xScale * 0.99
                        #print "xScale now %f" % self.xScale
                        adjusted = 1
                if not adjusted:
                    makeAdjustments = 0

            #print "after: xScale is %f" % self.xScale
            #import sys; sys.exit()
            #print "xScale is %f" % self.xScale
        else:
            #print gm[0]
            #print "xScale is set.  Fix me!"
            #import sys; sys.exit()
            pass



    def vectorString(self):

        self.yOrigin = -1.0 * (self.yScale / 2.0)

        stringList = []

        bboxX = self.width
        bboxY = (self.maxY * self.yScale)

        if self.svg:
            self.yOrigin = bboxY + (self.yScale / 2.0)
            self.yOrigin += 3 # fudge
        #print "xOrigin = %s, yOrigin = %s" % (self.xOrigin, self.yOrigin)

        if self.svg:
            stringList.append('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN"\n')
            stringList.append(' "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n')
            stringList.append('<svg xmlns="http://www.w3.org/2000/svg"\n')
            stringList.append(' xmlns:xlink="http://www.w3.org/1999/xlink"\n')
            stringList.append(' xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"')
            stringList.append(' width="%ipx" height="%ipx">\n' % (bboxX, bboxY))
            stringList.append('\n<title>An SVG tree from p4</title>\n')
            stringList.append('<defs>\n')
            stringList.append('<style type="text/css"><![CDATA[line {stroke:blue; stroke-width:5; stroke-linecap:round;}]]></style>\n')
            stringList.append('</defs>\n')
            stringList.append('<g inkscape:label="Layer 1" inkscape:groupmode="layer" id="layer1">\n')

        else:
            stringList.append("%!PS-Adobe-3.0 EPSF-3.0\n")
            stringList.append('%%%%BoundingBox: 0  0 %i %i\n' % (bboxX, bboxY))

        if not self.svg:
            stringList.append('\n 1 setlinecap\n')
            stringList.append(' 1 setlinejoin\n\n')
            #stringList.append(' 0.4 setlinewidth\n')
            stringList.append(' 1 setlinewidth\n')

            stringList.append('\n% Draw horizontal lines.\n')
        for n in self.tree.iterNodesNoRoot():
            if self.svg:
                #if ((n.parent.xPos - n.xPos) * self.xScale) > 1.0:
                stringList.append('<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f"/>\n' % (
                        self.rootNameContribution + self.xOrigin + (n.parent.xPos * self.xScale), self.yOrigin - (n.yPos * self.yScale),
                        self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale), self.yOrigin - (n.yPos * self.yScale)))
            else:
                stringList.append('%f %f moveto \n' % \
                        (self.rootNameContribution + self.xOrigin + (n.parent.xPos * self.xScale), self.yOrigin + (n.yPos * self.yScale)))
                stringList.append('%f %f lineto stroke\n' % \
                        (self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale), self.yOrigin + (n.yPos * self.yScale)))
        if not self.svg:
            stringList.append('\n% Draw vertical lines.\n')
        for n in self.tree.iterNodes():
            if not n.isLeaf:
                if 0:
                    if n == self.tree.root:
                        print "root leftChild = %i, rightmostChild = %i" % \
                              (n.leftChild.nodeNum, n.rightmostChild().nodeNum)
                        print "line from %f %f" % \
                              ((n.xPos), (n.leftChild.yPos))
                        print "       to %f %f" % \
                              ((n.xPos), (n.rightmostChild().yPos))
                    else:
                        print "nodeNum %i, rootNameContribution = %s, xOrigin= %s, n.xPos = %s, self.xScale=%s," % (
                            n.nodeNum, self.rootNameContribution, self.xOrigin, n.xPos, self.xScale)
                        print "self.yOrigin=%s, n.leftChild.yPos=%s, self.yScale=%s" % (
                            self.yOrigin, n.leftChild.yPos, self.yScale)
                if self.svg:
                    stringList.append('<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f"/>\n' % (
                        self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale), self.yOrigin - (n.leftChild.yPos * self.yScale),
                        self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale), self.yOrigin - (n.rightmostChild().yPos * self.yScale)))
                else:
                    stringList.append('%f %f moveto \n' % \
                            (self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale), self.yOrigin + (n.leftChild.yPos * self.yScale)))
                    stringList.append('%f %f lineto stroke\n' % \
                            (self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale), self.yOrigin + (n.rightmostChild().yPos * self.yScale)))

        # Terminal node names
        if not self.svg:
            stringList.append('\n% Node names of (non-root) terminal nodes.\n')
            #stringList.append('\n/Helvetica-Oblique findfont %i scalefont setfont\n' % textSize)
            #stringList.append('/Times-Italic findfont %i scalefont setfont\n' % self.textSize)
            stringList.append('/Times-Roman findfont %i scalefont setfont\n' % self.textSize)
        for n in self.tree.iterLeavesNoRoot():
            if n.name:
                if self.svg:
                    stringList.append('<text x="%.0f" y="%.0f">%s</text>\n' % (
                        self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale) + self.nameOffset,
                        self.yOrigin - (n.yPos * self.yScale) + self.nameDrop, n.name))
                else:
                    stringList.append('%f %f moveto \n' % \
                                      (self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale) + self.nameOffset,
                                       self.yOrigin + (n.yPos * self.yScale) - self.nameDrop))
                    stringList.append('(%s) show\n' % n.name)

        # Put on root.name
        if self.tree.root.name:
            if not self.svg:
                stringList.append('\n% Root name.  This may need adjustment.\n')
            n = self.tree.root
            # If it is a terminal root, then we want textSize.  If it is
            # an internal node, we want labelTextSize.  Also, the position is
            # slightly different.
            if n.isLeaf:
                if self.svg:
                    stringList.append('<text x="%.0f" y="%.0f" text-anchor="end">%s</text>\n' % (
                        self.rootNameContribution + self.xOrigin, self.yOrigin - (n.leftChild.yPos * self.yScale) + self.nameDrop, n.name))
                else:
                    stringList.append('%f %f moveto \n' % \
                                      (self.rootNameContribution + self.xOrigin, self.yOrigin + (n.leftChild.yPos * self.yScale) - self.nameDrop))
                    # we want textSize here, but that is currently what it is, so no need to specify
                    stringList.append('%% Terminal root; should already be textSize, ie %s scalefont.\n' % self.textSize)
            else:
                if self.svg:
                    stringList.append('<text x="%.0f" y="%.0f" text-anchor="end">%s</text>\n' % (
                        self.rootNameContribution + self.xOrigin, self.yOrigin - (n.yPos * self.yScale) + self.nameDrop, n.name))
                else:
                    stringList.append('%f %f moveto \n' %  (self.rootNameContribution + self.xOrigin, self.yOrigin + (n.yPos * self.yScale) - self.nameDrop))
                    stringList.append('/Times-Roman findfont %i scalefont setfont\n' % self.labelTextSize)
            # move to the left by the width of the name, and a bit more.
            # nameOffset is textSize/3.0, but it is a bit small, so use
            # textSize/2.0
            if self.svg:
                pass
            else:
                stringList.append('(%s) stringwidth pop -1.0 mul 0 rmoveto\n' % n.name)
                stringList.append('-%s 0 rmoveto\n' % (self.textSize * 0.5))
                stringList.append('(%s) show\n' % n.name)


        # Internal node names, and branch labels for all nodes.
        if not self.svg:
            stringList.append('\n% Node names of internal (non-root) nodes.\n')
            if self.putInternalNodeNamesOnBranches:
                stringList.append('% (Here, internal node names are placed on branches).\n')
            else:
                stringList.append('% (Here, internal node names are placed on nodes).\n')
            stringList.append('/Times-Roman findfont %i scalefont setfont\n' % self.labelTextSize)
        for n in self.tree.iterNodesNoRoot():
            if not n.isLeaf and n.name:
                if self.putInternalNodeNamesOnBranches:
                    if self.svg:
                        pass
                    else:
                        stringList.append('%f %f moveto \n' % \
                            (self.rootNameContribution + self.xOrigin + ((n.xPos + n.parent.xPos) * 0.5 * self.xScale),
                             self.yOrigin + (n.yPos * self.yScale) + 2 ))
                        #if n.name[0] == '\'' and n.name[-1] == '\'':  # Is this needed?
                        #    theName = n.name[1:-1]
                        #else:
                        #    theName = n.name
                        theName = n.name
                        stringList.append('(%s) stringwidth pop -0.5 mul 0 rmoveto\n' % theName)
                        stringList.append('(%s) show\n' % theName)
                if not self.putInternalNodeNamesOnBranches:
                    if self.svg:
                        pass
                    else:
                        stringList.append('%f %f moveto \n' % \
                                          (self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale) + (self.labelTextSize / 3.0),
                                           self.yOrigin + (n.yPos * self.yScale) - (self.labelTextSize /2.0) + (self.labelTextSize / 3.0)))
                        #if n.name[0] == '\'' and n.name[-1] == '\'':
                        #    theName = n.name[1:-1]
                        #else:
                        #    theName = n.name
                        stringList.append('(%s) show\n' % n.name)

        if self.putInternalNodeNamesOnBranches:
            # If node names are on branches, then branch names are ignored.
            pass
        else:
            if not self.svg:
                stringList.append('\n% Branch names.\n')
            for n in self.tree.iterNodesNoRoot():
                if hasattr(n.br, 'name') and n.br.name:
                    if self.svg:
                        pass
                    else:
                        stringList.append('%f %f moveto \n' % \
                            (self.rootNameContribution + self.xOrigin + ((n.xPos + n.parent.xPos) * 0.5 * self.xScale),
                             self.yOrigin + (n.yPos * self.yScale) + 2 ))
                        if n.br.name[0] == '\'' and n.br.name[-1] == '\'':
                            theName = n.br.name[1:-1]
                        else:
                            theName = n.br.name
                        stringList.append('(%s) stringwidth pop -0.5 mul 0 rmoveto\n' % theName)
                        stringList.append('(%s) show\n' % theName)

        if self.svg:
            stringList.append('</g>\n')
            stringList.append('<use xlink:href="#tree1" transform="scale(1.0)"/>\n')
            stringList.append('</svg>\n\n')
        else:
            stringList.append('\nshowpage \n\n')
        return string.join(stringList, '')


    def textString(self, returnAsList=False):

        self.width = int(round(self.width))
        self.yOrigin = 0.0
        if 0:
            print "xOrigin = %s, yOrigin = %s" % (self.xOrigin, self.yOrigin)
            print "xScale = %s, yScale = %s" % (self.xScale, self.yScale)

        nRows = 0
        for n in self.tree.iterNodes():
            n.yPos = int(round(2.0 * n.yPos)) - 2
            if n.yPos > nRows:
                nRows = n.yPos
        nRows += 1

        if 0:
            print "index   isLeaf        xPos         yPos"
            for n in self.tree.iterNodes():
                print "%4i" % n.nodeNum,
                print "%8i" % n.isLeaf,
                print "%12s" % n.xPos,
                print "%12s" % n.yPos
            print "%s rows" % nRows

        # Make a field of blanks on which to draw
        stringList = []
        for i in range(nRows):
            stringList.append([' '] * self.width)

        # Draw the horizontal lines
        #print 'textDrawModelThing = %s' % self.textDrawModelThing
        for n in self.tree.iterNodesNoRoot():
            # Choose a symbol
            if hasattr(n.br, 'textDrawSymbol') and n.br.textDrawSymbol:
                theSymbol = n.br.textDrawSymbol
            else:
                theSymbol = '-'
            if self.partNum > -1:
                tp = self.tree.model.parts[self.partNum]
                if self.textDrawModelThing == var.TEXTDRAW_COMP:
                    if n.parts[self.partNum].compNum >= 0 and tp.comps[n.parts[self.partNum].compNum].symbol:
                        theSymbol = tp.comps[n.parts[self.partNum].compNum].symbol
                elif self.textDrawModelThing == var.TEXTDRAW_RMATRIX:
                    if n.br.parts[self.partNum].rMatrixNum >= 0 and tp.rMatrices[n.br.parts[self.partNum].rMatrixNum].symbol:
                        theSymbol = tp.rMatrices[n.br.parts[self.partNum].rMatrixNum].symbol
                elif self.textDrawModelThing == var.TEXTDRAW_GDASRV:
                    if n.br.parts[self.partNum].gdasrvNum >= 0 and tp.gdasrvs[n.br.parts[self.partNum].gdasrvNum].symbol:
                        theSymbol = tp.gdasrvs[n.br.parts[self.partNum].gdasrvNum].symbol
                else:
                    #raise Glitch, "TreePicture.textString().  No good textDrawModelThing.  Fix me"
                    pass

            start = int(round(self.rootNameContribution + self.xOrigin + (n.parent.xPos * self.xScale)))
            end = int(round(self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale)))
            #print "start = %s, end = %s, width = %s" % (start, end, self.width)
            for i in range(start + 1, end):
                stringList[n.yPos][i] = theSymbol
            stringList[n.yPos][start] = '+'
            # When the root is a single monofurcating stem, replace the '+' with theSymbol
            if n.parent == self.tree.root and self.tree.root.getNChildren() == 1:
                stringList[n.yPos][start] = theSymbol

        # Draw the vertical lines
        for n in self.tree.iterNodes():
            if not n.isLeaf:
                theXPos = int(round(self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale)))
                top = n.leftChild.yPos
                bot = n.rightmostChild().yPos
                for i in range(bot + 1, top):
                    stringList[i][theXPos] = '|'

        # Put on terminal node names
        if 1:
            for n in self.tree.iterLeavesNoRoot():
                #if self.textShowNodeNums:
                #    theName = '%s:%s' % (n.nodeNum, n.name)
                #else:
                #    theName = '%s' % n.name
                if self.textShowNodeNums:
                    theName = '%s:' % n.nodeNum
                else:
                    theName = ''
                theName += ' ' * self.nameOffset + '%s' % n.name
                #print theName
                start = int(round(self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale)))
                #start = start + self.nameOffset
                j = 0
                if theName:
                    for i in range(start, start + len(theName)):
                        stringList[n.yPos][i] = theName[j]
                        j += 1

        # Put on internal node numbers and names
        #print "showInternalNodeNames = %s" % self.showInternalNodeNames
        #print "textShowNodeNums = %s" % self.textShowNodeNums
        for n in self.tree.iterInternalsNoRoot():
            theName = ''
            if self.textShowNodeNums:
                theName = `n.nodeNum`
            #print "x %s" % theName
            if self.showInternalNodeNames and n.name:
                if self.textShowNodeNums:
                    theName += ':%s' % n.name
                    #theName += ':' + ' ' * self.nameOffset + '%s' % n.name
                else:
                    theName += '%s' % n.name
                    #theName += ' ' * self.nameOffset + '%s' % n.name
            #print "y %s" % theName
            start = int(round(self.rootNameContribution + self.xOrigin + (n.xPos * self.xScale)))
            if self.textShowNodeNums:
                pass # We want to over-write the vertical line with the node number
            else:
                start = start + 1 # So we do not over-write a vertical line
            #start = start + self.nameOffset
            j = 0
            for i in range(start, start + len(theName)):
                stringList[n.yPos][i] = theName[j]
                j += 1

        # Put on root node name and number
        n = self.tree.root
        theName = None
        if not self.showInternalNodeNames and not n.isLeaf:
            if self.textShowNodeNums:
                theName = '%s' % n.nodeNum
            else:
                pass
        else:
            if self.textShowNodeNums:
                if n.name:
                    theName = '%s:%s' % (n.name, n.nodeNum)
                else:
                    theName = '%s' % n.nodeNum
            else:
                if n.name:
                    theName = n.name
        if theName:
            start = int(round(self.xOrigin))
            j = 0
            for i in range(start, start + len(theName)):
                stringList[n.yPos][i] = theName[j]
                j += 1

        # Make strings from the list elements
        stringList.reverse()
        for i in range(len(stringList)):
            #stringList[i] = string.join(stringList[i] + ['|'], '')
            stringList[i] = string.join(stringList[i], '')
            stringList[i] = string.rstrip(stringList[i])

        # Put a model key on the end
        if self.partNum > -1:
            tp = self.tree.model.parts[self.partNum]
            if self.textDrawModelThing == var.TEXTDRAW_COMP:
                stringList.append('The tree above shows part %s comps' % self.partNum)
                if self.tree.model.parts[self.partNum].nComps:
                    for i in range(self.tree.model.parts[self.partNum].nComps):
                        stringList.append('    %-3i %s' % (i, self.tree.model.parts[self.partNum].comps[i].symbol))
                    rt = self.tree.root
                    stringList.append('    root (node %s) has comp %s, symbol %s' % (
                        rt.nodeNum,
                        tp.comps[rt.parts[self.partNum].compNum].num,
                        tp.comps[rt.parts[self.partNum].compNum].symbol))
                else:
                    stringList.append('    No comps defined for part %s.' % self.partNum)
            elif self.textDrawModelThing == var.TEXTDRAW_RMATRIX:
                stringList.append('The tree above shows part %s rMatrices' % self.partNum)
                if self.tree.model.parts[self.partNum].nRMatrices:
                    for i in range(self.tree.model.parts[self.partNum].nRMatrices):
                        stringList.append('    %-3i %s' % (i, self.tree.model.parts[self.partNum].rMatrices[i].symbol))
                else:
                    stringList.append('    No rMatrices defined for part %s.' % self.partNum)
            elif self.textDrawModelThing == var.TEXTDRAW_GDASRV:
                stringList.append('The tree above shows part %s gdasrvs' % self.partNum)
                if self.tree.model.parts[self.partNum].nGdasrvs:
                    for i in range(self.tree.model.parts[self.partNum].nGdasrvs):
                        stringList.append('    %-3i %s' % (i, self.tree.model.parts[self.partNum].gdasrvs[i].symbol))
                else:
                    stringList.append("    No gdasrv's defined for part %s." % self.partNum)

        if returnAsList:
            return stringList
        else:
            s = string.join([''] + stringList, '\n')
            return s


