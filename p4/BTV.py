import time,os,sys,random,string,math
from p4.Tree import Tree
from Glitch import Glitch
try:
    from Tkinter import *
except ImportError:
    raise Glitch, "TV and BTV need Tkinter, and it does not seem to be installed."


def randomColour():
    # Colours are #xxyyzz, where each pair is a hex number.  Colours
    # with e and f as the first digit in all 3 pairs are sometimes too
    # faint.
    while 1:
        c = ['#']
        firsts = []
        for i in range(3):
            theChoice = random.choice(string.hexdigits)
            c.append(theChoice)
            firsts.append(theChoice)
            c.append(random.choice(string.hexdigits))
        lightCount = firsts.count('f') + firsts.count('e')
        if lightCount < 3:
            break
    return ''.join(c)


class CanvasA(Canvas):
    def __init__(self, *args, **kwargs):
        Canvas.__init__(self, *args, **kwargs)
        self.btv_tree = None
        self.longestNameLen = 0
        self.nameOffset = 4
        self.pointsPerLetter = 8
        self.lineWid = 3
        self.nameStuff = 0
        self.lineUpLeaves = True

    def btv_findLongestTaxName(self):
        longest = 0
        for n in self.btv_tree.iterLeavesNoRoot():
            if len(n.name) > longest:
                longest = len(n.name)
        self.longestNameLen = longest
        self.nameStuff = self.nameOffset + (self.longestNameLen * self.pointsPerLetter)
        

    def draw2(self):
        self.delete(ALL)
        if 0:
            bord = 5
            self.create_rectangle(bord, bord, self.curWinWidth - bord, self.curWinHeight - bord)

        
        #print self.box
        if 1:

            if self.box[2] != self.box[0]:
                self.xScale = (self.curWinWidth - (20. + self.nameStuff)) / (self.box[2] - self.box[0])
            else:
                self.xScale = (self.curWinWidth - (20. + self.nameStuff)) /0.001
            self.xOrig = 10. - (self.box[0] * self.xScale)

            if self.box[3] != self.box[1]:
                self.yScale = (self.curWinHeight - 20.) / (self.box[3] - self.box[1])
            else:
                self.yScale = (self.curWinHeight - 20.) / 2.
            self.yOrig = 10. - (self.box[1] * self.yScale)

            myCapStyle='round'
            for n in self.btv_tree.root.iterPreOrder():
                if n != self.btv_tree.root:
                    if n.yPos < self.box[3] and n.yPos > self.box[1]:
                        if self.lineUpLeaves:
                            nParentXPos = n.parent.xPosL
                            nXPos = n.xPosL
                        else:
                            nParentXPos = n.parent.xPos1
                            nXPos = n.xPos1
                        if n.isLeaf:
                            self.create_line((nParentXPos * self.xScale) + self.xOrig,
                                             (n.yPos * self.yScale) + self.yOrig,
                                             (nXPos * self.xScale) + self.xOrig,
                                             (n.yPos * self.yScale) + self.yOrig,
                                             fill=n.br.color, width=self.lineWid, capstyle=myCapStyle)
                            self.create_text((nXPos * self.xScale) + self.xOrig + self.nameOffset,
                                             (n.yPos * self.yScale) + self.yOrig,
                                             anchor='w',
                                             text=n.name, fill=n.nameColor)
                        else:
                            self.create_line((nParentXPos * self.xScale) + self.xOrig,
                                             (n.yPos * self.yScale) + self.yOrig,
                                             (nXPos * self.xScale) + self.xOrig,
                                             (n.yPos * self.yScale) + self.yOrig,
                                             fill=n.br.color, width=self.lineWid, capstyle=myCapStyle)
                            if n.name:
                                if 0:
                                    myAnchor = 'w'
                                    for ch in n.iterChildren():
                                        theDiff = math.fabs((ch.yPos - n.yPos) * self.yScale)
                                        #print "internal %s, diff is %f" % (n.name, theDiff)
                                        if theDiff < 7.5:
                                            if ch.yPos >= n.yPos:
                                                myAnchor = 'sw'
                                            else:
                                                myAnchor = 'nw'
                                    self.create_text((nXPos * self.xScale) + self.xOrig + self.nameOffset,
                                                 (n.yPos * self.yScale) + self.yOrig,
                                                 anchor=myAnchor,
                                                 text=n.name, fill='darkgreen')
                                else:
                                    self.create_text((nXPos * self.xScale) + self.xOrig - self.nameOffset,
                                                 (n.yPos * self.yScale) + self.yOrig,
                                                 anchor='se',
                                                 text=n.name, fill='darkgreen')

            # vertical lines
            for n in self.btv_tree.root.iterPreOrder():
                if not n.isLeaf:
                    if self.lineUpLeaves:
                        nXPos = n.xPosL
                    else:
                        nXPos = n.xPos1
                    self.create_line((nXPos * self.xScale) + self.xOrig,
                                     (n.leftChild.yPos * self.yScale) + self.yOrig,
                                     (nXPos * self.xScale) + self.xOrig,
                                     (n.rightmostChild().yPos * self.yScale) + self.yOrig,
                                     width=self.lineWid, capstyle='round')

class CanvasB(Canvas):
    def __init__(self, *args, **kwargs):
        Canvas.__init__(self, *args, **kwargs)
        self.btv_tree = None
        
        self.rect = None
        self.startx = 0.0
        self.starty = -1.5
        self.endx = 1.0
        self.endy = 20.5
        self.lineWid = 1
        self.xOrig = 10
        self.yOrig = 10
        self.lineUpLeaves = True
        self.btv = None
        self.setBinds()

    def setBinds(self):
        self.bind('<Button-1>', self.mouseDown)
        self.bind('<Button1-Motion>', self.mouseMotion)
        self.bind('<Button1-ButtonRelease>', self.mouseUp)
        self.bind_all('<Key>', self.btv_key)

    def draw1(self):

        self.delete(ALL)

        if 0:
            bord = 5
            self.create_rectangle(bord, bord,
                                  self.curWinWidth - bord,
                                  self.curWinHeight - bord, outline='red')

        self.xScale = self.curWinWidth - 20

        self.yScale = (self.curWinHeight - 20.) / self.btv_tree.nLeaves

        myCapStyle='round'
        #for n in self.btv_tree.root.iterPreOrder():
        for n in self.btv_tree.nodes:
            if n != self.btv_tree.root:
                if self.lineUpLeaves:
                    nParentXPos = n.parent.xPosL
                    nXPos = n.xPosL
                else:
                    nParentXPos = n.parent.xPos1
                    nXPos = n.xPos1
                if n.isLeaf:
                    self.create_line((nParentXPos * self.xScale) + self.xOrig,
                                     (n.yPos * self.yScale) + self.yOrig,
                                     (nXPos * self.xScale) + self.xOrig,
                                     (n.yPos * self.yScale) + self.yOrig,
                                     fill=n.br.color, width=self.lineWid, capstyle=myCapStyle)
                else:
                    self.create_line((nParentXPos * self.xScale) + self.xOrig,
                                     (n.yPos * self.yScale) + self.yOrig,
                                     (nXPos * self.xScale) + self.xOrig,
                                     (n.yPos * self.yScale) + self.yOrig,
                                     fill=n.br.color, width=self.lineWid, capstyle=myCapStyle)

        
        if self.btv.setsVar:
            theSet = self.btv.setsVar.get()
            #print "Got theSet %i" % theSet
            if theSet:
                # Re-draw the red lines, cuz they may have been covered in black.
                for n in self.btv_tree.nodes:
                    if n.isLeaf and n.br.color == 'red':
                        if self.lineUpLeaves:
                            nParentXPos = n.parent.xPosL
                            nXPos = n.xPosL
                        else:
                            nParentXPos = n.parent.xPos1
                            nXPos = n.xPos1
                        self.create_line((nParentXPos * self.xScale) + self.xOrig,
                                     (n.yPos * self.yScale) + self.yOrig,
                                     (nXPos * self.xScale) + self.xOrig,
                                     (n.yPos * self.yScale) + self.yOrig,
                                     fill=n.br.color, width=self.lineWid, capstyle=myCapStyle)
                        
                
        # vertical lines
        for n in self.btv_tree.root.iterPreOrder():
            if not n.isLeaf:
                if self.lineUpLeaves:
                    nXPos = n.xPosL
                else:
                    nXPos = n.xPos1
                self.create_line((nXPos * self.xScale) + self.xOrig,
                                 (n.leftChild.yPos * self.yScale) + self.yOrig,
                                 (nXPos * self.xScale) + self.xOrig,
                                 (n.rightmostChild().yPos * self.yScale) + self.yOrig,
                                 width=self.lineWid, capstyle='round')


    def mouseDown(self, event):
        if self.rect:
            self.delete(self.rect)
        self.startx = event.x
        self.starty = event.y
        self.bind_all('<Key>', self.btv_key)

    def mouseMotion(self, event):
        if 0:
            #self.aspectRatio == (self.endy - self.starty) / (self.endx - self.startx)
            x = event.x
            self.endy = event.y
            self.endx = ((self.endy - self.starty) / self.aspectRatio) + self.startx


            if self.endy > self.curWinHeight:
                self.endy = self.curWinHeight
                self.endx = ((self.endy - self.starty) / self.aspectRatio) + self.startx
            if self.endx > self.curWinWidth:
                self.endx = self.curWinWidth
                self.endy = (self.aspectRatio * (self.endx - self.startx)) + self.starty
            
        self.endx = event.x
        self.endy = event.y

        if self.endy > self.curWinHeight:
            self.endy = self.curWinHeight
        if self.endx > self.curWinWidth:
            self.endx = self.curWinWidth

        if (self.startx != event.x)  and (self.starty != event.y):
            if 0:
                if self.startx > self.endx:
                    temp = self.startx
                    self.startx = self.endx
                    self.endx = temp
                if self.starty > self.endy:
                    temp = self.starty
                    self.starty = self.endy
                    self.endy = temp
            self.drawRect()

    def drawFirstRect(self):
        self.startx  = (-0.012 * float(self.xScale)) + self.xOrig
        self.starty  = (-2.5 * float(self.yScale)) + self.yOrig
        self.endx  = (1.012 * float(self.xScale)) + self.xOrig

        theYNum = self.btv_tree.nLeaves /4
        if 23.5 < theYNum:
            theYNum = 23.5
        self.endy = (theYNum * float(self.yScale)) + self.yOrig
        #print "drawFirstRect here. ", self.startx, self.starty, self.endx, self.endy
        self.rect = self.create_rectangle(
            self.startx, self.starty, self.endx, self.endy, width=2, outline='red')
        self.update_idletasks()
        
    def drawRect(self):
        #print "drawRect here. ", self.startx, self.starty, self.endx, self.endy
        self.delete(self.rect)
        self.rect = self.create_rectangle(
            self.startx, self.starty, self.endx, self.endy, width=2, outline='red')
        # this flushes the output, making sure that
        # the rectangle makes it to the screen
        # before the next event is handled
        self.update_idletasks()

    def mouseUp(self, event):
        self.doBox()

    def doBox(self):
        self.box = [0.0] * 4
        smallestx = self.startx
        biggestx = self.endx
        if smallestx > biggestx:
            smallestx = self.endx
            biggestx = self.startx
        smallesty = self.starty
        biggesty = self.endy
        if smallesty > biggesty:
            smallesty = self.endy
            biggesty = self.starty
            
        self.box[0] = (smallestx - self.xOrig) / float(self.xScale)
        self.box[1] = (smallesty - self.yOrig) / float(self.yScale)
        self.box[2] = (biggestx - self.xOrig) / float(self.xScale)
        self.box[3] = (biggesty - self.yOrig) / float(self.yScale)
        self.btv_ca.box = self.box
        self.btv_ca.draw2()

    def btv_key(self, event):
        if event.keysym == 'Up':
            #print "Up"
            halfRectHeight = (self.endy - self.starty) / 2.
            self.endy -= halfRectHeight
            self.starty -= halfRectHeight
            self.drawRect()
            self.doBox()
            
        elif event.keysym == 'Right':
            #print "Right"
            pass
        elif event.keysym == 'Down':
            #print "Down"
            halfRectHeight = (self.endy - self.starty) / 2.
            self.endy += halfRectHeight
            self.starty += halfRectHeight
            self.drawRect()
            self.doBox()
        elif event.keysym == 'Left':
            #print "Left"
            pass

btvHelpText = """This is BTV, a Big Tree Viewer.

The full tree is shown in the panel on the right.  A selection from it
is shown in the panel on the left.  The selection is given by the red
rectangle in the right panel.  On startup, and whenever the window is
resized, a default selection is made at the top of the tree.  You can
make a new selection by mouse-dragging a new rectangle.

You can move the selection up and down using the up and down arrow
keys.


"""

class BTV(Canvas):
    """A big tree viewer.  Needs Tkinter.

    Instantiate it by passing it a p4 Tree object, eg

        BTV(myTree)

    """
    def __init__(self, tree):

        self.tk_root = Tk()
        self.tk_root.withdraw()
        master = Toplevel(self.tk_root)
        master.protocol("WM_DELETE_WINDOW", self.__close_help)
        Canvas.__init__(self, master)
        self.master.title("BTV")
        #master.resizable(1,1)
        self.closed = False

        #tree.lineUpLeaves(rootToLeaf=1.0)
        self.setsVar = None
        self.btv_tree = tree
        self.btv_ca = CanvasA(self, width=394, height=400, background='white')
        self.btv_ca.btv_tree = self.btv_tree
        self.btv_ca.btv_findLongestTaxName()
        self.btv_ca.pack(expand=1, fill=BOTH, side=LEFT)
        self.btv_cb = CanvasB(self, width=194, height=400, background='white')
        self.btv_cb.btv_tree = self.btv_tree
        self.btv_cb.btv_ca = self.btv_ca
        self.btv_ca.btv_cb = self.btv_cb
        self.btv_cb.btv = self
        self.btv_cb.pack(expand=1, fill=BOTH, side=RIGHT)
        # set xPos, which does not change
        biggestXPos = 0.0
        self.btv_tree.root.xPos = 0.0
        self.btv_tree.root.yPos = 0.0
        for n in self.btv_tree.iterPreOrder():
            if n != self.btv_tree.root:
                n.xPos = n.parent.xPos + n.br.len
                if n.xPos > biggestXPos:
                    biggestXPos = n.xPos
        self.btv_tree.root.xPos1 = 0.0
        for n in self.btv_tree.iterNodesNoRoot():
            n.xPos1 = n.xPos / biggestXPos

        self.btv_tree.lineUpLeaves(overWriteBrLens=False)
        self.btv_tree.root.xPosL = 0.0
        self.btv_tree.root.yPosL = 0.0
        for n in self.btv_tree.iterPreOrder():
            if n != self.btv_tree.root:
                n.xPosL = n.parent.xPosL + n.br.lenL

        # set yPos
        counter = 0
        for n in self.btv_tree.root.iterPostOrder():
            if n.isLeaf:
                if n == self.btv_tree.root:
                    pass
                else:
                    n.yPos = float(counter)
                    counter += 1
            else:
                n.yPos = (n.leftChild.yPos + n.rightmostChild().yPos) / 2.0

        # Count leaves in the tree
        self.btv_tree.nLeaves = 0
        for n in self.btv_tree.nodes:
            if n.isLeaf:
                self.btv_tree.nLeaves += 1
                n.nameColor = 'navy'
            if n.br:
                n.br.color1 = randomColour()
                n.br.color = n.br.color1

        self.configure(height=400,width=600, background='white')
        self.pack(expand=True, fill=BOTH)

        # These lines following don't work on my linux box.
        #self.btv_ca.curWinHeight = self.btv_ca.winfo_height()
        #self.btv_ca.curWinWidth = self.btv_ca.winfo_width()
        #self.btv_cb.curWinHeight = self.btv_cb.winfo_height()
        #self.btv_cb.curWinWidth = self.btv_cb.winfo_width()

        self.btv_ca.curWinHeight = 394
        self.btv_ca.curWinWidth = 394
        self.btv_cb.curWinHeight = 394
        self.btv_cb.curWinWidth = 194

        #self.counter = 0
        self.bind('<Configure>', self.btv_configure) # notices window size changes
        self.bind("<Activate>", self._onActivate)

        self.btv_menu = Menu(self)
        master.config(menu=self.btv_menu)
        self.btv_helpMenu = Menu(self.btv_menu)
        self.btv_menu.add_cascade(label="Help", menu=self.btv_helpMenu)
        self.btv_helpMenu.add_command(label="Help", command=self.btv_doHelp)

        self.btv_gramMenu = Menu(self.btv_menu)
        self.btv_menu.add_cascade(label="Gram", menu=self.btv_gramMenu)
        self.btv_gramMenu.add_command(label="Cladogram", command=self.btv_doLineUpLeaves)
        self.btv_gramMenu.add_command(label="Phylogram", command=self.btv_doOriginalBranchScale)

        if self.btv_tree.nexusSets and self.btv_tree.nexusSets.taxSets:
            self.setsVar = IntVar(master=master)
            self.btv_setsMenu = Menu(self.btv_menu)
            self.btv_menu.add_cascade(label="TaxSets", menu=self.btv_setsMenu)
            self.btv_setsMenu.add_radiobutton(label="None", variable=self.setsVar, value=0,
                                              command=self.btv_doColorTaxSets)
            txSetIndx = 1
            for txSet in self.btv_tree.nexusSets.taxSets:
                self.btv_setsMenu.add_radiobutton(label=txSet.name, variable=self.setsVar,
                                                 value=txSetIndx, command=self.btv_doColorTaxSets)
                txSetIndx += 1
        #root.mainloop()

    def btv_doColorTaxSets(self):
        #print "BTV doColorTaxSets here!"
        theSet = self.setsVar.get()
        #print "sets is %s\n" % theSet
        if theSet == 0:
            for n in self.btv_tree.nodes:
                if n.isLeaf:
                    n.nameColor = 'navy'
                if n.br:
                    n.br.color = n.br.color1
        else:
            theTaxSet = self.btv_tree.nexusSets.taxSets[theSet - 1]
            for n in self.btv_tree.nodes:
                if n.isLeaf:
                    theIndx = self.btv_tree.taxNames.index(n.name)
                    if theTaxSet.mask[theIndx] == '1':
                        n.nameColor = 'red'
                        if n.br:
                            n.br.color = 'red'
                    else:
                        n.nameColor = 'navy'
                        if n.br:
                            n.br.color = 'black'
                else:
                    if n.br:
                        n.br.color = 'black'
                    
        self.btv_cb.draw1()
        self.btv_cb.drawFirstRect()
        self.btv_cb.doBox()

    def _onActivate(self, e):
        #print "BTV onActivate.  self=%s" % self
        self.btv_cb.setBinds()


    def __close_help(self):
        """Close the window"""
        self.closed = True
        self.master.destroy()
        self.tk_root.update()

    def btv_doHelp(self):
        tl = Toplevel()
        tl.label = Label(tl, text=btvHelpText)
        tl.label.pack()

    def btv_doLineUpLeaves(self):
        self.btv_ca.lineUpLeaves = True
        self.btv_cb.lineUpLeaves = True
        self.btv_cb.draw1()
        self.btv_cb.drawRect()
        self.btv_cb.doBox()
        
    def btv_doOriginalBranchScale(self):
        self.btv_ca.lineUpLeaves = False
        self.btv_cb.lineUpLeaves = False
        self.btv_cb.draw1()
        self.btv_cb.drawRect()
        self.btv_cb.doBox()
        

        
    def btv_configure(self, event):
        #print "btv_configure() called.  event.type=%s" % event.type
        #print dir(event)
        ew = event.width
        eh = event.height

        #print "btv_configure:  event.width=%s, event.height=%s, ca.width=%s, cb.width=%s" % (
        #    ew, eh, self.btv_ca.curWinWidth, self.btv_cb.curWinWidth)

        if 1:
            if (ew - 12) > (self.btv_ca.curWinWidth + self.btv_cb.curWinWidth):
                self.btv_ca.configure(width=(2 * (ew/3.)) - 6)
                self.btv_cb.configure(width=(ew/3.) - 6)
            elif (ew + 12) < (self.btv_ca.curWinWidth + self.btv_cb.curWinWidth):
                self.btv_ca.configure(width=(2 * (ew/3.)) - 6)
                self.btv_cb.configure(width=(ew/3.) - 6)
        if (eh - 6) > (self.btv_ca.winfo_height()):
            self.btv_ca.configure(height=eh - 6)
            self.btv_cb.configure(height=eh - 6)

        self.update()

        self.btv_ca.curWinHeight = self.btv_ca.winfo_height()
        self.btv_ca.curWinWidth = self.btv_ca.winfo_width()
        self.btv_cb.curWinHeight = self.btv_cb.winfo_height()
        self.btv_cb.curWinWidth = self.btv_cb.winfo_width()

        self.btv_cb.draw1()
        self.btv_cb.drawFirstRect()
        self.btv_cb.doBox()
        



tvHelpText = """If your tree is too big,
you may want to try the Big Tree Viewer.
"""


class TV(Canvas):
    """A tree viewer.  Needs Tkinter.

    Instantiate it by passing it a p4 Tree object, eg

        TV(myTree)

    """

    def __init__(self, tree, title='TV'):

        self.tk_root = Tk()
        self.tk_root.withdraw()
        master = Toplevel(self.tk_root)
        master.protocol("WM_DELETE_WINDOW", self.__close_help)
        Canvas.__init__(self, master) # returns None
        # self.master is now an instance of Tkinter.Toplevel
        self.master.title(title)
        #master.resizable(1,1)
        self.closed = False

        self.longestNameLen = 0
        self.nameOffset = 4
        self.pointsPerLetter = 8
        self.lineWid = 3
        self.nameStuff = 0
        self.lineUpLeaves = True
        self.xOrig = 10
        self.yOrig = 10

        self.tv_tree = tree

        # Count leaves in the tree
        self.tv_tree.nLeaves = 0
        for n in self.tv_tree.nodes:
            if n.isLeaf:
                self.tv_tree.nLeaves += 1
                n.nameColor = 'navy'
            if n.br:
                n.br.color = 'black'


        # Find the longest name
        longest = 0
        for n in self.tv_tree.iterLeavesNoRoot():
            if len(n.name) > longest:
                longest = len(n.name)
        self.longestNameLen = longest
        self.nameStuff = self.nameOffset + (self.longestNameLen * self.pointsPerLetter)

        # set xPos, which does not change
        biggestXPos = 0.0
        self.tv_tree.root.xPos = 0.0
        self.tv_tree.root.yPos = 0.0
        for n in self.tv_tree.iterPreOrder():
            if n != self.tv_tree.root:
                n.xPos = n.parent.xPos + n.br.len
                if n.xPos > biggestXPos:
                    biggestXPos = n.xPos
        self.tv_tree.root.xPos1 = 0.0
        for n in self.tv_tree.iterNodesNoRoot():
            n.xPos1 = n.xPos / biggestXPos

        self.tv_tree.lineUpLeaves(overWriteBrLens=False)
        self.tv_tree.root.xPosL = 0.0
        self.tv_tree.root.yPosL = 0.0
        for n in self.tv_tree.iterPreOrder():
            if n != self.tv_tree.root:
                n.xPosL = n.parent.xPosL + n.br.lenL

        # set yPos
        counter = 0
        for n in self.tv_tree.root.iterPostOrder():
            if n.isLeaf:
                if n == self.tv_tree.root:
                    pass
                else:
                    n.yPos = float(counter)
                    counter += 1
            else:
                n.yPos = (n.leftChild.yPos + n.rightmostChild().yPos) / 2.0

        #self.pack()
        self.configure(height=400,width=400, background='white')
        self.pack(expand=True, fill=BOTH)

        # These lines following don't work on my linux box.
        #self.tv_ca.curWinHeight = self.tv_ca.winfo_height()
        #self.tv_ca.curWinWidth = self.tv_ca.winfo_width()
        #self.tv_cb.curWinHeight = self.tv_cb.winfo_height()
        #self.tv_cb.curWinWidth = self.tv_cb.winfo_width()

        self.curWinHeight = 394
        self.curWinWidth = 394

        #self.counter = 0
        self.bind('<Configure>', self.tv_configure) # notices window size changes

        self.tv_menu = Menu(self)
        #print "self.tv_menu is %s" % self.tv_menu
        master.config(menu=self.tv_menu)
        self.tv_helpMenu = Menu(self.tv_menu)
        self.tv_menu.add_cascade(label="Help", menu=self.tv_helpMenu)
        self.tv_helpMenu.add_command(label="Help", command=self.tv_doHelp)

        self.tv_gramMenu = Menu(self.tv_menu)
        self.tv_menu.add_cascade(label="Gram", menu=self.tv_gramMenu)
        self.tv_gramMenu.add_command(label="Cladogram", command=self.tv_doLineUpLeaves)
        self.tv_gramMenu.add_command(label="Phylogram", command=self.tv_doOriginalBranchScale)

        if self.tv_tree.nexusSets and self.tv_tree.nexusSets.taxSets:
            self.setsVar = IntVar(master=master)
            self.tv_setsMenu = Menu(self.tv_menu)
            self.tv_menu.add_cascade(label="TaxSets", menu=self.tv_setsMenu)
            self.tv_setsMenu.add_radiobutton(label="None", variable=self.setsVar, value=0, command=self.tv_doColorTaxSets)
            txSetIndx = 1
            for txSet in self.tv_tree.nexusSets.taxSets:
                self.tv_setsMenu.add_radiobutton(label=txSet.name, variable=self.setsVar,
                                                 value=txSetIndx, command=self.tv_doColorTaxSets)
                txSetIndx += 1

        self.bind("<1>", self.mouseDown)
        self.bind("<B1-Motion>", self.mouseMove)
        

    def tv_doColorTaxSets(self):
        #print "doColorTaxSets here!"
        theSet = self.setsVar.get()
        #print "sets is %s\n" % theSet
        if theSet == 0:
            for n in self.tv_tree.iterLeavesNoRoot():
                n.nameColor = 'navy'
        else:
            theTaxSet = self.tv_tree.nexusSets.taxSets[theSet - 1]
            for n in self.tv_tree.iterLeavesNoRoot():
                theIndx = self.tv_tree.taxNames.index(n.name)
                if theTaxSet.mask[theIndx] == '1':
                    n.nameColor = 'red'
                else:
                    n.nameColor = 'navy'
        self.draw1()
        
        

    def __close_help(self):
        """Close the window"""
        self.closed = True
        self.master.destroy()
        self.tk_root.update()

    def tv_doHelp(self):
        #tl = Toplevel()
        tl = Toplevel(self.tk_root)
        tl.label = Label(tl, text=tvHelpText)
        tl.label.pack()

    def tv_doLineUpLeaves(self):
        self.lineUpLeaves = True
        self.draw1()
        
    def tv_doOriginalBranchScale(self):
        self.lineUpLeaves = False
        self.draw1()

    def tv_configure(self, event):
        #print "tv_configure() called.  event.type=%s" % event.type
        #print dir(event)
        ew = event.width
        eh = event.height

        self.update()

        self.curWinHeight = self.winfo_height()
        self.curWinWidth = self.winfo_width()
        
        #print "tv_configure:  event.width=%s, event.height=%s, self.height=%s, self.width=%s" % (
        #    ew, eh, self.curWinHeight, self.curWinWidth)

        self.draw1()


    def draw1(self):
        #print "draw1() here!"
        self.delete(ALL)
        #self.xScale = (self.curWinWidth - (20. + self.nameStuff)) /0.001

        self.xScale = self.curWinWidth - (20. + self.nameStuff)
        self.yScale = (self.curWinHeight - 20.) / self.tv_tree.nLeaves

        if 1:
            myCapStyle='round'
            for n in self.tv_tree.root.iterPreOrder():
                if n != self.tv_tree.root:
                    if self.lineUpLeaves:
                        nParentXPos = n.parent.xPosL
                        nXPos = n.xPosL
                    else:
                        nParentXPos = n.parent.xPos1
                        nXPos = n.xPos1
                    if n.isLeaf:
                        self.create_line((nParentXPos * self.xScale) + self.xOrig,
                                         (n.yPos * self.yScale) + self.yOrig,
                                         (nXPos * self.xScale) + self.xOrig,
                                         (n.yPos * self.yScale) + self.yOrig,
                                         fill=n.br.color, width=self.lineWid, capstyle=myCapStyle)
                        self.create_text((nXPos * self.xScale) + self.xOrig + self.nameOffset,
                                         (n.yPos * self.yScale) + self.yOrig,
                                         anchor='w',
                                         text=n.name, fill=n.nameColor)
                    else:
                        self.create_line((nParentXPos * self.xScale) + self.xOrig,
                                         (n.yPos * self.yScale) + self.yOrig,
                                         (nXPos * self.xScale) + self.xOrig,
                                         (n.yPos * self.yScale) + self.yOrig,
                                         fill=n.br.color, width=self.lineWid, capstyle=myCapStyle)
                        if n.name:
                            if 0:
                                myAnchor = 'w'
                                for ch in n.iterChildren():
                                    theDiff = math.fabs((ch.yPos - n.yPos) * self.yScale)
                                    #print "internal %s, diff is %f" % (n.name, theDiff)
                                    if theDiff < 7.5:
                                        if ch.yPos >= n.yPos:
                                            myAnchor = 'sw'
                                        else:
                                            myAnchor = 'nw'
                                self.create_text((nXPos * self.xScale) + self.xOrig + self.nameOffset,
                                             (n.yPos * self.yScale) + self.yOrig,
                                             anchor=myAnchor,
                                             text=n.name, fill='darkgreen')
                            else:
                                fred = self.create_text((nXPos * self.xScale) + self.xOrig - self.nameOffset,
                                             (n.yPos * self.yScale) + self.yOrig,
                                             anchor='se',
                                             text=n.name, fill='darkgreen')
                                self.tag_bind(fred, "<Any-Enter>", self.mouseEnter)
                                self.tag_bind(fred, "<Any-Leave>", self.mouseLeave)
                                

            # vertical lines
            for n in self.tv_tree.root.iterPreOrder():
                if not n.isLeaf:
                    if self.lineUpLeaves:
                        nXPos = n.xPosL
                    else:
                        nXPos = n.xPos1
                    self.create_line((nXPos * self.xScale) + self.xOrig,
                                     (n.leftChild.yPos * self.yScale) + self.yOrig,
                                     (nXPos * self.xScale) + self.xOrig,
                                     (n.rightmostChild().yPos * self.yScale) + self.yOrig,
                                     width=self.lineWid, capstyle='round')


    # The stuff below from canvas-moving-w-mouse.py, by "matt"

    ###################################################################
    ###### Event callbacks for THE CANVAS (not the stuff drawn on it)
    ###################################################################
    def mouseDown(self, event):
        # remember where the mouse went down
        self.lastx = event.x
        self.lasty = event.y

    def mouseMove(self, event):
        # whatever the mouse is over gets tagged as CURRENT for free by tk.
        self.move(CURRENT, event.x - self.lastx, event.y - self.lasty)
        self.lastx = event.x
        self.lasty = event.y

    ###################################################################
    ###### Event callbacks for canvas ITEMS (stuff drawn on the canvas)
    ###################################################################
    def mouseEnter(self, event):
        # the CURRENT tag is applied to the object the cursor is over.
        # this happens automatically.
        self.itemconfig(CURRENT, fill="red")

    def mouseLeave(self, event):
        # the CURRENT tag is applied to the object the cursor is over.
        # this happens automatically.
        self.itemconfig(CURRENT, fill="darkgreen")


