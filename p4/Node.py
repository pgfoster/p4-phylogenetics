import func
from Var import var

class NodeBranchPart(object):
    def __init__(self):
        self.rMatrixNum = -1
        self.gdasrvNum = -1
        #self.bigP = None


class NodeBranch(object):
    def __init__(self):
        self.len = 0.1
        #self.textDrawSymbol = '-'  # See var.modelSymbols for some alternative symbols
        self.rawSplitKey = None    # Odd or even
        self.splitKey = None       # Only even
        #self.name = None
        #self.uName = None          # under-name
        #self.color =  None         # US spelling.
        #self.support = None        # A float, so that it can preserve all its significant
                                   # digits, yet can be formatted flexibly for output.
        #self.biRootCount = None  # For cons trees, where the input trees are
                                 # bi-Rooted, ie have bifurcating roots.  This
                                 # is the number of compatible input trees that
                                 # were rooted on this branch.
        self.parts = []    # NodeBranchPart() objects
        self.lenChanged = False

class NodePart(object):
    def __init__(self):
        #self.pats = None
        #self.nPats = 0
        self.compNum = -1
        #self.cl = None
        #self.cl2 = None



class Node(object):
    """A Node is a vertex in a Tree.  All but the root have a branch.

    A Node has pointers to its parent, leftChild, and sibling, any of which may be None.
    """

    def __init__(self):
        self.name = None
        self.nodeNum = -1
        self.parent = None
        self.leftChild = None
        self.sibling = None
        self.isLeaf = 0
        self.cNode = None   # Pointer to a c-struct
        self.seqNum = -1    # Zero-based seq numbering of course, so -1 means no sequence.
        self.br = NodeBranch()
        #self.rootCount = None           # For cons trees, where the input trees do not
                                        # have bifurcating roots.  This is the number of
                                        # compatible input trees that were rooted on this node.
        self.parts = []     # NodePart objects
        self.doDataPart = 0
        self.flag = 0




    ##Ignore
    def wipe(self):
        """Set the pointers parent, leftChild, and sibling to None"""
        
        self.parent = None
        self.leftChild = None
        self.sibling = None

    def rightmostChild(self):
        """Find and return the rightmostChild of self.

        If self has no children, return None.
        """
        n = self.leftChild
        if not n:
            return None
        while n.sibling:
            n = n.sibling
        return n

    def leftSibling(self):
        """Find and return the sibling on the left.

        A node has a pointer to its sibling, but that is the sibling
        on the right.  It is a bit awkward to find the sibling on the
        left, as you need to go via the parent and the leftChild of
        the parent.

        If there is no parent, return None.  If there is no
        leftSibling, return None.
        """
        if not self.parent:
            #print 'leftSibling(%i).  No parent.  returning None.' % self.nodeNum
            return None
        lsib = self.parent.leftChild
        if lsib == self:
            #print 'leftSibling(%i).  self is self.parent.leftChild.  returning None.' % self.nodeNum
            return None
        while lsib:
            if lsib.sibling == self:
                #print 'leftSibling(%i): returning node %i' % (self.nodeNum, lsib.nodeNum)
                return lsib
            lsib = lsib.sibling


    # These next 3 were suggestions from Rick Ree.  Thanks, Rick!
    # Then I added a couple more.  Note that all of these use
    # recursion, and so could bump into the recursion limit, and might
    # fail on large trees.  However, I tried iterPreOrder() on a
    # random tree of 10,000 taxa, and it was fine.

    # You can temporarily set a different recursion limit with the sys module.
    # oldlimit = sys.getrecursionlimit()
    # sys.setrecursionlimit(newLimit)

    # See also Tree.iterNodesNoRoot()
    
    def iterChildren(self):
        n = self.leftChild
        while n:
            yield n
            n = n.sibling

    def iterPostOrder(self):
        for c in self.iterChildren():
            for n in c.iterPostOrder():
                yield n
        yield self

    def iterPreOrder(self):
        yield self
        for c in self.iterChildren():
            for n in c.iterPreOrder():
                yield n

    def iterLeaves(self):
        for n in self.iterPreOrder():
            if n.isLeaf:
                yield n

    def iterInternals(self):
        for n in self.iterPreOrder():
            if not n.isLeaf:
                yield n

    def iterDown(self, showDown=False):
        """Iterates over all the nodes below self (including self)

        Starts by returning self.  And then iterates over all nodes below self.

        It does so by a combination of Node.iterPreOrder() and
        Node.iterDown() (ie recursively).  Now sometimes we want to
        know if the nodes that are returned come from iterDown()
        (strictly) or not (ie from iterPreOrder()).  If that bit of
        info is needed, then you can turn on the arg ``showDown``.
        (The following is probably bad Python practice!) When that is done, whenever
        iterDown() is called the first node that is returned will have
        the attribute ``down`` set to True.  But after it is returned,
        that ``down`` attribute is zapped (to try to keep the bloat
        down ...).  So you need to test ``if hasattr(yourNode,
        'down'):`` before you actually use it.
        
        """

        if showDown:
            self.down = True
        yield self
        if showDown:
            del(self.down)
        if self.parent:
            for c in self.parent.iterChildren():
                if c == self:
                    for n in c.parent.iterDown(showDown):
                        yield n
                else:
                    for n in c.iterPreOrder():
                        yield n
    


    # ###############################
    def getNChildren(self):
        """Returns the number of children that the node has."""
        if not self.leftChild:
            return 0
        c = self.leftChild
        counter = 0
        while c:
            c = c.sibling
            counter += 1
        return counter

    def isAncestorOf(self, otherNode):
        """Asks whether self is an an ancestor of otherNode."""
        n = otherNode
        while 1:
            n = n.parent
            if not n:
                return False
            elif n == self:
                return True
            

    def _ladderize(self, biggerGroupsOnBottom):
        """This is only used by Tree.ladderize()."""
        
        #print '====Node %i' % self.nodeNum
        if not self.leftChild:
            pass
        else:
            nLeaves = []
            children = []
            ch = self.leftChild
            while ch:
                nL = len([n2 for n2 in ch.iterLeaves()])
                nLeaves.append(nL)
                ch.nLeaves = nL
                children.append(ch)
                ch = ch.sibling
            #print '  nLeaves = %s' % nLeaves
            allOnes = True
            for ch in children:
                if ch.nLeaves > 1:
                    allOnes = False
                    break
            if not allOnes:
                children = func.sortListOfObjectsOnAttribute(children, 'nLeaves')
                if not biggerGroupsOnBottom:
                    children.reverse()
                #print '\n    Children\n    ------------'
                #for ch in children:
                #    print '    node=%i, nLeaves=%i' % (ch.nodeNum, ch.nLeaves)
                self.leftChild = children[0]
                theLeftChild = self.leftChild
                theLeftChild.sibling = None
                for ch in children[1:]:
                    theLeftChild.sibling = ch
                    theLeftChild = ch
                    theLeftChild.sibling = None
                for ch in children:
                    del(ch.nLeaves)
                for ch in self.iterChildren():
                    ch._ladderize(biggerGroupsOnBottom)

        

if var.usePfAndNumpy:
    import sys
    import pf
    #def __del__(self, freeNode=pf.p4_freeNode, dp_freeNode=pf.dp_freeNode, mysys=sys):
    def __del__(self, freeNode=pf.p4_freeNode, mysys=sys):
    #def __del__(self, freeNode=pf.p4_freeNode, dp_freeNode=pf.dp_freeNode):
    #def __del__(self, freeNode=pf.p4_freeNode):
        #if self.nodeNum == 0:
        #mysys.stdout.write('Node.__del__()   deleting node %i\n' % self.nodeNum)
        #mysys.stdout.flush()
        if self.cNode:  # Generally, cNodes are deleted before the cTree is freed.  freeNode requires the cTree!
            mysys.stdout.write('Node.__del__()  node %i (%s) has a cNode (%s).  How?!?\n' % (
                self.nodeNum, self, self.cNode))
            if self.doDataPart:
                dp_freeNode(self.cNode)
            else:
                freeNode(self.cNode)
            self.cNode = None

    Node.__del__ = __del__
    del(__del__)
    
