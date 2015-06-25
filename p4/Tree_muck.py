from Var import var
import string,math,random,copy,os
import types
import func
from Glitch import Glitch
from Node import Node,NodeBranch
if var.usePfAndNumpy:
    import numpy

def node(self, specifier):
    """Get a node based on a specifier.

    The *specifier* can be a nodeNum, name, or node object.
    """

    nodeNum = None

    gm = ['Tree.node()']

    if type(specifier) == types.IntType:
        nodeNum = specifier
    elif var.usePfAndNumpy and type(specifier) == numpy.int32:
        nodeNum = specifier
    elif isinstance(specifier, Node):   # if its a node object
        if specifier in self.nodes:
            return specifier
        else:
            gm.append("The specifier is a node object, but is not part of self.")
            raise Glitch, gm
    elif type(specifier) == types.StringType:   # if its a string
        for n in self.iterNodes():
            if n.name == specifier:
                return n
        if nodeNum == None:    # if we haven't found a node matching the specier...
            gm.append("Specifier string '%s' is not a node name.  What gives?" % specifier)
            raise Glitch, gm

    else:
        gm.append("I don't understand the specifier '%s', type '%s'." % (specifier, type(specifier)))
        raise Glitch, gm

    if nodeNum < 0 or nodeNum >= len(self.nodes):
        gm.append("The node number is out of range.")
        raise Glitch, gm

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
    # we need to know the node second from the right, which will become the newRightmostChild
    newRightmostChild = rotateNode.leftChild # as a first guess...
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
    So by default I checkBiRoot=True, and throw a Glitch if there is
    one.  If you want to draw such a pathological tree with a node
    with a single child, set checkBiRoot=False, and it will allow it.

    """

    gm = ['Tree.reRoot()']

    # The user probably does not want to reRoot if the current root is bifurcating.  So check that.
    if checkBiRoot:
        if self.root.getNChildren() == 2:
            gm.append("The tree has a bifurcating root, so you probably do not")
            gm.append("want to reRoot() it.  You can remove the bifurcating root")
            gm.append("with yourTree.removeRoot().  If you really want to reRoot()")
            gm.append("with a bifurcating root, set checkBiRoot=False in the reRoot() args.")
            raise Glitch, gm
    if self.root.isLeaf and not self.root.name:
        gm.append("The root is a leaf, but has no name.")
        gm.append("So when you reRoot() it, some other leaf will have no name.")
        gm.append("That is a recipe for trouble, and is not allowed.")
        raise Glitch, gm

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
                gm.append("Setting 'moveInternalName' implies keeping node names with their branches.")
                gm.append("The root in this tree has a name, but has no branch.")
                gm.append("So that does not work.")
                gm.append("Set arg stompRootName to work around this.")
                raise Glitch, gm
        if self.root.isLeaf and self.root.leftChild.name:
            assert self.root.name
            gm.append("The current root is a leaf, with a name.")
            gm.append("Its sole child has a name also, which is an 'internal' node name.")
            gm.append("Arg moveInternalName is turned on.")
            gm.append("So when the tree gets re-rooted, that internal node name  should stay with its branch, not its node.")
            gm.append("But the rerooted branch will already have a name-- the current root name.")
            gm.append("So that does not work.")
            raise Glitch, gm

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

        ##               +----------2:A
        ##    +----------1
        ##    |          +----------3:B
        ##    0
        ##    +----------4:C
        ##    |
        ##    +----------5:D


        ##    +----------2:A
        ##    |
        ##    1----------3:B
        ##    |
        ##    |          +----------4:C
        ##    +----------0
        ##               +----------5:D



        oldRoot = self.root   # in case this is not the first time around this loop...
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
                oldRoot.br.rawSplitKey = 1L << self.taxNames.index(oldRoot.name)  # "<<" is left-shift

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
        raise Glitch, gm
    
    if var.usePfAndNumpy:
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
        newRoot = oldRoot.leftChild # Try this one first:
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
                break # its good enough, use it
            newRoot = newRoot.sibling # its not good, try its sib

        # If we got this far an failed to find a good root, then
        # newRoot is None.  We are dealing with a tree where all the
        # root children are either leaves or have only one child.  Ok,
        # its unusual, but we do need a new root, even if it is a
        # terminal node.  So arbitrarily choose the leftChild to root
        # on.
        if not newRoot:
            newRoot = oldRoot.leftChild
            #print "x Re-rooting on node number %i" % newRoot.nodeNum
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
            #print "y reRooting on nodeNum %i" % newRoot.nodeNum
            if newRoot == oldRoot.leftChild:
                newRoot.sibling.br.len = newRoot.sibling.br.len + newRoot.br.len
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
        gm.append("Removing the root with more than two children is not implemented.")
        gm.append("Are you even sure you want to do that?")
        raise Glitch, gm


    oldRoot.wipe()
    self.nodes.remove(oldRoot)
    del oldRoot
    if newRoot:
        self.root = newRoot
    else:
        gm.append("No newRoot?   Programming error.")
        raise Glitch, gm
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
        #self.draw()
        #print "the specifier was %s" % specifier
        #print "the specified node was node number %i" % rNode.nodeNum
        #raise Glitch
        return
    rNodeParnt = rNode.parent
    
    if var.usePfAndNumpy:
        self.deleteCStuff()

    # For cases where the tree is originally with a single child root
    # -- we don't want to then delete that root below.
    assert self.root.leftChild
    isOriginallySingleChildRoot = False
    if not self.root.leftChild.sibling:
        isOriginallySingleChildRoot = True

    # For cases where we originally have a bi-Root, that we would like to keep.
    isOriginallyBiRoot = False
    if self.root.getNChildren() == 2:
        isOriginallyBiRoot = True
        

    #print "Removing node number", rNode.nodeNum
    #hitList = []
    #self.recursivelyListNodeIndicesDownTo(hitList, rNode)

    # getNodeNumsAbove does not require setting self.preAndPostOrderAreValid = 0.  Its irrelevant.
    hitList = self.getNodeNumsAbove(rNode.nodeNum) # does not include rNode
    hitList.append(rNode.nodeNum)
    #hitList.sort()
    #hitList.reverse()
    #print "Hit list is", hitList
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
                #print "leftSib is node %i" % leftSib.nodeNum
                leftSib.sibling = rNode.sibling
            else:
                gm.append("leftSib is None.  This shouldn't happen")
                gm.append("Programming error?")
                raise Glitch, gm

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
        #self.dump(node=1)
        self.preOrder = None
        self.postOrder = None
        self.preAndPostOrderAreValid = 0
        #self.draw()  # This won't work unless preAndPostOrderAreValid set to 0




    # the parent of the removed node may now only have one child,
    # in which case it (the parent) should be removed

    if alsoRemoveSingleChildParentNode or alsoRemoveBiRoot or haveRemovedSingleChildRoot:

        ignoreBrLens = True
        if var.usePfAndNumpy:
            self.preOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
            self.postOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
        else:
            self.preOrder = [var.NO_ORDER] * len(self.nodes)
            self.postOrder = [var.NO_ORDER] * len(self.nodes)

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
                    #print 'rNodeParnt is a left child'
                    rNodeParnt.parent.leftChild = rNodeParnt.leftChild
                else:
                    #print 'rNodeParnt is a middle or right child'
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
    if var.usePfAndNumpy:
        self.preOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
        self.postOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
    else:
        self.preOrder = [var.NO_ORDER] * len(self.nodes)
        self.postOrder = [var.NO_ORDER] * len(self.nodes)

    if len(self.nodes) > 1:
        self.setPreAndPostOrder()

def removeAboveNode(self, specifier, newName):
    """Remove everything above an internal node, making it a leaf, and so needing a new name.
    """
    rNode = self.node(specifier)
    assert rNode != self.root
    assert not rNode.isLeaf
    toDelete = [n for n in rNode.iterPostOrder() if n != rNode]
    #print [n.nodeNum for n in toDelete]
    for n in toDelete:
        #print 'deleting node %i' % n.nodeNum
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
    #print "Collapsing node %i" % theNode.nodeNum

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
    #self.nodes.remove(theNode)
    #del(theNode)
    
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
        raise Glitch, gm
    rNodeParnt = rNode.parent
    if not allowSingleChildNode:
        if rNodeParnt.getNChildren() < 3:
            #self.draw()
            gm.append("The arg allowSingleChildNode is turned off.")
            gm.append("This would be for those cases where the parent of the subTree has more than 2 children.")
            raise Glitch, gm

    #self.deleteCStuff()

    #print "rNode is node %i" % rNode.nodeNum
    #print "rNodeParnt is node %i" % rNodeParnt.nodeNum
    
    # Disconnect it from the tree.
    # the rNode is the left, middle, or right child of the parent
    if rNodeParnt.leftChild == rNode:
        #print "its the left child"
        rNodeParnt.leftChild = rNode.sibling
        rNode.sibling = None
        rNode.parent = None
    else:
        #print "its a middle or a right child"
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
        raise Glitch, gm
    stNode.parent = newParent
    if beforeNode == None:  # easy -- just add it to the rightmost child.
        theRMChildOfNewParent = newParent.rightmostChild()
        theRMChildOfNewParent.sibling = stNode
    else:
        bNode = self.node(beforeNode)
        if bNode.parent != newParent:
            gm.append("The parent of the 'beforeNode' should be the newParent.")
            raise Glitch, gm
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
        raise Glitch, gm
    elif aNode2 == aNode1.parent:
        temp = aNode1
        aNode1 = aNode2
        aNode2 = temp
    else:
        gm.append("The 2 specified nodes should have a parent-child relationship")
        raise Glitch, gm

    if var.usePfAndNumpy:
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
        if var.usePfAndNumpy:
            self.preOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
            self.postOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
        else:
            self.preOrder = [var.NO_ORDER] * len(self.nodes)
            self.postOrder = [var.NO_ORDER] * len(self.nodes)
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
        raise Glitch, gm
    if not self.taxNames:
        gm.append("Self (ie the tree) needs to have taxNames set.")
        raise Glitch, gm

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

    from Trees import Trees
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

    #self.draw()
    self.root._ladderize(biggerGroupsOnBottom)
    self.preAndPostOrderAreValid = 0
    #self.draw()



def randomizeTopology(self, randomBrLens=True):

    gm = ["Tree.randomizeTopology()"]
    if self.root.getNChildren() != 3 or not self.isFullyBifurcating():
        gm.append("Should be a fully bifurcating tree, this week.  Fix me?")
        raise Glitch, gm
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

        #print "lChild = node %i" % lChild.nodeNum

        ##    +----------1:oldLeftSib
        ##    |
        ##    +----------2:lChild
        ##    0
        ##    +----------3:lChildSib
        ##    |
        ##    +----------4:oldLChildSibSib


        ##    +----------1:oldLeftSib
        ##    |
        ##    |          +----------3:lChild
        ##    0----------2(n)
        ##    |          +----------4:lChildSib
        ##    |
        ##    +----------5:oldLChildSibSib

        n = internals[iIndx]
        iIndx += 1
        n.parent = None
        n.sibling = None
        n.leftChild = None
        n.nodeNum = nodeNum
        nodeNum += 1
        lChildSib = lChild.sibling  # guarranteed to have one
        oldLChildSibSib = lChildSib.sibling # ditto
        oldLeftSib = lChild.parent.leftChild # first guess ...
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


    #self.dump(all=True)
    #self.draw()

    # The way it is now, the root rightmost child is always a
    # leaf.  Not really random, then, right?  So choose a random
    # internal node, and re-root it there.
    #print "nTax=%i, len(t.nodes)=%i" % (nTax, len(t.nodes))
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
    #self.dump(all=True)
    
    
    



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
        raise Glitch, gm
    self.checkTaxNames()

    f = open(thePaupLogFileName, 'r')
    while 1:
        aLine = f.readline()
        if not aLine:
            gm.append("No Bipartitions line in %s?" % thePaupLogFileName)
            raise Glitch, gm
        if aLine.startswith('Bipartitions found'):
            #print aLine
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
        #print "a Got Line ==>%s<==" % aLine
        if not aLine:
            gm.append('Unexpected end of file.')
            raise Glitch, gm
        aLine = aLine.rstrip()
        #print "b Got Line ==>%s<==" % aLine
        if len(aLine):
            if aLine[0] == '-':
                #print prevLine
                #print aLine
                if prevLine.endswith('Freq') or prevLine.endswith('%'):
                    if prevLine.endswith('Freq'):
                        hasFreqOnly = 1
                    elif prevLine.endswith('%'):
                        hasFreqOnly = 0
                        
                    cycleType = LAST_CYCLE
                    #print "We are now in the last cycle.  hasFreqOnly=%s" % hasFreqOnly
                    splitLine = string.split(prevLine)
                    thisKeyLength = len(splitLine[0])
                    keyCounter = 0
                elif cycleType == MIDDLE_CYCLE:
                    #print "We are now in a middle cycle"
                    thisKeyLength = len(prevLine)
                    keyCounter = 0
                elif cycleType == FIRST_CYCLE:
                    gm.append('This should never happen.')
                    raise Glitch, gm
                else:
                    cycleType = FIRST_CYCLE
                    #print "We are now in the first cycle"
                    thisKeyLength = len(prevLine)


            elif aLine[0] in ['.', '*']:
                if cycleType in [FIRST_CYCLE, MIDDLE_CYCLE] and len(aLine) != thisKeyLength:
                    gm.append("Unequal key lengths?!?")
                    gm.append("thisKeyLength = %i" % thisKeyLength)
                    gm.append("%s" % aLine)
                    raise Glitch, gm

                if cycleType == FIRST_CYCLE:
                    theKeys.append(aLine)
                elif cycleType == MIDDLE_CYCLE:
                    theKeys[keyCounter] = theKeys[keyCounter] + aLine
                    keyCounter = keyCounter + 1
                    if keyCounter > nKeys:
                        gm.append("Too many keys.")
                        raise Glitch, gm
                elif cycleType == LAST_CYCLE:
                    # It will usually be a line like: ...*....*.     67.83  67.9%
                    # But it will sometimes be a line like: ...*....*.        68
                    splitLine = string.split(aLine)
                    theKey = splitLine[0]
                    if hasFreqOnly:
                        theSupport = float(splitLine[-1])
                        theSupport /= 100.0
                    else:
                        theSupport = float(splitLine[-1][:-1]) # don't read the % sign
                        theSupport /= 100.0

                    if len(theKey) != thisKeyLength:
                        gm.append("Unequal key lengths?!?")
                        gm.append("thisKeyLength = %i" % thisKeyLength)
                        gm.append("%s" % theKey)
                        raise Glitch, gm
                    if accumulatedKeyLength:
                        theKeys[keyCounter] = theKeys[keyCounter] + theKey
                        theHash[theKeys[keyCounter]] = theSupport
                        keyCounter = keyCounter + 1
                        if keyCounter > nKeys:
                            gm.append("Too many keys.")
                            raise Glitch, gm
                    else:  # LAST_CYCLE is also the FIRST_CYCLE, ie the table is in one section.
                        theKeys.append(theKey)
                        theHash[theKey] = theSupport
                else:
                    gm.append("This should never happen.")
                    raise Glitch, gm
            else:
                pass # Skip lines with numbers
        else: # a blank line
            if cycleType == FIRST_CYCLE:
                cycleType = MIDDLE_CYCLE
                nKeys = len(theKeys)
                #print "At the end of the first cycle, got %i keys" % nKeys
                accumulatedKeyLength = thisKeyLength
            elif cycleType == MIDDLE_CYCLE:
                accumulatedKeyLength += thisKeyLength
            elif cycleType == LAST_CYCLE:
                accumulatedKeyLength += thisKeyLength
                nKeys = len(theKeys)
                #print "finished"
                #print theKeys[0]
                break
            else:
                pass
    f.close()

    if 0:
        print "Finished getting split strings, and supports"
        print "Got %i items in the hash" % len(theHash)
        print "accumulatedKeyLength = %i" % accumulatedKeyLength
        print "nKeys = %i" % nKeys
        #self.draw()

    # I will need to make splitStrings (in dot-star notation) from
    # splitKeys.  To do that, I can use the
    # func.getSplitStringFromKey(theKey, nTax) function.

    self.makeSplitKeys()
    for n in self.nodes:
        if n != self.root:
            if not n.isLeaf:
                theNodeSplitString = func.getSplitStringFromKey(n.br.splitKey, self.nTax)
                if theHash.has_key(theNodeSplitString):
                    if hasattr(n.br, 'support') and n.br.support is not None:
                        gm.append("Node %i already has a br.support." % n.nodeNum)
                        gm.append("I am refusing to clobber it with the split support.")
                        gm.append("Either fix the tree or fix this method.")
                        raise Glitch, gm
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
        raise Glitch, gm
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
    f.write("p4_renameForPhylip_originalNames = %s\np4_renameForPhylip_dict = %s\n" % (originalNames,d))
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
        from __main__ import p4_renameForPhylip_dict,p4_renameForPhylip_originalNames
    else:
        gm.append("The dictionary file '%s' can't be found." % dictFName)
        raise Glitch, gm
    for n in self.iterNodes():
        if n.isLeaf:
            if p4_renameForPhylip_dict.has_key(n.name):
                n.name = p4_renameForPhylip_dict[n.name]
            else:
                gm.append("The dictionary does not contain a key for '%s'." % n.name)
                raise Glitch, gm
    if p4_renameForPhylip_originalNames:
        self.taxNames = p4_renameForPhylip_originalNames
    else:
        if self.taxNames:
            gm.append("self.taxNames is set, and should be replaced, but")
            gm.append("p4_renameForPhylip_originalNames is None. ?!?")
            raise Glitch, gm
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
        raise Glitch, gm
    loc = {}
    execfile(dictFileName, {}, loc)
    try:
        p4DupeSeqRenameDict = loc['p4DupeSeqRenameDict']
    except KeyError:
        gm.append("Can't get the dictionary named 'p4DupeSeqRenameDict' from the dict file.")
    #print p4DupeSeqRenameDict

    kk = p4DupeSeqRenameDict.keys()
    #print kk

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
    #print "levels = %i" % levels
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
            n.br.lenL = float(n.lul_pos)/float(levels) * rootToLeaf
        else:
            n.br.lenL = float(n.lul_pos - p.lul_pos)/float(levels) * rootToLeaf

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
        self.removeNode(ch, alsoRemoveSingleChildParentNode=False, alsoRemoveBiRoot=False)


        
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
    from Tree import Tree
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
                    pass # skip parent and sibling
                else:
                    # parents and siblings.  There will always be a parent.
                    st.nodes[nodeNumDict[selfNode.nodeNum]].parent = \
                           st.nodes[nodeNumDict[selfNode.parent.nodeNum]]
                    if selfNode.sibling:
                        st.nodes[nodeNumDict[selfNode.nodeNum]].sibling = \
                           st.nodes[nodeNumDict[selfNode.sibling.nodeNum]]
            st.root.leftChild = st.nodes[nodeNumDict[dupeNode.nodeNum]]
            st.root.leftChild.parent = st.root
    else: # down
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
    assert theSubTree.root.leftChild and not theSubTree.root.leftChild.sibling # its a root on a stick
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
    #print 
    #for n in self.nodes:
    #    print n.nodeNum
    #print "self.taxNames is %s" % self.taxNames
    #print "subTreeTaxNames %s" % subTreeTaxNames
    if self.taxNames:
        self._taxNames += subTreeTaxNames
    
    if self._nTax:
        self._nTax += len(subTreeTaxNames)
    self.preAndPostOrderAreValid=False
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
    #self.draw()
    #print "attachmentNode is %i" % aNode.nodeNum
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
    #self.draw()
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
    self.preAndPostOrderAreValid=False
    self.preOrder = None
    self.postOrder = None
    self.getPreAndPostOrderAboveRoot()
    #self.dump(node=True)
    #self.draw()
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
    #self.draw()
    #print "attachmentNode is %i" % aNode.nodeNum
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
    self.preAndPostOrderAreValid=False
    self.preOrder = None
    self.postOrder = None
    self.getPreAndPostOrderAboveRoot()
    #self.dump(node=True)
    #self.draw()
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
        for n in theNode.iterInternals(): # Includes theNode, which we want
            #if n.getNChildren() != 2:
            #    return False
            if n.leftChild and n.leftChild.sibling:
                if n.leftChild.sibling.sibling:
                    return False
            else:
                return False
        return True
    else:
        for n in theNode.iterDown():  # Includes theNode, which we do not want
            if n != theNode:
                if not n.isLeaf:
                    if n == self.root:
                        #if n.getNChildren() != 3:
                        #    return False
                        if n.leftChild and n.leftChild.sibling and n.leftChild.sibling.sibling:
                            if n.leftChild.sibling.sibling.sibling:
                                return False
                        else:
                            return False
                    else:
                        #if n.getNChildren() != 2:
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
        upperNode = self.node(upperNodeSpec) # This makes sure that upperNode is part of self.
    else:
        candidates = [n for n in self.iterInternalsNoRoot()]
        upperNode = random.choice(candidates)


    # Want the upperNode to have at least 2 children
    upperChildren = [n for n in upperNode.iterChildren()]
    if len(upperChildren) < 2:
        gm.append("upperNode needs to have at least 2 children.")
        raise Glitch, gm
    if upperNode.parent:
        lowerNode = upperNode.parent
    else:
        gm.append("upperNode needs to have a parent node.")
        raise Glitch, gm
    lowerChildren = [n for n in lowerNode.iterChildren() if n != upperNode]
    if 0:
        if lowerNode.parent:
            if len(lowerChildren) < 1:
                gm.append("The lower node has a parent.")
                gm.append("It needs at least one more child besides the upperNode.")
                raise Glitch, gm
        else:
            if len(lowerChildren) < 2:
                gm.append("The lower node does not have a parent.")
                gm.append("It needs at least 2 children besides the upperNode.")
                raise Glitch, gm
    if len(lowerChildren) < 1:
        gm.append("The lower node needs at least one more child besides the upperNode.")
        raise Glitch, gm
        
    upperSubTreeNode = random.choice(upperChildren)
    lowerSubTreeNode = random.choice(lowerChildren)

    upperSubTree = self.pruneSubTreeWithoutParent(upperSubTreeNode, allowSingleChildNode=True)
    lowerSubTree = self.pruneSubTreeWithoutParent(lowerSubTreeNode, allowSingleChildNode=True)

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
        if var.usePfAndNumpy:
            self.preOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
            self.postOrder = numpy.array([var.NO_ORDER] * len(self.nodes), numpy.int32)
        else:
            self.preOrder = [var.NO_ORDER] * len(self.nodes)
            self.postOrder = [var.NO_ORDER] * len(self.nodes)

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
        raise Glitch, gm

    pnNode = self.node(pruneNode)
    grNode = self.node(graftNode)

    #self.draw()

    if var.usePfAndNumpy:
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
    #print [n.nodeNum for n in subTreeNodes]
    if grNode in subTreeNodes:
        if above:
            gm.append("grNode %i is part of the pruned subtree from %i-%i up.  No workee!" % (
                grNode.nodeNum, pnNodeParnt.nodeNum, pnNode.nodeNum))
        else:
            gm.append("grNode %i is part of the pruned subtree below %i-%i.  No workee!" % (
                grNode.nodeNum, pnNodeParnt.nodeNum, pnNode.nodeNum))
            
        raise Glitch, gm

    # Prune it from the tree.
    if pnNodeParnt == self.root:
        if grNode.parent == self.root: # as well,
            gm.append("prune node and graft node both have root as parent -- ie same origin and destination.")
            raise Glitch, gm
        # Check if removal of the subtree will result in only 2 taxa
        singles = 0
        for ch in self.root.iterChildren():
            if ch != pnNode:
                if not ch.leftChild:
                    singles += 1
        if singles == 2:
            gm.append("Removing subtree at %i will leave only 2 taxa." % pnNode.nodeNum)
            raise Glitch, gm

        newRoot = None
        for ch in self.root.iterChildren():
            if ch != pnNode:
                if ch.leftChild:
                    newRoot = ch
                    break
        #print "rerooting to node %i" % newRoot.nodeNum
        self.reRoot(newRoot)
        #self.draw()
            
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
        else: # if pnNodeParnt.parent is the root, it can have 3 children, and maybe this ...
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
        else: # if pnNodeParnt.parent is the root, it can have 3 children, and maybe this ...
            pnNodeParnt.parent.leftChild.sibling.sibling = pnNodeLeftSib

    # To look at the pruned tree, before grafting ...
    if 0:
        print "removal of subtree at %i-%i gives .." % (pnNodeParnt.nodeNum, pnNode.nodeNum)
        self._nTax = 0
        self.preAndPostOrderAreValid = 0
        self.draw()  # This won't work unless preAndPostOrderAreValid set to 0

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
        self.draw()  # This won't work unless preAndPostOrderAreValid set to 0

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
            # iterPostOrder() ends with the pNode.  We need at least the pNode and 3 others, so total 3 is too few.
            nodesUp = [n2 for n2 in pNode.iterPostOrder()]
            #self.draw()
            #print "--", pNode.nodeNum, [n2.nodeNum for n2 in nodesUp]
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
        #print "===", "pNode=%i" % pNode.nodeNum, "gNNum=%i" % gNNum, subtreeNodeNums
        if pNode.parent==self.root and self.node(gNNum).parent==self.root:
            continue
        if not myAbove and self.node(gNNum).parent == pNode:
            continue
        break
    if safety >= giveUp:
        print "pNode is %i, above=%s" % (pNode.nodeNum, myAbove)
        self.draw()
        raise Glitch
    #print "spr()  pruneNum %i, above=%s, graftNum %i" % (pNode.nodeNum, myAbove, gNNum)
    
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
        #inTree.makeSplitKeys()
        
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
            rfDist = sDupe.topologyDistance(inTree, metric='sd', resetSplitKeySet=True)
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
    
    
        

