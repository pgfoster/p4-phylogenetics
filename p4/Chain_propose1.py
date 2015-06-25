import random,math
import pf,func
from Var import var
from Glitch import Glitch
import sys

#localCalls = 0


def proposeRoot3(self, theProposal):
    """For non-biRooted trees.  Root on another internal node."""
    internalsNoRoot = [n for n in self.propTree.iterInternalsNoRoot()]
    if len(internalsNoRoot):
        newRoot = random.choice(internalsNoRoot)
        self.propTree.reRoot(newRoot, moveInternalName=False, fixRawSplitKeys=self.mcmc.constraints)
    else:
        print "Chain.proposeRoot3().  No other internal nodes.  Fix me."
    self.logProposalRatio = 0.0
    self.logPriorRatio = 0.0
    #if self.mcmc.constraints:
    #    print "checkSplitKeys() at the end of root3"
    #    self.propTree.checkSplitKeys()

def proposeBrLen(self, theProposal):
    #gm = ['Chain.proposeBrLen']

    # Choose a node.
    nodesNoRoot = [n for n in self.propTree.iterNodesNoRoot()]
    theNode = random.choice(nodesNoRoot)
    #theNode = self.propTree.nodes[1]
    oldBrLen = theNode.br.len

    if 1: # "Multiplier" proposal
        newBrLen = oldBrLen * math.exp(theProposal.tuning * (random.random() - 0.5))

        # Logarithmic reflect if needed
        while (newBrLen < var.BRLEN_MIN) or (newBrLen > var.BRLEN_MAX):
            if newBrLen < var.BRLEN_MIN:
                newBrLen = var.BRLEN_MIN * var.BRLEN_MIN / newBrLen
            elif newBrLen > var.BRLEN_MAX:
                newBrLen = var.BRLEN_MAX * var.BRLEN_MAX / newBrLen
        theNode.br.len = newBrLen
        self.logProposalRatio = math.log(newBrLen/oldBrLen)
    else:  # Sliding window.  
        newBrLen = oldBrLen + (theProposal.tuning * (random.random() - 0.5))
        #newBrLen = oldBrLen + (2.0 * (random.random() - 0.5))

        # Linear reflect if needed
        while (newBrLen < var.BRLEN_MIN) or (newBrLen > var.BRLEN_MAX):
            if newBrLen < var.BRLEN_MIN:
                newBrLen = (var.BRLEN_MIN - newBrLen) + var.BRLEN_MIN
            elif newBrLen > var.BRLEN_MAX:
                newBrLen = var.BRLEN_MAX - (newBrLen - var.BRLEN_MAX)
        theNode.br.len = newBrLen
        self.logProposalRatio = 0.0

    if hasattr(self.mcmc.tunings, 'doInternalBrLenPrior') and self.mcmc.tunings.doInternalBrLenPrior:
        if theNode.isLeaf:
            self.logPriorRatio = self.mcmc.tunings.brLenPriorLambda * (oldBrLen - newBrLen)
        else:
            self.logPriorRatio = self.mcmc.tunings.brLenPriorLambdaForInternals * (oldBrLen - newBrLen)
    else:
        if self.mcmc.tunings.brLenPriorType == 'exponential':
            self.logPriorRatio = self.mcmc.tunings.brLenPriorLambda * (oldBrLen - newBrLen)
        else:
            self.logPriorRatio = 0.

    if var.doMcmcSp:
        theNode.br.lenChanged = True
        
    


def proposeLocal(self, theProposal):  # from BAMBE and MrBayes.

    # doAbort is set if brLens are too long or too short, or if a constraint is violated.
    
    #global localCalls
    #localCalls += 1
    #print 'localCalls %i' % localCalls

    gm = ['Chain.proposeLocal()']
    theProposal.topologyChanged = 0
    theProposal.doAbort = False
    pTree = self.propTree
    dbug = False
    if 0 and self.mcmc.gen == 0:
        dbug = True

    if dbug:
        for n in pTree.iterNodesNoRoot():
            n.br.oldLen = n.br.len
            n.br.oldNode = n


    if dbug:
        print "=" * 80
        #print "proposeLocal() starting with this tree ..."
        #pTree.draw(width=80, addToBrLen=0.0)
        #for n in pTree.iterInternalsNoRoot():
        #    print n.nodeNum, n.br.splitKey
        for n in pTree.iterNodes():
            n.oldName = n.name
        for n in pTree.iterLeavesNoRoot():
            n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
        for n in pTree.iterInternalsNoRoot():
            n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
        print "proposeLocal() starting with this tree ..."
        pTree.draw(width=100, addToBrLen=0.2, model=True)
        if self.mcmc.constraints:
            pTree.checkSplitKeys(useOldName=True)

    if pTree.root.getNChildren() == 2:
        isBiRoot = True
        gm.append("This method is not working for biRoot'd trees yet.")
        raise Glitch, gm
    else:
        isBiRoot = False

    assert pTree.nInternalNodes > 1, "For local, we need trees with more than 1 internal node."

    # We want to find a node with a great-grandparent, but that might be
    # impossible!, eg if it is a 4-taxon tree rooted on an internal
    # node.  So at least find a node with a grandparent.
    usedNodes = [n for n in pTree.iterNodes()]
    candidateC = random.choice(usedNodes)
    safety = 0
    while 1:
        if candidateC.parent and candidateC.parent.parent:
            break
        else:
            safety += 1
            if safety > 100:
                pTree.draw()
                gm.append("Unable to find a node with a grandparent after 100 tries.")
                gm.append("The propTree has %i internal nodes." % pTree.nInternalNodes)
                del(pTree.nInternalNodes)
                gm.append("Recalculated: the propTree has %i internal nodes." % pTree.nInternalNodes)
                raise Glitch, gm
            candidateC = random.choice(usedNodes)

    # Check whether candidateC has a great-grandparent.  (We know it
    # has a grandparent).  If so, then fine.  If not, we need to
    # re-root.
    oldRoot = None
    if candidateC.parent.parent.parent:
        pass
    else:
        oldRoot = pTree.root
        possibleRoots = []
        for n in oldRoot.iterChildren():
            if n != candidateC.parent:
                possibleRoots.append(n)
        if not possibleRoots:
            print "=" * 50
            pTree.draw()
            gm.append("Programming error. Could not find any possibleRoots")
            raise Glitch, gm
            
        newRoot = random.choice(possibleRoots)
        pTree.reRoot(newRoot, moveInternalName=False, fixRawSplitKeys=self.mcmc.constraints)
        
    if 0 and dbug:
        print "candidateC is node %i" % candidateC.nodeNum
        if oldRoot:
            print "I had to reRoot to node %i." % newRoot.nodeNum
        pTree.draw(width=80, showInternalNodeNames=1, addToBrLen=0.0)


    ## We want a tree like this:
    ##                   +------c
    ##           +-------|(v)
    ##    +------|(u)    +------d
    ##    |      |
    ##    |(a)   +-------b
    ##    |
    ##    +------X (which might be None)

    # set up the nodes as in Larget and Simon MBE, pg 754, fig 4
    c = candidateC
    v = c.parent
    safety = 0
    while c != v.leftChild:
        pTree.rotateAround(v)
        safety += 1
        if safety > 100:
            gm.append("Unable to make c as v's leftChild.")
            raise Glitch, gm
    assert c.sibling
    d = c.sibling
    u = v.parent
    safety = 0
    while v != u.leftChild:
        pTree.rotateAround(u)
        safety += 1
        if safety > 100:
            gm.append("Unable to make v as u's leftChild.")
            raise Glitch, gm
    assert v.sibling
    b = v.sibling
    a = u.parent
    safety = 0
    while a.leftChild != u:
        pTree.rotateAround(a)
        safety += 1
        if safety > 100:
            gm.append("Unable to make u as a's leftChild.")
            raise Glitch, gm

    if dbug:
        #v.oldName = v.name
        #u.oldName = u.name
        #c.oldName = c.name
        #d.oldName = d.name
        #a.oldName = a.name
        #b.oldName = b.name
        for n in pTree.iterNodes():
            n.name = n.oldName
        for n in pTree.iterLeavesNoRoot():
            n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
        for n in pTree.iterInternalsNoRoot():
            n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
        
        if v.name:
            v.name += '(v)'
        else:
            v.name = 'v'
        if u.name:
            u.name += '(u)'
        else:
            u.name = 'u'
        if c.name:
            c.name += '(c)'
        else:
            c.name = 'c'
        if d.name:
            d.name += '(d)'
        else:
            d.name = 'd'
        if a.name:
            a.name += '(a)'
        else:
            a.name = 'a'
        if b.name:
            b.name += '(b)'
        else:
            b.name = 'b'
        

        
        print "Label the nodes a,b,c,d,u,v, and arrange into a 'standard form' ..."
        pTree.draw(width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)
        
    ## At this point, the tree should look like this:
    ##                   +------c
    ##           +-------|(v)
    ##    +------|(u)    +------d
    ##    |      |
    ##    |(a)   +-------b
    ##    |
    ##    +------X (which might be None)

    m = c.br.len + v.br.len + u.br.len
    x = u.br.len
    y = x + v.br.len
    newMRatio = math.exp(theProposal.tuning * (random.random() - 0.5))  # by default, 0.909 to 1.1
    newM = m * newMRatio

    # Hopefully these checks will not be needed forever.
    ##tooShort = False
    ##if c.br.len < var.BRLEN_MIN:
    ##    gm.append("c.br.len (%g) is too short" % c.br.len)
    ##    tooShort = True
    ##elif v.br.len < var.BRLEN_MIN:
    ##    gm.append("v.br.len (%g) is too short" % v.br.len)
    ##    tooShort = True
    ##elif u.br.len < var.BRLEN_MIN:
    ##    gm.append("u.br.len (%g) is too short" % u.br.len)
    ##    tooShort = True
    ##elif newM < (3.0 * var.BRLEN_MIN):
    ##    gm.append("newM (%g) is too short." % newM)
    ##    tooShort = True
    ##if tooShort:
    ##    if self.lastProposal: # a tuple, see a few lines below
    ##        gm.append("Last proposal: %s, accepted=%s, topologyChanged=%s" % self.lastProposal)
    ##    else:
    ##        gm.append("This appears to be the first proposal.")
    ##    raise Glitch, gm

    if 0 and dbug:
        print
        print "m, the sum of brLens from a up to c, is %f" % m
        print "x, from a to u, is %f" % x
        print "y, from a to v, is %f" % y
        print "newMRatio is %f" % newMRatio
        print "newM is %f" % newM
    
    #################################################################
    # Detach either u or v, then re-attach somewhere between a and c.
    #################################################################

    if random.random() < 0.5:
        ## detach u
        ##                   +------c
        ##           +-------|(v)
        ##    +------|(u)    +------d
        ##    |      |
        ##    |(a)   +-------b
        ##    |
        ##    +------X
        ##
        ##               +----------c
        ##    +----------|(v)
        ##    |(a)       +----------d
        ##    |
        ##    +----------X

        newY = y * newMRatio
        newX = random.random() * newM

        # newX should be at least var.BRLEN_MIN away from newY
        safety = 0
        while math.fabs(newY - newX) < var.BRLEN_MIN:
            newX = random.random() * newM
            safety += 1
            if safety > 100:
                if dbug:
                    print "Unable to place newX sufficiently far away from newY"
                theProposal.doAbort = True
                return

        if 0 and dbug:
            print "Choose to detach node u (not v)"
            print "newY is (%f * %f =) %f" % (y, newMRatio, newY)
            print "newX, a random spot along newM, is %f" % newX
            if newX < newY:
                print "-> Since newX is still less than newY, there will be no topology change."
            else:
                print "-> Since newX is now more than newY, there will be a topology change."

        a.leftChild = v
        v.parent = a
        v.sibling = u.sibling # which might be None

        # now re-attach at newX
        if newX < newY:
            # no topology change, set up the same as above
            ##               +----------c
            ##    +----------|(v)
            ##    |(a)       +----------d
            ##    |
            ##    +----------X
            ##
            ##           newX    newY   newM
            ##           +       +      +
            ##
            ##                   +------c
            ##           +-------|(v)
            ##    +------|(u)    +------d
            ##    |      |
            ##    |(a)   +-------b
            ##    |
            ##    +------X
            ##
            a.leftChild = u
            u.parent = a
            u.leftChild = v
            u.sibling = v.sibling
            v.parent = u
            v.sibling = b
            u.br.len = newX
            v.br.len = newY - newX
            c.br.len = newM - newY
            if 0 and dbug:
                print "-> detach u, reattach between a and v, so no topology change"
                pTree.draw(width=80, showInternalNodeNames=1, addToBrLen=0.0)


        else:
            ## a topology change
            ##
            ##               +----------c
            ##    +----------|(v)
            ##    |(a)       +----------d
            ##    |
            ##    +----------X
            ##
            ##
            ##           newY    newX   newM
            ##           +       +      +
            ##
            ##                   +------c
            ##           +-------|(u)
            ##    +------|(v)    +------b
            ##    |      |
            ##    |(a)   +-------d
            ##    |
            ##    +------X
            v.leftChild = u
            u.parent = v
            u.leftChild = c
            u.sibling = d
            c.parent = u
            c.sibling = b
            v.br.len = newY
            u.br.len = newX - newY
            c.br.len = newM - newX
            if self.mcmc.constraints:
                tempRawSplitKey = v.br.rawSplitKey
                tempSplitKey = v.br.splitKey
                v.br.rawSplitKey = u.br.rawSplitKey
                v.br.splitKey = u.br.splitKey
                u.br.rawSplitKey = tempRawSplitKey
                u.br.splitKey = tempSplitKey
            pTree.preAndPostOrderAreValid = 0
            theProposal.topologyChanged = 1
            if dbug:
                print "-> detach u, re-attach between v and c, so there is a topology change"
                for n in pTree.iterNodes():
                    n.name = n.oldName
                for n in pTree.iterLeavesNoRoot():
                    n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
                for n in pTree.iterInternalsNoRoot():
                    n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)

                if v.name:
                    v.name += '(v)'
                else:
                    v.name = 'v'
                if u.name:
                    u.name += '(u)'
                else:
                    u.name = 'u'
                if c.name:
                    c.name += '(c)'
                else:
                    c.name = 'c'
                if d.name:
                    d.name += '(d)'
                else:
                    d.name = 'd'
                if a.name:
                    a.name += '(a)'
                else:
                    a.name = 'a'
                if b.name:
                    b.name += '(b)'
                else:
                    b.name = 'b'
                pTree.draw(width=100, showInternalNodeNames=1, addToBrLen=0.2)

    else:
        ## detach v
        ##                   +------c
        ##           +-------|(v)
        ##    +------|(u)    +------d
        ##    |      |
        ##    |(a)   +-------b
        ##    |
        ##    +------X
        ##
        ##
        ##               +----------c
        ##    +----------|(u)
        ##    |(a)       +----------b
        ##    |
        ##    +----------X
        ##

        newX = x * newMRatio
        newY = random.random() * newM

        # newY should be at least var.BRLEN_MIN away from newX
        safety = 0
        while math.fabs(newY - newX) < var.BRLEN_MIN:
            newY = random.random() * newM
            safety += 1
            if safety > 100:
                if dbug:
                    print "Unable to place newY sufficiently far away from newX"
                theProposal.doAbort = True
                return

        if 0 and dbug:
            print "Choose to detach node v (not u)"
            print "newX is (%f * %f =) %f" % (x, newMRatio, newX)
            print "newY, a random spot along newM, is %f" % newY
            if newY < newX:
                print "-> Since newY is now less than newX, there will be a topology change."
            else:
                print "-> Since newY is still more than newX, there will not be a topology change."

        u.leftChild = c
        c.parent = u
        c.sibling = b

        # now reattach at newY
        if newX < newY:
            # no topology change
            ##
            ##               +----------c
            ##    +----------|(u)
            ##    |(a)       +----------b
            ##    |
            ##    +----------X
            ##
            ##           newX    newY   newM
            ##           +       +      +
            ##
            ##                   +------c
            ##           +-------|(v)
            ##    +------|(u)    +------d
            ##    |      |
            ##    |(a)   +-------b
            ##    |
            ##    +------X
            ##
            u.leftChild = v
            v.parent = u
            v.leftChild = c
            v.sibling = b
            c.parent = v
            c.sibling = d
            u.br.len = newX
            v.br.len = newY - newX
            c.br.len = newM - newY
            if 0 and dbug:
                print "-> detach v, reattach between u and c, so no topology change"
                pTree.draw(width=80, showInternalNodeNames=1, addToBrLen=0.0)
        else:
            ## with a topology change
            ##
            ##               +----------c
            ##    +----------|(u)
            ##    |(a)       +----------b
            ##    |
            ##    +----------X
            ##
            ##
            ##           newY    newX   newM
            ##           +       +      +
            ##
            ##                   +------c
            ##           +-------|(u)
            ##    +------|(v)    +------b
            ##    |      |
            ##    |(a)   +-------d
            ##    |
            ##    +------X
            a.leftChild = v
            v.parent = a
            v.leftChild = u
            v.sibling = u.sibling
            u.parent = v
            u.sibling = d
            v.br.len = newY
            u.br.len = newX - newY
            c.br.len = newM - newX
            if self.mcmc.constraints:
                tempRawSplitKey = v.br.rawSplitKey
                tempSplitKey = v.br.splitKey
                v.br.rawSplitKey = u.br.rawSplitKey
                v.br.splitKey = u.br.splitKey
                u.br.rawSplitKey = tempRawSplitKey
                u.br.splitKey = tempSplitKey
            pTree.preAndPostOrderAreValid = 0
            theProposal.topologyChanged = 1
            if dbug:
                print "-> detach v, re-attach between a and u, so there is a topology change"
                print "   (splitKeys are wrong on nodes v and u, in the figure below)"
                for n in pTree.iterNodes():
                    n.name = n.oldName
                for n in pTree.iterLeavesNoRoot():
                    n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
                for n in pTree.iterInternalsNoRoot():
                    n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)

                if v.name:
                    v.name += '(v)'
                else:
                    v.name = 'v'
                if u.name:
                    u.name += '(u)'
                else:
                    u.name = 'u'
                if c.name:
                    c.name += '(c)'
                else:
                    c.name = 'c'
                if d.name:
                    d.name += '(d)'
                else:
                    d.name = 'd'
                if a.name:
                    a.name += '(a)'
                else:
                    a.name = 'a'
                if b.name:
                    b.name += '(b)'
                else:
                    b.name = 'b'
                pTree.draw(width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)
                if self.mcmc.constraints:
                    pTree.checkSplitKeys(useOldName=True, glitch=False)
                

    
    # Check that new brLens are not too short or too long.  Ronquist
    # suggests that if that is the case, then just abort, rather than
    # fussing with reflections.
    if 0:
        if c.br.len < var.BRLEN_MIN or v.br.len < var.BRLEN_MIN or u.br.len < var.BRLEN_MIN:
            #if dbug:
            print "At least 1 brLen is too short."
            theProposal.doAbort = True
            return
        elif c.br.len > var.BRLEN_MAX or v.br.len > var.BRLEN_MAX or u.br.len > var.BRLEN_MAX:
            #if dbug: 
            print "At least 1 brLen is too long.  Aborting. (No big deal ...)"
            theProposal.doAbort = True
            return

    if 1:
        complain = True
        if c.br.len < var.BRLEN_MIN:
            if complain:
                print "c  %i  too short" % self.mcmc.gen
            theProposal.doAbort = True
            return
        if v.br.len < var.BRLEN_MIN:
            if complain:
                print "v  %i  too short" % self.mcmc.gen
            theProposal.doAbort = True
            return
        if u.br.len < var.BRLEN_MIN:
            if complain:
                print "u  %i  too short" % self.mcmc.gen
            theProposal.doAbort = True
            return
        if c.br.len > var.BRLEN_MAX:
            if complain:
                print "c  %i  too long" % self.mcmc.gen
            theProposal.doAbort = True
            return
        if v.br.len > var.BRLEN_MAX:
            if complain:
                print "v  %i  too long" % self.mcmc.gen
            theProposal.doAbort = True
            return
        if u.br.len > var.BRLEN_MAX:
            if complain:
                print "u  %i  too long" % self.mcmc.gen
            theProposal.doAbort = True
            return





    #print self.mcmc.constraints, theProposal.topologyChanged
    #raise Glitch

    if self.mcmc.constraints and theProposal.topologyChanged:
        ## Check whether any constraints have been involved, and if so abort.
        ##
        ## It was like this:
        ##                   +------c
        ##           +-------|(v)
        ##    +------|(u)    +------d
        ##    |      |
        ##    |(a)   +-------b
        ##    |
        ##    +------X
        ##
        ## And now its like this:
        ##                   +------c
        ##           +-------|(u)
        ##    +------|(v)    +------b
        ##    |      |
        ##    |(a)   +-------d
        ##    |
        ##    +------X
        #
        # So splits on branches on all nodes except u are unaffected.
        # The rawSplitKey on node v might be affected, if there was a
        # re-rooting, but the splitKey is still ok.  So the
        # rawSplitKey and the splitKey for node u in its new position
        # needs to be re-calculated.

        oldUSplitKey = u.br.splitKey
        pTree.recalculateSplitKeysOfNodeFromChildren(u, self.mcmc.constraints.allOnes)
        #print "u, node %i recalculated br.rawSplitKey=%s, br.splitKey = %s" % (
        #    u.nodeNum, u.br.rawSplitKey, u.br.splitKey)
        
        if dbug:
            for n in pTree.iterNodes():
                n.name = n.oldName
            pTree.checkTaxNames()
            for n in pTree.iterLeavesNoRoot():
                n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
            for n in pTree.iterInternalsNoRoot():
                n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
            pTree.draw(width=100, showInternalNodeNames=1, addToBrLen=0.2)
            if self.mcmc.constraints:
                pTree.checkSplitKeys(useOldName=True, glitch=False)  # If v is wrong, it gets righted later.

        for sk in self.mcmc.constraints.constraints:
            if oldUSplitKey == sk:
                #print "oldUSplitKey was constrained to %i, and now it is gone. --> abort!" % sk
                theProposal.doAbort = True
                if dbug:
                    for n in pTree.iterNodes():
                        n.name = n.oldName
                return        

    if var.doMcmcSp:
        c.br.lenChanged = True
        u.br.lenChanged = True
        v.br.lenChanged = True

            

    # Now evaluate the prior ratio.  If doInternalBrLenPrior is set,
    # then we do it like Yang and Rannala Syst Biol 54(3):455-470,
    # 2005.  Do it before reRooting (does it make any difference?)
    #print self.mcmc.tunings.doInternalBrLenPrior
    if hasattr(self.mcmc.tunings, 'doInternalBrLenPrior') and self.mcmc.tunings.doInternalBrLenPrior:
        # originally this
        ##                   +------c
        ##           +-------|(v)
        ##    +------|(u)    +------d
        ##    |      |
        ##    |(a)   +-------b
        ##    |
        ##    +------X (which might be None)
        # Either u or v was detached, and reattached somewhere between a and c.

        # See the "long way" calculation below, which shows how it can
        # be done with a single prior.  The complication here is that
        # there are two priors, depending on whether it is a leaf or
        # not.
        theSum = 0.0
        if newX < newY: # no topology change
            # Nodes a or c might be leaf nodes.  Nodes u and v cannot be.
            # Do the 3 edges in turn.  First the a edge, then the internal edge, then the c edge.
            if a.isLeaf:
                theSum += (self.mcmc.tunings.brLenPriorLambda * x) - (self.mcmc.tunings.brLenPriorLambda * newX)
            else:
                theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * x) - (
                    self.mcmc.tunings.brLenPriorLambdaForInternals * newX)
            theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * (y - x)) - (
                self.mcmc.tunings.brLenPriorLambdaForInternals * (newY - newX))
            if c.isLeaf:
                theSum += (self.mcmc.tunings.brLenPriorLambda * (m - y)) - (self.mcmc.tunings.brLenPriorLambda * (newM - newY))
            else:
                theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * (m - y)) - (
                    self.mcmc.tunings.brLenPriorLambdaForInternals * (newM - newY))
        else: # with topology change
            if a.isLeaf:
                theSum += (self.mcmc.tunings.brLenPriorLambda * x) - (self.mcmc.tunings.brLenPriorLambda * newY)
            else:
                theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * x) - (
                    self.mcmc.tunings.brLenPriorLambdaForInternals * newY)
            theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * (y - x)) - (
                self.mcmc.tunings.brLenPriorLambdaForInternals * (newX - newY))
            if c.isLeaf:
                theSum += (self.mcmc.tunings.brLenPriorLambda * (m - y)) - (self.mcmc.tunings.brLenPriorLambda * (newM - newX))
            else:
                theSum += (self.mcmc.tunings.brLenPriorLambdaForInternals * (m - y)) - (
                    self.mcmc.tunings.brLenPriorLambdaForInternals * (newM - newX))

        if 0:
            # Slow check, via priorDensities.
            theta1 = self.mcmc.tunings.brLenPriorLambda
            theta2 = self.mcmc.tunings.brLenPriorLambdaForInternals

            if newX < newY:
                if a.isLeaf:
                    theta = self.mcmc.tunings.brLenPriorLambda
                else:
                    theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                prDensNu1 = theta * math.exp(-theta * x)
                prDensNuStar1 = theta * math.exp(-theta * newX)

                theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                prDensNu2 = theta * math.exp(-theta * (y - x))
                prDensNuStar2 = theta * math.exp(-theta * (newY - newX))
                if c.isLeaf:
                    theta = self.mcmc.tunings.brLenPriorLambda
                else:
                    theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                prDensNu3 = theta * math.exp(-theta * (m - y))
                prDensNuStar3 = theta * math.exp(-theta * (newM - newY))
            else:
                if a.isLeaf:
                    theta = self.mcmc.tunings.brLenPriorLambda
                else:
                    theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                prDensNu1 = theta * math.exp(-theta * x)
                prDensNuStar1 = theta * math.exp(-theta * newY)

                theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                prDensNu2 = theta * math.exp(-theta * (y - x))
                prDensNuStar2 = theta * math.exp(-theta * (newX - newY))
                if c.isLeaf:
                    theta = self.mcmc.tunings.brLenPriorLambda
                else:
                    theta = self.mcmc.tunings.brLenPriorLambdaForInternals
                prDensNu3 = theta * math.exp(-theta * (m - y))
                prDensNuStar3 = theta * math.exp(-theta * (newM - newX))
            prRat = (prDensNuStar1 * prDensNuStar2 * prDensNuStar3) / (prDensNu1 * prDensNu2 * prDensNu3)
            logPrRat = math.log(prRat)

            if math.fabs(logPrRat - theSum) > 1.e-10:
                print "xxzz differs.  logPrRat=%g, theSum=%g" % (logPrRat, theSum)
            #else:
            #    print "s",
            

        self.logPriorRatio = theSum
            
        
    else:  # Do not doInternalBrLenPrior
        if self.mcmc.tunings.brLenPriorType == 'uniform':
            self.logPriorRatio = 0.0
        elif self.mcmc.tunings.brLenPriorType == 'exponential':
            self.logPriorRatio = self.mcmc.tunings.brLenPriorLambda * (m - newM)
            if 0:  # Do the same calculation the long way, edge by edge.
                #print "logPriorRatio = %+.4f" % self.logPriorRatio,
                foo0 = (self.mcmc.tunings.brLenPriorLambda * m) - (self.mcmc.tunings.brLenPriorLambda * newM)
                #print "%+.4f" % foo0,

                foo = 0.0
                if newX < newY: # no topology change
                    ##           newX    newY   newM
                    ##           +       +      +
                    ##
                    ##                   +------c
                    ##           +-------|(v)
                    ##    +------|(u)    +------d
                    ##    |      |
                    ##    |(a)   +-------b
                    ##    |
                    ##    +------X
                    foo += (self.mcmc.tunings.brLenPriorLambda * x) - (self.mcmc.tunings.brLenPriorLambda * newX)
                    foo += (self.mcmc.tunings.brLenPriorLambda * (y - x)) - (self.mcmc.tunings.brLenPriorLambda * (newY - newX))
                    foo += (self.mcmc.tunings.brLenPriorLambda * (m - y)) - (self.mcmc.tunings.brLenPriorLambda * (newM - newY))
                else: # with topology change
                    ##           newY    newX   newM
                    ##           +       +      +
                    ##
                    ##                   +------c
                    ##           +-------|(u)
                    ##    +------|(v)    +------b
                    ##    |      |
                    ##    |(a)   +-------d
                    ##    |
                    ##    +------X
                    foo += (self.mcmc.tunings.brLenPriorLambda * x) - (self.mcmc.tunings.brLenPriorLambda * newY)
                    foo += (self.mcmc.tunings.brLenPriorLambda * (y - x)) - (self.mcmc.tunings.brLenPriorLambda * (newX - newY))
                    foo += (self.mcmc.tunings.brLenPriorLambda * (m - y)) - (self.mcmc.tunings.brLenPriorLambda * (newM - newX))
                #print "%+.4f" % foo
                if (math.fabs(self.logPriorRatio - foo0) > 1.e-10):
                    print "differs-- foo0, %g %g" % (self.logPriorRatio, foo0)
                if (math.fabs(self.logPriorRatio - foo) > 1.e-10):
                    print "differs-- foo, %g %g" % (self.logPriorRatio, foo)
        else:
            raise Glitch, "This should not happen."

    if oldRoot:
        if dbug:
            print '-------------------- about to reRoot -----------'
            pTree.draw(width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)
            if self.mcmc.constraints:
                pTree.checkSplitKeys(useOldName=True, glitch=False)
            for n in pTree.iterNodes():
                n.name = n.oldName
            
        pTree.reRoot(oldRoot, moveInternalName=False, fixRawSplitKeys=self.mcmc.constraints)

        if dbug:
            print '--------------after reRoot --------------'
            for n in pTree.iterLeavesNoRoot():
                n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
            for n in pTree.iterInternalsNoRoot():
                n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
            pTree.draw(width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)
            if self.mcmc.constraints:
                pTree.checkSplitKeys(useOldName=True)
            

        if self.mcmc.constraints and theProposal.topologyChanged:
            # Node v in its new position might now need to have its
            # rawSplitKey updated.  The splitKey should still be ok.
            children = [n for n in v.iterChildren()]
            #print "leaves above v, node %i: %s" % (v.nodeNum, [n.nodeNum for n in children])
            #for n in children:
            #    print "    %2i  %4i  %4i" % (n.nodeNum, n.br.rawSplitKey, n.br.splitKey)
            x = children[0].br.rawSplitKey
            #print "rawSplitKeys: ", x,
            for n in children[1:]:
                y = n.br.rawSplitKey
                #print y,
                x = x | y  # '|' is bitwise "OR".
            v.br.rawSplitKey = x
            #print '=>', x
            #print "v, node %i recalculated br.rawSplitKey=%s" % (v.nodeNum, v.br.rawSplitKey)
            
    if dbug:
        if theProposal.topologyChanged:
            print "The topology CHANGED"
        else:
            print "Topology -- no change."
        #pTree.draw(width=80)
        
        for n in pTree.iterNodes():
            n.name = n.oldName
        pTree.checkTaxNames()
        for n in pTree.iterLeavesNoRoot():
            n.name += '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
        for n in pTree.iterInternalsNoRoot():
            n.name = '[%s,%s]' % (n.br.rawSplitKey, n.br.splitKey)
        pTree.draw(width=100, showInternalNodeNames=1, addToBrLen=0.2, model=True)
        if self.mcmc.constraints:
            pTree.checkSplitKeys(useOldName=True)
        #c.name = c.oldName
        #u.name = u.oldName
        #v.name = v.oldName
        #d.name = d.oldName
        #a.name = a.oldName
        #b.name = b.oldName
        if 0:
            del(c.oldName)
            del(u.oldName)
            del(v.oldName)
            del(d.oldName)
            del(a.oldName)
            del(b.oldName)
        for n in pTree.iterNodes():
            n.name = n.oldName
        


    self.logProposalRatio = 3.0 * math.log(newMRatio)

    if 0:
        for n in pTree.iterNodesNoRoot():
            if n.br.lenChanged:
                print "l node %2i br.lenChanged" % n.nodeNum
    if dbug:
        if self.mcmc.constraints:
            pTree.checkSplitKeys()
    #if self.mcmc.constraints:
    #    print "checkSplitKeys() at the end of local"
    #    pTree.checkSplitKeys()

    # Check if we have a new combo of comp and rMatrix
    if theProposal.topologyChanged:
        for pNum in range(pTree.model.nParts):
            if 0 and self.mcmc.gen == 14:
                print
                print pTree.model.parts[pNum].bQETneedsReset
                print "a is node %i" % a.nodeNum
                print "u is node %i" % u.nodeNum
                print "v is node %i" % v.nodeNum
            for n in [a, u, v]:
                if n.br:
                    if pTree.model.parts[pNum].bQETneedsReset[n.parts[pNum].compNum][n.br.parts[pNum].rMatrixNum]:
                        #print "bQETneedsReset is set for %i, %i." % (
                        #    n.parts[pNum].compNum, n.br.parts[pNum].rMatrixNum)
                        pf.p4_resetBQET(pTree.model.cModel, pNum, n.parts[pNum].compNum, n.br.parts[pNum].rMatrixNum)
            if 0 and self.mcmc.gen == 14:
                print
                print pTree.model.parts[pNum].bQETneedsReset
                

    if dbug:
        for n in pTree.iterNodesNoRoot():
            if math.fabs(n.br.len - n.br.oldLen) > 0.0000001:
                #print "Node %2i br len changed" % n.nodeNum
                assert n.br.lenChanged
            else:
                if n.br.lenChanged:
                    print "Node %2i lenChanged set, but its the same length." % n.nodeNum
                    raise Glitch

            # If the branch has been inverted, we will want to recalculate
            # the bigPDecks, even if the length has not really changed.
            # Trigger that intent by setting lenChanged.
            if n.br.oldNode != n:
                #print "Node %2i branch: oldNode %2i" % (n.nodeNum, n.br.oldNode.nodeNum)
                n.br.lenChanged = 1



def proposeETBR_Blaise(self, theProposal):
    """Adapted and modified from Jason Evans' excellent Crux v 1.2.0

    Many thanks to JE, who wrote such clear code.  The Lakner et al
    version was modified by JE for Crux so that it works on
    polytomies, and so that it works on leaf nodes.

    Many thanks to Blaise Li who convincingly pointed out that LOCAL
    was not enough, and then pointed out that even the Crux version of
    eTBR did not appear to be reversible when there were polytomies.
    This version has a suggested modification from Blaise as described
    in his poster 'An eTBR proposal for non-binary trees in MCMC
    Bayesian phylogeny', B Li and P Foster, presented (by BL) at the
    15th Evolutionary Biology Meeting, September 27-30, 2011
    Marseilles.

    """

    # doAbort is set if brLens are too long or too short, or if a constraint is violated.
    gm = ['Chain.proposeETBR_Blaise()']

    # Now works with constraints.
    theProposal.topologyChanged = 0
    theProposal.doAbort = False
    pTree = self.propTree
    dbug = False

    for n in pTree.iterNodesNoRoot():
        n.br.oldLen = n.br.len
        n.br.oldNode = n
        
    if 0 and self.mcmc.gen == 217:
        dbug = True
        print "------------- eTBR gen %i -----------------" % self.mcmc.gen
        if 0:
            currentLogLike = self.propTree.logLike
            self.propTree.calcLogLike(verbose=0)  # with _commonCStuff()
            theDiff = math.fabs(currentLogLike - self.propTree.logLike)
            if theDiff > 1.e-9:
                gm.append("propTree like diff %f (%g)" % (theDiff, theDiff))
                raise Glitch, gm

    oldRoot = pTree.root
    eA = None
    eX = None
    eY = None

    etbrLambda = theProposal.tuning
    etbrPExt = self.mcmc.tunings.etbrPExt

    if 0 and dbug:
        print "=" * 50
        pTree.draw()
        print "starting with the tree above."

    # Choose a node, not the root.  It will have edge eA in Jason's diagram.
    # y0 will be the asterisk node in Jason's diagram.  It may be extended, below.
    y0 = None
    while not y0:
        nNum = random.choice(pTree.preOrder)
        if nNum != var.NO_ORDER and nNum != pTree.root.nodeNum:
            y0 = pTree.node(nNum)
    x0 = y0.parent
    if dbug:
        print "y0 is node %i" % y0.nodeNum
        y0.br.textDrawSymbol = '='
        if y0.name:
            y0.name += '_y0'
        else:
            y0.name = 'y0'
        if x0.name:
            x0.name += '_x0'
        else:
            x0.name = 'x0'

    # Name the edge here, used below when we modify the br.len
    eA = y0.br
            

    ################## Extend x
    # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
    r1 = x0  

    #  If x0 is not a leaf, it is 'unconstrained', and since x0 is the
    #  parent of y0, it will always be so in p4.
    x0Degree = pTree.getDegree(x0)
    x0Uncon = (x0Degree > 1)  # I think this will always be so.
    if x0Uncon:
        myRan = random.randrange(x0Degree - 1)
        r0 = pTree.nextNode(y0, x0)
        for i in range(myRan):
            r0 = pTree.nextNode(r0, x0)
        if r0 == x0:
            r0 = r0.parent
        x1 = r0
        if dbug:
            if x1.name:
                x1.name += '_x1'
            else:
                x1.name = 'x1'

        # So we are set up like this ...
        #   nR1      nR0
        #   nX0==eX==nX1

        # We name the edge here, and use it below when we modify the br.len
        if x0.parent == x1:
            eX = x0.br
        elif x1.parent == x0:
            eX = x1.br
        else:
            raise Glitch, "This should not happen"

        #  if r0 is not a leaf, it is 'unconstrained'
        r0Degree = pTree.getDegree(r0)
        r0Uncon = (r0Degree > 1)
        while r0Uncon:
            # randomly determine whether to extend
            myRan = random.random()
            if etbrPExt < myRan:
                break
            # choose an edge (direction) to go, and step to the next node
            myRan = random.randrange(r0Degree - 1)
            # We now call nextNode(spoke, hub) where hub is r0.  If r0
            # is the parent of r1, then the spoke of nextNode(spoke,
            # hub) is r1.  But if its the other way around, if r1 is
            # the parent of r0, then spoke is r0.
            if r1.parent == r0:
                r0new = pTree.nextNode(r1, r0)
            elif r0.parent == r1:
                r0new = pTree.nextNode(r0, r0)
            else:
                gm.append("This shouldn't happen.")
                raise Glitch, gm
            
            for i in range(myRan):
                r0new = pTree.nextNode(r0new, r0)
            if r0new == r0:
                r0new = r0.parent
            r1 = r0
            r0 = r0new
            r0Degree = pTree.getDegree(r0)
            r0Uncon = (r0Degree > 1)
        if dbug:
            if r1.name:
                r1.name += '_r1'
            else:
                r1.name = 'r1'
            if r0.name:
                r0.name += '_r0'
            else:
                r0.name = 'r0'            
            #pTree.draw()
            

        # Perform rearrangement unless it would be a no-op.  It would be a
        # no-op if r0 was still x1.
        if r0 == x1:
            # We did not extend, at all.
            if dbug:
                print "No extension from x1 was done (because x1=r0), so no rearrangement on the x side."
            pass
        else:
            # Do the rearrangement.  We need to have x0 above x1, and
            # r0 above r1.  We could guarrantee that by re-rooting on
            # r1, but that may not be needed.
            if x1.parent == x0:
                pTree.reRoot(x1, moveInternalName=False)
            elif x0.parent == x1:
                pass
            else:
                gm.append("This shouldn't happen.")
                raise Glitch, gm
            if dbug:
                x0.br.textDrawSymbol = 'X'
                
            # The r1-r0 branch might be pointing up or down, but we need
            # it to be pointing up, such that r0 is the child of r1 and so
            # r0.parent = r1.  If it is not, then we need to re-root.

            if r0.parent == r1:
                # no need to re-root, the r branch points up
                pass
            elif r1.parent == r0:
                pTree.reRoot(r1, moveInternalName=False)
            else:
                gm.append("This shouldn't happen.")
                raise Glitch, gm

            assert r0.parent == r1
            assert x0.parent == x1

            if dbug:
                pTree.draw()
                print "The drawing above is just before the rearrangement."

            # Get children of x0 that are not y0
            ch = [n for n in x0.iterChildren() if n != y0]
            assert ch
            if len(ch) == 1:
                myChoice = ch[0]
            else:
                myChoice = random.choice(ch)

            # Three spr moves.
            pTree.pruneSubTreeWithoutParent(myChoice, allowSingleChildNode=True)
            pTree.reconnectSubTreeWithoutParent(myChoice, x1)

            if x0 != r1:
                pTree.pruneSubTreeWithoutParent(x0) 
                pTree.reconnectSubTreeWithoutParent(x0, r1)

            pTree.pruneSubTreeWithoutParent(r0) 
            pTree.reconnectSubTreeWithoutParent(r0, x0)

    if 1 and dbug:
        pTree.setPreAndPostOrder()
        pTree.draw()
        if x0Uncon and r0 != x1:
            print "The drawing above shows that X extended"
        else:
            print "The drawing above shows that X did not extend."

    ################# Extend y
    # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
    s1 = y0  

    #  If y0 is not a leaf, it is 'unconstrained'
    y0Degree = pTree.getDegree(y0)
    y0Uncon = (y0Degree > 1)  # y0 will sometimes be a leaf
    if y0Uncon:
        myRan = random.randrange(y0Degree - 1)
        s0 = pTree.nextNode(y0, y0)
        for i in range(myRan):
            s0 = pTree.nextNode(s0, y0)
        assert s0.parent == y0
        y1 = s0
        if dbug:
            if y1.name:
                y1.name += '_y1'
            else:
                y1.name = 'y1'

        # name the edge here, to be used below when we modify the br.len
        eY = y1.br

        # So we are set up like this ...
        #   nS1      nS0
        #   nY0==eY==nY1
        
        #  if s0 is not a leaf, it is 'unconstrained'
        s0Degree = pTree.getDegree(s0)
        s0Uncon = (s0Degree > 1)
        while s0Uncon:
            # randomly determine whether to extend
            myRan = random.random()
            if etbrPExt < myRan:
                break
            # choose an edge (direction) to go, and step to the next node
            myRan = random.randrange(s0Degree - 1)
            # We now call nextNode(spoke, hub) where hub is s0.  If s0
            # is the parent of s1, then the spoke of nextNode(spoke,
            # hub) is s1.  But if its the other way around, if s1 is
            # the parent of s0, then spoke is s0.

            # Since we are extending Y here, and this entire subtree goes up, it should always be that s0.parent is s1.
            if s1.parent == s0:
                gm.append("s1.parent is s0.  This should not happen.")
                raise Glitch, gm
            elif s0.parent == s1:
                s0new = pTree.nextNode(s0, s0)
            else:
                gm.append("This shouldn't happen.")
                raise Glitch, gm
            
            for i in range(myRan):
                s0new = pTree.nextNode(s0new, s0)

            assert s0new != s0
            #if s0new == s0:
            #    s0new = s0.parent
            s1 = s0
            s0 = s0new
            s0Degree = pTree.getDegree(s0)
            s0Uncon = (s0Degree > 1)
        if dbug:
            if s1.name:
                s1.name += '_s1'
            else:
                s1.name = 's1'
            if s0.name:
                s0.name += '_s0'
            else:
                s0.name = 's0'
            y1.br.textDrawSymbol = 'Y'
            #pTree.draw()

        # Perform rearrangement unless it would be a no-op.  It would be a
        # no-op if s0 was still y1.
        if s0 == y1:
            # We did not extend, at all.
            if dbug:
                print "No Y extension was made (because y1=s0), so nothing to do."
        else:
            # Do the rearrangement.  Since it is Y that is being
            # extended, it should always be pointing such that
            # y1.parent is y0.
            assert y1.parent == y0

            # Because we are extending Y, the s1-s0 branch should always be
            # be pointing such that s0 is the child of s1 and so
            # s0.parent = s1.
            assert s0.parent == s1

            # We want y0.parent to be y1, and we want s0.parent to be
            # s1; and only the latter is currently true.  So reRoot to
            # y1.

            pTree.reRoot(y1, moveInternalName=False)
            assert y0.parent == y1
            assert s0.parent == s1

            if dbug:
                pTree.draw()
                print "The drawing above is the tree just before rearrangement on the Y side."

            # Get children of y0 that are not x0
            ch = [n for n in y0.iterChildren() if n != x0]
            assert ch
            if len(ch) == 1:
                myChoice = ch[0]
            else:
                myChoice = random.choice(ch)

            # Three spr moves.
            pTree.pruneSubTreeWithoutParent(myChoice, allowSingleChildNode=True) # myChoice is returned, also
            pTree.reconnectSubTreeWithoutParent(myChoice, y1)

            if y0 != s1:
                pTree.pruneSubTreeWithoutParent(y0) 
                pTree.reconnectSubTreeWithoutParent(y0, s1)
                
            pTree.pruneSubTreeWithoutParent(s0) 
            pTree.reconnectSubTreeWithoutParent(s0, y0)


    if dbug:
        pTree.setPreAndPostOrder()
        pTree.draw()
        if y0Uncon and s0 != y1:
            print "The drawing above shows that Y extended"
        else:
            print "The drawing above shows that Y did not extend."

    if oldRoot != pTree.root:
        pTree.reRoot(oldRoot, moveInternalName=False)
    if dbug:
        pTree.draw()
        print "The above is back to the original root."

    if dbug:
        for n in pTree.nodes:
            if n.isLeaf:
                while '_' in n.name:
                    n.name = n.name[:-3]
            elif not n.isLeaf:
                n.name = None
            if n.br:
                n.br.textDrawSymbol = '-'

    xRearranged = False
    if x0Uncon:
        if r0 != x1:
            xRearranged = True
    yRearranged = False
    if y0Uncon:
        if s0 != y1:
            yRearranged = True
    if xRearranged or yRearranged:
        theProposal.topologyChanged = True

    if dbug:
        print "-" * 20
        print "xRearranged = %s" % xRearranged
        print "yRearranged = %s" % yRearranged
        if eA:
            for n in pTree.iterNodesNoRoot():
                if n.br == eA:
                    print "eA is from node", n.nodeNum
        else:
            print "eA is None"
        if eX:
            for n in pTree.iterNodesNoRoot():
                if n.br == eX:
                    print "eX is from node", n.nodeNum
        else:
            print "eX is None"
        if eY:
            for n in pTree.iterNodesNoRoot():
                if n.br == eY:
                    print "eY is from node", n.nodeNum
        else:
            print "eY is None"
        print "x0Uncon is", x0Uncon
        print "y0Uncon is", y0Uncon

    # Are we violating constraints?
    if 1:
        if self.mcmc.constraints and theProposal.topologyChanged:
            pTree.makeSplitKeys()
            pTreeSKSet = set([n.br.splitKey for n in pTree.iterInternalsNoRoot()])
            isViolating = False
            for sk in self.mcmc.constraints.constraints:
                if sk not in pTreeSKSet:
                    isViolating = True
                    break
            if isViolating:
                theProposal.doAbort = True
                return

    # n.flag is set if the condLikes need recalculating.  Edges eA,
    # eX, and eY will have their bigPDecks recalculated, and all the
    # nodes below those will have their flags set (in
    # Chain.proposeSp()).  However, there are additional nodes that
    # need the flag set.
    if xRearranged:
        p = x1
        while 1:
            p.flag = 1
            p = p.parent
            if not p:
                break
    if yRearranged:
        p = y1
        while 1:
            p.flag = 1
            p = p.parent
            if not p:
                break

    
    # Do the branch length changes.  Taken nearly verbatim from Crux.
    # Thanks Jason!  The 3 edges to modify are eA, eX, and eY, which
    # were named above, as NodeBranch objects.
    
    lnProp = 0.0

    # Generate branch length multipliers and set new branch lengths.
    # eA.
    u = random.random()
    lnMA = etbrLambda * (u - 0.5)
    lnProp += lnMA
    mA = math.exp(lnMA)
    vA0 = eA.len
    vA1 = vA0 * mA
    eA.len = vA1
    eA.lenChanged = True
    # eX.
    if x0Uncon:
        u = random.random()
        lnMX = etbrLambda * (u - 0.5)
        lnProp += lnMX
        mX = math.exp(lnMX)
        vX0 = eX.len
        vX1 = vX0 * mX
        eX.len = vX1
        eX.lenChanged = True
    # eY.
    if y0Uncon:
        u = random.random()
        lnMY = etbrLambda * (u - 0.5)
        lnProp += lnMY
        mY = math.exp(lnMY)
        vY0 = eY.len
        vY1 = vY0 * mY
        eY.len = vY1
        eY.lenChanged = True
   
    # More from Crux ...
    # The prior ratio is the product of the prior ratios for each modified
    # branch length.  The number of internal branches does not change
    # (though the number of polytomies may change), so the topology prior
    # ratio is always 1.
    if self.mcmc.tunings.brLenPriorType == 'exponential':
        lnPrior = -self.mcmc.tunings.brLenPriorLambda * (vA1-vA0)
        if x0Uncon:
            lnPrior += -self.mcmc.tunings.brLenPriorLambda * (vX1-vX0)
        if y0Uncon:
            lnPrior += -self.mcmc.tunings.brLenPriorLambda * (vY1-vY0)
    elif self.mcmc.tunings.brLenPriorType == 'uniform':
        lnPrior = 0.0
        
    # The proposal ratio is the product of the proposal ratios for
    # extension of each end of eA, as well as the branch multipliers.  The
    # ratio is 1 for the constrained/constrained and
    # unconstrained/unconstrained extension cases.
    #
    # nY0/nR0.
    if x0Uncon and r0 is not x1:
        if y0Uncon:
            if not r0Uncon:
                lnProp += math.log(1.0 - etbrPExt)
        elif r0Uncon:
            lnProp += math.log(1.0 / (1.0 - etbrPExt))
    # nX0/nS0.
    if y0Uncon and s0 is not y1:
        if x0Uncon:
            if not s0Uncon:
                lnProp += math.log(1.0 - etbrPExt)
        elif s0Uncon:
            lnProp += log(1.0 / (1.0 - etbrPExt))

    self.propTree.preAndPostOrderAreValid = False

    self.logPriorRatio = lnPrior
    self.logProposalRatio = lnProp
    
    # Check if we have a new combo of comp and rMatrix.  This might be
    # more efficient if I cleverly only look at the affected nodes.
    if 1:
        if theProposal.topologyChanged:
            for n in pTree.iterNodesNoRoot():
                for pNum in range(pTree.model.nParts):
                    theCompNum = n.parts[pNum].compNum
                    theRMatrixNum = n.br.parts[pNum].rMatrixNum
                    if pTree.model.parts[pNum].bQETneedsReset[theCompNum][theRMatrixNum]:
                        pf.p4_resetBQET(pTree.model.cModel, pNum, theCompNum, theRMatrixNum)


    # This stuff below could probably be done more cleverly, but this works.
    for n in pTree.iterNodesNoRoot():
        if math.fabs(n.br.len - n.br.oldLen) > 0.0000001:
            #print "Node %2i br len changed" % n.nodeNum
            assert n.br.lenChanged
        #else:
        #    if n.br.lenChanged:
        #        print "Node %2i lenChanged set, but its the same length." % n.nodeNum
        #        raise Glitch

        # If the branch has been inverted, we will want to recalculate
        # the bigPDecks, even if the length has not really changed.
        # Trigger that intent by setting lenChanged.
        if n.br.oldNode != n:
            #print "Node %2i branch: oldNode %2i" % (n.nodeNum, n.br.oldNode.nodeNum)
            n.br.lenChanged = 1

def proposeETBR(self, theProposal):
    """Adapted from Jason Evans' excellent Crux v 1.2.0

    Many thanks to JE, who wrote such clear code.  The Crux version
    came from Lakner et al, but was modified by JE so that it works on
    polytomies, and so that it works on leaf nodes.  Also many thanks
    to Blaise Li who convincingly pointed out that LOCAL was not enough.

    It does not work with constraints yet.
    """

    # doAbort is set if brLens are too long or too short, or if a constraint is violated.
    gm = ['Chain.proposeETBR()']

    if self.mcmc.constraints:
        gm.append("Sorry, due to lazy programming, proposeETBR() does not work with constraints yet.")
        raise Glitch, gm
    theProposal.topologyChanged = 0
    theProposal.doAbort = False
    pTree = self.propTree
    dbug = False

    if 0 and self.mcmc.gen == 404:
        dbug = True
        if 0:
            currentLogLike = self.propTree.logLike
            self.propTree.calcLogLike(verbose=0)  # with _commonCStuff()
            theDiff = math.fabs(currentLogLike - self.propTree.logLike)
            if theDiff > 1.e-9:
                gm.append("propTree like diff %f (%g)" % (theDiff, theDiff))
                raise Glitch, gm

    oldRoot = pTree.root
    eA = None
    eX = None
    eY = None

    etbrLambda = theProposal.tuning
    etbrPExt = self.mcmc.tunings.etbrPExt

    if 1 and dbug:
        print "=" * 50
        #pTree.draw()
        #print "starting with the tree above."

    # Choose a node, not the root.  It will have edge eA in Jason's diagram.
    # y0 will be the asterisk node in Jason's diagram.  It may be extended, below.
    y0 = None
    while not y0:
        nNum = random.choice(pTree.preOrder)
        if nNum != var.NO_ORDER and nNum != pTree.root.nodeNum:
            y0 = pTree.node(nNum)
    x0 = y0.parent
    if dbug:
        print "y0 is node %i" % y0.nodeNum
        y0.br.textDrawSymbol = '='
        if y0.name:
            y0.name += '_y0'
        else:
            y0.name = 'y0'
        if x0.name:
            x0.name += '_x0'
        else:
            x0.name = 'x0'

    # Name the edge here, used below when we modify the br.len
    eA = y0.br
            

    ################## Extend x
    if dbug:
        xUp = None
    # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
    r1 = x0  

    #  If x0 is not a leaf, it is 'unconstrained', and since x0 is the
    #  parent of y0, it will always be so in p4.
    x0Degree = pTree.getDegree(x0)
    x0Uncon = (x0Degree > 1)  # I think this will always be so.
    if x0Uncon:
        myRan = random.randrange(x0Degree - 1)
        r0 = pTree.nextNode(y0, x0)
        for i in range(myRan):
            r0 = pTree.nextNode(r0, x0)
        if r0 == x0:
            r0 = r0.parent
        x1 = r0
        if dbug:
            if x1.name:
                x1.name += '_x1'
            else:
                x1.name = 'x1'

        # So we are set up like this ...
        #   nR1      nR0
        #   nX0==eX==nX1

        # We name the edge here, and used it below when we modify the br.len
        if x0.parent == x1:
            eX = x0.br
        elif x1.parent == x0:
            eX = x1.br
        else:
            raise Glitch, "This should not happen"

        #  if r0 is not a leaf, it is 'unconstrained'
        r0Degree = pTree.getDegree(r0)
        r0Uncon = (r0Degree > 1)
        while r0Uncon:
            # randomly determine whether to extend
            myRan = random.random()
            if etbrPExt < myRan:
                break
            # choose an edge (direction) to go, and step to the next node
            myRan = random.randrange(r0Degree - 1)
            # We now call nextNode(spoke, hub) where hub is r0.  If r0
            # is the parent of r1, then the spoke of nextNode(spoke,
            # hub) is r1.  But if its the other way around, if r1 is
            # the parent of r0, then spoke is r0.
            if r1.parent == r0:
                r0new = pTree.nextNode(r1, r0)
            elif r0.parent == r1:
                r0new = pTree.nextNode(r0, r0)
            else:
                gm.append("This shouldn't happen.")
                raise Glitch, gm
            
            for i in range(myRan):
                r0new = pTree.nextNode(r0new, r0)
            if r0new == r0:
                r0new = r0.parent
            r1 = r0
            r0 = r0new
            r0Degree = pTree.getDegree(r0)
            r0Uncon = (r0Degree > 1)
        if dbug:
            if r1.name:
                r1.name += '_r1'
            else:
                r1.name = 'r1'
            if r0.name:
                r0.name += '_r0'
            else:
                r0.name = 'r0'            
        if dbug:
            pTree.draw()

        # Perform rearrangement unless it would be a no-op.  It would be a
        # no-op if r0 was still x1.
        if r0 == x1:
            # We did not extend, at all.
            #print "No extension from x1 was done, so no rearrangement on the x side."
            pass
        else:
            # Do the rearrangement.  We want to collapse edge x1-x0, and so we
            # need to know if it is pointing up or down.
            if x1.parent == x0:
                xUp = x1
                xDn = x0
            elif x0.parent == x1:
                xUp = x0
                xDn = x1
            else:
                gm.append("This shouldn't happen.")
                raise Glitch, gm
            assert xUp.parent == xDn
            if dbug:
                #xUp.name += '_xUp'
                #xDn.name += '_xDn'
                xUp.br.textDrawSymbol = 'X'
                
            # The r1-r0 branch might be pointing up or down, but we need
            # it to be pointing up, such that r0 is the child of r1 and so
            # r0.parent = r1.  If it is not, then we need to re-root.

            if r0.parent == r1:
                # no need to re-root, the r branch points up
                pass
            elif r1.parent == r0:
                pTree.reRoot(r1, moveInternalName=False)
            else:
                gm.append("This shouldn't happen.")
                raise Glitch, gm
            assert r0.parent == r1
            assert xUp.parent == xDn

            # Collapse the node between xUp and xDn, ...
            if r1 == xUp:
                r1 = xDn
                if dbug:
                    # fix the internal node names
                    x1.name = x1.name[:-3]
                    x0.name += "_r1"
            theRightmostChild = xUp.rightmostChild()
            theLeftSib = xUp.leftSibling()
            if theLeftSib:
                theLeftSib.sibling = xUp.leftChild
            else:
                xDn.leftChild = xUp.leftChild
            for n in xUp.iterChildren():
                n.parent = xDn
            theRightmostChild.sibling = xUp.sibling
            xUp.wipe()  # needed?
            assert r0.parent == r1

            # ... and insert that node (xUp) between r0 and r1.
            xUp.parent = r1
            xUp.leftChild = r0
            r0.parent = xUp
            if r1.leftChild == r0:
                r1.leftChild = xUp
            else:
                oldCh = r1.leftChild
                while oldCh.sibling != r0:
                    oldCh = oldCh.sibling
                oldCh.sibling = xUp
            if r0.sibling:
                xUp.sibling = r0.sibling
                r0.sibling = None

            if dbug:
                pTree.setPreAndPostOrder()
                pTree.draw()

            pTree.pruneSubTreeWithoutParent(y0) # this method returns y0 as well.
            pTree.reconnectSubTreeWithoutParent(y0, xUp)

            #if oldRoot != pTree.root:
            #    pTree.reRoot(oldRoot, moveInternalName=False)

    if 1 and dbug:
        pTree.setPreAndPostOrder()
        pTree.draw()
        if x0Uncon and r0 != x1:
            print "The drawing above shows that X extended"
        else:
            print "The drawing above shows that X did not extend."

    ################# Extend y
    # x0, x1, y0, and y1 do not change, but r0, r1, s0, and s1 change.
    s1 = y0  

    #  If y0 is not a leaf, it is 'unconstrained'
    y0Degree = pTree.getDegree(y0)
    y0Uncon = (y0Degree > 1)  # y0 will sometimes be a leaf
    if y0Uncon:
        myRan = random.randrange(y0Degree - 1)
        s0 = pTree.nextNode(y0, y0)
        for i in range(myRan):
            s0 = pTree.nextNode(s0, y0)
        assert s0.parent == y0
        y1 = s0
        if dbug:
            if y1.name:
                y1.name += '_y1'
            else:
                y1.name = 'y1'

        # name the edge here, to be used below when we modify the br.len
        eY = y1.br

        # So we are set up like this ...
        #   nS1      nS0
        #   nY0==eY==nY1
        #  if s0 is not a leaf, it is 'unconstrained'
        s0Degree = pTree.getDegree(s0)
        s0Uncon = (s0Degree > 1)
        while s0Uncon:
            # randomly determine whether to extend
            myRan = random.random()
            if etbrPExt < myRan:
                break
            # choose an edge (direction) to go, and step to the next node
            myRan = random.randrange(s0Degree - 1)
            # We now call nextNode(spoke, hub) where hub is s0.  If s0
            # is the parent of s1, then the spoke of nextNode(spoke,
            # hub) is s1.  But if its the other way around, if s1 is
            # the parent of s0, then spoke is s0.

            # Since we are extending Y here, and this entire subtree goes up, it should always be that s0.parent is s1.
            if s1.parent == s0:
                gm.append("s1.parent is s0.  This should not happen.")
                raise Glitch, gm
            elif s0.parent == s1:
                s0new = pTree.nextNode(s0, s0)
            else:
                gm.append("This shouldn't happen.")
                raise Glitch, gm
            
            for i in range(myRan):
                s0new = pTree.nextNode(s0new, s0)

            assert s0new != s0
            #if s0new == s0:
            #    s0new = s0.parent
            s1 = s0
            s0 = s0new
            s0Degree = pTree.getDegree(s0)
            s0Uncon = (s0Degree > 1)
        if dbug:
            if s1.name:
                s1.name += '_s1'
            else:
                s1.name = 's1'
            if s0.name:
                s0.name += '_s0'
            else:
                s0.name = 's0'
            y1.br.textDrawSymbol = 'Y'
            pTree.draw()

        # Perform rearrangement unless it would be a no-op.  It would be a
        # no-op if s0 was still y1.
        if s0 == y1:
            # We did not extend, at all.
            pass
        else:
            # Do the rearrangement.  We want to collapse edge y1-y0.
            # Since it is Y that is being extended, it should always
            # be pointing such that y1.parent is y0.
            assert y1.parent == y0

            # Because we are extending Y, the s1-s0 branch should always be
            # be pointing such that s0 is the child of s1 and so
            # s0.parent = s1.
            assert s0.parent == s1

            # We can't easily prune off the subtree downwards, so reRoot to y0
            pTree.reRoot(y0, moveInternalName=False)
            if dbug:
                pTree.draw()

            # Collapse the node between y1 and y0, ...
            if s1 == y1:
                s1 = y0
                if dbug:
                    # fix the internal node names
                    y1.name = y1.name[:-3]
                    y0.name += "_s1"
            theRightmostChild = y1.rightmostChild()
            theLeftSib = y1.leftSibling()
            if theLeftSib:
                theLeftSib.sibling = y1.leftChild
            else:
                y0.leftChild = y1.leftChild
            for n in y1.iterChildren():
                n.parent = y0
            theRightmostChild.sibling = y1.sibling
            y1.wipe()  # needed?
            assert s0.parent == s1

            # ... and insert that node (y1) between s0 and s1.
            y1.parent = s1
            y1.leftChild = s0
            s0.parent = y1
            if s1.leftChild == s0:
                s1.leftChild = y1
            else:
                oldCh = s1.leftChild
                while oldCh.sibling != s0:
                    oldCh = oldCh.sibling
                oldCh.sibling = y1
            if s0.sibling:
                y1.sibling = s0.sibling
                s0.sibling = None

            if dbug:
                pTree.setPreAndPostOrder()
                pTree.draw()

            # Since the tree is rooted on y0 at the moment, we prune
            # off the x-subtree.  Usually this will be x0, but it can
            # happen that due to a rearrangement above that x1 is now
            # between y0 and x0.

            theX = None
            if x0.parent == y0:
                theX = x0
            elif x1.parent == y0:
                theX = x1
                #assert x0.parent == x1  Nope, not always.
            else:
                gm.append("Fix me.")
                raise Glitch, gm
            pTree.pruneSubTreeWithoutParent(theX) # this method returns x0 as well.
            pTree.reconnectSubTreeWithoutParent(theX, y1)

            if dbug:
                pTree.setPreAndPostOrder()
                pTree.draw()


    if oldRoot != pTree.root:
        pTree.reRoot(oldRoot, moveInternalName=False)

    if dbug:
        pTree.setPreAndPostOrder()
        pTree.draw()
        if y0Uncon and s0 != y1:
            print "The drawing above shows that Y extended"
        else:
            print "The drawing above shows that Y did not extend."

    if dbug:
        for n in pTree.nodes:
            if n.isLeaf:
                while '_' in n.name:
                    n.name = n.name[:-3]
                #if n.name.endswith('_xUp'):
                #    n.name = n.name[:-4]
                #elif n.name.endswith('_xDn'):
                #    n.name = n.name[:-4]
                #else:
                #    n.name = n.name[:-3]
            elif not n.isLeaf:
                n.name = None
            if n.br:
                n.br.textDrawSymbol = '-'
        #if xUp:
        #    xUp.br.textDrawSymbol = '-'
        #pTree.draw()

    xRearranged = False
    if x0Uncon:
        if r0 != x1:
            xRearranged = True
    yRearranged = False
    if y0Uncon:
        if s0 != y1:
            yRearranged = True
    if xRearranged or yRearranged:
        theProposal.topologyChanged = True

    if dbug:
        print "-" * 20
        print "xRearranged = %s" % xRearranged
        print "yRearranged = %s" % yRearranged
        if eA:
            for n in pTree.iterNodesNoRoot():
                if n.br == eA:
                    print "eA is from node", n.nodeNum
        else:
            print "eA is None"
        if eX:
            for n in pTree.iterNodesNoRoot():
                if n.br == eX:
                    print "eX is from node", n.nodeNum
        else:
            print "eX is None"
        if eY:
            for n in pTree.iterNodesNoRoot():
                if n.br == eY:
                    print "eY is from node", n.nodeNum
        else:
            print "eY is None"
        print "x0Uncon is", x0Uncon
        print "y0Uncon is", y0Uncon

    # n.flag is set if the condLikes need recalculating.  Edges eA,
    # eX, and eY will have their bigPDecks recalculated, and all the
    # nodes below those will have their flags set (in
    # Chain.proposeSp()).  However, there are additional nodes that
    # need the flag set.
    if xRearranged:
        p = x1
        while 1:
            p.flag = 1
            p = p.parent
            if not p:
                break
    if yRearranged:
        p = y0
        while 1:
            p.flag = 1
            p = p.parent
            if not p:
                break

    # if r1 is above r0, unusually, then eX points 'up' towards r1.
    # In that case, the branch on node x0 needs to have its lenChanged
    # set.
    if xRearranged:
        if r1.parent == x0:
            if x0.parent == r0:
                x0.br.lenChanged = True # by virtue of it being upside down.

    # if s0 is above y1, then there was a rearrangement due to the
    # re-rooting, so we need to set node.br.lenChanged from y0 down to s1.

    if yRearranged:
        if y1.isAncestorOf(s0):
            if s1.isAncestorOf(y0):
                p = y0
                while 1:
                    p.br.lenChanged = True
                    p = p.parent
                    if p == s1:
                        break
    

    
    # Do the branch length changes.  Taken nearly verbatim from Crux.
    # Thanks Jason!  The 3 edges to modify are eA, eX, and eY, which
    # were named above, as NodeBranch objects.
    
    lnProp = 0.0

    # Generate branch length multipliers and set new branch lengths.
    # eA.
    u = random.random()
    lnMA = etbrLambda * (u - 0.5)
    lnProp += lnMA
    mA = math.exp(lnMA)
    vA0 = eA.len
    vA1 = vA0 * mA
    eA.len = vA1
    eA.lenChanged = True
    # eX.
    if x0Uncon:
        u = random.random()
        lnMX = etbrLambda * (u - 0.5)
        lnProp += lnMX
        mX = math.exp(lnMX)
        vX0 = eX.len
        vX1 = vX0 * mX
        eX.len = vX1
        eX.lenChanged = True
    # eY.
    if y0Uncon:
        u = random.random()
        lnMY = etbrLambda * (u - 0.5)
        lnProp += lnMY
        mY = math.exp(lnMY)
        vY0 = eY.len
        vY1 = vY0 * mY
        eY.len = vY1
        eY.lenChanged = True
   
    # More from Crux ...
    # The prior ratio is the product of the prior ratios for each modified
    # branch length.  The number of internal branches does not change
    # (though the number of polytomies may change), so the topology prior
    # ratio is always 1.
    lnPrior = -self.mcmc.tunings.brLenPriorLambda * (vA1-vA0)
    if x0Uncon:
        lnPrior += -self.mcmc.tunings.brLenPriorLambda * (vX1-vX0)
    if y0Uncon:
        lnPrior += -self.mcmc.tunings.brLenPriorLambda * (vY1-vY0)
        
    # The proposal ratio is the product of the proposal ratios for
    # extension of each end of eA, as well as the branch multipliers.  The
    # ratio is 1 for the constrained/constrained and
    # unconstrained/unconstrained extension cases.
    #
    # nY0/nR0.
    if x0Uncon and r0 is not x1:
        if y0Uncon:
            if not r0Uncon:
                lnProp += math.log(1.0 - etbrPExt)
        elif r0Uncon:
            lnProp += math.log(1.0 / (1.0 - etbrPExt))
    # nX0/nS0.
    if y0Uncon and s0 is not y1:
        if x0Uncon:
            if not s0Uncon:
                lnProp += math.log(1.0 - etbrPExt)
        elif s0Uncon:
            lnProp += log(1.0 / (1.0 - etbrPExt))

    self.propTree.preAndPostOrderAreValid = False

    self.logPriorRatio = lnPrior
    self.logProposalRatio = lnProp
    
    # Check if we have a new combo of comp and rMatrix.  This might be
    # more efficient if I cleverly only look at the affected nodes.
    if 1:
        if theProposal.topologyChanged:
            for n in pTree.iterNodesNoRoot():
                for pNum in range(pTree.model.nParts):
                    theCompNum = n.parts[pNum].compNum
                    theRMatrixNum = n.br.parts[pNum].rMatrixNum
                    if pTree.model.parts[pNum].bQETneedsReset[theCompNum][theRMatrixNum]:
                        pf.p4_resetBQET(pTree.model.cModel, pNum, theCompNum, theRMatrixNum)
    





def proposePolytomy(self, theProposal):
    theProposal.doAbort = False
    dbug = False
    if dbug:
        #print "proposePolytomy() starting with this tree ..."
        #self.propTree.draw(width=80, addToBrLen=0.2)
        print "j There are %i internal nodes." % self.propTree.nInternalNodes
        if self.propTree.nInternalNodes == 1:
            print "-> so its a star tree -> proposeDeleteEdge is not possible."
        elif self.propTree.nInternalNodes == self.propTree.nTax - 2:
            print "-> so its a fully-resolved tree, so proposeAddEdge is not possible."

    if self.propTree.nInternalNodes == 1: # a star tree
        self.proposeAddEdge(theProposal)
    elif self.propTree.nInternalNodes == self.propTree.nTax - 2:
        candidateNodes = self._getCandidateNodesForDeleteEdge()
        if candidateNodes:
            self.proposeDeleteEdge(theProposal, candidateNodes)
        else:
            #gm = ["proposePolytomy()"]
            #gm.append("The tree is fully resolved, so I can't proposeAddEdge()")
            #gm.append("But there are no suitable nodes to remove.")
            #raise Glitch, gm
            theProposal.doAbort = True
            self.curTree._nInternalNodes = self.propTree._nInternalNodes
            return
    else:
        r = random.random()
        #r = 0.4
        if r < 0.5:
            self.proposeAddEdge(theProposal)
        else:
            candidateNodes = self._getCandidateNodesForDeleteEdge()
            if candidateNodes:
                self.proposeDeleteEdge(theProposal, candidateNodes)
            else:
                self.proposeAddEdge(theProposal)
    #if self.mcmc.constraints:
    #    print "checkSplitKeys() at the end of polytomy"
    #    self.propTree.checkSplitKeys()
    

def proposeAddEdge(self, theProposal):
    gm = ["Chain.proposeAddEdge()"]
    #print "proposeAddEdge() here"
    dbug = False
    pTree = self.propTree
    if 0:
        print "proposeAddEdge(), starting with this tree ..."
        pTree.draw()
        print "k There are %i internal nodes." % pTree.nInternalNodes
        print "root is node %i" % pTree.root.nodeNum
    allPolytomies = []
    for n in pTree.iterInternalsNoRoot():
        if n.getNChildren() > 2:
            allPolytomies.append(n)
    if pTree.root.getNChildren() > 3:
        allPolytomies.append(pTree.root)
    
    theChosenPolytomy = random.choice(allPolytomies)

    # We want to choose one of the possible ways to add a node.  See
    # Lewis et al page 246, left top.  "The number of distinct ways of
    # dividing k edges into two groups, making sure that at least 3
    # edges are attached to each node afterwards, is 2^{k-1} - k - 1".
    # For non-root polytomies (with 3 or more children), it is
    # straightforward, but for root polytomies (ie with 4 or more
    # children) it is different.  I think in the case of root
    # polytomies that they will be equivalent to non-root polytomies
    # if I arbitrarily consider one randomly chosen child node to
    # take the role that the parent takes in the non-root-polytomies.
    # So a 4-child root will be considered to have a parent-like node
    # and 3 children.
    if theChosenPolytomy != pTree.root:
        nChildren = theChosenPolytomy.getNChildren()
        k = nChildren + 1
        childrenNodeNums = pTree.getChildrenNums(theChosenPolytomy) 
    else:
        # Its the root.  So we say that a random child takes the role
        # of the "parent", for purposes of these calculations.
        nChildren = theChosenPolytomy.getNChildren() - 1  # n - 1 children
        k = nChildren + 1
        childrenNodeNums = pTree.getChildrenNums(theChosenPolytomy) # Yes, all children.
        
    nPossibleWays = math.pow(2, k-1) - k - 1
    if dbug:
        print "These nodes are polytomies: %s" % [n.nodeNum for n in allPolytomies]
        print "We randomly choose to do node %i" % theChosenPolytomy.nodeNum
        print "It has %i children, so k=%i, so there are %i possible ways to add a node." % (
            nChildren, k, nPossibleWays)

    # We want to choose one of the possible ways to add a node, but we
    # want to choose it randomly.  I'll describe it for the case with
    # nChildren=5, so k is 6.  We know already that there are
    # nPossibleWays=25 different ways to add a node.  The complication
    # is that we could make a new group of 2, 3, or 4 nInNewGroup, and it will be
    # different numbers of possible ways in each.  The numbers of each are given by
    # func.nChoosek(), so there are 10 ways to make a group of 2 from 5
    # children, 10 ways to make a group of 3 from 5 children, and 5
    # ways to make a group of 4 from 5 children.  So thats [10, 10,
    # 5], which sums to 25 (nPossibleWays).  So we can make a
    # cumulative sum list ie [10, 20, 25], and use it to choose one
    # group randomly.
    nChooseKs = []
    for i in range(2, nChildren):
        nChooseKs.append(func.nChooseK(nChildren, i))
    cumSum = [nChooseKs[0]]
    for i in range(len(nChooseKs))[1:]:
        cumSum.append(nChooseKs[i] + cumSum[i-1])
    ran = random.randrange(nPossibleWays)
    for i in range(len(cumSum)):
        if ran < cumSum[i]:
            break
    nInNewGroup = i + 2
    # Ok, so we have decided that of the nChildren of
    # theChosenPolytomy, we will make a new node with a group of
    # nInNewGroup of them.  For that, we can use random.sample().
    newChildrenNodeNums = random.sample(childrenNodeNums, nInNewGroup)


    if dbug:
        print "The nChooseKs are %s" % nChooseKs
        print "The cumSum is %s" % cumSum
        print "Since there are nPossibleWays=%i, we choose a random number from 0-%i" % (
            nPossibleWays, nPossibleWays-1)
        print "->We chose a random number: %i" % ran
        print "So we choose the group at index %i, which means nInNewGroup=%i" % (i, nInNewGroup)
        print "So we make a new node with newChildrenNodeNums %s" % newChildrenNodeNums
        #sys.exit()

    # Choose to add a node between theChosenPolytomy and the first in
    # the list of newChildrenNodeNums.  The node that we add will be
    # chosen from pTree.nodes for the first node where both the parent
    # and the leftChild are None.
    firstNode = pTree.nodes[newChildrenNodeNums[0]]
    for newNode in pTree.nodes:
        if not newNode.parent and not newNode.leftChild:
            break
    #print "Got newNode = %i" % newNode.nodeNum

    # Add the newNode between theChosenPolytomy and firstNode
    newNode.parent = theChosenPolytomy
    newNode.leftChild = firstNode
    firstNode.parent = newNode
    if theChosenPolytomy.leftChild == firstNode:
        theChosenPolytomy.leftChild = newNode
    else:
        oldCh = theChosenPolytomy.leftChild
        while oldCh.sibling != firstNode:
            oldCh = oldCh.sibling
        oldCh.sibling = newNode
    if firstNode.sibling:
        newNode.sibling = firstNode.sibling
        firstNode.sibling = None
    pTree.setPreAndPostOrder()
    pTree._nInternalNodes += 1

    if 0:
        #pTree.setPreAndPostOrder()
        pTree.draw()
    
    for nodeNum in newChildrenNodeNums[1:]:
        n = pTree.pruneSubTreeWithoutParent(nodeNum)
        pTree.reconnectSubTreeWithoutParent(n, newNode)

    # Choose a branch length for newNode.  See LewisHolderHolsinger eqn 6.
    newNode.br.len = - (1.0/self.mcmc.tunings.brLenPriorLambda) * math.log(1. - random.random())
    if newNode.br.len < var.BRLEN_MIN or newNode.br.len > var.BRLEN_MAX:
        safety = 0
        while newNode.br.len < var.BRLEN_MIN or newNode.br.len > var.BRLEN_MAX:
            newNode.br.len = - (1.0/self.mcmc.tunings.brLenPriorLambda) * math.log(1. - random.random())
            safety += 1
            if safety > 20:
                gm.append("Unable to find a good branch length for the new edge.")
                gm.append("Probably a programming error.")
                raise Glitch, gm
    if var.doMcmcSp:
        newNode.br.lenChanged = True

    # Calculate the rawSplitKey and splitKey.
    if self.mcmc.constraints:
        children = [n for n in newNode.iterChildren()]
        x = children[0].br.rawSplitKey
        for n in children[1:]:
            y = n.br.rawSplitKey
            x = x | y  # '|' is bitwise "OR".
        newNode.br.rawSplitKey = x
        if 1 & newNode.br.rawSplitKey: # Ie "Does rawSplitKey contain a 1?" or "Is rawSplitKey odd?"
            if self.mcmc.constraints:
                newNode.br.splitKey = self.mcmc.constraints.allOnes ^ newNode.br.rawSplitKey # "^" is xor, a bit-flipper.
            else:
                allOnes = 2L**(self.propTree.nTax) - 1
                newNode.br.splitKey = allOnes ^ newNode.br.rawSplitKey
        else:
            newNode.br.splitKey = newNode.br.rawSplitKey
    
    # Its a newly-added node, possibly in a new context.  We need to
    # deal with model stuff if it isHet.  The model.isHet if any part
    # isHet.
    if pTree.model.isHet:
        for pNum in range(pTree.model.nParts):
            mp = pTree.model.parts[pNum]
            if mp.isHet:
                # We want to make the model prams appropriate, and not
                # just use what was left over from the last time the
                # node was in the tree.
                if mp.nComps > 1:
                    # Pick a child, and use the same comp.
                    oneChildNum = random.choice(newChildrenNodeNums)
                    newNode.parts[pNum].compNum = pTree.nodes[oneChildNum].parts[pNum].compNum
                    mp.comps[newNode.parts[pNum].compNum].nNodes += 1 
                if mp.nRMatrices > 1:
                    oneChildNum = random.choice(newChildrenNodeNums)
                    newNode.br.parts[pNum].rMatrixNum = pTree.nodes[oneChildNum].br.parts[pNum].rMatrixNum
                    mp.rMatrices[newNode.br.parts[pNum].rMatrixNum].nNodes += 1 
                if mp.nGdasrvs > 1:
                    oneChildNum = random.choice(newChildrenNodeNums)
                    newNode.br.parts[pNum].gdasrvNum = pTree.nodes[oneChildNum].br.parts[pNum].gdasrvNum
                    mp.gdasrvs[newNode.br.parts[pNum].gdasrvNum].nNodes += 1
                
    if dbug:
        pTree.setPreAndPostOrder()
        pTree.draw()
    
    # Now the Hastings ratio.  First calculate gamma_B.  If the
    # current tree is a star tree (nInternalNodes == 1) and the
    # proposed tree is not fully resolved (ie is less than
    # len(self.propTree.nodes) - 2), then gamma_B is 0.5.
    if (self.curTree.nInternalNodes == 1) and (pTree.nInternalNodes < (len(pTree.nodes) - 2)):
        gamma_B = 0.5
    # If the proposed tree is fully resolved and the current tree is not the star tree
    elif (pTree.nInternalNodes == (len(pTree.nodes) - 2)) and (self.curTree.nInternalNodes > 1):
        gamma_B = 2.0
    else:
        gamma_B = 1.0

    # n_e is number of internal edges present before the Add-edge move.  That would be self.curTree.nInternalNodes - 1
    n_e = float(self.curTree.nInternalNodes - 1)
    # n_p is the number of polytomies present before the move, len(allPolytomies)
    n_p = float(len(allPolytomies))
    hastingsRatio = (gamma_B * n_p * float(nPossibleWays)) / (1.0 + n_e)

    if dbug:
        print "The new node is given a random branch length of %f" % newNode.br.len
        print "For the Hastings ratio ..."
        print "gamma_B is %.1f" % gamma_B
        print "n_e is %.0f" % n_e
        print "k is (still) %i, and (2^{k-1} - k - 1) = nPossibleWays is still %i" % (k, nPossibleWays)
        print "n_p = %.0f is the number of polytomies present before the move." % n_p
        print "So the hastings ratio is %f" % hastingsRatio

    self.logProposalRatio = math.log(hastingsRatio)


    if 0:
        priorRatio = self.mcmc.tunings.brLenPriorLambda * math.exp(- self.mcmc.tunings.brLenPriorLambda * newNode.br.len)
        if dbug:
            print "The self.mcmc.tunings.brLenPriorLambda is %f" % self.mcmc.tunings.brLenPriorLambda
            print "So the prior ratio is %f" % priorRatio

        self.logPriorRatio = math.log(priorRatio)

        # The Jacobian
        jacobian = 1.0 / (self.mcmc.tunings.brLenPriorLambda * math.exp(- self.mcmc.tunings.brLenPriorLambda * newNode.br.len))
        self.logJacobian = math.log(jacobian)
        print "logPriorRatio = %f, logJacobian = %f" % (self.logPriorRatio, self.logJacobian)

    # Here I pull a fast one, as explained in Lewis et al.  The
    # priorRatio and the Jacobian terms cancel out.  So the logs might
    # as well be zeros.
    self.logPriorRatio = 0.0
    #self.logJacobian = 0.0
    # That was easy, wasn't it?
    if self.mcmc.tunings.doPolytomyResolutionClassPrior:
        # We are gaining a node.  So the prior ratio is T_{n,m + 1} /
        # (T_{n,m} * C) .  We have the logs, and the result is the
        # log.
        if 0:
            print "-" * 30
            print 'curTree.nInternalNodes', self.curTree.nInternalNodes
            print 'pTree.nInternalNodes', pTree.nInternalNodes
            print 'logBigT[curTree.nInternalNodes]', theProposal.logBigT[self.curTree.nInternalNodes]
            #print math.exp(theProposal.logBigT[self.curTree.nInternalNodes])
            print 'C ', self.mcmc.tunings.polytomyPriorLogBigC
            print 'logBigT[pTree.nInternalNodes]', theProposal.logBigT[pTree.nInternalNodes]
            #print math.exp(theProposal.logBigT[pTree.nInternalNodes])
            print "-" * 30
        self.logPriorRatio = (theProposal.logBigT[self.curTree.nInternalNodes] -
                              (self.mcmc.tunings.polytomyPriorLogBigC +
                              theProposal.logBigT[pTree.nInternalNodes]))
        
    else:
        if self.mcmc.tunings.polytomyPriorLogBigC:
            self.logPriorRatio =  -self.mcmc.tunings.polytomyPriorLogBigC
        else:
            self.logPriorRatio = 0.0
    #print "\ngaining a node, m %2i->%2i. logPriorRatio is %f" % (self.curTree.nInternalNodes,
    #                                                              pTree.nInternalNodes, self.logPriorRatio)
    
    
def _getCandidateNodesForDeleteEdge(self):
    pTree = self.propTree
    nodesWithInternalEdges = [n for n in pTree.iterInternalsNoRoot()]

    # Remove any that might violate constraints.
    if self.mcmc.constraints:
        nodesToRemove = []
        for n in nodesWithInternalEdges:
            if n.br.splitKey in self.mcmc.constraints.constraints:
                nodesToRemove.append(n)
        for n in nodesToRemove:
            nodesWithInternalEdges.remove(n)
        

    # We need to check that we will not be deleting modelThings that are only on one node.
    if pTree.model.isHet:
        nodesToRemove = []
        for pNum in range(pTree.model.nParts):
            if pTree.model.parts[pNum].isHet:
                mp = pTree.model.parts[pNum]
                if mp.nComps > 1:
                    for mtNum in range(mp.nComps):
                        mt = mp.comps[mtNum]
                        if mt.nNodes <= 1: # These modelThings are on only one node
                            for nNum in range(len(nodesWithInternalEdges)):
                                n = nodesWithInternalEdges[nNum]
                                if n.parts[pNum].compNum == mtNum:
                                    if n not in nodesToRemove:
                                        nodesToRemove.append(n)
                if mp.nRMatrices > 1:
                    for mtNum in range(mp.nRMatrices):
                        mt = mp.rMatrices[mtNum]
                        if mt.nNodes <= 1: # These modelThings are on only one node
                            for nNum in range(len(nodesWithInternalEdges)):
                                n = nodesWithInternalEdges[nNum]
                                if n.br.parts[pNum].rMatrixNum == mtNum:
                                    if n not in nodesToRemove:
                                        nodesToRemove.append(n)
                if mp.nGdasrvs > 1:
                    for mtNum in range(mp.nGdasrvs):
                        mt = mp.gdasrvs[mtNum]
                        if mt.nNodes <= 1: # These modelThings are on only one node
                            for nNum in range(len(nodesWithInternalEdges)):
                                n = nodesWithInternalEdges[nNum]
                                if n.br.parts[pNum].gdasrvNum == mtNum:
                                    if n not in nodesToRemove:
                                        nodesToRemove.append(n)
        #print "There are %i nodesWithInternalEdges, and I need to remove %i nodes" % (
        #    len(nodesWithInternalEdges) ,len(nodesToRemove))
        for n in nodesToRemove:
            nodesWithInternalEdges.remove(n)
    return nodesWithInternalEdges
    
def proposeDeleteEdge(self, theProposal, candidateNodes):
    
    dbug = False
    pTree = self.propTree
    #print "doing proposeDeleteEdge()"
    if 0:
        print "proposeDeleteEdge(), starting with this tree ..."
        pTree.draw()
        print "m There are %i internal nodes (before deleting the edge)." % pTree.nInternalNodes
                        
    if not candidateNodes:
        raise Glitch, "proposeDeleteEdge() could not find a good node to attempt to delete."
    
    theChosenNode = random.choice(candidateNodes)
    if dbug:
        print "There are %i candidateNodes." % len(candidateNodes)
        print "node nums %s" % [n.nodeNum for n in candidateNodes]
        print "Randomly choose node %s" % theChosenNode.nodeNum

    
    if pTree.model.isHet:
        for pNum in range(pTree.model.nParts):
            if pTree.model.parts[pNum].isHet:
                mp = pTree.model.parts[pNum]
                if mp.nComps > 1:
                    mp.comps[theChosenNode.parts[pNum].compNum].nNodes -= 1 
                if mp.nRMatrices > 1:
                    mp.rMatrices[theChosenNode.br.parts[pNum].rMatrixNum].nNodes -= 1 
                if mp.nGdasrvs > 1:
                    mp.gdasrvs[theChosenNode.br.parts[pNum].gdasrvNum].nNodes -= 1 

    theNewParent = theChosenNode.parent
    theRightmostChild = theChosenNode.rightmostChild()
    theLeftSib = theChosenNode.leftSibling()
    if theLeftSib:
        theLeftSib.sibling = theChosenNode.leftChild
    else:
        theNewParent.leftChild = theChosenNode.leftChild
    for n in theChosenNode.iterChildren():
        n.parent = theNewParent
    theRightmostChild.sibling = theChosenNode.sibling
    theChosenNode.wipe()
    pTree.setPreAndPostOrder()
    pTree._nInternalNodes -= 1
    #print pTree.preOrder
    #if dbug:
    #    pTree.draw()

    if var.doMcmcSp:
        p = theNewParent.leftChild
        while p != pTree.root:
            p = p.parent
            p.flag = 1
            #print "setting flag of node %i" % p.nodeNum
    
        

    # Hastings ratio.  First calculate the gamma_D.  If the current
    # tree is fully resolved and the proposed tree is not the star
    # tree, then gamma_D is 0.5
    if (self.curTree.nInternalNodes == len(pTree.nodes) - 2) and pTree.nInternalNodes != 1:
        gamma_D = 0.5
    # If the proposed tree is the star tree and the current tree is not fully resolved
    elif (self.curTree.nInternalNodes < len(pTree.nodes) - 2) and pTree.nInternalNodes == 1:
        gamma_D = 2.
    else:
        gamma_D = 1.

    # n_e is the number of internal edges in existence before the move, which would be nInternalNodes - 1
    n_e = float(self.curTree.nInternalNodes - 1)
    # nStar_p is the number of polytomies in the tree after the move.
    nStar_p = 0
    for n in pTree.iterInternalsNoRoot():
        if n.getNChildren() > 2:
            nStar_p += 1
    if pTree.root.getNChildren() > 3:
        nStar_p += 1
    nStar_p = float(nStar_p)
    # kStar is the number of edges emanating from the polytomy created (or enlarged) by the move.
    kStar = theNewParent.getNChildren()
    if theNewParent.parent:
        kStar += 1

    hastingsRatio = (gamma_D * n_e) / (nStar_p * (2**(kStar - 1) - kStar - 1))
    self.logProposalRatio = math.log(hastingsRatio)

    if 0:
        # Now the prior ratio.  The prior probability density f(nu) for a
        # branch length is lambda * exp(-lambda * nu).  To a first
        # approximation, with equal priors on topologies, the prior ratio
        # is 1/f(nu)
        priorRatio = 1.0/(self.mcmc.tunings.brLenPriorLambda * math.exp(- self.mcmc.tunings.brLenPriorLambda * theChosenNode.br.len))
        if dbug:
            print "The self.mcmc.tunings.brLenPriorLambda is %f" % self.mcmc.tunings.brLenPriorLambda
            print "So the prior ratio is %f" % priorRatio

        self.logPriorRatio = math.log(priorRatio)    

        # The Jacobian
        jacobian = self.mcmc.tunings.brLenPriorLambda * math.exp(- self.mcmc.tunings.brLenPriorLambda * theChosenNode.br.len)
        self.logJacobian = math.log(jacobian)
        print "logPriorRatio = %f, logJacobian = %f" % (self.logPriorRatio, self.logJacobian)

    # Here I pull a fast one, as explained in Lewis et al.  The
    # priorRatio and the Jacobian terms cancel out.  So the logs might
    # as well be zeros.
    self.logPriorRatio = 0.0
    #self.logJacobian = 0.0
    # That was easy, wasn't it?

    if self.mcmc.tunings.doPolytomyResolutionClassPrior:
        # We are losing a node.  So the prior ratio is (T_{n,m} * C) /
        # T_{n,m - 1}.  We have the logs, and the result is the log.
        if 0:
            print "-" * 30
            print 'curTree.nInternalNodes', self.curTree.nInternalNodes
            print 'pTree.nInternalNodes', pTree.nInternalNodes
            print 'logBigT[curTree.nInternalNodes]', theProposal.logBigT[self.curTree.nInternalNodes]
            #print math.exp(theProposal.logBigT[self.curTree.nInternalNodes])
            print 'C ', self.mcmc.tunings.polytomyPriorLogBigC
            print 'logBigT[pTree.nInternalNodes]', theProposal.logBigT[pTree.nInternalNodes]
            #print math.exp(theProposal.logBigT[pTree.nInternalNodes])
            print "-" * 30
        self.logPriorRatio = ((theProposal.logBigT[self.curTree.nInternalNodes] +
                               self.mcmc.tunings.polytomyPriorLogBigC) -
                              theProposal.logBigT[pTree.nInternalNodes])
        
    else:
        if self.mcmc.tunings.polytomyPriorLogBigC:
            self.logPriorRatio =  self.mcmc.tunings.polytomyPriorLogBigC
        else:
            self.logPriorRatio = 0.0

    #print " losing a node, m %2i->%2i. logPriorRatio is %f" % (self.curTree.nInternalNodes,
    #                                                           pTree.nInternalNodes, self.logPriorRatio)

