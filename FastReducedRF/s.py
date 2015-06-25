import fastReducedRF
import numpy
import pyublas # not explicitly used--but makes converters available

nTax = 7
import string
tNames = list(string.uppercase[:nTax])

# Here is a supertree, the first of many.
var.warnReadNoFile = False
read("(((D, E), A), B, (C, (F, G)));")
bigT = var.trees.pop()
bigT.taxNames = tNames
bigT.draw()
assert len(bigT.postOrder)

# These are the input trees.
read("(D, A, (C, (G, B)));")
read("(E, G, (D, F));")

for t in var.trees:
    t.draw()

# Start it up.  bigTr is a fastReducedRF.Tr object, a vaguely Tree-like object
# for the bigT
frrf = fastReducedRF.Frrf(len(tNames))
bigTr = frrf.setBigT(len(bigT.nodes), bigT.nTax, bigT.postOrder)

# Set the topology if bigTr.
for n in bigT.nodes:
    if n.parent:
        bigTr.setParent(n.nodeNum, n.parent.nodeNum)
    if n.leftChild:
        bigTr.setLeftChild(n.nodeNum, n.leftChild.nodeNum)
    else:
        bigTr.setNodeTaxNum(n.nodeNum, tNames.index(n.name))
    if n.sibling:
        bigTr.setSibling(n.nodeNum, n.sibling.nodeNum)

# Make Tr ojects for the input trees, and set their topologies.
for t in var.trees:
    #t.beta = 2.0  # STMcmc would do this
    tr = frrf.appendInTree(len(t.nodes), t.nTax, t.postOrder)
    for n in t.nodes:
        if n.parent:
            tr.setParent(n.nodeNum, n.parent.nodeNum)
        if n.leftChild:
            tr.setLeftChild(n.nodeNum, n.leftChild.nodeNum)
        else:
            tr.setNodeTaxNum(n.nodeNum, tNames.index(n.name))
        if n.sibling:
            tr.setSibling(n.nodeNum, n.sibling.nodeNum)

#tr.dump()
#frrf.dump()
#print "=" * 25
frrf.setInTreeTaxBits()
#frrf.dump()
frrf.setInTreeInternalBits()
frrf.maybeFlipInTreeBits()
print
frrf.setBigTInternalBits()
print "-" * 50
#frrf.dump()
sd = frrf.getLogLike(1.)
print sd
rfDist = bigT.inputTreesToSuperTreeDistances(var.trees, doSd=True, doScqdist=False)
print rfDist
print "====="
#frrf.dump()

if 0:
    for i in range(2000):
        bigT.randomSpr()
        bigT.nni()
        #bigT.write()
        if 1:
            #bigT.draw()
            frrf.wipeBigTPointers()
            for n in bigT.nodes:
                if n.parent:
                    bigTr.setParent(n.nodeNum, n.parent.nodeNum)
                if n.leftChild:
                    bigTr.setLeftChild(n.nodeNum, n.leftChild.nodeNum)
                #else:
                #    bigTr.setNodeTaxNum(n.nodeNum, tNames.index(n.name))
                if n.sibling:
                    bigTr.setSibling(n.nodeNum, n.sibling.nodeNum)
            frrf.setBigTInternalBits()
            sd = frrf.getSymmDiff()
            #print sd
            if 1:
                rfDist = bigT.inputTreesToSuperTreeDistances(var.trees, doSd=True, doScqdist=False)
                if sd == rfDist:
                    pass
                    #print "ok"
                else:
                    print "====================== differs ", sd, rfDist
                
        
    



