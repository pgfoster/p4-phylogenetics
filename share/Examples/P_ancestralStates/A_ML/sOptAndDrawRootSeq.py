var.verboseRead = 0
var.warnReadNoFile = 0
read('d.nex')
d = Data()
d.compoSummary()

read('simTree.nex')
t = var.trees[0]
t.data = d
c1 = t.newComp(free=1, spec='empirical')
c2 = t.newComp(free=1, spec='empirical')

# Put the c1 comp on all the nodes of the tree.  Then put c2 on the
# root, over-riding c1 that is already there.
t.setModelThing(c1, node=0, clade=1)
t.setModelThing(c2, node=0, clade=0)

t.newRMatrix(free=0, spec='ones')
t.setNGammaCat(nGammaCat=1)
t.setPInvar(free=0, val=0.0)

t.optLogLike()

print('\nAfter optimizing, the composition of the model for the non-root nodes is:') 
print(t.model.parts[0].comps[0].val)
print('...and the composition of the root model is:')
print(t.model.parts[0].comps[1].val)
t.write()

counts = [0] * 4
for rep in range(100):
    ancSt = t.ancestralStateDraw()
    for i in range(4):
        ch = 'acgt'[i]
        cnt = ancSt.count(ch)
        counts[i] += cnt
mySum = float(sum(counts))
print("\noptimized      draws")
for i in range(4):
    print("  %.5f     %.4f" % (t.model.parts[0].comps[1].val[i], counts[i]/mySum))
