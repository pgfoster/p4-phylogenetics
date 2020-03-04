var.verboseRead = 0
var.warnReadNoFile = 0
read('data.nex')
d = Data()
d.compoSummary()

read('(A,B,C,D,E);')
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
