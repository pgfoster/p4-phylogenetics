# Calculate a likelihood.
read("d.nex")
d = Data()

read("t.nex")
t = var.trees[0]

t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=0, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=1, val=0.2)

t.optLogLike()
t.writeNexus('optTree.nex')
t.model.dump()

