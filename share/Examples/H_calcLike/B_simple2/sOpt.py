var.verboseRead = 0
read('t.nex')
t = var.trees[0]
read('d.nex')

t.data = Data()
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setPInvar(free=1, val=0.0)
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=1.0)
t.optLogLike()
# If I wanted to save the tree (with branch lengths, of course), I could
#t.writeNexus('optTree.nex')
# If I wanted to see the optimized model params, I could
#t.model.dump()



