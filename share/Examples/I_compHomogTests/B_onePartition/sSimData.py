read('t.nex')
t = var.trees[0]
var.alignments.append(func.newEmptyAlignment(dataType='dna', symbols=None, taxNames=t.taxNames, length=300))
d = Data()
t.data = d

c1 = t.newComp(free=1, spec='specified', val=[0.0, 0.2, 0.3])
t.newRMatrix(free=1, spec='specified', val=[1.2, 6.5, 1.3, 9.8, 1.1, 1.0])
t.setPInvar(free=1, val=0.1)
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)

t.simulate()
t.data.writeNexus('d.nex')
t.tPickle('noOpt')
