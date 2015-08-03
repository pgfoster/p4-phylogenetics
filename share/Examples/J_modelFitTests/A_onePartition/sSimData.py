read('t.nex')
t = var.trees[0]
var.alignments.append(func.newEmptyAlignment(dataType='dna', symbols=None, taxNames=t.taxNames, length=1000))
d = Data()
t.data = d

c1 = t.newComp(free=1, spec='specified', val=[0.1, 0.2, 0.3])
c2 = t.newComp(free=1, spec='specified', val=[0.4, 0.3, 0.1])
t.setModelThing(c1, node=0, clade=1)
t.setModelThing(c2, node=1, clade=1)
t.newRMatrix(free=0, spec='specified', val=[1.2, 6.5, 1.3, 9.8, 1.1, 4.3])
t.setPInvar(free=0, val=0.0)
t.setNGammaCat(nGammaCat=1)

t.simulate()
t.data.writeNexus('d.nex')
t.tPickle('noOpt')
