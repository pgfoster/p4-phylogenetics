var.warnReadNoFile = 0
var.verboseRead = 0
func.reseedCRandomizer(os.getpid())

nTax = 5
taxNames = list(string.uppercase[:nTax])
a = func.newEmptyAlignment(dataType='dna', taxNames=taxNames, length=100)
d = Data([a])

read('(B:0.5, ((D:0.4, A:0.3), C:0.5), E:0.5);')
t = var.trees[0]
t.taxNames = taxNames

t.data = d
c1 = t.newComp(free=1, spec='specified', val=[0.1, 0.2, 0.3])
c2 = t.newComp(free=1, spec='specified', val=[0.2, 0.3, 0.4])
t.setModelThing(c1, node=t.root, clade=1)
t.setModelThing(c2, node=3, clade=1)
r1 = t.newRMatrix(free=1, spec='specified', val=[1.2, 3.1, 5.7, 2.2, 8.2, 0.78])
r2 = t.newRMatrix(free=1, spec='specified', val=[4.1, 1.3, 1.7, 12.2, 1.2, 3.9])
t.setModelThing(r1, node=2, clade=1)
t.setModelThing(r2, node=1, clade=1)
t.setModelThing(r2, node=7, clade=1)
t.setNGammaCat(nGammaCat=1)
t.setPInvar(free=0, val=0.0)
t.draw(model=True)

t.simulate()
t.tPickle('sim')
d.writeNexus('d.nex')
