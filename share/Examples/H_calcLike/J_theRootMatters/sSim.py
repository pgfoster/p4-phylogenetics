var.warnReadNoFile = 0
var.verboseRead = 0
func.reseedCRandomizer(os.getpid())

nTax = 4
taxNames = list(string.uppercase[:nTax])
a = func.newEmptyAlignment(dataType='dna', taxNames=taxNames, length=500)
d = Data([a])
read('((A:0.4, B:0.4), C:0.4, D:0.4);')
t = var.trees[0]
t.data = d
c1 = t.newComp(free=1, spec='specified', val=[0.25, 0.2, 0.3])
c2 = t.newComp(free=1, spec='specified', val=[0.2, 0.3, 0.25])
t.setModelThing(c1, node=t.root, clade=1) 
t.setModelThing(c2, node=1, clade=1)
t.newRMatrix(free=1, spec='specified', val=[1.3, 3.2, 1.7, 8.0, 2.1, 0.9])
t.setNGammaCat(nGammaCat=1)
t.setPInvar(free=0, val=0.0)
t.simulate()
#t.draw(model=True)
d.writeNexus('d.nex')
