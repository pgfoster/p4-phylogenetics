var.warnReadNoFile = 0
var.verboseRead = 0
func.reseedCRandomizer(os.getpid())

nTax = 5
taxNames = list(string.uppercase[:nTax])
a = func.newEmptyAlignment(dataType='dna', taxNames=taxNames, length=400)
d = Data([a])

t = func.randomTree(taxNames=taxNames)
t.data = d
c1 = t.newComp(free=1, spec='specified', val=[0.1, 0.2, 0.3])
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=1)
t.setPInvar(free=0, val=0.0)

t.simulate()
d.writeNexus('d.nex')
