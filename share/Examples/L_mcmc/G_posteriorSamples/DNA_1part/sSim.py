nTax = 11
taxNames = list(string.uppercase[:nTax])
a = func.newEmptyAlignment(dataType='dna', taxNames=taxNames, length=200)
d = Data([a])

t = func.randomTree(taxNames=taxNames)

t.data = d
t.newComp(free=0, spec='specified', val=[0.4, 0.3, 0.0])
t.newRMatrix(free=0, spec='specified', val=[2., 3., 4., 5., 6., 7.])
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=0, val=0.5)
t.setPInvar(free=0, val=0.0)

func.reseedCRandomizer(os.getpid())
t.simulate()

d.writeNexus('d.nex', writeDataBlock=True)

