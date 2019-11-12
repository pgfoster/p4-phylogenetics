nTax = 7
taxNames = list(string.ascii_uppercase[:nTax])
a = func.newEmptyAlignment(dataType='dna', taxNames=taxNames, length=511)
d = Data([a])

t = func.randomTree(taxNames=taxNames)

t.data = d
t.newComp(free=0, spec='specified', val=[0.1, 0.2, 0.3])
t.newRMatrix(free=0, spec='specified', val=[2., 3., 4., 5., 6., 7.])
t.setNGammaCat(nGammaCat=1)
#t.newGdasrv(free=0, val=0.5)
t.setPInvar(free=0, val=0.0)

func.reseedCRandomizer(os.getpid())
t.simulate()

d.writeNexus('d.nex', writeDataBlock=True)
t.writeNexus("simTree.nex")

