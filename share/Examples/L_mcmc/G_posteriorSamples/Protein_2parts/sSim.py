nTax = 5
taxNames = list(string.uppercase[:nTax])
a = func.newEmptyAlignment(dataType='protein', taxNames=taxNames, length=200)
b = func.newEmptyAlignment(dataType='protein', taxNames=taxNames, length=133)
d = Data([a,b])

t = func.randomTree(taxNames=taxNames)

t.data = d
pNum = 0
t.newComp(partNum=pNum, free=0, spec='wag')
t.newRMatrix(partNum=pNum, free=0, spec='wag')
t.setNGammaCat(partNum=pNum, nGammaCat=4)
t.newGdasrv(partNum=pNum, free=0, val=0.5)
t.setPInvar(partNum=pNum, free=0, val=0.0)
t.setRelRate(partNum=pNum, val=0.5)

pNum = 1
t.newComp(partNum=pNum, free=0, spec='jtt')
t.newRMatrix(partNum=pNum, free=0, spec='jtt')
t.setNGammaCat(partNum=pNum, nGammaCat=4)
t.newGdasrv(partNum=pNum, free=0, val=0.5)
t.setPInvar(partNum=pNum, free=0, val=0.2)
t.setRelRate(partNum=pNum, val=1.5)

t.model.relRatesAreFree = 1
func.reseedCRandomizer(os.getpid())
t.simulate()

a.bluntEndLigate(b)
a.writeNexus('d.nex', writeDataBlock=True)

