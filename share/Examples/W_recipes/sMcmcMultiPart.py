# Do an MCMC with more than one data partition.
read("d.nex")
a = var.alignments[0]
a.setCharPartition('p1')
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d

pNum = 0
t.newComp(partNum=pNum, free=1, spec='empirical')
t.newRMatrix(partNum=pNum, free=1, spec='ones')
t.setNGammaCat(partNum=pNum, nGammaCat=4)
t.newGdasrv(partNum=pNum, free=1, val=0.5)
t.setPInvar(partNum=pNum, free=1, val=0.1)
t.setRelRate(partNum=pNum, val=5.0)

pNum = 1
t.newComp(partNum=pNum, free=1, spec='empirical')
t.newRMatrix(partNum=pNum, free=1, spec='ones')
t.setNGammaCat(partNum=pNum, nGammaCat=1)
#t.newGdasrv(partNum=pNum, free=0, val=0.5)
t.setPInvar(partNum=pNum, free=1, val=0.1)
t.setRelRate(partNum=pNum, val=2.0)

t.model.relRatesAreFree = True

m = Mcmc(t, nChains=4, runNum=0, sampleInterval=10, checkPointInterval=2000)
m.run(4000)
func.summarizeMcmcPrams(skip=200)


