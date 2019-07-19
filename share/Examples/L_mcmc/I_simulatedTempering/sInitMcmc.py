read("d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=1)
#t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)

cpInterv = 100000

# We use only one chain.  The sampleInterval is ignored --- Whenever
# the chain has a temperature of zero, it is sampled.  Potentially a
# lot of samples!  Setting simTemp and simTempMax turns on simulated
# tempering.  simTemp is the number of temperatures, and needs to be 2
# or more.

m = Mcmc(t, nChains=1, runNum=0, checkPointInterval=cpInterv, simTemp=10, simTempMax=10.0)
m.run(1, writeSamples=False)
m.gen = -1
m.checkPoint()
