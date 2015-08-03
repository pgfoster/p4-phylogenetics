read("../../K_thermus/noTRuberNoGapsNoAmbiguities.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)
m = Mcmc(t, nChains=1, runNum=0, sampleInterval=100, checkPointInterval=None, simulate=2)
m.run(100000)
func.summarizeMcmcPrams(skip=500)


