theRunNum = int(var.argvAfterDoubleDash[0])

read("../d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=1)
#t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)
m = Mcmc(t, nChains=4, runNum=theRunNum, sampleInterval=10, checkPointInterval=2000, simulate=None)

m.run(2000)

if 1 and func.which("gnuplot"):
    h = Numbers('mcmc_likes_%i' % theRunNum, col=1)
    h.plot()

