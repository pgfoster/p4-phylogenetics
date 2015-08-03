read("../d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=1)
t.setPInvar(free=0, val=0.0)
m = Mcmc(t, nChains=1, runNum=0, sampleInterval=10, checkPointInterval=2000)

# Set up a strong prior on polytomies
if 1:
    m.prob.polytomy = 1.0
    m.prob.brLen = 0.001
    m.tunings.doPolytomyResolutionClassPrior = True  # False by default
    import math
    m.tunings.polytomyPriorLogBigC = math.log(10.)   # Default is zero

m.run(4000)



