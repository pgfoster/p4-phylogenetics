read('../../noTRuberNoGapsNoAmbiguities.nex')
d = Data()

t = func.randomTree(taxNames=d.taxNames)
t.data = d
t.newComp(free=1, spec='empirical')
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)
t.setModelComponentsOnNodesRandomly()
m = Mcmc(t, nChains=4, runNum=0, sampleInterval=20, checkPointInterval=20000, simulate=3)

m.run(40000)

if 0 and func.which("gnuplot"):
    h = Numbers('mcmc_likes_0', col=1)
    h.plot()

if 1:
    cpm = McmcCheckPointReader()
    cpm.writeProposalAcceptances()
    cpm.writeSwapVectors()
    #cpm.compareSplits(0, 1, verbose=True)
func.summarizeMcmcPrams()

