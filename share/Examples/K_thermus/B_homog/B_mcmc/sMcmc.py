theRunNum = 0

read('../../noTRuberNoGapsNoAmbiguities.nex')
d = Data()

t = func.randomTree(taxNames=d.taxNames)
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)
m = Mcmc(t, nChains=4, runNum=theRunNum, sampleInterval=10, checkPointInterval=5000, simulate=3)

if 0:
    m.tunings.chainTemp = 0.15
    m.tunings.parts[0].rMatrix = 0.3
    m.tunings.parts[0].comp = 0.2
    m.tunings.parts[0].gdasrv = 0.811

if 1:
    m.prob.comp = 1
    m.prob.rMatrix = 1
    m.prob.pInvar = 0
    m.prob.gdasrv = 1
    m.prob.local = 1
    m.prob.root3 = 0
    m.prob.compLocation = 0
    m.prob.rMatrixLocation = 0
    m.prob.relRate = 0

m.autoTune()

m.run(15000)

if 1 and func.which("gnuplot"):
    h = Numbers('mcmc_likes_%i' % theRunNum, col=1)
    h.plot()

if 1:
    cpm = McmcCheckPointReader()
    cpm.writeProposalAcceptances()
    cpm.writeSwapMatrices()
    cpm.compareSplits(1, 2, verbose=True)
func.summarizeMcmcPrams(skip=1000)
