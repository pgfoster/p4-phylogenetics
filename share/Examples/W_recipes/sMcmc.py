# Do an MCMC.

# The defaults for these values are much smaller.  These suggested
# changes may not be needed all the time, but seems to help for
# difficult data.
var.PIVEC_MIN = 1.e-6
var.RATE_MIN = 1.e-6
var.BRLEN_MIN = 1.e-5
var.GAMMA_SHAPE_MIN = 0.15

# Assuming more than one run, we set the run number by calling this script as
# p4 sMcmc.py -- 0     # 0, 1, 2, ...
rNum = int(var.argvAfterDoubleDash[0])

read("d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)

nGen = 1000000
nSamples = 2000
sInterv = int(nGen / nSamples)
nCheckPoints = 2
cpInterv = int(nGen / nCheckPoints)  # or None

m = Mcmc(t, nChains=4, runNum=rNum, sampleInterval=sInterv, checkPointInterval=cpInterv)

m.doCpo = True
m.cpo_startGen = 500000


m.run(nGen)

# func.summarizeMcmcPrams(skip=1000)


