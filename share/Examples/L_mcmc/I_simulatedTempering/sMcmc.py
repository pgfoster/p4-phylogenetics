# Hopefully these will not be needed, but in the event of 
# a bad matrix exponentiation these might be useful
# var.PIVEC_MIN = 1.e-7
# var.RATE_MIN = 1.e-7
# var.BRLEN_MIN = 1.e-6
# var.GAMMA_SHAPE_MIN = 0.1

read("d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d


c = t.newComp(free=1, spec='empirical', symbol='-')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=1)
#t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)

nGen = 20000
nCp = 1        # number of checkPoints
cpInterv = int(nGen / nCp)

# The tempCurveLogBase affects the spacing of the temperatures, which
# in turn affects the logPiDiffs and the acceptances.
var.mcmc_simTemp_tempCurveLogBase = 5.0

# We use only one chain.  The sampleInterval is ignored --- Whenever
# the chain has a temperature of zero, it is sampled.  -- Potentially a
# lot of samples!  Setting simTemp and simTempMax turns on simulated
# tempering.  simTemp is the number of temperatures, and needs to be 2
# or more.

m = Mcmc(t, nChains=1, runNum=0, checkPointInterval=cpInterv, simTemp=8, simTempMax=60.0)

# The trial and error part below might not tune well enough (as shown
# by badly unequal occupancies --- it does not need to be exact).  So
# you might need to do it again, and you can get a head start by using
# the previous longSampleTunings.  They are all 1.0 by default.
#m.simTemp_longSampleTunings[0] = 90.0
#m.simTemp_longSampleTunings[1] = 8.0
#m.simTemp_longSampleTunings[-1] = 1.3

# A pre-run, without writing samples or a checkPoint.  Then zero the gen number. 
m.run(20000, writeSamples=False)
m.gen = -1

# This adjusts the the longSampleTunings, hopefully so that
# occupancies are more or less equal in all the temperatures.
for rNum in range(10):
    print("-" * 50)
    print("trial and error %i" % rNum)
    m.simTemp_trialAndError(m.simTemp * 1000)

# The main run
m.run(nGen)


