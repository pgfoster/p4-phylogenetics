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
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=1)
#t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)

nGen = 100000
nCp = 2        # number of checkPoints
cpInterv = int(nGen / nCp)

# We use only one chain.  The sampleInterval is ignored --- Whenever
# the chain has a temperature of zero, it is sampled.  -- Potentially a
# lot of samples!  Setting simTemp and simTempMax turns on simulated
# tempering.  simTemp is the number of temperatures, and needs to be 2
# or more.

m = Mcmc(t, nChains=1, runNum=0, checkPointInterval=cpInterv, simTemp=6, simTempMax=10.0)

# Mystery hack.
m.simTemp_tunePPLong_tunings = [3., 1.]

# A pre-run, writing neither samples nor checkPoint.  Then zero the gen number. 
m.run(20000, writeSamples=False)
m.gen = -1

for n in range(10):
    m.simTemp_trialAndError(6000)

# The main run
m.run(nGen)


