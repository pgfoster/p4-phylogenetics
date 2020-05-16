# These are to prevent bad matrix exponentiations
var.PIVEC_MIN = 1.e-7
var.RATE_MIN = 1.e-7
var.BRLEN_MIN = 1.e-6
var.GAMMA_SHAPE_MIN = 0.1

# You may be doing more than one run, and if so you will want to set
# the runNum arg when you instantiate the Mcmc object below.  I have
# runNum=rNum below.  If you are using the p4 script you can set the
# variable rNum here with var.argvAfterDoubleDash, by

# rNum = int(var.argvAfterDoubleDash[0])

# otherwise ...
rNum = 0

read("d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d

# Set your model.  This one is a simple tree-homog model.
c = t.newComp(free=1, spec='empirical', symbol='-')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=1)
#t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)

nGen = 200000
nCp = 2        # number of checkPoints
cpInterv = int(nGen / nCp)

# The tempCurveLogBase affects the spacing of the temperatures, which
# in turn affects the logPiDiffs and the acceptances, and affects how
# easy it is to tune it so that you get more-or-less equal occupancies
# in all the temperatures.  Small numbers (eg 1.5) make it more
# linear, bigger numbers (eg 4.) make it more curved.

var.mcmc_simTemp_tempCurveLogBase = 2.0

# We use only one chain.  The sampleInterval is ignored --- Whenever
# the chain has a temperature of zero, it is potentially sampled.  --
# that could be a lot more samples than needed.  To "thin" the
# samples, we can use var.mcmc_simTemp_thinning, a list of digit
# strings to say when to sample.  By default it is ['0'], meaning that
# samples are taken only when the gen number ends with '0'.  To make
# it sample less often, you might set it like this ---

# var.mcmc_simTemp_thinning = ['00', '50']   

# Instantiate an Mcmc object.  Setting simTemp and simTempMax turns on
# simulated tempering.  simTemp is the number of temperatures, and
# needs to be 2 or more.

m = Mcmc(t, nChains=1, runNum=rNum, checkPointInterval=cpInterv, simTemp=8, simTempMax=60.0)

# The trial and error part below might not tune well enough (as shown
# by badly unequal occupancies --- it does not need to be exact).  So
# you might need to do it again, and you can get a head start by using
# the previous longSampleTunings.  They are all 1.0 by default.
# See the files mcmc_simTemp_0 etc for occupancies during the main run.

#m.simTemp_longSampleTunings[0] = 90.0
#m.simTemp_longSampleTunings[1] = 8.0
#m.simTemp_longSampleTunings[-1] = 1.3

# A optional pre-run, without writing samples or a checkPoint.  Then zero the gen number. 
m.run(20000, writeSamples=False)
m.gen = -1

# This adjusts the longSampleTunings by trial-and-error, hopefully so
# that occupancies are more or less equal in all the temperatures.
for rNum in range(20):
    print("-" * 50)
    print("trial and error %i" % rNum)
    m.simTemp_trialAndError(m.simTemp * 1000)

# The main run
m.run(nGen)


