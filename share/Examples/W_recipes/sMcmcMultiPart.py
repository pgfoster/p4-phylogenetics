# Do an MCMC with more than one data partition.

# The defaults for these values are much smaller.  These suggested
# changes may not be needed all the time, but seems to help for
# difficult data.
var.PIVEC_MIN = 1.e-6
var.RATE_MIN = 1.e-6
var.BRLEN_MIN = 1.e-5
var.GAMMA_SHAPE_MIN = 0.15

# Assuming more than one run, we set the run number by calling this script as
# p4 sMcmcMultiPart.py -- 0     # 0, 1, 2, ...
rNum = int(var.argvAfterDoubleDash[0])

read("d.nex")
a = var.alignments[0]
a.setCharPartition('p1')
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d

pNum = 0
t.newComp(partNum=pNum, free=1, spec='empirical')
t.newRMatrix(partNum=pNum, free=1, spec='ones')
t.setNGammaCat(partNum=pNum, nGammaCat=4)
t.newGdasrv(partNum=pNum, free=1, val=0.5)
t.setPInvar(partNum=pNum, free=1, val=0.1)
t.setRelRate(partNum=pNum, val=5.0)

pNum = 1
t.newComp(partNum=pNum, free=1, spec='empirical')
t.newRMatrix(partNum=pNum, free=1, spec='ones')
t.setNGammaCat(partNum=pNum, nGammaCat=1)
#t.newGdasrv(partNum=pNum, free=0, val=0.5)
t.setPInvar(partNum=pNum, free=1, val=0.1)
t.setRelRate(partNum=pNum, val=2.0)

t.model.relRatesAreFree = True

nGen = 1000000
nSamples = 2000
sInterv = int(nGen / nSamples)
nCheckPoints = 2
cpInterv = int(nGen / nCheckPoints)  # or None

m = Mcmc(t, nChains=4, runNum=rNum, sampleInterval=sInterv, checkPointInterval=cpInterv)
m.run(nGen)
# func.summarizeMcmcPrams(skip=200)


