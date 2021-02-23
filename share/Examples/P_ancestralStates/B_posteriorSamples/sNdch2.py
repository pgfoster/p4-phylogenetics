var.PIVEC_MIN = 1.e-6
var.RATE_MIN = 1.e-6
var.BRLEN_MIN = 1.e-5
var.GAMMA_SHAPE_MIN = 0.15

rNum = 0
read('d.nex')
#read('simTree.nex')
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d

pNum = 0
#  comp for each node
for n in t.iterNodes():
    c = t.newComp(partNum=pNum, free=1, spec='empirical', symbol="-")
    t.setModelComponentOnNode(c, node=n, clade=False)
t.newRMatrix(partNum=pNum, free=0, spec='ones')
t.setNGammaCat(partNum=pNum, nGammaCat=1)       # no asrv for this demo
#t.newGdasrv(free=1, val=0.5)
t.setPInvar(partNum=pNum, free=0, val=0.0)
t.model.parts[pNum].ndch2 = True
t.model.parts[pNum].ndch2_writeComps = True     # often not

myCmd = "rm -f mcmc*_%i*" % rNum
os.system(myCmd)

sInterv = 100 
cpInterv = None

nGen = 100000
nSamples = 200
sInterv = int(nGen / nSamples)
nCheckPoints = 2
cpInterv = int(nGen / nCheckPoints)  # or None

m = Mcmc(t, nChains=4, runNum=rNum, sampleInterval=sInterv, checkPointInterval=cpInterv)
m.run(nGen)


