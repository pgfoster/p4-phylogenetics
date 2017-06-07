rNum = 0
read('d.nex')
#read('simTree.nex')
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d
var.strictRunNumberChecking = False

pNum = 0
#  comp for each node
for n in t.iterNodes():
    c = t.newComp(partNum=pNum, free=1, spec='empirical', symbol="-")
    t.setModelThing(c, node=n, clade=False)
t.newRMatrix(partNum=pNum, free=0, spec='ones')
t.setNGammaCat(partNum=pNum, nGammaCat=1)
#t.newGdasrv(free=1, val=0.5)
t.setPInvar(partNum=pNum, free=0, val=0.0)
t.model.parts[pNum].ndch2 = True
t.model.parts[pNum].ndch2_writeComps = True

myCmd = "rm -f mcmc*_%i*" % rNum
os.system(myCmd)

sInterv = 100 
cpInterv = None

m = Mcmc(t, nChains=4, runNum=rNum, sampleInterval=sInterv, checkPointInterval=cpInterv)
m.prob.local = 2.
m.prob.eTBR = 2.
m.prob.root3 = 10.0
m.prob.brLen = 0
m.prob.allBrLens = 2.0
m.prob.rMatrix = 0
m.prob.rMatrixDir = 0.0
m.prob.comp = 0
m.prob.compDir = 0.0
m.prob.allCompsDir = 0.0
m.prob.ndch2_leafCompsDir = 0.5
m.prob.ndch2_internalCompsDir = 0.5
#m.prob.ndch2_internalCompsDirAlpha = 1.
#m.prob.ndch2_leafCompsDirAlpha = 1.
#m.prob.gdasrv = 5.0
m.tunings.allBrLens = 0.09
m.tunings.parts[0].ndch2_leafCompsDir = 50000.
m.tunings.parts[0].ndch2_internalCompsDir = 1500.
m.tunings.local = 0.001
m.tunings.chainTemp = 1.0


print(m.prob)
print(m.tunings)

if 1:
    # autoTune was difficult, so it was done "by hand"
    m.run(50000)
    os.system("rm -f mcmc*_%i*" % rNum)
    # m.autoTune(giveUpAfter=15, writeTunings=False)

m.gen = -1
nGen = 100000
sInterv = int(nGen / 200)
cpInterv = int(nGen / 2)
m.sampleInterval=sInterv
m.checkPointInterval=cpInterv
m.run(nGen)


tp = TreePartitions('mcmc_trees_%i.nex' % rNum, skip=100)
t = tp.consensus()
t.writeNexus("cons%i.nex" % rNum)

