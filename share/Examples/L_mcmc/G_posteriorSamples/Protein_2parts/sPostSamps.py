from p4.PosteriorSamples import PosteriorSamples

read("d.nex")
read('sets.nex')
a = var.alignments[0]
a.setCharPartition('p1')
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d

pNum=0
t.newComp(partNum=pNum, free=1, spec='wag')
t.newRMatrix(partNum=pNum, free=0, spec='wag')
t.setNGammaCat(partNum=pNum, nGammaCat=4)
t.newGdasrv(partNum=pNum, free=1, val=0.5)
t.setPInvar(partNum=pNum, free=0, val=0.0)
t.setRelRate(partNum=pNum, val=0.5)

pNum = 1
t.newComp(partNum=pNum, free=1, spec='jtt')
t.newRMatrix(partNum=pNum, free=0, spec='jtt')
t.setNGammaCat(partNum=pNum, nGammaCat=4)
t.newGdasrv(partNum=pNum, free=1, val=0.5)
t.setPInvar(partNum=pNum, free=1, val=0.1)
t.setRelRate(partNum=pNum, val=2.0)

t.model.relRatesAreFree = True
t.calcLogLike()

func.reseedCRandomizer(os.getpid())
t.calcLogLike()

#ps = PosteriorSamples(t, runNum=0, program='p4', verbose=3)
#ps = PosteriorSamples(t, runNum=1, program='mrbayes', mbBaseName='mbout', verbose=3)
ps = PosteriorSamples(t, runNum=1, program='mrbayes', mbBaseName='mbout32', verbose=3)

for sampNum in range(0,10):
    t2 = ps.getSample(sampNum)
    t2.data = d
    t2.simulate()
    ret = t2.data.simpleBigXSquared()
    print ret[0], ret[1]
    

