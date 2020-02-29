var.verboseRead = 0
read('t.nex')
t = var.trees[0]
read('d.nex')
var.alignments[0].setCharPartition('cp1')
t.data = Data()

# Model for Part 0
pNum = 0
t.newComp(partNum=pNum, free=1, spec='empirical')
t.newRMatrix(partNum=pNum, free=1, spec='ones')
t.setPInvar(partNum=pNum, free=1, val=0.1)
t.setNGammaCat(partNum=pNum, nGammaCat=1)
t.setRelRate(partNum=pNum, val=0.5)

# Model for Part 1
pNum = 1
t.newComp(partNum=pNum, free=1, spec='empirical')
t.newRMatrix(partNum=pNum, free=1, spec='ones')
t.setPInvar(partNum=pNum, free=0, val=0.0)
t.setNGammaCat(partNum=pNum, nGammaCat=4)
t.newGdasrv(partNum=pNum, free=1, val=2.0)
t.setRelRate(partNum=pNum, val=1.2)

# Model for Part 2
pNum = 2
t.newComp(partNum=pNum, free=1, spec='empirical')
t.newRMatrix(partNum=pNum, free=1, spec='ones')
t.setPInvar(partNum=pNum, free=1, val=0.1)
t.setNGammaCat(partNum=pNum, nGammaCat=1)
t.setRelRate(partNum=pNum, val=0.1)

t.model.relRatesAreFree = 1

print('Optimizing...') 
t.optLogLike(newtAndBrentPowell=0, allBrentPowell=1)
#t.optLogLike(newtAndBrentPowell=1, allBrentPowell=0)
t.tPickle('opt')


