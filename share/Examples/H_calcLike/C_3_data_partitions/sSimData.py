read('t.nex')
t = var.trees[0]
a = func.newEmptyAlignment(dataType='dna', symbols=None, taxNames=t.taxNames, length=1500)
b = func.newEmptyAlignment(dataType='standard', symbols='xyz', taxNames=t.taxNames, length=1000)
read('setsForSim.nex')
a.setCharPartition('cp1')

d = Data([a,b])
t.data = d

# Model for Part 0
pNum = 0
t.newComp(partNum=pNum, free=0, spec='specified', val=[0.1, 0.2, 0.3])
t.newRMatrix(partNum=pNum, free=0, spec='specified', val=[1.2, 6.5, 1.3, 9.8, 1.1, 5.4])
t.setPInvar(partNum=pNum, free=0, val=0.2)
t.setNGammaCat(partNum=pNum, nGammaCat=1)
t.setRelRate(partNum=pNum, val=1.2)

# Model for Part 1
pNum = 1
t.newComp(partNum=pNum, free=0, spec='specified', val=[0.2, 0.3, 0.4])
t.newRMatrix(partNum=pNum, free=0, spec='specified', val=[4.2, 1.5, 1.1, 2.3, 10.8, 0.3])
t.setPInvar(partNum=pNum, free=0, val=0.0)
t.setNGammaCat(partNum=pNum, nGammaCat=4)
t.newGdasrv(partNum=pNum, free=0, val=0.5)
t.setRelRate(partNum=pNum, val=5.0)

# Model for Part 2
pNum = 2
t.newComp(partNum=pNum, free=0, spec='specified', val=[0.4, 0.3])
t.newRMatrix(partNum=pNum, free=0, spec='specified', val=[1.2, 2.3, 0.9])
t.setPInvar(partNum=pNum, free=0, val=0.3)
t.setNGammaCat(partNum=pNum, nGammaCat=1)
t.setRelRate(partNum=pNum, val=0.2)

t.simulate()
d.writeNexus('d.nex')

