# Simulate data with more than one data partition.
nTax = 5
taxNames = list(string.uppercase[:nTax])

dnaAlign = func.newEmptyAlignment(dataType='dna', taxNames=taxNames, length=300)
read(r"#nexus begin sets; charpartition p1 = first:1-100, second:101-.; end;")
dnaAlign.setCharPartition('p1')
d = Data([dnaAlign])
#d.dump()
read("t.nex")
t = var.trees[0]
#t.taxNames = taxNames

t.data = d

pNum = 0
t.newComp(partNum=pNum, free=0, spec='specified', val=[0.2, 0.1, 0.3])
t.newRMatrix(partNum=pNum, free=0, spec='specified', val=[2., 3., 4., 5., 6., 7.])
t.setNGammaCat(partNum=pNum, nGammaCat=4)
t.newGdasrv(partNum=pNum, free=0, val=1.2)
t.setPInvar(partNum=pNum, free=0, val=0.2)
t.setRelRate(partNum=pNum, val=1.5)

pNum = 1
t.newComp(partNum=pNum, free=0, spec='specified', val=[0.45, 0.05, 0.45])
t.newRMatrix(partNum=pNum, free=0, spec='specified', val=[1.2, 3.4, 0.4, 5.9, 0.6, 1.0])
t.setNGammaCat(partNum=pNum, nGammaCat=1)
#t.newGdasrv(partNum=pNum, free=0, val=0.5)
t.setPInvar(partNum=pNum, free=0, val=0.35)
t.setRelRate(partNum=pNum, val=0.3)

func.reseedCRandomizer(os.getpid())
t.simulate()
d.writeNexus('d.nex')
