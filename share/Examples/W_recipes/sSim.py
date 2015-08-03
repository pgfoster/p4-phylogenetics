# Simulate data
if 0:
    read("d.nex")
    d = Data()
if 1:
    nTax = 5
    taxNames = list(string.uppercase[:nTax])
    a = func.newEmptyAlignment(dataType='dna', taxNames=taxNames, length=200)
    d = Data([a])

if 0:
    read("t.nex")
    t = var.trees[0]
    #t.taxNames = taxNames
if 0:
    read('(B:0.5, ((D:0.4, A:0.3), C:0.5), E:0.5);')
    t = var.trees[0]
    t.taxNames = taxNames
if 1:
    t = func.randomTree(taxNames=taxNames)

t.data = d
t.newComp(free=0, spec='specified', val=[0.1, 0.2, 0.3])
t.newRMatrix(free=0, spec='specified', val=[2., 3., 4., 5., 6., 7.])
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=0, val=0.5)
t.setPInvar(free=0, val=0.2)

func.reseedCRandomizer(os.getpid())
t.simulate()

if 1:
    d.writeNexus('d.nex', writeDataBlock=True)
if 0:
    d.alignments[0].writePhylip("d.phy")

