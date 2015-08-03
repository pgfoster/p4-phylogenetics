nTax = 10
taxNames = list(string.uppercase[:nTax])
a = func.newEmptyAlignment(dataType='dna', taxNames=taxNames, length=96)
d = Data([a])
t = func.randomTree(taxNames=taxNames)
t.data = d
t.newComp(free=0, spec='specified', val=[0.1, 0.2, 0.3])
t.newRMatrix(free=0, spec='specified', val=[2., 3., 4., 5., 6., 7.])
t.setNGammaCat(nGammaCat=1)
#t.newGdasrv(free=0, val=0.5)
t.setPInvar(free=0, val=0.0)

func.reseedCRandomizer(os.getpid())
t.simulate()

for i in [2,7]:
    s = a.sequences[0]
    s2 = a.sequences[i]
    s2.sequence = list(s2.sequence)
    for pos in range(len(a)):
        s2.sequence[pos] = s.sequence[pos]
    s2.sequence = ''.join(s2.sequence)

for i in [4,6]:
    s = a.sequences[1]
    s2 = a.sequences[i]
    s2.sequence = list(s2.sequence)
    for pos in range(len(a)):
        s2.sequence[pos] = s.sequence[pos]
    s2.sequence = ''.join(s2.sequence)

if 1:
    a.writeNexus('d3.nex', writeDataBlock=True)
if 0:
    d.alignments[0].writePhylip("d.phy")

