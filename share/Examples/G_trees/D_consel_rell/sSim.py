var.warnReadNoFile = 0
var.verboseRead = 0

read('((D:0.4, E:0.3):0.03, B:0.5, (C:0.3, A:0.4):0.03);')
taxNames = list(string.uppercase[:5])
var.alignments.append(func.newEmptyAlignment(dataType='dna', taxNames=taxNames, length=500))
d = Data()
t = var.trees[0]
t.data = d
t.taxNames = taxNames
t.newComp(partNum=0, free=0, spec='specified', val=(0.1, 0.2, 0.3))
t.newRMatrix(partNum=0, free=0, spec='specified', val=(1.4, 5.8, 2.1, 18.5, 1.4, 7.2))
t.setNGammaCat(partNum=0, nGammaCat=1)
t.setPInvar(partNum=0, free=0, val=0.1)
t.simulate()
d.writeNexus('d.nex')

