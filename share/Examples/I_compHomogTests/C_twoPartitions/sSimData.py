read('t.nex')
t = var.trees[0]
var.alignments.append(func.newEmptyAlignment(dataType='dna', symbols=None, taxNames=t.taxNames, length=1000))
var.alignments.append(func.newEmptyAlignment(dataType='dna', symbols=None, taxNames=t.taxNames, length=700))
d = Data()
t.data = d

t.newComp(partNum=0, free=1, spec='specified', val=[0.1, 0.2, 0.3])
t.newRMatrix(partNum=0, free=1, spec='specified', val=[1.2, 6.5, 1.3, 9.8, 1.1, 0.85])
t.setPInvar(partNum=0, free=1, val=0.1)
t.setNGammaCat(partNum=0, nGammaCat=4)
t.newGdasrv(partNum=0, free=1, val=0.5)

t.newComp(partNum=1, free=1, spec='specified', val=[0.5, 0.3, 0.2, 0.0])
t.newRMatrix(partNum=1, free=1, spec='specified', val=[4.3, 2.1, 18.7, 2.3, 1.2, 4.2])
t.setPInvar(partNum=1, free=0, val=0.0)
t.setNGammaCat(partNum=1, nGammaCat=1)

t.simulate()
t.data.alignments[0].sequences[2].sequence = '-' * 1000
#t.data.writeNexus('d.nex')
t.data.alignments[0].writeNexus('dA.nex')
t.data.alignments[1].writeNexus('dB.nex')

t.tPickle('noOpt')
