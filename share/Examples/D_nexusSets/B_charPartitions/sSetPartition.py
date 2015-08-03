var.doCheckForDuplicateSequences = False
var.verboseRead = 0

read('ds.nex')

b = var.alignments[1]
b.writePhylip()
b.excludeCharSet('cs2')
b.writePhylip()
b.setCharPartition('cD')
d = Data()
d.dump()

b.setCharPartition(None)
d = Data()
d.dump()
