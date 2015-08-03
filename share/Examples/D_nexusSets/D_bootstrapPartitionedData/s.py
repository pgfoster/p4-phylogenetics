var.doCheckForDuplicateSequences=False
read('d.nex')
a=var.alignments[0]
a.setCharPartition('cp1')
d = Data()
d.alignments[0].writePhylip()

oneBoot = d.bootstrap()
oneBoot.alignments[0].writePhylip()


