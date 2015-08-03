read('protein.nex')
a=var.alignments[0]
a.recodeDayhoff()
a.writeNexus('recoded.nex')
