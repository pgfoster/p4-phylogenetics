read('njTree.nex')
t = var.trees[0]
t.restoreDupeTaxa()
t.draw(addToBrLen=0.0)
t.writeNexus(fName='restoredNJTree.nex')
