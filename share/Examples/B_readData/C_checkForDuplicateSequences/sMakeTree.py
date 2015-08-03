read("d3_noDupes.phy")
a = var.alignments[0]
dm = a.pDistances()
t = dm.njUsingPaup()
t.draw(addToBrLen=0.0)
t.writeNexus('njTree.nex')
