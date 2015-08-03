read('d3.nex')
a = var.alignments[0]
a.checkForDuplicateSequences(removeDupes=True, makeDict=True)
a.writePhylip(fName='d3_noDupes.phy')


