var.verboseRead = 0
var.doCheckForDuplicateSequences=False

read('ds.nex')

if 0:
    var.nexusSets.dump()

if 1:
    a=var.alignments[0]
    a.writePhylip(None)
    c = a.subsetUsingCharSet('cs2', inverse=True)
    c.writePhylip(None)
if 1:
    b = var.alignments[1]
    b.writePhylip(None)
    c = b.subsetUsingCharSet('cs1')
    c.writePhylip(None)
