var.doCheckForDuplicateSequences = False
read('ds.nex')
a=var.alignments[0]
a.writePhylip(None)
b = a.subsetUsingCharSet('cs2')
b.writePhylip(None)

if 0:  # Turn on to see some details
    a.nexusSets.dump()
    
