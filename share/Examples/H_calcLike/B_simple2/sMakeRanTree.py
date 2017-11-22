taxNames = list(string.ascii_uppercase[:7])
t = func.randomTree(taxNames)
t.writeNexus('t.nex')
