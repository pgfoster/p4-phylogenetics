taxNames = list(string.ascii_uppercase[:5])
t = func.randomTree(taxNames)
t.writeNexus('t.nex')
