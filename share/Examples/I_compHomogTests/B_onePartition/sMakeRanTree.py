taxNames = list(string.ascii_uppercase[:4])
t = func.randomTree(taxNames)
t.writeNexus('t.nex')
