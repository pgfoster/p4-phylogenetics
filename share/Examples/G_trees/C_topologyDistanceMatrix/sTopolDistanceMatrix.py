taxNames = list(string.uppercase[:7])
for i in range(6):
    t = func.randomTree(taxNames)
    t.name = 't%i' % (i + 1)
    var.trees.append(t)

tt = Trees(taxNames=taxNames)
dm = tt.topologyDistanceMatrix('wrf')
dm.writeNexus()
