skip=200
tp = TreePartitions("mcmc_trees_0.nex", skip=skip)
tp.read("mcmc_trees_1.nex", skip=skip)
tp.dump()
t = tp.consensus()

# put support on node.name's, for the text drawing
for n in t.iterInternalsNoRoot():
    n.name = "%.0f" % (100. * n.br.support)
t.draw()

# Save it
t.writeNexus(fName='cons.nex')
