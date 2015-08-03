tp = TreePartitions("mcmc_trees_0.nex", skip=200)
t = tp.consensus()

# put support on node.name's, for the text drawing
for n in t.iterInternalsNoRoot():
    n.name = "%.0f" % (100. * n.br.support)
t.draw()

# Save it
t.writeNexus(fName='cons.nex')
