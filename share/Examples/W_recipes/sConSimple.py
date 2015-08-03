# Make a consensus tree, uncomplicated.
tp = TreePartitions("mcmc_trees_0.nex", skip=500)
t = tp.consensus()
for n in t.iterInternalsNoRoot():
    n.name = "%.0f" % (100. * n.br.support)
t.draw()
t.writeNexus('cons.nex')
