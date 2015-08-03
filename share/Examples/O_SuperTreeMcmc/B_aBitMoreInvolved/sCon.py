tp = TreePartitions("mcmc_trees_0.nex", skip=500)
tp.read("mcmc_trees_1.nex", skip=500)
t = tp.consensus(minimumProportion=0.5)
for n in t.iterInternalsNoRoot():
    n.name = "%.0f" % (100. * n.br.support)
t.draw()
t.name = 'stMcmc'
t.writeNexus('stMcmcCons.nex')
