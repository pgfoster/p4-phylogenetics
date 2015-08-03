read('inTrees.phy')

stm = STMcmc(var.trees, sampleInterval=10, beta=2.0)
stm.run(500)

tp = TreePartitions("mcmc_trees_0.nex", skip=25)
t = tp.consensus(minimumProportion=0.5)
for n in t.iterInternalsNoRoot():
    n.name = "%.0f" % (100. * n.br.support)
t.draw()
