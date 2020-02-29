rNum = 0
tp = TreePartitions('mcmc_trees_%i.nex' % rNum, skip=100)
t = tp.consensus()

# Transfer node.br.support to node.name as percent
for n in t.iterInternalsNoRoot():
    n.name = "%.0f" % (100. * n.br.support)

t.draw()
t.writeNexus("cons%i.nex" % rNum)

