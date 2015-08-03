tp = TreePartitions('mcmc_trees_0.nex', skip=1000)
t = tp.consensus()

for n in t.iterInternalsNoRoot():
    n.name = '%.0f' % (100.0 * n.br.support)

t.write()
t.draw(width=30, showNodeNums=0)

