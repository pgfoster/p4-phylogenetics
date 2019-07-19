# Read checkPoints from an MCMC.

# Read them in ...
cpr = McmcCheckPointReader(theGlob='*100000')

m = cpr.mm[0]
t = m.treePartitions.consensus()
for n in t.iterInternalsNoRoot():
    n.name = "%.0f" % (n.br.support * 100)
t.draw()


