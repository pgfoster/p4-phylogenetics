# Read in the last one ...
cpr = McmcCheckPointReader(theGlob='*', last=True)

m = cpr.mm[0]
t = m.treePartitions.consensus()
for n in t.iterInternalsNoRoot():
    n.name = "%.0f" % (n.br.support * 100)
t.draw()


