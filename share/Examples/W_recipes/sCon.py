# Make a consensus tree, with variations.
mySkip = 500
tp = TreePartitions("mcmc_trees_0.nex", skip=mySkip)

# read another file
tp.read("mcmc_trees_1.nex", skip=mySkip)

t = tp.consensus()

# Re-root
t.reRoot(t.node('myOutTaxon').parent)

# Collapse poorly supported nodes
toCollapse = [n for n in t.iterInternalsNoRoot() if n.br.support < 0.7]
for n in toCollapse:
    t.collapseNode(n)

# Transfer node.br.support to node.name as percent
for n in t.iterInternalsNoRoot():
    n.name = "%.0f" % (100. * n.br.support)

t.draw()
t.writeNexus('cons01_70.nex')
