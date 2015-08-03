read('paupConTree.nex')
t = var.trees[0]

# The only point of the data is to get a taxNames list.
read('d.nex')
d = Data()
t.taxNames = d.taxNames

t.readBipartitionsFromPaupLogFile('paupLog')

# The support gets put in node.br.support, as a float from 0-1.  To
# see it in a drawing or write it in newick format, we move the
# support to the node.name, formatting it as percent with no decimal
# places.
for n in t.root.iterInternals():
    if n != t.root:
        if n.br.support:
            n.name = '%i' % round(100. * float(n.br.support))

t.draw(width=80, showNodeNums=0)
t.write()
