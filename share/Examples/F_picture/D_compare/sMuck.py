read('master.nex')
t = var.trees[0]
for n in t.iterInternalsNoRoot():
    n.name = None

if 1:
    one = t.dupe()
    one.draw()
    one.collapseNode(one.node(36))

