To collapse nodes in a tree
---------------------------

To collapse nodes, don't do it this way::

    for n in t.iterInternalsNoRoot():
        t.collapseNode(n)

because it modifies the tree as it iterates, and so it does not
work. Do it this way instead::

    # First make a list of nodes that you want to collapse
    toCollapse = [n for n in t.iterInternalsNoRoot()]
    
    # Then collapse them
    for n in toCollapse:
        t.collapseNode(n)

Or perhaps more realistically::
    
    toCollapse = [n for n in t.iterInternalsNoRoot() if n.br.support < 0.7]
    for n in toCollapse:
        t.collapseNode(n)

