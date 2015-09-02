=====
Trees
=====


Tree, Node, and Trees classes
-----------------------------

    Trees in p4 are described by a list of nodes, one of which is the root.
    All trees have a root, but this has nothing to do with the number of
    children that the root has.  In phylogenetics when we talk about rooted
    trees we usually mean that we think that one node, the root, is the
    ancestor of all the other nodes; however the root in p4 carries no such
    implication.  The root in p4 is just a node that is a starting place for
    the rest of the nodes.  Can the root ever be a leaf? --- Not usually, but it could be.

    Nodes are the vertexes of the tree, where the branches meet.  (My
    dictionary tells me that both 'vertices' and 'vertexes' are correct plurals
    of vertex.)  P4 describes relationships among nodes with 3 pointers- we
    have a pointer to the parent node, to a child node (which I imagine to
    be always on the left, so it is the leftChild), and to a sibling node
    (which I imagine to always be on the right).  Any of these pointers
    might be None; for example the root node parent is None, the leftChild
    of any tip node is None, and the sibling of the rightmost child of a
    node is None ::

        #    leftChild
        #             \\
        #              \\
        #               Node --- sibling
        #               |
        #               |
        #               parent

    All the nodes except the root in a p4 Tree also have NodeBranch
    attributes, which embody the branches leading to the node from the
    parent.  It has information about for example the length of the branch,
    and so for a node n we have n.br.len.

    The usual way to get around trees is by recursive traversal, and p4
    provides a few avenues for that, both as Node methods and as Tree
    methods.  See Knuth's The Art of Computer Programming for a good
    explanation of preorder and postorder tree traversal.  (Using a previous
    version of p4 that used recursion to read trees, I was given a large
    tree to read, and was surprised to find that it caused p4 to bump into
    the recursion limit of Python.  That limit can be re-set, but instead I
    re-wrote tree reading in p4 in a stack-based approach, rather than a
    recursion-based approach.)

    Sometimes it is natural to deal with trees by the bunch, and so we have
    the Trees class.  It is this class that is able to write Nexus tree
    files that use translations, and it is this class that interfaces with
    consel to compare several trees.


    This week, p4 will read trees in Nexus format, or in Phylip or raw
    Newick format.  The Nexus format is described in the 1997 paper in
    Systematic Biology (MadSwofMad97).  I am not sure where the
    Phylip/Newick format is described.  Older Phylip tree files began with
    a number, indicating the number of trees in the file, but it appears
    that newer Phylip tree files do not have or need that number.  When did
    it change?

    P4 reads both Nexus and Phylip trees with Nexus-like rules, so if you
    have punctuation or spaces in your taxon names then you need to put the
    name in single quotes.  P4 fails to read multi-word Phylip taxon names,
    so for example p4 will gag on the tree::

         ((A, B), (C, Homo sapiens));

    because one of the taxa has a space in it.  That sort of taxon name is
    perfectly acceptable if put inside single quotes.

    The tree files output by Tree-Puzzle have modifications to the
    Phylip/Newick format, where [comments bounded by square brackets]
    precede the tree description.  P4 can handle those.

    Branch lengths go after colons in the tree description, as usual.  You
    can label internal nodes, including the root, immediately after
    unparens, as in the following::

         p4> read("((A, B)the_AB_clade:0.2, (C, 'Homo sapiens')98:0.3)theRoot;")
         p4> t=var.trees[0]
         p4> t.draw()
                          +------2:A
                 +--------1:the_AB_clade
                 |        +------3:B
         theRoot:0
                 |          +------5:C
                 +----------4:98
                            +------6:Homo sapiens

    P4 separates Node objects (vertices) and NodeBranch objects (edges).
    It is the branch that holds the length of the branch.  Branches and branch
    lengths can be accessed by::

      myNode.br
      myNode.br.len


    **After doing a consensus tree, it is the branch that holds the
    support values.**  That makes the branch supports immune to
    re-rooting.  But it also means that in order to see them and save them
    in newick or Nexus format you will need to convert them to
    node.names.  This has the advantage that you can format the conversion
    to your liking --- you can make it percent, or as a a value from zero
    to 1.0 with however many significant figures you wish::

        tp = TreePartitions('myTrees.nex')
        t = tp.consensus()

        # list the support values
        for n in t.iterInternalsNoRoot():
            print "node %3i branch support is %f" % (n.nodeNum, n.br.support)

        # perhaps collapse nodes where the support is too small for you
        toCollapse = [n for n in t.iterInternalsNoRoot() if n.br.support < 0.7]
        for n in toCollapse:
            t.collapseNode(n)

        # optionally re-root it
        t.reRoot(t.node('myOutgroupTaxon').parent)

        # Make the internal node names the percent support values
        for n in t.iterInternalsNoRoot():
            n.name = '%.0f' % (100. * n.br.support)

        # Drawings get a bit messy with both node nums and support
        t.draw(showNodeNums=False)
        t.writeNexus('consTreeWithSupport.nex')

    The Nexus/Newick format can hold node names and branch lengths, but it
    is awkward and non-standard to make it hold more information (eg split
    support, branch colour, model usage information, rooting information,
    and so on).  Perhaps an XML-based format would be useful here, but these
    are early days for XML in phylogenetics, and XML files can get very big
    very quickly.  Update: NeXML looks interesting (files are still big, tho).

    As an alternative, you can store information-rich trees by pickling them, a
    standard Python thing to do to archive objects to files.  Trees, nodes,
    and models may have pointers to c-structures; these should not be
    archived, nor should you archive data objects along with the trees.  To
    facilitate pickling trees you can use the Tree method ``tPickle()``,
    which strips out the c-pointers and data before pickling the tree.
    Trees so pickled (with the p4_tPickle suffix) are autorecognized by p4
    with the ``read()`` function or at the command line.  


Tree pictures and drawings
--------------------------

    You can make a text drawing of trees to the screen with the Tree.draw()
    method.  It provides some control over the presentation, for example the
    width, whether node numbers are displayed, and so on.

    You can make a nice encapsulated postscript drawing of a tree with the
    :meth:`Tree.Tree.eps` method, or an svg drawing with the :meth:`Tree.Tree.svg` method.  While
    they are nice vector graphics, these diagrams are fairly basic, and if
    you want a tree drawing program with more ability then you might
    consider using the Gram package (which uses p4, but Gram is not included in
    p4).  Gram is very flexible, and uses LaTeX for typesetting.

    There is a GUI tree viewer in p4, using the Tree method :meth:`Tree.Tree.tv()`, usable
    with interactive p4.

    Big trees, that are really too big to print, are a special problem for
    both paper and screen.  If you have a tree with 1000 taxa, and each
    taxon is only 1 mm high (too small to read) then the drawing will be 1 m
    on the page or the screen.  If you make the text big enough to read, say
    1cm, then it will be 10 m high!  One solution, that seems to work for
    trees up to about 5K or so taxa, uses the Tree method :meth:`Tree.Tree.btv()` (Big Tree
    Viewer).  This requires a python with Tkinter installed.  This viewer
    is in 2 parts, where in the right panel you can see the whole tree in
    outline with a viewport, and in the right panel you get to see what is
    in that viewport.

    See :ref:`drawing-trees-examples`

Topology distance
-----------------

    You can compare tree topologies such that the root of the tree and any
    rotations of clades do not matter.  For example, these 2 trees have
    identical topologies (by this definition), even though they do not look
    like each other::

         +--------1:A
         |
         |        +--------3:B
         0--------2
         |        +--------4:C
         |
         |        +--------6:D
         +--------5
                  +--------7:E

         +--------1:E
         |
         |--------2:D
         0
         |                  +--------5:C
         |        +---------4
         +--------3         +--------6:B
                  |
                  +---------7:A

    With the :meth:`Tree.Tree.topologyDistance` method you can compare
    topologies without taking branch length differences into account,
    or you can use metrics that do take branch lengths into account.
    The default metric is the symmetric difference, aka the unweighted
    Robinson Foulds distance, which ignores branch lengths.  The
    :meth:`Tree.Tree.topologyDistance` method also provides the
    weighted Robinson-Foulds distance, and the branch length distance,
    which take branch lengths into account.  These are described in
    Felsenstein's book.  To do several trees at once, you can use the
    :meth:`Trees.Trees.topologyDistanceMatrix` method.

Patristic distances
-------------------

    This is just the length along the tree path between all pairs of nodes.
    The method :meth:`Tree.Tree.patristicDistanceMatrix` returns a DistanceMatrix
    object, which you probably want to write to a file.  For example, you
    might say::

         t = var.trees[0]
         dm = t.patristicDistanceMatrix()
         dm.writeNexus('patristic.nex')



