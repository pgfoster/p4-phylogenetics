===============================
Making a tree from compositions
===============================

The composition of sequences can influence the resulting tree, perhaps
via an obvious compositional attraction between unrelated taxa, or
perhaps via a complex interplay of attractions and repulsions.  To
help to visualize this effect, you can make a tree based solely on the
composition of the sequences in an alignment.

For example, lets say that you do a phylogenetic analysis on your
alignment and get a tree with an unexpected grouping.  You speculate
that the unexpected grouping that you see in your tree is due to
compositional effects, and so you make a composition tree as described
here.  If you do not see the same grouping, then you can reject your
speculation that it was due to compositional effects.

We would want to have some indication of the confidence that we have
in the composition tree, and the bootstrap can do that for us.  So we
can bootstrap the original alignment, and for each pseudoreplicate
make a composition tree, and then make a consensus with support values::

    a = func.readAndPop('myData.phy')
    tList = []
    for bNum in range(200):
        b = a.bootstrap()
        d = b.compositionEuclideanDistanceMatrix()
        t = d.bionj()
        tList.append(t)

    tt = Trees(trees=tList, taxNames=a.taxNames)
    #tt.writeNexus('bnj200.nex')  
    tp = TreePartitions(tt)
    t = tp.consensus()
    for n in t.iterInternalsNoRoot():
        n.name = "%.0f" % (n.br.support * 100.) # convert to percent
    t.writeNexus("compTree.nex")

