Doing MRP, Matrix representation / parsimony
============================================

Lets say that you have a number of trees with overlapping taxon sets, and you would like to do an MRP analysis on them to find a supertree.  The first step would be to make an input matrix suitable for a parsimony program like PAUP or TNT.  For this you can use the function :func:`MRP.mrp`.  The input for that function is a list of p4 tree objects, so you might do something like::

  read('myTrees.phy')
  a = mrp(var.trees)
  a.writeNexus('mr.nex')

The file ``mr.nex`` is a matrix representation of the input trees, and is suitable for input to your parsimony search program.

Although it is not needed for the parsimony analysis, each representation of an input tree in ``mr.nex`` has a NEXUS charset showing which sites in the alignment correspond to which input trees.  This allows reconstruction of the input trees from the matrix, as described in :func:`MRP.reverseMrp`.  

