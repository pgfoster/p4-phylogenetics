# We need to have a tree with its data.  The tree needs a model, and
# it needs to have model parameters optimized.

# Read in the tree and give it a name
read("myTree.nex")
t = var.trees[0]

# Read in the data and give it a name
read("myData.nex")
d = Data()

# Attach the data and a model to the tree.
t.data = d
t.newComp(free=0, spec='empirical')
t.newRMatrix(free=0, spec='rtRev')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)

# optimize
t.optLogLike()

# Do the test on the optimized tree.  The compoTestUsingSimulations
# shows that the data fail the test (P=0.0), and so the comps are
# heterogeneous.  The chi-squared test thinks that the sequences are
# homogeneous (P=1.0)
t.compoTestUsingSimulations(nSims=100, doIndividualSequences=0, doChiSquare=1, verbose=1)


