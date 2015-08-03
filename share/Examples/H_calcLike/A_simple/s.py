# Read in the tree and the data
read('t.nex')
t = var.trees[0]
read('d.nex')

# Make a Data object and attach it to the tree.
t.data = Data()

# Define the model.  This is the Juke-Cantor model with no among-site
# rate variation.
t.newComp(free=0, spec='equal')
t.newRMatrix(free=0, spec='ones')
t.setPInvar(free=0, val=0.0)
t.setNGammaCat(nGammaCat=1)

# Finally,
t.optLogLike()


