read("mcmc_trees_0.nex")
tt = Trees()
# Use a previously-made cons tree 
var.trees = []
read("cons.nex")
t = var.trees[0]
tt.trackSplitsFromTree(t)
