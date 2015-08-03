read('inTrees.phy')
inTrees = var.trees
var.trees = []


read('stMcmcCons.nex')
read('mrpMajRuleConsTree.nex')
read('mrpStrictConsTree.nex')

var.trees[1].name = 'mrpMajRule'
var.trees[2].name = 'mrpStrict'

tt = Trees(taxNames=var.trees[0].taxNames)
tt.inputTreesToSuperTreeDistances(inTrees)    
