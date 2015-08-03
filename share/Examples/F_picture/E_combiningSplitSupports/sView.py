read('dB.nex')
a = var.alignments[0]
read('mbout.con')
tMB = var.trees[0]
var.trees.pop()   

read('paupBootTree.nex')
tPAUP = var.trees[1] 

tMB.taxNames = a.taxNames
tPAUP.taxNames = a.taxNames

tMB.tvTopologyCompare(tPAUP)

read('combinedSupportsTree.nex')
tCombined = var.trees[2]
tCombined.tv()



