read('inTrees.phy')
inTrees = var.trees
var.trees = []
read('stMcmcCons.nex')
read('mrpMajRuleConsTree.nex')
read('mrpStrictConsTree.nex')
var.trees[1].name = 'mrpMajRule'
var.trees[2].name = 'mrpStrict'


from p4.SuperTreeSupport import SuperTreeSupport
print "%20s  %6s  %6s  %6s  %6s  %6s" % (' ', 'S', 'P', 'Q', 'R', 'V')
for st in var.trees:
    
    sts = SuperTreeSupport(st, inTrees)
    
    sts.doSaveDecoratedTree = False
    sts.decoratedFilename='mytree.nex'
    sts.doSaveIndexTree=False
    sts.indexFilename='mytreeIndex.nex'
    sts.csvFilename='mytreeIndex.csv'
    sts.doDrawTree=False
    sts.verbose=0
    
    sts.superTreeSupport()
    print "%20s  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f" % (st.name, sts.S, sts.P, sts.Q, sts.R, sts.V)
    
