var.doCheckForBlankSequences=False
var.verboseRead = 0
read('d.nex')
d = Data()
read('opt.p4_tPickle')
t = var.trees[0]
t.data = d
t.simsForModelFitTests(reps=11)
t.modelFitTests()
