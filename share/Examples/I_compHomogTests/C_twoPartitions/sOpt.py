read('noOpt.p4_tPickle')
t = var.trees[0]
var.doCheckForBlankSequences = False
read('dA.nex')
read('dB.nex')
t.data = Data()
t.optLogLike(verbose=1)
t.tPickle('opt')

