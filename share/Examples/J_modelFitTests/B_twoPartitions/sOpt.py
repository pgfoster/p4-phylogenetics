var.doCheckForBlankSequences=False
read('noOpt.p4_tPickle')
t = var.trees[0]
read('d.nex')
t.data = Data()
t.optLogLike(verbose=1)
t.tPickle('opt')

