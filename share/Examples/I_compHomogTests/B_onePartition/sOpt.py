read('noOpt.p4_tPickle')
t = var.trees[0]
read('d.nex')
t.data = Data()
t.optLogLike(newtAndBrentPowell=1, allBrentPowell=0, verbose=1)
t.model.dump()
t.tPickle('opt')

