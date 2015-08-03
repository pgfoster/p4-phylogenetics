read('sim.p4_tPickle')
t = var.trees[0]
read('d.nex')
t.data = Data()
t.optLogLike(newtAndBrentPowell=True, allBrentPowell=False)
t.model.dump()
