read('sim.p4_tPickle')
t = var.trees[0]
read('d.nex')
t.data = Data()
t.optLogLike()
t.model.dump()
