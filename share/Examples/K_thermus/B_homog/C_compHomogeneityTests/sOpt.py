read('../../tt.nex')
t = var.trees[1] # Under this model, the attract tree is the ML tree
read('../../noTRuberNoGapsNoAmbiguities.nex')
t.data = Data()
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)
t.optLogLike(newtAndBrentPowell=1, allBrentPowell=0, verbose=1)
t.model.dump()
t.tPickle('opt')

