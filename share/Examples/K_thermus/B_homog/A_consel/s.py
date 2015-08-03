read('../../noTRuberNoGapsNoAmbiguities.nex')
d = Data()
read('../../tt.nex')

for t in var.trees:
    t.data = d
    t.newComp(free=1, spec='empirical')
    t.newRMatrix(free=1, spec='ones')
    t.setNGammaCat(nGammaCat=4)
    t.newGdasrv(free=1, val=0.5)
    t.setPInvar(free=0, val=0.0)
    t.optLogLike()

tt = Trees()
tt.data = d
tt.consel()


