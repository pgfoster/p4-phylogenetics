read('../../noTRuberNoGapsNoAmbiguities.nex')
d = Data()
read('../../tt.nex')

for t in var.trees:
    t.data = d
    c1 = t.newComp(free=1, spec='empirical')
    c2 = t.newComp(free=1, spec='empirical')
    t.newRMatrix(free=1, spec='ones')
    t.setNGammaCat(nGammaCat=4)
    t.newGdasrv(free=1, val=0.5)
    t.setPInvar(free=0, val=0.0)
    t.setModelThing(c1, node=t.root, clade=1)
    t.setModelThing(c2, node=t.node('Deinococcus'))
    t.setModelThing(c2, node=t.node('Bacillus'))
    t.optLogLike()

tt = Trees()
tt.data = d
tt.consel()


