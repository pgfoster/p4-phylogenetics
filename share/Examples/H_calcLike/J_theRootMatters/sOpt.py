var.warnReadNoFile = 0
var.verboseRead = 0
func.reseedCRandomizer(os.getpid())

read('d.nex')
d = Data()
read('((A:0.4, B:0.4), C:0.4, D:0.4);')
t = var.trees[0]
t.data = d

print '==== Homogeneous model ========================='
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=1)
t.setPInvar(free=0, val=0.0)

#t.draw(model=1)
t.optLogLike()

print '\n--- Re-rooting ------'
t.reRoot(1)
#t.draw(model=1)
t.optLogLike()

t.reRoot(0)
t.model = None

print '==== Hetero model ========================='
c1 = t.newComp(free=1, spec='empirical')
c2 = t.newComp(free=1, spec='empirical')
t.setModelThing(c1, node=t.root, clade=1) 
t.setModelThing(c2, node=1, clade=1)
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=1)
t.setPInvar(free=0, val=0.0)

t.draw(model=1)
t.optLogLike()

print '\n--- Re-rooting ------'
t.reRoot(1)
t.draw(model=1)
t.optLogLike()

