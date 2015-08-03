var.warnReadNoFile = 0
var.verboseRead = 0
read('d.nex')
d = Data()
taxNames = list(string.uppercase[:5])
read('((A, B), C, (D, E));')
read('((D, E), B, (C, A));')
read('((C, E), (B, D), A);')

for i in range(len(var.trees)):
    var.trees[i].name = 't%i' % (i + 1)

for t in var.trees:
    t.taxNames = taxNames
    t.data = d
    t.newComp(partNum=0, free=1, spec='empirical')
    t.newRMatrix(partNum=0, free=1, spec='ones')
    t.setNGammaCat(partNum=0, nGammaCat=1)
    t.setPInvar(partNum=0, free=1, val=0.1)
    t.optLogLike()
    t.tPickle()
