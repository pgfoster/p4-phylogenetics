var.warnReadNoFile = 0
var.verboseRead = 0
read('protein.nex')
a = var.alignments[0]
a.recodeDayhoff()
#a.writePhylip(None)
read('(((A:0.4, (B:0.4, C:0.4):0.05):0.05, D:0.4):0.1, E:0.4, F:0.4);')
t = var.trees[0]
t.data = Data()
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=1)
t.setPInvar(free=0, val=0)
t.optLogLike(newtAndBrentPowell=1, allBrentPowell=0)
t.write()
t.model.dump()
