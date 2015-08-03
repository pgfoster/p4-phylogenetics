read('../A_NavidiAlignFromYang/navidiSSRNA.nex')
read('t.nex')

d=Data()
d.calcUnconstrainedLogLikelihood1()
u = d.unconstrainedLogLikelihood


#fResults = open('results', 'w')
fResults = sys.stdout
fResults.write('\nTable 2\n=======\n\nHKY, no rates\n\n')
fResults.write(' tree     diff                             kappa  \n')


for t in var.trees:
    t.data = d
    t.newComp(free=1, spec='empirical')
    t.newRMatrix(free=1, spec='2p', val=2.0)
    t.setPInvar(free=0, val=0.0)
    t.optLogLike(verbose=0)
    fResults.write('%5s  ' % t.name)
    fResults.write('%7.2f   ' % (u - t.logLike))
    fResults.write('(logLike = %f)   ' % t.logLike)
    fResults.write('%4.2f   ' % t.model.parts[0].rMatrices[0].val)
    fResults.write('\n')


fResults.write('\n\nHKY, with dG rates\n\n')
fResults.write(' tree     diff                             kappa  shape\n')


for t in var.trees:
    t.model = None
    t.newComp(free=1, spec='empirical')
    t.newRMatrix(free=1, spec='2p', val=2.0)
    t.setNGammaCat(nGammaCat=4)
    t.newGdasrv(free=1, val=0.5)
    t.setPInvar(free=0, val=0.0)
    t.optLogLike(verbose=0)
    fResults.write('%5s  ' % t.name)
    fResults.write('%7.2f   ' % (u - t.logLike))
    fResults.write('(logLike = %f)   ' % t.logLike)
    fResults.write('%4.2f   ' % t.model.parts[0].rMatrices[0].val)
    fResults.write('%4.2f   ' % t.model.parts[0].gdasrvs[0].val)
    fResults.write('\n')

if fResults != sys.stdout:
    fResults.close()

