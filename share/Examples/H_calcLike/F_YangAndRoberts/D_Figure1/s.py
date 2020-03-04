read('../A_NavidiAlignFromYang/navidiSSRNA.nex')
d = Data()
read('t.nex')
#readFile('t31.WithRates.nex')
t=var.trees[0]

t.data = d
t.newComp(free=1, spec='empirical')
t.newComp(free=1, spec='empirical')
t.newComp(free=1, spec='empirical')
t.newComp(free=1, spec='empirical')
t.newComp(free=1, spec='empirical')
t.newComp(free=1, spec='empirical')
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='2p', val=2.0)
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)

for i in range(7):
    t.setModelThing(t.model.parts[0].comps[i], node=i, clade=0)

d.calcUnconstrainedLogLikelihood1()
u = d.unconstrainedLogLikelihood

print("Optimizing.  This will take a few minutes...")
t.optLogLike()
t.write()
t.draw(model=1)
t.model.dump()

#fName = '%s_mt.nex' % tName
#t.writeNexus(fName)
#t.p4ModelsBlock.writeNexusFile(fName, append=1)

#fResults = open('results', 'w')
fResults = sys.stdout
fResults.write('\n\n')
fResults.write(' tree    diff                             kappa  shape\n')
fResults.write('%5s  ' % t.name)
fResults.write('%6.2f   ' % (u - t.logLike))
fResults.write('(logLike = %f)   ' % t.logLike)
fResults.write('%4.2f   ' % t.model.parts[0].rMatrices[0].val)
fResults.write('%4.2f\n'  % t.model.parts[0].gdasrvs[0].val)
if fResults != sys.stdout:
    fResults.close()
