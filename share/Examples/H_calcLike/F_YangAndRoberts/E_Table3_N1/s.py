read('../A_NavidiAlignFromYang/navidiSSRNA.nex')
read('t.nex')

d = Data()
d.calcUnconstrainedLogLikelihood1()
u = d.unconstrainedLogLikelihood

#fResults = open('results', 'w')
fResults = sys.stdout
fResults.write('\nYang and Roberts Table 3, HKY + dG + N1\n\n')
fResults.write('tree   diff   kappa    shape\n')


for tNum in [0, 2, 14]:
#for tNum in range(len(var.trees)):  # if you are patient ...
    t=var.trees[tNum]
    #print 'Doing tree %s...' % t.name
    #t.draw()
    t.data = d

    if t.name != 't0':
        cInternals = t.newComp(free=1, spec='empirical')
    cRoot = t.newComp(free=1, spec='empirical')
    cS_solfat = t.newComp(free=1, spec='empirical')
    cH_salin = t.newComp(free=1, spec='empirical')
    cE_coli = t.newComp(free=1, spec='empirical')
    cH_sapiens = t.newComp(free=1, spec='empirical')

    if t.name != 't0':
        t.setModelThing(cInternals, node=t.root, clade=1)
    t.setModelThing(cRoot, node=t.root, clade=0)
    t.setModelThing(cS_solfat, node=t.node('S_solfat'))
    t.setModelThing(cH_salin, node=t.node('H_salin'))
    t.setModelThing(cE_coli, node=t.node('E_coli'))
    t.setModelThing(cH_sapiens, node=t.node('H_sapiens'))

    t.newRMatrix(free=1, spec='2p', val=2.0)
    t.setNGammaCat(nGammaCat=4)
    t.newGdasrv(free=1, val=0.5)
    t.setPInvar(free=0, val=0.0)    
    t.optLogLike(verbose=0)
    fResults.write('%4s  ' % t.name)
    fResults.write('%6.2f   ' % (u - t.logLike))
    #fResults.write('(logLike = %f)   ' % t.logLike)
    fResults.write('%4.2f    ' % t.model.parts[0].rMatrices[0].val)
    fResults.write('%4.2f\n'  % t.model.parts[0].gdasrvs[0].val)

if fResults != sys.stdout:
    fResults.close()
