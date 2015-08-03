from p4.PosteriorSamples import PosteriorSamples

read("d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='specified', val=[2., 3., 4., 5., 6., 7.])
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=0, val=0.5)
t.setPInvar(free=0, val=0.0)

func.reseedCRandomizer(os.getpid())
t.calcLogLike()

ps = PosteriorSamples(t, runNum=0, program='p4', mbBaseName='mbout', verbose=3)
#ps = PosteriorSamples(t, runNum=1, program='mrbayes', mbBaseName='mbout', verbose=3)
#ps = PosteriorSamples(t, runNum=1, program='mrbayes', mbBaseName='mbout32', verbose=3)

# Iterate over samples
for sampNum in range(0,10):
    t2 = ps.getSample(sampNum)
    t2.data = d
    t2.simulate()
    ret = t2.data.simpleBigXSquared()
    print ret[0]
    

