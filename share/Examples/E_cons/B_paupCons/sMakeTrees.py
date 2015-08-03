nTax = int(var.argvAfterDoubleDash[0])
nTrees = int(var.argvAfterDoubleDash[1])

#taxNames = list(string.uppercase[:20])
t = func.randomTree(nTax=nTax)
#t.draw()
a = func.newEmptyAlignment(dataType='dna', taxNames=t.taxNames, length=10)
t.data = Data([a])
t.newComp(spec='equal')
t.newRMatrix()
t.setPInvar()
t.setNGammaCat()
t.simulate()
t.data.writeNexus('d.nex') # For paup, and for the taxNames
t.calcLogLike()

# Subvert an mcmc object to generate related trees.
m = Mcmc(t, nChains=1, runNum=0, sampleInterval=1, checkPointInterval=1, simulate=None, verbose=True)
m.run(1)

# Find the 'local' proposal; give it a name.
for p in m.proposals:
    if p.name == 'local':
        break

tList = []
for i in range(nTrees):
    m.chains[0].proposeLocal(p)
    #m.chains[0].propTree.draw()
    t = m.chains[0].propTree.dupe()
    t.model = None
    t.data = None
    t.logLike = None
    t.recipWeight = 1
    tList.append(t)
    
tt = Trees(trees=tList)
tt.writeNexus('tt.nex', withTranslation=1)
