# Do an MCMC with the NDCH2 model
# Here I assume only one data partition, part 0

read("d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d

# Fully parameterize the comps
for n in t.nodes:
    c = t.newComp(free=1, spec='empirical', symbol='-')
    t.setModelThing(c, node=n, clade=0)

t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)
 
t.model.parts[0].ndch2 = True
t.model.parts[0].ndch2_writeComps = False # Usually too many



nGen = 1000000
nSamples = 2000
sInterv = int(nGen / nSamples)
cpInterv = int(nGen / 2)

m = Mcmc(t, nChains=4, runNum=0, sampleInterval=sInterv, checkPointInterval=cpInterv)

m.prob.allCompsDir = 0.0
m.prob.ndch2_internalCompsDir = 1.0
m.prob.ndch2_internalCompsDirAlpha = 1.0
m.prob.ndch2_leafCompsDir = 1.0
m.prob.ndch2_leafCompsDirAlpha = 1.0

print(m.prob)

m.run(nGen)


