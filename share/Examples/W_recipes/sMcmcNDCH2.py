# Do an MCMC with the NDCH2 model
# Here I assume only one data partition, part 0

# var.strictRunNumberChecking = False

read("d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d

# Fully parameterize the comps, including a comp for the root
for n in t.nodes:
    c = t.newComp(free=1, spec='empirical', symbol='-')
    t.setModelThing(c, node=n, clade=0)

# adjust these appropriately ...
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)

# Turn on NDCH2
t.model.parts[0].ndch2 = True
t.model.parts[0].ndch2_writeComps = False # Usually too many

# Set up the MCMC
nGen = 1000000
nSamples = 2000
sInterv = int(nGen / nSamples)
nCheckPoints = 2
cpInterv = int(nGen / nCheckPoints)

# Instantiate an Mcmc object.  Adjust nChains and runNum appropriately.
m = Mcmc(t, nChains=4, runNum=0, sampleInterval=sInterv, checkPointInterval=cpInterv)

# Check what proposals are turned on
# print(m.prob)

m.run(nGen)


