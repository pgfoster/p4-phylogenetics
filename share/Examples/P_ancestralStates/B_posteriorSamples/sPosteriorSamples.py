rNum = 0
read('d.nex')
#read('simTree.nex')
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d
var.strictRunNumberChecking = False

pNum = 0
#  comp for each node
for n in t.iterNodes():
    c = t.newComp(partNum=pNum, free=1, spec='empirical', symbol="-")
    t.setModelThing(c, node=n, clade=False)
t.newRMatrix(partNum=pNum, free=0, spec='ones')
t.setNGammaCat(partNum=pNum, nGammaCat=1)
#t.newGdasrv(free=1, val=0.5)
t.setPInvar(partNum=pNum, free=0, val=0.0)
t.model.parts[pNum].ndch2 = True
t.model.parts[pNum].ndch2_writeComps = True

t.calcLogLike(verbose=False)

# Instantiate
ps = PosteriorSamples(t, runNum=0, program='p4', verbose=0)

# If I were to do simulations, then it would be more complicated because I would
# need to have a separate Data object in which to keep the reference data.
# However, since I am not doing simulations, I am doing ancestral state
# reconstructions on the root node, it is straightforward.  

counts = [0] * 4

# In this example we have 2000 samples.
for sampNum in range(100,200):
    print(sampNum)
    t2 = ps.getSample(sampNum)
    t2.data = d
    asd = t2.ancestralStateDraw()

    # At this point we can do something useful with the ancestral state draw.
    # Or just look at the composition, as is done here.
    for i in range(4):
        ch = 'acgt'[i]
        cnt = asd.count(ch)
        counts[i] += cnt
mySum = float(sum(counts))
print()
for i in range(4):
    print("%.4f" % (counts[i]/mySum))


