# Simulate very hetero data

# Below you specify the dataType, length, and relRate of each data
# partition (and symbols, if it is standard dataType).  There is one
# part per alignment, in this example.  You also specify the number of
# taxa in the tree

import random
import math

def randomNumbersThatSumTo1(length, minimum):
    maximum = 1.0 - ((length - 1) * minimum)
    myRange = maximum - minimum
    rnn = [random.random() for i in range(length)]
    mySum = sum(rnn)
    rnn = [minimum + (rn / (mySum / myRange)) for rn in rnn]
    assert math.fabs(sum(rnn) - 1.0) < 1.0e-14
    assert min(rnn) >= minimum
    return rnn

class MySPart(object):
    def __init__(self, dataType, sLen, relRate, symbols=None):
        self.dataType = dataType
        self.sLen = sLen
        self.relRate = relRate
        self.symbols = None
        if dataType == 'standard':
            assert symbols
            self.symbols = symbols

myParts = []
myParts.append(MySPart('dna', 231, 4.6))
myParts.append(MySPart('protein', 411, 0.3))
myParts.append(MySPart('standard', 197, 2.7, symbols='123456'))

nTax = 17
taxNames = ['t%i' % i for i in range(nTax)]
aa = []
for mP in myParts:
    if mP.dataType != 'standard':
        a = func.newEmptyAlignment(dataType=mP.dataType, taxNames=taxNames, length=mP.sLen)
    else:
        a = func.newEmptyAlignment(dataType=mP.dataType, symbols=mP.symbols, taxNames=taxNames, length=mP.sLen)
    aa.append(a)
d = Data(aa)
t = func.randomTree(taxNames=taxNames)
t.data = d
for aNum in range(len(aa)):
    a = aa[aNum]
    mP = myParts[aNum]
    for myNode in t.iterNodes():
        newVal = randomNumbersThatSumTo1(a.dim, var.PIVEC_MIN)
        mt = t.newComp(partNum=aNum, free=1, spec='specified', val=newVal)
        t.setModelThing(mt, myNode, clade=False)

    if a.dataType != 'protein':
        for myNode in t.iterNodesNoRoot():
            nRates = ((a.dim * a.dim) - a.dim) / 2
            newVal = randomNumbersThatSumTo1(nRates, var.RATE_MIN)
            mt = t.newRMatrix(partNum=aNum, free=1, spec='specified', val=newVal)
            t.setModelThing(mt, myNode, clade=False)
    else:
        t.newRMatrix(partNum=aNum, free=0, spec='wag')
    t.setNGammaCat(partNum=aNum, nGammaCat=4)
    t.newGdasrv(partNum=aNum, free=1, val=0.5)
    t.setPInvar(partNum=aNum, free=0, val=0.0)
    t.setRelRate(partNum=aNum, val=mP.relRate)

t.model.relRatesAreFree = True
func.reseedCRandomizer(os.getpid())
t.simulate()
t.name = 'simTree'

d.writeNexus('d.nex', writeDataBlock=True)
t.tPickle()

