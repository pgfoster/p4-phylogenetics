# Get the test quantity, X^2, from the original data.
read("../../K_thermus/noTRuberNoGapsNoAmbiguities.nex")
d = Data()
ret = d.compoChiSquaredTest()
#print ret
originalStat = ret[0][0]

# Get the sim stats
n = Numbers('mcmc_sims_0', col=1, skip=500)

# Evaluate the tail area probability
n.tailAreaProbability(originalStat)

