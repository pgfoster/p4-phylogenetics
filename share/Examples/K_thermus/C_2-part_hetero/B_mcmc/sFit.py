read('../../noTRuberNoGapsNoAmbiguities.nex')
d = Data()

# Unconstrained likelihood
d.calcUnconstrainedLogLikelihood2() # Installs the result in d.unconstrainedLogLikelihood
unc = d.unconstrainedLogLikelihood
n = Numbers('mcmc_sims_0', col=1, skip=1000)
print('Original unconstrained log like is %s' % unc)
print('Here is the posterior predictive distribution:')
n.histo()
n.tailAreaProbability(unc)
print()

# Do bigXSq
ret = d.compoChiSquaredTest(verbose=0)
bigX2 = ret[0][0]
n = Numbers('mcmc_sims_0', col=2, skip=1000)
print('Original X2 is %s' % bigX2)
print('Here is the posterior predictive distribution:')
n.histo()
n.tailAreaProbability(bigX2)

