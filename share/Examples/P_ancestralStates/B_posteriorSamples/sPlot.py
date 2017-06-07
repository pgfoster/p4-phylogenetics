rNum = 0
sk = 100

n = Numbers('mcmc_likes_%i' % rNum, col=1, skip=sk)
n.plot()

n = Numbers('mcmc_hypers_%i' % rNum, col=2, skip=sk)
n.plot()

