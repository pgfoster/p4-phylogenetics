read("../d.nex")
d = Data()
m = func.unPickleMcmc(0, d)
m.run(2000)
m = func.unPickleMcmc(1, d)
m.run(2000)

#n = Numbers('mcmc_likes_0', col=1)
#n.plot()
