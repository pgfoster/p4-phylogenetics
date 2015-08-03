# Restart an MCMC.
read("d.nex")
d = Data()
m = func.unPickleMcmc(0, d)
if 0:
    m.tunings.chainTemp = 0.15
    m.tunings.relRate = 1.2
    #m.tunings.parts[0].rMatrix = 1000.0
    #m.tunings.parts[0].comp = 50.
    #m.prob.comp = 0
m.run(4000)


