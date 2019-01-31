# Restart an MCMC.
read("d.nex")
d = Data()
m = func.unPickleMcmc(0, d)

if 0:
    # You probably don't want to do this, but 
    # it is possible to change the mcmc here.
    # eg
    # m.sampleInterval = 200

# Restart the run with a set number of generations, eg
# m.run(100000)

# Or set the number of gens in terms of the checkPointInterval
print("The checkpoint interval is currently %i" % m.checkPointInterval)

# For example, say to do two more checkpoints' worth of gens.
# m.run(m.checkPointInterval * 2)


