# Restart an MCMC.
read("d.nex")
d = Data()
myRunNum = 0
m = func.unPickleMcmc(myRunNum, d)

# Settings made to the var module in the old MCMC are not present in
# the pickled file, so you may need to set them again ---
# var.PIVEC_MIN = 1.e-7
# var.RATE_MIN = 1.e-7
# var.BRLEN_MIN = 1.e-6
# var.GAMMA_SHAPE_MIN = 0.1

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


