# Restart an MCMC from a checkpoint.

# Settings made to the var module in the original MCMC are not in
# the pickled file, so you will probably need to set them again ---
var.PIVEC_MIN = 1.e-6
var.RATE_MIN = 1.e-6
var.BRLEN_MIN = 1.e-5
var.GAMMA_SHAPE_MIN = 0.15

read("d.nex")
d = Data()

# Assuming more than one run, we set the run number by calling this script as
# p4 sRestartMcmc.py -- 0     # eg 0, 1, 2, ...
rNum = int(var.argvAfterDoubleDash[0])

m = func.unPickleMcmc(rNum, d)

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


