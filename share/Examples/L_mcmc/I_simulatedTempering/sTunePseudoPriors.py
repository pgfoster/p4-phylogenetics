# How many reps to do ---
nReps = int(var.argvAfterDoubleDash[0])

read("d.nex")
d = Data()
m = func.unPickleMcmc(0, d)


# This is something like a burn-in, where the pseudo priors on the
# temperatures are tuned by trial and error.  We want the occupancy to
# be approximately equal over all temperatures.
for bNum in range(nReps):
    print("-" * 50)
    print("trial and error %i" % bNum)
    m.simTemp_trialAndError(m.simTemp * 5000, verbose=False)

# m.simTemp_trialAndError() resets m.gen to -1.
# This will over-write previous mcmc_checkPoint_0.0 files.
m.checkPoint()
