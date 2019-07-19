# Restart an MCMC.
read("d.nex")
d = Data()
m = func.unPickleMcmc(0, d)

# The pseudo prior tuning continues through the run.  Here, samples
# are taken and written.  But note that the samples might be very unevenly spaced.
m.run(100000)

