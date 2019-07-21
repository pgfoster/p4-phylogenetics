# Simulate some data
p4 sSim.py

# Set up the MCMC with simulated tempering.  Do a pre-run with
# writeSamples turned off.  Then do a long run, collecting samples.
# On-the-fly tuning of temperature pseudo priors continues during the
# run, writing tuning results to mcmc_simTemp_0
p4 sMcmc.py

# We can use gnuplot to plot the likelihoods vs generation.  They may
# be unevenly spaced.
bash gplot.sh

# We can plot the likelihoods again, but this time evenly spaced
p4 sPlotLikelihoods.py

# We can look at which temperature the chain was at for each
# generation.  It should be random, but it may be uneven.
p4 sPlotTemperature.py

# We can make a consensus tree from the last checkpoint.
p4 sCons.py


# Optionally compare with MCMCMC with MrBayes
# mb326 m.nex      # or whatever your mrbayes is called ...


# bash cleanup.sh

