# Simulate some data
p4 sSim.py

# Set up the MCMC with simulated tempering.
# This checkpoints the Mcmc without doing anything.
p4 sInitMcmc.py

# I am not sure how many rounds of trial and error tuning I will need.
# Hopefully this will be enough, I should check at the end before
# going on to the "run".
p4 sTunePseudoPriors.py -- 20

# The run.  Samples are written.  If more samples are needed, this
# script can be re-run.  On-the-fly tuning of temperature pseudo
# priors continues during the run, writing tuning results to
# mcmc_simTemp_0
p4 sMcmcRun.py

# We can use gnuplot to plot the likelihoods vs generation.  Note that
# they are unevenly spaced!
bash gplot.sh

# We can plot the likelihoods again, but this time evenly spaced
p4 sPlotLikelihoods.py

# We can look at which temperature the chain was at for each
# generation.  It should be random, but it may be uneven.
p4 sPlotTemperature.py

# We can make a consensus tree from the checkpoint.
p4 sCons.py


# Optionally ...
# mb326 m.nex      # or whatever your mrbayes is called ...
# bash cleanup.sh

