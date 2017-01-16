# Read checkPoints from an MCMC.

# Read them in ...
cpr = McmcCheckPointReader(theGlob='*')

# A table of acceptances
cpr.writeProposalAcceptances()

# Compare splits using average std deviation of split frequencies
#cpr.compareSplitsAll()
# or between only two checkpoints ...
#cpr.compareSplits(0, 1)

# How was swapping between heated chains?
#cpr.writeSwapMatrices()

#cpr.writeProposalProbs()

# Get a single Mcmc object ..
if 0:
    m = cpr.mm[0]
    print(m.tunings)


