cpr = McmcCheckPointReader()
cpr.writeProposalAcceptances()
cpr.writeSwapVectors()
#cpr.writeProposalProbs()
m = cpr.mm[0]
m.tunings.dump(advice=False)
cpr.compareSplits(2, 3)
print("\n\nComparing all splits from all pairs of checkPoints ...")
cpr.compareSplitsAll()

