cpr = McmcCheckPointReader(theGlob='*')
cpr.writeProposalAcceptances()
# How was swapping between heated chains?
cpr.writeSwapVectors()

if 0:
    m = cpr.mm[0]
    print(m.tunings)


