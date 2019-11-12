cpr = McmcCheckPointReader(theGlob='*', last=True)

m = cpr.mm[0]
n = Numbers(m.simTemp_tNums)
n.plot()


