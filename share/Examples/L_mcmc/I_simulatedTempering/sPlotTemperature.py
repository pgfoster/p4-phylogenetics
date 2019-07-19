cpr = McmcCheckPointReader(theGlob='*100000')

m = cpr.mm[0]
n = Numbers(m.simTemp_tNums)
n.plot()


