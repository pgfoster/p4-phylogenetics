ls = LeafSupport('../input.nex')

#Setting the proportion of quartets or triplets to sample, range 0.0 - 1.0
ls.useAllQuartets = False
ls.noQuartetsToUse=0.4

ls.leafSupport()

