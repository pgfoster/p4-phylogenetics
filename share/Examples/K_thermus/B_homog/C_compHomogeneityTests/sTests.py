var.verboseRead = 0
read('../../noTRuberNoGapsNoAmbiguities.nex')
d = Data()
read('opt.p4_tPickle')
t = var.trees[0]
t.data = d
t.compoTestUsingSimulations(nSims=100, doIndividualSequences=0, doChiSquare=1, verbose=1)
