var.verboseRead = 0
read('d.nex')
d = Data()
d.compoChiSquaredTest(verbose=1, skipColumnZeros=1, useConstantSites=1, skipTaxNums=None, getRows=0)

read('opt.p4_tPickle')
t = var.trees[0]
t.data = d
t.compoTestUsingSimulations(nSims=100, doIndividualSequences=0, doChiSquare=0, verbose=1)
