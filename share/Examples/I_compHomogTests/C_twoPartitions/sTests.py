var.verboseRead = 0
var.doCheckForBlankSequences = False
read('dA.nex')
read('dB.nex')
d = Data()
d.compoChiSquaredTest(verbose=1, skipColumnZeros=1, useConstantSites=1, skipTaxNums=[[2], []], getRows=1)

read('opt.p4_tPickle')
t = var.trees[0]
t.data = d
t.compoTestUsingSimulations(nSims=100, doIndividualSequences=1, doChiSquare=0, verbose=1)
