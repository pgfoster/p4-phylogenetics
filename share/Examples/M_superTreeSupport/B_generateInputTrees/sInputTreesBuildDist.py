from p4.SuperTreeSupport import SuperTreeInputTrees

stit = SuperTreeInputTrees('balanced64.nex', distributionTrees='FelidaeRVS.tre')
stit.writeInputTreesToFile = True
stit.outputFile = 'inputtreesBuiltDist.tre'
stit.noOutputTrees = 20
stit.generateInputTrees()


