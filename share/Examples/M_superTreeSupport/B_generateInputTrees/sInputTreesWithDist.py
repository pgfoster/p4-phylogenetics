from p4.SuperTreeSupport import SuperTreeInputTrees

stit = SuperTreeInputTrees('balanced64.nex')
stit.writeInputTreesToFile = True
stit.outputFile = 'inputtreesWithDist.tre'
stit.useTaxonDistribution = True 
stit.noOutputTrees = 20
stit.generateInputTrees()
