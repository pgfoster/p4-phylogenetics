from p4.LeafSupport import CherryRemover

#These three lines removes the cherries and saves the resulting trees in a file. 
rc = CherryRemover('../input.nex')
rc.removeCherries()
rc.saveTrees('input.cherries.removed.nex')

ls = LeafSupport('input.cherries.removed.nex')

ls.defineClade(['Pholiderpeton','Proterogyrinus', 'Baphetes',
                'Megalocephalus', 'Loxomma', 'Balanerpeton-Dendrerpeton' ,'Gephyrostegus',
                'Crassigyrinus' ,'Eucritta' ,'Whatcheeria'])


ls.defineGroup(['Pholiderpeton', 'Proterogyrinus', 'Baphetes',
                'Megalocephalus', 'Loxomma', 'Balanerpeton-Dendrerpeton' ,'Gephyrostegus',
                'Crassigyrinus' ,'Eucritta' ,'Whatcheeria'])

ls.defineTaxSet(['Eusthenoperon','Baphetes', 'Megalocephalus',
                 'Loxomma','Crassigyrinus' ,'Eucritta' ,'Whatcheeria'])


#Set this to true if you want use unresolved quartets or triplets when calculating the stabilities. 
ls.treatUnresolvedAsValid = False

#The score of unresolved quartets or triplets will be divided among the three resolved versions. 
#Setting this option to true will divide the unresolved score equally among the resolutions. 
#Setting it to false will divide the score proportionally according to the score of the resolutions.  
ls.equalDistUnresolved=True


ls.useAllQuartets = False
ls.noQuartetsToUse=0.9

ls.writeCsv = True
ls.csvFilename='myCSV.csv'

ls.leafSupport()




