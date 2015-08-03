ls = LeafSupport('../input.nex')

# Defining a tax set will create a list of stabilities using only
# quartets or triplets defined by the members of set. This allows the
# user to investigate any set of taxa and there relative stabilities.

ls.defineTaxSet(['Eusthenoperon','Baphetes', 'Megalocephalus',
                 'Loxomma','Crassigyrinus' ,'Eucritta' ,'Whatcheeria'])

ls.leafSupport()




