ls = LeafSupport('../input.nex')

# A clade can be specified using the taxon names in the following
# way. This will produce a list of stabilities using only the quartets
# or triplets that are defined by the taxa in the clade.
ls.defineClade(['Pholiderpeton', 'Proterogyrinus', 'Baphetes',
                'Megalocephalus', 'Loxomma', 'Balanerpeton'
                ,'Dendrerpeton' ,'Gephyrostegus', 'Crassigyrinus'
                ,'Eucritta' ,'Whatcheeria'])


# A broader way of examining clades is to use a definition based on
# the number of input trees the clade appears in.  This way all clades
# that are present in more than 50 % of the input trees will have a
# separate stability list produced.

if 0:
    ls.exploreClades=True
    ls.cladePercentage=50


ls.leafSupport()




