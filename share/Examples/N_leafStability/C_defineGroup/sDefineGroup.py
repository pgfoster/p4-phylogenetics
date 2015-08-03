ls = LeafSupport('../input.nex')

# A group of taxa is defined like so. Defining a group will create a
# list of stabilities based on the quartets or triplets defined by the
# boundary of the group, i.e. only quartets that have one side in the
# group and the other outside the group will be included. This will
# make it easy to identify any taxa that is has poor membership
# support with the group.

ls.defineGroup(['Pholiderpeton', 'Proterogyrinus', 'Baphetes',
                'Megalocephalus', 'Loxomma', 'Balanerpeton'
                ,'Dendrerpeton' ,'Gephyrostegus', 'Crassigyrinus'
                ,'Eucritta' ,'Whatcheeria'])


ls.leafSupport()




