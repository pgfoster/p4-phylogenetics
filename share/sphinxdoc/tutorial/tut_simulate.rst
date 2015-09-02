===============
Simulating data
===============

The ability to simulate data is useful in many contexts.  You can only
simulate into an existing Data object, and so if you want to simulate a
new data set, you need to make a blank Data object in which to simulate.
You can make one or more new empty alignments and make a Data object
like this::

     a = func.newEmptyAlignment(dataType='dna', taxNames=myTaxNames, length=200)
     d = Data([a])

Then you can attach that Data object to a tree, define a model, and
simulate with the Tree ``simulate()`` method.

This week, it uses the random() function in stdlib.  You can reseed it
by, for example::

     func.reseedCRandomizer(os.getpid())
