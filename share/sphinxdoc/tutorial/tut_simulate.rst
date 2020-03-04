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

Generation of random numbers uses the GSL random number
generator.  The state is held in var.gsl_rng, which is None by
default.  If you do a simulation using this method, it will
use ``var.gsl_rng`` if it exists, or make it if it does not exist
yet.  When it makes it, it seeds the state based on the
current time.  That should give you lots of variation in the
simulations.

If on the other hand you want to make simulations that are the
same you can reseed the randomizer with the same seed whenever
you do it, like this::

    if not var.gsl_rng:
        var.gsl_rng = pf.gsl_rng_get()
    # unusually, set the seed with each simulation
    mySeed = 23    # your chosen int seed
    pf.gsl_rng_set(var.gsl_rng, mySeed)
