
===============================
Compositional homogeneity tests
===============================

The usual, poor, way to do a compositional homogeneity test is to do a
Chi-square test on the data.  Paup does this on the data overall, and
Tree-Puzzle does it on individual sequences.  In the latter, since there
are many simultaneous comparisons, the power of the conclusions might be
compromized.  P4 does the test both ways, using the method
:meth:`Data.Data.compoChiSquaredTest`.  

This test uses the X^2 statistic.  It is well-known that this sort of
test suffers from a high probability of type II error, because the
chi-square curve is not an appropriate null distribution by which to
assess significance of the X^2 statistic.

A better null distribution can be obtained by simulating data on the
tree and model in question, and using X^2 statistics extracted from
those simulations to make a null distribution appropriate to the problem
at hand.  This is done using the Tree method
:meth:`Tree.Tree.compoTestUsingSimulations`.

===============
Model fit tests
===============

For Bayesian model fit assessment you can do posterior predictive
simulations during the MCMC.  See :ref:`post-pred-sims-label`.

P4 implements two ML-based model fit tests-- the tree- and model-based
composition fit test, and Goldman's version of the Cox test.  Both
rely on simulations that generally need to be followed by
optimizations, so they are expensive tests.  In both tests the
simulation/optimization part is the time-consuming part, and since
both tests can use the same simulations, in p4 they are done together.
First, the simulations are done using the Tree method
:meth:`Tree.Tree.simsForModelFitTests`.  After that, the simulation
results are digested with the Tree method
:meth:`Tree.Tree.modelFitTests`.


