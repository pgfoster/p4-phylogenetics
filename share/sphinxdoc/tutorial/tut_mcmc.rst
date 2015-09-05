===========================
Bayesian analysis with MCMC
===========================

P4 has a basic MCMC for doing Bayesian analyses.  It uses the
Metropolis-coupled MCMC, or MCMCMC, that has worked so well in MrBayes.
To start an MCMC, you first set up your data, tree, and model as usual.
Then you instantiate an Mcmc object, and run the chain, like this::

     # t is the tree, with model and data attached
     m = Mcmc(t, nChains=4, runNum=0, sampleInterval=1000, checkPointInterval=250000, simulate=2)
     m.run(1000000) # one million generations

Here nChains says to make 4 MCMC chains in parallel in an MCMCMC, where
1 chain is the cold chain, and 3 are heated chains.  See the
documentation for MrBayes for more explanation.

As the MCMC runs, it saves trees to a tree file, saves likelihoods to
another file, and saves model parameters to another file.  These samples
need to be digested or summarized later.

The Mcmc class has a big doc string.  Say help(Mcmc).


A simple MCMC
-------------

::

     # You will usually start with a data file, and a random tree.
     read("d.nex")
     d = Data()
     t = func.randomTree(taxNames=d.taxNames)
     t.data = d

     # Now define a model, as usual.  This one is GTR+IG.
     t.newComp(free=1, spec='empirical')
     t.newRMatrix(free=1, spec='ones')
     t.setNGammaCat(nGammaCat=4)
     t.newGdasrv(free=1, val=0.5)
     t.setPInvar(free=1, val=0.2)

     # Now instantiate an Mcmc object.  This one has 4 MCMCMC chains, is set
     # to do a million generations (below, in run()), collect 2000 samples
     # (every 500 generations) and is set to take 4 checkPoints.

     m = Mcmc(t, nChains=4, runNum=0, sampleInterval=500, checkPointInterval=250000)

     # Here you can adjust the proposal probabilities, and tuning parameters.
     m.prob.polytomy = 1.0       # Default is zero
     m.prob.brLen = 0.001        # Default is zero
     m.tunings.chainTemp = 0.12  # Default is 0.15, this week

     # This starts the run, specifying the number of gens you want.  You can
     # restart it from a checkPoint if you find you need more.
     m.run(1000000)


Proposal probabilities
----------------------

Proposals to change the model parameters or tree topology are made with
default probabilities.  These default probabilities are calculated based
on the number of parameters involved.  You can adjust these proposal
probabilities, or perhaps turn a proposal type off altogether, by
setting attributes of Mcmc.prob.  These are user-specifiable relative
probabilities, with default being 1; they are not proper probabilities,
and do not sum to 1. 

For an Mcmc object ``m``, I can show these current relative probs by
saying::

     print m.prob

To show the real intended rates (the kind that add up to 1.0) that are
calculated from these relative probs, you can do a::

     m.writeProposalIntendedProbs()

or, after a run (or on checkpoints, see below)::

     m.writeProposalProbs()

For my Mcmc object ``m``, I could turn off the tree topology changes by
setting local to zero.  I can still allow branch lengths to be free by
turning on the brLens proposal, which is off by default. ::

     m.prob.local = 0
     m.prob.brLen = 1.0

These relative m.prob values affect proposal probabilities, but not in a
linear way.  For changes with inherently low proposal probabilities,
such as ``relRate``, doubling the relative prob approximately doubles the
final proposal probability.  As I mentioned, you can check the final
calculated proposal probabilities with ``m.writeProposalProbs()``.


The polytomy proposal
---------------------

This is an implementation of the move described in Lewis et al (2005)
"Polytomies and Bayesian phylogenetic inference" Syst Biol 54:241-253.
This move is turned off by default, but you can turn it on for an Mcmc
object m by saying::

     m.prob.polytomy = 1.0
     m.prob.brLen = 0.001

(If you have local turned on, which you probably do, you will also want
to turn on the brLen proposal in addition, for those cases when the
current tree is a star tree; the local proposal needs 3 contiguous
edges, and so the local proposal does not work on a star tree.  In that
case, the chain uses the brLen proposal instead of local.)

This move tends to decrease spurious high support for short internal
nodes.  The prior probability favouring increased multifurcations or
increased resolution can be adjusted.  Actually, as in the paper, there
are two ways to do this; you can use the polytomy prior, or you can use
the resolution class prior.  If you want to use the resolution class
prior, you can turn it on, for an Mcmc object m, by::

     m.tunings.doPolytomyResolutionClassPrior = True  # False by default

but if it is left as False then you use the polytomy prior.  In both
cases, there is a 'C' constant to set, that determines how much
multifurcations are favoured.  You can set it via its log, for example
for an Mcmc object m ::

     m.tunings.polytomyPriorLogBigC = 1.0   # Default is zero

By setting the logBigC to 1.0, above, I am setting the C to e as in the
paper.  Assuming the default of using polytomy priors, leaving logBigC
set to zero (ie C=1), the default, means that every tree (of any
resolution) has the same prior probability.  Setting it to higher
numbers favours polytomies.

As in the paper, you can use the resolution class priors, and set the C
constant to 1.1 (actually 10/9 = 1.111111), 2, or 10, for an Mcmc
object m, like this::

     import math
     m.tunings.doPolytomyResolutionClassPrior = True
     m.tunings.polytomyPriorLogBigC = math.log(10./9.)   # Or 2, or 10, or ...


Topology constraints
--------------------

You can specify constraints on the tree topology with a Constraints
object.  To make a Constraints object, you need a list of taxNames (in
the same order as your data, trees, and so on), and you also need a
(usually partly-resolved) tree object::

     myTaxNames = myData.taxNames
     read('myConstraintsTreeFile.phy')
     myConstraintsTree = var.trees.pop()
     myConstraints = Constraints(myTaxNames, myConstraintsTree)

You can pass a Constraints object to func.randomTree() and Mcmc() to
enforce constraints.  If you are starting a Mcmc with a randomTree,
then it should have the same constraints as you pass to the Mcmc::

     t = func.randomTree(taxNames=myTaxNames, constraints=myConstraints)
     m = Mcmc(t, ..., constraints=myConstraints)


Tuning the MCMC
---------------

Consider a proposal to change a parameter in a model.  There is a
current state, and a new state is proposed somewhere in a window
centered on the current state.  If the proposed state gives a better
likelihood then the proposal is accepted (or based on the Metropolis
Hastings algorithm, if the likelihood is only a little worse it will
also be accepted sometimes).  Now if the proposal window is narrow, then
proposals will be accepted a lot because each proposal is near the
current state.  However, the values of the parameter
explored/proposed/sampled by the Mcmc will be in a narrow range, which
does not promote good 'mixing'.  If however the window is wide, then
proposed values will usually be so far away that they will not be
accepted; in this case the parameter stays in its current state, which
is also not good for mixing.  So proposals need to be 'tuned' so that
proposal acceptance rates are not too big and not too small.  This is
done by setting tunings of the Mcmc object, as for example in::

     m.tunings.relRate = 0.5          # Not part-specific
     m.tunings.parts[0].comp = 0.2    # Part-specific

Tunings have default values, but they may be inappropriate for your
analysis, and will need adjustment.  I can do this by trial and error,
on short trial runs.  (See below concerning autoTune(), which does this
for you.)  To skip burnin, I start with the last tree and model from a
previous trial.  It sometimes takes me a few trials to get the
acceptance probabilities within a good range.  (What is a good range?
-- Neither too big nor too small.  The authors of the MrBayes program
suggest that good acceptance rates fall in the range 10-70 per cent;
that is probably ok.)

Tunings for the same proposal can be different in different data
partitions (parts).  Tunings for comp, rMatrix, gdasrv, and pInvar are
part-specific.

Tunings for the proposals in p4 are like window sizes-- if you increase
the number, the window gets bigger, and the acceptance probability goes
down.  Tunings for gdasrv, pInvar, relRate, comp, rMatrix, and local are
like that.  There are no tunings for re-rooting.  The tuning for for
moving model things around the tree is a bit of a hack -- it is the log
prior ratio, and by default is zero.

The exchange between pairs of chains in the MCMCMC is affected by the
chainTemp tuning.  Bigger numbers means less acceptance of proposed
exchanges.

You can use the Mcmc.autoTune() method to tune automatically::

     m = Mcmc(t, ...)
     m.autoTune()
     m.run(...)


Assessing the MCMC run
----------------------

After an MCMC run, you can use various diagnostics to assess whether
the run was good; or rather to assess whether the run was not bad (not
quite the same thing).  Probably the best way to do that would be to
query checkpoints, and that is described in the next section.  Here I
describe some things that you can ask an Mcmc object.

You will probably want to know about how often proposals were accepted.
You can demand::

     m.writeProposalAcceptances()

That gives a table showing how many proposals were made, and the
acceptance rate.  The acceptance rate should not be too big or too
small, as explained above (see `Tuning`).  If the acceptance rates are
bad, you can change the tuning for that proposal.  Recall that to see
un-normalized proposal probs, and tunings, you can ask to::

     print m.probs
     print m.tunings

Another table, given by::

     m.writeProposalIntendedProbs()

shows how often different proposals are set to be proposed.  These
rates add up to 1, and are affected by m.probs.  This table is printed
at the beginning of a run.

When using an MCMCMC, you can see how often the cold and heated chains
exchanged with one another by::

     m.writeSwapMatrix()

If there are too many or too few swaps, you can adjust the temperature
``note Tuning``


Checkpoints
-----------

When running an MCMC, you can write checkpoint files from time to time.
These files are the state of the Mcmc at that time.  Mcmc runs can be
restarted using checkpoints.  Also, you can do diagnostics on
checkpoints like you can on live Mcmc objects, but with checkpoints you
can do the diagnostics after the run has finished (and the Mcmc object
no longer exists) or during the run (querying finished checkpoints), and
you can do the diagnostics on several checkpoints to see how the
diagnostics change over time.

To tell your Mcmc to make checkpoints, say (among other args when you
start up an Mcmc) for example::

     m = Mcmc(t, ..., checkPointInterval=250000, ...)

To tell it not to do checkpoints, set it to zero or None::

     m = Mcmc(t, ..., checkPointInterval=None, ...)

The checkpoint interval should divide the number of generations
requested by run() evenly, so for example::

     m = Mcmc(t, ..., checkPointInterval=250000, ...)
     m.run(1000000)

will work, but ::

     m = Mcmc(t, ..., checkPointInterval=250000, ...)
     m.run(900000)

will not work.

I generally aim to collect perhaps 4 or 5 checkpoints in my runs.  If
you collect more checkpoints (by collecting them more often) then they
will each contain fewer samples, and so the estimates might be a bit
more noisy than if you take checkpoints at bigger intervals representing
more samples.

There is a class, McmcCheckPointReader(), that is good for reading and
digesting checkpoints.  When you start it up like this::

     cpr = McmcCheckPointReader()

then it will slurp in all the checkpoint files in the current directory.
There are other ways to start it up - read the class doc string.
Having got a McmcCheckPointReader object, then you can ask it to ::

     cpr.writeProposalAcceptances()

or ::

     cpr.writeSwapMatrices()

which calls those methods on all the checkpoints.  When you read in
checkpoints, they are full Mcmc objects (except that they do not have
data attached), and so you can ask questions of them as you would ask of
an Mcmc object, as::

     m = cpr.mm[0]    # get the first Mcmc object
     m.writeProposalAcceptances()

Using the McmcCheckPointReader is useful for seeing how things, for
example acceptance rates, change over the run.  Perhaps the most useful
thing that a McmcCheckPointReader can do is compare the split supports
between runs, using the average standard deviation of split support (or
split frequency).  When you do an Mcmc analysis it is good practice to
do more than one run, and comparing split supports between runs
measures of topological agreement between runs.

You can compare split supports between two different checkpoints, or
between all pairs, as::

     cpr.compareSplits(0,1)    # Compare checkpoint 0 with 1
     cpr.compareSplitsAll()    # Compare all checkpoints to each other


Restarting from checkpoints
---------------------------

You can restart an Mcmc from a checkpoint.  If the previous run
finished normally (ie was not a crash) then it is easy.  A checkpoint
is a full Mcmc object without the data - so you need to give it the
data to get going.  You can use the function::

     func.unPickleMcmc(runNum, theData, verbose=True)

to get the Mcmc object from the checkpoint.  That function will get the
last checkpoint from the specified run (runNum) and return an Mcmc
object.  So for example you might restart runNum 0 by::

     read("../d.nex")
     d = Data()
     m = func.unPickleMcmc(0, d)
     m.run(1000000)

which tells it to continue on in the same way as it was, and do another
million generations.  It would be possible to make changes to the Mcmc
before continuing the run, as::

     m = func.unPickleMcmc(0, d)
     < make changes here ...>
     m.run(1000000)

If the Mcmc crashed, you can restart from a checkpoint as above, but
first you will want to repair the output files to get rid of output
lines that were written after the last checkpoint but before the crash.


Output
------

Trees sampled in the MCMC are placed in a file, and the log likelihoods
of those sample trees are found in another file.  Model parameters are
put in another file.

You can make a consensus tree like this::

     tp = TreePartitions('trees.nex', skip=1000)  # skip burnin
     t = tp.consensus()

You will often want to transfer the node support values to internal node
names -- see the doc string for TreePartitions.

You can do a quick-and-dirty convergence test by plotting the log
likelihood values found in the output file.  If it reaches a plateau,
then it is assumed to have converged.  However, this method is
unreliable.  Certainly if the log likes have not reached a plateau then
it has not converged, but the reverse cannot be assumed.  Comparing
split supports as described above (see Checkpoints) offers better
convergence diagnostics.

You can also look at sampled model parameters, found in the mcmc_prams_N
(N = 0, 1, ...) using::

     func.summarizeMcmcPrams()

You can set it to skip a burnin.

.. _post-pred-sims-label:

Posterior predictive simulations
--------------------------------

To help assess model fit in a MCMC you can set up the MCMC to simulate
data sets based on the current tree and model every writeInterval.  When
a simulation is made, a test quantity is extracted from it and written
to a file.  

You can turn on simulations during an MCMC when you instantiate the
Mcmc object.  You set the *simulate* arg to a number from 1-31.  For
example, if you set it to 1, you get the unconstrained or multinomial
likelihood, if 2 you get X^2, and if 3 you get both.  The available
test quantities are in the class doc string for Mcmc here :class:`Mcmc.Mcmc`.

If you want to do the simulations after the MCMC is finished (or
perhaps from a MrBayes run), see :class:`PosteriorSamples.PosteriorSamples`

The idea is that you have a single test quantity from your data (or
data partition) and you want to compare it to the range that is
generated from the posterior distribution, to see whether the model
that you are using might have plausibly generated your original data.
A usual way to do that comparison is to use
:func:`func.tailAreaProbability` (or the same within the Numbers class
:meth:`Numbers.Numbers.tailAreaProbability`).  Here is an example from the
source code:: 

    # Get the test quantity, X^2, from the original data.
    read("../../K_thermus/noTRuberNoGapsNoAmbiguities.nex")
    d = Data()
    ret = d.compoChiSquaredTest()
    #print ret
    originalStat = ret[0][0]

    # Get the sim stats
    n = Numbers('mcmc_sims_0', col=1, skip=500)

    # Evaluate the tail area probability
    n.tailAreaProbability(originalStat)

From which the output is something like::

    Part 0: Chi-square = 47.836914, (dof=12) P = 0.000003
    # The stat is 47.8369140221
    # The distribution has 500 items
    # The distribution goes from 0.598552 to 10.158933
    # Items in the distribution were >= theStat 0 times.
    # The tail-area probability is 0.000000

In this example, the model does not fit, and could not have plausibly
generated the data from which the original test quantity was 47.8.
