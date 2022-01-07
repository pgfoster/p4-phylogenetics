=====================
Models and likelihood
=====================

In order to calculate the likelihood of a tree, you need to specify a
model.  Since the model depends on the datatype and data partition
structure, before you can specify a model you need to specify the data.

Having specified a model, you can calculate the likelihood without
optimization with :meth:`tree.Tree.calcLogLike`;  You can optimize
any free parameters with the :meth:`p4.tree.Tree.optLogLike` method. ::

     t = var.trees[0]
     t.data = Data()   # specify the data
     <... specify the model ...>
     t.calcLogLike()



Describing models
=================

You need to specify 

* a composition, 
* a rate matrix, and 
* the proportion of invariant sites.  
* Optionally you can specify gamma distributed among-site rate
  variation (gdasrv), with the number of discrete gamma classes, nGammaCat.


Here is a simple one, for an F81 model, with no among-site
rate variation::

    t.newComp(free=1, spec='empirical')
    t.newRMatrix(free=0, spec='ones')
    t.setNGammaCat(nGammaCat=1)       # no 'gamma'
    t.setPInvar(free=0, val=0.0)      # no proportion of invariant sites

(Note that while the F81 is a DNA model, the model above would also
work for other datatypes.)

Parameters of the model may or may not be *free*
(ie adjustable or optimizable), and so that needs to specified.
Depending on the *spec* (ie model specification), some model parameter
numerical values may need to be specified as a *val* argument.  For
example, to define the a composition for the third data partition
(partNum=2) for which you want the model composition to be fixed to
certain values (ie not free), you might say::

     t.newComp(partNum=2, free=0, spec='specified', val=[0.4, 0.3, 0.2])

A simple Poisson or Jukes-Cantor-like model can be described for single
partition data by the following::

     t.newComp(free=0, spec='equal')
     t.newRMatrix(free=0, spec='ones')

Here the *spec* for the composition is ``equal``, meaning that the
character state frequencies are all equal; since they are equal you do
not need to specify the numbers.  In specifying the rate matrix, the
``spec='ones'`` arg means that the rates of change from one base to another
are all the same.  Here the model parameters are not free, and so they
would not change in a likelihood optimization or in a Bayesian MCMC
analysis.

If all the parameters on the model above were set to be free (ie by
setting free=1) then it would become a GTR model.  If the comp was set
free but the rMatrix remained fixed then it would be an F81 model.

K2P and HKY models are specified by setting the rate matrix spec to
'2p'.  Empirical protein models are specified by setting the rate matrix
spec to 'd78', 'jtt', 'wag', or 'mtrev24'.

Multi-partition models for mult-partition data
----------------------------------------------

The model can differ over the data, and if you have more than one data
partition you need to specify the ``partNum`` when you specify the
components of a model::

    pNum = 0    # First partition, F81
    t.newComp(partNum=pNum, free=1, spec='empirical')
    t.newRMatrix(partNum=pNum, free=0, spec='ones')
    t.setNGammaCat(partNum=pNum, nGammaCat=1)
    t.setPInvar(partNum=pNum, free=0, val=0.0)
    t.setRelRate(partNum=pNum, val=0.9)

    pNum = 1   # second partition, F81+G
    t.newComp(partNum=pNum, free=1, spec='empirical')
    t.newRMatrix(partNum=pNum, free=0, spec='ones')
    t.setNGammaCat(partNum=pNum, nGammaCat=4)
    t.newGdasrv(partNum=pNum, free=1, val=0.5)
    t.setPInvar(partNum=pNum, free=0, val=0.0)
    t.setRelRate(partNum=pNum, val=1.1)


If your model is heterogeneous over the data, you can optionally specify
the relative rates of the different partitions, for example::

     t.setRelRate(partNum=0, val=0.5)
     t.setRelRate(partNum=1, val=1.5)

If you want to make these free parameters, you can say::

     t.model.relRatesAreFree = 1

Tree-heterogeneous models
-------------------------

The composition or rate matrix of the model can differ over the tree.
(The gdasrv can also differ over the tree, but this seems to be less
useful.)  

To specify it exactly (which is sometimes useful, but it would not
generally be useful at the start of an MCMC), first you specify the model attribute, and then
you put the model attribute on the tree, for example, if we start with this
tree::

              +--------2:one
     +--------1
     |        +--------3:two
     0
     |--------4:three
     |
     +--------5:four

and put 2 compositions on it, like this::

     A = t.newComp(spec='empirical', symbol='A')
     B = t.newComp(spec='empirical', symbol='B')
     t.setModelComponentOnNode(A, node=t.root, clade=1)
     t.setModelComponentOnNode(B, node=1, clade=1)
     t.draw(model=1)

then we end up with a tree like this::

              +BBBBBBBB2:one
     +BBBBBBBB1
     |        +BBBBBBBB3:two
     0
     |AAAAAAAA4:three
     |
     +AAAAAAAA5:four
     Part 0 comps
         0   A
         1   B
         root (node 0) has comp 0, symbol A

Here I have specified 2 compositions, A and B.  We place A on the root
node, but because we specify clade=1 that composition is applied over
the entire tree.  Then we place composition B on node 1, also
clade-wise, and in that part of the tree B displaces (ie over-rides) A.

An alternative to placing model components on the tree explicitly as
above, you can also :meth:`Tree.Tree.setModelComponentsOnNodesRandomly`.

The multinomial or unconstrained likelihood is a property of the data
only, and does not need a tree or model.  It can only be calculated if
there are no gaps or ambiguities in the data.  There are 2 ways to
calculate it-- you can either take the data partitions into account, or
not.  The former uses the Data method
:meth:`Data.Data.calcUnconstrainedLogLikelihood1`.  The result is placed in the data
attribute ``unconstrainedLogLikelihood``.  The Data method
:meth:`Data.Data.calcUnconstrainedLogLikelihood2` calculates the unconstrained log
like of each data partition.  Note that the unconstrained log like of
the combined data is not the sum of the unconstrained log likes of the
separate partitions.

