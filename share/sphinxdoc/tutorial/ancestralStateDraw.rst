========================================
Ancestral character state reconstruction
========================================

I think that ancestral state reconstruction should be done in a probabilistic way. I can’t see that it would be good enough to simply ask for the most probable ancestral state; it would be better to make draws from the posterior probability, and that is the way that ``p4`` does it.

Ziheng Yang has a very good description of ancestral state reconstruction in his "Molecular Evolution — A statistical approach" (2014)

The user interface method is :meth:`p4.tree.Tree.ancestralStateDraw`, which needs a tree with attached model and data.  It returns a string, which is a single draw from the inferred ancestral states.

There are examples on how to use it included in the ``share/Examples`` directory.  The first example uses an ML optimized model, and the second example makes draws from the posterior distribution of a Bayesian analysis.


 
