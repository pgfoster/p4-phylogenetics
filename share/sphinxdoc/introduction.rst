Introduction
============


What's it good for?
-------------------

**Phylogenetics toolkit**

* P4 can be used as a phylogenetic toolkit, the elements of which you can string together in different ways depending on the job at hand. It is useful for programmatic manipulation of phylogenetic data and trees. If you want to do something interesting with your trees or data, p4 might have at least some of what you want to do already in place.

* P4 will read data in a few of the common phylogenetic formats (eg
  Nexus, Phylip, clustalw, fasta, pir/nbrf), but does not read other formats in bioinformatics (eg EMBL, genbank). P4 will read in trees in Nexus or Phylip format.

* P4 will do some elementary data manipulation, eg extracting a Nexus-defined charset from an alignment, or converting data from one format to another. P4 will also do tree manipulation, and tree drawing.

* It has a big tree viewer, to be able to view big trees (eg up to 5000 taxa) on the screen.

* P4 is meant to be easily extensible, so if you want to do something that it cannot do, it is often easy to add that functionality.

**Heterogeneous models**

Heterogeneity in the process of evolution should be reflected in
heterogeneity in the models that we use. Different genes behave
differently in evolution, and it so we can have data-heterogeneous
models, where things like the rate matrix, composition, and even data
type can differ in different data partitions. However, the process of
evolution can also differ over the tree, and tree-heterogeneous models
can be used to model that. P4 allows maximum likelihood determination
of individual trees under flexible combinations of these
models. (However, p4 does not do tree searching in ML).


**Bayesian analysis**

P4 implements a MCMC for doing Bayesian phylogenetic analysis. It uses
NNI (actually Larget and Simon's LOCAL) and eTBR for tree
rearrangement. It can use the data- and tree-heterogeneous models
described above. It has a full implementation of Paul Lewis et al's
polytomy proposal (P. O. Lewis, M. T. Holder
and K. E. Holsinger, 2005. Polytomies and Bayesian phylogenetic
inference, Syst Biol 54(2): 241--253). P4 writes checkpoint files
during the MCMC run, so that if you want to continue a run after it
has ended, or if you want to recover from a crash or power outage, it
is possible to restart from a checkpoint. You can specify topological
constraints in the trees.

**A citation**

For the moment, until there is a better one, the best citation is 

* Foster, P.G. 2004.  Modeling compositional heterogeneity.  Syst. Biol. 53: 485-495.


