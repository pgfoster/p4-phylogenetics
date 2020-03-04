===============================
Simulations with reference data
===============================


.. _simsWithRefData:

- 2016-05-03

The ususal way to simulate data based on a tree and a model
-----------------------------------------------------------

The way that p4 has been doing simulations, and I belive the way most other simulation programs do it, is to generate a random sequence for the root and evolve that sequence over the given tree using the given model.  The root sequence is random, but has the model composition.  For ``+G+I`` models, the discrete gamma category and whether it is an invariant site is also chosen for each site in the root sequence, and that choice is maintained throughout the tree.  

Here is an example simulation using p4.  It uses a GTR+G model; the pInvar model is not used here because some programs used below (PAML/evolver and PhyloBayes) do not use it.

.. code-block:: python

    from p4 import *
    nTax=5
    taxNames = list(string.uppercase[:nTax])
    a = func.newEmptyAlignment(dataType='dna', taxNames=taxNames, length=80)
    d = Data([a])
    read('(B:0.5, ((D:0.4, A:0.3):0.1, C:0.5):0.1, E:0.5);')
    t = var.trees[0]
    t.taxNames = taxNames

    t.data = d
    t.newComp(free=0, spec='specified', val=[0.1, 0.2, 0.3])
    t.newRMatrix(free=0, spec='specified', val=[2., 3., 4., 5., 6., 7.])
    t.setNGammaCat(nGammaCat=4)
    t.newGdasrv(free=0, val=0.5)
    t.setPInvar(free=0, val=0.0)

    t.simulate()

    d.writeNexus('d.nex', writeDataBlock=True)
    d.alignments[0].writePhylip("d.phy", interleave=False, flat=True)

    # phylobayes will use this tree, and it does not like spaces in the tree file,
    # so use spaceAfterComma=False
    t.writeNewick(fName="simTree.phy", spaceAfterComma=False)

and here is the resulting file, ``d.phy`` ---

::

     5  80
    A          tcggctggtcatttggcgcccttgtttcggcggttgtgtatctggttactcgtctgttttcccttagctgggtctttgtt
    B          tcgagtgagatttttgcttccttgtttctgcgggttcgtccgggggcactcggcagatctccgttcgttgcatctttttt
    C          ttgcgtgatatttttcgtgcgttgttccctggtctctgtacctgagtactcggcagcgtgcttttgggtggttcgttatt
    D          tcgtttgagcttttttcagcctggttccggcgtttttgtagctgactactgggcagtttgccttaggtactttcgttatg
    E          tcggctgatatttgtcgtgccttgttcctgcgggtttgtgtgtgtcaatgcgcccgctttcggtttgtcggtggttgttt

Seq-Gen
~~~~~~~

A commonly used simulation program is the venerable Seq-Gen, written by Andrew Rambaut, and
available at  `http://tree.bio.ed.ac.uk/software/seqgen <http://tree.bio.ed.ac.uk/software/seqgen>`_.
It is described in the paper available at `http://bioinformatics.oxfordjournals.org/content/13/3/235.short <http://bioinformatics.oxfordjournals.org/content/13/3/235.short>`_.
In this example I am using ``v1.3.3``, with the command

::

    seq-gen -mGTR -r2.,3.,4.,5.,6.,7. -f0.1,0.2,0.3,0.4 -l80 -n1 -a0.5 -g4 -i0.0 simTree.phy > seq-gen_out.phy 

Here is the resulting simulation --- 

::

     5 80
    B         TTTCTGAGGCGCTCGTTCGGGCCCGTATTCTGTTGCTCGTGCTTTTTTTTCTCCGTAACGTTGGGGTTCGTCCGATTTCT
    D         CTGCGGCGGCGCTCTTTAGTGTCCGGCGTGTTTTTCTCGCTCGTCTTTTGTTCCGTAGCTTCGGTGTTTGTTCTAGTTGT
    A         TTTCGGTGGCGCTCCTATGTGCCCTTCGTTTTTTGCTAGCTCTTGTGTCTATCCGTGTTTTTGGTCTTAGTGCTATGTCT
    C         TTTTCGGTGCGCTCCTTTGGTCCCGTATTGGGTTTCGTGTTCTTTTGTTGCTCCGTGATTTTGGTGTTTGTGCTACATCT
    E         GTTACGTGACTCTCCTTGGCTTCCGCGTTGTTTGGCTGCATCTTGTGTCTCCCCGTGTGTTGAGTTTTTGTCGTAGTTCT

(With Seq-Gen you can specify an ancestral sequence, but I do not see an option where you can, for each site, set the gamma category or whether it is an invariant site.)

INDELible
~~~~~~~~~

INDELible is a very capable simulation program written by William Fletcher when he was in Ziheng Yang's group, and described in this paper `http://mbe.oxfordjournals.org/content/26/8/1879.full <http://mbe.oxfordjournals.org/content/26/8/1879.full>`_.
The manual, tutorial, and source code can be found at `http://abacus.gene.ucl.ac.uk/software/indelible/ <http://abacus.gene.ucl.ac.uk/software/indelible/>`_.

As in PAML, the composition (state frequencies) is given in the order ``TCAG``, and the rate matrix parameters are ``TC TA TG CA CG`` with  ``AG=1``.
So my composition of ``ACGT = [0.1, 0.2, 0.3, 0.4]`` would be ``[0.4, 0.2, 0.1, 0.3]``, and 
my RMatrix ``AC AG AT CG CT GT [2., 3., 4., 5., 6., 7.]`` would be reordered as ``[6., 4., 7., 2., 5., 3.]``, and then scaled as ``2.000 1.333 2.333 0.667 1.667 1.000``.

I used this control file to run INDELible.

::

    [TYPE] NUCLEOTIDE 1
    [MODEL]    myGTR
      [submodel]  GTR 2.000  1.333  2.333  0.667  1.667
      [statefreq] 0.4 0.2 0.1 0.3
      [rates]     0.0 0.5 4
    [TREE] myTree  (B:0.5,((D:0.4,A:0.3):0.1,C:0.5):0.1,E:0.5);
    [PARTITIONS] myPartition
      [myTree myGTR 80]
    [EVOLVE] myPartition 1 indelible_out


The simulation resulted in this alignment ---

::

    5  80
    B     CCTTTATGCGCTGGTGTGGTTGTAGGATTGTCCGTAGCTTGCTGGGACTATTTTGGTGGGTAGCAGGTGGGTGGGCAGCC     
    D     CCTTTATGCTCTGCTGTGGTTGTAGGTTTTTTTGTATCTTGCTCGGACTATGTTGTGGGGTTCCAGGTGGGAGTTCGGCC     
    A     CCTTTATGCGCCGCTGTGGTTGTAGGTTTGTTTGTATGTCGCTGGGACGATGTTTTGGGGGGACATTTCGGTGTTCAGCC     
    C     CGTTCTTGCGCCGGTGTGGTTGTAGGTTTTGTTGTATCTTGCTGCTACTATTTCGGGGGGTGTCATTTTGGTGATCAGCC     
    E     CCTTGTTGCTCGGGTGTGGTTGTAGGTTGGTTTGTAGCTGGCTGGGACCATTTCGCGGGGGGGCAGGTCGTTGTTCAGCC     

(In the 2009 Fletcher and Yang paper referred to above, the authors say that before INDELible  "... only MySSP (Rosenberg 2005) can simulate under nonstationary and nonhomogenous models."  It appears that they were not aware of my description of my ``p4`` software to do non-stationary, non-homogeneous simulations published in my 2004 paper `http://sysbio.oxfordjournals.org/content/53/3/485.full <http://sysbio.oxfordjournals.org/content/53/3/485.full>`_.)

Evolver from PAML
~~~~~~~~~~~~~~~~~

`PAML <http://abacus.gene.ucl.ac.uk/software/paml.html>`_ is a suite of programs, including ``evolver``, which does a few things including simulating sequences on a tree and model.
I used PAML v 4.9a.

In PAML, DNA bases are in the order TCAG, and presumably the rate matrix parameters are the same as in INDELible (I'm not sure, and it is not obvious from the documentation).  I will use the following control file.  I don't see a way to use the pInvar model, but I am not surprised as I know Ziheng is not keen on it.

::

     0     * 0,1:seqs or patterns in paml format (mc.paml); 2:paup format (mc.nex); 3: paup JC69 format
     -1234567   * random number seed (odd number)

    5 80 1  * <# seqs>  <# nucleotide sites>  <# replicates>
    -1         * <tree length, use -1 if tree below has absolute branch lengths>

    (B:0.5, ((D:0.4, A:0.3):0.1, C:0.5):0.1, E:0.5);

    7          * model: 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV
    2.000  1.333  2.333  0.667  1.667   * kappa or rate parameters in model
    0.5  4     * <alpha>  <#categories for discrete gamma>

    0.4 0.2 0.1 0.3    * base frequencies
      T        C        A        G

And here is the resulting simulated data ---

::

    5 80 

    B           TACCGCTGTT TTTTTGTGTG TAGTGTGTTT GCTGGCGGTT GACTCTCAAT TTCGAGGGAT GTCGTGGCTC GGCGTTCCTT 
    D           GACTGCGGTT TGTTTGTGTT TATGGGTTCT GGTCGCAGTT GTCGCTGAAA GGCGATGCTT GGCGTCTCGT CTTTCTCCCT 
    A           GACTCCGGTT TCTTTGTTTT TAGGGCTGGG GCTCACTGTG GACACTGAAC GTCGTTGGTT TGTGTGCCGT TTCACTCCGT 
    C           TTCTGCGTTT TTTCTGTGTG TAGGCCTTTT GTTGCTTGTC GACTGTTAAA TTCGGTGTTT TGGGTATCGT TTCACGCCGC 
    E           GACAGCGGTT GTTGTGTGTT TAGGGGTTGT GTTGTAGGTC GACGCTCAAT CTCGTCGGGC GTTGGTGCTT ATCGCTCCGT 

Ziheng Yang comments in the manual ---

"Some people wanted to specify the sequence at the root rather than letting the program generate a random sequence. This can be achieved by putting a sequence in the file RootSeq.txt. The sequence cannot have ambiguities or gaps or stop codons. In almost all simulations, it is simply wrong to fix the root sequence, so you should resist the temptation of making the mistake. If you want the simulation to reflect your particular gene, you may estimate parameters under a model from that gene and then simulate data sets using the parameter estimates."


... and in the PAML google groups, Ziheng comments in answer to a question about ``evolver`` ---

"I also seem to remember writing some warning notes against using fixed sequences at the root, but can't remember where it was.  from my experience, most users' justification of using fixed root sequence (like making my data look more similar to my real observed data) is not sensible, and you should be wary of the problems."

PhyloBayes has an unusual way to simulate data --- using a tree, model, and reference data
------------------------------------------------------------------------------------------

`PhyloBayes <http://megasun.bch.umontreal.ca/People/lartillot/www/index.htm>`_ is an extraordinary Bayesian phylogenetics program written by Nicolas Lartillot, and first described in the paper `http://mbe.oxfordjournals.org/content/21/6/1095.short <http://mbe.oxfordjournals.org/content/21/6/1095.short>`_.  It uses posterior predictive simulations to assess fit of the model to the data.  The simulations in PhyloBayes are unusual.  To show this, I will use PhyloBayes version 4.1c (not, in this case, the mpi version), and the data generated above by ``p4`` to run a short MCMC using the GTR+G model using a fixed topology --- the simulation tree above.  After that I tell the ancillary program ``ppred`` to simulate three sets of data from the last three posterior samples.

I think the following runs a GTR+G model ---

.. code-block:: sh

    pb -s -f -gtr -ncat 1 -dgam 4 -d d.phy -T simTree.phy -x 5 103 r1 


After the MCMC I used ``ppred`` to do the simulations, making three simulated alignments, here ---

::

    5       80
    A     GCTACCGATTGTTGTCGCTCCTTTTTGCCCGGTTTATGTGGGTGGCCAACCTCGCGAGCCCGAGGTGGCCTCGGCTTACT
    B     GCTTCTGACTCTGTATTGTGCTTATTTCTGGCCTGCTGTGCGTGTTTAACCTCGCGATCGCTACCTGCCGTTGGTTTAGA
    C     TCTACTGATTATTCCCGGGCCTTTTTCCCGCGGTTTTGTGCGTGCGGAACCTCGTGTTCTCGATCTGGCCTTGGTTTTTT
    D     ACTTCTCATAATTGCTCGTCCTTGTCTCTCGGCCTCTGTGAGTGCTTAATTCCGGGACCGCGATGTGCGGGTTGTGTGGT
    E     TCGGCTGATATTTGTCGTGCCTTGTTCCTGCGGGTTTGTGTGTGTCAATGCGCCCGCTTTCGGTTTGTCGGTGGTTGTTT
    5       80
    A     TTGGTTGTTATTTATAGTTTTTTCTTTCCGCGGCTTTGTGGGGGGTTATCCCGCCTTTTTCTTTTTGGAGGTGTTTGTGG
    B     TTGTCTGAGATTTATTGTTCTTTTTTTCTGCGGGTGTGTGTGCGATTATTCGCCTGGTTTCTTTTGGGTGTTCGGTTCGG
    C     TGGTTTGATATTTTTGGTTTCTTGTTCCTGCGGCTACGTGTGGGTCTATCGCTCCTGTTTCAGTGTGTCAGGGTTTTTGG
    D     TTGGTTGAATTTTATAGTTTCTTGTTGCAGTGTTTTTGTCTGGGTTTATCCCGCCTCTTTCATTTTGTGGCCGTTTTTGT
    E     TCGGCTGATATTTGTCGTGCCTTGTTCCTGCGGGTTTGTGTGTGTCAATGCGCCCGCTTTCGGTTTGTCGGTGGTTGTTT
    5       80
    A     TGGGTTGATTTTTCTCTCGCTTGTTTCCTTGGCGCTTGTGTGTGGGGATGGTTCTGCTGACGTTGGTTTGTCTGGTGGCT
    B     TCGCCTGACTTTTCTGGTGCCTCTTTCCTCGGCGGTTTTTGGTGTGGATGGTTCCGACGTCGTTTTTTCGGGGGTTGTCT
    C     TCGGTTGATTTTTTTCGTGTTTACTTCCGCGGCCGTTTTTCGGGGGTATTGTGCTGATGTCAATTTTTCGGTCGGTGTCT
    D     TGGCGTGATTTTTCTCGTGCCTATGTCCCCGGCGTATGTACGTGTGGATGGTTCGGCTGCCGTTGTTTCGTTGGCTGTCC
    E     TCGGCTGATATTTGTCGTGCCTTGTTCCTGCGGGTTTGTGTGTGTCAATGCGCCCGCTTTCGGTTTGTCGGTGGTTGTTT

For reference, here I repeat ``d.phy``, made above by ``p4``, and used as the "original data" in the PhyloBayes analysis.

::

     5  80
    A          tcggctggtcatttggcgcccttgtttcggcggttgtgtatctggttactcgtctgttttcccttagctgggtctttgtt
    B          tcgagtgagatttttgcttccttgtttctgcgggttcgtccgggggcactcggcagatctccgttcgttgcatctttttt
    C          ttgcgtgatatttttcgtgcgttgttccctggtctctgtacctgagtactcggcagcgtgcttttgggtggttcgttatt
    D          tcgtttgagcttttttcagcctggttccggcgtttttgtagctgactactgggcagtttgccttaggtactttcgttatg
    E          tcggctgatatttgtcgtgccttgttcctgcgggtttgtgtgtgtcaatgcgcccgctttcggtttgtcggtggttgttt

Notice that the PhyloBayes simulations are more-or-less similar to the original data ``d.phy``.  Although this example does not show it, PhyloBayes will even match the positions of original alignment gaps in the simulated data.

How similar are the simulations to the original data?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simulations were made from Seq-Gen as above, and from PhyloBayes, collecting 100 simulations for each.  For each, the number of differences between the original datset and the simulation was counted up, sequence by sequence, and position by position.  Since there are five sequences each 80 characters long the maximum difference is 400.  (And yes, that the Seq-Gen sequences rearranged the order of the taxa was taken into account in doing this measurement ...!).    :numref:`fig-simsDiffsA-label` shows that most of the sites differed in the Seq-Gen simulations, but most of the sites were the same in the PhyloBayes simulations.

.. _fig-simsDiffsA-label:

.. figure:: SimsWithRefDataFigs/simsDiffsA.svg

    Differences between the original data and the simulated data.  Black bars are from Seq-Gen, and white bars are from PhyloBayes.

An implementation in p4
-----------------------

I think the way PhyloBayes is doing the simulations is to use a posterior sample of the tree and model, and then use that together with the original data to make draws from probabilistic estimates of the root character states.  This simulated root is then evolved as usual on the tree using the model, to make the simulated data.
The root sequence simulations are as usual based on the posterior sample of the model parameters and branch lengths of the tree, but additionally and unusually the simulations are also based on the orginal data.  I implemented such a strategy, inspired by PhyloBayes, in p4.

The root sequence simulation is based on the *conditional likelihoods*  at the root, themselves dependent on the model prameters and sampled topology.  Conditional likelihoods are used in likelihood calculations and so are available if a likelihood based on the posterior sample has been calculated using the original data.  The root sequence simulation is tantamount to a *sampled probabilistic ancestral state reconstruction*.  The root state is drawn from the character states, and if the +G model is used a draw is made from the gamma category, and if the pInvar model is used and the original data site is constant a probabilistic choice is made about whether the site is invariant.  If we only look at (not-CAT, eg GTR-like) models that do not have among-site rate variation, then the probability :math:`P` of the root state being :math:`j` given leaf data :math:`X`, tree :math:`T`, and model parameters :math:`\theta` with character state frequencies :math:`\pi_i`, is 



.. math::

    P(j|X,T,\theta) = \frac{\pi_j L(j)}{\sum_i \pi_i L(i)}

where :math:`L(j)` is the conditional likelihood at the root for character state :math:`j`, and :math:`\sum_i \pi_i L(i)` is the *site likelihood*.  Among-site rate variation (gamma and pInvar) are similar.  This way to draw ancestral states would apply to any model in ``p4``, including the tree-heterogeneous models. 


I ran an MCMC using p4 with the same data ``d.phy`` as above, with the GTR+G model.  After, I used the   :class:`~p4.posteriorsamples.PosteriorSamples` class to get the samples (tree+model) from which I generated simulations.  I did simulations both with and without the refTree+model+refData.  Differences are plotted below.  The refTree simulations appear to be similar to the PhyloBayes simulations above, so perhaps I implemented it correctly.  The :math:`X^2` values from the two simulation sets, with and without the refData, were similar (:numref:`fig-simsX2A-label`).

.. code-block:: python

    def aligDifference(a, b):
        diffs = 0
        for i in range(a.nTax):
            sA = a.sequences[i]
            sB = b.sequences[i]
            for j in range(a.nChar):
                if sA.sequence[j] != sB.sequence[j]:
                    diffs += 1
        return diffs

    read('d.phy')

    # The Data d will be attached to the sim tree, and so the data contents will
    # change with each simulation.  The refData stays the same.
    d = Data()
    refData = d.dupe()
    t = func.randomTree(taxNames=d.taxNames)
    t.data = d

    pNum = 0
    t.newComp(partNum=pNum, free=1, spec='empirical')
    t.newRMatrix(partNum=pNum, free=1, spec='ones')
    t.setNGammaCat(partNum=pNum, nGammaCat=4)
    t.newGdasrv(partNum=pNum, free=1, val=0.5)
    t.setPInvar(partNum=pNum, free=0, val=0.0)

    # Check to make sure its all good to go.
    t.calcLogLike(verbose=False)

    # Instantiate
    ps = PosteriorSamples(t, runNum=0, program='p4', verbose=0)

    myDiffs = []
    myDiffsWithRefTree = []
    bigXSq = []
    bigXSqRT = []

    for sampNum in range(100,200):
        t2 = ps.getSample(sampNum)
        t2.data = d
        t2.simulate(refTree=None)
        bigXSq.append(t2.data.simpleBigXSquared()[0])
        myDiffs.append(aligDifference(refData.alignments[0], t2.data.alignments[0]))


        # Now do the sims with a refTree+model+data
        refTree = t2.dupe()
        refTree.data = refData
        t2.simulate(refTree=refTree)
        bigXSqRT.append(t2.data.simpleBigXSquared()[0])
        myDiffsWithRefTree.append(aligDifference(refData.alignments[0], t2.data.alignments[0]))


.. _fig-simsDiffsB-label:

.. figure:: SimsWithRefDataFigs/simsDiffsB.svg

    Differences between the original data and the simulated data, using p4.  White bars use a refTree+model+refData, and black bars are without.


.. _fig-simsX2A-label:

.. figure:: SimsWithRefDataFigs/simsX2A.svg

    X\ :sup:`2`\ values from the simulations using refTree+model+refData (white bars), and without (black bars).
