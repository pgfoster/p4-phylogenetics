.. _kosiol_ais:

===================================
Kosiol's AIS, almost invariant sets
===================================

This is for grouping amino acids into sets where there is a high probability of exchange within the set but exchange between sets has a lower probability.
The method is described in Kosiol et al (2004) J. Theoret. Biol. 228: 97--106.
The program is available `here <http://www.ebi.ac.uk/goldman/AIS/>`_.

The program wants the equilibrium frequencies for the model, the q-matrix (bigQ), and the eigenvectors.  Here is an example using the wag model, with wag frequencies::

    # The easiest way to get the bigQ and such from p4 is to just set up
    # the appropriate protein model as usual, and then calculate a
    # likelihood.

    # Read in some data.  You could use your own data if you wanted to use
    # empirical comps.
    seqs = """ 2 20
    one       arndcqeghilkmfpstwyv
    two       rndcqeghilkmfpstwyva
    """

    # a tree, with the corresponding taxa.
    ts = "(one, two);"

    read(seqs)
    t = func.readAndPop(ts)
    d = Data()
    t.data = d
    t.newComp(partNum=0, free=0, spec='wag')    # or maybe you want empirical for your own data?
    t.newRMatrix(partNum=0, free=0, spec='wag')
    t.setNGammaCat(partNum=0, nGammaCat=1)
    t.setPInvar(partNum=0, free=0, val=0.0)
    t.calcLogLike()
    #t.model.dump()

    # Write the aa freqs
    f = file('equi', 'w')
    f.write('20\n')
    for i in range(20):
        f.write("%f\n" % t.model.parts[0].comps[0].val[i])
    f.close()

    bigQ = t.model.getBigQ()

    # Write the bigQ
    f = file('q', 'w')
    f.write('20\n')
    for i in range(20):
        for j in range(20):
            f.write("%5g  " % bigQ[i][j])
        f.write('\n')
    f.close()

    # Get the eigensystem
    import numpy
    import numpy.linalg
    evals,evecs = numpy.linalg.eig(bigQ)

    # look, if you want
    if 0:
        numpy.set_printoptions(precision=4, linewidth=300)
        print evals
        print evecs

    # According to the web page, "The right eigenvectors should be ordered
    # according to the absolute value of their eigenvalues."  Well, the
    # output from numpy, which uses lapack, is not so ordered.  So do it.
    sorter = numpy.argsort(evals)
    sorter = sorter[::-1]   # reverse
    #print sorter

    f = file('evec', 'w')
    f.write('20\n')
    f.write('20\n')
    for colNum in sorter:
        for rowNum in range(20):
            f.write("%5g\n" % evecs[rowNum][colNum])
    f.close()

Kosiol's program ``ais`` asks for these 3 files made above, and the
number of groups that you want, and it suggests a grouping.

When I did that with Dayhoff 78, with D78 composition, I got these
groups --- 

- Set 0 = { R N D Q E H K T }
- Set 1 = { C V }
- Set 2 = { A G P S }
- Set 3 = { I L M }
- Set 4 = { W }
- Set 5 = { F Y }

While the groups from the log odds table are 

1. c
2. stpag 
3. ndeq 
4. hrk 
5. milv 
6. fyw 

So both similar and different.  Maybe it differs because the groups from the log odds are
from a PAM250 matrix, which would correspond to highly diverged
sequences?


