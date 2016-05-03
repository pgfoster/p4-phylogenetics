To combine tree supports from two analyses onto one tree
--------------------------------------------------------

You will often analyse the same data with two different
methods, each of which gives you a tree with support values, but you
may want to summarize the results of both trees using a single tree
figure.  The trees may differ in topology, or they may be the same.
Even if they are the same topology, the branch lengths may differ.
One approach is to consider one tree to be the "master" tree, giving
the figure its topology and branch lengths, but you then somehow
combine the split supports from the other tree onto the master.  

Consider these two trees::

    #          +-------------- E
    #   +------|0.84
    #   |      +----------- B
    #   |
    #   |-------- G
    #   |
    #   |                      +-------- C
    #   |                      |
    #   |               +------|0.90          +-------------- D
    #   |               |      |      +-------|0.95
    #   |               |      |      |       +----------- F
    #   |               |      +------|0.87
    #   |       +-------|0.95         |       +------------ H
    #   |       |       |             +-------|0.88
    #   |       |       |                     +-------------------- A
    #   +-------|0.59   |
    #           |       +-------------- I
    #           |
    #           +--------------- J
    #   
    #           +--------------- D
    #   +-------|73
    #   |       +------------- F
    #   |
    #   |       +-------------- H
    #   |-------|77
    #   |       +---------------------- A
    #   |
    #   |               +-------- C
    #   |       +-------|86
    #   |       |       +---------------- I
    #   |       |
    #   +-------|88      +---------------- J
    #           |        |
    #           +--------|58             +--------------- E
    #                    |       +-------|98
    #                    +-------|99     +------------ B
    #                            |
    #                            +-------- G

 
These two trees are both fully resolved, and mostly have the same
splits -- there is only one NNI difference (can you see that?).  If we make the
first tree the master tree, it is a fairly easy job to merge the split
supports from the second tree onto that first tree, as shown below.

Doing the combining "by hand" is an option, but for bigger trees
(which may share only some of the same splits) the task becomes
onerous, time-consuming, and error-prone.  The process can be
automated, and one way to do that is shown here.  The
process involves assigning numbers to each split in both trees;
those numbers being a numerical equivalent of the "dot-star"
notation for splits --- for example ``.**.**`` might be binary ``011011``.
Since the taxa are the same in both trees, and in the same order,
those numbers can be used by each split in the master tree to find
the same split (if it exists) in the other tree.  Then the supports
can be combined::


    #          +-------------- E
    #   +------|0.84/98
    #   |      +----------- B
    #   |
    #   |-------- G
    #   |
    #   |                      +-------- C
    #   |                      |
    #   |               +------|0.90/-        +-------------- D
    #   |               |      |      +-------|0.95/73
    #   |               |      |      |       +----------- F
    #   |               |      +------|0.87/88
    #   |       +-------|0.95/58      |       +------------ H
    #   |       |       |             +-------|0.88/77
    #   |       |       |                     +-------------------- A
    #   +-------|0.59/99|
    #           |       +-------------- I
    #           |
    #           +--------------- J

To start, we need a valid list of taxnames, which we might get from an alignment::

    a = func.readAndPop('myData.nex')
    # Now a.taxNames has the list we want

We read in and name our two trees, and make sure they both have the same taxNames::

    tMaster = func.readAndPop('myMainTree.nex')
    tSecondary = func.readAndPop('mySecondaryTree.nex')

    tMaster.taxNames = a.taxNames
    tSecondary.taxNames = a.taxNames

At this point we are set up and good to go::
    
    # Split keys are numerical versions of the 'dot-star' split notation.
    # The same split on the two trees would have the same split key.
    tMaster.makeSplitKeys()
    tSecondary.makeSplitKeys()
    
    # Make a dictionary, so that we can fish out nodes in the secondary tree
    # given a split key.  Split keys are found on node branches, here
    # n.br.
    myDict = {}
    for n in tSecondary.iterInternalsNoRoot():
        myDict[n.br.splitKey] = n
    
    for nM in tMaster.iterInternalsNoRoot():
        # Given a split key in the master tree, we can find the
        # corresponding node in the secondary tree, using the split key with
        # the dictionary.
        nS = myDict.get(nM.br.splitKey)
        # If there was none, then nS is None
        if nS:
            nM.name = '%s/%s' % (nM.name, nS.name)
        else:
            nM.name = '%s/-' % nM.name
        #print nM.name
    tMaster.writeNexus('combinedSupportsTree.nex')


There is an example of this in ``share/Examples/F_picture/E_combiningSplitSupports``.
