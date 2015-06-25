# A couple of Trees methods.  The first one works.

from Glitch import Glitch

def trackSplitsFromTree(self, theTree, windowSize=200, stride=100, fName='trackSplitsOut.py'):
    """See how slits from theTree changes over the trees in self.

    This looks at how some splits change over the trees in self (self
    is a Trees object).  The splits that are tracked are the ones in
    the arg theTree.  If you only wanted to track one split, you could
    supply such a tree with only that one split.

    It uses a moving window, given by the windowSize arg.  The stride
    is the distance between centers of the windows.  So in the
    default, with windowSize=200 and stride=100, the first trees
    looked at will be from 0--199, and the second window will be from
    100 -- 299.  Note that because the stride is less than the
    windowSize, there are overlapping windows.

    The output is to a python file, which you can examine or plot later.
    """
    
    theTree.makeSplitKeys()

    # Decorate the internal nodes of theTree with splitKeys
    for n in theTree.iterInternalsNoRoot():
        if n.name:
            n.name += " (%s)" % n.br.splitKey
        else:
            n.name = "(%s)" % n.br.splitKey
    theTree.draw()
    print "The drawing above shows the splitKeys that are being tracked."

    if fName:
        f = file(fName, 'w')
        textDrawList = theTree.textDrawList()
        for l in textDrawList:
            f.write("#  %s\n" % l)
        f.write("#\n# The drawing above shows the splitKeys that are being tracked.\n")
        f.write("#\n# windowSize=%s, stride=%s.\n" % (windowSize, stride))

    # Remove the decoration from above.
    for n in theTree.iterInternalsNoRoot():
        splitKeyStringWithBlank = " (%s)" % n.br.splitKey
        if n.name.endswith(splitKeyStringWithBlank):
            n.name = n.name[: -(len(splitKeyStringWithBlank))]
        else:
            n.name = None
    #theTree.draw()
    
    
    # Determine whether we need to makeSplitKeys()
    t = self.trees[0]
    if hasattr(t, "splitKeys") and t.splitKeys:
        pass
    else:
        for t in self.trees:
            t.makeSplitKeys()
            t.splitKeys = [n.br.splitKey for n in t.iterNodesNoRoot()]
            #print '\nsplitKeys = %s' % t.splitKeys
    tracks = {}
    kk = []
    for n in theTree.iterInternalsNoRoot():
        theSplitKey = n.br.splitKey
        print "=" * 50
        print "Looking at split %s" % theSplitKey
        kk.append(theSplitKey)
        tracks[theSplitKey] = []
        startTNum = 0
        while len(self.trees) - startTNum >= windowSize:
            print 'trees %4i to %4i: ' % (startTNum, (startTNum + windowSize) - 1),
            theTrees = self.trees[startTNum:(startTNum + windowSize)]

            nTrees = len(theTrees)
            splitCount = 0
            for t in theTrees:
                if theSplitKey in t.splitKeys:
                    splitCount += 1
            print ' nTrees=%3i, splitCount= %3i' % (nTrees, splitCount)
            tracks[theSplitKey].append([startTNum + (0.5 * windowSize), (float(splitCount)/nTrees)])
            startTNum += stride
    if fName:
        f.write("kk = %s\n" % kk)
        f.write("tracks = %s\n" % tracks)
        f.close()

##Ignore
def trackModelThingsFromTree(self, theTree, windowSize=200, stride=100, fName='trackModelThings.py'):

    complaintHead = '\nTrees.trackModelThingsFromTree()'
    gm = [complaintHead]

    gm.append("This method is not working yet.")
    raise Glitch, gm
    
    theTree.makeSplitKeys()

    # Decorate theTree with splitKeys for all nodes.
    for n in theTree.iterNodesNoRoot():
        if n.name:
            n.name += " (%s)" % n.br.splitKey
        else:
            n.name = "(%s)" % n.br.splitKey
    theTree.splitKeys = [n.br.splitKey for n in theTree.iterNodesNoRoot()]
    theTree.draw()
    print "The drawing above shows the splitKeys that are being tracked."

    f = file(fName, 'w')
    textDrawList = theTree.textDrawList()
    for l in textDrawList:
        f.write("#  %s\n" % l)
    f.write("#\n# The drawing above shows the splitKeys that are being tracked.\n")
    f.write("#\n# windowSize=%s, stride=%s.\n" % (windowSize, stride))

    # Remove the decoration from above.
    for n in theTree.iterNodesNoRoot():
        splitKeyStringWithBlank = " (%s)" % n.br.splitKey
        if n.name.endswith(splitKeyStringWithBlank):
            n.name = n.name[: -(len(splitKeyStringWithBlank))]
        else:
            n.name = None
    #theTree.draw()

    # We need to have read in the model comments when we read in the trees for self.
    if not hasattr(self.trees[0], "modelInfo"):
        gm.append("The first tree has no modelInfo.")
        gm.append("The trees should have been read in with")
        gm.append("var.doTreeReadMcmcModelUsageComments turned on.")
        raise Glitch, gm

    mi = self.trees[0].modelInfo


    # Determine whether we need to makeSplitKeys()
    t = self.trees[0]
    if hasattr(t, "splitKeys") and t.splitKeys:
        pass
    else:
        for t in self.trees:
            t.makeSplitKeys()
            t.splitKeys = [n.br.splitKey for n in t.iterNodesNoRoot()]
            #print '\nsplitKeys = %s' % t.splitKeys

    from Trees import Trees  # I need to make Trees objects below.  This very method
                             # is a Trees method, so importing Trees should
                             # generally not be needed.  But I can't do this earlier.
    from TreePartitions import TreePartitions

    # The root buisiness is not implemented yet.  The way I do it
    # should be guided by theTree, the reference tree.  It is rooted
    # on a certain node.  I should print out how many of the input
    # trees were rooted on that node, and what the compCounts were.
    startTNum = 0
    while len(self.trees) - startTNum >= windowSize:
        print "Trees %6i -- %6i" % (startTNum, (startTNum + windowSize))
        tt2 = Trees(self.trees[startTNum:(startTNum + windowSize)], taxNames=self.taxNames)
        tp = TreePartitions(tt2)
        for pNum in range(mi.nParts):
            print "  Part %i, nComps=%i" % (pNum, mi.parts[pNum].nComps)
            if mi.parts[pNum].nComps > 1:
                #for s in tp.splits:
                #    if s.key in theTree.splitKeys:
                #        print "    Split key: %12s, compCounts=%s" % (s.key, s.modelUsage.parts[pNum].compCounts)
                for k in theTree.splitKeys:
                    if tp.splitsHash.has_key(k):
                        s = tp.splitsHash[k]
                        print "    Split key: %12s, compCounts=%s" % (k, s.modelUsage.parts[pNum].compCounts)
        startTNum += stride
        
    f.close()

