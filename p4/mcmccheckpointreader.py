import os
import p4.func
import pickle
import math
import numpy
import glob
from p4.p4exceptions import P4Error


class McmcCheckPointReader(object):

    """Read in and display mcmc_checkPoint files.

    Three options--

    To read in a specific checkpoint file, specify the file name by
    fName=whatever

    To read in the most recent (by os.path.getmtime()) checkpoint
    file, say last=True

    If you specify neither of the above, it will read in all the
    checkPoint files that it finds.

    Where it looks is determined by theGlob, which by default is '*',
    ie everything in the current directory.  If you want to look
    somewhere else, you can specify eg::

        theGlob='SomeWhereElse/*' 

    or, if it is unambiguous, just::

        theGlob='S*/*' 

    So you might say::

        cpr = McmcCheckPointReader(theGlob='*_0.*')

    to get all the checkpoints from the first run, run 0.  Then, you
    can tell the cpr object to do various things.  Eg::

        cpr.writeProposalAcceptances()

    But perhaps the most powerful thing about it is that it allows
    easy access to the checkpointed Mcmc objects, in the list mm.  Eg
    to get the first one, ask for::

        m = cpr.mm[0]

    and m is an Mcmc object, complete with all its records of
    proposals and acceptances and so on.  And the TreePartitions
    object.  No data, tho, of course.

    (Sorry!  -- Lazy documentation.  See the source code for more that it can do.)
    """

    def __init__(self, fName=None, theGlob='*', last=False, verbose=True):
        self.mm = []
        if not fName:
            #fList = [fName for fName in os.listdir(os.getcwd()) if fName.startswith("mcmc_checkPoint")]
            #fList = glob.glob(theGlob)
            # print "Full glob = %s" % fList
            fList = [fName for fName in glob.glob(theGlob) if
                     os.path.basename(fName).startswith("mcmc_checkPoint")]
            # print fList
            if not fList:
                raise P4Error("No checkpoints found in this directory.")
            if last:
                # Find the most recent
                mostRecent = os.path.getmtime(fList[0])
                mostRecentFileName = fList[0]
                if len(fList) > 1:
                    for fName in fList[1:]:
                        mtime = os.path.getmtime(fName)
                        if mtime > mostRecent:
                            mostRecent = mtime
                            mostRecentFileName = fName
                f = open(mostRecentFileName, 'rb')
                m = pickle.load(f)
                f.close()
                self.mm.append(m)

            else:
                # get all the files
                for fName in fList:
                    f = open(fName, 'rb')
                    m = pickle.load(f)
                    f.close()
                    self.mm.append(m)

                self.mm = p4.func.sortListOfObjectsOn2Attributes(
                    self.mm, "gen", 'runNum')
        else:
            # get the file by name
            f = open(fName, 'rb')
            m = pickle.load(f)
            f.close()
            self.mm.append(m)
        if verbose:
            self.dump()

    def dump(self, extras=False):
        print("McmcCheckPoints (%i checkPoints read)" % len(self.mm))
        if extras:
            print("%12s %12s %12s %12s %12s %12s %12s" % (
                " ", "index", "run", "gen+1", "cpInterval", "sampInterv", "nSamps"))
            print("%12s %12s %12s %12s %12s %12s %12s" % (
                " ", "-----", "---", "-----", "----------", "----------", "------"))
            for i in range(len(self.mm)):
                m = self.mm[i]
                assert m.checkPointInterval % m.sampleInterval == 0
                if m.simTemp:
                    thisNSamps = m.treePartitions.nTrees
                else:
                    thisNSamps = int(m.checkPointInterval /  m.sampleInterval)
                    assert thisNSamps == m.treePartitions.nTrees
                # print "    %2i    run %2i,  gen+1 %11i" % (i, m.runNum, m.gen+1)
                print("%12s %12s %12s %12s %12s %12s %12s" % (
                    " ", i, m.runNum, m.gen + 1, m.checkPointInterval, m.sampleInterval, thisNSamps))
        else:
            print("%12s %12s %12s %12s %12s" % (
                " ", "index", "run", "gen+1", "nSamps"))
            print("%12s %12s %12s %12s %12s" % (
                " ", "-----", "---", "-----", "------"))
            for i in range(len(self.mm)):
                m = self.mm[i]
                assert m.checkPointInterval % m.sampleInterval == 0
                if hasattr(m, "simTemp") and m.simTemp:
                    thisNSamps = m.treePartitions.nTrees
                else:
                    thisNSamps = int(m.checkPointInterval /  m.sampleInterval)
                    assert thisNSamps == m.treePartitions.nTrees
                # print "    %2i    run %2i,  gen+1 %11i" % (i, m.runNum, m.gen+1)
                print("%12s %12s %12s %12s %12s" % (
                    " ", i, m.runNum, m.gen + 1, thisNSamps))
            

    def compareSplits(self, mNum1, mNum2, verbose=True, minimumProportion=0.1):
        """Do the TreePartitions.compareSplits() method between two checkpoints 

        Args:
            mNum1, mNum2 (int): indices to Mcmc checkpoints in self

        Returns:
            a tuple of asdoss and the maximum difference in split supports

        """

        # Should we be only looking at splits within the 95% ci of the topologies?
        m1 = self.mm[mNum1]
        m2 = self.mm[mNum2]
        tp1 = m1.treePartitions
        tp2 = m2.treePartitions

        if verbose:
            print("\nMcmcCheckPointReader.compareSplits(%i,%i)" % (mNum1, mNum2))
            print("%12s %12s %12s %12s %12s" % ("mNum", "runNum", "start", "gen+1", "nTrees"))
            for i in range(5):
                print("   ---------", end=' ')
            print()
            for mNum in [mNum1, mNum2]:
                print(" %10i " % mNum, end=' ')
                m = self.mm[mNum]
                print(" %10i " % m.runNum, end=' ')
                print(" %10i " % (m.startMinusOne + 1), end=' ')
                print(" %10i " % (m.gen + 1), end=' ')
                # for i in m.splitCompares:
                #    print i
                print(" %10i " % m.treePartitions.nTrees)

        asdos, maxDiff, meanDiff = p4.func._compareSplitsBetweenTwoTreePartitions(
            tp1, tp2, minimumProportion, verbose=verbose)
        asdos2, maxDiff2, meanDiff2= p4.func._compareSplitsBetweenTwoTreePartitions(
            tp2, tp1, minimumProportion, verbose=verbose)
        if math.fabs(asdos - asdos2) > 0.000001:
            print("Reciprocal assdos differs:  %s  %s" % (asdos, asdos2))

        if asdos == None and verbose:
            print("No splits > %s" % minimumProportion)
        return asdos, maxDiff, meanDiff


    def compareSplitsAll(self, precision=3, linewidth=120):
        """Do func.compareSplitsBetweenTreePartitions() for all pairs

        Output is verbose.  Shows 
        - average standard deviation of split frequencies (or supports), like MrBayes
        - maximum difference between split supports from each pair of checkpoints, like PhyloBayes

        Returns:
            None

        """
        tpp = [m.treePartitions for m in self.mm]
        p4.func.compareSplitsBetweenTreePartitions(tpp)

    def writeProposalAcceptances(self):
        for m in self.mm:
            m.writeProposalAcceptances()

    def writeSwapMatrices(self):
        for m in self.mm:
            if m.nChains > 1:
                m.writeSwapMatrix()

    def writeSwapVectors(self):
        for m in self.mm:
            if m.nChains > 1:
                m.writeSwapVector()

    def writeProposalProbs(self):
        for m in self.mm:
            m.writeProposalProbs()
