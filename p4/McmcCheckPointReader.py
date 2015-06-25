import os
import func
import cPickle
import math
import numpy
import glob
from Glitch import Glitch

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
            #print "Full glob = %s" % fList
            fList = [fName for fName in glob.glob(theGlob) if 
                     os.path.basename(fName).startswith("mcmc_checkPoint")]
            #print fList
            if not fList:
                raise Glitch, "No checkpoints found in this directory."
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
                f = file(mostRecentFileName)
                m = cPickle.load(f)
                f.close()
                self.mm.append(m)
                
            else:
                # get all the files
                for fName in fList:
                    f = file(fName)
                    m = cPickle.load(f)
                    f.close()
                    self.mm.append(m)

                self.mm = func.sortListOfObjectsOn2Attributes(self.mm, "gen", 'runNum')
        else:
            # get the file by name
            f = file(fName)
            m = cPickle.load(f)
            f.close()
            self.mm.append(m)
        if verbose:
            self.dump()

    def dump(self):
        print "McmcCheckPoints (%i checkPoints read)" % len(self.mm)
        print "%12s %12s %12s %12s" % (" ", "index", "run", "gen+1")
        print "%12s %12s %12s %12s" % (" ", "-----", "---", "-----")
        for i in range(len(self.mm)):
            m = self.mm[i]
            #print "    %2i    run %2i,  gen+1 %11i" % (i, m.runNum, m.gen+1)
            print "%12s %12s %12s %12s" % (" ", i, m.runNum, m.gen+1)
            

    def compareSplits(self, mNum1, mNum2, verbose=True, minimumProportion=0.1):
        """Should we be only looking at splits within the 95% ci of the topologies?"""
        m1 = self.mm[mNum1]
        m2 = self.mm[mNum2]
        tp1 = m1.treePartitions
        tp2 = m2.treePartitions
        
        if verbose:
            print "\nMcmcCheckPointReader.compareSplits(%i,%i)" % (mNum1, mNum2)
            print "%12s %12s %12s %12s %12s" % ("mNum", "runNum", "start", "gen+1", "nTrees")
            for i in range(5):
                print "   ---------",
            print
            for mNum in [mNum1, mNum2]:
                print " %10i " % mNum,
                m = self.mm[mNum]            
                print " %10i " % m.runNum,
                print " %10i " % (m.startMinusOne + 1),
                print " %10i " % (m.gen + 1),
                #for i in m.splitCompares:
                #    print i
                print " %10i " % m.treePartitions.nTrees

        asdos = self.compareSplitsBetweenTwoTreePartitions(tp1, tp2, minimumProportion, verbose=verbose)
        asdos2 = self.compareSplitsBetweenTwoTreePartitions(tp2, tp1, minimumProportion, verbose=verbose)
        if math.fabs(asdos - asdos2) > 0.000001:
            print "Reciprocal assdos differs:  %s  %s" % (asdos, asdos2)
            
        if asdos == None and verbose:
                print "No splits > %s" % minimumProportion
        return asdos

        
    def compareSplitsBetweenTwoTreePartitions(tp1, tp2, minimumProportion, verbose=False):
        ret = tp1.compareSplits(tp2, minimumProportion=minimumProportion)
        #print ret
        if ret != []:
            sumOfStdDevs = 0.0
            nSplits = len(ret)
            diffs = []
            for i in ret:
                #print "            %.3f  %.3f    " % (i[2][0], i[2][1]),
                stdDev = math.sqrt(func.variance(i[2]))
                #print "%.5f" % stdDev
                sumOfStdDevs += stdDev
                diffs.append(math.fabs(i[2][0] - i[2][1]))
            quot = sumOfStdDevs/nSplits
            if verbose:
                #print "  %f " % sumOfStdDevs,
                print "     nSplits=%i, average of std devs of splits %.4f " % (nSplits, quot)
                print "     max difference %f, mean difference %f" % (max(diffs), sum(diffs)/nSplits)
            return quot
        else:
            return None

    compareSplitsBetweenTwoTreePartitions = staticmethod(compareSplitsBetweenTwoTreePartitions)

    def compareSplitsAll(self, precision=3, linewidth=120):
        nM = len(self.mm)
        nItems = ((nM * nM) - nM)/2
        results = numpy.zeros((nM, nM), numpy.float)
        vect = numpy.zeros(nItems, numpy.float)
        vCounter = 0
        for mNum1 in range(1, nM):
            for mNum2 in range(mNum1):
                ret = self.compareSplits(mNum1, mNum2, verbose=False)
                #print "+++ ret = %s" % ret
                if ret == None:
                    ret = 0.0
                results[mNum1][mNum2] = ret
                results[mNum2][mNum1] = ret
                vect[vCounter] = ret
                vCounter += 1
                if 0:
                    print " %10i " % mNum1,
                    print " %10i " % mNum2,
                    print "%.3f" % ret

        # Save current numpy printoptions, and restore, below.
        curr = numpy.get_printoptions()
        numpy.set_printoptions(precision=precision, linewidth=linewidth)
        print results
        numpy.set_printoptions(precision=curr['precision'], linewidth=curr['linewidth'])
        
        print "For the %i values in one triangle," % nItems
        print "max =  ", vect.max()
        print "min =  ", vect.min()
        print "mean = ", vect.mean()
        print "var =  ", vect.var()
            
        
    def writeProposalAcceptances(self):
        for m in self.mm:
            m.writeProposalAcceptances()

    def writeSwapMatrices(self):
        for m in self.mm:
            if m.nChains > 1:
                m.writeSwapMatrix()
            
    def writeProposalProbs(self):
        for m in self.mm:
            m.writeProposalProbs()
