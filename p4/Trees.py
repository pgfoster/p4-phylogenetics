import sys,os,math
import func
from Var import var
from Tree import Tree
from Glitch import Glitch


class Trees(object):
    """A bunch of trees, all with the same taxNames.

    This class would be good for doing things that you might want to
    do with a bunch of trees rather than just a single tree.  For
    example, you might want to write out the trees to a nexus file but
    use a translation to keep the size of the file down -- this class
    knows how to do that.

    You can start this up in two ways::

        tt = Trees()           # All the trees in var.trees
        tt = Trees(myTreeList) # All the trees in myTreeList

    This class insists on knowing the taxNames.  If the list of
    taxNames is provided when you instantiate a Trees object, then it
    uses that list.  (A good place to get such a list might be an
    alignment object, as myAlignment.taxNames).

      tt = Trees(taxNames=myTaxNamesList)

    Otherwise it looks for a taxNames attribute in the first tree.
    """

    from Trees_converge import trackSplitsFromTree,trackModelThingsFromTree

    def __init__(self, trees=None, taxNames=[]):
        gm = ['Trees()']
        if trees == None:
            if not var.trees:
                gm.append("Arg trees is not given or is empty, and var.trees is empty.  No trees.")
                raise Glitch, gm
            else:
                self.trees = var.trees
        else:
            if type(trees) != type([]):
                gm.append("If provided, the trees arg should be a list.")
                raise Glitch, gm
            if not len(trees):
                gm.append("The list of trees appears to be empty.")
                raise Glitch, gm
            if not isinstance(trees[0], Tree):
                gm.append("The first item in the input list is not a Tree object.")
                raise Glitch, gm
            self.trees = trees
        if taxNames:
            if type(taxNames) == type([]) and len(taxNames) and type(taxNames[0]) == type('string'):
                pass
            else:
                gm.append("If provided, arg taxNames should be a list of (at least one) string(s).")
        else:
            if trees:
                if trees[0].taxNames:
                    taxNames = trees[0].taxNames
            else:
                if var.trees[0].taxNames:
                    taxNames = var.trees[0].taxNames
        if not taxNames:
            gm.append("I can't find a taxNames list in the input trees. ")
            gm.append("In this case you need to feed this a taxNames list when you instantiate.")
            gm.append("Often you can get a good taxNames from yourAlignment.taxNames.")
            gm.append("(The order often matters, and generally should be same for all your analyses.)")
            raise Glitch, gm
        self.taxNames = None
        self.setTaxNames(taxNames)
        self.data = None
        self.writeBranchLengths = True

    nTax = property(lambda self: len(self.taxNames))

    def setTaxNames(self, theTaxNames=None):
        """Sets and checks taxNames for self and all trees in self.

        You can provide the taxNames as the arg theTaxNames, or if you
        do not then it looks for a taxNames list in the first tree.
        If it does not find one there, it makes a list from the first
        tree, and sorts it.

        This propagates self.taxNames to all trees in self.trees, and
        then calls self.checkTaxNames()."""

        gm = ['Trees.setTaxNames()']
        if 0:
            print gm[0]
            print '    theTaxNames=%s' % theTaxNames

        if theTaxNames:
            self.taxNames = theTaxNames
        else:
            # Get a list from the first tree.
            t = None
            if len(self.trees):
                t = self.trees[0]
            if not t:
                print gm[0]
                print "    No trees?"
            if t.taxNames:
                self.taxNames = t.taxNames
            else:
                # Make up a list from node names.
                self.taxNames = []
                for n in t.nodes:
                    if n.isLeaf and n.name:
                        self.taxNames.append(n.name)
                self.taxNames.sort()

        # propagate and check
        for t in self.trees:
            t.taxNames = self.taxNames
        self.checkTaxNames()

    def checkTaxNames(self):
        """Check that all trees have all taxNames."""

        gm = ['Trees.checkTaxNames()']
        if not self.taxNames:
            gm.append("No taxNames.")
            raise Glitch, gm
        for t in self.trees:
            if t.taxNames != self.taxNames:
                if t.name:
                    gm.append("Tree %s taxNames is not the same object as self.taxNames." % t.name)
                else:
                    gm.append("Tree taxNames is not the same object as self.taxNames.")
                raise Glitch, gm
        for t in self.trees:
            t.checkTaxNames()


    def dump(self):
        """Print summary info about self."""
        
        print "Trees dump."
        print "  There are %i trees." % len(self.trees)
        print "  nTax is %i" % self.nTax
        print "  taxNames"
        if self.nTax < 11:
            for i in range(self.nTax):
                print "    %s" % self.taxNames[i]
        else:
            for i in range(10):
                print "    %s" % self.taxNames[i]
            print "    ... (and more)"
        if self.data:
            print "  There is a data object attached."
        else:
            print "  There is no data object attached."

    def writeNexus(self, fName=None, append=0, withTranslation=0, writeTaxaBlock=1, likeMcmc=0):
        """Write the trees out in NEXUS format, in a trees block.

        If *fName* is None, the default, it is written to `sys.stdout`.

        #NEXUS is written unless we are *append*-ing.
        """

        gm = ['Trees.writeNexus()']

        if withTranslation:
            if not self.taxNames:
                gm.append("No taxNames.  Set taxNames if you want withTranslation.")
                raise Glitch, gm

            translationHash = {}
            i = 1
            for tName in self.taxNames:
                translationHash[tName] = i
                i += 1
        else:
            translationHash = None

        if fName == None or fName == sys.stdout:
            f = sys.stdout
            if not append:
                f.write('#NEXUS\n\n')
        else:
            if append:
                import os
                if os.path.isfile(fName):
                    try:
                        f = open(fName, 'a')
                    except IOError:
                        gm.append("Can't open %s for appending." % fName)
                        raise Glitch, gm
                else:
                    if 0:
                        print gm[0]
                        print "    'append' is requested,"
                        print "    but '%s' is not a regular file (doesn't exist?)." \
                              % fName
                        print "    Writing to a new file instead."
                    try:
                        f = open(fName, 'w')
                        f.write('#NEXUS\n\n')
                    except IOError:
                        gm.append("Can't open %s for writing." % fName)
                        raise Glitch, gm

            else:
                try:
                    f = open(fName, 'w')
                    f.write('#NEXUS\n\n')
                except IOError:
                    gm.append("Can't open %s for writing." % fName)
                    raise Glitch, gm

        if writeTaxaBlock:
            if self.taxNames:
                f.write('begin taxa;\n')
                f.write('  dimensions ntax=%s;\n' % self.nTax)
                f.write('  taxlabels')
                for tN in self.taxNames:
                    f.write(' %s' % func.nexusFixNameIfQuotesAreNeeded(tN))
                f.write(';\nend;\n\n')
            else:
                gm.append("writeTaxaBlock is set, but there is no taxNames.  How did *that* happen?!?")
                raise Glitch, gm




        f.write('begin trees;\n')

        # write the "translate" command
        if withTranslation:
            f.write('    translate\n')
            for i in range(self.nTax - 1):
                f.write('        %3i %s,\n' % (i + 1, func.nexusFixNameIfQuotesAreNeeded(self.taxNames[i])))
            f.write('        %3i %s\n' % (self.nTax, func.nexusFixNameIfQuotesAreNeeded(self.taxNames[-1])))
            f.write('    ;\n')

        # write the models comment
        if self.trees:
            first = self.trees[0]
            if likeMcmc:
                if not withTranslation:
                    gm.append("withTranslation is turned off, but likeMcmc is turned on.")
                    gm.append("This will cause grief in TreePartitions, so is prohibited.")
                    gm.append("Both on or both off, please.")
                    raise Glitch, gm
                f.write('    [&&p4 models p%i' % first.model.nParts)
                for pNum in range(first.model.nParts):
                    f.write(' c%i.%i' % (pNum, first.model.parts[pNum].nComps))
                    f.write(' r%i.%i' % (pNum, first.model.parts[pNum].nRMatrices))
                    f.write(' g%i.%i' % (pNum, first.model.parts[pNum].nGdasrvs))
                f.write(']\n')


        # write each tree
        for t in self.trees:
            if not self.writeBranchLengths:
                t.stripBrLens()
            
            if t.logLike:
                f.write('    [logLike for tree %s is %f]\n' % (t.name, t.logLike))

            f.write('    tree %s = [&U] ' % func.nexusFixNameIfQuotesAreNeeded(t.name))
            if t.recipWeight:
                #if t.recipWeight == 1:
                #    f.write('[&W 1] ')
                #else:
                f.write('[&W 1/%i] ' % t.recipWeight)
            if hasattr(t, 'weight'):
                f.write('[&W %f] ' % t.weight)

            t.writeNewick(f, withTranslation=withTranslation, translationHash=translationHash, doMcmcCommandComments=likeMcmc)

        f.write('end;\n\n')
        if f != sys.stdout:
            f.close()


    def writeNewick(self, fName='intree'):
        """Write out the list of trees in Newick format.

        Good for phylip and puzzle.  Long names are not modified,
        which might get you into trouble.
        """

        gm = ['Trees.writeNewick()']
        if fName == None or fName == sys.stdout:
            flob = sys.stdout
        else:
            try:
                flob = open(fName, 'w')
            except IOError:
                gm.append("Can't open %s for writing." % fName)
                raise Glitch, gm

        #flob.write(' %i\n' % len(self.trees))
        for t in self.trees:
            t.writeNewick(flob)
        #flob.write('\n')
        if fName != sys.stdout:
            flob.close()

    def getTreesWithSplit(self, someTaxa):
        """Return a list of trees that have a split.

        The arg *someTaxa* is a list of taxNames that define the split
        that you are looking for.  It returns a list, not a Trees
        object.  (If no trees are found, it returns an empty list.)"""

        sk = func.getSplitKeyFromTaxNames(self.taxNames, someTaxa)
        foundTrees = []
        for t in self.trees:
            t.makeSplitKeys()
            for n in t.nodes:
                if n.br:
                    if n.br.splitKey:
                        if n.br.splitKey == sk:
                            foundTrees.append(t)
                            break
                    else:
                        gm = ['Trees.getTreesWithSplit()']
                        gm.append("No node.br.splitKey?")
                        raise Glitch, gm
        return foundTrees


    def topologyDistanceMatrix(self, metric='sd'):
        """Returns a DistanceMatrix object showing topology distances.

        Uses the :meth:`Tree.Tree.topologyDistance` method to compare
        trees.  That method returns distances between pairs of trees,
        and this method simply collates those distances into a
        :class:`DistanceMatrix.DistanceMatrix` object, which is
        returned.

        See :meth:`Tree.Tree.topologyDistance` for an explanation of the
        different metrics.  The metrics are given in
        ``var.topologyDistanceMetrics``"""

        from DistanceMatrix import DistanceMatrix
        d = DistanceMatrix()
        d.names = []
        for t in self.trees:
            d.names.append(t.name)
        d.dim = len(d.names)
        #print d.names
        #sys.exit()
        d.matrix = []
        for i in range(d.dim):
            d.matrix.append([0.0] * d.dim)

        for i in range(d.dim):
            t1 = self.trees[i]
            for j in range(i + 1, d.dim):
                t2 = self.trees[j]
                theDist = t1.topologyDistance(t2, metric=metric)
                d.matrix[i][j] = theDist
                d.matrix[j][i] = theDist
        return d

    def inputTreesToSuperTreeDistances(self, inputTrees, latex=False):
        """Make a table of supertree distances.

        See the method
        :meth:`Tree.Tree.inputTreesToSuperTreeDistances`, which gives
        supertree distances from a set of input trees to a single
        Tree.  This Trees method does the same for a bunch of trees,
        outputting the results in a table.

        And it returns the results as a list of lists as well.
        """
        nn = []
        sdd = []
        qdd = []
        for t in self.trees:
            sd,qd = t.inputTreesToSuperTreeDistances(inputTrees)
            nn.append(t.name)
            sdd.append(sd)
            qdd.append(qd)
        longest = 0
        for n in nn:
            if len(n) > longest:
                longest = len(n)
        name_sig = " %%-%is" % longest

        longest = 0
        for n in sdd:
            s = "%s" % n
            if len(s) > longest:
                longest = len(s)
        sd_sig = " %%-%is" % longest

        longest = 0
        for n in qdd:
            s = "%s" % n
            if len(s) > longest:
                longest = len(s)
        qd_sig = " %%-%is" % longest

        if latex:
            print "\\begin{center}"
            print "\\begin{tabular}{lrr} \\toprule"
            sig = "%s & %s & %s \\\\" % (name_sig, sd_sig, qd_sig)
        else:
            sig = "%s %s  %s" % (name_sig, sd_sig, qd_sig)
        results = []
        for tNum in range(len(self.trees)):
            nm = nn[tNum]
            sd = sdd[tNum]
            qd = qdd[tNum]
            print sig % (nm, sd, qd)
            results.append([nm, sd, qd])
        if latex:
            print "\\bottomrule"
            print "\\end{tabular}"
            print "\\end{center}"
        return results

    def treeProbabilities(self, writeNexus=False):
        """Order the trees in self by frequency of occurrence.

        This provides a posterior probability for the tree topology if
        trees are sampled from an MCMC. The tree with the highest number
        of occurrences is the MAP tree.

        If *writeNexus* is True a Nexus tree file called
        ``treeProbabilities.nex`` is written with tree posterior
        probabilities as weights.  If that file exists, it will be
        silently over-written!

        The tree probs list of lists is returned.

        (Thanks to Cymon Cox for this method.)

        """

        gm = ["Trees.treeProbabilities()"]

        # Make the splitKey dictionary
        skd = {}
        for t in self.trees:
            t.makeSplitKeys()
            skk = [n.br.splitKey for n in t.iterNodesNoRoot()]
            skk.sort()
            skk = tuple(skk)
            it = skd.get(skk)
            if it:
                it[0] += 1
            else:
                dTree = t.dupe()     # Don't mess with the trees in self
                dTree.stripBrLens()
                skd[skk] = [1, dTree]

        #for v in skd.itervalues():
        #    print v

        #Flatten to a list, order by counts
        skl = skd.values()
        skl = func.sortListOfListsOnListElementNumber(skl, 0)
        skl.reverse()

        # Check that we have all the trees
        nTrees = sum([i[0] for i in skl])
        if nTrees != len(self.trees):
            gm.append("Programming error: all trees are not accounted for...")
            raise Glitch, gm

        # Calculate frequency from counts
        trprobs = [(float(count)/float(nTrees), tree) for count, tree in skl]
        if writeNexus:
            count = 1
            for prob, tree in trprobs:
                tree.weight = prob
                tree.name = "t_%i" % count
                count += 1
            tt = Trees([i[1] for i in trprobs])
            tt.writeNexus(fName="treeProbabilities.nex", withTranslation=True)

        return trprobs
        
                       

    def consel(self, rankByInputOrder=False, clobber=False, quiet=1, tidy=1, returnResults=False):
        """Use Shimo's consel programs to compare trees.

        The trees in self should all be optimized, and all should have
        models attached.

        Self, the Trees object, should have a Data object attached,
        as self.data.

        The default is to rank the output by support.  You
        might rather have the output be in the same order as the
        input.  Which output format is determined by the
        'rankByInputOrder' arg.

        The analysis produces various files.  Whether they are
        overwritten if they already exist is controlled by the arg
        'clobber'.  The default is to refuse to overwrite.

        The consel programs produce lots of informative rubbish on the
        screen.  Setting 'quiet=1' or 'quiet=True', (the default)
        sends that output to the bit bucket.  The output of consel is
        a table, and that table is output to the screen, unless you
        set 'quiet=2', in which case it is not output to the screen.

        The 'tidy' arg says whether to delete intermediate files,
        including consel_out.ci, consel_out.pv, consel_out.rmt,
        consel_out.vt, and the site likelihoods 'sitelh' file.  By
        default it is 'tidy=1', (equivalent to 'tidy=True'), meaning
        delete them.  If you set 'tidy=2' then even the ouput file
        'conselOut' is deleted.

        Usually the results table will be printed to the screen
        (unless quiet=2), and also put in a file 'conselOut' (unless
        tidy=2), and nothing is returned.  If you want the numbers,
        set returnResults=True, and then this method will return a
        list of tuples of strings of the results table contents.

        """

        gm = ['Trees.consel()']

        if not self.trees or len(self.trees) == 0:
            gm.append("No trees?")
            raise Glitch, gm

        # Check if consel is installed
        progs = ['makermt', 'consel', 'catpv']
        for progName in progs:
            if func.which(progName):
                pass
            else:
                gm.append("The programs")
                for p in progs:
                    gm.append("        %s" % p)
                gm.append("need to be in your path.")
                gm.append("Can't find %s" % progName)
                raise Glitch, gm
        
        # Check for bad arg vals
        # False equates to zero, and True equates to 1.  True does not equate to 2
        if quiet not in [0, 1, 2]:
            gm.append("arg 'quiet' should be set to one of 0 (or False), 1 (or True), or 2.  Got '%s'" % quiet)
            raise Glitch, gm
        if tidy not in [0, 1, 2]:
            gm.append("arg 'tidy' should be set to one of 0 (or False), 1 (or True), or 2.  Got '%s'" % tidy)
            raise Glitch, gm

        # Can we make siteLikes?
        allTreesHaveData = True
        for t in self.trees:
            if not t.data:
                allTreesHaveData = False
                break
        if not allTreesHaveData:
            if not self.data:
                gm.append("You need to 'myTreesObject.data = myDataObject'")
                raise Glitch, gm
        for t in self.trees:
            if not t.model:
                gm.append("Tree %s has no model." % t.name)
                raise Glitch, gm

        # Attach self.data to the trees, if needed
        if not allTreesHaveData:
            for t in self.trees:
                if not t.data:
                    t.data = self.data

        # Are we about to clobber?
        #for fName in ['siteLikes.mt', 'consel_out.ci', 'consel_out.pv', 'consel_out.rmt', 'consel_out.vt']:
        for fName in ['siteLikes.sitelh', 'consel_out.ci', 'consel_out.pv', 'consel_out.rmt', 'consel_out.vt']:
            if os.path.exists(fName):
                if clobber:
                    os.remove(fName)
                else:
                    gm.append("Refusing to overwrite file %s" % fName)
                    raise Glitch, gm

        #print "Calculating siteLikes ..."
        for t in self.trees:
            t.getSiteLikes()
            t.deleteCStuff()    # Be memory efficient, but there is still a lot of inefficient re-malloc'ing.
        nSiteLikes = len(self.trees[0].siteLikes)

        # Write the site likes to a file in xx.sitelh, ie puzzle format
        #print "Writing siteLikes ..."
        f = file('siteLikes.sitelh', 'w')
        f.write('%i %i\n' % (len(self.trees), nSiteLikes))
        for i in range(len(self.trees)):
            f.write('t%i\t' % i)
            t = self.trees[i]
            f.write('%12.8f' % math.log(t.siteLikes[0]))
            for j in range(1, nSiteLikes):
                f.write(' %12.8f' % math.log(t.siteLikes[j]))
            f.write('\n')
        f.close()

        #print "Invoking makermt, consel, and catpv ..."
        if quiet:
            os.system('makermt --puzzle siteLikes consel_out > /dev/null')
            os.system('consel consel_out > /dev/null')
        else:
            os.system('makermt --puzzle siteLikes consel_out')
            os.system('consel consel_out')
        if rankByInputOrder:
            if quiet in [0, 1]:
                os.system('catpv -s 1 consel_out |tee conselOut')
            else:  # quiet=2:
                os.system('catpv -s 1 consel_out > conselOut')
        else:
            if quiet in [0, 1]:
                os.system('catpv consel_out |tee conselOut')
            else:  # quiet=2
                os.system('catpv consel_out > conselOut')

        if returnResults:
            fh = open("conselOut", "r")
            lines = fh.readlines()
            #for l in lines:
            #    print "x ", l,
            fh.close()

            # The output from catpv has both leading and trailing blank lines.

            # The header differs depending on whether rankByInputOrder is
            # turned on -- if so then there is an extra header line, as
            # shown below ("# 0 1+ 2 9 ...") that is not there when
            # rankByInputOrder is not turned on.

            """

            # reading consel_out.pv
            #   0    1+     2      9     10  |     3      4      5      6      7      8 
            # rank item    obs     au     np |     bp     pp     kh     sh    wkh    wsh |
            #    1    1 -1594.0  1.000  1.000 |  1.000  1.000  1.000  1.000  1.000  1.000 |
            #    2    2 1594.0  2e-06  2e-06 |      0      0      0      0      0      0 |

            """
            for line in lines:
                if not line.startswith("#"):
                    lines.remove(line)
            linesToSkip = 2
            if rankByInputOrder:
                linesToSkip += 1
            results = []
            for line in lines[linesToSkip:]: 
                a, b = line.split("|")[:2]
                results.append(tuple(a.split()[1:] + b.split()))
            
        
        if tidy:
            #os.system('rm consel_out.ci consel_out.pv consel_out.rmt consel_out.vt siteLikes.mt')
            os.system('rm consel_out.ci consel_out.pv consel_out.rmt consel_out.vt siteLikes.sitelh')
            if tidy == 2:
                os.system('rm conselOut')
        if returnResults:
            return results


if var.usePfAndNumpy:
    import pf
    def rell(self, bootCount=10000, seedIsPid=1):
        """This compares several trees by the RELL method.

        The trees in self should all be optimized, and all should have
        models attached.

        Self, the Trees object, should have a Data object attached,
        as self.data.

        Calculations are done in C, so it is fast.  It uses the gsl random
        number generator.  Setting the *seedIsPid* (ie turning it on) is
        appropriate, and is the default.  """

        gm = ['Trees.rell()']

        if not self.trees or len(self.trees) == 0:
            gm.append("No trees?")
            raise Glitch, gm

        # Can we make siteLikes?
        if not self.data:
            gm.append("No data.  You need to 'myTreesObject.data = myDataObject'")
            raise Glitch, gm
        for t in self.trees:
            if not t.model:
                gm.append("Tree %s has no model." % t.name)
                raise Glitch, gm

        # Attach self.data to the trees, if needed
        for t in self.trees:
            if not t.data:
                t.data = self.data

        #print "Calculating siteLikes ..."
        for t in self.trees:
            t.getSiteLikes()
            t.deleteCStuff()    # Be memory efficient, but there is still a lot of inefficient re-malloc'ing.

        nTrees = len(self.trees)
        for t in self.trees:
            t.siteLikes = map(math.log, t.siteLikes)
        nChar = len(self.trees[0].siteLikes)

        theSeed = 0
        if seedIsPid:
            theSeed = os.getpid()

        #print "Setting C-memory..."
        rStuff = pf.setRellMemory(nTrees, nChar, theSeed)
        for i in range(nTrees):
            t = self.trees[i]
            for j in range(nChar):
                pf.pokeRellMemory(i, j, t.siteLikes[j], rStuff)

        #print "Doing bootstrap ..."
        winners = pf.rell(bootCount, rStuff)  # see data.c
        pf.freeRellMemory(rStuff)

        #print winners

        print "\nRELL bootstrap results"
        print "======================"
        for i in range(nTrees):
            t = self.trees[i]
            print "%3i   %20s  %1.3f" % (i, t.name, (float(winners[i])/float(bootCount)))


    Trees.rell = rell
    del(rell)
    
           
