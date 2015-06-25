import time,os,string,glob,sys,math
import pf,func
from Glitch import Glitch

##Ignore
def _fixFileName(fName):
    if fName.count('.') or fName.count(' '):
        fName = list(fName)
        for i in range(len(fName)):
            theChar = fName[i]
            if theChar == '.' or theChar == ' ':
                fName[i] = '_'
        fName = string.join(fName, '')
    return fName

longMessage1 = """

(Boring old) Chi-square test for compositional homogeneity
==========================================================

    The statistic is Sum[(Obs - Exp)^2 / Exp],
        where Exp comes from the data.
    The statistic is X^2 (X squared), in the sense used by
        Sokal & Rohlf in 'Biometry', to distinguish it from
        'chi square', which is a distribution, not a statistic.
    This is (mostly) the same as the test in PAUP, using the
        basefreq command.  There are small differences having
        to do with calculation of degrees of freedom.
    Significance is assessed by Chi-square.
    
"""

longMessage2 = """

Tree- and model-based composition fit test
==========================================

    The statistic is Sum[(Obs - Exp)^2 / Exp],
        like the statistic used in the Chi-squared test,
        except that Exp comes from the model, not from the data.
    I call the statistic X^2_m (X squared sub m) since it is like
        the X^2 statistic used in the chi-square test above, but uses
        the m subscript to say that the expected values come from the model.
    Significance is assessed by simulations on the tree and model.
    A critical point of 95% is used to decide if the data fits or not.
"""




def simsForModelFitTests(self, reps=10, seed=None):
    """Do simulations for model fit tests.

    The model fit tests are the Goldman-Cox test, and the tree- and
    model-based composition fit test.  Both of those tests require
    simulations, optimization of the tree and model parameters on the
    simulated data, and extraction of statistics for use in the null
    distribution.  So might as well do them together.  The Goldman-Cox
    test is not possible if there are any gaps or ambiguities, and in
    that case Goldman-Cox simulation stats are not collected.

    Doing the simulations is therefore the time-consuming part, and so
    this method facilitates doing that job in sections.  If you do
    that, set the random number seed to different numbers.  If the
    seed is not set, the process id is used.  (So obviously you should
    explicitly set the seed if you are doing several runs in the same
    process.)  Perhaps you may want to do the simulations on different
    machines in a cluster.  The stats are saved to files.  The output
    files have the seed number attached to the end, so that different
    runs of this method will have different output file names.
    Hopefully.

    When your model uses empirical comps, simulation uses the
    empirical comp of the original data for simulation (good), then
    the optimization part uses the empirical comp of the
    newly-simulated data (also good, I think).  In that case, if it is
    tree-homogeneous, the X^2_m statistic would be identical to the
    X^2 statistic.

    You would follow this method with the modelFitTests() method,
    which uses all the stats files to make null distributions to
    assess significance of the same stats from self."""

    #gm = ['Tree.simsForModelFitTests()']

    # Make a new data object in which to do the sims, so we do not over-write self
    #print "a self.data = %s" % self.data
    #self.data.dump()
    savedData = self.data
    self.data = None # This triggers self.deleteCStuff()
    self.data = savedData.dupe()

    # We need 2 trees, one for sims, and one for evaluations.  We can
    # use self for sims.  Make a copy for evaluations.
    evalTree = self.dupe()
    evalTree.data = self.data

    # make sure all the memory works ...
    self.calcLogLike(verbose=0)
    evalTree.calcLogLike(verbose=0)

    # We can't do the Goldman-Cox test if there are any gaps or
    # ambiguities.
    doGoldmanCox = True
    for a in self.data.alignments:
        if a.hasGapsOrAmbiguities():
            doGoldmanCox = False
            break
    #print "sims doGoldmanCox = %s" % doGoldmanCox

    # Collect info about the observed data
    statsHashList = []  # one for each data part
    for pNum in range(self.data.nParts):
        h = {}
        statsHashList.append(h)
        h['individualNSites'] = []
        h['observedIndividualCounts'] = []
        for j in range(self.data.nTax):
            h['individualNSites'].append(pf.partSequenceSitesCount(self.data.parts[pNum].cPart, j)) # no gaps or qmarks
            h['observedIndividualCounts'].append(self.data.parts[pNum].composition([j]))  
            # (In the line above, its temporarily composition, not counts)
            #print "got seq %i comp = %s' % (j, h['observedIndividualCounts"][-1])
        # At the moment, h['observedIndividualCounts'] has composition,
        # not counts.  So multiply by h['individualNSites']
        for i in range(self.data.nTax):
            for j in range(self.data.parts[pNum].dim):
                h['observedIndividualCounts'][i][j] *= h['individualNSites'][i]

    # We will want to skip any sequences composed of all gaps
    skipTaxNums = []
    for pNum in range(self.data.nParts):
        stn = []
        for tNum in range(self.data.nTax):
            if not statsHashList[pNum]['individualNSites'][tNum]:
                stn.append(tNum)
        skipTaxNums.append(stn)
    #print "skipTaxNums = %s" % skipTaxNums

    if seed == None:
        seed = os.getpid()
    pf.reseedCRandomizer(int(seed))

    # Open up some output files in which to put the sim data
    outfileBaseName = 'sims' # Could be an argument, user-assignable.
    if doGoldmanCox:
        f2Name = outfileBaseName + '_GoldmanStats_%s' % seed
        f2 = open(f2Name, 'w')
        f2.write('# part\tunconstr L\t log like \tGoldman-Cox stat\n')

    f3Name = outfileBaseName + '_CompStats_%s' % seed
    f3 = open(f3Name, 'w')

    # When sims are done when the comp is empirical (whether or not
    # free) we need to re-set the comps based on the newly-simulated
    # data.  So first find out if any comps are empirical.
    hasEmpiricalComps = 0
    for mp in self.model.parts:
        for c in mp.comps:
            if c.spec == 'empirical':
                hasEmpiricalComps = 1
                break
    #print "hasEmpiricalComps=%s" % hasEmpiricalComps

    # Do the sims
    for i in range(reps):
        self.simulate()
        if hasEmpiricalComps:
            evalTree.setEmpiricalComps()   # Set empirical comps based on newly-simulated data
        evalTree.optLogLike(verbose=0)     # The time-consuming part

        if doGoldmanCox:
            if self.data.nParts > 1:
                self.data.calcUnconstrainedLogLikelihood2()
                diff = self.data.unconstrainedLogLikelihood - evalTree.logLike
                f2.write('-1\t%f\t%f\t%f\n' % (self.data.unconstrainedLogLikelihood, evalTree.logLike, diff))
            for pNum in range(self.data.nParts):
                unc = pf.getUnconstrainedLogLike(self.data.parts[pNum].cPart)
                like = pf.p4_partLogLike(evalTree.cTree, self.data.parts[pNum].cPart, pNum, 0) # 0 for getSiteLikes
                diff = unc - like
                f2.write('%i\t%f\t%f\t%f\n' % (pNum, unc, like, diff))

        for pNum in range(self.data.nParts):
            h = statsHashList[pNum]
            # pf.p4_expectedCompositionCounts returns a tuple of tuples
            # representing the counts of the nodes in proper alignment order.
            ret = pf.p4_expectedCompositionCounts(evalTree.cTree, pNum)
            h['expectedIndividualCounts'] = list(ret) # alignment order
            #print h['expectedIndividualCounts']
            h['overallSimStat'] = 0.0
            h['individualSimStats'] = [0.0] * self.data.nTax
            for seqNum in range(self.data.nTax):
                if seqNum in skipTaxNums[pNum]:
                    pass
                else:
                    obsCounts = list(pf.singleSequenceBaseCounts(self.data.parts[pNum].cPart, seqNum))
                    # obsCounts is the counts observed in the simulation.
                    # It assumes that there are no gaps.  If there are gaps, adjust the counts.
                    if h['individualNSites'][seqNum] != self.data.parts[pNum].nChar:
                        factor = float(h['individualNSites'][seqNum]) / self.data.parts[pNum].nChar
                        #print "factor = %s" % factor
                        for j in range(self.data.parts[pNum].dim):
                            obsCounts[j] = float(obsCounts[j]) * factor
                    #print obsCounts
                    for j in range(self.data.parts[pNum].dim):
                        # Avoid dividing by Zero.
                        if h['expectedIndividualCounts'][seqNum][j]:
                            dif = obsCounts[j] - h['expectedIndividualCounts'][seqNum][j]
                            h['individualSimStats'][seqNum] += ((dif * dif) / h['expectedIndividualCounts'][seqNum][j])
                    h['overallSimStat'] += h['individualSimStats'][seqNum]

            f3.write('%i\t' % pNum)
            for seqNum in range(self.data.nTax):
                f3.write('%f\t' % h['individualSimStats'][seqNum])
            f3.write('%f\n' % h['overallSimStat'])
            #print h['overallSimStat']

    if doGoldmanCox:
        f2.close()
    f3.close()

    # Replace the saved data
    self.data = savedData # Since we are replacing an exisiting data, this triggers self.deleteCStuff()



def modelFitTests(self, fName = 'model_fit_tests_out', writeRawStats=0):
    """Do model fit tests on the data.

    The two tests are the Goldman-Cox test, and the tree- and model-
    based composition fit test.  Both require simulations with
    optimizations in order to get a null distribution, and those
    simulations need to be done before this method.  The simulations
    should be done with the simsForModelFitTests() method.

    Self should have a data and a model attached, and be optimized.

    The Goldman-Cox test (Goldman 1993.  Statistical tests of models
    of DNA substitution.  J Mol Evol 36: 182-198.) is a test for
    overall fit of the model to the data.  It does not work if the
    data have gaps or ambiguities.

    The tree- and model-based composition test asks the question:
    'Does the composition implied by the model fit the data?'  If the
    model is homogeneous and empirical comp is used, then this is the
    same as the chi-square test except that the null distribution
    comes from simulations, not from the chi-square distribution.  In
    that case only the question is, additionally, 'Are the data
    homogeneous in composition?', ie the same question asked by the
    chi-square test.  However, the data might be heterogeneous, and
    the model might be heterogeneous over the tree; the tree- and
    model-based composition fit test can ask whether the heterogeneous
    model fits the heterogeneous data.  The composition is tested in
    each data partition, separately.  The test is done both overall,
    ie for all the sequences together, and for individual sequences.

    If you just want a compo homogeneity test with empirical
    homogeneous comp, try the compoTestUsingSimulations() method-- its
    way faster, because there are not optimizations in the sims part.

    Output is verbose, to a file."""

    gm = ['Tree.modelFitTests()']
    self.calcLogLike(verbose=0)
    doOut = True # Usually True.  Set to False for debugging, experimentation, getting individual stats, etc

    # We can't do the Goldman-Cox test if there are any gaps or
    # ambiguities.
    doGoldmanCox = True
    for a in self.data.alignments:
        if a.hasGapsOrAmbiguities():
            doGoldmanCox = False
            break
    #print "test doGoldmanCox = %s" % doGoldmanCox


    rawFName = '%s_raw.py' % fName
    #flob = sys.stderr
    #fRaw = sys.stderr
    if doOut:
        flob = file(fName, 'w')
    else:
        flob = None
    if writeRawStats:
        fRaw = file(rawFName, 'w')
    else:
        fRaw = None
    
    #######################
    # Goldman-Cox stats
    #######################

    # For a two-part data analysis, the first few lines of the
    # sims_GoldmanStats_* file will be like the following.  Its in
    # groups of 3-- the first one for all parts together (part number
    # -1), and the next lines for separate parts.

    ##    # part  unconstr L       log like       Goldman-Cox stat
    ##    -1      -921.888705     -1085.696919    163.808215
    ##    0       -357.089057     -430.941958     73.852901
    ##    1       -564.799648     -654.754962     89.955314
    ##    -1      -952.063037     -1130.195799    178.132761
    ##    0       -362.164119     -439.709824     77.545705
    ##   ... and so on.

    # For a one-part analysis, it will be the same except that one sim
    # gets only one line, starting with zero.
    
    if doGoldmanCox:
        goldmanOverallSimStats = []
        if self.data.nParts > 1:
            goldmanIndividualSimStats = []
            for partNum in range(self.data.nParts):
                goldmanIndividualSimStats.append([])

        import glob
        goldmanFNames = glob.glob('sims_GoldmanStats_*')
        #print "nParts=%s" % self.data.nParts
        #print goldmanFNames
        for fName1 in goldmanFNames:
            f2 = open(fName1)
            aLine = f2.readline()
            if not aLine:
                gm.append("Empty file %s" % fName1)
                raise Glitch, gm
            if aLine[0] != '#':
                gm.append("Expecting a '#' as the first character in file %s" % fName1)
                raise Glitch, gm
            aLine = f2.readline()
            #print "a got line %s" % aLine,
            while aLine:
                if self.data.nParts > 1:
                    splitLine = aLine.split()
                    if len(splitLine) != 4:
                        gm.append("Bad line in Goldman stats file %s" % fName1)
                        gm.append("'%s'" % aLine)
                        raise Glitch, gm
                    if int(splitLine[0]) != -1:
                        gm.append("Bad line in Goldman stats file %s" % fName1)
                        gm.append("First item should be -1")
                        gm.append("'%s'" % aLine)
                        raise Glitch, gm
                    #print splitLine[-1]
                    goldmanOverallSimStats.append(float(splitLine[-1]))

                    aLine = f2.readline()
                    #print "b got line %s" % aLine,
                    if not aLine:
                        gm.append("Premature end to file %s" % fName1)
                        raise Glitch, gm

                for partNum in range(self.data.nParts):
                    splitLine = aLine.split()
                    #print "partNum %i, splitLine=%s" % (partNum, splitLine)
                    if len(splitLine) != 4:
                        gm.append("Bad line in Goldman stats file %s" % fName1)
                        gm.append("'%s'" % aLine)
                        raise Glitch, gm
                    try:
                        splitLine[0] = int(splitLine[0])
                    except ValueError:
                        gm.append("Bad line in Goldman stats file %s" % fName1)
                        gm.append("First item should be the partNum %i" % partNum)
                        gm.append("'%s'" % aLine)
                        raise Glitch, gm
                    if splitLine[0] != partNum:
                        gm.append("Bad line in Goldman stats file %s" % fName1)
                        gm.append("First item should be the partNum %i" % partNum)
                        gm.append("'%s'" % aLine)
                        raise Glitch, gm
                    #for taxNum in range(self.data.nTax):
                    #    print splitLine[taxNum + 1]
                    #print splitLine[-1]
                    if self.data.nParts == 1:
                        goldmanOverallSimStats.append(float(splitLine[-1]))
                    else:
                        goldmanIndividualSimStats[partNum].append(float(splitLine[-1]))

                    aLine = f2.readline()
                    #print "c got line %s" % aLine,
            f2.close()


        #print "goldmanOverallSimStats =", goldmanOverallSimStats
        #print "goldmanIndividualSimStats =", goldmanIndividualSimStats
        #sys.exit()

        if doOut:
            flob.write('Model fit tests\n===============\n\n')
            flob.write('The data that we are testing have %i taxa,\n' % self.data.nTax)

            if len(self.data.alignments) == 1:
                flob.write('1 alignment, ')
            else:
                flob.write('%i alignments, ' % len(self.data.alignments))
            if self.data.nParts == 1:
                flob.write('and 1 data partition.\n')
            else:
                flob.write('and %i data partitions.\n' % self.data.nParts)

            flob.write('The lengths of those partitions are as follows:\n')
            flob.write('                  partNum    nChar \n')
            for i in range(self.data.nParts):
                flob.write('                      %3i   %5i\n' % (i, self.data.parts[i].nChar))
        self.data.calcUnconstrainedLogLikelihood2()
        if doOut:
            flob.write("\nThe unconstrained likelihood is %f\n" % self.data.unconstrainedLogLikelihood)
            flob.write('(This is the partition-by-partition unconstrained log likelihood, \n')
            flob.write('ie the sum of the unconstrained log likes from each partition separately, \n')
            flob.write('and so will not be the same as that given by PAUP, if the data are partitioned.)\n')


            flob.write('\n\nGoldman-Cox test for overall model fit\n')
            flob.write    ('======================================\n')
            flob.write('The log likelihood for these data for this tree is %f\n' % self.logLike)
            flob.write('The unconstrained log likelihood for these data is %f\n' % self.data.unconstrainedLogLikelihood)
        originalGoldmanCoxStat = self.data.unconstrainedLogLikelihood - self.logLike
        if doOut:
            flob.write('The Goldman-Cox statistic for the original data is the difference, %f\n' % originalGoldmanCoxStat)
            if self.data.nParts > 1:
                flob.write('(The unconstrained log likelihood for these data is calculated partition by partition.)\n')
            flob.write('\n')

        if self.data.nParts > 1:
            originalGoldmanCoxStatsByPart = []
            if doOut:
                flob.write('Stats by partition.\n')
                flob.write('part\t unconstrLogL\t log like \tGoldman-Cox stat\n')
                flob.write('----\t ----------\t -------- \t----------------\n')
            for partNum in range(self.data.nParts):
                unc = pf.getUnconstrainedLogLike(self.data.parts[partNum].cPart)
                like = pf.p4_partLogLike(self.cTree, self.data.parts[partNum].cPart, partNum, 0)
                diff = unc - like
                if doOut:
                    flob.write('  %i\t%f\t%f\t   %f\n' % (partNum, unc, like, diff))
                originalGoldmanCoxStatsByPart.append(diff)

        # Do the overall stat
        nSims = len(goldmanOverallSimStats)
        if doOut: flob.write('\nThere were %i simulations.\n\n' % nSims)

        if writeRawStats:
            fRaw.write('# Goldman-Cox null distributions.\n')
            if self.data.nParts > 1:
                fRaw.write('# Simulation stats for overall data, ie for all data partitions combined.\n')
            else:
                fRaw.write('# Simulation stats.\n')
            fRaw.write('goldman_cox_overall = %s\n' % goldmanOverallSimStats)
            if self.data.nParts > 1:
                for partNum in range(self.data.nParts):
                    fRaw.write('# Simulation stats for data partition %i\n' % partNum)
                    fRaw.write('goldman_cox_part%i = %s\n' % (partNum, goldmanIndividualSimStats[partNum]))


        prob =  func.tailAreaProbability(originalGoldmanCoxStat, goldmanOverallSimStats, verbose=0)
        if doOut:
            flob.write( '\n              Overall Goldman-Cox test: ')
            if prob <= 0.05:
                flob.write('%13s' % "Doesn't fit.")
            else:
                flob.write('%13s' %  'Fits.')
            flob.write('    P = %5.3f\n' % prob)

        if self.data.nParts > 1:
            if doOut: flob.write('  Tests for individual data partitions:\n')
            for partNum in range(self.data.nParts):
                prob =  func.tailAreaProbability(originalGoldmanCoxStatsByPart[partNum],
                                                 goldmanIndividualSimStats[partNum], verbose=0)
                if doOut:
                    flob.write( '                               Part %-2i: ' % partNum)
                    if prob <= 0.05:
                        flob.write('%13s' % 'Doesn\'t fit.')
                    else:
                        flob.write('%13s' %  'Fits.')
                    flob.write('    P = %5.3f\n' % prob)


    #########################
    # COMPOSITION
    #########################

    statsHashList = []
    for pNum in range(self.data.nParts):
        h = {}
        statsHashList.append(h)
        h['individualNSites'] = []
        h['observedIndividualCounts'] = []
        for j in range(self.data.nTax):
            #print pf.partSequenceSitesCount(self.data.parts[pNum].cPart, j)
            h['individualNSites'].append(pf.partSequenceSitesCount(self.data.parts[pNum].cPart, j)) # no gaps or qmarks
            #print self.data.parts[pNum].composition([j])
            h['observedIndividualCounts'].append(self.data.parts[pNum].composition([j]))
            # The line above is temporarily composition, not counts
        # pf.expectedCompositionCounts returns a tuple of tuples
        # representing the counts of the nodes in proper alignment order.
        h['expectedIndividualCounts'] = list(pf.p4_expectedCompositionCounts(self.cTree, pNum)) # alignment order

        # At the moment, h['observedIndividualCounts'] has composition,
        # not counts.  So multiply by h['individualNSites']
        for i in range(self.data.nTax):
            for j in range(self.data.parts[pNum].dim):
                h['observedIndividualCounts'][i][j] *= h['individualNSites'][i]

    # We will want to skip any sequences composed of all gaps
    skipTaxNums = []
    for pNum in range(self.data.nParts):
        stn = []
        for tNum in range(self.data.nTax):
            if not statsHashList[pNum]['individualNSites'][tNum]:
                stn.append(tNum)
        skipTaxNums.append(stn)
    #print "skipTaxNums = %s" % skipTaxNums

    # Do the boring old compo chi square test.
    if doOut: flob.write(longMessage1) # explanation ...
    for pNum in range(self.data.nParts):
        h = statsHashList[pNum]
        # Can't use func.xSquared(), because there might be column
        # zeros.
        #print "observedIndividualCounts = %s' % h['observedIndividualCounts"]
        nRows = len(h['observedIndividualCounts'])
        nCols = len(h['observedIndividualCounts'][0])
        theSumOfRows = func._sumOfRows(h['observedIndividualCounts']) # I could have just used nSites, above
        theSumOfCols = func._sumOfColumns(h['observedIndividualCounts'])
        #print theSumOfCols
        isOk = 1
        columnZeros = []
        #for j in range(len(theSumOfRows)):
        #    if theSumOfRows[j] == 0.0:
        #        gm.append("Zero in a row sum.  Programming error.")
        #        raise Glitch, gm
        for j in range(len(theSumOfCols)):
            if theSumOfCols[j] <= 0.0:
                columnZeros.append(j)
        theExpected = func._expected(theSumOfRows, theSumOfCols)
        #print "theExpected = %s" % theExpected
        #print "columnZeros = %s" % columnZeros
        xSq = 0.0
        for rowNum in range(nRows):
            if rowNum in skipTaxNums[pNum]:
                pass
            else:
                xSq_row = 0.0
                for colNum in range(nCols):
                    if colNum in columnZeros:
                        pass
                    else:
                        theDiff = h['observedIndividualCounts'][rowNum][colNum] - theExpected[rowNum][colNum]
                        xSq_row += (theDiff * theDiff) / theExpected[rowNum][colNum]
                xSq += xSq_row
        dof = (nCols - len(columnZeros) - 1) * (nRows - len(skipTaxNums[pNum]) - 1)
        prob = func.chiSquaredProb(xSq, dof)
        if doOut: flob.write('        Part %i: Chi-square = %f, (dof=%i) P = %f\n' % (pNum, xSq, dof, prob))

    for pNum in range(self.data.nParts):
        h = statsHashList[pNum]
        h['overallStat'] = 0.0
        h['individualStats'] = [0.0] * self.data.nTax
        for i in range(self.data.nTax):
            if i in skipTaxNums[pNum]:
                pass # h['individualStats'] stays at zeros
            else:
                for j in range(self.data.parts[pNum].dim):
                    # Avoid dividing by Zero.
                    if h['expectedIndividualCounts'][i][j]:
                        dif = h['observedIndividualCounts'][i][j] - h['expectedIndividualCounts'][i][j]
                        h['individualStats'][i] += ((dif * dif) /h['expectedIndividualCounts'][i][j])
                h['overallStat'] += h['individualStats'][i]

        h['overallSimStats'] = []
        h['individualSimStats'] = []
        for i in range(self.data.nTax):
            h['individualSimStats'].append([])

        if 0:
            print "h['individualNSites'] = %s" % h['individualNSites']
            print "h['observedIndividualCounts'] = %s" % h['observedIndividualCounts']
            print "h['expectedIndividualCounts'] = %s" % h['expectedIndividualCounts']
            print "h['overallStat'] = %s" % h['overallStat']
            print "h['individualStats'] = %s" % h['individualStats']
            raise Glitch, gm



    import glob
    compoFNames = glob.glob('sims_CompStats_*')
    #print compoFNames
    for fName1 in compoFNames:
        f2 = open(fName1)
        aLine = f2.readline()
        if not aLine:
            gm.append("Empty file %s" % fName1)
            raise Glitch, gm
        #print "a got line %s" % aLine,
        while aLine:
            for partNum in range(self.data.nParts):
                h = statsHashList[partNum]
                splitLine = aLine.split()
                if len(splitLine) != (self.data.nTax + 2):
                    gm.append("Bad line in composition stats file %s" % fName1)
                    gm.append("'%s'" % aLine)
                    raise Glitch, gm
                if int(splitLine[0]) != partNum:
                    gm.append("Bad line in composition stats file %s" % fName1)
                    gm.append("First item should be the partNum %i" % partNum)
                    gm.append("'%s'" % aLine)
                    raise Glitch, gm
                #for taxNum in range(self.data.nTax):
                #    print splitLine[taxNum + 1]
                #print splitLine[-1]
                h['overallSimStats'].append(float(splitLine[-1]))
                for i in range(self.data.nTax):
                    h['individualSimStats'][i].append(float(splitLine[i + 1]))
                #raise Glitch, gm

                aLine = f2.readline()
                if not aLine:
                    break
                #print "b got line %s" % aLine,
        f2.close()

    nSims = len(statsHashList[0]['overallSimStats'])
    if doOut:
        flob.write(longMessage2) # Explain tree- and model-based compo fit stat, X^2_m
        flob.write( '    %i simulation reps were used.\n\n' % nSims)

    spacer1 = ' ' * 10
    for partNum in range(self.data.nParts):
        h = statsHashList[partNum]
        if doOut:
            flob.write('Part %-2i:\n-------\n\n' % partNum)
            flob.write('Statistics from the original data\n')
            flob.write('%s%30s: %f\n' % (spacer1, 'Overall observed stat', h['overallStat']))
            flob.write('%s%30s:\n' %  (spacer1, 'Stats for individual taxa'))
            for taxNum in range(self.data.nTax):
                if taxNum not in skipTaxNums[partNum]:
                    flob.write('%s%30s: %f\n' % (spacer1, self.data.taxNames[taxNum], h['individualStats'][taxNum]))
                else:
                    flob.write('%s%30s: skipped\n' % (spacer1, self.data.taxNames[taxNum]))

            flob.write('\nAssessment of fit from null distribution from %i simulations\n' % nSims)
            flob.write('%s%30s:  ' % (spacer1, 'Overall'))
        prob =  func.tailAreaProbability(h['overallStat'], h['overallSimStats'], verbose=0)
        if doOut:
            if prob <= 0.05:
                flob.write('%13s' % 'Doesn\'t fit.')
            else:
                flob.write('%13s' %  'Fits.')
            flob.write('    P = %5.3f\n' % prob)
        #############
        theRet= prob
        #############
        for taxNum in range(self.data.nTax):
            if doOut: flob.write('%s%30s:  ' % (spacer1, self.data.taxNames[taxNum]))
            if taxNum in skipTaxNums[partNum]:
                if doOut: flob.write('%13s\n' % 'skipped.')
            else:
                prob =  func.tailAreaProbability(h['individualStats'][taxNum],
                                                 h['individualSimStats'][taxNum], verbose=0)
                if doOut:
                    if prob <= 0.05:
                        flob.write('%13s' % "Doesn't fit.")
                    else:
                        flob.write('%13s' %  'Fits.')
                    flob.write('    P = %5.3f\n' % prob)

        if writeRawStats:
            fRaw.write('#\n# Tree and model based composition fit test\n')
            fRaw.write('# =========================================\n')
            fRaw.write('# Simulation statistics, ie the null distributions\n\n')
            fRaw.write('# Part %i:\n' % partNum)
            fRaw.write('part%i_overall_compo_null = %s\n' % (partNum, h['overallSimStats']))
            for taxNum in range(self.data.nTax):
                fRaw.write('part%i_%s_compo_null = %s\n' % (partNum,
                                                             _fixFileName(self.data.taxNames[taxNum]),
                                                             h['individualSimStats'][taxNum]))


    if flob and flob != sys.stdout: # Yes, it is possible to close sys.stdout
        flob.close()
    if fRaw and fRaw != sys.stdout:
        fRaw.close()
    return theRet


def compoTestUsingSimulations(self, nSims=100, doIndividualSequences=0, doChiSquare=0, verbose=1):
    """Compositional homogeneity test using a null distribution from simulations.

    This does a compositional homogeneity test on each data partition.
    The statistic used here is X^2, obtained via
    Data.compoChiSquaredTest().

    The null distribution of the stat is made using simulations, so of
    course you need to provide a tree with a model, with optimized
    branch lengths and model parameters.  This is a comp homogeneity
    test, so the model should be tree-homogeneous.
    
    The analysis usually tests all sequences in the data partition
    together (like paup), but you can also 'doIndividualSequences'
    (like puzzle).  Beware that the latter is a multiple simultaneous
    stats test, and so the power may be compromized.

    For purposes of comparison, this test can also do compo tests in
    the style of PAUP and puzzle, using chi-square to assess
    significance.  Do this by turning 'doChiSquare' on.  The compo test
    in PAUP tests all sequences together, while the compo test in
    puzzle tests all sequences separately.  There are advantages and
    disadvantages to the latter-- doing all sequences separately
    allows you to identify the worst offenders, but suffers due to the
    problems of multiple simultaneous stats tests.  There are slight
    differences between the computation of the Chi-square in PAUP and
    puzzle and the p4 version.  The compo test in PAUP (basefreq) does
    the chi-squared test, but if sequences are blank it still counts
    them in the degrees of freedom; p4 does not count blank sequences
    in the degrees of freedom.  Puzzle simply uses the row sums, ie
    the contributions of each sequence to the total X-squared, and
    assesses significance with chi-squared using the number of symbols
    minus 1 as the degrees of freedom.  Ie for DNA dof=3, for protein
    dof=19.  Puzzle correctly gets the composition from sequences with
    gaps, but does not do the right thing for sequences with
    ambiguities like r, y, and so on.  P4 does calculate the
    composition correctly when there are such ambiguities.  So p4 will
    give you the same numbers as paup and puzzle for the chi-squared
    part as long as you don't have blank sequences or ambiguities like
    r and y.

    This uses the Data.compoChiSquaredTest() method to get the
    stats. See the doc string for that method, where it describes how
    zero column sums (ie some character is absent) can be dealt with.
    Here, when that method is invoked, 'skipColumnZeros' is turned on,
    so that the analysis is robust against data with zero or low
    values for some characters.

    For example::

        # First, do a homog opt, and pickle the optimized tree.
        # Here I use a bionj tree, but you could use whatever.
        read('d.nex')
        a = var.alignments[0]
        dm = a.pDistances()
        t = dm.bionj()
        d = Data()
        t.data = d
        t.newComp(free=1, spec='empirical')
        t.newRMatrix(free=1, spec='ones')
        t.setNGammaCat(nGammaCat=4)
        t.newGdasrv(free=1, val=0.5)
        t.setPInvar(free=0, val=0.0)
        t.optLogLike()
        t.name = 'homogOpt'
        t.tPickle()

        # Then, do the test ...
        read('homogOpt.p4_tPickle')
        t = var.trees[0]
        read('d.nex')
        d = Data()
        t.data = d
        t.compoTestUsingSimulations()

        # Output would be something like ...
        # Composition homogeneity test using simulations.
        # P-values are shown.

        #             Part Num       0     
        #            Part Name      all    
        # --------------------    -------- 
        #        All Sequences     0.0000  

        # Or using more sims for more precision, and also doing the
        # Chi-square test for contrast ...
        
        t.compoTestUsingSimulations(nSims=1000, doChiSquare=True)

        # Output might be something like ...
        # Composition homogeneity test using simulations.
        # P-values are shown.
        # (P-values from Chi-Square are shown in parens.)

        #             Part Num       0     
        #            Part Name      all    
        # --------------------    -------- 
        #        All Sequences     0.0140  
        #   (Chi-Squared Prob)    (0.9933) 

    It is often the case, as above, that this test will show
    significance while the Chi-square test does not.
    
    """

    gm = ['Tree.compoTestUsingSimulations()']

    #print "inComp = %s" % self.model.parts[0].comps[0].val

    if not self.data:
        gm.append("No data.  Set the data first.")
        raise Glitch, gm
    if not self.model:
        gm.append("No model.  You need to set the model first.")
        raise Glitch, gm
    self.modelSanityCheck()
    if self.model.isHet:
        gm.append("The model for this tree is tree-heterogeneous.")
        gm.append("This test is not valid for tree-hetero models.")
        raise Glitch, gm

    # Make a new data object in which to do the sims, so we do not over-write self
    #print "a self.data = %s" % self.data
    #self.data.dump()
    savedData = self.data
    self.data = None # This triggers self.deleteCStuff()
    self.data = savedData.dupe()

    #print "b self.data = %s" % self.data
    #self.data.dump()
    #raise Glitch, gm

    # Check for missing sequences in any of the parts.  Missing seq
    # nums go in skips, a list of lists.
    skips = []
    for pNum in range(self.data.nParts):
        skips.append([])
    for pNum in range(self.data.nParts):
        for tNum in range(self.data.nTax):
            nSites = pf.partSequenceSitesCount(self.data.parts[pNum].cPart, tNum) # no gaps, no missings
            if not nSites:
                skips[pNum].append(tNum)

    # Get the original stats from self.data.
    # compoChiSquaredTest(self, verbose=1, skipColumnZeros=0, useConstantSites=1, skipTaxNums=None, getRows=0)
    original = self.data.compoChiSquaredTest(verbose=0,
                                             skipColumnZeros=1,
                                             skipTaxNums=skips,
                                             getRows=doIndividualSequences)
    #print "original =", original

    # Make some empty lists in which to put our stats
    full = []
    if doIndividualSequences:
        rows = []
    for pNum in range(self.data.nParts):
        full.append([])
        if doIndividualSequences:
            onePartRows = []
            for i in range(self.data.nTax):
                onePartRows.append([])
            rows.append(onePartRows)

    # Do the sims
    for i in range(nSims):
        #if i < 5:
        #    print "%i simComp = %s" % (i, self.model.parts[0].comps[0].val)
        self.simulate()
        ret = self.data.compoChiSquaredTest(skipColumnZeros=1,
                                            skipTaxNums=skips, getRows=doIndividualSequences, verbose=0)
        #print "%i ret=%s" % (i, ret)
        for pNum in range(self.data.nParts):
            full[pNum].append(ret[pNum][0])
            if doIndividualSequences:
                for tNum in range(self.data.nTax):
                    if tNum not in skips[pNum]:
                        rows[pNum][tNum].append(ret[pNum][3][tNum])


    # Find the longest part name length, and heading width, so the output looks nice.
    partWid = 8
    for p in self.data.parts:
        if len(p.name) > partWid:
            partWid = len(p.name)
    partWid += 2

    headWid = 20
    for tN in self.data.taxNames:
        if len(tN) > headWid:
            headWid = len(tN)
    headWid += 2
    #headSig = '%-' + `headWid` + 's'
    headSig = '%' + `headWid - 2` + 's  '

    # Get the all-sequences tail area probs
    partTaps = []
    for pNum in range(self.data.nParts):
        partTaps.append(func.tailAreaProbability(original[pNum][0], full[pNum], verbose = 0))

    # Intro
    if verbose:
        print "Composition homogeneity test using simulations."
        print "P-values are shown."
        if doChiSquare:
            print "(P-values from Chi-Square are shown in parens.)"
        print

    # Print the Part Nums and Part Names
    if verbose:
        print headSig % 'Part Num',
        for pNum in range(self.data.nParts):
            print string.center('%i' % pNum, partWid),
        print
        print headSig % 'Part Name',
        for pNum in range(self.data.nParts):
            print string.center('%s' % self.data.parts[pNum].name, partWid),
        print
        print headSig % ('-' * (headWid - 2)),
        for pNum in range(self.data.nParts):
            print string.center('%s' % ('-' * (partWid - 2)), partWid),
        print

    # Print the all-sequences results
    if verbose:
        print headSig % 'All Sequences',
        for pNum in range(self.data.nParts):
            print string.center('%6.4f' % partTaps[pNum], partWid),
        print
        if doChiSquare:
            print headSig % '(Chi-Squared Prob)',
            for pNum in range(self.data.nParts):
                print string.center('(%6.4f)' % original[pNum][2], partWid),
            print


    if doIndividualSequences and verbose:
        print
        #print "Individual sequences"
        #print "--------------------"

        for tNum in range(self.data.nTax):
            print headSig % self.data.taxNames[tNum],
            for pNum in range(self.data.nParts):
                if tNum not in skips[pNum]:
                    ret = func.tailAreaProbability(original[pNum][3][tNum], rows[pNum][tNum], verbose = 0)
                    print string.center('%6.4f' % ret, partWid),
                else:
                    print string.center('%s' % ('-' * 4), partWid),
            print
            if doChiSquare:
                print headSig % ' ',
                for pNum in range(self.data.nParts):
                    dof = self.data.parts[pNum].dim - 1 # degrees of freedom
                    if tNum not in skips[pNum]:
                        ret = func.chiSquaredProb(original[pNum][3][tNum], dof)
                        print string.center('(%6.4f)' % ret, partWid),
                    else:
                        print string.center('%s' % ('-' * 4), partWid),
                print

    # Replace the saved data
    self.data = savedData # Since we are replacing an exisiting data, this triggers self.deleteCStuff()

    return partTaps[0]





def bigXSquaredSubM(self, verbose=False):
    """Calculate the X^2_m stat

    This can handle gaps and ambiguities.

    Column zeros in the observed is not a problem with this stat,
    as we are dividing by the expected composition, and that comes
    from the model, which does not allow compositions with values
    of zero.  """

    if not self.cTree:
        self._commonCStuff(resetEmpiricalComps=True)
    l = []
    for pNum in range(self.data.nParts):
        if verbose:
            print "Part %i" % pNum
            print "======"
        obs = []
        nSites = [] # no gaps or ?
        for taxNum in range(self.nTax):
            thisNSites = pf.partSequenceSitesCount(self.data.parts[pNum].cPart, taxNum)
            comp = self.data.parts[pNum].composition([taxNum])
            for symbNum in range(self.data.parts[pNum].dim):
                comp[symbNum] *= thisNSites
            nSites.append(thisNSites)
            obs.append(comp)
        if verbose:
            print "\n  Observed"
            print " " * 10,
            for symbNum in range(self.data.parts[pNum].dim):
                print "%8s" % self.data.parts[pNum].symbols[symbNum],
            print
            for taxNum in range(self.nTax):
                print "%10s" % self.taxNames[taxNum],
                for symbNum in range(self.data.parts[pNum].dim):
                    print "%8.4f" % obs[taxNum][symbNum],
                print "   n=%i" % nSites[taxNum]
                

        # pf.p4_expectedCompositionCounts returns a tuple of tuples
        # representing the counts of the nodes in proper alignment order.
        exp = list(pf.p4_expectedCompositionCounts(self.cTree, pNum))
        if verbose:
            print "\n  Expected"
            print " " * 10,
            for symbNum in range(self.data.parts[pNum].dim):
                print "%8s" % self.data.parts[pNum].symbols[symbNum],
            print
            for taxNum in range(self.nTax):
                print "%10s" % self.taxNames[taxNum],
                for symbNum in range(self.data.parts[pNum].dim):
                    print "%8.4f" % exp[taxNum][symbNum],
                print "   n=%i" % nSites[taxNum]

        # do the summation
        theSum = 0.0
        for taxNum in range(self.nTax):
            for symbNum in range(self.data.parts[pNum].dim):
                x = obs[taxNum][symbNum] - exp[taxNum][symbNum]
                theSum += (x * x) / exp[taxNum][symbNum]

        l.append(theSum)
        if verbose:
            print "The bigXSquaredSubM stat for this part is %.5f" % theSum
    return l


def compStatFromCharFreqs(self, verbose=False):
    """Calculate a statistic from observed and model character frequencies.

    Call it c_m, little c sub m.

    It is calculated from observed character frequencies and character
    frequencies expected from the (possibly tree-heterogeneous) model.

    It would be the sum of abs(obs-exp)/exp
    """

    if not self.cTree:
        self._commonCStuff(resetEmpiricalComps=True)
    # pf.p4_expectedComposition returns a tuple of tuples of tuples
    # representing the counts of the nodes in proper alignment order.
    exp = pf.p4_expectedComposition(self.cTree)
    l = []
    for pNum in range(self.data.nParts):
        if verbose:
            print "Part %i" % pNum
            print "======"
        obs = []
        for taxNum in range(self.nTax):
            comp = self.data.parts[pNum].composition([taxNum])
            obs.append(comp)
        if verbose:
            print "\n  Observed"
            print " " * 10,
            for symbNum in range(self.data.parts[pNum].dim):
                print "%8s" % self.data.parts[pNum].symbols[symbNum],
            print
            for taxNum in range(self.nTax):
                print "%10s" % self.taxNames[taxNum],
                for symbNum in range(self.data.parts[pNum].dim):
                    print "%8.4f" % obs[taxNum][symbNum],
                print
                

        if verbose:
            print "\n  Expected"
            print " " * 10,
            for symbNum in range(self.data.parts[pNum].dim):
                print "%8s" % self.data.parts[pNum].symbols[symbNum],
            print
            for taxNum in range(self.nTax):
                print "%10s" % self.taxNames[taxNum],
                for symbNum in range(self.data.parts[pNum].dim):
                    print "%8.4f" % exp[pNum][taxNum][symbNum],
                print

        # do the summation
        theSum = 0.0
        for taxNum in range(self.nTax):
            for symbNum in range(self.data.parts[pNum].dim):
                theSum += math.fabs(obs[taxNum][symbNum] - exp[pNum][taxNum][symbNum]) / exp[pNum][taxNum][symbNum]

        l.append(theSum)
        if verbose:
            print "The c_m stat for this part is %.5f" % theSum
    return l
    
    
