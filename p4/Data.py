from Alignment import Alignment
import sys,time,os
import pf,func
from Var import var
from Glitch import Glitch


class Data:
    """All the alignments that you want to work with, in one place.

    Initialize this with one of
      - nothing (or None),
      - a list of Alignment objects, or
      - a single Alignment object.

    If you initialize with nothing (or None), then all alignments in
    var.alignments are used.  If you initialize with a list of
    alignments, then that is used.  You can initialize with an empty
    list to get an empty Data object.

    """
    
    #def __del__(self, freeData=pf.freeData, dp_freeData=pf.dp_freeData):
    def __del__(self, freeData=pf.freeData):
        if self.alignments:
            for a in self.alignments:
                a.parts = []
        self.alignments = None
        if self.parts:
            #print len(self.parts)
            for p in self.parts:
                #if p.cPart:
                #    freePart(p.cPart)  # this is not as good as p.__del__(), as it leaves some memory un-freed
                p.__del__()
        #else:
        #    print 0
        if self.cData:
            if self.doDataPart:
                dp_freeData(self.cData)
            else:
                freeData(self.cData)
            self.cData = None
        self.parts = None
        self.taxNames = None

##    ##Ignore
##    def wipe(self):
##        if self.cData:
##            freeData(self.cData)
##            self.cData = None
##        for p in self.parts:
##            if p.cPart:
##                pf.freePart(p.cPart)
##                p.cPart = None
##            del(p)
##        self.taxNames = None
##        for a in self.alignments:
##            del(a)
##        self.alignments = None


    def __init__(self, alignments=None):
        
        gm = ['Data.__init__()']
        self.nParts = 0
        self.parts = []
        self.alignments = []
        self.nTax = 0
        self.taxNames = []
        self.cData = None
        self.unconstrainedLogLikelihood = None
        if alignments:
            if isinstance(alignments, Alignment):
                #Passed in a single alignment object not a list
                alignments = [alignments]
            else:
                if type(alignments) != type([]):
                    gm.append("The 'alignments' arg should be a list or a single Alignment object.")
                    raise Glitch, gm
                for a in alignments:
                    if isinstance(a, Alignment):
                        pass
                    else:
                        gm.append("Something in the 'alignments' arg was not an Alignment.")
                        raise Glitch, gm
            self._fill(alignments)
        elif alignments == []:
            pass
        elif var.alignments:
            self._fill(var.alignments)

        # temporary!  Only used in __del__()
        self.doDataPart = var.doDataPart


    def dump(self):
        """Print rubbish about self."""
        
        print "Data dump"
        if self.nParts == 1:
            if var.doDataPart:
                print "    There is 1 dataPart"
            else:
                print "    There is 1 part"
        else:
            if var.doDataPart:
                print "    There are %i dataParts" % self.nParts
            else:
                print "    There are %i parts" % self.nParts

        for p in self.parts:
            print "        name=%s, nChar %i, dataType %s, cPart %s" % \
                  (p.name, p.nChar, p.dataType, p.cPart)

        print "    There are %i taxa" % self.nTax

        if len(self.alignments) == 1:
            print "    There is 1 alignment"
        else:
            print "    There are %i alignments" % len(self.alignments)


        if self.cData:
            print "    The cData is %s" % self.cData
        else:
            print "    There is no cData"

        if self.unconstrainedLogLikelihood:
            print "    The unconstrainedLogLikelihood is %s" % self.unconstrainedLogLikelihood
        else:
            pass

    
    def _fill(self, alignments):
        
        # Fill self with Parts from all alignments.
        #
        # This method is called from __init__(), and it is generally
        # not needed on its own.  If we get here, we can be fairly sure
        # that arg alignments is a non-empty list of Alignment
        # objects. This method calls the Alignment method _initParts()

        gm = ["Data._fill()"]
        
        self.alignments = alignments

        # Make a part out of the first alignment.
        if not len(self.alignments):
            gm.append("There are no alignments")
            raise Glitch, gm
        a = self.alignments[0]
        if var.doDataPart:
            a.initDataParts()
        else:
            a._initParts()
        if not len(a.parts):
            gm.append("First alignment failed to make a part")
            raise Glitch, gm
        self.taxNames = a.taxNames
        self.nTax = len(self.taxNames)
        for p in a.parts:
            self.parts.append(p)
        self.nParts = len(self.parts)

        # Now do subsequent alignments ...
        for aligNum in range(len(self.alignments))[1:]:
            a = self.alignments[aligNum]
            if self.nTax != len(a.sequences):
                gm.append("Additional alignment is not the same size as the first alignment.")
                if a.fName:
                    gm.append('(New alignment from file %s.)' % a.fName)
                gm.append("From the first alignment, nTax is %s." % self.nTax)
                gm.append("However, (zero-based) alignment %i has %i sequences." % (aligNum, len(a.sequences)))
                raise Glitch, gm
            if self.nTax != len(a.taxNames):
                gm.append("Additional alignment appears to be not the same size as the first alignment.")
                if a.fName:
                    gm.append('(New alignment from file %s.)' % a.fName)
                gm.append("From the first alignment, nTax is %s." % self.nTax)
                gm.append("However, (zero-based) alignment %i has %i taxNames." % (aligNum, len(a.taxNames)))
                raise Glitch, gm
            for i in range(self.nTax):
                if self.taxNames[i] != a.taxNames[i]:
                    gm.append("Name mis-match in (zero-based) taxon number %i," % i)
                    gm.append("in (zero-based) alignment %i." % aligNum)
                    if a.fName:
                        gm.append('(New alignment from file %s.)' % a.fName)
                    gm.append("Newly-added alignment taxname %s is not the" % a.taxNames[i])
                    gm.append("    same as first alignment taxname %s" % self.taxNames[i])
                    raise Glitch, gm
            if var.doDataPart:
                a.initDataParts()
            else:
                a._initParts()
            if not len(a.parts):
                gm.append("Additional alignment failed to make a part.")
                if a.fName:
                    gm.append('(New alignment from file %s.)' % a.fName)
                raise Glitch, gm
            for p in a.parts:
                self.parts.append(p)
            self.nParts = len(self.parts)




    def calcUnconstrainedLogLikelihood1(self):
        """Calculate likelihood under the multinomial model.

        This calculates the unconstrained (multinomial) log like
        without regard to character partitions.  The result is placed
        in the data variable unconstrainedLogLikelihood.  If there is
        more than one partition, it makes a new temporary alignment
        and puts all the sequences in one part in that alignment.  So
        it ultimately only works on one data partition.  If there is
        more than one alignment, there is possibly more than one
        datatype, and so this method will refuse to do it.  Note that
        the unconstrained log like of the combined data is not the sum
        of the unconstrained log likes of the separate partitions.

        See also calcUnconstrainedLogLikelihood2

        """

        if len(self.alignments) > 1:
            gm = ["Data.calcUnconstrainedLogLikelihood()"]
            gm.append("This method is not implemented for more than one alignment.")
            raise Glitch, gm
        if self.nParts == 1:  # no problem
            self.unconstrainedLogLikelihood = pf.getUnconstrainedLogLike(self.parts[0].cPart)
        else:
            a = self.alignments[0]
            import copy
            newAlig = Alignment()
            newAlig.dataType = a.dataType
            newAlig.symbols = a.symbols
            newAlig.dim = a.dim
            newAlig.equates = a.equates
            newAlig.taxNames = a.taxNames
            for s in a.sequences:
                newAlig.sequences.append(copy.deepcopy(s))
            newAlig.checkLengthsAndTypes()
            newAlig._initParts()
            #newAlig.dump()
            self.unconstrainedLogLikelihood = pf.getUnconstrainedLogLike(newAlig.parts[0].cPart)
            del(newAlig)


    def calcUnconstrainedLogLikelihood2(self):
        """Calculate likelihood under the multinomial model.

        This calculates the unconstrained log like of each data
        partition and places the sum in the Data (self) variable
        unconstrainedLogLikelihood.  Note that the unconstrained log
        like of the combined data is not the sum of the unconstrained
        log likes of the separate partitions.  See also
        calcUnconstrainedLogLikelihood1

        """
        uncon = 0.0
        for p in self.parts:
            #print "            %i    %f" % (p.cPart, pf.getUnconstrainedLogLike(p.cPart))
            uncon = uncon + pf.getUnconstrainedLogLike(p.cPart)
        self.unconstrainedLogLikelihood = uncon


    def _setCStuff(self):
        if self.cData:
            gm = ["Data._setCStuff()"]
            gm.append("This should only be called if self.cData does not exist!")
            raise Glitch, gm
        else:
            if var.doDataPart:
                self.cData = pf.dp_newData(self.nTax, self.nParts)
                for i in range(self.nParts):
                    p = self.parts[i]
                    pf.dp_pokeDataPartInData(p.cPart, self.cData, i)
            else:
                self.cData = pf.newData(self.nTax, self.nParts)
                for i in range(self.nParts):
                    p = self.parts[i]
                    pf.pokePartInData(p.cPart, self.cData, i)
            #print "Made Data.cData = %s" % self.cData



    def writeNexus(self, fName=None, writeDataBlock=0, interleave=0, flat=0, append=0):
        """Write all the alignments in self to a Nexus file.

        If writeDataBlock=1, then taxa and characters are written to a
        'data' block, rather than the default, which is to write
        separate 'taxa' and 'characters' blocks.

        Arg 'flat' gives sequences all on one line.
        Arg 'append', if 0, writes #NEXUS first.  If 1, does not write #NEXUS.
        """

        # There may be more than one alignment, and we need to do the first
        # one first, because it may or may not be appended, while the remaining
        # alignments are appended for sure.
        if len(self.alignments):
            a = self.alignments[0]
            #if a.parts and len(a.parts):
            #    a.resetSequencesFromParts()        # simulate should be responsible for this
            a.writeNexus(fName, writeDataBlock, interleave, flat, append)
            for a in self.alignments[1:]:
                #if a.parts and len(a.parts):
                #    a.resetSequencesFromParts()
                a.writeNexus(fName, writeDataBlock, interleave, flat, append=1)


    def resetSequencesFromParts(self):
        for a in self.alignments:
            if a.parts:
                a.resetSequencesFromParts()
            else:
                raise Glitch, "Alignment has no parts."



    def compoSummary(self):
        """A verbose composition summary, one for each data partition."""

        print "\n\nData composition summary"
        print "========================\n"

        # Make a name format (eg '%12s') that is long enough for the longest name
        longestNameLen = 7 # to start
        for i in self.taxNames:
            if len(i) > longestNameLen:
                longestNameLen = len(i)
        nameFormat = '%' + '%i' % (longestNameLen + 1) + 's'

        for i in range(len(self.parts)):
            p = self.parts[i]
            print "Part %i" % i
            print "%s" % (' ' * (longestNameLen + 1)),
            for j in range(len(p.symbols)):
                print "%10s" % p.symbols[j],
            print "%10s" % 'nSites'
            #print ''
            #cumulativeComps = [0.0] * len(p.symbols)
            grandTotalNSites = 0
            for k in range(p.nTax):
                c = p.composition([k])
                #print "tax %s, part.composition() returns %s" % (k, c)
                nSites = pf.partSequenceSitesCount(p.cPart, k)
                grandTotalNSites = grandTotalNSites + nSites
                print nameFormat % self.taxNames[k],

                # Usually sum(c) will be 1.0, unless the sequence is
                # empty.  We don't want to test "if sum(c) == 0.0:" or
                # "if sum(c):" cuz of small numbers.
                if sum(c) > 0.99:
                    for j in range(len(p.symbols)):
                        print "%10.4f" % c[j],
                        #cumulativeComps[j] = cumulativeComps[j] + (c[j] * nSites)
                else: # Empty sequence, all zeros.  Write dashes.
                    for j in range(len(p.symbols)):
                        print "%10s" % '-',
                print "%10s" % nSites
            c = p.composition()
            print nameFormat % 'mean',
            for j in range(len(p.symbols)):
                print "%10.4f" % c[j],
            #print "%10s" % grandTotalNSites
            print "%10.4f" % (float(grandTotalNSites)/self.nTax)
            print "\n"



    def compoChiSquaredTest(self, verbose=1, skipColumnZeros=0, useConstantSites=1, skipTaxNums=None, getRows=0):
        """A chi square composition test for each data partition.

        So you could do, for example::

            read('myData.nex')

            # Calling Data() with no args tells it to make a Data object 
            # using all the alignments in var.alignments
            d = Data()

            # Do the test.  By default it is verbose, and prints results.
            # Additionally, a list of lists is returned
            ret = d.compoChiSquaredTest()

            # With verbose on, it might print something like ---
            # Part 0: Chi-square = 145.435278, (dof=170) P = 0.913995
            
            print ret
            # The list of lists that it returns might be something like ---
            # [[145.43527849758556, 170, 0.91399521077908041]]
            # which has the same numbers as above, with one 
            # inner list for each data partition.

        If your data has more than one partition::

            read('first.nex')
            read('second.nex')
            d = Data()
            d.compoChiSquaredTest()

            # Output something like ---
            # Part 0: Chi-square = 200.870463, (dof=48) P = 0.000000
            # Part 1: Chi-square = 57.794704, (dof=80) P = 0.971059
            # [[200.87046313430443, 48, 0.0], [57.794704451018163, 80, 0.97105866938683427]]

        where the last line is returned.  With *verbose* turned off,
        the ``Part N`` lines are not printed.

        This method returns a list of lists, one for each data
        partition.  If *getRows* is off, the default, then it is a
        list of 3-item lists, and if *getRows* is turned on then it is
        a list of 4-item lists.  In each inner list, the first is the
        X-squared statistic, the second is the degrees of freedom, and
        the third is the probability from chi-squared.  (The expected
        comes from the data.)  If *getRows* is turned on, the 4th item
        is a list of X-sq contributions from individual rows (ie
        individual taxa), that together sum to the X-sq for the whole
        partition as found in the first item.  This latter way is the
        way that Tree-Puzzle does it.

        Note that this ostensibly tests whether the data are
        homogeneous in composition, but it does not work on sequences
        that are related.  That is, testing whether the X^2 stat is
        significant using the chi^2 curve has a high probability of
        type II error for phylogenetic sequences.

        However, the X-squared stat can be used in valid ways.  You
        can simulate data under the tree and model, and so generate a
        valid null distribution of X^2 values from the simulations, by
        which to assess the significance of the original X^2.  You can
        use this method to generate X^2 values.

        A problem arises when a composition of a character is zero.
        If that happens, we can't calculate X-squared because there
        will be a division by zero.  If *skipColumnZeros* is set to 1,
        then those columns are simply skipped.  They are silently
        skipped unless verbose is turned on.

        So lets say that your original data have all characters, but
        one of them has a very low value.  That is reflected in the
        model, and when you do simulations based on the model you
        occasionally get zeros for that character.  Here it is up to
        you: you could say that the the data containing the zeros are
        validly part of the possibilities and so should be included,
        or you could say that the data containing the zeros are not
        valid and should be excluded.  You choose between these by
        setting *skipColumnZeros*.  Note that if you do not set
        *skipColumnZeros*, and then you analyse a partition that has
        column zeros, the result is None for that partition.

        Another problem occurs when a partition is completely missing
        a sequence.  Of course that sequence does not contribute to
        the stat.  However, in any simulations that you might do, that
        sequence *will* be there, and *will* contribute to the stat.
        So you will want to skip that sequence when you do your calcs
        from the simulation.  You can do that with the *skipTaxNums*
        arg, which is a list of lists.  The outer list is nParts long,
        and each inner list is a list of taxNums to exclude.

        """

        if not useConstantSites:
            newData = Data([])
            aligs = []
            for a in self.alignments:
                #aligs.append(a.removeConstantSites())
                aligs.append(a.subsetUsingMask(a.constantMask(), theMaskChar='1', inverse=1))
            newData._fill(aligs)
            theResult = newData.compoChiSquaredTest(verbose=verbose,
                                                    skipColumnZeros=skipColumnZeros,
                                                    useConstantSites=1, skipTaxNums=skipTaxNums,
                                                    getRows=getRows)
            del(newData)
            return theResult

        gm = ['Data.compoChiSquaredTest()']
        nColumnZeros = 0
        results = []

        # check skipTaxNums
        if skipTaxNums != None:
            if type(skipTaxNums) != type([]):
                gm.append("skipTaxNums should be a list of lists.")
                raise Glitch, gm
            if len(skipTaxNums) != self.nParts:
                gm.append("skipTaxNums should be a list of lists, nParts long.")
                raise Glitch, gm
            for s in skipTaxNums:
                if type(s) != type([]):
                    gm.append("skipTaxNums should be a list of lists.")
                    raise Glitch, gm
                for i in s:
                    if type(i) != type(1):
                        gm.append("skipTaxNums inner list items should be tax numbers.")
                        gm.append("Got %s" % i)
                        raise Glitch, gm

        # Check for blank sequences.  Its a pain to force the user to do this.
        hasBlanks = False
        blankSeqNums = []
        for partNum in range(self.nParts):
            p = self.parts[partNum]
            partBlankSeqNums = []
            for taxNum in range(self.nTax):
                if skipTaxNums and skipTaxNums[partNum] and taxNum in skipTaxNums[partNum]:
                    pass
                else:
                    nSites = pf.partSequenceSitesCount(p.cPart, taxNum) # no gaps, no missings
                    if not nSites:
                        partBlankSeqNums.append(taxNum)
            if partBlankSeqNums:
                hasBlanks = True
            blankSeqNums.append(partBlankSeqNums)
        if hasBlanks:
            gm.append("These sequence numbers were found to be blank. They should be excluded.")
            gm.append("%s" % blankSeqNums)
            gm.append("Set the arg skipTaxNums to this list.")
            raise Glitch, gm
                      

        for partNum in range(self.nParts):
            gm = ['Data.compoChiSquaredTest()  Part %i' % partNum]
            p = self.parts[partNum]
            comps = []
            for taxNum in range(self.nTax):
                if skipTaxNums and skipTaxNums[partNum] and taxNum in skipTaxNums[partNum]:
                    pass
                else:
                    oneComp = p.composition([taxNum])
                    nSites = pf.partSequenceSitesCount(p.cPart, taxNum) # no gaps, no missings
                    #print "tax %i, nSites=%i, oneComp=%s" % (taxNum, nSites, oneComp)
                    if nSites:
                        for k in range(len(oneComp)):
                            oneComp[k] = oneComp[k] * nSites
                        comps.append(oneComp)
                    else:
                        gm.append("(Zero-based) sequence %i is blank, and should be excluded." % taxNum)
                        gm.append("You need to add the number %i to the arg skipTaxNums list of lists." % taxNum)
                        gm.append("(I could do that automatically, but it is best if *you* do it, explicitly.)")
                        gm.append("You can use the Alignment method checkForBlankSequences(listSeqNumsOfBlanks=True)")
                        gm.append("to help you get those inner lists.")
                        raise Glitch, gm
            #print "comps=", comps


            # Here we calculate the X^2 stat.  But we want to check
            # for columns summing to zero.  So we can't use
            # func.xSquared()
            nRows = len(comps)
            nCols = len(comps[0])
            theSumOfRows = func._sumOfRows(comps) # I could have just kept nSites, above
            theSumOfCols = func._sumOfColumns(comps)
            #print theSumOfCols
            isOk = 1
            columnZeros = []
            for j in range(len(theSumOfRows)):
                if theSumOfRows[j] == 0.0:
                    gm.append("Zero in a row sum.  Programming error.")
                    raise Glitch, gm
            for j in range(len(theSumOfCols)):
                if theSumOfCols[j] == 0.0:
                    if skipColumnZeros:
                        columnZeros.append(j)
                    else:
                        if verbose:
                            print gm[0]
                            print "    Zero in a column sum."
                            print "    And skipColumnZeros is not set, so I am refusing to do it at all."
                        isOk = 0
                        nColumnZeros += 1

            theExpected = func._expected(theSumOfRows, theSumOfCols)
            #print "theExpected = ", theExpected
            #print "columnZeros = ", columnZeros
            if isOk:
                if getRows:
                    xSq_rows = []
                xSq = 0.0
                alreadyGivenZeroWarning = 0
                k = 0
                for taxNum in range(self.nTax):
                    if skipTaxNums and skipTaxNums[partNum] and taxNum in skipTaxNums[partNum]:
                        if getRows:
                            xSq_rows.append(0.0)  # this taxon is not in comps.  Add a placeholder
                    else:  # k is the counter for comps and theExpected, taxNum without the skips
                        xSq_row = 0.0
                        for j in range(nCols):
                            if j in columnZeros:
                                if skipColumnZeros:
                                    if verbose and not alreadyGivenZeroWarning:
                                        print gm[0]
                                        print "    Skipping (zero-based) column number(s) %s, which sum to zero." % columnZeros
                                        alreadyGivenZeroWarning = 1
                                else:
                                    gm.append("Programming error.")
                                    raise Glitch, gm
                            else:
                                theDiff = comps[k][j] - theExpected[k][j]
                                xSq_row += (theDiff * theDiff) / theExpected[k][j]
                        xSq += xSq_row
                        if getRows:
                            xSq_rows.append(xSq_row)
                        k += 1
                #print xSq_rows
                dof = (p.dim - len(columnZeros) - 1) * (len(comps) - 1)
                prob = pf.chiSquaredProb(xSq, dof)
                if verbose:
                    print "Part %i: Chi-square = %f, (dof=%i) P = %f" % (partNum, xSq, dof, prob)
                    if getRows:
                        #print "        rows = %s" % xSq_rows
                        print "%20s  %7s  %s" % ('taxName', 'xSq_row', 'P (like puzzle)')
                        for tNum in range(self.nTax):
                            if not skipTaxNums or tNum not in skipTaxNums[partNum]:
                                thisProb = pf.chiSquaredProb(xSq_rows[tNum], self.parts[partNum].dim - 1)
                                print "%20s  %7.5f  %7.5f" % (self.taxNames[tNum], xSq_rows[tNum], thisProb)
                            else:
                                print "%20s    ---      ---" % self.taxNames[tNum]
                if getRows:
                    results.append([xSq, dof, prob, xSq_rows])
                else:
                    results.append([xSq, dof, prob])
            else: # ie not isOk, ie there is a zero in a column sum
                results.append(None) # Maybe a bad idea.  Maybe it should just die, above.
        if nColumnZeros and verbose:
            print "There were %i column zeros." % nColumnZeros
        return results


    def simpleBigXSquared(self):
        """No frills calculation of bigXSquared.

        As in :meth:`Data.Data.compoChiSquaredTest`, but with no
        options, and hopefully faster.  It can't handle gaps or
        ambiguities.  It should be ok for simulations.  It returns a
        list of bigXSquared numbers, one for each data partition.

        If a character happens to not be there, then a column will be
        zero, and so it can't be calculated.  In that case -1.0 is
        returned for that part.
        """
        l = []
        for p in self.parts:
            l.append(pf.partBigXSquared(p.cPart))
        return l
        
    def simpleConstantSitesCount(self):
        """No frills constant sites count.

        It can't handle gaps or ambiguities.  It should be ok
        for simulations.  It returns a list of constant sites counts,
        one for each data partition.

        For each part, of the sites that are not all gaps+ambigs, if
        the sites that are not gaps or ambigs are all the same, then
        it is considered here to be a constant site.

        """
        l = []
        for p in self.parts:
            l.append(pf.partSimpleConstantSitesCount(p.cPart))
        return l
        

                
                
            


    def dupe(self):
        """Copy, making new cParts."""
        
        import copy
        aligListCopy = copy.deepcopy(self.alignments)
        for alig in aligListCopy:
            # We do not want the cPart's, but neither do we want to free the originals.
            for p in alig.parts:
                p.cPart = None
            del(alig.parts)
            alig.parts = []

        return  Data(aligListCopy)



    def bootstrap(self, seed=None):
        """Returns a new data object, filled with bootstrapped data.

        It is a non-parametric bootstrap.  Data partitions are handled
        properly, that is if your data has a charpartition, the
        bootstrap has the same charpartition, and sites are sampled
        only from the appropriate charpartition subset.  """

        gm = ['Data.bootstrap()']

        import copy
        aligListCopy = copy.deepcopy(self.alignments)
        for alig in aligListCopy:
            # We do not want the cPart's, but neither do we want to free the originals.
            for p in alig.parts:
                p.cPart = None
            del(alig.parts)
            alig.parts = []

        d = Data([])
        d._fill(aligListCopy)

        if not self.cData:
            self._setCStuff()
        d._setCStuff()

        if 0:
            print "\nSELF\n===="
            self.dump()
            print "\n\nNEW DATA\n========"
            d.dump()
            raise Glitch

        isNewGSL_RNG = 0
        if not var.gsl_rng:
            var.gsl_rng = pf.get_gsl_rng()
            isNewGSL_RNG = 1
            #print "got var.gsl_rng = %i" % var.gsl_rng

        # Set the GSL random number generator seed, only if it is a new GSL_RNG
        if isNewGSL_RNG:
            if seed != None:
                try:
                    newSeed = int(seed)
                    pf.gsl_rng_set(var.gsl_rng, newSeed)
                except ValueError:
                    print gm[0]
                    print "    The seed should be convertable to an integer"
                    print "    Using the process id instead."
                    pf.gsl_rng_set(var.gsl_rng,  os.getpid())
            else:
                pf.gsl_rng_set(var.gsl_rng,  os.getpid())

        pf.bootstrapData(self.cData, d.cData, var.gsl_rng)

        # Data.resetSequencesFromParts() uses
        # Alignment.resetSequencesFromParts(), which uses
        # partSeq = pf.symbolSequences(self.parts[i].cPart)
        # which uses thePart->sequences

        d.resetSequencesFromParts()
        return d


    def meanNCharsPerSite(self):
        """Mean number of different characters per site, of variable sites only.

        Constant sites are ignored.  Ambiguities and gaps are ignored.

        This is implemented in C, allowing multiple parts.  It is also
        implemented in pure Python in the Alignment class, for single
        parts (which also optionally gives you a distribution in
        addition to the mean); see
        :meth:`Alignment.Alignment.meanNCharsPerSite`.
        
        """
        l = []
        for p in self.parts:
            l.append(pf.partMeanNCharsPerSite(p.cPart))
        return l
    
