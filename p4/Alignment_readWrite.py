import string,sys,os
import func
from Var import var
from SequenceList import Sequence
from Glitch import Glitch


def readOpenPhylipFile(self, flob, nTax, nChar):
    """Read flob to get data in phylip format.

    The user would generally not need to call this method directly.
    It is called by read() etc.  """
    
    gm = ['Alignment.readOpenPhylipFile()']

    dbug = False

    if dbug:
        print "\nreadOpenPhylipFile here"
        if hasattr(flob, 'name'):
            print "    fileName is %s" % flob.name
        print "    nTax is", nTax
        print "    nChar is", nChar
        print "    var.phylipDataMaxNameLength is", var.phylipDataMaxNameLength

    # It is difficult to tell whether it is sequential or interleaved.
    # For example, is the following alignment sequential or interleaved?

    ##          2 14
    ##        acctg     aaaa
    ##        gaattc    aaaa
    ##        gaattc    cccc
    ##        cccccc    cccc

    ##        It might be sequential:
    ##            acctg    aaaagaattcaaaa
    ##           gaattc    cccccccccccccc

    ##        or it might be interleaved:
    ##            acctg    aaaagaattccccc
    ##           gaattc    aaaacccccccccc

    ##        And it is impossible to tell!

    # Here is the top of a file that I had to read.

    ##     36  209
    ##    Dvi        flnsfnakleqpvrqhlknvyaclamstmsaalgaaagflsaigalvfff
    ##    Gsp        finsfnskleqpvrqhlknvyacltmatmaaavgasagflsgigalvffg
    ##    Tca        flnsfsnsleapvrqhlknvyaclamstmaaaigasagflsgigaliffg
    ##    Ban        finsfqnrlespvrqhlknvygtlmmtcgaasagvyagilsaiagaalml
    ##    Bmo        fvnsfqnrleppvrqhlknvyatlmmtcvsasagvyagflsaivgaglml

    # After reading in Dvi, and then the next 3 lines, I had a
    # sequence name (Dvi) and exactly 209 characters, all of which
    # were aa characters.  So it was decided that it was Sequential!  Wrong!

    # Slurp the file into a list of lines, except for the header line.
    # Skip blank lines.
    # If a line startswith '//', stop (its a paml file)
    # If a line is composed wholly of numerals, don't collect it (paml again)
    # If a line starts with  '[', stop collecting until after a line that starts with ']' (paml again)

    # Sometimes paml files have stuff at the end, without any marking.
    # So later we need to stop collecting sequences after nTax.

    # First, get all the lines in a list.  The headerLine is separate.
    # Blank lines and lines that have no potential characters are ignored.
    myPotentialChars = string.letters + '?-'
    
    flob.seek(0, 0)
    headerLine = None
    theLines = []
    inComment = False
    while 1:
        aLine = flob.readline()
        if not aLine:
            break
        aLine = aLine.strip()
        if not aLine: # blank line
            pass
        elif aLine.startswith("//"):
            break
        else:
            if not headerLine:
                headerLine = aLine
            else:
                if aLine.startswith('['):
                    inComment = True
                elif inComment:
                    if aLine.startswith(']'):
                        inComment = False
                else:
                    # Check if it has letters.
                    hasLetters = False
                    for c in aLine:
                        if c in myPotentialChars:
                            hasLetters = True
                            break
                    if hasLetters:
                        theLines.append(aLine)
    assert theLines

    # Now we have the headerLine separate, and all the non-blank lines
    # in theLines list.  Check the headerLine.
    # If a paml file is interleaved, it seems to have an I in the header line.
    splFirstLine = headerLine.split()
    assert len(splFirstLine) >= 2   # More than 2 sometimes
    try:
        firstInt = int(splFirstLine[0])
        secondInt = int(splFirstLine[1])
    except ValueError:  # This should have been caught before, so this should never happen.
        gm.append('bad first line %s' % headerLine)
        raise Glitch, gm
    assert firstInt == nTax
    assert secondInt == nChar

    # Look for an I as the third symbol in the header line.
    # Not used.
    #gotPamlInterleavedI = False
    #if len(splFirstLine) > 2:
    #    thirdBit = splFirstLine[2]
    #    if thirdBit == 'I':
    #        gotPamlInterleavedI = True

    # There are 2 conflicting strategies.  The strict phylip file has
    # a set maximum tax name size, and all names must fit in that size
    # -- usually 10 (but it is a compile-time option).  If the tax
    # name is shorter, it must make up for it with blanks.  And there
    # is no requirement for a space between the tax name and the
    # sequence.  The tax name may have spaces.

    # The other strategy sets off the tax name from the sequence with
    # a space, and so that space is needed.  That implies, I think,
    # that the taxname may not have spaces.  The taxname is often
    # longer than 10.  There is no requirement for the taxname to fill
    # a certain width, and so short names do not need to be filled
    # with multiple blanks -- a single blank will do for both long and
    # short tax names.

    # I am calling the former 'strict' and the latter 'whitespaceSeparatesNames'

    # First I make a half-assed attempt to determine if it is
    # sequential or interleaved.  The decision could easily be wrong,
    # and thats ok.
    
    # The first line might be on a name by itself, which implies
    # sequential.  If it has no spaces, then that will be assumed.

    # The first line might be something like
    # SomeLongName    acgtacgtacgt    
    # which probably means interleaved, with a space.  Not strict
    # phylip format -- the name is longer than 10.

    # The first line might (rarely) be
    # SomeNameABacgtacgtacgtacgt
    # which probably means interleaved, proper phylip format

    # Or the first line might be
    # Short Nameacgtacgtacgt
    # which again implies proper phylip format.
    
    isSequential = None # undecided to start.
    moduloRemainderIsZero = None
    whitespaceSeparatesNames = True 

    # Check whether the number of lines is some multiple of nTax -- if
    # so then it is probably interleaved.  If not, it cannot be interleaved.
    remainder = len(theLines) % nTax
    if remainder == 0:
        moduloRemainderIsZero = True
    else:
        moduloRemainderIsZero = False
        isSequential = True

    if dbug:
        print "    moduloRemainderIsZero is %s" % moduloRemainderIsZero
            
    
    if isSequential == None:
        # Look at the first line of stuff, after the numbers.  
        # If the total length of the line is var.phylipDataMaxNameLength
        # (usually 10) or less, then I assume it is all name.
        firstLine = theLines[0]
        if len(firstLine) <= (var.phylipDataMaxNameLength):
            if dbug:
                print "Got the name %s by virtue of it being a short name on the first data line." % firstLine
                print "Setting 'isSequential' to True, based on that."
            isSequential = True
        elif not firstLine.count(' '):
            isSequential = True
            

    if isSequential == None:
        isSequential = True
    haveTriedSequential_strict = False
    haveTriedSequential = False
    haveTriedInterleaved_strict = False
    haveTriedInterleaved = False
    sequentialResult = None    # (None, True, or False, depending on don't know, success, or failure)
    interleavedResult = None
    gotIt = False

    ###############################################################################
    ###############################################################################
    while not gotIt:
        if dbug:
            if hasattr(flob, 'name'):
                print "  fileName is %s" % flob.name
            print "  isSequential is %s" % isSequential
            print "  whitespaceSeparatesNames is %s" % whitespaceSeparatesNames
            print "  haveTriedSequential = %s" % haveTriedSequential
            print "  haveTriedSequential_strict = %s" % haveTriedSequential_strict
            print "  seqentialResult is %s" % sequentialResult
            print "  haveTriedInterleaved is %s" % haveTriedInterleaved
            print "  haveTriedInterleaved_strict is %s" % haveTriedInterleaved_strict
            print "  interleavedResult is %s" % interleavedResult

        if 0:
            if len(theLines) >= 5:
                theRange = range(5)
            else:
                theRange = range(len(theLines))
            for lineNum in theRange:
                print theLines[lineNum]
        
        if whitespaceSeparatesNames:
            if isSequential and not haveTriedSequential:
                ret = self._readPhylipSequential(nTax, nChar, theLines)
                #print "a ret = %s" % ret
                if ret:
                    gotIt = True
                else:
                    haveTriedSequential = True
                    self.sequences = []
                    if not haveTriedInterleaved:
                        isSequential = False
                    else:
                        if not haveTriedSequential_strict:
                            whitespaceSeparatesNames = False
                        elif not haveTriedInterleaved_strict:
                            whitespaceSeparatesNames = False
                            isSequential = False
            elif not isSequential and not haveTriedInterleaved:
                ret = self._readPhylipInterleaved(nTax, nChar, theLines)
                #print "b ret = %s" % ret
                if ret:
                    gotIt = True
                else:
                    haveTriedInterleaved = True
                    self.sequences = []
                    if not haveTriedSequential:
                        isSequential = True
                    else:
                        if not haveTriedSequential_strict:
                            whitespaceSeparatesNames = False
                        elif not haveTriedInterleaved_strict:
                            whitespaceSeparatesNames = False
                            isSequential = False
        else:
            if isSequential:
                ret = self._readPhylipSequentialStrict(nTax, nChar, theLines)
                #print "c ret = %s" % ret
                if ret:
                    gotIt = True
                else:
                    haveTriedSequential_strict = True
                    self.sequences = []
                    if not haveTriedInterleaved_strict:
                        isSequential = True
                    else:
                        if not haveTriedSequential:
                            whitespaceSeparatesNames = True
                        elif not haveTriedInterleaved:
                            whitespaceSeparatesNames = True
                            isSequential = False
            else:
                ret = self._readPhylipInterleavedStrict(nTax, nChar, theLines)
                #print "d ret = %s" % ret
                if ret:
                    gotIt = True
                else:
                    haveTriedInterleaved_strict = True
                    self.sequences = []
                    if not haveTriedSequential_strict:
                        isSequential = True
                    else:
                        if not haveTriedSequential:
                            whitespaceSeparatesNames = True
                        elif not haveTriedInterleaved:
                            whitespaceSeparatesNames = True
                            isSequential = False
        #print "x gotIt is now %s" % gotIt
        if gotIt:
            break
        if haveTriedSequential and haveTriedInterleaved and haveTriedSequential_strict and haveTriedInterleaved_strict:
            gm.append("Failed to read the phylip or phylip-like file.")
            if not var.verboseRead:
                gm.append("(For more info, turn var.verboseRead on)")
            raise Glitch, gm
    ###########################################################################################
    ###########################################################################################
    
    # Paml interleaved files sometimes have dots after the first sequence.
    firstSequence = self.sequences[0]
    for s in self.sequences[1:]:
        if s.sequence.count('.'):
            s.sequence = list(s.sequence)
            for seqPos in range(nChar):
                if s.sequence[seqPos] == '.':
                    s.sequence[seqPos] = firstSequence.sequence[seqPos]
            s.sequence = ''.join(s.sequence)

    # Set the dataType
    for s in self.sequences:
        ret = None
        ret = func.isDnaRnaOrProtein(s.sequence) # returns 1,2 or 0, respectively
        if ret == 1:
            s.dataType = 'dna'
            s.symbols = 'acgt'
        elif ret == 0:
            s.dataType = 'protein'
            s.symbols = 'arndcqeghilkmfpstwyv'
        else:
            raise Glitch, "Got rna sequence.  Fix me."
    
    # Having read in all the sequences, check for valid characters.
    if len(self.sequences) > 0:
        bads = 0
        if self.sequences[0].dataType == 'dna':
            for s in self.sequences:
                j = 0
                while j < nChar:
                    if s.sequence[j] not in var.validDnaChars:
                        print "Got bad character '%s' in (zero-based) dna sequence %s " % \
                                                    (s.sequence[j], self.sequences.index(s))
                        print "          sequence name: %s" % s.name
                        print "          at (zero-based) position %s" % j
                        bads = bads + 1
                        if bads > 10:
                            print "...and possibly others"
                            break
                    j = j + 1
                if bads > 10:
                    break
            if bads:
                gm.append("See the list of bad chars above")
                raise Glitch, gm
        elif self.sequences[0].dataType == 'protein':
            for s in self.sequences:
                j = 0
                while j < nChar:
                    if s.sequence[j] not in var.validProteinChars:
                        print "Got bad character '%s' in (zero-based) protein sequence %s " % \
                                                    (s.sequence[j], self.sequences.index(s))
                        print "          sequence name: %s" % s.name
                        print "          at (zero-based) position %s" % j
                        bads = bads + 1
                        if bads > 10:
                            print "...and possibly others"
                            break
                    j = j + 1
                if bads > 10:
                    break
            if bads:
                gm.append("See the list of bad chars above")
                raise Glitch, gm

    #for s in self.sequences:
    #    print s.name
    #    print s.dataType
    #sys.exit()
    

def _readPhylipSequential(self, nTax, nChar, theLines):

    if var.verboseRead:
        print "Attempting to read the phylip file assuming"
        print "  sequential format, and that whitespaceSeparatesNames"

    lineNum = 0
    
    for taxNum in range(nTax):
        theSequenceBits = []
        seqLenSoFar = 0
        aLine = theLines[lineNum]
        lineNum += 1
        splitLine = aLine.split()
        theName = splitLine[0]
        if len(splitLine) > 1:
            for aBit in splitLine[1:]:
                bBit = func.stringZapWhitespaceAndDigits(aBit)
                theSequenceBits.append(bBit)
                seqLenSoFar += len(bBit)
        while seqLenSoFar < nChar:
            try:
                aLine = theLines[lineNum]
                lineNum += 1
            except IndexError:
                if var.verboseRead:
                    print "    Early termination"
                return False
            theSequenceBits.append(func.stringZapWhitespaceAndDigits(aLine))
            seqLenSoFar += len(theSequenceBits[-1])
            

        if seqLenSoFar != nChar:
            if var.verboseRead:
                print "    Did not get exactly nChar."
                print "    Expected: ", seqLenSoFar
                print "    Got: ", nChar
            return False

        s = Sequence()
        s.name = theName
        s.sequence = string.join(theSequenceBits, '')
        s.sequence = string.lower(s.sequence)
        self.sequences.append(s)

    if len(self.sequences) != nTax:
        if var.verboseRead:
            print "    Did not get correct nTax."
        return False
                      
    isBad = False
    for s in self.sequences:
        if len(s.sequence) != nChar:
            isBad = True
            break
    if isBad:
        if var.verboseRead:
            print "    At least one sequence has the wrong nChar."
        return False

    if var.verboseRead:
        print "    Got a candidate alignment."
    return True

def _readPhylipInterleaved(self, nTax, nChar, theLines):

    if var.verboseRead:
        print "Attempting to read the phylip file assuming"
        print "  interleaved format, and that whitespaceSeparatesNames"

    lineNum = 0
    for i in range(nTax):
        s = Sequence()
        aLine = theLines[lineNum]
        lineNum += 1
        splitLine = aLine.split()
        s.name = splitLine[0]
        s.theSequenceBits = []
        s.seqLenSoFar = 0
        if len(splitLine) > 1:
            for aBit in splitLine[1:]:
                bBit = func.stringZapWhitespaceAndDigits(aBit)
                s.theSequenceBits.append(bBit)
                s.seqLenSoFar += len(bBit)
        #print "%15s %s" % (s.name, bBit)
        self.sequences.append(s)
    if len(self.sequences) != nTax:
        if var.verboseRead:
            print "    Did not get exactly nTax sequences"
        return False

    # Subsequent cycles don't have names-- just sequence
    seqNum = 0
    isDone = False
    while not isDone:
        try:
            aLine = theLines[lineNum]
            lineNum += 1
        except IndexError:
            break
        #print "aLine a: %s" % aLine
        segment = func.stringZapWhitespaceAndDigits(aLine)
        if len(segment):
            self.sequences[seqNum].theSequenceBits.append(segment)
            self.sequences[seqNum].seqLenSoFar += len(segment)
            seqNum += 1
            if seqNum == nTax:
                seqNum = 0
                #print "aLine b: %s" % aLine
                #print "=" * 50
                # Finished a cycle.  If we have enough sequence, break
                if self.sequences[seqNum].seqLenSoFar >= nChar:
                    isDone = True
        else:
            print '_readPhylipInterleaved() bad line?: %s' % aLine


    for s in self.sequences:
        s.sequence = string.join(s.theSequenceBits, '')
        s.sequence = string.lower(s.sequence)
        del s.theSequenceBits
        del s.seqLenSoFar

    if len(self.sequences) != nTax:
        if var.verboseRead:
            print "    Did not get correct nTax."
        return False
                      
    isBad = False
    for s in self.sequences:
        #s.write()
        #print len(s.sequence), nChar
        if len(s.sequence) != nChar:
            isBad = True
            break
    if isBad:
        if var.verboseRead:
            print "    At least one sequence has the wrong nChar."
            print "    Expected: ", nChar
            print "    Got: ", len(s.sequence)
        return False

    if var.verboseRead:
        print "    Got a candidate alignment."
    return True
    

def _readPhylipSequentialStrict(self, nTax, nChar, theLines):

    if var.verboseRead:
        print "Attempting to read the phylip file assuming"
        print "  sequential format, with strict phylip format"

    lineNum = 0
    
    for taxNum in range(nTax):
        theSequenceBits = []
        seqLenSoFar = 0
        aLine = theLines[lineNum]
        lineNum += 1
        theName = aLine[:var.phylipDataMaxNameLength].strip()
        aBit = func.stringZapWhitespaceAndDigits(aLine[var.phylipDataMaxNameLength:])
        theSequenceBits.append(aBit)
        seqLenSoFar = len(aBit)
        while seqLenSoFar < nChar:
            try:
                aLine = theLines[lineNum]
                lineNum += 1
            except IndexError:
                if var.verboseRead:
                    print "    Early termination"
                return False
            aBit = func.stringZapWhitespaceAndDigits(aLine)
            theSequenceBits.append(aBit)
            seqLenSoFar += len(aBit)
            

        if seqLenSoFar != nChar:
            if var.verboseRead:
                print "    Did not get exactly nChar."
            return False

        s = Sequence()
        s.name = theName
        s.sequence = string.join(theSequenceBits, '')
        s.sequence = string.lower(s.sequence)
        self.sequences.append(s)

    if len(self.sequences) != nTax:
        if var.verboseRead:
            print "    Did not get correct nTax."
        return False
                      
    isBad = False
    for s in self.sequences:
        if len(s.sequence) != nChar:
            isBad = True
            break
    if isBad:
        if var.verboseRead:
            print "    At least one sequence has the wrong nChar."
        return False

    if var.verboseRead:
        print "    Got a candidate alignment."
    return True

def _readPhylipInterleavedStrict(self, nTax, nChar, theLines):
    if var.verboseRead:
        print "Attempting to read the phylip file assuming"
        print "  interleaved format, with strict phylip format"

    lineNum = 0
    for i in range(nTax):
        s = Sequence()
        aLine = theLines[lineNum]
        lineNum += 1
        
        s.name = aLine[:var.phylipDataMaxNameLength].strip()
        s.theSequenceBits = []
        aBit = func.stringZapWhitespaceAndDigits(aLine[var.phylipDataMaxNameLength:])
        s.theSequenceBits.append(aBit)
        s.seqLenSoFar = len(aBit)
        self.sequences.append(s)
    if len(self.sequences) != nTax:
        if var.verboseRead:
            print "    Did not get exactly nTax sequences"
        return False

    # Subsequent cycles don't have names-- just sequence
    seqNum = 0
    isDone = False
    while not isDone:
        try:
            aLine = theLines[lineNum]
            lineNum += 1
        except IndexError:
            break
        #print "aLine a: %s" % aLine
        segment = func.stringZapWhitespaceAndDigits(aLine)
        if len(segment):
            self.sequences[seqNum].theSequenceBits.append(segment)
            self.sequences[seqNum].seqLenSoFar += len(segment)
            seqNum += 1
            if seqNum == nTax:
                seqNum = 0
                # Finished a cycle.  If we have enough sequence, break
                if self.sequences[seqNum].seqLenSoFar >= nChar:
                    isDone = True
        else:
            print '_readPhylipInterleavedStrict() bad line?: %s' % aLine
            
    for s in self.sequences:
        s.sequence = string.join(s.theSequenceBits, '')
        s.sequence = string.lower(s.sequence)
        del s.theSequenceBits
        del s.seqLenSoFar

    if len(self.sequences) != nTax:
        if var.verboseRead:
            print "    Did not get correct nTax."
        return False
                      
    isBad = False
    for s in self.sequences:
        if len(s.sequence) != nChar:
            isBad = True
            break
    if isBad:
        if var.verboseRead:
            print "    At least one sequence has the wrong nChar."
        return False

    if var.verboseRead:
        print "    Got a candidate alignment."
    return True
    

    



##Ignore
def readOpenClustalwFile(self, flob):

    gm = ['Alignment.readOpenClustalwFile()']

    dbug = 0
    if dbug:
        print "\nreadOpenClustalwFile() here"

    # readline up to the first line of sequence
    aLine = flob.readline()
    if not aLine:
        gm.append("No sequence?")
        raise Glitch, gm
    while len(aLine) <= 1:
        aLine = flob.readline()
        if not aLine:
            gm.append("No sequence?")
            raise Glitch, gm
    if dbug:
        print "a- got aLine: '%s'" % aLine

    # Do the first cycle:
    while len(aLine) > 1 and aLine[0] not in string.whitespace:
        s = Sequence()
        splitLine = string.split(aLine)
        if len(splitLine) != 2:
            gm.append("Odd line:\n%s" % aLine)
            raise Glitch, gm
        s.name = splitLine[0]
        s.temp = []
        s.temp.append(string.lower(string.strip(splitLine[1])))
        if dbug:
            print "b got name: %s" % s.name
            print "b got a line of sequence:\n          %s" % s.temp
        self.sequences.append(s)
        aLine = flob.readline()
        if not aLine:
            break

    #nSeqs = len(self.sequences)

    nSeqs = len(self.sequences)
    if not nSeqs:
        gm.append("No sequences?")
        raise Glitch, gm
    if dbug:
        print "Got %i sequences (after first cycle)" % nSeqs

    # Do subsequent cycles
    while 1:
        while len(aLine) <= 1 or aLine[0] in string.whitespace:
            aLine = flob.readline()
            if not aLine:
                break
        if aLine:
            if dbug:
                print "begin subsequent cycle with aLine:\n    %s" % aLine
            for s in self.sequences:
                if aLine[0] in string.whitespace or len(aLine) <= 1:
                    gm.append("Bad line:\n    %s" % aLine)
                    #raise Glitch, gm
                splitLine = string.split(aLine)
                if len(splitLine) != 2:
                    gm.append("Odd line:\n%s" % aLine)
                    raise Glitch, gm
                if s.name != splitLine[0]:
                    gm.append("Problem: existing name %s does not match new name %s" % (s.name, splitLine[0]))
                    raise Glitch, gm
                s.temp.append(string.lower(string.strip(splitLine[1])))
                if dbug:
                    print "got name: %s" % s.name
                    print "got a line of sequence:\n          %s" % s.temp
                aLine = flob.readline()
                if not aLine:
                    break
        else:
            break

    #sys.exit()

    nChar = len(string.join(self.sequences[0].temp, ''))

    for s in self.sequences:
        if dbug:
            print "got name: %s" % s.name
            print "got sequence:\n          %s" % s.temp
        s.temp = string.join(s.temp, '')
        if len(s.temp) != nChar:
            gm.append("Something is wrong with the sequence length of the clustalw file.")
            gm.append("(zero-based) sequence %i" % self.sequences.index(s))
            gm.append("expected %i, got %i" % (nChar, len(s.temp)))
            raise Glitch, gm

        if func.isDnaRnaOrProtein(s.temp): # returns 1,2 or 0, respectively
            s.sequence = s.temp
            s.dataType = 'dna'
            s.symbols = 'acgt'
        else:
            s.sequence = s.temp
            s.dataType = 'protein'
            s.symbols = 'arndcqeghilkmfpstwyv'
        del(s.temp)

    # Having read in all the sequences, check for valid characters
    if len(self.sequences) > 0:
        bads = 0
        if self.sequences[0].dataType == 'dna':
            for s in self.sequences:
                j = 0
                while j < nChar:
                    if s.sequence[j] not in var.validDnaChars:
                        print "bad character '%s' in (zero-based) dna sequence %s " % \
                                                    (s.sequence[j], self.sequences.index(s))
                        print "          sequence name: %s" % s.name
                        print "          at (zero-based) position %s" % j
                        bads = bads + 1
                        if bads > 10:
                            print "...and possibly others"
                            break
                    j = j + 1
                if bads > 10:
                    break
            if bads:
                raise Glitch, gm
        elif self.sequences[0].dataType == 'protein':
            for s in self.sequences:
                j = 0
                while j < nChar:
                    if s.sequence[j] not in var.validProteinChars:
                        print "bad character '%s' in (zero-based) protein sequence %s " % \
                                                    (s.sequence[j], self.sequences.index(s))
                        print "          sequence name: %s" % s.name
                        print "          at (zero-based) position %s" % j
                        bads = bads + 1
                        if bads > 10:
                            print "...and possibly others"
                            break
                    j = j + 1
                if bads > 10:
                    break
            if bads:
                raise Glitch, gm

    if dbug:
        for t in self.sequences:
            print '%20s  %-30s' % ('name', t.name)
            print '%20s  %-30s' % ('comment', t.comment)
            print '%20s  %-30s' % ('sequence', t.sequence)
            print '%20s  %-30s' % ('type', t.dataType)
            print ''






def writeNexus(self, fName=None, writeDataBlock=0,  interleave=0, flat=0, append=0, userText=''):
    """Write self in Nexus format.

    If writeDataBlock=1, then a data block is written, rather than the
    default, which is to write a taxa and a characters block.

    Flat gives sequences all on one line.
    Append, if 0, writes #NEXUS first.  If 1, does not write #NEXUS.

    userText is anything, eg a comment or another Nexus block, that
    you might want to add.  It goes at the end. """

    gm = ["Alignment.writeNexus()"]
    if not fName:
        fName = sys.stdout
    if fName == sys.stdout:
        f = sys.stdout
        if append:
            pass
        else:
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
                print "Alignment: writeNexusFile() 'append' is requested,"
                print "    but '%s' is not a regular file (maybe it doesn't exist?)." % fName
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

    if not writeDataBlock:
        f.write('begin taxa;\n')
        f.write('  dimensions ntax=%s;\n' % len(self.sequences))
        f.write('  taxlabels')
        for i in range(len(self.sequences)):
            f.write(' %s' % func.nexusFixNameIfQuotesAreNeeded(self.sequences[i].name))
        f.write(';\n')
        f.write('end;\n\n')
    else: # ie writeDataBlock=1
        f.write('begin data;\n')
        if self.excludeDelete:
            if self.length < self.excludeDelete.length:
                f.write('  [%i characters have been excluded]\n' % (self.excludeDelete.length - self.length))
        f.write('  dimensions ntax=%s' % len(self.sequences))
    if not writeDataBlock:
        f.write('begin characters;\n')
        if self.excludeDelete:
            if self.length < self.excludeDelete.length:
                f.write('  [%i characters have been excluded]\n' % (self.excludeDelete.length - self.length))
        f.write('  dimensions')

    f.write(' nChar=%s;\n' % self.length)

    f.write('  format')
    if self.dataType == 'dna':
        f.write(' datatype=dna')
    elif self.dataType == 'protein':
        f.write(' datatype=protein')
    elif self.dataType == 'rna':
        f.write( ' datatype=rna')
    elif self.dataType == 'standard':
        f.write(' datatype=standard')
        f.write(' symbols=\"%s\"' % self.symbols)
    if self.equates:
        # We generally do not want to write the usual set of equates.
        if self.dataType == 'dna' or self.dataType == 'rna':
            usualset = set(list('nrykmswbdhv'))
        elif self.dataType == 'protein':
            usualset = set(['x', 'b', 'z'])
        else:
            usualset = set([]) # no usual equates for standard dataType

        # kset is the currently defined set. Now we want to get ks, the ones to write.
        kset = set(self.equates.keys())
        if usualset == kset:
            ks = []   # write none
        # Does the usualset have items not in the kset?  This would be odd, unexpected, but possible.
        elif usualset.difference(kset):
            ks = list(kset) # write all, as it is so unusual.
        elif kset.difference(usualset):
            ks = list(kset.difference(usualset)) # inefficient calculating it twice ...
        #print " [ks=%s]" % ''.join(ks),

        if ks:
            f.write(' equate=\"')
            for k in ks[:-1]:
                if len(self.equates[k]) == 1:
                    f.write('%s=%s ' % (k, self.equates[k]))
                else:
                    f.write('%s={%s} ' % (k, self.equates[k]))
            if len(self.equates[ks[-1]]) == 1:
                f.write('%s=%s\"' % (ks[-1], self.equates[ks[-1]]))
            else:
                f.write('%s={%s}\"' % (ks[-1], self.equates[ks[-1]]))
    if interleave:
        f.write(' interleave')
    f.write(' gap=-')
    f.write(' missing=?')
    f.write(';\n')
    f.write('  matrix\n')
    if interleave:
        if flat:
            gm.append("'interleave' option does not make sense with 'flat' option.")
            if f != sys.stdout:
                f.close()
            raise Glitch, gm
        else:
            # first, get the length of the longest name
            longest = 0
            for i in range(len(self.sequences)):
                s = self.sequences[i]
                if len(s.name) > longest:
                    longest = len(func.nexusFixNameIfQuotesAreNeeded(s.name))
            #formatString = '    %' + `-longest` + 's '  # boring left-justified
            formatString = '    %' + `longest` + 's '   # cool right-justified
            # print "format string is '%s'" % formatString
            if longest > 10:
                wid = 50
            else:
                wid = 60
            pos = 0
            left = len(self.sequences[0].sequence)
            while left > 0:
                for i in range(len(self.sequences)):
                    s = self.sequences[i]
                    f.write(formatString % func.nexusFixNameIfQuotesAreNeeded(s.name))
                    if left >= wid:
                        f.write('%s\n' % s.sequence[pos: pos + wid])
                    elif left > 0:
                        f.write('%s\n' % s.sequence[pos:])
                pos = pos + wid
                left = left - wid
                if left > 0:
                    f.write('\n')

    if not interleave:
        if flat:
            # first, get the length of the longest name
            longest = 0
            for i in range(len(self.sequences)):
                s = self.sequences[i]
                if len(s.name) > longest:
                    longest = len(func.nexusFixNameIfQuotesAreNeeded(s.name))
            #formatString = '    %' + `-longest` + 's '  # boring left-justified
            formatString = '    %' + `longest` + 's '   # cool right-justified
            # print "format string is '%s'" % formatString
            for i in range(len(self.sequences)):
                s = self.sequences[i]
                f.write(formatString % s.name)
                f.write('%s\n' % s.sequence)
        if not flat:
            wid = 60
            for i in range(len(self.sequences)):
                s = self.sequences[i]
                f.write('    %s\n' % func.nexusFixNameIfQuotesAreNeeded(s.name))
                left = len(s.sequence)
                pos = 0
                while left >= wid:
                    f.write('      %s\n' % s.sequence[pos: pos + wid])
                    pos = pos + wid
                    left = left - wid
                if left > 0:
                    f.write('      %s\n' % s.sequence[pos:])
    f.write('  ;\n')
    f.write('end;\n\n')

    if self.nexusSets:
        from NexusSets import NexusSets
        #print "self.nexusSets = %s" % self.nexusSets
        if self.excludeDelete:
            if self.length < self.excludeDelete.length:
                f.write('[skipping out-of-sync nexus sets block]\n')
            else:
                self.nexusSets.writeNexusToOpenFile(f)
        else:
            self.nexusSets.writeNexusToOpenFile(f)


    if userText:
        f.write(userText)
        f.write("\n")
        
    if f != sys.stdout:
        f.close()


def writePhylip(self, fName=None, interleave=True, whitespaceSeparatesNames=True, flat=False):
    """Write the alignment in Phylip format.
    
    If interleave is turned off, then sequences are
    written sequentially.

    Phylip and phylip-like formats are too varied.  The strict Phylip
    format has a set number of spaces for the taxon name, and there
    may not necessarily be a space between the name and the sequence.
    The name size is commonly 10 spaces, but it need not be -- it is a
    compile-time option in Phylip.

    Other programs, eg phyml and PAML, use a phylip-like format where
    the tax name is set off from the sequence by whitespace.  There is
    no set number of spaces that the sequence needs to occupy, and
    there may not be spaces in the tax name.

    This method used to write strict, real phylip format by default,
    where there is a set number of spaces for the taxon name, and
    where there may not necessarily be a space between the name and
    the sequence.  The name size is commonly 10 spaces, but it need
    not be-- it is set by var.phylipDataMaxNameLength (default 10).

    However, it no longer does that by default, as people use the
    format too loosely.  So now 'whitespaceSeparatesNames' is turned
    on by default.  It accommodates names longer than 10 chars.

    If you want to write strict Phylip format, turn
    'whitespaceSeparatesNames' off.  Note that in that format,
    described `here
    <http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles>`_,
    it is ok to have blanks in the sequence, and so p4 puts a blank at
    the beginning of the sequence when writing this way, even though
    the meaning of 'whitespaceSeparatesNames' would suggest not.  Here
    the main effect of whitespaceSeparatesNames=False is that it is
    assumed that whitespace will not be used to trigger the end of the
    name in the downstream program.  Rather it is assumed that in the
    downstream program the end of the name will triggered by reading
    in var.phylipDataMaxNameLength characters.  If I were to set
    whitespaceSeparatesNames=True when I write (because I expect that
    the downstream program that will read it expects that) then I
    check for internal blank spaces in the names (which would be
    disallowed), and also I then allow more than
    var.phylipDataMaxNameLength characters.  However, in both cases
    there is at least one blank after the name.

    Arg 'flat' puts it all on one line.  This is not compatible with
    'interleave', of course.
    
    """

    gm = ['Alignment.writePhylip(fName=%s, interleave=%s, whitespaceSeparatesNames=%s, flat=%s)' % (
        fName, interleave, whitespaceSeparatesNames, flat)]

    # Find the maxNameLen.
    maxNameLen = 0
    namesHaveSpaces = False
    for s in self.sequences:
        theNameLen = len(s.name)
        if theNameLen > maxNameLen:
            maxNameLen = theNameLen
        if s.name.count(' '):
            namesHaveSpaces = True
    #print 'The longest name length in this alignment is %i' % maxNameLen
            
    if whitespaceSeparatesNames and namesHaveSpaces:
        crimes = []
        for s in self.sequences:
            if s.name.count(' '):
                crimes.append("has space in the name: %s" % s.name)
        gm.append("whitespaceSeparatesNames is set, but some tax names")
        gm.append("have spaces.  That won't work. -- Fix it.")
        for crime in crimes:
            gm.append(crime)
        raise Glitch, gm

    doStrictOkAnyway = False
    if maxNameLen < var.phylipDataMaxNameLength:
        doStrictOkAnyway = True

    if interleave and flat:
        gm.append("Both 'interleave' and 'flat' are turned on -- does not work.")
        raise Glitch, gm
    
    # Check and complain if any taxNames will be truncated.
    if not whitespaceSeparatesNames and maxNameLen > var.phylipDataMaxNameLength:
        gm.append('The longest name length in this alignment is %i' % maxNameLen)
        gm.append('var.phylipDataMaxNameLength is now %i' % var.phylipDataMaxNameLength)
        gm.append('Sequence names will be truncated.  Fix it.')
        gm.append("You may want to use the 'renameForPhylip()' method.")
        raise Glitch, gm

    nameWid = var.phylipDataMaxNameLength + 1
    spacer1 = var.phylipDataMaxNameLength + 1
    if whitespaceSeparatesNames and (maxNameLen >= nameWid):
        nameWid = maxNameLen + 1
        spacer1 = 11
    #print 'The nameWid is %i' % nameWid
    
    
    if fName == None or fName == sys.stdout:
        f = sys.stdout
    else:
        try:
            f = open(fName, 'w')
        except IOError:
            gm.append("Can't open '%s' for writing." % fName)
            raise Glitch, gm
    f.write(' %i  %i\n' % (len(self.sequences), self.length))

    if interleave:
        #wid1 = 50
        if nameWid < 50:
            wid1 = 61 - nameWid
        else:
            wid1 = 0
            
        # do the first row
        if self.length >= wid1:
            upper = wid1
        else:
            upper = self.length
        for k in self.sequences:
            #theFormat = "%-" + "%is" % (var.phylipDataMaxNameLength + 1)
            #f.write(theFormat % k.name[0:var.phylipDataMaxNameLength])
            theFormat = "%-" + "%is" % (nameWid)
            #sys.stdout.write(theFormat % k.name[0:nameWid])
            #sys.stdout.write('%s\n' % k.sequence[0:upper])
            f.write(theFormat % k.name[0:nameWid])
            f.write('%s\n' % k.sequence[0:upper])
        lower = upper
        f.write('\n')
        # do subsequent rows
        wid2 = 50
        while lower < self.length:
            if self.length >= lower + wid2:
                upper = lower + wid2
            else:
                upper = self.length
            for k in self.sequences:
                #f.write('             ')
                #f.write(' ' * (var.phylipDataMaxNameLength + 1))
                f.write(' ' * spacer1)
                f.write('%s\n' % k.sequence[lower:upper])
            lower = upper
            f.write('\n')
    if not interleave:
        if flat:
            #print "nameWid = ", nameWid
            theFormat = "%-" + "%is" % nameWid
            for s in self.sequences:
                f.write(theFormat % s.name[:nameWid])
                f.write('%s\n' % s.sequence)
        else:
            #wid = 50
            if nameWid < 50:
                wid1 = 61 - nameWid
            else:
                wid1 = 0
            wid2 = 50
            theFormat = "%-" + "%is" % nameWid
            for s in self.sequences:
                #theFormat = "%-" + "%is   " % var.phylipDataMaxNameLength
                #f.write(theFormat % s.name[:var.phylipDataMaxNameLength])
                f.write(theFormat % s.name[:nameWid])
                left = len(s.sequence)
                pos = 0
                if left >= wid1:
                    f.write('%s\n' % s.sequence[pos: pos + wid1])
                    pos = pos + wid1
                    left = left - wid1
                else:
                    f.write('%s\n' % s.sequence[pos:])
                while left >= wid2:
                    f.write(' ' * spacer1)
                    f.write('%s\n' % s.sequence[pos: pos + wid2])
                    pos = pos + wid2
                    left = left - wid2
                if left > 0:
                    f.write(' ' * spacer1)
                    f.write('%s\n' % s.sequence[pos:])
                f.write('\n')

    if f != sys.stdout:
        f.close()


def writeMolphy(self, fName=None):
    """Another phylogenetics program, another format.

    Does anybody use this anymore?
    """

    gm = ["Alignment.writeMolphy()"]
    if fName == None or fName == sys.stdout:
        f = sys.stdout
    else:
        try:
            f = open(fName, 'w')
        except IOError:
            gm.append("Can't open %s for writing." % fName)
            raise Glitch, gm
    f.write(' %i  %i\n' % (len(self.sequences), self.length))

    wid = 50
    for i in range(len(self.sequences)):
        s = self.sequences[i]
        f.write('%s\n' % s.name)
        left = len(s.sequence)
        pos = 0
        if left >= wid:
            f.write('%s\n' % s.sequence[pos: pos + wid])
            pos = pos + wid
            left = left - wid
        else:
            f.write('%s\n' % s.sequence[pos:])
        while left >= wid:
            f.write('%s\n' % s.sequence[pos: pos + wid])
            pos = pos + wid
            left = left - wid
        if left > 0:
            f.write('%s\n' % s.sequence[pos:])
        f.write('\n')

    if f != sys.stdout:
        f.close()



##Ignore
def readOpenGdeFile(self, flob, inverseMasks=0):
    dbug = 0
    gm = ['Alignment: readOpenGdeFile()']
    sList = []
    maskList = []
    aLine = flob.readline()
    aLine = string.strip(aLine)
    if dbug:
        print aLine
    if aLine[0] != '{':
        gm.append("The first character is not '{' --- not a GDE file??")
        gm.append("(Note that p4 only reads proper gde files, not gde flat files.)")
        raise Glitch, gm
    while aLine:
        s = Sequence()
        if aLine[0] == '{':
            aLine = string.strip(flob.readline())
        else:
            gm.append("Problem with GDE file at line:\n\t%s" % aLine)
            raise Glitch, gm
        if aLine[:4] == 'name':
            # the name line might be like: name    "D.ethanoge"
            # or it might be like: name    "C.acetobA ", ie with a spurious space.
            if dbug:
                print "got name line: \n%s" % aLine
            splitLine = string.split(aLine, '\"')
            s.name = string.strip(splitLine[1])
            #print s.name
            if not func.nexusCheckName(s.name):
                gm.append("Bad name '%s'" % s.name)
                raise Glitch, gm
        else:
            gm.append("I was expecting a name line: instead I got line:\n\t%s" % aLine)
            raise Glitch, gm
        aLine = string.strip(flob.readline())
        if aLine[:4] == 'type':
            if dbug:
                print "get type line: \n%s" % aLine
            splitLine = string.split(aLine)
            type = splitLine[1][1:-1]
            #print type
            if type not in ['DNA', 'PROTEIN', 'MASK', 'TEXT', 'RNA']:
                gm.append("Unknown type '%s'" % type)
                raise Glitch, gm
        else:
            gm.append("I was expecting a type line: instead I got line:\n\t%s" % aLine)
            raise Glitch, gm

        if type == 'DNA' or type == 'RNA':
            s.dataType = 'dna'
        elif type == 'PROTEIN':
            s.dataType = 'protein'
        elif type == 'MASK':
            s.dataType = 'mask'   # a dodgy manoever, letting a Sequence object hold a mask or text.
        elif type == 'TEXT':
            s.dataType = 'text'

        if dbug:
            print "dataType set to: %s" % s.dataType

        # Next, we are on the lookout for an 'offset' line, which
        # may or may not be there.  If we come across the first
        # line of the sequence, break
        offset = 0
        while 1:
            aLine = flob.readline()
            # print "a Reading line: %s" % aLine
            if aLine[:6] == 'offset':
                aLine = string.strip(aLine)
                if dbug:
                    print "got offset line: \n%s" % aLine
                splitLine = string.split(aLine)
                try:
                    offset = int(splitLine[1])
                    if dbug:
                        print "get offset: %i" % offset
                except ValueError:
                    gm.append("Bad offset in line '%s'" % aLine)
                    raise Glitch, gm
            elif aLine[:8] == 'sequence' and aLine[:11] != 'sequence-ID':
                aLine = string.strip(aLine)
                break

        # Now aLine has the first line of the sequence
        if dbug:
            print "\nGot first line of sequence: %s" % aLine
        # Get the rest of the sequence.
        splitLine = string.split(aLine)
        thisSeqList = []
        thisSeqList.append(splitLine[1])
        aLine = string.strip(flob.readline())
        while aLine[0] != '}':
            thisSeqList.append(aLine)
            aLine = string.strip(flob.readline())
        #print "Got thisSeqList: %s" % thisSeqList
        s.sequence = (offset * '-') + string.join(thisSeqList, '')[1:-1]
        #print "Finished: name:'%s', dataType %s, sequence = '%s'" % (s.name, s.dataType, s.sequence)
        #print "Got last line of group: '%s'" % aLine\
        if s.dataType == 'mask':
            maskList.append(s)
        else:
            sList.append(s)

        while 1:
            aLine = flob.readline()
            if not aLine:
                break
            if aLine[0] == '{':
                break

    # At this point we have read in to the end of the file, and
    # have a sList containing Sequences, and a maskList containing
    # masks.

    # Make all the sequences the same length.  GDE certainly does not do this.
    if len(sList) == 0:
        gm.append("Finished reading file, but got no sequences.  What gives?")
        raise Glitch, gm
    maxSeqLen = 0
    for s in sList:
        if len(s.sequence) > maxSeqLen:
            maxSeqLen = len(s.sequence)
    for s in maskList:
        if len(s.sequence) > maxSeqLen:
            maxSeqLen = len(s.sequence)
    for s in sList:
        if len(s.sequence) != maxSeqLen:
            s.sequence = s.sequence + ((maxSeqLen - len(s.sequence)) * '-')
    for s in maskList:
        if len(s.sequence) != maxSeqLen:
            s.sequence = s.sequence + ((maxSeqLen - len(s.sequence)) * '-')

    # Now we have a bunch of sequences.  Fill in the alignment (ie self)
    self.sequences = sList
    self.checkNamesForDupes()
    if 0:
        self.checkLengthsAndTypes()
        self.nexusSets = NexusSets()
        self.nexusSets.nChar = self.length
        self.nexusSets.setPredefinedCharSets(self)
    else:
        # check lengths
        topLen = len(self.sequences[0].sequence)
        for s in self.sequences:
            #print "xxx dataType=%s" % s.dataType
            if len(s.sequence) != topLen:
                gm.append("Got sequences of unequal lengths")
                raise Glitch,gm
        self.length = topLen

        # check dataType's
        topType = self.sequences[0].dataType
        for s in self.sequences:
            if s.dataType != topType:
                gm.append("Mixed dataTypes.")
                raise Glitch, gm
        self.dataType = topType
        if self.dataType == 'dna':
            self.symbols = 'acgt'
            self.dim = 4
            if not self.equates:
                self.equates = {'n': 'acgt', 'm': 'ac', 'k': 'gt', # 'x': 'acgt', 
                                'h': 'act', 'y': 'ct', 'v': 'acg',
                                'w': 'at', 'd': 'agt', 'b': 'cgt',
                                'r': 'ag', 's': 'cg'}
        elif self.dataType == 'protein':
            self.symbols = 'arndcqeghilkmfpstwyv'
            self.dim = 20
            if not self.equates:
                self.equates = {'b': 'dn', 'x': 'arndcqeghilkmfpstwyv', 'z': 'eq'}
        else:
            gm.append("Can't deal with dataType '%s'." % self.dataType)
            raise Glitch, gm

    #self.dump()
    # Is the following really needed?
    #for s in self.sequences:
    #    if s.name not in self.taxNames:
    #        self.taxNames.append(s.name)
    #    else:
    #        gm.append("Duplicated name '%s'" % s.name)
    #        raise Glitch, gm

    if 1:
        if len(maskList):
            from NexusSets import NexusSets
            self.nexusSets = NexusSets()
        for s in maskList:
            from NexusSets import CharSet
            c = CharSet()
            c.nexusSets = self.nexusSets
            c.nChar = self.length
            c.name = s.name
            c.lowName = string.lower(s.name)
            c.format = 'vector'
            if inverseMasks:
                s.sequence = list(s.sequence)
                for i in range(len(s.sequence)):
                    if s.sequence[i] == '0':
                        s.sequence[i] = '1'
                    else:
                        s.sequence[i] = '0'
                s.sequence = string.join(s.sequence, '')
            c.mask = s.sequence
            for c1 in self.nexusSets.charSets:
                if c1.lowName == c.lowName:
                    gm.append("Duplicated charSet (ie mask) name %s" % c.name)
                    raise Glitch, gm
            self.nexusSets.charSets.append(c)


