import sys,re,string,os,cStringIO
import func
import copy
from Var import var
from Glitch import Glitch
from subprocess import Popen,PIPE

class Sequence(object):
    """A container for a single molecular sequence.

    We have

    - **sequence** a string, the molecular sequence
    - **name**     a string, the name
    - **dataType** either 'dna', or 'protein', or None, meaning 'standard'
    """

    def __init__(self):
        self.name = None
        self.comment = None
        self.sequence = ''
        self.dataType = None

    def _getNChar(self):
        if self.sequence:
            return len(self.sequence)
        else:
            return 0
        
    nChar = property(_getNChar)

    # If the self.sequence len is zero, and I don't have
    # __nonzero__(), then "assert self" will raise an AssertionError,
    # basing that response on the result of len(self).  Having
    # __nonzero__() makes "assert self" work.
    def __nonzero__(self):
        return True
    def __len__(self):
        if self.sequence:
            return len(self.sequence)
        else:
            return 0

    def dump(self):
        """Print rubbish about self."""
        print '%15s: %s' % ('name', self.name)
        if self.comment: print '%15s: %s' % ('comment', self.comment)
        #if self.dataType == 'dna':
            #if self.transl_table:
            #    print "%15s: %s' % ('transl_table", self.transl_table)
        if self.sequence:
            print '%15s: %s' % ('sequence', self.sequence[:25]),
            if len(self.sequence) > 25: print "..."
            else: print ''

    def dupe(self):
        """Return a duplicate of self."""
        
        return copy.deepcopy(self)
    
    def reverseComplement(self):
        """Convert self.sequence, a DNA sequence, to its reverse complement.

        Ambigs are handled correctly.  I think.
        """
        
        assert self.dataType == 'dna'
        self.sequence = list(self.sequence)
        self.sequence.reverse()
        # {'b': 'cgt', 'd': 'agt', 'h': 'act', 'k': 'gt', 'm': 'ac',
        #  'n': 'acgt', 's': 'cg', 'r': 'ag', 'w': 'at', 'v': 'acg', 'y': 'ct'}  # 'x': 'acgt', 
        for i in range(len(self.sequence)):
            c = self.sequence[i]
            if c == 'a':
                self.sequence[i] = 't'
            elif c == 't':
                self.sequence[i] = 'a'
            elif c == 'c':
                self.sequence[i] = 'g'
            elif c == 'g':
                self.sequence[i] = 'c'
            elif c == '-':
                pass
            elif c == 'n':
                pass
            #elif c == 'x':
            #    pass
            elif c == 'r':
                self.sequence[i] = 'y'
            elif c == 'y':
                self.sequence[i] = 'r'

            elif c == 'b':
                self.sequence[i] = 'v'
            elif c == 'd':
                self.sequence[i] = 'h'
            elif c == 'h':
                self.sequence[i] = 'd'
            elif c == 'k':
                self.sequence[i] = 'm'
            elif c == 'm':
                self.sequence[i] = 'k'
            elif c == 's':
                pass
                #self.sequence[i] = 's'
            elif c == 'w':
                pass
                #self.sequence[i] = 'w'
            elif c == 'v':
                self.sequence[i] = 'b'
            else:
                gm = ["Sequence.reverseComplement()"]
                if c in string.uppercase:
                    gm.append("Got uppercase '%s' How did that happen? -- can only handle lowercase." % c)
                else:
                    gm.append("Sequence.reverseComplement().  Got char '%s' What is it?" % c)
                raise Glitch, gm
            
        self.sequence = ''.join(self.sequence)
                
    def writeFastaToOpenFile(self, flob, width=60, doComment=True, writeExtraNewline=True):
        flob.write('>%s' % self.name)
        if doComment and self.comment:
            flob.write(' %s' % self.comment)
        flob.write('\n')
        left = len(self.sequence)
        pos = 0
        while left >= width:
            if var.writeFastaUppercase:
                flob.write('%s\n' % self.sequence[pos: pos + width].upper())
            else:
                flob.write('%s\n' % self.sequence[pos: pos + width])
            pos = pos + width
            left = left - width
        if left > 0:
            if var.writeFastaUppercase:
                flob.write('%s\n' % self.sequence[pos:].upper())
            else:
                flob.write('%s\n' % self.sequence[pos:])
        if writeExtraNewline:
            flob.write('\n')

    def write(self):
        self.writeFastaToOpenFile(sys.stdout)

    def writeFasta(self, fName=None, width=60, doComment=True):
        isFlob = False
        if not fName or fName == sys.stdout:
            f = sys.stdout
            isFlob = True
        elif hasattr(fName, 'write'): # an open file-like object
            f = fName
            isFlob = True
        else:
            f = file(fName, 'w')
        self.writeFastaToOpenFile(f, width=width, doComment=doComment)
        if not isFlob:
            f.close()
        
    def translate(self, transl_table=1, checkStarts=False, nnn_is_gap=False):
        """Returns a protein Sequence from self, a DNA sequence.

        Self is translated using
        :meth:`GeneticCode.GeneticCode.translate`, so it handles
        ambiguities.  At the moment, we can only do translations where the
        frame of the codon is 123, ie the first sequence position is the
        first position of the codon.  The default transl_table is the
        standard (or so-called universal) genetic code, but you can change
        it.

        Other available translation tables, this week::

            if transl_table == 1: # standard
            elif transl_table == 2: # vertebrate mito
            elif transl_table == 4: # Mold, Protozoan,
                                    # and Coelenterate Mitochondrial Code
                                    # and the Mycoplasma/Spiroplasma Code
            elif transl_table == 5: # invertebrate mito
            elif transl_table == 9: # echinoderm mito

            and now 6, 10, 11, 12, 13, 14, 21.

        (These are found in :class:`GeneticCode.GeneticCode`)

        See also :meth:`Alignment.Alignment.checkTranslation` and
        :meth:`Alignment.Alignment.checkTranslation`.

        If the arg *checkStarts* is turned on (by default it is not turned
        on) then this method checks whether the first codon is a start
        codon, and if it is then it uses it.

        Arg *nnn_is_gap* is for odd sequences where there are long
        stretches of 'nnn' codons, which probably should be gaps.
        Probably best to correct those elsewise.

        """

        gm = ['Sequence.translate()']
        if self.dataType != 'dna':
            gm.append("Self should be a DNA Sequence")
            raise Glitch, gm

        if self.nChar % 3 != 0:
            gm.append("The length of self should be a multiple of 3")
            raise Glitch, gm
        nTriplets = self.nChar / 3

        from GeneticCode import GeneticCode
        gc = GeneticCode(transl_table)

        prSeq = Sequence()
        prSeq.dataType = 'protein'
        prSeq.name = self.name
        prSeq.sequence = ['-'] * nTriplets
        
        dnaSeq = self.sequence
        protSeq = prSeq.sequence
        for j in range(nTriplets):
            theCodon = dnaSeq[(j * 3):(j * 3) + 3]
            #print theCodon
            if theCodon == '---':
                protSeq[j] = '-'
            elif theCodon.count('-'):
                print "    seq %s, position %4i, dnaSeq %4i, codon '%s' is incomplete" % (self.name, j, (j*3), theCodon)
            elif theCodon == 'nnn':
                if nnn_is_gap:
                    print "    seq %s, position %4i, dnaSeq %4i, codon '%s' translating to a gap ('-')" % (
                        self.name, j, (j*3), theCodon)
                    protSeq[j] = '-'
                else:
                    protSeq[j] = 'x'
            else:
                protSeq[j] = gc.translate(theCodon)
                if checkStarts and j == 0:
                    if theCodon in gc.startList:
                        sys.stderr.write("    Seq %s. The first codon, '%s', is a start codon -- making it m\n" % (
                            self.name, theCodon))
                        protSeq[j] = 'm'
                    else:
                        sys.stderr.write("    Seq %s. The first codon, '%s', is not a start codon\n" % (
                            self.name, theCodon))

        # Get rid of stop translation '*'
        if prSeq.sequence[-1] == '*':
            prSeq.sequence.pop()
        
        prSeq.sequence = string.join(prSeq.sequence, '')
        return prSeq

        


class SequenceList(object):
    """A container for a list of Sequence objects.

    The usual input would be a fasta file::

        read('sequences.fas')
        sl = var.sequenceLists.pop()

        # see what you have
        sl.dump()

        # look at the sequences
        for s in sl.sequences:
            print s.name, s.dataType

        # Get at sequences by name from a dictionary
        sl.makeSequenceForNameDict()
        s = sl.sequenceForNameDict['mammoth']

        # align them using muscle
        a = sl.muscle()

        
    """
        
    def __init__(self, flob = None):
        self.sequences = []
        self.fName = None
        self.sequenceForNameDict = None
        if flob:
            self._readFastaFile(flob)
            if hasattr(flob, 'name'):
                self.fName = flob.name

    def makeSequenceForNameDict(self):
        self.sequenceForNameDict = {}
        for s in self.sequences:
            assert not self.sequenceForNameDict.has_key(s.name), "duped name %s" % s.name
            self.sequenceForNameDict[s.name] = s

    def _readFastaFile(self, flob):
        flob.seek(0)
        gm = ['SequenceList._readFastaFile()']
        if hasattr(flob, 'name'):
            gm = ['SequenceList._readFastaFile(%s)' % flob.name]

        complaintAboutLength = """
        Lines should not be longer than 120 characters.
        This will be overlooked here, but other programs may gag.
        """
        alreadyComplainedAboutLength = 0

        # The first line might start with a ';'
        # Move the position to the first '>'
        aLine = flob.readline()
        while aLine[0] != '>':
            aLine = flob.readline()

        mySeq = Sequence()
        mySeq.sequence = []  # So I can append lines.  I'll change it back to a string later

        # Note that anything after the comment char ';' is ignored
        aLine = aLine[:string.find(aLine, ';')]

        # Get rid of any readseq-generated checksum and bases, ie >theName, 8 bases, A87261CF checksum.
        #print aLine[-9:]
        if aLine[-9:] == 'checksum.':
            i = string.rfind(aLine, ',') # find the last comma
            if i != -1:                  # -1 if a comma was not found
                aLine = aLine[:i]
                if aLine[-5:] == 'bases':
                    i = string.rfind(aLine, ',')
                    if i != -1:
                        aLine = aLine[:i]
        #print "aLine is now '%s'" % aLine

        # Split the first line into two parts, at the first space, if there is one.
        i = string.find(aLine, ' ')
        if i != -1:
            mySeq.headLineList = [string.strip(aLine[1:i]), string.strip(aLine[i:])]
        else:
            mySeq.headLineList = [string.strip(aLine[1:])]
        #print mySeq.headLineList

        # This may seem backwards, but it works...
        aLine = flob.readline()
        while aLine:
            #print "aLine, before find: %s   %s" % (aLine, aLine[-1])
            # There seems to be a bug with find.  If aLine does not end in '\n',
            # then it gets shortened by one.  Here is a workaround.
            if aLine[-1] == '\n' or aLine[-1] == '\r':
                aLine = aLine[:string.find(aLine, ';')]
            else:
                aLine = aLine + '\n'
                aLine = aLine[:string.find(aLine, ';')]
            #print "aLine, after find: %s" % aLine

            if len(aLine) and aLine[0] == '>':
                self.sequences.append(mySeq)
                mySeq = Sequence()
                mySeq.sequence = []

                # Get rid of any readseq-generated checksum
                #print aLine[-9:]
                if aLine[-9:] == 'checksum.':
                    i = string.rfind(aLine, ',') # find the last comma
                    if i != -1:                  # -1 if a comma was not found
                        aLine = aLine[:i]
                        if aLine[-5:] == 'bases':
                            i = string.rfind(aLine, ',')
                            if i != -1:
                                aLine = aLine[:i]
                i = string.find(aLine, ' ')
                if i != -1:
                    mySeq.headLineList = [string.strip(aLine[1:i]), string.strip(aLine[i:])]
                else:
                    mySeq.headLineList = [string.strip(aLine[1:])]
                #print mySeq.headLineList
                aLine = flob.readline()
            else:
                assert mySeq
                mySeq.sequence.append(string.lower(string.strip(aLine)))
                aLine = flob.readline()
        self.sequences.append(mySeq)

        # now fix the sequences
        toLowerTransTable = string.maketrans(string.uppercase[:26], string.lowercase[:26])
        for mySeq in self.sequences:
            mySeq.sequence = string.joinfields(mySeq.sequence, '')
            #print mySeq.sequence
            if string.count(mySeq.sequence, '.'):
                gm.append("Dots don't work in a fasta file, do they?")
                raise Glitch, gm
            mySeq.sequence = string.translate(mySeq.sequence, toLowerTransTable,
                                          string.digits + string.whitespace + '\0')
            dType = func.isDnaRnaOrProtein(mySeq.sequence)
            if dType == 1:
                #print "Its dna"
                mySeq.dataType = 'dna'
            elif dType == 2:
                #print "Its rna"
                #print "Converting RNA to DNA"
                mySeq.sequence = list(mySeq.sequence)
                for i in range(len(mySeq.sequence)):
                    if mySeq.sequence[i] == 'u':
                        mySeq.sequence[i] = 't'
                mySeq.sequence = string.join(mySeq.sequence, '')
                mySeq.dataType = 'dna'
            else:
                #print "Its protein"
                mySeq.dataType = 'protein'

            if mySeq.headLineList:
                if len(mySeq.headLineList) == 1:
                    mySeq.name = mySeq.headLineList[0]
                elif len(mySeq.headLineList) == 2:
                    mySeq.name = mySeq.headLineList[0]
                    mySeq.comment = mySeq.headLineList[1]
        if 0:
            for mySeq in self.sequences:
                print '%20s  %-30s' % ('name', mySeq.name)
                print '%20s  %-30s' % ('comment', mySeq.comment)
                print '%20s  %-30s' % ('sequence', mySeq.sequence)
                print '%20s  %-30s' % ('dataType', mySeq.dataType)
                print ''
                

        # check for invalid chars
        if len(self.sequences) > 0:
            bads = 0
            if self.sequences[0].dataType == 'dna':
                for s in self.sequences:
                    j = 0
                    while j < len(s.sequence):
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
                    gm.append("Got bad characters.")
                    raise Glitch, gm
            if self.sequences[0].dataType == 'protein':
                for s in self.sequences:
                    j = 0
                    while j < len(s.sequence):
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
                    gm.append("Got bad characters.")
                    raise Glitch, gm
        flob.close()


    def _readOpenPirFile(self, flob):

        # http://www.ebi.ac.uk/help/formats.html
        # NBRF/PIR Format:

        # * A sequence in PIR format consists of:
        #      1. One line starting with
        #            1. a ">" (greater-than) sign, followed by
        #            2. a two-letter code describing the sequence type (P1, F1, DL, DC, RL, RC, or XX), followed by
        #            3. a semicolon, followed by
        #            4. the sequence identification code (the database ID-code).
        #      2. One line containing a textual description of the sequence.
        #      3. One or more lines containing the sequence itself. The end of the sequence is marked
        #          by a "*" (asterisk) character.
        # * A file in PIR format may comprise more than one sequence.

        # Sequence type         Code
        # Protein (complete)    P1
        # Protein (fragment)    F1
        # DNA (linear)          DL
        # DNA (circular)        DC
        # RNA (linear)          RL
        # RNA (circular)        RC
        # tRNA                  N3
        # other functional RNA  N1

        flob.seek(0)
        gm = ['SequenceList._readOpenPirFile()']
        if hasattr(flob, 'name'):
            gm = ['SequenceList._readOpenPirFile(%s)' % flob.name]

        ll = [l.strip() for l in flob.readlines()]
        lNum = -1
        while 1:
            lNum += 1
            try:
                aLine = ll[lNum]
                #print "a %4i aLine: '%s'" % (lNum, aLine)
            except IndexError:
                break
            #print "a1 %4i aLine: '%s'" % (lNum, aLine)
            if not aLine or not aLine.startswith('>'):
                try:
                    lNum += 1
                    aLine = ll[lNum]
                    #print "b %4i aLine: '%s'" % (lNum, aLine)
                except IndexError:
                    break
            #print "c %4i aLine: %s" % (lNum, aLine)
            if aLine[3] != ';':
                gm.append("First line is: %s" % aLine.rstrip())
                gm.append("4th char should be ';'")
                raise Glitch, gm
            twoChars = aLine[1:3]
            if twoChars not in ['P1']:
                gm.append("First line is: %s" % aLine.rstrip())
                gm.append("Code characters '%s' are not recognized / implemented.  Fix me?" % twoChars)
                raise Glitch, gm
            
            seqObj = Sequence()
            if twoChars == 'P1':
                seqObj.dataType = 'protein'
            else:
                gm.append("Pir datatype code '%s' is not implemented.  Fix me." % twoChars)
                raise Glitch, gm
            seqObj.sequence = []  # So I can append lines.  I'll change it back to a string later
            splLine = aLine.split(';')
            seqObj.name = splLine[1]
            #print "got pir seq name %s" % seqObj.name

            # Get the comment line.
            lNum += 1
            try:
                aLine = ll[lNum]
                #print "d %4i aLine: %s" % (lNum, aLine)
                if aLine:
                    seqObj.comment = aLine
                    #print "got comment '%s' for pir seq %s" % (seqObj.comment, seqObj.name)
                #else:
                #    print "No comment line for %s" % seqObj.name
            except IndexError:
                gm.append("premature end to pir file, in sequence %s" % seqObj.name)
                raise Glitch, gm
            
            while 1:
                lNum += 1
                try:
                    aLine = ll[lNum]
                    #print "e %4i aLine: %s" % (lNum, aLine)
                except IndexError:
                    break
                if not aLine:
                    gm.append("Misplaced blank line in pir sequence %s" % seqObj.name)
                    raise Glitch, gm
                if aLine[0] == '>':
                    break
                seqObj.sequence.append(aLine)
                if aLine.endswith('*'):
                    break
            seqObj.sequence = ''.join(seqObj.sequence)
            assert seqObj.sequence.endswith('*')
            seqObj.sequence = seqObj.sequence[:-1]
            self.sequences.append(seqObj)

        # now fix the sequences
        toLowerTransTable = string.maketrans(string.uppercase[:26], string.lowercase[:26])
        for seqObj in self.sequences:
            if string.count(seqObj.sequence, '.'):
                gm.append("Dots don't work in a pir file, do they?")
                raise Glitch, gm
            seqObj.sequence = string.translate(seqObj.sequence, toLowerTransTable,
                                          string.digits + string.whitespace + '\0')
        if 0:
            for seqObj in self.sequences:
                print '%20s  %-30s' % ('name', seqObj.name)
                print '%20s  %-30s' % ('comment', seqObj.comment)
                print '%20s  %-30s' % ('sequence', seqObj.sequence)
                print '%20s  %-30s' % ('dataType', seqObj.dataType)
                print ''

        # check for invalid chars
        if len(self.sequences) > 0:
            bads = 0
            if self.sequences[0].dataType == 'dna':
                for s in self.sequences:
                    j = 0
                    while j < len(s.sequence):
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
                    gm.append("Got bad characters.")
                    raise Glitch, gm
            if self.sequences[0].dataType == 'protein':
                for s in self.sequences:
                    j = 0
                    while j < len(s.sequence):
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
                    gm.append("Got bad characters.")
                    raise Glitch, gm
        flob.close()
        return self  # ie success


##    def checkForEqualSequenceLengths(self):
##        if len(self.sequences) <= 1:
##            # print "SequenceList: checkForEqualSequenceLengths:"
##            # print "    The list has only one sequence.  Returning zero."
##            return None
##        else:
##            l0 = len(self.sequences[0].sequence)
##            for s in self.sequences[1:]:
##                if l0 != len(s.sequence):
##                    return False
##            return True

##    ##Ignore
##    def checkSequenceDataTypes(self):
##        if len(self.sequences) == 1:
##            # print "SequenceList: checkSequenceDataTypes:"
##            # print "    The list has only one sequence.  Returning zero."
##            return 0
##        elif len(self.sequences) > 1:
##            t0 = self.sequences[0].dataType
##            print "checkSequencDataTypes()  first sequence is %s" % t0
##            for s in self.sequences[1:]:
##                if t0 != s.dataType:
##                    print "checkSequencDataTypes()  sequence %s is %s" % (s.name, s.dataType)
##                    return 0
##            if t0 == 'dna':
##                return 1
##            elif t0 == 'protein':
##                return 2
##        else:
##            print "SequenceList: checkSequenceDataTypes:"
##            print "    No sequences.  Returning zero"
##            return 0


    def alignment(self):
        """Make self into an alignment, and return it.

        If all the sequences are the same length and type, then self,
        a sequenceList, could be an Alignment.  This method generates
        an Alignment instance, runs the Alignment method
        checkLengthsAndTypes(), and returns the Alignment.

        If you feed p4 a fasta sequence, it makes SequenceList object,
        and runs this method on it.  If it works then p4 puts the
        Alignment object in var.alignments, and if not it puts the
        SequenceList object in var.sequenceLists.

        It is possible that p4 might think that some short sequences
        are DNA when they are really protein.  In that case it will
        fail to make an alignment, because it will fail the types
        check.  So what you can do is something like this::

            sl = var.sequenceLists[0]
            for s in sl.sequences:
                s.dataType = 'protein'
            a = sl.alignment()

        """

        from Alignment import Alignment
        a = Alignment()
        a.fName = self.fName
        import copy
        a.sequences = copy.deepcopy(self.sequences)  # self will be deleted
        a.fName = self.fName
        a.checkLengthsAndTypes()
        return a


    def writeFasta(self, fName=None, comment=1, width=60, append=0, seqNum=None):
        """Write out the sequences in Fasta format.

        This will write to stdout by default, or a file name, or to an
        open file-like object, eg a StringIO object.

        The sequences may have comments, which are written by default.
        If you don't want comments, say comment=None

        If seqNum=None, the default, then all the sequences are
        written.  But you can also just write one sequence, given by
        its number.   Write out a bunch to the same file with 'append'.
        """

        complaintHead = '\nSequenceList.writeFasta()'

        isFlob = False
        originalTell = None
        if fName == None or fName == sys.stdout:
            f = sys.stdout
            isFlob = True
        elif hasattr(fName, 'write'):  # an open, file-like object
            f = fName
            originalTell = f.tell()
            isFlob = True
        else:
            if append:
                if os.path.isfile(fName):
                    try:
                        f = open(fName, 'a')
                    except IOError:
                        print complaintHead
                        print "    Can't open %s for appending." % fName
                        sys.exit()
                else:
                    if 0:
                        print complaintHead
                        print "    'append' is requested,"
                        print "    but '%s' is not a regular file (maybe it doesn't exist?)." \
                              % fName
                        print "    Writing to a new file instead."
                    try:
                        f = open(fName, 'w')
                    except IOError:
                        print complaintHead
                        print "    Can't open %s for writing." % fName
                        sys.exit()

            else:
                try:
                    f = open(fName, 'w')
                except IOError:
                    print complaintHead
                    print "    Can't open %s for writing." % fName
                    sys.exit()

        if seqNum == None:
            for i in range(len(self.sequences)):
                s = self.sequences[i]
                f.write('>%s' % s.name)
                if comment and s.comment:
                    f.write(' %s' % s.comment)
                f.write('\n')
                left = len(s.sequence)
                pos = 0
                while left >= width:
                    if var.writeFastaUppercase:
                        f.write('%s\n' % s.sequence[pos: pos + width].upper())
                    else:
                        f.write('%s\n' % s.sequence[pos: pos + width])
                    pos = pos + width
                    left = left - width
                if left > 0:
                    if var.writeFastaUppercase:
                        f.write('%s\n' % s.sequence[pos:].upper())
                    else:
                        f.write('%s\n' % s.sequence[pos:])
                f.write('\n')
        else:
            try:
                theInt = int(seqNum)
                if theInt < 0 or theInt >= len(self.sequences):
                    print complaintHead
                    print "    seqNum %i is out of range." % seqNum
                    sys.exit()
            except ValueError:
                print complaintHead
                print "    seqNum should be an integer."
                sys.exit()
            s = self.sequences[theInt]
            f.write('>%s' % s.name)
            if comment and s.comment:
                f.write(' %s' % s.comment)
            f.write('\n')
            left = len(s.sequence)
            pos = 0
            while left >= width:
                if var.writeFastaUppercase:
                    f.write('%s\n' % s.sequence[pos: pos + width].upper())
                else:
                    f.write('%s\n' % s.sequence[pos: pos + width])
                pos = pos + width
                left = left - width
            if left > 0:
                if var.writeFastaUppercase:
                    f.write('%s\n' % s.sequence[pos:].upper())
                else:
                    f.write('%s\n' % s.sequence[pos:])
            f.write('\n')

        #f.read()
        if isFlob and f != sys.stdout:
            if hasattr(f, 'seek'):
                f.seek(originalTell)
        else:
            f.close()
        
        


    def checkNamesForDupes(self):
        if not var.doCheckForDuplicateSequenceNames:
            return
        snDict = {}
        for s in self.sequences:
            ret = snDict.get(s.name)
            if ret:
                snDict[s.name] += 1
            else:
                snDict[s.name] = 1
        nDupes = 0
        for k,v in snDict.iteritems():
            if v > 1:
                print "Got %2i copies of sequence name %s" % (v, k)
                nDupes += 1
        if nDupes:
            gm = ["SequenceList.checkNamesForDupes()"]
            if self.fName:
                gm.append("File name %s" % self.fName)
            gm.append("Got %i duplicate sequence names." % nDupes)
            gm.append("(If you want to turn off checking, set ")
            gm.append("var.doCheckForDuplicateSequenceNames to False)")
            raise Glitch, gm


    def dump(self):
        if isinstance(self, SequenceList):
            print "\nSequenceList dump:"
            if self.fName:
                print "  File name is %s" % self.fName
        if len(self.sequences) == 1:
            print "  There is 1 sequence"
        else:
            print "  There are %s sequences" % len(self.sequences)

        nSeqsToDo = len(self.sequences)
        if nSeqsToDo > 12:
            nSeqsToDo = 10
        for i in range(nSeqsToDo):
            if isinstance(self, SequenceList):
                print "  %3i %5s %s" % (i, len(self.sequences[i].sequence), self.sequences[i].name)
            else: # Alignment, don't print sequence lengths
                print "  %3i   %s" % (i, self.sequences[i].name)
            # self.sequences[i].dump()
        if len(self.sequences) > nSeqsToDo:
            print "  ... and %i others..." % (len(self.sequences) - nSeqsToDo)
        print ''




    def renameForPhylip(self, dictFName='p4_renameForPhylip_dict.py'):
        """Rename with strict phylip-friendly short boring names.

        It saves the old names (together with the new) in a python
        dictionary, in a file, by default named
        p4_renameForPhylip_dict.py"""

        gm = ['SequenceList.renameForPhylip()']
        if os.path.exists(dictFName):
            gm.append("The dictionary file '%s' already exists." % dictFName)
            raise Glitch, gm
        if hasattr(self, 'taxNames'):
            originalNames = self.taxNames
        else:
            originalNames = [s.name for s in self.sequences]
        d = {}
        for i in range(len(self.sequences)):
            s = self.sequences[i]
            newName = 's%i' % i
            d[newName] = s.name
            s.name = newName
        f = file(dictFName, 'w')
        f.write("p4_renameForPhylip_originalNames = %s\np4_renameForPhylip_dict = %s\n" % (originalNames,d))
        f.close()

    def restoreNamesFromRenameForPhylip(self, dictFName='p4_renameForPhylip_dict.py'):
        """Given the dictionary file, restore proper names.

        The dictionary file is by default named p4_renameForPhylip_dict.py"""

        gm = ["SequenceLists.restoreNamesFromRenameForPhylip()"]
        if os.path.exists(dictFName):
            import __main__
            execfile(dictFName, __main__.__dict__,  __main__.__dict__)
            from __main__ import p4_renameForPhylip_dict
        else:
            gm.append("The dictionary file '%s' can't be found." % dictFName)
            raise Glitch, gm
        for s in self.sequences:
            if p4_renameForPhylip_dict.has_key(s.name):
                s.name = p4_renameForPhylip_dict[s.name]
            else:
                gm.append("The dictionary does not contain a key for '%s'." % s.name)
                raise Glitch, gm
        del(__main__.p4_renameForPhylip_dict)
        del(__main__.p4_renameForPhylip_originalNames)



    def muscle(self):
        """Do an alignment with muscle.

        Its all done in memory -- no files are written.

        An alignment object is returned.

        The order of the sequences in the new alignment is made to be
        the same as the order in self.

        """
        flob = cStringIO.StringIO()
        self.writeFasta(fName=flob)
        p = Popen(["muscle"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        ret = p.communicate(input=flob.getvalue())
        flob.close()
        try:
            a = func.readAndPop(ret[0])
        except Glitch:
            print ret
            raise Glitch, "Something didn't work ..."

        a.makeSequenceForNameDict()
        newSequenceList = []
        for sSelf in self.sequences:
            newSequenceList.append(a.sequenceForNameDict[sSelf.name])
        a.sequences = newSequenceList
        return a

    def clustalo(self):
        """Do an alignment with clustalo.

        Its all done in memory -- no files are written.

        An alignment object is returned.

        The order of the sequences in the new alignment is made to be
        the same as the order in self.

        """
        flob = cStringIO.StringIO()
        self.writeFasta(fName=flob)
        p = Popen(["clustalo", "-i", "-"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        ret = p.communicate(input=flob.getvalue())
        #ret = p.communicate()
        if ret[1]:
            print ret
            raise Glitch, "clustalo()  Something wrong here ..."
        flob.close()
        a = func.readAndPop(ret[0])
        a.makeSequenceForNameDict()
        newSequenceList = []
        for sSelf in self.sequences:
            newSequenceList.append(a.sequenceForNameDict[sSelf.name])
        a.sequences = newSequenceList
        return a
    
            

