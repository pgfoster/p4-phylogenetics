import sys
import re
import string
import os
import io
import p4.func
import copy
from p4.var import var
from p4.p4exceptions import P4Error
from subprocess import Popen, PIPE
from builtins import object       # For Py2/3 compatibility, needed for redefinition of __bool__() below in Py2
from p4.sequence import Sequence

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

    def __init__(self, flob=None):

        #: A list of Sequence objects
        self.sequences = []
        #: If it came from a file with a name, this is it.
        self.fName = None
        #: Allows you to find Sequence objects from their Sequence.name
        self.sequenceForNameDict = None
        if flob:
            self._readFastaFile(flob)
            if hasattr(flob, 'name'):
                self.fName = flob.name

    def makeSequenceForNameDict(self):
        self.sequenceForNameDict = {}
        for s in self.sequences:
            assert s.name not in self.sequenceForNameDict, "duped name %s" % s.name
            self.sequenceForNameDict[s.name] = s

    def _readFastaMakeSeq(self, splHeadLine, sLineList):
        gm = ['SequenceList._readFastaMakeSeq()']

        if not splHeadLine or not splHeadLine[0]:
            gm.append("No name for new fasta sequence.  This should not happen.")
            raise P4Error(gm)
        if not sLineList:
            gm.append("No sequence for %s" % splHeadLine)
            raise P4Error(gm)

        mySeq = Sequence()
        mySeq.name = splHeadLine[0]
        if len(splHeadLine) == 2:
            mySeq.comment = splHeadLine[1]
        mySeq.sequence = ''.join(sLineList).lower()
        return mySeq

        
    def _readFastaReadHeadLine(self, aLine):
        gm = ['SequenceList._readFastaReadHeadLine(%s)']
        assert aLine.startswith(">")  # or else we would not be here.
        # There should be no space after the ">"
        if aLine[1] in string.whitespace:
            gm.append("The '>' should not be followed by whitespace.")
            raise P4Error(gm)

        # In this next line, the comment, if it exists, picks up a newline.  Get
        # rid of it with strip().
        splHeadLine = [myWord.strip() for myWord in aLine[1:].split(None, 1)]
        return splHeadLine
        
        

    def _readFastaFile(self, flob):
        flob.seek(0)
        gm = ['SequenceList._readFastaFile()']
        if hasattr(flob, 'name'):
            gm = ['SequenceList._readFastaFile(%s)' % flob.name]

        complaintAboutLength = """
        Lines should not be longer than 120 characters.
        This will be overlooked here, but other programs may gag.
        """
        alreadyComplainedAboutLength = False

        # The first line might start with a ';'
        # Move the position to the first '>'
        aLine = flob.readline()
        while aLine[0] != '>':
            aLine = flob.readline()

        if not aLine:
            gm.append("Unable to find a line that starts with '>'")
            raise P4Error(gm)

        sList = []
        splHeadLine = self._readFastaReadHeadLine(aLine)
        # print(splHeadLine)

        # read the rest of the flob
        while 1:
            aLine = flob.readline()
            # print("read aLine: %s" % aLine, end='')

            # If we are at the end, stash the last sequence and break.
            if not aLine:
                if not splHeadLine:
                    break
                else:
                    if not sList:
                        gm.append("No sequence for '%s'?" % splHeadLine[0])
                        raise P4Error(gm)
                    mySeq = self._readFastaMakeSeq(splHeadLine, sList)
                    # print("Got seq, %s, %s, %s" % (mySeq, mySeq.name, mySeq.sequence))
                    self.sequences.append(mySeq)
                    del(sList)
                    break
                
            elif aLine[0] == '>':
                # Stash the previous sequence
                if not sList:
                    gm.append("No sequence for '%s'?" % splHeadLine[0])
                    raise P4Error(gm)
                mySeq = self._readFastaMakeSeq(splHeadLine, sList)
                self.sequences.append(mySeq)

                sList = []
                splHeadLine = self._readFastaReadHeadLine(aLine)
                
            else:
                sList.append(aLine.strip())


        # now fix the sequences
        for mySeq in self.sequences:
            dType = p4.func.isDnaRnaOrProtein(mySeq.sequence)
            if dType == 1:
                # print "Its dna"
                mySeq.dataType = 'dna'
            elif dType == 2:
                # print "Its rna"
                # print "Converting RNA to DNA"
                mySeq.sequence = list(mySeq.sequence)
                for i in range(len(mySeq.sequence)):
                    if mySeq.sequence[i] == 'u':
                        mySeq.sequence[i] = 't'
                mySeq.sequence = ''.join(mySeq.sequence)
                mySeq.dataType = 'dna'
            else:
                # print "Its protein"
                mySeq.dataType = 'protein'

        if 0:
            for mySeq in self.sequences:
                print('%20s  %-30s' % ('name', mySeq.name))
                print('%20s  %-30s' % ('comment', mySeq.comment))
                print('%20s  %-30s' % ('sequence', mySeq.sequence))
                print('%20s  %-30s' % ('dataType', mySeq.dataType))
                print('')

        # check for invalid chars
        if len(self.sequences) > 0:
            bads = 0
            if self.sequences[0].dataType == 'dna':
                for s in self.sequences:
                    j = 0
                    while j < len(s.sequence):
                        if s.sequence[j] not in var.validDnaChars:
                            print("bad character '%s' in (zero-based) dna sequence %s " % \
                                (s.sequence[j], self.sequences.index(s)))
                            print("          sequence name: %s" % s.name)
                            print("          at (zero-based) position %s" % j)
                            bads = bads + 1
                            if bads > 10:
                                print("...and possibly others")
                                break
                        j = j + 1
                    if bads > 10:
                        break
                if bads:
                    gm.append("Got bad characters.")
                    raise P4Error(gm)
            if self.sequences[0].dataType == 'protein':
                for s in self.sequences:
                    j = 0
                    while j < len(s.sequence):
                        if s.sequence[j] not in var.validProteinChars:
                            print("bad character '%s' in (zero-based) protein sequence %s " % \
                                (s.sequence[j], self.sequences.index(s)))
                            print("          sequence name: %s" % s.name)
                            print("          at (zero-based) position %s" % j)
                            bads = bads + 1
                            if bads > 10:
                                print("...and possibly others")
                                break
                        j = j + 1
                    if bads > 10:
                        break
                if bads:
                    gm.append("Got bad characters.")
                    raise P4Error(gm)
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
                # print "a %4i aLine: '%s'" % (lNum, aLine)
            except IndexError:
                break
            # print "a1 %4i aLine: '%s'" % (lNum, aLine)
            if not aLine or not aLine.startswith('>'):
                try:
                    lNum += 1
                    aLine = ll[lNum]
                    # print "b %4i aLine: '%s'" % (lNum, aLine)
                except IndexError:
                    break
            # print "c %4i aLine: %s" % (lNum, aLine)
            if aLine[3] != ';':
                gm.append("First line is: %s" % aLine.rstrip())
                gm.append("4th char should be ';'")
                raise P4Error(gm)
            twoChars = aLine[1:3]
            if twoChars not in ['P1']:
                gm.append("First line is: %s" % aLine.rstrip())
                gm.append(
                    "Code characters '%s' are not recognized / implemented.  Fix me?" % twoChars)
                raise P4Error(gm)

            seqObj = Sequence()
            if twoChars == 'P1':
                seqObj.dataType = 'protein'
            else:
                gm.append(
                    "Pir datatype code '%s' is not implemented.  Fix me." % twoChars)
                raise P4Error(gm)
            # So I can append lines.  I'll change it back to a string later
            seqObj.sequence = []
            splLine = aLine.split(';')
            seqObj.name = splLine[1]
            # print "got pir seq name %s" % seqObj.name

            # Get the comment line.
            lNum += 1
            try:
                aLine = ll[lNum]
                # print "d %4i aLine: %s" % (lNum, aLine)
                if aLine:
                    seqObj.comment = aLine
                    # print "got comment '%s' for pir seq %s" % (seqObj.comment, seqObj.name)
                # else:
                #    print "No comment line for %s" % seqObj.name
            except IndexError:
                gm.append(
                    "premature end to pir file, in sequence %s" % seqObj.name)
                raise P4Error(gm)

            while 1:
                lNum += 1
                try:
                    aLine = ll[lNum]
                    # print "e %4i aLine: %s" % (lNum, aLine)
                except IndexError:
                    break
                if not aLine:
                    gm.append(
                        "Misplaced blank line in pir sequence %s" % seqObj.name)
                    raise P4Error(gm)
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
        myZaps = string.digits + string.whitespace + '\0'
        for seqObj in self.sequences:
            if '.' in seqObj.sequence:
                gm.append("Dots don't work in a pir file, do they?")
                raise P4Error(gm)
            seqObj.sequence = seqObj.sequence.lower()
            seqObj.sequence = re.sub('['+myZaps+']', '', seqObj.sequence)
            
        if 0:
            for seqObj in self.sequences:
                print('%20s  %-30s' % ('name', seqObj.name))
                print('%20s  %-30s' % ('comment', seqObj.comment))
                print('%20s  %-30s' % ('sequence', seqObj.sequence))
                print('%20s  %-30s' % ('dataType', seqObj.dataType))
                print('')

        # check for invalid chars
        if len(self.sequences) > 0:
            bads = 0
            if self.sequences[0].dataType == 'dna':
                for s in self.sequences:
                    j = 0
                    while j < len(s.sequence):
                        if s.sequence[j] not in var.validDnaChars:
                            print("bad character '%s' in (zero-based) dna sequence %s " % \
                                (s.sequence[j], self.sequences.index(s)))
                            print("          sequence name: %s" % s.name)
                            print("          at (zero-based) position %s" % j)
                            bads = bads + 1
                            if bads > 10:
                                print("...and possibly others")
                                break
                        j = j + 1
                    if bads > 10:
                        break
                if bads:
                    gm.append("Got bad characters.")
                    raise P4Error(gm)
            if self.sequences[0].dataType == 'protein':
                for s in self.sequences:
                    j = 0
                    while j < len(s.sequence):
                        if s.sequence[j] not in var.validProteinChars:
                            print("bad character '%s' in (zero-based) protein sequence %s " % \
                                (s.sequence[j], self.sequences.index(s)))
                            print("          sequence name: %s" % s.name)
                            print("          at (zero-based) position %s" % j)
                            bads = bads + 1
                            if bads > 10:
                                print("...and possibly others")
                                break
                        j = j + 1
                    if bads > 10:
                        break
                if bads:
                    gm.append("Got bad characters.")
                    raise P4Error(gm)
        flob.close()
        return self  # ie success

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

        from p4.alignment import Alignment
        a = Alignment()
        a.fName = self.fName
        import copy
        a.sequences = copy.deepcopy(self.sequences)  # self will be deleted
        a.fName = self.fName
        a.checkLengthsAndTypes()
        return a

    def writeFasta(self, fName=None, comment=1, width=60, append=0, seqNum=None, writeExtraNewline=True):
        """Write out the sequences in Fasta format.

        This will write to stdout by default, or a file name, or to an
        open file-like object, eg a StringIO object.

        The sequences may have comments, which are written by default.
        If you don't want comments, say comment=None

        By default, sequences are wrapped when they are too long.
        You can set the length at which to wrap the sequences.
        Set width=0 if you want your sequences in one (long) line.

        If seqNum=None, the default, then all the sequences are
        written.  But you can also just write one sequence, given by
        its number.   Write out a bunch to the same file with 'append'.

        By default, a blank line will be written after each sequence.
        If you prefer your fasta without these extra lines, say
        writeExtraNewline=False.
        """

        complaintHead = '\nSequenceList.writeFasta()'

        isFlob = False
        #originalTell = None
        if fName == None or fName == sys.stdout:
            f = sys.stdout
            isFlob = True
        elif hasattr(fName, 'write'):  # an open, file-like object
            f = fName
            #originalTell = f.tell()
            isFlob = True
        else:
            assert isinstance(fName, str)  # a file name
            if append:
                if os.path.isfile(fName):
                    try:
                        f = open(fName, 'a')
                    except IOError:
                        print(complaintHead)
                        print("    Can't open %s for appending." % fName)
                        sys.exit()
                else:
                    if 0:
                        print(complaintHead)
                        print("    'append' is requested,")
                        print("    but '%s' is not a regular file (maybe it doesn't exist?)." \
                              % fName)
                        print("    Writing to a new file instead.")
                    try:
                        f = open(fName, 'w')
                    except IOError:
                        print(complaintHead)
                        print("    Can't open %s for writing." % fName)
                        sys.exit()

            else:
                try:
                    f = open(fName, 'w')
                except IOError:
                    print(complaintHead)
                    print("    Can't open %s for writing." % fName)
                    sys.exit()

        if seqNum == None:
            for i in range(len(self.sequences)):
                s = self.sequences[i]
                #print("SequenceList.writeFasta() s.name %s, type(s.name) %s" % (s.name, type(s.name)))
                f.write('>%s' % s.name)
                if comment and s.comment:
                    f.write(' %s' % s.comment)
                f.write('\n')
                left = len(s.sequence)
                pos = 0
                if width > 0:
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
                if writeExtraNewline:
                    f.write('\n')
        else:
            try:
                theInt = int(seqNum)
                if theInt < 0 or theInt >= len(self.sequences):
                    print(complaintHead)
                    print("    seqNum %i is out of range." % seqNum)
                    sys.exit()
            except ValueError:
                print(complaintHead)
                print("    seqNum should be an integer.")
                sys.exit()
            s = self.sequences[theInt]
            f.write('>%s' % s.name)
            if comment and s.comment:
                f.write(' %s' % s.comment)
            f.write('\n')
            left = len(s.sequence)
            pos = 0
            if width > 0:
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
            if writeExtraNewline:
                f.write('\n')

        if isFlob:
            if f != sys.stdout:
                f.close()
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
        for k, v in snDict.items():
            if v > 1:
                print("Got %2i copies of sequence name %s" % (v, k))
                nDupes += 1
        if nDupes:
            gm = ["SequenceList.checkNamesForDupes()"]
            if self.fName:
                gm.append("File name %s" % self.fName)
            gm.append("Got %i duplicate sequence names." % nDupes)
            gm.append("(If you want to turn off checking, set ")
            gm.append("var.doCheckForDuplicateSequenceNames to False)")
            raise P4Error(gm)

    def dump(self):
        if isinstance(self, SequenceList):
            print("\nSequenceList dump:")
            if self.fName:
                print("  File name is %s" % self.fName)
        if len(self.sequences) == 1:
            print("  There is 1 sequence")
        else:
            print("  There are %s sequences" % len(self.sequences))

        nSeqsToDo = len(self.sequences)
        if nSeqsToDo > 12:
            nSeqsToDo = 10
        for i in range(nSeqsToDo):
            if isinstance(self, SequenceList):
                print("  %3i %5s %s" % (i, len(self.sequences[i].sequence), self.sequences[i].name))
            else:  # Alignment, don't print sequence lengths
                print("  %3i   %s" % (i, self.sequences[i].name))
            # self.sequences[i].dump()
        if len(self.sequences) > nSeqsToDo:
            print("  ... and %i others..." % (len(self.sequences) - nSeqsToDo))
        print('')

    def renameForPhylip(self, dictFName='p4_renameForPhylip_dict.py'):
        """Rename with strict phylip-friendly short boring names.

        It saves the old names (together with the new) in a python
        dictionary, in a file, by default named
        p4_renameForPhylip_dict.py"""

        gm = ['SequenceList.renameForPhylip()']
        if os.path.exists(dictFName):
            gm.append("The dictionary file '%s' already exists." % dictFName)
            raise P4Error(gm)
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
        f = open(dictFName, 'w')
        f.write("p4_renameForPhylip_originalNames = %s\np4_renameForPhylip_dict = %s\n" % (
            originalNames, d))
        f.close()

    def restoreNamesFromRenameForPhylip(self, dictFName='p4_renameForPhylip_dict.py'):
        """Given the dictionary file, restore proper names.

        The dictionary file is by default named p4_renameForPhylip_dict.py"""

        gm = ["SequenceLists.restoreNamesFromRenameForPhylip()"]
        if os.path.exists(dictFName):
            import __main__
            exec(open(dictFName).read(), __main__.__dict__,  __main__.__dict__)
            from __main__ import p4_renameForPhylip_dict
        else:
            gm.append("The dictionary file '%s' can't be found." % dictFName)
            raise P4Error(gm)
        for s in self.sequences:
            if s.name in p4_renameForPhylip_dict:
                s.name = p4_renameForPhylip_dict[s.name]
            else:
                gm.append("The dictionary does not contain a key for '%s'." % s.name)
                raise P4Error(gm)
        del(__main__.p4_renameForPhylip_dict)
        del(__main__.p4_renameForPhylip_originalNames)

    def muscle(self):
        """Do an alignment with muscle.

        Its all done in memory -- no files are written.

        An alignment object is returned.

        The order of the sequences in the new alignment is made to be
        the same as the order in self.

        """
        flob = io.BytesIO()
        self.writeFastaToBytesFlob(flob)
        p = Popen(["muscle"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        ret = p.communicate(input=flob.getvalue())
        flob.close()
        try:
            a = p4.func.readAndPop(ret[0].decode())
        except P4Error:
            print(ret)
            raise P4Error("Something didn't work ...")

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
        flob = io.BytesIO()          # gotta be Bytes for subprocess
        self.writeFastaToBytesFlob(flob)
        p = Popen(["clustalo", "-i", "-"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        ret = p.communicate(input=flob.getvalue())
        #ret = p.communicate()
        if ret[1]:
            print(ret)
            raise P4Error("clustalo()  Something wrong here ...")
        flob.close()
        #print(ret)      # it is a bytes string
        a = p4.func.readAndPop(ret[0].decode())
        a.makeSequenceForNameDict()
        newSequenceList = []
        for sSelf in self.sequences:
            newSequenceList.append(a.sequenceForNameDict[sSelf.name])
        a.sequences = newSequenceList
        return a

    def writeFastaToBytesFlob(self, flob):
        """For subprocesses, eg muscle and clustalo

        No comment, no extra new line.  Width 60.
        """
        assert isinstance(flob, io.BytesIO)
        width = 60
        for s in self.sequences:
            st = '>%s\n' % s.name
            flob.write(str.encode(st))
            left = len(s.sequence)
            pos = 0
            while left >= width:
                if var.writeFastaUppercase:
                    st = '%s\n' % s.sequence[pos: pos + width].upper()
                    flob.write(str.encode(st))
                else:
                    st = '%s\n' % s.sequence[pos: pos + width]
                    flob.write(str.encode(st))
                pos = pos + width
                left = left - width
            if left > 0:
                if var.writeFastaUppercase:
                    st = '%s\n' % s.sequence[pos:].upper()
                    flob.write(str.encode(st))
                else:
                    st = '%s\n' % s.sequence[pos:]
                    flob.write(str.encode(st))


        
