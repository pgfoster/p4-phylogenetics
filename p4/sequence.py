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

class Sequence(object):

    """A container for a single molecular sequence.

    Data attributes

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
    """(property) return the length of the sequence, or zero"""

    # See the comment in alignment.py, for the same redefinition of __bool__().
    def __bool__(self):
        return True

    def __len__(self):
        if self.sequence:
            return len(self.sequence)
        else:
            return 0

    def dump(self):
        """Print rubbish about self."""
        print('%15s: %s' % ('name', self.name))
        if self.comment:
            print('%15s: %s' % ('comment', self.comment))
        # if self.dataType == 'dna':
           # if self.transl_table:
           #    print "%15s: %s' % ('transl_table", self.transl_table)
        if self.sequence:
            print('%15s: %s' % ('sequence', self.sequence[:25]), end=' ')
            if len(self.sequence) > 25:
                print("...")
            else:
                print('')

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
            # elif c == 'x':
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
                if c in string.ascii_uppercase:
                    gm.append("Got uppercase '%s' How did that happen? -- can only handle lowercase." % c)
                else:
                    gm.append("Sequence.reverseComplement().  Got char '%s' What is it?" % c)
                raise P4Error(gm)

        self.sequence = ''.join(self.sequence)

    def writeFastaToOpenFile(self, flob, width=60, doComment=True, writeExtraNewline=True):
        flob.write('>%s' % self.name)
        if doComment and self.comment:
            flob.write(' %s' % self.comment)
        flob.write('\n')
        left = len(self.sequence)
        pos = 0
        if width > 0:
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

    def writeFasta(self, fName=None, width=60, doComment=True, writeExtraNewline=True):
        isFlob = False
        if not fName or fName == sys.stdout:
            f = sys.stdout
            isFlob = True
        elif hasattr(fName, 'write'):  # an open file-like object
            f = fName
            isFlob = True
        else:
            f = open(fName, 'w')
        self.writeFastaToOpenFile(f, width=width, doComment=doComment, writeExtraNewline=writeExtraNewline)
        if not isFlob:
            f.close()

    def translate(self, transl_table=1, checkStarts=False, nnn_is_gap=False):
        """Returns a protein Sequence from self, a DNA sequence.

        Self is translated using
        :meth:`p4.geneticcode.GeneticCode.translate`, so it handles
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

        (These are found in :class:`p4.geneticcode.GeneticCode`)

        See also :meth:`p4.alignment.Alignment.checkTranslation` and
        :meth:`p4.alignment.Alignment.translate`.

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
            raise P4Error(gm)

        if self.nChar % 3 != 0:
            gm.append("The length of self should be a multiple of 3")
            raise P4Error(gm)
        nTriplets = self.nChar / 3

        from geneticcode import GeneticCode
        gc = GeneticCode(transl_table)

        prSeq = Sequence()
        prSeq.dataType = 'protein'
        prSeq.name = self.name
        prSeq.sequence = ['-'] * nTriplets

        dnaSeq = self.sequence
        protSeq = prSeq.sequence
        for j in range(nTriplets):
            theCodon = dnaSeq[(j * 3):(j * 3) + 3]
            # print theCodon
            if theCodon == '---':
                protSeq[j] = '-'
            elif theCodon.count('-'):
                print("    seq %s, position %4i, dnaSeq %4i, codon '%s' is incomplete" % (self.name, j, (j * 3), theCodon))
            elif theCodon == 'nnn':
                if nnn_is_gap:
                    print("    seq %s, position %4i, dnaSeq %4i, codon '%s' translating to a gap ('-')" % (
                        self.name, j, (j * 3), theCodon))
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

        prSeq.sequence = ''.join(prSeq.sequence)
        return prSeq

    def checkTranslation(self, theProteinSequence, transl_table=1, checkStarts=False):
        """Check that self translates to theProteinSequence

        Self should be a DNA sequence.  It is translated using
        :meth:`p4.geneticcode.GeneticCode.translate` (so it should handle
        ambiguities) and compared against theProteinSequence.  The
        theProteinSequence name and gap pattern should
        be the same as in the DNA sequence.  The default transl_table is
        the standard (or so-called universal) genetic code.

        Other available translation tables, this week::

            if transl_table == 1:   # standard
            elif transl_table == 2: # vertebrate mito
            elif transl_table == 4: # Mold, Protozoan,
                                    # and Coelenterate Mitochondrial Code
                                    # and the Mycoplasma/Spiroplasma Code
            elif transl_table == 5: # invertebrate mito
            elif transl_table == 9: # echinoderm mito

            # and now 6, 10, 11, 12, 13, 14, 21.

        (These are found in :class:`~p4.geneticcode.GeneticCode`)

        See also :meth:`p4.alignment.Alignment.translate` and 
        :meth:`p4.alignment.Alignment.checkTranslation`

        If the arg *checkStarts* is turned on (by default it is not turned
        on) then this method checks whether the first codon is a start
        codon.

        """

        gm = ['Sequence.checkTranslation()']
        if self.dataType != 'dna':
            gm.append("Self should be a DNA sequence.")
            raise P4Error(gm)
        if not theProteinSequence or \
                not isinstance(theProteinSequence, p4.sequencelist.Sequence) or \
                theProteinSequence.dataType != 'protein':
            gm.append("Something wrong with theProteinSequence")
            raise P4Error(gm)

        if self.name != theProteinSequence.name:
            gm.append("The sequence names of self and theProteinSequence are not the same")
            raise P4Error(gm)

        if self.nChar != (3 * theProteinSequence.nChar):
            gm.append("The length of the DNA sequence should be 3 times that of theProteinSequence")
            gm.append("DNA sequence (self):  %i" % self.nChar)
            gm.append("Protein sequence:     %i  ( * 3 = %i)" %
                      (theProteinSequence.nChar, (3 * theProteinSequence.nChar)))
            raise P4Error(gm)

        gc = p4.geneticcode.GeneticCode(transl_table)

        pLen = theProteinSequence.nChar
        crimes = 0
        for j in range(pLen):
            theCodon = self.sequence[(3 * j) + 0] + \
                self.sequence[(3 * j) + 1] + \
                self.sequence[(3 * j) + 2]
            if theCodon == '---':
                if theProteinSequence.sequence[j] != '-':
                    print("    position %4i, codon '---' is '%s', should be '-'" % (j, theProteinSequence.sequence[j]))
                    crimes += 1
            elif theCodon.count('-'):
                print("    position %4i, codon '%s' is incomplete" % (j, theCodon))
                crimes += 1
            # elif theCodon in gc.code:
            #     if gc.code[theCodon] != theProteinSequence.sequence[j]:
            #         print "    position %4i, codon '%s' is '%s', should be '%s'" % (
            #             j, theCodon, theProteinSequence.sequence[j], gc.code[theCodon])
            #         crimes += 1
            # else:
            #     print("    position %4i, codon '%s' is not a known codon" % (j, theCodon))
            #     crimes += 1
            else:
                tr = gc.translate(theCodon)
                if tr != theProteinSequence.sequence[j]:
                    print("    position %4i, codon '%s' is '%s', should be '%s'" % (
                        j, theCodon, theProteinSequence.sequence[j], gc.code[theCodon]))
                    crimes += 1

                # If arg checkStarts is turned on -- Is the first
                # codon a start?  -- if not, it is not a crime
                if checkStarts and j == 0:
                    if theCodon in gc.startList:
                        print("    Seq %i (%s). The first codon, '%s', is a start codon" % (i, self.name, theCodon))
                    else:
                        print("    Seq %i (%s). The first codon, '%s', is not a start codon" % (i, self.name, theCodon))
            if crimes > 6:
                break
        if crimes > 6:
            print("    ... and possibly others, skipped.")



