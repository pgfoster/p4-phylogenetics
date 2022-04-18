"""Various functions."""
import os
import sys
import re
import string
import math
import io
import random
import glob
import time
import pickle
import random
import inspect
import datetime
import subprocess

from p4.var import var
# from p4.sequencelist import Sequence, SequenceList
import p4.sequencelist
from p4.alignment import Alignment
from p4.nexus import Nexus
from p4.tree import Tree
from p4.node import Node
from p4.p4exceptions import P4Error
from p4.constraints import Constraints
import p4.pf as pf
import numpy
from p4.pnumbers import Numbers
import p4.version
from p4.nexustoken import nextTok



def nexusCheckName(theName):
    """Check to see if theName conforms to Nexus standards

    See page 597 of Maddison, Swofford, and Maddison 1997, where they say "Names
    are single NEXUS words; they cannot consist entirely of digits (e.g., a
    taxon called 123 is illegal).

    The all-digit name restriction can be relaxed in p4 by setting 
    var.nexus_allowAllDigitNames.

    Single digit names are prohibited, regardless.
    """
    if not isinstance(theName, str):
        print("func.nexusCheckName() %s is not str, is %s" % (theName, type(theName)))
        return 0
    if len(theName) == 1 and theName[0] not in string.ascii_letters:
        return 0
    if not var.nexus_allowAllDigitNames:
        try:
            int(theName)
            return 0  
        except ValueError:
            return 1
    else:
        return 1



def nexusUnquoteAndDeUnderscoreName(theName):
    """Deal with underscores and quotes.  Returns theName

    If theName is not quoted, convert any underscores to spaces.  If
    theName is quoted, remove the outside quotes, and convert any
    cases of 2 single quotes in a row to 1 single quote.

    This does not appear to be used in the rest of p4.
    """

    if theName[0] == "'":
        assert theName[-1] == "'", \
            "func.nexusUnquoteAndDeUnderscoreName().  First char is a single quote, but last char is not."
        theName = theName[1:-1]
        return theName.replace("''", "'")
    if '_' in theName:
        return theName.replace('_', ' ')
    else:
        return theName


def nexusUnquoteName(theName):
    """Deal with quotes.  Returns theName

    If theName is not quoted, just return it.  If
    theName is quoted, remove the outside quotes, and convert any
    cases of 2 single quotes in a row to 1 single quote.
    """

    if theName[0] == "'":
        if theName[-1] != "'":
            gm = ['func.nexusUnquoteName()']
            gm.append('the name is %s' % theName)
            gm.append("First char is a single quote, but last char is not.")
            raise P4Error(gm)
        theName = theName[1:-1]
        if theName.count("''"):
            return theName.replace("''", "'")
        else:
            return theName
    else:
        return theName


def nexusFixNameIfQuotesAreNeeded(theName, verbose=0):
    """Add quotes if needed, for writing to a Nexus file.

    Returns a possibly modified theName.  Usually it will not need
    quoting, so it just returns theName unchanged.

    If theName is None, or if it starts with a single quote, just return it.

    If it has (nexus-defined) punctuation, or spaces, then put quotes
    around it before returning it.  If there are internal single
    quotes, then double them, nexus-style.  Except if there are any
    already doubled single quotes, don't double them.  
    """

    if theName == None:
        return theName
    if theName.startswith("'"):
        return theName
    quotesAreNeeded = 0
    # for c in var.nexus_punctuation:
    for c in var.punctuation:
        if c in theName:
            quotesAreNeeded = 1
            break
    if not quotesAreNeeded:
        for c in string.whitespace:
            if c in theName:
                quotesAreNeeded = 1
                break
    if quotesAreNeeded:
        if verbose:
            oldName = theName
        if theName.count("'"):
            # If we have doubled quotes, don't re-double them
            if theName.count("''"):
                pass
            else:
                theName = theName.replace("'", "''")
        newName = "'" + theName + "'"
        if verbose:
            print("Warning. Nexus quoting |%s| to |%s|" % (oldName, newName))
        return newName
    else:
        return theName


def isDnaRnaOrProtein(aString):
    """Attempts to determinine the data type by the composition.

    Returns 1 for DNA, 2 for RNA, and 0 for protein.  Or so it thinks.

    It only works for lowercase symbol letters.
    """
    nGaps = aString.count('-')
    nQs = aString.count('?')
    strLenNoGaps = len(aString) - (nGaps + nQs)
    acgn = (aString.count('a') +
            aString.count('c') +
            aString.count('g') +
            aString.count('n'))
    t = aString.count('t')
    u = aString.count('u')
    if 0:
        print("stringLength(no gaps) = %s" % strLenNoGaps)
        print("acgn = %f" % acgn)
        print("t = %f" % t)
        print("u = %f" % u)

    # If it is 90% or better acgn +t or +u, then assume it is dna or rna
    threshold = 0.9 * strLenNoGaps
    if acgn + t >= threshold:
        return 1
    elif acgn + u >= threshold:
        return 2

    # At this point we still don't know.  It might be dna or rna with
    # lots of ambigs.  But we will only allow that if the content of
    # acgn +t or +u is more than threshold1.  That should
    # prevent eg 'rymkswhvd' from being assigned as dna.

    threshold1 = 0.75 * strLenNoGaps
    tally = (aString.count('r') +
             aString.count('y') +
             aString.count('m') +
             aString.count('k') +
             aString.count('s') +
             aString.count('w') +
             aString.count('h') +
             aString.count('b') +
             aString.count('v') +
             aString.count('d'))
    threshold2 = 0.99 * strLenNoGaps
    # print "  string length is %i, strLenNoGaps is %i" % (len(aString), strLenNoGaps)
    # print "  threshold2 is %3.1f, acgn + t + ambigs = %i" % (threshold, (acgn + t + tally))
    # print "  threshold1 is %3.1f, acgn + t = %i" % (threshold, (acgn + t))
    if (acgn + t + tally >= threshold2) and (acgn + t >= threshold1):
        return 1
    elif (acgn + u + tally >= threshold2) and (acgn + u >= threshold1):
        return 2
    else:
        return 0


def stringZapWhitespaceAndDigits(inString):
    out = list(inString)
    theRange = range(len(out))
    for i in reversed(theRange):
        if out[i] in string.whitespace or out[i] in string.digits:
            del out[i]
    return ''.join(out)


def dump():
    """A top-level dump of p4 trees, files, sequenceLists, and alignments."""

    print("\np4 dump")
    if len(var.fileNames) == 0:
        print("    Haven't read any files")
    elif len(var.fileNames) == 1:
        print("    Read one file: %s" % var.fileNames[0])
    else:
        print("    Read files:")
        for i in var.fileNames:
            print("                %s" % i)

    if len(var.alignments) == 0:
        # print "    There are no alignments."
        pass
    elif len(var.alignments) == 1:
        print("    There is 1 alignment.")
        # if recursive:
        #    var.alignments[0].dump()
    else:
        print("    There are %i alignments." % len(var.alignments))
        # if recursive:
        #    for i in var.alignments:
        #        i.dump()

    if len(var.sequenceLists) == 0:
        # print "    There are no sequenceLists."
        pass
    elif len(var.sequenceLists) == 1:
        print("    There is 1 sequenceList.")
        # if recursive:
        #    var.sequenceLists[0].dump()
    else:
        print("    There are %i sequenceLists." % len(var.sequenceLists))
        # if recursive:
        #    for i in var.sequenceLists:
        #        i.dump()

    if len(var.trees) == 0:
        # print "    There are no trees."
        pass
    elif len(var.trees) == 1:
        print("    There is 1 tree.")
        # if recursive:
        #    var.trees[0].dump()
    else:
        print("    There are %i trees." % len(var.trees))
        # if recursive:
        #    for i in var.trees:
        #        i.dump()


def read(stuff):
    """Read in data, trees, or python code from a file or from a string.

    For example::

        read('myTreeFile.nex')

    or ::

        read('*.nex')

    or ::

        read('((A, (B, C)), D, (E, F));')

    This is meant to be the main way to get phylogenetic stuff from
    files into p4.  The argument *stuff* can be a file name, or
    filenames described with wildcards (a 'glob', eg ``*.nex``), or it
    can be a string.  If you are specifying a literal file name or a
    glob, you will of course want to put it in quotes.

    If you want it to read a file, and you mis-specify the file name,
    it will try to read the bad file name as a string (and generally
    fail, of course).

    This method will recognize these data files --

    -   nexus
    -   fasta
    -   gde (not gde flat files)
    -   clustalw (``*.aln``)
    -   phylip (sequential or interleaved)

    and these tree files --

    -   nexus
    -   phylip/newick

    If there is a suffix, one of nex, nexus, aln, phy, phylip, gde,
    fasta, fas, fsa, or py, it will use that to make a choice of what
    kind of file it is.  (So don't give a fasta file the suffix 'phy'
    or it will fail for sure.)

    Fortunately the various kinds of files (including Python) are
    almost mutually exclusive.  A fasta file must start with either a
    '>' or a ';' as the first character.  A gde file has '{' as the
    first character.  A phylip data file must have the two integers as
    the first non-whitespace things on the first line.  In a nexus
    file, the first non-whitespace character must be either a '#' or a
    '['.  (Python files might start with a '#' also.)

    It is easy to fool.  For example, start any old nonsense file with
    two integers, and this function will think that it is a phylip
    data file!  So have a little bit of care what you feed this
    function.

    If you have trouble and want to know what it is thinking as it
    attempts to read your file, set var.verboseRead=True.

    Does anybody use GDE anymore?  This only reads parts of a GDE
    formatted alignment.  It reads the name, type, and sequence.  It
    does not read the Sequence-ID, creation date, direction,
    strandedness, nor the comments.  It correctly handles offset.  It
    fills out the ends of the sequences to make an alignment. In gde
    files, MASK sequences are changed into Nexus CharSet instances.
    The charsets made from the masks are of course printed out if you
    print out the alignment as a Nexus file.

    """

    if var.verboseRead:
        print("var.verboseRead is turned on")

    # Is this still a Bug?: single node phylip or raw newick trees must be specified as
    # ()A; they cannot be specified simply as A.

    gm = ['p4.read()']
    if isinstance(stuff, str):
        pass
    else:
        gm.append("I was expecting a string argument.")
        raise P4Error(gm)
    #nAlignments = len(var.alignments)

    if os.path.exists(stuff):
        readFile(stuff)
    else:
        # Is it a glob?
        myFlist = glob.glob(stuff)
        # print "read(). stuff=%s,  glob result: %s" % (stuff, myFlist)
        if myFlist:  # It appears to be a glob
            for fName in myFlist:
                readFile(fName)
        else:  # Nope, not a glob.  Is it a string, not a file name?
            if var.warnReadNoFile:
                print("\nread()")
                print("    A file by the name specified by the argument cannot be found.")
                print("    So I am assuming that it is to be taken as a string.")
                print("    Maybe it was a mis-specified file name?")
                print("    (You can turn off this warning by turning var.warnReadNoFile off.)\n")

            if 0:
                if sys.version_info < (3,):
                    stuff = unicode(stuff)
                    flob = io.StringIO(stuff)
            if 1:
                if sys.version_info < (3,):
                    flob = io.BytesIO(stuff)
                else:
                    flob = io.StringIO(stuff)
            
            _decideFromContent('<input string>', flob)


def readFile(fName):
    """If its a data or tree file, read it.  If its python code, exec it."""

    gm = ['func.readFile(%s)' % fName]
    # print(gm)

    # I should check if the file is a text file, an executable, or whatever.
    try:
        flob = open(fName)
    except IOError:
        gm.append("Can't open %s.  Are you sure you have the right name?" % fName)
        raise P4Error(gm)

    # print(flob, type(flob), flob.name)
    # print(dir(flob))

    # See if there is an informative suffix on the file name
    # If there is a suffix, but the file cannot be read,
    # it is a serious error, and death follows.
    result = re.search('(.+)\.(.+)', fName)
    if result:
        #baseName = result.group(1)
        # print("got result.group(2) = %s" % result.group(2))
        suffix = result.group(2).lower()
        #print("readFile: got suffix '%s'" % suffix)
        if suffix == 'py':
            flob.close()
            import __main__
            # print("__main__.__dict__ is %s" % __main__.__dict__)
            
            #execfile(fName, __main__.__dict__,  __main__.__dict__)
            #exec(open(fName).read(), __main__.__dict__,  __main__.__dict__)
            #exec(flob.read(), __main__.__dict__,  __main__.__dict__)

            # The following is better than simple exec(open(fName).read())
            # because in the event of a traceback the former (below) knows the
            # file name, but the simple version does not.  As explained on
            # stackoverflow, "(The compile call isn't strictly needed, but it
            # associates the filename with the code object making debugging a
            # little easier.)"
            with open(fName) as f:
                myCode = compile(f.read(), fName, 'exec')
                exec(myCode,  __main__.__dict__,  __main__.__dict__)

            if hasattr(__main__, 'pyFileCount'):
                __main__.pyFileCount += 1
            return
        elif suffix == 'nex' or suffix == 'nexus':
            _tryToReadNexusFile(fName, flob)
            return
        elif suffix in ['fasta', 'fas', 'fsa']:
            ret = _tryToReadFastaFile(fName, flob)
            if not ret:
                gm.append("Failed to read supposed fasta file '%s'" % fName)
                raise P4Error(gm)
            return
        elif suffix == 'gde':
            ret = _tryToReadGdeFile(fName, flob)
            if not ret:
                gm.append("Failed to read supposed gde file '%s'" % fName)
                raise P4Error(gm)
            return
        elif suffix in ['pir', 'nbrf']:  
            ret = _tryToReadPirFile(fName, flob)
            if not ret:
                gm.append("Failed to read supposed pir file '%s'" % fName)
                raise P4Error(gm)
            return
        elif suffix == 'phy' or suffix == 'phylip':
            ret = _tryToReadPhylipFile(fName, flob, None)
            if not ret:
                gm.append("Failed to read supposed phylip file '%s'" % fName)
                raise P4Error(gm)
            return
        elif suffix == 'aln':
            ret = _tryToReadClustalwFile(fName, flob)
            if not ret:
                gm.append("Failed to read supposed clustalw file '%s'" % fName)
                raise P4Error(gm)
            return
        elif result.group(2) in ['p4_tPickle']:  # preserve uppercase
            if var.verboseRead:
                print("Trying to read '%s' as a pickled Tree file..." % fName)
            # It should be a binary open
            flob.close()
            flob = open(fName, "rb")
            ret = pickle.load(flob)
            if not ret:
                gm.append("Failed to read supposed p4_tPickle file '%s'." % fName)
                raise P4Error(gm)
            else:
                if isinstance(ret, Tree):
                    ret.fName = fName
                    var.trees.append(ret)
                    var.fileNames.append(fName)
                    if var.verboseRead:
                        print("Got a tree from file '%s'." % fName)
                else:
                    gm.append("Failed to get a Tree from '%s'" % fName)
                    raise P4Error(gm)
            return
        else:
            _decideFromContent(fName, flob)
    else:
        _decideFromContent(fName, flob)

    # if var.verboseRead:
    # print "(You can turn off these messages by turning var.verboseRead
    # off.)\n"


def _decideFromContent(fName, flob):
    gm = ["func._decideFromContent()"]

    firstLine = False
    while not firstLine:
        firstLine = flob.readline()
        if not firstLine:  # end of the file
            break
        firstLine = firstLine.strip()
        if firstLine:  # got some content
            break
        else:
            # print "blank line"
            pass
    # if firstLine:
    #    print "got firstLine = %s" % firstLine
    # else:
    #    print "Got a blank file."
    if not firstLine:
        gm.append("Input '%s' is empty." % fName)
        raise P4Error(gm)
    else:
        # Fasta, phylip, and clustalw files have clues on the first line
        # If these files fail to be what they are first supposed to be,
        # then the program does not die, it just complains.
        # News, May 2001.  I will allow blank lines at the beginning.

        if firstLine[0] == '>' or firstLine[0] == ';':
            if var.verboseRead:
                print("Guessing that '%s' is a fasta file..." % fName)
            ret = _tryToReadFastaFile(fName, flob, firstLine)
            if not ret:
                if var.verboseRead:
                    print("Failed to read '%s' as a fasta file." % fName)
            if ret:
                return

        elif firstLine[0] == 'C':
            if var.verboseRead:
                print("First letter is 'C'-- guessing that '%s' is a clustalw file..." % fName)
            ret = _tryToReadClustalwFile(fName, flob, firstLine)
            if not ret:
                if var.verboseRead:
                    print("Failed to read '%s' as a clustalw file." % fName)
            if ret:
                return

        else:  # Maybe it is a phylip file?
            if var.verboseRead:
                print("Guessing that '%s' is a phylip file..." % fName)
            # either data or trees

            # Deal with what punctuation is considered to be, for the tokenizer
            # This week, the default is that punctuation is nexus_punctuation
            # but it could be set by the user
            punctuationWasNexusPunctuation = False
            if var.punctuation == var.nexus_punctuation:
                punctuationWasNexusPunctuation = True
                var.punctuation = var.phylip_punctuation
                
            ret = _tryToReadPhylipFile(fName, flob, firstLine)
            
            # Reset the punctuation if it was default, not if it was set by the user
            if punctuationWasNexusPunctuation:
                var.punctuation = var.nexus_punctuation

            if not ret:
                if var.verboseRead:
                    print("Failed to read '%s' as a phylip file." % fName)
            if ret:
                return

        # If we are here, then its not a fasta or phylip file
        # The first non-white char of a Nexus file must be a '[' or a '#'
        # So first, check for that.
        # Then confirm that it is a Nexus file, by trying to intantiate a NexusFile object
        # We have to do this because the possible beginnings of Nexus files are
        # too varied, and so we need the nextTok() function to sort it all out.
        # For example, it is perfectly valid to start a Nexus file with
        #     [a comment #nexus] #Ne[this is a comment]Xus

        flob.seek(0)
        firstChar = ' '
        while firstChar in string.whitespace:
            firstChar = flob.read(1)
        # print "Got firstChar candidate '%s'" % firstChar
        # Redundant: this problem seems to be caught above, in
        # _tryToReadPhylipFile()
        if firstChar == '':
            gm.append("This file seems to be composed only of whitespace.")
            raise P4Error(gm)

        # print "got firstChar = %s" % firstChar
        flob.seek(0)
        if var.verboseRead:
            print("Guessing that '%s' is a nexus or gde file..." % fName)
        if firstChar in ['[', '#']:
            # it might be a nexus file
            if var.verboseRead:
                print("Guessing that '%s' is a nexus file..." % fName)
            ret = _tryToReadNexusFile(fName, flob)
            if ret:
                return
        elif firstChar == '{':
            if var.verboseRead:
                print("Guessing that '%s' is a gde file..." % fName)
            ret = _tryToReadGdeFile(fName, flob)
            if ret:
                return

        # If we are here then it isn't a nexus or gde file.  Previously I then
        # tried to read it as a python script, even though it does not end in
        # '.py'.  Probably a bad idea.  So just give up.

        if var.verboseRead:
            print("Giving up on trying to read '%s'." % fName)

        if fName == '<input string>':
            flob.seek(0)
            first100 = flob.read(100)
            if first100:
                gm = ["Couldn't make sense out of the input '%s'" % first100]
            else:
                gm = ["Couldn't make sense out of the input '%s'" % fName]
            gm.append("It is not a file name, and I could not make sense out of it otherwise.")
        else:
            gm = ["Couldn't make sense of the input '%s'." % fName]
        flob.close()
        raise P4Error(gm)


def _tryToReadNexusFile(fName, flob):
    if var.verboseRead:
        print("Trying to read '%s' as a nexus file..." % fName)
    nf = Nexus()

    # nf.readNexusFile()
    #    returns -1 if it does not start with a #nexus token
    #    returns 1 otherwise
    ret = nf.readNexusFile(flob)
    # print "Nexus.readNexusFile() returned a %s" % ret
    if ret == -1:
        if var.verboseRead:
            print("Failed to get '%s' as a nexus file." % fName)
    else:
        if 0:
            print("\n******************************************\n")
            nf.dump()
            print("\n******************************************\n")
        if len(nf.alignments):
            for a in nf.alignments:
                a.checkLengthsAndTypes()
                if var.doCheckForAllGapColumns:
                    a.checkForAllGapColumns()
                if var.doCheckForBlankSequences:
                    a.checkForBlankSequences()
                if var.doCheckForDuplicateSequences:
                    a.checkForDuplicateSequences()
                var.alignments.append(a)
            nf.alignments = []
        if len(nf.trees):
            for t in nf.trees:
                t.checkDupedTaxonNames()
                if t.taxNames:
                    t.checkTaxNames()
            var.trees += nf.trees
            nf.trees = []

        # check for duplicate tree names (turned off)
        if 0:
            if len(var.trees):
                tnList = []
                for t in var.trees:
                    tnList.append(t.name.lower())
                for t in var.trees:
                    if tnList.count(t.name.lower()) > 1:
                        print("P4: Warning: got duplicated tree name '%s' (names are compared in lowercase)" \
                              % t.name)
                        # print "Lowercased tree names: %s" % tnList
                        # sys.exit()
        if hasattr(flob, 'name'):
            var.fileNames.append(flob.name)
        if var.verboseRead:
            print("Got nexus file '%s'" % fName)
        return 1


def _tryToReadFastaFile(fName, flob, firstLine=None):
    if not firstLine:
        firstLine = False
        while not firstLine:
            firstLine = flob.readline()
            if not firstLine:  # end of the file
                break
            firstLine = firstLine.strip()
            if firstLine:  # got some content
                break
            else:
                # print "blank line"
                pass
        if 0:
            if var.verboseRead:
                print("got firstLine = %s" % firstLine)
            else:
                print("Got a blank file.")
    if not firstLine:
        gm = ["_tryToReadFastaFile: the file '%s' is empty!" % fName]
        raise P4Error(gm)
    else:
        if var.verboseRead:
            print("Trying to read '%s' as a fasta file..." % fName)
        if len(firstLine) <= 1:
            if var.verboseRead:
                print("First line is blank--- not fasta")
            return
        if firstLine[0] not in ">;":
            if var.verboseRead:
                print("First char is neither '>' nor ';' ---not fasta")
            return

        flob.seek(0)
        sl = p4.sequencelist.SequenceList(flob)   # this parses the file contents
        if hasattr(flob, 'name'):
            sl.fName = flob.name
            var.fileNames.append(flob.name)
        sl.checkNamesForDupes()

        # If we have equal sequence lengths, then it might be an
        # alignment
        hasEqualSequenceLens = True
        if len(sl.sequences) <= 1:
            hasEqualSequenceLens = None  # ie not applicable
        else:
            len0 = len(sl.sequences[0].sequence)
            for s in sl.sequences[1:]:
                if len(s.sequence) != len0:
                    hasEqualSequenceLens = False

        if not hasEqualSequenceLens:
            if var.verboseRead:
                print("The sequences appear to be different lengths")
            var.sequenceLists.append(sl)
        else:
            if var.verboseRead:
                print("The sequences appear to be all the same length")
            try:
                # includes a call to checkLengthsAndTypes()
                a = sl.alignment()
                # a.checkLengthsAndTypes()
            except:
                if var.verboseRead:
                    print("Its not an alignment, even tho the sequences are all the same length.")
                    print("    Maybe p4 (erroneously?) thinks that the sequences are different dataTypes.")
                var.sequenceLists.append(sl)
                if var.verboseRead:
                    print("Got fasta file '%s'." % fName)
                return 1

            if var.verboseRead:
                print("The fasta file appears to be an alignment.")

            if var.doCheckForAllGapColumns:
                a.checkForAllGapColumns()
            if var.doCheckForBlankSequences:
                a.checkForBlankSequences()
            if var.doCheckForDuplicateSequences:
                a.checkForDuplicateSequences()
            var.alignments.append(a)

        if var.verboseRead:
            print("Got fasta file '%s'." % fName)
        return 1


def _tryToReadPhylipFile(fName, flob, firstLine):
    # print "tryToReadPhylipFile here"
    # print "firstLine is '%s'" % firstLine
    gm = ["func._tryToReadPhylipFile()"]
    if not firstLine:
        firstLine = flob.readline()
    # print "B firstLine is '%s'" % firstLine
    if not firstLine:
        gm.append("The file %s is empty." % fName)
        raise P4Error(gm)
    splitLine = firstLine.split()

    # If theres 2 numbers on the first line, it may be a phylip data file
    if len(splitLine) >= 2:
        try:
            firstNum = int(splitLine[0])
            secondNum = int(splitLine[1])

            if var.verboseRead:
                print("Trying to read '%s' as a phylip data file..." % fName)
            a = Alignment()
            if hasattr(flob, 'name'):
                a.fName = flob.name
            a.readOpenPhylipFile(flob, firstNum, secondNum)
            a.checkNamesForDupes()
            a.checkLengthsAndTypes()
            if var.doCheckForAllGapColumns:
                a.checkForAllGapColumns()
            if var.doCheckForBlankSequences:
                a.checkForBlankSequences()
            if var.doCheckForDuplicateSequences:
                a.checkForDuplicateSequences()
            var.alignments.append(a)
            var.fileNames.append(fName)
            if var.verboseRead:
                print("Got '%s' as a phylip-like file." % fName)
            return 1
        except ValueError:
            if var.verboseRead:
                print("Does not seem to be a phylip or phylip-like data file.")

    # Ok, so the file did not start with 2 integers.  It might still
    # be a Phylip tree file.
    flob.seek(0, 0)  # Go to the beginning of the flob
    firstChar = flob.read(1)
    while firstChar and firstChar in string.whitespace:
        firstChar = flob.read(1)
    if not firstChar:
        gm.append("No non-whitespace chars found.")
        raise P4Error(gm)
    # print "got firstChar '%s'" % firstChar
    if firstChar not in ['(', '[']:  # The '[' for puzzle output.
        if var.verboseRead:
            print("First char is not '(' or '['.  This does not seem to be a phylip or puzzle tree file.")
        return
    if var.verboseRead:
        print("Trying to read '%s' as a phylip tree file..." % fName)
    flob.seek(0, 0)
    theseTrees = []

    while 1:
        savedPosition = flob.tell()
        # Just to check whether there is a token that can be read ...
        tok = nextTok(flob)
        if not tok:
            break
        if tok != '(':
            if var.verboseRead:
                print("First char was '%s'," % firstChar, end=' ')
                print("so I thought it was a phylip or puzzle tree file.")
                print("However, after having read in %i trees," % len(theseTrees))
                print(" it confused me by starting a supposed new tree with a '%s'" % tok)
            return
        flob.seek(savedPosition, 0)  # Throw the token away.
        t = Tree()
        t.name = 't%i' % len(theseTrees)
        t.parseNewick(flob, None)  # None is the translationHash
        t._initFinish()
        theseTrees.append(t)
    if len(theseTrees) == 0:
        return
    else:
        for t in theseTrees:
            t.checkDupedTaxonNames()
            # if t.taxNames:  # Why would it?
            #    t.checkTaxNames()
        var.trees += theseTrees
    if hasattr(flob, 'name'):
        var.fileNames.append(flob.name)
    if var.verboseRead:
        print("Got %i trees from phylip tree file '%s'" % (len(theseTrees), fName))
    return 1


def _tryToReadClustalwFile(fName, flob, firstLine=None):
    if not firstLine:
        firstLine = flob.readline()
    if not firstLine:
        gm.append(
            "func. _tryToReadClustalwFile()  The file %s is empty." % fName)
        raise P4Error(gm)
    expectedFirstLine = 'CLUSTAL'
    if firstLine.startswith(expectedFirstLine):
        if var.verboseRead:
            print("Trying to read '%s' as a clustalw file..." % fName)
        a = Alignment()
        if hasattr(flob, 'name'):
            a.fName = flob.name
            var.fileNames.append(flob.name)
        a._readOpenClustalwFile(flob)
        a.checkNamesForDupes()
        a.checkLengthsAndTypes()
        if var.doCheckForAllGapColumns:
            a.checkForAllGapColumns()
        if var.doCheckForBlankSequences:
            a.checkForBlankSequences()
        if var.doCheckForDuplicateSequences:
            a.checkForDuplicateSequences()
        var.alignments.append(a)
        if var.verboseRead:
            print("Got '%s' as a clustalw file." % fName)
        return 1


def _tryToReadGdeFile(fName, flob):
    if var.verboseRead:
        print("Trying to read '%s' as a gde file..." % fName)
    a = Alignment()
    if hasattr(flob, 'name'):
        a.fName = flob.name
        var.fileNames.append(flob.name)
    a._readOpenGdeFile(flob)
    # a.writePhylip()
    if var.doCheckForAllGapColumns:
        a.checkForAllGapColumns()
    if var.doCheckForBlankSequences:
        a.checkForBlankSequences()
    if var.doCheckForDuplicateSequences:
        a.checkForDuplicateSequences()
    var.alignments.append(a)
    if var.verboseRead:
        print("Got '%s' as a gde file." % fName)
    return 1


def _tryToReadPirFile(fName, flob):
    if var.verboseRead:
        print("Trying to read '%s' as a pir file..." % fName)
    flob.seek(0)
    sl = p4.sequencelist.SequenceList()
    ret = sl._readOpenPirFile(flob)
    if not ret:
        if var.verboseRead:
            print("Reading it as a pir file didn't work.")
            return None
    else:
        # print "Got sl"
        if hasattr(flob, 'name'):
            sl.fName = flob.name
            var.fileNames.append(flob.name)
        sl.checkNamesForDupes()

        # If we have equal sequence lengths, then it might be an alignment
        hasEqualSequenceLens = True
        if len(sl.sequences) <= 1:
            hasEqualSequenceLens = None  # ie not applicable
        else:
            len0 = len(sl.sequences[0].sequence)
            for s in sl.sequences[1:]:
                if len(s.sequence) != len0:
                    hasEqualSequenceLens = False

        if not hasEqualSequenceLens:
            if var.verboseRead:
                print("The sequences appear to be different lengths")
            var.sequenceLists.append(sl)
        else:
            if var.verboseRead:
                print("The sequences appear to be all the same length")
            try:
                # includes a call to checkLengthsAndTypes()
                a = sl.alignment()
            except:
                if var.verboseRead:
                    print("Its not an alignment, even tho the sequences are all the same length.")
                    print("    Maybe p4 (erroneously?) thinks that the sequences are different dataTypes.")
                var.sequenceLists.append(sl)
                if var.verboseRead:
                    print("Got pir file '%s'." % fName)
                return 1

            if var.verboseRead:
                print("The pir file appears to be an alignment.")

            if var.doCheckForAllGapColumns:
                a.checkForAllGapColumns()
            if var.doCheckForBlankSequences:
                a.checkForBlankSequences()
            if var.doCheckForDuplicateSequences:
                a.checkForDuplicateSequences()
            var.alignments.append(a)

        if var.verboseRead:
            print("Got pir file '%s'." % fName)
        return 1


def splash():
    """Print a splash screen for p4."""
    print('')

    # from p4.version import versionString, dateString
    # print("p4 v %s, %s" % (versionString, dateString))
    
    print("""
usage:
    p4
 or
    p4 [-i] [-x] [-d] [yourScriptOrDataFile] [anotherScriptOrDataFile ...]
 or
    p4 --help

p4 is a Python package for phylogenetics.
p4 is also the name of a Python script that loads the p4 package.""")

    print("""
There is documentation at http://p4.nhm.ac.uk """)

    print("""
Using the p4 script, after reading in the (optional) files on the
command line, p4 goes interactive unless one of the files on the
command line is a Python script.  Use the -i option if you want to go
interactive even if you are running a script.  Use the -x option to
force exit, even if there was no Python script read.  If you use the
-d option, then p4 draws any trees that are read in on the command
line, and then exits.

Peter Foster
The Natural History Museum, London
p.foster@nhm.ac.uk""")

    if var.examplesDir:
        print("\nSee the examples in %s" % var.examplesDir)
    print('')
    print("(Control-d to quit.)\n")

def splash2(outFile=None, verbose=True):
    """Another splash, showing things like version, git hash, and date

    If verbose is set, it gets printed to sys.stdout.

    If you set an outFile, it will also be appended to that file.

    It also returns the info as a list of strings.
    """

    # Collect all the info in a list of strings
    stuff = []

    # Stolen from Cymon.  Thanks!
    #stuff.append("\nSummary from func.splash2()")
    #stuff.append("%16s: %s" % ("P4 version", p4.version.versionString))
    lp = os.path.dirname(inspect.getfile(p4))
    stuff.append("%16s: %s" % ("Library path", lp))

    # Get git version.
    if os.path.isdir(os.path.join(os.path.dirname(lp), '.git')):
        try:
            # I got these from https://stackoverflow.com/questions/14989858/get-the-current-git-hash-in-a-python-script
            # subprocess.check_output(['git', 'rev-parse', 'HEAD'])
            # ret = subprocess.check_output(['git', '-C', '%s' % lp, 'rev-parse', '--short', 'HEAD'])
            ret = subprocess.check_output(['git', '-C', '%s' % lp, 'log', '-1', '--date=short', '--pretty=format:"%h -- %cd -- %cr"'])
            #ret = ret.strip()    # get rid of newline, needed for rev-parse
            ret = ret[1:-1]       # get rid of quotes, needed for log
            stuff.append("%16s: %s" % ("git hash", ret))

        except subprocess.CalledProcessError:
            #print("%16s: %s" % ("git hash", "Not a git repo?"))
            pass
    else:
        stuff.append("%16s: %s" % ("git hash", "Not a git repo"))


    stuff.append("%16s: %s" % ("Python version", ".".join([str(i) for i in sys.version_info[:-2]])))
    #print("%16s: %s" % ("Date" , datetime.datetime.now().strftime("%d/%m/%Y")))
    stuff.append("%16s: %s" % ("Today's date" , datetime.datetime.now().strftime("%Y-%m-%d")))  # iso 8601 see https://xkcd.com/1179/
    host = os.uname()[1].split('.')[0]
    stuff.append("%16s: %s" % ("Host", host))
    #stuff.append("\n")

    if outFile:
        print("Appending splash2 info to file %s" % outFile) 
        fh = open(outFile, "a")
        for aLine in stuff:
            print(aLine, file=fh)
        fh.close()
    if verbose:
        for aLine in stuff:
            print(aLine)
    return stuff


def randomTree(taxNames=None, nTax=None, name='random', seed=None, biRoot=0, randomBrLens=1, constraints=None):
    """Make a simple random Tree.

    You can supply a list of taxNames, or simply specify nTax.  In the
    latter case the specified number of (boringly-named) leaves will
    be made.

    The default is to have 'randomBrLens', where internal nodes get
    brLens of 0.02 - 0.05, and terminal nodes get brLens of 0.02 -
    0.5. Branch lengths are all 0.1 if randomBrLens is turned off.

    This method starts with a star tree and keeps adding nodes until
    it is fully resolved.  If 'biRoot' is set, it adds one more node,
    and roots on that node, to make a bifurcating root.

    Repeated calls will give different random trees, without having to
    do any seed setting.  If for some reason you want to make
    identical random trees, set the seed to some positive integer, or
    zero.

    If you want the tree to have some topological constraints, say so
    with a Constraints object.  That can include a constrained root
    position.

    Returns a tree.
    """

    complaintHead = '\nrandomTree()'
    gm = [complaintHead]

    # we need either taxNames or nTax
    if not taxNames and not nTax:
        gm.append("You need to supply either taxNames or nTax.")
        raise P4Error(gm)
    if taxNames and nTax:
        if len(taxNames) != nTax:
            gm.append("You need not supply both taxNames and nTax,")
            gm.append("but if you do, at least they should match, ok?")
            raise P4Error(gm)
    if taxNames:  # implies not []
        nTax = len(taxNames)
    elif nTax:
        taxNames = []
        for i in range(nTax):
            taxNames.append('t%i' % i)

    if constraints:
        assert isinstance(constraints, Constraints)
        if constraints.rootConstraints:
            nCon = len(constraints.rootConstraints)
            if constraints.rTree.isBiRoot():
                if nCon != 1:
                    gm.append("biRoot trees with rootConstraints should")
                    gm.append("have exactly 1 rootConstraint split.")
                    gm.append(f"Got {nCon} rootConstraints {constraints.rootConstraints}")
                    raise P4Error(gm)
                if not biRoot:
                    gm.append("Arg biRoot is not turned on.")
                    gm.append("However, the constraint rTree is biRooted.")
                    raise P4Error(gm)
            elif constraints.rTree.isTriRoot: 
                if biRoot:
                    gm.append("Arg biRoot is turned on")
                    gm.append("However, the contraint rTree is tri-rooted.")
                    raise P4Error(gm)
            else:
                gm.append("The constraint rTree should be either biRoot or triRoot.")
                raise P4Error(gm)


    # Make a random list of indices for the taxNames
    if seed != None:  # it might be 0
        random.seed(seed)
    indcs = list(range(nTax))
    random.shuffle(indcs)

    # Instantiate the tree, and add a root node
    t = Tree()
    t._taxNames = taxNames
    n = Node()
    n.br = None
    t.name = name
    n.nodeNum = 0
    t.nodes.append(n)
    t.root = n

    # Add the left child of the root
    n = Node()
    n.isLeaf = 1
    t.root.leftChild = n
    t.nodes.append(n)
    n.nodeNum = 1
    n.name = taxNames[indcs[0]]
    n.parent = t.root
    previousNode = n

    # Add the rest of the terminal nodes
    nodeNum = 2
    for i in range(nTax)[1:]:
        n = Node()
        n.isLeaf = 1
        t.nodes.append(n)
        n.name = taxNames[indcs[nodeNum - 1]]
        previousNode.sibling = n
        previousNode = n
        n.parent = t.root
        n.nodeNum = nodeNum
        nodeNum += 1
    # t.dump(node=1)
    # t.draw()
    # constraints.dump()

    nNodesAddedForConstraints = 0
    if constraints:
        for aConstraint in constraints.constraints:
            # print("doing aConstraint %s  %i" % (getSplitStringFromKey(aConstraint, nTax), aConstraint))
            #t.dump(tree=0, node=1)
            t.setPreAndPostOrder()
            eTaxNames = []
            for i in range(nTax):
                tester = 1 << i
                # Does aConstraint contain the tester bit?
                if tester & aConstraint:
                    eTaxNames.append(taxNames[i])
            # print "aConstraint %s" % eTaxNames

            # check that they all share the same parent
            firstParent = t.node(eTaxNames[0]).parent
            for tN in eTaxNames[1:]:
                if t.node(tN).parent != firstParent:
                    gm.append("constraint %s" % getSplitStringFromKey(
                        aConstraint, constraints.tree.nTax))
                    gm.append("'%s' parent is not node %i" %
                              (tN, firstParent.nodeNum))
                    gm.append(
                        'It appears that there are incompatible constraints.')
                    raise P4Error(gm)

            n = Node()
            n.nodeNum = nodeNum
            nodeNum += 1
            chosenName = random.choice(eTaxNames)
            eTaxNames.remove(chosenName)
            # print 'adding a new parent for %s' % chosenName
            chosenNode = t.node(chosenName)
            chosenNodeOldSib = chosenNode.sibling
            chosenNodeOldLeftSib = chosenNode.leftSibling()  # could be None

            n.parent = firstParent
            n.leftChild = chosenNode
            n.sibling = chosenNodeOldSib
            chosenNode.parent = n
            if chosenNodeOldLeftSib:
                chosenNodeOldLeftSib.sibling = n
            else:
                firstParent.leftChild = n
            chosenNode.sibling = None
            t.nodes.append(n)
            nNodesAddedForConstraints += 1
            oldChosenNode = chosenNode

            if 0:
                t.preOrder = None
                t.postOrder = None
                t.preAndPostOrderAreValid = False
                t.draw()

            while eTaxNames:
                chosenName = random.choice(eTaxNames)
                # print "adding '%s'" % chosenName
                eTaxNames.remove(chosenName)
                chosenNode = t.node(chosenName)
                chosenNodeOldSib = chosenNode.sibling
                chosenNodeOldLeftSib = chosenNode.leftSibling()
                if 0:
                    if chosenNodeOldLeftSib:
                        print('chosenNodeOldLeftSib = %s' % chosenNodeOldLeftSib.nodeNum)
                    else:
                        print('chosenNodeOldLeftSib = None')
                    if chosenNodeOldSib:
                        print('chosenNodeOldSib = %s' % chosenNodeOldSib.nodeNum)
                    else:
                        print('chosenNodeOldSib = None')

                chosenNode.parent = n
                oldChosenNode.sibling = chosenNode
                if chosenNodeOldLeftSib:
                    chosenNodeOldLeftSib.sibling = chosenNodeOldSib
                else:
                    firstParent.leftChild = chosenNodeOldSib
                chosenNode.sibling = None
                oldChosenNode = chosenNode

                if 0:
                    t.preOrder = None
                    t.postOrder = None
                    t.preAndPostOrderAreValid = False
                    t.draw()

    if 0:
        t.preOrder = None
        t.postOrder = None
        t.preAndPostOrderAreValid = False
        t.draw()
    # sys.exit()

    # Now we have a star tree.  Now add internal nodes until it is all
    # resolved, which needs nTax - 3 nodes
    # print 'nNodesAddedForConstraints is %i' % nNodesAddedForConstraints
    for i in range(nTax - 3 - nNodesAddedForConstraints):
        # Pick a random node, that has a sibling, and a
        # sibling.sibling.  This week I am first making a list of
        # suitables and then choosing a random one, rather than first
        # choosing a random node and then asking whether it is
        # suitable.  It could be made more efficient by not re-making
        # the ssNodes list every time, but rather just keeping it up
        # to date.  But that is not so simple, so re-make the list
        # each time.
        ssNodes = []
        for n in t.nodes:
            if n.sibling and n.sibling.sibling:
                # This next thing can happen if there are constraints.
                # But it turns out that getNChildren() is very, very slow!  To be avoided!
                # if n.parent == t.root and t.root.getNChildren() == 3:
                #    pass
                if n.parent == t.root and t.root.leftChild.sibling and \
                        t.root.leftChild.sibling.sibling and not \
                        t.root.leftChild.sibling.sibling.sibling:
                    pass
                else:
                    ssNodes.append(n)
        lChild = random.choice(ssNodes)

        # print "lChild = node %i" % lChild.nodeNum

        # +----------1:oldLeftSib
        # |
        # +----------2:lChild
        # 0
        # +----------3:lChildSib
        # |
        # +----------4:oldLChildSibSib

        # +----------1:oldLeftSib
        # |
        # |          +----------3:lChild
        # 0----------2(n)
        # |          +----------4:lChildSib
        # |
        # +----------5:oldLChildSibSib

        n = Node()
        n.nodeNum = nodeNum
        nodeNum = nodeNum + 1
        lChildSib = lChild.sibling  # guarranteed to have one
        oldLChildSibSib = lChildSib.sibling  # ditto
        # oldLeftSib = lChild.parent.leftChild # first guess ...
        # if oldLeftSib != lChild:
        #    while oldLeftSib.sibling != lChild:
        #        oldLeftSib = oldLeftSib.sibling
        # else:
        #    oldLeftSib = None
        oldLeftSib = lChild.leftSibling()  # could be none
        if 0:
            if oldLeftSib:
                print("oldLeftSib = %i" % oldLeftSib.nodeNum)
            else:
                print("oldLeftSib = None")
            print("lChildSib = %i" % lChildSib.nodeNum)
            if oldLChildSibSib:
                print("oldLChildSibSib = %i" % oldLChildSibSib.nodeNum)
            else:
                print("oldLChildSibSib = None")

        if oldLeftSib:
            oldLeftSib.sibling = n
        else:
            lChild.parent.leftChild = n

        n.parent = lChild.parent
        lChild.parent = n
        n.leftChild = lChild
        lChildSib.parent = n
        n.sibling = oldLChildSibSib
        lChildSib.sibling = None
        t.nodes.append(n)

    if 0:
        print("before rooting", "=" * 50)
        t.preOrder = None
        t.postOrder = None
        t.preAndPostOrderAreValid = False
        t.makeSplitKeys()
        for n in t.iterInternalsNoRoot():
            n.name = f"{n.br.splitKey}"
        t.draw()

    if biRoot:
        # addNodeBetweenNodes() requires t.preOrder.
        t.preOrder = numpy.array([var.NO_ORDER] * len(t.nodes), numpy.int32)
        t.postOrder = numpy.array([var.NO_ORDER] * len(t.nodes), numpy.int32)
        if len(t.nodes) > 1:
            t.setPreAndPostOrder()
        if constraints and constraints.rootConstraints:
            # We have checked that we have a valid rootConstraints, above
            theKey = constraints.rootConstraints[0]
            t.makeSplitKeys()
            rNode = t.nodeForSplitKeyDict.get(theKey)
            assert rNode, "biRoot root constraint node not found"
            newNode = t.addNodeBetweenNodes(rNode, rNode.parent)
            t.reRoot(newNode, moveInternalName=False)
        else:
            # pick a random node, with a parent
            n = t.nodes[random.randrange(1, len(t.nodes))]
            newNode = t.addNodeBetweenNodes(n, n.parent)
            t.reRoot(newNode, moveInternalName=False)
    else:
        if constraints and constraints.rootConstraints:
            print(f"constraints.rootConstraints are {constraints.rootConstraints}")
            print("Below is the rootConstraints tree")
            constraints.rTree.draw()

            # Re-root to the first leaf, so we can post-order traverse
            # without bothering with the root
            t.reRoot(t.node(t.taxNames[0]))
            t.makeSplitKeys()
            # print("Below is the re-rooted tree I")
            # for n in t.iterInternalsNoRoot():
            #     n.name = f"{n.br.splitKey}"
            # t.draw()

            for n in t.iterPostOrder():
                # print(f"checking node {n.nodeNum}, with splitKey {n.br.splitKey}")
                if n.br.splitKey in constraints.rootConstraints:
                    t.reRoot(n.parent, moveInternalName=False)
                    break

            # print("Below is the re-rooted tree II")
            # t.draw()

            # Check
            rootSplits = [n.br.splitKey for n in t.root.iterChildren()]
            # print(f"rootSplits are {rootSplits}")
            isBad = False
            for sk in constraints.rootConstraints:
                if sk not in rootSplits:
                    gm.append(f"root constraint split {sk} is not a root split")
                    isBad = True
            if isBad:
                gm.append(f"Something is wrong. root splits are {rootSplits}")
                gm.append(f"root constraints are {constraints.rootConstraints}")
                gm.append(f"Maybe incompatible root contraints?")
                raise P4Error(gm)
        else:
            # The way it is now, the root rightmost child is always a
            # leaf.  Not really random, then, right?  So choose a random
            # internal node, and re-root it there.
            # print "nTax=%i, len(t.nodes)=%i" % (nTax, len(t.nodes))
            if nTax > 3:
                n = t.nodes[random.randrange(nTax + 1, len(t.nodes))]
                t.reRoot(n, moveInternalName=False)

    # The default is to have randomBrLens, where internal nodes get
    # brLens of 0.02 - 0.05, and terminal nodes get brLens of 0.2 -
    # 0.5.  Branch lengths are all 0.1 if randomBrLens is turned
    # off.
    if randomBrLens:
        for n in t.nodes:
            if n != t.root:
                if n.isLeaf:
                    n.br.len = 0.05 + (random.random() * 0.45)
                else:
                    n.br.len = 0.02 + (random.random() * 0.03)

    # _initFinish() checks that leaf nodes have names, and sets preAndPostOrder
    t._initFinish()

    # Added Sept 2020, because NDCH2 wants nodeNums in order.
    newNodes = []
    for i,n in enumerate(t.iterPreOrder()):
        n.nodeNum = i
        newNodes.append(n)
    t.nodes = newNodes
    t.setPreAndPostOrder()

    return t


def newEmptyAlignment(dataType=None, symbols=None, taxNames=None, length=None):
    """Make de novo and return an Alignment object, made of gaps.

    It is not placed in var.alignments.
    """

    complaintHead = ['\nnewEmptyAlignment()']
    gm = complaintHead
    # check for silliness
    if not dataType:
        gm.append(
            "No dataType. You need to specify at least the dataType, taxNames, and sequenceLength.")
        raise P4Error(gm)
    if not taxNames:
        gm.append(
            "No taxNames. You need to specify at least the dataType, taxNames, and sequenceLength.")
        raise P4Error(gm)
    if not length:
        gm.append(
            "No length.  You need to specify at least the dataType, taxNames, and sequenceLength.")
        raise P4Error(gm)
    goodDataTypes = ['dna', 'protein', 'standard']
    if dataType not in goodDataTypes:
        gm.append("dataType '%s' is not recognized.")
        gm.append("I only know about %s" % goodDataTypes)
        raise P4Error(gm)
    if dataType == 'standard':
        if not symbols:
            gm.append("For standard dataType you need to specify symbols.")
            raise P4Error(gm)
    else:
        if symbols:
            gm.append(
                "You should not specify symbols for %s dataType." % dataType)
            raise P4Error(gm)

    a = Alignment()
    a.length = length
    a.dataType = dataType

    # dataTypes, symbols, dim
    if a.dataType == 'dna':
        a.symbols = 'acgt'
        a.dim = 4
        a.equates = {'n': 'acgt', 'm': 'ac', 'k': 'gt',  # 'x': 'acgt',
                     'h': 'act', 'y': 'ct', 'v': 'acg',
                     'w': 'at', 'd': 'agt', 'b': 'cgt',
                     'r': 'ag', 's': 'cg'}
    elif a.dataType == 'protein':
        a.symbols = 'arndcqeghilkmfpstwyv'
        a.dim = 20
        a.equates = {'b': 'dn', 'x': 'arndcqeghilkmfpstwyv', 'z': 'eq'}
    elif a.dataType == 'standard':
        a.symbols = symbols
        a.dim = len(symbols)

    # Make sequences, composed of gaps.
    for n in taxNames:
        s = p4.sequencelist.Sequence()
        s.name = n
        s.dataType = a.dataType
        s.sequence = '-' * a.length
        a.sequences.append(s)

    return a


def getSplitStringFromKey(theKey, nTax, escaped=False):
    """Convert a long int binary split key to dot-star notation."""

    ss = ['.'] * nTax
    #ss = ['0'] * nTax
    for i in range(nTax):
        tester = 2 ** i
        if tester & theKey:
            ss[i] = '*'
            #ss[i] = '1'
    if escaped:
        return '\\' + '\\'.join(ss)
    else:
        return ''.join(ss)


def getSplitKeyFromTaxNames(allTaxNames, someTaxNames):
    """Make a long int binary split key from a list of taxNames.

    allTaxNames  -> an ordered list, nTax long
    someTaxNames -> list of taxNames on one side of the split.

    The split key that is returned will always be even.
    For example, assuming ::

        allTaxNames ['A', 'B', 'C', 'D']
        someTaxNames ['B', 'D']

    The bits for all the taxa are::

        A   B   C   D
        1   2   4   8

    So  the split key for ``['B', 'D']`` will be 2 + 8 = 10

    Another ::

        func.getSplitKeyFromTaxNames(['A', 'B', 'C', 'D'], ['B', 'D'])
        # returns 10

    However, if the splitKey is odd, it is bit-flipped.  So if
    ``someTaxNames = ['A', 'D']``, the raw split key for ``['A','D']``
    will be 1 + 8 = 9, binary '1001', which is then xor'd with
    1111, giving 6.

    Another ::

        getSplitKeyFromTaxNames(['A', 'B', 'C', 'D'], ['A', 'D'])
        returns 6

    """

    gm = ['func.getSplitKeyFromTaxNames()']
    if not len(allTaxNames) or not len(someTaxNames):
        gm.append("Got an empty arg?!?")
        raise P4Error(gm)
    theIndices = []
    for tn in someTaxNames:
        try:
            theIndex = allTaxNames.index(tn)
        except ValueError:
            gm.append("The taxName '%s' is not in allTaxNames." % tn)
            raise P4Error(gm)
        if theIndex not in theIndices:  # Duped indices would be a Bad Thing
            theIndices.append(theIndex)
    # print "theIndices = %s" % theIndices

    theRawSplitKey = 0
    for i in theIndices:
        theRawSplitKey += 1 << i  # "<<" is left-shift

    if 1 & theRawSplitKey:  # Is it odd?  or Does it contain a 1?
        allOnes = 2 ** (len(allTaxNames)) - 1
        theSplitKey = allOnes ^ theRawSplitKey  # "^" is xor, a bit-flipper.
        return theSplitKey
    else:
        return theRawSplitKey


def _sumOfRows(aList):
    """
    Adds up the rows of a 2d matrix, returning the vector.
    Eg _sumOfRows([[2,3], [6,13]]) returns [5, 19]
    """
    if not isinstance(aList[0], list):
        print("_sumOfRows: not a 2D array.  Assume its a row vector and return sum")
        return sum(aList)
    outList = []
    for i in aList:
        outList.append(sum(i))
    return outList


def _sumOfColumns(aList):
    """
    Adds up the rows of a 2d matrix, returning the vector.
    Eg _sumOfColumns([[2,3], [6,13]]) returns [8, 16]
    """
    if not isinstance(aList[0], list):
        print("_sumOfColumns: not a 2D array.  Assume its a column vector and return sum")
        return sum(aList)
    theLen = len(aList[0])
    for i in aList:
        if theLen != len(i):
            print("_sumOfColumns: unequal rows")
            return None
    outList = [0] * theLen
    for i in aList:
        for j in range(theLen):
            outList[j] = outList[j] + i[j]
    return outList


def _expected(sor, soc):  # sumOfRows, sumOfCols
    nRows = len(sor)
    nCols = len(soc)
    grand = float(sum(sor))
    expectedOut = []
    for i in range(nRows):
        outRow = []
        for j in range(nCols):
            outRow.append((sor[i] / grand) * soc[j])
        expectedOut.append(outRow)
    return expectedOut


def xSquared(observed):
    """Calculate the X^2 statistic from an R x C table.

    Arg observed is a 2D R x C table.  This stat is sometimes called
    Chi-squared.  """

    gm = ["func.xSquared()"]
    nRows = len(observed)
    nCols = len(observed[0])
    theSumOfRows = _sumOfRows(observed)
    theSumOfCols = _sumOfColumns(observed)
    theExpected = _expected(theSumOfRows, theSumOfCols)
    # print theExpected
    for i in theSumOfRows:
        if i == 0.0:
            gm.append(
                "Sum of rows includes a zero.  Can't calculate xSquared.")
            raise P4Error(gm)
    for i in theSumOfCols:
        if i == 0.0:
            gm.append(
                "Sum of cols includes a zero.  Can't calculate xSquared.")
            raise P4Error(gm)

    xSq = 0.0
    for i in range(nRows):
        for j in range(nCols):
            xSq = xSq + ((observed[i][j] - theExpected[i][j]) *
                         (observed[i][j] - theExpected[i][j]) / theExpected[i][j])
    return xSq


def variance(seq):
    """This would not be good for a lot of data. n - 1 weighted."""
    sumSeq = float(sum(seq))
    return (_sumOfSquares(seq) - ((sumSeq * sumSeq) / len(seq))) / (len(seq) - 1)
    # return (_sumOfSquares(seq) - ((sumSeq * sumSeq) / len(seq))) / len(seq)


def _stdErrorOfTheDifferenceBetweenTwoMeans(seq1, seq2):
    """This could use some re-coding to handle short (<30) n"""
    if len(seq1) == len(seq2) and len(seq1) > 30:
        return math.sqrt((variance(seq1) + variance(seq2)) / len(seq1))
    else:
        gm = ["_stdErrorOfTheDifferenceBetweenTwoMeans()"]
        gm.append(
            "I can only deal with sequences of equal length, each more than 30 long.")
        raise P4Error(gm)


def mean(seq):
    """Simple, pure-python mean.  For big lists, use something better."""

    return float(sum(seq)) / len(seq)


def studentsTStat(seq1, seq2):
    """Returns Student's t statistic for 2 lists or tuples.

    Mean of seq1 - mean of seq2, divided by the
    _stdErrorOfTheDifferenceBetweenTwoMeans(seq1, seq2)
    """
    return (mean(seq1) - mean(seq2)) / _stdErrorOfTheDifferenceBetweenTwoMeans(seq1, seq2)


def tailAreaProbability(theStat, theDistribution, verbose=1):
    """Calculate the tail area probability of theStat.
    
    That is the number of items in theDistribution that are greater
    than or equal to theStat.  

    Arg theDistribution need not be sorted.

    This function is used by  :meth:`~p4.pnumbers.Numbers.tailAreaProbability()`.

    Args:
    theStat (float): the test quantity
    theDistribution (list of floats): The null distribution
    verbose (Bool): Whether to print expanded results to stdout
    
    Returns:
    list: [theStat, (theMin, theMax), tailAreaProbability]
    
    """

    gm = ["tailAreaProbability()"]
    theLen = len(theDistribution)
    for i in theDistribution:
        try:
            float(i)
        except TypeError:
            gm.append("Item '%s' from theDistribution does not seem to be a float." % i)
            raise P4Error(gm)
    try:
        float(theStat)
    except TypeError:
        gm.append("theStat '%s' does not seem to be a float." % theStat)
        raise P4Error(gm)

    hits = 0
    theMax = theDistribution[0]
    theMin = theDistribution[0]
    for x in theDistribution:
        if x < theMin:
            theMin = x
        elif x > theMax:
            theMax = x
        if x >= theStat:
            hits = hits + 1
    tap = float(hits) / float(theLen)
    if verbose:
        print("# The stat is %s" % theStat)
        print("# The distribution has %i items" % len(theDistribution))
        print("# The distribution goes from %s to %s" % (theMin, theMax))
        print("# Items in the distribution were >= theStat %i times." % hits)
        print("# The tail-area probability is %f" % tap)
    return [theStat, (theMin, theMax), tap]


def ls():
    """Like the shell ls
    """
    fList = os.listdir('.')
    fList.sort()
    for f in fList:
        print(f)


def which(what, verbose=0):
    """Asks if an auxiliary program is available.

    This uses the shell command 'which' to find whether a program (as
    given by the argument 'what') is in the path.  It returns 0 or 1.
    If verbose is turned on, it speaks the path, if it exists."""

    if not isinstance(what, str):
        raise P4Error("function which().  I was expecting a string argument.")
    f = os.popen('which %s 2> /dev/null' % what, 'r')
    aLine = f.readline()
    f.close()
    if aLine:
        aLine = aLine[:-1]
        # tcsh does this, but I have not tested this part.
        if aLine.endswith('Command not found.'):
            aLine = None
    if aLine:
        if verbose:
            print(aLine)
        return 1
    else:
        return 0


def which2(program):
    """Find an executable

    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def _writeRMatrixTupleToOpenFile(theTuple, dim, flob, offset=23):
    gm = ["func._writeRMatrixTupleToOpenFile()"]
    if dim < 2:
        gm.append("dim must be two or more for this to work.")
        raise P4Error(gm)

    isShort = False  # For backward compatibility
    if var.rMatrixNormalizeTo1:
        isShort = False
        if len(theTuple) != (((dim * dim) - dim) / 2):
            if len(theTuple) == ((((dim * dim) - dim) / 2) - 1):
                isShort = True
            else:
                gm.append("var.rMatrixNormalizeTo1 is %i" %
                          var.rMatrixNormalizeTo1)
                gm.append("The length of the tuple (%i) is " % len(theTuple))
                gm.append("incommensurate with the dim (%i)" % dim)
                gm.append("(should be %i)" % (((dim * dim) - dim) / 2))
                raise P4Error(gm)
    else:
        isShort = True
        if len(theTuple) != (((dim * dim) - dim) / 2) - 1:
            gm.append("var.rMatrixNormalizeTo1 is %i" %
                      var.rMatrixNormalizeTo1)
            gm.append("The length of the tuple (%i) is " % len(theTuple))
            gm.append("incommensurate with the dim (%i)" % dim)
            gm.append("(should be %i)" % ((((dim * dim) - dim) / 2) - 1))
            raise P4Error(gm)

    if dim == 3:
        flob.write('%s\n' % repr(theTuple))

    else:
        # if var.rMatrixNormalizeTo1:
        if not isShort:
            # print "theTuple=%s" % theTuple
            tuplePos = 0
            row = 0
            flob.write('(')
            for i in range((dim - row) - 1):
                flob.write('%10.6f, ' % theTuple[tuplePos])
                tuplePos += 1
            flob.write('\n')
            row += 1
            formatString = '%' + '%is' % offset
            while row < dim - 2:
                flob.write(formatString % ' ')
                flob.write(' ')
                for i in range(row):
                    flob.write('            ')
                for i in range((dim - row) - 1):
                    flob.write('%10.6f, ' % theTuple[tuplePos])
                    tuplePos = tuplePos + 1
                flob.write('\n')
                row += 1
            flob.write(formatString % ' ')
            flob.write(' ')
            for i in range(row):
                flob.write('            ')
            # print "dim=%i, row=%i, range((dim-row) - 1) = %s, tuplePos=%i" %
            # (dim, row, range((dim-row)-1), tuplePos)
            flob.write('%10.6f )' % theTuple[tuplePos])
            flob.write('\n')
        else:
            tuplePos = 0
            row = 0
            flob.write('(')
            for i in range((dim - row) - 1):
                flob.write('%10.6f, ' % theTuple[tuplePos])
                tuplePos += 1
            flob.write('\n')
            row += 1
            formatString = '%' + '%is' % offset
            while row < dim - 3:
                flob.write(formatString % ' ')
                flob.write(' ')
                for i in range(row):
                    flob.write('            ')
                for i in range((dim - row) - 1):
                    flob.write('%10.6f, ' % theTuple[tuplePos])
                    tuplePos = tuplePos + 1
                flob.write('\n')
                row += 1
            flob.write(formatString % ' ')
            flob.write(' ')
            for i in range(row):
                flob.write('            ')
            for i in range((dim - row) - 2):
                flob.write('%10.6f, ' % theTuple[tuplePos])
                tuplePos += 1
            flob.write('%10.6f )' % theTuple[tuplePos])
            flob.write('\n')


def _writeCharFreqToOpenFile(theCharFreq, dim, symbols, flob, offset=23):
    formatString = '%' + '%is' % offset
    s = 0.0
    if symbols:
        if dim == 2:
            flob.write('(')
            flob.write('     [%c] %f,)\n' % (symbols[0], theCharFreq[0]))
            flob.write(formatString % ' ')
            flob.write('[and [%c] %f, ' % (symbols[1], 1.0 - theCharFreq[0]))
        else:
            flob.write('(')
            flob.write('     [%c] %f,\n' % (symbols[0], theCharFreq[0]))
            s = theCharFreq[0]
            for i in range(dim - 2)[1:]:
                flob.write(formatString % ' ')
                flob.write('     [%c] %f,\n' % (symbols[i], theCharFreq[i]))
                s = s + theCharFreq[i]
            flob.write(formatString % ' ')
            flob.write('     [%c] %f)\n' %
                       (symbols[dim - 2], theCharFreq[dim - 2]))
            s = s + theCharFreq[dim - 2]
            flob.write(formatString % ' ')
            flob.write('[and [%c] %f, ' % (symbols[dim - 1], 1.0 - s))
    else:
        if dim == 2:
            flob.write('(')
            flob.write('%f,)\n' % theCharFreq[0])
            flob.write(formatString % ' ')
            flob.write('[and %f, ' % (1.0 - theCharFreq[0]))
        else:
            flob.write('(')
            flob.write('%f,\n' % theCharFreq[0])
            s = theCharFreq[0]
            for i in range(dim - 2)[1:]:
                flob.write(formatString % ' ')
                flob.write('%f,\n' % theCharFreq[i])
                s = s + theCharFreq[i]
            flob.write(formatString % ' ')
            flob.write('%f)\n' % theCharFreq[dim - 2])
            s = s + theCharFreq[dim - 2]
            flob.write(formatString % ' ')
            flob.write('[and %f, ' % (1.0 - s))


def fixCharsForLatex(theString):
    if not theString:
        return theString
    l = list(theString)
    # print l
    inMath = False
    for i in range(len(l)):
        if l[i] in string.ascii_letters or l[i] in string.digits:
            pass
        elif l[i] == '$':
            if inMath:
                inMath = False
            else:
                inMath = True
        elif l[i] in ['%', '#']:
            l[i] = '\\' + l[i]
        elif l[i] in ['_']:
            if inMath:
                pass
            else:
                l[i] = '\\' + l[i]
        elif l[i] in ['|']:
            if inMath:
                pass
            else:
                l[i] = '$' + l[i] + '$'
        else:
            pass
    if l[0] == "'" and l[-1] == "'":
        del(l[-1])
        del(l[0])
        theRange = range(len(l))
        for i in reversed(theRange[:-1]):
            if l[i] == "'" and l[i - 1] == "'":
                del(l[i])
    return ''.join(l)


def maskFromNexusCharacterList(nexusCharListString, maskLength, invert=0):
    """Returns a mask string, converted from a Nexus char list.

    Convert a Nexus characters list to a mask string composed of zeros
    and 1's Eg char list ``r'1 2-4 6-10'`` (don't forget the 'r' if you
    have any backslash characters) becomes (for a maskLength of 16)
    ``1111011111000000`` And ``r'10-.'`` results in ``0000000001111111`` (for
    maskLength 16).  Note that Nexus char lists are 1-based.

    Invert inverts the zeros and 1's.
    """

    gm = ["maskFromNexusCharacterList()"]
    from p4.alignment import cListPat, cList2Pat, cListAllPat
    #cListPat = re.compile('(\d+)-?(.+)?')
    #cList2Pat = re.compile('(.+)\\\\(\d+)')
    #cListAllPat = re.compile('all\\\\?(\d+)?')

    # print("char list is: %s" % nexusCharListString)
    cList = nexusCharListString.split()
    # print("cList is %s" % cList)
    mask = ['0'] * maskLength
    for c in cList:        # eg 6-10\2
        first = None       # the first item eg 6
        second = None      # the second item eg 10
        third = None       # the third item eg 2
        result = cListPat.match(c)
        if result:
            # print("%s\t%s" % (result.group(1), result.group(2)))
            first = result.group(1)
            if result.group(2):
                r2 = cList2Pat.match(result.group(2))
                if r2:
                    # print("%s\t%s" % (r2.group(1), r2.group(2)))
                    second = r2.group(1)
                    third = r2.group(2)
                else:
                    second = result.group(2)
        else:
            result = cListAllPat.match(c)
            if result:
                first = '1'
                second = '.'
                third = result.group(1)
            else:
                gm.append(
                    "Can't parse maskFromNexusCharacterList '%s'" % nexusCharListString)
                raise P4Error(gm)
        # print("first = %s, second = %s, third = %s" % (first, second, third))
        if not first:
            gm.append("Can't parse maskFromNexusCharacterList '%s'" %
                      nexusCharListString)
            raise P4Error(gm)
        elif first and not second:  # its a single
            if first.lower() == 'all':
                for i in range(len(mask)):
                    mask[i] = '1'
            elif first == '.':
                mask[-1] = '1'
            else:
                try:
                    it = int(first)
                    mask[it - 1] = '1'
                except ValueError:
                    gm.append("Can't parse '%s' in maskFromNexusCharacterList '%s'"
                              % (first, nexusCharListString))
                    raise P4Error(gm)
        elif first and second:  # its a range
            try:
                start = int(first)
            except ValueError:
                gm.append("Can't parse '%s' in maskFromNexusCharacterList '%s'"
                          % (first, nexusCharListString))
                raise P4Error(gm)
            if second == '.':
                fin = len(mask)
            else:
                try:
                    fin = int(second)
                except ValueError:
                    gm.append("Can't parse '%s' in maskFromNexusCharacterList '%s'" %
                              (second, nexusCharListString))
                    raise P4Error(gm)
            if third:
                try:
                    bystep = int(third)
                except ValueError:
                    gm.append("Can't parse '%s' in maskFromNexusCharacterList '%s'" %
                              (third, nexusCharListString))
                    raise P4Error(gm)
                for spot in range(start - 1, fin, bystep):
                    mask[spot] = '1'
            else:
                for spot in range(start - 1, fin):
                    mask[spot] = '1'
    if invert:
        for i in range(len(mask)):
            if mask[i] == '0':
                mask[i] = '1'
            elif mask[i] == '1':
                mask[i] = '0'
    return ''.join(mask)


def polar2square(angleLenList):
    """Convert a coord in polar coords to usual (Cartesian? square?) coords.

    Input is a list composed of the angle in radians and the length
    (ie from the origin).  A list of [x,y] is returned.
    """
    if angleLenList[1] == 0.0:
        return [0, 0]
    elif angleLenList[1] < 0.0:
        raise P4Error("func.polar2square error: len is less than zero")
    if angleLenList[0] == (math.pi / 2.0) or angleLenList[0] == (math.pi / -2.0):
        adj = 0.0
    else:
        adj = math.cos(angleLenList[0]) * angleLenList[1]
    if angleLenList[0] == math.pi or angleLenList[0] == -math.pi:
        opp = 0.0
    else:
        opp = math.sin(angleLenList[0]) * angleLenList[1]
    return [adj, opp]


def square2polar(xyList):
    """Convert usual square coords to polar.

    Arg is [x,y], and [angle, length] is returned.
    """

    angle = math.atan2(xyList[1], xyList[0])
    len = math.hypot(xyList[0], xyList[1])
    return [angle, len]


def factorial(n):
    """Return n!

    Its fast for n up to 30, cuz it uses a dictionary, rather than doing the
    computations.  Not needed for newer Pythons -- its in the math module."""

    try:
        n = int(n)
    except:
        raise P4Error("n should be (at least convertible to) an int.")
    assert n >= 0, "n should be zero or positive."

    fact = {0: 1, 1: 1, 2: 2, 3: 6, 4: 24, 5: 120, 6: 720, 7: 5040, 8: 40320, 9: 362880,
            10: 3628800, 11: 39916800, 12: 479001600, 13: 6227020800, 14: 87178291200,
            15: 1307674368000, 16: 20922789888000, 17: 355687428096000, 18: 6402373705728000,
            19: 121645100408832000, 20: 2432902008176640000, 21: 51090942171709440000,
            22: 1124000727777607680000, 23: 25852016738884976640000, 24: 620448401733239439360000,
            25: 15511210043330985984000000, 26: 403291461126605635584000000,
            27: 10888869450418352160768000000, 28: 304888344611713860501504000000,
            29: 8841761993739701954543616000000, 30: 265252859812191058636308480000000}
    if n <= 30:
        return fact[n]
    else:
        total = 1
        while n > 1:
            total *= n
            n -= 1
        return total


def nChooseK(n, k):
    """Get the number of all possible len k subsets from range(n)."""
    
    assert isinstance(n, int)
    assert isinstance(k, int)

    assert n >= 0, "n should be zero or more."
    assert k <= n, "k should be less than or equal to n"
    assert k >= 0, "k should be zero or more."

    nFact = factorial(n)
    kFact = factorial(k)
    nMinusKFact = factorial(n - k)
    return nFact / (kFact * nMinusKFact)


def nUnrootedTrees(nTaxa):
    upper = (nTaxa * 2) - 5
    nTrees = 1
    i = 3
    while i <= upper:
        nTrees *= i
        i += 2
    return nTrees


def nRootedTrees(nTaxa):
    upper = (nTaxa * 2) - 3
    nTrees = 1
    i = 3
    while i <= upper:
        nTrees *= i
        i += 2
    return nTrees


def nRootedTreesWithMultifurcations(nTaxa):
    """Returns a list T(n,m) with m from zero to n-1

    n is the number of leaves
    m is the number of internal nodes

    The first number in the returned list will always be zero, and the
    second number will always be 1.  See Felsenstein, page 27.  So for
    example, for nTaxa = 8 (as in the example), this function returns
    [0, 1, 246, 6825, 56980, 190575, 270270, 135135].  """

    # Make a table.
    table = []
    for nInts in range(nTaxa):
        table.append([0] * (nTaxa + 1))
    for nT in range(2, nTaxa + 1):
        table[1][nT] = 1
    for nInt in range(2, nTaxa):
        for nTx in range(nInt + 1, nTaxa + 1):
            table[nInt][nTx] = (
                (nTx + nInt - 2) * table[nInt - 1][nTx - 1]) + (nInt * table[nInt][nTx - 1])

    if 0:
        # Print out the table.
        print("%-9s|" % "nTx ->", end=' ')
        for nTx in range(1, nTaxa + 1):
            print("%10i" % nTx, end=' ')
        print()
        for nTx in range(nTaxa + 1):
            print("%10s" % " ---------", end=' ')
        print()
        print("%8s |" % "nInt")

        for nInt in range(1, nTaxa):
            print("%8i |" % nInt, end=' ')
            for nTx in range(nInt):
                print("%10s" % "", end=' ')
            for nTx in range(nInt + 1, nTaxa + 1):
                print("%10i" % table[nInt][nTx], end=' ')
            print()
    results = []
    for nInt in range(nTaxa):
        results.append(table[nInt][nTaxa])
    return results


def nUnrootedTreesWithMultifurcations(nTaxa):
    """Returns a list T(n,m) with m from zero to n-2

    n is the number of leaves
    m is the number of internal nodes

    The first number in the returned list will always be zero, and the
    second number will always be 1.  See Felsenstein, page 27.  So for
    example, for nTaxa = 9 (as in the example), this function returns
    [0, 1, 246, 6825, 56980, 190575, 270270, 135135].  """

    # From Felsenstein, page 27, 28.
    return nRootedTreesWithMultifurcations(nTaxa - 1)


def dirichlet1(inSeq, alpha, theMin, theMax=None, u=None):
    """Modify inSeq with a draw from a Dirichlet distribution with a single alpha value.

    *inSeq* is a list, and a copy is made, modified, normalized so that
    it sums to 1.0, and returned.

    This function uses the function random.gammavariate(x, 1.0), which
    takes draws from a gamma distribution (not gamma function).  Now
    random.gammavariate(x, 1.0) can return zero for x less than about
    0.001. Eg here are results for 1000 draws from different x values

    x=1.0000000000  min= 0.000319714   mean=   0.96978
    x=0.1000000000  min= 1.65643e-38   mean=   0.10074
    x=0.0100000000  min=4.03328e-309   mean=  0.013386
    x=0.0010000000  min=           0   mean= 0.0026256
    x=0.0001000000  min=           0   mean= 1.625e-09
    x=0.0000100000  min=           0   mean=8.8655e-15
    x=0.0000010000  min=           0   mean=1.0435e-268
    x=0.0000001000  min=           0   mean=         0

    Assuming we do not want zeros, a big alpha would help, but only by
    however big alpha is, so for example if alpha is 100 then we can
    still usually get a non-zero from inSeq vals of 0.0001.  One hack
    that is adopted here is to possibly add a constant u that might be
    0.1 or more to x (depending on alpha), so that the arg x is always
    big enough to return a non-zero.  Stick that in the arg u, which
    is by default None

    """

    gm = ['func.dirichlet1()']

    kk = len(inSeq)
    if theMax == None:
        theMax = 1.0 - ((kk - 1) * theMin)
    safety = 0
    safetyLimit = 300
    while 1:
        #theSeq = inSeq[:]
        if u:
            theSeq = [
                random.gammavariate((inSeq[i] * alpha) + u, 1.0) for i in range(kk)]
        else:
            theSeq = [
                random.gammavariate(inSeq[i] * alpha, 1.0) for i in range(kk)]
        # print safety, theSeq
        theSum = sum(theSeq)
        theSeq = [v / theSum for v in theSeq]
        thisMin = min(theSeq)
        thisMax = max(theSeq)
        isOk = True
        if (thisMin < theMin) or (thisMax > theMax):
            isOk = False
        if isOk:
            return theSeq
        safety += 1
        if safety > safetyLimit:
            gm.append(
                "Tried more than %i times to get good dirichlet values, and failed.  Giving up." % safetyLimit)
            gm.append("inSeq: %s" % inSeq)
            gm.append("theMin: %s, theMax: %s, u=%s" % (theMin, theMax, u))
            raise P4Error(gm)


def unPickleMcmc(runNum, theData, verbose=True):
    """Unpickle a checkpoint, return an Mcmc ready to go."""

    gm = ["func.unPickleMcmc()"]
    try:
        runNum = int(runNum)
    except (ValueError, TypeError):
        gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
        raise P4Error(gm)
    if runNum < 0:
        gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
        raise P4Error(gm)
    baseName = "mcmc_checkPoint_%i." % runNum
    ff = glob.glob("%s*" % baseName)

    pickleNums = []
    for f in ff:
        genNum = int(f.split('.')[1])
        pickleNums.append(genNum)
    if not pickleNums:  # an empty sequence
        gm.append("Can't find any checkpoints for runNum %i." % runNum)
        gm.append("Got the right runNum?")
        raise P4Error(gm)
    theIndx = pickleNums.index(max(pickleNums))
    fName = ff[theIndx]
    if verbose:
        print("...unpickling Mcmc in %s" % fName)

    f = open(fName, 'rb')
    m = pickle.load(f)
    f.close()

    # Restore gsl_rng state
    if not var.gsl_rng:
        var.gsl_rng = pf.gsl_rng_get()

    # accommodate old checkpoints that do not have this stuff
    if hasattr(m, "gsl_rng_state_ndarray"):
        the_gsl_rng_size = pf.gsl_rng_size(var.gsl_rng) # size of the state
        assert m.gsl_rng_state_ndarray.shape[0] == the_gsl_rng_size
        pf.gsl_rng_setstate(var.gsl_rng, m.gsl_rng_state_ndarray)
    else:
        # an old checkpoint with no random state; these steps are also done in Mcmc.__init__()
        pf.gsl_rng_set(var.gsl_rng, int(time.time()))
        the_gsl_rng_size = pf.gsl_rng_size(var.gsl_rng) # size of the state
        # A place to store the state.  Empty to start.  It is stored during a checkpoint.
        m.gsl_rng_state_ndarray = numpy.array(['0'] * the_gsl_rng_size, numpy.dtype('B'))  # B is unsigned byte

    if hasattr(m, 'randomState'):
        # restore random (python module) state
        assert m.randomState
        random.setstate(m.randomState)
    else:
        # an old checkpoint, that does not have this yet.
        m.randomState = None

    m.tree.data = theData
    m.tree.calcLogLike(verbose=False, resetEmpiricalComps=False)
    if m.simulate:
        m.simTree.data = theData.dupe()
        m.simTree.calcLogLike(verbose=False, resetEmpiricalComps=False)
    for chNum in range(m.nChains):
        ch = m.chains[chNum]
        ch.curTree.data = theData
        ch.curTree.calcLogLike(verbose=False, resetEmpiricalComps=False)
        ch.propTree.data = theData
        ch.propTree.calcLogLike(verbose=False, resetEmpiricalComps=False)

    m._setLogger()

    return m


def unPickleSTMcmc(runNum, verbose=True):
    """Unpickle a STMcmc checkpoint, return an STMcmc ready to go."""

    gm = ["func.unPickleSTMcmc()"]
    try:
        runNum = int(runNum)
    except (ValueError, TypeError):
        gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
        raise P4Error(gm)
    if runNum < 0:
        gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
        raise P4Error(gm)
    baseName = "mcmc_checkPoint_%i." % runNum
    ff = glob.glob("%s*" % baseName)

    pickleNums = []
    for f in ff:
        genNum = int(f.split('.')[1])
        pickleNums.append(genNum)
    if not pickleNums:  # an empty sequence
        gm.append("Can't find any checkpoints for runNum %i." % runNum)
        gm.append("Got the right runNum?")
        raise P4Error(gm)
    theIndx = pickleNums.index(max(pickleNums))
    fName = ff[theIndx]
    if verbose:
        print("...unpickling Mcmc in %s" % fName)

    f = open(fName, 'rb')
    m = pickle.load(f)
    f.close()

    if m.stRFCalc == 'fastReducedRF':
        import pyublas      # needed
        for chNum in range(m.nChains):
            ch = m.chains[chNum]
            ch.startFrrf()

    if m.modelName == 'SPA' and var.stmcmc_useFastSpa:
        import p4.fastspa as fastspa
        m.fspa = fastspa.FastSpa(m.useSplitSupport)
        for tNum, t in enumerate(m.trees):
            m.fspa.setInTr(tNum, t.nTax, m.nTax, t.baTaxBits.to01(), t.firstTax)
            for n in t.internals:
                if n.br and hasattr(n.br, "support"):
                    support = n.br.support
                else:
                    support = -1.0
                m.fspa.setInTrNo(tNum, n.stSplitKey.to01(), support)
        #m.fspa.summarizeInTrs()

        for chNum in range(m.nChains):
            ch = m.chains[chNum]
            ch.setupBitarrayCalcs()




            ch.getTreeLogLike_spa_bitarray()
            if var.stmcmc_useFastSpa:
                #print("Here E.  bitarray propTree.logLike is %f" % self.propTree.logLike)
                fspaLike = ch.stMcmc.fspa.calcLogLike(ch.chNum)
                diff = math.fabs(ch.propTree.logLike - fspaLike)
                #print("Got fspaLike %f, diff %g" % (fspaLike, diff))
                if diff > 1e-13:
                    gm.append("bad fastspa likelihood calc, %f vs %f, diff %f" % (ch.propTree.logLike, fspaLike, diff))
                    raise P4Error(gm)

            ch.curTree.logLike = ch.propTree.logLike
        

    m._setLogger()
    return m


#################################################

def recipes(writeToFile=var.recipesWriteToFile):
    """Reminders, suggestions, multi-step methods...

    These should be in the Sphinx docs also.
    """

    gm = ["func.recipes()"]
    if not var.examplesDir:
        gm.append("Can't find the Examples directory.")
        raise P4Error(gm)
    recipesDir = os.path.join(var.examplesDir, 'W_recipes')

    # # One line blurb, and default file name.
    # recipesList = [
    #     ["Calculate likelihood.", 'sLike.py'],
    #     ["Calculate likelihood with more than one data partition.", 'sLikeMultiPart.py'],
    #     ["Simulate data.", 'sSim.py'],
    #     ["Simulate data with more than one data partition.", 'sSimMultiPart.py'],
    #     ["Do an MCMC.", 'sMcmc.py'],
    #     ["Do an MCMC with more than one data partition.", 'sMcmcMultiPart.py'],
    #     ["Read checkPoints from an MCMC.", 'sReadMcmcCheckPoints.py'],
    #     ["Restart an MCMC.", 'sRestartMcmc.py'],
    #     ["Make a consensus tree.", 'sCon.py'],
    #     ]

    fList = glob.glob("%s/*.py" % recipesDir)
    # print fList
    bNames = [os.path.basename(nm) for nm in fList]
    # print bNames
    firstLines = []
    for fN in fList:
        f = open(fN)
        firstLine = f.readline()
        f.close()
        if firstLine.startswith('#'):
            firstLine = firstLine[1:]
            firstLine = firstLine.strip()
        else:
            firstLine = "Fix me!"
        firstLines.append(firstLine)
    recipesList = []
    for fNum in range(len(fList)):
        recipesList.append([firstLines[fNum], bNames[fNum]])
    # print recipesList

    for recNum in range(len(recipesList)):
        rec = recipesList[recNum]
        print("%s.  %s" % (string.ascii_uppercase[recNum], rec[0]))

    ret = input('Tell me a letter: ')
    
    # print "Got %s" % ret
    if ret == '':
        return
    elif ret[0] in string.ascii_uppercase:
        recNum = string.ascii_uppercase.index(ret[0])
    elif ret[0] in string.ascii_lowercase:
        recNum = string.ascii_lowercase.index(ret[0])
    else:
        return
    # print "Got recNum %i" % recNum

    if writeToFile:
        if sys.version_info < (3,):
            ret = raw_input('Write it to file name [default %s] :' % recipesList[recNum][1])
        else:
            ret = input('Write it to file name [default %s] :' % recipesList[recNum][1])

        if ret == '':
            theFName = recipesList[recNum][1]
        else:
            theFName = ret.strip()
        if os.path.exists(theFName):
            "The file %s already exists.  I'm refusing to over-write it, so I'm doing nothing." % theFName
            return
        print("Writing to file '%s' ..." % theFName)
        os.system("cp %s %s" %
                  (os.path.join(recipesDir, recipesList[recNum][1]), theFName))
    else:
        print("\n")
        os.system("cat %s" % os.path.join(recipesDir, recipesList[recNum][1]))


def uninstall():
    """Uninstall the p4 package."""

    print("""
This function is for uninstalling an installed version of p4.  It uses
the p4.installation module to get the locations of files.  It may need
to be done as root, or using sudo.""")

    weAreInteractive = False
    if os.getenv('PYTHONINSPECT'):   # The p4 script sets this
        weAreInteractive = True
    elif len(sys.argv) == 1 and sys.argv[0] == '':  # Command 'p4' with no args
        weAreInteractive = True
    # How can I tell if python was invoked with the -i flag?
    if not weAreInteractive:
        return
    try:
        import p4.installation
    except ImportError:
        raise P4Error("Unable to import the p4.installation module.")
    print("""
    
This function will remove the p4 script file
  %s
the library files in the directory
  %s
and the documentation in the directory
  %s

""" % (installation.p4ScriptPath, installation.p4LibDir, installation.p4DocDir))

    ret = input('Ok to do this? [y/n]')

    ret = ret.lower()
    if ret not in ['y', 'yes']:
        return
    print("Ok, deleting ...")
    if os.path.exists(installation.p4LibDir):
        os.system("rm -fr %s" % installation.p4LibDir)
    else:
        raise P4Error("Could not find %s" % installation.p4LibDir)
    if os.path.exists(installation.p4DocDir):
        os.system("rm -fr %s" % installation.p4DocDir)
    else:
        raise P4Error("Could not find %s" % installation.p4DocDir)
    if os.path.exists(installation.p4ScriptPath):
        os.system("rm %s" % installation.p4ScriptPath)
    else:
        raise P4Error("Could not find %s" % installation.p4ScriptPath)


###########################################################
# Tkinter stuff
###########################################################

# def startTkThread():
# if var.tk_thread_running:
# print "Tk thread is already running, it appears."
# return
####    import thread
####    import atexit
####    from Queue import Queue

####    var.tk_request = Queue(0)
####    var.tk_result = Queue(1)

# thread.start_new_thread(_tk_thread,())
####    var.tk_thread_running = True
# atexit.register(_tkShutdown)

# def _tk_thread():
##    import Tkinter
# print "_tk_thread() here!"
##    var.tk_root = Tkinter.Tk()
# print "var.tk_root is %s" % var.tk_root
# var.tk_root.withdraw()
##    var.tk_root.after(var.tk_pollInterval, _tk_pump)
# var.tk_root.mainloop()

# def _tk_pump():
# global _thread_running
# while not var.tk_request.empty():
##        command,returns_value = var.tk_request.get()
# try:
##            result = command()
# if returns_value:
# var.tk_result.put(result)
# except:
##            var.tk_thread_running = False
# if returns_value:
# var.tk_result.put(None) # release client
# raise # re-raise the exception -- kills the thread
# if var.tk_thread_running:
##        var.tk_root.after(var.tk_pollInterval, _tk_pump)


# def _tkShutdown():
# shutdown the tk thread
# global _thread_running
# _tkExec(sys.exit)
##    var.tk_thread_running = False
# time.sleep(.5) # give tk thread time to quit


##############################################################
##############################################################
##############################################################

def spaceDelimitedToTabDelimited(fName, outFName=None):
    """Convert space-delimited data files to tab-delimited.

    The outfilename by default is the infilename with .tabbed stuck on
    the end.
    """

    assert os.path.isfile(fName)

    if outFName:
        oFName = outFName
    else:
        oFName = "%s.tabbed" % fName

    f = open(fName)
    ll = f.readlines()
    f.close()

    f = open(oFName, 'w')
    for l in ll:
        sl = l.split()
        jl = '\t'.join(sl)
        f.write("%s\n" % jl)
    f.close()


def uniqueFile(file_name):
    """Returns an open file and filename, modified from file_name.

    If the file_name is foo.bar, then the new filename will start with
    foo and end with .bar, with a bit of unique nonsense in between.
    With a simple file_name input the new file is made in current
    directory, but by supplying a file_name including a path, you can
    create the new file elsewhere.

    Don't forget to close the file.
    """
    # I got this from the web somewhere.
    import tempfile
    dirname, filename = os.path.split(file_name)
    prefix, suffix = os.path.splitext(filename)

    fd, filename = tempfile.mkstemp(suffix, prefix + "_", dirname)
    return os.fdopen(fd, 'w'), filename


def writeInColour(theString, colour='blue'):
    goodColours = [
        'red', 'RED', 'blue', 'BLUE', 'cyan', 'CYAN', 'violet', 'VIOLET']
    if colour not in goodColours:
        raise P4Error(
            "func.printColour().  The colour should be one of %s" % goodColours)
    codeDict = {
        'red': '\033[0;31m',
        'RED': '\033[1;31m',
        'blue': '\033[0;34m',
        'BLUE': '\033[1;34m',
        'cyan': '\033[0;36m',
        'CYAN': '\033[1;36m',
        'violet': '\033[0;35m',
        'VIOLET': '\033[1;35m',
    }
    backToBlackCode = '\033[m'
    sys.stdout.write("%s%s%s" % (codeDict[colour], theString, backToBlackCode))


def setTerminalColour(theColour):
    goodTerminalColours = [
        'red', 'RED', 'blue', 'BLUE', 'cyan', 'CYAN', 'violet', 'VIOLET']
    terminalColourCodeDict = {
        'red': '\033[0;31m',
        'RED': '\033[1;31m',
        'blue': '\033[0;34m',
        'BLUE': '\033[1;34m',
        'cyan': '\033[0;36m',
        'CYAN': '\033[1;36m',
        'violet': '\033[0;35m',
        'VIOLET': '\033[1;35m',
    }
    #self.terminalBackToBlackCode = '\033[m'

    assert theColour in goodTerminalColours, "The colour must be one of %s" % goodTerminalColours
    sys.stdout.write(terminalColourCodeDict[theColour])


def unsetTerminalColour():
    sys.stdout.write('\033[m')

# Color and colour.
writeInColor = writeInColour
setTerminalColor = setTerminalColour
unsetTerminalColor = unsetTerminalColour


def _sumOfSquares(seq):
    """Pure Python.  Converts to floats, returns a float."""

    if not seq:
        return 0.0
    mySum = 0.0
    for it in seq:
        mySum += (float(it) * float(it))
    return mySum


def sortListOfObjectsOnAttribute(aListOfObjects, attributeString):
    """Returns a new sorted list."""

    # Repaired for py3. 

    # In the original (that worked with py2 all these years) used pairs, where,
    # 'paired' was a list of 2-element tuples, (thingToSortOn, theOriginalObject)

    # The problem in py3 is in tie scores in thingToSortOn, where sort() will
    # then attempt to compare objects, and that will generally fail, cuz objA <
    # objB (in Py3) will throw a
    # TypeError: '<' not supported between instances of 'Foo' and 'Foo'

    # So I simply add a third element, the indx in the middle, that allows
    # reliable sorting when there is a tie score at the first position.

    def tripling(anObject, indx, a=attributeString):
        return (getattr(anObject, a), indx, anObject)
    triplets = []
    for indx, ob in enumerate(aListOfObjects):
        triplets.append(tripling(ob, indx, a=attributeString))
    #print(triplets)
    triplets.sort()
    def stripit(aTriplet):
        return aTriplet[2]
    return [stripit(tr) for tr in triplets]


def sortListOfObjectsOn2Attributes(aListOfObjects, attributeString1, attributeString2):
    """Returns a new sorted list."""

    # Fixed for py3 as in sortListOfObjectsOnAttribute, adding a 4th element
    def quadding(anObject, indx, a=attributeString1, b=attributeString2):
        return (getattr(anObject, a), getattr(anObject, b), indx, anObject)
    quads = []
    for indx, ob in  enumerate(aListOfObjects):
        quads.append(quadding(ob, indx, a=attributeString1, b=attributeString2))
    #print(quads)
    quads.sort()
    def stripit(quad):
        return quad[3]
    return [stripit(qu) for qu in quads]

def sortListOfListsOnListElementNumber(aListOfLists, elementNumber):
    """Returns a new sorted list."""

    def pairing(aList, e=elementNumber):
        return (aList[e], aList)
    paired = [pairing(aList) for aList in aListOfLists]
    # 'paired' is a list of 2-element tuples, (thingToSortOn, theOriginalElement)
    # print "paired = ", paired
    paired.sort()

    def stripit(pair):
        return pair[1]
    return [stripit(pr) for pr in paired]

def readAndPop(stuff):
    """Read in simple stuff, pop the single object from var lists, and return it.

    The stuff to be read in must be convertible to a single object,
    one of Alignment, SequenceList, or Tree.  When that is read, the
    stuff as usual goes into one of var.alignments, var.sequenceLists,
    or var.trees.  The single object is popped from where it ends up,
    and returned.

    """
    gm = ['func.readAndPop()']
    #assert os.path.isfile(fName)
    onSL = len(var.sequenceLists)
    onAlig = len(var.alignments)
    onTrees = len(var.trees)
    read(stuff)
    nnSL = len(var.sequenceLists) - onSL
    nnAlig = len(var.alignments) - onAlig
    nnTrees = len(var.trees) - onTrees
    mySum = nnSL + nnAlig + nnTrees
    if mySum != 1:
        if mySum < 1:
            gm.append("no appropriate objects were made.")
        else:
            gm.append("Got %i objects.  Only 1 allowed." % mySum)
        raise P4Error(gm)
    if nnSL:
        return var.sequenceLists.pop()
    elif nnAlig:
        return var.alignments.pop()
    elif nnTrees:
        return var.trees.pop()

##########################################################################
##########################################################################
##########################################################################


def charsets(names, lens, fName=None):
    """Write a nice NEXUS sets block, given partition names and lengths.

    For example::

        geneNames = 'coi cytb nad5'.split()
        lens = []
        for geneName in geneNames:
            read('%s.nex' % geneName)
            lens.append(var.alignments[-1].nChar)
        func.charsets(geneNames, lens)

    which writes::

        #nexus

        begin sets;
          charset coi = 1 - 131;  [nChar = 131]
          charset cytb = 132 - 352;  [nChar = 221]
          charset nad5 = 353 - 521;  [nChar = 169]
          charpartition p1 =  coi:coi, cytb:cytb, nad5:nad5 ;
          [partition p1 = 3:coi, cytb, nad5;]
        end;

    """

    gm = ['func.charsets()']
    if len(names) != len(lens):
        gm.append("len of names (%i) is not the same as the len of lengths (%i)" % (
            len(names), len(lens)))
        print("names: ", names)
        print("lens: ", lens)
        raise P4Error(gm)

    start = 1
    if fName:
        f = open(fName, 'w')
    else:
        f = sys.stdout

    f.write("#nexus\n\n")
    f.write("begin sets;\n")
    for cNum in range(len(names)):
        f.write("  charset %s = %i - %i;" %
                (names[cNum], start, start - 1 + lens[cNum]))
        f.write(" [nChar = %i]\n" % lens[cNum])
        start += lens[cNum]
    f.write("  charpartition p1 = ")
    pp = ["%s:%s" % (nm, nm) for nm in names]
    f.write(', '.join(pp))
    f.write(';\n')

    f.write("  [partition p1 = %i:" % len(names))
    f.write(', '.join(names))
    f.write(';]\n')

    f.write("end;\n")
    if fName:
        f.close()


def reseedCRandomizer(newSeed):
    """Set a new seed for the c-language random() function.

    For those things in the C-language that use random(), this
    re-seeds the randomizer.  Reseed to different integers to make
    duplicate runs of something come out different.  This is not for
    the GSL random stuff.  Confusing to have 2 systems, innit?  And
    then there is the 3rd system that Python uses in the random
    module.  Sorry!

    If you wanted for example to repeat an MCMC, there are four random
    number generators to contend with.  You could do something like
    this (where I am using a seed of zero for all -- you could use
    your own seed, and they do not need to be the same)::

        # for C random() function, used in Brent-Powell optimization
        func.reseedCRandomizer(0)
        
        # for gsl, used in simulations (amongst other places)
        var.gsl_rng = pf.gsl_rng_get()
        pf.gsl_rng_set(var.gsl_rng, 0)
        
        # for the Python random library, used a lot in python code
        random.seed(0)
        
        # for Numpy and Scipy.  Used in func.dirichlet2()
        numpy.random.seed(0)

    """

    pf.reseedCRandomizer(newSeed)


def gsl_meanVariance(seq, mean=None, variance=None):
    """Use gsl to compute both the mean and variance.

    Arg seq can be a list or a numpy array.

    Returns a 2-tuple of single-item NumPy arrays.  To save a little
    time, you can pass the mean and variance to this function, in
    the form of zero-dimensional, single item NumPy arrays
    (eg mean = numpy.array(0.0))

    The numpy built-in variance function does not use n-1 weighting.
    This one (from gsl) does use n-1 weighting.

    """

    if isinstance(seq, numpy.ndarray):
        mySeq = seq
    else:
        mySeq = numpy.array(seq, numpy.float)

    if mean is None:
        mean = numpy.array([0.0])
    if variance is None:
        variance = numpy.array([0.0])
    pf.gsl_meanVariance(mySeq, len(seq), mean, variance)
    if 0:
        print("slow p4: mean=%f, variance=%f" % (func.mean(list(seq)), func.variance(list(seq))))
        print("gsl: mean=%f, variance=%f" % (mean, variance))
        # different than gsl-- no n-1 weighting.
        print("numpy: mean=%f, variance=%f" % (mySeq.mean(), mySeq.var()))
    return (mean, variance)


def chiSquaredProb(xSquared, dof):
    """Returns the probability of observing X^2."""
    return pf.chiSquaredProb(xSquared, dof)


def gsl_ran_gamma(a, b, seed=None):
    """See also random.gammavariate()"""

    complaintHead = '\nfunc.gsl_ran_gamma()'
    gm = complaintHead
    try:
        a = float(a)
        b = float(b)
    except:
        gm.append("Both a and b should be floats.")
        raise P4Error(gm)

    isNewGSL_RNG = 0
    if not var.gsl_rng:
        var.gsl_rng = pf.gsl_rng_get()
        isNewGSL_RNG = 1
        # print "got var.gsl_rng = %i" % var.gsl_rng
        # sys.exit()

        # Set the GSL random number generator seed, only if it is a new GSL_RNG
        if isNewGSL_RNG:
            if seed != None:
                try:
                    newSeed = int(seed)
                    pf.gsl_rng_set(var.gsl_rng, newSeed)
                except ValueError:
                    print(complaintHead)
                    print("    The seed should be convertible to an integer")
                    print("    Using the process id instead.")
                    pf.gsl_rng_set(var.gsl_rng,  os.getpid())
            else:
                pf.gsl_rng_set(var.gsl_rng,  os.getpid())

    return pf.gsl_ran_gamma(var.gsl_rng, a, b)


def dirichlet2(inSeq, outSeq, alpha, theMin):
    """Modify inSeq with a draw from a Dirichlet distribution with a single alpha value.

    Args *inSeq* and *outSeq* are both numpy arrays.  Arg *inSeq*
    is not modified itself; the modification is placed in *outSeq*
    (and nothing is returned).  The result is normalized to 1.0

    This is about 20% slower than :func:`func.dirichlet1` -- not quite sure why.
    """

    gm = ['func.dirichlet2()']
    assert isinstance(inSeq, numpy.ndarray)
    assert isinstance(outSeq, numpy.ndarray)
    kk = len(inSeq)
    theMax = 1.0 - ((kk - 1) * theMin)
    safety = 0
    safetyLimit = 300
    while 1:
        for i in range(kk):
            #theSeq[i] = pf.gsl_ran_gamma(var.gsl_rng, theSeq[i] * alpha, 1.0)
            # outSeq[i] = random.gammavariate(inSeq[i] * alpha, 1.0)   --- this
            # is very slow!
            outSeq[i] = numpy.random.gamma(inSeq[i] * alpha)
            #outSeq[i] = inSeq[i] * numpy.random.gamma(alpha)
        outSeq /= outSeq.sum()
        thisMin = numpy.amin(outSeq)
        thisMax = numpy.amax(outSeq)
        isOk = True
        if (thisMin < theMin) or (thisMax > theMax):
            isOk = False
        if isOk:
            return
        safety += 1
        if safety > safetyLimit:
            gm.append(
                "Tried more than %i times to get good dirichlet values, and failed.  Giving up." % safetyLimit)
            raise P4Error(gm)


def gsl_ran_dirichlet(alpha, theta):
    """Make a random draw from a dirichlet distribution.

    Args:

        alpha and theta (numpy arrays): both the same length (more
            than 1).  The length is the dimension of the dirichlet.
            The contents of *theta* are over-written (without being
            used).  The draw ends up in *theta*.  It is normalized so
            that it sums to 1.0.


    """

    complaintHead = '\nfunc.gsl_ran_dirichlet()'
    gm = complaintHead
    try:
        assert isinstance(alpha, numpy.ndarray)
        assert isinstance(theta, numpy.ndarray)
    except AssertionError:
        gm.append(" alpha, theta, should be numpy arrays.")
        raise P4Error(gm)

    assert len(alpha) > 1
    assert len(theta) == len(alpha)

    if not var.gsl_rng:
        var.gsl_rng = pf.gsl_rng_get()
        pf.gsl_rng_set(var.gsl_rng, int(time.time()))

    pf.gsl_ran_dirichlet(var.gsl_rng, len(theta), alpha, theta)


def studentsTTest1(seq, mu=0.0, verbose=True):
    """Test whether a sample differs from mu.

    From wikipedia.

    Arg 'seq' is a list of numbers.  Internally it is converted to a
    numpy array of floats, so the input seq need not be floats, and
    need not be a numpy array, although it does not hurt to be either.

    Arg 'mu' is by default zero.

    Returns the p-value.
    """

    sq = numpy.array(seq, dtype=numpy.float)
    m, v = gsl_meanVariance(sq)
    s = numpy.sqrt(v)
    n = len(sq)
    sqN = numpy.sqrt(n)
    t = (m - mu) / (s / sqN)
    dof = n - 1
    p = pf.studentsTProb(t, dof)

    if verbose:
        print("mean =", m)
        print("std dev =", s)
        print("t statistic =", t)
        print("prob =", p)
    return p


def effectiveSampleSize(data, mean):
    """As done in Tracer v1.4, by Drummond and Rambaut.  Thanks guys!

    But see :func:`p4.func.summarizeMcmcPrams`, which gives ESSs.

    """

    nSamples = len(data)
    maxLag = int(nSamples / 3)
    if maxLag > 1000:
        maxLag = 1000

    gammaStatAtPreviousLag = numpy.array([0.0])
    assert isinstance(data, numpy.ndarray), "Arg 'data' should be a numpy.array.  Its %s" % type(data)
    assert isinstance(mean, numpy.ndarray), "Arg 'mean' should be a numpy.array.  Its %s" % type(mean)
    gammaStat = numpy.array([0.0])
    varStat = numpy.array([0.0])
    gammaStatAtLagZero = numpy.array([0.0])
    if 0:
        lag = 0
        while lag < maxLag:
            gammaStat[0] = 0.0
            for j in range(nSamples - lag):
                gammaStat[0] += (data[j] - mean) * (data[j + lag] - mean)

            # if lag == 0:
            #    print "mean is %f" % mean
            #    print "lag is 0, gammaStat = %f" % gammaStat[0]

            gammaStat[0] /= (nSamples - lag)

            if lag == 0:
                varStat[0] = gammaStat
                gammaStatAtLagZero[0] = gammaStat
                # print "got gammaStatAtLagZero = %f" % gammaStatAtLagZero[0]
            elif (lag % 2) == 0:
                if gammaStatAtPreviousLag + gammaStat > 0:
                    varStat[0] += 2.0 * (gammaStatAtPreviousLag + gammaStat)
                else:
                    break
            lag += 1
            gammaStatAtPreviousLag[0] = gammaStat
            #gammaStat[0] = 0.0

        # print gammaStatAtLagZero, gammaStat, varStat, lag
        # print "maxLag is %i" % maxLag
        # stdErrorOfMean = numpy.sqrt(varStat / nSamples)  ??!?
        #ACT = stepSize * varStat / gammaStatAtLagZero
        #ESS = (stepSize * nSamples) / ACT;
        # stepSize is not needed
        ESS1 = nSamples * (gammaStatAtLagZero / varStat)
        # print "got ESS1 %f" % ESS1

    if 1:
        pf.effectiveSampleSize(data, mean, nSamples, maxLag, gammaStatAtPreviousLag,
                               gammaStat, varStat, gammaStatAtLagZero)

        # print(f"nSamples {nSamples} type {type(nSamples)}")
        # print(f"gammaStatAtLagZero {gammaStatAtLagZero} type {type(gammaStatAtLagZero)}")
        # print(f"varStat {varStat} type {type(varStat)}")
        
        if varStat == 0:
            return 0.0  # Is this best for this case?
        else:
            ESS2 = nSamples * (gammaStatAtLagZero / varStat)
        #fabsDiff = numpy.fabs(ESS1 - ESS2)
        # print "ESS1 is %f, ESS2 is %f, diff is %g" % (ESS1, ESS2, fabsDiff)
    return ESS2[0]


def summarizeMcmcPrams(skip=0, run=-1, theDir='.', makeDict=False):
    """Find the mean, variance, and ess of mcmc parameters.

    Ess is effective sample size, as in Tracer by Drummond and Rambaut.

    The numbers are found in mcmc_prams_N, N=0, 1, etc.  If arg 'run' is
    set to -1, the default, then all runs are done.  Alternatively you
    can set the run to a specific run number, and that is the only one
    that is done.

    The 'profile', with the names of the parameters, and the number of each, is
    found in mcmc_pramsProfile_N.py (N=0,1,2 ...).  It is not essential, but it
    gives names to the parameters.

    """

    gm = ["func.summarizeMcmcPrams()"]
    nPrams = None
    pramsProfile = None

    numsList = None
    if run == -1:
        runNum = 0
        # read in the prams profile only once.  Assume it will apply to all runs.
        try:
            loc = {}
            theFName = "mcmc_pramsProfile_%i.py" % runNum
            exec(open(os.path.join(theDir, theFName)).read(), {}, loc)
            # loc =locals()  no workee.
            # print "loc = %s" % loc
            nPrams = loc['nPrams']
            pramsProfile = loc['pramsProfile']
        except IOError:
            print("The file '%s' cannot be found." % theFName)
            if makeDict:
                print("Cannot make dictionary without %s" % theFName)
                return None

    else:
        runNum = run
    totalLinesRead = 0

    while 1:

        if run != -1:
            # possibly different prams profiles for each run
            try:
                loc = {}
                theFName = "mcmc_pramsProfile_%i.py" % runNum
                exec(open(os.path.join(theDir, theFName)).read(), {}, loc)
                # loc =locals()  no workee.
                # print "loc = %s" % loc
                nPrams = loc['nPrams']
                pramsProfile = loc['pramsProfile']
            except IOError:
                print("The file '%s' cannot be found." % theFName)
                if makeDict:
                    print("Cannot make dictionary without %s" % theFName)
                    return None

        try:
            theFName = os.path.join(theDir, "mcmc_prams_%i" % runNum)
            flob = open(theFName)
            if not makeDict:
                print("Reading prams from file %s" % theFName)
        except IOError:
            break

        theLines = flob.readlines()
        flob.close()
        runNum += 1
        skipsDone = 0
        linesRead = 0
        for aLine in theLines:
            ll = aLine.lstrip()
            if ll.startswith("#"):
                pass
            elif not ll:
                pass
            elif ll[0] not in string.digits:
                pass
            else:
                if skipsDone < skip:
                    skipsDone += 1
                else:
                    splitLine = aLine.split()
                    # If it does not exist yet, then make it now.
                    if not numsList:
                        thisNPrams = len(splitLine) - 1
                        if nPrams:
                            if not thisNPrams == nPrams:
                                gm.append(
                                    "thisNPrams = %i, nPrams = %s" % (thisNPrams, nPrams))
                                raise P4Error(gm)
                        else:
                            nPrams = thisNPrams
                        numsList = []
                        for pramNum in range(nPrams):
                            numsList.append([])
                    for pramNum in range(nPrams):
                        try:
                            theOne = splitLine[pramNum + 1]
                        except IndexError:
                            gm.append("Line '%s'.  " % aLine.rstrip())
                            gm.append(
                                "Can't get parameter number %i  " % pramNum)
                            raise P4Error(gm)
                        try:
                            aFloat = float(theOne)
                            numsList[pramNum].append(aFloat)
                        except (ValueError, TypeError):
                            gm.append("Can't make sense of '%s'" % theOne)
                            raise P4Error(gm)
                    linesRead += 1
        if not makeDict:
            print("  skipped %i lines" % skipsDone)
            print("  read %i lines" % linesRead)
        totalLinesRead += linesRead
        if run != -1:
            break
    if not makeDict:
        print("Read %i pram lines in total." % totalLinesRead)

    # print numsList
    spacer1 = ' ' * 20
    pramsdict = {}
    if pramsProfile:
        if not makeDict:
            print("%s   %26s         mean      variance       ess  " % (spacer1, ' '))
            print("%s   %26s       --------    --------    --------" % (spacer1, ' '))
        pramCounter = 0
        for partNum in range(len(pramsProfile)):
            if len(pramsProfile) > 1:
                if not makeDict:
                    print("Data partition %i" % partNum)
            if len(pramsProfile[partNum]):
                pramsdict["part%i" % partNum] = []
                # print pramsProfile[partNum]
                for pramNum in range(len(pramsProfile[partNum])):
                    pString = pramsProfile[partNum][pramNum][0]
                    pramCounts = pramsProfile[partNum][pramNum][1]
                    for p in range(pramCounts):
                        if not makeDict:
                            print("%s%3i %22s[%2i]   " % (spacer1, pramCounter, pString, p), end=' ')
                        d = numpy.array(numsList[pramCounter], numpy.float)
                        m, v = gsl_meanVariance(d)
                        ess = effectiveSampleSize(d, m)
                        stats = []
                        stats.append("%s[%i]" % (pString, p))
                        if m == 0.0:
                            if not makeDict:
                                print("  0.0      ", end=' ')
                            else:
                                stats.append("0.0")
                        elif m < 0.00001:
                            if not makeDict:
                                print("%10.3g " % m, end=' ')
                            else:
                                stats.append("%.3g" % m)
                        elif m < 1.0:
                            if not makeDict:
                                print("%10.6f " % m, end=' ')
                            else:
                                stats.append("%.6f" % m)
                        else:
                            if not makeDict:
                                print("%10.4f " % m, end=' ')
                            else:
                                stats.append("%.4f" % m)
                        if v == 0.0:
                            if not makeDict:
                                print("  0.0      ", end=' ')
                            else:
                                stats.append("0.0")
                        elif v < 0.000001:
                            if not makeDict:
                                print("%10.3g " % v, end=' ')
                            else:
                                stats.append("%.3g" % v)
                        elif v < 1.0:
                            if not makeDict:
                                print("%10.6f " % v, end=' ')
                            else:
                                stats.append("%.6f " % v)
                        else:
                            if not makeDict:
                                print("%10.4f " % v, end=' ')
                            else:
                                stats.append("%.6f" % v)
                        if not makeDict:
                            print("%10.1f " % ess, end=' ')
                        else:
                            stats.append("%.1f" % ess)
                        if not makeDict:
                            print()
                        pramCounter += 1
                        pramsdict["part%i" % partNum].append(stats)

            else:
                if not makeDict:
                    print("        No parameters in this data partition.")
    else:  # no pramsProfile
        print("%9s  mean      variance       ess  " % ' ')
        print("%9s--------    --------    --------" % ' ')
        for pramNum in range(nPrams):
            print("  %2i  " % pramNum, end=' ')
            d = numpy.array(numsList[pramNum], numpy.float)
            m, v = gsl_meanVariance(d)
            ess = effectiveSampleSize(d, m)
            if m == 0.0:
                print("  0.0      ", end=' ')
            elif m < 0.00001:
                print("%10.3g " % m, end=' ')
            elif m < 1.0:
                print("%10.6f " % m, end=' ')
            else:
                print("%10.4f " % m, end=' ')

            if v == 0.0:
                print("  0.0      ", end=' ')
            elif v < 0.000001:
                print("%10.3g " % v, end=' ')
            elif v < 1.0:
                print("%10.6f " % v, end=' ')
            else:
                print("%10.4f " % v, end=' ')

            print("%10.1f " % ess, end=' ')
            print()

    if makeDict:
        return pramsdict


def newtonRaftery94_eqn16(logLikes, delta=0.1, verbose=False):
    """Importance sampling, as in Newton and Raftery 1994, equation 16"""

    lla = numpy.array(logLikes)

    # Start the iteration with the harm mean, so calculate it.

    # This first way is the way from Susko I think it was, via Jessica
    # Leigh.  It is not very good.  But it probably does not matter,
    # as it is only a starting point for an iteration.
    if 1:
        theMax = -numpy.min(lla)
        diff = 700. - theMax
        shifted = (-lla) + diff
        expd = numpy.exp(shifted)
        theSum = numpy.sum(expd)
        theHarmMean = float(-(numpy.log(theSum) - diff))
        # print theHarmMean
        # print type(theHarmMean)

        return pf.newtonRaftery94_eqn16(lla, len(logLikes), theHarmMean, delta, int(verbose))

    #n = Numbers(logLikes)
    # theHarmMean = n.harmonicMeanOfLogs()  # better
    # return pf.newtonRaftery94_eqn16(lla, len(logLikes), theHarmMean, delta,
    # int(verbose))

def readJplace(fName, verbose=False):
    """Read ``*.jplace`` files for phylogenetic placement of short reads

    It returns a tuple --- (tree, heatSum).
    The tree is decorated with a 'heat' attribute on each node.br.
    The heatSum is the sum of all the heats.
    """
    
    import json

    assert fName.endswith("jplace"), "Is %s a jplace file? --- it does not end with '.jplace'" % fName

    # read with json
    flob = open(fName)
    jobj = json.load(flob)
    flob.close()

    # Read the tree
    inBraces = False
    tLst2 = []
    tStr = jobj['tree']
    #print(tStr)
    tLst = list(tStr)
    for c in tLst:
        if inBraces:
            if c == '}':
                inBraces = False
        else:
            if c == '{':
                inBraces = True
            else:
                tLst2.append(c)
    tStr = ''.join(tLst2)
    #print(tStr)
    t = readAndPop(tStr)

    if verbose:
        print("Got tree with %i nodes" % len(t.nodes))

    # JPlace files use postorder node numbering. P4 uses preorder.
    # Po is PostOrder.
    nodeForPoNumDict = {}
    for poNum, nNum in enumerate(t.postOrder):
        #print("post-order =", poNum, "node num = ",nNum)
        nodeForPoNumDict[poNum] = t.nodes[nNum]

    placements = jobj['placements']
    #n_p is the number of placement lines
    n_p = len(placements)
    if verbose:
        print("There are %i placement lines (n_p)" % n_p)

    pResults = [0.0] * len(t.nodes)

    for pLine in placements:
        #print(type(pLine))     # a dict
        nn = pLine['n']
        assert len(nn) == 1

        pp = pLine['p']
        pLine['S_lwr'] = sum([it[2] for it in pp])
        #print("Got sum of lwr's %f" % pLine['S_lwr'])
        for pNode in pp:
            nNum = pNode[0]
            contrib = pNode[2] / pLine['S_lwr']
            pResults[nNum] += contrib

    heatSum = 0.0
    for poNum, pR in enumerate(pResults):
        n = nodeForPoNumDict[poNum]
        if n == t.root:
            assert math.fabs(pR) < 1.e-14, "Something wrong. Root has placement result %g" % pR
        else:
            n.br.heat = pR
            if verbose:
                print("Got node.br 'heat'", pR)
            heatSum += pR
    if verbose:
        print("Sum of heats", heatSum)
        print("returning a tuple of the heat-decorated Tree from %s, and heatSum..." % fName)
    return (t, heatSum)


def matrixLowerTriangleToUpperTriangle(ltList, dim):
    """Rearrange matrix from a lower triangle to an upper triangle list

    Args:
        ltList (list): a lower triangle of a matrix in the form of a list

        dim (int): the dimension of the (square) matrix

    Returns:
        the same items rearranged as the upper triangle, in the form of a list.

    If your matrix is like this::

        - - - -
        A - - -
        B D - -
        C E F -

    where the dim is 4, then the lower triangle is the list [A,B,D,C,E,F].
    This function rearranges that so that it is the upper triangle::

        - A B C
        - - D E
        - - - F
        - - - -

    and returns the list [A,B,C,D,E,F].

    This is useful for user-specified empirical protein rate matrices
    where the rate matrix is given as a PAML-style lower triangle, but
    you need a p4-style upper triangle::

        pamlStyleRates = [190 rates]
        upTriangle = func.matrixLowerTriangleToUpperTriangle(pamlStyleRates, 20)
        # then when you specify your model ...
        t.newRMatrix(free=0, spec='specified', val=upTriangle)

    """
    
    assert isinstance(dim, int)
    expectedLen = int(((dim * dim) - dim)/2)
    assert len(ltList) == expectedLen, f"Lower triangle len is {len(ltList)}, expected {expectedLen}"

    bigM = []
    for rNum in range(dim):
        bigM.append(['0.0'] * dim)

    i = 0
    for rNum in range(1,dim):
        for cNum in range(rNum):
            #print(ltList[i])
            bigM[cNum][rNum] = ltList[i]
            i += 1
    upper = []
    for rNum in range(dim-1):
        for cNum in range(rNum+1, dim):
            #print(rNum, cNum, bigM[rNum][cNum])
            upper.append(bigM[rNum][cNum])
    return upper


def _compareSplitsBetweenTwoTreePartitions(tp1, tp2, minimumProportion, verbose=False):
    """Returns a tuple of asdoss, maximum of the differences and mean of the differences

    This calls the method TreePartitions.compareSplits(), and digests the
    results returned from that.

    Args:
        tp1, tp2 (TreePartition): TreePartition objects
        minimumProportion (float): passed to TreePartitions.compareSplits()

    Returns:
        (asdoss, maxOfDiffs, meanOfDiffs)

    """

    ret = tp1.compareSplits(tp2, minimumProportion=minimumProportion)

    #print(ret)  # a list of 3-item lists
    #  1. The split key
    #  2. The split string
    #  3. A list of the 2 supports

    if not ret:
        return None

    sumOfStdDevs = 0.0
    nSplits = len(ret)
    diffs = []
    for i in ret:
        # print "            %.3f  %.3f    " % (i[2][0], i[2][1]),
        stdDev = math.sqrt(p4.func.variance(i[2]))
        # print "%.5f" % stdDev
        sumOfStdDevs += stdDev
        diffs.append(math.fabs(i[2][0] - i[2][1]))
    asdoss = sumOfStdDevs / nSplits
    maxOfDiffs = max(diffs)
    meanOfDiffs = sum(diffs)/nSplits
    if verbose:
        print("     nSplits=%i, average of std devs of split supports %.4f " % (nSplits, asdoss))
        print("     max of differences %f, mean of differences %f" % (maxOfDiffs, meanOfDiffs))
    return (asdoss, maxOfDiffs, meanOfDiffs)


def compareSplitsBetweenTreePartitions(treePartitionsList, precision=3, linewidth=120):
    """Pairwise ASDOSS (ASDSF) and MaxDiff calculations 

    Output is verbose.  Shows 
    - average standard deviation of split frequencies (or supports), like MrBayes
    - maximum difference between split supports from each pair of checkpoints, like PhyloBayes

    Returns:
        None

    """

    nM = len(treePartitionsList)
    supportMins = numpy.zeros(nM, dtype=numpy.float64)
    supportMaxs = numpy.zeros(nM, dtype=numpy.float64)
    nItems = int(((nM * nM) - nM) / 2)
    asdosses = numpy.zeros((nM, nM), dtype=numpy.float64)
    vect = numpy.zeros(nItems, dtype=numpy.float64)
    mdvect = numpy.zeros(nItems, dtype=numpy.float64)
    maxDiffs = numpy.zeros((nM, nM), dtype=numpy.float64)

    for mNum1 in range(nM):
        tp1 = treePartitionsList[mNum1]
        tp1.finishSplits()
        ret = tp1.getProportionRange(verbose=False)
        supportMins[mNum1] = ret[0]
        supportMaxs[mNum1] = ret[1]
    sRanges = supportMaxs - supportMins


    minimumProportion=0.1
    vCounter = 0
    for mNum1 in range(1, nM):
        tp1 = treePartitionsList[mNum1]
        for mNum2 in range(mNum1):
            tp2 = treePartitionsList[mNum2]
            thisAsdoss, thisMaxDiff, thisMeanDiff = p4.func._compareSplitsBetweenTwoTreePartitions(
                tp1, tp2, minimumProportion, verbose=False)
            #if thisAsdoss == None and verbose:
            #    print("No splits > %s" % minimumProportion)

            if thisAsdoss == None:
                thisAsdoss = 0.0
            asdosses[mNum1][mNum2] = thisAsdoss
            asdosses[mNum2][mNum1] = thisAsdoss
            vect[vCounter] = thisAsdoss
            maxDiffs[mNum1][mNum2] = thisMaxDiff
            maxDiffs[mNum2][mNum1] = thisMaxDiff
            mdvect[vCounter] = thisMaxDiff
            vCounter += 1

            if 0:
                print(" %10i " % mNum1, end=' ')
                print(" %10i " % mNum2, end=' ')
                print("%.3f" % thisAsdoss)

    # Save current numpy printoptions, and restore, below.
    curr = numpy.get_printoptions()
    numpy.set_printoptions(precision=precision, linewidth=linewidth)

    print("Support values within tree partitions. Maxs, mins, and ranges")
    print(supportMaxs)
    print(supportMins)
    print(sRanges)
    print()

    print("Pairwise average standard deviation of split frequency values ---")
    print(asdosses)
    print()
    print("For the %i values in one triangle," % nItems)
    print("max =  %.3f" % vect.max())
    print("min =  %.3f" % vect.min())
    print("mean = %.3f" % vect.mean())
    #print("var =  %.2g", vect.var())

    print()
    print("Pairwise maximum differences in split supports between the runs ---")
    print(maxDiffs)
    print()
    print("For the %i values in one triangle," % nItems)
    print("max =  %.3f" % mdvect.max())
    print("min =  %.3f" % mdvect.min())
    print("mean = %.3f" % mdvect.mean())
    #print("var =  ", mdvect.var())


    # Reset printoptions back to what it was
    numpy.set_printoptions(precision=curr['precision'], linewidth=curr['linewidth'])

def getProteinEmpiricalModelComp(spec):
    """Return comps as a list of floats, normalized to sum to 1.0"""

    assert spec in var.rMatrixProteinSpecs, "spec should be one of var.rMatrixProteinSpecs"

    # Several of these empirical protein comps are from the dat
    # files in PAML.  Thanks, Ziheng!

    if spec == 'cpREV':
        val = [0.0755, 0.0621, 0.0410, 0.0371, 0.0091,
               0.0382, 0.0495, 0.0838, 0.0246, 0.0806,
               0.1011, 0.0504, 0.0220, 0.0506, 0.0431,
               0.0622, 0.0543, 0.0181, 0.0307, 0.0660]

    elif spec == 'd78':
        # These first values have a couple more decimal places.  I
        # think I got these from some obscure code in a back corner of
        # the NCBI ftp site.  I believe it is obtainable by raising
        # the d78 matrix to a high power.  It is a more precise comp,
        # but is not used here because it is not standard.
        # val = [0.08713, 0.04090, 0.04043, 0.04687, 0.03347,
        #        0.03826, 0.04953, 0.08861, 0.03362, 0.03689,
        #        0.08536, 0.08048, 0.01475, 0.03977, 0.05068,
        #        0.06958, 0.05854, 0.01049, 0.02992, 0.06472]

        # This next set of values is from her paper, and is the set
        # that everybody uses.
        # val = [0.087, 0.041, 0.040, 0.047, 0.033,
        #        0.038, 0.05, 0.089, 0.034, 0.037,
        #        0.085, 0.08, 0.015, 0.04, 0.051,
        #        0.07, 0.058, 0.01, 0.03, 0.065]

        # These values are from Goldman's recommendations (Kosiol &
        # Goldman)
        val = [0.087127, 0.040904, 0.040432, 0.046872, 0.033474,
               0.038255, 0.049530, 0.088612, 0.033619, 0.036886,
               0.085357, 0.080481, 0.014753, 0.039772, 0.050680,
               0.069577, 0.058542, 0.010494, 0.029916, 0.064718]

    elif spec == 'jtt':
        # val = [0.077,0.051, 0.043, 0.052, 0.02,
        #        0.041, 0.062, 0.074, 0.023, 0.052,
        #        0.091, 0.059, 0.024, 0.04, 0.051,
        #        0.069, 0.059, 0.014, 0.032, 0.066]

        # again, a Goldman recommendation
        val = [0.076862, 0.051057, 0.042546, 0.051269, 0.020279,
               0.041061, 0.061820, 0.074714, 0.022983, 0.052569,
               0.091111, 0.059498, 0.023414, 0.040530, 0.050532,
               0.068225, 0.058518, 0.014336, 0.032303, 0.066374]

    elif spec == 'mtREV24':
        val = [0.072, 0.019, 0.039, 0.019, 0.006,
               0.025, 0.024, 0.056, 0.028, 0.088,
               0.168, 0.023, 0.054, 0.061, 0.054,
               0.072, 0.086, 0.029, 0.033, 0.043]

    elif spec == 'mtmam':
        val = [0.0692, 0.0184, 0.0400, 0.0186, 0.0065,
               0.0238, 0.0236, 0.0557, 0.0277, 0.0905,
               0.1675, 0.0221, 0.0561, 0.0611, 0.0536,
               0.0725, 0.0870, 0.0293, 0.0340, 0.0428]

    elif spec == 'wag':
        val = [0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078,
               0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466,
               0.086209,  0.0620286, 0.0195027, 0.0384319, 0.0457631,
               0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956]

    elif spec == 'rtRev':
        val = [0.0646, 0.0453, 0.0376, 0.0422, 0.0114, 0.0606,
               0.0607, 0.0639, 0.0273, 0.0679, 0.1018, 0.0751,
               0.0150, 0.0287, 0.0681, 0.0488, 0.0622, 0.0251,
               0.0318, 0.0619]

    elif spec == 'tmjtt94':
        val = [0.105068479, 0.015695291, 0.018494452, 0.008897331,
               0.021893432, 0.014095771, 0.009697091, 0.075777267,
               0.016794962, 0.118764371, 0.163450965, 0.011196641,
               0.033290013, 0.077676697, 0.025992202, 0.056782965,
               0.052284315, 0.022293312, 0.032390283, 0.119464161]

    elif spec == 'tmlg99':
        val = [0.100632, 0.014017, 0.014706, 0.010371, 0.030668,
               0.015152, 0.011343, 0.069235, 0.017501, 0.107722,
               0.155161, 0.009723, 0.038730, 0.086453, 0.031761,
               0.064333, 0.044847, 0.028277, 0.036988, 0.112380]

    elif spec == 'lg':
        val = [0.079066, 0.055941, 0.041977, 0.053052, 0.012937,
               0.040767, 0.071586, 0.057337, 0.022355, 0.062157,
               0.099081, 0.064600, 0.022951, 0.042302, 0.044040,
               0.061197, 0.053287, 0.012066, 0.034155, 0.069147]

    elif spec == 'blosum62':
        val = [0.074, 0.052, 0.045, 0.054, 0.025,
               0.034, 0.054, 0.074, 0.026, 0.068,
               0.099, 0.058, 0.025, 0.047, 0.039,
               0.057, 0.051, 0.013, 0.032, 0.073]

    elif spec == 'hivb':
        val = [0.060490222, 0.066039665, 0.044127815, 0.042109048, 0.020075899,
               0.053606488, 0.071567447, 0.072308239, 0.022293943, 0.069730629,
               0.098851122, 0.056968211, 0.019768318, 0.028809447, 0.046025282,
               0.05060433, 0.053636813, 0.033011601, 0.028350243, 0.061625237]

    elif spec == 'mtart':
        val = [0.054116, 0.018227, 0.039903, 0.020160, 0.009709,
               0.018781, 0.024289, 0.068183, 0.024518, 0.092638,
               0.148658, 0.021718, 0.061453, 0.088668, 0.041826,
               0.091030, 0.049194, 0.029786, 0.039443, 0.057700]

    elif spec == 'mtzoa':
        val = [0.068880,    0.021037,    0.030390,    0.020696,    0.009966,
               0.018623,    0.024989,    0.071968,    0.026814,    0.085072,
               0.156717,    0.019276,    0.050652,    0.081712,    0.044803,
               0.080535,    0.056386,    0.027998,    0.037404,    0.066083]

    elif spec == 'gcpREV':
        val = [0.079510, 0.056001, 0.040459, 0.033220, 0.009051,
               0.037505, 0.049675, 0.080233, 0.021880, 0.080496,
               0.107512, 0.049324, 0.020776, 0.047731, 0.039916,
               0.073820, 0.053615, 0.016705, 0.030790, 0.071781]

    elif spec == 'stmtREV':
        val = [0.046181, 0.053408, 0.036197, 0.023332, 0.023417,
               0.039040, 0.034128, 0.038916, 0.016464, 0.089153,
               0.161731, 0.055134, 0.023326, 0.091125, 0.034471,
               0.077108, 0.041860, 0.020078, 0.030542, 0.064385]

    elif spec == 'vt':  #from iqtree
        val = [0.0770764620135024, 0.0500819370772208, 0.0462377395993731, 0.0537929860758246, 
               0.0144533387583345, 0.0408923608974345, 0.0633579339160905, 0.0655672355884439, 
               0.0218802687005936, 0.0591969699027449, 0.0976461276528445, 0.0592079410822730, 
               0.0220695876653368, 0.0413508521834260, 0.0476871596856874, 0.0707295165111524, 
               0.0567759161524817, 0.0127019797647213, 0.0323746050281867, 0.0669190817443274]

    elif spec == 'pmb':  # from iqtree
        val = [0.076, 0.054, 0.038, 0.045, 0.028, 
               0.034, 0.053, 0.078, 0.030, 0.060, 
               0.096, 0.052, 0.022, 0.045, 0.042, 
               0.068, 0.056, 0.016, 0.036, 0.071]

    theSum = sum(val)
    for i in range(len(val)):
        val[i] /= theSum
    return val

def getProteinEmpiricalModelRMatrix(spec, upperTriangle=True):
    """Return RMatrix as a numpy array, or upper triangle as a list"""

    dim = 20
    assert spec in var.rMatrixProteinSpecs, "spec should be one of var.rMatrixProteinSpecs"
    specNum = var.rMatrixProteinSpecNumberForNameDict[spec]
    a = numpy.zeros((dim, dim), numpy.float)
    pf.getBigR(specNum, a)

    if upperTriangle:
        uT = []
        for rNum in range(dim - 1):
            for cNum in range(rNum + 1, dim):
                uT.append(a[rNum][cNum])
        return uT
    else:
        return a


