"""Various functions."""

import os
import sys
import re
import string
import math
import cStringIO
import random
import glob
import time
import types
from Var import var
from SequenceList import Sequence,SequenceList
from Alignment import Alignment
from Nexus import Nexus
from Tree import Tree
from Node import Node
from Glitch import Glitch
from Constraints import Constraints
if var.usePfAndNumpy:
    try:
        import pf
        import numpy
    except ImportError:
        raise Glitch, "Can't import at least one of pf or numpy."
    from Numbers import Numbers


def nexusCheckName(theName):
    """Check to see if theName conforms to Nexus standards."""
    if type(theName) != type('string'):
        return 0
    if var.nexus_allowAllDigitNames:
        return 1
    if len(theName) == 1 and theName[0] not in string.letters:
        return 0
    try:
        int(theName)
        return 0 # we don't allow all-digit names
    except ValueError:
        return 1

def nexusUnquoteAndDeUnderscoreName(theName):
    """Deal with underscores and quotes.  Returns theName

    If theName is not quoted, convert any underscores to spaces.  If
    theName is quoted, remove the outside quotes, and convert any
    cases of 2 single quotes in a row to 1 single quote.
    """

    if theName[0] == '\'':
        assert theName[-1] == "'", "func.nexusUnquoteAndDeUnderscoreName().  First char is a single quote, but last char is not."
        theName = theName[1:-1]
        if theName.count('\'\''):
            return theName.replace('\'\'', '\'')
        else:
            return theName
    if string.count(theName, '_'):
        return theName.replace('_', ' ')
    else:
        return theName

def nexusUnquoteName(theName):
    """Deal with quotes.  Returns theName

    If theName is not quoted, just return it.  If
    theName is quoted, remove the outside quotes, and convert any
    cases of 2 single quotes in a row to 1 single quote."""

    if theName[0] == '\'':
        if theName[-1] != "'":
            gm = ['func.nexusUnquotName()']
            gm.append('the name is %s' % theName)
            gm.append("First char is a single quote, but last char is not.")
            raise Glitch, gm
        theName = theName[1:-1]
        if theName.count('\'\''):
            return theName.replace('\'\'', '\'')
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
    already doubled single quotes, don't double them.  """

    if theName == None:
        return theName
    if theName.startswith('\''):
        return theName
    quotesAreNeeded = 0
    #for c in var.nexus_punctuation:
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
        if theName.count('\''):
            # If we have doubled quotes, don't re-double them
            if theName.count('\'\''):
                pass
            else:
                theName = theName.replace('\'', '\'\'')
        newName = '\'' + theName + '\''
        if verbose:
            print "Warning. Nexus quoting |%s| to |%s|" % (oldName, newName)
        return newName
    else:
        return theName
    
        

def isDnaRnaOrProtein(aString):
    """Attempts to determinine the data type by the composition.

    Returns 1 for DNA, 2 for RNA, and 0 for protein.  Or so it thinks.

    It only works for lowercase symbol letters.
    """
    nGaps = string.count(aString, '-')
    nQs = string.count(aString, '?')
    strLenNoGaps = len(aString) - (nGaps + nQs)
    acgn = (string.count(aString, 'a') +
           string.count(aString, 'c') +
           string.count(aString, 'g') +
           string.count(aString, 'n'))
    t = string.count(aString, 't')
    u = string.count(aString, 'u')
    if 0:
        print "stringLength(no gaps) = %s" % strLenNoGaps
        print "acgn = %f" % acgn
        print "t = %f" % t
        print "u = %f" % u

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
    tally = (string.count(aString, 'r') +
             string.count(aString, 'y') +
             string.count(aString, 'm') +
             string.count(aString, 'k') +
             string.count(aString, 's') +
             string.count(aString, 'w') +
             string.count(aString, 'h') +
             string.count(aString, 'b') +
             string.count(aString, 'v') +
             string.count(aString, 'd') )
    threshold2 = 0.99 * strLenNoGaps
    #print "  string length is %i, strLenNoGaps is %i" % (len(aString), strLenNoGaps)
    #print "  threshold2 is %3.1f, acgn + t + ambigs = %i" % (threshold, (acgn + t + tally))
    #print "  threshold1 is %3.1f, acgn + t = %i" % (threshold, (acgn + t))
    if (acgn + t + tally >= threshold2) and (acgn + t >= threshold1):
        return 1
    elif (acgn + u + tally >= threshold2) and (acgn + u >= threshold1):
        return 2
    else:
        return 0


def stringZapWhitespaceAndDigits(inString):
    out = list(inString)
    theRange = range(len(out))
    theRange.reverse()
    for i in theRange:
        if out[i] in string.whitespace or out[i] in string.digits:
            del out[i]
    return string.join(out, '')


def dump():
    """A top-level dump of p4 trees, files, sequenceLists, and alignments."""

    print "\np4 dump"
    if len(var.fileNames) == 0:
        print "    Haven't read any files"
    elif len(var.fileNames) == 1:
        print "    Read one file: %s" % var.fileNames[0]
    else:
        print "    Read files:"
        for i in var.fileNames:
            print "                %s" % i

    if len(var.alignments) == 0:
        #print "    There are no alignments."
        pass
    elif len(var.alignments) == 1:
        print "    There is 1 alignment."
        #if recursive:
        #    var.alignments[0].dump()
    else:
        print "    There are %i alignments." % len(var.alignments)
        #if recursive:
        #    for i in var.alignments:
        #        i.dump()

    if len(var.sequenceLists) == 0:
        #print "    There are no sequenceLists."
        pass
    elif len(var.sequenceLists) == 1:
        print "    There is 1 sequenceList."
        #if recursive:
        #    var.sequenceLists[0].dump()
    else:
        print "    There are %i sequenceLists." % len(var.sequenceLists)
        #if recursive:
        #    for i in var.sequenceLists:
        #        i.dump()

    if len(var.trees) == 0:
        # print "    There are no trees."
        pass
    elif len(var.trees) == 1:
        print "    There is 1 tree."
        #if recursive:
        #    var.trees[0].dump()
    else:
        print "    There are %i trees." % len(var.trees)
        #if recursive:
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

    # Is this still a Bug?: single node phylip or raw newick trees must be specified as
    # ()A; they cannot be specified simply as A.
    
    gm = ['p4.read()']
    if type(stuff) != type('string'):
        gm.append("I was expecting a string argument.")
        raise Glitch, gm
    #nAlignments = len(var.alignments)
    if os.path.exists(stuff):
        #var.nexus_doFastNextTok = False
        readFile(stuff)
    else:
        # Is it a glob?
        myFlist = glob.glob(stuff)
        #print "read(). stuff=%s,  glob result: %s" % (stuff, myFlist)
        if myFlist: # It appears to be a glob
            for fName in myFlist:
                readFile(fName)
        else: # Nope, not a glob.  Is it a string, not a file name?
            saved_nexus_doFastNextTok = var.nexus_doFastNextTok
            if var.warnReadNoFile:
                print "\nread()"
                print "    A file by the name specified by the argument cannot be found."
                print "    So I am assuming that it is to be taken as a string."
                print "    Maybe it was a mis-specified file name?"
                print "    (You can turn off this warning by turning var.warnReadNoFile off.)\n"
            var.nexus_doFastNextTok = False
            flob = cStringIO.StringIO(stuff)
            _decideFromContent('<input string>', flob)
            #_decideFromContent(stuff, flob)
            #if var.verboseRead:
            #    print "(You can turn off these messages by turning var.verboseRead off.)\n"
            var.nexus_doFastNextTok = saved_nexus_doFastNextTok



def readFile(fName):
    """If its a data or tree file, read it.  If its python code, execfile it."""

    gm = ['func.readFile(%s)' % fName]
    #print gm
    #print 'func.readFile().  nexus_doFastNextTok=%s' % var.nexus_doFastNextTok
    # I should check if the file is a text file, an executable, or whatever.
    try:
        flob = file(fName, "U") # Universal line endings.
    except IOError:
        gm.append("Can't open %s.  Are you sure you have the right name?" % fName)
        raise Glitch, gm

    
    # See if there is an informative suffix on the file name
    # If there is a suffix, but the file cannot be read,
    # it is a serious error, and death follows.
    result = re.search('(.+)\.(.+)', fName)
    if result:
        #baseName = result.group(1)
        #print "got result.group(2) = %s" % result.group(2)
        suffix = string.lower(result.group(2))
        #print "readFile: got suffix '%s'" % suffix
        if suffix == 'py':
            flob.close()
            import __main__
            execfile(fName, __main__.__dict__,  __main__.__dict__)
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
                raise Glitch, gm
            return
        elif suffix == 'gde':
            ret = _tryToReadGdeFile(fName, flob)
            if not ret:
                gm.append("Failed to read supposed gde file '%s'" % fName)
                raise Glitch, gm
            return
        elif suffix in ['pir', 'nbrf']:  # new, july 2010
            ret = _tryToReadPirFile(fName, flob)
            if not ret:
                gm.append("Failed to read supposed pir file '%s'" % fName)
                raise Glitch, gm
            return
        elif suffix == 'phy' or suffix == 'phylip':
            ret = _tryToReadPhylipFile(fName, flob, None)
            if not ret:
                gm.append("Failed to read supposed phylip file '%s'" % fName)
                raise Glitch, gm
            return
        elif suffix == 'aln':
            ret = _tryToReadClustalwFile(fName, flob)
            if not ret:
                gm.append("Failed to read supposed clustalw file '%s'" % fName)
                raise Glitch, gm
            return
        elif result.group(2) in ['p4_tPickle']: # preserve uppercase
            if var.verboseRead:
                print "Trying to read '%s' as a pickled Tree file..." % fName
            import cPickle
            ret = cPickle.load(flob)
            if not ret:
                gm.append("Failed to read supposed p4_tPickle file '%s'." % fName)
                raise Glitch, gm
            else:
                if isinstance(ret, Tree):
                    ret.fName = fName
                    var.trees.append(ret)
                    var.fileNames.append(fName)
                    if var.verboseRead:
                        print "Got a tree from file '%s'." % fName
                else:
                    gm.append("Failed to get a Tree from '%s'" % fName)
                    raise Glitch, gm
            return
        else:
            _decideFromContent(fName, flob)
    else:
        _decideFromContent(fName, flob)

    #if var.verboseRead:
    #    print "(You can turn off these messages by turning var.verboseRead off.)\n"




def _decideFromContent(fName, flob):
    gm = ["func._decideFromContent()"]

    firstLine = False
    while not firstLine:
        firstLine = flob.readline()
        if not firstLine: # end of the file
            break
        firstLine = firstLine.strip()
        if firstLine: # got some content
            break
        else:
            #print "blank line"
            pass
    #if firstLine:
    #    print "got firstLine = %s" % firstLine
    #else:
    #    print "Got a blank file."
    if not firstLine:
        gm.append("Input '%s' is empty." % fName)
        raise Glitch, gm
    else:
        # Fasta, phylip, and clustalw files have clues on the first line
        # If these files fail to be what they are first supposed to be,
        # then the program does not die, it just complains.
        # News, May 2001.  I will allow blank lines at the beginning.

        if firstLine[0] == '>' or firstLine[0] == ';':
            if var.verboseRead:
                print "Guessing that '%s' is a fasta file..." % fName
            ret = _tryToReadFastaFile(fName, flob, firstLine)
            if not ret:
                if var.verboseRead:
                    print "Failed to read '%s' as a fasta file." % fName
            if ret:
                return

        elif firstLine[0] == 'C':
            if var.verboseRead:
                print "First letter is 'C'-- guessing that '%s' is a clustalw file..." % fName
            ret = _tryToReadClustalwFile(fName, flob, firstLine)
            if not ret:
                if var.verboseRead:
                    print "Failed to read '%s' as a clustalw file." % fName
            if ret:
                return

        else: # Maybe its a phylip file?
            if var.verboseRead:
                print "Guessing that '%s' is a phylip file..." % fName
            # either data or trees
            ret = _tryToReadPhylipFile(fName, flob, firstLine) 
            if not ret:
                if var.verboseRead:
                    print "Failed to read '%s' as a phylip file." % fName
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
        #print "Got firstChar candidate '%s'" % firstChar
        if firstChar == '':  # Redundant: this problem seems to be caught above, in _tryToReadPhylipFile()
            gm.append("This file seems to be composed only of whitespace.")
            raise Glitch, gm

        #print "got firstChar = %s" % firstChar
        flob.seek(0)
        if firstChar in ['[', '#']:
            # it might be a nexus file
            if var.verboseRead: print "Guessing that '%s' is a nexus file..." % fName
            ret = _tryToReadNexusFile(fName, flob)
            if ret:
                return
        elif firstChar == '{':
            if var.verboseRead: print "Guessing that '%s' is a gde file..." % fName
            ret = _tryToReadGdeFile(fName, flob)
            if ret:
                return

        # If we are here then it isn't a nexus or gde file.
        # Last alternative: its assumed to be a python script
        if var.verboseRead:
            #print "As a last resort, guessing that '%s' is a python file..." % fName
            print "Trying to read '%s' as a python file..." % fName
        try:
            theString = flob.getvalue() # get it, for a possible error message below, while the flob is still open
        except AttributeError:
            theString = None
        flob.close()
        try:
            import __main__
            execfile(fName, __main__.__dict__, __main__.__dict__)
            if hasattr(__main__, 'pyFileCount'):
                __main__.pyFileCount += 1
        except:
            if var.verboseRead:
                gm = ["Failed to read '%s' as a python file..." % fName]
                gm.append("Giving up on trying to read '%s' at all." % fName)
            else:
                if type(flob) == type(cStringIO.StringIO('foo')): # probably a mis-spelled file name
                    if theString:
                        if len(theString) < 100:
                            gm = ["Couldn't make sense out of the input '%s'" % theString]
                        else:
                            gm = ["Couldn't make sense out of the input '%s ...'" % theString[100]]
                    gm.append("It is not a file name, and I could not make sense out of it otherwise.")
                else:
                    gm=["Couldn't make sense of the input '%s'." % fName]
            raise Glitch, gm



def _tryToReadNexusFile(fName, flob):
    if var.verboseRead:
        print "Trying to read '%s' as a nexus file..." % fName
    nf = Nexus()

    # nf.readNexusFile()
    #    returns -1 if it does not start with a #nexus token
    #    returns 1 otherwise
    ret = nf.readNexusFile(flob)
    #print "Nexus.readNexusFile() returned a %s" % ret
    if ret == -1:
        if var.verboseRead:
            print "Failed to get '%s' as a nexus file." % fName
    else:
        if 0:
            print "\n******************************************\n"
            nf.dump()
            print "\n******************************************\n"
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
                    tnList.append(string.lower(t.name))
                for t in var.trees:
                    if tnList.count(string.lower(t.name)) > 1:
                        print "P4: Warning: got duplicated tree name '%s' (names are compared in lowercase)" \
                              % t.name
                        #print "Lowercased tree names: %s" % tnList
                        #sys.exit()
        if hasattr(flob, 'name'):
            var.fileNames.append(flob.name)
        if var.verboseRead:
            print "Got nexus file '%s'" % fName
        return 1


def _tryToReadFastaFile(fName, flob, firstLine=None):
    if not firstLine:
        firstLine = False
        while not firstLine:
            firstLine = flob.readline()
            if not firstLine: # end of the file
                break
            firstLine = firstLine.strip()
            if firstLine: # got some content
                break
            else:
                #print "blank line"
                pass
        #if firstLine:
        #    print "got firstLine = %s" % firstLine
        #else:
        #    print "Got a blank file."
    if not firstLine:
        gm = ["_tryToReadFastaFile: the file '%s' is empty!" % fName]
        raise Glitch, gm
    else:
        if var.verboseRead:
            print "Trying to read '%s' as a fasta file..." % fName
        if len(firstLine) > 1:
            if firstLine[0] == '>' or firstLine[0] == ';':
                flob.seek(0)
                try:
                    # this code is useless -- it always succeeds, never excepts.
                    sl = SequenceList(flob)
                except:
                    if var.verboseRead:
                        print "Reading it as a fasta file didn't work (no 'SequenceList' object returned)"
                    return None
                else:
                    if hasattr(flob, 'name'):
                        sl.fName = flob.name
                        var.fileNames.append(flob.name)
                    sl.checkNamesForDupes()

                    # If we have equal sequence lengths, then it might be an alignment
                    hasEqualSequenceLens = True
                    if len(sl.sequences) <= 1:
                        hasEqualSequenceLens = None # ie not applicable
                    else:
                        len0 = len(sl.sequences[0].sequence)
                        for s in sl.sequences[1:]:
                            if len(s.sequence) != len0:
                                hasEqualSequenceLens = False
                                
                        
                    if not hasEqualSequenceLens:
                        if var.verboseRead:
                            print "The sequences appear to be different lengths"
                        var.sequenceLists.append(sl)
                    else:
                        if var.verboseRead:
                            print "The sequences appear to be all the same length"
                        try:
                            a = sl.alignment()    # includes a call to checkLengthsAndTypes()
                            #a.checkLengthsAndTypes()
                        except:
                            if var.verboseRead:
                                print "Its not an alignment, even tho the sequences are all the same length."
                                print "    Maybe p4 (erroneously?) thinks that the sequences are different dataTypes."
                            var.sequenceLists.append(sl)
                            if var.verboseRead:
                                print "Got fasta file '%s'."  % fName
                            return 1

                        if var.verboseRead:
                            print "The fasta file appears to be an alignment."

                        if var.doCheckForAllGapColumns:
                            a.checkForAllGapColumns()
                        if var.doCheckForBlankSequences:
                            a.checkForBlankSequences()
                        if var.doCheckForDuplicateSequences:
                            a.checkForDuplicateSequences()
                        var.alignments.append(a)
                        
                    if var.verboseRead:
                        print "Got fasta file '%s'."  % fName
                    return 1
            else:
                if var.verboseRead:
                    print "First char is neither '>' nor ';' ---not fasta"
        else:
            if var.verboseRead:
                print "First line is blank--- not fasta"



def _tryToReadPhylipFile(fName, flob, firstLine):
    #print "tryToReadPhylipFile here"
    #print "firstLine is '%s'" % firstLine
    gm = ["func._tryToReadPhylipFile()"]
    if not firstLine:
        firstLine = flob.readline()
    #print "B firstLine is '%s'" % firstLine
    if not firstLine:
        gm.append("The file %s is empty." % fName)
        raise Glitch, gm
    splitLine = firstLine.split()

    # If theres 2 numbers on the first line, it may be a phylip data file
    if len(splitLine) >= 2:
        try:
            firstNum = int(splitLine[0])
            secondNum = int(splitLine[1])
            
            if var.verboseRead:
                print "Trying to read '%s' as a phylip data file..." % fName
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
                print "Got '%s' as a phylip-like file." % fName
            return 1
        except ValueError:
            if var.verboseRead:
                print "Does not seem to be a phylip or phylip-like data file."

    # Ok, so the file did not start with 2 integers.  It might still
    # be a Phylip tree file.
    flob.seek(0,0) # Go to the beginning of the flob
    firstChar = flob.read(1)
    while firstChar and firstChar in string.whitespace:
        firstChar = flob.read(1)
    if not firstChar:
        gm.append("No non-whitespace chars found.")
        raise Glitch, gm
    #print "got firstChar '%s'" % firstChar
    if firstChar not in ['(', '[']:  # The '[' for puzzle output.
        if var.verboseRead:
            print "First char is not '(' or '['.  This does not seem to be a phylip or puzzle tree file."
        return
    if var.verboseRead:
        print "Trying to read '%s' as a phylip tree file..." % fName
    flob.seek(0,0)
    theseTrees = []

    # We need to import nextTok.
    if var.nexus_doFastNextTok:
        from NexusToken2 import nextTok
    else:
        from NexusToken import nextTok
    
    while 1:
        savedPosition = flob.tell()
        tok = nextTok(flob)  # Just to check whether there is a token that can be read ...
        if not tok:
            break
        if tok != '(':
            if var.verboseRead:
                print "First char was '%s'," % firstChar,
                print "so I thought it was a phylip or puzzle tree file."
                print "However, after having read in %i trees," % len(theseTrees)
                print " it confused me by starting a supposed new tree with a '%s'" % tok
            return
        flob.seek(savedPosition, 0) # Throw the token away.
        t = Tree()
        t.name = 't%i' % len(theseTrees)
        t.parseNewick(flob, None) # None is the translationHash
        t.initFinish()
        theseTrees.append(t)
    if len(theseTrees) == 0:
        return
    else:
        for t in theseTrees:
            t.checkDupedTaxonNames()
            #if t.taxNames:  # Why would it?
            #    t.checkTaxNames()
        var.trees += theseTrees
    if hasattr(flob, 'name'):
        var.fileNames.append(flob.name)
    if var.verboseRead:
        print "Got %i trees from phylip tree file '%s'" % (len(theseTrees), fName)
    return 1



def _tryToReadClustalwFile(fName, flob, firstLine = None):
    if not firstLine:
        firstLine = flob.readline()
    if not firstLine:
        gm.append("func. _tryToReadClustalwFile()  The file %s is empty." % fName)
        raise Glitch, gm
    expectedFirstLine = 'CLUSTAL'
    if firstLine.startswith(expectedFirstLine):
        if var.verboseRead:
            print "Trying to read '%s' as a clustalw file..." % fName
        a = Alignment()
        if hasattr(flob, 'name'):
            a.fName = flob.name
            var.fileNames.append(flob.name)
        a.readOpenClustalwFile(flob)
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
            print "Got '%s' as a clustalw file." % fName
        return 1


def _tryToReadGdeFile(fName, flob):
    if var.verboseRead:
        print "Trying to read '%s' as a gde file..." % fName
    a = Alignment()
    if hasattr(flob, 'name'):
        a.fName = flob.name
        var.fileNames.append(flob.name)
    a.readOpenGdeFile(flob)
    #a.writePhylip()
    if var.doCheckForAllGapColumns:
        a.checkForAllGapColumns()
    if var.doCheckForBlankSequences:
        a.checkForBlankSequences()
    if var.doCheckForDuplicateSequences:
        a.checkForDuplicateSequences()
    var.alignments.append(a)
    if var.verboseRead:
        print "Got '%s' as a gde file." % fName
    return 1

def _tryToReadPirFile(fName, flob):
    if var.verboseRead:
        print "Trying to read '%s' as a pir file..." % fName
    flob.seek(0)
    sl = SequenceList()
    ret = sl._readOpenPirFile(flob)
    if not ret:
        if var.verboseRead:
           print "Reading it as a pir file didn't work."
           return None
    else:
        #print "Got sl"
        if hasattr(flob, 'name'):
            sl.fName = flob.name
            var.fileNames.append(flob.name)
        sl.checkNamesForDupes()

        # If we have equal sequence lengths, then it might be an alignment
        hasEqualSequenceLens = True
        if len(sl.sequences) <= 1:
            hasEqualSequenceLens = None # ie not applicable
        else:
            len0 = len(sl.sequences[0].sequence)
            for s in sl.sequences[1:]:
                if len(s.sequence) != len0:
                    hasEqualSequenceLens = False


        if not hasEqualSequenceLens:
            if var.verboseRead:
                print "The sequences appear to be different lengths"
            var.sequenceLists.append(sl)
        else:
            if var.verboseRead:
                print "The sequences appear to be all the same length"
            try:
                a = sl.alignment()    # includes a call to checkLengthsAndTypes()
            except:
                if var.verboseRead:
                    print "Its not an alignment, even tho the sequences are all the same length."
                    print "    Maybe p4 (erroneously?) thinks that the sequences are different dataTypes."
                var.sequenceLists.append(sl)
                if var.verboseRead:
                    print "Got pir file '%s'."  % fName
                return 1

            if var.verboseRead:
                print "The pir file appears to be an alignment."

            if var.doCheckForAllGapColumns:
                a.checkForAllGapColumns()
            if var.doCheckForBlankSequences:
                a.checkForBlankSequences()
            if var.doCheckForDuplicateSequences:
                a.checkForDuplicateSequences()
            var.alignments.append(a)

        if var.verboseRead:
            print "Got pir file '%s'."  % fName
        return 1



def splash():
    """Print a splash screen for p4."""
    print ''

    from version import versionString, dateString
    print "p4 v %s, %s" % (versionString, dateString)
    print """
usage:
    p4
 or
    p4 [-i] [-x] [-d] [yourScriptOrDataFile] [anotherScriptOrDataFile ...]
 or
    p4 --help

p4 is a Python package for phylogenetics.
p4 is also the name of a Python script that loads the p4 package."""
    
    print """
There is documentation at http://p4.nhm.ac.uk """

    print """
Using the p4 script, after reading in the (optional) files on the
command line, p4 goes interactive unless one of the files on the
command line is a Python script.  Use the -i option if you want to go
interactive even if you are running a script.  Use the -x option to
force exit, even if there was no Python script read.  If you use the
-d option, then p4 draws any trees that are read in on the command
line, and then exits.

Peter Foster
The Natural History Museum, London
p.foster@nhm.ac.uk"""

    if var.examplesDir:
        print "\nSee the examples in %s" % var.examplesDir
    print ''
    print "(Control-d to quit.)\n"


        


def randomTree(taxNames=None, nTax=None, name='random', seed=None, biRoot=0, randomBrLens=1, constraints=None):
    """Make a simple random Tree.

    You can supply a list of taxNames, or simply specify nTax.  In the
    latter case the specified number of (boringly-named) leaves will
    be made.

    The default is to have 'randomBrLens', where internal nodes get
    brLens of 0.02 - 0.05, and terminal nodes get brLens of 0.02 -
    0.5.  Branch lengths are all 0.1 if randomBrLens is turned
    off.

    This method starts with a star tree and keeps adding nodes
    until it is fully resolved.  If 'biRoot' is set, it adds one
    more node, and roots on that node, to make a bifurcating root.

    Repeated calls will give different random trees, without
    having to do any seed setting.  If for some reason you want to
    make identical random trees, set the seed to some positive
    integer, or zero.

    If you want the tree to have some topological constraints, say so
    with a Constraints object.

    Returns a tree."""

    complaintHead = '\nrandomTree()'
    gm = [complaintHead]

    # we need either taxNames or nTax
    if not taxNames and not nTax:
        gm.append("You need to supply either taxNames or nTax.")
        raise Glitch, gm
    if taxNames and nTax:
        if len(taxNames) != nTax:
            gm.append("You need not supply both taxNames and nTax,")
            gm.append("but if you do, at least they should match, ok?")
            raise Glitch, gm
    if taxNames:  # implies not []
        nTax = len(taxNames)
    elif nTax:
        taxNames = []
        for i in range(nTax):
            taxNames.append('t%i' % i)

    if constraints:
        assert isinstance(constraints, Constraints)

    # Make a random list of indices for the taxNames
    import random
    if seed != None: # it might be 0
        random.seed(seed)
    indcs = range(nTax)
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
    #t.dump(node=1)
    #t.draw()
    #constraints.dump()

    nNodesAddedForConstraints = 0
    if constraints:
        for aConstraint in constraints.constraints:
            #print "doing aConstraint %s  %i" % (getSplitStringFromKey(aConstraint, nTax), aConstraint)
            #t.dump(tree=0, node=1)
            t.setPreAndPostOrder()
            eTaxNames = []
            for i in range(nTax):
                tester = 1L << i
                if tester & aConstraint: # Does aConstraint contain the tester bit?
                    eTaxNames.append(taxNames[i])
            #print "aConstraint %s" % eTaxNames
            
            # check that they all share the same parent
            firstParent = t.node(eTaxNames[0]).parent
            for tN in eTaxNames[1:]:
                if t.node(tN).parent != firstParent:
                    gm.append("constraint %s" % getSplitStringFromKey(aConstraint, constraints.tree.nTax))
                    gm.append("'%s' parent is not node %i" % (tN, firstParent.nodeNum))
                    gm.append('It appears that there are incompatible constraints.')
                    raise Glitch, gm
            
            n = Node()
            n.nodeNum = nodeNum
            nodeNum += 1
            chosenName = random.choice(eTaxNames)
            eTaxNames.remove(chosenName)
            #print 'adding a new parent for %s' % chosenName
            chosenNode = t.node(chosenName)
            chosenNodeOldSib = chosenNode.sibling
            chosenNodeOldLeftSib = chosenNode.leftSibling() # could be None
            
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
                #print "adding '%s'" % chosenName
                eTaxNames.remove(chosenName)
                chosenNode = t.node(chosenName)
                chosenNodeOldSib = chosenNode.sibling
                chosenNodeOldLeftSib = chosenNode.leftSibling()
                if 0:
                    if chosenNodeOldLeftSib:
                        print 'chosenNodeOldLeftSib = %s' % chosenNodeOldLeftSib.nodeNum
                    else:
                        print 'chosenNodeOldLeftSib = None'
                    if chosenNodeOldSib:
                        print 'chosenNodeOldSib = %s' % chosenNodeOldSib.nodeNum
                    else:
                        print 'chosenNodeOldSib = None'
                
                
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
    #sys.exit()

    # Now we have a star tree.  Now add internal nodes until it is all
    # resolved, which needs nTax - 3 nodes
    #print 'nNodesAddedForConstraints is %i' % nNodesAddedForConstraints
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
                #if n.parent == t.root and t.root.getNChildren() == 3:
                #    pass
                if n.parent == t.root and t.root.leftChild.sibling and \
                       t.root.leftChild.sibling.sibling and not \
                       t.root.leftChild.sibling.sibling.sibling:
                    pass
                else:
                    ssNodes.append(n)
        lChild = random.choice(ssNodes)

        #print "lChild = node %i" % lChild.nodeNum

        ##    +----------1:oldLeftSib
        ##    |
        ##    +----------2:lChild
        ##    0
        ##    +----------3:lChildSib
        ##    |
        ##    +----------4:oldLChildSibSib


        ##    +----------1:oldLeftSib
        ##    |
        ##    |          +----------3:lChild
        ##    0----------2(n)
        ##    |          +----------4:lChildSib
        ##    |
        ##    +----------5:oldLChildSibSib

        n = Node()
        n.nodeNum = nodeNum
        nodeNum = nodeNum + 1
        lChildSib = lChild.sibling  # guarranteed to have one
        oldLChildSibSib = lChildSib.sibling # ditto
        #oldLeftSib = lChild.parent.leftChild # first guess ...
        #if oldLeftSib != lChild:
        #    while oldLeftSib.sibling != lChild:
        #        oldLeftSib = oldLeftSib.sibling
        #else:
        #    oldLeftSib = None
        oldLeftSib = lChild.leftSibling() # could be none
        if 0:
            if oldLeftSib:
                print "oldLeftSib = %i" % oldLeftSib.nodeNum
            else:
                print "oldLeftSib = None"
            print "lChildSib = %i" % lChildSib.nodeNum
            if oldLChildSibSib:
                print "oldLChildSibSib = %i" % oldLChildSibSib.nodeNum
            else:
                print "oldLChildSibSib = None"

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
        t.preOrder = None
        t.postOrder = None
        t.preAndPostOrderAreValid = False
        t.draw()
        #t.summarizeNodes()

    if biRoot:
        # pick a random node, with a parent
        n = t.nodes[random.randrange(1, len(t.nodes))]

        # addNodeBetweenNodes() requires t.preOrder.
        if 1:
            if var.usePfAndNumpy:
                t.preOrder = numpy.array([var.NO_ORDER] * len(t.nodes), numpy.int32)
                t.postOrder = numpy.array([var.NO_ORDER] * len(t.nodes), numpy.int32)
            else:
                t.preOrder = [var.NO_ORDER] * len(t.nodes)
                t.postOrder = [var.NO_ORDER] * len(t.nodes)
            if len(t.nodes) > 1:
                t.setPreAndPostOrder()
        nodeNum = t.addNodeBetweenNodes(n, n.parent)
        t.reRoot(nodeNum, moveInternalName=False)
    else:
        # The way it is now, the root rightmost child is always a
        # leaf.  Not really random, then, right?  So choose a random
        # internal node, and re-root it there.
        #print "nTax=%i, len(t.nodes)=%i" % (nTax, len(t.nodes))
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

    t.initFinish()
    return t




def newEmptyAlignment(dataType=None, symbols=None, taxNames=None, length=None):
    """Make de novo and return an Alignment object, made of gaps.

    It is not placed in var.alignments.
    """

    complaintHead = ['\nnewEmptyAlignment()']
    gm = complaintHead
    # check for silliness
    if not dataType:
        gm.append("No dataType. You need to specify at least the dataType, taxNames, and sequenceLength.")
        raise Glitch, gm
    if not taxNames:
        gm.append("No taxNames. You need to specify at least the dataType, taxNames, and sequenceLength.")
        raise Glitch, gm
    if not length:
        gm.append("No length.  You need to specify at least the dataType, taxNames, and sequenceLength.")
        raise Glitch, gm
    goodDataTypes = ['dna', 'protein', 'standard']
    if dataType not in goodDataTypes:
        gm.append("dataType '%s' is not recognized.")
        gm.append("I only know about %s" % goodDataTypes)
        raise Glitch, gm
    if dataType == 'standard':
        if not symbols:
            gm.append("For standard dataType you need to specify symbols.")
            raise Glitch, gm
    else:
        if symbols:
            gm.append("You should not specify symbols for %s dataType." % dataType)
            raise Glitch, gm

    from Alignment import Alignment
    from SequenceList import Sequence
    a = Alignment()
    a.length = length
    a.dataType = dataType

    # dataTypes, symbols, dim
    if a.dataType == 'dna':
        a.symbols = 'acgt'
        a.dim = 4
        a.equates = {'n': 'acgt', 'm': 'ac', 'k': 'gt', # 'x': 'acgt', 
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
        s = Sequence()
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
        tester = 2L ** i
        if tester & theKey:
            ss[i] = '*'
            #ss[i] = '1'
    if escaped:
        return '\\' + string.join(ss, '\\')
    else:
        return string.join(ss, '')


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
        
    So  the split key for ``['B', 'D']`` will be 2 + 8 = 10L

    Another ::
    
        func.getSplitKeyFromTaxNames(['A', 'B', 'C', 'D'], ['B', 'D'])
        # returns 10L

    However, if the splitKey is odd, it is bit-flipped.  So if
    ``someTaxNames = ['A', 'D']``, the raw split key for ``['A','D']``
    will be 1 + 8 = 9L, binary '1001', which is then xor'd with
    1111, giving 6L.

    Another ::
    
        getSplitKeyFromTaxNames(['A', 'B', 'C', 'D'], ['A', 'D'])
        returns 6L

    """

    gm = ['func.getSplitKeyFromTaxNames()']
    if not len(allTaxNames) or not len(someTaxNames):
        gm.append("Got an empty arg?!?")
        raise Glitch, gm
    theIndices = []
    for tn in someTaxNames:
        try:
            theIndex = allTaxNames.index(tn)
        except ValueError:
            gm.append("The taxName '%s' is not in allTaxNames." % tn)
            raise Glitch, gm
        if theIndex not in theIndices:  # Duped indices would be a Bad Thing
            theIndices.append(theIndex)
    #print "theIndices = %s" % theIndices

    theRawSplitKey  = 0L
    for i in theIndices:
        theRawSplitKey += 1L << i # "<<" is left-shift

    if 1 & theRawSplitKey:  # Is it odd?  or Does it contain a 1?
        allOnes = 2L**(len(allTaxNames)) - 1
        theSplitKey = allOnes ^ theRawSplitKey  # "^" is xor, a bit-flipper.
        return theSplitKey
    else:
        return theRawSplitKey


def _sumOfRows(aList):
    """
    Adds up the rows of a 2d matrix, returning the vector.
    Eg _sumOfRows([[2,3], [6,13]]) returns [5, 19]
    """
    if type(aList[0]) != type([]):
        print "_sumOfRows: not a 2D array.  Assume its a row vector and return sum"
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
    if type(aList[0]) != type([]):
        print "_sumOfColumns: not a 2D array.  Assume its a column vector and return sum"
        return sum(aList)
    theLen = len(aList[0])
    for i in aList:
        if theLen != len(i):
            print "_sumOfColumns: unequal rows"
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
    #print theExpected
    for i in theSumOfRows:
        if i == 0.0:
            gm.append("Sum of rows includes a zero.  Can't calculate xSquared.")
            raise Glitch, gm
    for i in theSumOfCols:
        if i == 0.0:
            gm.append("Sum of cols includes a zero.  Can't calculate xSquared.")
            raise Glitch, gm

    xSq = 0.0
    for i in range(nRows):
        for j in range(nCols):
            xSq = xSq + ((observed[i][j] - theExpected[i][j]) * \
                        (observed[i][j] - theExpected[i][j]) / theExpected[i][j])
    return xSq


def variance(seq):
    """This would not be good for a lot of data. n - 1 weighted."""
    sumSeq = float(sum(seq))
    return (_sumOfSquares(seq) - ((sumSeq * sumSeq) / len(seq))) / (len(seq) - 1)
    #return (_sumOfSquares(seq) - ((sumSeq * sumSeq) / len(seq))) / len(seq)

def _stdErrorOfTheDifferenceBetweenTwoMeans(seq1, seq2):
    """This could use some re-coding to handle short (<30) n"""
    if len(seq1) == len(seq2) and len(seq1) > 30:
        return math.sqrt((variance(seq1) + variance(seq2)) / len(seq1))
    else:
        gm = ["_stdErrorOfTheDifferenceBetweenTwoMeans()"]
        gm.append("I can only deal with sequences of equal length, each more than 30 long.")
        raise Glitch, gm

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
    than or equal to theStat.  theDistribution need not be sorted."""

    gm = ["tailAreaProbability()"]
    theLen = len(theDistribution)
    for i in theDistribution:
        try:
            float(i)
        except TypeError:
            gm.append("Item '%s' from theDistribution does not seem to be a float." % i)
            raise Glitch, gm
    try:
        float(theStat)
    except TypeError:
        gm.append("theStat '%s' does not seem to be a float." % theStat)
        raise Glitch, gm
        
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
    tap = float(hits)/float(theLen)
    if verbose:
        print "# The stat is %s" % theStat
        print "# The distribution has %i items" % len(theDistribution)
        print "# The distribution goes from %s to %s" % (theMin, theMax)
        print "# Items in the distribution were >= theStat %i times." % hits
        print "# The tail-area probability is %f" % tap
    return tap


def ls():
    """Like the shell ls
    """
    import os
    fList = os.listdir('.')
    fList.sort()
    for f in fList:
        print f

def which(what, verbose=0):
    """Asks if an auxiliary program is available.

    This uses the shell command 'which' to find whether a program (as
    given by the argument 'what') is in the path.  It returns 0 or 1.
    If verbose is turned on, it speaks the path, if it exists."""

    if type(what) != type('aString'):
        raise Glitch, "function which().  I was expecting a string argument."
    import os
    f = os.popen('which %s 2> /dev/null' % what, 'r')
    aLine = f.readline()
    f.close()
    if aLine:
        aLine = aLine[:-1]
        if aLine.endswith('Command not found.'): # tcsh does this, but I have not tested this part.
            aLine = None
    if aLine:
        if verbose:
            print aLine
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


##Ignore
def writeRMatrixTupleToOpenFile(theTuple, dim, flob, offset=23):
    gm = ["func.writeRMatrixTupleToOpenFile()"]
    if dim < 2:
        gm.append("dim must be two or more for this to work.")
        raise Glitch, gm

    isShort = False # For backward compatibility
    if var.rMatrixNormalizeTo1:
        isShort = False
        if len(theTuple) != (((dim * dim) - dim) / 2):
            if len(theTuple) == ((((dim * dim) - dim) / 2) - 1):
                isShort = True
            else:
                gm.append("var.rMatrixNormalizeTo1 is %i" % var.rMatrixNormalizeTo1)
                gm.append("The length of the tuple (%i) is " % len(theTuple))
                gm.append("incommensurate with the dim (%i)" % dim)
                gm.append("(should be %i)" % (((dim * dim) - dim) / 2))
                raise Glitch, gm
    else:
        isShort = True
        if len(theTuple) != (((dim * dim) - dim) / 2) - 1:
            gm.append("var.rMatrixNormalizeTo1 is %i" % var.rMatrixNormalizeTo1)
            gm.append("The length of the tuple (%i) is " % len(theTuple))
            gm.append("incommensurate with the dim (%i)" % dim)
            gm.append("(should be %i)" % ((((dim * dim) - dim) / 2) - 1))
            raise Glitch, gm
        
    if dim == 3:
        flob.write('%s\n' % repr(theTuple))
        
    else:
        #if var.rMatrixNormalizeTo1:
        if not isShort:  
            #print "theTuple=%s" % theTuple
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
            #print "dim=%i, row=%i, range((dim-row) - 1) = %s, tuplePos=%i" % (dim, row, range((dim-row)-1), tuplePos)
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

##Ignore
def writeCharFreqToOpenFile(theCharFreq, dim, symbols, flob, offset=23):  ##Ignore
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
            flob.write('     [%c] %f)\n' % (symbols[dim - 2], theCharFreq[dim - 2]))
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
    #print l
    inMath = False
    for i in range(len(l)):
        if l[i] in string.letters or l[i] in string.digits:
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
    if l[0] == '\'' and l[-1] == '\'':
        del(l[-1])
        del(l[0])
        theRange = range(len(l))
        theRange.reverse()
        for i in theRange[:-1]:
            if l[i] == '\'' and l[i - 1] == '\'':
                del(l[i])
    return string.join(l, '')



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
    from Alignment_muck import cListPat, cList2Pat, cListAllPat
    #cListPat = re.compile('(\d+)-?(.+)?')
    #cList2Pat = re.compile('(.+)\\\\(\d+)')
    #cListAllPat = re.compile('all\\\\?(\d+)?')

    #print "char list is: %s" % nexusCharListString
    cList = string.split(nexusCharListString)
    #print "cList is %s" % cList
    import array
    mask = array.array('c', maskLength * '0')
    for c in cList:        # eg 6-10\2
        first = None       # the first item eg 6
        second = None      # the second item eg 10
        third = None       # the third item eg 2
        result = cListPat.match(c)
        if result:
            # print "%s\t%s" % (result.group(1), result.group(2))
            first = result.group(1)
            if result.group(2):
                r2 = cList2Pat.match(result.group(2))
                if r2:
                    # print "%s\t%s" % (r2.group(1), r2.group(2))
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
                gm.append("Can't parse maskFromNexusCharacterList '%s'" % nexusCharListString)
                raise Glitch, gm
        # print "first = %s, second = %s, third = %s" % (first, second, third)
        if not first:
            gm.append("Can't parse maskFromNexusCharacterList '%s'" % nexusCharListString)
            raise Glitch, gm
        elif first and not second: # its a single
            if string.lower(first) == 'all':
                for i in range(len(mask)):
                    mask[i] = '1'
            elif first == '.':
                mask[-1] = '1'
            else:
                try:
                    it = int(first)
                    mask[it - 1] = '1'
                except ValueError:
                    gm.append("Can't parse '%s' in maskFromNexusCharacterList '%s'" \
                          % (first, nexusCharListString))
                    raise Glitch, gm
        elif first and second:  # its a range
            try:
                start = int(first)
            except ValueError:
                gm.append("Can't parse '%s' in maskFromNexusCharacterList '%s'" \
                      % (first, nexusCharListString))
                raise Glitch, gm
            if second == '.':
                fin = len(mask)
            else:
                try:
                    fin = int(second)
                except ValueError:
                    gm.append("Can't parse '%s' in maskFromNexusCharacterList '%s'" % \
                                    (second, nexusCharListString))
                    raise Glitch, gm
            if third:
                try:
                    bystep = int(third)
                except ValueError:
                    gm.append("Can't parse '%s' in maskFromNexusCharacterList '%s'" % \
                          (third, nexusCharListString))
                    raise Glitch, gm
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
    return mask.tostring()


def polar2square(angleLenList):
    """Convert a coord in polar coords to usual (Cartesian? square?) coords.

    Input is a list composed of the angle in radians and the length
    (ie from the origin).  A list of [x,y] is returned.
    """
    if angleLenList[1] == 0.0:
        return [0, 0]
    elif angleLenList[1] < 0.0:
        raise Glitch, 'func.polar2square error: len is less than zero'
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
        raise Glitch, "n should be (at least convertible to) an int."
    assert n >= 0, "n should be zero or positive."

    fact = {0: 1, 1: 1, 2: 2, 3: 6, 4: 24, 5: 120, 6: 720, 7: 5040, 8: 40320, 9: 362880,
            10: 3628800, 11: 39916800, 12: 479001600, 13: 6227020800L, 14: 87178291200L,
            15: 1307674368000L, 16: 20922789888000L, 17: 355687428096000L, 18: 6402373705728000L,
            19: 121645100408832000L, 20: 2432902008176640000L, 21: 51090942171709440000L,
            22: 1124000727777607680000L, 23: 25852016738884976640000L, 24: 620448401733239439360000L,
            25: 15511210043330985984000000L, 26: 403291461126605635584000000L,
            27: 10888869450418352160768000000L, 28: 304888344611713860501504000000L,
            29: 8841761993739701954543616000000L, 30: 265252859812191058636308480000000L}
    if n <= 30:
        return fact[n]
    else:
        total = 1L
        while n > 1:
            total *= n
            n -= 1
        return total
        
def nChooseK(n, k):
    """Get the number of all possible len k subsets from range(n)."""
    try:
        n = int(n)
        k = int(k)
    except ValueError, TypeError:
        raise Glitch, "n and k should be (at least convertible to) ints."
    
    assert n >= 0, "n should be zero or more."
    assert k <= n, "k should be less than or equal to n"
    assert k >= 0, "k should be zero or more."

    nFact = factorial(n)
    kFact = factorial(k)
    nMinusKFact = factorial(n-k)
    return nFact / (kFact * nMinusKFact)

def nUnrootedTrees(nTaxa):
    upper = (nTaxa * 2) - 5
    nTrees = 1L
    i = 3
    while i <= upper:
        nTrees *= i
        i += 2
    return nTrees

def nRootedTrees(nTaxa):
    upper = (nTaxa * 2) - 3
    nTrees = 1L
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
    [0L, 1L, 246L, 6825L, 56980L, 190575L, 270270L, 135135L].  """
    
    # Make a table.
    table = []
    for nInts in range(nTaxa):
        table.append([0L] * (nTaxa + 1))
    for nT in range(2, nTaxa + 1):
        table[1][nT] = 1L
    for nInt in range(2, nTaxa):
        for nTx in range(nInt + 1, nTaxa + 1):
            table[nInt][nTx] = ((nTx + nInt - 2) * table[nInt - 1][nTx -1]) + (nInt * table[nInt][nTx - 1])

    if 0:
        # Print out the table.
        print "%-9s|" % "nTx ->",
        for nTx in range(1,nTaxa + 1):
            print "%10i" % nTx,
        print
        for nTx in range(nTaxa + 1):
            print "%10s" % " ---------",
        print
        print "%8s |" % "nInt"


        for nInt in range(1, nTaxa):
            print "%8i |" % nInt,
            for nTx in range(nInt):
                print "%10s" % "",
            for nTx in range(nInt + 1, nTaxa + 1):
                print "%10i" % table[nInt][nTx],
            print
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
    [0L, 1L, 246L, 6825L, 56980L, 190575L, 270270L, 135135L].  """

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
            theSeq = [random.gammavariate((inSeq[i] * alpha) + u, 1.0) for i in range(kk)]
        else:
            theSeq = [random.gammavariate(inSeq[i] * alpha, 1.0) for i in range(kk)]
        #print safety, theSeq
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
            gm.append("Tried more than %i times to get good dirichlet values, and failed.  Giving up." % safetyLimit)
            gm.append("inSeq: %s" % inSeq)
            gm.append("theMin: %s, theMax: %s, u=%s" % (theMin, theMax, u))
            raise Glitch, gm




    

    


def unPickleMcmc(runNum, theData, verbose=True):
    """Unpickle a checkpoint, return an Mcmc ready to go."""
    
    gm = ["func.unPickleMcmc()"]
    try:
        runNum = int(runNum)
    except (ValueError, TypeError):
        gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
        raise Glitch, gm
    if runNum < 0:
        gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
        raise Glitch, gm
    baseName = "mcmc_checkPoint_%i." % runNum
    ff = glob.glob("%s*" % baseName)

    pickleNums = []
    for f in ff:
        genNum = int(f.split('.')[1])
        pickleNums.append(genNum)
    if not pickleNums: # an empty sequence
        gm.append("Can't find any checkpoints for runNum %i." % runNum)
        gm.append("Got the right runNum?")
        raise Glitch, gm
    theIndx = pickleNums.index(max(pickleNums))
    fName = ff[theIndx]
    if verbose:
        print "...unpickling Mcmc in %s" %  fName

    import cPickle
    f = file(fName)
    m = cPickle.load(f)
    f.close()

    # 30 Aug 2011.  Fix for backward compatibility.  Need to add
    # modelPart.rjComp_k, comp.rj_f, and comp.rj_isInPool for pickles
    # from before those attributes were added.
    for myPart in m.tree.model.parts:
        if not hasattr(myPart, 'rjComp_k'):
            myPart.rjComp_k = 1
        for myComp in myPart.comps:
            if not hasattr(myComp, 'rj_f'):
                myComp.rj_f = 0.0
            if not hasattr(myComp, 'rj_isInPool'):
                myComp.rj_isInPool = False
    for myChain in m.chains:
        for myTree in [myChain.curTree, myChain.propTree]:
            for  myPart in myTree.model.parts:
                if not hasattr(myPart, 'rjComp_k'):
                    myPart.rjComp_k = 1
                for myComp in myPart.comps:
                    if not hasattr(myComp, 'rj_f'):
                        myComp.rj_f = 0.0
                    if not hasattr(myComp, 'rj_isInPool'):
                        myComp.rj_isInPool = False

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
    
    return m

def unPickleSTMcmc(runNum, verbose=True):
    """Unpickle a STMcmc checkpoint, return an STMcmc ready to go."""
    
    gm = ["func.unPickleSTMcmc()"]
    try:
        runNum = int(runNum)
    except (ValueError, TypeError):
        gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
        raise Glitch, gm
    if runNum < 0:
        gm.append("runNum should be an int, 0 or more.  Got %s" % runNum)
        raise Glitch, gm
    baseName = "mcmc_checkPoint_%i." % runNum
    ff = glob.glob("%s*" % baseName)

    pickleNums = []
    for f in ff:
        genNum = int(f.split('.')[1])
        pickleNums.append(genNum)
    if not pickleNums: # an empty sequence
        gm.append("Can't find any checkpoints for runNum %i." % runNum)
        gm.append("Got the right runNum?")
        raise Glitch, gm
    theIndx = pickleNums.index(max(pickleNums))
    fName = ff[theIndx]
    if verbose:
        print "...unpickling Mcmc in %s" %  fName

    import cPickle
    f = file(fName)
    m = cPickle.load(f)
    f.close()

    if m.stRFCalc == 'fastReducedRF':
        import pyublas      # needed
        for chNum in range(m.nChains):
            ch = m.chains[chNum]
            ch.startFrrf()
        
    return m


#################################################

def recipes(writeToFile=var.recipesWriteToFile):
    """Reminders, suggestions, multi-step methods...

    These should be in the Sphinx docs also.
    """
    
    gm = ["func.recipes()"]
    if not var.examplesDir:
        gm.append("Can't find the Examples directory.")
        raise Glitch, gm
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
    #print fList
    bNames = [os.path.basename(nm) for nm in fList]
    #print bNames
    firstLines = []
    for fN in fList:
        f = file(fN)
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
    #print recipesList
    
    for recNum in range(len(recipesList)):
        rec = recipesList[recNum]
        print "%s.  %s" % (string.uppercase[recNum], rec[0])
    ret = raw_input('Tell me a letter: ')
    #print "Got %s" % ret
    if ret == '':
        return
    elif ret[0] in string.uppercase:
        recNum = string.uppercase.index(ret[0])
    elif ret[0] in string.lowercase:
        recNum = string.lowercase.index(ret[0])
    else:
        return
    #print "Got recNum %i" % recNum

    if writeToFile:
        ret = raw_input('Write it to file name [default %s] :' % recipesList[recNum][1])
        if ret == '':
            theFName = recipesList[recNum][1]
        else:
            theFName = string.strip(ret)
        if os.path.exists(theFName):
            "The file %s already exists.  I'm refusing to over-write it, so I'm doing nothing." % theFName
            return
        print "Writing to file '%s' ..." % theFName
        os.system("cp %s %s" % (os.path.join(recipesDir, recipesList[recNum][1]), theFName))
    else:
        print "\n"
        os.system("cat %s" % os.path.join(recipesDir, recipesList[recNum][1]))
    
        
        
    
    
def uninstall():
    """Uninstall the p4 package."""

    print """
This function is for uninstalling an installed version of p4.  It uses
the p4.installation module to get the locations of files.  It may need
to be done as root, or using sudo."""

    weAreInteractive = False
    if os.getenv('PYTHONINSPECT'):   # The p4 script sets this
        weAreInteractive = True
    elif len(sys.argv) == 1 and sys.argv[0] == '':  # Command 'p4' with no args
        weAreInteractive = True
    # How can I tell if python was invoked with the -i flag?
    if not weAreInteractive:
        return
    try:
        import installation
    except ImportError:
        raise Glitch, "Unable to import the p4.installation module."
    print """
    
This function will remove the p4 script file
  %s
the library files in the directory
  %s
and the documentation in the directory
  %s

""" % (installation.p4ScriptPath, installation.p4LibDir, installation.p4DocDir)
    
    ret = raw_input('Ok to do this? [y/n]')
    ret = string.lower(ret)
    if ret not in ['y', 'yes']:
        return
    print "Ok, deleting ..."
    if os.path.exists(installation.p4LibDir):
        os.system("rm -fr %s" % installation.p4LibDir)
    else:
        raise Glitch, "Could not find %s" % installation.p4LibDir
    if os.path.exists(installation.p4DocDir):
        os.system("rm -fr %s" % installation.p4DocDir)
    else:
        raise Glitch, "Could not find %s" % installation.p4DocDir
    if os.path.exists(installation.p4ScriptPath):
        os.system("rm %s" % installation.p4ScriptPath)
    else:
        raise Glitch, "Could not find %s" % installation.p4ScriptPath
    
    

###########################################################
### Tkinter stuff
###########################################################

####def startTkThread():
####    if var.tk_thread_running:
####        print "Tk thread is already running, it appears."
####        return
####    import thread
####    import atexit
####    from Queue import Queue

####    var.tk_request = Queue(0)
####    var.tk_result = Queue(1)
    
####    thread.start_new_thread(_tk_thread,())
####    var.tk_thread_running = True
####    atexit.register(_tkShutdown)

##def _tk_thread():
##    import Tkinter
##    #print "_tk_thread() here!"
##    var.tk_root = Tkinter.Tk()
##    #print "var.tk_root is %s" % var.tk_root
##    var.tk_root.withdraw()
##    var.tk_root.after(var.tk_pollInterval, _tk_pump)
##    var.tk_root.mainloop()

##def _tk_pump():
##    #global _thread_running
##    while not var.tk_request.empty():
##        command,returns_value = var.tk_request.get()
##        try:
##            result = command()
##            if returns_value:
##                var.tk_result.put(result)
##        except:
##            var.tk_thread_running = False
##            if returns_value:
##                var.tk_result.put(None) # release client
##            raise # re-raise the exception -- kills the thread
##    if var.tk_thread_running:
##        var.tk_root.after(var.tk_pollInterval, _tk_pump)


##def _tkShutdown():
##    # shutdown the tk thread
##    #global _thread_running
##    #_tkExec(sys.exit)
##    var.tk_thread_running = False
##    time.sleep(.5) # give tk thread time to quit


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
    
    f = file(fName, "U") # Universal line endings.
    ll = f.readlines()
    f.close()

    f = file(oFName, 'w')
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

    fd, filename = tempfile.mkstemp(suffix, prefix+"_", dirname)
    return os.fdopen(fd, 'w'), filename

def writeInColour(theString, colour='blue'):
    goodColours = ['red', 'RED', 'blue', 'BLUE', 'cyan', 'CYAN', 'violet', 'VIOLET']
    if colour not in goodColours:
        raise Glitch, "func.printColour().  The colour should be one of %s" % goodColours
    codeDict = {
        'red': '\033[0;31m', 
        'RED':'\033[1;31m',
        'blue':'\033[0;34m',
        'BLUE':'\033[1;34m',
        'cyan':'\033[0;36m',
        'CYAN':'\033[1;36m',
        'violet':'\033[0;35m',
        'VIOLET':'\033[1;35m',
        }
    backToBlackCode = '\033[m'
    sys.stdout.write("%s%s%s" % (codeDict[colour], theString, backToBlackCode))

def setTerminalColour(theColour):
    goodTerminalColours = ['red', 'RED', 'blue', 'BLUE', 'cyan', 'CYAN', 'violet', 'VIOLET']
    terminalColourCodeDict = {
        'red': '\033[0;31m', 
        'RED':'\033[1;31m',
        'blue':'\033[0;34m',
        'BLUE':'\033[1;34m',
        'cyan':'\033[0;36m',
        'CYAN':'\033[1;36m',
        'violet':'\033[0;35m',
        'VIOLET':'\033[1;35m',
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
    """Pure Python, using reduce."""
    
    def addSq(x, y): return float(x) + (float(y) * float(y))
    return reduce(addSq, seq, 0)



def sortListOfObjectsOnAttribute(aListOfObjects, attributeString):
    """Returns a new sorted list."""
    
    def pairing(anObject, a=attributeString):
        return (getattr(anObject, a), anObject)
    paired = map(pairing, aListOfObjects)
    #print paired
    # 'paired' is a list of 2-element tuples, (thingToSortOn, theOriginalObject)
    paired.sort()
    def stripit(pair):
        return pair[1]
    return map(stripit, paired)

def sortListOfObjectsOn2Attributes(aListOfObjects, attributeString1, attributeString2):
    """Returns a new sorted list."""
    
    def tripling(anObject, a=attributeString1, b=attributeString2):
        return (getattr(anObject, a), getattr(anObject, b), anObject)
    tripled = map(tripling, aListOfObjects)
    # 'tripled' is a list of 3-element tuples, (thingToSortOn, secondThingToSortOn, theOriginalObject)
    tripled.sort()
    def stripit(triple):
        return triple[2]
    return map(stripit, tripled)


def sortListOfListsOnListElementNumber(aListOfLists, elementNumber):
    """Returns a new sorted list."""
    
    def pairing(aList, e = elementNumber):
        return (aList[e], aList)
    paired = map(pairing, aListOfLists)
    # 'paired' is a list of 2-element tuples, (thingToSortOn, theOriginalElement)
    # print "paired = ", paired
    paired.sort()
    def stripit(pair):
        return pair[1]
    return map(stripit, paired)


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
        raise Glitch, gm
    if nnSL:
        return var.sequenceLists.pop()
    elif nnAlig:
        return var.alignments.pop()
    elif nnTrees:
        return var.trees.pop()

###############################################################################################################
###############################################################################################################
###############################################################################################################
    

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
        gm.append("len of names (%i) is not the same as the len of lengths (%i)" % (len(names), len(lens)))
        print "names: ", names
        print "lens: ", lens
        raise Glitch, gm
                  
    start = 1
    if fName:
        f = file(fName, 'w')
    else:
        f = sys.stdout

    f.write("#nexus\n\n")
    f.write("begin sets;\n")
    for cNum in range(len(names)):
        f.write("  charset %s = %i - %i;" % (names[cNum], start, start - 1 + lens[cNum]))
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


if var.usePfAndNumpy:
    
    def reseedCRandomizer(newSeed):
        """Set a new seed for the c-language random() function.

        For those things in the C-language that use random(), this
        re-seeds the randomizer.  Reseed to different integers to make
        duplicate runs of something come out different.  This is not for
        the GSL random stuff.  Confusing to have 2 systems, innit?  And
        then there is the 3rd system that Python uses in the random
        module.  Sorry!"""

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

        if type(seq) == numpy.ndarray:
            mySeq = seq
        else:
            mySeq = numpy.array(seq, numpy.float)

        if type(mean) == types.NoneType:
            mean = numpy.array([0.0])
        if type(variance) == types.NoneType:
            variance = numpy.array([0.0])
        pf.gsl_meanVariance(mySeq, len(seq), mean, variance)
        if 0:
            from p4 import func
            print "slow p4: mean=%f, variance=%f" % (func.mean(list(seq)), func.variance(list(seq)))
            print "gsl: mean=%f, variance=%f" % (mean, variance)
            print "numpy: mean=%f, variance=%f" % (mySeq.mean(), mySeq.var())  # different than gsl-- no n-1 weighting.
        return (mean, variance)

    def chiSquaredProb(xSquared, dof):
        """Returns the probability of observing X^2."""
        import pf
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
            raise Glitch, gm

        isNewGSL_RNG = 0
        if not var.gsl_rng:
            var.gsl_rng = pf.get_gsl_rng()
            isNewGSL_RNG = 1
            #print "got var.gsl_rng = %i" % var.gsl_rng
            #sys.exit()

            # Set the GSL random number generator seed, only if it is a new GSL_RNG
            if isNewGSL_RNG:
                if seed != None:
                    try:
                        newSeed = int(seed)
                        pf.gsl_rng_set(var.gsl_rng, newSeed)
                    except ValueError:
                        print complaintHead
                        print "    The seed should be convertible to an integer"
                        print "    Using the process id instead."
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
        assert type(inSeq) == numpy.ndarray
        assert type(outSeq) == numpy.ndarray
        kk = len(inSeq)
        theMax = 1.0 - ((kk - 1) * theMin)
        safety = 0
        safetyLimit = 300
        while 1:
            for i in range(kk):
                #theSeq[i] = pf.gsl_ran_gamma(var.gsl_rng, theSeq[i] * alpha, 1.0)
                #outSeq[i] = random.gammavariate(inSeq[i] * alpha, 1.0)   --- this is very slow!
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
                gm.append("Tried more than %i times to get good dirichlet values, and failed.  Giving up." % safetyLimit)
                raise Glitch, gm
    


    def gsl_ran_dirichlet(alpha, theta, seed=None):
        """Make a random draw from a dirichlet distribution.

        Args *alpha* and *theta* should be NumPy arrays, both the same
        length (more than 1).  The length is the dimension of the
        dirichlet.  The contents of *theta* are over-written (without
        being used).  The draw ends up in *theta*.  It is normalized
        so that it sums to 1.0.

        This handles making the GSL random number generator, or
        re-using it if it was made before.  If it is newly made, you
        can optionally set its *seed*; otherwise the pid is used.
        """

        complaintHead = '\nfunc.gsl_ran_dirichlet()'
        gm = complaintHead
        try:
            assert type(alpha) == numpy.ndarray
            assert type(theta) == numpy.ndarray
        except AssertionError:
            gm.append(" alpha, theta, should be numpy arrays.")
            raise Glitch, gm

        assert len(alpha) > 1
        assert len(theta) == len(alpha)

        isNewGSL_RNG = 0
        if not var.gsl_rng:
            var.gsl_rng = pf.get_gsl_rng()
            isNewGSL_RNG = 1
            #print "got var.gsl_rng = %i" % var.gsl_rng
            #sys.exit()

            # Set the GSL random number generator seed, only if it is a new GSL_RNG
            if isNewGSL_RNG:
                if seed != None:
                    try:
                        newSeed = int(seed)
                        pf.gsl_rng_set(var.gsl_rng, newSeed)
                    except ValueError:
                        print complaintHead
                        print "    The seed should be convertable to an integer"
                        print "    Using the process id instead."
                        pf.gsl_rng_set(var.gsl_rng,  os.getpid())
                else:
                    pf.gsl_rng_set(var.gsl_rng,  os.getpid())

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
        m,v = gsl_meanVariance(sq)
        s = numpy.sqrt(v)
        n = len(sq)
        sqN = numpy.sqrt(n)
        t = (m - mu) / (s / sqN)
        dof = n - 1
        p = pf.studentsTProb(t, dof)

        if verbose:
            print "mean =", m
            print "std dev =", s
            print "t statistic =", t
            print "prob =", p
        return p


    def effectiveSampleSize(data, mean):
        """As done in Tracer v1.4, by Drummond and Rambaut.  Thanks guys!

        But see :func:`func.summarizeMcmcPrams`, which gives ESSs.
        
        """

        nSamples = len(data)
        maxLag = int(nSamples / 3)
        if maxLag > 1000:
            maxLag = 1000

        gammaStatAtPreviousLag = numpy.array([0.0])
        assert type(data) == type(gammaStatAtPreviousLag), "Arg 'data' should be a numpy.array.  Its %s" % type(data)
        assert type(mean) == type(gammaStatAtPreviousLag), "Arg 'mean' should be a numpy.array.  Its %s" % type(mean)
        gammaStat = numpy.array([0.0])
        varStat = numpy.array([0.0])
        gammaStatAtLagZero = numpy.array([0.0])
        if 0:
            lag = 0
            while lag < maxLag:
                gammaStat[0] = 0.0
                for j in range(nSamples - lag):
                    gammaStat[0] += (data[j] - mean) * (data[j + lag] - mean)

                #if lag == 0:
                #    print "mean is %f" % mean
                #    print "lag is 0, gammaStat = %f" % gammaStat[0]

                gammaStat[0] /= (nSamples - lag)

                if lag == 0:
                    varStat[0] = gammaStat
                    gammaStatAtLagZero[0] = gammaStat
                    #print "got gammaStatAtLagZero = %f" % gammaStatAtLagZero[0]
                elif (lag % 2) == 0:
                    if gammaStatAtPreviousLag + gammaStat > 0:
                        varStat[0] += 2.0 * (gammaStatAtPreviousLag + gammaStat)
                    else:
                        break
                lag += 1
                gammaStatAtPreviousLag[0] = gammaStat
                #gammaStat[0] = 0.0


            #print gammaStatAtLagZero, gammaStat, varStat, lag
            #print "maxLag is %i" % maxLag
            #stdErrorOfMean = numpy.sqrt(varStat / nSamples)  ??!?
            #ACT = stepSize * varStat / gammaStatAtLagZero
            #ESS = (stepSize * nSamples) / ACT;
            ESS1 = nSamples * (gammaStatAtLagZero / varStat)   # stepSize is not needed
            #print "got ESS1 %f" % ESS1

        if 1:
            pf.effectiveSampleSize(data, mean, nSamples, maxLag, gammaStatAtPreviousLag,
                                   gammaStat, varStat, gammaStatAtLagZero)
            ESS2 = nSamples * (gammaStatAtLagZero / varStat)
            #fabsDiff = numpy.fabs(ESS1 - ESS2)
            #print "ESS1 is %f, ESS2 is %f, diff is %g" % (ESS1, ESS2, fabsDiff)
        return ESS2[0]


    def summarizeMcmcPrams(skip=0, run=-1, theDir='.'):
        """Find the mean, variance, and ess of mcmc parameters.

        Ess is effective sample size, as in Tracer by Drummond and Rambaut.

        The numbers are found in mcmc_prams_N, N=0, 1, etc.  If arg 'run' is
        set to -1, the default, then all runs are done.  Alternatively you
        can set the run to a specific run number, and that is the only one
        that is done.

        The 'profile', with the names of the parameters, and the number of
        each, is found in mcmc_pramsProfile.py.  It is not essential, but
        it gives names to the parameters.

        """

        gm = ["func.summarizeMcmcPrams()"]
        nPrams = None
        pramsProfile = None
        try:
            loc = {}
            execfile(os.path.join(theDir, "mcmc_pramsProfile.py"), {}, loc)
            #loc =locals()  no workee.
            #print "loc = %s" % loc
            nPrams = loc['nPrams']
            pramsProfile = loc['pramsProfile']
        except IOError:
            print "The file 'mcmc_pramsProfile.py' cannot be found."

        numsList = None
        if run == -1:
            runNum = 0
        else:
            runNum = run
        totalLinesRead = 0
        while 1:
            try:
                theFName = os.path.join(theDir, "mcmc_prams_%i" % runNum)
                flob = file(theFName)
                print "Reading prams from file %s" % theFName 
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
                        if not numsList:  # If it does not exist yet, then make it now.
                            thisNPrams = len(splitLine) - 1
                            if nPrams:
                                if not thisNPrams == nPrams:
                                    gm.append("thisNPrams = %i, nPrams = %s" % (thisNPrams, nPrams))
                                    raise Glitch, gm
                            else:
                                nPrams = thisNPrams
                            numsList = []
                            for pramNum in range(nPrams):
                                numsList.append([])
                        for pramNum in range(nPrams):
                            try:
                                theOne = splitLine[pramNum + 1]
                            except IndexError:
                                gm.append("Line '%s'.  " % string.rstrip(aLine))
                                gm.append("Can't get parameter number %i  " % pramNum) 
                                raise Glitch, gm
                            try:
                                aFloat = float(theOne)
                                numsList[pramNum].append(aFloat)
                            except (ValueError, TypeError):
                                gm.append("Can't make sense of '%s'" % theOne)
                                raise Glitch, gm
                        linesRead += 1
            print "  skipped %i lines" % skipsDone
            print "  read %i lines" % linesRead
            totalLinesRead += linesRead
            if run != -1:
                break

        print "Read %i pram lines in total." % totalLinesRead

        #print numsList
        spacer1 = ' ' * 20
        if pramsProfile:
            print "%s   %16s         mean      variance       ess  " % (spacer1, ' ')
            print "%s   %16s       --------    --------    --------" % (spacer1, ' ')
            pramCounter = 0
            for partNum in range(len(pramsProfile)):
                if len(pramsProfile) > 1:
                    print "Data partition %i" % partNum
                if len(pramsProfile[partNum]):
                    #print pramsProfile[partNum]
                    for pramNum in range(len(pramsProfile[partNum])):
                        pString = pramsProfile[partNum][pramNum][0]
                        pramCounts = pramsProfile[partNum][pramNum][1]
                        for p in range(pramCounts):
                            print "%s%3i %12s[%2i]   " % (spacer1, pramCounter, pString, p),
                            d = numpy.array(numsList[pramCounter], numpy.float)
                            m,v = gsl_meanVariance(d)
                            ess = effectiveSampleSize(d, m)
                            if m == 0.0:
                                print "  0.0      ",
                            elif m < 0.00001:
                                print "%10.3g " % m,
                            elif m < 1.0:
                                print "%10.6f " % m,
                            else:
                                print "%10.4f " % m,
                            if v == 0.0:
                                print "  0.0      ",
                            elif v < 0.000001:
                                print "%10.3g " % v,
                            elif v < 1.0:
                                print "%10.6f " % v,
                            else:
                                print "%10.4f " % v,
                            print "%10.1f " % ess,
                            print
                            pramCounter += 1

                else:
                    print "        No parameters in this data partition."
        else: # no pramsProfile
            print "%9s  mean      variance       ess  " % ' '
            print "%9s--------    --------    --------" % ' '
            for pramNum in range(nPrams):
                print "  %2i  " % pramNum,
                d = numpy.array(numsList[pramNum], numpy.float)
                m,v = gsl_meanVariance(d)
                ess = effectiveSampleSize(d, m)
                if m == 0.0:
                    print "  0.0      ",
                elif m < 0.00001:
                    print "%10.3g " % m,
                elif m < 1.0:
                    print "%10.6f " % m,
                else:
                    print "%10.4f " % m,

                if v == 0.0:
                    print "  0.0      ",
                elif v < 0.000001:
                    print "%10.3g " % v,
                elif v < 1.0:
                    print "%10.6f " % v,
                else:
                    print "%10.4f " % v,

                print "%10.1f " % ess,
                print



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
            #print theHarmMean
            #print type(theHarmMean)

            return pf.newtonRaftery94_eqn16(lla, len(logLikes), theHarmMean, delta, int(verbose))

        #n = Numbers(logLikes)
        #theHarmMean = n.harmonicMeanOfLogs()  # better
        #return pf.newtonRaftery94_eqn16(lla, len(logLikes), theHarmMean, delta, int(verbose))
