import string,array
import func
from Var import var
from Alignment import Alignment
#from NexusToken import * # nextTok() et al
from NexusSets import NexusSets
from SequenceList import Sequence
from Glitch import Glitch

# Some definitions from the MadSwofMad Syst Biol Nexus format paper (MSM97).
#
# Punctuation: '\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>', and notice
#               it does not contain the period.
# Token: a nexus word or punctuation
# Symbol: a single dark character

# Word: Except for special cases involving quotes or comments, a NEXUS
# word is an string of text characters that is bounded by whitespace
# or punctuation and that does not contain whitespace or punctuation.
# If the first character of a word is a single quote, then the word
# ends with the next single quote (unless that single quote is in a
# pair of consecutive single quotes; if so, then the word ends at
# the first unpaired single quote).  Any character, including
# punctuation and whitespace, may be contained within a quoted word.  A
# word cannot consist of only whitespace and punctuation.
#

# Underscores are are considered equivalent to blank spaces, except
# that underscores are dark characters and blank spaces are
# whitespace.  Any doubled single quotes within a quoted word should
# be converted to single quotes within the program.

# Names, for nexus objects, are single nexus words; they cannot
# consist entirely of digits. (Page 597, at the top).

# SYMBOLS.  MSM97 Page 599.  "For standard datatypes, a symbols
# subcommand will replace the default symbols list of "0 1".  For DNA,
# RNA, nucleotide, and protein datatypes, a symbols subcommand will
# not replace the default symbols list, but will add character-state
# symbols to the symbols list.  

# Page 595 (point 5) "New commands in public blocks can be added ..."
# (And Mesquite has done this, adding "TITLE" commands to taxa and
# trees blocks (at least).  And LINK command to trees.

# Comments and Command comments.  Page 618.  There is an inconsistency
# in MSM97 in the definition of command comments.  In the section
# "Command comment", it says that the only symbols that start a
# command comment are '\' and '&'.  However, in the section "Comment",
# it says that command comments start with any of '&', '%', '/', '\',
# or '@'.  For simplicity, I will adopt the former.

# Bug:  Visible comments (and all comments) are ignored in blocks that
# are skipped, ie unknown blocks.

# Bug: Multi-line comments are not returned.
# Bug: Quoted tokens that are longer than the longest line may fail, so line endings within quoted tokens are not allowed.
# Bug: When var.nexus_getLineEndingsAsTokens, if line endings are \r, fast nextTok() returns \n

class Nexus:

    def __init__(self):

        # self.nexusData is a NexusData() instance.  It is a
        # container for nexus data, taxa, and characters blocks.
        self.nexusData = None
        self.nexusSets = []
        self.alignments = []
        self.trees = []

    def readNexusFile(self, flob):

        if hasattr(flob, 'name'):
            gm = ['Nexus.readNexusFile()  fName=%s' % flob.name]
        else:
            gm = ['Nexus.readNexusFile()']

        # We need to import nextTok.
        #print 'var.nexus_doFastNextTok = %s' % var.nexus_doFastNextTok
        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok
            from NexusToken2 import checkLineLengths
            checkLineLengths(flob)
        else:
            from NexusToken import nextTok
        
        # Go to the beginning and check if it starts with '#NEXUS'
        flob.seek(0)
        #print "Nexus.readNexusFile()  about to nextTok() a.  "
        tok = nextTok(flob)
        #print 'xyz tok = %s' % tok
        #sys.exit()
        if tok:
            lowTok = string.lower(tok)
        else:
            gm.append("Can't find the first token.  Is it an empty file?!?")
            raise Glitch, gm
        if lowTok != '#nexus':
            gm.append("Got first token: %s" % tok)
            gm.append("Hmmm..., this doesn't appear to be a nexus file. ")
            gm.append("The first token is not '#NEXUS'.")
            # -1 is the signal to func._tryToReadNexusFile()
            # that it is not a nexus file.  It should be non-fatal.
            # No, I changed my mind-- it should raise an error.
            #return -1
            raise Glitch, gm
        while 1:
            #print "Nexus.readNexusFile()  about to nextTok() b"
            tok = nextTok(flob)
            #print "readNexusFile: got '%s'" % tok
            if not tok:
                break
            lowTok = string.lower(tok)
            if lowTok == 'begin':
                self.readBlock(flob)
                # We may have enough for an Alignment by now.
                # If so, make the alignment and then wipe out self.nexusData
                if self.nexusData:
                    if len(self.nexusData.sequences):
                        a = self.nexusData.alignment()
                        if hasattr(flob, 'name'):
                            a.fName = flob.name
                        #a.nexusSets = NexusSets()
                        #a.nexusSets.nChar = a.length
                        #a.nexusSets.setPredefinedCharSets(a)
                        self.alignments.append(a)
                        self.nexusData = None

            # This next bit is counter to MSM97, but is commonly seen
            # and not much of a crime, so don't die.  But do complain.
            elif lowTok == '#nexus':
                print gm[0]
                print "    Skipping spurious '%s'" % tok
            else:
                gm.append("I was expecting a 'begin', to start a block.")
                gm.append("Token '%s' is not recognized." % tok)
                raise Glitch, gm

        # Thats it, the file has been read in to the end.
        # We may have enough for an Alignment by now.
        # If so, make the alignment and then wipe out self.nexusData

        # I suppose it is possible to read a Taxa block in one file
        # and a characters block in a subsequent file.  Will this do
        # it?  Will that ever happen?
        if self.nexusData:
            if len(self.nexusData.sequences):
                a = self.nexusData.alignment()
                if hasattr(flob, 'name'):
                    a.fName = flob.name
                #a.nexusSets = NexusSets()
                #a.nexusSets.nChar = a.length
                #a.nexusSets.setPredefinedCharSets(a)
                self.alignments.append(a)
                self.nexusData = None
        flob.close()
        return 1

    def readBlock(self, flob):

        if hasattr(flob, 'name'):
            gm = ['Nexus.readBlock() from file %s' % flob.name]
        else:
            gm = ['Nexus.readBlock()']
        # We start here having read the tok "begin".  So get the next
        # token to find out what kind of block.

        # We need to import nextTok.
        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok,nexusSkipPastNextSemiColon,nexusSkipPastBlockEnd
        else:
            from NexusToken import nextTok,nexusSkipPastNextSemiColon,nexusSkipPastBlockEnd

        tok = nextTok(flob)
        if not tok:
            gm.append("Failed to get block type.  Premature end of file?")
            raise Glitch, gm
        blockType = string.lower(tok)
        if var.verboseRead:
            print "    Entering '%s' block..." % tok

        if blockType in ['data', 'taxa', 'characters']:
            if blockType == 'taxa':
                if self.nexusData:  
                    # We already have a nexusData, but we have just come
                    # on another taxa block.  So zap the current
                    # nexusData.
                    #if var.verboseRead:
                    #    print "        Already have a previous nexusData.  The new one will be the current one."
                    del(self.nexusData)
                    self.nexusData = None
                                    
            if not self.nexusData:
                self.nexusData = NexusData()
            self.nexusData.readBlock(flob, blockType)
            nexusSkipPastNextSemiColon(flob)
            if var.verboseRead:
                print "    Finished '%s' block." % tok
        elif blockType == 'trees':
            self.readTreesBlock(flob)
            nexusSkipPastNextSemiColon(flob)
            if var.verboseRead:
                print "    Finished '%s' block." % tok
            # print "got trees block"
            if len(self.trees):
                if self.nexusData and self.nexusData.taxNames:
                    for t in self.trees:
                        # We will have a nexusData object even if we have only read in a taxa block.

                        # Should I dupe the taxNames??  ie t.taxNames =
                        # self.nexusData.taxNames[:]?  I don't, but not
                        # doing so might be a source of bugs if the tree
                        # taxNames ever change.                    
                        #t.taxNames = self.nexusData.taxNames[:]

                        #print 'self.nexusData.taxNames = %s' % self.nexusData.taxNames
                        t.taxNames = self.nexusData.taxNames
                        #if var.doRepairDupedTaxonNames:
                        #    t.checkDupedTaxonNames()
                        #t.checkTaxNames()
                    # I don't know why I do (did) this.  This is done by func._tryToReadNexusFile()
                    #var.trees += self.trees
                    #self.trees = []
                #else:
                #    for t in self.trees:
                #        if var.doRepairDupedTaxonNames:
                #            t.checkDupedTaxonNames()
                    
            
        elif blockType == 'sets':

            if not var.nexusSets:
                var.nexusSets = NexusSets()
            var.nexusSets._continueReadingFromNexusFile(flob)
            # When the above finishes, we have just read the 'end' of the block, but not the semi-colon


            nexusSkipPastNextSemiColon(flob)
            if var.verboseRead:
                print "    Finished '%s' block." % tok

            #self.nexusSets.dump()
        else:
            if var.nexus_warnSkipUnknownBlock:
                print "    Skipping unknown Nexus block '%s'" % blockType
            nexusSkipPastNextSemiColon(flob)
            nexusSkipPastBlockEnd(flob)


    def readTreesBlock(self, flob):
        from Tree import Tree
        if hasattr(flob, 'name'):
            gm = ['Nexus.readTreesBlock() from file %s' % flob.name]
        else:
            gm = ['Nexus.readTreesBlock()']

        # Maybe I need to save these, and restore them later??
        ##    nexus_writeVisibleComments = 0      # write, but do not necessarily get, all [!...]
        ##    nexus_getP4CommandComments = 0      # all [&&p4 ...]
        ##    nexus_getWeightCommandComments = 0  # all [&w ...]
        ##    nexus_getAllCommandComments = 0     # all [&...]   (should include all [\...] also...)

        if 0:
            print gm[0]
            print "    var.nexus_writeVisibleComments = %s" % var.nexus_writeVisibleComments
            print "    var.nexus_getP4CommandComments = %s" % var.nexus_getP4CommandComments
            print "    var.nexus_getWeightCommandComments = %s" % var.nexus_getWeightCommandComments
            print "    var.nexus_getAllCommandComments = %s" % var.nexus_getAllCommandComments

        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok,nexusSkipPastNextSemiColon,safeNextTok
        else:
            from NexusToken import nextTok,nexusSkipPastNextSemiColon,safeNextTok

        # We have read the word 'trees' in 'begin trees;', but have
        # not read the semicolon yet.  So the first thing to do is ...
        nexusSkipPastNextSemiColon(flob)
        commandName = nextTok(flob)
        if commandName:
            lowCommandName = string.lower(commandName)
        else:
            gm.append("Expecting a tree definition")
            raise Glitch, gm
        translationHash = None
        hasDoneATree = None
        doMcmcModelComments = None
        theModelInfo = None
        while lowCommandName and lowCommandName != 'end' and lowCommandName != 'endblock':
            if lowCommandName == 'translate':
                if hasDoneATree:
                    gm.append("The 'translation' command should come before any 'tree' commands, ok?")
                    raise Glitch, gm
                elif translationHash:
                    gm.append("You can't have more than one 'translation' command in a trees block, ok?")
                    raise Glitch, gm
                translationHash = self.readTranslateCommand(flob)

                if var.doTreeReadMcmcModelUsageComments:
                    # Read comments like [&&p4 models p1 c0.2 r0.1 g0.0] following the translate command
                    # We will make var.doTreeReadMcmcModelUsageComments be nParts, eg for p1, nParts=1
                    savedState = var.nexus_getP4CommandComments
                    var.nexus_getP4CommandComments = 1
                    lowTok = lowCommandName
                    while lowTok != 'tree':
                        tok = safeNextTok(flob)
                        if tok[0] == '[':
                            #print "xxyyx got comment: %s" % tok
                            from TreePartitions import _getModelInfo
                            theModelInfo = _getModelInfo(tok)
                            if theModelInfo:
                                # theModelInfo.check() returns
                                #   1 for complete info, but no heterogeneity over the tree
                                #   2 for complete info, and heterogeneity over the tree
                                ret = theModelInfo.check()
                                #print 'got ret=%s' % ret
                                if ret == 1:
                                    pass # doMcmcModelComments remains None
                                elif ret == 2:
                                    doMcmcModelComments = 1
                        else:
                            lowTok = string.lower(tok)
                    lowCommandName = lowTok
                    var.nexus_getP4CommandComments = savedState
                    if var.doTreeReadMcmcModelUsageComments and not theModelInfo:
                        gm.append("""You have var.doTreeReadMcmcModelUsageComments set, but I was
                        not able to find a model info command comment between the
                        translate command and the first tree.  Should you really have
                        var.doTreeReadMcmcModelUsageComments set?  Or is this not a proper
                        p4 mcmc trees output file with model comments?""")
                        raise Glitch, gm
                    continue


            elif lowCommandName == 'tree':
                newTree = Tree()
                if doMcmcModelComments and theModelInfo:
                    newTree.parseNexus(flob, translationHash, doModelComments=theModelInfo.nParts)
                    newTree.modelInfo = theModelInfo
                else:
                    newTree.parseNexus(flob, translationHash)
                self.trees.append(newTree)
                hasDoneATree = 1
            elif lowCommandName == 'utree':
                if var.nexus_allowUTREE:
                    newTree = Tree()
                    if doMcmcModelComments and theModelInfo:
                        newTree.parseNexus(flob, translationHash, doModelComments=theModelInfo.nParts)
                        newTree.modelInfo = theModelInfo
                    else:
                        newTree.parseNexus(flob, translationHash)
                    self.trees.append(newTree)
                    hasDoneATree = 1
                else:
                    gm.append("Use of '%s' is deprecated." % commandName)
                    gm.append("Use 'tree treeName = [&U] ...' instead,")
                    gm.append("    (although p4 does not care if you use the [&U] or not.)")
                    gm.append("(You can force reading of UTREE commands with var.nexus_allowUTREE.)")
                    raise Glitch, gm
            elif lowCommandName[0] not in string.lowercase:
                gm.append("Got spurious '%s' (...exiting)" % commandName)
                raise Glitch, gm
            else:
                if var.verboseRead:
                    print "        skipping unknown trees block command '%s'" % commandName
                nexusSkipPastNextSemiColon(flob)
            commandName = nextTok(flob)
            if commandName:
                lowCommandName = string.lower(commandName)
            else:
                gm.append("Expecting a trees block command (eg 'tree'), or 'end'")
                raise Glitch, gm


    def readTranslateCommand(self, flob):

        if hasattr(flob, 'name'):
            gm = ['Nexus.readTranslateCommand() (in a trees block) from file %s' % flob.name]
        else:
            gm = ['Nexus.readTranslateCommand() (in a trees block)']

        ##        We have two references: Maddison et al, Syst Biol. 46:590-621,
        ##        1997., and the PAUP 4.0b8 manual, in pdf.

        ##        Maddison et al says, on page 613:

        ##        The syntax for the TREEES block is
        ##        BEGIN TREES;
        ##            [TRANSLATE arbitrary-token-used-in-tree-description
        ##             valid-taxon-name
        ##             [, arbitrary-token-used-in-tree-description
        ##             valid-taxon-name...];]
        ##            [TREE [*] tree-name=tree-specification;]
        ##        END;

        ##        The PAUP manual says this:

        ##        The syntax for the TREES block follows:
        ##        BEGIN TREES [ block-name ] ;
        ##            [ TRANSLATE token taxon-name [ , token taxon-name  ]
        ##               ... ; ]
        ##            [ TREE [*] name = tree-specification; ]
        ##        END;

        translationHash = {}

        # We need to import safeNextTok.
        if var.nexus_doFastNextTok:
            from NexusToken2 import safeNextTok
        else:
            from NexusToken import safeNextTok


        while 1:
            keyTok = func.nexusUnquoteName(safeNextTok(flob, 'Nexus: readTranslateCommand'))
            #print "x got keyTok '%s'" % keyTok
            if keyTok == None or keyTok == ';':
                break
            valueTok = safeNextTok(flob, 'Nexus: readTranslateCommand')
            valueTok = func.nexusUnquoteName(valueTok)
            #print "  got valueTok '%s'" % valueTok

            if valueTok == None or valueTok == ';':
                gm.append("Translate items should be pairs.")
                gm.append("Got nothing to pair with '%s'" % keyTok)
                raise Glitch, gm

            # It should be that the keyTok is the integer (usually,
            # but not necessarily) and the valueTok is the
            # valid-taxon-name.  For older MrBayes (before v3) it was
            # backwards, and PAUP could handle it.  But here I will
            # assume that it is not backwards.

            # The valueTok should not be all digits, unless
            # var.nexus_allowAllDigitNames is turned on.
            if not var.nexus_allowAllDigitNames:
                ret = min([c in string.digits for c in valueTok])
                if ret == True:  # meaning that valueTok is all digits
                    gm.append("Got all-digit name '%s'" % valueTok)
                    gm.append("But var.nexus_allowAllDigitNames is currently False.")
                    gm.append("If you want to allow this, set it to True")
                    raise Glitch, gm
                
            translationHash[keyTok] = valueTok
            commaTok = safeNextTok(flob, 'Nexus: readTranslateCommand')
            #print "  got commaTok '%s'" % commaTok
            if commaTok == ',':
                pass
            elif commaTok == ';':
                break
            else:
                gm.append("Expecting a comma or semi-colon, but got '%s'" % commaTok)
                gm.append("Just after '%s'" % valueTok)
                raise Glitch, gm

        return translationHash


class NexusData:
    def __init__(self):
        self.nTax = None
        self.nChar = None
        self.dataType = 'standard'
        self.formatCommandSymbols = None
        self.symbols = None
        self.nexus_gap = None #'-'
        self.nexus_missing = None # '?'
        self.nexus_matchchar = None
        self.formatCommandEquates = {}
        self.equates = {}
        self.interleave = 0
        self.taxNames = []
        self.lowTaxNames = []
        self.sequences = []

    def readBlock(self, flob, lowerBlockType):
        if hasattr(flob, 'name'):
            gm = ['NexusData.readBlock() from file %s' % flob.name]
        else:
            gm = ['NexusData.readBlock()']
        if lowerBlockType == 'taxa':
            self.readTaxaBlock(flob, lowerBlockType)
            if len(self.taxNames) !=self.nTax:
                gm.append("The number of taxNames, %i, must be equal to nTax, %i." % (len(self.taxNames), self.nTax))
                raise Glitch, gm
        elif lowerBlockType == 'data' or lowerBlockType == 'characters':
            if lowerBlockType == 'characters':
                if len(self.taxNames) == 0 or not self.nTax:
                    gm.append("Reading characters block")
                    gm.append("ntax and taxNames must be defined before reading characters")
                    raise Glitch, gm
            self.readDataBlock(flob, lowerBlockType)
            if len(self.sequences):


                if self.dataType == 'dna':
                    if self.formatCommandSymbols:
                        if var.nexus_ignoreFormatCommandSymbols == True:
                            print gm[0]
                            print "Ignoring extra symbols '%s' from the format command." % self.formatCommandSymbols
                        else:
                            gm.append("Got extra symbols '%s' from the format command." % self.formatCommandSymbols)
                            raise Glitch, gm
                    self.symbols = 'acgt'
                    self.equates = {}
                    self.equates['r'] = 'ag'
                    self.equates['y'] = 'ct'
                    self.equates['m'] = 'ac'
                    self.equates['k'] = 'gt'
                    self.equates['s'] = 'cg'
                    self.equates['w'] = 'at'
                    self.equates['h'] = 'act'
                    self.equates['b'] = 'cgt'
                    self.equates['v'] = 'acg'
                    self.equates['d'] = 'agt'
                    self.equates['n'] = 'acgt'
                    #self.equates['x'] = 'acgt'
                    if self.formatCommandEquates:
                        for k,v in self.formatCommandEquates.items():
                            if k in self.symbols:
                                gm.append("Equate key '%s' is one of the character symbols.  Bad." % k)
                                raise Glitch, gm
                            if self.equates.has_key(k):
                                gm.append("Equates are currently %s" % self.equates)
                                gm.append("Equate key '%s' from the format command is already in the equates." % k)
                                raise Glitch, gm
                            for c in v:
                                if c not in self.symbols:
                                    gm.append("Equate %s:%s" % (k,v))
                                    gm.append("%s is not in symbols." % c)
                                    raise Glitch(gm, 'nexus_equateIsNotInSymbols')
                            self.equates[k] = v
                elif self.dataType == 'protein':
                    if self.formatCommandSymbols:
                        if var.nexus_ignoreFormatCommandSymbols == True:
                            print gm[0]
                            print "Ignoring extra symbols '%s' from the format command." % self.formatCommandSymbols
                        else:
                            gm.append("Got extra symbols '%s' from the format command." % self.formatCommandSymbols)
                            raise Glitch, gm
                    self.symbols = 'arndcqeghilkmfpstwyv'
                    self.equates = {}
                    self.equates['b'] = 'dn'
                    self.equates['z'] = 'eq'
                    self.equates['x'] = 'arndcqeghilkmfpstwyv'
                    if self.formatCommandEquates:
                        for k,v in self.formatCommandEquates.items():
                            if k in self.symbols:
                                gm.append("Equate key '%s' is one of the character symbols.  Bad." % k)
                                raise Glitch, gm
                            if self.equates.has_key(k):
                                gm.append("Equates are currently %s" % self.equates)
                                gm.append("Equate key '%s' from the format command is already in the equates." % k)
                                raise Glitch, gm
                            for c in v:
                                if c not in self.symbols:
                                    gm.append("Equate %s:%s" % (k,v))
                                    gm.append("%s is not in symbols." % c)
                                    raise Glitch(gm, 'nexus_equateIsNotInSymbols')
                            self.equates[k] = v
                elif self.dataType == 'rna':
                    if self.formatCommandSymbols:
                        if var.nexus_ignoreFormatCommandSymbols == True:
                            print gm[0]
                            print "Ignoring extra symbols '%s' from the format command." % self.formatCommandSymbols
                        else:
                            gm.append("Got extra symbols '%s' from the format command." % self.formatCommandSymbols)
                            raise Glitch, gm
                    self.symbols = 'acgu'
                    self.equates = {}
                    self.equates['r'] = 'ag'
                    self.equates['y'] = 'cu'
                    self.equates['m'] = 'ac'
                    self.equates['k'] = 'gu'
                    self.equates['s'] = 'cg'
                    self.equates['w'] = 'au'
                    self.equates['h'] = 'acu'
                    self.equates['b'] = 'cgu'
                    self.equates['v'] = 'acg'
                    self.equates['d'] = 'agu'
                    self.equates['n'] = 'acgu'
                    #self.equates['x'] = 'acgu'
                    if self.formatCommandEquates:
                        for k,v in self.formatCommandEquates.items():
                            if k in self.symbols:
                                gm.append("Equate key '%s' is one of the character symbols.  Bad." % k)
                                raise Glitch, gm
                            if self.equates.has_key(k):
                                gm.append("Equates are currently %s" % self.equates)
                                gm.append("Equate key '%s' from the format command is already in the equates." % k)
                                raise Glitch, gm
                            for c in v:
                                if c not in self.symbols:
                                    gm.append("Equate %s:%s" % (k,v))
                                    gm.append("%s is not in symbols." % c)
                                    raise Glitch(gm, 'nexus_equateIsNotInSymbols')
                            self.equates[k] = v
                elif self.dataType == 'standard':
                    if self.formatCommandSymbols:
                        self.symbols = self.formatCommandSymbols
                    else:
                        self.symbols = '01'
                    if self.formatCommandEquates:
                        for k,v in self.formatCommandEquates.items():
                            if k in self.symbols:
                                gm.append("Equate key '%s' is one of the character symbols.  Bad." % k)
                                raise Glitch, gm
                            for c in v:
                                if c not in self.symbols:
                                    gm.append("Equate %s:%s" % (k,v))
                                    gm.append("%s is not in symbols." % c)
                                    raise Glitch, gm
                        self.equates = self.formatCommandEquates

                eKeys = ''.join(self.equates.keys())
                charsUsedSoFar = self.symbols + eKeys
                if self.nexus_gap:
                    charsUsedSoFar += self.nexus_gap
                if self.nexus_missing:
                    charsUsedSoFar += self.nexus_missing
                if self.nexus_matchchar:
                    charsUsedSoFar += self.nexus_matchchar
                    
                if not self.nexus_missing and '?' not in charsUsedSoFar:
                    self.nexus_missing = '?'
                    charsUsedSoFar += '?'

                self.propagateMatchchars()

                if self.dataType == 'dna' and self.symbols != 'acgt':
                    gm.append("Got dna dataType, but symbols were re-defined to '%s'" % self.symbols)
                    raise Glitch, gm
                elif self.dataType == 'rna' and self.symbols != 'acgu':
                    gm.append("Got rna dataType, but symbols were re-defined to '%s'" % self.symbols)
                    raise Glitch, gm
                elif self.dataType == 'protein' and self.symbols != 'arndcqeghilkmfpstwyv':
                    gm.append("Got protein dataType, but symbols were re-defined to '%s'" % self.symbols)
                    raise Glitch, gm

                # In some pathological treebase aligns, the gapchar is
                # a letter N.  We can't allow that sort of nonsense.
                eKeys = ''.join(self.equates.keys())
                if self.nexus_gap:
                    if self.nexus_gap in self.symbols:
                        gm.append("The gap symbol '%s' is one of the character state symbols '%s'.  Bad." % (
                            self.nexus_gap, self.symbols))
                        raise Glitch, gm
                    if self.nexus_gap in eKeys:
                        gm.append("The gap symbol '%s' is one of the equate symbols '%s'.  Bad." % (
                            self.nexus_gap, self.eKeys))
                        raise Glitch(gm, 'nexus_badSymbolForGap')
                        #del(self.equates[self.nexus_gap])
                        #eKeys = ''.join(self.equates.keys())

                if self.nexus_missing:
                    if self.nexus_missing in self.symbols:
                        gm.append("The symbol for missing, '%s', is one of the character state symbols '%s'.  Bad." % (
                            self.nexus_missing, self.symbols))
                        raise Glitch, gm
                    if self.nexus_missing in eKeys:
                        gm.append("The symbol for missing, '%s', is one of the equate symbols '%s'.  Bad." % (
                            self.nexus_missing, eKeys))
                        raise Glitch(gm, 'nexus_badSymbolForMisssing')
                        #del(self.equates[self.nexus_missing])
                        #eKeys = ''.join(self.equates.keys())

                if self.nexus_gap and self.nexus_missing:
                    if self.nexus_gap == self.nexus_missing:
                        gm.append("The gap char and the missing char are the same.")
                        raise Glitch(gm, 'nexus_badGapMissing')
                        
                # put the sequences in lowercase
                if self.nTax:
                    for i in range(self.nTax):
                        self.sequences[i] = string.lower(self.sequences[i])

                
                # Ask whether any of the charsUsedSoFar is represented more than once in that string.
                for symb in charsUsedSoFar:
                    if charsUsedSoFar.count(symb) != 1:
                        gm.append("current valid symbols include '%s'" % charsUsedSoFar) 
                        gm.append('Symbol %s is used more than once.' % symb)
                        raise Glitch, gm
                
                #print "dataType=%s, charsUsedSoFar=%s" % (self.dataType, charsUsedSoFar)
                valids = charsUsedSoFar
                if self.dataType == 'protein':
                    valids += '*'
                i = 0
                bads = 0
                if self.nTax:
                    for i in range(self.nTax):
                        j = 0
                        while j < self.nChar:
                            if self.sequences[i][j] not in valids:
                                theComplaint = "Bad character '%s' in (zero-based) sequence " % self.sequences[i][j]
                                theComplaint += "%s at (zero-based) position %s" % (i, j)
                                gm.append(theComplaint)
                                bads += 1
                                if bads > 10:
                                    gm.append("...and possibly others")
                                    break
                            j = j + 1
                        if bads > 10:
                            break
                    if bads:
                        raise Glitch(gm, 'nexus_badCharacter')
            else:
                gm.append("NexusData  read %s block" % lowerBlockType)
                gm.append("The point of reading this kind of block is to")
                gm.append("read in some sequences, but none were found.")
                raise Glitch, gm



    def readTaxaBlock(self, flob, blockType):
        # See page 597 in MSM97.
        # There are only two commands possible-- dimensions and
        # taxlabels.  Dimensions must appear before taxlabels.  Only
        # one of each is allowed per block.

        # Correction: Although page 597 in MSM97 implies that there
        # are only 2 commands possible, but on page 595 it says that
        # new commands in public blocks can be added.

        if hasattr(flob, 'name'):
            gm = ['NexusData.readTaxaBlock() from file %s' % flob.name]
        else:
            gm = ['NexusData.readTaxaBlock()']

        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok,nexusSkipPastNextSemiColon
        else:
            from NexusToken import nextTok,nexusSkipPastNextSemiColon

        nexusSkipPastNextSemiColon(flob)   # to get to the end of 'begin taxa maybe with other stuff;'
        commandName = nextTok(flob)
        if commandName:
            lowCommandName = string.lower(commandName)
        else:
            gm.append("Failed to read any commands.")
            raise Glitch, gm

        # it is required that the dimensions command come first, then the taxlabels command.
        hasDoneDimensionsCommand = False
        hasDoneTaxlabelsCommand = False
        while lowCommandName != None and lowCommandName != 'end' and lowCommandName != 'endblock':
            if lowCommandName == 'dimensions':
                if hasDoneDimensionsCommand:
                    gm.append("Found a second dimensions command.")
                    raise Glitch, gm
                else:
                    self.readDimensionsCommand(flob, blockType)
                    if not self.nTax:
                        gm.append("ntax must be defined in the taxa block.")
                        raise Glitch, gm
                    hasDoneDimensionsCommand = True
            elif lowCommandName == 'taxlabels':
                if not hasDoneDimensionsCommand:
                    gm.append("The 'dimensions' command must precede the 'taxlabels' command.")
                    raise Glitch, gm
                self.readTaxlabelsCommand(flob)
                hasDoneTaxlabelsCommand = True
                # print "got taxa %s" % self.taxNames

            else:
                if var.verboseRead:
                    print "        skipping unknown taxa block command '%s'" % commandName
                nexusSkipPastNextSemiColon(flob)

            commandName = nextTok(flob)
            if commandName:
                lowCommandName = string.lower(commandName)
            else:
                gm.append("Premature end of file.")
                raise Glitch, gm
        if not hasDoneDimensionsCommand:
            gm.append("Failed to find a 'dimensions' command.")
            raise Glitch, gm
        if not hasDoneTaxlabelsCommand:
            gm.append("Failed to find a 'taxlabels' command.")
            raise Glitch, gm


    def readDataBlock(self, flob, blockType):  # either 'data' or 'characters' blocks
        if hasattr(flob, 'name'):
            gm = ['NexusData.readDataBlock() from file %s' % flob.name]
        else:
            gm = ['NexusData.readDataBlock()']
        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok,nexusSkipPastNextSemiColon
        else:
            from NexusToken import nextTok,nexusSkipPastNextSemiColon

        nexusSkipPastNextSemiColon(flob)
        notImplemented = ['eliminate', 'charstatelabels', 'charlabels', 'statelabels', 'taxlabels']
        # implemented: dimensions, format,  matrix, and thats all.
        commandName = nextTok(flob)
        if commandName:
            lowCommandName = string.lower(commandName)
        else:
            return
        hasDoneDimensionsCommand = None
        hasDoneFormatCommand = None
        while lowCommandName != None and lowCommandName != 'end' and lowCommandName != 'endblock':
            # print "got command '%s', %s" % (commandName, lowCommandName)
            if lowCommandName == 'dimensions':
                if hasDoneDimensionsCommand:
                    gm.append("NexusData: read %s block" % blockType)
                    gm.append("Attempting a second dimensions command.  Not allowed.")
                    raise Glitch, gm
                self.readDimensionsCommand(flob, blockType)
                hasDoneDimensionsCommand = 1
            elif lowCommandName == 'format':
                if hasDoneFormatCommand:
                    gm.append("NexusData: read %s block" % blockType)
                    gm.append("Attempting a second format command.  Not allowed.")
                    raise Glitch, gm
                self.readFormatCommand(flob, blockType)
                hasDoneFormatCommand = 1
            elif lowCommandName == 'matrix':
                if not self.nTax:
                    gm.append("Read %s block" % lowCommandName)
                    gm.append("ntax must be defined before reading characters")
                    raise Glitch, gm
                if not self.nChar:
                    gm.append("Read %s block" % lowCommandName)
                    gm.append("nChar must be defined before reading characters")
                    raise Glitch, gm
                if self.interleave:
                    self.readInterleaveMatrix(flob, blockType)
                else:
                    self.readNonInterleaveMatrix(flob, blockType)
            elif lowCommandName in notImplemented:
                gm.append("%s block command '%s' is not implemented" %  (blockType, commandName))
                raise Glitch(gm, 'nexus_commandNotImplemented')
            elif lowCommandName[0] not in string.lowercase:
                gm.append("Got spurious '%s' (...exiting)" % commandName)
                raise Glitch, gm
            else:
                if var.verboseRead:
                    print "        skipping unknown %s block command '%s'" % (blockType, commandName)
                nexusSkipPastNextSemiColon(flob)
            commandName = nextTok(flob)
            if commandName:
                lowCommandName = string.lower(commandName)
            else:
                return


    def parseSubcommandEqualsArg(self, flob, blockType, command, sub):
        """SomeCommand someSubCommand=someArg"""
        
        if hasattr(flob, 'name'):
            gm = ['NexusData.parseSubcommandEqualsArg(blockType=%s,command=%s, sub=%s) from file %s' % (blockType, command, sub, flob.name)]
        else:
            gm = ['NexusData.parseSubcommandEqualsArg(blockType=%s,command=%s, sub=%s)' % (blockType, command, sub)]
        
        #print "parseSubcommandEqualsArg sub=%s, command=%s" % (sub,command)

        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok,safeNextTok
        else:
            from NexusToken import nextTok,safeNextTok

        tok = nextTok(flob)
        if tok:
            if tok == '=':
                tok = nextTok(flob)
                #print "parseSubcommandEqualsArg, sub=%s, got token '%s'" % (sub, tok)
                if tok:
                    if tok == ';':
                        gm.append("%s block %s subcommand '%s': premature command end." %  (blockType, command, sub))
                        raise Glitch, gm
                    elif sub == 'nchar':
                        try:
                            self.nChar = int(tok)
                            # print "got nChar = %i" % self.nChar
                        except ValueError:
                            gm.append("Bad arg for %s block %s subcommand '%s'"  % (blockType, command, sub))
                            raise Glitch, gm
                    elif sub == 'ntax':
                        try:
                            self.nTax = int(tok)
                            # print "got nTax = %i" % self.nTax
                        except ValueError:
                            gm.append("Bad arg for %s block %s subcommand '%s'"  % (blockType, command, sub))
                            raise Glitch, gm
                    elif sub == 'datatype':
                        lowTok = string.lower(tok)
                        assumedToBeDNA = False
                        if lowTok == 'nucleotide':
                            # assume it is dna
                            assumedToBeDNA = True
                            print "Got datatype '%s', assuming that it is DNA." % tok
                        if assumedToBeDNA or lowTok == 'dna':
                            self.dataType = 'dna'

                        elif lowTok == 'protein':
                            self.dataType = 'protein'
                            
                        elif lowTok == 'standard':
                            self.dataType = 'standard'
                            
                        elif lowTok == 'rna':
                            self.dataType = 'rna'

                        else:
                            gm.append("%s block format subcommand '%s': " % (blockType, sub))
                            gm.append("I only do dna, protein, and standard.")
                            gm.append("I cannot do '%s'." % tok)
                            raise Glitch, gm

                    elif sub == 'gap':
                        lowTok = string.lower(tok)
                        self.nexus_gap = lowTok
                        #if lowTok != var.nexus_gap:
                            #gm.append("%s block format subcommand '%s': " % (blockType, sub))
                            #gm.append("The '-' character is the only one allowed for gaps.")
                            #gm.append("Deal with it.  (exiting...)")
                            #raise Glitch, gm
                    elif sub == 'missing':
                        lowTok = string.lower(tok)
                        self.nexus_missing = lowTok
                        #if lowTok != '?':
                        #    gm.append("%s block format subcommand '%s': " % (blockType, sub))
                        #    gm.append("The '?' character is the only one allowed for missing.")
                        #    gm.append("Deal with it.  (exiting...)")
                        #    raise Glitch, gm
                    elif sub == 'matchchar':
                        lowTok = string.lower(tok)
                        self.nexus_matchchar = lowTok
                        #if lowTok != '.':
                        #    gm.append("%s block 'format' subcommand '%s': " % (blockType, sub))
                        #    gm.append("The '.' character is the only one allowed for matchchar.")
                        #    gm.append("Deal with it.  (exiting...)")
                        #    raise Glitch, gm
                    elif sub == 'equate':
                        #if not self.symbols:
                        #    gm.append("NexusData, 'format' subcommand 'equate'")
                        #    gm.append("should be preceded by dataType, and symbols if needed.")
                        #    raise Glitch, gm
                        if tok == '\"':
                            # print "equate starts well, at least"
                            tok = nextTok(flob)
                            a = None
                            eq = None
                            b = None
                            while tok != '\"':
                                if not a:
                                    a = tok
                                    if len(a) != 1:
                                        gm.append("Bad equate arg %s.  Single characters only, please" % tok)
                                        raise Glitch, gm
                                    #if a in self.symbols:
                                    #    gm.append("Nexus format command. Bad equate key %s" % a)
                                    #    gm.append("should not be one of the symbols")
                                elif not eq:
                                    if tok != '=':
                                        gm.append("Bad equate arg near '%s'.  No '='?" % tok)
                                        raise Glitch, gm
                                    eq = 1
                                elif not b:
                                    b = tok
                                    #print "Got equate third token '%s', length %i" % (tok, len(tok))
                                    if b == '{':
                                        bList = []
                                        tok = safeNextTok(flob, 'format subcommand equate')
                                        #print "got tok '%s'" % tok
                                        while tok != '}':
                                            bList.append(string.lower(tok))
                                            tok = safeNextTok(flob, 'format subcommand equate')
                                            #print "got tok '%s'" % tok
                                        b = string.join(bList, '')
                                    elif b == '(':
                                        bList = []
                                        tok = safeNextTok(flob, 'format subcommand equate')
                                        #print "got tok '%s'" % tok
                                        while tok != ')':
                                            bList.append(string.lower(tok))
                                            tok = safeNextTok(flob, 'format subcommand equate')
                                            #print "got tok '%s'" % tok
                                        b = string.join(bList, '')
                                    #print "Got b %s" % b
                                    self.formatCommandEquates[a] = b.lower()
                                    #print "got equate %s = %s" % (a, b)
                                    a = None
                                    eq = None
                                    b = None
                                tok = nextTok(flob)
                            if a or eq or b:
                                gm.append("Bad equate arg. Leftovers?")
                                raise Glitch, gm
                        else:
                            gm.append("%s block format subcommand '%s': must be bounded by double quotes." % (blockType, sub))
                            raise Glitch, gm
                    elif sub == 'symbols':
                        if tok == '\"':
                            symbolsList = []
                            tok = nextTok(flob)
                            while tok != '\"':
                                symbolsList.append(string.lower(tok))
                                tok = safeNextTok(flob, 'format subcommand symbols')
                            newSymbols = string.join(symbolsList, '')
                            #if len(newSymbols) > 20:
                            #    gm.append("Got %i new symbols.  Rather a lot, yes?" % len(newSymbols))
                            #    gm.append("Got new symbols = %s" % newSymbols)
                            #    raise Glitch, gm

                            # See notes on symbols, at the top of this
                            # file, in the MSM97 notes.  If its not
                            # standard datatype, I should *add* the
                            # newly-defined symbols, like this:
                            #if self.datatype in ['dna', 'protein']:
                            #    self.symbols += newSymbols
                            # Except that I don't like that.  This could be considered a bug.  I deal with this later.
                            self.formatCommandSymbols = newSymbols
                            #print "wxyz Got symbols: %s" % self.symbols
                        else:
                            gm.append("%s block format subcommand '%s': arg must be bounded by double quotes." % (blockType, sub))
                            raise Glitch, gm
                else:
                    gm.append("%s block %s subcommand '%s' needs an arg: premature end" % (blockType, command, sub))
                    raise Glitch, gm
            else:
                gm.append("%s block %s subcommand '%s' needs an arg: no '='" % (blockType, command, sub))
                raise Glitch, gm
        else:
            gm.append("%s block %s subcommand '%s' needs an arg: did not even get a token." % (blockType, command, sub))
            raise Glitch, gm

    def readTaxlabelsCommand(self, flob):
        if hasattr(flob, 'name'):
            gm = ['NexusData.readTaxlabelsCommand() from file %s' % flob.name]
        else:
            gm = ['NexusData.readTaxlabelsCommand()']
        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok
        else:
            from NexusToken import nextTok

        tok = nextTok(flob)
        while tok and tok != ';':
            theName = func.nexusUnquoteName(tok)
            if not func.nexusCheckName(theName):
                gm.append("Bad nexus name '%s'" % theName)
                raise Glitch(gm, 'nexus_badName')
            lowName = string.lower(theName)
            if lowName in self.lowTaxNames:
                gm.append("Got duplicated (lowercased) taxname %s" % theName)
                raise Glitch(gm, 'nexus_duplicatedTaxnames')
            self.lowTaxNames.append(lowName)
            self.taxNames.append(theName)
            tok = nextTok(flob)
        if len(self.taxNames) != self.nTax:
            gm.append("The number of taxon names in the 'taxlabels' command (%s)" % len(self.taxNames))
            gm.append("is not the same as ntax (%s)" % self.nTax)
            raise Glitch(gm, 'nexus_nTaxProblem')

    def readDimensionsCommand(self, flob, blockType):
        if hasattr(flob, 'name'):
            gm = ['NexusData.readDimensionsCommand() from file %s' % flob.name]
        else:
            gm = ['NexusData.readDimensionsCommand()']
        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok
        else:
            from NexusToken import nextTok

        sub = nextTok(flob)  # get subcommand
        if sub:
            lowSub = string.lower(sub)
        else:
            return
        while lowSub != None and lowSub != ';':
            # print "got subCommand '%s', %s" % (sub, lowSub)
            # print "blockType = %s" % blockType
            if blockType == 'taxa' and lowSub == 'ntax':
                self.parseSubcommandEqualsArg(flob, blockType, 'dimensions', lowSub)
            elif blockType == 'characters' and lowSub == 'nchar':
                self.parseSubcommandEqualsArg(flob, blockType, 'dimensions', lowSub)
            elif blockType == 'data' and lowSub in ['nchar', 'ntax']:
                self.parseSubcommandEqualsArg(flob, blockType, 'dimensions', lowSub)
            elif lowSub == '=':
                gm.append("%s block dimensions command: spurious '%s'" % (blockType, sub))
                raise Glitch, gm
            elif lowSub == 'newtaxa':
                gm.append("%s block dimensions subcommand '%s' not implemented" % (blockType, sub))
                raise Glitch, gm
            else:
                print "skipping unknown %s block dimensions subcommand '%s'" % (blockType, sub)
                sub = nextTok(flob)
                if sub:
                    # print "after unknown: got '%s'" % sub
                    if sub == '=':
                        sub = nextTok(flob)
                        if sub:
                            # print "after equals after unknown: got '%s'" % sub
                            if sub == ';':
                                gm.append("%s block dimensions command: spurious '='" % blockType)
                                raise Glitch, gm
                            sub = nextTok(flob)
                            if sub:
                                lowSub = string.lower(sub)
                                continue
                            else:
                                return
                        else:
                            gm.append("%s block dimensions command: spurious '='" % blockType)
                            raise Glitch, gm
                    else:
                        lowSub = string.lower(sub)
                        continue
                else:
                    return

            sub = nextTok(flob)
            if sub:
                lowSub = string.lower(sub)
                continue
            else:
                return

    def readFormatCommand(self, flob, blockType):
        if hasattr(flob, 'name'):
            gm = ['NexusData.readFormatCommand() from file %s' % flob.name]
        else:
            gm = ['NexusData.readFormatCommand()']
        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok
        else:
            from NexusToken import nextTok

        notImplemented = ['respectcase',
                        'labels', 'nolabels', 'transpose', 'items', 'statesformat',
                        'tokens', 'notokens']
        sub = nextTok(flob)
        if sub:
            lowSub = string.lower(sub)
        else:
            return
        lastSubcommand = None
        while lowSub != None and lowSub != ';':
            # print "got subCommand '%s', %s" % (sub, lowSub)
            #print "q1 got subCommand '%s', lowSub=%s, self.symbols=%s" % (sub, lowSub,self.symbols)
            if lowSub in ['datatype', 'equate', 'gap', 'missing', 'matchchar', 'symbols']:
                self.parseSubcommandEqualsArg(flob, blockType, 'format', lowSub)
                lastSubcommand = lowSub
            elif lowSub == 'interleave':
                self.interleave = 1
                lastSubcommand = lowSub
            elif lowSub == '=':
                gm.append("%s block format command: spurious '%s'" % (blockType, sub))
                if lastSubcommand == 'interleave':
                    gm.append("(Note that the Nexus standard does not allow eg 'interleave=yes' -- it is just 'interleave')")
                #elif lastSubcommand == 'respectcase':
                #    gm.append("(Note that the Nexus standard does not allow eg 'respectcase=yes' -- it is just 'respectcase')")
                raise Glitch, gm
            elif lowSub in notImplemented:
                #lastSubcommand = lowSub
                gm.append("%s block format subcommand '%s' not implemented" % (blockType, sub))
                raise Glitch, gm
            else:
                print "skipping unknown %s block format subcommand '%s'" % (blockType, sub)
                #raise Glitch
                sub = nextTok(flob)
                if sub:
                    # print "after unknown: got '%s'" % sub
                    if sub == '=':
                        sub = nextTok(flob)
                        if sub:
                            # print "after equals after unknown: got '%s'" % sub
                            if sub == ';':
                                gm.append("%s block format command: spurious '='" % blockType)
                                raise Glitch, gm
                            sub = nextTok(flob)
                            if sub:
                                lowSub = string.lower(sub)
                                continue
                            else:
                                return
                        else:
                            gm.append("%s block format command: spurious '='" % blockType)
                            raise Glitch, gm
                    else:
                        lowSub = string.lower(sub)
                        continue
                else:
                    return

            sub = nextTok(flob)
            if sub:
                lowSub = string.lower(sub)
                continue
            else:
                print "m2 self.symbols=%s" % self.symbols
                return

    def readInterleaveMatrix(self, flob, blockType):
        if hasattr(flob, 'name'):
            gm = ['NexusData.readInterleaveMatrix() from file %s' % flob.name]
        else:
            gm = ['NexusData.readInterleaveMatrix()']
        dbug = 0
        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok
        else:
            from NexusToken import nextTok

        if dbug:
            print "Nexus: readInterleaveMatrix here, blockType = %s" % blockType
        for i in range(self.nTax):
            self.sequences.append([])
        var.nexus_getLineEndingsAsTokens = 1
        tok = func.nexusUnquoteName(nextTok(flob))
        #if dbug:
        #    if tok == '\n' or tok == '\r':
        #        print '%10s: %s' % ('firsttok', '\\n')
        #    else:
        #        print '%10s: %s' % ('firsttok', tok)
        tokens = []
        counter = 0
        while tok != None and tok != ';':
            if dbug:
                if tok == '\n' or tok == '\r':
                    print '%10s: %s' % ('tok', 'line ending')
                else:
                    print '%10s: %s' % ('tok', tok)
            if (tok == '\n' or tok == '\r') and len(tokens):
                sn = counter % self.nTax
                if dbug:
                    print "        got one line for sequence number %i" % sn
                if blockType == 'data':
                    if len(self.taxNames) < self.nTax:
                        self.taxNames.append(tokens[0])
                    elif len(self.taxNames) >= self.nTax:
                        if self.taxNames[sn] != tokens[0]:
                            gm.append("interleaved matrix name mismatch: %s does not match %s" % \
                                     (self.taxNames[sn], tokens[0]))
                            raise Glitch, gm
                elif blockType == 'characters':
                    if self.taxNames[sn] != tokens[0]:
                        gm.append("Name mismatch: %s does not match (previous) %s" % \
                                 (self.taxNames[sn], tokens[0]))
                        raise Glitch, gm
                self.sequences[sn].append(string.join(tokens[1:], ''))
                # print "%s: %s" % (tokens[0], self.sequences[sn])
                tokens = []
                counter = counter + 1
            elif tok == '\n' or tok == '\r':  #  and len(tokens) == 0
                pass
            else:
                if len(tokens) == 0:
                    # Its the first token, which would be the tax name.  Check it.
                    if not func.nexusCheckName(tok):
                        gm.append("Problem with nexus name '%s': it does not appear to be nexus-compliant." % tok)
                        raise Glitch(gm, 'nexus_badName')
                tokens.append(tok)
                #if dbug:
                #    print "       tokensLength is now %i" % len(tokens)
            tok = func.nexusUnquoteName(nextTok(flob))
        for i in range(self.nTax):
            self.sequences[i] = string.join(self.sequences[i], '')
        var.nexus_getLineEndingsAsTokens = 0 # back to normal

    def readNonInterleaveMatrix(self, flob, blockType):
        if hasattr(flob, 'name'):
            gm = ['NexusData.readNonInterleaveMatrix() from file %s' % flob.name]
        else:
            gm = ['NexusData.readNonInterleaveMatrix()']
        if var.nexus_doFastNextTok:
            from NexusToken2 import nextTok
        else:
            from NexusToken import nextTok

        tok = func.nexusUnquoteName(nextTok(flob))
        tokens = []
        tokensLen = 0
        counter = 0
        while tok != None and tok != ';':
            #print "readNonInterleaveMatrix().  Got tok %s" % tok
            if not func.nexusCheckName(tok):
                gm.append("Problem with nexus name '%s': it does not appear to be nexus-compliant." % tok)
                raise Glitch(gm, 'nexus_badName')
            # get taxname
            if blockType == 'data':
                self.taxNames.append(tok)
            elif blockType == 'characters':
                if self.taxNames[counter]:
                    if self.taxNames[counter] != tok:
                        gm.append("Name mismatch in characters block matrix:")
                        gm.append("\t%s does not match %s" % (self.taxNames[counter], tok))
                        raise Glitch, gm
            # print "got taxname %s" % tok
            # get sequence
            tok = func.nexusUnquoteName(nextTok(flob))
            # print "got token: '%s'" % tok
            while tok != ';':
                tokensLen += len(tok)
                tokens.append(tok)
                if tokensLen == self.nChar:
                    self.sequences.append(string.join(tokens, ''))
                    tokens = []
                    tokensLen = 0
                    counter = counter + 1
                    tok = func.nexusUnquoteName(nextTok(flob))
                    break
                elif tokensLen > self.nChar:
                    gm.append("Sequence for taxon %s appears to be too long" % self.taxNames[counter])
                    gm.append("%s" % tokens)
                    raise Glitch(gm, 'nexus_badSequenceLength')
                else:
                    tok = func.nexusUnquoteName(nextTok(flob))
                    # print "got token: '%s'" % tok
                    if not tok:
                        gm.append("End of file reached while reading data.")
                        raise Glitch, gm
        if not counter:
            gm.append("No sequences?")
            raise Glitch, gm

        elif counter != self.nTax:
            gm.append("The number of sequences (%i) does not match ntax (%i)" % (counter, self.nTax))
            raise Glitch, gm
        #for i in range(self.nTax):
        #    print "%s: " % self.taxNames[i],
        #    print "%s" % self.sequences[i]

    def propagateMatchchars(self):
        if not self.nexus_matchchar:
            return
        gm = ['NexusData.propagateMatchchars()']
        # are there any matchchars?
        if string.count(self.sequences[0], self.nexus_matchchar):
            gm.append("can't have matchchars in the first sequences")
            raise Glitch, gm
        haveToDoIt = False
        # I can't remember why I had to do this next bit, so I am turning it off.
        #if len(self.equates):
        #    haveToDoIt = True # due to equates      Why is this?
        #else:    # see if there are any matchchars
        
        for aSeq in self.sequences[1:]:
            if string.count(aSeq, self.nexus_matchchar):
                haveToDoIt = True
                break
        if haveToDoIt:
            tempSeqs = []
            for aSeq in self.sequences:
                tempSeqs.append(array.array('c', aSeq))
            for aSeq in tempSeqs[1:]:
                for j in range(len(aSeq)):
                    if aSeq[j] == self.nexus_matchchar:
                        aSeq[j] = tempSeqs[0][j]
            # put them back into strings
            self.sequences = []
            for aSeq in tempSeqs:
                self.sequences.append(aSeq.tostring())

    def alignment(self):
        gm = ['NexusData.alignment()']

        # This method is used in a couple of places in Nexus.readNexusFile(), where fName is set.

        # At this point, matchchars have been propagated, and we have
        # checked that no sequences contain bad chars.  And checked
        # that no symbol occurs more than once (eg some char is both
        # an equate and a gapchar.)  And that neither gapchar nor
        # missingchar is one of the symbols or equate keys.
        
        if len(self.sequences):
            eKeys = ''.join(self.equates.keys())
            needToSwitchGaps = False
            needToSwitchMissing = False
            if self.nexus_gap:
                if self.nexus_gap != '-':
                    for aSeq in self.sequences:
                        if aSeq.count(self.nexus_gap):
                            needToSwitchGaps = True
                            break
            if self.nexus_missing:
                if self.nexus_missing != '?':
                    for aSeq in self.sequences:
                        if aSeq.count(self.nexus_missing):
                            needToSwitchMissing = True
                            break
                    
                
            
            #print "needToSwitchGaps = %s, needToSwitchMissing = %s" % (needToSwitchGaps, needToSwitchMissing)
            #print "self.nexus_gap=%s, self.nexus_missing=%s" % (self.nexus_gap, self.nexus_missing)
            if needToSwitchGaps and var.nexus_warnSwitchGapChar:
                print "Warning: nexus gap char '%s' will be changed to the more usual '-'" % self.nexus_gap
                print "(To turn off this warning, set var.nexus_warnSwitchGapChar = False.)"
            if needToSwitchMissing and var.nexus_warnSwitchMissingChar:
                print "Warning: nexus missing char '%s' will be changed to the more usual '?'" % self.nexus_missing
                print "(To turn off this warning, set var.nexus_warnSwitchMissingChar = False.)"

            if needToSwitchGaps and needToSwitchMissing:
                temp = '\t'
                for seqNum in range(len(self.sequences)):
                    aSeq = self.sequences[seqNum]
                    #print "before: %s" % aSeq
                    if aSeq.count(self.nexus_gap) or aSeq.count(self.nexus_missing):
                        aSeqList = list(aSeq)
                        for seqPos in range(len(aSeqList)):
                            if aSeqList[seqPos] == self.nexus_missing:
                                aSeqList[seqPos] = temp
                            if aSeqList[seqPos] == self.nexus_gap:
                                aSeqList[seqPos] = '-'
                            if aSeqList[seqPos] == temp:
                                aSeqList[seqPos] = '?'
                        aSeq = ''.join(aSeqList)
                        self.sequences[seqNum] = aSeq
                    #print " after: %s" % aSeq
            elif needToSwitchGaps:
                for seqNum in range(len(self.sequences)):
                    aSeq = self.sequences[seqNum]
                    if aSeq.count(self.nexus_gap):
                        aSeq = aSeq.replace(self.nexus_gap, '-')
                        self.sequences[seqNum] = aSeq
            elif needToSwitchMissing:
                for seqNum in range(len(self.sequences)):
                    aSeq = self.sequences[seqNum]
                    if aSeq.count(self.nexus_missing):
                        aSeq = aSeq.replace(self.nexus_missing, '?')
                        self.sequences[seqNum] = aSeq
                    
            a = Alignment()
            for i in range(self.nTax):
                s = Sequence()
                s.name = self.taxNames[i]
                s.sequence = self.sequences[i]
                s.dataType = self.dataType
                a.sequences.append(s)
            a.symbols = self.symbols
            a.dim = len(self.symbols)
            a.equates = self.equates
            a.checkLengthsAndTypes()
            a.checkNamesForDupes() # a SequenceList method, inherited

            return a
        else:
            print "\nNexusData.alignment()"
            print "    an alignment has been requested but the nexus "
            print "    file has no sequences.  Returning None."
            return None


