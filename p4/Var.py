# -*- coding: latin-1 -*-
import sys,os
#import numpy
from Glitch import Glitch

##          A           Ala            Alanine
##          R           Arg            Arginine
##          N           Asn            Asparagine
##          D           Asp            Aspartic acid
##          C           Cys            Cysteine
##          Q           Gln            Glutamine
##          E           Glu            Glutamic acid
##          G           Gly            Glycine
##          H           His            Histidine
##          I           Ile            Isoleucine
##          L           Leu            Leucine
##          K           Lys            Lysine
##          M           Met            Methionine
##          F           Phe            Phenylalanine
##          P           Pro            Proline
##                  O           Pyl            Pyrrolysine (22nd amino acid)
##                  U           Sec            Selenocysteine (21st amino acid)
##          S           Ser            Serine
##          T           Thr            Threonine
##          W           Trp            Tryptophan
##          Y           Tyr            Tyrosine
##          V           Val            Valine
##          B           Asx            Aspartic acid or Asparagine
##          Z           Glx            Glutamine or Glutamic acid
##                  J           Xle            Isoleucine or Valine (mass spec ambiguity)
##          X           Xaa            Any or unknown amino acid


class Var(object):
    """An instance of this class is the global bucket.

    An instance of this class is made, called 'var'.  It can be read
    and changed by the user, and it is imported by other classes that
    need access to 'global' variables contained in it.
    """

    def __init__(self):

        
        self._examplesDir = None  # self.examplesDir is a property, see below

        #: Files read in from the command line.
        self.fileNames = []       

        #: A list of alignments.  Only those alignments read from the command line and from read()
        self.alignments = []
        
        #: A list of SequenceList objects from files from the command line or from read() 
        self.sequenceLists = []   

        #: Tree objects from files from the command line or from read() 
        self.trees = []

        #: A single NexusSets object (not a list) from files from the command line
        #: or from read().  If more sets info is read in, it all gets
        #: put in the single NexusSets object
        self.nexusSets = None     

        #: Whether to do the splash screen on startup with no args.
        #: When you get tired of it, turn it off in a customization
        #: file.  You can always see it again with :func:`func.splash`.
        self.doSplash = 1

        #: Running commentary on p4's attempts to read in files.
        self.verboseRead = 0            

        #: If you give read() a string that is not a file, p4 points
        #: it out, in case its a mistake.
        self.warnReadNoFile = 1         
                                        
        #: In the p4 script, non-filename command-line args are accommodated,
        #: awkwardly.  They need to be put after                   
        #: a double dash.  They can then be found                                                     
        #: here.
        self.argvAfterDoubleDash = None
        
        self.allowDupedTaxonNames = False  # Dangerous to change this.  Not widely implemented.
                                           # Experimental hack.  It does not work on taxNames
                                           # -- just on the tree nodes.  Setting it to
                                           # True is verbose, and at least tells you that you
                                           # have dupes.  Setting it to 2 is silent --
                                           # especially dangerous.
        self.doRepairDupedTaxonNames = 0   # 0 means don't do it.
                                           # 1 means do it verbosely,
                                           # 2 means do it silently.

        self.warnAboutTerminalRootWithNoName = True  # It is possible to have root-on-a-stick style trees,
                                                     # where the root may or may not be a leaf, and may or
                                                     # may not have a taxon name.  If it is a root on a stick,
                                                     # it is by default considered a leaf, and if it has no
                                                     # name then a warning is given, if this is set.
                                                     
        self.doTreeReadMcmcModelUsageComments = 0  # These comments are always read in by TreePartitions
                                                   # objects.  This says whether they are read
                                                   # in by Tree objects as well.

        #: Require that MCMC run numbers are in order.
        self.strictRunNumberChecking = True #Ensure chains are init'ed sequentially

        self.NO_ORDER = -10000  # pre and postOrder for unused nodes.

        #: Beast trees use command comments eg [&stuff=0.12345] a lot.  To
        #: read them, turn var.nexus_getAllCommandComments on, and then
        #: turn this on.
        self.nexus_readBeastTreeCommandComments = False

        #: The NEXUS spec does not allow names that are all numerals.
        #: Temporarily override by True-ing this.
        self.nexus_allowAllDigitNames = False

        #: In Nexus data, the characters for gap and missing can be set
        #: to almost anything you like, but I like them to be '-' and
        #: '?', which is more usual.  So if needed, I switch them.
        #: Dodgy! -- But at least I warn the user, unless that warning
        #: is turned off.
        self.nexus_warnSwitchGapChar = True
        self.nexus_warnSwitchMissingChar = True
        """See :attr:`Var.nexus_warnSwitchGapChar`"""

        #: In nexus, if you have dna, rna, or protein datatype, you
        #: can, according to the nexus standard, add extra symbols,
        #: defined in the format command.  I don't like that, so I
        #: don't allow it.  However, at least I warn the user, either
        #: as a warning (if True) or a glitch (if False).
        self.nexus_ignoreFormatCommandSymbols = True

        self.nexus_warnSkipUnknownBlock = True
        """Print a little message if p4 encounters a NEXUS block that it does not know what to do with."""
        
        self.nexus_allowUTREE = False
        """UTREE has been deprecated since at least 1997.

        But if you need to read UTREE-containing tree files, you can
        turn it on with this.
        """

        self.newick_allowSpacesInNames = False
        """Allow unquoted taxon names in Newick and Nexus trees to have spaces.

        This will allow reading trees like (A, B C, (D E, F)G H);
        Multiple spaces are collapsed into single spaces.
        """

        self.topologyDistanceMetrics = ['sd', 'wrf', 'bld', 'diffs', 'scqdist']  # 'triplet'? 'thquartet',
        """A list of metrics.

        Used by :meth:`Tree.Tree.topologyDistance` and :meth:`Trees.Trees.topologyDistanceMatrix`
        """

        # These are used by TreePicture
        self.TEXTDRAW_NONE = 0
        self.TEXTDRAW_COMP = 1
        self.TEXTDRAW_RMATRIX = 2
        self.TEXTDRAW_GDASRV = 3
        self.TEXTDRAW_NTHINGS = 4 # ie the max number

        self.recipesWriteToFile=True  # See func.recipes()
        
        self.nexus_punctuation = '\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>'  # no period, '!', #, @, _, &, ^, %, $
        self.phylip_punctuation = '\(\)\,\;\:'
        self.punctuation = self.nexus_punctuation
        
        self.validDnaChars = 'acgt-?nrykmswbdhv'              # no u
        self.validRnaChars = 'acgu-?nrykmswbdhv'              # no t
        self.validNucleotideChars = 'acgt-?nurykmswbdhv'      # u and t, both
        self.validProteinChars = 'acdefghiklmnpqrstvwy-?xbzju*' # "*" (stop) for Genbank and gde.
        
        self.phylipDataMaxNameLength = 10
        """The length of the taxon name in phylip data.

        Unfortunately the phylip data format has a lot of variations.
        Even in the 'official' phylip data format, the name length is
        a compile-time variable.  Its usually 10, tho.  You can change
        it with this."""

        self.writeFastaUppercase = False
        """Whether fasta sequences are written in uppercase."""


        # Check sequence files that are read in and become Alignment objects.
        self.doCheckForAllGapColumns = True  
        self.doCheckForBlankSequences = True 
        self.doCheckForDuplicateSequences = True
        self.doCheckForDuplicateSequenceNames = True


        

        self.SAME = 0
        self.DIFFERENT = 1
        self.OK = 0
        self.NOT_OK = 1

       

        # The python that shiped with MacOS was usually old and was
        # not built with readline.  Thats ok, 'cuz it was easy enough
        # to install a new python that had readline.  The one from
        # pythonmac.org was fine.  Now the python that comes with
        # MacOS 10.5 is fairly up-to-date, and it is built with
        # readline -- sort of.  The problem is that the readline
        # python module does not use GNU libreadline, it uses the
        # editline library (libedit ?), and it has a different syntax
        # for parse_and_bind().  And it is partly broken, in that it
        # does not pass a 'text' arg to the drop-in 'complete'
        # function.
        #
        # It hurts.
        #
        # I suspect that the pythonmac.org python is still built with
        # gnu readline, and is ok, but I have not confirmed.  If you
        # are using the Mac python that uses editline, then change the
        # following to 'True'.  Completion will be slightly broken, in
        # that getting the argspec inside functions and methods will
        # not work, but the rest of it seems to work.
        #
        # Update on the above.  I got my python2.5, and now python2.6
        # from MacPorts, and they were both fine.  And 2.7 from HomeBrew.
        self.readlineUsesEditline = False


        try:
            import pf
            import numpy
            self.usePfAndNumpy = True
        except ImportError:
            self.usePfAndNumpy = False
        #self.usePfAndNumpy = False
        
        if self.usePfAndNumpy:
            self.gsl_rng = None # A pointer to a gsl random number generator   SEEMS TO BE A LITTLE MEMORY LEAK!
            self.doDataPart = 0        # Experimental
            
            # Modify behavior of NexusToken.nextTok() function.
            self._nexus_writeVisibleComments = numpy.array([0], numpy.int32)   # write, but do not get, all [!...]
            self._nexus_getP4CommandComments = numpy.array([0], numpy.int32)             # all [&&p4 ...]
            self._nexus_getWeightCommandComments = numpy.array([1], numpy.int32)         # all [&w ...]
            self._nexus_getAllCommandComments = numpy.array([0], numpy.int32)            # all [&...]
            self._nexus_getLineEndingsAsTokens = numpy.array([0], numpy.int32)
            self.nexus_doFastNextTok = True      # nextTok in C, from NexusToken2.  Does not work for CStrings.

            self.rMatrixProteinSpecs = ['cpREV', 'd78', 'jtt', 'mtREV24', 'mtmam',
                                        'wag', 'rtRev', 'tmjtt94', 'tmlg99', 'lg',
                                        'blosum62', 'hivb', 'mtart', 'mtzoa', 'gcpREV']
            """A list of the currently available protein models."""
            
            self.rMatrixSpecs = ['ones', '2p', 'specified', 'optimized'] + self.rMatrixProteinSpecs
            self.compSpecs = ['equal', 'empirical', 'specified'] + self.rMatrixProteinSpecs

            self.modelSymbols = ['=', '@', '#', '$', '%', '&', '*', '+', 'a', 'b',
                                 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
                                 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
                                 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D', 'E', 'F',
                                 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
                                 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
                                 '1', '2', '3', '4', '5', '6', '7', '8', '9', '<',
                                 '!', '>', '£', '^', '[', ']', ':', '(', ')', '~',
                                 '"', '|'
                                ] # no '-', '?'
            """A list of symbols used in text tree drawings, showing model disposition on the tree."""
            

            self.PIVEC_MIN = 1e-18 # see Pf/defines.h
            self.RATE_MIN = 1.0e-8 # ie for rMatrices
            self.RATE_MAX = 1.0e8
            self.GAMMA_SHAPE_MIN = 0.000001
            self.GAMMA_SHAPE_MAX = 300.0
            self.PINVAR_MIN = 0.0
            self.PINVAR_MAX = 0.99
            self.RELRATE_MIN = 1.0e-8
            self.RELRATE_MAX = 1.0e8
            self.KAPPA_MIN = 0.000001
            self.KAPPA_MAX = 100.0
            self.BRLEN_MIN = 1.0e-8
            self.BRLEN_MAX = 3.0


            self.GAP_CODE = -1
            self.QMARK_CODE = -2
            self.N_LIKE = -3
            self.EQUATES_BASE = -64

            self.doMcmcSp = True # speedier like calcs in mcmc

            # This next variable (var.rMatrixNormalizeTo1 is a property)
            # determines whether the elements of a GTR-like rMatrix are
            # normalized to sum to 1, or else whether the bottom corner
            # rate (G->T in DNA GTR) is set to 1 and the other rate
            # elements are relative to it.
            self._rMatrixNormalizeTo1 = numpy.array([1], numpy.int32)

            self.mcmcMaxCompAndRMatrixTuning = 0.9

            #self.allowUnusedComps = False
            self.rjCompUniformAllocationPrior = True
            self.rjRMatrixUniformAllocationPrior = True

        else:
            # Modify behavior of NexusToken.nextTok() function.
            self._nexus_writeVisibleComments = 0   # write, but do not get, all [!...]
            self._nexus_getP4CommandComments = 0             # all [&&p4 ...]
            self._nexus_getWeightCommandComments = 1         # all [&w ...]
            self._nexus_getAllCommandComments = 0            # all [&...]
            self._nexus_getLineEndingsAsTokens = 0
            self.nexus_doFastNextTok = False      # nextTok in C, from NexusToken2.  Does not work for CStrings.
            self._rMatrixNormalizeTo1 = 1
            



        

    def _del_nothing(self):
        gm = ["Don't/Can't delete this property."]
        raise Glitch, gm

    ##    self._nexus_writeVisibleComments = 0             # write, but do not get, all [!...]
    ##    self._nexus_getP4CommandComments = 0             # all [&&p4 ...]
    ##    self._nexus_getWeightCommandComments = 1         # all [&w ...]
    ##    self._nexus_getAllCommandComments = 0            # all [&...]

    def _get_nexus_writeVisibleComments(self):
        if self.usePfAndNumpy:
            return self._nexus_writeVisibleComments[0]
        else:
            return self._nexus_writeVisibleComments

    def _set_nexus_writeVisibleComments(self, newVal):
        try:
            newVal = int(newVal)
        except:
            gm = ['This property should be set to an int.']
            raise Glitch, gm
        if self.usePfAndNumpy:
            self._nexus_writeVisibleComments[0] = newVal
        else:
            self._nexus_writeVisibleComments = newVal

    nexus_writeVisibleComments = property(_get_nexus_writeVisibleComments,
                                          _set_nexus_writeVisibleComments, _del_nothing)

    def _get_nexus_getP4CommandComments(self):
        if self.usePfAndNumpy:
            return self._nexus_getP4CommandComments[0]
        else:
            return self._nexus_getP4CommandComments

    def _set_nexus_getP4CommandComments(self, newVal):
        try:
            newVal = int(newVal)
        except:
            gm = ['This property should be set to an int.']
            raise Glitch, gm
        if self.usePfAndNumpy:
            self._nexus_getP4CommandComments[0] = newVal
        else:
            self._nexus_getP4CommandComments[0] = newVal

    nexus_getP4CommandComments = property(_get_nexus_getP4CommandComments,
                                          _set_nexus_getP4CommandComments, _del_nothing)
    """Whether to get p4-specific command comments in NEXUS trees."""

    def _get_nexus_getWeightCommandComments(self):
        if self.usePfAndNumpy:
            return self._nexus_getWeightCommandComments[0]
        else:
            return self._nexus_getWeightCommandComments

    def _set_nexus_getWeightCommandComments(self, newVal):
        try:
            newVal = int(newVal)
        except:
            gm = ['This property should be set to an int.']
            raise Glitch, gm
        if self.usePfAndNumpy:
            self._nexus_getWeightCommandComments[0] = newVal
        else:
            self._nexus_getWeightCommandComments = newVal

    nexus_getWeightCommandComments = property(_get_nexus_getWeightCommandComments,
                                              _set_nexus_getWeightCommandComments, _del_nothing)
    """Whether to get the 'weight' command comment in a NEXUS tree."""

    def _get_nexus_getAllCommandComments(self):
        if self.usePfAndNumpy:
            return self._nexus_getAllCommandComments[0]
        else:
            return self._nexus_getAllCommandComments

    def _set_nexus_getAllCommandComments(self, newVal):
        try:
            newVal = int(newVal)
        except:
            gm = ['This property should be set to an int.']
            raise Glitch, gm
        if self.usePfAndNumpy:
            self._nexus_getAllCommandComments[0] = newVal
        else:
            self._nexus_getAllCommandComments = newVal

    nexus_getAllCommandComments = property(_get_nexus_getAllCommandComments,
                                           _set_nexus_getAllCommandComments, _del_nothing)
    """Whether to get all command comments in NEXUS tree files."""

    def _get_nexus_getLineEndingsAsTokens(self):
        if self.usePfAndNumpy:
            return self._nexus_getLineEndingsAsTokens[0]
        else:
            return self._nexus_getLineEndingsAsTokens

    def _set_nexus_getLineEndingsAsTokens(self, newVal):
        try:
            newVal = int(newVal)
        except:
            gm = ['This property should be set to an int.']
            raise Glitch, gm
        if self.usePfAndNumpy:
            self._nexus_getLineEndingsAsTokens[0] = newVal
        else:
            self._nexus_getLineEndingsAsTokens = newVal

    nexus_getLineEndingsAsTokens = property(_get_nexus_getLineEndingsAsTokens,
                                            _set_nexus_getLineEndingsAsTokens, _del_nothing)

    def _get_rMatrixNormalizeTo1(self):
        if self.usePfAndNumpy:
            return self._rMatrixNormalizeTo1[0]
        else:
            return self._rMatrixNormalizeTo1

    def _set_rMatrixNormalizeTo1(self, newVal):
        #try:
        #    newVal = int(newVal)
        #    self._rMatrixNormalizeTo1[0] = newVal
        #except:
        #    gm = ['This property should be set to an int.']
        #    raise Glitch, gm
        gm = ["Var._set_rMatrixNormalizeTo1()"]
        gm.append("This fundamental variable affects array lengths, and so should not be changed during a run--")
        gm.append("It should only be set in Var.py, only to be read at start-up.")

        raise Glitch, gm

    rMatrixNormalizeTo1 = property(_get_rMatrixNormalizeTo1,
                                   _set_rMatrixNormalizeTo1, _del_nothing)

    def _getExamplesDir(self):
        if self._examplesDir:
            return self._examplesDir
        else:
            try:
                from installation import p4ExamplesDir
                if os.path.exists(os.path.join(p4ExamplesDir, 'A_quickStart')) and \
                       os.path.exists(os.path.join(p4ExamplesDir, 'W_recipes')):
                    self._examplesDir = p4ExamplesDir
                    return p4ExamplesDir
            except ImportError:
                # It was not installed.  It is running 'in-place' -- no installation.py file.
                import p4
                pth = p4.__file__
                #print pth
                pth = os.path.split(pth)[0]
                pth = os.path.split(pth)[0]
                #print pth
                if os.path.exists(pth):
                    pth = os.path.join(pth, 'share')
                    if os.path.exists(pth):
                        pth = os.path.join(pth, 'Examples')
                        if os.path.exists(os.path.join(pth, 'A_quickStart')) and \
                               os.path.exists(os.path.join(pth, 'W_recipes')):
                            self._examplesDir = pth
                            return pth
            
    def _setExamplesDir(self, theDir):
        gm = ["Var._setExamplesDir()"]
        # No whitespace, no trailing slash
        theDir = theDir.strip()
        if theDir.endswith("/"):
            theDir = theDir[:-1]
        # if theDir does not exist, the following throws an OSError 
        theDirList = os.listdir(theDir)
        if 'A_quickStart' in theDirList and 'W_recipes' in theDirList:
            self._examplesDir = theDir
        else:
            gm.append("The directory '%s'" % theDir)
            gm.append("does not appear to contain the p4 Examples.")
            raise Glitch, gm

    examplesDir = property(_getExamplesDir, _setExamplesDir, _del_nothing)
    """The directory where the p4 examples can be found.

    Assuming that p4 can find them.
    """



# Make a single instance.
var = Var()

