import sys
import os
import numpy
from p4.p4exceptions import P4Error

# A           Ala            Alanine
# R           Arg            Arginine
# N           Asn            Asparagine
# D           Asp            Aspartic acid
# C           Cys            Cysteine
# Q           Gln            Glutamine
# E           Glu            Glutamic acid
# G           Gly            Glycine
# H           His            Histidine
# I           Ile            Isoleucine
# L           Leu            Leucine
# K           Lys            Lysine
# M           Met            Methionine
# F           Phe            Phenylalanine
# P           Pro            Proline
# O           Pyl            Pyrrolysine (22nd amino acid)
# U           Sec            Selenocysteine (21st amino acid)
# S           Ser            Serine
# T           Thr            Threonine
# W           Trp            Tryptophan
# Y           Tyr            Tyrosine
# V           Val            Valine
# B           Asx            Aspartic acid or Asparagine
# Z           Glx            Glutamine or Glutamic acid
# J           Xle            Isoleucine or Valine (mass spec ambiguity)
# X           Xaa            Any or unknown amino acid


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
        #: file.  You can always see it again with :func:`p4.func.splash`.
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

        # Dangerous to change this.  Not widely implemented.
        self.allowDupedTaxonNames = False
        # Experimental hack.  It does not work on taxNames
        # -- just on the tree nodes.  Setting it to
        # True is verbose, and at least tells you that you
        # have dupes.  Setting it to 2 is silent --
        # especially dangerous.
        self.doRepairDupedTaxonNames = 0   # 0 means don't do it.
        # 1 means do it verbosely,
        # 2 means do it silently.

        # It is possible to have root-on-a-stick style trees,
        self.warnAboutTerminalRootWithNoName = True
        # where the root may or may not be a leaf, and may or
        # may not have a taxon name.  If it is a root on a stick,
        # it is by default considered a leaf, and if it has no
        # name then a warning is given, if this is set.

        # These comments are always read in by TreePartitions
        self.doTreeReadMcmcModelUsageComments = 0
        # objects.  This says whether they are read
        # in by Tree objects as well.

        #: Require that MCMC run numbers are in order.
        # Ensure chains are init'ed sequentially
        self.strictRunNumberChecking = False   # changed from True April 2019

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

        self.topologyDistanceMetrics = [
            'sd', 'wrf', 'bld', 'diffs', 'scqdist']  # 'triplet'? 'thquartet',
        """A list of metrics.

        Used by :meth:`Tree.Tree.topologyDistance` and :meth:`Trees.Trees.topologyDistanceMatrix`
        """

        # These are used by TreePicture
        self.TEXTDRAW_NONE = 0
        self.TEXTDRAW_COMP = 1
        self.TEXTDRAW_RMATRIX = 2
        self.TEXTDRAW_GDASRV = 3
        self.TEXTDRAW_NTHINGS = 4  # ie the max number

        self.recipesWriteToFile = True  # See func.recipes()

        # # no period, '!', #, @, _, &, ^, %, $
        self.nexus_punctuation = """()[]{}\\/,;:=*'"`+-<>"""  # no period, '!', #, @, _, &, ^, %, $
        self.nexus_safeChars = "!#@_&^%$"
        #self.phylip_punctuation = '\(\)\,\;\:'
        self.phylip_punctuation = '(),;:'
        self.punctuation = self.nexus_punctuation

        self.validDnaChars = 'acgt-?nrykmswbdhv'              # no u
        self.validRnaChars = 'acgu-?nrykmswbdhv'              # no t
        self.validNucleotideChars = 'acgt-?nurykmswbdhv'      # u and t, both
        # "*" (stop) for Genbank and gde.
        self.validProteinChars = 'acdefghiklmnpqrstvwy-?xbzju*'

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

        # A pointer to a gsl random number generator   SEEMS TO BE A LITTLE
        # MEMORY LEAK!
        self.gsl_rng = None
        self.doDataPart = 0        # Experimental

        # Modify behavior of NexusToken.nextTok() function.
        self._nexus_writeVisibleComments = numpy.array([0], numpy.int32)   # write, but do not get, all [!...]
        self._nexus_getP4CommandComments = numpy.array([0], numpy.int32)             # all [&&p4 ...]
        self._nexus_getWeightCommandComments = numpy.array([1], numpy.int32)         # all [&w ...]
        self._nexus_getAllCommandComments = numpy.array([0], numpy.int32)            # all [&...]
        self._nexus_getLineEndingsAsTokens = numpy.array([0], numpy.int32)

        self.rMatrixProteinSpecs = ['cpREV', 'd78', 'jtt', 'mtREV24', 'mtmam',
                                    'wag', 'rtRev', 'tmjtt94', 'tmlg99', 'lg',
                                    'blosum62', 'hivb', 'mtart', 'mtzoa',
                                    'gcpREV', 'stmtREV', 'vt']
        """A list of the currently available protein models."""

        self.rMatrixSpecs = [
            'ones', '2p', 'specified', 'optimized'] + self.rMatrixProteinSpecs
        self.compSpecs = [
            'equal', 'empirical', 'specified'] + self.rMatrixProteinSpecs

        self.modelSymbols = ['=', '@', '#', '$', '%', '&', '*', '+', 'a', 'b',
                             'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
                             'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
                             'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D', 'E', 'F',
                             'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
                             'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
                             '1', '2', '3', '4', '5', '6', '7', '8', '9', '<',
                             '!', '>', '^', '[', ']', ':', '(', ')', '~',
                             '"', '|'
                             ]  # no '-', '?'
        """A list of symbols used in text tree drawings, showing model disposition on the tree."""

        # This next variable (var.rMatrixNormalizeTo1 is a property)
        # determines whether the elements of a GTR-like rMatrix are
        # normalized to sum to 1, or else whether the bottom corner
        # rate (G->T in DNA GTR) is set to 1 and the other rate
        # elements are relative to it.
        self._rMatrixNormalizeTo1 = numpy.array([1], numpy.int32)

        self._PIVEC_MIN = numpy.array([1.0e-13], numpy.float64)      # Changed from 1e-18 March 2019, to avoid getting negative bigP values
        self._PIVEC_MAX = numpy.array([0.999], numpy.float64)      
        #self._RATE_MIN = 1.0e-14                                 # ie for rMatrices, changed from 1.0e-8 nov 2016, because self._rMatrixNormalizeTo1 is set
        self._RATE_MIN = numpy.array([1.0e-13], numpy.float64)    # changed from 1.e-14 March 2019, to avoid getting negative bigP values
        self._RATE_MAX = numpy.array([0.9999999], numpy.float64)  # changed from 1.0e8 nov 2016, also because self._rMatrixNormalizeTo1 is set
        self._GAMMA_SHAPE_MIN = numpy.array([0.1], numpy.float64)  # changed from 0.000001 april 2019, with alpha=0.1, slowest cat is 5.e-7
        self._GAMMA_SHAPE_MAX = numpy.array([300.0], numpy.float64)
        self._PINVAR_MIN = numpy.array([0.0], numpy.float64)
        self._PINVAR_MAX = numpy.array([0.99], numpy.float64)
        self._RELRATE_MIN = numpy.array([1.0e-8], numpy.float64)
        self._RELRATE_MAX = numpy.array([1.0e8], numpy.float64)
        self._KAPPA_MIN = numpy.array([0.000001], numpy.float64)
        self._KAPPA_MAX = numpy.array([100.0], numpy.float64)
        self._BRLEN_MIN = numpy.array([1.0e-8], numpy.float64)
        self._BRLEN_MAX = numpy.array([3.0], numpy.float64)

        self._newtAndBrentPowellOptPassLimit = numpy.array([50], numpy.int32)

        self.GAP_CODE = -1
        self.QMARK_CODE = -2
        self.N_LIKE = -3
        self.EQUATES_BASE = -64

        self.doMcmcSp = True  # speedier like calcs in mcmc

        # # Modify behavior of NexusToken.nextTok() function.
        # # write, but do not get, all [!...]
        # self._nexus_writeVisibleComments = 0
        # self._nexus_getP4CommandComments = 0             # all [&&p4 ...]
        # self._nexus_getWeightCommandComments = 1         # all [&w ...]
        # self._nexus_getAllCommandComments = 0            # all [&...]
        # self._nexus_getLineEndingsAsTokens = 0
        # self._rMatrixNormalizeTo1 = 1
        self._interactiveHelper = None
        self._excepthookEditor = None
        self.allowEmptyCharSetsAndTaxSets = False

        #self.mcmc_swapVector = False  # (old) matrix or (new) vector, both mcmc and stmcmc
        self.mcmc_swapTunerSampleSize = 250  # mcmc and stmcmc
        self.stmcmc_useFastSpa = False
        self.mcmc_sameBigTToStartOnAllChains = False # mcmc and stmcmc, for debugging, and fixed toplogy runs
        #self.mcmc_doTuneChainTemp = False
        self.mcmc_allowUnresolvedStartingTree = False
        self.mcmc_simTemp_tempCurveLogBase = 2.8    # higher -> more curvey; lower -> more linear
        self.mcmc_logTunings = False  # verbose logging of on-the-fly changes to tunings 

    def _del_nothing(self):
        gm = ["Don't/Can't delete this property."]
        raise P4Error(gm)

    def _getPIVEC_MIN(self):
        return self._PIVEC_MIN[0]
    def _setPIVEC_MIN(self, newValue):
        self._PIVEC_MIN[0] = newValue
    PIVEC_MIN = property(_getPIVEC_MIN, _setPIVEC_MIN, _del_nothing, "(property) PIVEC_MIN (float)")

    def _getPIVEC_MAX(self):
        return self._PIVEC_MAX[0]
    def _setPIVEC_MAX(self, newValue):
        self._PIVEC_MAX[0] = newValue
    PIVEC_MAX = property(_getPIVEC_MAX, _setPIVEC_MAX, _del_nothing, "(property) PIVEC_MAX (float)")

    def _getRATE_MIN(self):
        return self._RATE_MIN[0]
    def _setRATE_MIN(self, newValue):
        self._RATE_MIN[0] = newValue
    RATE_MIN = property(_getRATE_MIN, _setRATE_MIN, _del_nothing, "(property) RATE_MIN (float)")

    def _getRATE_MAX(self):
        return self._RATE_MAX[0]
    def _setRATE_MAX(self, newValue):
        self._RATE_MAX[0] = newValue
    RATE_MAX = property(_getRATE_MAX, _setRATE_MAX, _del_nothing, "(property) RATE_MAX (float)")

    def _getGAMMA_SHAPE_MIN(self):
        return self._GAMMA_SHAPE_MIN[0]
    def _setGAMMA_SHAPE_MIN(self, newValue):
        self._GAMMA_SHAPE_MIN[0] = newValue
    GAMMA_SHAPE_MIN = property(_getGAMMA_SHAPE_MIN, _setGAMMA_SHAPE_MIN, _del_nothing, "(property) GAMMA_SHAPE_MIN (float)")

    def _getGAMMA_SHAPE_MAX(self):
        return self._GAMMA_SHAPE_MAX[0]
    def _setGAMMA_SHAPE_MAX(self, newValue):
        self._GAMMA_SHAPE_MAX[0] = newValue
    GAMMA_SHAPE_MAX = property(_getGAMMA_SHAPE_MAX, _setGAMMA_SHAPE_MAX, _del_nothing, "(property) GAMMA_SHAPE_MAX (float)")

    def _getPINVAR_MIN(self):
        return self._PINVAR_MIN[0]
    def _setPINVAR_MIN(self, newValue):
        self._PINVAR_MIN[0] = newValue
    PINVAR_MIN = property(_getPINVAR_MIN, _setPINVAR_MIN, _del_nothing, "(property) PINVAR_MIN (float)")

    def _getPINVAR_MAX(self):
        return self._PINVAR_MAX[0]
    def _setPINVAR_MAX(self, newValue):
        self._PINVAR_MAX[0] = newValue
    PINVAR_MAX = property(_getPINVAR_MAX, _setPINVAR_MAX, _del_nothing, "(property) PINVAR_MAX (float)")

    def _getRELRATE_MIN(self):
        return self._RELRATE_MIN[0]
    def _setRELRATE_MIN(self, newValue):
        self._RELRATE_MIN[0] = newValue
    RELRATE_MIN = property(_getRELRATE_MIN, _setRELRATE_MIN, _del_nothing, "(property) RELRATE_MIN (float)")

    def _getRELRATE_MAX(self):
        return self._RELRATE_MAX[0]
    def _setRELRATE_MAX(self, newValue):
        self._RELRATE_MAX[0] = newValue
    RELRATE_MAX = property(_getRELRATE_MAX, _setRELRATE_MAX, _del_nothing, "(property) RELRATE_MAX (float)")

    def _getKAPPA_MIN(self):
        return self._KAPPA_MIN[0]
    def _setKAPPA_MIN(self, newValue):
        self._KAPPA_MIN[0] = newValue
    KAPPA_MIN = property(_getKAPPA_MIN, _setKAPPA_MIN, _del_nothing, "(property) KAPPA_MIN (float)")

    def _getKAPPA_MAX(self):
        return self._KAPPA_MAX[0]
    def _setKAPPA_MAX(self, newValue):
        self._KAPPA_MAX[0] = newValue
    KAPPA_MAX = property(_getKAPPA_MAX, _setKAPPA_MAX, _del_nothing, "(property) KAPPA_MAX (float)")

    def _getBRLEN_MIN(self):
        return self._BRLEN_MIN[0]
    def _setBRLEN_MIN(self, newValue):
        self._BRLEN_MIN[0] = newValue
    BRLEN_MIN = property(_getBRLEN_MIN, _setBRLEN_MIN, _del_nothing, "(property) BRLEN_MIN (float)")

    def _getBRLEN_MAX(self):
        return self._BRLEN_MAX[0]
    def _setBRLEN_MAX(self, newValue):
        self._BRLEN_MAX[0] = newValue
    BRLEN_MAX = property(_getBRLEN_MAX, _setBRLEN_MAX, _del_nothing, "(property) BRLEN_MAX (float)")

    def _get_newtAndBrentPowellOptPassLimit(self):
        return self._newtAndBrentPowellOptPassLimit[0]
    def _set_newtAndBrentPowellOptPassLimit(self, newValue):
        self._newtAndBrentPowellOptPassLimit[0] = newValue
    newtAndBrentPowellOptPassLimit = property(_get_newtAndBrentPowellOptPassLimit, 
                                              _set_newtAndBrentPowellOptPassLimit, 
                                              _del_nothing, "(property) newtAndBrentPowellOptPassLimit (int)")


    # self._nexus_writeVisibleComments = 0             # write, but do not get, all [!...]
    # self._nexus_getP4CommandComments = 0             # all [&&p4 ...]
    # self._nexus_getWeightCommandComments = 1         # all [&w ...]
    # self._nexus_getAllCommandComments = 0            # all [&...]

    def _get_nexus_writeVisibleComments(self):
        return self._nexus_writeVisibleComments[0]

    def _set_nexus_writeVisibleComments(self, newVal):
        try:
            newVal = int(newVal)
        except:
            gm = ['This property should be set to an int.']
            raise P4Error(gm)
        self._nexus_writeVisibleComments[0] = newVal

    nexus_writeVisibleComments = property(_get_nexus_writeVisibleComments,
                                          _set_nexus_writeVisibleComments, _del_nothing)
    """(property) nexus_writeVisibleComments (int)"""

    def _get_nexus_getP4CommandComments(self):
        return self._nexus_getP4CommandComments[0]

    def _set_nexus_getP4CommandComments(self, newVal):
        try:
            newVal = int(newVal)
        except:
            gm = ['This property should be set to an int.']
            raise P4Error(gm)
        self._nexus_getP4CommandComments[0] = newVal

    nexus_getP4CommandComments = property(_get_nexus_getP4CommandComments,
                                          _set_nexus_getP4CommandComments, _del_nothing)
    """Whether to get p4-specific command comments in NEXUS trees."""

    def _get_nexus_getWeightCommandComments(self):
        return self._nexus_getWeightCommandComments[0]

    def _set_nexus_getWeightCommandComments(self, newVal):
        try:
            newVal = int(newVal)
        except:
            gm = ['This property should be set to an int.']
            raise P4Error(gm)
        self._nexus_getWeightCommandComments[0] = newVal

    nexus_getWeightCommandComments = property(_get_nexus_getWeightCommandComments,
                                              _set_nexus_getWeightCommandComments, _del_nothing)
    """Whether to get the 'weight' command comment in a NEXUS tree."""

    def _get_nexus_getAllCommandComments(self):
        return self._nexus_getAllCommandComments[0]

    def _set_nexus_getAllCommandComments(self, newVal):
        try:
            newVal = int(newVal)
        except:
            gm = ['This property should be set to an int.']
            raise P4Error(gm)
        self._nexus_getAllCommandComments[0] = newVal

    nexus_getAllCommandComments = property(_get_nexus_getAllCommandComments,
                                           _set_nexus_getAllCommandComments, _del_nothing)
    """Whether to get all command comments in NEXUS tree files."""

    def _get_nexus_getLineEndingsAsTokens(self):
        return self._nexus_getLineEndingsAsTokens[0]

    def _set_nexus_getLineEndingsAsTokens(self, newVal):
        try:
            newVal = int(newVal)
        except:
            gm = ['This property should be set to an int.']
            raise P4Error(gm)
        self._nexus_getLineEndingsAsTokens[0] = newVal

    nexus_getLineEndingsAsTokens = property(_get_nexus_getLineEndingsAsTokens,
                                            _set_nexus_getLineEndingsAsTokens, _del_nothing)
    """(property) Whether line endings get returned as nexus tokens (int)"""

    def _get_rMatrixNormalizeTo1(self):
        return self._rMatrixNormalizeTo1[0]

    def _set_rMatrixNormalizeTo1(self, newVal):
        # try:
        #    newVal = int(newVal)
        #    self._rMatrixNormalizeTo1[0] = newVal
        # except:
        #    gm = ['This property should be set to an int.']
        #    raise P4Error(gm)
        gm = ["Var._set_rMatrixNormalizeTo1()"]
        gm.append(
            "This fundamental variable affects array lengths, and so should not be changed during a run--")
        gm.append(
            "It should only be set in Var.py, only to be read at start-up.")

        raise P4Error(gm)

    rMatrixNormalizeTo1 = property(_get_rMatrixNormalizeTo1,
                                   _set_rMatrixNormalizeTo1, _del_nothing)
    """(property) not user-settable"""

    def _get_interactiveHelper(self):
        return self._interactiveHelper

    def _set_interactiveHelper(self, newVal):
        goodValues = [None, 'bpython', 'ipython']
        if newVal in goodValues:
            self._interactiveHelper = newVal
        else:
            gm = ['This property should be set to one of %s' % goodValues]
            raise P4Error(gm)

    interactiveHelper = property(_get_interactiveHelper,
                                 _set_interactiveHelper, _del_nothing)
    """For interactive use, set the helper.

    Set to bpython, or ipython.  Default is None.
    """

    def _get_excepthookEditor(self):
        return self._excepthookEditor

    def _set_excepthookEditor(self, newVal):
        assert newVal == None or isinstance(newVal,str)
        self._excepthookEditor = newVal

    excepthookEditor = property(_get_excepthookEditor,
                                 _set_excepthookEditor, _del_nothing)
    """The editor that is called by excepthook, or None

    Setting this to None or to the name of an editor acts as a Boolean to say
    whether to follow a traceback with a hook to call your editor.

    It is for using p4 at the terminal, writing scripts and source in an editor
    such as emacs or vi.  The various files and lines of the traceback are given
    numbers that you can type in to get to the source.  Handy!  

    The editor needs to be clever enough to use the "+N" command line option to
    be able to go to a particular line (N).  Emacs and vi can both do this.

    Set to, for example 'emacsclient -n', or 'vim'.  Default is None.

    """


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
                # It was not installed.  It is running 'in-place' -- no
                # installation.py file.
                import p4
                pth = p4.__file__
                # print pth
                pth = os.path.split(pth)[0]
                pth = os.path.split(pth)[0]
                # print pth
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
            raise P4Error(gm)

    examplesDir = property(_getExamplesDir, _setExamplesDir, _del_nothing)
    """The directory where the p4 examples can be found.

    Assuming that p4 can find them.
    """


# Make a single instance.
var = Var()
