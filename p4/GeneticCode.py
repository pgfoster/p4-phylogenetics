import string
import sys


class GeneticCode:
    """A container for NCBI translation tables.

    See the ncbi translation tables, which this week are at
    http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c

    (If they move, poke around the 'taxonomy browser' area.)

    This week we have

    - **1** standard
    - **2** vertebrate mito
    - **3** yeast mito
    - **4** Mold, Protozoan
    - **5** invertbrate mito
    - **6** The Ciliate, Dasycladacean and Hexamita Nuclear Code
    - **9** echinoderm and flatworm mito
    - **10**  Euplotid Nuclear Code
    - **11**  Bacterial and Plant Plastid Code
    -  **12** Alternative Yeast Nuclear Code
    -  **13** Ascidian Mitochondrial Code
    -  **14** Alternative Flatworm Mitochondrial Code
    -  **21** Trematode Mitochondrial Code
    -  **24** Pterobranchia mito

    If more transl_tables are needed, you should be able to just drop
    them in, with a little tweaking.

    This provides
    
    - **code**  A dictionary.  So you can ask for eg myGC.code['ggg']
    - **codonsForAA**  Another dictionary, where you can ask for eg myGC.codonsForAA['v']
    - **startList**  A list of start codons

    **Methods**
     .. autosummary::

        GeneticCode.translate
        GeneticCode.wise2Table

    Wikipedia says: The joint nomenclature committee of the
    IUPAC/IUBMB has officially recommended the three-letter symbol Sec
    and the one-letter symbol U for selenocysteine.  The UGA codon is
    made to encode selenocysteine by the presence of a SECIS element
    (SElenoCysteine Insertion Sequence) in the mRNA.
    
    """
    
    def __init__(self, transl_table=1):
        self.transl_table = transl_table
        self.code = {}
        self.codonsForAA = {}
        self.startList = []

        if transl_table == 1: # standard
            AAs    = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '---M---------------M---------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 2: # vertebrate mito
            AAs      = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG'
            Starts   = '--------------------------------MMMM---------------M------------'
            Base1    = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2    = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3    = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 3:  #3. The Yeast Mitochondrial Code (transl_table=3)
            AAs    = 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '----------------------------------MM----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 4: # Mold, Protozoan,
                                # and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '--MM---------------M------------MMMM---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 5: # invertebrate mito
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG'
            Starts = '---M----------------------------MMMM---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 6: # The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
            AAs    = 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '-----------------------------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        # tables 7 and 8 have been deleted from NCBI.
        
        elif transl_table == 9: # echinoderm and flatworm mito
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG'
            Starts = '-----------------------------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 10: # The Euplotid Nuclear Code (transl_table=10)
            AAs    = 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '-----------------------------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 11: # The Bacterial and Plant Plastid Code (transl_table=11)
            AAs    = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '---M---------------M------------MMMM---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 12: # The Alternative Yeast Nuclear Code (transl_table=12)
            AAs    = 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            Starts = '-------------------M---------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 13: # The Ascidian Mitochondrial Code (transl_table=13)
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG'
            Starts = '---M------------------------------MM---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 14: # The Alternative Flatworm Mitochondrial Code (transl_table=14)
            AAs    = 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG'
            Starts = '-----------------------------------M----------------------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
            
        elif transl_table == 21: # Trematode Mitochondrial Code (transl_table=21)
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG'
            Starts = '-----------------------------------M---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

        elif transl_table == 24: # Pterobranchia mitochondrial code (transl_table=24)
            AAs    = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG'
            Starts = '---M---------------M---------------M---------------M------------'
            Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
            Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
            Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'


        else:
            print "GeneticCode: I don't know transl_table %i.  Get it from NCBI and add it!" % transl_table
            sys.exit()

        for i in range(64):
            theCodon = string.lower(Base1[i] + Base2[i] + Base3[i])
            theAA = string.lower(AAs[i])
            self.code[theCodon] = theAA
            if self.codonsForAA.has_key(theAA):
                self.codonsForAA[theAA].append(theCodon)
            else:
                self.codonsForAA[theAA] = [theCodon]
            if Starts[i] == 'M':
                self.startList.append(theCodon)

        if 1:
            self.codonsForAA['b'] = self.codonsForAA['n'] + self.codonsForAA['d']
            self.codonsForAA['z'] = self.codonsForAA['q'] + self.codonsForAA['e']

        if 1:
            self.codonsForAA['u'] = ['tga']  # selenocysteine

        if 0:
            k = self.code.keys()
            k.sort()
            for aKey in k:
                print "%20s  %-30s" % (aKey, self.code[aKey])

    def wise2Table(self):
        """Output in a form suitable to replace ``codon.table`` in ``genewise`` in Wise2.

        See the `Wise2 web site. <http://www.ebi.ac.uk/Tools/Wise2/doc_wise2.html>`_

        By default, ``genewise`` from the Wise2 package uses the
        standard genetic code, defined in the file ``codon.table``.
        However, you can supply your own table, and you can use this
        method to make a suitable file.  Due to lazy programming, this
        method prints to stdout, so you will need to put the output in
        a file yourself.  Put that file, suitably named (eg
        ``codon.table5`` or whatever) in the ``wisecfg`` directory
        (where the original ``codon.table`` resides), which might be
        ``/usr/local/src/wise2.2.0/wisecfg`` or some such location.

        Then, when you call ``genewise`` you can use the ``-codon``
        option to set the codon table file that you want to use, eg::

            genewise -genes -cdna -trans -pep -pretty -silent -codon codon.table5 guideFileName dnaFileName
            
        """

        print "! this is a codon table"
        print "! by p4, for ncbi %i" % self.transl_table
        for first in "tcag":
            for second in "tcag":
                for third in "tcag":
                    lcod = "%s%s%s" % (first,second,third)
                    ret = self.code[lcod]
                    if ret == '*':
                        ret = 'X'
                    print lcod.upper(), ret.upper()
                    
                    
    def translate(self, theCodon, verbose=1):
        """Translate a codon, handling ambiguities.

        This method will translate a codon, depending of course on the
        transl_table (which is specific to *self*), correctly handling
        ambiguities appropriate to the transl_table.  It does not give
        info about whether the codon is potentially a start codon.

        This method is used by the methods
        :meth:`Alignment.Alignment.translate` and
        :meth:`Alignment.Alignment.checkTranslation`.
        
        A translation like that from codon ``ggg`` to amino acid ``g``
        is direct and easy.  However, this method will also translate
        ambiguous codons ``ggy`` or ``ggs`` (and so on) to ``g``,
        unequivocally, because the four codons that start with ``gg``
        all code for ``g`` (which appears to be true for all
        translation tables, but that is not a requirement here).

        If all possible disambiguations of an ambiguous codon code for
        a particular amino acid then this method will return that
        amino acid; otherwise the translation is ambiguous and this
        method returns ``x``. 

        So for example the codon ``tgr`` will translate to ``x`` using
        transl_table = 1 (standard) because ``tga`` is a stop codon
        and ``tgg`` codes for ``w``.  However, using transl_table = 2
        (vertebrate mito), codon ``tgr`` will translate to ``w``
        because both ``tga`` and ``tgg`` code for ``w``::

            >>> gc = GeneticCode(transl_table=1)
            >>> gc.translate('tgr')
                #   codon 'tgr' translates to ['*', 'w'] -- ambiguous -- returning 'x'
            >>> gc = GeneticCode(transl_table=2)
            >>> gc.translate('tgr')
                #   codon 'tgr' translates to 'w'

        Exceptions to the latter rule are codons that ambiguously code
        for either ``d`` or ``n``, which return ambiguous amino acid
        ``b``, and codons that ambigously code for either ``q`` or
        ``e``, which return ambiguous amino acid ``z``.  See the
        example below.
        
        If arg *verbose* is 0, it does not speak (except for errors,
        of course).  If its 1, it speaks for ambiguous translations.
        If its 2, it speaks for all translations.  The default is 1::

            gc = GeneticCode(transl_table=1)

            for cdn in 'gga ggy ray ggn aam sar ccm'.split():
                gc.translate(cdn, verbose=2)

        prints::
        
            codon 'gga' translates to 'g'
            codon 'ggy' translates to 'g'
            codon 'ray' translates to ['d', 'n'] -- ambiguous aa 'b'
            codon 'ggn' translates to 'g'
            codon 'aam' translates to ['k', 'n'] -- ambiguous -- returning 'x'
            codon 'sar' translates to ['q', 'e'] -- ambiguous aa 'z'
            codon 'ccm' translates to 'p'

        """

        gm = ["GeneticCode.translate(), for codon '%s'" % theCodon]
        #print theCodon
        assert len(theCodon) == 3
        if verbose not in [0, 1, 2]:
            gm.append("Arg verbose should be one of 0, 1, or 2.")
            raise Glitch, gm
        
        dnaEquateKeys =  ['b', 'd', 'h', 'k', 'm', 'n', 's', 'r', 'w', 'v', 'y']
        equates = {'b': 'cgt', 'd': 'agt', 'h': 'act', 'k': 'gt', 'm': 'ac',
        'n': 'acgt', 's': 'cg', 'r': 'ag', 'w': 'at', 'v': 'acg', 'y': 'ct'}
        # {'x': 'arndcqeghilkmfpstwyv', 'b': 'dn', 'z': 'eq'}

        # If its in self.code, its easy...
        ret = self.code.get(theCodon)
        if ret:
            if verbose >= 2:
                print "    codon '%s' translates to '%s'" % (theCodon, ret)
            return ret

        # Ok, so not easy.  Check for valid characters.
        if theCodon[0] not in 'acgtbdhkmnsrwvy':
            gm.append("The first position of codon '%s' is not a lowercase DNA character." % theCodon)
            raise Glitch, gm
        if theCodon[1] not in 'acgtbdhkmnsrwvy':
            gm.append("The second position of codon '%s' is not a lowercase DNA character." % theCodon)
            raise Glitch, gm
        if theCodon[2] not in 'acgtbdhkmnsrwvy':
            gm.append("The third position of codon '%s' is not a lowercase DNA character." % theCodon)
            raise Glitch, gm

        # Expand the ambiguities.  Eg 'gcy' becomes ['gcc', 'gct']
        expanded = []
        expanded2 = []
        dnaEquateKeys =  'bdhkmnsrwvy'
        equates = {'b': 'cgt', 'd': 'agt', 'h': 'act', 'k': 'gt', 'm': 'ac',
        'n': 'acgt', 's': 'cg', 'r': 'ag', 'w': 'at', 'v': 'acg', 'y': 'ct'}
        bSet = set(['d', 'n'])
        zSet = set(['e', 'q'])

        if theCodon[0] in dnaEquateKeys:
            c = theCodon[0]
            vv = equates[c]
            for v in vv:
                expanded.append("%s%s%s" % (v, theCodon[1], theCodon[2]))
        else:
            expanded = [theCodon]
        if theCodon[1] in dnaEquateKeys:
            c = theCodon[1]
            vv = equates[c]
            for cd in expanded:
                for v in vv:
                    expanded2.append("%s%s%s" % (cd[0], v, cd[2]))
            expanded = expanded2
            expanded2 = []
        if theCodon[2] in dnaEquateKeys:
            c = theCodon[2]
            vv = equates[c]
            for cd in expanded:
                for v in vv:
                    expanded2.append("%s%s%s" % (cd[0], cd[1], v))
            expanded = expanded2

        #print expanded
        translations = []
        for cd in expanded:
            tr = self.code.get(cd)
            if not tr:
                gm.append("Could not translate expanded codon '%s'" % cd)
                raise Glitch, gm
            translations.append(tr)
        #print translations
        if not translations:
            gm.append("Did not get any translations from expanded codons.")
            raise Glitch, gm
        tSet = set(translations)
        if len(tSet) == 1:
            # Easy again
            if verbose >= 1:
                print "    codon '%s' translates to '%s'" % (theCodon, translations[0])
            return translations[0]
        else:
            tList = list(tSet)

            if len(tSet) == 2:
                if tSet == bSet:
                    if verbose >= 1:
                        print "    codon '%s' translates to %s -- ambiguous aa 'b'" % (theCodon, tList)
                    return 'b'
                elif tSet == zSet:
                    if verbose >= 1:
                        print "    codon '%s' translates to %s -- ambiguous aa 'z'" % (theCodon, tList)
                    return 'z'
                else:
                    if verbose >= 1:
                        print "    codon '%s' translates to %s -- ambiguous -- returning 'x'" % (theCodon, tList)
                    return 'x'
            else: # more than 2
                if verbose:
                    print "    codon '%s' translates to %s -- ambiguous -- returning 'x'" % (theCodon, tList)
                return 'x'

