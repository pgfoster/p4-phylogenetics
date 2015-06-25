import os,sys,string,array,types
import copy
from Var import var
# Don't bother with NexusToken2, cuz sets blocks are small
from NexusToken import nexusSkipPastNextSemiColon,safeNextTok 
import func
from Glitch import Glitch


##    [Examples from the paup manual,
##    but note the bad charpartition subset names '1' and '2'.  P4 would not allow those names.]
##    charset coding = 2-457 660-896;
##    charset noncoding = 1 458-659 897-898;
##    charpartition gfunc = 1:coding, 2:noncoding;

# Notes from MadSwofMad97.
# TaxSet taxset-name [({Standard | Vector})] = taxon-set;  # standard is default
# TaxPartition partition-name [([{[No]Tokens}]       # tokens is default
#                             [{standard|vector}])]  # standard is default
#                              = subset-name:taxon-set [, subset-name:taxon-set...];
# eg TaxSet outgroup=1-4;
#    TaxSet beetles=Omma-.;
#
# taxpartition populations=1:1-3, 2:4-6, 3:7 8;  # note bad taxpartition names 1, 2, 3
# taxpartition populations (vector notokens) = 11122233;
# 

class CaseInsensitiveDict(dict):
    """A dictionary that is case insensitive, for Nexus"""

    def __init__(self, default=None):
        dict.__init__(self)
        self.default = default
        #self.keyDict = {}


    def __setitem__(self, key, val):
        if type(key) != types.StringType:
            gm = ["CaseInsensitiveDict()"]
            gm.append("The key must be a string.  Got '%s'" % key)
            raise Glitch, gm
        lowKey = string.lower(key)
        dict.__setitem__(self, lowKey, val)
        #self.keyDict[string.lower(key)] = key

    def __getitem__(self, key):
        if type(key) != types.StringType:
            gm = ["CaseInsensitiveDict()"]
            gm.append("The key must be a string.  Got '%s'" % key)
            raise Glitch, gm
        lowKey = string.lower(key)
        try:
            return dict.__getitem__(self, lowKey)
        except KeyError:
            return self.default

    def get(self, key, *args):
        if not args:
            args = (self.default,)
        return dict.get(self, key, *args)



#########################################################################
# CLASS    NexusSets
#########################################################################


class NexusSets(object):
    """A container for Nexus CharSet, CharPartition, and TaxSet objects.

    When the first Nexus sets block is read, a NexusSets object is
    made and saved as ``var.nexusSets``.  ``CharSet``, ``TaxSet``, and
    ``CharPartition`` objects are placed in it, as they are
    read/created.  TaxPartition commands are not implemented.  Here is
    a simple nexus sets block that only has charsets::

        #nexus

        begin sets;
          charset pos1 = 1-.\\3;
          charset pos2 = 2-.\\3;
          charset pos3 = 3-.\\3;
        end;

    To get the third positions only, you could say::

        read('myAlignment.phy')
        a = var.alignments[0]
        read('mySets.nex')       # the sets block above
        b = a.subsetUsingCharSet('pos3')

    What happens above when the mySets.nex file is read is that a
    NexusSets object is created as ``var.nexusSets`` and populated
    with the three charsets as CharSet objects.  Then when you asked
    for a subset, a copy of that NexusSets object was made and applied
    to the alignment.

    Notice that the length of the alignment is not part of the
    information in the sets block, and so things remain undefined
    in ``var.nexusSets`` until the nexus sets are applied to a
    particular alignment.  One consequence of this somewhat awkward
    system is that the same charsets could then be applied to another
    alignment of a different size::

        read('myAlignment.phy')
        aA = var.alignments[0]
        read('anotherAlignment.nex')
        aB = var.alignments[1]
        read('mySets.nex')       # the sets block above
        bA = aA.subsetUsingCharSet('pos3')
        bB = aB.subsetUsingCharSet('pos3')

    In the above example, ``bA.nexusSets`` and ``bB.nexusSets`` are
    both derived from ``var.nexusSets`` but are independent of it, and
    different from each other.
        
    So when an Alignment (or Tree object) wants to use ``var.nexusSets``, it
    makes a copy of it, and attaches the copy as
    theAlignment.nexusSets or theTree.nexusSets

    Here is another example, including a ``charPartition`` definition::
    
        begin sets;
          charset gene1 = 1-213;
          charset gene2 = 214-497;
          charPartition cpName = gene1:gene1, gene2:gene2;
        end;

    For an alignment, you can then set a **character partition** by ::
    
        a.setCharPartition(cpName)

    Do this *before* you make a Data object, to partition the alignment.

    You can also use charsets to extract subsets, eg via::

        b = a.subsetUsingCharSet(csName)


    Setting a charPartition or asking for a subset will trigger
    applying ``var.nexusSets`` to the alignment, but you can also do
    it explicitly, by::

        myTree.setNexusSets()

    NexusSets knows about predefined 'constant', 'gapped', and
    'remainder' charsets.  It does not know about 'missambig' or
    'uninf' charsets.

    NexusSets can either be in the default standard format or in
    vector format -- you can change them to vector format with the ::

        mySet.vectorize()

    method, and you can change them to standard format with the ::

        mySet.standardize()

    method.  For taxSets, you can use actual tax names (rather than
    numbers or ranges) by invoking the method::

        myTaxSet.setUseTaxNames()

    which sets the attribute 'useTaxNames' to True, and puts the
    taxNames for the taxSet in the ::

        taxSet.taxNames

    list, which might be handy.

    You can see the current state of a NexusSets object using ::

       myNexusSets.dump()

    It can also be written out as a nexus sets block.  If an Alignment object
    has a ``nexusSets`` attribute then if you ask the alignment to write
    itself to a nexus file then the Alignment.nexusSets is also
    written.  If you would rather it not be written, delete it first.
    If you would rather it be written to a separate file, do that
    first and then delete it.

    One nice thing about taxsets is that :meth:`Tree.Tree.tv` and
    :meth:`Tree.Tree.btv` know about them and can display them.

    """
    def __init__(self):

        self.charSets = []
        self.charSetsDict = CaseInsensitiveDict()
        self.charSetLowNames = []
        self.taxSets = []
        self.taxSetsDict = CaseInsensitiveDict()
        self.taxSetLowNames = []
        self.charPartitions = []
        self.charPartitionsDict = CaseInsensitiveDict()
        self.charPartitionLowNames = []
        self.charPartition = None
        #self.alignment = None
        self.aligNChar = None
        self.taxNames = []
        self.nTax = None
        self.predefinedCharSetLowNames = ['constant', 'gapped']

        # The nexus format defines several "predefined" charSets.
        # For all datatypes:
        #      constant
        #      gapped
        #      missambig
        #      remainder
        #      uninf
        # I only have implemented 2-- constant and gapped.  The
        # 'remainder' charSet is handled by p4, but not as a CharSet
        # object, since its content depends on the context.
        
        cS = CharSet(self)
        cS.num = -1
        cS.name = 'constant'
        cS.lowName = 'constant'
        cS.format = 'vector'
        #self.charSets.append(cS)
        self.constant = cS
        self.charSetsDict['constant'] = self.constant

        cS = CharSet(self)
        cS.num = -1
        cS.name = 'gapped'
        cS.lowName = 'gapped'
        cS.format = 'vector'
        #self.charSets.append(cS)
        self.gapped = cS
        self.charSetsDict['gapped'] = self.gapped

    def _continueReadingFromNexusFile(self, flob):
        gm = ['NexusSets._continueReadingFromNexusFile()']
        if hasattr(flob, 'name') and flob.name:
            gm.append("file name %s" % flob.name)
        if 0:
            print gm[0]
            print '    var.nexus_doFastNextTok = %s' % var.nexus_doFastNextTok
        nexusSkipPastNextSemiColon(flob)
        commandName = safeNextTok(flob, gm[0])
        lowCommandName = string.lower(commandName)
        #print 'got lowCommandName = %s' % lowCommandName
        while lowCommandName not in [None, 'end', 'endblock']:
            #print "Got lowCommandName '%s'" % lowCommandName
            if lowCommandName == 'charset':
                self._readCharSetCommand(flob)
            elif lowCommandName == 'charpartition':
                self._readCharPartitionCommand(flob)
            elif lowCommandName == 'taxset':
                self._readTaxSetCommand(flob)
            elif lowCommandName == 'taxpartition':
                print
                print gm[0]
                if len(gm) > 1:
                    print gm[1]
                print "    Sorry-- taxpartition is not implemented."
                nexusSkipPastNextSemiColon(flob)
            else:
                gm.append("Got unrecognized sets block command '%s'" %  commandName)
                raise Glitch, gm
            commandName = safeNextTok(flob, 'NexusSets.continueReadingFromNexusFile()')
            lowCommandName = string.lower(commandName)

    def _readCharSetCommand(self, flob):
        # We have just read 'charset'.  The next thing we expect is the charset name.
        gm = ['NexusSets._readCharSetCommand()']
        if hasattr(flob, 'name') and flob.name:
            gm.append("file name %s" % flob.name)
        name = func.nexusUnquoteName(safeNextTok(flob, 'NexusSets: _readCharSetCommand'))
        #print "readCharSetCommand: got name '%s'" % name
        lowName = string.lower(name)
        if not func.nexusCheckName(lowName):
            gm.append("Bad charSet name '%s'" % name)
            raise Glitch, gm

        # Check for duped names
        if lowName in self.charSetLowNames:
            gm.append("Duplicated charSet name '%s'" % name)
            raise Glitch, gm
        elif lowName in self.predefinedCharSetLowNames:
            gm.append("You cannot use the name '%s' -- it is predefined." % name)
            raise Glitch, gm

        cs = CharSet(self)
        cs.name = name
        cs.lowName = lowName

        cs.readTaxOrCharSetDefinition(flob)
        cs.num = len(self.charSets)
        self.charSets.append(cs)
        self.charSetsDict[name] = cs
        self.charSetLowNames.append(cs.lowName)

    def _readTaxSetCommand(self, flob):
        # We have just read 'taxset'.  The next thing we expect is the taxset name.
        gm = ['NexusSets._readTaxSetCommand()']
        if hasattr(flob, 'name') and flob.name:
            gm.append("file name %s" % flob.name)
        name = func.nexusUnquoteName(safeNextTok(flob, 'NexusSets: readTaxSetCommand'))
        #print "readTaxSetCommand: got name '%s'" % name
        lowName = string.lower(name)
        if not func.nexusCheckName(lowName):
            gm.append("Bad taxSet name '%s'" % name)
            raise Glitch, gm

        # Check for duped names
        if lowName in self.taxSetLowNames:
            gm.append("Duplicated taxSet name '%s'" % name)
            raise Glitch, gm

        ts = TaxSet(self)
        ts.name = name
        ts.lowName = lowName

        ts.readTaxOrCharSetDefinition(flob)
        ts.num = len(self.taxSets)
        self.taxSets.append(ts)
        self.taxSetsDict[name] = ts
        self.taxSetLowNames.append(ts.lowName)

    def _readCharPartitionCommand(self, flob):
        gm = ['NexusSets._readCharPartitionCommand()']
        if hasattr(flob, 'name') and flob.name:
            gm.append("file name %s" % flob.name)
        name = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
        #print "readCharPartitionCommand: got name '%s'" % name
        lowName = string.lower(name)
        if not func.nexusCheckName(lowName):
            gm.append("Bad charPartition name '%s'" % name)

        if lowName in self.charPartitionLowNames:
            gm.append("Duplicated charPartition name '%s'" % name)
            raise Glitch, gm

        cp = CharPartition(self)
        cp.name = name
        cp.lowName = lowName

        cp.readCharPartitionDefinition(flob)
        self.charPartitions.append(cp)
        self.charPartitionsDict[name] = cp
        self.charPartitionLowNames.append(cp.lowName)


    def dump(self):
        print "        NexusSets dump"
        if self.constant:
            print "            Predefined char set 'constant'"
            self.constant.dump()
        if self.gapped:
            print "            Predefined char set 'gapped'"
            self.gapped.dump()
        print "            There are %i non-predefined char sets" % len(self.charSets)
        for cs in self.charSets:
            cs.dump()
        print "            There are %i tax sets" % len(self.taxSets)
        for ts in self.taxSets:
            ts.dump()
        print "            There are %i char partitions" % len(self.charPartitions)
        for cp in self.charPartitions:
            cp.dump()
        if self.charPartition:
            print "            self.charPartition.name is %s" % func.nexusFixNameIfQuotesAreNeeded(self.charPartition.name)
        else:
            print "            There is no self.charPartition"

    def write(self):
        """Write self in Nexus format to stdout."""
        self.writeNexusToOpenFile(sys.stdout)

    def writeNexus(self, fName=None):
        """Write self in Nexus format to stdout or a file."""

        if fName:
            f = file(fName, 'w')
        else:
            f = sys.stdout
        f.write('#nexus\n\n')
        self.writeNexusToOpenFile(f)
        if fName:
            f.close()


    def writeNexusToOpenFile(self, flob):
        """This only writes non-trivial stuff.

        Ie if self has only constant and gapped charsets, then it does
        not write anything."""
        
        if self.charSets or self.charPartitions or self.taxSets:
            flob.write('begin sets;\n')
            for cs in self.charSets:
                cs.writeNexusToOpenFile(flob)
            for cp in self.charPartitions:
                cp.writeNexusToOpenFile(flob)
            for ts in self.taxSets:
                ts.writeNexusToOpenFile(flob)
            flob.write('end;\n\n')

    def newCharSet(self, name, mask=None):
        cs = CharSet(self)
        cs.name = name
        cs.name = name.lower()
        cs.num = len(self.charSets)
        if mask:
            cs.format = 'vector'
            cs.mask = mask
        else:
            pass
        self.charSets.append(cs)
        self.charSetsDict[cs.name] = cs

    def dupeCharSet(self, existingCharSetName, newName):
        theCS = self.charSetsDict.get(existingCharSetName)
        if not theCS:
            raise Glitch, "NexusSets.dupeCharSet() -- can't find char set '%s'" % existingCharSetName
        
        cs = CharSet(self)
        cs.name = newName
        cs.name = newName.lower()
        cs.num = len(self.charSets)
        self.charSets.append(cs)
        self.charSetsDict[cs.name] = cs

        cs.format = theCS.format
        cs.triplets = copy.deepcopy(theCS.triplets) # its a list of lists
        cs.tokens = theCS.tokens[:]
        cs.mask = theCS.mask
        cs.aligNChar = theCS.aligNChar

class TaxOrCharSet(object):
    def __init__(self, theNexusSets):
        self.nexusSets = theNexusSets
        self.num = -1
        self.name = None
        self.lowName = None
        self._format = 'standard'  # or 'vector'  So it should be a property.
        self.triplets = []
        self.tokens = []
        self.mask = None
        self.className = 'TaxOrCharSet'
        self.lowTaxNames = []
        self.taxNames = []
        self.useTaxNames = None # undecided
        

    def _getFormat(self):
        return self._format
    def _setFormat(self, newFormat):
        assert newFormat in ['standard', 'vector']
        self._format = newFormat
    format = property(_getFormat, _setFormat)


    def dump(self):
        print "                   %s %i" % (self.className, self.num)
        print "                                   name: %s" % self.name
        if hasattr(self, 'aligNChar'):
            print "                              aligNChar: %s" % self.aligNChar
        print "                                 format: %s" % self.format
        if hasattr(self, 'useTaxNames'):
            print "                            useTaxNames: %s" % self.useTaxNames
        print "                               triplets: "
        for t in self.triplets:
            print "                                         %s" % t
        if hasattr(self, 'numberTriplets'):
            print "                         numberTriplets: "
            for t in self.numberTriplets:
                print "                                         %s" % t
        print "                                 tokens: %s" % self.tokens
        print "                                   mask: %s" % self.mask
        if self.mask:
            print "                          mask 1s-count: %s" % self.mask.count('1')



    def readTaxOrCharSetDefinition(self, flob):
        gm = ['%s.readTaxSetDefinition()' % self.className]
        if hasattr(flob, 'name') and flob.name:
            gm.append("file name %s" % flob.name)
        tok = safeNextTok(flob, gm[0])
        lowTok = string.lower(tok)
        #print "readTaxSetDefinition: get tok '%s'" % tok
        if lowTok == '=':
            pass
        elif lowTok == '(':
            #['standard', 'vector']:
            tok = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
            lowTok = string.lower(tok)
            if lowTok == 'standard':
                pass
            elif lowTok == 'vector':
                self.format = 'vector'
            else:
                gm.append("Unexpected '%s'" % tok)
                gm.append("(I was expecting either 'standard' or")
                gm.append("'vector' following the parenthesis.)")
                raise Glitch, gm
            tok = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
            if tok == ')':
                pass
            else:
                gm.append("Unexpected '%s'" % tok)
                gm.append("(I was expecting an unparentheis after '%s')" % self.format)
                raise Glitch, gm
            tok = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
            if tok != '=':
                gm.append("Unexpected '%s'" % tok)
                gm.append("I was expecting an '=' after '(%s)'" % self.format)
                raise Glitch, gm
        else:
            gm.append("Unexpected '%s'" % tok)
            raise Glitch, gm

        # Now we are on the other side of the '='
        tok = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
        lowTok = string.lower(tok)
        while lowTok not in [None, ';', 'end', 'endblock']:
            self.tokens.append(tok)
            tok = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
            lowTok = string.lower(tok)

        if self.format == 'vector':
            self.mask = string.join(self.tokens, '')
            self.tokens = []
            for i in range(len(self.mask)):
                if self.mask[i] not in ['0', '1']:
                    gm.append("%s '%s', vector format" % (self.className, self.name))
                    gm.append("The vector must be all zeros or ones.")
                    raise Glitch, gm
            #print self.mask

        # do a once-over sanity check, and convert integer strings to ints
        #print "xx1 self.tokens is now %s" % self.tokens
        for tokNum in range(len(self.tokens)):
            tok = self.tokens[tokNum]
            lowTok = string.lower(tok)
            if lowTok in ['.', 'all', '-', '\\']:
                pass
            elif self.className == 'CharSet' and lowTok in self.nexusSets.charSetLowNames:
                #print "    xx3 %s is an existing charSet" % tok
                pass
            elif self.className == 'CharSet' and lowTok in self.nexusSets.predefinedCharSetLowNames:
                #print "    xx3 %s is a pre-defined charSet" % tok
                pass
            elif self.className == 'TaxSet' and lowTok in self.nexusSets.taxSetLowNames:
                #print "    xx4 %s is an existing taxSet" % tok
                pass
            else:
                #print "    xx5"
                try:
                    intTok = int(tok)
                    self.tokens[tokNum] = intTok
                except ValueError:
                    if self.className == 'TaxSet':
                        pass
                    elif self.className == 'CharSet':
                        gm.append("I don't understand the token '%s'" % tok)
                        raise Glitch, gm
                
        # Now I want to make a list of triplets representing eg 23-87\3
        # first item = 23, second item = 87, third = 3
        # not all will exist for each part of the char definition.
        tokNum = 0
        self.triplets = []
        while tokNum < len(self.tokens):
            tok = self.tokens[tokNum]
            #print "Considering tok[%i]  '%s'" % (tokNum, tok)
            if type(tok) == type('str'):
                lowTok = string.lower(tok)
            else:
                lowTok = None

            if self.className == 'TaxSet' and lowTok in self.nexusSets.taxSetLowNames or \
                   self.className == 'charSet' and lowTok in self.nexusSets.charSetLowNames:
                aTriplet = [tok, None, None]
                self.triplets.append(aTriplet)
                tokNum += 1
                if tokNum < len(self.tokens):
                    if self.tokens[tokNum] == '-':
                        gm.append("%s '%s' definition" % (self.className, self.name))
                        gm.append("An existing tax or char set may not be followed by a '-'")
                        raise Glitch, gm
                    if self.tokens[tokNum] == '\\':
                        gm.append("%s '%s' definition" % (self.className, self.name))
                        gm.append("An existing tax or char set may not be followed by a '\\'")
                        raise Glitch, gm

            elif tok == 'all':
                aTriplet = [tok, None, None]
                self.triplets.append(aTriplet)
                tokNum += 1
                if tokNum < len(self.tokens):
                    if self.tokens[tokNum] == '-':
                        gm.append("%s '%s' definition" % (self.className, self.name))
                        gm.append("Tax or char set 'all' may not be followed by a '-'")
                        raise Glitch, gm
                    if self.tokens[tokNum] == '\\':
                        gm.append("%s '%s' definition" % (self.className, self.name))
                        gm.append("Tax or char set 'all' may not be followed by a '\\'")
                        raise Glitch, gm

            elif tok == '-':
                gm.append("%s '%s' definition" % (self.className, self.name))
                gm.append("Out of place '-'")
                raise Glitch, gm

            elif tok == '\\':
                gm.append("%s '%s' definition" % (self.className, self.name))
                gm.append("Out of place '\\'")
                raise Glitch, gm

            elif tok == '.':
                aTriplet = [tok, None, None]
                self.triplets.append(aTriplet)
                tokNum += 1
                if tokNum < len(self.tokens):
                    if self.tokens[tokNum] == '-':
                        gm.append("%s '%s' definition" % (self.className, self.name))
                        gm.append("Tax or char set '.' may not be followed by a '-'")
                        raise Glitch, gm
                    if self.tokens[tokNum] == '\\':
                        gm.append("%s '%s' definition" % (self.className, self.name))
                        gm.append("Tax or char set '.' may not be followed by a '\\'")
                        raise Glitch, gm

            elif type(tok) == type(1) or type(tok) == type('str'):
                aTriplet = [tok, None, None]
                tokNum += 1
                if tokNum < len(self.tokens):
                    if self.tokens[tokNum] == '-':
                        tokNum += 1
                        if tokNum < len(self.tokens):
                            if type(self.tokens[tokNum]) == type('str'): # maybe '.'
                                aTriplet[1] = self.tokens[tokNum]
                            elif type(self.tokens[tokNum]) == type(1):
                                if type(aTriplet[0]) == type(1):
                                    if self.tokens[tokNum] > aTriplet[0]:
                                        aTriplet[1] = self.tokens[tokNum]
                                    else:
                                        gm.append("%s '%s' definition" % (self.className, self.name))
                                        gm.append("If a range is defined by two numbers,")
                                        #gm.append("(as it appears to be -- %s %s %s)" % (
                                        #    aTriplet[0], aTriplet[1], aTriplet[2]))
                                        gm.append("the second number of a range must be bigger than")
                                        gm.append("the first.")
                                        raise Glitch, gm
                                else:
                                    aTriplet[1] = self.tokens[tokNum]

                            else:
                                raise Glitch, gm

                            tokNum += 1

                            if tokNum < len(self.tokens):
                                if self.tokens[tokNum] == '\\':
                                    tokNum += 1
                                    if tokNum < len(self.tokens):
                                        if type(self.tokens[tokNum]) == type(1):
                                            aTriplet[2] = self.tokens[tokNum]
                                        else:
                                            gm.append("%s '%s' definition" % (self.className, self.name))
                                            gm.append("Step value of a range must be a number")
                                            gm.append("(Got '%s')" % self.tokens[tokNum])
                                            raise Glitch, gm
                                        tokNum += 1
                
                self.triplets.append(aTriplet)
        #print "xxy self.mask = %s" % self.mask
        if not self.triplets and not self.mask:
            gm.append("%s '%s' definition" % (self.className, self.name))
            gm.append("Got no definition (no triplets or mask)")
            raise Glitch, gm

        if 0:
            print gm[0]
            print "    Got self.triplets %s" % self.triplets



    def setMask(self):
        """Set self.mask."""

        gm = ["%s.setMask()  name='%s'" % (self.className, self.name)]

        if self.format == 'vector':
            if self.mask:
                pass
            else:
                gm.append("vector format, but no mask?")
                raise Glitch, gm
        elif self.format == 'standard':
            if 0:
                print gm[0]
                self.dump()
                
            if not len(self.triplets):
                gm.append("standard format, but we have no triplets? - no definition?")
                raise Glitch, gm

            if self.className == 'CharSet':
                thisMaskLen = self.aligNChar
                existingSetNames = self.nexusSets.charSetLowNames
                existingSets = self.nexusSets.charSets
                theTriplets = self.triplets
            elif self.className == 'TaxSet':
                thisMaskLen = self.nexusSets.nTax
                existingSetNames = self.nexusSets.taxSetLowNames
                existingSets = self.nexusSets.taxSets
                theTriplets = self.numberTriplets
            mask = array.array('c', thisMaskLen * '0')

            for aTriplet in theTriplets:
                if 0:
                    print gm[0]
                    print "        '%s' aTriplet=%s" % (self.name, aTriplet)
                first = aTriplet[0]
                second = aTriplet[1]
                third = aTriplet[2]

                lowFirst = None
                lowSecond = None
                if type(first) == type('str'):
                    lowFirst = string.lower(first)
                if type(second) == type('str'):
                    lowSecond = string.lower(second)
                    
                if first and not second: # its a single, or an existing set, not a range
                    if lowFirst:
                        if lowFirst == 'all':
                            for i in range(thisMaskLen):
                                mask[i] = '1'
                        if lowFirst in existingSetNames:
                            for aSet in existingSets:
                                if lowFirst == aSet.lowName:
                                    if not aSet.mask:
                                        aSet.setMask()
                                    for j in range(thisMaskLen):
                                        if aSet.mask[j] == '1':
                                            mask[j] = '1'
                        # Maybe its a predefined charset --- constant or gapped
                        elif self.className == 'CharSet' and lowFirst in self.nexusSets.predefinedCharSetLowNames:
                            aSet = None
                            if lowFirst == 'constant':
                                aSet = self.nexusSets.constant
                            elif lowFirst == 'gapped':
                                aSet = self.nexusSets.gapped
                            assert aSet
                            for j in range(thisMaskLen):
                                if aSet.mask[j] == '1':
                                    mask[j] = '1'
                        else:
                            gm.append("I don't know '%s'" % first)
                            raise Glitch, gm
                                
                    elif first == '.':
                        mask[-1] = '1'
                    elif type(first) == type(1):
                        if first > 0 and first <= thisMaskLen:
                            mask[first - 1] = '1'
                        else:
                            # This will have been checked before.
                            gm.append("Component '%s' is out of range of mask len (%s)" % (first, thisMask))
                            raise Glitch, gm
                elif first and second:
                    # Its a range.
                    start = int(first)
                    if second == '.':
                        fin = len(mask)
                    else:
                        fin = int(second)
                    if third:
                        bystep = int(third)
                        #print "mask len %i, start-1 %i, fin %i, bystep %i" % (len(mask), (start-1), fin, bystep)
                        for spot in range(start - 1, fin, bystep):
                            mask[spot] = '1'
                    else:
                        for spot in range(start - 1, fin):
                            mask[spot] = '1'
                #print "            finished incorporating triplet %s into '%s' mask." % (aTriplet, self.name)
            mask = mask.tostring()
            # print "Got char set '%s' mask '%s'" % (self.name, mask)
            self.mask = mask


    def invertMask(self):
        """Change zeros to ones, and non-zeros to zero."""

        gm = ['%s.invertMask()' % self.className]
        if not self.mask:
            self.dump()
            gm.append("The charset has no mask")
            raise Glitch, gm
        self.mask = list(self.mask)
        for i in range(len(self.mask)):
            if self.mask[i] == '0':
                self.mask[i] = '1'
            else:
                self.mask[i] = '0'
        self.mask = string.join(self.mask, '')

    def write(self):
        """Write self in Nexus format to stdout."""
        self.writeNexusToOpenFile(sys.stdout)

    def writeNexus(self):
        """Write self in Nexus format to stdout."""
        self.writeNexusToOpenFile(sys.stdout)

    def writeNexusToOpenFile(self, flob):
        if self.className == 'CharSet':
            theSetName = 'charSet'
        else:
            theSetName = 'taxSet'
            
        if self.format == 'standard':
            flob.write('  %s %s =' % (theSetName, self.name))
            if self.useTaxNames:
                for tN in self.taxNames:
                    flob.write(" %s" % func.nexusFixNameIfQuotesAreNeeded(tN))
            else:
                #for i in self.tokens:
                #    flob.write(' %s' % i)
                previousTok = None
                for theTok in self.tokens:
                    if type(theTok) == types.StringType:
                        if theTok not in ['-', '\\']:
                            tok = func.nexusFixNameIfQuotesAreNeeded(theTok)
                        else:
                            tok = theTok
                    else:
                        tok = theTok
                    if previousTok != None:
                        # tokens will be either ints or strings
                        previousType = type(previousTok)
                        #print "previousTok = %s, previousType = %s" % (previousTok, previousType)
                        if type(tok) == previousType:    # usually put in a space
                            if tok in ['-'] or previousTok in ['-']: # except in this case
                                flob.write('%s' % tok)
                            else:
                                flob.write(' %s' % tok)
                        else:                          # usually no space
                            if tok in ['-'] or previousTok in ['-']:
                                flob.write('%s' % tok)
                            else:                       # except in this case
                                flob.write(' %s' % tok)
                        previousTok = tok
                        #print "previousTok = %s, previousType = %s" % (previousTok, previousType)

                    else:
                        flob.write(' %s' % tok)
                        previousTok = tok

            flob.write(';\n')
        elif self.format == 'vector':
            flob.write('  %s %s (vector) = ' % (theSetName, self.name))
            flob.write('%s;\n' % self.mask)


    def vectorize(self):
        if self.format == 'vector':
            return
        if not self.mask:
            self.setMask()
        #self.triplets = []
        #self.tokens = []
        self.format = 'vector'

    def standardize(self):
        if self.format == 'standard':
            return
        self.triplets = []
        self.tokens = []
        thisTriplet = []
        for mPos in range(len(self.mask)):
            #print "mPos=%i  mask=%s  thisTriplet=%s" % (mPos, self.mask[mPos], thisTriplet)
            if self.mask[mPos] == '0':
                if thisTriplet:
                    if thisTriplet[0] == mPos:
                        thisTriplet.append(None)
                        thisTriplet.append(None)
                    else:
                        thisTriplet.append(mPos)
                        thisTriplet.append(None)
                    #print "   finished triplet -- %s" % thisTriplet
                    self.triplets.append(thisTriplet)
                    thisTriplet = []
                
            else:
                if thisTriplet:
                    pass
                else:
                    thisTriplet.append(mPos + 1)
                    #print "   started triplet -- %s" % thisTriplet
        if thisTriplet:
            if thisTriplet[0] == len(self.mask):
                thisTriplet.append(None)
                thisTriplet.append(None)
            else:
                thisTriplet.append(mPos + 1)
                thisTriplet.append(None)
            #print "   finished last triplet -- %s" % thisTriplet
            self.triplets.append(thisTriplet)
        #print self.triplets

        for triplet in self.triplets:
            if triplet[1] == None:
                self.tokens.append(triplet[0])
            else:
                self.tokens.append(triplet[0])
                self.tokens.append('-')
                self.tokens.append(triplet[1])
        self.format = 'standard'
        #self.dump()


class CharSet(TaxOrCharSet):
    def __init__(self, theNexusSets):
        TaxOrCharSet.__init__(self, theNexusSets)
        self.className = 'CharSet'
        self.aligNChar = None

    def getNChar(self):
        self.setMask()
        return self.mask.count('1')

    def setAligNChar(self, aligNChar):
        gm =['CharSet.setAligNChar()']
        #print "CharSet name=%s, format=%s, aligNChar=%i" % (self.name, self.format, aligNChar)
        self.aligNChar = aligNChar
        if self.format == 'standard':
            for aTriplet in self.triplets:
                first = aTriplet[0]
                second = aTriplet[1]
                third = aTriplet[2]
                if first and not second: # its a single
                    if type(first) == type(1):
                        if first > 0 and first <= self.aligNChar:
                            pass
                        else:
                            gm.append("Charset '%s' definition" % self.name)
                            gm.append("Charset definition element '%s' is out of range" % first)
                            gm.append("(aligNChar = %i)" % self.aligNChar)
                            raise Glitch, gm
                        pass
                elif first and second:  # its a range
                    try:
                        start = int(first)
                    except ValueError:
                        gm.append("Charset '%s' definition" % self.name)
                        gm.append("Can't parse definition element '%s'" % first)
                        raise Glitch, gm
                    if second == '.':
                        fin = self.aligNChar
                    else:
                        try:
                            fin = int(second)
                        except ValueError:
                            gm.append("Charset '%s' definition" % self.name)
                            gm.append("Can't parse definition element '%s'" % second)
                            raise Glitch, gm
                    if third:
                        try:
                            bystep = int(third)
                        except ValueError:
                            gm.append("Charset '%s' definition" % self.name)
                            gm.append("Can't parse definition element '%s'" % third)
                            raise Glitch, gm
        elif self.format == 'vector':
            #print "charset %s, vector format %s, mask %s" % (self.name, self.format, self.mask)
            if self.mask:
                if len(self.mask) == self.aligNChar:
                    pass
                else:
                    gm.append("len(self.mask) is %i, but aligNChar is %i" % (len(self.mask), self.aligNChar))
                    raise Glitch, gm
        else:
            gm.append("bad format %s" % self.format)
            raise Glitch, gm

class TaxSet(TaxOrCharSet):
    def __init__(self, theNexusSets):
        TaxOrCharSet.__init__(self, theNexusSets)
        self.className = 'TaxSet'
        self.numberTriplets = []



    def setNumberTriplets(self):
        gm = ['TaxSet.setNumberTriplets()']
        if not self.nexusSets.lowTaxNames:
            self.nexusSets.lowTaxNames = [string.lower(txName) for txName in self.nexusSets.taxNames]
        self.numberTriplets = []
        #print "self.triplets = %s" % self.triplets
        
        for tr in self.triplets:
            #print "setNumberTriplets() tr=%s" % tr
            numTr = []
            for itemNum in range(2):
                trItem = tr[itemNum]
                #print " considering '%s'" % trItem
                if trItem == None:
                    numTr.append(trItem)
                elif type(trItem) == type(1):
                    numTr.append(trItem)
                elif trItem == '.':
                    numTr.append(self.nexusSets.nTax)
                else:
                    assert type(trItem) == type('str')
                    lowTrItem = string.lower(trItem)
                    if lowTrItem in self.nexusSets.taxSetLowNames:
                        numTr.append(trItem)
                    else:
                        if lowTrItem not in self.nexusSets.lowTaxNames:
                            gm.append("Triplet %s" % tr)
                            gm.append("'%s' is a string, but not in the taxNames." % trItem)
                            raise Glitch, gm
                        theIndx = self.nexusSets.lowTaxNames.index(lowTrItem)
                        theIndx += 1
                        numTr.append(theIndx)
            trItem = tr[2]
            if trItem == None:
                numTr.append(None)
            else:
                assert type(trItem) == type(1)
                numTr.append(trItem)
            assert len(numTr) == 3
            #print numTr

            first = numTr[0]
            # first might be a pre-existing taxSet name
            if type(first) == type('str'):
                pass
            else:
                second = numTr[1]
                assert type(first) == type(1) and first != 0
                if type(second) == type(1):
                    assert second != 0
                    if second <= first:
                        gm.append("Triplet %s" % tr)
                        gm.append("Triplet expressed as numbers. %s" % numTr)
                        gm.append("This appears to be a range, but the second number")
                        gm.append("is not bigger than the first.")
                        raise Glitch, gm
                    assert second <= self.nexusSets.nTax
                assert first <= self.nexusSets.nTax
                              
            self.numberTriplets.append(numTr)


    def setUseTaxNames(self):
        if self.useTaxNames:
            return
        #if not self.mask:
        #    self.setMask()
        if not self.taxNames:
            for pos in range(len(self.mask)):
                c = self.mask[pos]
                if c == '1':
                    self.taxNames.append(self.nexusSets.taxNames[pos])
        self.useTaxNames = True
        


class CharPartitionSubset(object):
    def __init__(self):
        self.name = None
        self.lowName = None
        self.tokens = []
        self.mask = None
        self.triplets = []

    def dump(self):
        print "                              -- CharPartitionSubset"
        print "                                         name: %s" % func.nexusFixNameIfQuotesAreNeeded(self.name)
        print "                                     triplets: "
        for t in self.triplets:
            print "                                               %s" % t
        print "                                       tokens: %s" % self.tokens
        #for t in self.tokens:
        #    print "                                               %s" % t
        print "                                         mask: %s" % self.mask

    def writeNexusToOpenFile(self, flob):   ##Ignore
        flob.write('%s:' % self.name)
        #print self.tokens
        #for i in self.tokens:
        #    flob.write(' %s' % i)
        previousTok = None
        for i in self.tokens:
            if previousTok != None:
                # tokens will be either ints or strings
                previousType = type(previousTok)
                #print "previousTok = %s, previousType = %s" % (previousTok, previousType)
                if type(i) == previousType:    # put in a space
                    flob.write(' %s' % i)
                else:                          # no space
                    flob.write('%s' % i)
                previousTok = i
            else:
                flob.write(' %s' % i)
                previousTok = i



class CharPartition(object):
    def __init__(self, theNexusSets):
        self.nexusSets = theNexusSets
        self.name = None
        self.lowName = None
        self.tokens = []
        self.subsets = []

    ##Ignore
    def readCharPartitionDefinition(self, flob):
        gm = ['CharPartition.readCharPartitionDefinition()']
        if hasattr(flob, 'name') and flob.name:
            gm.append("file name %s" % flob.name)
        tok = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
        lowTok = string.lower(tok)
        while lowTok != '=':
            if lowTok == '(':
                tok = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
                lowTok = string.lower(tok)
                while lowTok != ')':
                    if lowTok in ['notokens', 'vector']:
                        gm.append("Got charpartition modifier: '%s'" % tok)
                        gm.append("It is not implemented.")
                        gm.append("Only 'tokens' and 'standard' are implemented.")
                        raise Glitch, gm
                    elif lowTok in ['tokens', 'standard']:
                        pass
                    else:
                        gm.append("Got charpartition modifier: '%s'" % tok)
                        gm.append("This is not understood.")
                        gm.append("(Only 'tokens' and 'standard' are implemented.)")
                        raise Glitch, gm
                    tok = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
                    lowTok = string.lower(tok)
            else:
                gm.append("Got unexpected token: '%s'" % tok)
                gm.append("I was expecting either an '=' or something in parentheses.")
                raise Glitch, gm

        tok = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
        lowTok = string.lower(tok)
        while lowTok not in [None, ';', 'end', 'endblock']:
            self.tokens.append(tok)
            tok = func.nexusUnquoteName(safeNextTok(flob, gm[0]))
            lowTok = string.lower(tok)

        #print "readCharPartitionDefinition: tokens %s" % self.tokens



        # Divide into CharPartitionSubset instances
        i = 0
        while i< len(self.tokens):
            aSubset = CharPartitionSubset()
            aSubset.name = self.tokens[i]
            if not func.nexusCheckName(aSubset.name):
                gm.append("CharPartition '%s' definition:" % self.name)
                gm.append("Bad subset name (%s, I think)" %  aSubset.name)
                raise Glitch, gm
            aSubset.lowName = string.lower(aSubset.name)
            i += 1
            if i >= len(self.tokens):
                gm.append("CharPartition '%s' definition:" % self.name)
                gm.append("Subset name (%s) should be followed by a colon" %  aSubset.name)
                raise Glitch, gm
            if self.tokens[i] != ':':
                gm.append("CharPartition '%s' definition:" % self.name)
                gm.append("Subset name (%s) should be followed by a colon" %  aSubset.name)
                raise Glitch, gm
            i += 1
            if i >= len(self.tokens):
                gm.append("CharPartition '%s' definition:" % self.name)
                gm.append("Subset name (%s) and colon should be followed" %  aSubset.name)
                gm.append("by a subset definition (charSet or charSet definition)")
                raise Glitch, gm
            while i < len(self.tokens) and self.tokens[i] != ',':
                aSubset.tokens.append(self.tokens[i])
                i += 1
            i += 1
            self.subsets.append(aSubset)

        # do a once-over sanity check,
        # check for duplicated names
        # and convert integer strings to ints
        existingPartNames = []
        for aSubset in self.subsets:
            #print "Checking charPartitionPart '%s'" % aSubset.name
            #print "    existingPartNames '%s'" % existingPartNames
            if aSubset.lowName in existingPartNames:
                gm.append("CharPartition '%s' definition:" % self.name)
                gm.append("Duplicated subset name (%s, I think)" %  aSubset.name)
                raise Glitch, gm
            existingPartNames.append(aSubset.lowName)
            for i in range(len(aSubset.tokens)):
                tok = aSubset.tokens[i]
                lowTok = string.lower(tok)
                #print "considering '%s', ord(lowTok[0])=%i" % (lowTok, ord(lowTok[0]))
                if lowTok in ['.', 'all', '-', '\\', 'remainder']:  # Does not pick up '.'!!!!
                    pass
                elif lowTok in self.nexusSets.charSetLowNames:
                    pass
                elif lowTok in self.nexusSets.predefinedCharSetLowNames:
                    pass
                else:
                    #print "             lowTok=%s, ord(lowTok[0])=%s, ord('.')=%s" % (
                    #    lowTok, ord(lowTok[0]), ord('.'))
                    try:
                        intTok = int(tok)
                        aSubset.tokens[i] = intTok
                    except ValueError:
                        gm.append("CharPartition '%s' definition:" % self.name)
                        gm.append("Can't understand '%s' in subset '%s' definition" % \
                              (tok, aSubset.name))
                        gm.append("(If you are using read('whatever'), and there are backslashes,")
                        gm.append("are you using raw strings, ie read(r'whatever')?)")
                        raise Glitch, gm


    def setSubsetMasks(self):
        """Make charParititionSubset.mask's appropriate to the Alignment.

        This is called by theAlignment.setCharPartition().
        """

        gm = ['CharPartition.setSubsetMasks()']

        assert self.nexusSets.aligNChar
        

        # Make a list of triplets representing eg 23-87\3
        # first item = 23, second item = 87, third = 3
        # Not all will exist for each part of the char definition.
        for aSubset in self.subsets:
            i = 0
            aSubset.triplets = []
            while i < len(aSubset.tokens):
                tok = aSubset.tokens[i]
                if type(tok) == type('string'):
                    lowTok = string.lower(tok)
                else:
                    lowTok = None
                #print "Doing triplets: looking at tok '%s'" % tok
                if lowTok and lowTok in self.nexusSets.charSetLowNames or \
                       lowTok in self.nexusSets.predefinedCharSetLowNames:
                    aTriplet = [lowTok, None, None]
                    aSubset.triplets.append(aTriplet)
                    i += 1
                    if i < len(aSubset.tokens):
                        if aSubset.tokens[i] == '-':
                            gm.append("CharPartition '%s' definition" % self.name)
                            gm.append("Subset '%s' definition" % aSubset.name)
                            gm.append("An existing char set may not be followed by a '-'")
                            raise Glitch, gm
                        if aSubset.tokens[i] == '\\':
                            gm.append("CharPartition '%s' definition" % self.name)
                            gm.append("Subset '%s' definition" % aSubset.name)
                            gm.append("An existing char set may not be followed by a '\\'")
                            raise Glitch, gm

                elif lowTok in ['all', 'remainder']:
                    aTriplet = [lowTok, None, None]
                    aSubset.triplets.append(aTriplet)
                    i += 1
                    if lowTok == 'remainder' and i < len(aSubset.tokens):
                        gm.append("CharPartition '%s' definition" % self.name)
                        gm.append("Subset '%s' definition" % aSubset.name)
                        gm.append("Char set 'remainder' must be the last one in the charPartition definition")
                        raise Glitch, gm

                    if i < len(aSubset.tokens):
                        if aSubset.tokens[i] == '-':
                            gm.append("CharPartition '%s' definition" % self.name)
                            gm.append("Subset '%s' definition" % aSubset.name)
                            gm.append("Char set '%s' may not be followed by a '-'" % lowTok)
                            raise Glitch, gm
                        if aSubset.tokens[i] == '\\':
                            gm.append("CharPartition '%s' definition" % self.name)
                            gm.append("Subset '%s' definition" % aSubset.name)
                            gm.append("Char set '%s' may not be followed by a '\\'" % lowTok)
                            raise Glitch, gm
                elif tok == '-':
                    gm.append("CharPartition '%s' definition" % self.name)
                    gm.append("Subset '%s' definition" % aSubset.name)
                    gm.append("Out of place '-'")
                    raise Glitch, gm

                elif tok == '\\':
                    gm.append("CharPartition '%s' definition" % self.name)
                    gm.append("Subset '%s' definition" % aSubset.name)
                    gm.append("Out of place '\\'")
                    raise Glitch, gm

                elif tok == '.':
                    aTriplet = [tok, None, None]
                    aSubset.triplets.append(aTriplet)
                    i += 1
                    if i < len(aSubset.tokens):
                        if aSubset.tokens[i] == '-':
                            gm.append("CharPartition '%s' definition" % self.name)
                            gm.append("Subset '%s' definition" % aSubset.name)
                            gm.append("Char set '.' may not be followed by a '-'")
                            raise Glitch, gm
                        if aSubset.tokens[i] == '\\':
                            gm.append("CharPartition '%s' definition" % self.name)
                            gm.append("Subset '%s' definition" % aSubset.name)
                            gm.append("Char set '.' may not be followed by a '\\'")
                            raise Glitch, gm

                elif type(tok) == type(1):
                    aTriplet = [tok, None, None]
                    i = i + 1
                    if i < len(aSubset.tokens):
                        if aSubset.tokens[i] == '-':
                            i = i + 1
                            if i < len(aSubset.tokens):
                                if aSubset.tokens[i] == '.':
                                    aTriplet[1] = aSubset.tokens[i]
                                elif type(aSubset.tokens[i]) == type(1):
                                    if aSubset.tokens[i] > aTriplet[0]:
                                        aTriplet[1] = aSubset.tokens[i]
                                    else:
                                        gm.append("CharPartition '%s' definition" % self.name)
                                        gm.append("Subset '%s' definition" % aSubset.name)
                                        gm.append("Second number of a character range must be bigger than")
                                        gm.append("the first.")
                                        raise Glitch, gm

                                else:
                                    gm.append("CharPartition '%s' definition" % self.name)
                                    gm.append("Subset '%s' definition" % aSubset.name)
                                    gm.append("Second item of a character range must be either a")
                                    gm.append("number or a '.'.  I got '%s'" % aSubset.tokens[i])
                                    raise Glitch, gm

                                i = i + 1
                                if i < len(aSubset.tokens):
                                    if aSubset.tokens[i] == '\\':
                                        i = i + 1
                                        if i < len(aSubset.tokens):
                                            if type(aSubset.tokens[i]) == type(1):
                                                aTriplet[2] = aSubset.tokens[i]
                                            else:
                                                gm.append("CharPartition '%s' definition" % self.name)
                                                gm.append("Subset '%s' definition" % aSubset.name)
                                                gm.append("Step value of a range must be a number")
                                                gm.append("(Got '%s')" % aSubset.tokens[i])
                                                raise Glitch, gm

                                            i = i + 1
                    aSubset.triplets.append(aTriplet)
                else:
                    gm.append("CharPartition '%s' definition" % self.name)
                    gm.append("Subset '%s' definition" % aSubset.name)
                    gm.append("token '%s' is not understood." % tok)
                    raise Glitch, gm

            if 0:
                print gm[0]
                print "Got aSubset (%s) triplets %s" % (aSubset.name, aSubset.triplets)
                #sys.exit()


            aSubset.mask = array.array('c', self.nexusSets.aligNChar * '0')

            for aTriplet in aSubset.triplets:
                #print "setSubsetMasks()  Looking at triplet '%s'" % aTriplet
                first = aTriplet[0]
                second = aTriplet[1]
                third = aTriplet[2]
                lowFirst = None
                lowSecond = None
                if type(first) == type('str'):
                    lowFirst = string.lower(first)
                if type(second) == type('str'):
                    lowSecond = string.lower(second)
                    
                if first and not second: # its a single
                    #print "Got single: %s" % first
                    if lowFirst == 'all':
                        for i in range(self.nexusSets.aligNChar):
                            aSubset.mask[i] = '1'
                    elif lowFirst in self.nexusSets.predefinedCharSetLowNames:
                        theCS = None
                        if lowFirst == 'constant':
                            theCS = self.nexusSets.constant
                        elif lowFirst == 'gapped':
                            theCS = self.nexusSets.gapped
                        assert theCS
                        assert theCS.mask
                        for j in range(self.nexusSets.aligNChar):
                            if theCS.mask[j] == '1':
                                aSubset.mask[j] = '1'
                    elif lowFirst in self.nexusSets.charSetLowNames:
                        theCS = None
                        for cs in self.nexusSets.charSets:
                            if lowFirst == cs.lowName:
                                theCS = cs
                                break
                        assert theCS
                        assert theCS.mask
                        for j in range(self.nexusSets.aligNChar):
                            if theCS.mask[j] == '1':
                                aSubset.mask[j] = '1'
                    elif first == '.': # Its legit to use this as a single char.
                        aSubset.mask[-1] = '1'
                    elif type(first) == type(1):
                        if first > 0 and first <= self.nexusSets.aligNChar:
                            aSubset.mask[first - 1] = '1'
                        else:
                            gm.append("CharPartition '%s' definition" % self.name)
                            gm.append("Subset '%s' definition" % aSubset.name)
                            gm.append("Charset definition element '%s' is out of range" % first)
                            gm.append("(aligNChar = %i)" % self.nexusSets.aligNChar)
                            raise Glitch, gm
                    elif lowFirst == 'remainder':
                        #print "Got first == remainder"
                        for i in range(self.nexusSets.aligNChar):
                            aSubset.mask[i] = '1'
                        #print "Got new aSubset.mask = %s" % aSubset.mask
                        for ss in self.subsets[:-1]:
                            if ss.mask:
                                #print "Previous mask: %s" % ss.mask
                                for j in range(self.nexusSets.aligNChar):
                                    if ss.mask[j] == '1':
                                        aSubset.mask[j] = '0'
                            else:
                                gm.append("CharPartition '%s' definition" % self.name)
                                gm.append("Subset '%s' definition" % aSubset.name)
                                gm.append("When implementing 'remainder' charset")
                                gm.append("Found that subset '%s' had no mask" % ss)
                                raise Glitch, gm
                    else:
                        gm.append("CharPartition '%s' definition" % self.name)
                        gm.append("Subset '%s' definition" % aSubset.name)
                        gm.append("Charset definition element '%s' is not understood" % first)
                        raise Glitch, gm

                elif first and second:  # its a range
                    try:
                        start = int(first)
                    except ValueError:
                        gm.append("CharPartition '%s' definition" % self.name)
                        gm.append("Subset '%s' definition" % aSubset.name)
                        gm.append("Can't parse definition element '%s'" % first)
                        raise Glitch, gm
                    if second == '.':
                        fin = len(aSubset.mask)
                    else:
                        try:
                            fin = int(second)
                        except ValueError:
                            gm.append("CharPartition '%s' definition" % self.name)
                            gm.append("Subset '%s' definition" % aSubset.name)
                            gm.append("Can't parse definition element '%s'" % second)
                            raise Glitch, gm
                    if third:
                        try:
                            bystep = int(third)
                        except ValueError:
                            gm.append("CharPartition '%s' definition" % self.name)
                            gm.append("Subset '%s' definition" % aSubset.name)
                            gm.append("Can't parse definition element '%s'" % third)
                        for spot in range(start - 1, fin, bystep):
                            aSubset.mask[spot] = '1'
                    else:
                        for spot in range(start - 1, fin):
                            aSubset.mask[spot] = '1'
            aSubset.mask = aSubset.mask.tostring()
            #print "Got char subset '%s' mask '%s'" % (aSubset.name, aSubset.mask)
            if aSubset.mask.count('1') == 0:
                gm.append("The mask for charPartitionSubset '%s' is empty." % aSubset.name)
                raise Glitch, gm


    ##Ignore
    def checkForOverlaps(self):
        gm = ['CharParitition.checkForOverlaps()']
        unspanned = 0
        for i in range(self.nexusSets.aligNChar):
            sum = 0
            for aSubset in self.subsets:
                if aSubset.mask[i] == '1':
                    sum += 1
            if sum > 1:
                gm.append("Char partition '%s'" % self.name)
                gm.append("The problem is that there are overlapping subsets in this")
                gm.append("charpartition.  The same position is in more than one subset.")
                gm.append("Zero-based position %i, one-based position %i." % (i, i + 1))
                raise Glitch, gm
            if sum < 1:
                unspanned = 1
        if unspanned:
            gm.append("Char partition '%s'" % self.name)
            gm.append("You should be aware that this partition does not span")
            gm.append("the entire sequence.  Hopefully that is intentional.")


    def dump(self):
        print "                CharPartition:     name: %s" % func.nexusFixNameIfQuotesAreNeeded(self.name)
        print "                                 tokens: %s" % self.tokens #string.join(self.tokens)
        #for t in self.tokens:
        #    print "                                         %s" % t
        print "                      number of subsets: %s" % len(self.subsets)
        for aSubset in self.subsets:
            aSubset.dump()

    ##Ignore
    def writeNexusToOpenFile(self, flob):
        flob.write('  charPartition %s = ' % self.name)
        #print " [ %s subsets ] " % len(self.subsets)
        for aSubset in self.subsets[:-1]:
            aSubset.writeNexusToOpenFile(flob)
            flob.write(', ')
        self.subsets[-1].writeNexusToOpenFile(flob)
        flob.write(';\n')

    def mask(self):
        if not self.nexusSets.aligNChar:
            self.nexusSets.aligNChar = self.theNexusSets.aligNChar
        self.setSubsetMasks()
        import array
        m = array.array('c', self.nexusSets.aligNChar * '0')
        for i in range(self.nexusSets.aligNChar):
            for aSubset in self.subsets:
                if aSubset.mask[i] == '1':
                    m[i] = '1'
        return m.tostring()










