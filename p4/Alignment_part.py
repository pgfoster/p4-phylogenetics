import string
import pf
from Glitch import Glitch
from Var import var
from Part import Part

def _initParts(self):
    gm = ['Alignment._initParts()']

    if len(self.parts):
        for p in self.parts:
            del(p)
    self.parts = []
    if self.equates:
        eqSymb = self.equates.keys()
        eqSymb.sort()
        eqSymb = string.join(eqSymb, '')
    else:
        eqSymb = ''

    if len(self.sequences) and self.length and  self.symbols and self.dim:
        pass
    else:
        gm.append("Can't allocate part.")
        if not len(self.sequences):
            gm.append("-no sequences.")
        elif not self.length:
            gm.append("-the sequences have no length")
        elif not self.symbols:
            gm.append("-no symbols")
        elif not self.dim:
            gm.append("-dim not set")
        raise Glitch, gm

    if not self.nexusSets or not self.nexusSets.charPartition: # its all one part
        aPart = Part()
        aPart.alignment = self
        aPart.name = 'all'
        aPart.lowName = 'all'
        aPart.dataType = self.dataType
        aPart.dim = self.dim
        aPart.symbols = self.symbols
        aPart.equates = self.equates
        aPart.nTax = len(self.sequences)
        aPart.nChar = self.length
        assert aPart.nChar

        if 0:
            print gm[0]
            print "    symbols=%s" % self.symbols

        aPart.cPart = pf.newPart(len(self.sequences), self.length,
                                     eqSymb, self.symbols)
        if not aPart or not aPart.cPart:
            gm.append("Failed to get memory for part.")
            raise Glitch, gm

        # Make the equates table
        verbose = 0
        equatesTable = []
        if verbose:
            print "equates is %s" % self.equates
            print "eqSymb is %s" % eqSymb # the keys
            print "symbols is %s" % self.symbols
        for i in range(len(eqSymb)):
            if verbose: print "%3s: " % eqSymb[i],
            e = self.equates[eqSymb[i]]
            if verbose: print "%8s : " % e,
            for s in self.symbols:
                if s in e:
                    if verbose: print "%1i" % 1,
                    equatesTable.append('1')
                else:
                    if verbose: print "%1i" % 0,
                    equatesTable.append('0')
            if verbose: print ''
        equatesTable = string.join(equatesTable, '')
        if verbose:
            print "\n\nequatesTable:"
            print equatesTable
        pf.pokeEquatesTable(aPart.cPart, equatesTable)

        sList = []
        for s in self.sequences:
            sList.append(s.sequence)
        if 0:
            print gm[0]
            print "sList = %s" % sList
            print "joined = %s" % string.join(sList, '')
        pf.pokeSequences(aPart.cPart, string.join(sList, ''))
        #print "about to makePatterns ..."
        pf.makePatterns(aPart.cPart)
        #print "about to setInvar"
        pf.setGlobalInvarSitesVec(aPart.cPart)

        #pf.dumpPart(aPart.cPart)
        self.parts.append(aPart)

    elif self.nexusSets.charPartition:
        for cpp in self.nexusSets.charPartition.subsets:
            #print "Doing subset '%s', mask: %s" % (cpp.name, cpp.mask)
            #print "About to subsetUsingMask (self length is %i)" % self.length
            b = self.subsetUsingMask(cpp.mask)
            b._initParts()    # This very method, but now there are no charPartitions in b.
            b.parts[0].name = cpp.name
            b.parts[0].lowName = string.lower(cpp.name)
            self.parts.append(b.parts[0])
            b.parts = [] # so we don't try free-ing it twice



def initDataParts(self):
    gm = ['Alignment.initDataParts()']

    if len(self.parts):
        for p in self.parts:
            del(p)
    self.parts = []

    if len(self.sequences) and self.length and  self.symbols and self.dim:
        pass
    else:
        gm.append("Can't allocate part.")
        if not len(self.sequences):
            gm.append("-no sequences.")
        elif not self.length:
            gm.append("-the sequences have no length")
        elif not self.symbols:
            gm.append("-no symbols")
        elif not self.dim:
            gm.append("-dim not set")
        raise Glitch, gm


    if not self.nexusSets or not self.nexusSets.charPartition: # its all one part
        aPart = DataPart(self)
        self.parts.append(aPart)

    elif self.nexusSets.charPartition:
        for cpp in self.nexusSets.charPartition.subsets:
            #print "Doing subset '%s', mask: %s" % (cpp.name, cpp.mask)
            #print "About to subsetUsingMask (self length is %i)" % self.length
            b = self.subsetUsingMask(cpp.mask)
            b.initDataParts()    # This very method, but now there are no charPartitions in b.
            b.parts[0].alignment  = self
            b.parts[0].name = cpp.name
            b.parts[0].lowName = string.lower(cpp.name)
            self.parts.append(b.parts[0])
            b.parts = [] # so we don't try free-ing the new part twice




def resetSequencesFromParts(self):
    """Gets the sequences from Part.cPart, and installs them in self."""

    #print "Alignment.resetSequencesFromParts() here."
    if (not self.parts) or len(self.parts) == 0:
        gm = ["Alignment.resetSequencesFromParts()"]
        gm.append("No parts.")
        raise Glitch, gm

    if not var.doDataPart:
        if len(self.parts) == 1 and self.parts[0].name == 'all':
            allSeq = pf.symbolSequences(self.parts[0].cPart)
            #print "allSeq[0:20] = %s" % allSeq[0:20]
            for i in range(len(self.sequences)):
                self.sequences[i].sequence = allSeq[(i * self.length): ((i + 1) * self.length)]
        else:
            for i in range(len(self.sequences)):
                self.sequences[i].sequence = list(self.sequences[i].sequence)
            for i in range(len(self.parts)):
                partSeq = pf.symbolSequences(self.parts[i].cPart)
                #print partSeq
                spot = 0
                m = self.nexusSets.charPartition.subsets[i].mask
                for s in self.sequences:
                    for k in range(self.length):
                        if m[k] == '1':
                            s.sequence[k] = partSeq[spot]
                            spot += 1
            for i in range(len(self.sequences)):
                self.sequences[i].sequence = string.join(self.sequences[i].sequence, '')
    else:
        if len(self.parts) == 1:
            for i in range(len(self.sequences)):
                self.sequences[i].sequence = self.parts[0].sequenceString(i)
        else:
            for i in range(len(self.sequences)):
                self.sequences[i].sequence = list(self.sequences[i].sequence)
            for pNum in range(len(self.parts)):
                for sNum in range(len(self.sequences)):
                    partSeq = self.parts[pNum].sequenceString(sNum)
                    print partSeq
                    spot = 0
                    m = self.nexusSets.charPartition.subsets[pNum].mask
                    s = self.sequences[sNum]
                    for k in range(self.length):
                        if m[k] == '1':
                            s.sequence[k] = partSeq[spot]
                            spot += 1
            for i in range(len(self.sequences)):
                self.sequences[i].sequence = string.join(self.sequences[i].sequence, '')



def resetPartsContentFromSequences(self):
    """Reset Part.cPart sequences from self.sequences.

    It then makes patterns, and sets the global invariant sites
    array.  """

    gm = ['Alignment.resetPartsContentFromSequences()']
    if len(self.parts) == 1: # its all one part
        aPart = self.parts[0]
        if not var.doDataPart:
            sList = []
            for s in self.sequences:
                sList.append(s.sequence)
            pf.pokeSequences(aPart.cPart, string.join(sList, ''))
            # are the following necessary?
            pf.makePatterns(aPart.cPart)
            pf.setGlobalInvarSitesVec(aPart.cPart)
        else:
            for sNum in range(len(self.sequences)):
                s = self.sequences[sNum]
                for cNum in range(self.length):
                    theChar = s.sequence[cNum]
                    if theChar == '-':
                        aPart.seq[sNum, cNum] = var.GAP_CODE
                    if theChar == '?':
                        aPart.seq[sNum, cNum] = var.QMARK_CODE
                    elif aPart.equateSymbols and theChar in aPart.equateSymbols:
                        aPart.seq[sNum, cNum] = var.EQUATES_BASE + aPart.equateSymbols.index(theChar)
                    else:
                        aPart.seq[sNum, cNum] = aPart.symbols.index(theChar)

    elif len(self.parts) > 1:
        # the number of parts is also the length of the subsets list
        if self.nexusSets and self.nexusSets.charPartition and \
               self.nexusSets.charPartition.subsets and \
               len(self.nexusSets.charPartition.subsets) == len(self.parts):
            pass
        else:
            gm.append('Something is wrong with the nexusSets or its charPartition')
            raise Glitch, gm
        for i in range(len(self.parts)):
            cpSubset = self.nexusSets.charPartition.subsets[i]
            aPart = self.parts[i]
            b = self.subsetUsingMask(cpSubset.mask)
            if not var.doDataPart:
                sList = []
                for s in b.sequences:
                    sList.append(s.sequence)
                pf.pokeSequences(aPart.cPart, string.join(sList, ''))
                # are the following necessary?
                pf.makePatterns(aPart.cPart)
                pf.setGlobalInvarSitesVec(aPart.cPart)
            else:
                for sNum in range(len(b.sequences)):
                    s = b.sequences[sNum]
                    for cNum in range(b.length):
                        theChar = s.sequence[cNum]
                        if theChar == '-':
                            aPart.seq[sNum, cNum] = var.GAP_CODE
                        if theChar == '?':
                            aPart.seq[sNum, cNum] = var.QMARK_CODE
                        elif aPart.equateSymbols and theChar in aPart.equateSymbols:
                            aPart.seq[sNum, cNum] = var.EQUATES_BASE + aPart.equateSymbols.index(theChar)
                        else:
                            aPart.seq[sNum, cNum] = aPart.symbols.index(theChar)

    else:
        gm.append("No parts.")
        raise Glitch, gm
