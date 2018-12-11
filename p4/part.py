from p4.p4exceptions import P4Error
import p4.pf as pf
import p4.func
from p4.var import var


class Part:

    def __del__(self, freePart=pf.freePart):
        self.alignment = None
        # print "Part.__del__() here.  cPart=%s" % self.cPart
        if self.cPart:
            # print "Part.__del__()   about to free part %i" % self.cPart
            freePart(self.cPart)
            self.cPart = None

    def __init__(self):
        self.alignment = None
        self.name = None             # not lowercased
        self.lowName = None
        self.dataType = 'standard'
        self.dim = None
        self.symbols = '01'
        self.equates = None
        self.nTax = None
        self.nChar = None
        self.cPart = None  # the pointer to the c-structure
        self.seq = None
        #self.invarVec = None
        #self.invarArray = None

    def dump(self):
        print("        Part dump:")
        print("            alignment = %s" % self.alignment)
        print("            nChar = %s" % self.nChar)
        print("            cPart = %s" % self.cPart)
        print("            name = %s" % self.name)
        if self.cPart:
            pf.dumpPart(self.cPart)

    def composition(self, sequenceNumberList=None):
        """Like Alignment.composition(), but for the part, only."""

        # print "About to start part composition()"
        gm = ['Part: composition()']
        if not sequenceNumberList:
            sequenceNumberList = range(self.nTax)
        else:
            if not isinstance(sequenceNumberList, list):
                gm.append("The sequenceNumberList should be a list, ok?")
                raise P4Error(gm)
            for i in sequenceNumberList:
                if not isinstance(i, int):
                    gm.append("The sequenceNumberList should be integers, ok?")
                    raise P4Error(gm)
                if i < 0 or i > self.nTax - 1:
                    gm.append(
                        "Item '%i' in sequenceNumberList is out of range" % i)
                    raise P4Error(gm)

        if not self.cPart:
            self.alignment._initParts()

        for i in range(self.nTax):
            if i in sequenceNumberList:
                # print "self.cPart = %s, i = %s" % (self.cPart, i)
                pf.pokePartTaxListAtIndex(self.cPart, 1, i)
            else:
                pf.pokePartTaxListAtIndex(self.cPart, 0, i)
        return pf.partComposition(self.cPart)
