import sys
import os
import func
from Var import var
from Glitch import Glitch
from subprocess import Popen,PIPE

class DistanceMatrix:
    """A container for distances between sequences (usually).

    The numbers are in self.matrix, a self.dim * self.dim list of lists.
    
    There is also a self.names attribute, which is usually for sequence names.
    """

    def __init__(self):
        self.dim = None
        self.names = None
        self.matrix = None
        self.message = None

    def setDim(self, dim):
        #print "setDim here"
        self.dim = dim
        self.matrix = []
        for i in range(dim):
            l = [0.0] * dim
            self.matrix.append(l)

    def writeNexus(self, fName=None, writeTaxaBlock=1, append=0, digitsAfterDecimal=6):
        """Write out self in Nexus format.
        
        If writeTaxaBlock=1, then a Taxa block is written before the
        Distances block. Append, if 0, writes #NEXUS first.  If 1,
        does not write #NEXUS.  """

        gm = ["DistanceMatrix.writeNexus()"]
        import string
        assert self.dim, "Distance Matrix.writeNexus() no dim"

        if fName == None or fName == sys.stdout:
            f = sys.stdout
        else:
            if append:
                import os
                if os.path.isfile(fName):
                    try:
                        f = open(fName, 'a')
                    except IOError:
                        gm.append("Can't open %s for appending." % fName)
                        raise Glitch, gm
                else:
                    print gm[0]
                    assert not os.path.lexists()
                    print "    'append' is requested,"
                    print "    but '%s' is not a regular file (maybe it doesn't exist?)." % fName
                    print "    Writing to a new file instead."
                    try:
                        f = open(fName, 'w')
                        f.write('#NEXUS\n\n')
                    except IOError:
                        gm.append("Can't open %s for writing." % fName)
                        raise Glitch, gm

            else:
                try:
                    f = open(fName, 'w')
                except IOError:
                    gm.append("Can't open %s for writing." % fName)
                    raise Glitch, gm
        self.writeNexusToOpenFile(f, writeTaxaBlock, append, digitsAfterDecimal)
        if f != sys.stdout:
            f.close()
        

    def writeNexusToOpenFile(self, flob, writeTaxaBlock, append, digitsAfterDecimal):

        gm = ["DistanceMatrix.writeNexusToOpenFile()"]
        import string
        assert self.dim, "Distance Matrix.writeNexusToOpenFile() no dim"
        f = flob

        if append:
            pass
        else:
            f.write('#NEXUS\n\n')

        if writeTaxaBlock:
            f.write('begin taxa;\n')
            f.write('  dimensions ntax=%s;\n' % self.dim)
            f.write('  taxlabels')
            for i in range(self.dim):
                f.write(' %s' % func.nexusFixNameIfQuotesAreNeeded(self.names[i]))
            f.write(';\n')
            f.write('end;\n\n')


        f.write('begin distances;\n')
        if self.message:
            f.write("\n  [%s\n  ]\n" % self.message)
            
        f.write('  format triangle=both;\n')
        f.write('  matrix\n')

        # Make the format strings.
        if digitsAfterDecimal > 8:
            totWid = digitsAfterDecimal + 5
        else:
            totWid = 11
        #colNameFormat = '%' + '-%is' % totWid
        if digitsAfterDecimal == 0:
            numberFormat = ' %' + '%ii   ' % (totWid - 4)
        elif digitsAfterDecimal > 0:
            numberFormat = '%' + '%i.%if' % (totWid, digitsAfterDecimal)
        else:
            gm.append("digitsAfterDecimal may not be below zero.")
            raise Glitch, gm

        if self.names:
            longestNameLength = 0
            for i in self.names:
                if len(i) > longestNameLength:
                    longestNameLength = len(i)
            nameFormat = '%' + '%is ' % longestNameLength
            f.write('  [')
            f.write(nameFormat % ' ')
            for i in self.names:
                #f.write('%-10s' % i[:9])
                colName = '%s' % i[:totWid - 1]
                f.write(string.center(colName, totWid))
            f.write(']\n')

        for i in range(self.dim):
            if self.names:
                f.write('  ')
                f.write(nameFormat % func.nexusFixNameIfQuotesAreNeeded(self.names[i]))
            for j in range(self.dim):
                #f.write('%10.6f' % self.matrix[i][j])
                f.write(numberFormat % self.matrix[i][j])
            f.write('\n')

        f.write('  ;\n')
        f.write('end;\n')

    def writePhylip(self, fName=None, append=0, digitsAfterDecimal=6):
        """Write out self in Phylip format"""

        gm = ["DistanceMatrix.writePhylip()"]
        assert self.dim, "Distance Matrix.writePhylip() no dim"

        if fName == None or fName == sys.stdout:
            f = sys.stdout
        else:
            if append:
                import os
                if os.path.isfile(fName):
                    try:
                        f = open(fName, 'a')
                    except IOError:
                        gm.append("Can't open %s for appending." % fName)
                        raise Glitch, gm
                else:
                    print gm[0]
                    assert not os.path.lexists()
                    print "    'append' is requested,"
                    print "    but '%s' is not a regular file (maybe it doesn't exist?)." % fName
                    print "    Writing to a new file instead."
                    try:
                        f = open(fName, 'w')
                    except IOError:
                        gm.append("Can't open %s for writing." % fName)
                        raise Glitch, gm

            else:
                try:
                    f = open(fName, 'w')
                except IOError:
                    gm.append("Can't open %s for writing." % fName)
                    raise Glitch, gm
        self.writePhylipToOpenFile(f, digitsAfterDecimal)
        if f != sys.stdout:
            f.close()
        
    def writePhylipToOpenFile(self, flob, digitsAfterDecimal):

        gm = ["DistanceMatrix.writePhylipToOpenFile()"]
        assert self.dim, "Distance Matrix.writePhylipToOpenFile() no dim"
        f = flob

        f.write('%i\n' % self.dim)

        # Make the format strings.
        if digitsAfterDecimal > 8:
            totWid = digitsAfterDecimal + 5
        else:
            totWid = 11
        #colNameFormat = '%' + '-%is' % totWid
        if digitsAfterDecimal == 0:
            numberFormat = ' %' + '%ii   ' % (totWid - 4)
        elif digitsAfterDecimal > 0:
            numberFormat = '%' + '%i.%if' % (totWid, digitsAfterDecimal)
        else:
            gm.append("digitsAfterDecimal may not be below zero.")
            raise Glitch, gm

        longestNameLength = 0
        for i in self.names:
            if len(i) > longestNameLength:
                longestNameLength = len(i)
        nameFormat = '%' + '%is ' % longestNameLength

        for i in range(self.dim):
            if self.names:
                f.write('  ')
                f.write(nameFormat % func.nexusFixNameIfQuotesAreNeeded(self.names[i]))
            for j in range(self.dim):
                #f.write('%10.6f' % self.matrix[i][j])
                f.write(numberFormat % self.matrix[i][j])
            f.write('\n')

        f.write('\n')

    def readPhylipFile(self, theFileName):
        """Read a distance matrix in phylip format.
        """
        gm = ['DistanceMatrix.readPhylipFile()']
        import string
        self.names = []
        f = open(theFileName, 'r')
        splitString = string.split(f.readline())
        if len(splitString) != 1:
            gm.append("The first line should have the number of taxa, and thats all.")
            raise Glitch, gm
        try:
            theDim = int(splitString[0])
        except ValueError:
            gm.append("Could not get an integer from the first line.")
            raise Glitch, gm
        self.setDim(theDim)
        for i in range(theDim):
            aLine = f.readline()
            if not aLine:
                gm.append("File too short?")
                raise Glitch, gm
            splitLine = string.split(aLine)
            if len(splitLine) == 0:
                gm.append("Empty line?")
                gm.append("Got: '%s'" % aLine)
                raise Glitch, gm
            if len(splitLine) <= 1:
                gm.append("Line too short.")
                gm.append("Got: '%s'" % aLine)
                raise Glitch, gm
            self.names.append(splitLine[0])
            #print "got name %s" % splitLine[0]
            j = 0
            for k in splitLine[1:]:
                try:
                    self.matrix[i][j] = float(k)
                    #print "i=%i, j=%i, got dist %f" % (i, j, self.matrix[i][j])
                    j = j + 1
                    if j >= theDim:
                        break
                except ValueError:
                    gm.append("    Could not convert %s to a float." % k)
                    raise Glitch, gm
            while j < theDim:
                aLine = f.readline()
                if not aLine:
                    gm.append("File too short?")
                    raise Glitch, gm
                splitLine = string.split(aLine)
                if len(splitLine) == 0:
                    gm.append("Empty line?")
                    gm.append("Got: '%s'" % aLine)
                    raise Glitch, gm
                for k in splitLine:
                    try:
                        self.matrix[i][j] = float(k)
                        #print "i=%i, j=%i, got dist %f" % (i, j, self.matrix[i][j])
                        j = j + 1
                        if j >= theDim:
                            break
                    except ValueError:
                        gm.append("Could not convert %s to a float." % k)
                        raise Glitch, gm
        f.close()


    def njUsingPaup(self, paupPath='paup'):
        """Use paup to make a neighbor-joining tree, which is returned.

        The resulting tree is read in by p4, and is returned.

        We interact with paup by writing files, but care is taken that
        existing files are not overwritten, because new file names are
        made to be unique. 

        If this does not work well, try setting the paupPath arg.
        """

        gm = ["DistanceMatrix.njUsingPaup()"]

        #filename    = sha.new(str(os.getpid())).hexdigest()[-10:]
        #dmFName     = os.path.join(pathPrefix, "%s.dmat" % filename)
        #treeFName   = os.path.join(pathPrefix, "%s.tree" % filename)
        #pFName      = os.path.join(pathPrefix, "%s.cmds" % filename)

        #tempfile.mkstemp(suffix='', prefix='tmp', dir=None, text=False)
        #if pathPrefix:
        #    theDir = pathPrefix
        #else:
        #    theDir = None
        flob_dm, dmFName_fq = func.uniqueFile('tmp.dm')
        flob_tf, treeFName_fq = func.uniqueFile('tmp.tree') #tempfile.mkstemp(suffix='tree', dir=theDir)
        flob_tf.close()
        flob_pf, pFName = func.uniqueFile('tmp.cmds') # tempfile.mkstemp(suffix='cmds', dir=theDir)

        # Throw the dir and dirname away. 
        dirname, dmFName = os.path.split(dmFName_fq)
        dirname, treeFName = os.path.split(treeFName_fq)

        # Make the paup commands
        paupCommandString = """#nexus
        begin paup;
          execute %s;
          set crit=dist;
          dset negbrlen=setzero;
          nj;
          savetrees file=%s format=altnex brlens=yes taxablk=yes replace=yes;
          quit;
        end;
        
        """ % (dmFName, treeFName)

        #print paupCommandString

        # Write the files, do the analysis
        #writeNexusToOpenFile(self, flob, writeTaxaBlock, append, digitsAfterDecimal)
        self.writeNexusToOpenFile(flob_dm, True, False, 6)
        flob_dm.close()
        flob_pf.write(paupCommandString)
        flob_pf.close()

        os.system("%s -n %s > /dev/null" % (paupPath, pFName))

        # This is the result.  The tree, if it exists, is read in by p4.
        oldLen = len(var.trees)
        func.read(treeFName)
        newLen = len(var.trees)
        if newLen == oldLen + 1:
            pass
        else:
            gm.append("I was expecting exactly one tree.  Got %i" % (oldLen - newLen))
            raise Glitch, gm
        t = var.trees.pop()

        # Tidy up.
        os.remove(treeFName)
        os.remove(pFName)
        os.remove(dmFName)

        for n in t.iterNodesNoRoot():
            if n.br.len < 0.0:
                n.br.len = 0.0

        return t

    def bionj(self):
        """Use bionj to make a neighbor-joining tree, which is returned.

        The resulting tree is read in by p4, and is returned.

        We interact with bionj by writing files, but care is taken that
        existing files are not overwritten, because new file names are
        made to be unique. 

        If the branch lengths are less than zero, they are made to be zero.
        """

        gm = ["DistanceMatrix.bionj()"]

        flob_dm, dmFName_fq = func.uniqueFile('tmp.dm')
        flob_tf, treeFName_fq = func.uniqueFile('tmp.tree') #tempfile.mkstemp(suffix='tree', dir=theDir)
        flob_tf.close()

        # Throw the dir and dirname away. 
        dirname, dmFName = os.path.split(dmFName_fq)
        dirname, treeFName = os.path.split(treeFName_fq)

        # Write the files, do the analysis
        self.writePhylipToOpenFile(flob_dm, 6)
        flob_dm.close()

        os.system("bionj %s %s >/dev/null" % (dmFName, treeFName))

        # This is the result.  The tree, if it exists, is read in by p4.
        oldLen = len(var.trees)
        func.read(treeFName)
        newLen = len(var.trees)
        if newLen == oldLen + 1:
            pass
        else:
            gm.append("I was expecting exactly one tree.  Got %i" % (oldLen - newLen))
            raise Glitch, gm
        t = var.trees.pop()

        # Tidy up.
        os.remove(treeFName)
        os.remove(dmFName)

        for n in t.iterNodesNoRoot():
            if n.br.len < 0.0:
                n.br.len = 0.0

        return t
    
    def fastme(self):
        """Use fastme to make a minimum-evolution tree, which is returned.

        The resulting tree is read in by p4, and is returned.

        We interact with fastme by writing files, but care is taken that
        existing files are not overwritten, because new file names are
        made to be unique. 

        If the branch lengths are less than zero, they are made to be zero.
        """

        gm = ["DistanceMatrix.fastme()"]

        flob_dm, dmFName_fq = func.uniqueFile('tmp.dm')
        flob_tf, treeFName_fq = func.uniqueFile('tmp.tree') #tempfile.mkstemp(suffix='tree', dir=theDir)
        flob_tf.close()

        # Throw the dir and dirname away. 
        dirname, dmFName = os.path.split(dmFName_fq)
        dirname, treeFName = os.path.split(treeFName_fq)

        # Write the files, do the analysis
        self.writePhylipToOpenFile(flob_dm, 6)
        flob_dm.close()

        os.system("fastme -i %s -o %s" % (dmFName, treeFName))

        # This is the result.  The tree, if it exists, is read in by p4.
        oldLen = len(var.trees)
        func.read(treeFName)
        newLen = len(var.trees)
        if newLen == oldLen + 1:
            pass
        else:
            gm.append("I was expecting exactly one tree.  Got %i" % (oldLen - newLen))
            raise Glitch, gm
        t = var.trees.pop()

        # Tidy up.
        os.remove(treeFName)
        os.remove(dmFName)

        for n in t.iterNodesNoRoot():
            if n.br.len < 0.0:
                n.br.len = 0.0
            if not n.isLeaf:
                n.name = None

        return t
    
        
