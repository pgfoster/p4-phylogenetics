from Glitch import Glitch
from Tree import Tree
import Nexus
from Var import var
import sys

import os,string,cStringIO,copy

class TreeFileLite(object):
    """Get trees in big files without reading the lot into memory.

    P4 Tree objects are a little obese, and large tree files will
    flood your RAM.  This class addresses that problem by reading in
    the file as text, and then creating Tree objects only on demand.
    If the trees are not saved then there should not be such a problem
    with memory.

    Instantiate with a file name.  It can handle mcmc output from p4,
    Beast, and MrBayes, and phylip format.

    This can handle tree descriptions with line breaks.  However, it
    does not know about nexus-style 'commenting-out' [ie with square
    brackets, like this].  Also, it is not particularly robust with
    regard to being case-insensitive.  So while the usual way of
    reading in tree files via the read() command will handle nexus
    tree lines that start with tReE or trEe, TreeFileLite cannot, due
    to lazy programming.  So unless your file conforms to the
    expectations of TreeFileLite, it would be best to use read().

    To decrease bloat, it is not loaded by default when you start up
    p4.  To access it, you need to do::

      from p4.TreeFileLite import TreeFileLite

    The only method is getTree(), although you can get the tLines if
    you want.

    Eg to just get a few Tree objects::

      from p4.TreeFileLite import TreeFileLite 
      tfl = TreeFileLite('mcmc_trees_0.nex')
      for i in [23, 45, 67]:
          t = tfl.getTree(i)
          t.draw()

    or, to write some trees, as text (not as Tree objects), to a new
    file::
    
      from p4.TreeFileLite import TreeFileLite
      tfl = TreeFileLite('myBigFile.nex')
      f = file('mySmallerFile.nex', 'w')
      f.write(tfl.header)
      for i in range(24000,25000):
          f.write('tree %s\\n' % tfl.tLines[i])
      f.write('end;')
      f.close()


    """

    def __init__(self, fName=None, verbose=1):
        

        gm = ["TreeFileLite()  init"]
        self.fName = fName
        self.verbose = verbose
        self.translationHash = None
        self.tLines = []
        self.header = None
        self._readTreeFile()
        #self._readMrBayesFile()
        
        self.nSamples = len(self.tLines)
        if self.nSamples:
            if self.verbose >= 1:
                print "Got %i samples." % self.nSamples
        else:
            gm.append("Got 0 tree samples.")
            raise Glitch, gm


    def getTree(self, treeNum):
        savedDoFastNextTok = var.nexus_doFastNextTok
        var.nexus_doFastNextTok = False
        tLine = self.tLines[treeNum]
        if self.verbose >= 3:
            print tLine
        f = cStringIO.StringIO(tLine)
        t = Tree()
        if tLine.startswith("("):
            t.parseNewick(f, translationHash=self.translationHash)
            t.setPreAndPostOrder()
        else:
            t.parseNexus(f, translationHash=self.translationHash)
        var.nexus_doFastNextTok = savedDoFastNextTok
        return t
    
    def _readTreeFile(self):
        gm = ["TreeFileLite._readTreeFile()"]
        # Read in the trees
        try:
            f = file(self.fName, "U")
        except IOError:
            gm.append("Can't find tree file '%s'" % self.fName)
            raise Glitch, gm
        fLines = f.readlines()
        f.close()

        # We cant use doFastNextTok
        savedDoFastNextTok = var.nexus_doFastNextTok
        var.nexus_doFastNextTok = False

        # If it is not a nexus file, it must be a phylip file, so we
        # are done.
        lNum = 0
        aLine = fLines[0].strip()
        if not aLine.startswith("#"):
            self.tLines = fLines
            return

        # So assume it is nexus.  Get the 'header', which might be
        # useful.  Its everything up to the first tree line.
        headerLines = []
        lNum = 0
        aLine = fLines[0]
        sLine = aLine.lstrip()
        lowLine = string.lower(sLine)
        while 1:
            if lowLine.startswith("tree"):
                break
            headerLines.append(aLine)
            lNum += 1
            try:
                aLine = fLines[lNum]
                sLine = aLine.lstrip()
                lowLine = string.lower(sLine)
            except IndexError:
                headerLines = [] # something went wrong ...
                break

        if headerLines:
            self.header = ''.join(headerLines)
        
        
        
        # Get the translate command, if it exists
        translateLines = []
        lNum = 0
        aLine = fLines[0].strip()
        lowLine = string.lower(aLine)
        #print "a aLine: %s" % aLine
        try:
            while not lowLine.startswith("translate"):
                lNum += 1
                aLine = fLines[lNum].strip()
                lowLine = string.lower(aLine)
                if lowLine.startswith('tree'): # then we have gone too far
                    lNum = 0
                    aLine = fLines[0].strip()
                    lowLine = string.lower(aLine)
                    break
        except IndexError:
            # no translate line, so go back to the beginning
            lNum = 0
            aLine = fLines[0].strip()
            lowLine = string.lower(aLine)

        #print "b lowLine: %s" % lowLine
        
        # If we got a translate line, then parse the translate command.
        assert lowLine
        if lowLine.startswith("translate"):
            lNum += 1
            aLine = fLines[lNum].strip()
            while not aLine.endswith(";"):
                translateLines.append(aLine)
                lNum += 1
                aLine = fLines[lNum].strip()
            translateLines.append(aLine)
            translateFlob = cStringIO.StringIO(' '.join(translateLines))
            nx = Nexus.Nexus()
            self.translationHash = nx.readTranslateCommand(translateFlob)
            #print self.translationHash
        var.nexus_doFastNextTok = savedDoFastNextTok

        while not aLine.startswith("tree ") and not aLine.startswith("TREE "):
            lNum += 1
            aLine = fLines[lNum].strip()

        # Get the tree lines.
        self.tLines = []
        while 1:
            if aLine.startswith("tree ") or aLine.startswith("TREE "):
                tempLine = aLine
                # accommodate trees with line breaks.
                while aLine.find(";") < 0:
                    lNum += 1
                    aLine = fLines[lNum].strip()
                    tempLine += aLine
                self.tLines.append(tempLine[5:])
            lNum += 1
            aLine = fLines[lNum].strip()
            if aLine.startswith("end;") or aLine.startswith("End;") or aLine.startswith("ENDBLOCK;") or aLine.startswith('END'):
                break
        


