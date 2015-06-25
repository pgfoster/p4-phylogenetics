from Glitch import Glitch
from Tree import Tree
import Nexus
from Var import var

import os,string,cStringIO,copy

class PosteriorSamples(object):
    """A container for mcmc samples from files.

    This would be useful if you wanted to do eg posterior predictive
    simulations based on the posterior distribution from an already
    completed run.  This reads in the tree and parameter files, and
    puts the two together to make p4 Tree objects with a model
    attached.  For each, you would then attach a data object and do
    your simulation.

    After getting it instantiated, the only method is
    getSample(sampNum), which gets the zero-based sampNum'th sample.

    This will handle output from p4 or MrBayes.  It will only do the
    MrBayes models that p4 can do as well.  It will only do one 'run'
    at a time.  It will do models on partitioned data, and p4
    tree-hetero models.

    **Args**
    

    tree
             A p4 tree, with model attached.  It must also have
             taxNames attached, in their proper order.
    runNum
             The run that you want to use.  It would be zero-based for
             p4 files, and 1-based for MrBayes files.  Eg for the p4
             file mcmc_trees_0.nex you would say runNum=0, and for the
             MrBayes file foo.run2.p you would say runNum=2.
    program
             'p4' or 'mrbayes'
    mbBaseName
             If you are trying to read MrBayes files, you need to
             supply the 'baseName' for the files.  Eg if the run1 tree
             file is foo.run1.t and the run1 parameters file is
             foo.run1.p, you would supply mbBaseName='foo'.
    directory
             The directory where the tree files and parameter files
             are located.  By default, it is '.', which is the current
             directory.
    verbose
             An integer, from 0-3.  Zero is quiet.  The bigger the
             number, the chattier it gets.


    **Example**::

        # We need a tree with an attached model.  Lets say that we have a one
        # partition model, with a GTR+G model.  Do the usual p4 model setup.

        read('myTree.nex')
        t = var.trees[0]

        # Read in some data.
        read('myData.nex')
        d = Data()

        t.data = d
        t.newComp(free=1, spec='empirical')
        t.newRMatrix(free=1, spec='ones')
        t.setNGammaCat(nGammaCat=4)
        t.newGdasrv(free=1, val=0.5)
        t.setPInvar(free=0, val=0.0)

        # Check to make sure its all good to go.
        # (This step is needed for empirical comps, at least)
        t.data = d
        t.calcLogLike()

        # Instantiate
        ps = PosteriorSamples(t, runNum=1, program='mrbayes', mbBaseName='mbout', verbose=2)

        # Iterate over samples, saving your favourite test quantity from simulations as you go.
        myStats = []
        for sampNum in range(500,1000):
            t2 = ps.getSample(sampNum)
            t2.data = d
            t2.simulate()
            myStats.append(t2.data.simpleBigXSquared()[0])
    
    """
    def __init__(self, tree, runNum, program='p4', mbBaseName=None, directory='.', verbose=1):
        

        gm = ["PosteriorSamples()  init"]
        if not isinstance(tree, Tree) or not tree.model:
            gm.append("Instantiate with a p4 tree with a model attached.")
            raise Glitch, gm
        if not tree.taxNames:
            gm.append("The tree should have taxNames (in proper order!) attached.")
            raise Glitch, gm
                      
        self.tree = tree

        # Check the tree, by calculating the likelihood
        #self.tree.calcLogLike(verbose=False)
        self.model = copy.deepcopy(self.tree.model)
        for pNum in range(self.model.nParts):
            for compNum in range(self.model.parts[pNum].nComps):
                if self.model.parts[pNum].comps[compNum].spec == 'empirical':
                    self.model.parts[pNum].comps[compNum].spec = 'specified'
                if not self.model.parts[pNum].comps[compNum].val:
                    gm.append("Comp %i in partition %i has no val set." % (compNum, pNum))
                    gm.append("Maybe fix by calculating a likelihood before?")
                    raise Glitch, gm
        self.model.cModel = None
        self.runNum = int(runNum)
        self.goodPrograms = ['p4', 'mrbayes']
        lowProgram = string.lower(program)
        if program not in self.goodPrograms:
            gm.append("The program generating the files should be one of %s" % self.goodPrograms)
            raise Glitch, gm
        self.program = lowProgram
        self.verbose = verbose
        assert os.path.isdir(directory)
        self.directory = directory
        
        if self.program == 'p4':
            self._readP4Files()
        elif self.program == 'mrbayes':
            self.mbBaseName = mbBaseName
            self._readMrBayesFiles()

        self.nSamples = len(self.tLines)
        if self.tree.model.nFreePrams:
            nPLines = len(self.pLines)
            if self.nSamples and self.nSamples == nPLines:
                if self.verbose >= 1:
                    print "Got %i samples." % self.nSamples
            else:
                gm.append("Got %i tree samples, but %i parameter samples." % (self.nSamples, nPLines))
                raise Glitch, gm
        else:
            #print "Got %i samples. (no free parameters)" % self.nSamples
            pass
        

    def getSample(self, sampNum):
        if self.program == 'p4':
            return self._getP4SampleTree(sampNum)
        elif self.program == 'mrbayes':
            return self._getMrBayesSampleTree(sampNum)

    def _readP4Files(self):
        gm = ["PosteriorSamples._readP4Files()"]
        # Read in the trees
        fName = "mcmc_trees_%i.nex" % self.runNum
        if self.directory != '.':
            fName = os.path.join(self.directory, fName)
        try:
            f = file(fName)
        except IOError:
            gm.append("Can't find tree file '%s'" % fName)
            raise Glitch, gm
        fLines = f.readlines()
        f.close()

        # Get the translate command
        savedDoFastNextTok = var.nexus_doFastNextTok
        var.nexus_doFastNextTok = False
        lNum = 0
        aLine = fLines[0].strip()
        translateLines = []
        while not aLine.startswith("translate"):
            lNum += 1
            aLine = fLines[lNum].strip()
        lNum += 1
        aLine = fLines[lNum].strip()
        while not aLine.startswith(";"):
            translateLines.append(aLine)
            lNum += 1
            aLine = fLines[lNum].strip()
        translateLines.append(aLine)
        translateFlob = cStringIO.StringIO(' '.join(translateLines))
        nx = Nexus.Nexus()
        self.translationHash = nx.readTranslateCommand(translateFlob)
        #print self.translationHash
        var.nexus_doFastNextTok = savedDoFastNextTok

        # Get the models definition, if it exists.  Move to the first tree line.
        while not aLine.startswith("tree t_"):
            if aLine.startswith("[&&p4 models"):
                self.modelLine = aLine
            lNum += 1
            aLine = fLines[lNum].strip()

        # Get the tree lines.
        self.tLines = []
        while not aLine.startswith("end;"):
            self.tLines.append(aLine[5:])
            lNum += 1
            aLine = fLines[lNum].strip()

        if self.tree.model.nFreePrams:
            # Read in the prams
            fName = "mcmc_prams_%i" % self.runNum
            if self.directory != '.':
                fName = os.path.join(self.directory, fName)
            try:
                f = file(fName)
            except IOError:
                gm.append("Can't find prams file '%s'" % fName)
                raise Glitch, gm
            fLines = f.readlines()
            f.close()
            lNum = 0
            aLine = fLines[0].strip()
            while aLine[0] not in string.digits:
                lNum += 1
                aLine = fLines[lNum].strip()
            self.pLines = []
            while aLine:
                self.pLines.append(aLine)
                lNum += 1
                try:
                    aLine = fLines[lNum].strip()
                except IndexError:
                    break

            # Read in the pramsProfile
            self.nPrams = None
            self.pramsProfile = None
            fName = "mcmc_pramsProfile.py"
            if self.directory != '.':
                fName = os.path.join(self.directory, fName)
            try:
                loc = {}
                execfile(fName, {}, loc)
                #loc =locals()  no workee.
                #print "loc = %s" % loc
                self.nPrams = loc['nPrams']
                self.pramsProfile = loc['pramsProfile']
            except IOError:
                print "The file '%s' cannot be found." % fName
            

    def _getP4SampleTree(self, sampNum):
        savedDoFastNextTok = var.nexus_doFastNextTok
        var.nexus_doFastNextTok = False
        tLine = self.tLines[sampNum]
        if self.verbose >= 3:
            print tLine
        f = cStringIO.StringIO(tLine)
        t = Tree()
        t.parseNexus(f, translationHash=self.translationHash, doModelComments=self.tree.model.nParts)
        var.nexus_doFastNextTok = savedDoFastNextTok
        t.taxNames = self.tree.taxNames

        for n in t.iterLeavesNoRoot():
            n.seqNum = t.taxNames.index(n.name)
        
        t.model = copy.deepcopy(self.model)

        if self.tree.model.nFreePrams:
            pLine = self.pLines[sampNum]
            if self.verbose >= 3:
                print pLine
            splitPLine = pLine.split()

            pGenNum = int(splitPLine[0])
            splitTName = t.name.split('_')
            tGenNum = int(splitTName[1])
            if tGenNum != pGenNum:
                raise Glitch, "something wrong. tGenNum=%i, pGenNum=%i" % (tGenNum, pGenNum)
            if self.verbose >= 2:
                print "(zero-based) sample %i is gen %i" % (sampNum, tGenNum)

            #t.model.dump()

            splIndx = 1
            for pNum in range(len(self.pramsProfile)):
                compNum = 0
                rMatrixNum = 0
                gdasrvNum = 0
                for desc in self.pramsProfile[pNum]:
                    if desc[0] == 'relRate':
                        t.model.parts[pNum].relRate = float(splitPLine[splIndx])
                        splIndx += 1                    
                    elif desc[0] == 'comp':
                        vv = []
                        for i in range(desc[1]):
                            vv.append(float(splitPLine[splIndx]))
                            splIndx += 1
                        for i in range(desc[1]):
                            if vv[i] < var.PIVEC_MIN:
                                vv[i] = var.PIVEC_MIN * 1.1
                        thisSum = sum(vv)
                        factor = 1.0 / thisSum # must sum to one
                        for i in range(desc[1]):
                            t.model.parts[pNum].comps[compNum].val[i] = vv[i] * factor
                        compNum += 1
                    elif desc[0] == 'rMatrix':
                        vv = []
                        for i in range(desc[1]):
                            vv.append(float(splitPLine[splIndx]))
                            splIndx += 1
                        if len(vv) == 1:
                            # its a '2p' model, with a kappa
                            t.model.parts[pNum].rMatrices[rMatrixNum].val[0] = vv[0]
                        else:
                            # gtr
                            for i in range(desc[1]):
                                if vv[i] < var.RATE_MIN:
                                    vv[i] = var.RATE_MIN * 1.1
                            thisSum = sum(vv)
                            factor = 1.0 / thisSum # must sum to one
                            for i in range(desc[1]):
                                t.model.parts[pNum].rMatrices[rMatrixNum].val[i] = vv[i] * factor
                        rMatrixNum += 1
                    elif desc[0] == 'gdasrv':
                        t.model.parts[pNum].gdasrvs[gdasrvNum].val[0] = float(splitPLine[splIndx])
                        splIndx += 1
                        gdasrvNum += 1
                    elif desc[0] == 'pInvar':
                        t.model.parts[pNum].pInvar.val = float(splitPLine[splIndx])
                        splIndx += 1

            if splIndx != len(splitPLine):
                raise Glitch, "Something is wrong.  After reading, splIndx=%i, but len split pram line=%i" % (
                    splIndx, len(splitPLine))
        return t
        
        
    def _readMrBayesFiles(self):
        gm = ["PosteriorSamples._readMrBayesFiles()"]
        # Read in the trees
        fName = "%s.run%i.t" % (self.mbBaseName, self.runNum)
        if self.directory != '.':
            fName = os.path.join(self.directory, fName)
        try:
            f = file(fName)
        except IOError:
            gm.append("Can't find tree file '%s'" % fName)
            raise Glitch, gm
        fLines = f.readlines()
        f.close()

        # Get the translate command
        savedDoFastNextTok = var.nexus_doFastNextTok
        var.nexus_doFastNextTok = False
        lNum = 0
        aLine = fLines[0].strip()
        translateLines = []
        while not aLine.startswith("translate"):
            lNum += 1
            aLine = fLines[lNum].strip()
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

        # Get the models definition, if it exists.  Move to the first tree line.
        # MrBayes3.2 uses 'gen', 3.1.2 uses 'rep'.
        while not aLine.startswith("tree gen.") and not aLine.startswith("tree rep."):
            lNum += 1
            aLine = fLines[lNum].strip()

        # Get the tree lines.
        self.tLines = []
        while not aLine.startswith("end;"):
            self.tLines.append(aLine[5:])
            lNum += 1
            aLine = fLines[lNum].strip()
        
        # Read in the prams
        fName = "%s.run%i.p" % (self.mbBaseName, self.runNum)
        if self.directory != '.':
            fName = os.path.join(self.directory, fName)
        try:
            f = file(fName)
        except IOError:
            gm.append("Can't find prams file '%s'" % fName)
            raise Glitch, gm
        fLines = f.readlines()
        f.close()

        # Get the header line, starts with Gen
        lNum = 1
        aLine = fLines[lNum].strip()
        self.pramsHeader = aLine.split()
        if self.verbose >= 2:
            print "pramsHeader: %s" % self.pramsHeader
        assert self.pramsHeader[0] == 'Gen'

        # Collect pram lines
        lNum += 1
        aLine = fLines[lNum].strip()
        self.pLines = []
        while aLine:
            self.pLines.append(aLine)
            lNum += 1
            try:
                aLine = fLines[lNum].strip()
            except IndexError:
                break
        #print self.pLines
        self.nPrams = len(self.pramsHeader)
        if self.verbose >= 2:
            print "pram line length is %i" % self.nPrams

    def _getMrBayesSampleTree(self, sampNum):
        savedDoFastNextTok = var.nexus_doFastNextTok
        var.nexus_doFastNextTok = False
        tLine = self.tLines[sampNum]
        if self.verbose >= 3:
            print tLine
        f = cStringIO.StringIO(tLine)
        t = Tree()
        t.parseNexus(f, translationHash=self.translationHash, doModelComments=self.tree.model.nParts) # doModelComments is nParts
        var.nexus_doFastNextTok = savedDoFastNextTok
        t.taxNames = self.tree.taxNames

        for n in t.iterLeavesNoRoot():
            n.seqNum = t.taxNames.index(n.name)
        
        t.model = copy.deepcopy(self.model)

        pLine = self.pLines[sampNum]
        if self.verbose >= 3:
            print pLine
        splitPLine = pLine.split()

        pGenNum = int(splitPLine[0])
        splitTName = t.name.split('.')
        tGenNum = int(splitTName[1])
        if tGenNum != pGenNum:
            raise Glitch, "something wrong. tGenNum=%i, pGenNum=%i" % (tGenNum, pGenNum)
        if self.verbose >= 2:
            print "(zero-based) sample %i is gen %i" % (sampNum, tGenNum)

        #t.model.dump()

        splIndx = 3
        while splIndx < self.nPrams:
            pNum = 0
            #print "splIndx = %i, pramsHeader = %s" % (splIndx, self.pramsHeader[splIndx])
            if self.pramsHeader[splIndx].startswith('r(A<->C)'):
                if self.tree.model.nParts > 1:
                    try:
                        splitPramHeader = self.pramsHeader[splIndx].split('{')[1][:-1]
                        pNum = int(splitPramHeader)
                        pNum -= 1
                    except:
                        raise Glitch, "could not get the part number"
                thisSum = 0.0
                for i in range(6):
                    theFloat = float(splitPLine[splIndx])
                    t.model.parts[pNum].rMatrices[0].val[i] = theFloat
                    thisSum += theFloat
                    splIndx += 1
                factor = 1.0 / thisSum  # must sum to one
                for i in range(6):
                    t.model.parts[pNum].rMatrices[0].val[i] *= factor
            elif self.pramsHeader[splIndx].startswith('pi(A)'):
                if self.tree.model.nParts > 1:
                    try:
                        splitPramHeader = self.pramsHeader[splIndx].split('{')[1][:-1]
                        pNum = int(splitPramHeader)
                        pNum -= 1
                    except:
                        raise Glitch, "could not get the part number"
                thisSum = 0.0
                for i in range(4):
                    theFloat = float(splitPLine[splIndx])
                    t.model.parts[pNum].comps[0].val[i] = theFloat
                    thisSum += theFloat
                    splIndx += 1
                factor = 1.0 / thisSum # must sum to one
                for i in range(4):
                    t.model.parts[pNum].comps[0].val[i] *= factor
            elif self.pramsHeader[splIndx].startswith('alpha'):
                if self.tree.model.nParts > 1:
                    try:
                        splitPramHeader = self.pramsHeader[splIndx].split('{')[1][:-1]
                        pNum = int(splitPramHeader)
                        pNum -= 1
                    except:
                        raise Glitch, "could not get the part number"
                #print "got pNum = %i" % pNum
                #print "got splitPLine[%i] = %s" % (splIndx, splitPLine[splIndx])
                t.model.parts[pNum].gdasrvs[0].val[0] = float(splitPLine[splIndx])
                splIndx += 1
            elif self.pramsHeader[splIndx].startswith('pinvar'):
                if self.tree.model.nParts > 1:
                    try:
                        splitPramHeader = self.pramsHeader[splIndx].split('{')[1][:-1]
                        pNum = int(splitPramHeader)
                        pNum -= 1
                    except:
                        raise Glitch, "could not get the part number"
                t.model.parts[pNum].pInvar.val = float(splitPLine[splIndx])
                splIndx += 1
            elif self.pramsHeader[splIndx].startswith('m'):
                if self.tree.model.nParts > 1:
                    try:
                        splitPramHeader = self.pramsHeader[splIndx].split('{')[1][:-1]
                        pNum = int(splitPramHeader)
                        pNum -= 1
                    except:
                        raise Glitch, "could not get the part number"
                t.model.parts[pNum].relRate = float(splitPLine[splIndx])
                splIndx += 1
            else:
                print "splIndx=%i.  Got unknown pram %s.  Fix me!" % (splIndx, self.pramsHeader[splIndx])
                splIndx += 1
                

        if splIndx != len(splitPLine):
            raise Glitch, "Something is wrong.  After reading, splIndx=%i, but len split pram line=%i" % (
                splIndx, len(splitPLine))
        return t
    
        
