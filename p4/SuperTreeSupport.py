import sys, csv, random
from Tree import Tree
from func import read
from Var import var
from Glitch import Glitch
from p4.ReducedStrictConsensus import Intersection, TreeBuilderFromSplits

class SuperTreeInputTrees(object):

    def __init__(self, inputTree, distributionTrees=None):
        """
        SuperTreeInputTrees is a utility to create sets of input trees. 
        The input trees are primarily to be used to evaluate super tree
        construction methods. 
        
        Invocation removing a fixed number of taxa from each prospective input tree:
        
            stit = SuperTreeInputTrees(inputTree)
            stit.writeInputTreesToFile = True
            stit.outputFile = 'myInputtrees.tre'
            stit.noTaxaToRemove = 32 
            stit.noOutputTrees = 10
            stit.generateInputTrees()
        
        
        Invocation using built in distribution gathered from real world super tree cases::
        
            stit = SuperTreeInputTrees(inputTree)
            stit.writeInputTreesToFile = True
            stit.outputFile = 'myInputtrees.tre'
            stit.useTaxonDistribution = True
            stit.generateInputTrees()
        
        The user can generate a distribution of their own by supplying a list of p4 trees or a tree file. 
        The order of the trees is important, supertree and then all other trees. This goes for both list and 
        file. Like so::
        
            stit = SuperTreeInputTrees(inputTree, distributionTrees='myTreefile.nex')
            stit.writeInputTreesToFile = True
            stit.outputFile = 'myInputtrees.tre'
            stit.useTaxonDistribution = True
            stit.generateInputTrees()
        
        Placeholders which allow access to data after completed computations::
         
            stit.outputTrees 
            stit.dist
        
        """
        
        self.writeInputTreesToFile = False
        self.outputFile = 'inputtrees.tre'
        
        self.useTaxonDistribution = False #Set to False if you want to have a set number of taxa in the output trees 
        self.noTaxaToRemove = 32 #Only meaningful if setting useTaxonDistribution = False 
        self.noOutputTrees = 10
        
        gm = ['SuperTreeInputTrees()']
        
        if isinstance(inputTree, Tree):
            self.inputTree = inputTree            # not a list.
        elif type(inputTree) == type(""):
            var.trees = []
            read(inputTree)
            if len(var.trees) > 1:
                gm.append('Sorry, supply only one tree as supertree')
                raise Glitch, gm
            self.inputTree = var.trees.pop()       # this was originally a list, ie [var.trees.pop()]
        else:
            gm.append("Input tree was neither a p4 Tree nor a valid filename")
            gm.append("Got %s" % inputTree)
            raise Glitch, gm
        
        if not self.inputTree._taxNames:
            self.inputTree._setTaxNamesFromLeaves()
        
        self.outputTrees = []
        
        self.normalizedDist = []

        #Distributions gathered from real world supertree input
        #The dists are first a list of input tree taxon set sizes and the supertree taxon set size
        #Using this data we can normalize the dists to fit the size of trees we want

        #BunnyRSVNormal set from Wilkinson et al 2005, Syst Biol 54:823
#        self.dist = [[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 13, 13, 13, 13, 13, 14, 14, 15, 17, 17, 18, 18, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 22, 22, 23, 24, 25, 25, 25, 25, 25, 25, 26, 27, 28, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 37, 38, 38, 40, 40, 41, 47, 51, 51, 52, 52, 52, 68, 70, 78, 78, 79, 80, 80], 80]

        #CanidaeRVS set from Wilkinson et al 2005, Syst Biol 54:823
        #self.dist = [[3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 6, 7, 8, 8, 9, 10, 11, 11, 11, 12, 16, 16, 20, 23, 24, 30, 30, 33, 34, 34, 34, 34, 34], 34]

        #CarnivoraRVS set from Wilkinson et al 2005, Syst Biol 54:823
        #self.dist = [[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12], 12]

        #DavideDinoMRP set from Wilkinson et al 2005, Syst Biol 54:823
        #self.dist = [[4, 4, 4, 5, 6, 6, 6, 7, 8, 8, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 17, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 20, 20, 20, 22, 23, 23, 24, 24, 25, 26, 27, 27, 28, 28, 29, 29, 29, 29, 29, 30, 30, 30, 31, 31, 31, 31, 33, 33, 33, 33, 36, 37, 37, 38, 38, 39, 42, 45, 47, 48, 50, 53, 53, 66, 70, 71, 74, 74, 75, 75, 76, 78, 78, 80, 86, 86, 92, 94, 96, 100, 101, 102, 102, 103, 105, 110, 111, 111, 139, 148, 149, 153, 173, 199, 204, 217, 240, 269, 270, 271, 272, 272, 273, 273, 273, 273, 274, 275], 277]

        #FelidaeRVS set from Wilkinson et al 2005, Syst Biol 54:823
        self.dist = [[3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 7, 9, 9, 10, 10, 14, 16, 17, 24, 25, 28, 29, 29, 30, 30, 32, 34, 36, 36, 36, 36, 36, 36, 36, 36], 36]

        #KennedyPageData set from Wilkinson et al 2005, Syst Biol 54:823
        #self.dist = [[14, 16, 17, 20, 30, 30, 90], 122]

        #ViverridaeRVS set from Wilkinson et al 2005, Syst Biol 54:823
        #self.dist = [[4, 5, 10, 16, 19, 33, 34, 34, 34], 34]
        
        if distributionTrees:
            self.useTaxonDistribution = True
            if type(distributionTrees) == type([]):        
                for t in distributionTrees:
                    if not isinstance(t, Tree):
                        gm.append("Input trees should be a list of p4 Tree objects. Got %s" % t)
                        raise Glitch, gm
                superTree = distributionTrees.pop(0)
                inputTrees = distributionTrees
            elif type(distributionTrees) == type(""):
                var.trees = []
                read(distributionTrees)
                if len(var.trees) < 1:
                    gm.append('Sorry, at least one tree must be supplied as input tree')
                    raise Glitch, gm
                superTree = var.trees.pop(0)
                inputTrees = var.trees
            self._generateDistribution(superTree, inputTrees)

    def _generateDistribution(self, superTree, inputTrees):

        if not superTree._taxNames:
            superTree._setTaxNamesFromLeaves()

        distribution = []
        for tree in inputTrees:
            if not tree._taxNames:
                tree._setTaxNamesFromLeaves()
            distribution.append(len(tree.taxNames))
    
        distribution.sort()

        self.dist = [distribution,len(superTree.taxNames)]

#        print [distribution,len(superTree.taxNames)]

    #Prepares the dist by normalizing and changing it to reflect instead the no of taxa to remove.
    #This change is done simply to ease the other calculations
    def _prepareDistribution(self, noOfTaxaInSuperTree):
        #    print 'Supertree taxa:       ',noOfTaxaInSuperTree
        #    print 'Dist supertree taxa:  ',dist[1]
        normalizationIndex = self.dist[1]/float(noOfTaxaInSuperTree)
        #    print 'Normalization index:  ',normalizationIndex
        #    print 'Dist supertree taxa N:',dist[1]/normalizationIndex
        for size in self.dist[0]:
            normalized = size/normalizationIndex
            if int(noOfTaxaInSuperTree-normalized) > 2:
#                print '%s, %s ' % (size, int(noOfTaxaInSuperTree-normalized))
                self.normalizedDist.append(int(noOfTaxaInSuperTree-normalized))
    

    #Determines how we deal with the no of output trees
    #In the simple case this is a set number, but we can also use a distribution
    def _noOutputTrees(self):
        return self.noOutputTrees

    #Determines how we deal with the no of taxa to remove
    #In the simple case this is a set number, but we can also use a distribution
    def _noTaxaToRemove(self):
        if self.useTaxonDistribution:
            return random.choice(self.normalizedDist)
        else:
            return self.noTaxaToRemove 

    #Randomizes a taxon for removal, does so without regard for the shape of the tree
    def _getTaxaToRemove(self, nodes):
        node = nodes[random.randint(0, len(nodes)-1)]
        if node.isLeaf:
            return node 
        else:
            while not node.isLeaf:
                node = nodes[random.randint(0, len(nodes)-1)]
        return node

    def generateInputTrees(self):
        gm = ['SuperTreeInputTrees.generateInputTrees()']
    #Check if input values are valid
        if not self.useTaxonDistribution:
            if self.noTaxaToRemove >= len(self.inputTree.taxNames) - 2:
                gm.append('The number of taxa to remove would leave less than 3 taxa in the tree, quite uninformative')
                raise Glitch, gm

    #Prepare the distribution by normalizing it to the size of the input tree
        if self.useTaxonDistribution:
            self._prepareDistribution(len(self.inputTree.taxNames))

    #Creates the output trees and removes taxa from them accoring to the settings
    #Checks if the output tree has the correct number of taxa as a precation
        for i in range(self._noOutputTrees()):
            tree = self.inputTree.dupe()
            tree.name = 'inputtree' + str(i+1)
            taxa2Remove = self._noTaxaToRemove()
            for j in range(taxa2Remove):
                tree.removeNode(self._getTaxaToRemove(tree.nodes))
            
            tree._setTaxNamesFromLeaves()
            if len(tree.taxNames) == len(self.inputTree.taxNames) - taxa2Remove:
                self.outputTrees.append(tree)
            else:
                print 'Bugger, the correct number of taxa were not removed, taxa remaining: ', len(tree.taxNames)
                print 'Expected: ', len(self.inputTree.taxNames) - taxa2Remove
    
        #Writes the trees to file
        if self.writeInputTreesToFile:
            for tree in self.outputTrees:
            #    tree.draw()
                tree. writeNewick(fName=self.outputFile, withTranslation=0, translationHash=None, doMcmcCommandComments=0, toString=False, append=True)


class SuperTreeSupport(object):
    """Supertree support measures
    
    Super tree support can be used to calculate a number of support measures for a set of trees and 
    a supertree. The measures can be at split level and placed on the supertree for image production or
    at tree level with a number of summary measures. 
    
    The support of the input trees for a supertree is measured by counting the number of
    input trees that support(S), conflict(Q), permits(P) or are relevant(R) with the splits in the supertree. 
    
    Supply a supertree and the input trees used to create it. Filenames or trees will do. 
    A single supertree and a list of input trees. 
    
    For example::

        sts = SuperTreeSupport('supertree.nex', 'input.nex')

    or::
    
        read('input.nex')
        inputTrees = var.trees
        sts = SuperTreeSupport('supertree.nex', inputTrees)

        sts.doSaveDecoratedTree = True
        sts.decoratedFilename='mytree.nex'
        sts.doSaveIndexTree=False
        sts.indexFilename='mytreeIndex.nex'
        sts.csvFilename='mytreeIndex.csv'
        sts.doDrawTree=True
        sts.verbose=1

        sts.superTreeSupport()
    
    After completing the analysis there are a number of placeholders that allows access to the resulting data::
    
        sts.decoratedSuperTree
        sts.indexSuperTree
        sts.csvList
        sts.T, no. of input trees; 
        sts.L, no. of leaves; 
        sts.C, coverage (average proportion of leaves in the input tree); 
        sts.mean, mean taxon overlap among input trees
        sts.median, median taxon overlap among input trees
        sts.SC, number of supertree clades; 
        sts.U, no. of unsupported supertree clades; 
        sts.UC, no. of unsupported supertree clades that conflict with at least one input tree; 
        sts.UCC, no. of unsupported clades conflicting with all relevant input trees; 
        sts.QS, average qualitative support for supertree clades. Figures in parentheses are ranges.
        sts.S, average support
        sts.P, average permitted
        sts.Q, average conflict
        sts.R, average relevance
        sts.V, average V for supertree clades V = (s minus q)/(s + q)
        sts.VV, V+ = (s - q + p)/(s + q + p)
        sts.Vv, V minus = (s - q - p)/(s + q + p)
        sts.wV, wV = (ws minus q)/(ws + q)
        sts.wVV, wVV = (ws minus q +wp)/(ws + q + wp)
        sts.wVv, wVv = (ws minus q minus wp)/(ws + q + wp)
    
    These can be used for further analysis.
    
    Examples of support, conflict, relevance and permission: 
    
    **Support**::
    
        supertree:

                 +--------1:A
        +--------5:100
        |        +--------2:B
        0
        |--------3:C
        |
        +--------4:D

        input tree:

                 +--------1:A
        +--------5:100
        |        +--------2:B
        0
        |--------3:C


    **Conflict**::

        supertree: 

                 +--------1:A
        +--------5:100
        |        +--------2:B
        0
        |--------3:C
        |
        +--------4:D


        input tree: 

                 +--------1:A
        +--------5:100
        |        +--------2:C
        0
        |--------3:B
        |
        +--------4:D


    **Relevance**, i.e an input tree split that is irrelevant to the supertree split::

        supertree: 

                 +--------1:A
        +--------7:50
        |        +--------2:B
        |
        |--------3:C
        0
        |--------4:D
        |
        |--------5:E
        |
        +--------6:F

        input tree: 

                 +--------1:E
        +--------7:50
        |        +--------2:F
        |
        |--------3:C
        0
        |--------4:D
        |
        |--------5:A


    **Permission**, i.e. an input trees split that permits the supertree split but does not support it::

        supertree: 

                 +--------1:A
        +--------7:50
        |        +--------2:B
        |
        |--------3:C
        0
        |--------4:D
        |
        |--------5:E    
        |
        +--------6:F

        input tree: 

                 +--------1:A
        +--------7:50
        |        +--------2:B
        |        |
        0        |--------3:C
        |
        |--------4:D
        |
        |--------5:E    
   
    
    The analysed supertree can be decorated in three different ways. 
    
    I.   V = (s minus q)/(s + q) in standard mode, i.e. values between 0 - 100. 
    
    II.  V = (s minus q)/(s + q) in supertree mode, i.e. values between -1 - 1
    
    III. With indices that can be correlated to a csv file holding the support values. 
         This allows for further analysis of the support values and post mapping to the tree.
    """
    
    def __init__(self, supertree, inputTrees):
        
#        There are two ways of decorating the supertree with the support values. 
#        Standard conforms to the consensus tree tradition, i.e. values are presented between 
#        0 to 100 percent. Non standard adhears to the few supertree papers regarding support values
#        i.e -1 to 1. 
        self.doStandardDecoration=True
        
#        The decorated supertree can be saved to file
        self.doSaveDecoratedTree=False
        self.decoratedFilename='superTreeSupport.nex'
        
#        There is a option to save a supertree decorated with index values instead of support values. 
#        This can then be used with a csv file containing the support values for each index. 
#        Further analysis of the support values can be performed and then matched to the indecies in the
#        decorated supertree
        self.doSaveIndexTree=False
        self.indexFilename='supertreeIndex.nex'
        self.csvFilename='supertreeIndex.csv'
        
#        Draws the decorated supertree to screen
        self.doDrawTree=False
        
#        Produces output to screen
        self.verbose=1
        
#        Placeholders that allows access to the data after completing calculations
        self.decoratedSuperTree = None
        self.indexSuperTree = None
        self.csvList = None
        
#       Keeps track of splits for producing output  
        self.indexIntersections = []
        self.csvValues = []
        self.intersections = []
        
#        Let t be the number of input trees, 
#        s the number of input trees supporting a supertree clade, 
#        r the number of input trees that are irrelevant to the supertree clade, 
#        q the number of input trees that conflict with the supertree clade, 
#        p the number of input trees that permit the supertree clade, 
#        so that t = p + q + r + s.

        self.T = 0 #no. of input trees; 
        self.L = 0 #no. of leaves; 
        self.C = 0.0 #coverage (average proportion of leaves in the input tree); 
        self.SC = 0 #number of supertree clades; 
        self.U = 0 #no. of unsupported supertree clades; 
        self.UC = 0 #no. of unsupported supertree clades that conflict with at least one input tree; 
        self.UCC = 0 #no. of unsupported clades conflicting with all relevant input trees; 
        self.QS = 0.0 #average qualitative support for supertree clades. Figures in parentheses are ranges.
        self.S = 0.0 # average support
        self.P = 0.0 # average permitted
        self.Q = 0.0 # average conflict
        self.R = 0.0 # average relevance
        self.wS = 0.0 # average weighted support 
        self.wP = 0.0 # average weighted permitance
        self.V = 0.0 # average V for supertree cladesV = (s minus q)/(s + q)
        self.VV = 0.0 # V+ = (s minus q +p)/(s + q + p)
        self.Vv = 0.0 # V minus = (s minus q minus p)/(s + q + p)
        self.wV = 0.0 # wV = (ws minus q)/(ws + q)
        self.wVV = 0.0 # wVV = (ws minus q +wp)/(ws + q + wp)
        self.wVv = 0.0 # wVv = (ws minus q minus wp)/(ws + q + wp)
        
        
        gm = ['SuperTreeSupport()']
        
        var.warnReadNoFile = False
        
        if type(inputTrees) == type([]):
            for t in inputTrees:
                if not isinstance(t, Tree):
                    gm.append("Input trees should be a list of p4 Tree objects. Got %s" % t)
                    raise Glitch, gm
            self.inputTrees = inputTrees
        elif type(inputTrees) == type(""):
            var.trees = []
            read(inputTrees)
            if len(var.trees) < 1:
                gm.append('Sorry, at least one tree must be supplied as input tree')
                raise Glitch, gm
            self.inputTrees = var.trees
        else:
            gm.append("Input trees are neither a list of p4 Tree objects nor a valid filename.")
            raise Glitch, gm
        
        if isinstance(supertree, Tree):
            self.supertree = supertree            # not a list.
        elif type(supertree) == type(""):
            var.trees = []
            read(supertree)
            if len(var.trees) > 1:
                gm.append('Sorry, supply only one tree as supertree')
                raise Glitch, gm
            self.supertree = var.trees.pop()       # this was originally a list, ie [var.trees.pop()]
        else:
            gm.append("Supertree was neither a p4 Tree nor a valid filename")
            gm.append("Got %s" % supertree)
            raise Glitch, gm
        
        for tree in self.inputTrees:
            if not tree._taxNames:
                tree._setTaxNamesFromLeaves()

        #Mean and median overlap of the input trees
        overlapList = []
        meanOverlap = 0.0
        index = 0
        for i in range(0,len(self.inputTrees)-1):
            for j in range(i+1,len(self.inputTrees)):
                overlap = len(set(self.inputTrees[i].taxNames).intersection(set(self.inputTrees[j].taxNames)))
                overlapList.append(overlap)
                meanOverlap += overlap
                index += 1
        
        if index == 0:
            self.mean = 0
            self.median = 0
        else:
            self.mean = meanOverlap/index
            overlapList.sort()
            self.median = overlapList[len(overlapList)/2]  
        
        commonLeafSet = CommonLeafSet()
        self.splits = commonLeafSet.updateTreesToCommonLeafSet([self.inputTrees, [self.supertree]])
        self.bitkeys = commonLeafSet.getCommonBitkeys()
        self.taxnames = commonLeafSet.getCommonTaxNames()
        self.taxa2Bitkey = commonLeafSet.getCommonTaxa2Bitkey()

    def superTreeSupport(self):
        
        """
        Perform analyses on a number of input trees and the resulting supertree. 
        
        """
        gm = ['superTreeSupport()']

        allOnes = 2L**(len(self.bitkeys)) - 1

        inputSplits = self.splits[0]
        self.T = len(inputSplits)
        self.L = len(self.bitkeys)
        #if 'ROOT' in self.taxnames: 
        #   self.L -= 1
        # if 'mrp_outgroup' in self.taxnames:
        #    self.L -= 1
        
        coverage = 0.0
        for splits in inputSplits:
            t = self.L - self.popcount(splits[0][1], self.bitkeys)
            if 'ROOT' in self.taxnames and splits[0][1] & self.taxa2Bitkey['ROOT']:
                t += 1
            coverage += t
        
        self.C = (coverage/self.T)/self.L
        #print "self.splits = %s" % self.splits
        if self.splits[1]:
            supertreeSplits = self.splits[1][0]
        else:
            supertreeSplits = []
        self.SC = len(supertreeSplits)
        
        if self.verbose > 1:
            print 'Identifying supporting, conflicting and neutral input trees,'
            sys.stdout.flush()

        """
        A dictionary supertreeDict = {} will keep track of the supertree splits and their support
        The indeces from the supertreeSplits list will be used to keep track of the splits
        i.e. supertreeSplits[1] info will be in supertreeDict[1]
        Each entry in the dict will be a list formated like so
        [input trees hard support, input trees soft support(entailed), input trees neutral, input trees conflict]
        """

        def isPowerOf2(number): return (not(number&(number-1)))
        
        supertreeDict = {}
        for i in range(self.SC):
#           Used to store each supertree clades support values 
#                              [S,  P,  R,  Q,  WS, WP, V,  V+, V-, wV, wV+,wV-]
            supertreeDict[i] = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        
        supertreeSplitDict = {}
        for j in range(self.SC):
            supertreeSplitDict[j] = -1
        
        for splitSet in inputSplits:
            supportedSupertreeSplits = []
            permittedSupertreeSplits = []
            taxSplit = splitSet[0][1] ^ allOnes
            
            allSuperTreeSplitsResovled = False
            for split in splitSet:

                if allSuperTreeSplitsResovled:
                    break
            
                allSuperTreeSplitsResovled = True
                for i in range(self.SC):
                    if supertreeSplitDict[i] < 0:
                        allSuperTreeSplitsResovled = False
                        workingSplit = (supertreeSplits[i][0] & splitSet[0][1]) ^ supertreeSplits[i][0]
                        if not self.popcountLimit(taxSplit & workingSplit, 2) or not self.popcountLimit(taxSplit ^ workingSplit, 2):
                            supertreeSplitDict[i] = 2
                            
                        else:
                            
                            ones = split[0] & workingSplit
                            one1 = 0
                            if ones > 0:
                                if not isPowerOf2(ones):
                                    one1 = 2
                                else:
                                    one1 = 1
#                        print 'one1: ', one1
                        
                            zeroOnes = (split[0] ^ workingSplit) & workingSplit
                            zero1 = 0
                            if zeroOnes > 0:
                                if not isPowerOf2(zeroOnes):
                                    zero1 = 2
                                else:
                                    zero1 = 1
#                        print 'zero1: ', zero1
                        
                            oneZeroes = (split[0] ^ workingSplit) & split[0]
                            one0 = 0
                            if oneZeroes > 0:
                                if not isPowerOf2(oneZeroes):
                                    one0 = 2
                                else:
                                    one0 = 1
#                        print 'one0: ', one0
  
                            zeroZeroes = ((workingSplit ^ allOnes) & (split[0] ^ allOnes)) ^ split[1]
                            zero0 = 0 
                            if zeroZeroes > 0:
                                if not isPowerOf2(zeroZeroes):
                                    zero0 = 2
                                else:
                                    zero0 = 1
#                        print 'zero0: ', zero0
                                 
                        
                            if one1 and zero1 and one0 and zero0:
                                supertreeSplitDict[i] = 3
                        
                            if one1 and zero0 and (not zero1 and not one0):
                                #self.printSplitList([split, supertreeSplits[i]], self.bitkeys)
                                supportedSupertreeSplits.append(i)
                                supertreeSplitDict[i] = 0
                        
                            if zero1 and one0 and (not one1 and not zero0):
                                #self.printSplitList([split, supertreeSplits[i]], self.bitkeys)
                                supportedSupertreeSplits.append(i)
                                supertreeSplitDict[i] = 0

            
            for i in range(self.SC):
                if supertreeSplitDict[i] < 0:
                    supertreeDict[i][1] = supertreeDict[i][1] + 1
                    permittedSupertreeSplits.append(i)
                supertreeDict[i][supertreeSplitDict[i]] = supertreeDict[i][supertreeSplitDict[i]] + 1 
                supertreeSplitDict[i] = -1
                    
            if len(permittedSupertreeSplits) > 0:
                wp = 1.0/len(permittedSupertreeSplits)
                for index in permittedSupertreeSplits:
                    supertreeDict[index][5] = supertreeDict[index][5] + wp
                
            if len(supportedSupertreeSplits) > 0:
                ws = 1.0/len(supportedSupertreeSplits)
                for index in supportedSupertreeSplits:
                    supertreeDict[index][4] = supertreeDict[index][4] + ws
            
        if self.verbose > 1:
            print 'done.'
                    
        """
        How many input trees containing the taxa in the supertree split, conflict, agree, neutral. 
        
        """
        if self.verbose > 1:
            print 'Creating statistics and supertree with support values,'
            
        self._calcSupport(supertreeDict)
        
        if self.doSaveIndexTree:
            self.saveIndexTree(supertreeSplits, supertreeDict)
        
        if self.doDrawTree or self.doSaveDecoratedTree:
            self.buildDecoratedTree(supertreeSplits, supertreeDict)
        
        if self.doDrawTree:
            self.decoratedSuperTree.draw()
            print ''
        
        if self.doSaveDecoratedTree:
            self.decoratedSuperTree.writeNexus(self.decoratedFilename)        
        
        if self.verbose == 1:
            self.printReducedOutput()
        elif self.verbose > 1:
            self.printReducedOutput(printInstructions=True)
            
            

    def buildDecoratedTree(self, supertreeSplits, supertreeDict):
            
        for key, value in supertreeDict.items():
            inter = Intersection(supertreeSplits[key][0], supertreeSplits[key][1], 0, 0)
            name = ''
            if self.doStandardDecoration:
                name = str(int(round(100.0*(1.0 + value[9])/2.0)))
            else:
                name = str(value[9])
            inter.name = name
            self.intersections.append(inter)
            
        treeBuilderFromSplits = TreeBuilderFromSplits(self.bitkeys, self.taxnames)
        
        self.decoratedSuperTree = treeBuilderFromSplits.buildTreeFromInformativeList(self.intersections, treeName='decoratedSupertree')  

    def saveIndexTree(self, supertreeSplits, supertreeDict):
        for key, value in supertreeDict.items():
            inter = Intersection(supertreeSplits[key][0], supertreeSplits[key][1], 0, 0)
            if self.doStandardDecoration:
                self.csvValues.append([key,int(round(100.0*(value[0]/self.T))),int(round(100.0*(value[1]/self.T))), 
                                           int(round(100.0*(value[2]/self.T))), int(round(100.0*(value[3]/self.T))), 
                                           int(round(100.0*(value[4]/self.T))), int(round(100.0*(value[5]/self.T))),
                                           value[6], value[7], value[8], value[9], value[10], value[11]])
            else:
                self.csvValues.append([key, value[0]+value[1],value[0],value[1],value[2],value[3], value[4],value[5], value[6], value[7], value[8], value[9], value[10], value[11]]) 
            inter.name = str(key)
            self.indexIntersections.append(inter)
            
        treeBuilderFromSplits = TreeBuilderFromSplits(self.bitkeys, self.taxnames)
        
        self.indexSuperTree = treeBuilderFromSplits.buildTreeFromInformativeList(self.indexIntersections, treeName='indexSupertree')  
        
        self.indexSuperTree.writeNexus(self.indexFilename)
                
        writer = csv.writer(open(self.csvFilename, "wb"))
                
        writer.writerows(self.csvValues)
    
    def _calcSupport(self, supertreeDict):
        
        for value in supertreeDict.values():
            S = value[0]
            P = value[1]
            Q = value[3]
            R = value[2]
            
            #  print 'S+P+Q+R: %s, %s, %s, %s' % (S, P, Q, R)
            
            wS = value[4]
            wP = value[5]    
            self.S += S
            self.P += P
            self.Q += Q
            self.R += R
            self.wS += wS
            self.wP += wP
            
            if S < 1:
                self.U += 1
            if S < 1 and Q > 0:
                self.UC += 1
            if S < 1 and Q > 0 and Q == self.T - R:
                self.UCC += 1
            
            if S + Q > 0:
                value[6] = (S - Q)/(S + Q)
                self.V += value[6]
            if wS + Q > 0:
                value[9] = (wS - Q)/(wS + Q)
                self.wV += value[9]
            if S + Q + P > 0:
                value[7] = (S - Q + P)/(S + Q + P)
                self.VV += value[7]
                value[8] = (S - Q - P)/(S + Q + P)
                self.Vv += value[8]
            if wS + Q + wP > 0:
                value[10] = (wS - Q + wP)/(wS + Q + wP)
                self.wVV += value[10]
                value[11] = (wS - Q - wP)/(wS + Q + wP)
                self.wVv += value[11]
        if 0:
            print "self.SC =", self.SC
            print "self.S =", self.S
            print "self.P =", self.P
            print "self.Q =", self.Q
            print "self.R =", self.R
            print "self.V =", self.V
            print "self.VV =", self.VV
            print "self.Vv =", self.Vv
            print "self.wV =", self.wV
            print "self.wVV =", self.wVV
            print "self.wVv =", self.wVv

        if self.SC:
            self.S = self.S/self.SC
            self.P = self.P/self.SC
            self.Q = self.Q/self.SC
            self.R = self.R/self.SC
            self.V = self.V/self.SC
            self.VV = self.VV/self.SC
            self.Vv = self.Vv/self.SC
            self.wV = self.wV/self.SC
            self.wVV = self.wVV/self.SC
            self.wVv = self.wVv/self.SC
        

    def printReducedOutput(self, printInstructions=False):
        if printInstructions: 
            print 'I = no. of input trees;'
            print 'L = no. of leaves;'
            print 'Mean = mean taxon overlap among input trees'
            print 'Medi = median taxon overlap among input trees'
            print 'C = coverage (average proportion of leaves in the input tree);'
            print 'SC = number of supertree clades;'
            print 'U = no. of unsupported supertree clades;'
            print 'U* = no. of unsupported supertree clades that conflict with at least one input tree;'
            print 'U** = no. of unsupported clades conflicting with all relevant input trees;' 
            print 'V = avg of support-conflict / support+conflict;'
            print 'V+ = avg of support-conflict+permitting / support+conflict+permitting ;' 
            print 'V- = avg of support-conflict-permitting / support+conflict+permitting ;' 
            print 'wV = avg of weighted support-conflict / weighted support+conflict;'
            print 'wV+ = avg of weighted support-conflict+weighted permitting / weighted support+conflict+weighted permitting ;' 
            print 'wV- = avg of weighted support-conflict-weighted permitting / weighted support+conflict+weighted permitting ;'
            print 'I     L     Mean  Medi  C     SC    U     U*    U**   V     V+    V-    wV    wV+   wV-'

        print repr(self.T).ljust(6)[:5],
        print repr(self.L).ljust(6)[:5],
        print repr(self.mean).ljust(6)[:5],
        print repr(self.median).ljust(6)[:5],
        print repr(self.C).ljust(6)[:5],
        print repr(self.SC).ljust(6)[:5],
        print repr(self.U).ljust(6)[:5],
        print repr(self.UC).ljust(6)[:5],
        print repr(self.UCC).ljust(6)[:5],
        print repr(self.V).ljust(6)[:5],
        print repr(self.VV).ljust(6)[:5],
        print repr(self.Vv).ljust(6)[:5],
        print repr(self.wV).ljust(6)[:5],
        print repr(self.wVV).ljust(6)[:5],
        print repr(self.wVv).ljust(6)[:5]
                          
    def popcount(self, n, bitkeys):
        count = 0
        for bk in bitkeys:
            if n < bk:
                return count
            if bk & n:
                count += 1
        return count
    
    def popcountLimit(self, n, limit):
        count = 0
        for bk in self.bitkeys:
            if bk & n:
                count += 1
            if n >= limit:
                return True
        return False
    
    def printSplitList(self, list, bitkeys):
        for split in list:
            intersection = ''
            for bk in bitkeys:
                if bk & split[0]:
                    intersection = intersection + '1'
                elif bk & split[1]:
                    intersection = intersection + '?'
                else:
                    intersection = intersection + '*'
            print 'Split: %s' % (intersection)
        print ' '
    

class CommonLeafSet(object):
               
    def updateToCommonLeafSetWithTaxa(self, tfl, dict, bitkeys, taxnames):
        gm = ['updateToCommonLeafSetWithTaxa(tfl, dict, bitkeys, taxnames)']
        
        uniqueSet = set()
        for i in range(tfl.nSamples):
            t = tfl.getTree(i)
            if not t._taxNames:
                t._setTaxNamesFromLeaves()
            for name in t.taxNames:
                uniqueSet.add(name)
        list = []
        for name in uniqueSet:
            list.append(name)
               
        for name in list:
            if not name in taxnames:
                gm.append('Found taxa in that does not appear in the supplied taxa list: ' + name)
                raise Glitch, gm
               
        for name in taxnames:
            if not name in list:
                gm.append('Found taxa in supplied taxa list that does not appear in the supertree: '+ name)
                raise Glitch, gm
        
        self.weights = []
        splits = []
        for i in range(tfl.nSamples):
            t = tfl.getTree(i)
            if not t._taxNames:
                t._setTaxNamesFromLeaves()
            t.missingTaxa = 0L
            for name in taxnames:
                if not t.taxNames.count(name):
                    t.missingTaxa = t.missingTaxa | dict[name]
            for n in t.iterLeavesNoRoot():
                n.br.rc = dict[n.name]
            t._makeRCSplitKeys()
            weight = 1.0
            if hasattr(t, "weight"):
                if t.weight != None:
                    weight = t.weight
            elif hasattr(t, "recipWeight"):
                if t.recipWeight != None:
                    weight = 1.0/int(t.recipWeight)
            self.weights.append(weight)
            t.splits =[]
            for n in t.iterInternalsNoRoot():
                for m in n.br.rcList:
                    t.splits.append([m,t.missingTaxa])
            splits.append(t.splits)
        
        return splits    
        
    def updateTreesToCommonLeafSet(self, trees):
        """
        Creates common leafset 
        Requires input to be a list like so: [[inputtrees][supertree]]
        """
        uniqueSet = set()
    
        sampler = float(len(trees))/10
        
        index = 0.0
        
        for treeList in trees:
            for t in treeList:
                if t.root.getNChildren() <= 2:
                    t.addSibLeaf(t.root, "ROOT")
                index += 0.5
                
                if not t._taxNames:
                    t._setTaxNamesFromLeaves()
                for name in t.taxNames:
                    uniqueSet.add(name)
                if index % sampler == 0:
                    print '%s trees processed ' % (int (index))
                    sys.stdout.flush()
                
        self.list = []
        for name in uniqueSet:
            self.list.append(name)
        
        self.list.sort()
        
        self.dict = {}
        self.bitkeys = []
        for i in range(len(self.list)):
            self.bitkeys.append(1L << i)
            self.dict[self.list[i]] = 1L << i 
        
        self.weights = []
        self.treeNames = []
        tflsSplits = []
        for treeList in trees:
            splits = []
            for t in treeList:
                index += 0.5
                self.treeNames.append(t.name)
                if not t._taxNames:
                    t._setTaxNamesFromLeaves()
                t.missingTaxa = 0L
                for name in self.list:
                    if not t.taxNames.count(name):
                        t.missingTaxa = t.missingTaxa | self.dict[name]
                for n in t.iterLeavesNoRoot():
                    n.br.rc = self.dict[n.name]
                t._makeRCSplitKeys()
                weight = 1.0
                if hasattr(t, "weight"):
                    if t.weight != None:
                        weight = t.weight
                elif hasattr(t, "recipWeight"):
                    if t.recipWeight != None:
                        weight = 1.0/int(t.recipWeight)
                self.weights.append(weight)
                t.splits =[]
                for n in t.iterInternalsNoRoot():
                    for m in n.br.rcList:
                        t.splits.append([m,t.missingTaxa])
                if len(t.splits) > 0:
#                    for split in t.splits:
#                        intersection = ''
#                        for bk in self.bitkeys:
#                            if bk & split[0]:
#                                intersection = intersection + '1'
#                            elif bk & split[1]:
#                                intersection = intersection + '?'
#                            else:
#                                intersection = intersection + '*'
#                        print 'Split: %s' % (intersection)
#                    print ''
            
                    splits.append(t.splits)
                if index % sampler == 0:
                    print '%s trees processed ' % (int (index))
                    sys.stdout.flush()
            tflsSplits.append(splits)
            
        return tflsSplits
    
               
    def updateToCommonLeafSet(self, tfls):
        uniqueSet = set()
        
        trees = 0
        for tfl in tfls:
            trees += tfl.nSamples
        
        sampler = float(trees)/10
        
        index = 0.0
        for tfl in tfls:
            for i in range(tfl.nSamples):
                index += 0.5
                t = tfl.getTree(i)
                if not t._taxNames:
                    t._setTaxNamesFromLeaves()
                for name in t.taxNames:
                    uniqueSet.add(name)
                if index % sampler == 0:
                    print '%s trees processed ' % (int (index))
                    sys.stdout.flush()
                
        self.list = []
        for name in uniqueSet:
            self.list.append(name)
        
        self.list.sort()
        
        self.dict = {}
        self.bitkeys = []
        for i in range(len(self.list)):
            self.bitkeys.append(1L << i)
            self.dict[self.list[i]] = 1L << i 
        
        self.weights = []
        self.treeNames = []
        tflsSplits = []
        for tfl in tfls:
            splits = []
            for i in range(tfl.nSamples):
                index += 0.5
                t = tfl.getTree(i)
                self.treeNames.append(t.name)
                if not t._taxNames:
                    t._setTaxNamesFromLeaves()
                t.missingTaxa = 0L
                for name in self.list:
                    if not t.taxNames.count(name):
                        t.missingTaxa = t.missingTaxa | self.dict[name]
                for n in t.iterLeavesNoRoot():
                    n.br.rc = self.dict[n.name]
                t._makeRCSplitKeys()
                weight = 1.0
                if hasattr(t, "weight"):
                    if t.weight != None:
                        weight = t.weight
                elif hasattr(t, "recipWeight"):
                    if t.recipWeight != None:
                        weight = 1.0/int(t.recipWeight)
                self.weights.append(weight)
                t.splits =[]
                for n in t.iterInternalsNoRoot():
                    for m in n.br.rcList:
                        t.splits.append([m,t.missingTaxa])
                splits.append(t.splits)
                if index % sampler == 0:
                    print '%s trees processed ' % (int (index))
                    sys.stdout.flush()
            tflsSplits.append(splits)
            
        return tflsSplits
        
    def getCommonTaxa2Bitkey(self):
        return self.dict
        
    def getCommonTaxNames(self):
        return self.list 
    
    def getCommonBitkeys(self):
        return self.bitkeys
    
    def getWeights(self):
        return self.weights
    
    def getTreeNames(self):
        return self.treeNames
        
    
