from Glitch import Glitch
from TreePartitions import TreePartitions
from Tree import Tree
from Trees import Trees
from func import read, var
import sys, csv, random
from math import log, factorial, floor

ROOT_NODE_NAME = 'A_NAME_NOT_EASILY_FOUND_IN_A_TREE'

TAXON_SEPARATOR = ':'

class TreeSubsets(object):
    
        def __init__(self, inputTrees):
            gm = ['TreeSubsets()']
            if type(inputTrees) == type([]):
                trees = [] 
                for t in inputTrees:
                    if not isinstance(t, Tree):
                        gm.append("Input trees should be a list of p4 Tree objects. Got %s" % t)
                        raise Glitch, gm
                    else:
                        trees.append(t)
                if len(inputTrees) < 1:
                    gm.append('Sorry, at least one tree must be supplied as input tree')
                    raise Glitch, gm
                self.inputTrees = trees
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
        
            self.nonUnrootedTrees = []
            for t in self.inputTrees:
                if not len(t.taxNames):           
                    t._setTaxNamesFromLeaves()
#                    self.nonUnrootedTrees.append(t.dupe())
                if t.root.getNChildren() <= 2:
                    t.removeRoot()
        
            self.treeSets = []
            self.subtrees = {}
            self._buildConsensus()
            self.noTaxonSets = 0

        def _buildConsensus(self):
            trees = Trees(self.inputTrees)
        
            tp = TreePartitions(trees)
            self.consensusTree = tp.consensus(minimumProportion=0.49)
            self.consensusTree.setPreAndPostOrder()
            for n in self.consensusTree.root.iterInternals():
                if n != self.consensusTree.root:
                    if n.br.support:
                        n.name = '%i' % floor(100. * float(n.br.support))
        
        def getConsensusTree(self):
            return self.consensusTree
        
        def _printConsensus(self):
            self.consensusTree.draw(width=80, showNodeNums=0)            
        
        
        
        def _getTaxonSets(self):
            taxonSets = []
            taxNames = []
            t = self.consensusTree.dupe()
            while self.isSubtreeDivideable(t.root):
                
                for n in t.iterInternalsNoRootPostOrder():
                    if n.name == '100':
                        taxa = t.getAllLeafNames(n)                       
                        taxname = ''
                        for taxon in taxa: 
                            taxname = taxname + taxon + TAXON_SEPARATOR
                        taxname = taxname[:len(taxname)-1]
                        list = []
                        taxonSets.append(taxa)
                        taxNames.append(taxname)
                        t.addSibLeaf(n.parent, taxname)
                        t.removeNode(n, alsoRemoveSingleChildParentNode=False)
#                        t.draw(width=80, showNodeNums=0)
                        break
            
#            print 'TaxonSets'
#            for taxonSet in taxonSets:
#                print taxonSet
#            print 
            return taxonSets, taxNames
        
        def getSubTreesAndTaxonSetsFromInputTrees(self, verbose=False):

            taxonSets, taxNames = self._getTaxonSets()
            cherries = []
            for set in taxonSets:
                if len(set) == 2:
                    cherries.append(set)
            self.noTaxonSets = len(taxonSets)
            subtreeDicts = []
            
            if self.noTaxonSets == 0:
                for t in self.inputTrees:
                    subTrees = {}
                    subTrees['tree'] = t
                    subtreeDicts.append(subTrees)
                return subtreeDicts, ['tree']
            
            if verbose:
                index = 1
                sampler = 1
                if len(self.inputTrees) > 20:
                    sampler = int(float(len(self.inputTrees)/20))
                sys.stdout.write('0%                    100%\n')
                sys.stdout.write('  ')
                sys.stdout.flush()
            
            for t in self.inputTrees:
                if verbose:
                    index += 1
                    if (index + 1) % sampler == 0:
                        sys.stdout.write('.')
                        sys.stdout.flush()
                subtreeDicts.append(self._getSubTreesFromTree(t, taxonSets))
                t = None
                
            if verbose:
                sys.stdout.write('\n')
                sys.stdout.flush()
#            print 'SubtreeSets:'
#            for dict in subtreeDicts:
#                for name, subtree in dict.iteritems():
#                    print name
#                    subtree.draw(width=80, showNodeNums=0)
#                print '********************************************************'

            return subtreeDicts, taxNames, taxonSets, cherries
        
        def _getSubTreesFromTree(self, t, taxonSets):
            
            subTrees = {}
            
            self._getSubTreesFromSubTree(t.root, t, taxonSets, subTrees)
            
            return subTrees
            
        def getNamesInOrder(self):
            
            taxa = ''
            for leaf in self.consensusTree.iterLeavesPostOrder():
                taxa += leaf.name
            return taxa
            
        def _getSubTreesFromSubTree(self, node, t, taxonSets, subTrees):
            
            if node.isLeaf:
                return [node.name]
            
            names = []
            names.extend(self._getSubTreesFromSubTree(node.leftChild, t, taxonSets, subTrees))
            sibling = node.leftChild.sibling
            while sibling:
                if sibling.isLeaf:
                    names.append(sibling.name)
                else:
                    names.extend(self._getSubTreesFromSubTree(sibling, t, taxonSets, subTrees))
                sibling = sibling.sibling
            
            match = False
            taxset = None
            for taxonSet in taxonSets:
                match = True
                if len(taxonSet) == len(names):
                    
                    for taxname in taxonSet: 
                        if not taxname in names:
                            match = False
                            break
                else:
                    match = False      
                if match:
                    taxset = taxonSet
                    break
            
            if match:
                name = ''
                for taxon in taxset: 
                    name = name + taxon + TAXON_SEPARATOR
                name = name[:len(name)-1]

                if len(taxonSet) > 2:
                    subtree = t.dupeSubTree(node, True)
                    subtree.removeRoot()
                    subTrees[name]=subtree

                node.name = name
                node.isLeaf = True
                node.leftChild = None
                return [name]
            
            return names
            

        def isSubtreeDivideable(self, node):
            for n in node.iterInternals():
                if n.name == '100':
                    return True
            return False
            
class CherryRemover(object):
    
    def __init__(self, inputTrees):
        
        gm = ['CherryRemover()']
        if type(inputTrees) == type([]):
            trees = [] 
            for t in inputTrees:
                if not isinstance(t, Tree):
                    gm.append("Input trees should be a list of p4 Tree objects. Got %s" % t)
                    raise Glitch, gm
                else:
                    trees.append(t.dupe())
            if len(inputTrees) < 1:
                gm.append('Sorry, at least one tree must be supplied as input tree')
                raise Glitch, gm
            self.inputTrees = trees
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
    
        #print 'Got %s trees' % (len(self.inputTrees))
    
        self.taxname2Subtree = {}
        
        self.nonUnrootedTrees = []
        for t in self.inputTrees:
            if not len(t.taxNames):           
                t._setTaxNamesFromLeaves()
            self.nonUnrootedTrees.append(t.dupe())
            if t.root.getNChildren() <= 2:
                t.removeRoot()
       
        self.reducedTrees = []
        
        self.removeInfo = [] #list of lists like so [[name2ReplaceWith, taxa2remove],[name2ReplaceWith, taxa2remove],,,,,]
        
        self.removableNodes = []
        
        self.consensusTree = None
        
        self.allowMultiFurcatingCherries = False
        
        self.cherryCutOff = 1.0
        
    def _buildConsensus(self):
        trees = Trees(self.inputTrees)
        
        tp = TreePartitions(trees)
        self.consensusTree = tp.consensus(minimumProportion=0.49)
        self.consensusTree.setPreAndPostOrder()
        for n in self.consensusTree.root.iterInternals():
            if n != self.consensusTree.root:
                if n.br.support:
                    n.name = '%i' % round(100. * float(n.br.support))

    def _identifySubtreesForRemoval(self):
        node = None
        for n in self.consensusTree.iterInternalsNoRootPostOrder():
            n.removable = False
            if n.br.support >= self.cherryCutOff:

                if self.allowMultiFurcatingCherries:
                    n.removable = True
                    for l in n.iterChildren():
                        if not l.isLeaf:
                            if not l.removable:
                                n.removable = False
                    
                else:
                    if n.getNChildren() == 2:
                        n.removable = True
                        for l in n.iterChildren():
                            if not l.isLeaf:
                                if not l.removable:
                                    n.removable = False
                                
                if n.removable:
                    node = n
            elif node:
                self.removableNodes.append(node)
                node = None
        if node:
            self.removableNodes.append(node)
                
    def _buildRemoveInfo(self):
        
        for n in self.removableNodes:
            name = ''
            list = []
            for l in n.iterLeaves():
                name += l.name + '-'
                list.append(l.name)
            name = name[:len(name)-1]
            self.removeInfo.append([name,list])
            
        
    def _replaceAndRemove(self):
        
        first = True
        for t in self.nonUnrootedTrees:
            for info in self.removeInfo:
                for n in t.root.iterPostOrder():
                    if n.name == info[1][0]:
                        n.name = info[0]
                        break
                for name in info[1][1:]:
                    t.removeNode(name, alsoRemoveSingleChildParentNode=True, alsoRemoveBiRoot=False, alsoRemoveSingleChildRoot=True)
            if first:
                first = False
                t._setTaxNamesFromLeaves()
                taxnames = t.taxNames
            else:
                t.taxNames = taxnames
    
    def saveTrees(self, filename):
        self.nonUnrootedTrees[0]._setTaxNamesFromLeaves()
        trees = Trees(self.nonUnrootedTrees)
        trees.writeBranchLengths = False
        trees.writeNexus(filename, withTranslation=1)
        
    
    def getTrees(self):
        return self.nonUnrootedTrees
    
    def _printTrees(self):
        for t in self.nonUnrootedTrees:
            t.draw()
    
    def printConsensus(self):
        self._buildConsensus()
        self.consensusTree.draw(width=80, showNodeNums=0)
    
    def removeCherries(self):
        
        self._buildConsensus()
        
        self._identifySubtreesForRemoval()
        
        self._buildRemoveInfo()
        
        self._replaceAndRemove()
  
        
class LeafSupport(object):
            
    def __init__(self, inputTrees, identicalTaxaSets=True, removeCherries=False):
        """Calculates three different leaf support values for the taxa.
	
        The inputTrees arg should be a tree file or a list of p4 Tree objects. 

        
        rooted=True  - if the trees are to be treated as rooted 
        rooted=False - if the trees are to be treated as unrooted
	(True by default)
        
        identicalTaxaSets=True - if all trees have identical sets of taxa  
        identicalTaxaSets=False- if sets of taxa differ among trees
        
        equalDistUnresolved=True - if any unresolved quartets or
                                   triplets should be evenly divided
                                   among the resolved
        equalDistUnresolved=False - if distributed according the the
	                            distribution of the resolved
	                            quartets/triplets
        (True by default)
        
        treatUnresolvedAsValid=True - if unresolved quartet/triplets
	                              should be treated as a valid
	                              option, ie polytomies are
	                              accepted
        treatUnresolvedAsValid=False - polytomies should be treated as
	                              unwanted and included taxa
	                              punished
        (True by default)
	
        exploreClades=True - will find all clades that are present in
	                              all input trees. Each clade will
	                              be analysed separately and
	                              thereby possibly giving more
	                              information on clades that are
	                              omnipresent but not identical in
	                              all input trees.
        (False by default)
        
        eg
        
            ls = LeafSupport(filename.nex)
            ls.rooted=False
            ls.leafSupport()
        
        for unrooted trees
        and
        
            ls = LeafSupport(filename.nex)
            ls.leafSupport()
        
        for rooted
        
        The resulting leaf support will be presented in a matrix where each taxa have three support values.
        Maximum - the average maximum support for the taxa. Higher is better
        Difference - the average difference betwen the maximum support and the second highest support.  Higher is better
        Entropy - the entropy among the highest, medium and lowest support for the taxa. Lower is better.
        """
        self.rooted=True
        self.exploreClades=False
        self.cladePercentage=100
        self.identicalTaxaSets=identicalTaxaSets
        self.equalDistUnresolved=True
        self.treatUnresolvedAsValid=True
        self.writeCsv=False
        self.csvFilename='leafSupport.csv'
        self.useAllQuartets=True
        self.noQuartetsToUse=0.1 #Proportion of the total amount of possible quartets
        
        # Tobias -- Why is this here?
        #if not self.identicalTaxaSets:
        #    var.allowTreesWithDifferingTaxonSets = True
            
        gm = ['LeafSupport.__init__(filename)']
        self.inputTrees = []
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
        
        self.groups = []
        self.clades = []
        self.taxSets = []
        self.init = False
        
    def _init(self):
        if self.init:
            return
        
        self.init = True

        uniqueSet = set()
        if self.identicalTaxaSets:
            t = self.inputTrees[0]   
            if not t._taxNames:
                t._setTaxNamesFromLeaves() 
            taxnames = t.taxNames
        else:
            taxnames = []
            for t in self.inputTrees:
                if not t._taxNames:
                    t._setTaxNamesFromLeaves()
                for name in t.taxNames:
                    uniqueSet.add(name)
                
            for taxa in uniqueSet:
                taxnames.append(taxa)
            
        taxnames.sort()
        self.taxnames = taxnames

        self.sorted2QuartetType = {}

        self.taxa2IndexDict = {}
        self.index2TaxaDict = {}
        self.bitkeys = []
        self.bitkey2Index = {}
        self.index2TaxaDict[-1] = ROOT_NODE_NAME
        self.taxa2IndexDict[ROOT_NODE_NAME] = -1
        for i in range(len(taxnames)):
            self.taxa2IndexDict[taxnames[i]] = i 
            self.index2TaxaDict[i] = taxnames[i]
            self.bitkeys.append(1L << i)

            self.bitkey2Index[1L << i] = i

        st = SplitStripper()
        for t in self.inputTrees:
            st.buildInformativeSplitsFromTree(t, self.taxnames)
            st.reset()

        self.buildSortedQuartetList()
        
    def printTaxnames(self):
        
        index = 0 
        print 'Index, taxname'
        for name in self.taxnames:
            print '%s, %s' % (index, name)
            index = index + 1
            
    
    def defineGroup(self, list):
        if not self.init:
            self._init()
        gm = ['defineGroup.(list))']
        if type(list) != type([]):
            gm.append('Input must be a list of taxnames and/or taxnumbers')
            raise Glitch, gm
        group = []
        for i in list:
            if self.taxa2IndexDict.has_key(i):
        ##        print 'Taxname: %s, %s' % (i, self.taxa2IndexDict[i])
                group.append(self.taxa2IndexDict[i])
            elif self.index2TaxaDict.has_key(i):
        #        print 'Index: %s, %s' % (self.index2TaxaDict[i], i)
                group.append(i)
            else:
                gm.append('%s is not a tax name or number.' % (i))
                raise Glitch, gm

        #    print ''
        #print group
        group.sort()
        tuple = ()
        for i in group:
            tuple = tuple + (i, )
        self.groups.append(tuple)
      
    def defineTaxSet(self, list):
        if not self.init:
            self._init()
        gm = ['defineTaxSet.(list))']
        if type(list) != type([]):
            gm.append('Input must be a list of taxnames and/or taxnumbers')
            raise Glitch, gm
        taxSet = []
        for i in list:
            if self.taxa2IndexDict.has_key(i):
        #        print 'Taxname: %s, %s' % (i, self.taxa2IndexDict[i])
                taxSet.append(self.taxa2IndexDict[i])
            elif self.index2TaxaDict.has_key(i):
        #        print 'Index: %s, %s' % (self.index2TaxaDict[i], i)
                taxSet.append(i)
            else:
                gm.append('%s is not a tax name or number' % (i))
                raise Glitch, gm
        #    print ''
        #print group
        taxSet.sort()
        tuple = ()
        for i in taxSet:
            tuple = tuple + (i, )
        self.taxSets.append(tuple)
        
    def defineClade(self, list):
        if not self.init:
            self._init()
        gm = ['defineClade.(list))']
        if type(list) != type([]):
            gm.append('Input must be a list of taxnames and/or taxnumbers')
            raise Glitch, gm
        clade = []
        for i in list:
            if self.taxa2IndexDict.has_key(i):
        #        print 'Taxname: %s, %s' % (i, self.taxa2IndexDict[i])
                clade.append(i)
            elif self.index2TaxaDict.has_key(i):
        #        print 'Index: %s, %s' % (self.index2TaxaDict[i], i)
                clade.append(self.index2TaxaDict[i])
            else:
                gm.append('%s is not a tax name or number' % (i))
                raise Glitch, gm

        #    print ''
        #print clade
        clade.sort()
        tuple = ()
        for i in clade:
            tuple = tuple + (i, )
        self.clades.append(tuple)  
    
    def buildSortedQuartetList(self):
        
        noTaxnames = len(self.taxnames)
        
        self.sortedQuartets = []
        
        if self.useAllQuartets: 
            if self.rooted:
                for index in range(noTaxnames-2):
                    for j in range(index+1, noTaxnames-1):
                        for k in range(j+1, noTaxnames):
                            self.sortedQuartets.append((-1,self.bitkeys[index], self.bitkeys[j], self.bitkeys[k]))
            else:
                for index in range(noTaxnames-3):
                    for j in range(index+1, noTaxnames-2):
                        for k in range(j+1, noTaxnames-1):
                            for l in range(k+1, noTaxnames):
                                self.sortedQuartets.append((self.bitkeys[index], self.bitkeys[j], self.bitkeys[k], self.bitkeys[l]))
        else:
            taxon = 4
            if self.rooted:
                taxon = 3
            noOfPossibleQuartets = factorial(noTaxnames) / ( factorial(taxon) * factorial(noTaxnames - taxon) )
            
            quartetsToUse = int(noOfPossibleQuartets*self.noQuartetsToUse)
            
            print 'Using %s quartets out of %s possible' % (quartetsToUse, noOfPossibleQuartets)
#            print 'Quartets to use: ', quartetsToUse
#            print 'Number of taxa:  ',noTaxnames
            
            quartetDict = {}
            if self.rooted:
                i = 0
                while i < quartetsToUse:
                    for index in range(noTaxnames-2):
                        tuple = (-1,self.bitkeys[index])
                
                        taxa = []
                        taxa.append(random.randint(index+1, noTaxnames-1))
                        taxa.append(random.randint(index+1, noTaxnames-1))
            
                        while taxa[1] == taxa[0]:
                            taxa[1] = random.randint(index+1, noTaxnames-1)
                        taxa.sort()
                        tuple += (self.bitkeys[taxa[0]],self.bitkeys[taxa[1]])
                        if not quartetDict.has_key(tuple):
                            quartetDict[tuple] = 1
                            self.sortedQuartets.append(tuple)
                            i += 1
            else:
                i = 0
                while i < quartetsToUse:
                    for index in range(noTaxnames-3):
                        tuple = (self.bitkeys[index],)
                
                        taxa = []
                        taxa.append(random.randint(index+1, noTaxnames-1))
                        taxa.append(random.randint(index+1, noTaxnames-1))
            
                        while taxa[1] == taxa[0]:
                            taxa[1] = random.randint(index+1, noTaxnames-1)
            
                        taxa.append(random.randint(index+1, noTaxnames-1))
                        while taxa[2]== taxa[1] or taxa[2] == taxa[0]:        
                            taxa[2] = random.randint(index+1, noTaxnames-1)
                        taxa.sort()
                        tuple += (self.bitkeys[taxa[0]],self.bitkeys[taxa[1]], self.bitkeys[taxa[2]])
                        if not quartetDict.has_key(tuple):
                            self.sortedQuartets.append(tuple)
                            i += 1

            
    def updateQuartetDicts(self, t, weight, taxa2index):

        for quartet in self.sortedQuartets:
            
            resolutionFound = False
            for split in t.splits:
                
                if self.rooted:
                    hits = 0
                    one = False
                    two = False
                    three = False
                    if quartet[1] & split:
                        hits += 1
                        one = True
                    if quartet[2] & split:
                        hits += 1
                        two = True
                    if quartet[3] & split:
                        hits += 1
                        three = True
                        
                    if hits == 2:
                        if one and two:
                            q = ((quartet[0],self.bitkey2Index[quartet[3]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]]))
                            self.updateDicts(quartet, q, 1, weight, taxa2index)
                            resolutionFound =True
                            break
                        
                        if one and three:
                            q = ((quartet[0],self.bitkey2Index[quartet[2]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[3]]))
                            self.updateDicts(quartet, q, 2, weight, taxa2index)
                            resolutionFound = True
                            break
                        if two and three:
                            q = ((quartet[0],self.bitkey2Index[quartet[1]]),(self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]]))
                            self.updateDicts(quartet, q, 3, weight, taxa2index)
                            resolutionFound = True
                            break
                
                
                else:
                    hits = 0
                    one = False
                    two = False
                    three = False
                    four = False
                    if quartet[0] & split:
                        hits += 1
                        one = True
                    if quartet[1] & split:
                        hits += 1
                        two = True
                    if quartet[2] & split:
                        hits += 1
                        three = True
                    if quartet[3] & split:
                        hits += 1
                        four = True
                    
                    if hits == 2:
                        
                        if one and two or three and four:
                            q = ((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[1]]),(self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]]))
                            self.updateDicts(quartet, q, 1, weight, taxa2index)
                            resolutionFound =True
                            break
                        
                        if one and three or two and four:
                            q = ((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[2]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[3]]))
                            self.updateDicts(quartet, q, 2, weight, taxa2index)
                            resolutionFound = True
                            break
                        
                        if one and four or two and three:
                            q = ((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[3]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]]))
                            self.updateDicts(quartet, q, 3, weight, taxa2index)
                            resolutionFound = True
                            break
                        
            if not resolutionFound:
                if self.rooted:
                    q = ((quartet[0],self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]]),)
                else:
                    q = ((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]]),)
                    
                self.updateDicts(quartet, q, 4, weight, taxa2index)
        
        
            
    def getQuartetsFromTree(self, t, weight, taxa2index):
        
        quartets = []
        for quartet in self.sortedQuartets:
            
            resolutionFound = False
            for split in t.splits:
                
                if self.rooted:
                    hits = 0
                    one = False
                    two = False
                    three = False
                    if quartet[1] & split:
                        hits += 1
                        one = True
                    if quartet[2] & split:
                        hits += 1
                        two = True
                    if quartet[3] & split:
                        hits += 1
                        three = True
                        
                    if hits == 2:
                        if one and two:
                            q = ((quartet[0],self.bitkey2Index[quartet[3]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]]))
                            self.updateDicts(quartet, q, 1, weight, taxa2index)
                            quartets.append(((quartet[0],self.bitkey2Index[quartet[3]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]])))
                            resolutionFound =True
                            break
                        
                        if one and three:
                            q = ((quartet[0],self.bitkey2Index[quartet[2]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[3]]))
                            self.updateDicts(quartet, q, 2, weight, taxa2index)
                            quartets.append(((quartet[0],self.bitkey2Index[quartet[2]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[3]])))
                            resolutionFound = True
                            break
                        if two and three:
                            q = ((quartet[0],self.bitkey2Index[quartet[1]]),(self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]]))
                            self.updateDicts(quartet, q, 3, weight, taxa2index)
                            quartets.append(((quartet[0],self.bitkey2Index[quartet[1]]),(self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]])))
                            resolutionFound = True
                            break
                
                
                else:
                    hits = 0
                    one = False
                    two = False
                    three = False
                    four = False
                    if quartet[0] & split:
                        hits += 1
                        one = True
                    if quartet[1] & split:
                        hits += 1
                        two = True
                    if quartet[2] & split:
                        hits += 1
                        three = True
                    if quartet[3] & split:
                        hits += 1
                        four = True
                    
                    if hits == 2:
                        
                        if one and two or three and four:
                            q = ((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[1]]),(self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]]))
                            self.updateDicts(quartet, q, 1, weight, taxa2index)
                            quartets.append(((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[1]]),(self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]])))
                            resolutionFound =True
                            break
                        
                        if one and three or two and four:
                            q = ((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[2]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[3]]))
                            self.updateDicts(quartet, q, 2, weight, taxa2index)
                            quartets.append(((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[2]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[3]])))
                            resolutionFound = True
                            break
                        
                        if one and four or two and three:
                            q = ((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[3]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]]))
                            self.updateDicts(quartet, q, 3, weight, taxa2index)
                            quartets.append(((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[3]]),(self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]])))
                            resolutionFound = True
                            break
                        
            if not resolutionFound:
                if self.rooted:
                    q = ((quartet[0],self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]]),)
                    quartets.append(((quartet[0],self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]]),))
                else:
                    q = ((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]]),)
                    quartets.append(((self.bitkey2Index[quartet[0]],self.bitkey2Index[quartet[1]],self.bitkey2Index[quartet[2]],self.bitkey2Index[quartet[3]]),))
                
                self.updateDicts(quartet, q, 4, weight, taxa2index)
        
#        noTaxnames = len(self.taxnames)
#        taxon = 4
#        if self.rooted:
#            taxon = 3
#        noOfPossibleQuartets = factorial(noTaxnames) / ( factorial(taxon) * factorial(noTaxnames - taxon) )
#        if len(quartets) != noOfPossibleQuartets:
#            print 'Quartets: ', len(quartets)
#            print 'Possible quartets: ', noOfPossibleQuartets
        
#        print '------------------------------------------------------------------------------------------'
        return quartets
    
    def updateDicts(self, sortedQuartet, quartet, quartetType, weight, taxa2index):
        
#        if self.sorted2QuartetType.has_key(sortedQuartet):
#            print self.sorted2QuartetType[sortedQuartet]
#        else:
#            print 'Sorted quartet not present'
#        print sortedQuartet
#        print quartet
#        print quartetType
        
        if self.sorted2QuartetType.has_key(sortedQuartet):
            if self.sorted2QuartetType[sortedQuartet].has_key(quartetType):
#                print '   1'
                self.sorted2QuartetType[sortedQuartet][quartetType] = self.sorted2QuartetType[sortedQuartet][quartetType] + weight
            else:
#                print '   2'
                self.sorted2QuartetType[sortedQuartet][quartetType] = weight 
        else:
#            print '   3'
            dict = {}
            dict[quartetType] = weight
            self.sorted2QuartetType[sortedQuartet] = dict
        
  
    def leafSupport(self):
        
        if not self.init:
            self._init()
        
    
        if self.exploreClades or len(self.clades) > 0:
            cladeStripper = CladeStripper()
        
        
        #print 'Leaf stabilites for trees from %s' % (filename) 
        #print 'Processing trees'
        print 'Calculating leaf support'
        sys.stdout.flush()
        
        quartet2index = {}
        index2Quartet = {}
        taxa2index = {}
        
        self.sumOfWeights = 0.0
        cladeDict = {}
        
        index = 0

        treeIndex = 0
        sampler = 1
        if len(self.inputTrees) > 20:
            sampler = int(float(len(self.inputTrees)/20))
        sys.stdout.write('0%                    100%\n')
        sys.stdout.write('  ')
        sys.stdout.flush()
        
        for t in self.inputTrees:
            treeIndex += 1
            if (treeIndex + 1) % sampler == 0:
                sys.stdout.write('.')
                sys.stdout.flush()

            t.setPreAndPostOrder()
            weight = 1.0
            if hasattr(t, "weight"):
                if t.weight != None:
                    weight = t.weight
            elif hasattr(t, "recipWeight"):
                if t.recipWeight != None:
                    weight = 1.0/int(t.recipWeight)
            self.sumOfWeights += weight
            #print 'TreeIndex: ', treeIndex
            #print 'Quartets: ', len(quartetStripper.getQuartetSet(False, t, self.taxa2IndexDict))
            #quartetStripper.reset()
            
            if self.exploreClades or len(self.clades) > 0:
                cladeStripper.reset()
                for clade in cladeStripper.getCladesFromTree(t):
                    if len(clade) > 3:
                        if cladeDict.has_key(clade):
                            cladeDict[clade] = cladeDict[clade] + weight
                        else:
                            cladeDict[clade] = weight


            self.updateQuartetDicts(t, weight, taxa2index)
            #self.getQuartetsFromTree(t, weight, taxa2index)                         
                
        sys.stdout.write('\n')
        sys.stdout.flush()
        
        

        scoreDict = self.calcQuickScores(self.sorted2QuartetType, self.treatUnresolvedAsValid, self.equalDistUnresolved)
       
        
        resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning = self.formatResults(scoreDict, self.index2TaxaDict)
        
        results = [[resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning]]
         
        groupResults = []
     
        if self.exploreClades:
            omnipresentClades = []
            for clade, hits in cladeDict.items():
                if hits >= self.sumOfWeights * self.cladePercentage/100.0 :
                    tuple = ()
                    for name in clade:
                        tuple = tuple +(self.taxa2IndexDict[name],)
                        
                    omnipresentClades.append([tuple,(hits/float(self.sumOfWeights))*100])
        
            if len(omnipresentClades) == 0:
                print 'Sorry, no clades matching the %s proportion found in the input trees' % (self.cladePercentage)
      
            omnipresentClades = sorted(omnipresentClades, self._hitsCmp)
        
            for clade in omnipresentClades:
                translatedClade = ()
                for i in clade[0]:
                    translatedClade = translatedClade + (self.index2TaxaDict[i],)
                scoreDict = self.calcLimitedSetScores(self.sorted2QuartetType, index2Quartet, self.treatUnresolvedAsValid, self.equalDistUnresolved, clade[0])
                resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning = self.formatResults(scoreDict, self.index2TaxaDict)  
                results.append([resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning, translatedClade, clade[1]])
        
        if len(self.clades) > 0:
            for clade in self.clades:
                hits = 0.0
                if cladeDict.has_key(clade):
                    hits = (cladeDict[clade]/float(self.sumOfWeights))*100
                translatedClade = ()
                for name in clade:
                    translatedClade = translatedClade +(self.taxa2IndexDict[name],)
                
                scoreDict = self.calcLimitedSetScores(self.sorted2QuartetType, index2Quartet, self.treatUnresolvedAsValid, self.equalDistUnresolved, translatedClade)
                if len(scoreDict.keys()) > 0:
                    resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning = self.formatResults(scoreDict, self.index2TaxaDict)  
                    results.append([resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning, translatedClade, hits])
        
        if len(self.taxSets) > 0:
            for taxSet in self.taxSets:
                scoreDict = self.calcLimitedSetScores(self.sorted2QuartetType, index2Quartet, self.treatUnresolvedAsValid, self.equalDistUnresolved, taxSet)
#                scoreDict = self.calcMemberScores(taxa2index, index2Quartet, taxSet)
            
                if len(scoreDict.keys()) > 0:
                    resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning = self.formatResults(scoreDict, self.index2TaxaDict)  
                    results.append([resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning, taxSet, 0.0])
        

        if len(self.groups) > 0:
            for group in self.groups:
                scoreDict = self.calcMemberScores(self.sorted2QuartetType, index2Quartet, group)
                
                if len(scoreDict.keys()) > 0:
                    member = {}
                    nonMember = {}
                    for i,j in scoreDict.items():
#                        print '%s, %s' % (self.index2TaxaDict[i],j)
                        if i in group:
                            member[i] = j
                        else:
                            nonMember[i] = j
                    
                    resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning = self.formatResults(member, self.index2TaxaDict)
                    nMresultList, nMmaxAverage, nMdiffAverage, nMentAverage, nMdiffWarning, nMentWarning = self.formatResults(nonMember, self.index2TaxaDict)
                    
#                    resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning = self.formatResults(scoreDict, self.index2TaxaDict)  
                    groupResults.append([[resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning, group, 0.0],[nMresultList, nMmaxAverage, nMdiffAverage, nMentAverage, nMdiffWarning, nMentWarning, 0.0]])
#                    results.append([resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning, group, 0.0])
                else:
                    print 'Zero length'
        if self.writeCsv:
            writer = csv.writer(open(self.csvFilename, "wb"))
            list = []
            for result in results:
                if len(result) == 6:
                    list.append(['All taxa'])
                    list.append(['Percentage: ', 100.0])
                    self.csvListFromResult(result, list)
                    list.append([])
                else:
                    if result[7] > 0.0:
                        clade = ['Clade: ']
                        clade.extend(result[6])
                        list.append(clade)
                        list.append(['Percentage: ', result[7]])
                    else:
                        list.append(['Taxon set'])
                    
                    self.csvListFromResult(result, list)
                    list.append([])
                    
            for result in groupResults:
                list.append(['Groupmembership'])
                list.append(['Members'])
                self.csvListFromResult(result[0], list)
                list.append(['Nonmembers'])
                self.csvListFromResult(result[1], list)
                list.append([])
                
            writer.writerows(list)
        
        for result in results:        
##            self.printSupportTable(result[0],result[1], result[2], result[3], result[4], result[5])
            if len(result) == 6:
                print 'All taxa'
                self.printSupportTable(result[0],result[1], result[2], result[3], result[4], result[5])
                print ''
            else:
                if result[7] > 0.0:
                    print 'Clade appears in %s percent of the trees ' % (result[7])
                else:
                    print 'Taxon set'
                self.printSupportTable(result[0],result[1], result[2], result[3], result[4], result[5])
                print ''
        
        for result in groupResults:
            print 'Group membership'
            self.printGroupSupportTable(result)
            print ''
            
            
    def _hitsCmp(self, one, other):
        return int(other[1]) - int(one[1]) 
        
    def calcLimitedSetScores(self, taxa2index, index2Quartet, treatUnresolvedAsValid, equalDistUnresolved, clade):
        scoreDict = {}
        for sortedQuartet, dict in taxa2index.items():
            firstTaxon = -1
            if not self.rooted:
                firstTaxon = self.bitkey2Index[sortedQuartet[0]]
            if self._isCladeContainingQuartet(clade, (firstTaxon,self.bitkey2Index[sortedQuartet[1]],
                                                      self.bitkey2Index[sortedQuartet[2]],self.bitkey2Index[sortedQuartet[3]])):
                quartets = []
                for index, no in dict.items():
                    quartets.append(self.buildQuartet(index, no, sortedQuartet))

                if len(quartets) > 0:

                    max, diff, ent = self.calcAvgDiffEnt(quartets, treatUnresolvedAsValid, equalDistUnresolved)

                    for taxon in sortedQuartet:
                        if taxon >= 0:
                            taxa = self.bitkey2Index[taxon]
                            if scoreDict.has_key(taxa):
                                list = scoreDict[taxa]
                                scoreDict[taxa] = [list[0]+max, list[1]+diff, list[2]+ent, list[3]+1]
                            else:
                                scoreDict[taxa] = [max, diff, ent, 1]
        return scoreDict
            
    def _isCladeContainingQuartet(self, clade, quartet):

        for taxa in quartet:
            if taxa != -1:
                damned = True
                for i in clade:
                    if i == taxa:
                        damned = False
                if damned:
                    return False
        return True
            
    def calcMemberScores(self, taxa2index, index2Quartet, group):
        scoreDict = {}
        for sortedQuartet, dict in self.sorted2QuartetType.items():
            firstTaxon = -1
            if not self.rooted:
                firstTaxon = self.bitkey2Index[sortedQuartet[0]]
            if self._isCladeContainingHalfQuartet(group,(firstTaxon,self.bitkey2Index[sortedQuartet[1]],
                                                         self.bitkey2Index[sortedQuartet[2]],self.bitkey2Index[sortedQuartet[3]])):
                quartets = []
                for index, no in dict.items():
                    q = self.buildQuartet(index, no, sortedQuartet)
                    if self._isGroupHalfQuartet(group, q[0]):
                        quartets.append(q)
                if len(quartets) > 0:
                    max, diff, ent = self.calcAvgDiffEnt(quartets, self.treatUnresolvedAsValid, self.equalDistUnresolved)

                    for taxon in sortedQuartet:
                        if taxon >= 0:
                            taxa = self.bitkey2Index[taxon]
                            if scoreDict.has_key(taxa):
                                list = scoreDict[taxa]
                                scoreDict[taxa] = [list[0]+max, list[1]+diff, list[2]+ent, list[3]+1]
                            else:
                                scoreDict[taxa] = [max, diff, ent, 1]
        return scoreDict
    
    def _isGroupHalfQuartet(self, group, quartet):
        if len(quartet[0]) > 2:
            return False
        first = False
        second = False
        third = False
        fourth = False
        
        for i in group:
            if i == quartet[0][0]:
                first = True
            if i == quartet[0][1]:
                second = True
            if i == quartet[1][0]:
                third = True
            if i == quartet[1][1]:
                fourth = True
                
        if first and second and not (third or fourth):
            return True
        if third and fourth and not (first and second):
            return True
        
        return False
        
            
    def _isCladeContainingHalfQuartet(self, group, quartet):
        
#        print 'group:   ', group
#        print 'quartet: ', quartet
        
        truths = []
        for i in group:
            if i == quartet[0]:
                truths.append(True)
            if i == quartet[1]:
                truths.append(True)
            if i == quartet[2]:
                truths.append(True)
            if i == quartet[3]:
                truths.append(True)
                
        if len(truths) == 2:
            return True

        return False
    
    def buildQuartet(self, quartetType, weight, sortedQuartet):
        q = 0
        firstTaxon = -1
        if not self.rooted:
            firstTaxon = self.bitkey2Index[sortedQuartet[0]]
                
        if quartetType == 1:
            q = (((firstTaxon,self.bitkey2Index[sortedQuartet[1]]),
                              (self.bitkey2Index[sortedQuartet[2]],self.bitkey2Index[sortedQuartet[3]])),weight)
                
        if quartetType == 2:
            q = (((firstTaxon,self.bitkey2Index[sortedQuartet[2]]),
                              (self.bitkey2Index[sortedQuartet[1]],self.bitkey2Index[sortedQuartet[3]])),weight)
                
        if quartetType == 3:
            q = (((firstTaxon,self.bitkey2Index[sortedQuartet[3]]),
                              (self.bitkey2Index[sortedQuartet[1]],self.bitkey2Index[sortedQuartet[2]])),weight)
        if quartetType == 4:
            q= (((firstTaxon,self.bitkey2Index[sortedQuartet[1]],
                               self.bitkey2Index[sortedQuartet[2]],self.bitkey2Index[sortedQuartet[3]]),),weight)
        return q
    
    def calcQuickScores(self, taxa2index, treatUnresolvedAsValid, equalDistUnresolved):
        
        scoreDict = {}
        for sortedQuartet, dict in taxa2index.items():
#            print ((self.index2TaxaDict[self.bitkey2Index[sortedQuartet[0]]],self.index2TaxaDict[self.bitkey2Index[sortedQuartet[1]]]),
#                                      (self.index2TaxaDict[self.bitkey2Index[sortedQuartet[2]]],self.index2TaxaDict[self.bitkey2Index[sortedQuartet[3]]]))
#            print dict
            quartets = []
            for index, no in dict.items():
                
                quartets.append(self.buildQuartet(index, no, sortedQuartet))
                
            if len(quartets) > 0:

                max, diff, ent = self.calcAvgDiffEnt(quartets, treatUnresolvedAsValid, equalDistUnresolved)
                for taxon in sortedQuartet:
                    if taxon >= 0:
                        taxa = self.bitkey2Index[taxon]
                        if scoreDict.has_key(taxa):
                            list = scoreDict[taxa]
                            scoreDict[taxa] = [list[0]+max, list[1]+diff, list[2]+ent, list[3]+1]
                        else:
                            scoreDict[taxa] = [max, diff, ent, 1]
        return scoreDict
    
    def calcScores(self, taxa2index, index2Quartet, treatUnresolvedAsValid, equalDistUnresolved):
        scoreDict = {}
        for sortedQuartet, dict in taxa2index.items():
#            print sortedQuartet
#            print dict
            quartets = []
            for index, no in dict.items():

                quartets.append((index2Quartet[index],no))

            if len(quartets) > 0:

                max, diff, ent = self.calcAvgDiffEnt(quartets, treatUnresolvedAsValid, equalDistUnresolved)

                for taxa in sortedQuartet:
                    if taxa >= 0:
                        if scoreDict.has_key(taxa):
                            list = scoreDict[taxa]
                            scoreDict[taxa] = [list[0]+max, list[1]+diff, list[2]+ent, list[3]+1]
                        else:
                            scoreDict[taxa] = [max, diff, ent, 1]

        return scoreDict
        
    def formatResults(self, scoreDict, index2TaxaDict):
        resultList = []
        
#       Retrive the results from scoreDict 
        for taxa in scoreDict.keys():
            list = scoreDict[taxa]
            list = [round(list[0]/list[3],8),round(list[1]/list[3],8),round(list[2]/list[3],8),list[3]]
            list.append(index2TaxaDict[taxa])
            resultList.append(list)
        
#       Sort results based on Maximum column 
        for i in range(0,len(resultList)):
            for j in range(i+1,len(resultList)):
                if resultList[i][0] < resultList[j][0]:
                    temp = resultList[j]
                    resultList[j] = resultList[i]
                    resultList[i] = temp
                elif resultList[i][0] == resultList[j][0]:
                    if resultList[i][4] > resultList[j][4]:
                        temp = resultList[j]
                        resultList[j] = resultList[i]
                        resultList[i] = temp
        
#       Calc average result based on Maximum
        maxAverage = 0.0
        diffAverage = 0.0
        entAverage = 0.0
        for result in resultList:
            maxAverage += result[0]
            diffAverage += result[1]
            entAverage += result[2]
             
        maxAverage = maxAverage/len(resultList)
        diffAverage = diffAverage/len(resultList)
        entAverage = entAverage/len(resultList)
        
#        Add ranks to the resultlist
#        

        maxList = []
        diffList = []
        entList = []
        
        for result in resultList:
            maxList.append(result[0])
            diffList.append(result[1])
            entList.append(result[2])
        
        maxDict = {}
        maxList.sort(reverse=True)
        index = 0
        maxScore = 1.0
        for score in maxList:
            if score == maxScore:
                if score == 1.0 and index == 0:
                    index += 1
                maxDict[score] = index
            else:
                maxScore = score
                index += 1
                maxDict[score] = index
        
        diffDict = {}
        diffList.sort(reverse=True)
        index = 0
        diffScore = 1.0
        for score in diffList:
            if score == diffScore:
                if score == 1.0 and index == 0:
                    index += 1
                diffDict[score] = index
            else:
                diffScore = score
                index += 1
                diffDict[score] = index
        
        entDict = {}
        entList.sort(reverse=True)
        index = 0
        entScore = 1.0
        for score in entList:
            if score == entScore:
                if score == 1.0 and index == 0:
                    index += 1
                entDict[score] = index
            else:
                entScore = score
                index += 1
                entDict[score] = index
        
        for result in resultList:
            result.append(maxDict[result[0]])
            result.append(diffDict[result[1]])
            result.append(entDict[result[2]])
            
#       Checks if results from Difference and entropy disagre with the Maximum results
#       i.e. if the resulting sorting would differ if based on Difference or Entropy 
#       instead of Maximum
        diff = 10.0
        diffWarning = False
        ent = 10.0
        entWarning = False
        for result in resultList:
            if result[1] > diff + 0.000000000001:
#                print result[1]
#                print diff
#                print ''
                diffWarning = True
            else:
                diff = result[1]
            if result[2] > ent + 0.000000000001:
                entWarning = True
            else:
                ent = result[2]
                
        return resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning
    
    def csvListFromResult(self, result, csvList):
        csvList.append(['','Maximum','','Difference','','Entropy',''])
        csvList.append(['Taxa', 'LSMax', 'Rank', 'LSDiff','Rank', 'LSEnt','Rank'])
        printedAverage = False
        for r in result[0]:
            if result[1] > r[0] and not printedAverage:
                printedAverage = True
                csvList.append(['Average', result[1] ,'-', result[2], '-', result[3], '-'])
            csvList.append([r[4],r[0],r[5],r[1],r[6],r[2],r[7]])
        
        
    def printGroupSupportTable(self, groupResult):
        resultList = groupResult[0][0]
        maxAverage = groupResult[0][1]
        diffAverage = groupResult[0][2]
        entAverage  = groupResult[0][3]
        
        nMresultList = groupResult[1][0]
        nMmaxAverage = groupResult[1][1]
        nMdiffAverage = groupResult[1][2]
        nMentAverage  = groupResult[1][3]
        max = 0
        for result in resultList:
            if len(result[4]) > max:
                max = len(result[4])
        for result in nMresultList:
            if len(result[4]) > max:
                max = len(result[4])
        if max < 10:
            max = 10
        line = ''
        for i in range(max+64):
            line += '-'
        printedAverage = False
        print line
        print '| Members'.ljust(max+2),
        print '| Maximum'.ljust(12),
        print '| Rank'.ljust(6),
        print '| Difference'.ljust(12),
        print '| Rank'.ljust(6),
        print '| Entropy'.ljust(12),
        print '| Rank'.ljust(6),
        print '|'
        print line
        for result in resultList:
            if maxAverage > result[0] and not printedAverage:
                printedAverage = True
                print line
                print '|',
                print 'Average'.ljust(max),
                print '|',
                print str(maxAverage).ljust(10)[:10],
                print '|      |',
                print str(diffAverage).ljust(10)[:10],
                print '|      |',
                print str(entAverage).ljust(10)[:10],
                print '|      |'
                print line
            print '|',
            print str(result[4]).ljust(max),
            print '|',
            print str(result[0]).ljust(10)[:10],
            print '|',
            print str(result[5]).ljust(4)[:4],
            print '|',
            print str(result[1]).ljust(10)[:10],
            print '|',
            print str(result[5]).ljust(4)[:4],
            print '|',
            print str(result[2]).ljust(10)[:10],
            print '|',
            print str(result[5]).ljust(4)[:4],
            print '|'
        print line    
        print '| Non members'.ljust(max+62),
        print '|'
        print line
        
        printedAverage = False
        for result in nMresultList:
            if nMmaxAverage > result[0] and not printedAverage:
                printedAverage = True
                print line
                print '|',
                print 'Average'.ljust(max),
                print '|',
                print str(nMmaxAverage).ljust(10)[:10],
                print '|      |',
                print str(nMdiffAverage).ljust(10)[:10],
                print '|      |',
                print str(nMentAverage).ljust(10)[:10],
                print '|      |'
                print line
            print '|',
            print str(result[4]).ljust(max),
            print '|',
            print str(result[0]).ljust(10)[:10],
            print '|',
            print str(result[5]).ljust(4)[:4],
            print '|',
            print str(result[1]).ljust(10)[:10],
            print '|',
            print str(result[5]).ljust(4)[:4],
            print '|',
            print str(result[2]).ljust(10)[:10],
            print '|',
            print str(result[5]).ljust(4)[:4],
            print '|'
            
        print line
    
    def printSupportTable(self, resultList, maxAverage, diffAverage, entAverage, diffWarning, entWarning):
        max = 0
        for result in resultList:
            if len(result[4]) > max:
                max = len(result[4])
        if max < 10:
            max = 10
        line = ''
        for i in range(max+64):
            line += '-'
        printedAverage = False
        print line
        print '| Taxa'.ljust(max+2),
        print '| Maximum'.ljust(12),
        print '| Rank'.ljust(6),
        print '| Difference'.ljust(12),
        print '| Rank'.ljust(6),
        print '| Entropy'.ljust(12),
        print '| Rank'.ljust(6),
        print '|'
        print line
        for result in resultList:
            if maxAverage > result[0] and not printedAverage:
                printedAverage = True
                print line
                print '|',
                print 'Average'.ljust(max),
                print '|',
                print str(maxAverage).ljust(10)[:10],
                print '|      |',
                print str(diffAverage).ljust(10)[:10],
                print '|      |',
                print str(entAverage).ljust(10)[:10],
                print '|      |'
                print line
            print '|',
            print str(result[4]).ljust(max),
            print '|',
            print str(result[0]).ljust(10)[:10],
            print '|',
            print str(result[5]).ljust(4)[:4],
            print '|',
            print str(result[1]).ljust(10)[:10],
            print '|',
            print str(result[6]).ljust(4)[:4],
            print '|',
            print str(result[2]).ljust(10)[:10],
            print '|',
            print str(result[7]).ljust(4)[:4],
            print '|'
        print line
        print 'Results sorted by Maximum column'
        if diffWarning:
            print '* Note that the sorting in the Difference column differs from the Maximum column'
        if entWarning:
            print '** Note that the sorting in the Entropy column differs from the Maximum column'
       
    def getSortedQuartetTuple(self, quartet):
        list = []
        if len(quartet[0]) > 2:
            for taxa in quartet[0]:
                list.append(taxa)
            list.sort()
            return (list[0], list[1], list[2], list[3])
        else:
            for taxa in quartet[0]:
                list.append(taxa)
            for taxa in quartet[1]:
                list.append(taxa)
            list.sort()
            return (list[0], list[1], list[2], list[3])
        
        
    def calcAvgDiffEnt(self, quartetList, treatUnresolvedAsValid, equalDist=True):
        quartetTypes = {}
#        quartetTypes[quartetList[0]] = 1.0
        weight = 0.0
        found = 0.0
        for qW in quartetList:
#            print qW
            weight += qW[1]
            if len(qW[0]) > 1:
                if quartetTypes.has_key(qW[0]):
                    quartetTypes[qW[0]] = quartetTypes[qW[0]] + qW[1]
                else:
                    quartetTypes[qW[0]] = qW[1]
            else:
#                print '         ',qW
                found += 1.0

        unresolved = found
#        print 'Unresolved: ', unresolved
#        print 'Found: ', found
#        print ''
        
        if len(quartetTypes.values()) > 4:
            print 'Found more than four kinds of quartets for 4 taxa, not good'
#            print quartetTypes.values()
#            print quartetTypes.keys()
            return 0.0, 0.0, 0.0
        
        if unresolved == 0 and len(quartetTypes.values()) == 1:
#            print 'One kind of quartet'
#            print quartetTypes.values()
#            print ''
#            print 'len(quartetList) ', len(quartetList)
#            print 'len(self.inputTrees)', len(self.inputTrees)
            if self.sumOfWeights == weight: 
#                print '%s,%s,%s' % ('1.0', '1.0', '1.0')
#                print ''
                return 1.0, 1.0, 1.0
            else:
#                print '1'
                max = 1.0/(self.sumOfWeights/weight)
                med = 1.0/(self.sumOfWeights/weight) - (1.0 - 1.0/(self.sumOfWeights/(weight)))
                if med <= 0:
                    ent = 1 - (-1 * (max*log(max, 2)))/log(3,2)
                else:
                    ent = 1 - (-1 * (max*log(max, 2) + med*log(med, 2)))/log(3, 2)
                
#                print '%s, %s, %s' % (max, med, ent)
#                print ''
                return max, med, ent
        
        if unresolved > 0 and len(quartetTypes.values()) == 0:
#            print 'Only unresolved quartets'
#            for qW in quartetList:
#                print qW
            if treatUnresolvedAsValid:
                if self.sumOfWeights == weight:
#                    print '%s,%s,%s' % ('1.0', '1.0', '1.0')
#                    print ''
                    return 1.0,1.0,1.0
                else:
#                    print '2'
                    max = 1.0/(self.sumOfWeights/weight)
                    med = 1.0/(self.sumOfWeights/weight) - (1.0 - 1.0/(self.sumOfWeights/(weight)))
                    if med <= 0:
                        ent = 1 - (-1 * (max*log(max, 2)))/log(3,2)
                    else:
                        ent = 1 - (-1 * (max*log(max, 2) + med*log(med, 2)))/log(3, 2)
#                    print '%s, %s, %s' % (max, med, ent)
#                    print ''
                    return max, med, ent
            else:
#                print '3'
#                print '%s,%s,%s' % ('0.333333333333', '0.0', '0.0')
#                print ''
                return 0.333333333333, 0.0, 0.0
            
        
#        print 'More than one value'
        
        val = []
        index = 0
        for i in quartetTypes.values():
            index += 1
            val.append(i)
        for i in range(index, 3):
            val.append(0.0)
        
#        print 'QuartetTypes: '

#        if weight > self.sumOfWeights or weight < self.sumOfWeights:
#            print val
        quartets = self.sumOfWeights
#            print 'weight ', weight
#            print 'self.sumofweights', self.sumOfWeights
#            print 'unresolved: ', unresolved
#            print ''
        
#        if unresolved > 0:
#            print '                       UNRESOlVED'
                
        if unresolved > 0 and not treatUnresolvedAsValid:
            if equalDist:
                dist = unresolved/3.0
#                print 'Dist: %s' % (str(dist))
                for i in range(0, 3):
                    val[i] = val[i] + dist
            else:
                for i in range(0, 3):
                    val[i] = val[i] + unresolved*(val[i]/(quartets-unresolved))
#                    print 'Disto: %s' % (str(unresolved*(val[i]/quartets)))
            val.sort()
            val.append(0.0)
        elif unresolved > 0:
            val.sort()
            val.append(unresolved)
        else:
            val.sort()
            val.append(0.0)

        max = val[2]/quartets
        avg = (val[2] - val[1])/quartets
        med = val[1]/quartets
        min = val[0]/quartets

        if val[2] == 0:
            t1 = 0
        else:
            t1 = max*log(max, 2)
        
        if val[1] == 0:
            t2 = 0
        else:
            t2 = med*log(med, 2)
            
        if val[0] == 0:
            t3 = 0
        else:
            t3 = min*log(min, 2)
        
        if not treatUnresolvedAsValid:
            ent = 1 - (-1 * (t1 + t2 + t3))/log(3, 2)
        else:
            if val[3] == 0:
                t4 = 0
            else:
                t4 = min*log(val[3], 2)
            ent = 1 - (-1 * (t1 + t2 + t3 + t4))/log(4, 2)
#        print '4'
#        print '%s, %s, %s' % (max, med, ent)
#        print ''
        return max, avg, ent
                        
    def matchingQuartets(self, quartet, possible):
        for taxa in possible:
            found = False
            for tuple in quartet:
                if taxa in tuple:
                    found = True
            if not found:
                return False
        return True
        
        
    def allPossibleQuartets(self, taxaList, taxa2IndexDict, rooted):
        possibles = []
        for i in range(0, len(taxaList)-2):
            for j in range(i+1,len(taxaList)-1):
                for k in range(j+1, len(taxaList)):
                    if rooted:
                        possibles.append((taxa2IndexDict[taxaList[i]], taxa2IndexDict[taxaList[j]], taxa2IndexDict[taxaList[k]]))
                    else:
                        for l in range(k+1, len(taxaList)):
                            possibles.append((taxa2IndexDict[taxaList[i]], taxa2IndexDict[taxaList[j]], taxa2IndexDict[taxaList[k]], taxa2IndexDict[taxaList[l]]))
            
        return possibles
        
        
        
        
class TreeStripper(object):
    
    def __init__(self):
        self.set = set()
        self.translatedSet = None
        self.taxa2id = {}
        self.added = 0
        self.index = 0
    
    def reset(self):
        self.set = set()
        self.translatedSet = None
        self.taxa2id = {}
        self.added = 0
        self.index = 0
    
    def printSet(self):
        print 'Added: %s' % (self.added)
        print 'Set.len: %s' % (len(self.set))
        print 'Set'
        for t in self.set:
            print t
        print 'Translation table'
        for i,j in self.taxa2id.items():
            print 'Item: %s, Value: %s' % (i,j)
        if self.translatedSet:
            print 'Translated set'
            for t in self.translatedSet:
                print t    
                
    def getTranslatedSet(self):
        if self.translatedSet == None:
            self.translate()
        return self.translatedSet
        
    def translate(self):
        self.translatedSet = set()
        for t in self.set:
            t = ()
            for i in t:
                t + (self.taxa2id[i])
            self.translatedSet.add(t)
    
    def getSet(self):
        return self.set
    
    def getTranslationDict(self):
        return self.taxa2id
    
    def allPossiblePairs(self, list):
        possibles = []
        for i in range(0,len(list)-1):
            for j in range(i+1, len(list)):
                
                possibles.append([list[i],list[j]])
#        print 'len(possibles): %s' % (len(possibles))
        return possibles
    
    def getSubtreeTaxNames(self, tree, node):
        names = []
        unresolved = []
        if node.isLeaf:
#            print '11'
#            print node.name
            unresolved.append([node.name])
        else:
            list = []
#            if node.leftChild.isLeaf:
#                list.append(node.leftChild.name)
#                list.extend(tree.getAllLeafNames(node))
#            else:
            list.extend(tree.getAllLeafNames(node)) 
            if len(list) > 0:
#                print '121'
#                print list
                names.append(list)
        
        sibling = node.sibling
#        print '222222'
        while sibling:
            if sibling.isLeaf:
#                print '22'
#                print sibling.name
                unresolved.append([sibling.name])
#                print '-22'
            else:
#                print '33'
#                print tree.getAllLeafNames(sibling)
                names.append(tree.getAllLeafNames(sibling))
#                print '-33'
            sibling = sibling.sibling
#        print '3333333'
#        print 'Names: '
#        print names
#        print 'Unresolved: '
#        print unresolved
        return names, unresolved
                   
    def getSubtreeTaxNo(self, tree, node):
        leftChild = node.leftChild
#        list1 = tree.getNodeNumsAbove(leftChild, 1)
        siblings = []
        list1 = []
        if node.isLeaf:
            siblings.append(node.name)
        else:
            if leftChild.isLeaf:
                list1.append(leftChild.name)
                list1.extend(tree.getAllLeafNames(leftChild))
            else:
                list1 = tree.getAllLeafNames(leftChild)
 
        list2 = []
        sibling = node.sibling
        while sibling:
#            list2.extend(tree.getNodeNumsAbove(sibling, 1))
            if sibling.isLeaf:
                print '22'
                print sibling.name
                siblings.append(sibling.name)
            else:
                print '33'
                print tree.getAllLeafNames(sibling)
                list2.extend(tree.getAllLeafNames(sibling))
            sibling = sibling.sibling           
        return list1, list2, siblings
    
    def getSiblingTaxNames(self, tree, node):
        return self.getSubtreeTaxNames(tree, node.leftChild)

class TripletStripper(TreeStripper):
  
    def getTripletSetFromTree(self, tree):
        
        for node in tree.iterLeavesNoRoot():
            self.taxa2id[node.name] = self.index
            self.index += 1
        
        for node in tree.iterInternalsNoRoot():
            if node.parent.leftChild == node:
                names = self.getSubtreeTaxNames(tree, node)
            else:
                names = self.getSiblingTaxNames(tree, node)

            tuples = self.buildTripletsFromLists(names)
            self.added += len(tuples)
            for t in tuples:
                self.set.add(t)
            
            return self.set, self.taxa2id, self.index
    
    
    def getQuartetSet(self,rooted, tree, dict=None):
        self.reset()
#        tree.draw()
        for node in tree.iterNodes():
            if node == tree.root:
#                print node.name
                leafs = []
                for leaf in tree.iterLeavesNoRoot():
                    leafs.append(leaf.name)
#                    print leaf.name
                tuples = self.buildTripletsFromLists([[node.name],leafs], dict, rooted)
#            print 'node.nodeNum: %s' % (node.nodeNum)
            elif node.parent.leftChild == node:
#                print '1A'

                names, unresolved = self.getSubtreeTaxNames(tree, node)

                if len(unresolved) <= 2:
                    names.extend(unresolved)
                    tuples = self.buildTripletsFromLists(names, dict, rooted)
                else:
                    names.extend(unresolved)
                    tuples = self.buildTripletsFromLists(names, dict, rooted)
#                    print 'Building triplets from unresolved'
                    tuples.extend(self.buildTripletsFromUnresolved(unresolved, dict, rooted))
            else:
#                print '2A'
                if not node.isLeaf:
#                    print '2AA'
                    names, unresolved = self.getSiblingTaxNames(tree, node)
#                    print 'len(unresolved) 2'
#                    print len(unresolved)
#                    print unresolved
                    if len(unresolved) <= 2:
                        names.extend(unresolved)
                        tuples = self.buildTripletsFromLists(names, dict, rooted)
                    else:
                        names.extend(unresolved)
                        tuples = self.buildTripletsFromLists(names, dict, rooted)
#                        print 'Building triplets from unresolved'
                        tuples.extend(self.buildTripletsFromUnresolved(unresolved, dict, rooted))
#                tuples = self.buildTripletsFromLists(names)
#                else:
#                    print 'OOOOOOPPPPPS not covered'
            self.added += len(tuples)
            for t in tuples:
                self.set.add(t)
#        list = []
#        for t in self.set:
#            list.append(t)
    
        return self.set
        
    def stripTripletsFromTree(self, tree):
        
        for node in tree.iterLeavesNoRoot():
            self.taxa2id[node.name] = self.index
            self.index += 1
        
        for node in tree.iterNodesNoRoot():
#            print 'node.nodeNum: %s' % (node.nodeNum)
            if node.parent.leftChild == node:
                names, unresolved = self.getSubtreeTaxNames(tree, node)
#                print 'len(unresolved) 1'
#                print len(unresolved)
#                print unresolved
                if len(unresolved) <= 2:
                    names.extend(unresolved)
                    tuples = self.buildTripletsFromLists(names)
                else:
                    tuples = self.buildTripletsFromLists(names)
#                    print 'Building triplets from unresolved'
                    tuples.extend(self.buildTripletsFromUnresolved(unresolved))
            else:
                if not node.isLeaf:
                    names, unresolved = self.getSiblingTaxNames(tree, node)
#                    print 'len(unresolved) 2'
#                    print len(unresolved)
#                    print unresolved
                    if len(unresolved) <= 2:
                        names.extend(unresolved)
                        tuples = self.buildTripletsFromLists(names, dict)
                    else:
                        tuples = self.buildTripletsFromLists(names, dict)
#                        print 'Building triplets from unresolved'
                        tuples.extend(self.buildTripletsFromUnresolved(unresolved))
                tuples = self.buildTripletsFromLists(names)
            self.added += len(tuples)
            for t in tuples:
                self.set.add(t)
    
    def buildTripletsFromUnresolved(self, list, dict=None, rooted=True):
        triplets = []
        if dict:
            dict[ROOT_NODE_NAME] =-1
        if len(list) < 3:
            print 'Unresolved list to short'
        for i in range(0, len(list)-2):
            for j in range(i+1,len(list)-1):
                for k in range(j+1, len(list)):
                    if not rooted:
                        if dict:
                            list1 = []
                            list1.append(dict[list[i][0]])
                            list1.append(dict[list[j][0]])
                            list1.append(dict[list[k][0]])
                            list1.append(dict[ROOT_NODE_NAME])
                            list1.sort()
                            triplets.append(((list1[0], list1[1], list1[2], list1[3]),))
                        else:
#                            print '1'
                            triplets.append(((list[i][0], list[j][0], list[k][0]),))
                    else:
#                        print '2'
                        triplets.append(((list[i][0], list[j][0], list[k][0]),))
#        print triplets
        return triplets
          
    def buildTripletsFromLists(self, list, dict=None, rooted=True):

#        print 'Building triplets from list'
#        print list
        
        tuples = []
        for i in range(0, len(list)-1):
            if len(list[i]) > 1:
                possibles1 = self.allPossiblePairs(list[i])
            else:
                possibles1 = []
            for j in range(i+1, len(list)):
                if len(list[j]) > 1:
                    possibles2 = self.allPossiblePairs(list[j])
                else:
                    possibles2 = []
                
                
                if list[j] != [None]:
#                    print '1'
                    tuples.extend(self.combineLists(possibles1, list[j], dict, rooted))
                
                if list[i] != [None]:
#                    print '2'
                    tuples.extend(self.combineLists(possibles2, list[i], dict, rooted))
                      
        if len(list) >= 3:
#            print 'Building new triplets'
            
            for i in range(0, len(list)-2):
#                print '11'
                for j in range(i+1, len(list)-1):
#                    print '22'
                    for k in range(j+1, len(list)):
#                        print '33'
                        for ii in list[i]:
#                            print '44'
                            for jj in list[j]:
#                                print '55'
                                for kk in list[k]:
#                                    print 'Appending new triplet'
                                    tuples.append(self.buildTripletsFromThree(ii, jj, kk, dict, rooted))
                                     
        return tuples
    
    def buildTripletsFromThree(self, tax1, tax2, tax3, dict, rooted):
        
#        print 'In build from three'
#        print tax1
##        print tax2
#        print tax3
        
        if dict:
            dict[ROOT_NODE_NAME] = -1
        if rooted:
#            print 'Rooted'
            return ((tax1,tax2,tax3),)
        
        if not dict:
#            print 'No dict'
            return ((tax1, tax2, tax3, ROOT_NODE_NAME),)
        
        list = []
        list.append(dict[tax1])
        list.append(dict[tax2])
        list.append(dict[tax3])
        list.append(dict[ROOT_NODE_NAME])
        list.sort()
        
        return ((list[0], list[1], list[2], list[3]),)
    
                            
    def combineLists(self, list1, list2, dict, rooted):
        tuples = []
        if dict:
            dict[ROOT_NODE_NAME] =-1
#        print list1
#        print list2
        for i in list1:
            for j in list2:
                if not rooted:
                    if dict:
                        list = []
                        list.append(dict[i[0]])
                        list.append(dict[i[1]])
                        list.sort()
                        list1 = []
#                        print 'J: ', j
                        
                        list1.append(dict[j])
                        list1.append(dict[ROOT_NODE_NAME])
                        list1.sort()
                        tuples.append(((list[0],list[1]),(list1[0],list1[1])))
                    else:
#                        print '3'
                        tuples.append(((i[0],i[1]),(j[0],ROOT_NODE_NAME)))
                else:
#                    print '4'
                    tuples.append(((i[0],i[1]),j[0]))
#        print 'len(combined): %s' % (len(tuples))
        return tuples
        
class QuartetStripper(TreeStripper):
    
    def getQuartetSetFromTree(self, tree):
        ss = SplitStripper()
        splitSet, self.taxa2id, self.index = ss.getSplitsFromTree(tree)
        for s in splitSet:
            possibles1 = self._allPossiblePairs(s[0])
            possibles2 = self._allPossiblePairs(s[1])
            list = self._combineLists(possibles1, possibles2)
            for q in list:
                self.set.add(q)
        return self.set
    
    
    def getQuartetSet(self, extra, tree, dict=None):
#        set = set('')
        self.reset()
        quartets = []
        taxaNames = []
        for node in tree.iterLeavesNoRoot():
            taxaNames.append(node.name)
        
        for node in tree.iterNodesNoRoot():
            if node.parent.leftChild == node:
                names, unresolved = self.getSubtreeTaxNames(tree, node)
#                print '0'
#                print 'names: ', names
#                print 'unresolved: ', unresolved
                
                unresolvedLen = len(unresolved)
                if unresolvedLen <= 2 and unresolvedLen > 0:
                    names.extend(unresolved)
                    unresolved = []
#                    print 'Added unresolved to names'

#                print 'Names: ', names
#                print 'Unresolved: ', unresolved
                
                added = []
                
#                if unresolvedLen == 3:
#                    print '1'
#                    added = self.buildQuartetsFromUnresolvedTrio(unresolved, taxaNames)
#                    quartets.extend(added)
#                    print 'Quartets: ', quartets
#                    print 'Combine with all other taxaNames to make unresolved quartets'

                if unresolvedLen >= 3:
#                    print '2'
                    added = self.buildQuartetsFromUnresolvedTaxa(unresolved, taxaNames)
                    quartets.extend(added)
#                    print 'Quartets: ', quartets
#                    print 'Create all possible trios and combine them with all other taxaN'
#                    print 'to make unresolved quartets'
#                    print 'Then create all possible unresolved quartets using the taxaNames in'
#                    print 'unresolved i.e. all possible quartets'
                
#                if len(names) <= 1:
#                    print 'len(names) <= 1 so we do nothing about names'
                    
                if len(names) == 2:
                    tempNames = []
                    tempNames.extend(names[0])
                    tempNames.extend(names[1])
                    otherTaxa = []
                    for taxaname in taxaNames:
                        include = True
                        for name in tempNames:
                            if taxaname == name:
                                include = False
                                break
                        if include:
                            otherTaxa.append(taxaname)
#                    print '3.0'
#                    print 'tempNames: ',tempNames
#                    print 'otherTaxa: ',otherTaxa
                    added = self.buildQuartetsFrom2Lists(tempNames, otherTaxa)
                    quartets.extend(added)
#                    print 'Quartets: ', quartets
                    
                    if len(names[0]) > 1 or len(names[1]) > 1:
                        for i in range(0,len(names)):
                            if len(names[i]) > 1:
                                otherTaxa = []
                                for j in range(0,len(names)):
                                    if j != i:
                                        for taxNames in names[j]:
                                            otherTaxa.append(taxNames)
#                            print '3'
#                            print 'names[i]: ',names[i]
#                            print 'otherTaxa: ',otherTaxa
                            added = self.buildQuartetsFrom2Lists(names[i], otherTaxa)
                            quartets.extend(added)
                        
                if len(names) > 2:
                    tempNames = []
                    tempNames.extend(names[0])
                    tempNames.extend(names[1])
                    otherTaxa = []
                    for taxaname in taxaNames:
                        include = True
                        for name in tempNames:
                            if taxaname == name:
                                include = False
                                break
                        if include:
                            otherTaxa.append(taxaname)
#                    print '3.0'
#                    print tempNames
#                    print otherTaxa
                    added = self.buildQuartetsFrom2Lists(tempNames, otherTaxa)
                    quartets.extend(added)
#                    print 'Quartets: ', quartets
                    for i in range(0,len(names)):
                        if len(names[i]) > 1:
                            otherTaxa = []
                            for j in range(0,len(names)):
                                if j != i:
                                    for taxNames in names[j]:
                                        otherTaxa.append(taxNames)
#                            print '4'
                            added = self.buildQuartetsFrom2Lists(names[i], otherTaxa)
                            quartets.extend(added)
                
#                print ''
#            print added
        
        quartets = self.translateAndSortQuartets(quartets, dict)
        for quartet in quartets:
            self.set.add(quartet)
#        print set
        return self.set

    def translateAndSortQuartets(self, quartets, dict):
        translatedQuartets = []
        if dict:
            for quartet in quartets:
                if len(quartet[0]) == 2:
                    translatedQuartets.append(((dict[quartet[0][0]],dict[quartet[0][1]]),(dict[quartet[1][0]],dict[quartet[1][1]])))
                else:
                    translatedQuartets.append(([dict[quartet[0][0]],dict[quartet[0][1]],dict[quartet[0][2]],dict[quartet[0][3]]],))
        else:
            translatedQuartets = quartets
            
        for i in range(0, len(translatedQuartets)):
            translatedQuartets[i] = self.sortQuartet(translatedQuartets[i])
        return translatedQuartets
            
    def sortQuartet(self, quartet):
        q1 = quartet[0]
        if len(quartet[0]) == 2:
            q2 = quartet[1]
            if q1[0] > q1[1]:
                q1 = (q1[1],q1[0])
        
            if q2[0] > q2[1]:
                q2 = (q2[1], q2[0])
        
            if (q2[0] > q1[0] and q2[0] > q1[1]) or (q2[1] > q1[0] and q2[1] > q1[1]):
                return (q1, q2)
            else:
                return (q2, q1)
        else:
            q1.sort()
            return (((q1[0],q1[1],q1[2],q1[3]),))
    
    def buildQuartetsFromUnresolvedTaxa(self, unresolved, allTaxa):
        if len(unresolved) < 3:
            print 'Not a happy camper, less than four taxa in unresolved'
        
#        print unresolved
        unique = []
        for taxa in allTaxa:
            add = True
            for t in unresolved: 
                if taxa == t[0]:
                    add = False
            if add:
                unique.append(taxa)
            
        quartets = []
#        for i in range(0, len(unresolved)):
#            for j in range(i+1, len(unresolved)):
#                for k in range(j+1, len(unresolved)):
#                    quartets.extend(self.buildQuartetsFromUnresolvedTrio([unresolved[i],unresolved[j],unresolved[k]], allTaxa))
        for i in range(0, len(unresolved)):
            for j in range(i+1, len(unresolved)):
                for k in range(j+1, len(unresolved)):
                    for l in unique:
                        quartets.append(([unresolved[i][0],unresolved[j][0],unresolved[k][0],l],))
#        print quartets
#        print ''
        return quartets
    
    def buildQuartetsFromUnresolvedTrio(self, trio, allTaxa):
        t1 = trio[0][0]
        t2 = trio[1][0]
        t3 = trio[2][0]
        newTrio = (t1,t2,t3)
        quartets = []
        for taxa in allTaxa:
            if not taxa in newTrio:
                quartets.append(([t1,t2,t3,taxa],))
        return quartets
        

    def buildQuartetsFrom2Lists(self, taxa, otherTaxa):
        taxaPairs = []
        for i in range(0, len(taxa)):
            for j in range(i+1, len(taxa)):
                taxaPairs.append((taxa[i],taxa[j]))
        otherTaxaPairs = []
        for i in range(0, len(otherTaxa)):
            for j in range(i+1, len(otherTaxa)):
                otherTaxaPairs.append((otherTaxa[i],otherTaxa[j]))
        quartets = []
        for pair in taxaPairs:
            for otherPair in otherTaxaPairs:
                quartets.append((pair, otherPair))
        
        return quartets
        

    def getQuartetSetTest(self, dict, rooted, tree):
        for node in tree.iterNodesNoRoot():
#            print 'node.nodeNum: %s' % (node.nodeNum)
            if node.parent.leftChild == node:
                names, unresolved = self.getSubtreeTaxNames(tree, node)
#                print 'len(unresolved) 1'
#                print len(unresolved)
#                print unresolved
                if len(unresolved) <= 2:
                    names.extend(unresolved)
                    print 'Names: '
                    print names
                    tuples = self.buildQuartetsFromLists(names, dict)
                else:
                    print 'Names: '
                    print names
                    tuples = self.buildQuartetsFromLists(names, dict)
#                    print 'Building triplets from unresolved'
                    print 'Unresolved: '
                    print unresolved
                    tuples.extend(self.buildQuartetsFromUnresolved(unresolved, dict))
            else:
                if not node.isLeaf:
                    names, unresolved = self.getSiblingTaxNames(tree, node)
#                    print 'len(unresolved) 2'
#                    print len(unresolved)
#                    print unresolved
                    if len(unresolved) <= 2:
                        names.extend(unresolved)
                        print 'Names: '
                        print names
                        tuples = self.buildQuartetsFromLists(names, dict)
                    else:
                        print 'Names: '
                        print names
                        tuples = self.buildQuartetsFromLists(names, dict)
#                        print 'Building triplets from unresolved'
                        print 'Unresolved: '
                        print unresolved
                        tuples.extend(self.buildQuartetsFromUnresolved(unresolved, dict))
#                tuples = self.buildTripletsFromLists(names)
            self.added += len(tuples)
            for t in tuples:
                self.set.add(t)
        list = []
        for t in self.set:
            list.append(t)
    
        return list
    
    def buildQuartetsFromLists(self, list, dict):
        quartets = []
        for i in range(0, len(list)):
            for j in range(i+1, len(list)):
                possibles1 = self._allPossiblePairs(list[i])
                possibles2 = self._allPossiblePairs(list[j])
                quartets.extend(self._combineLists(possibles1, possibles2))
        print quartets
        return quartets
    
    def buildQuartetsFromUnresolved(self, list, dict):
        quartets = []
        if len(list) < 4:
            print 'Unresolved list to short'
        for i in range(0, len(list)-2):
            for j in range(i+1,len(list)-1):
                for k in range(j+1, len(list)):
                    for l in range(k+1, len(list)):
                        list = []
                        list.append(dict[list[i]])
                        list.append(dict[list[j]])
                        list.append(dict[list[k]])
                        list.append(dict[list[l]])
                        list.sort()
                        quartets.append(((list[0], list[1], list[2], list[3]),))
        print quartets
        return quartets
    
    def combineLists(self, list1, list2, dict):
        tuples = []
        for i in list1:
            for j in list2:

                list = []
                list.append(dict[i[0]])
                list.append(dict[i[1]])
                list.sort()
                list1 = []
#                        print j[0]
                list1.append(dict[j][0])
                list1.append(dict[j][1])
                list1.sort()
                tuples.append(((list[0],list[1]),(list1[0],list1[1])))
#        print 'len(combined): %s' % (len(tuples))
        return tuples
    
    def stripQuartetSetFromTree(self, tree):
        ss = SplitStripper()
        splitSet, self.taxa2id, self.index = ss.getSplitsFromTree(tree)
        for s in splitSet:
            possibles1 = self._allPossiblePairs(s[0])
            possibles2 = self._allPossiblePairs(s[1])
            list = self._combineLists(possibles1, possibles2)
            for q in list:
                self.set.add(q)
    
    def getQuartetSetFromSplits(self, splits):
        for s in splits:
            possibles1 = self._allPossiblePairs(s[0])
            possibles2 = self._allPossiblePairs(s[1])
            list = self._combineLists(possibles1, possibles2)
            for q in list:
                self.set.add(q)
        return self.set
        
    def _allPossiblePairs(self, list):
        if len(list) < 2:
            print 'No possible pairs, uninformative split'
        if len(list) ==2:
            return [(list[0], list[1])]
        possibles = []
        for i in range(0,len(list)-1):
            for j in range(i+1, len(list)):
                possibles.append((list[i],list[j]))
#        print 'len(possibles): %s' % (len(possibles))
        return possibles
    
    def _combineLists(self, list1, list2):
        combined = []
        for i in list1:
            if len(i) == 2:
                for j in list2:
                    if len(j) == 2:
                        if i < j:
                            combined.append((i,j))
                        else:
                            combined.append((j,i))
        return combined

class CladeStripper(TreeStripper):
    
    def getCladesFromTree(self, tree):
        for node in tree.iterNodesNoRoot():
            if node.parent.leftChild == node:
                names, unresolved = self.getSubtreeTaxNames(tree, node)
                names.extend(unresolved)
#                print 'names'
#                print names
#                print 'unresolved'
#                print unresolved
#                
                tuples = []
                for i in names:
                    if len(i) > 1:
                        list = []
                        for name in i:
                            list.append(name)
                        list.sort()
                        tuple = ()
                        for name in list:
                            tuple = tuple + (name,)
                        tuples.append(tuple)

            else:
#                print '2A'
                if not node.isLeaf:
#                    print '2AA'
                    names, unresolved = self.getSiblingTaxNames(tree, node)
                    names.extend(unresolved)
#                    print 'names'
#                    print names
#                    print 'unresolved'
#                    print unresolved
                    tuples = []
                    for i in names:
                        if len(i) > 1:
                            list = []
                            for name in i:
                                list.append(name)
                            list.sort()
                            tuple = ()
                            for name in list:
                                tuple = tuple + (name,)
                            tuples.append(tuple)
                
            for t in tuples:
                self.set.add(t)
        list = []
        for t in self.set:
            list.append(t)
    
        return list
        
         

    
class SplitStripper(TreeStripper):
    
    def getSplitsFromTree(self, tree):
        
        for node in tree.iterLeavesNoRoot():
            self.taxa2id[node.name] = self.index
            self.index += 1
        
        seT = set()
        if not tree._taxNames:
            tree._setTaxNamesFromLeaves()
        for name in tree.taxNames:
                seT.add(name)
        list = []
        for name in seT:
            list.append(name)
        list.sort()
        dict = {}
        rdict = {}
        self.bitkeys = []
        for i in range(len(list)):
            self.bitkeys.append(1L << i)
            dict[list[i]] = 1L << i 
            rdict[1L << i] = list[i]
            
        for n in tree.iterLeavesNoRoot():
            n.br.rc = dict[n.name]
        
        tree._makeRCSplitKeys()
        tree.splits =[]
        for n in tree.iterInternalsNoRoot():
            tree.splits.append(n.br.rc)
           
        for s in tree.splits:
            tuple1 = ()
            tuple2 = ()
            for i in self.bitkeys:
                if i & s:
                    tuple1 = tuple1 + (rdict[i],)
                else:
                    tuple2 = tuple2 + (rdict[i],)
                    
            list = []
            list.append(tuple1)
            list.append(tuple2)
            list.sort()
            
            self.set.add((list[0], list[1]))
        
        return self.set, self.taxa2id, self.index
        
        
    def buildInformativeSplitsFromTree(self, tree, taxnames):
        
        dict = {}
        rdict = {}
        self.bitkeys = []
        for i in range(len(taxnames)):
            self.bitkeys.append(1L << i)
            dict[taxnames[i]] = 1L << i 
            rdict[1L << i] = taxnames[i]
            
        for n in tree.iterLeavesNoRoot():
            n.br.rc = dict[n.name]
        
        tree._makeRCSplitKeys()
        tree.splits =[]
        for n in tree.iterInternalsNoRoot():
            tree.splits.append(n.br.rc)
        
        return self.bitkeys
        
        if 0:   
            splits = []
            for s in tree.splits:
                tuple1 = ()
                for i in self.bitkeys:
                    if i & s:
                        tuple1 = tuple1 + (rdict[i],)

                splits.append(tuple1)
        
            for i in splits:
                print i

        
            return splits
        
        
