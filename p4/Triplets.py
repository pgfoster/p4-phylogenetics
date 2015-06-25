import string
from Tree import Tree
from Node import Node


class Aho(object):


    def isTripletCompatibleWithSet(self, set, triplet):
        list = []
        list.append(triplet)
        for t in set:
            list.append(t)
        return self.isListCompatible(list)

    def isSetCompatible(self, set):
        list = []
        for t in set:
            list.append(t)
        return self.isListCompatible(list)
    
    def isListCompatible(self, list):
        if len(list) <= 1:
            return True
        
        components = self._buildComponents(list)
        
        if len(components) == 1 and components[0][0] <= 2:
            return True
        elif len(components) <= 1 and components[0][0] > 2:
            return False
        
        for c in components:
            if not self.isListCompatible(c[1]):
                return False
        
        return True
    
        
    def _buildComponents(self, list):
        
        components = [Component(list[0])]
        
        for i in range(1, len(list)):
            join = []
            next = []
            for c in components:
                if c.contains(list[i]):
                    c.add(list[i])
                    join.append(c)
                else:
                    next.append(c)
                    
            if len(join) >= 1:
                joined = join[0]
                for i in range(1,len(join)):
                    joined.join(join[i])
                components = []
                components.extend(next)
                components.append(joined)
            else:
                components.append(Component(list[i]))
        
        tripletLists = []
        for c in components:
            tripletLists.append([len(c.dict.keys()),c.getTriplets()])
#            print c.getTriplets()
                  
        return tripletLists
        
        
class Component(object):
    
    def __init__(self, triplet):
        self.dict = {}
        self.list = []
        self.add(triplet)
        
    def add(self, triplet):
        self.list.append(triplet)
        self.dict[triplet[0]] = 1
        self.dict[triplet[1]] = 1
        
    def join(self, component):
        self.list.extend(component.list)
        for k in component.dict.keys():
            self.dict[k] = 1
             
    def contains(self, triplet):
        if self.dict.has_key(triplet[0]):
            return True
        if self.dict.has_key(triplet[1]):
            return True
        return False
    
    def getTriplets(self):
        filteredList = []
        for t in self.list:
            if t[2] in self.dict:
                filteredList.append(t)
        return filteredList

class TripletStripper(object):
    
    def __init__(self):
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
        
    def getTripletSet(self):
        return self.set
    
    def getTranslationDict(self):
        return self.taxa2id
    
    def getTranslatedTripletSet(self):
        if self.translatedSet == None:
            self.translateTriplets()
        return self.translatedSet
        
    def translateTriplets(self):
        self.translatedSet = set()
        for t in self.set:
            self.translatedSet.add((self.taxa2id[t[0]], self.taxa2id[t[1]], self.taxa2id[t[2]]))
        
    def stripTripletsFromTree(self, tree):
        
        for node in tree.iterLeavesNoRoot():
            self.taxa2id[node.name] = self.index
            self.index += 1
        
        for node in tree.iterInternalsNoRoot():
#            print 'node.nodeNum: %s' % (node.nodeNum)
            if node.parent.leftChild == node:
                names = self.getSubtreeTaxNames(tree, node)
            else:
                names = self.getSiblingTaxNames(tree, node)
            
#            prinT = False
#            for list in names:
#                if len(list) >= 2:
#                    prinT = True
#            if prinT:
#                print names
            tuples = self.buildTripletsFromLists(names)
            self.added += len(tuples)
            for t in tuples:
                self.set.add(t)
          
    def buildTripletsFromLists(self, list):

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
                tuples.extend(self.combineLists(possibles1, list[j]))
                tuples.extend(self.combineLists(possibles2, list[i]))
                
                    
        return tuples
                            
    def allPossiblePairs(self, list):
        possibles = []
        for i in range(0,len(list)-1):
            for j in range(i+1, len(list)):
                possibles.append([list[i],list[j]])
#        print 'len(possibles): %s' % (len(possibles))
        return possibles
        
    def combineLists(self, list1, list2):
        tuples = []
        for i in list1:
            for j in list2:
                tuples.append((i[0],i[1],j[0]))
#        print 'len(combined): %s' % (len(tuples))
        return tuples
        
    def getSubtreeTaxNames(self, tree, node):
        names = []
        if node.isLeaf:
            names.append([node.name])
        else:
            list = []
            if node.leftChild.isLeaf:
                list.append(node.leftChild.name)
                list.extend(tree.getAllLeafNames(node.leftChild))
            else:
                list.extend(tree.getAllLeafNames(node.leftChild))
            if len(list) > 0:
                names.append(list)
        sibling = node.sibling
        while sibling:
            if sibling.isLeaf:
                names.append([sibling.name])
            else:
                names.append(tree.getAllLeafNames(sibling))
            sibling = sibling.sibling
        return names
                   
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
                siblings.append(sibling.name)
            else:
                list2.extend(tree.getAllLeafNames(sibling))
            sibling = sibling.sibling           
        return list1, list2, siblings
    
    def getSiblingTaxNames(self, tree, node):
        return self.getSubtreeTaxNames(tree, node.leftChild)
                             