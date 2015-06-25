
class QuartetSet(object):
       
    def __init__(self):
        self.quartetLogic = QuartetLogic()
        self.quartetSet = DualIndexNoDupes()
        self.closureSet = DualIndexNoDupes()
        self.conflictingSet = DualIndexNoDupes()
        
        self.compatibleSet = DualIndexNoDupes()
        self.compatibleClosure = DualIndexNoDupes()
        self.incompatibleSet = DualIndexNoDupes()
        
        self.taxNames = []
        self.bitKeys = []
    
    def addQuartet(self, l1, l2, r1, r2):
        if not l1 in self.taxNames:
            self.taxNames.append(l1)
            self.bitKeys.append(1L << self.taxNames.index(l1))
        if not l2 in self.taxNames:
            self.taxNames.append(l2)
            self.bitKeys.append(1L << self.taxNames.index(l2))
        if not r1 in self.taxNames:
            self.taxNames.append(r1)
            self.bitKeys.append(1L << self.taxNames.index(r1))
        if not r2 in self.taxNames:
            self.taxNames.append(r2)
            self.bitKeys.append(1L << self.taxNames.index(r2))
        
#        print 'l1: %s, l2: %s, r1: %s, r2: %s' % (self.taxNames.index(l1), self.taxNames.index(l2), self.taxNames.index(r1), self.taxNames.index(r2))
        
        quartet = Quartet((1L << self.taxNames.index(l1) | 1L << self.taxNames.index(l2)),
                          (1L << self.taxNames.index(r1) | 1L << self.taxNames.index(r2)))
        
        if quartet.isValidQuartet(self.bitKeys):
            compatible, q = self.checkCompatibility(quartet)
            if self.checkCompatibility(quartet):
                self.quartetSet.addQuartet(quartet)
                self.closureSet.addQuartet(quartet)
                self._closure(self.closureSet)       
                return
            else:
                print 'Quartet incompatible with set'
                return
        
        print 'Not a valid quartet, it was not added to the set'
    
    def checkCompatibility(self, quartet):
        list = self.quartetSet.list()
        for q in list:
            if self.quartetLogic.doesQuartetsContradict(quartet, q, self.bitKeys):
                return False, 
        list = self.closureSet.list()
        for q in list:
            if self.quartetLogic.doesQuartetsContradict(quartet, q, self.bitKeys):
                return False
        return True
    
    def checkCompatibleSet(self, quartet):
        list = self.compatibleSet.list()
        compatible = True
        offending = []
        for q in list:
            if self.quartetLogic.doesQuartetsContradict(quartet, q, self.bitKeys):
                compatible = False
                offending.append(q)
        return compatible, offending

    def checkIncompatibleSet(self, closure):
        list = self.closure.list()
        compatible = True
        offending = []
        for q in list:
            if self.incompatibleSet.contains(q.left, q.right):
                compatible = False
                offending.append(q)
        return compatible, offending
    
    def removeFromCompatibleSet(self, list):
        for q in list:
            return
            
    
    def sortQuartetSet(self):
        sortedTaxaNames = sorted(self.taxNames)
        sortedBitkeys = []
        for taxaName in sortedTaxaNames:
            sortedBitkeys.append(1L << sortedTaxaNames.index(taxaName))
        quartetList = self.quartetSet.list()
        sortedQuartetSet = DualIndexNoDupes()
        for q in quartetList:
            quartetTaxNames = self._getListOfTaxNames(q.left, q.right)
#            print 'List: %s' % (quartetTaxNames)
            sortedQuartetSet.addQuartet(Quartet(1L << sortedTaxaNames.index(quartetTaxNames[0]) | 
                                                1L << sortedTaxaNames.index(quartetTaxNames[1]),
                                                1L << sortedTaxaNames.index(quartetTaxNames[2]) |
                                                1L << sortedTaxaNames.index(quartetTaxNames[3])))
        self.taxNames = sortedTaxaNames
        self.bitKeys = sortedBitkeys
        self.quartetSet = sortedQuartetSet
    
    
    def _getListOfTaxNames(self, left, right):
        indeces = []
        indeces.append(self.taxNames[self._getIndex1(left)])
        indeces.append(self.taxNames[self._getIndex2(left)])
        indeces.append(self.taxNames[self._getIndex1(right)])
        indeces.append(self.taxNames[self._getIndex2(right)])
        return indeces
    
    def _getIndex1(self, number):
        for n in self.bitKeys:
            if number & n != 0:
                return self.bitKeys.index(n)
        
    def _getIndex2(self, number):
        first = True
        for n in self.bitKeys:
            if number & n != 0:
                if not first:
                    return self.bitKeys.index(n)
                if first:
                    first = False
    
    def rule1Closure(self):
        lenght = 0
        while lenght < len(self.quartetSet.list()):
            qList = self.quartetSet.list()
            lenght = len(qList)
            for i in range(0, len(qList)-1):
                for j in range(i+1,len(qList)):
                    quartet = self.quartetLogic.rule1a(qList[i], qList[j], self.bitKeys)
#                    quartet = qList[i].rule1a(qList[j], self.bitKeys)
                    if isinstance(quartet, Quartet):
                        self.quartetSet.addQuartet(quartet)
                    quartet = self.quartetLogic.rule1b(qList[i], qList[j], self.bitKeys)
#                    quartet = qList[i].rule1b(qList[j], self.bitKeys)
                    if isinstance(quartet, Quartet):
                        self.quartetSet.addQuartet(quartet)
                    

    def rule2Closure(self):
        lenght = 0
        while lenght < len(self.quartetSet.list()):
            qList = self.quartetSet.list()
            lenght = len(qList)
            for i in range(0, len(qList)-1):
                for j in range(i+1,len(qList)):
                    quartet = self.quartetLogic.rule2(qList[i], qList[j], self.bitKeys)
#                    quartet = qList[i].rule2(qList[j], self.bitKeys)
                    if isinstance(quartet, Quartet):
                        self.quartetSet.addQuartet(quartet)
                                           
    def closure(self):
        lenght = 0
        while lenght < len(self.quartetSet.list()):
            qList = self.quartetSet.list()
            lenght = len(qList)
            for i in range(0, len(qList)-1):
                for j in range(i+1,len(qList)):
                    quartet = self.quartetLogic.rule1a(qList[i], qList[j], self.bitKeys)
#                    quartet = qList[i].rule1a(qList[j], self.bitKeys)
                    if isinstance(quartet, Quartet):
                        self.quartetSet.addQuartet(quartet)
                    quartet = self.quartetLogic.rule1b(qList[i], qList[j], self.bitKeys)
#                    quartet = qList[i].rule1b(qList[j], self.bitKeys)
                    if isinstance(quartet, Quartet):
                        self.quartetSet.addQuartet(quartet)    
                    quartet = self.quartetLogic.rule2(qList[i], qList[j], self.bitKeys)
#                    quartet = qList[i].rule2(qList[j], self.bitKeys)
                    if isinstance(quartet, Quartet):
                        self.quartetSet.addQuartet(quartet)
    
    def _closure(self, set):
        lenght = 0
        while lenght < len(set.list()):
            qList = set.list()
            lenght = len(qList)
            for i in range(0, len(qList)-1):
                for j in range(i+1,len(qList)):
                    quartet = self.quartetLogic.rule1a(qList[i], qList[j], self.bitKeys)
                    if isinstance(quartet, Quartet):
                        set.addQuartet(quartet)
                    quartet = self.quartetLogic.rule1b(qList[i], qList[j], self.bitKeys)
                    if isinstance(quartet, Quartet):
                        set.addQuartet(quartet)    
                    quartet = self.quartetLogic.rule2(qList[i], qList[j], self.bitKeys)
                    if isinstance(quartet, Quartet):
                        set.addQuartet(quartet)
    
    def printQuartets(self):
        index = ''
        for n in self.taxNames:
            index += n
        print index
        for q in self.quartetSet.list():
            q.printQ(self.bitKeys)
    
    
#    The DualIndexNoDupes class holds objects based on a dual index system, ie two indeces are required 
#    for each object. If another object is stored using the same two indeces the one priorly stored
#    will be replaced. 
#    In the QuartetSet case this class is used to avoid the problem of hashing the long ints which is required 
#    if we want to store the quartets in a python set. For the set to contain only unique quartets there would have
#    to be a hash that could turn two long ints into a normal int. 
    
class DualIndexNoDupes(object):
    
    def __init__(self):
        self.dict = {}

#    Method specific for the Quartet class case
    def addQuartet(self, quartet):
        if isinstance(quartet, Quartet):
            self.add(quartet.left, quartet.right, quartet)
    
    def add(self, i1, i2, object):
        if self.dict.has_key(i1):
            self.dict[i1][i2] = object
        else:
            d = {}
            d[i2] = object
            self.dict[i1] = d
            
    def remove(self, i1, i2):
        if i1 in self.dict:
            d = dict[i1]
            if i2 in d:
                del d[i2]
    
    def contains(self, i1, i2):
        if i1 in self.dict:
            d = dict[i1]
            if i2 in d:
                return True
        return False
            
#    Returns a list of all objects stored 
    def list(self):
        list = []
        for d in self.dict.values():
            list.extend(d.values())
        return list
    
    def printD(self):
        list = self.list()
        for d in list:
            print '%s' % (d)



class QuartetLogic(object):
    
    def doesQuartetsContradict(self, one, other, bitkeys):
        if self.ruleOf3(one, other, bitkeys):
            if one.left & other.left != 0 and one.left & other.right != 0:
                if other.left & one.right != 0 and other.right & one.right != 0:
                    return True
                
            if one.right & other.left != 0 and one.right & other.right != 0:
                if other.left & one.left != 0 and other.right & one.left != 0:
                    return True
                
            if other.left & one.left != 0 and other.left & one.right != 0:
                if one.left & other.right != 0 and one.right & other.right != 0:
                    return True
                
            if other.right & one.left != 0 and other.right & one.right != 0:
                if one.left & other.left != 0 and one.right & other.left != 0:
                    return True
        return False
        
        
        
    
    def doesRule2Apply(self, one, other):
        if one.left == other.left:
            if one.right & other.right != 0:
                return True
        if one.left == other.right:
            if one.right & other.left != 0:
                return True
        if one.right == other.left:
            if one.left & other.right != 0:
                return True
        if one.left == other.right:
            if one.right & other.left != 0:
                return True
#        print 'did not find Quartet in doesRule2Apply'
        return False    

    def rule2(self, one, other, bitkeys):
        if one.left == other.left:
            if one.right & other.right != 0:
                quartet = Quartet(one.left, one.right ^ other.right, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet
#                print 'invalid quartet in if 1 applyRule2 '
        if one.left == other.right:
            if one.right & other.left != 0:
                quartet = Quartet(one.left, one.right ^ other.left, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet
#                print 'invalid quartet in if 2 applyRule2 '
        if one.right == other.left:
            if one.left & other.right != 0:
                quartet = Quartet(one.right, one.left ^ other.right, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet
#                print 'invalid quartet in if 3 applyRule2 '
        if one.left == other.right:
            if one.right & other.left != 0:
                quartet = Quartet(one.left, one.right ^ other.left, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet
#                print 'invalid quartet in if 4 applyRule2 '
#        print 'Did not find quartet in applyRule2'

    def doesRule1Apply(self, one, other, bitkeys):
        if one.ruleOf3(other, bitkeys):
            if (one.left & other.left) ^ (one.left & other.right) == one.left:
                return True
            if (one.right & other.left) ^ (one.right & other.right) == one.right:
                return True
            if (other.left & one.left) ^ (other.left & one.right) == other.left:
                return True
            if (other.right & one.left) ^ (other.right & one.right) == other.right:
                return True
#        print 'Could not find a quartet in doesRule1Apply'
        return False
    
    def rule1a(self, one, other, bitkeys):
        if self.ruleOf3(one, other, bitkeys):
            if one.left & other.left == 0:
                quartet = Quartet(one.left, other.left, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet    
#                one.printQ(bitkeys)
#                other.printQ(bitkeys)
#                quartet.printQ(bitkeys)
#                print 'invalid quartet in if 1 rule1a '
            if one.left & other.right == 0:
                quartet = Quartet(one.left, other.right, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet    
#                one.printQ(bitkeys)
#                other.printQ(bitkeys)
#                quartet.printQ(bitkeys)
#                print 'invalid quartet in if 1 rule1a '
            if one.right & other.left == 0:
                quartet = Quartet(one.right, other.left, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet    
#                one.printQ(bitkeys)
#                other.printQ(bitkeys)
#                quartet.printQ(bitkeys)
#                print 'invalid quartet in if 1 rule1a '
            if one.right & other.right == 0:
                quartet = Quartet(one.right, other.right, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet    
#                one.printQ(bitkeys)
#                other.printQ(bitkeys)
#                quartet.printQ(bitkeys)
#                print 'invalid quartet in if 1 rule1a '
    
    
    def rule1b(self, one, other, bitkeys):
        if self.ruleOf3(one, other, bitkeys):
            if one.left & other.left != 0:
                nleft = (one.left | other.left) ^ (one.left & other.left)
                if nleft & one.right == 0:
                    quartet = Quartet(nleft, one.right, one, other)
                else:
                    quartet = Quartet(nleft, other.right, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet    
##                one.printQ(bitkeys)
#                other.printQ(bitkeys)
#                quartet.printQ(bitkeys)
#                print 'invalid quartet in if 1 rule1b '
            if one.left & other.right != 0:
                nleft = (one.left | other.right) ^ (one.left & other.right)
                if nleft & one.right == 0:
                    quartet = Quartet(nleft, one.right, one, other)
                else:
                    quartet = Quartet(nleft, other.left, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet    
#                one.printQ(bitkeys)
#                other.printQ(bitkeys)
#                quartet.printQ(bitkeys)
#                print 'invalid quartet in if 2 rule1b '
            if one.right & other.left != 0:
                nleft = (one.right | other.left) ^ (one.right & other.left)
                if nleft & one.left == 0:
                    quartet = Quartet(nleft, one.left)
                else:
                    quartet = Quartet(nleft, other.right)
                if quartet.isValidQuartet(bitkeys):
                    return quartet    
#                one.printQ(bitkeys)
#                other.printQ(bitkeys)
#                quartet.printQ(bitkeys)
#                print 'invalid quartet in if 3 rule1b '
            if one.right & other.right != 0:
                nleft = (one.right | other.right ) ^ (one.right & other.right)
                if nleft & one.left == 0:
                    quartet = Quartet(nleft, one.left, one, other)
                else:
                    quartet = Quartet(nleft, other.left, one, other)
                if quartet.isValidQuartet(bitkeys):
                    return quartet    
#                one.printQ(bitkeys)
#                other.printQ(bitkeys)
#                quartet.printQ(bitkeys)
#                print 'invalid quartet in if 4 rule1b '

    def ruleOf3(self, one, other, bitkeys):
        if one.left == other.left or one.left == other.right or one.right == other.left or one.right == other.right:
            return False
        if self._popcountTo3(one, other, bitkeys):
            return True
        return False
    
    def _popcountTo3(self, one, other, bitkeys):
        count = 0
        countie = (one.left | one.right) & (other.left | other.right)
        for bk in bitkeys:
            if bk & countie:
                count += 1
            if count == 3:
                return True
        return False    
    
    
#    The Quartet class is a container class for the two sides of a quartet, ie left and right
#    The class can perform some simple evaluations to determine whether it is infact a valid quartet
#    It can also prin itself. 
#    The bitkeys object passed to isValidQuartet and printQ should be a sorted list of the binary numbers 
#    that make up the quartets ie if there are 4 taxa a,b,c,d the list should be 1000, 0100, 0010 and 0001 
#    corresponding to the taxa. The trackkeeping is dealt with if you use the QuartetSet class which will
#    take care of all your indexing needs 
    
class Quartet(object):

    def __init__(self, _left = 0L, _right = 0L, _parent1=None, _parent2=None):
        self.parent1 = _parent1
        self.parent2 = _parent2
        if _left < _right:
            self.left = _left
            self.right = _right
        else:
            self.left = _right
            self.right = _left
    
    def __eq__(self, other):
        return isinstance(other, Quartet) and self.left == other.left and self.right == other.right

    def getParents(self):
        return self.parent1, self.parent2

    def isValidQuartet(self, bitkeys):
        if self.left & self.right != 0 or self._popcount(bitkeys) != 4:
            return False
        return True
    
    def _popcount(self, bitkeys):
        count = 0
        countie = self.left | self.right
        for bk in bitkeys:
            if countie < bk:
                return count
            if bk & countie:
                count += 1
        return count
            
    def printQ(self, bitkeys):
        intersection = ''
        for bk in bitkeys:
            if bk & self.left:
                intersection = intersection + '*'
            elif bk & self.right:
                intersection = intersection + '0'
            else:
                intersection = intersection + '.'
        print '%s' % (intersection)
#        print '%s, q.l: %s, q.r: %s' % (intersection, self.left, self.right)
