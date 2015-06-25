from Tree import Tree
from Node import Node
from func import read
from Var import var
from Glitch import Glitch

import sys, csv, operator, time

multiProcessing = False
try:
    
    from multiprocessing import Process as Process
    from multiprocessing import Queue, cpu_count
    cpu_count = cpu_count()
    multiProcessing = True
except ImportError:
    from threading import Thread as Process
    from Queue import Queue
    

def printVerticalNumbers(noOfNumbers):
    if noOfNumbers > 99999:
        print 'Ops, to many taxa in the tree, more than 99999 of them'
        return
    Tens = Hundreds = Thousands = Tenthousands = False
    base = ''
    if noOfNumbers >= 10:
        Tens = True
        tens = ''
        tensIndex = 0
    if noOfNumbers >= 100:
        Hundreds = True
        hundreds = ''
        hundredsIndex = 0
    if noOfNumbers >= 1000:
        Thousands = True
        thousands = ''
        thousandsIndex = 0
    if noOfNumbers >= 10000:
        Tenthousands = True
        tenthousands = ''
        tenthousandsIndex = 0
            
    for i in range(1,noOfNumbers+1):
        base += str(i % 10)
        if Tens:
            if i < 10:
                tens += ' '
            else:
                if i % 10 == 0:
                    tensIndex += 1
                    if tensIndex % 10 == 0:
                        tensIndex = 0
                tens += str(tensIndex)
        if Hundreds:
            if i < 100:
                hundreds += ' '
            else:
                if i % 100 == 0:
                    hundredsIndex += 1
                    if hundredsIndex % 10 == 0:
                        hundredsIndex = 0
                hundreds += str(hundredsIndex)
        if Thousands:
            if i < 1000:
                thousands += ' '
            else:
                if i % 1000 == 0:
                    thousandsIndex += 1
                    if thousandsIndex % 10 == 0:
                        thousandsIndex = 0
                thousands += str(thousandsIndex)
        if Tenthousands:
            if i < 10000:
                tenthousands += ' '
            else:
                if i % 10000 == 0:
                    tenthousandsIndex += 1
                    if tenthousandsIndex % 10 == 0:
                        tenthousandsIndex = 0
                tenthousands += str(tenthousandsIndex) 
    if Tenthousands:
        print tenthousands
    if Thousands:
        print thousands
    if Hundreds:
        print hundreds
    if Tens:
        print tens                        
    print base
        

#
#minorVersion = int(sys.version_info[1])
#
#if minorVersion < 6:
#    from threading import Thread as Process
#else:
#    from multiprocessing import Process as Process, Queue

class ConcurrentIntersections(Process):
    
    def __init__(self, bitKeys, minimumProportion, minNoOfTaxa, weights, listOfSplits, rooted, queue, verbose=1):
        Process.__init__(self)
        self.weights = weights
        self.minProp = minimumProportion
        self.minNoOfTaxa = minNoOfTaxa
        self.bitkeys = bitKeys
        self.bits = len(self.bitkeys)
        self.dict = {}
#        print 'Thread len(): ', len(listOfSplits)
        self.listOfSplits = listOfSplits
        self.rooted = rooted
        self.noTrees = float(len(listOfSplits))
        self.cutoff = sum(self.weights) * self.minProp
        self.weightedNoTrees = sum(self.weights)
        self.weight = 0
        self.queue = queue
        self.verbose = verbose
        
    def run(self):
        if self.verbose:
            sampler = int( float(len(self.listOfSplits))/20)
            if sampler < 1:
                sampler = 1
            sys.stdout.write('0% ')
            sys.stdout.flush()
        weightSoFar = 0.0
        for q in range(len(self.listOfSplits)):
            completeList = self.filteredList(self.weightedNoTrees - weightSoFar)
            weightSoFar += self.weights[q]
            dict = {}
            for s in self.listOfSplits[q]:
                for i in completeList:    
                    t = i[0] & s[0], (i[0] ^ s[0]) | i[1] | s[1]
                    dict[t] = i[2]
                    if not self.rooted:
                        t = (i[0] ^ s[0]) & i[0], (i[0] ^ (~s[0] &((self.bitkeys[-1] << 1) -1))) | i[0] | s[0]
                        dict[t] = i[2]
            
            for i,j in dict.items():
                self.addAndUpdate(long(i[0]),long(i[1]), self.weights[q], j)
            
            for s in self.listOfSplits[q]:
                self.add(s[0], s[1], self.weights[q]) 
            
            if self.verbose:
                if (q + 1) % sampler == 0:
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    
        if self.verbose:
            sys.stdout.write(' 100%\n')
            sys.stdout.flush()
        self.filteredList(self.weightedNoTrees - weightSoFar)
        self.weight = weightSoFar
        
        self.queue.put((self.dict, self.weight))
    
    def filteredList(self, remaining):
        list = []
        newDict = {}
        for key,value in self.dict.items():
            dict = {}
            for k,v in value.items():
                if remaining + v >= self.cutoff and self.enoughTaxa(k):
                    list.append([key,k,v])
                    dict[k] = v
            if len(dict.values()) > 0:
                newDict[key] = dict
        self.dict = newDict
        return list
    
    def add(self, index1, index2, weight):
        if not self._isTrivial(index1, index2):
            if index1 in self.dict:
                if index2 in self.dict[index1]:
#                    self.dict[index1][index2] = weight + self.dict[index1][index2]
                    return
                else:
                    self.dict[index1][index2] = weight
                    return
            else:
                self.dict[index1] = {}
                self.dict[index1][index2] = weight  
                return
    
    def addAndUpdate(self, index1, index2, weight1, weight2):
        if not self._isTrivial(index1, index2):
            if index1 in self.dict:
                if index2 in self.dict[index1]:
                    self.dict[index1][index2] = weight1 + self.dict[index1][index2]
                    return
                else:
                    self.dict[index1][index2] =  weight1 + weight2 
                    return
            else:
                self.dict[index1] = {}
                self.dict[index1][index2] =  weight1 + weight2
                return
    
    def _isTrivial(self, index1, index2):
        if index1 == 0:
            return True
        if self._popcount(index1) <= 1:
            return True
        if self._popcount(index1 | index2) == len(self.bitkeys):
            return True
        if index1 | index2 == 0:
            return True
        return False
    
    def enoughTaxa(self, countie):
#        print 'self.minNoOfTaxa: ', self.minNoOfTaxa
        count = 0
        for bk in self.bitkeys:
            if not bk & countie:
                count += 1
            if count == self.minNoOfTaxa:
#                print 'True'
                return True
#        print 'False'
        return False
    
    def _popcount(self, countie):
        count = 0
        for bk in self.bitkeys:
            if countie < bk:
                return count
            if bk & countie:
                count += 1
        return count   
    
    def getWeight(self):
        return self.weight
    
    def getDict(self):
        return self.dict
    
    def getResults(self):
        return self.results

class ConcurrentCombineIntersectionDicts(Process):
    
    def __init__(self, bitKeys, minimumProportion, dict1, dict2, rooted, weights, weight, queue):
        Process.__init__(self)
        self.weights = weights
        self.weight = weight
        self.minProp = minimumProportion
        self.bitkeys = bitKeys
        self.dict = {}
        self.dict1 = dict1
        self.dict2 = dict2
        self.rooted = rooted
        self.cutoff = round(sum(self.weights) * self.minProp)
        self.weightedNoTrees = sum(self.weights)
        self.queue = queue
        
    def run(self):
#        print 'Started thread: ', self.getName()
        list = []
        for key, value in self.dict2.items():
            for k,v in value.items():
                list.append([key,k,v])
                
        self.dict = self.dict1
         
#        list = self.combine()    
        sampler = float(len(list))/10
        weightSoFar = self.weight
        dict = {}
        for q in range(len(list)):
            completeList = self.filteredList(self.weightedNoTrees - weightSoFar)
            for i in completeList:    
                t = i[0] & list[q][0], (i[0] ^ list[q][0]) | i[1] | list[q][1]
                if dict.has_key(t):
                    if dict[t] < i[2]:
                        dict[t] = i[2]
                else:
                    dict[t] = i[2]
                if not self.rooted:
                    t = (i[0] ^ list[q][0]) & i[0], (i[0] ^ (~list[q][0] &((self.bitkeys[-1] << 1) -1))) | i[0] | list[q][0]
                    if dict.has_key(t):
                        if dict[t] < i[2]:
                            dict[t] =  i[2]
                    else:
                        dict[t] =  i[2]
            if q + 1 % sampler == 0:
                print 'Thread %s processed %s combinations' % (self.getName(), q +1)
                sys.stdout.flush()
                
        for i,j in dict.iteritems():
            self.addAndUpdate(long(i[0]),long(i[1]), j, list[q][2])
        
        for q in range(len(list)):
            self.add(list[q][0], list[q][1], list[q][2])
    
            
                
        self.printDict(self.dict)
        self.queue.put(self.dict)
            #======================================================================
            # completeList = self.filteredList(self.weightedNoTrees - weightSoFar)
            #======================================================================
        
    def add(self, index1, index2, weight):
        if not self._isTrivial(index1, index2):
            if index1 in self.dict:
                if index2 in self.dict[index1]:
#                    self.dict[index1][index2] = weight + self.dict[index1][index2] 
                    return
                else:
                    self.dict[index1][index2] = weight
                    return
            else:
                self.dict[index1] = {}
                self.dict[index1][index2] = weight  
                return
    
    def addAndUpdate(self, index1, index2, weight1, weight2):
        if not self._isTrivial(index1, index2):
            if index1 in self.dict:
                if index2 in self.dict[index1]:
                    self.dict[index1][index2] = weight2 + self.dict[index1][index2]
                    return
                else:
                    self.dict[index1][index2] = weight1 + weight2
                    return
            else:
                self.dict[index1] = {}
                self.dict[index1][index2] = weight1 + weight2
                return
    
    def _isTrivial(self, index1, index2):
        if index1 == 0:
            return True
        if self._popcount(index1) <= 1:
            return True
        if self._popcount(index1 | index2) == len(self.bitkeys):
            return True
        if index1 | index2 == 0:
            return True
        return False
    
    def _popcount(self, countie):
        count = 0
        for bk in self.bitkeys:
            if countie < bk:
                return count
            if bk & countie:
                count += 1
        return count   
       
    def filteredList(self, remaining):
        list = []
        newDict = {}
        for key,value in self.dict.items():
            dict = {}
            for k,v in value.items():
#                if remaining + v >= self.cutoff:
                list.append([key,k, v])
                dict[k] = v
#                else:
#                    print 'Cut: ', self.cutoff
#                    print 'Rem + v: ', remaining +v
#                    print 'DOOOOOOOGH'
            if len(dict.values()) > 0:
                newDict[key] = dict
        self.dict = newDict
        return list
    
    def printDict(self, dict):
        
        for k, v in dict.items():
            for key, value in v.items():
                print '%s, %s, %s ' % (k, key, value)
                
    
    def combine(self):
        
        self.dict = self.dict1
        
        list = []
        for k, v in self.dict2.items():
            for k2, v2 in v.items():
                list.append([k,k2,v2])
                
        return list
        
#        list = []
#        for k1, v1 in self.dict1.items():
#            if k1 in self.dict2:
#                dict = {}
#                for k2, v2 in v1.items():
#                    if k2 in self.dict2[k1]:
#                        dict[k2] = self.dict2[k1][k2] + v2
#                    else:
#                        list.append([k1,k2,v2])
#                if len(dict.keys()) > 0:
#                    self.dict[k1] = dict
#            else:
#                for k2, v2 in v1.items():
#                    list.append([k1,k2,v2])
#        print list
#        for k1, v1 in self.dict2.items():
#            if k1 in self.dict1:
#                for k2, v2 in v1.items():
#                    if not k2 in self.dict1[k1]:
#                        list.append([k1,k2,v2])
#            else:
#                for k2, v2 in v1.items():
#                    list.append([k1,k2,v2])
#        print list
#        return list
        
    def getCombinedDict(self):
        return self.dict

class SimpleIntersections(object):
    
    def __init__(self, bitKeys, minimumProportion, minNoOfTaxa, weights, listOfSplits, rooted):
        self.weights = weights
        self.minProp = minimumProportion
        self.minNoOfTaxa = minNoOfTaxa
        self.bitkeys = bitKeys
        self.bits = len(self.bitkeys)
        self.dict = {}
        self.inputSplits = listOfSplits
        self.splits = []
        self.rooted = rooted
        self.cutoff = sum(self.weights) * self.minProp
        self.weightedNoTrees = sum(self.weights)
        self.allOnes = 2L**(len(self.bitkeys)) - 1
        
        
    def createIntersections(self):
        sampler = int( float(len(self.listOfSplits))/20)
        if sampler < 1:
            sampler = 1
        weightSoFar = 0.0
        sys.stdout.write('0% ')
        sys.stdout.flush()
        
        for inputList in self.inputSplits:
            
            for split in inputList:
                self.splits.append()
        
        for q in range(len(self.listOfSplits)):
            completeList = self.filteredList(self.weightedNoTrees - weightSoFar)
            weightSoFar += self.weights[q]
            dict = {}
            for s in self.listOfSplits[q]:
                for i in completeList:    
                    t = i[0] & s[0], (i[0] ^ s[0]) | i[1] | s[1]
                    dict[t] = i[2]
                    if not self.rooted:
                        t = (i[0] ^ s[0]) & i[0], (i[0] ^ (~s[0] &((self.bitkeys[-1] << 1) -1))) | i[0] | s[0]
                        dict[t] = i[2]
            
            for i,j in dict.items():
                self.addAndUpdate(long(i[0]),long(i[1]), self.weights[q], j)
            
            for s in self.listOfSplits[q]:
                self.add(s[0], s[1], self.weights[q]) 
            
            if (q + 1) % sampler == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
        
        sys.stdout.write(' 100%\n')
        sys.stdout.flush()
        self.filteredList(self.weightedNoTrees - weightSoFar)
        self.weight = weightSoFar
        
        self.queue.put((self.dict, self.weight))
        

class BuildIntersections(object):
    
    def __init__(self, bitKeys, minimumProportion, minNoOfTaxa, weights, listOfSplits, rooted, sortByNoSplits, verbose=1):
        self.verbose=verbose
        self.sortByNoSplits = sortByNoSplits
        self.weights = weights
        self.minProp = minimumProportion
        self.minNoOfTaxa = minNoOfTaxa
        self.bitkeys = bitKeys
        self.bits = len(self.bitkeys)
        self.dict = {}
        self.listOfSplits = listOfSplits
        self.rooted = rooted
        self.cutoff = sum(self.weights) * self.minProp
        self.weightedNoTrees = sum(self.weights)
        self.allOnes = 2L**(len(self.bitkeys)) - 1
#        print 'self.weightedNoTrees:', self.weightedNoTrees
        
    def buildIntersections(self, verbose=1):
        
        if multiProcessing and False:
            print 'Starting dual cpu computations'
            queue = Queue()
#            print 'CPUS: ', cpu_count
            threads = []
#            cpus = cpu_count
#            divider = len(self.listOfSplits)/cpus*4
        
#            print 'Len(): ', len(self.listOfSplits)
        
#            print 'Divider: ', divider
        
#            for i in range(cpus):
#                start = i*divider
#                stop = (i+1)*divider
#                print 'List of splits %s:%s ' % (start, stop)
#                print 'Len(): ', len(self.listOfSplits[start:stop])
#                thread = ConcurrentIntersections(self.bitkeys, self.minProp, self.minNoOfTaxa, self.weights[start:stop], self.listOfSplits[start:stop], self.rooted, queue)
#                thread.start()
#                threads.append(thread)

#        results = []
#        start = 0
#        stop = 25
#        thread = ConcurrentIntersections(self.bitkeys, self.minProp, self.weights[start:stop], self.listOfSplits[start:stop], self.rooted, queue)
#        thread.start()
#        threads.append(thread)
#        results.append(queue.get())
#        thread.join()
#        print 'List of splits %s:%s ' % (start,stop)
#        print 'Len(): ', len(self.listOfSplits[start:stop])
#        print self.listOfSplits[start:stop]
#        
#        start = 25
#        stop = 50
#        thread = ConcurrentIntersections(self.bitkeys, self.minProp, self.weights[start:stop], self.listOfSplits[start:stop], self.rooted, queue)
#        thread.start()
#        threads.append(thread)
#        results.append(queue.get())
#        thread.join()
#        print 'List of splits %s:%s ' % (start, stop)
#        print 'Len(): ', len(self.listOfSplits[start:stop])
#        print self.listOfSplits[start:stop]
#        

#            results = []
#            for i in range(cpus):
#                results.append(queue.get())
        
#        for t in threads:
#            t.join()
            thread1 = ConcurrentIntersections(self.bitkeys, self.minProp, self.minNoOfTaxa, self.weights, self.listOfSplits[len(self.listOfSplits)/2:], self.rooted, queue)
            thread2 = ConcurrentIntersections(self.bitkeys, self.minProp, self.minNoOfTaxa, self.weights, self.listOfSplits[:len(self.listOfSplits)/2], self.rooted, queue)
        
#            thread1.setName('Worker 1')
#            thread2.setName('Worker 2')

            results = []
            thread1.start()
            results.append(queue.get())
            thread1.join()
            dict1 = results[0][0]
            print dict1
            weight1 = results[0][1]     
            print 'weight1: ', weight1
            self.dict = dict1 
            
            
            thread2.start()
        
            results.append(queue.get())
            
            thread2.join()
            
            dict2 = results[0][0]
            print dict2
            weight2 = results[1][1]
            print 'weight2: ', weight2
            self.dict = dict2
#        
            combine = ConcurrentCombineIntersectionDicts(self.bitkeys, self.minProp, dict1, dict2, self.rooted, self.weights, weight1+weight2, queue)
#        
#            combine.setName('Combinator')
#
            combine.start()
#        
            combine.join()
#        
#            self.dict = combine.getCombinedDict()
            
            self.dict = queue.get()
#            self.dict = dict2
            
            print 'self.dict: ', self.dict        
        
        else:
#            print 'Starting single cpu computations'
            queue = Queue()
            threads = []
            thread = ConcurrentIntersections(self.bitkeys, self.minProp, self.minNoOfTaxa, self.weights, self.listOfSplits, self.rooted, queue, verbose=0)
            thread.start()
            threads.append(thread)
                   
            results = []
            results.append(queue.get())
        
            dict1 = results[0][0]
            weight1 = results[0][1]  
            self.dict = dict1   
        
        self.finalList = self.informativeList()
        
        return self.finalList

    def buildExtendedProfile(self, generations=2):
        
        alreadyAdded = {}
        newlyAdded = {}
        completeList = []
        for list in self.finalList:
            for i in list:
                if i[2] == 1.0 and self._addToDict(alreadyAdded, i[0], i[1]):
                    completeList.append(i)

        for counter in range(generations):
            list = []
            for i in range(0,len(completeList)-1):
                for j in range(i+1, len(completeList)):
                    extend = self._extend(completeList[i], completeList[j])
                    if len(extend) > 0:
                        if self._addToDict(alreadyAdded, extend[0][0], extend[0][1]):
                            self._addToDict(newlyAdded, extend[0][0], extend[0][1])
                            list.append(extend[0])
                        if self._addToDict(alreadyAdded, extend[1][0], extend[1][1]):
                            self._addToDict(newlyAdded, extend[1][0], extend[1][1])
                            list.append(extend[1])
            completeList.extend(list)
        
        list = []
        for key in newlyAdded.keys():
            list.extend(self._getInformativeList(newlyAdded[key],key))
     
        list = self._sortByLeafSet(list)

        tempList = []
        for l in list:
            tempList.extend(self._buildUnconflictingLists(l))
        
        return tempList
    
    def _addToDict(self, dict, i,j):
        if i in dict:
            if j not in dict[i]:
                dict[i][j] = self.weightedNoTrees
                return True
            else:
                return False
        else:
            dict[i] = {}
            dict[i][j] = self.weightedNoTrees
            
        return True
    
    def _extend(self,i,j):
 
        if i[0] & j[0] == 0:
            if i[0] & j[1] == 0:
                if i[1] & j[0] == 0:
                    return [[i[0],i[1] | j[1], self.weightedNoTrees],[j[0],i[1] | j[1], self.weightedNoTrees]]
        
        if i[0] | j[0] == i[0] or i[0] | j[0] == j[0]:
            if i[0] & j[1] == 0 and j[0] & i[1] == 0:
                return [[i[0],i[1] | j[1], self.weightedNoTrees],[j[0],i[1] | j[1], self.weightedNoTrees]]
        
        return []
    
    def popcountLimit(self, n, limit):
        count = 0
        for bk in self.bitkeys:
            if bk & n:
                count += 1
            if n >= limit:
                return True
        return False
    
    
    
    def isCompatible(self, i, j, excluded):
        
        def isPowerOf2(number): return (not(number&(number-1)))
        
        if not self.popcountLimit(i & j, 2) or not self.popcountLimit(i ^ j, 2):
            return True 
        else:
            ones = i & j
            one1 = 0
            if ones > 0:
                if not isPowerOf2(ones):
                    one1 = 2
                else:
                    one1 = 1
#                        print 'one1: ', one1
                        
            zeroOnes = (i ^ j) & j
            zero1 = 0
            if zeroOnes > 0:
                if not isPowerOf2(zeroOnes):
                    zero1 = 2
                else:
                    zero1 = 1
#                        print 'zero1: ', zero1
                        
            oneZeroes = (i ^ j) & i
            one0 = 0
            if oneZeroes > 0:
                if not isPowerOf2(oneZeroes):
                    one0 = 2
                else:
                    one0 = 1
#                        print 'one0: ', one0
  
            zeroZeroes = ((j ^ self.allOnes) & (i ^ self.allOnes)) ^ excluded
            zero0 = 0 
            if zeroZeroes > 0:
                if not isPowerOf2(zeroZeroes):
                    zero0 = 2
                else:
                    zero0 = 1
#                        print 'zero0: ', zero0
                                 
                        
            if one1 and zero1 and one0 and zero0:
                return False

        return True
    
    def quickList(self, list, bitkeys):
        
#        print 'Initial list'
#        for i in list:
#            i.printIntersection(bitkeys)
        
#        print 
        
        for i in range(0,len(list)-1):
            for j in range(i+1, len(list)):
                if list[i].frequency < list[j].frequency:
                    temp = list[i]
                    list[i] = list[j]
                    list[j] = temp
                elif list[i].frequency == list[j].frequency:
                    pop1 = self._pop(list[i].right, bitkeys)
                    pop2 = self._pop(list[j].right, bitkeys)
                    if pop1 > pop2:
                        temp = list[i]
                        list[i] = list[j]
                        list[j] = temp
                    elif pop1 == pop2:
                        if list[i].firstHit(bitkeys) > list[j].firstHit(bitkeys):
                            temp = list[i]
                            list[i] = list[j]
                            list[j] = temp
                                
                        elif list[i].firstHit(bitkeys) == list[j].firstHit(bitkeys):
                            if self._pop(list[i].left, bitkeys) < self._pop(list[j].left, bitkeys):
                                temp = list[i]
                                list[i] = list[j]
                                list[j] = temp

#        print 'Sorted'
#        for i in list:
#            i.printIntersection(bitkeys)        
        
        nonRedundantList = []
        
        for i in range(0,len(list)-1):
            for j in range(i+1, len(list)):
                if list[i].left | list[j].right == list[j].left | list[j].right:
                    if list[i].right | list[j].right == list[j].right:
                        list[j].hits = 0
        
        if self.verbose:
            print 'Splits'
        for intersection in list:
            if intersection.hits > 0:
                nonRedundantList.append(intersection)
                if self.verbose:
                    intersection.printIntersection(bitkeys, self.weightedNoTrees) 
        
        list = nonRedundantList
        
        self.fullList = []
        for i in nonRedundantList:
            copy = i.copy()
            copy.frequency = copy.frequency/self.weightedNoTrees
            self.fullList.append(copy)
            
        lists = []
        
        while len(list) > 0 and list[0].frequency >= self.cutoff:
            split = list.pop(0)
            exclude = split.right
            compatibleList = [split]
            for i in range(len(list)):
                added = False
                wannabe = list.pop(0)
                if wannabe.right | exclude == exclude: 
                    temp = wannabe.copy()
                    if wannabe.right != exclude:
                        added = True
                        list.append(temp)
                        wannabe.right = wannabe.right | exclude
                        wannabe.left = wannabe.left ^ (wannabe.left & exclude)
                    compatible = True
                    for j in range(len(compatibleList)):
                        if compatibleList[j].left == wannabe.left:
                            compatible = False
                            break
                        if not self.isCompatible(compatibleList[j].left, wannabe.left, exclude):
                            compatible = False
                            break
                    if compatible:
                        if self._pop(wannabe.left, bitkeys) > 1:
                            compatibleList.append(wannabe)
                        else:
                            if not added:
                                list.append(temp)
                    else:
                        if not added:
                            list.append(temp)
                    
                else:
                    list.append(wannabe)
            if len(compatibleList) > 0:
                lists.append(compatibleList)
        
#        print 'Lists'
#        for list in lists:
#            for i in list:
#                i.printIntersection(bitkeys)
#            print

#        print 'len(lists): ', len(lists) 
        
        for index in range(1,len(lists)):
            excluded = lists[index][0].right
            for i in lists[0]:
                j = i.copy()
                j.right = excluded
                j.left = j.left ^ (j.left & excluded)
                if j.left == 0 or j.popcount(bitkeys) <= 1:
                    pass
                else:
                    append = True
                    for k in lists[index]:
                        if k.left == j.left:
                            append = False
                            if k.frequency < j.frequency:
                                k.frequency = j.frequency
                            break
                        
                    if append:
                        lists[index].append(j)
        
        
        if self.sortByNoSplits:
            for i in range(0,len(lists)-1):
                for j in range(i, len(lists)):
                    if len(lists[i]) > len(lists[j]):
                        temp = lists[i]
                        lists[i] = lists[j]
                        lists[j] = temp
                    elif len(lists[i]) == len(lists[j]):
                        if lists[i][0].popcountExcluded(bitkeys) > lists[j][0].popcountExcluded(bitkeys):
                            temp = lists[i]
                            lists[i] = lists[j]
                            lists[j] = temp
            
        
            
#        print 'Partitioned splits'
        for list in lists:
            for i in list:
                i.setExcluded(bitkeys)
                i.frequency = i.frequency/self.weightedNoTrees
#                i.printIntersection(bitkeys) 
                
#            print ''
    
        return lists
        
    
    def informativeList(self):
        
        self._remove()
        
        list = []
        for key, value in self.dict.items():
            for k,v in value.items():
                if v >= self.cutoff:
                    list.append(Intersection(key,k,v))
        
#        for i in list:
#            i.printIntersection(self.bitkeys)
#        print ''

        for i in range(0,len(list)-1):
            for j in range(i+1, len(list)):
                if list[i].frequency < list[j].frequency:
                    temp = list[i]
                    list[i] = list[j]
                    list[j] = temp
                elif list[i].frequency == list[j].frequency:
                    if self._popcount(list[i].right) > self._popcount(list[j].right):
                        temp = list[i]
                        list[i] = list[j]
                        list[j] = temp
                    elif self._popcount(list[i].right) == self._popcount(list[j].right):
                        if list[i].firstHit(self.bitkeys) > list[j].firstHit(self.bitkeys):
                            temp = list[i]
                            list[i] = list[j]
                            list[j] = temp
                                
                        elif list[i].firstHit(self.bitkeys) == list[j].firstHit(self.bitkeys):
                            if self._popcount(list[i].left) < self._popcount(list[j].left):
                                temp = list[i]
                                list[i] = list[j]
                                list[j] = temp
                            

#        print 'Sorted'
#        for i in list:
#            i.printIntersection(self.bitkeys)        
        
        nonRedundantList = []
        
        for i in range(0,len(list)-1):
            for j in range(i+1, len(list)):
                if list[i].left | list[j].right == list[j].left | list[j].right:
                    if list[i].right | list[j].right == list[j].right:
                        list[j].hits = 0
        
        if self.verbose:
            print 'Splits'
        for intersection in list:
            if intersection.hits > 0:
                nonRedundantList.append(intersection)
                if self.verbose:
                    intersection.printIntersection(self.bitkeys, self.weightedNoTrees) 
        
        list = nonRedundantList
        
        self.fullList = []
        for i in nonRedundantList:
            copy = i.copy()
            copy.frequency = copy.frequency/self.weightedNoTrees
            self.fullList.append(copy)
            
        lists = []
        
        while len(list) > 0 and list[0].frequency >= self.cutoff:
            split = list.pop(0)
            exclude = split.right
            compatibleList = [split]
            for i in range(len(list)):
                added = False
                wannabe = list.pop(0)
                if wannabe.right | exclude == exclude: 
                    temp = wannabe.copy()
                    if wannabe.right != exclude:
                        added = True
                        list.append(temp)
                        wannabe.right = wannabe.right | exclude
                        wannabe.left = wannabe.left ^ (wannabe.left & exclude)
                    compatible = True
                    for j in range(len(compatibleList)):
                        if compatibleList[j].left == wannabe.left:
                            compatible = False
                            break
                        if not self.isCompatible(compatibleList[j].left, wannabe.left, exclude):
                            compatible = False
                            break
                    if compatible:
                        if self._popcount(wannabe.left) > 1:
                            compatibleList.append(wannabe)
                        else:
                            if not added:
                                list.append(temp)
                    else:
                        if not added:
                            list.append(temp)
                    
                else:
                    list.append(wannabe)
            if len(compatibleList) > 0:
                lists.append(compatibleList)
        
#        print 'len(lists): ', len(lists) 
        
        for index in range(1,len(lists)):
            excluded = lists[index][0].right
            for i in lists[0]:
                j = i.copy()
                j.right = excluded
                j.left = j.left ^ (j.left & excluded)
                if j.left == 0 or j.popcount(self.bitkeys) <= 1:
                    pass
                else:
                    append = True
                    for k in lists[index]:
                        if k.left == j.left:
                            append = False
                            if k.frequency < j.frequency:
                                k.frequency = j.frequency
                            break
                        
                    if append:
                        lists[index].append(j)
        
        
        if self.sortByNoSplits:
            for i in range(0,len(lists)-1):
                for j in range(i, len(lists)):
                    if len(lists[i]) > len(lists[j]):
                        temp = lists[i]
                        lists[i] = lists[j]
                        lists[j] = temp
                    elif len(lists[i]) == len(lists[j]):
                        if lists[i][0].popcountExcluded(self.bitkeys) > lists[j][0].popcountExcluded(self.bitkeys):
                            temp = lists[i]
                            lists[i] = lists[j]
                            lists[j] = temp
            
        
            
#        print 'Partitioned splits'
        for list in lists:
            for i in list:
                i.setExcluded(self.bitkeys)
                i.frequency = i.frequency/self.weightedNoTrees
#                i.printIntersection(self.bitkeys) 
                
#            print ''
    
        return lists
        
    def _informativeList(self):
        
        self._remove()
        
        list = []
        for key, value in self.dict.items():
            for k,v in value.items():
                list.append(Intersection(key,k,v))
        
        for i in range(0,len(list)-1):
                for j in range(i+1, len(list)):
                    if list[i].frequency < list[j].frequency:
                        temp = list[i]
                        list[i] = list[j]
                        list[j] = temp
                    elif list[i].frequency == list[j].frequency and self._popcount(list[i].right) < self._popcount(list[j].right):
                        temp = list[i]
                        list[i] = list[j]
                        list[j] = temp
        
#        for i in list:
#            i.printIntersection(self.bitkeys)
        
#        print 'Consolidate within: '
        list = self._consolidateWithinFrequency(list)
#        for i in list:
#            i.printIntersection(self.bitkeys)
        
        list = self._consolidateBetweenFrequencies(list)
#        print 'Consolidate between:'
#        for i in list:
#            i.printIntersection(self.bitkeys)
            
        for i in range(0, len(list)-1):
            for j in range(i+1, len(list)):
                if list[i].frequency == list[j].frequency:
                    if list[i].popcountExcluded(self.bitkeys) > list[j].popcountExcluded(self.bitkeys):
                        temp = list[i]
                        list[i] = list[j]
                        list[j] = temp
                    elif list[i].firstHit(self.bitkeys) > list[j].firstHit(self.bitkeys):
                        temp = list[i]
                        list[i] = list[j]
                        list[j] = temp
                    elif list[i].firstHit(self.bitkeys) == list[j].firstHit(self.bitkeys):
                        if list[i].popcount(self.bitkeys) > list[j].popcount(self.bitkeys):
                            temp = list[i]
                            list[i] = list[j]
                            list[j] = temp
                        
        
        print 'Splits'
        printVerticalNumbers(len(self.bitkeys))
        for i in list:
            i.printIntersection(self.bitkeys)
                    
        smallestList, remainder = self._smallestList(list)
        
#        for i in remainder:
#            i.printIntersection(self.bitkeys)
        
        lists = []
        lists.append(smallestList)
        oldList = smallestList
        added = {}
        for intersection in remainder:
            if not added.has_key(intersection.right):
                temp = []
                addable = self._addable2List(intersection, oldList)
                if addable == 2:
                    temp = self._add2List(intersection, oldList)
                    for i in remainder:
                        addableWithRestraint = self._addable2ListWithRestraint(i, temp)
                        if addableWithRestraint == 2:
                            added[i.right] = 1
                            self._add2ListWithRestraint(i, temp)
#                            temp.append(i)
                        elif addableWithRestraint == 1:
                            added[i.right] = 1
                            self._replaceWithRestraint(i, temp)
                
                elif addable == 1:
                    temp = self._replace(intersection, oldList)
                    for i in remainder: 
                        addableWithRestraint = self._addable2ListWithRestraint(i, temp)
                        if addableWithRestraint == 2:
                            self._add2ListWithRestraint(i, temp)
#                            temp.append(i)
                            added[i.right] = 1
                        elif addableWithRestraint == 1:
                            self._replaceWithRestraint(i, temp)
                            added[i.right] = 1
                        
                if len(temp) > 0:
                    lists.append(temp)
                    
        
        for i in range(0, len(lists)-1):
            for j in range(i+1, len(lists)):
                if lists[i][0].popcountExcluded(self.bitkeys) > lists[j][0].popcountExcluded(self.bitkeys):
                    temp = lists[i]
                    lists[i] = lists[j]
                    lists[j] = temp
        
                    
        for list in lists:
            for i in range(0,len(list)-1):
                for j in range(i+1, len(list)):
                    
                    if list[i].popcount(self.bitkeys) > list[j].popcount(self.bitkeys):
                        temp = list[i]
                        list[i] = list[j]
                        list[j] = temp
                    elif list[i].popcount(self.bitkeys) == list[j].popcount(self.bitkeys):
                        if list[i].firstHit(self.bitkeys) > list[j].firstHit(self.bitkeys):
                            temp = list[i]
                            list[i] = list[j]
                            list[j] = temp
                            
            for i in list:
                i.setExcluded(self.bitkeys)
                i.frequency = i.frequency/self.weightedNoTrees
        
        for list in lists:
            for i in list:
                i.printIntersection(self.bitkeys)
            print ''
        
        return lists
        
    def _replaceWithRestraint(self, intersection, intersections):
        counter = 0
        for i in intersections:
            if intersection.left == i.left:
               # print 'replaced %s with %s' % (i.frequency, intersection.frequency)
                i.frequency = intersection.frequency
                counter += 1
                
        if counter > 1:
            intersection.printIntersection(self.bitkeys)
                
    def _replace(self, intersection, intersections):
        copy = intersection.copy()
        jointExcluded = intersection.right | intersections[0].right
        copy.left = (copy.left & jointExcluded) ^ copy.left
        copy.right = jointExcluded
        first = True
        dict = {}
        for i in intersections:
            c = i.copy()
            c.left = (c.left & jointExcluded) ^ c.left
            if self._popcount(c.left) > 1:
                c.right = jointExcluded
                if copy.left == c.left:
                    if copy.frequency > c.frequency:
                        c.frequency = copy.frequency
#                    if first:
#                        first = False
                    if dict.has_key(c.left):
                        if c.frequency > dict[c.left].frequency:
                            dict[c.left] = c
                    else:
                        dict[c.left] = c
                else:
                    if dict.has_key(c.left):
                        if c.frequency > dict[c.left].frequency:
                            dict[c.left] = c
                    else:
                        dict[c.left] = c
                
        return dict.values()
    
    def _add2ListWithRestraint(self, intersection, intersections):
        for i in intersections:
            if intersection.left == i.left:
                print 'Adding identical intersection to list'
        intersections.append(intersection)
    
    def _add2List(self, intersection, intersections):
        jointExcluded = intersection.right | intersections[0].right
        copy = intersection.copy()
        copy.left = (copy.left & jointExcluded) ^ copy.left
        copy.right = jointExcluded
        add = True
        dict = {}
        for i in intersections:
            temp = i.copy()
            temp.left = (i.left & jointExcluded) ^ i.left
            temp.right = jointExcluded
            if self._popcount(temp.left) > 1:
                if dict.has_key(temp.left):
                    if dict[temp.left].frequency < temp.frequency:
                        if temp.left == copy.left and temp.frequency >= copy.frequency:
                            add = False
#                    list.append(temp)
                        dict[temp.left] = temp
                else:
                    if temp.left == copy.left and temp.frequency >= copy.frequency:
                        add = False
#                    list.append(temp)
                    dict[temp.left] = temp
        list = dict.values()
        if add:
            list.append(copy)
        return list
    
    def _addable2ListWithRestraint(self, proposed, intersections):
        if not proposed.right == intersections[0].right:
            return 0
        higher = False
        for intersection in intersections:
#            temp = proposed.right & intersection.left
            
#            if proposed.left == intersection.left ^ temp:
            if proposed.left == intersection.left:
                if proposed.frequency <= intersection.frequency:
                    return 0
                else:
                    higher = True
        if higher:
            return 1
        return 2
    
    def _addable2List(self, proposed, intersections):
        higher = False
        for intersection in intersections:
            temp = proposed.right & intersection.left
            if proposed.left == intersection.left ^ temp:
                if proposed.frequency <= intersection.frequency:
                    return 0
                else:
                    higher = True
        if higher:
            return 1
        return 2
            
    
    def _consolidateBetweenFrequencies(self, intersections):
        list = []
        
        for index in range(len(intersections)-1,-1, -1):
            if index == 0:
                list.append(intersections[index])
                list.reverse()
                return list
#            intersections[index].printIntersection(self.bitkeys)
            pop = self._popcount(intersections[index].right)
            left = intersections[index].left
            frequency = intersections[index].frequency
            replace = False
            for jndex in range(index-1, -1, -1):
#                intersections[jndex].printIntersection(self.bitkeys)
                if intersections[jndex].left == left:
#                    print 'true' 
                    if frequency < intersections[jndex].frequency:
#                        print 'True'
                        if pop > self._popcount(intersections[jndex].right):
#                            print 'TRUE'
                            replace = True
#            print ''
            if not replace:
                list.append(intersections[index])
                
        return list
    
    def _consolidateWithinFrequency(self,intersections):
        list = []
        for index in range(len(intersections)):
            if index == len(intersections)-1:
                list.append(intersections[index])
                return list    
            pop = self._popcount(intersections[index].right)
            if pop > 0:
                replace = False
                makeup = intersections[index].left | intersections[index].right
                hits = intersections[index].frequency
#                print '&', intersections[index].left | intersections[index].right
                for jndex in range(index+1, len(intersections)):
#                    print '.', intersections[jndex].left
                    if intersections[jndex].left == makeup:
#                        print 'true'
                        if intersections[jndex].frequency == hits:
                            replace = True
                if not replace:
                    list.append(intersections[index])
            else:
                list.append(intersections[index])
                
        return list
        
    def _smallestList(self, list):
        smallest = sys.maxint
        smallestList = []
        remainder = []
        for intersection in list:
            pop = self._popcount(intersection.right)
            if pop < smallest:
                smallest = pop
                remainder.extend(smallestList)
                smallestList = [intersection]
            elif pop == smallest:
                smallestList.append(intersection)
            else:
                remainder.append(intersection)
                
        return smallestList, remainder
    
    def _remove(self):
        
        dict = {}
        for key, value in self.dict.items():
            temp = {}        
            for v in value.values():
                temp[v] = 1
            tempDict = {}
            for k in temp.keys():
                list = []
                smallest = sys.maxint
                for ke, v in value.items():
                    if v == k:
                        list.append((ke,v))
                        pop = self._popcount(ke)
                        if pop < smallest:
                            smallest = pop
                
                for tuple in list:
                    pop = self._popcount(tuple[0])
                    if pop == smallest:
                        tempDict[tuple[0]] = tuple[1]
            
            dict[key] = tempDict
        
        self.dict = dict
                
    def removeRedundant(self):
        
        leafsets = {}
        
        for d1 in self.dict.keys():
            for k,v in self.dict[d1].items():
                leafset = d1 | k
                if leafsets.has_key(leafset):
                    if leafsets[leafset].has_key(v):
                        leafsets[leafset][v].append((d1,k))
                    else:
                        leafsets[leafset][v] = [(d1,k)]
                else:
                    leafsets[leafset] = {}
                    leafsets[leafset][v] = [(d1,k)]
 
        dict = {}
        for k,v in leafsets.items():
            for hits, list in v.items():
                biggestList = []
                biggest = 0
                for i in range(len(list)):
                    pop = self._popcount(list[i][0])
                    if pop > biggest:
                         biggestList = [list[i]]
                         biggest = pop
                    elif pop == biggest:
                        biggestList.append(list[i])
            
                for i in biggestList:
                    if i[0] in dict:
                        if i[1] not in dict[i[0]]:
                            dict[i[0]][i[1]] = hits
                    else:
                        dict[i[0]] = {}
                        dict[i[0]][i[1]] = hits   
        self.dict = dict
    
    def _getInformativeList(self, dict, left):
        
        freqDict = {}
        for v in dict.keys():
            if dict[v] >= self.cutoff:
                if freqDict.has_key(dict[v]):
                    freqDict[dict[v]].append(v)
                else:
                    freqDict[dict[v]] = [v]
        
        informativeList = []
        for v in freqDict.keys():
            for i in self._getMostInformative(freqDict[v]):
#                print 'v:', v
#                print 'self.weightedNoTrees: ', self.weightedNoTrees
                informativeList.append([left, i,v/self.weightedNoTrees])
            
        return informativeList

    def _getMostInformative(self, list):
        smallest = self._popcount(list[0])
        index = [list[0]]
        for i in range(1,len(list)):
            pop = self._popcount(list[i])
            if pop == smallest:
                index.append(list[i])
            elif pop < smallest:
                smallest = pop
                index = [list[i]]
        return index
    
    def _pop(self, countie, bitkeys):
        count = 0
        for bk in bitkeys:
            if countie < bk:
                return count
            if bk & countie:
                count += 1
        return count 
    
    def _popcount(self, countie):
        count = 0
        for bk in self.bitkeys:
            if countie < bk:
                return count
            if bk & countie:
                count += 1
        return count   
    
    def _sortByLeafSet(self, list):
        leafDict = {}
        for i in list:
            if leafDict.has_key(i[1]):
                leafDict[i[1]].append(i)
            else:
                leafDict[i[1]] = [i]
        list = []
        for v in leafDict.values():
            list.append(v)
        return list

    def _buildUnconflictingLists(self, list, sort=True):
        
        if sort:
            self._sortListByHits(list)
        
        unconflictingLists, unresolved = self._unConflictingList(list)

        if len(unresolved) > 0:
            lenght = len(unresolved) +1
            while len(unresolved) < lenght:
                unconflictingLists, unresolved = self._unconflictingList(unconflictingLists, unresolved, False)
                lenght = len(unresolved)
                
            if lenght > 0:
                unconflictingLists, unresolved = self._unconflictingList(unconflictingLists, unresolved, True)
                if len(unresolved) > 0:
                    print 'Damn, unconlictingList creation gone wrong'
    
        return unconflictingLists
    
    def _unConflictingList(self, list):
        return self._unconflictingList([[list.pop(0)]], list, False)
        
    def _unconflictingList(self, unconflictingLists, list, addUnresolvedToMultiple=False):
        unresolved = []
        for i in range(0,len(list)):
            temp = []
            for ii in unconflictingLists:
                temp.append(ii[:])
            
            insert = []
            for j in range(0,len(temp)):
                added = False
                conflict = False
                for k in temp[j]: 
                    if list[i][0] & k[0] != 0:
                        if list[i][0] & k[0] != list[i][0] and list[i][0] & k[0] != k[0]:
                            conflict = True
                if not conflict:
                    insert.append(j)

            if len(insert) > 0:
                if len(insert) > 1:
                    if addUnresolvedToMultiple:
                        for index in insert:
                            unconflictingLists[index].append(list[i])
                    else:
                        unresolved.append(list[i])
                else:
                    for index in insert:
                        unconflictingLists[index].append(list[i])
            else:
                unconflictingLists.append([list[i]])
        return unconflictingLists, unresolved
    
    def _sortListByHits(self, list):
        for i in range(0, len(list) -1):
            for j in range(i+1, len(list)):
                if list[i][2] < list[j][2]:
                    temp = list[i]
                    list[i] = list[j]
                    list[j] = temp
  
class TreeHandler(object):
    
    def __init__(self, inThing=None, skip=0, max=None):
        self.trees = []
        self.intersections = []
        self.nTrees = 0
        self.taxNames = None
        self.nTax = 0        # Merely len(taxNames)
        self.splits = []
        self.bitkeys = []
        self.splitsHash = {}
        # biSplits and biSplitsHash for the "even" side of the biRoot bifurcation.
        self.biSplits = []
        self.biSplitsHash = {}
        self.isBiRoot = None      # Is it a bifurcating root?
        self.doModelComments = 0
        self.modelInfo = None
        self.treesBuilt = False
        self.tfl = None
        self.verbose = 1
        if inThing:
            self.read(inThing, skip, max)
            
    def read(self, inThing, skip=0, max=None):
        """Read in a tree, or some trees, or a file of trees.

        Arg inThing can be a file name, a Tree instance, or a
        Trees instance.  Args skip and max only come into play if
        inThing is a file or a Trees object."""

        gm = ['TreePartitions.read()']

        if not inThing:
            gm.append("No input?")
            raise Glitch, gm
        #print "inThing = %s, type %s" % (inThing, type(inThing))
        if type(inThing) == type('string'):
            var.trees = []
            read(inThing)
            if len(var.trees) < 1:
                gm.append('Sorry, at least one tree must be supplied as input tree')
                raise Glitch, gm
            self.trees = var.trees
#            self.tfl = TreeFileLite(inThing)
        elif type(inThing) == type([]):
            for t in inThing:
                if not isinstance(t, Tree):
                    gm.append("Input trees should be a list of p4 Tree objects. Got %s" % t)
                    raise Glitch, gm
            self.trees = inThing
        else:
            gm.append("Sorry-- I can't grok your input.  What is '%s'?" % inThing)
            raise Glitch, gm
  
    def _printSplitList(self, list):
        for split in list:
#            for split in list:
            intersection = ''
            for bk in self.bitkeys:
                if bk & split[0]:
                    intersection = intersection + '1'
                elif bk & split[1]:
                    intersection = intersection + '?'
                else:
                    intersection = intersection + '*'
            print 'Split: %s' % (intersection)
        print ' '
 
    
    def buildTreesFromIntersections(self, intersections,minimumSplitsToBuildTree, treeLabel):
    
        print treeLabel
#        print 'Taxnumber to taxname translation: '
#        for i in range(0, len(self.taxNames)):
#            print '%s,%s' % ( i+1, self.taxNames[i])
            
        treeBuilder = TreeBuilderFromSplits(self.bitkeys, self.taxNames)
        
        conTrees = []
        for l in intersections:
            if len(l) >= minimumSplitsToBuildTree:
                l.sort()
                t = treeBuilder.buildTreeFromInformativeList(l)
                ex = 'Excluded taxa: '
                for j in l[0].excluded:
                    ex += self.taxNames[j] +', '
                if len(l[0].excluded) > 0:
                    print ex[0:len(ex)-2]
                else:
                    print ex + 'None'
#                sp = 'Included splits: '
                
                    
                t.draw(showInternalNodeNames=1, addToBrLen=0.2, width=None, showNodeNums=0, partNum=0, model=None)  
                print ''  
                conTrees.append(t)
        
#        for i in range(len(conTrees)):
#            if isinstance(conTrees[i], Tree):
#                ex = 'Exluded taxa: '
#                print 'len(conTrees[i].excluded): ', len(conTrees[i].excluded)
#                for j in conTrees[i].excluded:
#                    ex += self.taxNames[j] +', '
#                print ex[0:len(ex)-1]
#                self.printVerticalNumbers(len(self.taxNames))
#                for j in intersections[i]:
#                    j.printIntersection(self.bitkeys)
#                conTrees[i].draw(showInternalNodeNames=1, addToBrLen=0.2, width=None, showNodeNums=0, partNum=0, model=None)  
#                print ''     
 
        if len(conTrees) == 0:
            print 'No trees where produced'
            return []
        else:
            return conTrees
    
    def buildConsensusTreeHash(self, contree):
        if not contree._taxNames:
            contree._setTaxNamesFromLeaves()

#        t.taxNames.sort()
        dict = {}
        self.bitkeys = []
        for i in range(len(contree.taxNames)):
            self.bitkeys.append(1L << i)
            dict[contree.taxNames[i]] = 1L << i 
        
        self.taxNames = contree.taxNames
        
        
            
    
    def commonLeafSetTrees(self, trees):
        list = []
        
        if self.verbose:
            print 'Processing input trees'
            sys.stdout.flush()
        
            sampler = int(float(len(trees)/20))
            if sampler < 1:
                sampler = 1
        t = trees[0]
        if not t._taxNames:
            t._setTaxNamesFromLeaves()

        t.taxNames.sort()
        dict = {}
        self.bitkeys = []
        for i in range(len(t.taxNames)):
            self.bitkeys.append(1L << i)
            dict[t.taxNames[i]] = 1L << i 
        
        self.taxNames = t.taxNames
        
        self.weights = []
        splits = []
        treeNames = []
        index = 0
        if self.verbose:
            sys.stdout.write('0% ')
            sys.stdout.flush()
        for t in trees:
            index += 1
            if self.verbose:
                if index % sampler == 0:
                    sys.stdout.write('.')
                    sys.stdout.flush()
            treeNames.append(t.name)
            for n in t.iterLeavesNoRoot():
                n.br.rc = dict[n.name]
            t.splits = []
            t._makeRCSplitKeys(t.splits)
            weight = 1.0
            if hasattr(t, "weight"):
                if t.weight != None:
                    weight = t.weight
            elif hasattr(t, "recipWeight"):
                if t.recipWeight != None:
                    weight = 1.0/int(t.recipWeight)
            self.weights.append(weight)
            splits.append(t.splits)
        if self.verbose:
            sys.stdout.write(' 100%\n')
            sys.stdout.flush()
        
        return splits, treeNames
        

    
    def commonLeafSet(self, tfl):
        list = []
        
        print 'Processing input trees'
        sys.stdout.flush()
        
        sampler = int(float(tfl.nSamples)/20)
        if sampler < 1:
            sampler = 1
        t = tfl.getTree(1)
        if not t._taxNames:
            t._setTaxNamesFromLeaves()

        t.taxNames.sort()
        dict = {}
        self.bitkeys = []
        for i in range(len(t.taxNames)):
            self.bitkeys.append(1L << i)
            dict[t.taxNames[i]] = 1L << i 
        
        self.taxNames = t.taxNames
        
        self.weights = []
        splits = []
        treeNames = []
        index = 0
        sys.stdout.write('0% ')
        sys.stdout.flush()
        for i in range(tfl.nSamples):
            index += 1
            if index % sampler == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
            t = tfl.getTree(i)
            treeNames.append(t.name)
            for n in t.iterLeavesNoRoot():
                n.br.rc = dict[n.name]
            t.splits = []
            t._makeRCSplitKeys(t.splits)
            weight = 1.0
            if hasattr(t, "weight"):
                if t.weight != None:
                    weight = t.weight
            elif hasattr(t, "recipWeight"):
                if t.recipWeight != None:
                    weight = 1.0/int(t.recipWeight)
            self.weights.append(weight)
            splits.append(t.splits)
        
        sys.stdout.write(' 100%\n')
        sys.stdout.flush()
        
        return splits, treeNames
        
        
        
    def updateToCommonLeafSet(self, tfl):
        uniqueSet = set()
        
        print 'Processing input trees'
        
        sampler = int(float(tfl.nSamples)/20)
        if sampler < 1:
            sampler = 1
        
        index = 0.0
        sys.stdout.write('0% ')
        for i in range(tfl.nSamples):
            index += 0.5
            if index % sampler == 0:
                sys.stdout.write('.')
#                print '%s trees processed ' % (int (index))
                sys.stdout.flush()
            t = tfl.getTree(i)
#            t.setPreAndPostOrder()
            if not t._taxNames:
                t._setTaxNamesFromLeaves()
            for name in t.taxNames:
                uniqueSet.add(name)
        self.list = []
        for name in uniqueSet:
            self.list.append(name)
        
        self.list.sort()
        
        dict = {}
        self.bitkeys = []
        for i in range(len(self.list)):
            self.bitkeys.append(1L << i)
            dict[self.list[i]] = 1L << i 
        
        self.taxNames = self.list

        self.weights = []
        splits = []
        treeNames = []
        for i in range(tfl.nSamples):
            
            index += 0.5
            if index % sampler == 0:
                sys.stdout.write('.')
#                print '%s trees processed' % (int(index))
                sys.stdout.flush()
            t = tfl.getTree(i)
            treeNames.append(t.name)
            if not t._taxNames:
                t._setTaxNamesFromLeaves()
            t.missingTaxa = 0L
            for name in self.list:
                if not t.taxNames.count(name):
                    t.missingTaxa = t.missingTaxa | dict[name]
            for n in t.iterLeavesNoRoot():
                n.br.rc = dict[n.name]
#            t.setPreAndPostOrder()
            t._makeRCSplitKeys()
            weight = 1.0
            if hasattr(t, "weight"):
                if t.weight != None:
                    weight = t.weight
            elif hasattr(t, "recipWeight"):
                if t.recipWeight != None:
                    weight = 1.0/int(t.recipWeight)
            self.weights.append(weight)
            t.splits = []
            for n in t.iterInternalsNoRoot():
                for m in n.br.rcList:
                    t.splits.append([m,t.missingTaxa])
            splits.append(t.splits)
        
        sys.stdout.write(' 100%\n')
        
        return splits, treeNames
            
    def writeNexus(self, trees,  filename='rcTrees'):
        """
        This method will simply write the trees that were produced by calling reduced()
        If a filename is provided the trees will be written to a number of files starting
        with filename1.nex and then filename2.nex and so on and so forth. This is done to avoid 
        mixing trees with different leafsets in the same nexus file as this will make it hard
        to open the file in other software.
        Beware that if no filename is provided the default filename will be used and might therefore overwrite
        already completed analysis.
        
        """
        
        if len(trees) == 0:
            print 'No trees created, can not write to file'
            return
        
        for i in range(len(trees)):
            trees[i].writeNexus('%s%d%s' % (filename,i+1,'.nex'))
            
    def first(self, split):
        first = 0.0
        for bk in self.bitkeys:
            if bk & split:
                return first
            first += 1
        return first
    
    def popcount(self, n):
        count = 0
        for bk in self.bitkeys:
            if n < bk:
                return count
            if bk & n:
                count += 1
        return count

    def _reduceSplitsPair(self, reduceRules, splits1, splits2):
        newSplits1 = []
        newSplits1.extend(splits1)
        newSplits2 = []
        newSplits2.extend(splits2)

        splitList = []
        splitList.append(newSplits1)
        splitList.append(newSplits2)

        excludedTaxa = []
        for row in csvReader:
#            print row
            for name in row:
                if not taxa2bitkey.has_key(name):
                    print 'Taxname not in any tree: ', name
                    return
            if len(row) > 1:
                reduce2 = taxa2bitkey[row[0]]
                substitute = 0L
                for i in range(1,len(row)):
                    excludedTaxa.append(taxa2bitkey[row[i]])
                    substitute = substitute | taxa2bitkey[row[i]]
                
                for i in range(len(splitList)):
                    
#                    if i == 239:
#                        print 'Before'
#                        for split in splitList[i]:
#                            Intersection(split[0], split[1]).printIntersection(self.bitkeys)
                    for split in splitList[i]:
                        if split[0] & substitute != 0:
                            split[0] = (split[0] & substitute) ^ split[0]
                            split[0] = split[0] | reduce2
                        if split[1] & reduce2 != 0:
                            split[1] = split[1] ^ reduce2
                        split[1] = split[1] | substitute
#                    if i == 239:
#                        print 'After'
#                        for split in splitList[i]:
#                            Intersection(split[0], split[1]).printIntersection(self.bitkeys)
#                        sys.exit()
                            
        excluded = 0L
        for bitkey in excludedTaxa:
            excluded = excluded | bitkey
        return excluded
        

        return newSplits1, newSplits2

    def _reduceSplits(self, csvReader, splitList):
        
        bitkey2taxa = {}
        taxa2bitkey = {}
        for i in range(0,len(self.taxNames)):
            taxa2bitkey[self.taxNames[i]] = self.bitkeys[i]
            bitkey2taxa[self.bitkeys[i]] = self.taxNames[i]
   
        excludedTaxa = []
        for row in csvReader:
#            print row
            for name in row:
                if not taxa2bitkey.has_key(name):
                    print 'Taxname not in any tree: ', name
                    return
            if len(row) > 1:
                reduce2 = taxa2bitkey[row[0]]
                substitute = 0L
                for i in range(1,len(row)):
                    excludedTaxa.append(taxa2bitkey[row[i]])
                    substitute = substitute | taxa2bitkey[row[i]]
                
                for i in range(len(splitList)):
                    
#                    if i == 239:
#                        print 'Before'
#                        for split in splitList[i]:
#                            Intersection(split[0], split[1]).printIntersection(self.bitkeys)
                    for split in splitList[i]:
                        if split[0] & substitute != 0:
                            split[0] = (split[0] & substitute) ^ split[0]
                            split[0] = split[0] | reduce2
                        if split[1] & reduce2 != 0:
                            split[1] = split[1] ^ reduce2
                        split[1] = split[1] | substitute
#                    if i == 239:
#                        print 'After'
#                        for split in splitList[i]:
#                            Intersection(split[0], split[1]).printIntersection(self.bitkeys)
#                        sys.exit()
                            
        excluded = 0L
        for bitkey in excludedTaxa:
            excluded = excluded | bitkey
        return excluded
    
    def _sortCompatibleSplits(self, compatibleSplits):
        
        pops = []
        i = 0
        for split in compatibleSplits:
            pop = self.popcount(split[0])
            if pop > 1:
                first = self.first(split[0])/len(self.bitkeys)
                pop = pop - first
                pops.append([pop,i])
            i += 1
        
        
#        print 'unsorted', pops
        pops = sorted(pops, key=operator.itemgetter(0), reverse=True)
#        print 'sorted', pops
        
        sortedList = []
#        print 'Sorting indices'
        for i in range(0, len(pops)):
#            print pops[i][1]
            sortedList.append(compatibleSplits[pops[i][1]])
        
        return sortedList
        
class SemiStrictSuperTree(TreeHandler):

    def __init__(self, inThing=None, skip=0, max=None):
        TreeHandler.__init__(self, inThing, skip, max)    
    
    def semiStrictSuperTree(self, semiStrictOnly=True, treatMultifurcatingRootsAsUnrooted=0, readCSV=False, csvFile='clades.csv'):
        """
        Creates a semistrict supertree from a set of trees. Setting the semiStricOnly to False will
        cause the algorithm to try and resolve any unresolved nodes by using majority rule consensus. 
        
        Usage: 
        
        rd = Reduced('filename')
        rd.semiStrictSuperTree()
        
        to safe supertree to file:
        
        rd.writeNexus('outFilename')
        
        """
        self.conTrees = []
        
        self.splits, treeNames = self.updateToCommonLeafSet(self.tfl)
        
        if readCSV:
            csvReader = csv.reader(open(csvFile), delimiter=',', quotechar='"')
        
            excludedTaxa = self._reduceSplits(csvReader, self.splits)
        else:
            excludedTaxa = 0L
#        print 'Excluded taxa: ', excludedTaxa
        
        noSplitSets = len(self.splits)
        
#        Found no splits, exiting    
        if len(self.splits) < 1:
            print 'No splits found'
            return
        
#       Found only one set of splits, as the set is from a single tree it is assumed to be compatible 
        elif len(self.splits) == 1:
            print 'Only one set of splits'
            treeBuilder = TreeBuilderFromSplits(self.bitkeys, self.taxNames)
            intersections = []
            self._printSplitList(self.splits[0])
            for split in self.splits[0]:
                intersections.append(Intersection(split[0],split[1], 1, 1))
#            print self.splits
#            print self.splits[0]
            mrc = treeBuilder.buildTreeFromInformativeList(intersections)
            mrc.draw()
            print ' '
            return
        
        print 'Creating compatible split set, collapsing where necessary'
        
        set = []
        for t in range(0,len(self.splits)):
            for split in self.splits[t]:
                set.append(split)
        print 'Initial splits'
        self._printSplitList(set)
        
        overlappingLists = self._buildOverlappingLists(set)
        
        print 'Overlapping lists: ', len(overlappingLists)
        for list in overlappingLists:
            self._printSplitList(list)
        
        combinedLists = []
        for list in overlappingLists:
            combinedLists.append(self._combineList(list))
        
        print 'Combined overlapping lists'
        for list in combinedLists:
            self._printSplitList(list)    
        
        
#        print 'Creating compatible splits list'
        compatibleSplits = []
        for list in combinedLists:
            for split in list:
                compatibleSplits.append(split)
        
        uniqueSet = {}
        for split in compatibleSplits:
            if not uniqueSet.has_key((split[0],split[1])):
                uniqueSet[(split[0],split[1])] = 1
        
#        self._isListCompatible(uniqueSet.keys())
        
        compatibleSplits = []
        for split in uniqueSet.keys():
            compatibleSplits.append([split[0], split[1]])
        
        compatibleSplits = self._sortCompatibleSplits(compatibleSplits)
        
        print 'Sorted compatible splits'
        self._printSplitList(compatibleSplits)
        
                
        if semiStrictOnly:
#        To conform to the format used by treeBuilderFromSplits we need to convert splits into intersections
#        We use a set to remove any duplicate splits created during the refinement
        
#            print '1. Looking for incompatibility in fixed compatible splits'
            
                
            compatibleIntersections = []
            for s in compatibleSplits:
                i = Intersection(s[0],s[1], 1, 1)
                compatibleIntersections.append(i)
                
            if False:
                print 'List of unique splits'
                for s in compatibleSplits:
                    self._printSplitList([s])
            
            print 'Building tree'
            treeBuilder = TreeBuilderFromSplits(self.bitkeys, self.taxNames)
            mrc = treeBuilder.buildTreeFromInformativeList(compatibleIntersections, excluded=excludedTaxa)
#            mrc.dump()
            mrc.draw()
            print ' '
            self.treesBuilt = True
            self.conTrees.append(mrc)
            return mrc

    def _compatible(self, split1, split2):

        if split1[0] & split2[0] == 0:
#            print '1' 
#            self._printSplitList([split1])
#            self._printSplitList([split2])
            return True
        
        if split1[0] | split2[0] == split1[0] or split1[0] | split2[0] == split2[0]:
#            print '2' 
#            self._printSplitList([split1])
#            self._printSplitList([split2])
            return True
        
        if (split1[0] | split1[1]) | split2[0] == (split1[0] | split1[1]):
#            print '3' 
#            self._printSplitList([split1])
#            self._printSplitList([split2])
            return True
        
        if (split2[0] | split2[1]) | split1[0] == (split2[0] | split2[1]):
#            print '4' 
#            self._printSplitList([split1])
#            self._printSplitList([split2])
            return True
        
        if (split1[0] | split1[1]) | (split2[0] | split2[1]) == (split1[0] | split1[1]):
#            print '5' 
#            self._printSplitList([split1])
#            self._printSplitList([split2])
            return True
        
        if (split1[0] | split1[1]) | (split2[0] | split2[1]) == (split2[0] | split2[1]):
#            print '6' 
#            self._printSplitList([split1])
#            self._printSplitList([split2])
            return True
        
        if (split1[0] | (split2[0] | split2[1]) != (split1[0] | split1[1]) and split1[0] | (split2[0] | split2[1]) != (split2[0] | split2[1])):
#            print '7' 
#            self._printSplitList([split1])
#            self._printSplitList([split2])
            return False
        
        print 'Defaults to true, not good'
        print split1
        print split2
        print ' ' 
        return True 
    
    def _isListCompatible(self, list):
        for i in range(0, len(list)-1):
            for j in range(i+1, len(list)):
                if not self._compatible(list[i], list[j]):
                    print 'incompatible'
                    self._printSplitList([list[i]])
                    self._printSplitList([list[j]])
    

    def _buildOverlappingLists(self, list):
        
        overlappingList = [[list.pop()]]
        firstTime = True
        while len(list) > 0:
#            print 'while 1'
            i = 0
            while i <= 1:
                for oList in overlappingList:
                    remaining = []
                    for split in list:
                        add = False
                        for s in oList:
                            if s[0] & split[0] != 0:
                                add = True
                                break
                            
                        if add:
                            oList.append(split)
                        else:
                            remaining.append(split)
                    list = remaining
                i += 1
                
            if len(list) > 0:
                overlappingList.append([list.pop()])
#                
#            if len(list) == 1 and not firstTime:
#                overlappingList.append([list.pop()])
#            if len(list) == 1 and firstTime:
#                firstTime = False
#            
        return overlappingList
        
    def _combineList(self, list):
        
        allOnes = 2L**(len(self.bitkeys)) - 1    
        
        
#        zero will hold the unknown taxa off all splits in list
        zero = 0L
        for split in list:
            zero = zero | split[1]
#        print 'Zero: '
        self._printSplitList([[zero,0L]])
            
        xor = zero ^ allOnes
#        print 'Xor: '
        self._printSplitList([[xor, 0L]])
        
        includeSplit = []
        if len(list) > 2:
            for split in list:
#                self._printSplitList([split])
#            split[1] = split[1] & xor
#            self._printSplitList([split])
#            print '*********************************'
                includeSplit.append(False)
#            print '*********************************'    
        else:
            for split in list:
                includeSplit.append(True)
                
                
        includesCombined = {}    
            
        for i in range(0, len(list)-1):
            for j in range(i+1, len(list)):
                
#                self._printSplitList([list[i]])
#                self._printSplitList([list[j]])
                
                if list[i][0] & list[j][0] != 0 and list[i][0] | list[j][0] != list[i][0] and list[i][0] | list[j][0] != list[j][0] :
                    combined = list[i][0] | list[j][0]
                    if self.popcount(list[i][0] | list[i][1]) >= self.popcount(list[j][0] | list[j][1]):
#                        list[i][0] = combined
#                        list[i][1] = list[i][1] & xor
                        if includesCombined.has_key(list[i][0]):
                            includesCombined[list[i][0]] = -1
                        includesCombined[combined] = list[i][1] & xor
                    else:
#                        list[j][0] = combined
#                        list[j][1] = list[j][1] & xor
                        if includesCombined.has_key(list[j][0]):
                            includesCombined[list[j][0]] = -1
                        
                        includesCombined[combined] = list[j][1] & xor
#                    self._printSplitList([list[i]])
#                    self._printSplitList([list[j]])
#                    print '------------------------------------'
                else:
#                    if (list[i][0] | list[i][1]) | (list[j][0]) == list[j][0] or (list[j][0] | list[j][1]) | list[i][0] == list[i][0]:
                    includeSplit[i] = True
                    includeSplit[j] = True
              
#                if not self._compatible(list[i], list[j]):
#                    self._printSplitList([list[i]])
#                    self._printSplitList([list[j]])
#                    print 'incompatible****************************************************************************'

        newList = []
        for i in range(0,len(includeSplit)):
            if includeSplit[i]:
                list[i][1] = list[i][1] & xor
                newList.append(list[i])
        
#        print includesCombined
        for i in includesCombined.keys():
            if includesCombined[i] >= 0:
                newList.append([i, includesCombined[i]])
        
        return newList
            
class CompatibleTreePairs(TreeHandler):
    
    def __init__(self, inThing=None, skip=0, max=None):
        TreeHandler.__init__(self, inThing, skip, max)    
            
    def slidingWindow(self, csvInfile, csvOutFile, startYear, stopYear, maxDistance):
        
        csvReader = csv.reader(open(csvInfile), delimiter=',', quotechar='"')
        csvWriter = csv.writer(open(csvOutFile, "wb"))
        
        yearDict = {}
        first = True
        index = 0
        added = 0
        rows = []
        for row in csvReader:
            rows.append(row)
            if first:
                first = False
            else:
                index += 1
                name = str(row[0])
                for i in range(startYear, stopYear+1):
                    if name.rfind(str(i)) >=0 and name.rfind('SCC') < 0:
                        if yearDict.has_key(i):
                            added += 1
                            yearDict[i].append([index, row])
                        else:
                            added += 1
                            yearDict[i] = [[index,row]]
                        break
                        
        
#        print 'Counted rows: ', index
#        print 'Added rows: ', added
#        print ''
        
        ranges = []
        for index in range(startYear, stopYear-4):
#            print 'Range: %s, %s ' % (index, index+4)
            studies = []
            for i in range(index, index+5):
                if yearDict.has_key(i):
                    studies.extend(yearDict[i])
                
#            for studie in studies:
#                print studie[1][0]
#                print studie
                
#            print ''
#            print '[%s, %s, %'
            ranges.append([index,index+4,studies])
        
#        print ''
        outputRows = []
        for dataset in ranges:
            row = [dataset[0], dataset[1]]
            print '%s, %s' % (dataset[0], dataset[1])
            indecies = []
            print dataset
            for data in dataset[2]:
                indecies.append(data[0])
            print indecies
           
            
            sum = 0.0
            index = 0
            for data in dataset[2]: 
                for data1 in dataset[2]:
                    if data[0] != data1[0]:
                        if data[1][data1[0]] != '*':
                            index += 1
#                            print data[1][data1[0]]
                            sum += float(data[1][data1[0]])/maxDistance
            if sum == 0:
                row.append(0)
                print 'Mean: ', 0
                print ''
            else:
                row.append(sum/index)
                print 'Mean: ', sum/index
                print ''
                

            outputRows.append(row)
        
        csvWriter.writerows(outputRows)
    
    def formatNexusFileFromCSV(self, csvInFile, nexusOutFile, maxDistance):
        
        csvDataReader = csv.reader(open(csvInFile), delimiter=',', quotechar='"')
        data = []
        nTax = -1
        for row in csvDataReader:
            nTax += 1
            data.append(row)


        writer = open(nexusOutFile, 'w')
        writer.write('#NEXUS\n')
        writer.write('[!user defined distances]\n')
        writer.write('Begin distances;\n')
        dim = 'Dimensions ntax=%s;\n' % (str(nTax))
        writer.write(dim)
        writer.write('format nodiagonal;\n')
        writer.write('matrix\n')
        first = True
        for index in range(0,len(data)):
            if first:
                first = False
            else:
                print data[index][0:index]
                line = ''
                for datapoint in data[index][0:index]:
                    if datapoint == '*':
                        line += '0.0'
                    else:
                        line += datapoint
                    line += ' '
                if index == len(data)-1:
                    line += ';'
                line += '\n'
                print line
                writer.write(line)
                
        writer.write('end;\n')
        writer.write('[! nj with user defined distances]\n')
        writer.write('Begin paup;\n')
        writer.write('dset distance=user;\n')
        writer.write('nj;\n')
        writer.write('end;\n')

        writer.close() 
                
         
        
            
    def calcStatsFromCSV(self, csvInFileIntraGroups, csvInFileInterGroups, csvInFileData, csvOutFile, maxDistance):
        
        csvInterGroupReader = csv.reader(open(csvInFileInterGroups), delimiter=',', quotechar='"')
        
        csvIntraGroupReader = csv.reader(open(csvInFileIntraGroups), delimiter=',', quotechar='"')
        
        csvDataReader = csv.reader(open(csvInFileData), delimiter=',', quotechar='"')
        
        data = []
        for row in csvDataReader:
            if len(row) > 0:
                data.append(row)
        
        length = 0
        if len(data) > 0:
            length = len(data[0])
        
        intraGroups = {}
        for row in csvIntraGroupReader:
            if len(row) > 0:
#                print 'Adding:', row[0]
                intraGroups[row[0]] = [length,0]    
#        print intraGroups.keys()
        
        interGroups = {}
        for row in csvInterGroupReader:
            if len(row) > 0:
                interGroups[row[0]] = [0]        
#        print interGroups.keys()

        for index in range(1,len(data)):
            if interGroups.has_key(data[index][0]):
                interGroups[data[index][0]] = index
#                print 'Inter: %s, %s ' % (data[index][0], index)
            else:
                if intraGroups.has_key(data[index][0]):
#                    print 'data[index][0]:', data[index][0]
#                    print 'data[index][1]:', data[index][1]
#                    print 'intraGroups[data[index][0]]:',intraGroups[data[index][0]] 
#                    print '%s : %s,%s' % (data[index][0], index, index)
                    intraGroups[data[index][0]] = [index,index]
#                    if index < intraGroups[data[index][0]]:
                        
#                        intraGroups[data[index][0]] = [index,intraGroups[data[index][1]]]
#                    if index > intraGroups[data[index][1]]:
#                        intraGroups[data[index][0]] = [intraGroups[data[index][0],index]]
                else:
                    start = data[index][0].find('mpt')
                    substring = ''
                    if start > 0:
                        substring = data[index][0].partition('mpt')[0]
                    else:
                        start = data[index][0].find('SCC')
                        if start > 0:
                            substring = data[index][0].partition('SCC')[0]
#                    print 'Created substring:', substring
                    if intraGroups.has_key(substring):
                        if index < intraGroups[substring][0]:
                            intraGroups[substring] = [index,intraGroups[substring][1]]
#                        print 'Intra1: %s, %s - %s ' % (substring, index, intraGroups[substring][1])
                        elif index > intraGroups[substring][1]:
                            intraGroups[substring] = [intraGroups[substring][0],index]
#                        print 'Intra2: %s, %s - %s ' % (substring, intraGroups[substring][0], index ) 
                    else:
                        print 'Found nowhere to put: ', data[index][0]
                    
        writer = csv.writer(open(csvOutFile, "wb"))
        outputRows = []
        outputRows.append(['Real trees intra group distances, normalized, mean distance'])
        outputRows.append(['Study','Mean'])
#        print 'Study: max, min, mean'
        def compareStudyYear(x, y):
#            print '%s, %s' % (x[0], x[0][len(x[0])-4:len(x[0])])
#            print '%s, %s' % (y[0], y[0][len(y[0])-4:len(y[0])])
#            print ''
            return int(x[0][len(x[0])-4:len(x[0])]) - int(y[0][len(y[0])-4:len(y[0])]) 

        items = intraGroups.items()
        items.sort(cmp = compareStudyYear)
#        print items
        for k,v in items:
#            max, min, mean = self.calcIntraMinMaxMean(data, v[0], v[1])
            outputRows.append([k,self.calcIntraMinMaxMean(data, v[0], v[1], maxDistance)])
#            print '%s: %s, %s, %s' % (k, max, min, mean)
        
        outputRows.append([' '])
        outputRows.append(['Real trees inter group distances, normalized, mean distance'])
        items = intraGroups.items()
        items.sort(cmp = compareStudyYear)
#        print 'Study: max, min, mean'
        row = ['']
        for i in items:
            row.append(i[0])
        outputRows.append(row)
        
        for i in range(0, len(items)):
            row = [items[i][0]]
            for j in range(0, len(items)):
                row.append(self.calcIntraInterMinMaxMean(data, items[i][1][0], items[i][1][1], items[j][1][0], items[j][1][1], maxDistance))
#            max, min, mean = self.calcIntraMinMaxMean(data, v[0], v[1])
            
            outputRows.append(row)
#            print '%s: %s, %s, %s' % (k, max, min, mean)
        
        
        items = interGroups.items()
        items.sort(cmp = compareStudyYear)
        outputRows.append([' '])
        outputRows.append(['Expertograms intra tree distances, normalized, distance'])
        row = ['']
        for i in items:
            row.append(i[0])
        outputRows.append(row)
        
        for i in range(0, len(items)):
            row = [items[i][0]]
            for j in range(0, len(items)):
#                print items[i][1]
#                print items[j][1]
                if data[items[i][1]][items[j][1]] == '*':
                    row.append(data[items[i][1]][items[j][1]])
                else:
                    row.append(float(data[items[i][1]][items[j][1]])/maxDistance)
            outputRows.append(row)
        
        
        intraItems = intraGroups.items()
        intraItems.sort(cmp = compareStudyYear)
        outputRows.append([' '])
        outputRows.append(['Expertograms to real trees inter group distances, normalized, mean distance'])
        row = ['']
        for i in intraItems:
            row.append(i[0])
        outputRows.append(row)
        
        for key, value in items:
            row = []
            row.append(key)
            for d,v in intraItems:
                row.append(self.calcInterMinMaxMean(data, value, v[0], v[1], maxDistance))
            outputRows.append(row)
    
        writer.writerows(outputRows)
        
    def median(self, numbers):
#        "Return the median of the list of numbers."
        # Sort the list and take the middle element.
        n = len(numbers)
        copy = numbers[:] # So that "numbers" keeps its original order
        copy.sort()
        if n & 1:         # There is an odd number of elements
            return copy[n // 2]
        else:
            return (copy[n // 2 - 1] + copy[n // 2]) / 2

        
    def calcInterMinMaxMean(self, data, inter, start, stop, maxDistance):
        min = 99999999999
        max = 0
        sum = 0.0
        index = 0
        
#        print 'range: ', range(start, stop+1)
        
        for i in range(start, stop+1):
            if data[inter][i] == '*':
#                return '*/*/*'
                return '*'
#            print '%s,%s' % (inter,i)
            dist = float(data[inter][i])/maxDistance 
            if  dist < min:
                min = dist
            if dist > max:
                max = dist
            sum += dist
            index += 1
        
#        return str(max) +'/' + str(min) +'/' + str(sum/index)
        return sum/index
        
        
    def calcIntraInterMinMaxMean(self, data, start1, stop1, start2, stop2, maxDistance):
        min = 99999999999.0
        max = 0.0
        sum = 0.0
        index = 0.0
    
        for i in range(start1, stop1+1):
            for j in range(start2, stop2+1):
                if i != j:
#                print '%s, %s' % (data[i][0], data[j][0])
                    if data[i][j] == '*':
#                    return '*','*','*'
                        return '*'
#                print '%s,%s' % (i,j)
                    dist = float(data[i][j])/maxDistance
                    if dist < min:
                        min = dist
                    
                    if dist > max:
                        max = dist
                    sum += dist
                    index += 1

#        return max, min, sum/index
        if sum != 0:
            return sum/index
        return 0
        
    def calcIntraMinMaxMean(self, data, start, stop, maxDistance):
        min = 99999999999.0
        max = 0.0
        sum = 0.0
        index = 0.0
        
#        print range(start, stop+1)
        
        for i in range(start, stop+1):
            for j in range(i+1, stop+1):
#                print '%s, %s' % (data[i][0], data[j][0])
                if data[i][j] == '*':
#                    return '*','*','*'
                    return '*'
#                print '%s,%s' % (i,j)
                dist = float(data[i][j])/maxDistance
                if dist < min:
                    #===========================================================
                    # print 'MIN: %s < %s' % (float(data[i][j]) ,min)
                    #===========================================================
                    min = dist
                    
                if dist > max:
                    #===========================================================
                    # print 'MAX: %s > %s' % (float(data[i][j]) ,max)
                    #===========================================================
                    max = dist
                sum += dist
                index += 1
                
#        print 'Max: ', max
#        print 'Min: ', min
#        print 'Sum: ', sum
#        print 'Index: ', index
#        print 'Mean: ', sum/index
        
#        return max, min, sum/index
        if sum != 0:
            return sum/index
        return 0
            
    def createCompatibleTreePairs(self, readCSV=True, csvInFile='clades.csv', writeCSV=False, csvOutFile='treeMatrix.csv', writeNexus=True, nexusFilename='compatiblePair', calcDistances = False, distanceMeasure = 'sd'):
        
        from Trees import Trees
        
        self.splits, treeNames = self.updateToCommonLeafSet(self.tfl)
        
        self.bitkey2taxa = {}
        self.taxa2bitkey = {}
        for i in range(0,len(self.taxNames)):
            self.taxa2bitkey[self.taxNames[i]] = self.bitkeys[i]
            self.bitkey2taxa[self.bitkeys[i]] = self.taxNames[i]
        
        if readCSV:
            
            csvReader = csv.reader(open(csvInFile), delimiter=',', quotechar='"')
            self.replace = []
            self.reduceRules = []
            
            for row in csvReader:
                self.reduceRules.append(row)
                self.replace.append(self.taxa2bitkey[row[0]])
   
        if writeCSV:
            writer = csv.writer(open(csvOutFile, "wb"))
        
        treeBuilder = TreeBuilderFromSplits(self.bitkeys, self.taxNames, internalNames=False)
        
#        bitkey2taxa = {}
#        taxa2bitkey = {}
#        for i in range(0,len(self.taxNames)):
#            taxa2bitkey[self.taxNames[i]] = self.bitkeys[i]
#            bitkey2taxa[self.bitkeys[i]] = self.taxNames[i]
        
        index = 0
        
        csvOut = []
        if writeCSV:
            names = []
            names.append('')
            for i in range(0,len(treeNames)):
                names.append(treeNames[i])
            csvOut.append(names)

        if calcDistances: 
            max = 0.0
            min = 99999999999999.0
            sum = 0.0
            distances = 0
            if distanceMeasure not in ['sd', 'triplet', 'quartet']:
                print 'No such distance measure: ', distanceMeasure
                print 'Defaulting to Symmetric Difference'
                distanceMeasure = 'sd'
                
            print 'Distance measure: ', distanceMeasure

        dict = {}
        nonUnique = {}

        for i in range(0,len(self.splits)):
            row = []
            if writeCSV:
                row.append(treeNames[i])
            for j in range(0,len(self.splits)):
                if j >= i:
                    if i != j:
                        
                        set1, set2, taxnames, pairExcluded = self._commonTaxaFromTreePair(self.splits[i], self.splits[j])
                
                        if len(set1) > 0 and len(set2) > 0:
                            
                            uniqueSet = {}
                            for split in set1:
                                if not uniqueSet.has_key((split[0],split[1])):
                                    uniqueSet[(split[0],split[1])] = 1
                    
                            list = []
                            for s in uniqueSet.keys():
                                list.append([s[0],s[1]])
                        
                            list = self._sortCompatibleSplits(list)
                    
                            intersections1 = []
                            for s in list:
                                intersections1.append(Intersection(s[0],s[1], 1, 1))
                    
                            mrc = treeBuilder.buildTreeFromInformativeList(intersections1, treeNames[i])

                            mrc.taxNames = taxnames
                            uniqueSet = {}
                            for split in set2:
                                if not uniqueSet.has_key((split[0],split[1])):
                                    uniqueSet[(split[0],split[1])] = 1
                    
                            list = []
                            for s in uniqueSet.keys():
                                list.append([s[0],s[1]])
                            
                            list = self._sortCompatibleSplits(list)
                            
                            intersections2 = []
                            for s in list:
                                intersections2.append(Intersection(s[0],s[1], 1, 1))
                    
                            mrc2 = treeBuilder.buildTreeFromInformativeList(intersections2, treeNames[j])

                            mrc2.taxNames = taxnames
                            
                            if calcDistances:
                                if mrc.root.getNChildren() <= 2:
                                    mrc.removeRoot()
                                if mrc2.root.getNChildren() <= 2:
                                    mrc2.removeRoot()
                                
                                dist= mrc.topologyDistance(mrc2, distanceMeasure)
#                                dist = float(dist)/setSize
                                dict[(j,i)] = float(dist)
                                
                                if dist > max:
                                    max = dist
                                
                                if dist < min:
                                    min = dist
                                    
                                sum += dist
                                
                                
#                                if dist > 0:
                                distances += 1
                                
                                if dist == 0.0:
                                    study = 'Clark2004'
                                    if mrc.name.startswith(study) and mrc2.name.startswith(study): 
                                        print '%s, %s: %s' % (mrc.name, mrc2.name, dist)
                                        print
                                
                                if writeCSV:
                                    row.append(dist)
                                    
                                if dist == 0:
                                    nonUnique[mrc2.name] = 1
                                    sys.stdout.flush()
                        
                            if writeNexus:
                                trees = Trees([mrc,mrc2], taxnames)
                                trees.writeNexus(fName='%s%d%s' % (nexusFilename,index+1,'.nex'), append=0, withTranslation=1, writeTaxaBlock=1, likeMcmc=0)

                            index += 1
                    
                        else:
                            dict[(j,i)] = '*'
                            if writeCSV:
                                row.append('*')
                    else:
                        row.append('0')
                else:
                    row.append(dict[(i,j)])
            if writeCSV:
                csvOut.append(row)
                        
        print ' '
        
        if calcDistances: 
            print 'Max: ', max
            print 'Min: ', min
            print 'Mean: ', sum/distances
            
        if writeCSV:
            writer.writerows(csvOut)
            print 'Wrote csv-data to file: ', csvOutFile
            
        
        if writeNexus: 
            print 'Wrote tree data to nexus files on format %s ' % (nexusFilename + '#' + '.nex')
        
            
        nonUniquetrees = nonUnique.keys()
        
#        print sorted(nonUniquetrees)
        
    def _commonTaxaFromTreePair(self, splits1, splits2):
        
        replace = False
        for i in self.replace:
#            for split in splits1:
            if not splits1[0][1] & i:
                replace = True
                break
#            if replace:
#                break
        replace2 = False
        for i in self.replace:
#            for split in splits2:
            if not splits2[0][1] & i:
                replace2 = True
                break
#            if replace2:
#                break
        
        if (replace and not replace2) or (not replace and replace2):
#            print '=========================================================='
            splits1, splits2 = self._reduceSplitsPair(splits1, splits2)
        
        
        included1 = 2L**(len(self.bitkeys)) - 1
        included1 = included1 ^ splits1[0][1]
            
        included2 = 2L**(len(self.bitkeys)) - 1
        included2 = included2 ^ splits2[0][1]

        combined = included1 & included2
        excluded = 2L**(len(self.bitkeys)) - 1 ^ combined
        
        set1 = []
        set2 = []
        taxnames = []
        
        if self.popcount(combined) >= 3:
            
            for i in self.bitkeys:
                if i & combined:
                    taxnames.append(self.bitkey2taxa[i])
            
            temp = {}
            for split in splits1:
                if self.popcount(split[0] & combined) > 1 and ((split[0] & combined) ^ (split[1] | excluded)) ^ (2L**(len(self.bitkeys)) - 1) !=0:
                    temp[(split[0] & combined, split[1] | excluded)] = 1

            for i in temp.keys():
                set1.append([i[0],i[1]])
            temp = {}
            for split in splits2:
                if self.popcount(split[0] & combined) > 1 and ((split[0] & combined) ^ (split[1] | excluded)) ^ (2L**(len(self.bitkeys)) - 1) !=0:
                    temp[(split[0] & combined, split[1] | excluded)] = 1

            for i in temp.keys():    
                set2.append([i[0], i[1]])
            
        return set1,set2, taxnames, excluded
        
    def _reduceSplitsPair(self, splits1, splits2):
        newSplits1 = []
        for i in splits1:
            newSplits1.append([i[0], i[1]])
        
        newSplits2 = []
        for i in splits2:
            newSplits2.append([i[0], i[1]])
        
        splitList = []
        splitList.append(newSplits1)
        splitList.append(newSplits2)

        excludedTaxa = []
        for row in self.reduceRules:
            present = []
            for name in row:
                if self.taxa2bitkey.has_key(name):
                    present.append(name)
#                else:
#                    print 'Taxname not in any tree: ', name

            if len(present) > 1:
                reduce2 = self.taxa2bitkey[present[0]]
                substitute = 0L
                for i in range(1,len(present)):
                    excludedTaxa.append(self.taxa2bitkey[present[i]])
                    substitute = substitute | self.taxa2bitkey[present[i]]
                
                for i in range(len(splitList)):
                    for split in splitList[i]:
                        if split[0] & substitute != 0:
                            split[0] = (split[0] & substitute) ^ split[0]
                            split[0] = split[0] | reduce2
                        if split[1] & reduce2 != 0:
                            split[1] = split[1] ^ reduce2
                        split[1] = split[1] | substitute
                            
        excluded = 0L
        for bitkey in excludedTaxa:
            excluded = excluded | bitkey
#        return excluded

        return splitList[0], splitList[1]
            
    def _buildAndDisplayTree(self, set1):
        treeBuilder = TreeBuilderFromSplits(self.bitkeys, self.taxNames, internalNames=False)
        uniqueSet = {}
        for split in set1:
            if not uniqueSet.has_key((split[0],split[1])):
                uniqueSet[(split[0],split[1])] = 1
                    
        list = []
        for s in uniqueSet.keys():
            list.append([s[0],s[1]])
                        
        list = self._sortCompatibleSplits(list)
                    
        intersections1 = []
        for s in list:
            intersections1.append(Intersection(s[0],s[1], 1, 1))
                    
        mrc = treeBuilder.buildTreeFromInformativeList(intersections1, 'name')
        
        mrc.draw()
            

             
class Reduced(TreeHandler):
    
    def __init__(self, inThing=None, skip=0, max=None):
        self.minimumProportion=0.5
        self.rooted=True
        self.minNoOfTaxa = 3
        self.treatMultifurcatingRootsAsUnrooted=0
        self.minimumSplitsToBuildTree=1
        self.extendProfile=False
        self.saveTreesToFile=False
        self.saveSplitsToFile=False
        self.extendProfile=False
        self.saveExtendedTreesToFile=False
        self.saveExtendedToFile=False
        self.treeFileName='trees'
        self.splitFileName='splits.txt'
        self.extendedTreeFileName='extendedTrees'
        self.extendedFileName='extended.txt'
        self.sortByNoSplits=True
        

        TreeHandler.__init__(self, inThing, skip, max)
                
    """
    A container for 
    * reduced strict consensus 
    * reduced majority consensus, .
    * reduced strict supertree
    * reduced majority supertree
    
    Start it up like this:

        rd = Reduced(inThing)

    where inThing can be a file name

    If you are reading from a file (generally a bootstrap or mcmc
    output), you can skip some trees at the beginning, and optionally
    after that read only a maximum number of trees, like this::

        rd = Reduced('myFile', skip=1000, max=500)

    Then you can do::

        reduced strict consensus
        rd.reducedStrictConsensus()
        
        
        majority rule consensus
        rd.majorityRuleConsensus()
        
        
        write trees to file
        
        rd.writeNexus(filename)

    The input trees do not need to have the same taxa.

    These are the default settings:: 
        
      rd.minimumProportion=0.5
      rd.rooted=True
      rd.minNoOfTaxa = 3
      rd.treatMultifurcatingRootsAsUnrooted=0
      rd.minimumSplitsToBuildTree=1
      rd.extendProfile=False
      rd.saveTreesToFile=False
      rd.saveSplitsToFile=False
      rd.extendProfile=False
      rd.saveExtendedTreesToFile=False
      rd.saveExtendedToFile=False
      rd.treeFileName='trees'
      rd.splitFileName='splits.txt'
      rd.extendedTreeFileName='extendedTrees'
      rd.extendedFileName='extended.txt'
      self.sortByNoSplits=True
        
    To change one of them do this::
    
      rd = Reduced(inThing)
      rd.rooted = False
        
        
    """

    def _getTreeSplits(self):
        
        for t in self.trees:
            self._determineRooting(t)
        
        for t in self.trees:
            t.setPreAndPostOrder()
            t._makeRCSplitKeys()
            t.splits =[]
            for n in t.iterInternalsNoRoot():
                t.splits.append([n.br.rc,t.missingTaxa])
#            print t.splits

    def _determineRooting(self, theTree):
        
        gm = ['ReducedStrictConsensus._determineRooting()']
        
        nRootChildren = theTree.root.getNChildren()
        if not nRootChildren:
            gm.append("Root has no children.")
            raise Glitch, gm
        elif nRootChildren == 1:
            gm.append("Tree is rooted on a terminal node.  No workee.")
            raise Glitch, gm
        elif nRootChildren == 2:
            if self.isBiRoot == None:
                print 'Found bifurcating root'
                self.isBiRoot = 1
        else:
            if self.isBiRoot == None:
                print 'Found multifurcating root'
                self.isBiRoot = 0
            elif self.isBiRoot == 0:
                pass
            else:
                gm.append("Self.isBiRoot has been previously turned on, but now we have a non-biRoot tree.")
                raise Glitch, gm

    def writeIntersections(self, filename='rcIntersections'):
        
        if len(self.intersections) == 0:
            print 'No splits created, can not write splits to file'
            return
        
        fileHandle = open ( filename, 'w' )
        for name in self.taxNames:
            fileHandle.write(name)
            fileHandle.write(',')
        fileHandle.write('\n')
        
        list = []
        for il in self.intersections:
            list.extend(il)
        
#        print 'list1:', list
        for i in range(0,len(list)-1):
            for j in range(i, len(list)):
                if list[i].frequency < list[j].frequency:
                    temp = list[i]
                    list[i] = list[j]
                    list[j] = temp
                        
#        print 'list2:',list
        fileHandle.write('Taxnumber to taxname translation: \n')
        for i in range(0, len(self.taxNames)):
            fileHandle.write('%s,%s \n' % ( i+1, self.taxNames[i]))
#        self.printVerticalNumbers(len(self.taxNames))
        for i in list:
            fileHandle.write (i.intersectionToString(self.bitkeys) )
            fileHandle.write ('\n')
#            fileHandle.write('\n')

        fileHandle.close()
        
    def writeExtended(self, extendedFileName='rcExtended'):
        
        if len(self.extendedIntersections) > 0:
            list = []
            for il in self.extendedIntersections:
                list.extend(il)
                
#            print 'list1:', list
            for i in range(0,len(list)-1):
                for j in range(i, len(list)):
                    if list[i].frequency < list[i].frequency:
                        temp = list[i]
                        list[i] = list[j]
                        list[j] = temp
                        


            fileHandle = open ( extendedFileName, 'w')
            for name in self.taxNames:
                fileHandle.write(name)
                fileHandle.write(',')
            fileHandle.write('\n')
            for i in list:
                fileHandle.write(i.intersectionToString(self.bitkeys))
                fileHandle.write('\n')
                
            fileHandle.close()
        else:
            print 'No splits in the extended profile'
        
        
        

    


    
          
    def _rule1(self, compatibleSplits):
        highestbitkey = (self.bitkeys[-1] << 1) -1
        
        print 'HighestBitkey: %s' % (str(highestbitkey))
        
        split2Parents = {}
    
        parents = {}
        newSplits = []
        for i in range(0, len(compatibleSplits)-1):
            for j in range(i+1, len(compatibleSplits)):
#               Rule 1 
                if compatibleSplits[i][0] & compatibleSplits[j][0] != 0 and compatibleSplits[i][0] & (compatibleSplits[j][0] | compatibleSplits[j][1]) == 0 and compatibleSplits[j][0] & (compatibleSplits[i][0] | compatibleSplits[i][1]) == 0:
                    split1 = [compatibleSplits[i][0] & compatibleSplits[j][0],((compatibleSplits[j][0] | compatibleSplits[j][1]) | highestbitkey) | ((compatibleSplits[i][0] | compatibleSplits[i][1]) | highestbitkey)]
                    if split2Parent.has_key((split1[0],split1[1])):
                        split2Parent[(split1[0],split1[1])][i] = 1
                        split2Parent[(split1[0],split1[1])][j] = 1
                    else:
                        dict = {}
                        dict[i] = 1
                        split12Parent[(split1[0],split1[1])] = dict
                    newSplits.append(split1)
                    parents[(split1[0],split1[1])]
                    split2 = [compatibleSplits[i][0] | compatibleSplits[j][0],((compatibleSplits[j][0] | compatibleSplits[j][1]) | highestbitkey) & ((compatibleSplits[i][0] | compatibleSplits[i][1]) | highestbitkey)]
                    newSplits.append(split2)


    def _printMajorityRuleSplits(self, list, noSplitSets):
        for split in list:
#            for split in list:
            intersection = ''
            for bk in self.bitkeys:
                if bk & split[0]:
                    intersection = intersection + '1'
                elif bk & split[1]:
                    intersection = intersection + '?'
                else:
                    intersection = intersection + '*'
            print 'Split: %s, %s' % (intersection, str(float(split[2])/noSplitSets))
        print ' '
            
    def _findSmallestSet(self, bk, compatibleSplits):
        smallest = 0
        minPopcount = len(compatibleSplits) + 1
        for split in compatibleSplits:
            if split[0] & bk:
                popCount = self.popcount(split[0])
                if popCount < minPopcount:
                    minPopcount = popCount
                    smallest = split[0]
                elif popCount == minPopcount:
                    smallest = smallest | split[0]
        
        return smallest
    
    
    def _applySmallest(self,unknown, smallest, compatibleSplits):
#        print 'Unknown'
#        print str(unknown)
        for split in compatibleSplits:
            if split[0] & smallest != 0:
                print '11'
                self._printSplitList([split])
                split[0] = split[0] | unknown
                split[1] = split[1] ^ unknown
                self._printSplitList([split])
            if split[1] & unknown:
                print '22'
                self._printSplitList([split])
                split[1] = split[1] ^ unknown
                self._printSplitList([split])
##                print ''
                
                
    def _isOverlapping(self, split1, split2):
        if split1[0] | split2[0] == split1[0] or split1[0] | split2[0] == split2[0]:
#            print 'True 1'
            return True
        
        
        if (split1[0] | split1[1]) | (split2[0] | split2[1]) == (split1[0] | split1[1]):
#            print 'True 2'
            return True
        
        if (split1[0] | split1[1]) | (split2[0] | split2[1]) == (split2[0] | split2[1]):
#            print 'True 3'
            return True
        
        return False
        
    def _combineSet(self, set):
        newSet = set()
        used = {}
        for i in range(0, len(set)-1):
            for j in range(i, len(set)):
                if set[i][0] & set[j][0] != 0:
                    and1 = set[i][0] & set[j][1]
                    and2 = set[j][0] & set[i][1]
                    or1 = and1 | and2
                    
                    newSet.add((set[i][0] | or1, (set[i][1] &or1) ^ set[i][1]))
                    newSet.add((set[j][0] | or1, (set[j][1] &or1) ^ set[j][1]))
                    
                    used[set[i]] = 1
                    used[set[j]] = 1
                    
                else:
                    if not used.has_key(set[i]):
                        newSet.add(set[i])
                    if not used.has_key(set[j]):
                        newSet.add(set[j])
                    
                    
##                    set[i][0] = set[i][0] | or1
#                    set[i][1] = (set[i][1] &or1) ^ set[i][1]
                    
#                    set[j][0] = set[j][0] | or1
#                    set[j][1] = (set[j][1] &or1) ^ set[j][1]
#                     
#                    newSet.add((set[i][0],set[i][1]))
#                    newSet.add((set[j][0],set[j][1]))
                    
        return newSet
        
          
    def _combineSplits(self, split1, split2):
        
        if split1[0] & split2[0] != 0:
            and1 = split1[0] & split2[1]
            and2 = split2[0] & split1[1]
            or1 = and1 | and2
            
            split1[0] = split1[0] | or1
            split1[1] = (split1[1] & or1) ^ split1[1]
            split2[0] = split2[0] | or1
            split2[1] = (split2[1] & or1) ^ split2[1]        
          
    def _makeCompatible(self, split1, split2):
        if split1[0] & split2[0] == 0:
            return
        if split1[0] | split2[0] == split1[0] or split1[0] | split2[0] == split2[0]:
            return
        if (split1[0] | split1[1]) | split2[0] == (split1[0] | split1[1]):
            return
        if (split2[0] | split2[1]) | split1[0] == (split2[0] | split2[1]):
            return
        if (split1[0] | split1[1]) | (split2[0] | split2[1]) == (split1[0] | split1[1]):
            return
        if (split1[0] | split1[1]) | (split2[0] | split2[1]) == (split2[0] | split2[1]):
            return
        
        if (split1[0] | (split2[0] | split2[1]) != (split1[0] | split1[1]) and 
            split1[0] | (split2[0] | split2[1]) != (split2[0] | split2[1])):
#            print 'incompatible' 
#            self._printSplitList([split1])
#            self._printSplitList([split2])
            
            xor1 = split1[0] ^ split2[0]
            
            and1 = xor1 & split1[1]
            
            xor2 = xor1 ^ and1
            
            and2 = xor2 & split2[1]
            
            xor3 = and2 ^ xor2
            
            split1[0] = split1[0] ^ (split1[0] & xor3)
            
            split1[1] = split1[1] | xor3
            
            split2[0] = split2[0] ^ (split2[0] & xor3)
            
            split2[1] = split2[1] | xor3
            
#            print 'made compatible'
#            self._printSplitList([split1])
#            self._printSplitList([split2])
#            
            return
        
        print 'default, not good at all'
#        self._printSplitList([split1])
#        self._printSplitList([split2])

          

    
 
    
     
        


        

    
    def listOfSplits(self, treatMultifurcatingRootsAsUnrooted, minimumProportion=1.0):
        """
        listOfSplits will return a list of splits (maybe zero lenght) that fulfill the minimumProportion criterion
        """
        
        self.splits, treeNames = self.updateToCommonLeafSet(self.tfl)
        
        print 'Building intersectionlist'
        print time.ctime(time.time())
        
        intersectionList = BuildIntersections(self.bitkeys, minimumProportion, self.weights, self.splits, not treatMultifurcatingRootsAsUnrooted)
        splitData = intersectionList.buildIntersections()

        print 'Build intersectionlist'
        print time.ctime(time.time())
        
        intersections = []
        for list in splitData:
            tempList = []
            for sd in list:
                tempList.append(Intersection(sd[0],sd[1],sd[2]))
            intersections.append(tempList)
            
        for il in intersections:
            print ''
            for i in il:
                i.printIntersection(self.bitkeys)
        return intersections
        
    def reducedStrictConsensus(self):
        self.minimumProportion = 1.0
        trees, extendedTrees = self.reducedConsensusProfile()
        
        if self.saveTreesToFile:
            self.writeNexus(trees, self.treeFileName)
     
        if self.saveSplitsToFile:
            self.writeIntersections(self.splitFileName)
            
        if self.saveExtendedTreesToFile:
            self.writeNexus(extendedTrees, self.extendedTreeFileName)
            
        if self.saveExtendedToFile:
            self.writeExtended(self.extendedFileName)
            
    def majorityRuleConsensus(self):
        
        self.minimumProportion = 0.5
        trees, extendedTrees = self.reducedConsensusProfile()
        
        if self.saveTreesToFile:
            self.writeNexus(trees, self.treeFileName)
     
        if self.saveSplitsToFile:
            self.writeIntersections(self.splitFileName)
     
        if self.saveExtendedTreesToFile:
            self.writeNexus(extendedTrees, self.extendedTreeFileName)
            
        if self.saveExtendedToFile:
            self.writeExtended(self.extendedFileName)
            
    def reducedMemQuick(self):
        
        from p4.LeafSupport import TreeSubsets
        print 'Running divide and conquer to speed things up'
        sys.stdout.flush()
        ts = TreeSubsets(self.trees)
        self.trees = None
        self.buildConsensusTreeHash(ts.getConsensusTree())
        
        subtreeDicts, taxNames, taxonSets, cherries = ts.getSubTreesAndTaxonSetsFromInputTrees(verbose=True)
        
        ts = None
        if len(subtreeDicts[0].keys()) == 1: 
            print 'No subdivisions possible, solving whole tree'
        else:
            print 'Got %s subproblems to solve' % (len(cherries) + len(subtreeDicts[0].keys()))
        sys.stdout.flush()
        
        allTaxa = self.taxNames[:]
        allBitkeys = self.bitkeys[:]
        taxa2bitkey = {}
        for index in range(0,len(allTaxa)):
            taxa2bitkey[allTaxa[index]] = allBitkeys[index] 

        self.verbose = 0
        print 'Building intersectionlists, divide and conquer style'
        print time.ctime(time.time())
        sys.stdout.flush()
        
        weight = 0.0
        
        splits = []
        for i in range(0,len(taxNames)):
            
            if len(taxonSets[i]) > 2:
                trees = []
                for dict in subtreeDicts:
                    trees.append(dict[taxNames[i]])
                
                self.splits, treeNames = self.commonLeafSetTrees(trees)
                
                if not weight:
                    weight = sum(self.weights)
                    
                names = []
                bitkeys = []
                for index in range(0,len(self.taxNames)):
                    n = self.taxNames[index].split(':')
                    for t in n:
                        names.append(t)
                        bitkeys.append(self.bitkeys[index])
                        
                left = 0L
                right = 0L
                for i in range(0,len(allTaxa)):
                    if allTaxa[i] in names:
                        left = left ^ allBitkeys[i]
                    
                i = Intersection(left, right, weight)    
#                i.printIntersection(allBitkeys)
                splits.append(i)
            
                taxon2bitkey = {}
                for index2 in range(0,len(names)):
                    taxon2bitkey[names[index2]] = bitkeys[index2]
                
                intersectionList = BuildIntersections(self.bitkeys, self.minimumProportion, self.minNoOfTaxa, self.weights, self.splits, not self.treatMultifurcatingRootsAsUnrooted, self.sortByNoSplits, self.verbose)
        
                intersectionList.buildIntersections()
            
                for split in intersectionList.fullList:
                    left = 0L
                    right = 0L
                    for i in range(0,len(allTaxa)):
                        if taxon2bitkey.has_key(allTaxa[i]):
                            if taxon2bitkey[allTaxa[i]] & split.left:
                                left = left ^ allBitkeys[i]
                            elif taxon2bitkey[allTaxa[i]] & split.right:
                                right = right ^ allBitkeys[i]
                    i = Intersection(left, right, split.hits)
#                    i.printIntersection(allBitkeys)
                    splits.append(i)
                   
        for cherry in cherries:
            names = {}
            for set in cherry:
                for name in set.split(':'):
                    names[name] = 1
                    
            left = 0L
            right = 0L
            for i in range(0,len(allTaxa)):
                if names.has_key(allTaxa[i]):
                    left = left ^ allBitkeys[i]
                    
            i = Intersection(left, right, weight)    
#            i.printIntersection(allBitkeys)
            splits.append(i)

        print 'Built intersectionlists'
        print time.ctime(time.time())
        sys.stdout.flush()
        
        print 'Taxnumber to taxname translation: '
        for i in range(0, len(allTaxa)):
            print('%s,%s' % ( i+1, allTaxa[i]))
        
        intersectionList.verbose = True
        
        self.intersections = intersectionList.quickList(splits, allBitkeys)
        
        self.taxNames = allTaxa 
        self.bitkeys = allBitkeys
        
        self.extendedIntersections = []
        if self.extendProfile:
            extendedProfile = intersectionList.buildExtendedProfile(1)
            for list in extendedProfile:
                tempList = []
                for sd in list:
                    tempList.append(Intersection(sd[0],sd[1],sd[2]))
                    self.extendedIntersections.append(tempList)
        
        self.conTrees = self.buildTreesFromIntersections(self.intersections,self.minimumSplitsToBuildTree, 'Reduced strict consensus trees')
        
        if self.extendProfile and len(extendedProfile) > 0:
            self.extendedProfileTrees = self.buildTreesFromIntersections(self.extendedIntersections, self.minimumSplitsToBuildTree,  'Extended profile trees')
            return self.conTrees, self.extendedProfileTrees
            
        return self.conTrees, []        
            
    def reducedQuick(self):
        
#        list of taxanames in order from consensus tree, get from TreeSubsets
#        build taxname2bitkey with order from consensus tree
    

        
#        translate partial splits into full splits using taxonSet2taxonNames and taxname2bitkey

        from p4.LeafSupport import TreeSubsets
        print 'Running divide and conquer to speed things up'
        sys.stdout.flush()
        ts = TreeSubsets(self.trees)
        subtreeDicts, taxNames, taxonSets = ts.getSubTreesAndTaxonSetsFromInputTrees(verbose=False)
        self.buildConsensusTreeHash(ts.getConsensusTree())
        ts = None
        if len(subtreeDicts[0].keys()) == 1: 
            print 'No subdivisions possible, solving whole tree'
        else:
            print 'Got %s subproblems to solve' % (len(subtreeDicts[0].keys()))
        sys.stdout.flush()
        
        allTaxa = self.taxNames[:]
        allBitkeys = self.bitkeys[:]
        taxa2bitkey = {}
        for index in range(0,len(allTaxa)):
            taxa2bitkey[allTaxa[index]] = allBitkeys[index] 
        
#        print 'All taxa: ', allTaxa
        
#        taxname2splitList i.e. A:B:C -> [split, split,split]
#        taxname2taxonNames i.e. A:B:C -> ['A','B','C']
        taxname2splitList = {}
        taxname2taxonNames = {}
        
        self.verbose = 0
        print 'Building intersectionlists, divide and conquer style'
        print time.ctime(time.time())
        sys.stdout.flush()
        
        for i in range(0,len(taxNames)):

            trees = []
            for dict in subtreeDicts:
                trees.append(dict[taxNames[i]])
        
            self.splits, treeNames = self.commonLeafSetTrees(trees)
            taxname2taxonNames[taxNames[i]] = [self.taxNames[:], self.bitkeys[:]]
            
            if len(self.taxNames) > 2:
            
                intersectionList = BuildIntersections(self.bitkeys, self.minimumProportion, self.minNoOfTaxa, self.weights, self.splits, not self.treatMultifurcatingRootsAsUnrooted, self.sortByNoSplits, self.verbose)
        
                intersectionList.buildIntersections()
            
                taxname2splitList[taxNames[i]] = intersectionList.fullList

        
        weight = sum(self.weights)
        print 'Built intersectionlists'
        print time.ctime(time.time())
        sys.stdout.flush()
        
        splits = []
        for name in taxNames:
            temp = taxname2taxonNames[name][0]
            temp2 = taxname2taxonNames[name][1]
            names = []
            bitkeys = []
            for index in range(0,len(temp)):
                n = temp[index].split(':')
                for t in n:
                    names.append(t)
                    bitkeys.append(temp2[index])
#            print 'names:',names
            
            taxon2bitkey = {}
            for index2 in range(0,len(names)):
                taxon2bitkey[names[index2]] = bitkeys[index2]
            
            left = 0L
            right = 0L
            for i in range(0,len(allTaxa)):
                if taxon2bitkey.has_key(allTaxa[i]):
                    left = left ^ allBitkeys[i]

            i = Intersection(left, right, weight)
#            i.printIntersection(allBitkeys)          
            splits.append(i)  
                
            if index <= 1:
                pass
            else:
                for split in taxname2splitList[name]:
#                    split.printIntersection(bitkeys)
                    left = 0L
                    right = 0L
                    for i in range(0,len(allTaxa)):
#                        print allTaxa[i]
                        if taxon2bitkey.has_key(allTaxa[i]):
#                            print '1'
                            if taxon2bitkey[allTaxa[i]] & split.left:
#                                print '2'
                                left = left ^ allBitkeys[i]
                            elif taxon2bitkey[allTaxa[i]] & split.right:
#                                print '3'
                                right = right ^ allBitkeys[i]
#                        else:
#                            left = left ^ allBitkeys[i]
#                        right = right ^ allBitkeys[i]
#                    print 'Freq:', split.hits
                    i = Intersection(left, right, split.hits)
#                    i.printIntersection(allBitkeys)          
                    splits.append(i)
        
        intersectionList.verbose = True
        
        self.intersections = intersectionList.quickList(splits, allBitkeys)
        
        self.taxNames = allTaxa 
        self.bitkeys = allBitkeys
        
        
        self.extendedIntersections = []
        if self.extendProfile:
            extendedProfile = intersectionList.buildExtendedProfile(1)
            for list in extendedProfile:
                tempList = []
                for sd in list:
                    tempList.append(Intersection(sd[0],sd[1],sd[2]))
                    self.extendedIntersections.append(tempList)
        
        self.conTrees = self.buildTreesFromIntersections(self.intersections,self.minimumSplitsToBuildTree, 'Reduced strict consensus trees')
        
        if self.extendProfile and len(extendedProfile) > 0:
            self.extendedProfileTrees = self.buildTreesFromIntersections(self.extendedIntersections, self.minimumSplitsToBuildTree,  'Extended profile trees')
            return self.conTrees, self.extendedProfileTrees
            
        return self.conTrees, []
            
    def reducedConsensusProfile(self):
        """
        reducedConsensusProfile will depending on the minimumProportion and the input data produce 4 different outputs.
        * Given a set of trees with identical leafsets and a minimumProportion=1.0 a reduced strict consensus tree 
        will be produced, i.e. a tree where all splits can be found in all input trees
        * Given a set of trees with identical leafsets and a 1.0 < minimumProportion >= 0.5 a reduced majority profile
        will be produced, i.e. a set of trees defined on differing leafsets with differing support, this profile will
        include the reduced strict consensus tree if it exists
        * Given a set of trees with differing leafsets and a minimumProportion= 1.0 a reduced strict super tree will
        be produced, i.e. a tree where all splits can be found in all input trees
        * Given a set of trees with differing leafsets and a 1.0 < minimimProportion >= 0.5 a reduced majority supertree
        profile will be produced, i.e. a set of trees defined on differing leafsets with differing support, this profile will
        include the reduced strict supertree tree if it exists
        """
        
        self.splits, treeNames = self.commonLeafSetTrees(self.trees)
        
#        self.splits, treeNames = self.updateToCommonLeafSet(self.tfl)

        print 'Building intersectionlist'
        print time.ctime(time.time())
        sys.stdout.flush()
        intersectionList = BuildIntersections(self.bitkeys, self.minimumProportion, self.minNoOfTaxa, self.weights, self.splits, not self.treatMultifurcatingRootsAsUnrooted, self.sortByNoSplits)
        
        self.intersections = intersectionList.buildIntersections()

        print 'Done'
        print time.ctime(time.time())
        sys.stdout.flush()
        
        self.extendedIntersections = []
        if self.extendProfile:
            extendedProfile = intersectionList.buildExtendedProfile(1)
            for list in extendedProfile:
                tempList = []
                for sd in list:
                    tempList.append(Intersection(sd[0],sd[1],sd[2]))
                    self.extendedIntersections.append(tempList)
         
        self.conTrees = self.buildTreesFromIntersections(self.intersections,self.minimumSplitsToBuildTree, 'Reduced strict consensus trees')
        
        if self.extendProfile and len(extendedProfile) > 0:
            self.extendedProfileTrees = self.buildTreesFromIntersections(self.extendedIntersections, self.minimumSplitsToBuildTree,  'Extended profile trees')
            return self.conTrees, self.extendedProfileTrees
        
        
        return self.conTrees, []
        

        
class TreeBuilderFromSplits(object):
    
    def __init__(self, bitKeys, taxnames, internalNames=True):
        self.taxNames = taxnames
        self.bitkeys = bitKeys
        self.internalNames = internalNames
    
    def buildTreeFromInformativeList(self, list, treeName='conTreeName', excluded=0L):
        informativeBits = ((self.bitkeys[-1] << 1) -1) ^ list[0].right
##        print 'informative: ',informativeBits
#        print 'Excluded: ', excluded
        informativeBits = informativeBits ^ excluded
        
#        print 'informative: '
#        Intersection(informativeBits).printIntersection(self.bitkeys)
        
        conTree = Tree()
        conTree.name = treeName

        conTree.root = Node()
        conTree.root.nodeNum = 0
        conTree.root.isLeaf = 0
        conTree.nodes.append(conTree.root)

        index = 0
        leafNo = 1
        first = True
        rootBitkey = 0L
        n = Node()
        previous = n
        for bk in self.bitkeys:
            if bk & informativeBits:
                if first:
                    n.nodeNum = leafNo
                    conTree.nodes.append(n)
                    conTree.root.leftChild = n
                    n.parent = conTree.root
                    n.isLeaf = 1
                    n.name = self.taxNames[index]
                    n.br.splitKey = bk
                    rootBitkey = rootBitkey | bk
                    first = False
                else: 
                    n = Node()
                    n.nodeNum = leafNo
                    conTree.nodes.append(n)
                    previous.sibling = n
                    n.parent = conTree.root
                    n.isLeaf = 1
                    n.name = self.taxNames[index]
                    n.br.splitKey = bk
                    rootBitkey = rootBitkey | bk
                    previous = n
                leafNo += 1                   
            index += 1
        conTree.root.br.splitKey = rootBitkey  
        for s in list:
            self.applySplitToTree(s, conTree)
        
        conTree.setPreAndPostOrder()
        conTree._setTaxNamesFromLeaves()
        
        return conTree
    
    def popcount(self, n):
        count = 0
        for bk in self.bitkeys:
            if n < bk:
                return count
            if bk & n:
                count += 1
        return count
    
    def applySplitToTree(self, split, tree):
        
        insertNode = self.findInsertNode(split, tree)
        
        self.insertSplit(split, insertNode, tree)
        
        return 1 
            
    def findInsertNode(self, split, tree):
        
        node = tree.root
        andScore = 0
        nodeSplitScore = 9999999999
           
        for n in tree.nodes:
            if not n.isLeaf:
                tScore = self.popcount(n.br.splitKey & split.left)
                if tScore > andScore:
                    andScore = tScore
                    node = n
                    nodeSplitScore = self.popcount(node.br.splitKey)
                elif tScore == andScore and self.popcount(n.br.splitKey) < nodeSplitScore:
                    nodeSplitScore = self.popcount(n.br.splitKey)
                    node = n
        
        return node

    def insertSplit(self, split, node, tree):   
        
        n = Node()
        n.isLeaf = 0
        n.parent = node
        n.br.splitKey = split.left
        n.nodeNum = len(tree.nodes)
        tree.nodes.append(n)
        
   
        n.sibling = node.leftChild
        node.leftChild = n
        if split.name == None:
            if self.internalNames:
                n.name = '%.0f' % (100. * split.frequency)
        else:
            n.name = split.name
              
        addNode = self.getSiblingInSplit(node, split)
 
        while addNode:
#            print 'while addNode'
            self.addNodeToNode(n, addNode)
            addNode = self.getSiblingInSplit(node, split)
    
    def getSiblingInSplit(self, node, split):

        previous = node.leftChild
        n = previous.sibling
        
        if not n:
#            print 'no n'
            return None
        while n.sibling:
            if self.isNodeInSplit(n, split):
#                print 'found node in split'
                previous.sibling = n.sibling
                n.sibling = None
                return n
            if not n.sibling:
#                print 'reached end of siblings'
                return None
            previous = n
            n = n.sibling
        if self.isNodeInSplit(n, split):
#            print 'found node in split outside loop'
            previous.sibling = None
#           n.wipe()
            n.sibling = None
#            n.parent = None
            return n
#        print 'reached end of loop'
        return None
       
    def addNodeToNode(self, node, addNode):
        
        addNode.parent = node
        
        if node.getNChildren() == 0:
#            print 'node.getNChildren == 0'
            node.leftChild = addNode
        elif node.getNChildren() == 1:
#            print 'node.getNChildren == 1'
            leftChild = node.leftChild
            leftChild.sibling = addNode
        else:
#            print 'node.getNChildren == @@@'
            rightChild = node.rightmostChild()
            rightChild.sibling = addNode
        
    def isNodeInSplit(self, node, split):
        if node.br.splitKey | split.left == split.left:
            return True
        else:
            return False

    
class Intersection(object):
    
    def __init__(self, bitKey, bright=0L, bhits=0, bfrequency=0):
#        print 'Frequency: ', bhits
        self.left = bitKey
        self.right = bright
        self.hits = bhits
        self.frequency = bhits
        self.name = None
        self.excluded = []
        
    def __cmp__(self, other):
        return cmp(other.frequency, self.frequency)
        
    def __str__(self):
        intersection = ''
        for bk in bitkeys:
            if bk & self.left:
                intersection = intersection + '*'
            elif bk & self.right:
                intersection = intersection + '?'
            else:
                intersection = intersection + '.'
        return intersection +', ' +str(self.frequency*100)
    
    def copy(self):
        return Intersection(self.left, self.right, self.frequency)

    def firstHit(self, bitkeys):
        hit = 0
        for bk in bitkeys:
            if bk & self.left:
                return hit
            hit += 1
        return hit
    
    def setExcluded(self, bitkeys):
        self.excluded = []
        index = 0
        for bk in bitkeys:
            if self.right & bk:
                self.excluded.append(index)
            index += 1
    
    def popcount(self, bitkeys):
        count = 0
        for bk in bitkeys:
            if self.left < bk:
                return count
            if bk & self.left:
                count += 1
        return count     
    
    def getInformativeBits(self, bitkeys):
        highestNo = (bitkeys[-1] << 1) -1
        return highestNo ^ self.right
     
    def popcountLR(self, bitkeys):
        countie = self.left | self.right
        count = 0
        for bk in bitkeys:
            if countie < bk:
                return count
            if bk & countie:
                count += 1
        return count 
    
    def popcountExcluded(self, bitkeys):
        index = 0
        for bk in bitkeys:
            if self.right < bk:
                return index
            if self.right & bk:
                index += 1
        return index
    
    def isIntersectionTrivial(self, bitkeyList):
        if self.left == 0:
            return 1
        if self.popcount(bitkeyList) <= 1 or self.popcountLR(bitkeyList) == len(bitkeyList):
            return 1
        if self.left | self.right == 0:
            return 1
        return 0
        
    def dump(self):
        print "%s %s, %6s %s, %6s %s, %6s %s  " % ('L:', self.left, 'R:', self.right, 'F:', self.frequency, 'H:', self.hits)

    def printIntersection(self, bitkeys, weight=1.0):
        intersection = ''
        for bk in bitkeys:
            if bk & self.left:
                intersection = intersection + '*'
            elif bk & self.right:
                intersection = intersection + '?'
            else:
                intersection = intersection + '.'
        if self.name != None:
            print '%s  %s' % (intersection, self.name)
        else:
            f = self.frequency/weight
            if f <= 1:
                f = f *100
            f = int(f)
            print '%s  %s' % (intersection, f)

    def intersectionToString(self, bitkeys):
        intersection = ''
        for bk in bitkeys:
            if bk & self.left:
                intersection = intersection + '*'
            elif bk & self.right:
                intersection = intersection + '?'
            else:
                intersection = intersection + '.'
        return intersection +', ' +str(self.frequency*100)

    def addSplit(self, bitkey):
        newSplit = Intersection(0L)
        newSplit.left = self.left & bitkey
        newSplit.right = (self.left ^ bitkey) | self.right
#        print 'Creating new intersection: '
#        newSplit.dump()
        return newSplit 
    
    def addReverseSplit(self, bitkey, highestBitkey):    
        newSplit = Intersection(0L)
        newSplit.left = (self.left ^ bitkey) & self.left
        newSplit.right = (self.left ^ (~bitkey & ((highestBitkey << 1) -1))) | self.right
        return newSplit    
    
    def intersection(self, other):
        newSplit = Intersection(0L)
        newSplit.left = self.left & other.left
        newSplit.right = (self.left ^ other.left) | self.right | other.right
        return newSplit
    
    def reverseIntersection(self, other, highestBitkey):
        newSplit = Intersection(0L)
        newSplit.left = (self.left ^ other.left) & self.left
        newSplit.right = (self.left ^ (~other.left & ((highestBitkey << 1) -1))) | self.right | other.right
        return newSplit
    
    def matchingSplit(self, split):
        if self.left == split:
#            if self.left == 768:
#                print '%s i.l:%s i.r:%s == s:%s' % ('Found matching split', self.left, self.right, split)
            return 1
        if ((self.left ^ split) | self.right) == self.right:
#            if self.left == 768:
#                print '%s i.left:%s, i.right:%s, s:%s' % ('Found matching split XOR', self.left, self.right, split)
            return 1
#        if self.left == 768:
#            print '%s i.l:%s  i.r%s != s:%s' % ('Found NONmatching split', self.left, self.right, split)
        return 0
    
    def overlapping(self, other):
        if self.right & other.right != 0:
            return True
        return False
    
    def smaller(self, other):
        if (self.right & other.right) == self.right:
#            print 'True'
            return True
        return False
    
    def definesNonOverlappingSets(self, other):       
        if self.right == 0 or other.right == 0:
            return False
        if self.right & other.right == 0:
            return True
        return False
    
    def definesSmallerAndOverlappingSet(self, other, minimumProportion):     

        if self.right & other.right == self.right:
            return True
        return False  
    
    def isMoreInformative(self, list, minimumProportion):
        for other in list:
#            self.printIntersection()
#            other.sprintIntersection()
            if self.frequency > other.frequency and other.frequency < minimumProportion:
#                print '1'
                return 1
            if other.frequency > self.frequency and self.frequency < minimumProportion:
#                print '0'
                return 0
            selfORother = self.right | other.right 
            if selfORother == self.right:
#                print '0'
                return 0
            if selfORother == other.right:
#                print '1'
                return 1
        return 2

    def isSameSet(self, other):
        
#        if self.left | self.right == other.left | other.right:
#            return True
        
        if self.right == other.right:
            return True

    def isSameLooseSet(self, other):
        if self.left | self.right == other.left | other.right:
            return True
        return False

    def isOverlapping(self, other):
        
        if self.right == 0 or other.right == 0:
            return True
        
        if self.right == other.right:
            return True
        
#        if self.left & other.right != 0 or other.left & self.right != 0:
#            return 0
#        if self.left & other.left != 0:
#            print 'Is overlapping: %s:%s ' % (self.left, other.left)
#            return 1
#        print 'Is not overlapping: %s:%s ' % (self.left, other.left)
        
        
#        if self.left & other.left == 0 and self.right ^ other.right == 0:
#            return 1
        return 0
    
    def isOtherRedundant(self, other):
        if self.left | self.right == other.left | other.right:
#            if self.left | other.left == self.left:
            if self.right & other.right == self.right:
#            if self.left | self.right == other.left | other.right:
#                print 'Redundant o.l: %s, o.r: %s, i.l: %s' % (other.left, other.right, self.left)
                return 1
#        print 'NonRedundant o.l: %s, o.r: %s, i.l: %s' % (other.left, other.right, self.left)
        return 0 
    
    
    def isCorrectIntersection(self):
        if self.left & self.right != 0:
            return 0
        return 1
    
    def isIdentical(self, other):
        if self.left == other.left and self.right == other.right:
            return 1
        return 0

 
