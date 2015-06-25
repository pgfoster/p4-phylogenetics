# Matrix representation / parsimony.

from Tree import Tree
from Alignment import Alignment
from Glitch import Glitch
import func
from NexusSets import CharSet
from Node import Node
from TreePartitions import TreePartitions



def mrpSlice(self, pos, zeroBasedNumbering=True):
    """Pretty-print a mrp site, with no '?' positions.

    Zero-based numbering, unless arg zeroBasedNumbering is set to False.
    """
    if zeroBasedNumbering:
        ss = self.sequenceSlice(pos)
        print "mrp matrix position (zero-based) %i" % pos
    else:
        pos2 = pos - 1
        ss = self.sequenceSlice(pos2)
        print "mrp matrix position (1-based) %i" % pos2
    for txNum in range(self.nTax):
        if ss[txNum] == '?':
            pass
        else:
            print "  %25s %s" % (self.taxNames[txNum], ss[txNum])

Alignment.mrpSlice = mrpSlice
del(mrpSlice)

def mrp(trees, taxNames=None):
    """Code a list of trees with matrix representation.

    The input should be a list of p4 Tree objects.

    The argument 'taxNames' need not be supplied, but you can if you
    want to.

    This returns an alignment, with a character set for each input tree.

    For example, you might say::
      
      read('myTrees.phy')
      a = mrp(var.trees)
      a.writeNexus('a.nex')
    """

    gm = ['mrp()']
    if type(trees) != type([]):
        gm.append("The 'trees' arg should be a list of p4 tree objects.")
        raise Glitch, gm
    for t in trees:
        if not isinstance(t, Tree):
            gm.append("The 'trees' arg should be a list of p4 tree objects.")
            raise Glitch, gm
    
    myTaxNames = []
    for t in trees:
        for n in t.iterLeavesNoRoot():
            if n.name not in myTaxNames:
                myTaxNames.append(n.name)
    if taxNames:
        suppliedTaxNamesSet = set(taxNames)
        myTaxNamesSet = set(myTaxNames)
        if suppliedTaxNamesSet != myTaxNamesSet:
            print suppliedTaxNamesSet
            print myTaxNamesSet
            symDiff = myTaxNamesSet.symmetric_difference(suppliedTaxNamesSet)
            gm.append("The taxNames list supplied does not represent the taxa in the input trees.")
            gm.append("The symmetric difference is:")
            gm.append(symDiff)
            raise Glitch, gm
    else:
        taxNames = myTaxNames

    # make bitKey's for taxNames, in a dictionary
    txBkDict = {}
    for tNum in range(len(taxNames)):
        tx = taxNames[tNum]
        bk = 1L << tNum
        txBkDict[tx] = bk

    # Decorate trees with BitKeys, and count the number of splits.
    nSplits = 0
    for t in trees:
        tNSplits = 0
        for n in t.iterPostOrder():
            if not n == t.root:
                if n.isLeaf:
                    # order comes from taxNames, not from the tree
                    n.br.bitKey = 1L << taxNames.index(n.name)
                else:
                    nSplits += 1
                    tNSplits += 1
                    childrenNums = t.getChildrenNums(n)
                    try:
                        x = t.nodes[childrenNums[0]].br.bitKey
                        for i in childrenNums[1:]:
                            y = t.nodes[i].br.bitKey
                            x = x | y
                    except AttributeError:
                        print "t.preAndPostOrderAreValid = %s" % t.preAndPostOrderAreValid
                        #t.draw()
                        print "n is nodeNum %i" % n.nodeNum
                        print "childrenNums = %s" % childrenNums
                        raise AttributeError
                    n.br.bitKey = x

        t.nSplits = tNSplits
        n = t.root
        assert not n.isLeaf
        childrenNums = t.getChildrenNums(n)
        x = t.nodes[childrenNums[0]].br.bitKey
        for i in childrenNums[1:]:
            y = t.nodes[i].br.bitKey
            x = x | y
        t.taxBits = x

    if nSplits == 0:
        for t in trees:
            t.write()
        gm = ["mrp().  No splits were found in the input trees."]
        gm.append("That does not work.")
        raise Glitch, gm
    a = func.newEmptyAlignment(dataType='standard', symbols='01', taxNames=taxNames, length=nSplits)
    a.setNexusSets()
    for s in a.sequences:
        s.sequence = list(s.sequence)
    siteNum = 0
    tRange = range(len(taxNames))
    for tNum in range(len(trees)):
        t = trees[tNum]
        if t.nSplits:
            csName = 'cs%i' % tNum
            cs = CharSet(a.nexusSets)
            cs.nChar = nSplits
            cs.name = csName
            cs.num = tNum
            cs.lowName = csName
            cs.format = 'vector'
            cs.start = siteNum
            for n in t.iterPostOrder():
                if n != t.root:
                    if not n.isLeaf:
                        assert n.br.bitKey
                        for tNum in tRange:
                            tx = taxNames[tNum]
                            bk = txBkDict[tx]
                            s = a.sequences[tNum].sequence
                            if bk & n.br.bitKey:
                                s[siteNum] = '1'
                            elif bk & t.taxBits:
                                s[siteNum] = '0'
                            else:
                                s[siteNum] = '?'
                        siteNum += 1
            cs.mask = ['0'] * nSplits
            for cPos in range(cs.start, siteNum):
                cs.mask[cPos] = '1'
            cs.mask = ''.join(cs.mask)
            cs.standardize()
            a.nexusSets.charSets.append(cs)
            a.nexusSets.charSetsDict[cs.name] = cs
    for s in a.sequences:
        s.sequence = ''.join(s.sequence)
    return a

                    
def reverseMrp(alignment):
    """Reconstruct trees from a matrix representation.

    This needs character sets, one for each tree.

    You might say::
    
      read('a.nex')          # read the matrix representation in
      a = var.alignments[0]  # give the alignment a name
      a.setNexusSets()       # apply var.nexusSets to the alignment
      tt = reverseMrp(a)     # the function returns a list of tree objects
      for t in tt:
          t.write()

    """

    a = alignment
    assert a.nexusSets
    assert a.nexusSets.charSets
    tRange = range(len(a.taxNames))
    tt = []
    for cs in a.nexusSets.charSets:
        #print cs.name
        #cs.dump()
        #cs.vectorize()
        vPos = 0
        while cs.mask[vPos] == '0':
            vPos += 1
        firstVPos = vPos
        firstSite = a.sequenceSlice(vPos)
        #print firstSite
        thisTNames = []
        tNums = []
        for tPos in tRange:
            if firstSite[tPos] != '?':
                thisTNames.append(a.taxNames[tPos])
                tNums.append(tPos)
        #print thisTNames

        csMaskChar = cs.mask[vPos]
        while 1:
            if csMaskChar == '1':
                st = a.sequenceSlice(vPos)
                for tPos in tRange:
                    if st[tPos] == '?':
                        assert tPos not in tNums, "bad site %i, taxon %s" % (vPos, a.taxNames[tPos])
                    else:
                        assert tPos in tNums, "bad site %i, taxon %s" % (vPos, a.taxNames[tPos])
            vPos += 1
            try:
                csMaskChar = cs.mask[vPos]
            except IndexError:
                break

        vPos = firstVPos
        csMaskChar = cs.mask[vPos]
        partitions = []
        while 1:
            if csMaskChar == '1':
                st = a.sequenceSlice(vPos)
                txPos = 0
                aPart = []
                for tPos in tRange:
                    if st[tPos] != '?':
                        if st[tPos] == '1':
                            aPart.append(txPos)
                        txPos += 1
                partitions.append(aPart)
            vPos += 1
            try:
                csMaskChar = cs.mask[vPos]
            except IndexError:
                break
        #print partitions
        tp = TreePartitions()
        t = tp.makeTreeFromPartitions(partitions, taxNames=thisTNames)
        t.name = cs.name 
        tt.append(t)
    return tt
        

