# Use the mb cons tree as the master, and add supports from the paup boot tree.

# We need an ordered list of taxNames.
read('dB.nex')
a = var.alignments[0]
# a.taxNames is the list we want.

# Read in the two trees.  The mrbayes con file has 2 trees; we only
# want the first one.
read('mbout.con')
tMB = var.trees[0] # name the first
var.trees.pop()    # discard the second

read('paupBootTree.nex')
tPAUP = var.trees[1] # name it

tMB.taxNames = a.taxNames
tPAUP.taxNames = a.taxNames

# Split keys are numerical versions of the 'dot-star' split notation.
# The same split on the two trees would have the same split key.
tMB.makeSplitKeys()
tPAUP.makeSplitKeys()

# Make a dictionary, so that we can fish out nodes in the paup tree
# given a split key.  Split keys are found on node branches, here
# n.br.
myDict = {}
for n in tPAUP.iterInternalsNoRoot():
    myDict[n.br.splitKey] = n

for nMB in tMB.iterInternalsNoRoot():
    # Given a split key in the mrbayes tree, we can find the
    # corresponding node in the paup tree, using the split key with
    # the dictionary.
    nPAUP = myDict.get(nMB.br.splitKey)
    # If there was none, then nPAUP is None
    if nPAUP:
        nMB.name = '%s/%s' % (nMB.name, nPAUP.name)
    else:
        nMB.name = '%s/-' % nMB.name
    #print nMB.name
tMB.writeNexus('combinedSupportsTree.nex')
