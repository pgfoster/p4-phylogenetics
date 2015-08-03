var.verboseRead = 0
var.warnReadNoFile = 0

if 1:
    read('((A, B), C, (D, E));')
    read('((A, B), D, (E, C));')
    read('((C, A), (D, B), E);')
    taxNames = list(string.uppercase[:5])
    theSplitTax = ['A', 'C']
else:
    taxNames = list(string.uppercase[:7])
    for i in range(30):
        var.trees.append(func.randomTree(taxNames))
    theSplitTax = ['A', 'C', 'F']

tt = Trees(taxNames=taxNames)
treesWithSplit = tt.getTreesWithSplit(theSplitTax)

print 'These trees have the split containing %s' % theSplitTax
for t in treesWithSplit:
    t.draw()


