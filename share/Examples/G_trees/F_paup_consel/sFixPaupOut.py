import sys
nTrees = 15
nSites = 400

fIn = file('siteLikes')
fOut = file('siteLikes.txt', 'w')
fIn.readline() # skip the first line
fOut.write('Tree\t-lnL\tSite\t-lnL\n')

for i in range(nTrees):
    for j in range(nSites):
        aLine = fIn.readline()
        if not aLine:
            print('no workee!  Ran out of lines.')
            sys.exit()
        splitLine = aLine.split()
        if len(splitLine) != 2:
            print('no workee!  Bad line %s' % aLine)
            sys.exit()
        fOut.write('\t\t%s\t%s\n' % (splitLine[0], splitLine[1]))
    aLine = fIn.readline()
    fOut.write(aLine)

fIn.close()
fOut.close()


