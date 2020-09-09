import sys
import os
sys.path.append(os.path.abspath("tqDist-1.0.2/tqDist"))
import pytqdist

if 1:
    # This does not need p4
    # spaces in the string are not allowed.
    # eg taStr = "(A, B, (C, (D, E)));"
    taStr = "(A:0.1,B:0.3,(C:0.1,(D:0.1,E:0.1)thisName:0.1,F:0.1):0.1)theRoot;"
    tbStr = "(B:0.2,A:0.55,(D:0.1,(C:0.1,E:0.1)thisName:0.1,F:0.1):0.1)theRoot;"

    ret = pytqdist.qdistFromStrings(taStr, tbStr)
    print("Got qdist = %s" % ret)

if 1:
    # This needs p4
    nTx = 311
    tA = func.randomTree(nTax=nTx)
    toCollapse = []
    for n in tA.iterInternalsNoRoot():
        if random.random() < 0.2:
            toCollapse.append(n)
    for n in toCollapse:
        tA.collapseNode(n)

    tB = func.randomTree(nTax=nTx)
    toCollapse = []
    for n in tB.iterInternalsNoRoot():
        if random.random() < 0.2:
            toCollapse.append(n)
    for n in toCollapse:
        tB.collapseNode(n)
    tAstr = tA.writeNewick(toString=True, spaceAfterComma=False)
    tBstr = tB.writeNewick(toString=True, spaceAfterComma=False)

    ret = pytqdist.qdistFromStringsVerbose(tAstr, tBstr)
    print("Got qdist = %s" % ret)
