import sys
import os
sys.path.append(os.path.abspath("tqDist-1.0.2/tqDist"))
import tqdist

# spaces in the string are not allowed.
# taStr = "(A, B, (C, (D, E)));"
taStr = "(A:0.1,B:0.3,(C:0.1,(D:0.1,E:0.1)thisName:0.1,F:0.1):0.1)theRoot;"

tbStr = "(B:0.1,A:0.3,(C:0.1,(D:0.1,E:0.1)thisName:0.1,F:0.1):0.1)theRoot;"

ret = tqdist.qdist(taStr, tbStr)

print("Got qdist = %s" % ret)

