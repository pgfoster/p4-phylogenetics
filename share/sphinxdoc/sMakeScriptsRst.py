# Make the rst file scripts.rst
from p4.P4Error import P4Error

start = """
========================
Scripts for common tasks
========================

This is a cook book of reminders, suggestions, and points of departure
for customization.  You can copy and paste these into your own files.
Or you can use the p4 recipes function, ``func.recipes()``, to write new
files.

In these scripts I suggest variations, often by enclosing the variation
in an ``if 1:`` block to turn it on, or an ``if 0:`` block to turn it off.
Hopefully the meaning is clear and the (sometimes mutually exclusive)
possibilities are obvious.

"""

gm = ["func.recipes()"]
if not var.examplesDir:
    gm.append("Can't find the Examples directory.")
    raise Glitch, gm
recipesDir = os.path.join(var.examplesDir, 'W_recipes')

myRelPath = os.path.relpath(recipesDir)

fList = glob.glob("%s/*.py" % myRelPath)
#print fList
bNames = [os.path.basename(nm) for nm in fList]
#print bNames
firstLines = []
for fN in fList:
    f = file(fN)
    firstLine = f.readline()
    f.close()
    if firstLine.startswith('#'):
        firstLine = firstLine[1:]
        firstLine = firstLine.strip()
    else:
        firstLine = "Fix me!"
    firstLines.append(firstLine)
recipesList = []
for fNum in range(len(fList)):
    recipesList.append([firstLines[fNum], bNames[fNum]])
#print recipesList
    
#.. literalinclude:: pyplots/ellipses.py

start = start[1:]
print start
print
#print ".. literalinclude:: %s/sCon.py" % pre

for rNum in range(len(recipesList)):
    rec = recipesList[rNum]
    heading = rec[0]
    if heading.endswith('.'):
        heading = heading[:-1]
    uLine = '-' * len(heading)
    print heading
    print uLine
    print
    print ".. literalinclude:: %s" % fList[rNum]
    print
    
