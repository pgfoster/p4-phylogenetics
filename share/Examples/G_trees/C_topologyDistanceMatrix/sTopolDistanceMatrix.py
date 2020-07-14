taxNames = list(string.ascii_uppercase[:7])
for i in range(6):
    t = func.randomTree(taxNames)
    t.name = 't%i' % (i + 1)
    var.trees.append(t)

tt = Trees(taxNames=taxNames)
dm = tt.topologyDistanceMatrix('wrf')
dm.writeNexus()

try:
    import p4.scqdist
    print("\nQuartet distances from scqdist")
    dm = tt.topologyDistanceMatrix('scqdist')
    dm.writePhylip(digitsAfterDecimal=0)
except ImportError:
    print("\nUnable to import the scqdist module from p4")

try:
    import p4.pytqdist
    print("\nQuartet distances from tqdist")
    dm = tt.topologyDistanceMatrix('tqdist')
    dm.writePhylip(digitsAfterDecimal=0)
except ImportError:
    print("\nUnable to import the pytqdist module from p4")
