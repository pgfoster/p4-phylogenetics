read('../noTRuberNoGapsNoAmbiguities.nex')
a=var.alignments[0]
tList = []
print("Doing bootstrap ...")
for i in range(200):
    b = a.bootstrap()
    dm = b.compositionEuclideanDistanceMatrix()
    t = dm.njUsingPaup()
    tList.append(t)

tt = Trees(tList, taxNames=a.taxNames)
tp = TreePartitions(tt)
t = tp.consensus()
for n in t.iterInternalsNoRoot():
    n.name = "%.0f" % (100. * n.br.support)
t.draw()
