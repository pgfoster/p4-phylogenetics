read('t.nex')
t1 = var.trees[0]
t2 = var.trees[1]

# See page 532 in Felsenstein

print "The 'symmetric distance' = ", t1.topologyDistance(t2, metric='sd')
print "The 'weighted Robinson Foulds distance'  = ", t1.topologyDistance(t2, metric='wrf')
ret = t1.topologyDistance(t2, metric='bld')
print "The 'branch score' = %.4f" % (ret * ret)
print "The 'branch-length distance' = %.4f" % ret
