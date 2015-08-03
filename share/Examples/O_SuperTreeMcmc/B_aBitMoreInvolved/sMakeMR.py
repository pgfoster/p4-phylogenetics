from p4.mrp import mrp
read('inTrees.phy')
a = mrp(var.trees)
a.writeNexus('mr.nex')
