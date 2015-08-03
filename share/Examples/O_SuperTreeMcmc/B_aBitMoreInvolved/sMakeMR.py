from p4.MRP import mrp
read('inTrees.phy')
a = mrp(var.trees)
a.writeNexus('mr.nex')
