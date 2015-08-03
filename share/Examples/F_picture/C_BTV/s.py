read('t.1111.nex')
read('t555.nex')
read('t43.nex')
read('t39.nex')
read('sets1.nex')

if 1:
    t = var.trees[0]
    t.setNexusSets()
    t.btv()

    t = var.trees[1]
    t.btv()

# Out with the old nexusSets, in with the new.
var.nexusSets = None
read('sets2.nex')

t = var.trees[2]
t.setNexusSets()
t.tv()

t = var.trees[3]
t.setNexusSets()
t.tv()



