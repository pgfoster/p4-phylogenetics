from p4.LeafSupport import CherryRemover

#These three lines removes the cherries and saves the resulting trees in a file. 
rc = CherryRemover('../input.nex')
rc.removeCherries()
rc.saveTrees('input.cherries.removed.nex')

#The file is then used to calculate the leaf stabilities
ls = LeafSupport('input.cherries.removed.nex')

ls.removeCherries=True

ls.leafSupport()



