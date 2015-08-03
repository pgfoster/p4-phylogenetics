read("d.nex")
a = var.alignments[0]
dm = a.logDet(correction='L94', doPInvarOfConstants=True,
              pInvar=None, pInvarOfConstants=None,
              missingCharacterStrategy='fudge', minCompCount=1,
              nonPositiveDetStrategy='invert')
dm.writeNexus("dm.nex")
t = dm.bionj()
#t.taxNames = a.taxNames
t.writeNexus('tLogDet.nex')

