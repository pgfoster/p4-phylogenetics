read('d.nex')
a=var.alignments[0]
a.writePhylip()
m = func.maskFromNexusCharacterList("1 3-5 7 11", a.length, invert=0)
print(m)
b = a.subsetUsingMask(m)
b.writePhylip()
