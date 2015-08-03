read('../A_NavidiAlignFromYang/navidiSSRNA.nex')
a = var.alignments[0]

print '%20s   %7s%7s%7s%7s%7s' % ('species', 'T', 'C', 'A', 'G', 'G + C')
for i in range(len(a.sequences)):
    s = a.sequences[i]
    c = a.composition([i]) # comp of an individual sequence
    print '%20s   ' % s.name,
    for b in [3, 1, 0, 2]:   # p4 order is acgt, yang and roberts order is tcag
        print '%3.4f' % c[b],
    print '%3.4f' % (c[1] + c[2]),
    print ''

print '%20s   ' % 'average',
c = a.composition()  # overall comp
for b in [3, 1, 0, 2]:
    print '%3.4f' % c[b],
print '%3.4f' % (c[1] + c[2]),
print ''
