read('d.nex')
d = Data()
for i in range(3):
    read('t%i.p4_tPickle' % (i + 1))

tt = Trees()
tt.data = d
if 1:
    tt.consel()
else:
    tt.rell()
