tp = TreePartitions('tt.nex')
tp.writeSplits()  # If you like this sort of thing ...
t = tp.consensus()
t.draw()
