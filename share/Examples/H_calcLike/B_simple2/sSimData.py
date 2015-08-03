def pepper(self, proportion=0.1, andAmbiguities=1):
    """Pepper the alignment (self) with random gaps and ambiguities."""
    import random
    import math
    ambigs = 'rymkswhbvdnnn'
    #ambigs = 'n'
    for s in self.sequences:
        s.sequence = list(s.sequence)
        for i in range(self.length):
            if andAmbiguities:
                if random.random() < proportion:
                    s.sequence[i] = '-'
                elif random.random() < proportion:
                    s.sequence[i] = ambigs[int(math.floor(len(ambigs) * random.random()))]
            else:
                if random.random() < proportion:
                    s.sequence[i] = '-'
        s.sequence = string.join(s.sequence, '')

Alignment.pepper = pepper
del(pepper)

read('t.nex')
t = var.trees[0]
var.alignments.append(func.newEmptyAlignment(dataType='dna', symbols=None, taxNames=t.taxNames, length=500))
d = Data()
t.data = d
t.newComp(free=0, spec='specified', val=[0.1, 0.2, 0.3])
t.newRMatrix(free=0, spec='specified', val=[1.2, 6.5, 1.3, 9.8, 1.1, 4.5])
t.setPInvar(free=0, val=0.2)
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=0, val=0.5)
t.simulate()
t.data.alignments[0].pepper()
t.data.writeNexus('d.nex')
