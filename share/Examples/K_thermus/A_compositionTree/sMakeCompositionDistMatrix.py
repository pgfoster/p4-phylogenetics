read('../noTRuberNoGapsNoAmbiguities.nex')
a=var.alignments[0]
dm = a.compositionEuclideanDistanceMatrix()
dm.writeNexus('compDistMatrix.nex')
