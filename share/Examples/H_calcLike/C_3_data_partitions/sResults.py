read('opt.p4_tPickle')
t = var.trees[0]

print('The rate matrix for the second partition is:')
print(t.model.parts[1].rMatrices[0].val)

print('The composition of the model in the last partition is:')
print(t.model.parts[2].comps[0].val)

print('The relative rate of the first partition is:')
print(t.model.parts[0].relRate)

# To get the lot, say
#t.model.dump()
