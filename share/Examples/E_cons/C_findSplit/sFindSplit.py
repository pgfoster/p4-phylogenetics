allTaxNames = ['t0', 't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 't10', 't11', 't12', 't13', 't14', 't15', 't16', 't17', 't18', 't19', 't20', 't21', 't22', 't23', 't24', 't25', 't26', 't27', 't28', 't29', 't30', 't31', 't32', 't33', 't34', 't35', 't36', 't37', 't38', 't39', 't40', 't41', 't42', 't43', 't44']

some = ['t20', 't23']
splitKey = func.getSplitKeyFromTaxNames(allTaxNames, some)
splitString = func.getSplitStringFromKey(splitKey, len(allTaxNames), escaped=True)
print(some)
os.system('grep "%s" splits' % splitString)

some = ['t27', 't36']
splitKey = func.getSplitKeyFromTaxNames(allTaxNames, some)
splitString = func.getSplitStringFromKey(splitKey, len(allTaxNames), escaped=True)
print(some)
os.system('grep "%s" splits' % splitString)
