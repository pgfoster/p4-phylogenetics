var.verboseRead = 0
read('d.nex')
d = Data()
print("\n") # put in a couple of spaces
ret = d.compoChiSquaredTest(verbose=1)
print("""
The p4 output above is there because 'verbose' was turned on.  For
programmatic usage, you would want to turn it off.  In that case you
can get the numbers from the list of lists that the method returns.
In this case, it was

%s
""" % ret)
