
def hook(self, theTree):
    """Write stuff during the MCMC.

    You get the stuff from the Mcmc (which is self) or from the
    current tree.  Write whatever you like, however you like it.  You
    will need to know how to get it, tho ...

    In this demo, I write out the gen number+1, the composition, and
    the pInvar.

    You could put this hook in your ~/.p4 directory, if you like.
    """

    global hookFile
    hookFile.write("%10i " % (self.gen + 1))
    theCompVal = theTree.model.parts[0].comps[0].val
    hookFile.write("  %9.6f%9.6f%9.6f%9.6f" % (theCompVal[0], theCompVal[1], theCompVal[2], theCompVal[3]))
    hookFile.write("     %9.6f\n" % theTree.model.parts[0].pInvar.val)

Mcmc.hook = hook
del(hook)

read("../d.nex")
d = Data()
t = func.randomTree(taxNames=d.taxNames)
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=1, val=0.2)
m = Mcmc(t, nChains=1, runNum=0, sampleInterval=100, checkPointInterval=None)
hookFile = file("mcmc_hookFile", "w")
m.run(2000)
hookFile.close()



