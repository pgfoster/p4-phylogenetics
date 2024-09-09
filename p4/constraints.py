import p4.func
from p4.p4exceptions import P4Error


class Constraints(object):

    """A container for tree topology constraints.

    Constraints are defined with a partially-resolved tree.  We can
    then use a constraints object to ask whether a query tree has
    those constraints.  In addition a constrained rooting may be
    defined, and we can then ask whether a query tree is
    consistent with that rooting.  So for example we can use this
    to enforce a bifurcating root in a fixed position, where the
    content of the two clades off the root is fixed, but where the
    topology within the two clades can vary. (Also possible with a 
    biRoot'ed tree.)
    
    Init args ---

    taxNames
              A list of taxNames in the same order as in
              the data or alignment, and the same order as
              in other tree or trees objects.

    constraintTree
              A partially resolved tree object that
              describes the constraints.  You need to include
              all the taxNames, but they do not need to be in 
              a particular order.

    rooting
              (Boolean -- False by default) A constraint tree 
              can define the root (possibly with other constraints 
              as well).  If rooting is set, a biRooted
              constraint tree should be biRooted, and a triRooted
              constraint tree should be triRooted.  You need to
              include all the taxNames.  

    For example, this will constrain A+C, without rooting, in a 
    random biRoot'ed tree::

      tNames = list("ABCDE")
      cTree = func.readAndPop("((A, C), B, D, E);")
      constr = Constraints(tNames, constraintTree=cTree, rooting=False)
      t = func.randomTree(taxNames=tNames, constraints=constr, biRoot=True)

    A biRooted cTree might be ((A,B,C),(D,E,F)); while a triRooted 
    cTree might be (A,B,(C,D,E)); or ((A,B,C),(D,E,F),(G,H,I));  If 
    you turn on rooting, then those rootings will be enforced.

    Constraints can be nested, done with the constraintTree description, 
    as for example (A,B,(C,D,(E,F),G)). You can have nested constraints 
    when rooting is set.

    You can pass a Constraints object to p4.func.randomTree() and to
    Mcmc() to enforce constraints.

    """

    def __init__(self, taxNames, constraintTree=None, rooting=False):

        gm = ["Constraints init()"]
        self.taxNames = taxNames
        self.cTree = constraintTree
        self.rooting = rooting
        self.constraints = []
        self.allOnes = 2 ** (len(self.taxNames)) - 1

        if self.cTree:
            self.cTree.taxNames = taxNames
            self.cTree.makeSplitKeys()
            for n in self.cTree.iterInternalsNoRoot():
                # n.name = n.br.splitKey
                if n.br.splitKey not in self.constraints:
                    self.constraints.append(n.br.splitKey)

        if self.rooting:
            isCTreeBiRoot = self.cTree.isBiRoot()
            if not isCTreeBiRoot:
                isCTreeTriRoot = self.cTree.isTriRoot()
                if not isCTreeTriRoot:
                    gm.append("Constraints cTree is neither biRoot nor triRoot")
                    gm.append("When rooting is on, the tree should be one or the other")
                    raise P4Error(gm)

        assert self.constraints, "No constraints?"

    def areConsistentWithTree(self, aTree):
        """Ask whether aTree is consistent with the constraints in self

        It is the responsibility of the caller to aTree.makeSplitKeys()

        Returns True or False
        """
        assert aTree.nodeForSplitKeyDict
        areConsistent = True
        for sk in self.constraints:
            ret = aTree.nodeForSplitKeyDict.get(sk)
            if not ret:
                areConsistent = False
                break
        return areConsistent

    def areConsistentWithTreeRoot(self, aTree):
        """Ask if aTree root is consistent with constraints

        Additionally checks whether the degree of the root (ie biRoot,
        triRoot, etc) of aTree is the same as self.cTree.

        Returns True or False
        """

        assert self.rooting, "Constraints.areConsistentWithTreeRoot() rooting is not turned on."

        aTreeRootSplits = [n.br.splitKey for n in aTree.root.iterChildren()]
        areConsistent = True

        # First check splits in self.cTre
        selfRootSplits = [n.br.splitKey for n in self.cTree.root.iterChildren()]
        for sk in selfRootSplits:
            if not sk in aTreeRootSplits:
                areConsistent = False
                break
        if not areConsistent:
            return False

        # Check degree of the root (biRoot, triRoot, etc)
        selfNroot = len(selfRootSplits)
        aTreeNroot = len(aTreeRootSplits)
        if selfNroot != aTreeNroot:
            print(f"Constraints.cTree root has degree {selfNroot}")
            print(f"aTree root has degree {aTreeNroot}")
            areConsistent = False
        return areConsistent    

    def dump(self):
        """A summary"""

        print('Constraints.dump()')
        print(f'rooting: {self.rooting}')
        print('taxNames:')
        for i,txN in enumerate(self.taxNames):
            print('    %3i  %s' % (i, txN))
        print('constraints:')
        for cn in self.constraints:
            print(cn)
        if self.cTree:
            self.cTree.draw()
