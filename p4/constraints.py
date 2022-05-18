import p4.func
from p4.p4exceptions import P4Error


class Constraints(object):

    """A container for tree topology constraints.

    We may have usual constraints, and we are able ask whether a
    query tree has those constraints.  In addition we may have
    constraints defined on the root clades, and we are able to
    ask if a query tree is consistent with the number and content
    of those clades.  So for example we can use this to enforce a
    bifurcating root in a fixed position, where the content of
    the two clades off the root is fixed, but where the topology
    within the two clades can vary.
    
    Init args ---

    taxNames
              A list of taxNames in the same order as in
              the data or alignment, and the same order as
              in other tree or trees objects.

    constraintTree
              A partially resolved tree object that
              describes the constraints.  You need to include
              all the taxNames.

    rootConstraintTree
              A tree with splits that define the root (possibly
              with other constraints as well).  A biRooted
              constraint tree should be biRooted, and a triRooted
              constraint tree should be triRooted.  You need to
              include all the taxNames.  Constraints in the
              rootConstraintTree are added to the constraints
              from the constraintTree, so you might be able to
              define all the constraints you need with a
              rootConstraintTree.  


    For example::

      tNames = list("ABCDE")
      cTree = func.readAndPop('(A, B, (E, D), C);')
      rTree = func.readAndPop("((A, B), (C, D, E));")
      constr = Constraints(tNames, constraintTree=cTree, rootConstraintTree=rTree)
      t = func.randomTree(taxNames=tNames, constraints=constr, biRoot=True)

    A biRooted rTree (rootConstrainTree) might be (A,B),(C,D)); while a triRooted 
    rTree might be (A,B,(C,D)); or ((A,B),(C,D),(E,F));

    You can pass a Constraints object to p4.func.randomTree() and
    Mcmc() to enforce constraints.

    """

    def __init__(self, taxNames, constraintTree=None, rootConstraintTree=None):

        self.taxNames = taxNames
        self.cTree = constraintTree
        self.rTree = rootConstraintTree
        self.constraints = []
        self.rootConstraints = []
        self.allOnes = 2 ** (len(self.taxNames)) - 1

        if self.cTree:
            self.cTree.taxNames = taxNames
        if self.rTree:
            self.rTree.taxNames = taxNames
            
        # If there are two trees, are they compatible?
        if self.cTree and self.rTree:
            isCompatible = self.cTree.isCompatibleWith(self.rTree)
            assert isCompatible, "cTree and rTree are not compatible"
                
        if self.cTree:
            self.cTree.makeSplitKeys()
            for n in self.cTree.iterInternalsNoRoot():
                n.name = n.br.splitKey
                if n.br.splitKey not in self.constraints:
                    self.constraints.append(n.br.splitKey)
        if self.rTree:
            self.rTree.makeSplitKeys()
            for n in self.rTree.iterInternalsNoRoot():
                n.name = n.br.splitKey
                if n.parent == self.rTree.root:
                    if n.br.splitKey not in self.rootConstraints:
                        self.rootConstraints.append(n.br.splitKey)
                # Root constraints are also constraints
                if n.br.splitKey not in self.constraints:
                    self.constraints.append(n.br.splitKey)

        assert self.constraints or self.rootConstraints, "No constraints?"

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
        """Ask if aTree root is consistent with self.rootConstraints

        Additionally checks whether the degree of the root (ie biRoot,
        triRoot, etc) of aTree is the same as self.rTree.

        Returns True or False
        """
        aTreeRootSplits = [n.br.splitKey for n in aTree.root.iterChildren()]
        areConsistent = True

        # First check splits in self.rootConstraints
        for sk in self.rootConstraints:
            if not sk in aTreeRootSplits:
                areConsistent = False
                break
        if not areConsistent:
            return False

        # Check degree of the root (biRoot, triRoot, etc)
        selfNroot = self.rTree.root.getNChildren()
        aTreeNroot = len(aTreeRootSplits)
        if selfNroot != aTreeNroot:
            print(f"Constraints.rTree root has degree {selfNroot}")
            print(f"aTree root has degree {aTreeNroot}")
            areConsistent = False
        return areConsistent    

    def dump(self):
        """A summary"""

        print('Constraints.dump()')
        print('taxNames:')
        for i,txN in enumerate(self.taxNames):
            print('    %3i  %s' % (i, txN))
        print('constraints (includes root contraints if present):')
        for cn in self.constraints:
            print(cn)
        if self.cTree:
            self.cTree.draw()
        print('rootConstraints:')
        for i in self.rootConstraints:
            print(p4.func.getSplitStringFromKey(i, self.rTree.nTax), i)
        if self.rTree:
            self.rTree.draw()
