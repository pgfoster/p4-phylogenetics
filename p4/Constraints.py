import func
from Glitch import Glitch

class Constraints(object):
    """A container for tree topology constraints.

    taxNames
              A list of taxNames in the same order as in
              the data or alignment, and the same order as
              in other tree or trees objects.

    constraintTree
              A partially resolved tree object that
              describes the constraints.  You need to include
              all the taxNames.


    For example::
    
      tNames = ['A', 'B', 'C', 'D', 'E']
      read('(A, B, (E, D), C);')
      constTree = var.trees.pop()
      c = Constraints(tNames, constTree)
      t = func.randomTree(taxNames=tName, constraints=c)

    The constraint tree should not have a bifurcating root.

    You can pass a Constraints object to func.randomTree() and
    Mcmc() to enforce constraints.

    """
    def __init__(self, taxNames, constraintTree):
        
        if constraintTree.root.getNChildren() == 2:
            raise Glitch, 'The constraint tree should not have a bifurcating root.'
        self.tree = constraintTree
        self.tree.taxNames = taxNames
        self.tree.reRoot(self.tree.taxNames[0], moveInternalName=False)
        self.allOnes = 2L**(self.tree.nTax) - 1

        self.tree.makeSplitKeys()
        self.constraints = []
        internalsExceptTheFirst = [n for n in self.tree.iterInternalsNoRoot()][1:]
        for n in internalsExceptTheFirst:
            n.name = n.br.splitKey
            self.constraints.append(n.br.splitKey)
        assert self.constraints, "No constraints?"
        

    def dump(self):
        print 'Constraints.dump()'
        print 'taxNames:'
        for i in range(self.tree.nTax):
            print '    %3i  %s' % (i, self.tree.taxNames[i])
        print 'constraints:'
        for i in self.constraints:
            print func.getSplitStringFromKey(i, self.tree.nTax)
        self.tree.draw()

