import p4.func
from p4.p4exceptions import P4Error


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
      t = p4.func.randomTree(taxNames=tName, constraints=c)

    You can pass a Constraints object to p4.func.randomTree() and
    Mcmc() to enforce constraints.

    """

    def __init__(self, taxNames, constraintTree):

        self.tree = constraintTree
        self.tree.taxNames = taxNames
        self.allOnes = 2 ** (self.tree.nTax) - 1

        self.tree.makeSplitKeys()
        self.constraints = []
        #self.tree.draw()
        for n in self.tree.iterInternalsNoRoot():
            n.name = n.br.splitKey
            if n.br.splitKey not in self.constraints:
                self.constraints.append(n.br.splitKey)
        assert self.constraints, "No constraints?"

    def dump(self):
        print('Constraints.dump()')
        print('taxNames:')
        for i in range(self.tree.nTax):
            print('    %3i  %s' % (i, self.tree.taxNames[i]))
        print('constraints:')
        for i in self.constraints:
            print(p4.func.getSplitStringFromKey(i, self.tree.nTax), i)
        self.tree.draw()
