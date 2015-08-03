read(""" 2 2
one
ac
two
gt
""")
read('(one,two);')
t = var.trees[0]
t.data = Data()
t.newComp()
t.newRMatrix()
t.setPInvar()
t.calcLogLike()
