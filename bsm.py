from neuron import h
from time import time

h.register_BooleanSparseMatrix()

rank = 100000
ms = [h.BooleanSparseMatrix(rank, rank) for _ in range(2)]

r = h.Random()
r.Random123(1, 2, 3)
frac = 0.001
# average interval between true in a row is 1/frac
r.negexp(1 / frac)


def rpick():
    return int(r.repick())


def set(m):
    for i in range(int(m.nrow())):
        j = rpick()
        while j < int(m.ncol()):
            m.setelm(i, j, True)
            j += rpick()


for m in ms:
    start = time()
    set(m)
    m.pr()
    print("set time %g" % (time() - start))


start = time()
ms[0].elmmul(ms[1]).pr()
print("elmmul time %g" % (time() - start))
