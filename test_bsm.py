from neuron import h
from time import time

h.register_BooleanSparseMatrix()

rank = 10
ms = [h.BooleanSparseMatrix(rank, rank) for _ in range(2)]

r = h.Random()
r.Random123(1, 2, 3)
frac = 0.2
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
    set(m)
    m.pr()


def test_clone():
    m1 = ms[0].clone()
    assert m1.equal(ms[0])
    m2 = m1.clone()
    m2.setelm(3, 4, m2.getelm(3, 4) == 0)
    assert not m2.equal(m1)


def test_elmmul():
    m1 = ms[0].clone()
    m1.elmmul(ms[1])
    print("Test element multiplication")
    m1.pr()


def test_add():
    m1 = ms[0].clone()
    m1.add(ms[1])
    print("Test element addition")
    m1.pr()


def test_vector_getset():
    print("test_vector_getset")
    m1 = ms[0]
    m2 = h.BooleanSparseMatrix(m1.nrow(), m1.ncol())
    nrow = int(m1.nrow())
    for i in range(nrow):
        m2.setrow(i, m1.getrow(i))
    assert m2.equal(m1)


if __name__ == "__main__":
    test_clone()
    test_elmmul()
    test_add()
    test_vector_getset()
