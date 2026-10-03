"""Checks for oct/nonassociativity.md.

Verifies:
1. Index correspondence (Graves and Cayley): i, j, k, l -> e_1, e_2, e_3, e_4,
   il, jl, kl -> e_5, e_6, e_7, and (il)^2 = (jl)^2 = (kl)^2 = -1.
2. The alternative definition (cyclic triads (a, a+1, a+3) mod 7):
   i, j, k, l -> e_1, e_2, e_4, e_7, li = e_3, lj = e_6, lk = e_5, and that it is
   an octonion algebra (alternative, multiplicative norm).
3. Non-associativity and anti-associativity chains for (e_1 e_2) e_4.
4. The zero divisor (e_3 - e_1)(e_4 + e_6) = 0 when e_1 e_6 = e_7 is imposed.
5. Examples of alternativity, the 3-factor remark, the triad example,
   e_5 e_6 = -e_3 and (il)(jl) = -k.
6. Summary: for basis units, (e_a e_b) e_c is associative iff two indices coincide
   or {a, b, c} is a triad, and anti-associative otherwise; for general octonions,
   (xy)z is neither equal to x(yz) nor to -x(yz) in general.
"""

import itertools
import random

import sympy as sp

from common.octonion import basis_mul, triples


def make_mul(table):
    def mul(x, y):
        z = [0] * 8
        for i in range(8):
            if x[i] == 0:
                continue
            for j in range(8):
                if y[j] == 0:
                    continue
                s, k = table(i, j)
                z[k] += s * x[i] * y[j]
        return z

    return mul


def table_from_triples(ts):
    mul = {}
    for a, b, c in ts:
        for x, y, z in [(a, b, c), (b, c, a), (c, a, b)]:
            mul[(x, y)] = (1, z)
            mul[(y, x)] = (-1, z)

    def table(i, j):
        if i == 0:
            return (1, j)
        if j == 0:
            return (1, i)
        if i == j:
            return (-1, 0)
        return mul[(i, j)]

    return table


m = make_mul(basis_mul)


def e(i):
    v = [0] * 8
    v[i] = 1
    return v


def neg(v):
    return [-a for a in v]


def check_gc_index():
    print("=== 1. グレイブスとケイリーの添字対応 ===")
    i, j, k, l = e(1), e(2), e(3), e(4)
    assert m(i, j) == k
    assert m(i, l) == e(5) and m(j, l) == e(6) and m(k, l) == e(7)
    for x in (m(i, l), m(j, l), m(k, l)):
        assert m(x, x) == neg(e(0))
    print("il=e5, jl=e6, kl=e7, (il)^2=(jl)^2=(kl)^2=-1: OK")


def check_alt_definition():
    print("=== 2. 別の定義（巡回的な三つ組） ===")
    ts = [(a, (a % 7) + 1, ((a + 2) % 7) + 1) for a in range(1, 8)]
    print("三つ組:", ts)
    mm = make_mul(table_from_triples(ts))
    i, j, k, l = e(1), e(2), e(4), e(7)
    assert mm(i, j) == k
    assert mm(l, i) == e(3) and mm(l, j) == e(6) and mm(l, k) == e(5)
    assert mm(e(7), e(1)) == e(3)
    assert mm(e(6), e(7)) == e(2)
    assert mm(e(4), e(5)) == e(7)
    x = list(sp.symbols("x0:8", real=True))
    y = list(sp.symbols("y0:8", real=True))
    for lhs, rhs in [(mm(mm(x, x), y), mm(x, mm(x, y))), (mm(mm(y, x), x), mm(y, mm(x, x)))]:
        assert all(sp.expand(a - b) == 0 for a, b in zip(lhs, rhs))
    xy = mm(x, y)
    assert sp.expand(sum(c**2 for c in xy) - sum(c**2 for c in x) * sum(c**2 for c in y)) == 0
    print("i,j,k,l -> e1,e2,e4,e7, li=e3, lj=e6, lk=e5, 交代性・ノルムの乗法性: OK")


def check_chains():
    print("=== 3. 非結合性・反結合性 ===")
    lhs = m(m(e(1), e(2)), e(4))
    assert lhs == neg(m(e(4), m(e(1), e(2))))  # (e1e2)e4 = -e4(e1e2)
    assert lhs == e(7)
    assert neg(m(e(1), m(e(2), e(4)))) == e(7)
    # anti-associativity steps
    assert lhs == neg(m(e(1), m(e(2), e(4))))
    assert m(e(1), m(e(4), e(2))) == neg(m(m(e(1), e(4)), e(2)))
    assert m(m(e(4), e(1)), e(2)) == neg(m(e(4), m(e(1), e(2))))
    print("(e1e2)e4 = e7 = -e1(e2e4) = -e4(e1e2): OK")


def check_zero_divisor():
    print("=== 4. 零因子 ===")
    ts = [t for t in triples if set(t) != {1, 7, 6}] + [(1, 6, 7)]
    mm = make_mul(table_from_triples(ts))
    a = [0] * 8
    a[3], a[1] = 1, -1
    b = [0] * 8
    b[4], b[6] = 1, 1
    assert mm(e(1), e(6)) == e(7)
    assert mm(a, b) == [0] * 8
    assert m(a, b) != [0] * 8
    print("e1e6=e7 とすると (e3-e1)(e4+e6)=0: OK")


def check_examples():
    print("=== 5. 記事中の例 ===")
    # alternativity example
    l1 = m(m(e(1), e(2)), m(e(3), e(1)))
    assert l1 == m(e(3), m(e(3), e(1))) == m(m(e(3), e(3)), e(1))
    assert m(m(e(3), e(3)), e(1)) == m(m(m(e(1), e(2)), e(3)), e(1))
    # remark: 3 distinct factors
    assert m(m(e(1), e(2)), m(e(2), e(4))) == neg(m(m(e(3), e(2)), e(4)))
    # triad
    assert m(m(e(1), e(2)), e(3)) == m(e(1), m(e(2), e(3)))
    # e5 e6 chain
    assert m(e(5), e(6)) == neg(e(3))
    steps = [
        m(e(5), m(e(2), e(4))),
        neg(m(m(e(5), e(2)), e(4))),
        neg(m(m(m(e(1), e(4)), e(2)), e(4))),
        m(m(e(1), m(e(4), e(2))), e(4)),
        neg(m(m(e(1), m(e(2), e(4))), e(4))),
        m(m(m(e(1), e(2)), e(4)), e(4)),
        m(m(e(1), e(2)), m(e(4), e(4))),
        neg(e(3)),
    ]
    assert all(s == steps[0] for s in steps)
    i, j, l = e(1), e(2), e(4)
    il, jl = m(i, l), m(j, l)
    chain = [
        m(il, jl),
        neg(m(m(il, j), l)),
        m(m(i, m(l, j)), l),
        neg(m(m(i, jl), l)),
        m(m(m(i, j), l), l),
        m(m(i, j), m(l, l)),
        neg(e(3)),
    ]
    assert all(c == chain[0] for c in chain)
    print("交代性の例, 3因子の rem, 三つ組, e5e6=-e3, (il)(jl)=-k: OK")


def check_summary():
    print("=== 6. まとめの分類 ===")
    lines = [set(t) for t in triples]
    for a, b, c in itertools.product(range(1, 8), repeat=3):
        lhs = m(m(e(a), e(b)), e(c))
        rhs = m(e(a), m(e(b), e(c)))
        if len({a, b, c}) < 3 or {a, b, c} in lines:
            assert lhs == rhs
        else:
            assert lhs == neg(rhs)
    assert m(e(3), e(2)) == neg(e(1))
    x = [1, 1, 0, 0, 0, 0, 0, 0]  # 1 + e1
    assert m(m(x, e(2)), e(4)) == [0, 0, 0, 0, 0, 0, 1, 1]  # e6 + e7
    assert m(x, m(e(2), e(4))) == [0, 0, 0, 0, 0, 0, 1, -1]  # e6 - e7
    random.seed(1)
    x, y, z = ([random.randint(-3, 3) for _ in range(8)] for _ in range(3))
    lhs = m(m(x, y), z)
    rhs = m(x, m(y, z))
    assert lhs != rhs and lhs != neg(rhs)
    print("基底の虚数単位では結合的/反結合的の二択: OK")
    print("一般の八元数の例:", x, y, z, "では (xy)z ≠ ±x(yz)")


if __name__ == "__main__":
    check_gc_index()
    check_alt_definition()
    check_chains()
    check_zero_divisor()
    check_examples()
    check_summary()
