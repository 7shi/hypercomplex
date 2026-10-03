"""Checks for oct/01-octonion.md.

Verifies:
1. Basis rules: e_i^2 = -1, e_i e_j = -e_j e_i (i != j), e_i e_j = ±e_k with k distinct.
2. The 7 triads of Graves and Cayley: each unit lies in exactly 3 triads,
   each pair of distinct units lies in exactly one triad, and the first triad
   satisfies the quaternion relations.
3. Fano plane layout: e_1, e_2, e_3 on the circle, e_4 at the center, and
   e_5, e_6, e_7 obtained by adding indices (1+4, 2+4, 3+4).
4. Product examples in the article.
5. Associativity: (e_i e_j) e_k = e_i (e_j e_k) holds when {i, j, k} lies in one
   triad (or has repeated indices), and fails (anti-associative) otherwise.
6. Conjugate and norm: o o^* = sum of squares, and |op| = |o||p| symbolically.
7. Alternativity (xx)y = x(xy), (yx)x = y(xx) and flexibility (xy)x = x(yx)
   for general octonions, symbolically.
"""

import itertools

import sympy as sp

from common.octonion import basis_mul, triples


def mul_basis_vec(x, y):
    z = [0] * 8
    for i in range(8):
        if x[i] == 0:
            continue
        for j in range(8):
            if y[j] == 0:
                continue
            s, k = basis_mul(i, j)
            z[k] += s * x[i] * y[j]
    return z


def e(i):
    v = [0] * 8
    v[i] = 1
    return v


def neg(v):
    return [-a for a in v]


def check_basis():
    print("=== 1. 虚数単位の性質 ===")
    for i in range(1, 8):
        assert mul_basis_vec(e(i), e(i)) == neg(e(0))
        for j in range(1, 8):
            if i == j:
                continue
            s, k = basis_mul(i, j)
            assert k not in (0, i, j)
            assert mul_basis_vec(e(i), e(j)) == neg(mul_basis_vec(e(j), e(i)))
    print("e_i^2 = -1, 反交換性, e_i e_j = ±e_k: OK")


def check_triads():
    print("=== 2. 三つ組 ===")
    assert triples == [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
    for i in range(1, 8):
        assert sum(i in t for t in triples) == 3
    for i, j in itertools.combinations(range(1, 8), 2):
        assert sum(i in t and j in t for t in triples) == 1
    for a, b, c in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]:
        assert mul_basis_vec(e(a), e(b)) == e(c)
        assert mul_basis_vec(e(b), e(a)) == neg(e(c))
    print("各単位は3つの三つ組に含まれる, 各ペアは1つの三つ組に含まれる, 四元数の関係: OK")


def check_fano():
    print("=== 3. ファノ平面の配置 ===")
    lines = {frozenset(t) for t in triples}
    assert frozenset({1, 2, 3}) in lines  # circle
    for i in (1, 2, 3):
        assert (i, 4, i + 4) in triples  # lines through the center, oriented i -> 4 -> i+4
    # Remaining 3 lines (triangle sides) do not contain e_4
    rest = [t for t in triples if 4 not in t and set(t) != {1, 2, 3}]
    assert len(rest) == 3
    print("円(1,2,3), 中心e_4を通る線(i,4,i+4), 三角形の辺:", rest, "OK")


def check_examples():
    print("=== 4. 記事中の例 ===")
    m = mul_basis_vec
    assert m(e(1), e(2)) == e(3)
    assert m(e(2), e(1)) == neg(e(3))
    assert m(e(1), e(4)) == e(5)
    assert m(e(3), e(4)) == e(7)
    assert m(e(2), e(4)) == e(6)
    assert m(e(1), e(6)) == neg(e(7))
    assert m(m(e(1), e(2)), e(4)) == e(7)
    assert m(e(1), m(e(2), e(4))) == neg(e(7))
    assert m(m(e(1), e(2)), e(3)) == neg(e(0))
    assert m(e(1), m(e(2), e(3))) == neg(e(0))
    print("(e1e2)e4 = e7, e1(e2e4) = -e7, (e1e2)e3 = e1(e2e3) = -1: OK")


def check_assoc():
    print("=== 5. 結合性と三つ組 ===")
    lines = [set(t) for t in triples]
    for i, j, k in itertools.product(range(1, 8), repeat=3):
        lhs = mul_basis_vec(mul_basis_vec(e(i), e(j)), e(k))
        rhs = mul_basis_vec(e(i), mul_basis_vec(e(j), e(k)))
        s = {i, j, k}
        closed = len(s) < 3 or s in lines
        if closed:
            assert lhs == rhs, (i, j, k)
        else:
            assert lhs == neg(rhs), (i, j, k)
    print("三つ組内（または添字の重複あり）なら結合的, それ以外は反結合的: OK")


def check_norm():
    print("=== 6. 共役とノルム ===")
    o = sp.symbols("o0:8", real=True)
    p = sp.symbols("p0:8", real=True)
    oc = [o[0]] + [-x for x in o[1:]]
    oo = mul_basis_vec(list(o), oc)
    assert sp.expand(oo[0] - sum(x**2 for x in o)) == 0
    assert all(sp.expand(x) == 0 for x in oo[1:])
    oco = mul_basis_vec(oc, list(o))
    assert all(sp.expand(a - b) == 0 for a, b in zip(oo, oco))
    op = mul_basis_vec(list(o), list(p))
    lhs = sum(x**2 for x in op)
    rhs = sum(x**2 for x in o) * sum(x**2 for x in p)
    assert sp.expand(lhs - rhs) == 0
    print("oo* = o*o = Σo_i^2, |op|^2 = |o|^2|p|^2: OK")


def check_alternative():
    print("=== 7. 交代性 ===")
    x = list(sp.symbols("x0:8", real=True))
    y = list(sp.symbols("y0:8", real=True))
    m = mul_basis_vec
    pairs = [
        (m(m(x, x), y), m(x, m(x, y))),
        (m(m(y, x), x), m(y, m(x, x))),
        (m(m(x, y), x), m(x, m(y, x))),
    ]
    for lhs, rhs in pairs:
        assert all(sp.expand(a - b) == 0 for a, b in zip(lhs, rhs))
    print("(xx)y = x(xy), (yx)x = y(xx), (xy)x = x(yx): OK")


if __name__ == "__main__":
    check_basis()
    check_triads()
    check_fano()
    check_examples()
    check_assoc()
    check_norm()
    check_alternative()
