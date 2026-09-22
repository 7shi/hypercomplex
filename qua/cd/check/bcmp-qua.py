"""bcmp-qua.md の数式の数値検証。"""

import numpy as np

rng = np.random.default_rng(0)


def rc():
    return complex(*rng.normal(size=2))


# 複素数を対にした元 (a, b) = a + bj の積
def mul_bc(p, q):
    a, b = p
    c, d = q
    return (a * c - b * d, a * d + b * c)


def mul_qu(p, q):
    a, b = p
    c, d = q
    return (a * c - b * d.conjugate(), a * d + b * c.conjugate())


def conj_bc(p):  # i と j を反転（k は不変）
    a, b = p
    return (a.conjugate(), -b.conjugate())


def conj_qu(p):  # i, j, k をすべて反転
    a, b = p
    return (a.conjugate(), -b)


def comps(p):  # x0 + x1 i + x2 j + x3 k（k = ij）
    a, b = p
    return np.array([a.real, a.imag, b.real, b.imag])


def M_bc(p):
    a, b = p
    return np.array([[a, -b], [b, a]])


def M_qu(p):
    a, b = p
    return np.array([[a, b], [-b.conjugate(), a.conjugate()]])


def close(x, y):
    return np.allclose(np.asarray(x, dtype=complex), np.asarray(y, dtype=complex))


def rand():
    return (rc(), rc())


# 複素数の共役積
a1, a2, b1, b2 = rng.normal(size=4)
al, be = complex(a1, a2), complex(b1, b2)
assert close(al.conjugate() * be, complex(a1 * b1 + a2 * b2, a1 * b2 - a2 * b1))
assert close(be.conjugate() * al, (al.conjugate() * be).conjugate())
A = np.array([[a1, -a2], [a2, a1]])
assert close(A.T @ [b1, b2], [a1 * b1 + a2 * b2, a1 * b2 - a2 * b1])

# 追い越し規則: 四元数では ji = -ij、双複素数では ji = ij
one, i, j = (1 + 0j, 0j), (1j, 0j), (0j, 1 + 0j)
k = mul_qu(i, j)
assert close(mul_qu(j, i), (-k[0], -k[1]))
assert close(mul_bc(j, i), mul_bc(i, j))
assert close(mul_qu(j, j), (-1, 0)) and close(mul_bc(j, j), (-1, 0))

for _ in range(20):
    p, q, r = rand(), rand(), rand()
    # 結合法則
    assert close(mul_qu(mul_qu(p, q), r), mul_qu(p, mul_qu(q, r)))
    assert close(mul_bc(mul_bc(p, q), r), mul_bc(p, mul_bc(q, r)))
    # 行列表現が積を保つ
    assert close(M_qu(p) @ M_qu(q), M_qu(mul_qu(p, q)))
    assert close(M_bc(p) @ M_bc(q), M_bc(mul_bc(p, q)))
    # M(p) v(q) = v(pq)、v(q) は M(q) の第1列
    c, d = q
    e, f = mul_qu(p, q)
    assert close(M_qu(p) @ [c, -d.conjugate()], [e, -f.conjugate()])
    assert close(M_qu(q)[:, 0], [c, -d.conjugate()])
    # 行列式
    a, b = p
    assert close(np.linalg.det(M_qu(p)), abs(a) ** 2 + abs(b) ** 2)
    assert close(np.linalg.det(M_bc(p)), a * a + b * b)
    # 共役とエルミート共役
    assert close(M_qu(conj_qu(p)), M_qu(p).conj().T)
    assert close(M_bc(conj_bc(p)), M_bc(p).conj().T)
    # 四元数の共役積
    assert close(mul_qu(conj_qu(p), q), (a.conjugate() * c + b * d.conjugate(),
                                         a.conjugate() * d - b * c.conjugate()))
    x, y = comps(p), comps(q)
    s = comps(mul_qu(conj_qu(p), q))
    t = comps(mul_qu(conj_qu(q), p))
    assert close(s[0], x @ y) and close(t[0], s[0]) and close(t[1:], -s[1:])
    assert close(mul_qu(conj_qu(mul_qu(p, q)), one),
                 mul_qu(conj_qu(q), conj_qu(p)))
    xv, yv = x[1:], y[1:]
    assert close(s[1:], x[0] * yv - y[0] * xv - np.cross(xv, yv))
    assert close(comps(mul_qu(conj_qu(p), p)), [x @ x, 0, 0, 0])
    # 双複素数の共役積: 実部は内積、k 成分は対称、i, j 成分は反対称
    s = comps(mul_bc(conj_bc(p), q))
    t = comps(mul_bc(conj_bc(q), p))
    assert close(mul_bc(conj_bc(p), q), (a.conjugate() * c + b.conjugate() * d,
                                         a.conjugate() * d - b.conjugate() * c))
    assert close(s[0], x @ y)
    assert close(t[[0, 3]], s[[0, 3]]) and close(t[1:3], -s[1:3])
    # p*p の k 成分は 2 Im(a* b)
    assert close(comps(mul_bc(conj_bc(p), p)),
                 [x @ x, 0, 0, 2 * (a.conjugate() * b).imag])

# 基底とパウリ行列
s1 = np.array([[0, 1], [1, 0]])
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]])
assert close(M_qu(i), 1j * s3)
assert close(M_qu(j), 1j * s2)
assert close(M_qu(k), 1j * s1)

# 双複素数の零因子と共役積の例
ij = mul_bc(i, j)
assert close(mul_bc((1, 1j), (1, -1j)), (0, 0))
assert close(np.linalg.det(M_bc((1, 1j))), 0)
p = (1j, 1 + 0j)  # i + j
assert close(comps(mul_bc(conj_bc(p), p)), [2, 0, 0, -2])
assert close(ij, (0, 1j))

print("OK")
