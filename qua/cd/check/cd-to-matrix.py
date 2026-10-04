"""cd-to-matrix.md の数式の数値検証。"""

import numpy as np

rng = np.random.default_rng(0)


def rc():
    return complex(*rng.normal(size=2))


def close(x, y):
    return np.allclose(np.asarray(x, dtype=complex), np.asarray(y, dtype=complex))


# 四元数を 4x4 実行列（左乗算）で表し、p + qj の積を直接計算する
def quat(x0, x1, x2, x3):
    return np.array([x0, x1, x2, x3], dtype=float)


def qmul(x, y):
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y
    return quat(
        a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
        a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
        a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
        a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
    )


def from_pair(p, q):  # p + qj（k = ij なので qj = c j + d k）
    return quat(p.real, p.imag, q.real, q.imag)


def to_pair(x):
    return complex(x[0], x[1]), complex(x[2], x[3])


def cd(pq, rs):  # ケイリー＝ディクソンの構成法
    p, q = pq
    r, s = rs
    return (p * r - q * s.conjugate(), p * s + q * r.conjugate())


def M(pq):
    p, q = pq
    return np.array([[p.conjugate(), -q.conjugate()], [q, p]])


def v(pq):
    p, q = pq
    return np.array([p.conjugate(), q])


ok = True


def check(name, cond):
    global ok
    print(("OK  " if cond else "NG  ") + name)
    ok &= bool(cond)


# jp = p*j
p = rc()
j = quat(0, 0, 1, 0)
check("jp = p*j", close(qmul(j, from_pair(p, 0)), from_pair(0, p.conjugate())))

for _ in range(5):
    x, y = (rc(), rc()), (rc(), rc())
    # 四元数の積とケイリー＝ディクソンの構成法の一致
    check("(p+qj)(r+sj) = CD", close(to_pair(qmul(from_pair(*x), from_pair(*y))), cd(x, y)))
    # 行列とベクトルの計算
    pr = cd(x, y)
    check("M(p,q) v(r,s) = v(CD)", close(M(x) @ v(y), v(pr)))
    # 行列の積が積を保つ（記事では未主張）
    check("M(x)M(y) = M(xy)", close(M(x) @ M(y), M(pr)))
    # ベクトルは行列の第1列
    check("v = M の第1列", close(M(x)[:, 0], v(x)))

# 仮計算の段階では一致しないこと（p, -q; q, p で (r, s)）
x, y = (rc(), rc()), (rc(), rc())
p, q = x
r, s = y
naive = np.array([[p, -q], [q, p]]) @ np.array([r, s])
check("素朴な行列では一致しない", not close(naive, cd(x, y)))

# a+bi, c+di での展開と heuristic-cmatrix の行列
a, b, c, d = rng.normal(size=4)
Mh = np.array([[a - b * 1j, -(c - d * 1j)], [c + d * 1j, a + b * 1j]])
check("M = heuristic-cmatrix の行列", close(M((complex(a, b), complex(c, d))), Mh))
I = np.eye(2)
MI = np.array([[-1j, 0], [0, 1j]])
MJ = np.array([[0, -1], [1, 0]])
MK = np.array([[0, 1j], [1j, 0]])
check("M = aI+bM_I+cM_J+dM_K", close(Mh, a * I + b * MI + c * MJ + d * MK))

print("all OK" if ok else "FAILED")
