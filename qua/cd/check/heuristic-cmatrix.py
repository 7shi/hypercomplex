"""heuristic-cmatrix.md の数式の数値検証。"""

import numpy as np

rng = np.random.default_rng(0)

ok = True


def check(name, cond):
    global ok
    print(("OK  " if cond else "NG  ") + name)
    ok &= bool(cond)


def close(x, y):
    return np.allclose(np.asarray(x, dtype=complex), np.asarray(y, dtype=complex))


def qmul(x, y):  # x0 + x1 i + x2 j + x3 k
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y
    return np.array([
        a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
        a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
        a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
        a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
    ])


def naive(x):  # (a+bi) + (c+di)j -> (a+bi, c+di)
    a, b, c, d = x
    return np.array([a + b * 1j, c + d * 1j])


def twist(x):  # (a+bi) + (c+di)j -> (a-bi, c+di)
    a, b, c, d = x
    return np.array([a - b * 1j, c + d * 1j])


one, qi, qj, qk = np.eye(4)
x = rng.normal(size=4)
a, b, c, d = x

# 単純な対応での i, j の左乗算
check("i(a+bi+cj+dk) = (-b+ai)+(-d+ci)j", close(qmul(qi, x), [-b, a, -d, c]))
check("素朴: diag(i,i) で i を表せる", close(np.diag([1j, 1j]) @ naive(x), naive(qmul(qi, x))))
check("j(a+bi+cj+dk) = (-c+di)+(a-bi)j", close(qmul(qj, x), [-c, d, a, -b]))
# 素朴な対応では j の左乗算が複素線形でない：j(i x) と i(j x) の対応が食い違う
lhs = naive(qmul(qj, qmul(qi, x)))
rhs = 1j * naive(qmul(qj, x))
check("素朴: j の左乗算は複素線形でない", not close(lhs, rhs))

# ひねりを加えた表現
MI = np.array([[-1j, 0], [0, 1j]])
MJ = np.array([[0, -1], [1, 0]])
MK = np.array([[0, 1j], [1j, 0]])
I = np.eye(2)
check("k(a+bi+cj+dk) = (-d-ci)+(b+ai)j", close(qmul(qk, x), [-d, -c, b, a]))
for name, q, M in [("i", qi, MI), ("j", qj, MJ), ("k", qk, MK)]:
    check(f"M_{name.upper()} v(x) = v({name}x)", close(M @ twist(x), twist(qmul(q, x))))

# 乗算規則
check("M_I^2 = M_J^2 = M_K^2 = -I", all(close(M @ M, -I) for M in (MI, MJ, MK)))
check("M_I M_J = -M_J M_I = M_K", close(MI @ MJ, MK) and close(MJ @ MI, -MK))
check("M_J M_K = -M_K M_J = M_I", close(MJ @ MK, MI) and close(MK @ MJ, -MI))
check("M_K M_I = -M_I M_K = M_J", close(MK @ MI, MJ) and close(MI @ MK, -MJ))

# 線形結合と第1列
Mx = a * I + b * MI + c * MJ + d * MK
check("aI+bM_I+cM_J+dM_K の成分", close(Mx, [[a - b * 1j, -(c - d * 1j)], [c + d * 1j, a + b * 1j]]))
check("第1列 = v(x)", close(Mx[:, 0], twist(x)))


def M(x):
    a, b, c, d = x
    return a * I + b * MI + c * MJ + d * MK


for _ in range(5):
    x, y = rng.normal(size=4), rng.normal(size=4)
    check("M(x) v(y) = v(xy)", close(M(x) @ twist(y), twist(qmul(x, y))))
    check("M(x) M(y) = M(xy)", close(M(x) @ M(y), M(qmul(x, y))))

# パウリ行列への組み替え
IH, JH, KH = -MK, MJ, MI
check("I_H の成分", close(IH, [[0, -1j], [-1j, 0]]))
check("K_H = I_H J_H", close(IH @ JH, KH))
check("I_H, J_H, K_H の乗算規則",
      all(close(M @ M, -I) for M in (IH, JH, KH))
      and close(JH @ KH, IH) and close(KH @ IH, JH)
      and close(JH @ IH, -KH))
s1 = np.array([[0, 1], [1, 0]])
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]])
check("i I_H, i J_H, i K_H = σ1, σ2, σ3",
      close(1j * IH, s1) and close(1j * JH, s2) and close(1j * KH, s3))

print("all OK" if ok else "FAILED")
