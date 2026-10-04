"""real-to-cmatrix.md の数式の数値検証。"""

from itertools import permutations

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


E = np.eye(4)
one, qi, qj, qk = E


def left(q):  # 左乗算の 4x4 実行列（列 = 基底の像）
    return np.column_stack([qmul(q, e) for e in E])


Mi, Mj, Mk = left(qi), left(qj), left(qk)
x = rng.normal(size=4)
a, b, c, d = x

check("i x = (-b, a, -d, c)", close(qmul(qi, x), [-b, a, -d, c]))
check("j x = (-c, d, a, -b)", close(qmul(qj, x), [-c, d, a, -b]))
check("k x = (-d, -c, b, a)", close(qmul(qk, x), [-d, -c, b, a]))
check("M_i の成分", close(Mi, [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]]))
check("M_j の成分", close(Mj, [[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0]]))
check("M_k の成分", close(Mk, [[0, 0, 0, -1], [0, 0, -1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]))


def is_complex(M):  # 2x2 ブロックがすべて (α -β; β α) の形か
    for r in range(2):
        for s in range(2):
            B = M[2 * r:2 * r + 2, 2 * s:2 * s + 2]
            if not (np.isclose(B[0, 0], B[1, 1]) and np.isclose(B[0, 1], -B[1, 0])):
                return False
    return True


def to_complex(M):
    return np.array([[M[2 * r, 2 * s] + 1j * M[2 * r + 1, 2 * s] for s in range(2)] for r in range(2)])


check("(a,b,c,d) では M_i のみ複素行列に変換できる",
      is_complex(Mi) and not is_complex(Mj) and not is_complex(Mk))

names = "abcd"
good = []
for perm in permutations(range(4)):
    P = E[list(perm)]  # 並べ替えたベクトル = P x
    if all(is_complex(P @ M @ P.T) for M in (Mi, Mj, Mk)):
        good.append("(" + ",".join(names[p] for p in perm) + ")")
listed = "(a,b,d,c) (a,c,b,d) (a,d,c,b) (b,a,c,d) (b,c,d,a) (b,d,a,c) (c,a,d,b) (c,b,a,d) (c,d,b,a) (d,a,b,c) (d,b,c,a) (d,c,a,b)".split()
check("24通り中12通りが複素行列に変換できる", len(good) == 12)
check("12通りの一覧が記事と一致", sorted(good) == sorted(listed))

# (d,a,b,c) の並べ替え
P = E[[3, 0, 1, 2]]
check("P の成分", close(P, [[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]))
check("M_r' = P M_r P^{-1}", close(P @ Mi @ np.linalg.inv(P), P @ Mi @ P.T))
Mi2, Mj2, Mk2 = (P @ M @ P.T for M in (Mi, Mj, Mk))
check("M_i' の成分", close(Mi2, [[0, 0, 0, 1], [0, 0, -1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]]))
check("M_j' の成分", close(Mj2, [[0, 0, -1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, 1, 0, 0]]))
check("M_k' の成分", close(Mk2, [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]]))
check("並べ替えたベクトルへの作用",
      close(Mi2 @ P @ x, [c, -b, a, -d]) and close(Mj2 @ P @ x, [-b, -c, d, a]) and close(Mk2 @ P @ x, [a, -d, -c, b]))

IH, JH, KH = to_complex(Mi2), to_complex(Mj2), to_complex(Mk2)
check("I_H, J_H, K_H の成分",
      close(IH, [[0, -1j], [-1j, 0]]) and close(JH, [[0, -1], [1, 0]]) and close(KH, [[-1j, 0], [0, 1j]]))
v = np.array([d + 1j * a, b + 1j * c])
check("I_H v", close(IH @ v, [c - 1j * b, a - 1j * d]))
check("J_H v", close(JH @ v, [-b - 1j * c, d + 1j * a]))
check("K_H v", close(KH @ v, [a - 1j * d, -c + 1j * b]))

I = np.eye(2)
check("I_H J_H = K_H, J_H I_H = -K_H", close(IH @ JH, KH) and close(JH @ IH, -KH))
check("I_H^2 = J_H^2 = K_H^2 = -I", all(close(M @ M, -I) for M in (IH, JH, KH)))
check("J_H K_H = I_H, K_H I_H = J_H", close(JH @ KH, IH) and close(KH @ IH, JH))

# ベクトルの -i 倍と線形結合の第1列
check("-i (d+ia, b+ic) = (a-id, c-ib)", close(-1j * v, [a - 1j * d, c - 1j * b]))
Q = a * I + b * IH + c * JH + d * KH
check("線形結合の成分", close(Q, [[a - 1j * d, -c - 1j * b], [c - 1j * b, a + 1j * d]]))
check("第1列 = -i v", close(Q[:, 0], -1j * v))


def rep(x):
    a, b, c, d = x
    return a * I + b * IH + c * JH + d * KH


for _ in range(5):
    x, y = rng.normal(size=4), rng.normal(size=4)
    check("rep(x) rep(y) = rep(xy)", close(rep(x) @ rep(y), rep(qmul(x, y))))

# パウリ行列
s1 = np.array([[0, 1], [1, 0]])
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]])
check("i I_H = σ1, i J_H = σ2, i K_H = σ3", close(1j * IH, s1) and close(1j * JH, s2) and close(1j * KH, s3))
check("I_H = -iσ1 など", close(IH, -1j * s1) and close(JH, -1j * s2) and close(KH, -1j * s3))

print("all OK" if ok else "FAILED")
