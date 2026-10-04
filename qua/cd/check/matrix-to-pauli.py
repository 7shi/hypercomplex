"""matrix-to-pauli.md の数式の数値検証。"""

import numpy as np

rng = np.random.default_rng(0)

ok = True


def check(name, cond):
    global ok
    print(("OK  " if cond else "NG  ") + name)
    ok &= bool(cond)


def close(x, y):
    return np.allclose(np.asarray(x, dtype=complex), np.asarray(y, dtype=complex))


I = np.eye(2)
JH = np.array([[0, -1], [1, 0]])

# J_H は複素数 i の行列表現
x, y = rng.normal(size=2)
check("J_H (x,y) = (-y,x)", close(JH @ [x, y], [-y, x]))

# I_H J_H = -J_H I_H から b=c, d=-a
a, b, c, d = rng.normal(size=4)
X = np.array([[a, b], [c, d]])
check("XJ_H, -J_H X の成分", close(X @ JH, [[b, -a], [d, -c]]) and close(-JH @ X, [[c, d], [-a, -b]]))
A = np.array([[a, b], [b, -a]])
check("反交換する形", close(A @ JH, -JH @ A))
check("(a b; b -a)^2 = (a^2+b^2)I", close(A @ A, (a * a + b * b) * I))


def IH_(t):
    return 1j * np.array([[np.sin(t), np.cos(t)], [np.cos(t), -np.sin(t)]])


def KH_(t):
    return 1j * np.array([[np.cos(t), -np.sin(t)], [-np.sin(t), -np.cos(t)]])


for t in rng.uniform(0, 2 * np.pi, size=3):
    IH, KH = IH_(t), KH_(t)
    check("I_H^2 = -I（パラメータ表示）", close(IH @ IH, -I))
    check("K_H = I_H J_H（パラメータ表示）", close(IH @ JH, KH))
    check("K_H^2 = -I", close(KH @ KH, -I))
    check("J_H K_H = I_H, K_H I_H = J_H", close(JH @ KH, IH) and close(KH @ IH, JH))
    # 別のパラメータ表示 a=-i cos, b=i sin は θ のずれ
    alt = np.array([[-1j * np.cos(t), 1j * np.sin(t)], [1j * np.sin(t), 1j * np.cos(t)]])
    check("別のパラメータ表示 = θ-π/2", close(alt, IH_(t - np.pi / 2)))

IH, KH = IH_(np.pi), KH_(np.pi)
check("θ=π の I_H", close(IH, [[0, -1j], [-1j, 0]]))
check("θ=π の K_H", close(KH, [[-1j, 0], [0, 1j]]))

# 転置・共役・エルミート共役
check("J_H^T = -J_H", close(JH.T, -JH))
check("I_H^T = I_H, K_H^T = K_H", close(IH.T, IH) and close(KH.T, KH))
check("I_H^* = -I_H, K_H^* = -K_H", close(IH.conj(), -IH) and close(KH.conj(), -KH))
check("I_H, J_H, K_H は反エルミート", all(close(M.conj().T, -M) for M in (IH, JH, KH)))


# 実行列表現での転置
def real4(M):
    R = np.zeros((4, 4))
    for r in range(2):
        for s in range(2):
            z = M[r, s]
            R[2 * r:2 * r + 2, 2 * s:2 * s + 2] = [[z.real, -z.imag], [z.imag, z.real]]
    return R


M = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
(a, b), (c, d), (e, f), (g, h) = [(z.real, z.imag) for z in M.flatten()]
check("4x4 実行列の成分", close(real4(M), [[a, -b, c, -d], [b, a, d, c], [e, -f, g, -h], [f, e, h, g]]))
check("実行列の転置 = エルミート共役の実行列", close(real4(M).T, real4(M.conj().T)))

# パウリ行列
s1, s2, s3 = 1j * IH, 1j * JH, 1j * KH
check("σ1, σ2, σ3 の成分",
      close(s1, [[0, 1], [1, 0]]) and close(s2, [[0, -1j], [1j, 0]]) and close(s3, [[1, 0], [0, -1]]))
check("σ^2 = I", all(close(s @ s, I) for s in (s1, s2, s3)))
check("σ はエルミート", all(close(s.conj().T, s) for s in (s1, s2, s3)))
check("θ=0 では iI_H = -σ1", close(1j * IH_(0), -s1))
check("σ1σ2 = -σ2σ1 = iσ3 = -K_H", close(s1 @ s2, -s2 @ s1) and close(s1 @ s2, 1j * s3) and close(s1 @ s2, -KH))
check("σ2σ3 = -σ3σ2 = iσ1 = -I_H", close(s2 @ s3, -s3 @ s2) and close(s2 @ s3, 1j * s1) and close(s2 @ s3, -IH))
check("σ3σ1 = -σ1σ3 = iσ2 = -J_H", close(s3 @ s1, -s1 @ s3) and close(s3 @ s1, 1j * s2) and close(s3 @ s1, -JH))
check("σ3σ2 = I_H, σ1σ3 = J_H, σ2σ1 = K_H", close(s3 @ s2, IH) and close(s1 @ s3, JH) and close(s2 @ s1, KH))
check("σ1σ2σ3 = iI", close(s1 @ s2 @ s3, 1j * I))

# 双四元数の一般形（2行目と3行目の一致）
co = rng.normal(size=8)
lhs = (co[0] * I + co[1] * s1 + co[2] * s2 + co[3] * s3
       + co[4] * s3 @ s2 + co[5] * s1 @ s3 + co[6] * s2 @ s1 + co[7] * s1 @ s2 @ s3)
rhs = (co[0] * I + co[1] * s1 + co[2] * s2 + co[3] * s3
       - co[4] * 1j * s1 - co[5] * 1j * s2 - co[6] * 1j * s3 + co[7] * 1j * I)
check("双四元数の一般形の行列", close(lhs, rhs))
# 4行目以降の係数が四元数 i, j, k と h に対応すること
check("i ≅ I_H, j ≅ J_H, k ≅ K_H, h ≅ iI",
      close(s3 @ s2, IH) and close(s1 @ s3, JH) and close(s2 @ s1, KH) and close(s1 @ s2 @ s3, 1j * I))

print("all OK" if ok else "FAILED")
