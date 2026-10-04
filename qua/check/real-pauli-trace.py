"""real-pauli-trace.md の数式の数値検証。"""

from itertools import product

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
s1 = np.array([[0, 1], [1, 0]])
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]])

# パウリ行列の関係式
check("σ^2 = I", all(close(s @ s, I) for s in (s1, s2, s3)))
check("σ1σ2 = -σ2σ1 = iσ3", close(s1 @ s2, -s2 @ s1) and close(s1 @ s2, 1j * s3))
check("σ2σ3 = -σ3σ2 = iσ1", close(s2 @ s3, -s3 @ s2) and close(s2 @ s3, 1j * s1))
check("σ3σ1 = -σ1σ3 = iσ2", close(s3 @ s1, -s1 @ s3) and close(s3 @ s1, 1j * s2))

# 実パウリ行列
t1, t2, t3 = s1, -1j * s2, s3
check("-iσ2 = (0 -1; 1 0)", close(t2, [[0, -1], [1, 0]]))
check("σ1(-iσ2) = σ3", close(s1 @ t2, s3))
check("τ1^2 = I, τ2^2 = -I, τ3^2 = I", close(t1 @ t1, I) and close(t2 @ t2, -I) and close(t3 @ t3, I))
check("τ1τ2 = -τ2τ1 = τ3", close(t1 @ t2, t3) and close(t2 @ t1, -t3))
check("τ1τ3 = -τ3τ1 = τ2", close(t1 @ t3, t2) and close(t3 @ t1, -t2))
check("τ2τ3 = -τ3τ2 = τ1", close(t2 @ t3, t1) and close(t3 @ t2, -t1))
x, y = rng.normal(size=2)
check("τ2 (x,y) = (-y,x)", close(t2 @ [x, y], [-y, x]))

# 行列単位
E11 = np.array([[1, 0], [0, 0]])
E12 = np.array([[0, 1], [0, 0]])
E21 = np.array([[0, 0], [1, 0]])
E22 = np.array([[0, 0], [0, 1]])
check("E11 = (I+τ3)/2", close(E11, (I + t3) / 2))
check("E12 = (τ1-τ2)/2", close(E12, (t1 - t2) / 2))
check("E21 = (τ1+τ2)/2", close(E21, (t1 + t2) / 2))
check("E22 = (I-τ3)/2", close(E22, (I - t3) / 2))
p, q, r, s = rng.normal(size=4)
check("(p q; r s) の展開",
      close([[p, q], [r, s]], (p + s) / 2 * I + (q + r) / 2 * t1 + (-q + r) / 2 * t2 + (p - s) / 2 * t3))
B = np.array([m.flatten() for m in (I, t1, t2, t3)]).real
check("{I, τ1, τ2, τ3} は一次独立", np.linalg.matrix_rank(B) == 4)

# トレース
check("tr τ = 0", all(np.isclose(np.trace(t), 0) for t in (t1, t2, t3)))
basis = [I, t1, t2, t3]


def comb(c):
    return sum(ci * b for ci, b in zip(c, basis))


a = rng.normal(size=4)
b = rng.normal(size=4)
A, Bm = comb(a), comb(b)
check("tr A = 2a0", np.isclose(np.trace(A), 2 * a[0]))
check("A の成分", close(A, [[a[0] + a[3], a[1] - a[2]], [a[1] + a[2], a[0] - a[3]]]))
check("A^T = a0 I + a1 τ1 - a2 τ2 + a3 τ3", close(A.T, comb([a[0], a[1], -a[2], a[3]])))
check("tr(A^T B) = 2 a·b", np.isclose(np.trace(A.T @ Bm), 2 * a @ b))
# 交差項が消えること：異なる基底の積 X^T Y のトレースは 0
cross = all(np.isclose(np.trace(X.T @ Y), 0) for (m, X), (n, Y) in product(enumerate(basis), repeat=2) if m != n)
check("異なる基底の tr(X^T Y) = 0", cross)
check("同じ基底の tr(X^T X) = 2", all(np.isclose(np.trace(X.T @ X), 2) for X in basis))
check("a_r = tr(τ_r^T A)/2", all(np.isclose(np.trace(X.T @ A) / 2, a[n]) for n, X in enumerate(basis)))
check("a_2 = -tr(τ2 A)/2, τ2^T = -τ2", np.isclose(-np.trace(t2 @ A) / 2, a[2]) and close(t2.T, -t2))
check("tr(AB) = 2(a0b0+a1b1-a2b2+a3b3)", np.isclose(np.trace(A @ Bm), 2 * (a[0] * b[0] + a[1] * b[1] - a[2] * b[2] + a[3] * b[3])))
check("tr(A^T B) = 成分ごとの積の和（フロベニウス内積）", np.isclose(np.trace(A.T @ Bm), np.sum(A * Bm)))

print("all OK" if ok else "FAILED")
