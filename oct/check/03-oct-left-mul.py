"""Checks for oct/03-oct-left-mul.md (claims not covered by the other 03-*.py).

Verifies:
1. The matrix L_1 shown in the article, L_i^2 = -I, L_i L_j = -L_j L_i, and that
   L_i L_j (i != j) is not the left multiplication by any octonion (L_1 L_2 != L_3).
2. The sequence o_0, ..., o_7 = o_0 obtained by left-multiplying e_1, ..., e_7,
   L_7 L_6 ... L_1 = I, every step of the derivation, L_1 ... L_6 = L_7 and
   L_1 ... L_7 = -I.
3. The 64 increasing monomials in L_1, ..., L_6 are linearly independent (span M_8(R)),
   are mutually orthogonal under the Frobenius inner product, every non-identity
   monomial has trace 0 and anticommutes with some monomial C.
4. Transpose of a grade-k monomial is (-1)^k times the reversed product, B B^T = I,
   and c_k = Tr(B_k^T A)/8 recovers the coefficients.
5. The formulas R_i = (-L_i + three products)/2 for all i.
6. The expansion of P in the 7 triple products.
7. The first row of R_i^{-1} = -R_i has its only nonzero entry 1 at position i.
8. Proof ingredients added after review: L_u^2 = -|u|^2 I for pure imaginary u,
   the chain computing W(1) = 1, L_x L_x = L_{x^2}, L_x L_y P = L_{xy} P,
   and A P = L_{A(1)} P for arbitrary A.
"""

import itertools

import numpy as np

from common.octonion import L, R, P, basis_mul, triples

I = np.eye(8, dtype=int)


def prod(idx):
    M = I.copy()
    for i in idx:
        M = M @ L[i]
    return M


def check_l1():
    print("=== 1. L_1 と反交換性 ===")
    L1 = np.array(
        [
            [0, -1, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, -1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, -1, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 0, -1, 0],
        ]
    )
    assert (L[1] == L1).all()
    for i in range(1, 8):
        assert (L[i] @ L[i] == -I).all()
        for j in range(1, 8):
            if i != j:
                assert (L[i] @ L[j] == -L[j] @ L[i]).all()
    # L_i L_j is not L_a for any octonion a: it would have to be L_{(L_iL_j)(1)}
    for i, j in itertools.permutations(range(1, 8), 2):
        M = L[i] @ L[j]
        a = M[:, 0]
        La = sum(a[k] * L[k] for k in range(8))
        assert not (M == La).all()
    assert not (L[1] @ L[2] == L[3]).all()
    print("  OK: L_1 の行列, L_i^2=-I, 反交換, L_iL_j はどの左作用とも一致しない")


def oct_mul(x, y):
    z = np.zeros(8, dtype=int)
    for i in range(8):
        for j in range(8):
            if x[i] and y[j]:
                s, k = basis_mul(i, j)
                z[k] += s * x[i] * y[j]
    return z


def check_volume():
    print("=== 2. 7回の左乗算と体積要素 ===")
    expected = [
        [1, 1, 1, 1, 1, 1, 1, 1],
        [-1, 1, -1, 1, -1, 1, 1, -1],
        [1, 1, -1, -1, -1, 1, -1, 1],
        [1, 1, 1, 1, -1, -1, -1, -1],
        [1, -1, -1, -1, 1, -1, -1, -1],
        [1, -1, -1, 1, -1, 1, -1, 1],
        [1, -1, 1, 1, -1, -1, 1, -1],
        [1, 1, 1, 1, 1, 1, 1, 1],
    ]
    o = np.array(expected[0])
    for k in range(1, 8):
        e = np.zeros(8, dtype=int)
        e[k] = 1
        o = oct_mul(e, o)
        assert (o == np.array(expected[k])).all(), k
    assert (prod([7, 6, 5, 4, 3, 2, 1]) == I).all()
    # derivation steps: (-1)^m L_7...L_{m+1} = L_1...L_m
    for m in range(1, 8):
        lhs = (-1) ** m * prod(range(7, m, -1))
        assert (lhs == prod(range(1, m + 1))).all(), m
    assert (prod(range(1, 7)) == L[7]).all()
    assert (prod(range(1, 8)) == -I).all()
    print("  OK: o_1..o_7, L_7...L_1=I, 各段階, L_1...L_6=L_7, L_1...L_7=-I")


def monomials():
    for k in range(7):
        for c in itertools.combinations(range(1, 7), k):
            yield c


def check_basis():
    print("=== 3. 64個の基底とトレース ===")
    mons = list(monomials())
    Bs = [prod(c) for c in mons]
    A = np.array([B.flatten() for B in Bs])
    assert len(Bs) == 64 and np.linalg.matrix_rank(A) == 64
    G = A @ A.T
    assert (G == 8 * np.eye(64, dtype=int)).all()
    for c, B in zip(mons[1:], Bs[1:]):
        assert np.trace(B) == 0
        assert any((B @ C == -C @ B).all() for C in Bs), c
    print("  OK: 線形独立, フロベニウス内積で直交 (Tr(B_j^T B_k)=8δ), Tr=0, 反交換する C が存在")
    return mons, Bs


def check_transpose(mons, Bs):
    print("=== 4. 転置と成分抽出 ===")
    for c, B in zip(mons, Bs):
        k = len(c)
        assert (B.T == (-1) ** k * prod(reversed(c))).all()
        assert (B @ B.T == I).all()
    rng = np.random.default_rng(0)
    coef = rng.integers(-5, 6, size=64)
    A = sum(a * B for a, B in zip(coef, Bs))
    rec = [np.trace(B.T @ A) / 8 for B in Bs]
    assert np.allclose(rec, coef)
    print("  OK: B^T=(-1)^k 逆順, BB^T=I, c_k=Tr(B_k^T A)/8")


def check_right():
    print("=== 5. 右作用の表現 ===")
    for i in range(1, 8):
        pairs = []
        for a, b, c in triples:
            for x, y, z in [(a, b, c), (b, c, a), (c, a, b)]:
                if z == i:
                    pairs.append((x, y))
        assert len(pairs) == 3
        rhs = (-L[i] + sum(L[x] @ L[y] for x, y in pairs)) / 2
        assert np.allclose(R[i], rhs), i
    # the explicit forms in the article
    forms = {
        1: [(2, 3), (4, 5), (7, 6)],
        2: [(3, 1), (4, 6), (5, 7)],
        3: [(1, 2), (6, 5), (4, 7)],
        4: [(5, 1), (6, 2), (7, 3)],
        5: [(1, 4), (7, 2), (3, 6)],
        6: [(1, 7), (2, 4), (5, 3)],
        7: [(6, 1), (2, 5), (3, 4)],
    }
    for i, ps in forms.items():
        rhs = (-L[i] + sum(L[x] @ L[y] for x, y in ps)) / 2
        assert np.allclose(R[i], rhs), i
    print("  OK: R_i = (-L_i + e_i を生成する3対の積)/2 (i=1..7)")


def check_p():
    print("=== 6. 射影 P の展開 ===")
    terms = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
    rhs = (I - sum(prod(t) for t in terms)) / 8
    assert np.allclose(P, rhs)
    print("  OK: P = (I - Σ L_aL_bL_c)/8")


def check_rinv():
    print("=== 7. R_i^{-1} の第1行 ===")
    for i in range(1, 8):
        Rinv = -R[i]
        assert (Rinv @ R[i] == I).all()
        row = Rinv[0]
        assert row[i] == 1 and np.count_nonzero(row) == 1
    print("  OK: R_i^{-1} = -R_i の第1行は位置 i だけが 1")


def Lx(x):
    return sum(x[k] * L[k] for k in range(8))


def check_review_additions():
    print("=== 8. 証明に用いた関係 ===")
    rng = np.random.default_rng(1)
    e = np.eye(8, dtype=int)
    # W(1) chain: e2e1=-e3, e3(-e3)=1, e4*1=e4, e5e4=-e1, e6(-e1)=-e7, e7(-e7)=1
    v = e[0]
    expect = [e[1], -e[3], e[0], e[4], -e[1], -e[7], e[0]]
    for k, ex in zip(range(1, 8), expect):
        v = L[k] @ v
        assert (v == ex).all(), k
    for _ in range(20):
        u = np.zeros(8)
        u[1:] = rng.normal(size=7)
        assert np.allclose(Lx(u) @ Lx(u), -(u @ u) * np.eye(8))
        x, y = rng.normal(size=8), rng.normal(size=8)
        xx = Lx(x) @ x
        assert np.allclose(Lx(x) @ Lx(x), Lx(xx))
        xy = Lx(x) @ y
        assert np.allclose(Lx(x) @ Lx(y) @ P, Lx(xy) @ P)
        A = rng.normal(size=(8, 8))
        assert np.allclose(A @ P, Lx(A[:, 0]) @ P)
    print("  OK: W(1) の計算, L_u^2=-|u|^2 I, L_xL_x=L_{x^2}, L_xL_yP=L_{xy}P, AP=L_{A(1)}P")


if __name__ == "__main__":
    check_l1()
    check_volume()
    mons, Bs = check_basis()
    check_transpose(mons, Bs)
    check_right()
    check_p()
    check_rinv()
    check_review_additions()
