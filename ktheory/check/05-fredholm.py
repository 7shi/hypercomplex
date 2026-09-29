"""Checks for ktheory/05-fredholm.md (Fredholm index and Toeplitz operators).

Conventions: l^2 = sequences (a_0, a_1, ...) = Hardy space H^2 (Fourier coefficients of
nonnegative frequency). S is the right shift (S e_n = e_{n+1}), S^* the left shift.
For a Laurent polynomial f = sum_k c_k z^k the Toeplitz operator T_f = P M_f has matrix
(T_f)_{mn} = c_{m-n} (m, n >= 0), so T_z = S and T_{z^{-1}} = S^*.

Infinite operators are represented by M x M sections; since all operators here are banded,
the top-left N x N block of a product of sections is exact when N + (bandwidths) <= M.

1. Shifts: S^* S = I, S S^* = I - P_0, [S^*, S] = P_0 (trace 1). For the N x N section J
   of S, J^* J - J J^* = diag(1, 0, ..., 0, -1) (trace 0), ker J = e_{N-1}, ker J^* = e_0.
   The cyclic shift X of clif/03 has X^* X - X X^* = 0.
2. S - a = (I - a S^*) S, and (I - a S^*)^{-1} = sum a^n (S^*)^n for |a| < 1;
   S - b = -b (I - S / b) for |b| > 1.
3. Products: T_{conj(g) f h} = T_{conj(g)} T_f T_h for polynomials g, h (nonnegative
   frequencies), while T_z T_{z^{-1}} != T_1. Norm bound ||T_f|| <= max |f|.
4. Index formula: for f = z^{-m} p (p a polynomial without zeros on |z| = 1), the numerical
   winding number equals #{zeros of p in |z| < 1} - m, so -wind f = m - #{inside zeros}.
5. Finite sections: T_N(f) is square (index 0), but it has exactly |wind f| small singular
   values. If wind f < 0 their right singular vectors sit near e_0 (kernel of T_f) and the
   left ones near e_{N-1}; if wind f > 0 the roles swap (cokernel of T_f near e_0).
6. Additivity of the index for finite matrices (dim V1 - dim V3) and the kernel/cokernel
   bookkeeping used in the proof of the product lemma.
7. Sign convention: compressing M_f to the negative frequencies z^{-1}, z^{-2}, ... and
   identifying z^n <-> z^{-n-1} gives the Toeplitz operator with symbol f(z^{-1}), whose
   winding number is -wind f (so that index is +wind f).
"""

import numpy as np

rng = np.random.default_rng(5)


def toeplitz(coef, M):
    """M x M section of T_f for f = sum coef[k] z^k."""
    T = np.zeros((M, M), dtype=complex)
    for k, c in coef.items():
        T += c * np.eye(M, k=-k)
    return T


def mulpoly(a, b):
    c = {}
    for i, x in a.items():
        for j, y in b.items():
            c[i + j] = c.get(i + j, 0) + x * y
    return c


def conjpoly(a):
    return {-k: np.conj(c) for k, c in a.items()}


def evalpoly(a, z):
    return sum(c * z ** k for k, c in a.items())


def winding(f, n=20000):
    z = np.exp(2j * np.pi * np.arange(n + 1) / n)
    w = f(z)
    return np.sum(np.angle(w[1:] / w[:-1])) / (2 * np.pi)


M, N = 80, 40
I = np.eye(M)
S = toeplitz({1: 1}, M)
Ss = S.conj().T
blk = lambda A: A[:N, :N]

print("1. shifts")
assert np.allclose(blk(Ss @ S), np.eye(N))
P0 = np.zeros((N, N))
P0[0, 0] = 1
assert np.allclose(blk(S @ Ss), np.eye(N) - P0)
assert np.allclose(blk(Ss @ S - S @ Ss), P0)
print("  S^*S = I, SS^* = I - P_0, [S^*, S] = P_0: OK")
J = toeplitz({1: 1}, N)
D = J.conj().T @ J - J @ J.conj().T
assert np.allclose(D, np.diag([1] + [0] * (N - 2) + [-1]))
assert np.isclose(np.trace(D), 0)
assert np.linalg.matrix_rank(J) == N - 1
e = np.eye(N)
assert np.allclose(J @ e[N - 1], 0) and np.allclose(J.conj().T @ e[0], 0)
X = J.copy()
X[0, N - 1] = 1
assert np.allclose(X.conj().T @ X - X @ X.conj().T, 0)
assert np.linalg.matrix_rank(X - J) == 1
print("  J^*J - JJ^* = diag(1,0,...,0,-1), ker J = e_{N-1}, coker J = e_0: OK")
print("  cyclic X: X^*X = XX^*, X - J has rank 1: OK")

print("2. linear symbols")
for _ in range(5):
    a = 0.9 * rng.uniform() * np.exp(2j * np.pi * rng.uniform())
    assert np.allclose(blk((I - a * Ss) @ S), blk(S - a * I))
    inv = sum(a ** n * np.linalg.matrix_power(Ss, n) for n in range(400))
    assert np.allclose((I - a * Ss) @ inv, I)
    b = (1.1 + rng.uniform()) * np.exp(2j * np.pi * rng.uniform())
    assert np.allclose(S - b * I, -b * (I - S / b))
print("  S - a = (I - aS^*)S, Neumann series inverts I - aS^*: OK")

print("3. products")
for _ in range(5):
    g = {k: rng.normal() + 1j * rng.normal() for k in range(3)}
    h = {k: rng.normal() + 1j * rng.normal() for k in range(4)}
    f = {k: rng.normal() + 1j * rng.normal() for k in range(-3, 4)}
    lhs = toeplitz(mulpoly(mulpoly(conjpoly(g), f), h), M)
    rhs = toeplitz(conjpoly(g), M) @ toeplitz(f, M) @ toeplitz(h, M)
    assert np.allclose(blk(lhs), blk(rhs))
    th = np.linspace(0, 2 * np.pi, 4001)
    sup = np.max(np.abs(evalpoly(f, np.exp(1j * th))))
    assert np.linalg.norm(toeplitz(f, 200), 2) <= sup + 1e-9
assert not np.allclose(blk(S @ Ss), np.eye(N))
print("  T_{conj(g) f h} = T_conj(g) T_f T_h, T_z T_{1/z} != 1, ||T_f|| <= sup|f|: OK")

print("4. index formula via factorization")
for _ in range(20):
    m = int(rng.integers(0, 4))
    roots = []
    for _ in range(int(rng.integers(1, 6))):
        r = rng.choice([rng.uniform(0.1, 0.8), rng.uniform(1.3, 3)])
        roots.append(r * np.exp(2j * np.pi * rng.uniform()))
    c = rng.normal() + 1j * rng.normal()
    f = lambda z: c * z ** (-m) * np.prod([z - r for r in roots], axis=0)
    inside = sum(abs(r) < 1 for r in roots)
    assert np.isclose(winding(f), inside - m)
print("  wind(z^-m p) = #inside zeros - m: OK")

print("5. finite sections")
Nb = 300
cases = [
    {1: 1},
    {2: 1},
    {-1: 1},
    {-2: 1, 0: 0.3},
    {0: -0.5, 1: 1},
    mulpoly({-1: 1}, mulpoly({0: -0.3, 1: 1}, {0: -2, 1: 1})),
    mulpoly({0: -0.2, 1: 1}, {0: 0.4j, 1: 1}),
    mulpoly({-3: 1}, {0: -0.5, 1: 1}),
    {0: 1, 1: 0.3},
]
for coef in cases:
    w = round(winding(lambda z: evalpoly(coef, z)))
    T = toeplitz(coef, Nb)
    U, s, Vh = np.linalg.svd(T)
    small = s < 1e-6
    assert small.sum() == abs(w), (coef, w, s[-5:])
    head = lambda v: np.sum(np.abs(v[: Nb // 2]) ** 2)
    right = [head(Vh[i].conj()) for i in np.where(small)[0]]
    left = [head(U[:, i]) for i in np.where(small)[0]]
    if w < 0:
        assert all(r > 0.99 for r in right) and all(l < 0.01 for l in left)
    if w > 0:
        assert all(r < 0.01 for r in right) and all(l > 0.99 for l in left)
    print(f"  wind = {w:+d}: {small.sum()} small singular values, localized as expected: OK")

print("6. additivity for finite matrices")
for _ in range(20):
    n1, n2, n3 = rng.integers(1, 7, size=3)
    r1 = int(rng.integers(0, min(n1, n2) + 1))
    r2 = int(rng.integers(0, min(n2, n3) + 1))
    A = rng.normal(size=(n2, r1)) @ rng.normal(size=(r1, n1))
    B = rng.normal(size=(n3, r2)) @ rng.normal(size=(r2, n2))
    ind = lambda T: T.shape[1] - T.shape[0]  # dim ker - dim coker = dim V - dim W
    assert ind(B @ A) == ind(A) + ind(B)
    rk = np.linalg.matrix_rank
    kerA, cokA = n1 - rk(A), n2 - rk(A)
    kerB = n2 - rk(B)
    kerBA, cokBA = n1 - rk(B @ A), n3 - rk(B @ A)
    # dim ker BA = dim ker A + dim(ker B ∩ im A)
    KB = np.linalg.svd(B)[2][rk(B):].conj().T if kerB else np.zeros((n2, 0))
    capdim = KB.shape[1] + rk(A) - rk(np.hstack([KB, A])) if kerB else 0
    assert kerBA == kerA + capdim
    assert kerBA - cokBA == (kerA - cokA) + (kerB - (n3 - rk(B)))
print("  ind(BA) = ind A + ind B, dim ker BA = dim ker A + dim(ker B ∩ im A): OK")
print("7. negative frequencies")
for _ in range(5):
    f = {k: rng.normal() + 1j * rng.normal() for k in range(-3, 4)}
    # matrix of the compression on basis z^{-1-n}: <f z^{-1-n}, z^{-1-m}> = c_{n-m}
    Tneg = np.zeros((M, M), dtype=complex)
    for m_ in range(M):
        for n_ in range(M):
            Tneg[m_, n_] = f.get(n_ - m_, 0)
    finv = {-k: c for k, c in f.items()}
    assert np.allclose(Tneg, toeplitz(finv, M))
    if min(abs(evalpoly(f, np.exp(2j * np.pi * np.arange(1000) / 1000)))) > 1e-3:
        assert np.isclose(winding(lambda z: evalpoly(finv, z)), -winding(lambda z: evalpoly(f, z)))
print("  compression to negative frequencies = T_{f(1/z)}, wind f(1/z) = -wind f: OK")
print("All checks passed.")
