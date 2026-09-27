"""Checks for ktheory/02-vector-fields.md (vector fields on spheres and Clifford modules).

1. Real matrix representations of Cl_{0,k}(R) (generators squaring to -1) on R^{a_k}
   for k = 1..8: complex unit (k=1), quaternion left multiplications (k=2,3),
   octonion left multiplications (k=4..7), and a doubling to R^16 (k=8).
   Each generator is orthogonal and antisymmetric, J_i^2 = -1, J_i J_j = -J_j J_i.
2. Irreducibility and type: the commutant has dimension 2,4,4,4,2,1,1,1 (C, H, H, H, C,
   R, R, R), and the image algebra has dimension 2^k except for k = 3, 7 (2^{k-1}),
   where the pseudoscalar acts as a scalar (+-1), matching the direct-sum types 2H, 2R(8).
   A smaller representation space does not exist: a_k = 2,4,4,8,8,8,8,16 is the real
   dimension of the irreducible module read from the classification table.
3. Vector fields: x -> J_i x (i = 1..k) are tangent to S^{n-1} and orthonormal at every
   point (random points), including block-diagonal modules R^n = (R^{a_k})^m.
4. Linear orthonormal tangent fields force a Clifford module: for A_i with
   <A_i x, x> = 0 and <A_i x, A_j x> = delta_ij |x|^2 for all x, the polarization gives
   A_i^T = -A_i, A_i^T A_i = 1, A_i^T A_j + A_j^T A_i = 0 (checked symbolically for n = 2).
5. Radon-Hurwitz number: rho(n) - 1 = max{k : a_k | n} agrees with rho(n) = 8a + 2^b for
   n = 2^{4a+b} * odd, 0 <= b <= 3 (n = 1..1024); rho(n) = n only for n = 1, 2, 4, 8.
6. Hurwitz: a normed algebra with unit on R^n makes R^n a Cl_{0,n-1} module
   (left multiplication by imaginary units); a_{n-1} | n holds only for n = 1, 2, 4, 8.
   For the octonions, L_7 equals the pseudoscalar L_1...L_6 up to sign.
"""

import itertools

import numpy as np
import sympy as sp

from common.octonion import L as OL

rng = np.random.default_rng(1)


def a(k):
    """Real dimension of an irreducible Cl_{0,k}(R) module."""
    base = [1, 2, 4, 4, 8, 8, 8, 8]
    q, r = divmod(k, 8)
    return base[r] * 16**q


def reps():
    J = np.array([[0, -1], [1, 0]])
    quat = [OL[i][:4, :4] for i in (1, 2, 3)]  # triple (1,2,3) spans H
    oct_ = [OL[i] for i in range(1, 8)]
    Z = np.zeros((8, 8), dtype=int)
    I = np.eye(8, dtype=int)
    e8 = [np.block([[m, Z], [Z, -m]]) for m in oct_] + [np.block([[Z, -I], [I, Z]])]
    return {
        1: [J],
        2: quat[:2],
        3: quat,
        4: oct_[:4],
        5: oct_[:5],
        6: oct_[:6],
        7: oct_,
        8: e8,
    }


def span_dim(mats):
    return np.linalg.matrix_rank(np.array([m.flatten() for m in mats], dtype=float))


def commutant_dim(gens):
    n = gens[0].shape[0]
    rows = []
    for g in gens:  # X g - g X = 0 as a linear map on vec(X)
        rows.append(np.kron(np.eye(n), g.T) - np.kron(g, np.eye(n)))
    A = np.vstack(rows).astype(float)
    return n * n - np.linalg.matrix_rank(A)


def image_algebra(gens):
    n = gens[0].shape[0]
    mats = []
    for r in range(len(gens) + 1):
        for idx in itertools.combinations(range(len(gens)), r):
            m = np.eye(n, dtype=int)
            for i in idx:
                m = m @ gens[i]
            mats.append(m)
    return mats


print("1-2. Cl_{0,k} representations")
R = reps()
expected_comm = {1: 2, 2: 4, 3: 4, 4: 4, 5: 2, 6: 1, 7: 1, 8: 1}
for k, gens in R.items():
    n = gens[0].shape[0]
    I = np.eye(n, dtype=int)
    ok = all((g @ g == -I).all() and (g.T == -g).all() and (g.T @ g == I).all() for g in gens)
    ok &= all((gi @ gj == -gj @ gi).all() for gi, gj in itertools.combinations(gens, 2))
    cd = commutant_dim(gens)
    img = span_dim(image_algebra(gens))
    omega = np.linalg.multi_dot(gens) if len(gens) > 1 else gens[0]
    scalar = (omega == omega[0, 0] * I).all() and omega[0, 0] != 0
    print(f"  k={k}: n={n} (a_k={a(k)}) Clifford/orthogonal={ok} commutant={cd}"
          f" (expected {expected_comm[k]}) image dim={img} / 2^k={2**k}"
          f" pseudoscalar scalar={scalar}")
    assert n == a(k) and ok and cd == expected_comm[k]
    assert img == (2 ** (k - 1) if k in (3, 7) else 2**k)
    assert scalar == (k in (3, 7))

print("3. tangent orthonormal fields")
for k, gens in R.items():
    for m in (1, 2, 3):
        big = [np.kron(np.eye(m, dtype=int), g) for g in gens]
        n = big[0].shape[0]
        for _ in range(20):
            x = rng.normal(size=n)
            x /= np.linalg.norm(x)
            V = np.array([g @ x for g in big])
            assert np.allclose(V @ x, 0)
            assert np.allclose(V @ V.T, np.eye(k))
print("  x -> J_i x orthonormal and tangent: OK (k = 1..8, m = 1..3)")

print("4. linear orthonormal fields => Clifford relations (polarization, n = 2)")
n = 2
x = sp.Matrix(sp.symbols("x1:3"))
A = sp.Matrix(n, n, sp.symbols("a1:5"))
B = sp.Matrix(n, n, sp.symbols("b1:5"))
# <Ax, x> = x^T (A + A^T) x / 2 and <Ax, Bx> = x^T (A^T B + B^T A) x / 2, so vanishing
# for all x is equivalent to the symmetric parts vanishing
assert sp.expand((A * x).dot(x) - (x.T * (A + A.T) * x)[0] / 2) == 0
assert sp.expand((A * x).dot(B * x) - (x.T * (A.T * B + B.T * A) * x)[0] / 2) == 0
print("  quadratic forms <Ax,x> = x^T (A+A^T) x / 2, <Ax,Bx> = x^T (A^T B + B^T A) x / 2: OK")


def rho_formula(n):
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    q, b = divmod(c, 4)
    return 8 * q + 2**b


def rho_clifford(n):
    return 1 + max(k for k in range(0, 64) if n % a(k) == 0)


print("5. Radon-Hurwitz numbers")
for n in range(1, 1025):
    assert rho_formula(n) == rho_clifford(n)
print("  rho(n) =", [rho_formula(n) for n in range(1, 33)], "(n = 1..32)")
par = [n for n in range(1, 1025) if rho_formula(n) == n]
print("  rho(n) = n for n =", par)
assert par == [1, 2, 4, 8]

print("6. Hurwitz dimensions")
hur = [n for n in range(1, 1025) if n % a(n - 1) == 0]
print("  a_{n-1} | n for n =", hur)
assert hur == [1, 2, 4, 8]
# octonions: L_1..L_7 satisfy Clifford relations on R^8, L_7 = +-L_1...L_6
w6 = np.linalg.multi_dot(OL[1:7])
print("  L_1...L_6 = L_7:", (w6 == OL[7]).all(), "/ = -L_7:", (w6 == -OL[7]).all())
assert (w6 == OL[7]).all() or (w6 == -OL[7]).all()
print("All checks passed.")
