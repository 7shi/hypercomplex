"""Numerical checks for 08-g2.md.

Automorphisms of the octonions: phi(1) = 1, preservation of the
real/imaginary split and the norm, orthogonality with det 1, and the
equivalence with cross-product preservation on Im O; conjugation
rho_r being an automorphism exactly for triple-rotation angles
0, +-2pi/3 (r^3 real), and the diagonal-phase reading
theta_1 + theta_2 + theta_3 = 0 of the same condition; derivations
(Leibniz rule) forming the Lie algebra g2 with dim Der(O) = 14,
closed under commutator inside so(7); ad_u failing Leibniz by one
associator [u, zw] - [u,z]w - z[u,w] = -(u,z,w) and the corrected
D_{x,y}(z) = [[x,y],z] - 3(x,y,z) being a derivation, with such maps
spanning all of Der(O); the decomposition so(7) = g2 (+) {u cross .}
(21 = 14 + 7, direct sum); exp of a derivation being an automorphism
while a generic rotation exp of so(7) is not; basic triples
(u1, u2, u3) determining automorphisms (transitivity + freeness) and
the dimension count 14 = 6 + 5 + 3; the stabilizer chain: derivations
fixing e1 forming su(3) (dim 8, anti-Hermitian traceless on C^3 with
i = left multiplication by e1, exp unitary with det 1), fixing e1, e2
dim 3 (su(2)), fixing e1, e2, e4 dim 0.
"""

from itertools import combinations

import numpy as np

rng = np.random.default_rng(8)


def conj(a):
    c = -a.copy()
    c[0] = a[0]
    return c


def mul(a, b):
    """Cayley-Dickson product: (a1,a2)(b1,b2) = (a1 b1 - b2* a2, b2 a1 + a2 b1*)."""
    n = len(a)
    if n == 1:
        return a * b
    h = n // 2
    a1, a2, b1, b2 = a[:h], a[h:], b[:h], b[h:]
    return np.concatenate([
        mul(a1, b1) - mul(conj(b2), a2),
        mul(b2, a1) + mul(a2, conj(b1)),
    ])


def basis(n, i):
    v = np.zeros(n)
    v[i] = 1.0
    return v


def norm(a):
    return np.sqrt(a @ a)


def im(a):
    b = a.copy()
    b[0] = 0.0
    return b


def rand(n=8):
    return rng.standard_normal(n)


def rand_unit_im():
    v = im(rand())
    return v / norm(v)


def comm(a, b):
    return mul(a, b) - mul(b, a)


def assoc(a, b, c):
    return mul(mul(a, b), c) - mul(a, mul(b, c))


def cross(a, b):
    return im(mul(a, b))


def oexp(x):
    t, v = x[0], im(x)
    r = norm(v)
    u = v / r if r > 1e-12 else 0 * v
    e = basis(len(x), 0)
    return np.exp(t) * (np.cos(r) * e + np.sin(r) * u)


def mexp(A):
    """Matrix exponential by scaling and squaring with a Taylor series."""
    k = max(0, int(np.ceil(np.log2(max(np.abs(A).sum(), 1e-16)))) + 1)
    B = A / 2**k
    E, term = np.eye(len(A)), np.eye(len(A))
    for n in range(1, 30):
        term = term @ B / n
        E = E + term
    for _ in range(k):
        E = E @ E
    return E


E = [basis(8, i) for i in range(8)]
ONE = E[0]

# Left/right multiplication matrices on O
L = [np.column_stack([mul(E[i], E[j]) for j in range(8)]) for i in range(8)]
R = [np.column_stack([mul(E[j], E[i]) for j in range(8)]) for i in range(8)]


def is_auto(F, trials=20):
    ok = True
    for _ in range(trials):
        x, y = rand(), rand()
        ok &= np.allclose(F @ mul(x, y), mul(F @ x, F @ y))
    return ok


# Conjugation rho_r as an 8x8 matrix
def rho_matrix8(r):
    rinv = conj(r) / (r @ r)
    return np.column_stack([mul(mul(r, E[j]), rinv) for j in range(8)])


# Triple rotation by theta = automorphism iff theta = 0, +-2pi/3 (r^3 real)
ok = True
for theta, expect in [(0.0, True), (0.7, False), (np.pi / 3, False), (np.pi, False),
                      (2 * np.pi / 3, True), (4 * np.pi / 3, True)]:
    u = rand_unit_im()
    r = oexp(u * theta / 2)
    F = rho_matrix8(r)
    r3 = mul(mul(r, r), r)
    ok &= is_auto(F) == expect
    ok &= np.allclose(im(r3), 0) == expect
print("rho_r automorphism iff rotation angle in {0, +-2pi/3} iff r^3 real:", ok)

# Diagonal-phase reading: conjugation by exp(e1 t/2) has phases (t, t, t)
# on the complex basis (e2, e4, e6) with i = left mult by e1, and
# 3t = 0 mod 2pi is the automorphism condition
theta = 2 * np.pi / 3
F = rho_matrix8(oexp(E[1] * theta / 2))
ok = is_auto(F)
for j in (2, 4, 6):
    target = np.cos(theta) * E[j] + np.sin(theta) * mul(E[1], E[j])
    ok &= np.allclose(F @ E[j], target)
print("phases (t,t,t) on (e2,e4,e6), automorphism at 3t = 0 mod 2pi:", ok)

# Derivations: D(xy) = D(x)y + x D(y) as linear constraints on 8x8 D
def derivation_space(fixed=()):
    rows = []
    for i in range(8):
        for j in range(8):
            A = np.zeros((8, 64))
            w = mul(E[i], E[j])
            for k in range(8):
                for m in range(8):
                    A[k, k * 8 + m] += w[m]       # D(e_i e_j)
                    A[k, m * 8 + i] -= R[j][k, m]  # -(D e_i) e_j
                    A[k, m * 8 + j] -= L[i][k, m]  # -e_i (D e_j)
            rows.append(A)
    for f in fixed:
        A = np.zeros((8, 64))
        for k in range(8):
            A[k, k * 8 + f] = 1.0  # D e_f = 0
        rows.append(A)
    A = np.vstack(rows)
    _, s, Vt = np.linalg.svd(A)
    null = Vt[np.sum(s > 1e-8):]
    return [v.reshape(8, 8) for v in null]

DER = derivation_space()
print("dim Der(O) = 14 (= dim G2):", len(DER) == 14)
ok = all(np.allclose(D @ ONE, 0) and np.allclose(D, -D.T) for D in DER)
print("derivations kill 1 and are antisymmetric (in so(7)):", ok)

# Der(O) is closed under the matrix commutator (an honest Lie algebra)
flat = np.array([D.flatten() for D in DER])
ok = True
for _ in range(20):
    D1 = np.tensordot(rng.standard_normal(14), DER, 1)
    D2 = np.tensordot(rng.standard_normal(14), DER, 1)
    br = (D1 @ D2 - D2 @ D1).flatten()
    coef, *_ = np.linalg.lstsq(flat.T, br, rcond=None)
    ok &= np.allclose(flat.T @ coef, br)
print("[Der, Der] stays in Der (g2 closed under bracket):", ok)

# ad_u fails Leibniz by associators: [u,zw] - [u,z]w - z[u,w] = -3(u,z,w)
ok = True
for _ in range(20):
    u, z, w = rand(), rand(), rand()
    lhs = comm(u, mul(z, w)) - mul(comm(u, z), w) - mul(z, comm(u, w))
    ok &= np.allclose(lhs, -3 * assoc(u, z, w))
print("[u,zw] - [u,z]w - z[u,w] = -3(u,z,w):", ok)

# Corrected map D_{x,y}(z) = [[x,y],z] - 3(x,y,z) is a derivation
def D_xy(x, y):
    return np.column_stack([comm(comm(x, y), E[j]) - 3 * assoc(x, y, E[j])
                            for j in range(8)])

ok = True
for _ in range(10):
    x, y = rand(), rand()
    D = D_xy(x, y)
    for _ in range(5):
        z, w = rand(), rand()
        ok &= np.allclose(D @ mul(z, w), mul(D @ z, w) + mul(z, D @ w))
print("D_{x,y}(z) = [[x,y],z] - 3(x,y,z) satisfies the Leibniz rule:", ok)

mats = [D_xy(E[i], E[j]).flatten() for i, j in combinations(range(1, 8), 2)]
rank = np.linalg.matrix_rank(np.array(mats))
print("the maps D_{x,y} span all of Der(O): rank =", rank, "=", 14, ":", rank == 14)

# so(7) = Der(O) (+) cross maps: 21 = 14 + 7, direct sum
CROSS = []
for i in range(1, 8):
    M = np.zeros((8, 8))
    for j in range(1, 8):  # x -> e_i cross x on Im O (column 0 stays zero)
        M[:, j] = cross(E[i], E[j])
    CROSS.append(M)
ok = all(np.allclose(M, -M.T) and np.allclose(M @ ONE, 0) for M in CROSS)
print("cross maps kill 1 and are antisymmetric (in so(7)):", ok)
both = np.array([D.flatten() for D in DER] + [M.flatten() for M in CROSS])
print("so(7) = g2 (+) {u cross .}: dim =", np.linalg.matrix_rank(both),
      "= 21, direct sum:", np.linalg.matrix_rank(both) == 21)
crossflat = np.array([M.flatten() for M in CROSS])
ok = True
for M in CROSS:  # cross maps are not derivations (except 0)
    coef, *_ = np.linalg.lstsq(flat.T, M.flatten(), rcond=None)
    ok &= not np.allclose(flat.T @ coef, M.flatten())
print("no nonzero cross map is a derivation:", ok)

# exp(D) is an automorphism (in G2, orthogonal, det 1, preserves cross);
# exp of a generic so(7) element is a rotation but not an automorphism
D = np.tensordot(rng.standard_normal(14), DER, 1)
F = mexp(D)
ok = is_auto(F)
ok &= np.allclose(F.T @ F, np.eye(8)) and np.isclose(np.linalg.det(F), 1)
for _ in range(10):
    x, y = im(rand()), im(rand())
    ok &= np.allclose(F @ cross(x, y), cross(F @ x, F @ y))
print("exp(derivation): automorphism, orthogonal, det 1, preserves cross:", ok)
G = mexp(np.tensordot(rng.standard_normal(21), DER + CROSS, 1))
print("exp(generic so(7)): rotation but not an automorphism:",
      np.allclose(G.T @ G, np.eye(8)) and not is_auto(G))

# Basic triples: u1 in S^6, u2 in S^5 (perp u1), u3 in S^3 (perp 1,u1,u2,u1u2)
def basic_triple():
    u1 = rand_unit_im()
    v = im(rand())
    v -= (v @ u1) * u1
    u2 = v / norm(v)
    H = np.array([ONE, u1, u2, mul(u1, u2)])
    w = rand()
    for h in H:
        w -= (w @ h) * h
    u3 = w / norm(w)
    return u1, u2, u3

def triple_auto(u1, u2, u3):
    """The automorphism sending (e1, e2, e4) to (u1, u2, u3)."""
    cols = [ONE, u1, u2, mul(u1, u2), u3, mul(u1, u3), mul(u2, u3),
            mul(mul(u1, u2), u3)]
    return np.column_stack(cols)

ok = True
for _ in range(10):
    F = triple_auto(*basic_triple())
    ok &= is_auto(F)
    ok &= np.allclose(F.T @ F, np.eye(8)) and np.isclose(np.linalg.det(F), 1)
print("every basic triple defines an automorphism (transitivity), det 1:", ok)

# Freeness: an automorphism is determined by the images of e1, e2, e4
ok = True
for _ in range(10):
    D = np.tensordot(rng.standard_normal(14), DER, 1)
    F = mexp(D)
    F2 = triple_auto(F @ E[1], F @ E[2], F @ E[4])
    ok &= np.allclose(F, F2)
print("an automorphism is determined by the images of e1, e2, e4:", ok)

# Stabilizer chain: dims 14 -> 8 -> 3 -> 0, orbit dims 6, 5, 3
d0 = len(DER)
d1 = len(derivation_space(fixed=(1,)))
d2 = len(derivation_space(fixed=(1, 2)))
d3 = len(derivation_space(fixed=(1, 2, 4)))
print("stabilizer dims 14, 8, 3, 0:", (d0, d1, d2, d3) == (14, 8, 3, 0),
      "/ orbit dims (S^6, S^5, S^3):", (d0 - d1, d1 - d2, d2 - d3) == (6, 5, 3))

# Fixing e1, e2 also fixes e3 = e1 e2 (whole quaternion subalgebra pointwise),
# and the 3-dim stabilizer algebra acts on span(e4..e7) with su(2) structure
DER12 = derivation_space(fixed=(1, 2))
ok = all(np.allclose(D @ E[3], 0) for D in DER12)
X = [D for D in DER12]
G = np.array([[np.trace(A.T @ B) for B in X] for A in X])
Ch = np.linalg.cholesky(np.linalg.inv(G))  # orthonormalize
Y = [sum(Ch[i, j] * X[j] for j in range(3)) for i in range(3)]
brs = [(Y[0] @ Y[1] - Y[1] @ Y[0]), (Y[1] @ Y[2] - Y[2] @ Y[1]),
       (Y[2] @ Y[0] - Y[0] @ Y[2])]
rank = np.linalg.matrix_rank(np.array([b.flatten() for b in brs]))
print("fixing e1, e2 fixes e3; 3-dim algebra with full bracket (su(2)):",
      ok and rank == 3)

# Derivations fixing e1 = su(3): on C^3 = span(e2, e4, e6) with i = L_{e1},
# the matrices are anti-Hermitian and traceless; exp is unitary with det 1
DER1 = derivation_space(fixed=(1,))
IDX = {2: 0, 4: 1, 6: 2}
IM = {2: 3, 4: 5, 6: 7}   # e1 e2 = e3, e1 e4 = e5, e1 e6 = -e7
SGN = {2: 1.0, 4: 1.0, 6: -1.0}

def cmat(D):
    M = np.zeros((3, 3), dtype=complex)
    for j, cj in IDX.items():
        v = D @ E[j]
        for k, ck in IDX.items():
            M[ck, cj] = v[k] + 1j * SGN[k] * v[IM[k]]
    return M

ok = True
for D in DER1:
    ok &= np.allclose(D @ E[1], 0)
    ok &= np.allclose(D @ L[1], L[1] @ D)  # commutes with the complex structure
    M = cmat(D)
    ok &= np.allclose(M.conj().T, -M) and np.isclose(np.trace(M), 0)
print("derivations fixing e1: C-linear, anti-Hermitian, traceless (su(3)):", ok)
Ms = np.array([cmat(D).flatten() for D in DER1])
Mre = np.vstack([Ms.real, Ms.imag])
print("their C^3 matrices span all of su(3): rank =", np.linalg.matrix_rank(Ms), "= 8:",
      np.linalg.matrix_rank(Ms) == 8)
D = np.tensordot(rng.standard_normal(8), DER1, 1)
M = cmat(mexp(D))  # exp acts C-linearly; read off its C^3 matrix
print("exp of such a derivation acts on C^3 as SU(3):",
      np.allclose(M.conj().T @ M, np.eye(3)) and np.isclose(np.linalg.det(M), 1))
