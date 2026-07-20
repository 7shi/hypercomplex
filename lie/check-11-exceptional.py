"""Numerical checks for 11-exceptional.md.

Hermitian matrices over the octonions and the Jordan product
A o B = (AB + BA)/2: closure, commutativity, non-associativity;
the Jordan identity (A o B) o A^2 = A o (B o A^2) holding for
h_2(O) and h_3(O) but failing for h_4(O) (while h_4(H) is fine);
power associativity; the trace form tr(X o Y) and its associativity
tr((X o Y) o Z) = tr(X o (Y o Z)); Cayley-Hamilton for the cubic
Jordan algebra h_3(O) and the explicit determinant formula
det = xi1 xi2 xi3 + 2 Re(x1 x2 x3) - sum xi_i |x_i|^2;
derivations of h_3(O) (Leibniz rule) giving dim Der = 52 = f4,
killing the identity, preserving the trace, antisymmetric,
closed under bracket, with exp(D) an automorphism preserving
trace and det; rank-1 idempotents: the h_2(O) block P = v v*
always idempotent (Artin) reproducing the Hopf map S^15 -> S^8,
while a generic v in O^3 fails (OP^2 needs associative coordinates);
the stabilizer of E1 in f4 having dim 36 = so(9) (acting on a
9-dim vector block and a 16-dim spinor block), orbit dim 16 = OP^2;
the stabilizer of all diagonal idempotents having dim 28 = so(8),
preserving the three octonion slots with antisymmetric blocks,
each of the three block projections injective (triality) and the
three blocks tied by the invariance of Re((x1 x2) x3);
compact e6 = f4 (+) i L_{X} (X traceless): L_X symmetric,
[L_X, L_Y] a derivation with such brackets spanning all of f4,
[D, L_X] = L_{DX} (Leibniz), real dim 78, closed under bracket,
anti-Hermitian (exp unitary = compact) and det-preserving on the
complexified J, while exp(D + L_X) preserves the real det but not
the norm (the noncompact real form E6(-26)); the n = 2 row:
dim Der(h_2(K)) = 1, 3, 10, 36 = so(2), so(3), so(5), so(9) for
K = R, C, H, O; the Freudenthal magic square via the Tits formula
dim M(A, B) = dim Der(A) + dim Der(h_3(B)) + dim Im A * dim h_3(B)_0
with every entry computed from numerically obtained dimensions,
matching the known table and symmetric; E7 and E8 dimension
bookkeeping 133 = 78 + 27 + 27 + 1, 248 = 120 + 128 =
28 + 28 + 3*64 = 133 + 3 + 2*56.
"""

from itertools import combinations

import numpy as np

rng = np.random.default_rng(11)


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


def rand(n=8):
    return rng.standard_normal(n)


def mexp(A):
    """Matrix exponential by scaling and squaring with a Taylor series."""
    k = max(0, int(np.ceil(np.log2(max(np.abs(A).sum(), 1e-16)))) + 1)
    B = A / 2**k
    E, term = np.eye(len(A), dtype=A.dtype), np.eye(len(A), dtype=A.dtype)
    for n in range(1, 30):
        term = term @ B / n
        E = E + term
    for _ in range(k):
        E = E @ E
    return E


# --- matrices over a Cayley-Dickson algebra K (shape (n, n, d)) ---

def mmul(A, B):
    n, d = A.shape[0], A.shape[2]
    C = np.zeros((n, n, d), dtype=np.result_type(A, B))
    for i in range(n):
        for k in range(n):
            for j in range(n):
                C[i, k] += mul(A[i, j], B[j, k])
    return C


def hconj(A):
    n = A.shape[0]
    return np.stack([[conj(A[j, i]) for j in range(n)] for i in range(n)])


def jordan(A, B):
    return (mmul(A, B) + mmul(B, A)) / 2


def tr(A):
    return sum(A[i, i, 0] for i in range(A.shape[0]))


def herm_basis(n, d):
    """Hermitian basis, orthonormal for the trace form tr(X o Y):
    diagonal units first, then off-diagonal slots (i, j), i < j,
    scaled by 1/sqrt(2)."""
    bas = []
    for i in range(n):
        E = np.zeros((n, n, d))
        E[i, i, 0] = 1.0
        bas.append(E)
    s = 1 / np.sqrt(2)
    for i in range(n):
        for j in range(i + 1, n):
            for u in range(d):
                E = np.zeros((n, n, d))
                e = basis(d, u)
                E[i, j] = s * e
                E[j, i] = s * conj(e)
                bas.append(E)
    return bas


def rand_herm(n, d):
    bas = herm_basis(n, d)
    c = rng.standard_normal(len(bas))
    return sum(ci * B for ci, B in zip(c, bas))


# --- Jordan algebra closure, commutativity, Jordan identity ---

def jordan_identity_ok(n, d, trials=10):
    ok = True
    for _ in range(trials):
        A, B = rand_herm(n, d), rand_herm(n, d)
        A2 = jordan(A, A)
        ok &= np.allclose(jordan(jordan(A, B), A2), jordan(A, jordan(B, A2)))
    return ok


ok = True
for _ in range(10):
    A, B, C = (rand_herm(3, 8) for _ in range(3))
    ok &= np.allclose(jordan(A, B), hconj(jordan(A, B)))  # closure
    ok &= np.allclose(jordan(A, B), jordan(B, A))         # commutative
print("h3(O): A o B Hermitian and commutative:", ok)
A, B, C = (rand_herm(3, 8) for _ in range(3))
print("h3(O): o is not associative:",
      not np.allclose(jordan(jordan(A, B), C), jordan(A, jordan(B, C))))

print("Jordan identity holds in h2(O):", jordan_identity_ok(2, 8))
print("Jordan identity holds in h3(O):", jordan_identity_ok(3, 8))
print("Jordan identity FAILS in h4(O):", not jordan_identity_ok(4, 8))
print("Jordan identity holds in h4(H) (associative base):",
      jordan_identity_ok(4, 4))

# Power associativity in h3(O): A^2 o A^3 = A o (A o A^3) (= A^5 both ways)
A = rand_herm(3, 8)
A2 = jordan(A, A)
A3 = jordan(A, A2)
print("h3(O) power associativity (A^2 o A^3 = A o (A o A^3)):",
      np.allclose(jordan(A2, A3), jordan(A, jordan(A, A3))))

# Trace form associativity: tr((X o Y) o Z) = tr(X o (Y o Z))
ok = True
for _ in range(10):
    X, Y, Z = (rand_herm(3, 8) for _ in range(3))
    ok &= np.isclose(tr(jordan(jordan(X, Y), Z)), tr(jordan(X, jordan(Y, Z))))
print("trace form associativity tr((X o Y) o Z) = tr(X o (Y o Z)):", ok)

# Re((xy)z) = Re(x(yz)) for octonions: the associator is purely imaginary,
# so Re(xyz) needs no parentheses
ok = True
for _ in range(10):
    x, y, z = rand(), rand(), rand()
    ok &= np.isclose(mul(mul(x, y), z)[0], mul(x, mul(y, z))[0])
print("Re((xy)z) = Re(x(yz)) (Re(xyz) needs no parentheses):", ok)

# --- determinant of h3(O): Newton formula and explicit cubic form ---

def newton_det_traces(t1, t2, t3):
    return (t1**3 - 3 * t1 * t2 + 2 * t3) / 6


def ndet(A):
    A2 = jordan(A, A)
    A3 = jordan(A, A2)
    return newton_det_traces(tr(A), tr(A2), tr(A3))


def albert(xi, x):
    """[[xi1, x3, x2*], [x3*, xi2, x1], [x2, x1*, xi3]]"""
    A = np.zeros((3, 3, 8))
    for i in range(3):
        A[i, i, 0] = xi[i]
    A[1, 2], A[2, 1] = x[0], conj(x[0])
    A[2, 0], A[0, 2] = x[1], conj(x[1])
    A[0, 1], A[1, 0] = x[2], conj(x[2])
    return A


xi = rng.standard_normal(3)
x = [rand() for _ in range(3)]
A = albert(xi, x)
det1 = ndet(A)
det2 = (xi[0] * xi[1] * xi[2] + 2 * mul(mul(x[0], x[1]), x[2])[0]
        - sum(xi[i] * (x[i] @ x[i]) for i in range(3)))
print("explicit det = xi1 xi2 xi3 + 2Re(x1x2x3) - sum xi_i |x_i|^2:",
      np.isclose(det1, det2))

# Cayley-Hamilton: A^3 - t1 A^2 + sigma A - det I = 0
A = rand_herm(3, 8)
A2 = jordan(A, A)
A3 = jordan(A, A2)
t1, t2 = tr(A), tr(A2)
sigma = (t1**2 - t2) / 2
I3 = sum(herm_basis(3, 8)[:3])
print("Cayley-Hamilton A^3 - t1 A^2 + sigma A - det I = 0:",
      np.allclose(A3 - t1 * A2 + sigma * A - ndet(A) * I3, 0))

# --- structure tensor and derivations (generic Leibniz solver) ---

def jordan_structure(n, d):
    """S[i, j] = coordinates of B_i o B_j in the orthonormal basis."""
    bas = herm_basis(n, d)
    m = len(bas)
    flat = np.array([b.reshape(-1) for b in bas])
    S = np.zeros((m, m, m))
    for i in range(m):
        for j in range(i, m):
            c = flat @ jordan(bas[i], bas[j]).reshape(-1)
            S[i, j] = c
            S[j, i] = c
    return S


def derivation_space(S, fixed=(), sym=True, tol=1e-8):
    """Null space of the Leibniz constraints D(x_i x_j) = D(x_i)x_j + x_i D(x_j),
    plus D c = 0 for each coordinate vector c in fixed."""
    m = S.shape[0]
    idx = np.arange(m)
    rows = []
    for i in range(m):
        for j in range(i if sym else 0, m):
            blk = np.zeros((m, m, m))
            blk[idx, idx, :] = S[i, j]
            blk[:, :, i] -= S[:, j, :].T
            blk[:, :, j] -= S[i, :, :].T
            rows.append(blk.reshape(m, m * m))
    for c in fixed:
        blk = np.zeros((m, m, m))
        blk[idx, idx, :] = c
        rows.append(blk.reshape(m, m * m))
    Amat = np.vstack(rows)
    _, s, Vt = np.linalg.svd(Amat, full_matrices=True)
    null = Vt[np.sum(s > tol * s[0]):]
    return [v.reshape(m, m) for v in null]


S3 = jordan_structure(3, 8)
m3 = S3.shape[0]
print("dim h3(O) = 27 (= 3 + 3*8):", m3 == 27)

DER = derivation_space(S3)
print("dim Der(h3(O)) = 52 (= dim F4):", len(DER) == 52)

idc = np.zeros(m3)
idc[:3] = 1.0  # coordinates of the identity matrix
trv = idc      # tr(X) = <I, X> in the orthonormal basis
ok = all(np.allclose(D @ idc, 0) and np.allclose(trv @ D, 0)
         and np.allclose(D, -D.T) for D in DER)
print("derivations kill I, preserve trace, antisymmetric:", ok)

flat52 = np.array([D.reshape(-1) for D in DER])


def in_span(flat, v):
    coef, *_ = np.linalg.lstsq(flat.T, v, rcond=None)
    return np.allclose(flat.T @ coef, v)


ok = True
for _ in range(10):
    D1 = np.tensordot(rng.standard_normal(52), DER, 1)
    D2 = np.tensordot(rng.standard_normal(52), DER, 1)
    ok &= in_span(flat52, (D1 @ D2 - D2 @ D1).reshape(-1))
print("[Der, Der] stays in Der (f4 closed under bracket):", ok)


def jmul(xc, yc):
    """Jordan product in coordinates."""
    return np.einsum("ijk,i,j->k", S3, xc, yc)


def jdet(xc):
    x2 = jmul(xc, xc)
    x3 = jmul(xc, x2)
    return newton_det_traces(trv @ xc, trv @ x2, trv @ x3)


# exp(D) is an automorphism of the Jordan product, preserves trace and det
D = np.tensordot(rng.standard_normal(52), DER, 1)
F = mexp(D)
ok = True
for _ in range(10):
    xc, yc = rng.standard_normal(m3), rng.standard_normal(m3)
    ok &= np.allclose(F @ jmul(xc, yc), jmul(F @ xc, F @ yc))
    ok &= np.isclose(trv @ (F @ xc), trv @ xc)
    ok &= np.isclose(jdet(F @ xc), jdet(xc))
print("exp(derivation): Jordan automorphism, preserves trace and det:", ok)

# --- rank-1 idempotents: Hopf map and OP^2 ---

# h2(O): P = v v* is always idempotent (2 elements associate: Artin),
# reproducing the Hopf map S^15 -> S^8
o1, o2 = rand(), rand()
nrm = np.sqrt(o1 @ o1 + o2 @ o2)
o1, o2 = o1 / nrm, o2 / nrm
P = np.zeros((2, 2, 8))
P[0, 0] = (o1 @ o1) * basis(8, 0)
P[1, 1] = (o2 @ o2) * basis(8, 0)
P[0, 1] = mul(o1, conj(o2))
P[1, 0] = conj(P[0, 1])
print("h2(O): P = vv* (v in S^15) is a trace-1 idempotent (Hopf S^15 -> S^8):",
      np.allclose(jordan(P, P), P) and np.isclose(tr(P), 1))

# traceless part of h2(O): X o Y = (X . Y) I, the symmetrized Clifford relation
HB2 = herm_basis(2, 8)
I2 = HB2[0] + HB2[1]
ok = True
for _ in range(10):
    c1, c2 = rng.standard_normal(len(HB2)), rng.standard_normal(len(HB2))
    c1[1], c2[1] = -c1[0], -c2[0]  # traceless
    X = sum(c * B for c, B in zip(c1, HB2))
    Y = sum(c * B for c, B in zip(c2, HB2))
    XY = jordan(X, Y)
    ok &= np.allclose(XY, (tr(XY) / 2) * I2)
print("h2(O) traceless: X o Y = (X . Y) I (symmetrized Clifford relation):", ok)

# h3(O): generic v in O^3 does NOT give an idempotent, but v with
# components in a quaternion subalgebra does (OP^2 needs associativity)
def vvdag(vs):
    v = np.array(vs)
    v = v / np.sqrt(sum(vi @ vi for vi in v))
    P = np.zeros((3, 3, 8))
    for i in range(3):
        for j in range(3):
            P[i, j] = mul(v[i], conj(v[j]))
    return P


P = vvdag([rand(), rand(), rand()])
print("h3(O): generic v in O^3: vv* fails to be idempotent:",
      not np.allclose(jordan(P, P), P))
q = [np.concatenate([rand(4), np.zeros(4)]) for _ in range(3)]  # in H
P = vvdag(q)
print("h3(O): v with components in H: vv* is a trace-1 idempotent:",
      np.allclose(jordan(P, P), P) and np.isclose(tr(P), 1))

# --- stabilizers in f4: Spin(9) and the OP^2 orbit, Spin(8) and triality ---

e1c = np.zeros(m3)
e1c[0] = 1.0  # E1 = diag(1, 0, 0)
DER_E1 = derivation_space(S3, fixed=(e1c,))
print("stabilizer of E1 in f4: dim 36 (= so(9)), orbit dim 16 (= OP^2):",
      len(DER_E1) == 36 and 52 - len(DER_E1) == 16)

# slot index ranges in the orthonormal basis: diag 0..2, then
# (0,1): 3..10, (0,2): 11..18, (1,2): 19..26
SLOT = {(0, 1): np.arange(3, 11), (0, 2): np.arange(11, 19),
        (1, 2): np.arange(19, 27)}
# 9-dim vector block: span{E2 - E3, slot (1,2)}; build projector basis
V9vecs = [np.zeros(m3) for _ in range(9)]
V9vecs[0][1], V9vecs[0][2] = 1 / np.sqrt(2), -1 / np.sqrt(2)
for k, i in enumerate(SLOT[(1, 2)]):
    V9vecs[k + 1][i] = 1.0
S16vecs = [np.zeros(m3) for _ in range(16)]
for k, i in enumerate(np.concatenate([SLOT[(0, 1)], SLOT[(0, 2)]])):
    S16vecs[k][i] = 1.0
V9m = np.array(V9vecs)
S16m = np.array(S16vecs)
ok = True
for D in DER_E1:
    for v in V9vecs:  # D maps the 9-dim block to itself
        w = D @ v
        ok &= np.allclose(w, V9m.T @ (V9m @ w))
    for v in S16vecs:  # and the 16-dim block to itself
        w = D @ v
        ok &= np.allclose(w, S16m.T @ (S16m @ w))
print("stabilizer of E1 preserves the 9-dim and 16-dim blocks:", ok)
blocks = np.array([(V9m @ D @ V9m.T).reshape(-1) for D in DER_E1])
print("its action on the 9-dim block spans so(9): rank =",
      np.linalg.matrix_rank(blocks), "= 36:",
      np.linalg.matrix_rank(blocks) == 36)

# fixing all diagonal idempotents: dim 28 = so(8), slots preserved
e2c = np.zeros(m3)
e2c[1] = 1.0
DER0 = derivation_space(S3, fixed=(e1c, e2c))
print("stabilizer of E1, E2 (hence E3): dim 28 (= so(8)):", len(DER0) == 28)
ok = True
for D in DER0:
    for (i, j), sl in SLOT.items():
        for k in sl:
            v = np.zeros(m3)
            v[k] = 1.0
            w = D @ v
            keep = np.zeros(m3)
            keep[sl] = w[sl]
            ok &= np.allclose(w, keep)
print("it preserves each of the three octonion slots:", ok)


def slot_block(D, key):
    sl = SLOT[key]
    return D[np.ix_(sl, sl)]


ok = all(np.allclose(slot_block(D, k), -slot_block(D, k).T)
         for D in DER0 for k in SLOT)
print("the three 8x8 blocks are antisymmetric (three so(8) actions):", ok)
for key in SLOT:
    r = np.linalg.matrix_rank(np.array([slot_block(D, key).reshape(-1)
                                        for D in DER0]))
    assert r == 28
print("each block projection is injective (rank 28): one so(8) action",
      "determines the other two (triality):", True)

# the three blocks are tied by the invariance of Re((x1 x2) x3):
# with octonion actions B1 (slot (1,2)), B2 = K B(0,2) K, B3 (slot (0,1))
K8 = np.diag([1.0] + [-1.0] * 7)
ok = True
for _ in range(5):
    D = np.tensordot(rng.standard_normal(28), DER0, 1)
    B1 = slot_block(D, (1, 2))
    B2 = K8 @ slot_block(D, (0, 2)) @ K8
    B3 = slot_block(D, (0, 1))
    for _ in range(5):
        u, v, w = rand(), rand(), rand()
        val = (mul(mul(B1 @ u, v), w)[0] + mul(mul(u, B2 @ v), w)[0]
               + mul(mul(u, v), B3 @ w)[0])
        ok &= np.isclose(val, 0)
print("Re((B1 u)v w) + Re(u (B2 v) w) + Re(u v (B3 w)) = 0 (triality tie):", ok)

# the S3 permuting the diagonal entries (and the three slots) is realized
# inside F4: X -> P X P^T for a permutation matrix P is a Jordan automorphism
HB3 = herm_basis(3, 8)


def perm_op(p):
    """Coordinate matrix of X -> (X[p[i], p[j]])_{ij} in the orthonormal basis."""
    flat = np.array([b.reshape(-1) for b in HB3])
    cols = []
    for B in HB3:
        PB = np.stack([[B[p[i], p[j]] for j in range(3)] for i in range(3)])
        cols.append(flat @ PB.reshape(-1))
    return np.column_stack(cols)


ok = True
for p in ([1, 2, 0], [1, 0, 2]):  # a 3-cycle and a transposition
    F = perm_op(p)
    for _ in range(5):
        xc, yc = rng.standard_normal(m3), rng.standard_normal(m3)
        ok &= np.allclose(F @ jmul(xc, yc), jmul(F @ xc, F @ yc))
print("diagonal permutations (S3) are Jordan automorphisms (inside F4):", ok)


# G2 sits inside F4: applying a derivation of O entrywise is a derivation of
# h3(O) fixing the diagonal, and its three octonion blocks all equal D
# (cf. 09: G2 is where the three so(8) actions coincide)
def cd_structure(d):
    S = np.zeros((d, d, d))
    for i in range(d):
        for j in range(d):
            S[i, j] = mul(basis(d, i), basis(d, j))
    return S


DER_O = derivation_space(cd_structure(8), sym=False)
flat28 = np.array([D.reshape(-1) for D in DER0])
flat3 = np.array([b.reshape(-1) for b in HB3])
ok = len(DER_O) == 14
for _ in range(5):
    D8 = np.tensordot(rng.standard_normal(14), DER_O, 1)
    cols = []
    for B in HB3:
        DB = np.stack([[D8 @ B[i, j] for j in range(3)] for i in range(3)])
        cols.append(flat3 @ DB.reshape(-1))
    M = np.column_stack(cols)
    ok &= in_span(flat28, M.reshape(-1))
    ok &= np.allclose(slot_block(M, (1, 2)), D8)
    ok &= np.allclose(K8 @ slot_block(M, (0, 2)) @ K8, D8)
    ok &= np.allclose(slot_block(M, (0, 1)), D8)
print("Der(O) entrywise: a derivation of h3(O), all three blocks = D",
      "(G2 inside F4):", ok)

# --- compact e6 = f4 (+) i L_X (X traceless) on the complexified J ---

# multiplication operators L_X in the orthonormal basis are symmetric
Lb = [S3[i].T for i in range(m3)]  # L_{B_i}
ok = all(np.allclose(L, L.T) for L in Lb)
print("L_X is symmetric w.r.t. the trace form:", ok)

# traceless directions: 24 off-diagonal + 2 diagonal combinations
d1 = np.zeros(m3)
d1[0], d1[1] = 1 / np.sqrt(2), -1 / np.sqrt(2)
d2 = np.zeros(m3)
d2[0], d2[1], d2[2] = 1 / np.sqrt(6), 1 / np.sqrt(6), -2 / np.sqrt(6)
t26 = [d1, d2] + [np.eye(m3)[i] for i in range(3, m3)]


def Lop(xc):
    return np.tensordot(xc, Lb, 1)


# [L_X, L_Y] is a derivation, and such brackets span all of f4
brackets = []
ok = True
for xc in t26:
    for yc in t26[:6]:
        br = Lop(xc) @ Lop(yc) - Lop(yc) @ Lop(xc)
        ok &= in_span(flat52, br.reshape(-1))
for i, j in combinations(range(26), 2):
    br = Lop(t26[i]) @ Lop(t26[j]) - Lop(t26[j]) @ Lop(t26[i])
    brackets.append(br.reshape(-1))
rank = np.linalg.matrix_rank(np.array(brackets))
print("[L_X, L_Y] is a derivation; such brackets span f4: rank =", rank,
      "= 52:", ok and rank == 52)

# [D, L_X] = L_{DX} (the Leibniz rule in operator form)
ok = True
for _ in range(10):
    D = np.tensordot(rng.standard_normal(52), DER, 1)
    xc = rng.standard_normal(m3)
    ok &= np.allclose(D @ Lop(xc) - Lop(xc) @ D, Lop(D @ xc))
print("[D, L_X] = L_{DX}:", ok)

E6 = [D.astype(complex) for D in DER] + [1j * Lop(xc) for xc in t26]
flat78 = np.array([np.concatenate([M.reshape(-1).real, M.reshape(-1).imag])
                   for M in E6])
print("dim(f4 (+) i L) = 78 (= dim E6):",
      np.linalg.matrix_rank(flat78) == 78)
ok = True
for _ in range(10):
    M1 = np.tensordot(rng.standard_normal(78), E6, 1)
    M2 = np.tensordot(rng.standard_normal(78), E6, 1)
    br = M1 @ M2 - M2 @ M1
    ok &= in_span(flat78,
                  np.concatenate([br.reshape(-1).real, br.reshape(-1).imag]))
print("closed under bracket (a real Lie algebra, e6):", ok)
ok = all(np.allclose(M.conj().T, -M) for M in E6)
M = np.tensordot(rng.standard_normal(78), E6, 1)
U = mexp(M)
ok &= np.allclose(U.conj().T @ U, np.eye(m3))
print("generators anti-Hermitian, exp unitary (compact):", ok)
ok = True
for _ in range(5):
    zc = rng.standard_normal(m3) + 1j * rng.standard_normal(m3)
    ok &= np.isclose(jdet(U @ zc), jdet(zc))
xc = rng.standard_normal(m3)
haslop = not np.isclose(trv @ (U @ xc.astype(complex)), trv @ xc)
print("exp(e6) preserves the complexified det but not the trace:",
      ok and haslop)

# the real form D + L_X preserves the real det but not the norm (noncompact)
xt = rng.standard_normal(m3)
xt -= (trv @ xt) / 3 * idc
G = mexp(np.tensordot(rng.standard_normal(52), DER, 1) + Lop(xt))
ok = True
for _ in range(5):
    xc = rng.standard_normal(m3)
    ok &= np.isclose(jdet(G @ xc), jdet(xc))
grow = not np.isclose(np.linalg.norm(G @ xc), np.linalg.norm(xc))
print("exp(D + L_X) preserves the real det but not the norm (E6(-26)):",
      ok and grow)

# --- the n = 2 row: Der(h2(K)) = so(2), so(3), so(5), so(9) ---

dims2 = [len(derivation_space(jordan_structure(2, d))) for d in (1, 2, 4, 8)]
print("dim Der(h2(K)) for K = R, C, H, O:", dims2, "= so(2), so(3), so(5), so(9):",
      dims2 == [1, 3, 10, 36])

# --- Freudenthal magic square via the Tits dimension formula ---

der_K = [len(derivation_space(cd_structure(d), sym=False))
         for d in (1, 2, 4)] + [len(DER_O)]
print("dim Der(K) for K = R, C, H, O:", der_K, "= [0, 0, 3, 14]:",
      der_K == [0, 0, 3, 14])

der_J = [len(derivation_space(jordan_structure(3, d)))
         for d in (1, 2, 4)] + [len(DER)]
print("dim Der(h3(K)) for K = R, C, H, O:", der_J,
      "= [3, 8, 21, 52] (so(3), su(3), sp(3), f4):", der_J == [3, 8, 21, 52])

im_K = [0, 1, 3, 7]
j0 = [3 + 3 * d - 1 for d in (1, 2, 4, 8)]  # traceless part of h3(K)
print("dim h3(K)_0 = [5, 8, 14, 26]:", j0 == [5, 8, 14, 26])

magic = [[der_K[a] + der_J[b] + im_K[a] * j0[b] for b in range(4)]
         for a in range(4)]
target = [[3, 8, 21, 52], [8, 16, 35, 78], [21, 35, 66, 133],
          [52, 78, 133, 248]]
print("Tits formula magic square:")
for row in magic:
    print("   ", row)
print("matches [so3, su3, sp3, f4 / su3, su3^2, su6, e6 /",
      "sp3, su6, so12, e7 / f4, e6, e7, e8]:", magic == target)
print("the square is symmetric (the 'magic'):",
      all(magic[a][b] == magic[b][a] for a in range(4) for b in range(4)))

# --- E7 and E8 dimension bookkeeping ---

print("dim e7 = 133 = 78 + 27 + 27 + 1 (e6 (+) J (+) J' (+) R):",
      133 == 78 + 27 + 27 + 1)
print("Freudenthal triple: 56 = 1 + 27 + 27 + 1:", 56 == 1 + 27 + 27 + 1)
print("dim e8 = 248 = 120 + 128 (so(16) (+) half-spinor):",
      248 == 120 + 128)
print("dim e8 = 248 = 28 + 28 + 3*64 (triality form, cf. article 09):",
      248 == 28 + 28 + 3 * 64)
print("dim e8 = 248 = 133 + 3 + 2*56 (e7 (+) su(2) (+) 2 x 56):",
      248 == 133 + 3 + 2 * 56)
print("Tits M(O, O) = 14 + 52 + 7*26 = 248:", 14 + 52 + 7 * 26 == 248)
