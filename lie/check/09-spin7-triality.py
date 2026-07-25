"""Numerical checks for 09-spin7-triality.md.

Left multiplications L_i by the imaginary octonion units as 8x8 signed
permutation matrices: Clifford relations L_i^2 = -I, L_iL_j = -L_jL_i,
non-closure L_1L_2 != L_3, and the volume relations L_1...L_7 = -I,
L_7 = L_1...L_6 (Cl_{0,6} = M_8(R)); the antisymmetric split
so(8) = span{L_i} (+) span{L_iL_j} (28 = 7 + 21) and the bracket
closure of span{L_u} being all of so(8) (unit left multiplications
generate SO(8)); spin(7) = span{L_iL_j} closed under bracket, the
vector action [A, L_v] = L_{a v} with a in so(7), the double cover
(exp(t L_1L_2 / 2) rotates the (e_1, e_2) plane by t, turns the
spinor side by t/2 isoclinically, and gives -1 at t = 2pi); the
compatibility g(vz) = chi_g(v) g(z); the spinor-1 stabilizer
{A in spin(7) : A(1) = 0} = Der(O) = g_2 (dim 14, orbit dim 7) and
exp of it being an automorphism; spin(6) = {A in spin(7) :
[A, L_7] = 0} (dim 15) commuting with the complex structure J = L_7
(anti-Hermitian traceless on C^4: su(4)), its spinor-1 stabilizer
being su(3) (dim 8, fixes e_7, inside g_2); the doubled construction
gamma(v) (16x16) with gamma(v)^2 = |v|^2 I generating Cl_{8,0} =
M_16(R) (256-dim); spin(8) = span{gamma_i gamma_j} (dim 28) as
block-diagonal pairs (a_+, a_-), the vector action a_0 defined by
a_+ L_y - L_y a_- = L_{a_0 y}, the infinitesimal triality
a_+(yz) = (a_0 y)z + y(a_- z) and its group form
g_+(yz) = (g_0 y)(g_- z) (trivial for quaternions by associativity);
the three projections a -> a_+, a_-, a_0 being linear bijections
spin(8) -> so(8); triality being outer (a rank-2 plane rotation on
the vector side pulls back to isoclinic half-angle rotations
a_+^2 = -I/4 on the half-spinor sides); the three stabilizers of
"1" (vector, spinor S^+, spinor S^-) each of dim 21, with
{a_0(1) = 0} = {a_+ = a_-} = spin(7), all pairwise intersections
equal to the triple intersection g_2 (dim 14 = 21 + 21 - 28) whose
common action is by derivations; and the stabilizer chain
28 -> 21 -> 14 -> 8 -> 3 -> 0 with orbit dims 7, 7, 6, 5, 3
(spheres S^7, S^7, S^6, S^5, S^3), fixing the two "1"s, then the
vector units e_7, e_1, e_2 in turn.
"""

from itertools import combinations

import numpy as np

rng = np.random.default_rng(9)


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
I8 = np.eye(8)


def Lmat(v):
    """Left multiplication by a general octonion v (L is linear in v)."""
    return sum(v[i] * L[i] for i in range(8))


def is_auto(F, trials=20):
    ok = True
    for _ in range(trials):
        x, y = rand(), rand()
        ok &= np.allclose(F @ mul(x, y), mul(F @ x, F @ y))
    return ok


def span_rank(mats):
    return np.linalg.matrix_rank(np.array([M.flatten() for M in mats]))


def null_space(A, tol=1e-8):
    _, s, Vt = np.linalg.svd(A)
    return Vt[np.sum(s > tol):]


# --- Left multiplications as Clifford generators ---

ok = True
for i in range(1, 8):
    ok &= np.allclose(L[i] @ L[i], -I8) and np.allclose(L[i], -L[i].T)
    for j in range(1, 8):
        if i != j:
            ok &= np.allclose(L[i] @ L[j], -(L[j] @ L[i]))
print("L_i^2 = -I, antisymmetric, L_iL_j = -L_jL_i (i != j):", ok)

print("e_1e_2 = e_3 but L_1L_2 != L_3 (left actions do not close):",
      np.allclose(mul(E[1], E[2]), E[3]) and not np.allclose(L[1] @ L[2], L[3]))

vol7 = np.eye(8)
for i in range(1, 8):
    vol7 = vol7 @ L[i]
vol6 = np.eye(8)
for i in range(1, 7):
    vol6 = vol6 @ L[i]
print("L_1...L_7 = -I and L_7 = L_1...L_6 (volume element of Cl_{0,6}):",
      np.allclose(vol7, -I8) and np.allclose(vol6, L[7]))

# --- so(8) = span{L_i} (+) span{L_iL_j}, and left multiplications generate SO(8) ---

grade1 = [L[i] for i in range(1, 8)]
grade2 = [L[i] @ L[j] for i, j in combinations(range(1, 8), 2)]
print("dim span{L_i} = 7, dim span{L_iL_j} = 21, together 28 = dim so(8):",
      (span_rank(grade1), span_rank(grade2), span_rank(grade1 + grade2))
      == (7, 21, 28))

mats = [M.flatten() for M in grade1]
r_prev = 0
r = np.linalg.matrix_rank(np.array(mats))
while r > r_prev:
    r_prev = r
    base = [M.reshape(8, 8) for M in mats]
    for A in base:
        for B in base:
            mats.append((A @ B - B @ A).flatten())
    r = np.linalg.matrix_rank(np.array(mats))
print("bracket closure of span{L_u} has dim", r, "= 28 (generates so(8)):",
      r == 28)

# --- spin(7) = span{L_iL_j} and the double cover onto SO(7) ---

SPIN7 = grade2
ok = True
for _ in range(20):
    A = np.tensordot(rng.standard_normal(21), SPIN7, 1)
    B = np.tensordot(rng.standard_normal(21), SPIN7, 1)
    br = (A @ B - B @ A).flatten()
    flat21 = np.array([M.flatten() for M in SPIN7])
    coef, *_ = np.linalg.lstsq(flat21.T, br, rcond=None)
    ok &= np.allclose(flat21.T @ coef, br)
print("spin(7) = span{L_iL_j} is closed under the bracket (dim 21):", ok)

# vector action: [A, L_v] = L_{a v} with a in so(7)
ok = True
for _ in range(10):
    A = np.tensordot(rng.standard_normal(21), SPIN7, 1)
    a = np.column_stack([(A @ L[j] - L[j] @ A) @ ONE for j in range(8)])
    for _ in range(5):
        v = im(rand())
        M = A @ Lmat(v) - Lmat(v) @ A
        ok &= np.allclose(M, Lmat(a @ v))
    ok &= np.allclose(a, -a.T) and np.allclose(a @ ONE, 0)
print("[A, L_v] = L_{a v} with a antisymmetric fixing 1 (a in so(7)):", ok)

# double cover: g(t) = exp(t L_1L_2 / 2) rotates the (e_1, e_2) plane by t,
# while the spinor side turns by the half angle t/2 (isoclinic), and
# g(2pi) = -1
t = 0.8
g = mexp(t / 2 * (L[1] @ L[2]))
chi = np.column_stack([(g @ L[j] @ g.T) @ ONE for j in range(8)])
target = np.eye(8)
target[1, 1] = target[2, 2] = np.cos(t)
target[2, 1] = np.sin(t)
target[1, 2] = -np.sin(t)
ok = np.allclose(chi, target)
ok &= np.allclose(g @ ONE, np.cos(t / 2) * ONE + np.sin(t / 2) * E[3])
ok &= np.allclose(g + g.T, 2 * np.cos(t / 2) * I8)
ok &= np.allclose(mexp(2 * np.pi / 2 * (L[1] @ L[2])), -I8)
print("exp(t L_1L_2/2): vector angle t, spinor angle t/2 (1 -> e_3),",
      "isoclinic, -1 at t = 2pi:", ok)

# compatibility: g(v z) = chi_g(v) g(z) for g in Spin(7), v imaginary
ok = True
for _ in range(10):
    g = mexp(np.tensordot(rng.standard_normal(21), SPIN7, 1))
    chi = np.column_stack([(g @ L[j] @ g.T) @ ONE for j in range(8)])
    ok &= np.allclose(chi @ ONE, ONE)
    ok &= np.allclose(chi.T @ chi, I8) and np.isclose(np.linalg.det(chi), 1)
    for _ in range(5):
        v, z = im(rand()), rand()
        ok &= np.allclose(g @ mul(v, z), mul(chi @ v, g @ z))
print("g(vz) = chi_g(v) g(z), chi_g in SO(7) (fixes 1, orthogonal, det 1):", ok)


# --- the spinor-1 stabilizer in spin(7) is g_2 = Der(O) ---

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
    return [v.reshape(8, 8) for v in null_space(A)]


DER = derivation_space()
ev1 = np.array([(M @ ONE) for M in SPIN7])  # 21 x 8: A -> A(1)
stab = [np.tensordot(c, SPIN7, 1) for c in null_space(ev1.T)]
print("dim {A in spin(7) : A(1) = 0} =", len(stab), "= 14:", len(stab) == 14)
print("orbit dim {A(1)} =", np.linalg.matrix_rank(ev1), "= 7 (tangent of S^7):",
      np.linalg.matrix_rank(ev1) == 7)

both = np.array([M.flatten() for M in stab] + [D.flatten() for D in DER])
print("{A in spin(7) : A(1) = 0} = Der(O) = g_2 (span equality):",
      np.linalg.matrix_rank(both) == 14 and len(DER) == 14)

g = mexp(np.tensordot(rng.standard_normal(14), stab, 1))
print("exp of a spinor-1 stabilizer element is an automorphism (in G_2):",
      is_auto(g) and np.allclose(g @ ONE, ONE))

# --- spin(6) and su(4), su(3) ---

comm7 = np.array([((M @ L[7] - L[7] @ M)).flatten() for M in SPIN7])  # 21 x 64
spin6_coef = null_space(comm7.T)
SPIN6 = [np.tensordot(c, SPIN7, 1) for c in spin6_coef]
target6 = [L[i] @ L[j] for i, j in combinations(range(1, 7), 2)]
both6 = np.array([M.flatten() for M in SPIN6] + [M.flatten() for M in target6])
print("{A in spin(7) : [A, L_7] = 0} = span{L_iL_j (i,j <= 6)} = spin(6), dim",
      len(SPIN6), "= 15:", len(SPIN6) == 15
      and np.linalg.matrix_rank(both6) == 15)

J = L[7]
ok = np.allclose(J @ J, -I8)
for M in SPIN6:
    ok &= np.allclose(M @ J, J @ M)
    ok &= np.isclose(np.trace(M), 0) and np.isclose(np.trace(J @ M), 0)
print("J = L_7: J^2 = -I; spin(6) commutes with J, real and J-traces vanish",
      "(anti-Hermitian traceless on C^4 = su(4)):", ok)

ev6 = np.array([(M @ ONE) for M in SPIN6])
su3 = [np.tensordot(c, SPIN6, 1) for c in null_space(ev6.T)]
ok = len(su3) == 8
flat_der = np.array([D.flatten() for D in DER])
for M in su3:
    coef, *_ = np.linalg.lstsq(flat_der.T, M.flatten(), rcond=None)
    ok &= np.allclose(flat_der.T @ coef, M.flatten())
    ok &= np.allclose(M @ E[7], 0)
print("{A in spin(6) : A(1) = 0}: dim 8, inside g_2, fixes e_7 (su(3)):", ok)

orbit7 = np.array([np.column_stack([(M @ L[j] - L[j] @ M) @ ONE
                                    for j in range(8)]) @ E[7] for M in SPIN7])
print("orbit dim of the vector e_7 under spin(7) =",
      np.linalg.matrix_rank(orbit7), "= 6 (tangent of S^6):",
      np.linalg.matrix_rank(orbit7) == 6)

# --- the doubled construction: Cl_{8,0} and Spin(8) ---

GAM = []
for i in range(8):
    gmat = np.zeros((16, 16))
    gmat[:8, 8:] = Lmat(E[i])
    gmat[8:, :8] = Lmat(conj(E[i]))
    GAM.append(gmat)

I16 = np.eye(16)
ok = True
for _ in range(10):
    v = rand()
    gv = np.zeros((16, 16))
    gv[:8, 8:] = Lmat(v)
    gv[8:, :8] = Lmat(conj(v))
    ok &= np.allclose(gv @ gv, (v @ v) * I16)
print("gamma(v)^2 = |v|^2 I (via v(v* z) = |v|^2 z, two elements only):", ok)

prods = {0: [np.eye(16)]}
allp = [np.eye(16)]
for mask in range(1, 256):
    idx = [i for i in range(8) if mask & (1 << i)]
    M = np.eye(16)
    for i in idx:
        M = M @ GAM[i]
    allp.append(M)
print("the 256 products of gamma_0..gamma_7 span M_16(R) (Cl_{8,0}):",
      span_rank(allp) == 256)

SPIN8 = [GAM[i] @ GAM[j] for i, j in combinations(range(8), 2)]
ok = span_rank(SPIN8) == 28
for M in SPIN8:
    ok &= np.allclose(M[:8, 8:], 0) and np.allclose(M[8:, :8], 0)
print("spin(8) = span{gamma_i gamma_j}: dim 28, block-diagonal (a_+, a_-):", ok)


def parts(a):
    return a[:8, :8], a[8:, 8:]


def vec_action(ap, am):
    return np.column_stack([(ap @ Lmat(E[j]) - Lmat(E[j]) @ am) @ ONE
                            for j in range(8)])


ok = True
for _ in range(10):
    a = np.tensordot(rng.standard_normal(28), SPIN8, 1)
    ap, am = parts(a)
    a0 = vec_action(ap, am)
    ok &= np.allclose(a0, -a0.T)
    for _ in range(5):
        y = rand()
        ok &= np.allclose(ap @ Lmat(y) - Lmat(y) @ am, Lmat(a0 @ y))
        z = rand()
        ok &= np.allclose(ap @ mul(y, z), mul(a0 @ y, z) + mul(y, am @ z))
print("a_+ L_y - L_y a_- = L_{a_0 y}, a_0 in so(8);",
      "a_+(yz) = (a_0 y)z + y(a_- z):", ok)

# group-level triality: g_+(yz) = (g_0 y)(g_- z), g_0 in SO(8), sign-blind
ok = True
for _ in range(5):
    a = np.tensordot(rng.standard_normal(28), SPIN8, 1)
    ap, am = parts(a)
    gp, gm = mexp(ap), mexp(am)
    g0 = np.column_stack([gp @ mul(E[j], gm.T @ ONE) for j in range(8)])
    ok &= np.allclose(g0.T @ g0, I8) and np.isclose(np.linalg.det(g0), 1)
    for _ in range(5):
        y, z = rand(), rand()
        ok &= np.allclose(gp @ mul(y, z), mul(g0 @ y, gm @ z))
    g0_flip = np.column_stack([(-gp) @ mul(E[j], (-gm).T @ ONE)
                               for j in range(8)])
    ok &= np.allclose(g0, g0_flip)
print("g_+(yz) = (g_0 y)(g_- z) with g_0 in SO(8); (-g_+, -g_-) -> same g_0:",
      ok)

# for quaternions the same identity is trivial: p(yz) = (p y q^-1)(q z)
ok = True
for _ in range(10):
    p, q = rand(4), rand(4)
    p, q = p / norm(p), q / norm(q)
    y, z = rand(4), rand(4)
    ok &= np.allclose(mul(p, mul(y, z)), mul(mul(mul(p, y), conj(q)), mul(q, z)))
print("quaternions: p(yz) = (p y q^{-1})(q z) holds by associativity:", ok)

# the three projections a -> a_+, a_-, a_0 are bijections spin(8) -> so(8)
plus_flat = np.array([parts(M)[0].flatten() for M in SPIN8])
minus_flat = np.array([parts(M)[1].flatten() for M in SPIN8])
zero_flat = np.array([vec_action(*parts(M)).flatten() for M in SPIN8])
ranks = (np.linalg.matrix_rank(plus_flat), np.linalg.matrix_rank(minus_flat),
         np.linalg.matrix_rank(zero_flat))
print("projections a -> a_+, a_-, a_0 all have rank 28 (bijections onto so(8)):",
      ranks == (28, 28, 28))

# triality is outer: a plane rotation (rank 2) on the vector side pulls
# back to isoclinic half-angle rotations on the spinor sides
T = np.zeros((8, 8))
T[1, 0], T[0, 1] = 1.0, -1.0  # rotation generator in the (1, e_1) plane
coef, *_ = np.linalg.lstsq(zero_flat.T, T.flatten(), rcond=None)
a = np.tensordot(coef, SPIN8, 1)
ap, am = parts(a)
ok = np.allclose(vec_action(ap, am), T)
ok &= np.linalg.matrix_rank(T) == 2
ok &= np.allclose(ap @ ap, -I8 / 4) and np.allclose(am @ am, -I8 / 4)
print("vector plane rotation (rank 2) <-> isoclinic a_+^2 = a_-^2 = -I/4",
      "(rank 8): spectra differ, so triality is outer:", ok)

# the three stabilizers of "1" and their intersections
ev_plus = np.array([parts(M)[0] @ ONE for M in SPIN8])   # 28 x 8
ev_minus = np.array([parts(M)[1] @ ONE for M in SPIN8])
ev_zero = np.array([vec_action(*parts(M)) @ ONE for M in SPIN8])
d_plus = len(null_space(ev_plus.T))
d_minus = len(null_space(ev_minus.T))
d_zero = len(null_space(ev_zero.T))
print("stabilizers of 1 in V, S^+, S^-: dims", (d_zero, d_plus, d_minus),
      "= (21, 21, 21):", (d_zero, d_plus, d_minus) == (21, 21, 21))
print("orbit dims (transitive on the three unit spheres S^7):",
      (np.linalg.matrix_rank(ev_zero), np.linalg.matrix_rank(ev_plus),
       np.linalg.matrix_rank(ev_minus)) == (7, 7, 7))

# {a_0(1) = 0} = {a_+ = a_-}: the vector stabilizer is our diagonal Spin(7)
diff_flat = np.array([(parts(M)[0] - parts(M)[1]).flatten() for M in SPIN8])
vz = null_space(ev_zero.T)
ok = len(vz) == 21
for c in vz:
    ap, am = parts(np.tensordot(c, SPIN8, 1))
    ok &= np.allclose(ap, am)
print("{a in spin(8) : a_0(1) = 0} = {a_+ = a_-} (= spin(7), acting as",
      "L_iL_j on both halves):", ok)

pair_zp = len(null_space(np.vstack([ev_zero.T, ev_plus.T])))
pair_zm = len(null_space(np.vstack([ev_zero.T, ev_minus.T])))
pair_pm = len(null_space(np.vstack([ev_plus.T, ev_minus.T])))
triple = null_space(np.vstack([ev_zero.T, ev_plus.T, ev_minus.T]))
print("pairwise intersections of the three stabilizers: dims",
      (pair_zp, pair_zm, pair_pm), "= (14, 14, 14), triple:", len(triple),
      "= 14 (all equal g_2 = 21 + 21 - 28):",
      (pair_zp, pair_zm, pair_pm, len(triple)) == (14, 14, 14, 14))

ok = True
for c in triple:
    ap, am = parts(np.tensordot(c, SPIN8, 1))
    a0 = vec_action(ap, am)
    ok &= np.allclose(ap, am) and np.allclose(ap, a0)
    coef, *_ = np.linalg.lstsq(flat_der.T, ap.flatten(), rcond=None)
    ok &= np.allclose(flat_der.T @ coef, ap.flatten())
print("on the triple intersection a_+ = a_- = a_0 is a derivation (g_2):", ok)

# stabilizer chain 28 -> 21 -> 14 -> 8 (orbits S^7, S^7, S^6), then 08's
# 14 -> 8 -> 3 -> 0
ev_e7 = np.array([vec_action(*parts(M)) @ E[7] for M in SPIN8])
chain8 = len(null_space(np.vstack([ev_zero.T, ev_plus.T, ev_e7.T])))
print("adding the vector e_7 to the fixed set: dim", chain8, "= 8 (su(3));",
      "chain 28 -> 21 -> 14 -> 8 with orbits 7, 7, 6:",
      (28 - d_zero, d_zero - pair_zp, pair_zp - chain8) == (7, 7, 6))

# the chain continues inside g_2 by fixing e_7, then e_1, then e_2
c1 = len(derivation_space(fixed=(7,)))
c2 = len(derivation_space(fixed=(7, 1)))
c3 = len(derivation_space(fixed=(7, 1, 2)))
print("derivations fixing e_7 / e_7, e_1 / e_7, e_1, e_2: dims", (c1, c2, c3),
      "= (8, 3, 0), orbits 6, 5, 3 (S^6, S^5, S^3):",
      (c1, c2, c3) == (8, 3, 0))
