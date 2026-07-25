"""Numerical checks for 10-roots-dynkin.md.

Weights as rotation numbers: su(2) ladder [s3, s+-] = +-2 s+-, spin-j
irreps having integer weights -2j..2j in steps of 2 (quantization from
the closed torus circle exp(2 pi i H) = I), and the SO(3) condition
exp(i pi H) = I holding exactly for integer spin (even weights);
weights of the su(3) defining representation forming an equilateral
triangle with centroid 0; roots as weights of the adjoint
representation: the six su(3) roots (+-2,0), +-(1,sqrt3), +-(-1,sqrt3)
forming a regular hexagon with 6 + 2 = 8 = dim, E_12 = l1 + i l2 as a
ladder shifting weights by the root mu1 - mu2 = (2,0), and the braided
su(2) triples of 05 matching root pairs with coroot coefficients; the
integrality condition 2(b,a)/(a,a) in Z for roots and weights and the
resulting angle/length restrictions; rank-2 root systems: so(4) =
(+-1,+-1) splitting into two orthogonal A1 pairs (D2 = A1 x A1),
so(5) = B2 (8 roots, length ratio sqrt2, 45-degree grid), su(3) = A2,
and g2 = Der(O) (12 roots, ratio sqrt3, 30-degree grid) whose 6 long
roots + Cartan lie in the su(3) subalgebra fixing e1 (14 = 8 + 6);
simple roots and Cartan matrices computed from each algebra and
identified with A2, B2 = C2 (so(5) vs sp(2)), A3 = D3 (su(4) vs
so(6)), B3 (so(7)), D4 (so(8)) with its S3 diagram symmetry, and G2;
Weyl reflections in the simple roots regenerating the full root
systems; dim = #roots + rank for every algebra; and the exceptional
root systems from explicit coordinates: F4 (48 roots) and E8 (240
roots, all length sqrt2) with their Cartan matrices/diagram shapes,
and E7 (126) and E6 (72) cut out of E8 as the roots orthogonal to a
root / to an A2 pair, giving dim = #roots + rank = 52, 78, 133, 248.
"""

from itertools import combinations, permutations, product

import numpy as np

rng = np.random.default_rng(10)

TOL = 1e-8


def comm(A, B):
    return A @ B - B @ A


def expm(A, terms=200):
    S, T = np.eye(*A.shape, dtype=complex), np.eye(*A.shape, dtype=complex)
    for n in range(1, terms):
        T = T @ A / n
        S = S + T
    return S


# ---------------------------------------------------------------- su(2)

s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]], dtype=complex)
sp = (s1 + 1j * s2) / 2
sm = (s1 - 1j * s2) / 2
print("[s3, s+] = 2 s+, [s3, s-] = -2 s-:",
      np.allclose(comm(s3, sp), 2 * sp) and np.allclose(comm(s3, sm), -2 * sm))


def spin_rep(j):
    """Standard spin-j matrices; H has the integer weights 2m."""
    d = int(round(2 * j)) + 1
    m = j - np.arange(d)
    H = np.diag(2 * m).astype(complex)
    E = np.zeros((d, d), dtype=complex)
    for k in range(d - 1):
        E[k, k + 1] = np.sqrt(j * (j + 1) - m[k + 1] * (m[k + 1] + 1))
    F = E.conj().T
    return H, E, F


ok_lad = ok_int = ok_so3 = True
for j in (0.5, 1, 1.5, 2):
    H, E, F = spin_rep(j)
    ok_lad &= np.allclose(comm(H, E), 2 * E) and np.allclose(comm(E, F), H)
    w = np.diag(H).real
    ok_int &= np.allclose(w, np.round(w)) and np.allclose(np.sort(w), np.arange(-2 * j, 2 * j + 1, 2))
    # torus closes after one turn: exp(i theta H) at theta = 2 pi is I
    ok_int &= np.allclose(expm(1j * 2 * np.pi * H), np.eye(len(w)))
    # SO(3) circle closes at theta = pi: only for even weights (integer spin)
    is_id = np.allclose(expm(1j * np.pi * H), np.eye(len(w)))
    ok_so3 &= is_id == (float(j).is_integer())
print("spin-j ladder: [H,E] = 2E, [E,F] = H; weights = -2j..2j step 2:", ok_lad and ok_int)
print("exp(i pi H) = I exactly for integer spin (SO(3) reps = even weights):", ok_so3)

# ---------------------------------------------------------------- su(3)

l1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
l2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]])
l3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
l4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
l5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]])
l6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
l7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]])
l8 = np.diag([1, 1, -2]).astype(complex) / np.sqrt(3)
lam = [l1, l2, l3, l4, l5, l6, l7, l8]

# weights of the defining representation C^3, measured by (l3, l8)
mu = np.array([[l[k, k].real for l in (l3, l8)] for k in range(3)])
print("weights of C^3:", mu.round(3).tolist())
d01 = np.linalg.norm(mu[0] - mu[1])
print("equilateral triangle, centroid 0:",
      np.allclose([np.linalg.norm(mu[i] - mu[j]) for i, j in combinations(range(3), 2)], d01)
      and np.allclose(mu.sum(axis=0), 0))

# roots: eigenvalues of ad_{l3}, ad_{l8} on the complexified algebra
E12, E45, E67 = l1 + 1j * l2, l4 + 1j * l5, l6 + 1j * l7
su3_roots = []
ok = True
for E in (E12, E45, E67):
    a3 = comm(l3, E)[np.abs(E) > 0.5][0] / E[np.abs(E) > 0.5][0]
    a8 = comm(l8, E)[np.abs(E) > 0.5][0] / E[np.abs(E) > 0.5][0]
    ok &= np.allclose(comm(l3, E), a3 * E) and np.allclose(comm(l8, E), a8 * E)
    su3_roots += [np.array([a3.real, a8.real]), -np.array([a3.real, a8.real])]
su3_roots = np.array(su3_roots)
print("E_12, E_45, E_67 are simultaneous eigenvectors of ad_l3, ad_l8:", ok)
print("six roots:", su3_roots.round(3).tolist())
expected = np.array([[2, 0], [-2, 0], [1, np.sqrt(3)], [-1, -np.sqrt(3)],
                     [-1, np.sqrt(3)], [1, -np.sqrt(3)]])
match = all(min(np.linalg.norm(su3_roots - e, axis=1)) < TOL for e in expected)
print("roots = (+-2,0), +-(1,sqrt3), +-(-1,sqrt3):", match)
lens = np.linalg.norm(su3_roots, axis=1)
angs = sorted(round(np.degrees(np.arctan2(r[1], r[0]))) % 360 for r in su3_roots)
print("regular hexagon (all length 2, every 60 degrees):",
      np.allclose(lens, 2) and angs == [0, 60, 120, 180, 240, 300])
print("6 roots + rank 2 = dim su(3) = 8:", len(su3_roots) + 2 == 8)

# E_12 = 2 e_12 shifts weights by the root: mu1 - mu2 = (2, 0)
e12 = np.zeros((3, 3), dtype=complex)
e12[0, 1] = 1
print("E_12 = 2 e_12 (matrix unit), root = mu1 - mu2 = (2,0):",
      np.allclose(E12, 2 * e12)
      and np.allclose(mu[0] - mu[1], [2, 0]))
print("roots = weight differences of C^3:",
      all(min(np.linalg.norm(su3_roots - (mu[i] - mu[j]), axis=1)) < TOL
          for i in range(3) for j in range(3) if i != j))

# braided su(2) triples of 05: third element = (a3 l3 + a8 l8)/2 with (a3, a8) = root
ok = True
for (X, Y), r in (((l1, l2), (2, 0)), ((l4, l5), (1, np.sqrt(3))),
                  ((l6, l7), (-1, np.sqrt(3)))):
    ok &= np.allclose(comm(X, Y), 1j * (r[0] * l3 + r[1] * l8))  # [s1,s2] = 2i s3 pattern
print("braid triples: [Re E, Im E] = i (a3 l3 + a8 l8), coroot = root coords:", ok)

# integrality 2(b,a)/(a,a) in Z for all roots and weights (Euclidean coords)
def n_int(b, a):
    return 2 * (b @ a) / (a @ a)

ok = all(abs(n_int(b, a) - round(n_int(b, a))) < TOL for a in su3_roots for b in su3_roots)
ok &= all(abs(n_int(m, a) - round(n_int(m, a))) < TOL for a in su3_roots for m in mu)
print("2(b,a)/(a,a) integer for all root/root and weight/root pairs:", ok)

# ------------------------------------------------- generic root machinery


def coords(X, basis):
    M = np.stack([B.ravel() for B in basis]).T
    c, *_ = np.linalg.lstsq(M, X.ravel(), rcond=None)
    assert np.allclose(M @ c, X.ravel()), "element not in span of basis"
    return c


def ad_matrix(h, basis):
    return np.stack([coords(comm(h, B), basis) for B in basis]).T


def root_data(basis, cartan):
    """Roots of a compact algebra: imaginary parts of the simultaneous
    eigenvalues of ad_h for h in the Cartan, plus a Gram matrix for the
    root space from the (negative-definite) Killing form."""
    ads = [ad_matrix(h, basis) for h in cartan]
    M = sum(c * A for c, A in zip(rng.standard_normal(len(ads)), ads))
    vals, vecs = np.linalg.eig(M.astype(complex))
    roots, vectors = [], []
    for k in range(len(vals)):
        if abs(vals[k]) < 1e-6:
            continue
        v = vecs[:, k]
        r = [(v.conj() @ (A @ v)) / (v.conj() @ v) for A in ads]
        assert all(np.allclose(A @ v, rc * v, atol=1e-6) for A, rc in zip(ads, r))
        roots.append(np.array([rc.imag for rc in r]))
        vectors.append(v)
    K = np.array([[np.trace(A @ B).real for B in ads] for A in ads])
    G = np.linalg.inv(-K)  # Killing form is negative definite on a compact algebra
    return np.array(roots), vectors, G


def ip(a, b, G):
    return a @ G @ b


def simple_roots(roots, tol=1e-6):
    f = np.arange(1, len(roots[0]) + 1) * np.pi ** np.arange(len(roots[0]))  # generic
    pos = [r for r in roots if f @ r > 0]
    key = {tuple(np.round(r, 5)) for r in pos}
    simple = [r for r in pos
              if not any(tuple(np.round(r - a, 5)) in key for a in pos)]
    return pos, simple


def cartan_matrix(simple, G):
    n = len(simple)
    C = np.array([[2 * ip(simple[i], simple[j], G) / ip(simple[j], simple[j], G)
                   for j in range(n)] for i in range(n)])
    assert np.allclose(C, np.round(C), atol=1e-6)
    return np.round(C).astype(int)


def same_diagram(C1, C2):
    if C1.shape != C2.shape:
        return False
    n = len(C1)
    return any(np.array_equal(C1[np.ix_(p, p)], C2) for p in permutations(range(n)))


def weyl_closure(simple, G, tol=1e-4):
    roots = [np.array(r) for r in simple]
    queue = list(roots)
    while queue:
        r = queue.pop()
        for a in simple:
            s = r - 2 * ip(r, a, G) / ip(a, a, G) * a
            if all(np.linalg.norm(s - x) > tol for x in roots):
                roots.append(s)
                queue.append(s)
    return roots


def diagram_arms(C):
    """(degree list, arm lengths from the unique branch node) of a tree diagram."""
    n = len(C)
    adj = [[j for j in range(n) if j != i and C[i][j] != 0] for i in range(n)]
    deg = sorted(len(a) for a in adj)
    branch = [i for i in range(n) if len(adj[i]) == 3]
    arms = []
    for b in branch:
        for start in adj[b]:
            length, prev, cur = 1, b, start
            while True:
                nxt = [j for j in adj[cur] if j != prev]
                if not nxt:
                    break
                prev, cur = cur, nxt[0]
                length += 1
            arms.append(length)
    return deg, sorted(arms)


# expected Cartan matrices
def C_A(n):
    C = 2 * np.eye(n, dtype=int)
    for i in range(n - 1):
        C[i, i + 1] = C[i + 1, i] = -1
    return C

C_B2 = np.array([[2, -2], [-1, 2]])
C_B3 = np.array([[2, -1, 0], [-1, 2, -2], [0, -1, 2]])
C_D4 = np.array([[2, -1, -1, -1], [-1, 2, 0, 0], [-1, 0, 2, 0], [-1, 0, 0, 2]])
C_G2 = np.array([[2, -1], [-3, 2]])
C_F4 = np.array([[2, -1, 0, 0], [-1, 2, -2, 0], [0, -1, 2, -1], [0, 0, -1, 2]])

# --------------------------------------------- classical algebras as matrices


def so_basis(n):
    basis = []
    for i, j in combinations(range(n), 2):
        B = np.zeros((n, n), dtype=complex)
        B[i, j], B[j, i] = -1, 1
        basis.append(B)
    return basis


def so_cartan(n):
    hs = []
    for k in range(n // 2):
        h = np.zeros((n, n), dtype=complex)
        h[2 * k, 2 * k + 1], h[2 * k + 1, 2 * k] = -1, 1
        hs.append(h)
    return hs


def su_basis(n):
    basis = []
    for i, j in combinations(range(n), 2):
        B = np.zeros((n, n), dtype=complex)
        B[i, j], B[j, i] = 1, -1
        basis.append(B)
        B = np.zeros((n, n), dtype=complex)
        B[i, j] = B[j, i] = 1j
        basis.append(B)
    for k in range(n - 1):
        d = np.zeros(n)
        d[k], d[k + 1] = 1, -1
        basis.append(1j * np.diag(d).astype(complex))
    return basis


def su_cartan(n):
    return [1j * np.diag(np.eye(n)[k] - np.eye(n)[k + 1]).astype(complex)
            for k in range(n - 1)]


def sp2_basis():
    """Compact sp(2): anti-Hermitian 4x4 X with X^T O + O X = 0."""
    ah = []
    for i, j in combinations(range(4), 2):
        B = np.zeros((4, 4), dtype=complex)
        B[i, j], B[j, i] = 1, -1
        ah.append(B)
        B = np.zeros((4, 4), dtype=complex)
        B[i, j] = B[j, i] = 1j
        ah.append(B)
    for k in range(4):
        ah.append(1j * np.diag(np.eye(4)[k]).astype(complex))
    O = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
    rows = np.stack([np.concatenate([(B.T @ O + O @ B).real.ravel(),
                                     (B.T @ O + O @ B).imag.ravel()]) for B in ah])
    _, s, Vt = np.linalg.svd(rows.T)
    null = Vt[np.sum(s > TOL):]
    return [sum(c * B for c, B in zip(v, ah)) for v in null]


def check_algebra(name, basis, cartan, n_roots, expect_C=None, arms=None):
    roots, vecs, G = root_data(basis, cartan)
    pos, simple = simple_roots(roots)
    C = cartan_matrix(simple, G)
    okC = same_diagram(C, expect_C) if expect_C is not None else True
    if arms is not None:
        okC = diagram_arms(C) == arms
    closure = weyl_closure(simple, G)
    ok_weyl = len(closure) == len(roots) and all(
        min(np.linalg.norm(roots - np.array(r), axis=1)) < 1e-5 for r in closure)
    ok = (len(roots) == n_roots and len(basis) == n_roots + len(cartan)
          and len(simple) == len(cartan) and okC and ok_weyl)
    print(f"{name}: {len(roots)} roots, dim {len(basis)} = {n_roots} + rank "
          f"{len(cartan)}, diagram ok, Weyl closure ok:", ok)
    return roots, vecs, G, C


# so(4) = D2 = A1 x A1: roots (+-1,+-1) in two orthogonal pairs
roots4, _, G4, C4 = check_algebra("so(4)", so_basis(4), so_cartan(4), 4,
                                  expect_C=np.array([[2, 0], [0, 2]]))
print("so(4) roots (+-1,+-1), orthogonal pairs (D2 = A1 x A1):",
      all(min(np.linalg.norm(roots4 - np.array(e), axis=1)) < 1e-6
          for e in product([1, -1], repeat=2))
      and abs(ip(np.array([1, 1]), np.array([1, -1]), G4)) < TOL)

# so(5) = B2: 8 roots, two lengths ratio sqrt2, 45-degree fan
roots5, _, G5, _ = check_algebra("so(5)", so_basis(5), so_cartan(5), 8, expect_C=C_B2)
L = sorted(np.sqrt(ip(r, r, G5)) for r in roots5)
print("so(5) = B2: length ratio sqrt2:", np.isclose(L[-1] / L[0], np.sqrt(2)))

# sp(2) = C2: same diagram as B2 (Spin(5) = Sp(2))
sp2 = sp2_basis()
h_sp = [np.diag([1j, 0, -1j, 0]), np.diag([0, 1j, 0, -1j])]
print("sp(2) basis dim 10, Cartan inside:", len(sp2) == 10
      and all(np.allclose(sum(c * B for c, B in zip(coords(h, sp2), sp2)), h) for h in h_sp))
_, _, _, C_sp = check_algebra("sp(2)", sp2, h_sp, 8, expect_C=C_B2)

# su(4) and so(6): both A3 (Spin(6) = SU(4))
_, _, _, Csu4 = check_algebra("su(4)", su_basis(4), su_cartan(4), 12, expect_C=C_A(3))
_, _, _, Cso6 = check_algebra("so(6)", so_basis(6), so_cartan(6), 12, expect_C=C_A(3))
print("A3 = D3: su(4) and so(6) give the same diagram:", same_diagram(Csu4, Cso6))

# su(3) through the same machinery: A2
check_algebra("su(3)", su_basis(3), su_cartan(3), 6, expect_C=C_A(2))

# so(7) = B3, so(8) = D4
check_algebra("so(7)", so_basis(7), so_cartan(7), 18, expect_C=C_B3)
_, _, _, C8 = check_algebra("so(8)", so_basis(8), so_cartan(8), 24, expect_C=C_D4)

# D4 diagram S3 symmetry (triality): permuting the three outer nodes fixes C
syms = [p for p in permutations(range(4)) if np.array_equal(C_D4[np.ix_(p, p)], C_D4)]
print("D4 Cartan matrix: symmetry group of the diagram = S3 (6 permutations):",
      len(syms) == 6)

# --------------------------------------------------------- g2 from Der(O)


def conj_o(a):
    c = -a.copy()
    c[0] = a[0]
    return c


def mul_o(a, b):
    n = len(a)
    if n == 1:
        return a * b
    h = n // 2
    a1, a2, b1, b2 = a[:h], a[h:], b[:h], b[h:]
    return np.concatenate([
        mul_o(a1, b1) - mul_o(conj_o(b2), a2),
        mul_o(b2, a1) + mul_o(a2, conj_o(b1)),
    ])


E8b = [np.eye(8)[i] for i in range(8)]
Lo = [np.stack([mul_o(E8b[i], E8b[k]) for k in range(8)]).T for i in range(8)]
Ro = [np.stack([mul_o(E8b[k], E8b[i]) for k in range(8)]).T for i in range(8)]


def derivation_space(fixed=()):
    rows = []
    for i in range(8):
        for j in range(8):
            A = np.zeros((8, 64))
            w = mul_o(E8b[i], E8b[j])
            for k in range(8):
                for m in range(8):
                    A[k, k * 8 + m] += w[m]
                    A[k, m * 8 + i] -= Ro[j][k, m]
                    A[k, m * 8 + j] -= Lo[i][k, m]
            rows.append(A)
    for f in fixed:
        A = np.zeros((8, 64))
        for k in range(8):
            A[k, k * 8 + f] = 1.0
        rows.append(A)
    A = np.vstack(rows)
    _, s, Vt = np.linalg.svd(A)
    return [v.reshape(8, 8) for v in Vt[np.sum(s > TOL):]]


DER = [D.astype(complex) for D in derivation_space()]
STAB = [D.astype(complex) for D in derivation_space(fixed=(1,))]  # su(3) fixing e1
print("dim Der(O) = 14, stabilizer of e1 = 8 (su(3)):", len(DER) == 14 and len(STAB) == 8)

# Cartan of g2 inside the su(3): centralizer of a generic stabilizer element
X = sum(c * D for c, D in zip(rng.standard_normal(8), STAB))
adX = np.stack([coords(comm(X, D), DER) for D in DER]).T
_, s, Vt = np.linalg.svd(adX)
cart = [sum(c * D for c, D in zip(v, DER)) for v in Vt[np.sum(s > 1e-6):]]
print("centralizer of a generic element is 2-dimensional (rank g2 = 2):", len(cart) == 2)

g2_roots, g2_vecs, Gg, Cg2 = check_algebra("g2 = Der(O)", DER, cart, 12, expect_C=C_G2)
Lg = np.array(sorted(np.sqrt(ip(r, r, Gg)) for r in g2_roots))
print("g2: 6 short + 6 long, length ratio sqrt3:",
      np.allclose(Lg[:6], Lg[0]) and np.allclose(Lg[6:], Lg[0] * np.sqrt(3)))
cosang = sorted({round(ip(a, b, Gg) / np.sqrt(ip(a, a, Gg) * ip(b, b, Gg)), 6)
                 for a in g2_roots for b in g2_roots})
allowed = sorted({round(np.cos(np.radians(30 * k)), 6) for k in range(13)})
print("g2: all angles multiples of 30 degrees:", all(c in allowed for c in cosang))

# long roots of g2: their root vectors lie in the complexified su(3) = stab(e1)
Mstab = np.stack([D.ravel() for D in STAB]).T
long_ok = short_ok = True
for r, v in zip(g2_roots, g2_vecs):
    Xv = sum(c * D for c, D in zip(v, DER)).ravel()
    c, *_ = np.linalg.lstsq(Mstab, Xv, rcond=None)
    inside = np.allclose(Mstab @ c, Xv, atol=1e-6)
    if np.isclose(np.sqrt(ip(r, r, Gg)), Lg[-1]):
        long_ok &= inside
    else:
        short_ok &= not inside
print("6 long roots span the su(3) fixing e1 (with the Cartan: 8 = 2 + 6),"
      " short roots outside:", long_ok and short_ok)
long6 = np.array([r for r in g2_roots if np.isclose(np.sqrt(ip(r, r, Gg)), Lg[-1])])
cos6 = sorted(round(ip(a, b, Gg) / Lg[-1] ** 2, 6) for a, b in combinations(long6, 2))
print("long roots form a regular hexagon (A2):",
      cos6 == sorted([0.5] * 6 + [-0.5] * 6 + [-1.0] * 3))

# ----------------------------------- exceptional root systems by coordinates


def euclid_checks(name, roots, rank, arms=None, expect_C=None):
    roots = np.array(roots)
    G = np.eye(len(roots[0]))
    okZ = all(abs(n_int(b, a) - round(n_int(b, a))) < TOL for a in roots for b in roots)
    span = np.linalg.matrix_rank(roots)
    # project onto the span to run the simple-root machinery in full rank
    _, _, Vt = np.linalg.svd(roots)
    P = Vt[:span]
    pr = roots @ P.T
    pos, simple = simple_roots(pr)
    C = cartan_matrix(simple, np.eye(span))
    ok = okZ and span == rank and len(simple) == rank
    if expect_C is not None:
        ok &= same_diagram(C, expect_C)
    if arms is not None:
        ok &= diagram_arms(C) == arms
    closure = weyl_closure(simple, np.eye(span))
    ok &= len(closure) == len(roots)
    print(f"{name}: {len(roots)} roots, rank {rank}, integrality, diagram, "
          f"Weyl closure -> dim = {len(roots)} + {rank} = {len(roots) + rank}:", ok)
    return roots


e = np.eye(8)
E8_roots = []
for i, j in combinations(range(8), 2):
    for si, sj in product([1, -1], repeat=2):
        E8_roots.append(si * e[i] + sj * e[j])
for signs in product([0.5, -0.5], repeat=8):
    if sum(s < 0 for s in signs) % 2 == 0:
        E8_roots.append(np.array(signs))
E8_roots = np.array(E8_roots)
print("E8: 112 + 128 = 240 roots, all length sqrt2:",
      len(E8_roots) == 240 and np.allclose(np.linalg.norm(E8_roots, axis=1), np.sqrt(2)))
euclid_checks("E8", E8_roots, 8, arms=([1, 1, 1, 2, 2, 2, 2, 3], [1, 2, 4]))

a0 = e[0] + e[1]
E7_roots = [r for r in E8_roots if abs(r @ a0) < TOL]
euclid_checks("E7 = E8 roots orthogonal to one root", E7_roots, 7,
              arms=([1, 1, 1, 2, 2, 2, 3], [1, 2, 3]))

b0 = e[2] - e[1]
print("(a0, b0) form an A2 pair (angle 120):",
      np.isclose(a0 @ b0, -1) and np.isclose(b0 @ b0, 2))
E6_roots = [r for r in E8_roots if abs(r @ a0) < TOL and abs(r @ b0) < TOL]
euclid_checks("E6 = E8 roots orthogonal to an A2 pair", E6_roots, 6,
              arms=([1, 1, 1, 2, 2, 3], [1, 2, 2]))

e4 = np.eye(4)
F4_roots = []
for i, j in combinations(range(4), 2):
    for si, sj in product([1, -1], repeat=2):
        F4_roots.append(si * e4[i] + sj * e4[j])
for i in range(4):
    F4_roots += [e4[i], -e4[i]]
for signs in product([0.5, -0.5], repeat=4):
    F4_roots.append(np.array(signs))
print("F4: 24 + 8 + 16 = 48 roots:", len(F4_roots) == 48)
euclid_checks("F4", F4_roots, 4, expect_C=C_F4)

print("dim = #roots + rank: G2 14, F4 52, E6 78, E7 133, E8 248:",
      [12 + 2, 48 + 4, 72 + 6, 126 + 7, 240 + 8] == [14, 52, 78, 133, 248]
      and len(E7_roots) == 126 and len(E6_roots) == 72)
