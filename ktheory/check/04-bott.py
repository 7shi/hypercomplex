"""Checks for ktheory/04-bott.md (K-groups of spheres and Bott periodicity).

Conventions as in 03: E_g is glued by identifying (x, v) over D_- with (x, g(x) v) over
D_+. For a Cl_{0,k-1}(R) module W with generators acting as J_1, ..., J_{k-1}, the
clutching function on S^{k-1} is g_W(x) = x_0 + x_1 J_1 + ... + x_{k-1} J_{k-1}.

1. Rotation homotopy: diag(a, 1) R(t) diag(1, b) R(t)^{-1} (t in [0, pi/2]) goes from
   diag(a, b) to diag(ab, 1) inside GL(2, C). With a = b = z it deforms diag(z, z) (H + H)
   to diag(z^2, 1) (H^2 + 1); with a = z, b = z^{-1} it deforms diag(z, z^{-1}) to I.
   The determinant winding is preserved along the way.
2. Extension kills the bundle: if W extends to a Cl_{0,k} module (J_k added), then
   g_t(x) = cos t g_W(x) + sin t J_k is orthogonal for all t, so g_W is homotopic to the
   constant J_k. Checked for the irreducible modules of 02 (k = 1..8).
3. Restriction table (non-graded ABS): each irreducible Cl_{0,k} module restricted to
   Cl_{0,k-1} splits into irreducibles; for k-1 = 3, 7 the two irreducibles are told
   apart by the sign of the pseudoscalar omega = J_1 ... J_{k-1} (omega^2 = 1). The
   cokernel of the restriction map Z^{#irr(k)} -> Z^{#irr(k-1)} is
   Z_2, Z_2, 0, Z, 0, 0, 0, Z for k = 1..8 (KO~(S^k)). J_k anticommutes with omega when
   k-1 is odd, so it swaps the two eigenspaces (multiplicities (1, 1)).
4. Complex: Cl_{k-1}(C) modules modulo restrictions of Cl_k(C) modules give Z, 0
   alternately (K~(S^k) = Z for k even). For k = 2 the two irreducible Cl_1(C) modules
   (e_1 = +-i) give the clutching functions z and z^{-1} (H and H^{-1}).
5. The second irreducible of 2H (generators -i, -j, -k) gives the clutching function
   c -> L_{c^*}, and L_c 1 = c, L_{c^*} 1 = c^* (degree 1 and -1 on S^3).
6. Periodicity on the algebra side: Cl_{0,k+8} = Cl_{0,k} (x) Cl_{0,8} with
   e_i -> e_i (x) omega_8 (i <= k), e_{k+j} -> 1 (x) f_j acts on R^{16 a_k}, so
   a_{k+8} = 16 a_k (checked for k = 1..8 by building the matrices).
"""

import itertools

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from common.octonion import L as OL

rng = np.random.default_rng(4)


def rot(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s], [s, c]])


def winding(vals):
    ph = np.unwrap(np.angle(vals))
    return round((ph[-1] - ph[0]) / (2 * np.pi))


print("1. rotation homotopy diag(a,b) ~ diag(ab,1)")
th = np.linspace(0, 2 * np.pi, 401)
zs = np.exp(1j * th)
for name, bf in (("z, z", lambda z: z), ("z, z^-1", lambda z: 1 / z)):
    mins = []
    for t in np.linspace(0, np.pi / 2, 51):
        dets = []
        for z in zs:
            a, b = z, bf(z)
            g = np.diag([a, 1]) @ rot(t) @ np.diag([1, b]) @ rot(t).T
            mins.append(np.linalg.svd(g, compute_uv=False).min())
            dets.append(np.linalg.det(g))
        assert winding(np.array(dets)) == winding(zs * bf(zs))
    z = zs[37]
    g0 = np.diag([z, 1]) @ np.diag([1, bf(z)])
    g1 = np.diag([z, 1]) @ rot(np.pi / 2) @ np.diag([1, bf(z)]) @ rot(np.pi / 2).T
    assert np.allclose(g0, np.diag([z, bf(z)])) and np.allclose(g1, np.diag([z * bf(z), 1]))
    print(f"  a, b = {name}: endpoints OK, min singular value {min(mins):.3f} > 0,"
          f" det winding {winding(zs * bf(zs))} constant")
    assert min(mins) > 0.1


def reps():
    """Irreducible Cl_{0,k}(R) modules (k = 0..8), as in 02."""
    J = np.array([[0, -1], [1, 0]])
    quat = [OL[i][:4, :4] for i in (1, 2, 3)]
    oct_ = [OL[i] for i in range(1, 8)]
    Z = np.zeros((8, 8), dtype=int)
    I = np.eye(8, dtype=int)
    e8 = [np.block([[m, Z], [Z, -m]]) for m in oct_] + [np.block([[Z, -I], [I, Z]])]
    return {0: [], 1: [J], 2: quat[:2], 3: quat, 4: oct_[:4], 5: oct_[:5],
            6: oct_[:6], 7: oct_, 8: e8}


R = reps()
dim = {k: (g[0].shape[0] if g else 1) for k, g in R.items()}

print("2. an extended module gives a null-homotopic clutching function")
for k in range(1, 9):
    gens = R[k]
    n = dim[k]
    W, Jk = gens[:k - 1], gens[k - 1]
    for _ in range(10):
        x = rng.normal(size=k)
        x /= np.linalg.norm(x)
        gW = x[0] * np.eye(n) + sum(x[i] * W[i - 1] for i in range(1, k))
        for t in np.linspace(0, np.pi / 2, 11):
            gt = np.cos(t) * gW + np.sin(t) * Jk
            assert np.allclose(gt.T @ gt, np.eye(n))
print("  cos t g_W + sin t J_k is orthogonal (k = 1..8): OK")


def irreducibles(k):
    """List of irreducible Cl_{0,k} modules (generator lists)."""
    g = R[k]
    if k in (3, 7):
        return [g, [-m for m in g]]
    return [g]


def pseudoscalar(gens, n):
    return np.linalg.multi_dot(gens) if len(gens) > 1 else (gens[0] if gens else np.eye(n))


def restrict(gens, km1):
    """Multiplicities of the irreducible Cl_{0,km1} modules in the restriction."""
    n = gens[0].shape[0]
    sub = gens[:km1]
    if km1 in (3, 7):
        w = pseudoscalar(sub, n)
        assert np.allclose(w @ w, np.eye(n))
        ref = pseudoscalar(irreducibles(km1)[0], dim[km1])[0, 0]  # sign on the first irr.
        ev = np.round(np.linalg.eigvalsh((w + w.T) / 2)).astype(int)
        plus = int((ev == ref).sum()) // dim[km1]
        minus = int((ev == -ref).sum()) // dim[km1]
        # J_k anticommutes with omega: it swaps the eigenspaces
        assert np.allclose(gens[km1] @ w, -w @ gens[km1])
        return [plus, minus]
    return [n // dim[km1]]


print("3. restriction Cl_{0,k} -> Cl_{0,k-1} and its cokernel")
expected = ["Z_2", "Z_2", "0", "Z", "0", "0", "0", "Z"]


def coker(M):
    """Describe Z^rows / (column span of M)."""
    S = smith_normal_form(sp.Matrix(M), domain=sp.ZZ)
    rows = S.shape[0]
    diag = [abs(S[i, i]) for i in range(min(S.shape))]
    parts = [f"Z_{d}" for d in diag if d > 1]
    free = rows - sum(1 for d in diag if d != 0)
    parts += ["Z"] * free
    return " + ".join(parts) if parts else "0"


for k in range(1, 9):
    cols = [restrict(g, k - 1) for g in irreducibles(k)]
    M = [[c[i] for c in cols] for i in range(len(cols[0]))]
    res = coker(M)
    print(f"  k={k}: Cl_{{0,{k-1}}} irr dims {[dim[k-1]] * len(cols[0])},"
          f" Cl_{{0,{k}}} irr dims {[dim[k]] * len(cols)}, restriction {M} -> {res}")
    assert res == expected[k - 1]


print("4. complex Clifford algebras Cl_k(C)")
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Zp = np.diag([1, -1]).astype(complex)


def kron_all(ms):
    out = np.eye(1, dtype=complex)
    for m in ms:
        out = np.kron(out, m)
    return out


def gammas(m):
    """2m Hermitian anticommuting matrices squaring to 1 on C^{2^m} (Jordan-Wigner)."""
    I2 = np.eye(2, dtype=complex)
    out = []
    for j in range(m):
        out.append(kron_all([Zp] * j + [X] + [I2] * (m - j - 1)))
        out.append(kron_all([Zp] * j + [Y] + [I2] * (m - j - 1)))
    return out


def complex_irreducibles(k):
    """Irreducible Cl_k(C) modules, generators e_j = i gamma_j (e_j^2 = -1)."""
    m = k // 2
    gs = gammas(m)
    if k % 2 == 0:
        return [[1j * g for g in gs]]
    chi = kron_all([Zp] * m)  # chirality: anticommutes with all gammas, squares to 1
    return [[1j * g for g in gs] + [s * 1j * chi] for s in (1, -1)]


def complex_pseudo(gens, n):
    """Pseudoscalar normalised to square 1."""
    w = np.eye(n, dtype=complex)
    for g in gens:
        w = w @ g
    sq = (w @ w)[0, 0]
    return w / np.sqrt(sq)


for k in range(0, 9):
    for gens in complex_irreducibles(k):
        n = gens[0].shape[0] if gens else 1
        for a_, b_ in itertools.combinations(gens, 2):
            assert np.allclose(a_ @ b_, -b_ @ a_)
        for a_ in gens:
            assert np.allclose(a_ @ a_, -np.eye(n))

res_c = []
for k in range(1, 9):
    sub_irr = complex_irreducibles(k - 1)
    n_sub = sub_irr[0][0].shape[0] if k > 1 else 1
    cols = []
    for gens in complex_irreducibles(k):
        n = gens[0].shape[0]
        sub = gens[:k - 1]
        if len(sub_irr) == 2:
            w = complex_pseudo(sub, n)
            ref = complex_pseudo(sub_irr[0], n_sub)[0, 0]
            ev = np.round(np.linalg.eigvals(w).real).astype(int)
            cols.append([int((ev == round(ref.real)).sum()) // n_sub,
                         int((ev == -round(ref.real)).sum()) // n_sub])
        else:
            cols.append([n // n_sub])
    M = [[c[i] for c in cols] for i in range(len(cols[0]))]
    res_c.append(coker(M))
print("  K~(S^k), k = 1..8:", res_c)
assert res_c == ["0", "Z"] * 4
for s, gens in zip((1, -1), complex_irreducibles(1)):
    e1 = gens[0][0, 0]
    g = np.array([np.cos(t) + np.sin(t) * e1 for t in th])
    print(f"  Cl_1(C) irreducible e_1 = {e1:+.0f}: clutching winding {winding(g)}")
    assert winding(g) == (1 if e1 == 1j else -1)

print("5. the two irreducibles of Cl_{0,3} = 2H")
Hm = [OL[i][:4, :4] for i in range(4)]
for _ in range(10):
    c = rng.normal(size=4)
    c /= np.linalg.norm(c)
    gp = c[0] * Hm[0] + sum(c[i] * Hm[i] for i in (1, 2, 3))
    gm = c[0] * Hm[0] + sum(c[i] * (-Hm[i]) for i in (1, 2, 3))
    cc = c * np.array([1, -1, -1, -1])
    Lcc = sum(cc[i] * Hm[i] for i in range(4))
    assert np.allclose(gm, Lcc)
    assert np.allclose(gp @ np.eye(4)[0], c) and np.allclose(gm @ np.eye(4)[0], cc)
w3 = [int(pseudoscalar(g, 4)[0, 0]) for g in irreducibles(3)]
print(f"  pseudoscalar signs {w3}; second module gives L_(c*): OK")

print("6. periodicity: Cl_{0,k+8} = Cl_{0,k} (x) Cl_{0,8}")
w8 = pseudoscalar(R[8], 16)
assert np.allclose(w8 @ w8, np.eye(16))
for g in R[8]:
    assert np.allclose(g @ w8, -w8 @ g)
for k in range(1, 9):
    big = [np.kron(g, w8) for g in R[k]] + [np.kron(np.eye(dim[k]), f) for f in R[8]]
    n = big[0].shape[0]
    for a_, b_ in itertools.combinations(big, 2):
        assert np.allclose(a_ @ b_, -b_ @ a_)
    for a_ in big:
        assert np.allclose(a_ @ a_, -np.eye(n))
    assert n == 16 * dim[k]
print("  Cl_{0,k+8} acts on R^{16 a_k} (k = 1..8), a_{k+8} = 16 a_k: OK")
print("All checks passed.")
