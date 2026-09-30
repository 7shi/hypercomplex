"""Checks for ktheory/06-sphere-dirac.md (the Dirac operator on the sphere).

Conventions (clif-analysis/05): Cl_{n,0}(R), e_a^2 = +1, indices 0..n-1,
x = sum x_a e_a, D = sum e_a d_a, D^2 = Delta. E = sum x_a d_a (Euler operator),
L_ab = x_a d_b - x_b d_a (rotation fields, tangent to every sphere),
Gamma = sum_{a<b} e_a e_b L_ab. Inner product <F, G> = scalar part of rev(F) G.

1. Decomposition: x D = E + Gamma, and sum_{a<b} L_ab^2 = r^2 Delta - E^2 - (n-2) E
   (the spherical Laplacian Delta_S), and Delta_S = (n-2) Gamma - Gamma^2 on
   Cl-valued functions.
2. Left multiplication by a bivector e_a e_b is antisymmetric for <,>, and
   int_{S^{n-1}} L_ab f = 0 (rotation invariance), so Gamma is symmetric.
3. Spherical monogenics: the Cauchy-Kovalevskaya extension
   P = sum_j x_0^j / j! (-e_0 D')^j p of a degree-k polynomial p(x_1..x_{n-1})
   is monogenic; Gamma P_k = -k P_k, D(x P_k) = (n + 2k) P_k,
   Gamma (x P_k) = (k + n - 1) x P_k. dim M_k = 2^n C(k+n-2, n-2).
4. Spinor connection: nabla_X psi = d_X psi + c X x psi is compatible with the
   tangent Clifford multiplication X . psi = X x psi exactly for c = -1/2;
   (X x)^2 = -|X|^2 on the sphere (Cl_{0,n-1} convention = even subalgebra).
5. Dirac operator: sum over the tight frame v_ab = x_a e_b - x_b e_a of
   (v_ab x) nabla_{v_ab} equals (n-1)/2 - Gamma on the unit sphere.
   The rotation fields satisfy sum_{a<b} P(d_{v_ab} v_ab) = 0 (P = tangential projection),
   so the connection Laplacian needs no frame correction term.
6. Lichnerowicz: -sum nabla_{v_ab}^2 = -Delta_S - Gamma + (n-1)/4, hence
   Dirac^2 = nabla^* nabla + R/4 with R = (n-1)(n-2).
7. Spectrum: +(k + (n-1)/2) on P_k, -(k + (n-1)/2) on x P_k; |lambda| >= (n-1)/2 and
   the Friedrich bound m R / (4 (m-1)) with m = n-1 equals (n-1)^2/4.
8. n = 3: with the rotor R(theta, phi) sending (e_1, e_2, e_3) to the frame
   (theta-hat, phi-hat, x), psi = R psi_f turns Dirac into
   h1 (d_theta + cot(theta)/2) + h2 d_phi / sin(theta), h1 = e1 e3, h2 = e2 e3,
   which becomes lie/MEMO's -i(s1 (d_theta + cot/2) + s2 d_phi / sin) under
   h1 -> -i s1, h2 -> -i s2. R(phi + 2 pi) = -R. Multiplicity 2(k+1) of k+1 on
   complex spinors C^2 (j = k + 1/2).
"""

import itertools
import random
from math import comb

import sympy as sp

from common.clifford import MV, Alg, eq, mv


def peq(F, G):
    """Equality of polynomial multivectors (fast: expand only)."""
    return all(sp.expand(v) == 0 for v in (mv(F) - mv(G)).d.values())


def euler(A, F):
    return sum((A.X[a] * A.d(F, a) for a in range(A.n)), A.zero())


def L(A, F, a, b):
    return A.X[a] * A.d(F, b) - A.X[b] * A.d(F, a)


def gamma(A, F):
    return sum((A.e[a] * A.e[b] * L(A, F, a, b)
                for a in range(A.n) for b in range(a + 1, A.n)), A.zero())


def lap_s(A, F):
    return sum((L(A, L(A, F, a, b), a, b)
                for a in range(A.n) for b in range(a + 1, A.n)), A.zero())


def vfield(A, a, b):
    return A.X[a] * A.e[b] - A.X[b] * A.e[a]


def sphere_points(n, count, seed):
    """Rational points on S^{n-1} via inverse stereographic projection."""
    rnd = random.Random(seed)
    pts = []
    for _ in range(count):
        u = [sp.Rational(rnd.randint(-5, 5), rnd.randint(1, 4)) for _ in range(n - 1)]
        s = sum(t ** 2 for t in u)
        pts.append([2 * t / (s + 1) for t in u] + [(s - 1) / (s + 1)])
    return pts


def eq_on_sphere(A, F, G, seed=0, count=4):
    for p in sphere_points(A.n, count, seed):
        sub = dict(zip(A.X, p))
        H = (mv(F, A.neg) - mv(G, A.neg)).map(lambda v: sp.expand(v).subs(sub))
        if any(v != 0 for v in H.d.values()):
            return False
    return True


# --- 1. decomposition ---------------------------------------------------------

def check_decomposition():
    for n in (2, 3, 4, 5):
        A = Alg(n)
        F = A.rnd(n, deg=2, grades=(0, 1, 2) if n == 5 else None)
        assert peq(A.x * A.D(F), euler(A, F) + gamma(A, F))
        r2 = A.r2
        rhs = r2 * A.lap(F) - euler(A, euler(A, F)) - (n - 2) * euler(A, F)
        assert peq(lap_s(A, F), rhs)
        assert peq(lap_s(A, F), (n - 2) * gamma(A, F) - gamma(A, gamma(A, F)))
    print("1. x D = E + Gamma, Delta_S = sum L_ab^2 = (n-2) Gamma - Gamma^2")


# --- 2. symmetry of Gamma -----------------------------------------------------

def sphere_integral(poly, X):
    """Integral over the unit sphere S^{n-1} of a polynomial (exact)."""
    n = len(X)
    total = 0
    for mon, c in sp.Poly(sp.expand(poly), *X).terms():
        if any(k % 2 for k in mon):
            continue
        num = 2 * sp.prod([sp.gamma(sp.Rational(k + 1, 2)) for k in mon])
        total += c * num / sp.gamma(sp.Rational(sum(mon) + n, 2))
    return sp.simplify(total)


def inner(F, G):
    return (F.rev() * G).scalar()


def check_symmetry():
    for n in (2, 3, 4):
        A = Alg(n)
        F, G = A.rnd(10 + n, deg=2), A.rnd(20 + n, deg=2)
        for a in range(n):
            for b in range(a + 1, n):
                B = A.e[a] * A.e[b]
                assert sp.expand(inner(B * F, G) + inner(F, B * G)) == 0
                assert sphere_integral(L(A, mv(inner(F, G)), a, b).scalar(), A.X) == 0
        lhs = sphere_integral(inner(F, gamma(A, G)), A.X)
        rhs = sphere_integral(inner(gamma(A, F), G), A.X)
        assert sp.simplify(lhs - rhs) == 0
    print("2. bivector multiplication antisymmetric, int L_ab f = 0, Gamma symmetric")


# --- 3. spherical monogenics --------------------------------------------------

def ck_extension(A, p):
    """Monogenic extension of p(x_1..x_{n-1}) (Cl-valued) off x_0 = 0."""
    def Dp(F):
        return sum((A.e[a] * A.d(F, a) for a in range(1, A.n)), A.zero())
    P, term, j = A.zero(), p, 0
    while term.d:
        P = P + term.map(lambda v: v * A.X[0] ** j / sp.factorial(j))
        term = -(A.e[0] * Dp(term))
        j += 1
    return P


def random_poly(A, k, seed, even=None):
    rnd = random.Random(seed)
    mons = [m for m in itertools.product(range(k + 1), repeat=A.n - 1) if sum(m) == k]
    blades = [b for b in range(1 << A.n)
              if even is None or (bin(b).count("1") % 2 == 0) == even]
    return MV({b: sum(rnd.randint(-3, 3) * sp.prod([A.X[i + 1] ** e for i, e in enumerate(m)])
                      for m in mons) for b in blades}, A.neg)


def monogenic_dim(A, k, even=None):
    """Real dimension of Cl-valued (or even/odd-valued) homogeneous monogenic
    polynomials of degree k."""
    mons = [m for m in itertools.product(range(k + 1), repeat=A.n) if sum(m) == k]
    unknowns = []
    F = A.zero()
    for b in range(1 << A.n):
        if even is not None and (bin(b).count("1") % 2 == 0) != even:
            continue
        for m in mons:
            c = sp.Symbol(f"c{len(unknowns)}")
            unknowns.append(c)
            F = F + MV({b: c * sp.prod([A.X[i] ** e for i, e in enumerate(m)])}, A.neg)
    eqs = []
    for v in A.D(F).d.values():
        eqs += sp.Poly(sp.expand(v), *A.X).coeffs()
    M = sp.Matrix([[sp.diff(e, c) for c in unknowns] for e in eqs]) if eqs else sp.zeros(0, len(unknowns))
    return len(unknowns) - M.rank()


def check_monogenics():
    for n in (2, 3, 4):
        A = Alg(n)
        for k in range(4):
            P = ck_extension(A, random_poly(A, k, 100 * n + k))
            assert peq(A.D(P), 0)
            assert peq(euler(A, P), k * P)
            assert peq(gamma(A, P), -k * P)
            xP = A.x * P
            assert peq(A.D(xP), (n + 2 * k) * P)
            assert peq(gamma(A, xP), (k + n - 1) * xP)
        for k in range(3 if n < 4 else 2):
            assert monogenic_dim(A, k) == 2 ** n * comb(k + n - 2, n - 2), (n, k)
            for ev in (True, False):
                assert monogenic_dim(A, k, ev) == 2 ** (n - 1) * comb(k + n - 2, n - 2)
    print("3. Gamma P_k = -k P_k, Gamma(x P_k) = (k+n-1) x P_k, dim M_k = 2^n C(k+n-2,n-2)")


# --- 4-6. connection, Dirac operator, Lichnerowicz ---------------------------

def nabla(A, v, F, c=sp.Rational(-1, 2)):
    """Covariant derivative along a rotation field v = v_ab (given as (a, b))."""
    a, b = v
    return L(A, F, a, b) + c * (vfield(A, a, b) * A.x * F)


def check_connection():
    c = sp.Symbol("c")
    for n in (3, 4):
        A = Alg(n)
        psi = A.rnd(30 + n, deg=1, grades=range(0, n + 1, 2))
        nonzero = False
        for (a, b), (p, q) in itertools.product(
                [(0, 1), (1, 2)], [(0, 2), (1, 2), (0, 1)]):
            X, Y = vfield(A, a, b), vfield(A, p, q)
            # (X x)^2 = -|X|^2 on the sphere
            assert eq_on_sphere(A, X * A.x * X * A.x, -(X * X).scalar())
            lhs = nabla(A, (a, b), Y * A.x * psi, c)
            dXY = L(A, Y, a, b)                        # d_X Y
            normal = (dXY * A.x + A.x * dXY).scalar() / 2
            tang = dXY - normal * A.x                   # P d_X Y (on the sphere)
            rhs = tang * A.x * psi + Y * A.x * nabla(A, (a, b), psi, c)
            diff = mv(lhs) - mv(rhs)
            for pnt in sphere_points(n, 2, 7 + a + p):
                sub = dict(zip(A.X, pnt))
                vals = [sp.expand(v.subs(sub)) for v in diff.d.values()]
                assert all(v.subs(c, sp.Rational(-1, 2)) == 0 for v in vals)
                nonzero |= any(v.subs(c, 0) != 0 for v in vals)
        assert nonzero
    print("4. nabla_X = d_X - (1/2) X x is compatible with X . psi = X x psi; (X x)^2 = -|X|^2")


def dirac_s(A, F):
    return sum((vfield(A, a, b) * A.x * nabla(A, (a, b), F)
                for a in range(A.n) for b in range(a + 1, A.n)), A.zero())


def check_dirac_lichnerowicz():
    for n in (2, 3, 4, 5):
        A = Alg(n)
        acc = A.zero()
        for a in range(n):
            for b in range(a + 1, n):
                w = L(A, vfield(A, a, b), a, b)          # d_v v
                nrm = (w * A.x + A.x * w).scalar() / 2
                acc = acc + w - nrm * A.x                  # tangential part
        assert eq_on_sphere(A, acc, 0)
        F = A.rnd(40 + n, deg=2 if n < 5 else 1, grades=range(0, n + 1, 2))
        half = sp.Rational(n - 1, 2)
        target = half * F - gamma(A, F)
        assert eq_on_sphere(A, dirac_s(A, F), target)
        rough = -sum((nabla(A, (a, b), nabla(A, (a, b), F))
                      for a in range(n) for b in range(a + 1, n)), A.zero())
        assert eq_on_sphere(A, rough, -lap_s(A, F) - gamma(A, F) + sp.Rational(n - 1, 4) * F)
        # Dirac^2 = rough + R/4, computed with the extrinsic Dirac on the sphere
        # (Gamma preserves the degree, so apply (n-1)/2 - Gamma twice)
        D2 = half * target - gamma(A, target)
        R = (n - 1) * (n - 2)
        assert eq_on_sphere(A, D2, rough + sp.Rational(R, 4) * F)
    print("5. Dirac = (n-1)/2 - Gamma on the sphere")
    print("5'. sum_{a<b} nabla_{v_ab} v_ab = 0")
    print("6. -sum nabla^2 = -Delta_S - Gamma + (n-1)/4, Dirac^2 = nabla*nabla + (n-1)(n-2)/4")


def check_spectrum():
    for n in range(2, 12):
        m = n - 1
        low = sp.Rational(n - 1, 2) ** 2
        R = (n - 1) * (n - 2)
        assert low >= sp.Rational(R, 4)
        if m >= 2:
            assert sp.Rational(m * R, 4 * (m - 1)) == low
    print("7. min lambda^2 = (n-1)^2/4 > R/4, equal to the Friedrich bound")


# --- 8. n = 3: frame form -----------------------------------------------------

def num_zero(F, th, ph, count=5):
    rnd = random.Random(3)
    for _ in range(count):
        sub = {th: sp.Float(rnd.uniform(0.2, 2.9), 30), ph: sp.Float(rnd.uniform(0, 6.2), 30)}
        if any(abs(sp.N(v.subs(sub), 30)) > 1e-20 for v in F.d.values()):
            return False
    return True


def check_n3():
    A = Alg(3)
    th, ph = sp.symbols("theta phi", real=True)
    e1, e2, e3 = A.e
    e12, e31 = e1 * e2, e3 * e1
    Rphi = sp.cos(ph / 2) - e12 * sp.sin(ph / 2)
    Rth = sp.cos(th / 2) - e31 * sp.sin(th / 2)
    R = Rphi * Rth
    Rr = R.rev()
    sub = {A.X[0]: sp.sin(th) * sp.cos(ph), A.X[1]: sp.sin(th) * sp.sin(ph),
           A.X[2]: sp.cos(th)}
    xs = A.x.map(lambda v: v.subs(sub))
    assert eq((R * e3 * Rr - xs).simp(), 0)
    that = (R * e1 * Rr).simp()
    phat = (R * e2 * Rr).simp()
    assert eq(that, A.x.map(lambda v: sp.diff(v.subs(sub), th)))
    assert eq((phat * sp.sin(th)).simp(), A.x.map(lambda v: sp.diff(v.subs(sub), ph)))
    assert eq(R.map(lambda v: v.subs(ph, ph + 2 * sp.pi)).simp(), -R)

    # Dirac in spherical coordinates: tangent Clifford multiplication by the
    # orthonormal frame and nabla along it, acting on psi = R psi_f.
    f = [sp.cos(th) * sp.sin(ph / 2), sp.sin(th) ** 2 * sp.cos(3 * ph / 2),
         sp.cos(2 * th) + sp.sin(ph / 2), sp.sin(th) * sp.cos(th) * sp.cos(ph / 2)]
    psi_f = MV({0: f[0], 0b011: f[1], 0b101: f[2], 0b110: f[3]}, A.neg)
    psi = R * psi_f

    def dth(F):
        return F.map(lambda v: sp.diff(v, th))

    def dph(F):
        return F.map(lambda v: sp.diff(v, ph))

    nab_th = dth(psi) - sp.Rational(1, 2) * that * xs * psi
    nab_ph = dph(psi) / sp.sin(th) - sp.Rational(1, 2) * phat * xs * psi
    Dpsi = that * xs * nab_th + phat * xs * nab_ph
    h1, h2 = e1 * e3, e2 * e3
    expect = h1 * (dth(psi_f) + sp.cot(th) / 2 * psi_f) + h2 * dph(psi_f) / sp.sin(th)
    assert num_zero(Rr * Dpsi - expect, th, ph)

    # (theta-hat, phi-hat) Clifford multiplication agrees with the ambient formula
    # (n-1)/2 - Gamma: compare on a polynomial restricted to the sphere.
    F = A.rnd(77, deg=2, grades=(0, 2))
    lhs = (sp.Integer(1) * F - gamma(A, F)).map(lambda v: v.subs(sub))
    Fs = F.map(lambda v: v.subs(sub))
    n_th = dth(Fs) - sp.Rational(1, 2) * that * xs * Fs
    n_ph = dph(Fs) / sp.sin(th) - sp.Rational(1, 2) * phat * xs * Fs
    rhs = that * xs * n_th + phat * xs * n_ph
    assert num_zero(lhs - rhs, th, ph)

    # representation h1 -> -i s1, h2 -> -i s2 of Cl_{0,2} (= even part of Cl_{3,0})
    I = sp.I
    s1 = sp.Matrix([[0, 1], [1, 0]])
    s2 = sp.Matrix([[0, -I], [I, 0]])
    H1, H2 = -I * s1, -I * s2
    assert H1 ** 2 == -sp.eye(2) and H2 ** 2 == -sp.eye(2)
    assert H1 * H2 + H2 * H1 == sp.zeros(2)
    assert eq(h1 * h1, -1) and eq(h2 * h2, -1) and eq(h1 * h2 + h2 * h1, 0)

    # multiplicity on complex spinors: even-valued monogenics of degree k have real
    # dimension 4 (k+1); as a complex space of C^2-valued functions this is 2(k+1).
    for k in range(3):
        assert monogenic_dim(A, k, True) // 2 == 2 * (k + 1)
        j = sp.Rational(2 * k + 1, 2)
        assert (k + 1) ** 2 == j * (j + 1) + sp.Rational(1, 4)
        assert 2 * j + 1 == 2 * (k + 1)
    print("8. n=3: psi = R psi_f gives h1(d_th + cot/2) + h2 d_ph/sin, R(phi+2pi) = -R")


if __name__ == "__main__":
    check_decomposition()
    check_symmetry()
    check_monogenics()
    check_connection()
    check_dirac_lichnerowicz()
    check_spectrum()
    check_n3()
    print("06-sphere-dirac: all checks passed")
