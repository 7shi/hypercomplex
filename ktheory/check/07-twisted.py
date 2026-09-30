"""Checks for ktheory/07-twisted.md (the Dirac operator on S^2 twisted by H^n).

Conventions (ktheory/06, n = 3): Cl_{3,0}(R), e_a^2 = +1, x = x0 e0 + x1 e1 + x2 e2
on the unit sphere, spinors take values in the even subalgebra, X . psi = X x psi,
nabla_X = d_X - (1/2) X x, D_S = sum t_i . nabla_{t_i}. J = e0 e1 (the complex
structure is right multiplication by J), omega = e0 e1 e2, f1 = e0 e2, f2 = e1 e2,
U = exp(-J phi/2) exp(-e2 e0 theta/2) (Ũ x U = e2).

1. Chirality gamma psi = omega x psi J: gamma^2 = 1, anticommutes with X ., commutes
   with nabla_X and with right multiplication by J, is symmetric; D_S gamma = -gamma D_S.
   In the frame: Ũ (omega x) U = J, so gamma_f psi_f = J psi_f J; S^+ = U span(f1, f2).
2. Phi(alpha, beta) = f1 alpha - beta (alpha, beta in C_J = span(1, J)): for h(alpha, beta)
   = (2 alpha conj(beta), |alpha|^2 - |beta|^2) = x (unit vector), gamma Phi = Phi, and
   Phi(-conj beta, conj alpha) has gamma = -1. Local sections v_+ = (cos, sin e^{-J phi}),
   v_- = (cos e^{J phi}, sin) (half angles) with h = x, v_+ = v_- e^{-J phi}, and their
   expressions in x. Transitions: S^+ ~ z, S^- ~ z^{-1}.
3. Twisted connection nabla^a = nabla + (a . X) R_J. Gluing psi_+ = psi_- e^{n J phi}
   requires a_+ - a_- = -n grad(phi), grad(phi) = v01/(x0^2 + x1^2). For a = f(x2) v01 the
   curvature F = t2 . d_{t1} a - t1 . d_{t2} a equals -d/dx2[(1 - x2^2) f]; constant F = c
   gives f_+ = c/(1+x2), f_- = -c/(1-x2), and gluing forces c = -n/2.
   Stokes: int F dOmega = oint (a_+ - a_-) . dx = -2 pi n (numerically).
   D_n = D_S + a x R_J, symmetric, anticommutes with gamma; psi -> psi K (K = f2)
   swaps chirality and D_{-a}(psi K) = (D_a psi) K (H^n -> H^{-n}).
   n = 1: a_- . d_phi x = cos^2(theta/2) = <e, d_phi e> for the unit frame e = v_-
   of the Hopf line (projection connection); <e, d_theta e> = 0.
4. Lichnerowicz: for any tangent a, D^2 = nabla^a* nabla^a + 1/2 - F gamma;
   for the constant connection D_n^2 = nabla* nabla + 1/2 + (n/2) gamma.
5. Frame form on U_+: D_n = f1 (d_theta + cot/2) + f2/sin (d_phi + A R_J),
   A = -(n/2)(1 - cos theta); modes f1 g e^{J mu phi} with
   g_mu = sin^{-mu-1/2}(theta/2) cos^{mu-1/2-n}(theta/2) solve D_n = 0.
   Phi(v_+) = U f1 e^{-J phi/2}. Kernel elements Phi(v) P(v) (P homogeneous of degree
   m - 1, n = -m) satisfy D_n = 0 in both gauges and glue; |Phi(v)| = 1 and it is parallel
   for n = -1. Counting: dim ker D^+ = max(-n, 0), dim ker D^- = max(n, 0), index -n.
"""

import random

import sympy as sp

from common.clifford import MV, Alg, mv

A = Alg(3)
e0, e1, e2 = A.e
X0, X1, X2 = A.X
x = A.x
w = A.I
J = e0 * e1
f1, f2 = e0 * e2, e1 * e2
half = sp.Rational(1, 2)
th, ph = sp.symbols("theta phi", real=True)
SPH = {X0: sp.sin(th) * sp.cos(ph), X1: sp.sin(th) * sp.sin(ph), X2: sp.cos(th)}


def L(F, a, b):
    return A.X[a] * A.d(F, b) - A.X[b] * A.d(F, a)


def vfield(a, b):
    return A.X[a] * A.e[b] - A.X[b] * A.e[a]


PAIRS = [(0, 1), (0, 2), (1, 2)]


def dot(u, v):
    return (u * v + v * u).scalar() / 2


def ex(t, B):
    """exp(B t) for B^2 = -1."""
    return sp.cos(t) + B * sp.sin(t)


def cj(c):
    """complex conjugate in C_J = span(1, J)."""
    return MV({k: (-v if k == 0b011 else v) for k, v in mv(c, A.neg).d.items()}, A.neg)


def nab(F, a, b, al=None):
    """Covariant derivative along the rotation field v_ab (twisted by al)."""
    v = vfield(a, b)
    r = L(F, a, b) - half * v * x * F
    if al is not None:
        r = r + dot(al, v) * F * J
    return r


def dirac(F, al=None):
    return sum((vfield(a, b) * x * nab(F, a, b, al) for a, b in PAIRS), A.zero())


def rough(F, al=None):
    return -sum((nab(nab(F, a, b, al), a, b, al) for a, b in PAIRS), A.zero())


def gam(F):
    return w * x * F * J


def sphere_points(count, seed):
    rnd = random.Random(seed)
    pts = []
    for _ in range(count):
        u = [sp.Rational(rnd.randint(-5, 5), rnd.randint(1, 4)) for _ in range(2)]
        s = sum(t ** 2 for t in u)
        pts.append([2 * u[0] / (s + 1), 2 * u[1] / (s + 1), (s - 1) / (s + 1)])
    return pts


def zero_on_sphere(F, count=4, seed=0, numeric=False):
    for p in sphere_points(count, seed):
        sub = dict(zip(A.X, p))
        for v in mv(F, A.neg).d.values():
            val = v.subs(sub)
            if numeric:
                if abs(sp.N(val, 40)) > 1e-25:
                    return False
            elif sp.simplify(val) != 0:
                return False
    return True


def zero_num(F, count=5, seed=3, lo=0.2, hi=2.9):
    """F(theta, phi) = 0 at random numeric points."""
    rnd = random.Random(seed)
    for _ in range(count):
        s = {th: sp.Float(rnd.uniform(lo, hi), 40), ph: sp.Float(rnd.uniform(0, 6.2), 40)}
        if any(abs(sp.N(v.subs(s), 40)) > 1e-25 for v in mv(F, A.neg).d.values()):
            return False
    return True


def to_sph(F):
    return mv(F, A.neg).map(lambda v: v.subs(SPH))


def nth_ph(F):
    return F.map(lambda v: sp.diff(v, th)), F.map(lambda v: sp.diff(v, ph))


# --- 1. chirality -------------------------------------------------------------

def check_chirality():
    F = A.rnd(5, deg=2, grades=(0, 2))
    G = A.rnd(6, deg=1, grades=(0, 2))
    assert zero_on_sphere(gam(gam(F)) - F)
    for a, b in PAIRS:
        X = vfield(a, b)
        assert zero_on_sphere(X * x * gam(F) + gam(X * x * F))
        assert zero_on_sphere(nab(gam(F), a, b) - gam(nab(F, a, b)))
    assert zero_on_sphere(gam(F * J) - gam(F) * J)
    assert zero_on_sphere((gam(F).rev() * G).scalar() - (F.rev() * gam(G)).scalar())
    assert zero_on_sphere(dirac(gam(F)) + gam(dirac(F)))
    U = ex(-ph / 2, J) * ex(-th / 2, e2 * e0)
    Ur = U.rev()
    assert zero_num(Ur * to_sph(w * x) * U - J)
    assert zero_num(Ur * to_sph(x) * U - e2)
    for B, s in ((1, -1), (J, -1), (f1, 1), (f2, 1)):
        assert (J * mv(B, A.neg) * J - s * mv(B, A.neg)).d == {}
    # t1 . t2 . psi = -omega x psi for (t1, t2, x) with t1 t2 x = omega
    t1 = to_sph(x).map(lambda v: sp.diff(v, th))
    t2 = to_sph(x).map(lambda v: sp.diff(v, ph) / sp.sin(th))
    xs = to_sph(x)
    assert zero_num(t1 * t2 * xs - w)
    assert zero_num(t1 * xs * t2 * xs + w * xs)
    print("1. gamma = omega x (.) J: gamma^2 = 1, parallel, anticommutes with X . and D_S;"
          " frame: J psi_f J, S^+ = span(f1, f2)")


# --- 2. Hopf bundle and chirality ---------------------------------------------

def phi_map(al, be):
    return f1 * al - be


def check_hopf():
    rnd = random.Random(11)
    for _ in range(4):
        a0, a1, b0, b1 = (sp.Rational(rnd.randint(-6, 6), rnd.randint(1, 5)) for _ in range(4))
        al, be = a0 + a1 * J, b0 + b1 * J
        nrm = (al * cj(al) + be * cj(be)).scalar()
        c = 2 * al * cj(be) / nrm                      # x0 + x1 J
        r = (al * cj(al) - be * cj(be)).scalar() / nrm
        xv = c.scalar() * e0 + c.d.get(0b011, 0) * e1 + r * e2
        P = phi_map(al, be)
        assert (w * xv * P * J - P).d == {}
        Q = phi_map(-cj(be), cj(al))                  # the orthogonal line, over -x
        assert (w * xv * Q * J + Q).d == {}
    # local sections
    c, s = sp.cos(th / 2), sp.sin(th / 2)
    vm = (c * ex(ph, J), mv(s, A.neg))
    vp = (mv(c, A.neg), s * ex(-ph, J))
    for al, be in (vm, vp):
        assert zero_num(2 * al * cj(be) - to_sph(X0 + X1 * J))
        assert zero_num(al * cj(al) - be * cj(be) - to_sph(mv(X2, A.neg)))
    assert zero_num(vp[0] - vm[0] * ex(-ph, J)) and zero_num(vp[1] - vm[1] * ex(-ph, J))
    # expressions in x
    assert zero_num(vm[0] - to_sph((X0 + X1 * J) / sp.sqrt(2 * (1 - X2))), lo=0.3)
    assert zero_num(vm[1] - to_sph(mv(sp.sqrt((1 - X2) / 2), A.neg)), lo=0.3)
    assert zero_num(vp[0] - to_sph(mv(sp.sqrt((1 + X2) / 2), A.neg)), hi=2.8)
    assert zero_num(vp[1] - to_sph((X0 - X1 * J) / sp.sqrt(2 * (1 + X2))), hi=2.8)
    # transitions: S^+ spanned by Phi(v), S^- by Phi(-conj b, conj a)
    Pp, Pm = phi_map(*vp), phi_map(*vm)
    assert zero_num(Pp - Pm * ex(-ph, J))          # chi = Pp c+ = Pm c-  =>  c+ = e^{J phi} c-
    Qp, Qm = phi_map(-cj(vp[1]), cj(vp[0])), phi_map(-cj(vm[1]), cj(vm[0]))
    assert zero_num(Qp - Qm * ex(ph, J))           # c+ = e^{-J phi} c-
    print("2. Phi(alpha, beta) = f1 alpha - beta maps the Hopf line over x into S^+ (and the"
          " orthogonal line into S^-); S^+ ~ z, S^- ~ z^{-1}; local sections v_+-")


# --- 3. twisted connection and curvature --------------------------------------

n = sp.Symbol("n")
V01 = vfield(0, 1)
AP = (-n / 2) * V01 / (1 + X2)
AM = (n / 2) * V01 / (1 - X2)


def curvature(al):
    """F with sum_i t_i ^ nabla_{t_i} a = F omega x (via rotation fields)."""
    curl = A.zero()
    for a, b in PAIRS:
        v = vfield(a, b)
        dv = L(al, a, b)
        dv = dv - dot(dv, x) * x
        curl = curl + (v * dv - dv * v) / 2
    return curl, -(curl * w * x).scalar()


def check_connection():
    # gluing: d_X phi = grad(phi) . X, grad(phi) = v01/(x0^2 + x1^2)
    phi_fn = sp.atan2(X1, X0)
    for a, b in PAIRS:
        v = vfield(a, b)
        lhs = sp.simplify(L(mv(phi_fn, A.neg), a, b).scalar())
        assert sp.simplify(lhs - dot(V01, v) / (X0 ** 2 + X1 ** 2)) == 0
    assert zero_on_sphere(AP - AM + n * V01 / (X0 ** 2 + X1 ** 2))
    # F for a = f(x2) v01
    f = sp.Function("f")
    curl, F = curvature(f(X2) * V01)
    assert zero_on_sphere((curl - F * w * x).map(lambda t: t.subs(f(X2), X2 ** 2 + 3)
                                                   .doit()))
    Fexp = -sp.diff((1 - X2 ** 2) * f(X2), X2)
    for fx in (X2 ** 3 + 2, 1 / (2 + X2)):
        assert zero_on_sphere(mv((F - Fexp).subs(f(X2), fx).doit(), A.neg))
    c = sp.Symbol("c")
    for fx in (c / (1 + X2), -c / (1 - X2)):
        assert sp.simplify(Fexp.subs(f(X2), fx).doit() - c) == 0
    assert sp.solve(sp.Eq(c / (1 + X2) + c / (1 - X2), -n / (1 - X2 ** 2)), c) == [-n / 2]
    assert zero_on_sphere(curvature(AP)[1] + n / 2) and zero_on_sphere(curvature(AM)[1] + n / 2)
    # frame coefficient A(theta) = a_+ . d_phi x
    assert sp.simplify(to_sph(mv(dot(AP, V01), A.neg)).scalar()
                       + n / 2 * (1 - sp.cos(th))) == 0
    # Stokes: oint (a_+ - a_-) . dx along the equator (phi increasing) = -2 pi n
    eq_sub = {X0: sp.cos(ph), X1: sp.sin(ph), X2: 0}
    integrand = sp.simplify(dot(AP - AM, V01).subs(eq_sub))
    assert sp.simplify(sp.integrate(integrand, (ph, 0, 2 * sp.pi)) + 2 * sp.pi * n) == 0
    # hemisphere integrals of the constant curvature: -pi n each
    assert sp.integrate(-n / 2 * sp.sin(th), (th, 0, sp.pi / 2), (ph, 0, 2 * sp.pi)) == -sp.pi * n
    # D_n = D_S + a x R_J, gamma, symmetry of a x R_J
    F0 = A.rnd(7, deg=2, grades=(0, 2))
    for al in (AP, AM):
        assert zero_on_sphere(dirac(F0, al) - dirac(F0) - al * x * F0 * J)
        assert zero_on_sphere(gam(al * x * F0 * J) + al * x * gam(F0) * J)
    G0 = A.rnd(8, deg=1, grades=(0, 2))
    assert zero_on_sphere(((AP * x * F0 * J).rev() * G0).scalar()
                          - (F0.rev() * (AP * x * G0 * J)).scalar())
    # psi -> psi K: gluing, chirality, D_n(psi K) = (D_{-n} psi) K
    K = f2
    assert (K * J + J * K).d == {}
    assert (ex(n * ph, J) * K - K * ex(-n * ph, J)).map(sp.expand).d == {}
    assert zero_on_sphere(gam(F0 * K) + gam(F0) * K)
    assert zero_on_sphere(dirac(F0 * K, AP) - dirac(F0, AP.map(lambda t: t.subs(n, -n))) * K)
    assert zero_on_sphere(dirac(F0 * K, AP.map(lambda t: t.subs(n, -n))) - dirac(F0, AP) * K)
    # n = 1 projection connection: e = v_- (unit vector of the Hopf line)
    cc, ss = sp.cos(th / 2), sp.sin(th / 2)
    al, be = cc * ex(ph, J), mv(ss, A.neg)
    dth = lambda F: F.map(lambda t: sp.diff(t, th))
    dph = lambda F: F.map(lambda t: sp.diff(t, ph))
    assert zero_num(cj(al) * dph(al) + cj(be) * dph(be) - cc ** 2 * J)
    assert zero_num(cj(al) * dth(al) + cj(be) * dth(be))
    am1 = to_sph(AM.map(lambda t: t.subs(n, 1)))
    assert zero_num(mv(dot(am1, to_sph(V01)) - cc ** 2, A.neg))
    print("3. gluing a_+ - a_- = -n grad(phi); F = -d/dx2[(1-x2^2) f] for a = f(x2) v01;"
          " constant F forces a_+- and F = -n/2; oint = -2 pi n; D_n = D_S + a x R_J;"
          " psi K: n <-> -n; n = 1 is the projection connection")


# --- 4. Lichnerowicz ------------------------------------------------------------

def check_lichnerowicz():
    beta = (X0 * X1 + 2) * e0 + (X2 - X0 ** 2) * e1 + (X1 * X2 + X0) * e2
    al = beta - dot(beta, x) * x
    F0 = A.rnd(9, deg=1, grades=(0, 2))
    _, Fs = curvature(al)
    assert zero_on_sphere(dirac(dirac(F0, al), al) - rough(F0, al) - half * F0 + Fs * gam(F0))
    # rough Laplacian expansion: nabla^a* nabla^a = nabla* nabla - div(a) R_J - 2 nabla_a R_J + |a|^2
    div = sum((dot(L(al, a, b), vfield(a, b)) for a, b in PAIRS), 0)
    na = sum((dot(al, vfield(a, b)) * nab(F0, a, b) for a, b in PAIRS), A.zero())
    assert zero_on_sphere(rough(F0, al) - rough(F0) + div * F0 * J + 2 * na * J
                          - dot(al, al) * F0)
    for al in (AP, AM):
        for nn in (-2, 1, 3):
            a_n = al.map(lambda t: t.subs(n, nn))
            D2 = dirac(dirac(F0, a_n), a_n)
            assert zero_on_sphere(D2 - rough(F0, a_n) - half * F0
                                  - sp.Rational(nn, 2) * gam(F0))
    print("4. D^2 = nabla^a* nabla^a + 1/2 - F gamma; constant connection: + (n/2) gamma")


# --- 5. kernel ------------------------------------------------------------------

def frame_dirac(psi_f, nn):
    A_th = -sp.Rational(nn, 2) * (1 - sp.cos(th))
    d_t, d_p = nth_ph(psi_f)
    return f1 * (d_t + sp.cot(th) / 2 * psi_f) + f2 / sp.sin(th) * (d_p + A_th * psi_f * J)


def ambient_dirac_sph(psi, al):
    """D_n via the tangent frame (theta-hat, phi-hat) on functions of (theta, phi)."""
    xs = to_sph(x)
    that = xs.map(lambda v: sp.diff(v, th))
    phat = xs.map(lambda v: sp.diff(v, ph) / sp.sin(th))
    als = to_sph(al)
    d_t, d_p = nth_ph(psi)
    n_t = d_t - half * that * xs * psi + dot(als, that) * psi * J
    n_p = d_p / sp.sin(th) - half * phat * xs * psi + dot(als, phat) * psi * J
    return that * xs * n_t + phat * xs * n_p, (n_t, n_p)


def check_kernel():
    U = ex(-ph / 2, J) * ex(-th / 2, e2 * e0)
    # frame form of D_n (twisted) against the ambient definition
    psi_f = MV({0: sp.cos(th) * sp.sin(ph / 2), 0b011: sp.sin(th) ** 2 * sp.cos(3 * ph / 2),
                0b101: sp.cos(2 * th) + sp.sin(ph / 2), 0b110: sp.sin(th) * sp.cos(ph / 2)},
               A.neg)
    for nn in (-2, 3):
        Dn, _ = ambient_dirac_sph(U * psi_f, AP.map(lambda t: t.subs(n, nn)))
        assert zero_num(U.rev() * Dn - frame_dirac(psi_f, nn))
    # modes: f1 g_mu e^{J mu phi} solves D_n = 0
    for nn in (-3, -1, 0, 2):
        for mu in (-sp.Rational(5, 2), -half, half, sp.Rational(3, 2)):
            g = sp.sin(th / 2) ** (-mu - half) * sp.cos(th / 2) ** (mu - half - nn)
            assert zero_num(frame_dirac(f1 * g * ex(mu * ph, J), nn))
    # Phi(v_+) = U f1 e^{-J phi/2}
    cc, ss = sp.cos(th / 2), sp.sin(th / 2)
    vm = (cc * ex(ph, J), mv(ss, A.neg))
    vp = (mv(cc, A.neg), ss * ex(-ph, J))
    assert zero_num(phi_map(*vp) - U * f1 * ex(-ph / 2, J))
    # kernel elements Phi(v) P(v) in both gauges
    for m in (1, 2, 3):
        nn = -m
        ap, am = (a.map(lambda t: t.subs(n, nn)) for a in (AP, AM))
        for p in range(m):
            psis = []
            for (al, be), a_g, lo, hi in ((vp, ap, 0.1, 2.8), (vm, am, 0.3, 3.0)):
                mono = mv(1, A.neg)
                for _ in range(m - 1 - p):
                    mono = mono * al
                for _ in range(p):
                    mono = mono * be
                psi = phi_map(al, be) * mono
                Dn, nabs = ambient_dirac_sph(psi, a_g)
                assert zero_num(Dn, lo=lo, hi=hi)
                assert zero_num(to_sph(w * x) * psi * J - psi)
                if m == 1:
                    assert all(zero_num(t, lo=lo, hi=hi) for t in nabs)
                psis.append(psi)
            assert zero_num(psis[0] - psis[1] * ex(nn * ph, J))
            # frame mode: mu = -p - 1/2
            g = ss ** p * cc ** (m - 1 - p)
            assert zero_num(psis[0] - U * f1 * g * ex(-(p + half) * ph, J))
    assert zero_num(phi_map(*vp).rev() * phi_map(*vp) - 1)
    # counting bounded modes: n + 1/2 <= mu <= -1/2
    for nn in range(-6, 7):
        plus = sum(1 for k in range(-20, 20)
                   if (-(k + half) - half >= 0 and (k + half) - half - nn >= 0))
        minus = max(nn, 0)
        assert plus == max(-nn, 0) and plus - minus == -nn
    print("5. frame form, modes g_mu, Phi(v_+) = U f1 e^{-J phi/2}; Phi(v) P(v) in the kernel"
          " (both gauges, glued), parallel for n = -1; dim ker D^+ = max(-n,0), index -n")


def check_toeplitz_kernel():
    """(S*)^m on the Hardy space: kernel spanned by 1, z, ..., z^{m-1}."""
    N = 12
    Ss = sp.zeros(N, N)
    for k in range(N - 1):
        Ss[k, k + 1] = 1          # S* e_{k+1} = e_k
    for m in (1, 2, 3):
        T = Ss ** m
        ker = T.nullspace()
        assert len(ker) == m
        assert all(all(v[k] == 0 for k in range(m, N)) for v in ker)
    print("6. ker T_{z^{-m}} = ker (S*)^m = span(1, ..., z^{m-1})")


if __name__ == "__main__":
    check_chirality()
    check_hopf()
    check_connection()
    check_lichnerowicz()
    check_kernel()
    check_toeplitz_kernel()
