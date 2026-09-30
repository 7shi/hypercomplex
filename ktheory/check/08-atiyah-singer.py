"""Checks for ktheory/08-atiyah-singer.md (the Atiyah-Singer index theorem).

Conventions as in ktheory/06 and ktheory/07: Cl_{n,0}(R), e_a^2 = +1, spinors and
forms on S^{n-1} take values in the even subalgebra, X . psi = X x psi. On S^2
(n = 3): J = e0 e1, omega = e0 e1 e2, f1 = e0 e2, gamma psi = omega x psi J,
Phi(alpha, beta) = f1 alpha - beta.

1. Orientation and c1: the complex coordinate z = cot(theta/2) e^{i phi} (= alpha/beta)
   is orientation-reversing with respect to the outward orientation (theta, phi).
   With the Fubini-Study metric (1 + |z|^2)^n on the section (z, 1)^n of H^n,
   the integral of c1 = (i/2pi) dbar d log h in the z-orientation is -n; this equals
   (1/2pi) int F dOmega = -n of ktheory/07 (F = -n/2, outward dOmega).
   gamma = i t2. t1. for t1 t2 x = omega, i.e. gamma is the standard chirality of the
   z-orientation.
2. Principal symbol: at x = e2, xi = cos t e0 + sin t e1, xi . (f1 v) = -e^{-Jt} v.
   At random points, xi . maps S^+ to S^- and, in the unit bases Phi(alpha, beta) and
   Phi(-conj beta, conj alpha), is multiplication by a unit complex number of winding
   -1 in t (outward orientation). For a graded Cl_{0,m} module, e1^{-1} xi restricted
   to S^+ is the paravector xi_1 + sum_l xi_l h_l, h_l = e1^{-1} e_l, h_l^2 = -1,
   anticommuting (the clutching function g_W of ktheory/04).
3. dbar: on U_- write psi_- = Phi(alpha_-, beta_-) (1 + |z|^2)^{(n+1)/2} F. Then the
   frame form of D_n on positive chirality is a nonzero multiple of dF/dzbar.
   For F = z^k the length (1 + |z|^2)^{(n+1)/2} |z|^k stays bounded iff n = -m and
   k <= m - 1 (m polynomials; none for n >= 0). Riemann-Roch on S^2:
   deg + 1 with deg H^{n+1} = -(n+1) gives -n.
4. Forms on S^2 (vector calculus): d f = grad_S f, d a = F(a) = (curl a) . x,
   the adjoint of a -> F(a) is g -> -x cross grad g, the adjoint of grad is -div.
   F(grad f) = 0, div grad f = Delta_S f, div(x cross grad g) = 0,
   |x cross grad g| = |grad g| (exact integrals of polynomials over S^2).
5. Forms on S^4 in the even subalgebra of Cl_{5,0} at x = e0: forms
   t_{i1} ^ ... ^ t_{ik} <-> (t_{i1} x) ... (t_{ik} x); left multiplication by xi x is
   xi ^ - iota_xi (the symbol of d + delta); epsilon psi = x psi x is (-1)^k;
   B = i epsilon R_{v x} = i x psi v = -i (v^ + iota_v) is odd, anticommutes with xi x,
   B^2 = -|v|^2, so
   (xi x + t B)^2 = -(|xi|^2 + t^2 |v|^2); cos s xi x + sin s B has square -1
   for unit xi, v (the extension of ktheory/04).
6. Characteristic numbers: signature = p1/3, A-hat = -p1/24 in dimension 4:
   CP^2 (p1 = 3): sigma = 1, A-hat = -1/8; K3 (sigma = -16): A-hat = 2.
   Todd on S^2: ind dbar_L = c1(L) + 1 (c1(TS^2) = 2 in the z-orientation);
   D_n = dbar on H^{n+1}: c1(H^{n+1}) + 1 = -n = c1(H^n).
"""

import itertools
import math
import random

import sympy as sp

from common.clifford import MV, Alg, eq, mv

half = sp.Rational(1, 2)
th, ph = sp.symbols("theta phi", real=True)


# --- 1. orientation and c1 -----------------------------------------------------

def check_orientation_c1():
    r = sp.cot(th / 2)
    X, Y = r * sp.cos(ph), r * sp.sin(ph)
    jac = sp.Matrix([[X.diff(th), X.diff(ph)], [Y.diff(th), Y.diff(ph)]]).det()
    jac = sp.simplify(jac)
    assert all(jac.subs({th: t, ph: 0.3}) < 0 for t in (0.3, 1.2, 2.5))
    n = sp.symbols("n", integer=True)
    x, y = sp.symbols("x y", real=True)
    logh = n * sp.log(1 + x ** 2 + y ** 2)
    lap = sp.simplify(logh.diff(x, 2) + logh.diff(y, 2))
    # dbar d log h = (1/4) lap dzbar^dz, i dzbar^dz = -2 dx^dy
    dens = sp.simplify(-lap / (4 * sp.pi))
    rho = sp.symbols("rho", positive=True)
    tot = sp.integrate(sp.simplify(dens.subs({x: rho, y: 0})) * 2 * sp.pi * rho, (rho, 0, sp.oo))
    assert sp.simplify(tot + n) == 0
    # 07: F = -n/2 constant, area 4 pi
    assert sp.simplify(-n / 2 * 4 * sp.pi / (2 * sp.pi) + n) == 0
    # gamma = i t2. t1. (i = right J) for t1 t2 x = omega
    A = Alg(3)
    e0, e1, e2 = A.e
    J, w = e0 * e1, A.I
    xx, t1, t2 = e2, e0, e1
    assert eq(t1 * t2 * xx, w)
    psi = A.rnd(1, deg=0, grades=(0, 2))
    assert eq(t2 * xx * t1 * xx * psi * J, w * xx * psi * J)
    assert eq(t1 * xx * t2 * xx * psi * J, -(w * xx * psi * J))
    print("1. z-orientation is inward; c1(H^n) = -n there = (1/2pi) int F dOmega;"
          " gamma = i t2. t1.")


# --- 2. principal symbol ---------------------------------------------------------

def check_symbol():
    A = Alg(3)
    e0, e1, e2 = A.e
    J, w, f1 = e0 * e1, A.I, e0 * e2
    t = sp.symbols("t", real=True)
    xi = sp.cos(t) * e0 + sp.sin(t) * e1
    assert eq(xi * e2 * f1, -(sp.cos(t) - J * sp.sin(t)))

    def cplx(c):  # element of C_J -> complex number
        return complex(sp.N(c.d.get(0, 0))) + 1j * complex(sp.N(c.d.get(0b011, 0)))

    def phi(al, be):
        return f1 * (al[0] + al[1] * J) - (be[0] + be[1] * J)

    rnd = random.Random(5)
    for _ in range(4):
        # random unit (alpha, beta) in C^2
        v = [rnd.uniform(-1, 1) for _ in range(4)]
        s = math.sqrt(sum(u * u for u in v))
        a, b = complex(v[0], v[1]) / s, complex(v[2], v[3]) / s
        c = 2 * a * b.conjugate()
        r = abs(a) ** 2 - abs(b) ** 2
        xx = c.real * e0 + c.imag * e1 + r * e2
        xv = [c.real, c.imag, r]
        # orthonormal tangent t1, t2 with t1 t2 x = omega
        u = [1.0, 0.0, 0.0] if abs(xv[0]) < 0.9 else [0.0, 1.0, 0.0]
        d = sum(p * q for p, q in zip(u, xv))
        t1v = [p - d * q for p, q in zip(u, xv)]
        nn = math.sqrt(sum(p * p for p in t1v))
        t1v = [p / nn for p in t1v]
        t2v = [xv[1] * t1v[2] - xv[2] * t1v[1], xv[2] * t1v[0] - xv[0] * t1v[2],
               xv[0] * t1v[1] - xv[1] * t1v[0]]
        T1 = sum((p * E for p, E in zip(t1v, A.e)), A.zero())
        T2 = sum((p * E for p, E in zip(t2v, A.e)), A.zero())
        dev = (T1 * T2 * xx - w).d.values()
        assert all(abs(complex(sp.N(q))) < 1e-9 for q in dev)
        P = phi((a.real, a.imag), (b.real, b.imag))
        Q = phi((-b.real, b.imag), (a.real, -a.imag))
        vals = []
        N = 64
        for k in range(N):
            tt = 2 * math.pi * k / N
            X = math.cos(tt) * T1 + math.sin(tt) * T2
            img = X * xx * P
            cc = Q.rev() * img
            # the image lies in S^-: Q c
            assert all(abs(complex(sp.N(q))) < 1e-9 for q in (Q * cc - img).d.values())
            assert all(abs(complex(sp.N(q))) < 1e-9 for m_, q in cc.d.items() if m_ not in (0, 0b011))
            vals.append(cplx(cc))
        assert all(abs(abs(z) - 1) < 1e-9 for z in vals)
        wind = 0.0
        for k in range(N):
            z0, z1 = vals[k], vals[(k + 1) % N]
            wind += math.atan2((z1 / z0).imag, (z1 / z0).real)
        assert round(wind / (2 * math.pi)) == -1
    # graded Cl_{0,m} module: e1^{-1} xi is a paravector in h_l = e1^{-1} e_l
    for m in range(2, 7):
        B = Alg(0, m)
        e = B.e
        inv1 = -e[0]
        assert eq(inv1 * e[0], 1)
        h = [inv1 * e[l] for l in range(1, m)]
        for i in range(m - 1):
            assert eq(h[i] * h[i], -1)
            for j in range(i + 1, m - 1):
                assert eq(h[i] * h[j] + h[j] * h[i], 0)
        xs = sp.symbols(f"xi0:{m}", real=True)
        xi = sum((xs[a] * e[a] for a in range(m)), B.zero())
        assert eq(inv1 * xi, xs[0] + sum((xs[l] * h[l - 1] for l in range(1, m)), B.zero()))
    print("2. symbol at e2: xi.(f1 v) = -e^{-Jt} v; S^+ -> S^- with winding -1 (outward);"
          " e1^{-1} xi = paravector (m = 2..6)")


# --- 3. dbar -------------------------------------------------------------------

def check_dbar():
    n = sp.symbols("n", integer=True)
    I = sp.I
    F = sp.Function("F")(th, ph)
    r = sp.cot(th / 2)
    wgt = 1 + r ** 2
    # psi_+ = psi_- e^{n J phi}, Phi(v_-) = Phi(v_+) e^{J phi} = U f1 e^{J phi/2}
    v = sp.exp(I * (n + half) * ph) * wgt ** ((n + 1) / 2) * F
    A_ = -(n / 2) * (1 - sp.cos(th))
    Lv = -(v.diff(th) + sp.cot(th) / 2 * v) + (I / sp.sin(th)) * (v.diff(ph) + I * A_ * v)
    dr = -2 * sp.sin(th / 2) ** 2          # d theta / d r
    dzbar = (sp.exp(I * ph) / 2) * (dr * F.diff(th) + (I / r) * F.diff(ph))
    # Lv = c * dzbar with c free of F (checked for n = -3..3 at numeric points)
    Fth, Fph, F0 = sp.symbols("Fth Fph F0")
    rep = {F.diff(th): Fth, F.diff(ph): Fph}
    pts = [(0.4, 0.7), (1.3, 2.1), (2.6, 4.4)]
    for nn in range(-3, 4):
        Lv_s = sp.expand(Lv.subs(n, nn).subs(rep).subs(F, F0))
        dz_s = sp.expand(dzbar.subs(rep))
        for a_, b_ in pts:
            sub = {th: a_, ph: b_}
            assert abs(complex(sp.N(Lv_s.coeff(F0).subs(sub)))) < 1e-12
            c1 = complex(sp.N((Lv_s.coeff(Fth) / dz_s.coeff(Fth)).subs(sub)))
            c2 = complex(sp.N((Lv_s.coeff(Fph) / dz_s.coeff(Fph)).subs(sub)))
            assert abs(c1 - c2) < 1e-12 and abs(c1) > 1e-6
    # growth: (1+rho^2)^{(n+1)/2} rho^k bounded iff n = -m, k <= m-1
    rho = sp.symbols("rho", positive=True)
    for nn in range(-5, 4):
        good = [k for k in range(0, 10)
                if sp.limit((1 + rho ** 2) ** sp.Rational(nn + 1, 2) * rho ** k, rho, sp.oo) != sp.oo]
        assert len(good) == max(-nn, 0)
        # Riemann-Roch on S^2: deg + 1, deg H^{n+1} = -(n+1) in the z-orientation
        assert -(nn + 1) + 1 == -nn
    print("3. D_n^+ = (nonzero) dF/dzbar; bounded holomorphic F: polynomials of degree"
          " <= m-1; RR: deg H^{n+1} + 1 = -n")


# --- 4. forms on S^2 -------------------------------------------------------------

def sphere_integral(poly, X):
    n = len(X)
    total = 0
    for mon, c in sp.Poly(sp.expand(poly), *X).terms():
        if any(k % 2 for k in mon):
            continue
        num = 2 * sp.prod([sp.gamma(sp.Rational(k + 1, 2)) for k in mon])
        total += c * num / sp.gamma(sp.Rational(sum(mon) + n, 2))
    return sp.simplify(total)


def check_forms_s2():
    X = sp.symbols("x0:3", real=True)
    xv = sp.Matrix(X)
    r2 = sum(t ** 2 for t in X)

    def grad(f):
        return sp.Matrix([sp.diff(f, t) for t in X])

    def proj(b):
        return b - xv * (xv.dot(b))

    def grad_s(f):
        return proj(grad(f))

    def curl(a):
        return sp.Matrix([sp.diff(a[2], X[1]) - sp.diff(a[1], X[2]),
                          sp.diff(a[0], X[2]) - sp.diff(a[2], X[0]),
                          sp.diff(a[1], X[0]) - sp.diff(a[0], X[1])])

    def Fc(a):
        return curl(a).dot(xv)

    def div_s(a):
        jac = sp.Matrix(3, 3, lambda i, j: sp.diff(a[i], X[j]))
        return jac.trace() - (xv.T * jac * xv)[0]

    def L(f, a, b):
        return X[a] * sp.diff(f, X[b]) - X[b] * sp.diff(f, X[a])

    def lap_s(f):
        return sum(L(L(f, a, b), a, b) for a in range(3) for b in range(a + 1, 3))

    def on_sphere_zero(expr, seed=0):
        rnd = random.Random(seed)
        for _ in range(5):
            u = [sp.Rational(rnd.randint(-5, 5), rnd.randint(1, 4)) for _ in range(2)]
            s = sum(t ** 2 for t in u)
            p = [2 * u[0] / (s + 1), 2 * u[1] / (s + 1), (s - 1) / (s + 1)]
            if sp.simplify(sp.sympify(expr).subs(dict(zip(X, p)))) != 0:
                return False
        return True

    rnd = random.Random(3)

    def rpoly(deg=3):
        mons = [m for m in itertools.product(range(deg + 1), repeat=3) if sum(m) <= deg]
        return sum(rnd.randint(-3, 3) * X[0] ** m[0] * X[1] ** m[1] * X[2] ** m[2] for m in mons)

    f, g = rpoly(), rpoly()
    a = proj(sp.Matrix([rpoly(2), rpoly(2), rpoly(2)]))
    # F(grad f) = 0, div grad = Delta_S
    assert on_sphere_zero(Fc(grad_s(f)))
    assert on_sphere_zero(div_s(grad_s(f)) - lap_s(f))
    b = -xv.cross(grad(g))
    assert on_sphere_zero(xv.dot(b))
    assert on_sphere_zero(div_s(xv.cross(grad(g))))
    assert on_sphere_zero(b.dot(b) - grad_s(g).dot(grad_s(g)))
    # adjoints (exact integrals; restrict to r = 1 by homogenising is unnecessary:
    # integrands are evaluated on the unit sphere as polynomials)
    lhs = sphere_integral(g * Fc(a), X)
    rhs = sphere_integral(b.dot(a), X)
    assert sp.simplify(lhs - rhs) == 0
    lhs = sphere_integral(grad_s(f).dot(a), X)
    rhs = sphere_integral(-f * div_s(a), X)
    assert sp.simplify(lhs - rhs) == 0
    # integration by parts for Delta_S: int f Delta_S f = -int |grad_S f|^2
    assert sp.simplify(sphere_integral(f * lap_s(f), X)
                       + sphere_integral(grad_s(f).dot(grad_s(f)), X)) == 0
    # 1 and the area form are harmonic; Gauss-Bonnet: K = 1, area 4 pi
    assert sp.Rational(1, 2) / sp.pi * 4 * sp.pi == 2
    _ = r2
    print("4. S^2 forms: F(grad f) = 0, div grad = Delta_S, adjoint of F is -x cross grad,"
          " |x cross grad g| = |grad g|, chi = 2")


# --- 5. forms on S^4 -------------------------------------------------------------

def check_forms_s4():
    B = Alg(5)
    e = B.e
    x = e[0]
    tang = [1, 2, 3, 4]

    def form(idx):
        r = mv(1, B.neg)
        for i in idx:
            r = r * (e[i] * x)
        return r

    forms = {idx: form(idx) for k in range(5) for idx in itertools.combinations(tang, k)}

    def coords(p):
        """expand an even element in the form basis."""
        out = {}
        rest = p
        for idx, F in forms.items():
            c = (F.rev() * rest).scalar()
            if c != 0:
                out[idx] = c
        assert eq(sum((c * forms[idx] for idx, c in out.items()), B.zero()), p)
        return out

    # left multiplication by e_j x = e_j ^ - iota_{e_j}
    for j in tang:
        for idx, F in forms.items():
            got = coords(e[j] * x * F)
            if j in idx:
                pos = idx.index(j)
                rest = tuple(i for i in idx if i != j)
                exp = {rest: -(-1) ** pos}
            else:
                new = tuple(sorted(idx + (j,)))
                pos = new.index(j)
                exp = {new: (-1) ** pos}
            assert got == exp
        # epsilon = (-1)^k
    for idx, F in forms.items():
        assert eq(x * F * x, (-1) ** len(idx) * F)
    rnd = random.Random(2)
    xi = sum((sp.Rational(rnd.randint(-3, 3)) * e[i] for i in tang), B.zero())
    v = sum((sp.Rational(rnd.randint(-3, 3)) * e[i] for i in tang), B.zero())
    xi2, v2 = (xi * xi).scalar(), (v * v).scalar()
    t, s = sp.symbols("t s", real=True)

    def Lm(p):
        return xi * x * p

    def Bm(p):
        return sp.I * (x * (p * v * x) * x)

    for F in forms.values():
        assert eq(Lm(Bm(F)) + Bm(Lm(F)), 0)
        assert eq(Bm(Bm(F)), -v2 * F)
        assert eq(x * Bm(F) * x, -Bm(x * F * x))
        G = Lm(F) + t * Bm(F)
        assert eq(Lm(G) + t * Bm(G), -(xi2 + t ** 2 * v2) * F)
    # B_v = i x psi v = -i (v^ + iota_v) on forms (intrinsic form)
    for j in tang:
        for idx, F in forms.items():
            if j in idx:
                q = idx.index(j)
                exp = {tuple(i for i in idx if i != j): (-1) ** q}
            else:
                new = tuple(sorted(idx + (j,)))
                exp = {new: (-1) ** new.index(j)}
            rhs = sum((-sp.I * c * forms[k] for k, c in exp.items()), B.zero())
            assert eq(sp.I * (x * F * e[j]), rhs)
    # unit case: cos s xi x + sin s B has square -1
    xu = e[1]
    vu = e[2]
    for F in forms.values():
        def M(p):
            return sp.cos(s) * xu * x * p + sp.sin(s) * sp.I * (x * (p * vu * x) * x)
        assert eq(M(M(F)), -F)
    print("5. S^4 forms: xi x = xi^ - iota_xi, epsilon = x.x = (-1)^k, B = i eps R_{vx}"
          " anticommutes, (xi x + tB)^2 = -(|xi|^2 + t^2|v|^2)")


# --- 6. characteristic numbers ---------------------------------------------------

def check_numbers():
    p1 = sp.Integer(3)                        # CP^2
    assert p1 / 3 == 1 and -p1 / 24 == sp.Rational(-1, 8)
    sig = -16                                 # K3
    p1 = 3 * sig
    assert sp.Rational(-p1, 24) == 2
    for nn in range(-5, 6):
        c1_E = -nn                            # c1(H^n) in the z-orientation
        c1_L = -(nn + 1)                      # c1(H^{n+1}) = S^+ (x) H^n
        assert c1_L + 1 == c1_E == -nn        # Todd: 1 + c1(T)/2, c1(TS^2) = 2
    print("6. CP^2: sigma 1, A-hat -1/8; K3: A-hat 2; Todd on S^2 gives -n")


if __name__ == "__main__":
    check_orientation_c1()
    check_symbol()
    check_dbar()
    check_forms_s2()
    check_forms_s4()
    check_numbers()
