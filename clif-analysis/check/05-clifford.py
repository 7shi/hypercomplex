"""Checks for article 05 (Clifford analysis in general dimension).

Cl_{n,0}(R) with generators e_0..e_{n-1} (e_a^2 = 1) is implemented directly on
bitmask blades, as in 04. The Dirac operator is D = sum e_a d_a acting on
functions of x in R^n with values in Cl_{n,0}(R).

1. D^2 = Laplacian; sum e_a e_0 e_a = (2 - n) e_0, D x^{-1} = (n - 2)/|x|^2;
   the kernel x/|x|^n is left and right monogenic (n = 2..6);
   D |x|^{2-n} = (2 - n) x/|x|^n (n >= 3), D log|x| = x^{-1} (n = 2), and
   these potentials are harmonic.
2. |S^{n-1}| = 2 pi^{n/2}/Gamma(n/2): Gaussian integral in polar coordinates,
   the monomial formula on S^{n-1}, and n vol(B^n) = |S^{n-1}| with vol(B^n)
   computed by slicing (n = 2..6).
3. Kernel on a sphere: E(eps w) w = w E(eps w) = 1/eps^{n-1} for |w| = 1.
4. Oriented boundary element: for hyperspherical coordinates the tangent
   blade dX = r_{u1} ^ ... ^ r_{u_{n-1}} satisfies dX I^{-1} = J n and
   n ^ dX = J I (J > 0 the area element) for n = 2, 3, 4, so n dS = dX I^{-1};
   for n = 2 this is n ds = I dx, for n = 3 it is dX = I n dA (article 02).
   For generic vectors w and an (n-1)-blade dX, w dX commutes with I
   (n = 2..6), so E dX I^{-1} F = I^{-1} E dX F in every dimension.
5. Cauchy-Pompeiu on the unit ball at a = 0 (exact, left and right, n = 3, 5):
   |S^{n-1}| F(0) = oint E n F dS - int E DF dV with E = x/|x|^n.
6. Left monogenic polynomials and the mean value property with a symbolic center
   (exact, n = 5); Cauchy's formula at an off-center point (numeric, n = 3).
7. Residues: (1/|S^{n-1}|) oint n (E c) dS = c, and the residue of
   d_k E is 0 (n = 3, 5).
8. Fueter-Sce: in Cl_{n,0} with q = e_0 x, Fueter-type operator e_0 D,
   imaginary units h_l = e_0 e_l (l = 1..n-1):
   e_0 D (A + w B) = A_0 - B_r - (n-2) B/r + w (A_r + B_0) (both sides);
   axial Laplacian g_00 + g_rr + (n-2) g_r / r;
   Delta(v/r) = (v_00 + v_rr)/r + (n-4)(v_r/r^2 - v/r^3);
   ftilde(z^k) = q^k, ftilde(1/z) = q^{-1}, e_0 D ftilde = -(n-2) v/r;
   Delta^{(n-2)/2} ftilde is left and right regular for e_0 D (n = 4, 6);
   qbar/|q|^n = E e_0 is left and right regular for e_0 D, while
   (E e_0) D != 0 (right regularity differs from right monogenicity);
   Delta (x_j g(|x|^2)) = x_j (4 rho g'' + (2n+4) g'), which gives
   Delta^k q^{-1} = c_n qbar/|q|^n, k = (n-2)/2, c_n = (-4)^k (k!)^2
   (c_4 = -4, c_6 = 64; closed form checked for n = 4..12).
9. D = d - delta: with the identification of k-vectors and k-forms and
   *F = reverse(F) I, D ^ F is the exterior derivative and
   D . F = -delta F, delta = (-1)^{n(k+1)+1} * d * (n = 3, 4).
"""

import itertools
import math
import random

import numpy as np
import sympy as sp

# ---------------------------------------------------------------- Cl_{n,0}


def blade_mul(a, b):
    """Product of basis blades (bitmasks) in Cl_{n,0}: returns (sign, mask)."""
    s = 0
    t = a >> 1
    while t:
        s += bin(t & b).count("1")
        t >>= 1
    return (-1 if s & 1 else 1), a ^ b


def grade(m):
    return bin(m).count("1")


class MV:
    def __init__(self, d=None):
        self.d = {k: v for k, v in (d or {}).items() if v != 0}

    def __add__(self, o):
        o = mv(o)
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, 0) + v
        return MV(d)

    __radd__ = __add__

    def __neg__(self):
        return MV({k: -v for k, v in self.d.items()})

    def __sub__(self, o):
        return self + (-mv(o))

    def __rsub__(self, o):
        return mv(o) - self

    def __mul__(self, o):
        o = mv(o)
        d = {}
        for a, x in self.d.items():
            for b, y in o.d.items():
                s, m = blade_mul(a, b)
                d[m] = d.get(m, 0) + s * x * y
        return MV(d)

    def __rmul__(self, o):
        return mv(o) * self

    def __truediv__(self, c):
        return MV({k: v / c for k, v in self.d.items()})

    def map(self, f):
        return MV({k: f(v) for k, v in self.d.items()})

    def grade(self, k):
        return MV({m: v for m, v in self.d.items() if grade(m) == k})

    def rev(self):
        return MV({m: (-v if (grade(m) * (grade(m) - 1) // 2) % 2 else v) for m, v in self.d.items()})

    def scalar(self):
        return self.d.get(0, 0)


def mv(o):
    return o if isinstance(o, MV) else MV({0: sp.sympify(o)})


def zero_expr(v):
    v = sp.cancel(sp.together(sp.expand(v)))
    return v == 0 or sp.simplify(v) == 0


def eq(A, B):
    return all(zero_expr(v) for v in (mv(A) - mv(B)).d.values())


class Alg:
    """Cl_{n,0} with coordinates x_0..x_{n-1} and the Dirac operator."""

    def __init__(self, n, name="x"):
        self.n = n
        self.e = [MV({1 << a: 1}) for a in range(n)]
        self.X = sp.symbols(f"{name}0:{n}", real=True)
        self.x = sum((self.X[a] * self.e[a] for a in range(n)), MV())
        self.r2 = sum(v**2 for v in self.X)
        self.I = MV({(1 << n) - 1: 1})

    def d(self, F, a):
        return F.map(lambda v: sp.diff(v, self.X[a]))

    def D(self, F):
        return sum((self.e[a] * self.d(F, a) for a in range(self.n)), MV())

    def Dr(self, F):
        return sum((self.d(F, a) * self.e[a] for a in range(self.n)), MV())

    def lap(self, F):
        return sum((self.d(self.d(F, a), a) for a in range(self.n)), MV())

    def rnd(self, seed, deg=2, grades=None):
        rnd = random.Random(seed)
        mons = [m for m in itertools.product(range(deg + 1), repeat=self.n) if sum(m) <= deg]
        blades = [b for b in range(1 << self.n) if grades is None or grade(b) in grades]
        return MV({b: sum(rnd.randint(-3, 3) * sp.prod([self.X[a] ** p for a, p in enumerate(m)])
                          for m in mons) for b in blades})


# ---------------------------------------------------------------- 1.
for n in range(2, 7):
    C = Alg(n)
    if n <= 4:
        F = C.rnd(n, grades=None)
        assert eq(C.D(C.D(F)), C.lap(F))
    assert eq(sum((C.e[a] * C.e[0] * C.e[a] for a in range(n)), MV()), (2 - n) * C.e[0])
    assert eq(C.D(C.x / C.r2), mv((n - 2) / C.r2))
    kern = C.x / C.r2 ** sp.Rational(n, 2)
    assert eq(C.D(kern), 0) and eq(C.Dr(kern), 0)
    if n == 2:
        pot = mv(sp.log(C.r2) / 2)
        assert eq(C.D(pot), C.x / C.r2)
    else:
        pot = mv(C.r2 ** sp.Rational(2 - n, 2))
        assert eq(C.D(pot), (2 - n) * kern)
    assert eq(C.lap(pot), 0)
print("1. D^2 = Delta, sum e_a e0 e_a = (2-n) e0, D x^{-1} = (n-2)/|x|^2, x/|x|^n monogenic (both sides), "
      "D|x|^{2-n} = (2-n) x/|x|^n, D log|x| = x^{-1} (n = 2..6)")

# ---------------------------------------------------------------- 2.
r, t, ep = sp.symbols("r t epsilon", positive=True)


def area(n):
    return 2 * sp.pi ** sp.Rational(n, 2) / sp.gamma(sp.Rational(n, 2))


def int_sphere(c, W):
    """Exact integral over S^{n-1} of a polynomial in the unit normal W."""
    n = len(W)
    poly = sp.Poly(sp.expand(c), *W)
    tot = 0
    for exps, coef in poly.terms():
        if any(p % 2 for p in exps):
            continue
        tot += coef * 2 * sp.prod([sp.gamma(sp.Rational(p + 1, 2)) for p in exps]) \
            / sp.gamma(sp.Rational(sum(exps) + n, 2))
    return sp.simplify(tot)


vol = {1: sp.Integer(2)}
for n in range(2, 7):
    # Gaussian: pi^{n/2} = int e^{-|x|^2} dV = |S^{n-1}| int_0^oo r^{n-1} e^{-r^2} dr
    radial = sp.integrate(r ** (n - 1) * sp.exp(-r**2), (r, 0, sp.oo))
    assert sp.simplify(radial - sp.gamma(sp.Rational(n, 2)) / 2) == 0
    assert sp.simplify(sp.pi ** sp.Rational(n, 2) / radial - area(n)) == 0
    W = sp.symbols(f"w0:{n}", real=True)
    assert sp.simplify(int_sphere(1, W) - area(n)) == 0
    vol[n] = sp.simplify(sp.integrate(vol[n - 1] * (1 - t**2) ** sp.Rational(n - 1, 2), (t, -1, 1)))
    assert sp.simplify(n * vol[n] - area(n)) == 0
print("2. |S^{n-1}| = 2 pi^{n/2}/Gamma(n/2) = n vol(B^n):",
      ", ".join(f"n={n}: {sp.nsimplify(area(n))}" for n in range(2, 7)))

# ---------------------------------------------------------------- 3.
for n in range(2, 7):
    C = Alg(n)
    W = sp.symbols(f"w0:{n}", real=True)
    w = sum((W[a] * C.e[a] for a in range(n)), MV())
    Ew = ep * w / ep**n  # E(eps w) with |eps w| = eps
    s = sum(v**2 for v in W)
    assert eq(Ew * w, s / ep ** (n - 1)) and eq(w * Ew, s / ep ** (n - 1))
print("3. E(eps n) n = n E(eps n) = 1/eps^{n-1} (n = 2..6)")

# ---------------------------------------------------------------- 4.
# hyperspherical coordinates, oriented so that n ^ dX has the orientation of I
U = sp.symbols("u0:3", real=True)


def sphere_point(n):
    if n == 2:
        th = U[0]
        return [sp.cos(th), sp.sin(th)], [th]
    # x_{n-1} = cos u0, remaining coordinates = sin u0 * (point of S^{n-2})
    p, us = sphere_point(n - 1)
    pts = [sp.sin(U[n - 2]) * c for c in p] + [sp.cos(U[n - 2])]
    return pts, us + [U[n - 2]]


for n in (2, 3, 4):
    C = Alg(n)
    pts, us = sphere_point(n)
    nvec = sum((pts[a] * C.e[a] for a in range(n)), MV())
    tang = [sum((sp.diff(pts[a], u) * C.e[a] for a in range(n)), MV()) for u in us]
    prod_ = mv(1)
    for v in tang:
        prod_ = prod_ * v
    dX = prod_.grade(n - 1)
    # area element for u0 in (0, 2pi), other parameters in (0, pi)
    J = sp.prod([sp.sin(U[k]) ** k for k in range(1, n - 1)])
    lhs = dX * C.I.rev()  # I^{-1} = reverse(I)
    sgn = 1 if eq(lhs, J * nvec) else -1
    assert eq(lhs, sgn * J * nvec)
    assert eq((nvec * dX).grade(n), sgn * J * C.I)
    print(f"   n={n}: parameters {us}: dX I^-1 = {'+' if sgn > 0 else '-'}J n, n^dX = {'+' if sgn > 0 else '-'}J I")
C = Alg(2)
dx = sp.symbols("dx0:2", real=True)
dxv = dx[0] * C.e[0] + dx[1] * C.e[1]
assert eq(dxv * C.I.rev(), C.I * dxv)  # n ds = dx I^{-1} = I dx
C = Alg(3)
nn = sp.symbols("m0:3", real=True)
nv = sum((nn[a] * C.e[a] for a in range(3)), MV())
dXv = C.I * nv
assert eq(dXv * C.I.rev(), nv)  # dX = I n dA  <=>  n dA = dX I^{-1}
for n in range(2, 7):
    C = Alg(n)
    W = sp.symbols(f"w0:{n}", real=True)
    w = sum((W[a] * C.e[a] for a in range(n)), MV())
    blade = mv(1)
    for i in range(n - 1):
        blade = blade * sum((sp.Symbol(f"t{i}_{a}", real=True) * C.e[a] for a in range(n)), MV())
    dXg = blade.grade(n - 1)
    assert eq(w * dXg * C.I, C.I * w * dXg)
print("4. n ^ dX = J I  <=>  n dS = dX I^{-1}; n = 2: n ds = I dx, n = 3: dX = I n dA; "
      "E dX commutes with I (n = 2..6)")

# ---------------------------------------------------------------- 5.


def sphere_mv_int(F, W):
    return F.map(lambda v: int_sphere(v, W))


for n in (3, 5):
    C = Alg(n)
    W = sp.symbols(f"w0:{n}", real=True)
    w = sum((W[a] * C.e[a] for a in range(n)), MV())
    sub_b = {C.X[a]: W[a] for a in range(n)}
    sub_v = {C.X[a]: r * W[a] for a in range(n)}
    F = C.rnd(10 + n, deg=2) if n == 3 else C.rnd(10 + n, deg=1)
    F0 = F.map(lambda v: v.subs({x: 0 for x in C.X}))
    # boundary: on the unit sphere E n = 1, dS = dOmega
    bnd = sphere_mv_int(F.map(lambda v: v.subs(sub_b)), W)
    # volume: E(r w) DF(r w) r^{n-1} dr dOmega = w DF(r w) dr dOmega
    vl = (w * C.D(F).map(lambda v: v.subs(sub_v))).map(lambda v: sp.integrate(sp.expand(v), (r, 0, 1)))
    assert eq(bnd - sphere_mv_int(vl, W), area(n) * F0)
    vr = (C.Dr(F).map(lambda v: v.subs(sub_v)) * w).map(lambda v: sp.integrate(sp.expand(v), (r, 0, 1)))
    assert eq(bnd - sphere_mv_int(vr, W), area(n) * F0)
print("5. Cauchy-Pompeiu on the unit ball at 0, left and right (n = 3, 5)")

# ---------------------------------------------------------------- 6.
n = 5
C = Alg(n)
Z = [C.X[l] - C.X[0] * C.e[0] * C.e[l] for l in range(1, n)]
for z in Z:
    assert eq(C.D(z), 0)
P = Z[0] * Z[1] + Z[1] * Z[0]
assert eq(C.D(P), 0)
assert not eq(C.D(Z[0] * Z[1]), 0)
A = sp.symbols(f"a0:{n}", real=True)
R = sp.symbols("R", positive=True)
W = sp.symbols(f"w0:{n}", real=True)
for Fm in (Z[0], P, P * C.e[2] + Z[3]):
    shifted = Fm.map(lambda v: sp.expand(v.subs({C.X[a]: A[a] + R * W[a] for a in range(n)}, simultaneous=True)))
    avg = sphere_mv_int(shifted, W) / area(n)
    assert eq(avg, Fm.map(lambda v: v.subs({C.X[a]: A[a] for a in range(n)}, simultaneous=True)))
print("6. x_l - x_0 e_0 e_l and symmetrized products are left monogenic; mean value property (n = 5, exact)")

# numeric: Cauchy's formula at an off-center point in the unit ball of R^3


def np_mul(a, b, n):
    out = np.zeros(1 << n)
    for i in range(1 << n):
        if a[i] == 0:
            continue
        for j in range(1 << n):
            if b[j] == 0:
                continue
            s, m = blade_mul(i, j)
            out[m] += s * a[i] * b[j]
    return out


def np_vec(v, n):
    out = np.zeros(1 << n)
    for a in range(n):
        out[1 << a] = v[a]
    return out


n = 3
b = np.array([1.2, -0.9, 0.7])  # singularity outside the unit ball
c = np.array([0.3, -1.0, 0.5, 2.0, 0.0, 1.5, -0.7, 0.4])  # right coefficient


def F_num(x):
    y = x - b
    return np_mul(np_vec(y / np.linalg.norm(y) ** 3, n), c, n)


a = np.array([0.25, 0.1, -0.35])
Nt, Np = 60, 120
gl, glw = np.polynomial.legendre.leggauss(Nt)
tot = np.zeros(8)
for ct, wt in zip(gl, glw):
    st = math.sqrt(1 - ct * ct)
    for k in range(Np):
        ph = 2 * math.pi * k / Np
        nv = np.array([st * math.cos(ph), st * math.sin(ph), ct])
        y = nv - a
        E = np_vec(y / np.linalg.norm(y) ** 3, n)
        tot += wt * (2 * math.pi / Np) * np_mul(np_mul(E, np_vec(nv, n), n), F_num(nv), n)
assert np.allclose(tot / (4 * math.pi), F_num(a), atol=1e-10)
print("   Cauchy's formula at an off-center point of the unit ball (n = 3, numeric)")

# ---------------------------------------------------------------- 7.
for n in (3, 5):
    C = Alg(n)
    W = sp.symbols(f"w0:{n}", real=True)
    w = sum((W[a] * C.e[a] for a in range(n)), MV())
    cst = MV({m: sp.Integer(m % 5 - 2) for m in range(1 << n)})
    kern = C.x / C.r2 ** sp.Rational(n, 2)
    rho = sp.sqrt(C.r2)

    def on_sphere(F):
        # restrict to |x| = eps, x = eps w: substitute |x| first, then x
        return F.map(lambda v: sp.expand(v.subs(rho, ep).subs(C.r2, ep**2)
                                          .subs({C.X[a]: ep * W[a] for a in range(n)}, simultaneous=True)))

    def res(F):
        return sphere_mv_int(w * on_sphere(F), W).map(lambda v: sp.simplify(v * ep ** (n - 1))) / area(n)

    assert eq(res(kern * cst), cst)
    for k in range(n):
        assert eq(res(C.d(kern, k) * cst), 0)
        for l in range(n):
            assert eq(res(C.d(C.d(kern, k), l)), 0)
print("7. Res(E c) = c, Res(d_k E) = Res(d_k d_l E) = 0 (n = 3, 5)")

# ---------------------------------------------------------------- 8.
x0s, rs = sp.symbols("X R", positive=True)
zc = sp.Symbol("z")
cvals = {}
for n in (4, 6):
    C = Alg(n)
    x0 = C.X[0]
    h = [C.e[0] * C.e[l] for l in range(1, n)]
    q = C.e[0] * C.x
    qbar = C.x * C.e[0]
    assert eq(q, x0 + sum((C.X[l] * h[l - 1] for l in range(1, n)), MV()))
    rr = sp.sqrt(sum(C.X[l] ** 2 for l in range(1, n)))
    qvec = sum((C.X[l] * h[l - 1] for l in range(1, n)), MV())
    omega = qvec / rr
    FD = lambda F: C.e[0] * C.D(F)
    FDr = lambda F: C.Dr(F * C.e[0])  # right action of e0 D is (F e0) D

    def at(expr):
        return expr.subs({x0s: x0, rs: rr})

    # axial derivative formula, concrete A, B
    for A_, B_ in ((x0s**2 * rs + rs**3, x0s * rs**2 + x0s**3), (sp.exp(x0s) * rs**2, sp.sin(x0s) * rs)):
        F = mv(at(A_)) + omega * at(B_)
        rhs = mv(at(sp.diff(A_, x0s) - sp.diff(B_, rs) - (n - 2) * B_ / rs)) \
            + omega * at(sp.diff(A_, rs) + sp.diff(B_, x0s))
        assert eq(FD(F), rhs) and eq(FDr(F), rhs)
        g = at(A_)
        assert zero_expr(sum(sp.diff(g, v, 2) for v in C.X)
                         - at(sp.diff(A_, x0s, 2) + sp.diff(A_, rs, 2) + (n - 2) / rs * sp.diff(A_, rs)))
    # Delta(v/r) in axial form
    vf = sp.Function("v")(x0s, rs)
    g = vf / rs
    lap_ax = sp.diff(g, x0s, 2) + sp.diff(g, rs, 2) + (n - 2) / rs * sp.diff(g, rs)
    expect = (sp.diff(vf, x0s, 2) + sp.diff(vf, rs, 2)) / rs \
        + (n - 4) * (sp.diff(vf, rs) / rs**2 - vf / rs**3)
    assert sp.simplify(lap_ax - expect) == 0
    # ftilde and Sce's theorem
    k = (n - 2) // 2

    def lapk(F):
        for _ in range(k):
            F = C.lap(F)
        return F

    q2 = q * q
    q3 = q2 * q
    qinv = qbar / C.r2
    assert eq(qinv * q, 1)
    for f, qf in ((zc**2, q2), (zc**3, q3), (zc**4, q3 * q), (1 / zc, qinv)):
        fz = sp.expand_complex(f.subs(zc, x0s + sp.I * rs))
        u_, v_ = sp.re(fz), sp.im(fz)
        ft = mv(at(u_)) + omega * at(v_)
        assert eq(ft, qf), f
        assert eq(FD(ft), mv(at(-(n - 2) * v_ / rs))), f
        if f == 1 / zc:
            continue  # handled below through the radial formula
        lf = lapk(qf)
        assert eq(FD(lf), 0) and eq(FDr(lf), 0), f
        if n == 6 and f == zc**3:
            vr = at(v_ / rs)
            assert not zero_expr(sum(sp.diff(vr, v, 2) for v in C.X))
    # Delta (x_j g(rho)) = x_j (4 rho g'' + (2n + 4) g') with rho = |x|^2
    rho = sp.Symbol("rho", positive=True)
    for g in (rho ** sp.Rational(-3, 2), sp.exp(rho), rho**2):
        lhs_ = sum(sp.diff(C.X[1] * g.subs(rho, C.r2), v, 2) for v in C.X)
        Lg = 4 * rho * sp.diff(g, rho, 2) + (2 * n + 4) * sp.diff(g, rho)
        assert zero_expr(lhs_ - C.X[1] * Lg.subs(rho, C.r2))
    g = 1 / rho  # q^{-1} = qbar g(rho)
    for _ in range(k):
        g = sp.simplify(4 * rho * sp.diff(g, rho, 2) + (2 * n + 4) * sp.diff(g, rho))
    cn = sp.simplify(g * rho ** sp.Rational(n, 2))
    assert cn.is_number
    E = qbar / C.r2 ** sp.Rational(n, 2)  # left and right monogenic by 1.
    assert eq(FD(E), 0) and eq(FDr(E), 0)
    assert not eq(C.Dr(E), 0)  # right regular, but not right monogenic for D
    cvals[n] = cn
assert cvals[4] == -4 and cvals[6] == 64
rho = sp.Symbol("rho", positive=True)
for n in range(4, 13, 2):
    k = (n - 2) // 2
    g = 1 / rho
    for _ in range(k):
        g = sp.simplify(4 * rho * sp.diff(g, rho, 2) + (2 * n + 4) * sp.diff(g, rho))
    assert sp.simplify(g * rho ** sp.Rational(n, 2) - (-4) ** k * sp.factorial(k) ** 2) == 0
print("8. e0 D (A + wB) = A0 - B_r - (n-2)B/r + w(A_r + B0); Delta(v/r) = (v00+vrr)/r + (n-4)(v_r/r^2 - v/r^3); "
      "Delta^{(n-2)/2} ftilde left/right regular (e0 D) for z^2, z^3, z^4, 1/z; (E e0) D != 0;",
      ", ".join(f"Delta^{(n - 2) // 2} q^-1 = {c} qbar/|q|^{n}" for n, c in cvals.items()),
      "; c_n = (-4)^k (k!)^2 (n = 4..12)")

# ---------------------------------------------------------------- 9.
for n in (3, 4):
    C = Alg(n)
    star = lambda F: F.rev() * C.I
    for kg in range(n + 1):
        F = C.rnd(30 + kg, deg=2, grades={kg})
        DF = C.D(F)
        wedge = DF.grade(kg + 1)
        dot = DF.grade(kg - 1) if kg > 0 else MV()
        # exterior derivative: d(f dx^A) = sum d_a f dx^a ^ dx^A
        dF = MV()
        for m, v in F.d.items():
            for a in range(n):
                if not (m >> a) & 1:
                    s, mm = blade_mul(1 << a, m)
                    dF = dF + MV({mm: s * sp.diff(v, C.X[a])})
        assert eq(wedge, dF)
        # codifferential via Hodge star on (n-k)-vectors
        sF = star(F)
        dsF = MV()
        for m, v in sF.d.items():
            for a in range(n):
                if not (m >> a) & 1:
                    s, mm = blade_mul(1 << a, m)
                    dsF = dsF + MV({mm: s * sp.diff(v, C.X[a])})
        delta = (-1) ** (n * (kg + 1) + 1) * star(dsF)
        assert eq(dot, -delta)
        assert eq(DF, dF - delta)
print("9. D ^ F = dF, D . F = -delta F, D = d - delta (n = 3, 4)")
