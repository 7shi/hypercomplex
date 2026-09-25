"""Checks for article 04 (quaternionic analysis, Fueter regular functions).

Cl_{4,0}(R) with generators e0..e3 (e_a^2 = 1) uses the bitmask blades of
common.clifford. The quaternion units are i = e0 e1, j = e0 e2, k = e0 e3, the
quaternion variable is q = e0 x, and the Fueter operator is
e0 D = d0 + i d1 + j d2 + k d3 with the vector Dirac operator D = sum e_a d_a.
Quaternion-valued computations use sympy's Quaternion.

1. i, j, k: squares -1, pairwise anticommuting, ijk = -I4 (I4 = e0e1e2e3),
   I4^2 = 1, I4 central in the even part; on the half (1+I4)/2 the Hamilton
   relations ij = k etc. hold, on (1-I4)/2 ijk = +1.
2. e0 D = d0 + i d1 + j d2 + k d3, Dbar := D e0 = d0 - i d1 - j d2 - k d3,
   Dbar (e0 D) = (e0 D) Dbar = D^2 = Laplacian; (G e0) D = G (right Fueter).
3. Dimension counting: sum_a e_a e0 e_a = (2 - n) e0, D(e0 x) = (2 - n) e0 in
   Cl_{n,0} (n = 2, 3, 4), D x^{-1} = (n - 2)/|x|^2, D (x/|x|^n) = 0 (both sides).
   q^{-1} = x^{-1} e0 and E := qbar/|q|^4 = (x/|x|^4) e0.
4. Quaternion level: Dq = -2, D qbar = 4, D q^{-1} = 2/|q|^2, D q^2 = -4 x0,
   E left and right regular, Laplacian of components of regular functions = 0.
5. Products: zeta_l = x_l - x0 e_l regular (both sides); product rule
   D(fg) = (Df)g + f(Dg) + sum_a [eps_a, f] d_a g; D(zeta1 zeta2) = 2 x0 k,
   symmetrized product regular, zeta1^2 regular, zeta1 j left but not right
   regular, zeta1 o zeta2 = -x2 i with D = k.
6. Fueter's theorem: D(A + w B) = A0 - B_r - 2B/r + w (A_r + B0) for w = qvec/r,
   axial Laplacian of v/r = (v00 + vrr)/r; for f = z, z^2, z^3, 1/z, exp z:
   D ftilde = -2 v/r and Delta ftilde is left and right regular; ftilde(z^n)
   = q^n, ftilde(1/z) = q^{-1}, Delta q^{-1} = -4 E, Delta q^2 = -4,
   Delta q^3 = -4 (2q + qbar).
7. Two-sided fundamental theorem with quaternion coefficients on the unit
   4-cube (exact): int [(G D) F + G (D F)] dV = oint G n F dS.
8. Kernel on a sphere: E(eps n) n = n E(eps n) = 1/eps^3, |S^3| = 2 pi^2.
   Sphere integrals of polynomials in the unit normal use the exact monomial
   formula for S^3.
9. Cauchy-Pompeiu on the unit ball at a = 0 (exact, left and right):
   F(0) = (1/2pi^2)[oint E n F dS - int E DF dV].
10. Cauchy-Fueter formula at an off-center point (numeric) and the mean value
    property (exact) for regular functions.
"""

import itertools
import random

import numpy as np
import sympy as sp
from sympy.algebras.quaternion import Quaternion as Q

from common.clifford import MV, eq, mv


e = [MV({1 << a: 1}) for a in range(4)]
e0, e1, e2, e3 = e
i_, j_, k_ = e0 * e1, e0 * e2, e0 * e3
I4 = e0 * e1 * e2 * e3
X = sp.symbols("x0:4", real=True)
x0, x1, x2, x3 = X

# 1. quaternion units in Cl^0_{4,0}
for u in (i_, j_, k_):
    assert eq(u * u, -1)
for u, w in ((i_, j_), (j_, k_), (k_, i_)):
    assert eq(u * w, -(w * u))
assert eq(i_ * j_ * k_, -I4)
assert eq(I4 * I4, 1)
even = [MV({m: 1}) for m in range(16) if bin(m).count("1") % 2 == 0]
for b in even:
    assert eq(I4 * b, b * I4)
assert eq(e0 * I4, -(I4 * e0))
Pp, Pm = (1 + I4) / 2, (1 - I4) / 2
assert eq(Pp * Pp, Pp) and eq(Pp * Pm, 0)
assert eq(i_ * j_ * Pp, k_ * Pp) and eq(j_ * k_ * Pp, i_ * Pp) and eq(k_ * i_ * Pp, j_ * Pp)
assert eq(i_ * j_ * k_ * Pp, -Pp) and eq(i_ * j_ * k_ * Pm, Pm)
print("1. i=e0e1, j=e0e2, k=e0e3: ijk = -I4, I4 central in even part, Hamilton on (1+I4)/2")


# 2. operators
def rnd_even(seed, deg=2):
    rnd = random.Random(seed)
    mons = [m for m in itertools.product(range(deg + 1), repeat=4) if sum(m) <= deg]
    return MV({b: sum(rnd.randint(-3, 3) * sp.prod([X[a] ** p for a, p in enumerate(m)]) for m in mons)
               for b in range(16) if bin(b).count("1") % 2 == 0})


def dmv(F, a):
    return F.map(lambda v: sp.diff(v, X[a]))


def Dv(F):  # vector Dirac operator from the left
    return sum((e[a] * dmv(F, a) for a in range(4)), MV())


def Dv_r(F):  # from the right
    return sum((dmv(F, a) * e[a] for a in range(4)), MV())


def lap_mv(F):
    return sum((dmv(dmv(F, a), a) for a in range(4)), MV())


eps_mv = [mv(1), i_, j_, k_]
F = rnd_even(1)
G = rnd_even(2)
DF = sum((eps_mv[a] * dmv(F, a) for a in range(4)), MV())
DbF = dmv(F, 0) - sum((eps_mv[a] * dmv(F, a) for a in range(1, 4)), MV())
assert eq(e0 * Dv(F), DF)
assert eq(Dv(e0 * F), DbF)
assert eq(Dv(Dv(F)), lap_mv(F))
DbDF = dmv(DF, 0) - sum((eps_mv[a] * dmv(DF, a) for a in range(1, 4)), MV())
DDbF = dmv(DbF, 0) + sum((eps_mv[a] * dmv(DbF, a) for a in range(1, 4)), MV())
assert eq(DbDF, lap_mv(F)) and eq(DDbF, lap_mv(F))
GD = sum((dmv(G, a) * eps_mv[a] for a in range(4)), MV())
assert eq(Dv_r(G * e0), GD)
print("2. e0 D = d0 + i d1 + j d2 + k d3, Dbar = D e0, Dbar D = D Dbar = D^2 = Delta, (G e0) D = G D_Fueter")


# 3. dimension counting in Cl_{n,0}
def check_dim(n):
    E = [MV({1 << a: 1}) for a in range(n)]
    Y = sp.symbols(f"y0:{n}", real=True)
    x = sum((Y[a] * E[a] for a in range(n)), MV())
    r2 = sum(y**2 for y in Y)
    D_ = lambda F_: sum((E[a] * F_.map(lambda v: sp.diff(v, Y[a])) for a in range(n)), MV())
    Dr_ = lambda F_: sum((F_.map(lambda v: sp.diff(v, Y[a])) * E[a] for a in range(n)), MV())
    assert eq(sum((E[a] * E[0] * E[a] for a in range(n)), MV()), (2 - n) * E[0])
    assert eq(D_(E[0] * x), (2 - n) * E[0])
    assert eq(E[0] * D_(E[0] * x), 2 - n)
    assert eq(D_(x / r2), mv((n - 2) / r2))
    kern = x / r2 ** sp.Rational(n, 2)
    assert eq(D_(kern), 0) and eq(Dr_(kern), 0)


for n in (2, 3, 4):
    check_dim(n)
xv = sum((X[a] * e[a] for a in range(4)), MV())
r2 = sum(v**2 for v in X)
qv = e0 * xv
assert eq(qv, x0 + x1 * i_ + x2 * j_ + x3 * k_)
qbar = x0 - x1 * i_ - x2 * j_ - x3 * k_
assert eq(xv * e0, qbar)  # reversion of e0 x
assert eq(qv * qbar, r2)
assert eq((xv / r2) * e0 * qv, 1)  # q^{-1} = x^{-1} e0
assert eq((xv / r2**2) * e0, qbar / r2**2)  # E = (x/|x|^4) e0
print("3. sum e_a e0 e_a = (2-n) e0, D(e0 x) = (2-n) e0, D x^{-1} = (n-2)/|x|^2, D(x/|x|^n) = 0 (n = 2,3,4)")

# ---------------------------------------------------------------- quaternions
EPS = [Q(1, 0, 0, 0), Q(0, 1, 0, 0), Q(0, 0, 1, 0), Q(0, 0, 0, 1)]
ZQ = Q(0, 0, 0, 0)


def comps(f):
    return (f.a, f.b, f.c, f.d)


def qmap(f, g):
    return Q(*[g(c) for c in comps(f)])


def d(f, a):
    return qmap(f, lambda c: sp.diff(c, X[a]))


def DL(f):
    return sum((EPS[a] * d(f, a) for a in range(4)), ZQ)


def DR(f):
    return sum((d(f, a) * EPS[a] for a in range(4)), ZQ)


def lap(f):
    return qmap(f, lambda c: sum(sp.diff(c, v, 2) for v in X))


def qzero(f):
    return all(sp.simplify(sp.expand(c)) == 0 for c in comps(f))


def qeq(f, g):
    return qzero(f - g)


def scal(s):
    return Q(s, 0, 0, 0)


q = Q(*X)
qb = Q(x0, -x1, -x2, -x3)
rq2 = x0**2 + x1**2 + x2**2 + x3**2
qinv = qmap(qb, lambda c: c / rq2)
Ek = qmap(qb, lambda c: c / rq2**2)

# 4.
assert qeq(DL(q), scal(-2))
assert qeq(DL(qb), scal(4))
assert qeq(DL(qinv), scal(2 / rq2))
assert qeq(DL(q * q), scal(-4 * x0))
assert qzero(DL(Ek)) and qzero(DR(Ek))
print("4. Dq = -2, D qbar = 4, D q^{-1} = 2/|q|^2, D q^2 = -4 x0, E = qbar/|q|^4 left and right regular")

# 5. products
zeta = [Q(X[l], 0, 0, 0) - x0 * EPS[l] for l in (1, 2, 3)]
for z_ in zeta:
    assert qzero(DL(z_)) and qzero(DR(z_))


def rnd_q(seed, deg=2):
    rnd = random.Random(seed)
    mons = [m for m in itertools.product(range(deg + 1), repeat=4) if sum(m) <= deg]
    return Q(*[sum(rnd.randint(-3, 3) * sp.prod([X[a] ** p for a, p in enumerate(m)]) for m in mons) for _ in range(4)])


f, g = rnd_q(3), rnd_q(4)
rhs = DL(f) * g + f * DL(g) + sum(((EPS[a] * f - f * EPS[a]) * d(g, a) for a in range(4)), ZQ)
assert qeq(DL(f * g), rhs)
z1, z2 = zeta[0], zeta[1]
assert qeq(DL(z1 * z2), Q(0, 0, 0, 2 * x0))
assert qzero(DL(z1 * z2 + z2 * z1)) and qzero(DR(z1 * z2 + z2 * z1))
assert qzero(DL(z1 * z1))
assert qzero(DL(z1 * EPS[2])) and qeq(DR(z1 * EPS[2]), Q(0, 0, 0, -2))
# composition zeta1 o zeta2: substitute the components of zeta2 for x0..x3
comp = qmap(z1, lambda c: c.subs(dict(zip(X, comps(z2))), simultaneous=True))
assert qeq(comp, Q(0, -x2, 0, 0)) and qeq(DL(comp), Q(0, 0, 0, 1))
for h in (z1 * z2 + z2 * z1, Ek, z1 * EPS[2]):
    assert qzero(lap(h))
print("5. zeta_l regular; D(fg) = (Df)g + f(Dg) + sum [eps_a, f] d_a g; D(zeta1 zeta2) = 2 x0 k; "
      "symmetrized regular; zeta1 j left only; zeta1 o zeta2 = -x2 i not regular")

# 6. Fueter's theorem
rs = sp.sqrt(x1**2 + x2**2 + x3**2)
w = [x1 / rs, x2 / rs, x3 / rs]
t0, tr = sp.symbols("t0 tr", real=True)


def axial(A, B):
    """A(x0, r) + w B(x0, r) as a quaternion function of x0..x3."""
    s = {t0: x0, tr: rs}
    Ax, Bx = A.subs(s), B.subs(s)
    return Q(Ax, w[0] * Bx, w[1] * Bx, w[2] * Bx)


A_ = t0**2 * tr**3 + sp.sin(tr) * t0 + sp.exp(t0) * tr**2
B_ = tr**2 * t0 + tr**5 - sp.cos(t0) * tr
lhs = DL(axial(A_, B_))
expect = axial(A_.diff(t0) - B_.diff(tr) - 2 * B_ / tr, A_.diff(tr) + B_.diff(t0))
assert qeq(lhs, expect) and qeq(DR(axial(A_, B_)), expect)
vv = sp.Function("v")(t0, tr)
g_ = vv / tr
axlap = g_.diff(t0, 2) + g_.diff(tr, 2) + 2 / tr * g_.diff(tr)
assert sp.simplify(axlap - (vv.diff(t0, 2) + vv.diff(tr, 2)) / tr) == 0

Xc, Yc = sp.symbols("X Y", real=True)
zc = Xc + sp.I * Yc
cases = {
    "z": (zc, q),
    "z^2": (zc**2, q * q),
    "z^3": (zc**3, q * q * q),
    "1/z": (1 / zc, qinv),
    "exp z": (sp.exp(zc), None),
}
for name, (fz, qf) in cases.items():
    fz = sp.expand_complex(fz)
    u, v = sp.re(fz), sp.im(fz)
    assert sp.simplify(sp.diff(u, Xc) - sp.diff(v, Yc)) == 0
    ft = axial(u.subs({Xc: t0, Yc: tr}), v.subs({Xc: t0, Yc: tr}))
    if qf is not None:
        assert qeq(ft, qf), name
    assert qeq(DL(ft), scal(-2 * v.subs({Xc: x0, Yc: rs}) / rs)), name
    lf = lap(ft)
    assert qzero(DL(lf)) and qzero(DR(lf)), name
assert qeq(lap(qinv), qmap(Ek, lambda c: -4 * c))
assert qeq(lap(q * q), scal(-4))
assert qeq(lap(q * q * q), qmap(2 * q + qb, lambda c: -4 * c))
assert qzero(DL(2 * q + qb))
print("6. D(A + wB) = A0 - B_r - 2B/r + w(A_r + B0); Delta(v/r) = (v00+vrr)/r; "
      "Fueter for z, z^2, z^3, 1/z, exp z; Delta q^{-1} = -4E, Delta q^3 = -4(2q + qbar)")

# 7. two-sided fundamental theorem on the unit 4-cube
F4, G4 = rnd_q(5), rnd_q(6)
vol = G4 * DL(F4) + DR(G4) * F4 + ZQ


def integrate_cube(f, skip=None):
    for a in range(4):
        if a != skip:
            f = qmap(f, lambda c: sp.integrate(c, (X[a], 0, 1)))
    return f


lhs = integrate_cube(vol)
rhs = ZQ
for a in range(4):
    for val, sgn in ((1, 1), (0, -1)):
        face = qmap(G4 * (sgn * EPS[a]) * F4, lambda c: c.subs(X[a], val))
        rhs += integrate_cube(face, skip=a)
assert qeq(lhs, rhs)
print("7. int [(G D) F + G (D F)] dV = oint G n F dS on the unit 4-cube")

# 8. kernel on a sphere; integrals over S^3 of polynomials in the unit normal
# n = (n0, n1, n2, n3) use the exact monomial formula
#   int_{S^3} n0^a n1^b n2^c n3^d dOmega
#     = 2 prod Gamma((p+1)/2) / Gamma((a+b+c+d+4)/2)   (all exponents even)
r, ep = sp.symbols("r epsilon", positive=True)
N = sp.symbols("n0:4", real=True)
nq = Q(*N)


def int_S3(c):
    poly = sp.Poly(sp.expand(c), *N)
    tot = 0
    for exps, coef in poly.terms():
        if any(p % 2 for p in exps):
            continue
        tot += coef * 2 * sp.prod([sp.gamma(sp.Rational(p + 1, 2)) for p in exps]) / sp.gamma(sp.Rational(sum(exps) + 4, 2))
    return sp.simplify(tot)


assert sp.simplify(int_S3(1) - 2 * sp.pi**2) == 0
assert sp.simplify(int_S3(N[0] ** 2 + N[1] ** 2 + N[2] ** 2 + N[3] ** 2) - 2 * sp.pi**2) == 0
sph = lambda f, rad, c0=(0, 0, 0, 0): qmap(f, lambda c: c.subs({X[a]: c0[a] + rad * N[a] for a in range(4)}, simultaneous=True))
nbar = Q(N[0], -N[1], -N[2], -N[3])
nn1 = N[0] ** 2 + N[1] ** 2 + N[2] ** 2 + N[3] ** 2
# E(eps n) = eps nbar / eps^4 when |n| = 1
En = qmap(sph(Ek, ep), lambda c: sp.simplify(c.subs(N[0] ** 2, 1 - N[1] ** 2 - N[2] ** 2 - N[3] ** 2)))
assert qeq(En, qmap(nbar, lambda c: c / ep**3))
assert qeq(nbar * nq, scal(nn1)) and qeq(nq * nbar, scal(nn1))
print("8. E(eps n) = nbar/eps^3, so E n = n E = 1/eps^3 on the sphere; |S^3| = 2 pi^2")

# 9. Cauchy-Pompeiu on the unit ball at a = 0 (E(r n) = nbar / r^3, dV = r^3 dr dOmega)
F9 = rnd_q(7)
F0 = qmap(F9, lambda c: c.subs({v: 0 for v in X}))


def qint_S3(f):
    return qmap(f, int_S3)


bd = qint_S3(nbar * nq * sph(F9, 1))
vl = qint_S3(qmap(nbar * sph(DL(F9), r), lambda c: sp.integrate(sp.expand(c), (r, 0, 1))))
assert qeq(qmap(bd - vl, lambda c: c / (2 * sp.pi**2)), F0)
bd = qint_S3(sph(F9, 1) * nq * nbar)
vl = qint_S3(qmap(sph(DR(F9), r) * nbar, lambda c: sp.integrate(sp.expand(c), (r, 0, 1))))
assert qeq(qmap(bd - vl, lambda c: c / (2 * sp.pi**2)), F0)
print("9. Cauchy-Pompeiu (left and right) on the unit ball at a = 0")

# 10. Cauchy-Fueter at an off-center point (numeric) and mean value (exact)
H = lap(q * q * q) + (z1 * z2 + z2 * z1) * EPS[2] + Q(1, 2, -1, 3)
assert qzero(DL(H))
fH = sp.lambdify(X, list(comps(H)), "numpy")


def qmul_np(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
            a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)


a_pt = np.array([0.2, -0.1, 0.3, 0.1])
NQ = 48
gc, wc = np.polynomial.legendre.leggauss(NQ)
CH = (gc + 1) * np.pi / 2
TH = CH
PH = np.arange(2 * NQ) * np.pi / NQ
C, T, P = np.meshgrid(CH, TH, PH, indexing="ij")
W = np.einsum("i,j->ij", wc, wc)[:, :, None] * (np.pi / 2) ** 2 * (np.pi / NQ)
nn = (np.cos(C), np.sin(C) * np.cos(T), np.sin(C) * np.sin(T) * np.cos(P), np.sin(C) * np.sin(T) * np.sin(P))
dq = tuple(nn[a] - a_pt[a] for a in range(4))
rho2 = sum(c**2 for c in dq)
Enp = (dq[0] / rho2**2, -dq[1] / rho2**2, -dq[2] / rho2**2, -dq[3] / rho2**2)
Fv = tuple(np.broadcast_to(c, C.shape) for c in fH(*nn))
integrand = qmul_np(qmul_np(Enp, nn), Fv)
dO = np.sin(C) ** 2 * np.sin(T)
val = np.array([np.sum(c * dO * W) for c in integrand]) / (2 * np.pi**2)
assert np.allclose(val, fH(*a_pt), atol=1e-11), (val, fH(*a_pt))
# mean value around an arbitrary center c with radius R
cc = sp.symbols("c0:4", real=True)
R = sp.Symbol("R", positive=True)
mean = qmap(qint_S3(sph(H, R, cc)), lambda c: sp.simplify(c / (2 * sp.pi**2)))
assert qeq(mean, qmap(H, lambda c: c.subs(dict(zip(X, cc)))))
print("10. Cauchy-Fueter formula at a = (0.2,-0.1,0.3,0.1) (numeric); mean value on spheres (exact)")

print("all checks passed")
