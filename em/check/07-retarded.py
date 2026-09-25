"""Checks for em/07-retarded.md (hyperbolic fundamental solution and the sorting of properties).

1. Retarded solution: for r > 0, (d_0^2 - Lap)(f(x_0 - r)/r) = 0 for any f; the
   flux of -grad G, G = f(x_0 - r)/(4 pi r), through the sphere |x| = eps tends to
   f(x_0), which is the normalization (d_0^2 - Lap) G = delta^3(x) f(x_0).
2. D of the kernel (Cl_{1,3}): D(f(x_0 - r)/(4 pi r)) = f'/(4 pi r) (gamma_0 + n)
   + f/(4 pi r^2) n with n = sum (x_k/r) gamma_k; (gamma_0 + n)^2 = 0.
3. Retarded potential and the Lorenz gauge: for R = |x - y| and J(s, y),
   div_x(J(x_0 - R, y)/R) = -div_y(J(x_0 - R, y)/R) + (div J)(x_0 - R, y)/R and
   d_0(rho(x_0 - R, y)/R) = (d_s rho)(x_0 - R, y)/R, so continuity gives
   d_0 phi + c div A = 0 after integrating the total y-divergence away.
4. Space dimension 2: G_2 = 1/(2 pi sqrt(x_0^2 - r^2)) satisfies d_0^2 - Lap_2 = 0 inside
   the cone, so the support is not concentrated on the cone.
5. Kirchhoff/Poisson formula for polynomial data: with the spherical mean
   M_R[p] = sum_k R^{2k} Lap^k p / (2k+1)! (n = 3),
   u = d_t(t M_{ct}[g]) + t M_{ct}[h] solves u_tt = c^2 Lap u, u(0) = g, u_t(0) = h.
   The spherical mean formula is checked against direct integration for monomials.
6. Lost properties for the plane wave u = f(x_0 - x_3): the spherical mean of
   (x_0 - x_3)^2 over |x| = R is x_0^2 + R^2/3 (no mean value property); the C^infty
   profile f(s) = exp(-1/s^2) (s > 0), 0 (s <= 0) is a solution vanishing on a
   half-space without being zero (no unique continuation); cos(x_0 - x_3) is
   bounded and nonconstant (no Liouville).
7. Radiation: (f'/r)^2 4 pi r^2 is independent of r, (f/r^2)^2 4 pi r^2 -> 0.
8. Radiation: D(x_0 - r) = gamma_0 + n = k (null); for vectors (k ^ a)^2 = (k.a)^2 - k^2 a^2, so
   k.a = 0 gives (k ^ a)^2 = 0 (the 1/r part of the retarded field is null).
"""

import sympy as sp

from common.clifford import Alg, eq, mv

x0 = sp.Symbol("x0", real=True)
Xs = sp.symbols("x1:4", real=True)
r = sp.sqrt(sum(v**2 for v in Xs))
f = sp.Function("f")
s = sp.Symbol("s", real=True)


def box(u, xs=Xs):
    return sp.diff(u, x0, 2) - sum(sp.diff(u, v, 2) for v in xs)


# ---------------------------------------------------------------- 1.
G = f(x0 - r) / r
assert sp.simplify(box(G)) == 0
assert sp.simplify(box(f(x0 + r) / r)) == 0
eps = sp.Symbol("epsilon", positive=True)
R_ = sp.Symbol("R", positive=True)
Gr = f(x0 - R_) / (4 * sp.pi * R_)
flux = -sp.diff(Gr, R_) * 4 * sp.pi * R_**2
fpR = sp.Subs(sp.Derivative(f(s), s), s, x0 - R_).doit()
assert sp.simplify(flux - (f(x0 - R_) + R_ * fpR)) == 0  # -> f(x0) as R -> 0
fc = sp.exp(-s**2) * sp.cos(3 * s)  # a concrete profile
fluxc = flux.subs(f(x0 - R_), fc.subs(s, x0 - R_)).subs(fpR, sp.diff(fc, s).subs(s, x0 - R_))
assert sp.simplify(sp.limit(fluxc.subs(R_, eps), eps, 0) - fc.subs(s, x0)) == 0
print("1. box(f(x0 - r)/r) = 0 (r > 0); flux of -grad G through small spheres -> f(x0)")

# ---------------------------------------------------------------- 2.
S = Alg(1, 3)
g = S.e
Y = S.X
rY = sp.sqrt(Y[1]**2 + Y[2]**2 + Y[3]**2)
GY = f(Y[0] - rY) / (4 * sp.pi * rY)
DG = S.D(mv(GY, S.neg))
n = sum((Y[k] / rY * g[k] for k in (1, 2, 3)), S.zero())
fp = sp.Subs(sp.Derivative(f(s), s), s, Y[0] - rY).doit()
expect = fp / (4 * sp.pi * rY) * (g[0] + n) + f(Y[0] - rY) / (4 * sp.pi * rY**2) * n
assert eq(DG, expect)
assert eq((g[0] + n) * (g[0] + n), 0)
print("2. D G = f'/(4 pi r)(gamma_0 + n) + f/(4 pi r^2) n, (gamma_0 + n)^2 = 0")

# ---------------------------------------------------------------- 3.
ys = sp.symbols("y1:4", real=True)
Rxy = sp.sqrt(sum((Xs[k] - ys[k])**2 for k in range(3)))
# concrete test functions J_k(s, y), rho(s, y) (generic Function objects leave Subs
# expressions that sympy does not identify)
t_ = sp.Symbol("t_")
Jc = [t_ * ys[1]**2 + sp.sin(t_) * ys[0] * ys[2],
      sp.exp(-t_**2) * ys[2] + t_**3 * ys[0]**2,
      sp.cos(t_ * ys[1]) + ys[0] * ys[1] * ys[2]]
rhoc = sp.exp(t_) * ys[0] * ys[1] + t_ * ys[2]**2
ret = x0 - Rxy
Jret = [v.subs(t_, ret) for v in Jc]
divx = sum(sp.diff(Jret[k] / Rxy, Xs[k]) for k in range(3))
divy_total = sum(sp.diff(Jret[k] / Rxy, ys[k]) for k in range(3))
divJ_partial = sum(sp.diff(Jc[k], ys[k]) for k in range(3)).subs(t_, ret)
assert sp.simplify(divx - (-divy_total + divJ_partial / Rxy)) == 0
drho = sp.diff(rhoc.subs(t_, ret) / Rxy, x0)
assert sp.simplify(drho - sp.diff(rhoc, t_).subs(t_, ret) / Rxy) == 0
print("3. div_x(J(x0 - R, y)/R) = -div_y(...) + (div J)_ret/R, d0(rho_ret/R) = (d_s rho)_ret/R")

# ---------------------------------------------------------------- 4.
x1, x2 = Xs[0], Xs[1]
G2 = 1 / (2 * sp.pi * sp.sqrt(x0**2 - x1**2 - x2**2))
assert sp.simplify(box(G2, (x1, x2))) == 0
print("4. 2 space dimensions: 1/(2 pi sqrt(x0^2 - r^2)) solves the wave equation inside the cone")

# ---------------------------------------------------------------- 5.
lap = lambda p: sum(sp.diff(p, v, 2) for v in Xs)


def smean(p, R):
    tot, term, k = 0, p, 0
    while term != 0:
        tot += R**(2 * k) * term / sp.factorial(2 * k + 1)
        term = sp.expand(lap(term))
        k += 1
    return sp.expand(tot)


th, ph = sp.symbols("theta phi", real=True)
Rs = sp.Symbol("R", positive=True)
sph = [Rs * sp.sin(th) * sp.cos(ph), Rs * sp.sin(th) * sp.sin(ph), Rs * sp.cos(th)]
for mono in (Xs[0]**2, Xs[0]**2 * Xs[1]**2, Xs[2]**4, Xs[0] * Xs[1] * Xs[2]**2):
    shifted = mono.subs({Xs[k]: Xs[k] + sph[k] for k in range(3)}, simultaneous=True)
    avg = sp.integrate(sp.integrate(sp.expand(shifted) * sp.sin(th), (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)) / (4 * sp.pi)
    assert sp.expand(avg - smean(mono, Rs)) == 0
t, c = sp.Symbol("t", positive=True), sp.Symbol("c", positive=True)
gdat = Xs[0]**4 + 3 * Xs[0] * Xs[1]**2 * Xs[2] - Xs[2]**3
hdat = Xs[1]**2 * Xs[2]**2 + Xs[0]
u = sp.diff(t * smean(gdat, c * t), t) + t * smean(hdat, c * t)
assert sp.expand(sp.diff(u, t, 2) - c**2 * lap(u)) == 0
assert sp.expand(u.subs(t, 0) - gdat) == 0
assert sp.expand(sp.diff(u, t).subs(t, 0) - hdat) == 0
print("5. spherical mean sum R^{2k} Lap^k p/(2k+1)!; Kirchhoff u = d_t(t M[g]) + t M[h] solves the wave equation")

# ---------------------------------------------------------------- 6.
pw = (x0 - Xs[2])**2
shifted = pw.subs({Xs[k]: sph[k] for k in range(3)}, simultaneous=True)
avg = sp.integrate(sp.integrate(sp.expand(shifted) * sp.sin(th), (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)) / (4 * sp.pi)
assert sp.expand(avg - (x0**2 + Rs**2 / 3)) == 0
assert sp.simplify(box((x0 - Xs[2])**2)) == 0
bump = sp.exp(-1 / s**2)
for k in range(4):
    assert sp.limit(sp.diff(bump, s, k), s, 0, "+") == 0  # all derivatives vanish at 0
w = sp.Function("w")
assert sp.simplify(box(w(x0 - Xs[2]))) == 0
assert sp.simplify(box(sp.cos(x0 - Xs[2]))) == 0
print("6. mean of (x0 - x3)^2 over |x| = R is x0^2 + R^2/3; f(x0 - x3) solves for any f; "
      "exp(-1/s^2) is flat at 0 (solution vanishing on a half-space); cos(x0 - x3) bounded")

# ---------------------------------------------------------------- 7.
# Radiation: through the sphere of radius r, a field ~ f'/r gives energy flow ~ |f'|^2/r^2
# times the area 4 pi r^2, independent of r; a field ~ f/r^2 gives ~ 1/r^2 -> 0.
rr_, fp_, f_ = sp.symbols("r fp f", positive=True)
assert sp.simplify((fp_ / rr_)**2 * 4 * sp.pi * rr_**2 - 4 * sp.pi * fp_**2) == 0
assert sp.limit((f_ / rr_**2)**2 * 4 * sp.pi * rr_**2, rr_, sp.oo) == 0
print("7. 1/r part: flux through the sphere independent of r; 1/r^2 part: flux -> 0")

# ---------------------------------------------------------------- 8.
kk = S.D(mv(Y[0] - rY, S.neg))
assert eq(kk, g[0] + n)
Av = sp.symbols("A0:4", real=True)
av = sum((Av[m_] * g[m_] for m_ in range(4)), S.zero())
wedge = (kk * av - av * kk) / 2
dotk = ((kk * av + av * kk) / 2)
assert eq(wedge * wedge, dotk * dotk - (kk * kk) * (av * av))
print("8. D(x0 - r) = gamma_0 + n null; (k ^ a)^2 = (k.a)^2 - k^2 a^2, so the radiation part is null")
