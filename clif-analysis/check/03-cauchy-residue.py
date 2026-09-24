"""Checks for article 03 (Cauchy's integral theorem, integral formula, residues).

Cl_{2,0}(R) is realized as M_2(R) as in 02. Polar coordinates about a point
a: x - a = r n with n = cos(th) e1 + sin(th) e2, and t = n I is the unit
tangent (counterclockwise). z = e1 x is the even element x + I y.

1. Polar form of D: D = n (d_r + (I/r) d_th).
2. Annulus identity: d_r (r n F) + d_th (n I F) = r D F, so for monogenic F
   oint_{|x-a|=r} n F ds does not depend on r.
3. Two-sided annulus identity:
   d_r (r G n F) + d_th (G t F) = r [(G D) F + G (D F)],
   with G D = (d_r G) n + (1/r)(d_th G) t. Mean value: int H dth = 2 pi H(0)
   for monogenic H.
4. Kernel on a circle: (x - a)^{-1} n = 1/r, and
   (x - a)^{-1} n ds = I^{-1} dz / (z - a) with dz = e1 dx (a even-ified).
5. Cauchy-Pompeiu for a generic (non-monogenic) F:
   F(a) = (1/2pi) [ oint (x-a)^{-1} n F ds - int (x-a)^{-1} D F dA ]
   on the unit disk (a = 0, exact) and the unit square (a off-center, numeric).
   Also the right version F(a) = (1/2pi)[oint F n (x-a)^{-1} ds - int (F D)(x-a)^{-1} dA].
6. Residues, Res_a F = (1/2pi) oint_{small circle} n F ds:
   Res of x^{-1} c is c, Res of d_x x^{-1}, d_y x^{-1} (times c) is 0.
7. Example F = 1/(1 + z^2) (even): poles at x = +-e2, GA residues
   e1 Res_classical, and the residue theorem on a big circle.
"""

import mpmath as mp
import sympy as sp

x, y, r, th, eps = sp.symbols("x y r theta epsilon", real=True)

one = sp.eye(2)
e1 = sp.Matrix([[1, 0], [0, -1]])
e2 = sp.Matrix([[0, 1], [1, 0]])
I = e1 * e2
assert I**2 == -one
Z = sp.zeros(2, 2)


def simp(M):
    return M.applyfunc(lambda e: sp.simplify(sp.expand(e)))


def eq(A, B):
    return simp(A - B) == Z


def D(F):
    return e1 * F.diff(x) + e2 * F.diff(y)


def Dr(F):
    return F.diff(x) * e1 + F.diff(y) * e2


def poly(seed):
    import random

    rnd = random.Random(seed)
    mons = [x**i * y**j for i in range(4) for j in range(4 - i)]
    return sp.Matrix(2, 2, lambda i, j: sum(rnd.randint(-3, 3) * m for m in mons))


def vec(px, py):
    return px * e1 + py * e2


def inv_vec(px, py):
    return vec(px, py) / (px**2 + py**2)


F = poly(1)
G = poly(2)
nvec = vec(sp.cos(th), sp.sin(th))
tvec = vec(-sp.sin(th), sp.cos(th))
assert eq(nvec * I, tvec)
to_polar = {x: r * sp.cos(th), y: r * sp.sin(th)}


def polar(M):
    return M.subs(to_polar, simultaneous=True)


# 1. polar form of D
Fp = polar(F)
assert eq(polar(D(F)), nvec * (Fp.diff(r) + I * Fp.diff(th) / r))
print("1. D = n (d_r + (I/r) d_th)")

# 2. annulus identity
lhs = (r * nvec * Fp).diff(r) + (nvec * I * Fp).diff(th)
assert eq(lhs, r * polar(D(F)))
print("2. d_r (r n F) + d_th (n I F) = r D F")

# 3. two-sided annulus identity
Gp = polar(G)
assert eq(polar(Dr(G)), Gp.diff(r) * nvec + Gp.diff(th) * tvec / r)
lhs = (r * Gp * nvec * Fp).diff(r) + (Gp * tvec * Fp).diff(th)
assert eq(lhs, r * polar(Dr(G) * F + G * D(F)))
zz = x * one + y * I
H = simp(zz**3 - 2 * zz + 5 * one)
assert D(H) == Z
avg = simp(polar(H).applyfunc(lambda e: sp.integrate(sp.expand(e), (th, 0, 2 * sp.pi))))
assert eq(avg, 2 * sp.pi * H.subs({x: 0, y: 0}))
print("3. d_r (r G n F) + d_th (G t F) = r [(G D) F + G (D F)]; int H dth = 2 pi H(0)")

# 4. kernel on a circle
a1, a2 = sp.symbols("a1 a2", real=True)
px, py = a1 + eps * sp.cos(th), a2 + eps * sp.sin(th)
K = inv_vec(px - a1, py - a2)
assert eq(K * nvec, one / eps)
assert eq(nvec * K, one / eps)
# (x-a)^{-1} n ds = I^{-1} dz/(z-a): ds = eps dth, dx = t eps dth, dz = e1 dx
zma = e1 * vec(px - a1, py - a2)  # z - a, even
dz = e1 * tvec * eps
assert eq(K * nvec * eps, I.inv() * dz * zma.inv())
print("4. (x-a)^{-1} n = 1/eps on the circle; (x-a)^{-1} n ds = I^{-1} dz/(z-a)")

# 5. Cauchy-Pompeiu
# unit disk, a = 0, exact
Kp = inv_vec(r * sp.cos(th), r * sp.sin(th))
bd = simp((Kp * nvec * Fp).subs(r, 1))
bd = bd.applyfunc(lambda e: sp.integrate(sp.expand(e), (th, 0, 2 * sp.pi)))
ar = simp(Kp * polar(D(F)) * r)
ar = ar.applyfunc(lambda e: sp.integrate(sp.integrate(sp.expand(e), (r, 0, 1)), (th, 0, 2 * sp.pi)))
assert eq((bd - ar) / (2 * sp.pi), F.subs({x: 0, y: 0}))
# right version
bd = simp((Fp * nvec * Kp).subs(r, 1))
bd = bd.applyfunc(lambda e: sp.integrate(sp.expand(e), (th, 0, 2 * sp.pi)))
ar = simp(polar(Dr(F)) * Kp * r)
ar = ar.applyfunc(lambda e: sp.integrate(sp.integrate(sp.expand(e), (r, 0, 1)), (th, 0, 2 * sp.pi)))
assert eq((bd - ar) / (2 * sp.pi), F.subs({x: 0, y: 0}))
print("5a. Cauchy-Pompeiu (left and right) on the unit disk at a = 0")

# unit square, a = (1/3, 1/5), numeric
mp.mp.dps = 20
ax, ay = mp.mpf(1) / 3, mp.mpf(1) / 5
fF = sp.lambdify((x, y), F, "mpmath")
fDF = sp.lambdify((x, y), D(F), "mpmath")
E1, E2 = mp.matrix([[1, 0], [0, -1]]), mp.matrix([[0, 1], [1, 0]])


def kern(X, Y):
    dx, dy = X - ax, Y - ay
    return (dx * E1 + dy * E2) / (dx**2 + dy**2)


def mquad(f, *lims):
    return mp.matrix([[mp.quad(lambda *s: f(*s)[i, j], *lims) for j in range(2)] for i in range(2)])


edges = [  # (point(t), outward normal)
    (lambda s: (s, 0), -E2),
    (lambda s: (1, s), E1),
    (lambda s: (s, 1), E2),
    (lambda s: (0, s), -E1),
]
bdn = mp.matrix(2, 2)
for pt, nn in edges:
    bdn += mquad(lambda s: kern(*pt(s)) * nn * mp.matrix(fF(*pt(s))), [0, 1])
arn = mquad(lambda X, Y: kern(X, Y) * mp.matrix(fDF(X, Y)), [0, ax, 1], [0, ay, 1])
val = (bdn - arn) / (2 * mp.pi)
assert mp.mnorm(val - mp.matrix(fF(ax, ay)), 1) < mp.mpf(10) ** -12
print("5b. Cauchy-Pompeiu on the unit square at a = (1/3, 1/5) (numeric)")

# 6. residues of the kernel and its derivatives
c = sp.Matrix(2, 2, sp.symbols("c0:4", real=True))
Kxy = inv_vec(x, y)
assert eq(D(Kxy), Z)
circ = {x: eps * sp.cos(th), y: eps * sp.sin(th)}


def residue(Fx):
    integrand = simp(nvec * Fx.subs(circ, simultaneous=True) * eps)
    return simp(integrand.applyfunc(lambda e: sp.integrate(e, (th, 0, 2 * sp.pi))) / (2 * sp.pi))


assert eq(residue(Kxy * c), c)
for dk in (Kxy.diff(x), Kxy.diff(y), Kxy.diff(x, 2)):
    assert eq(D(dk), Z)
    assert eq(residue(dk * c), Z)
# 1/z = x^{-1} e1: residue e1
assert eq(e1 * vec(x, y) * inv_vec(x, y) * e1, one)
assert eq(residue(Kxy * e1), e1)
print("6. Res x^{-1} c = c, Res (d x^{-1}) c = 0, Res 1/z = e1")

# 7. F = 1/(1 + z^2)
Fz = simp((one + zz**2).inv())
assert eq(D(Fz), Z)
fz = sp.lambdify((x, y), Fz, "mpmath")
En = lambda t_: mp.cos(t_) * E1 + mp.sin(t_) * E2


def circle_res(cx, cy, rad):
    return mquad(lambda t_: En(t_) * mp.matrix(fz(cx + rad * mp.cos(t_), cy + rad * mp.sin(t_))) * rad, [0, mp.pi, 2 * mp.pi]) / (2 * mp.pi)


Ic = mp.matrix([[0, 1], [-1, 0]])
for sgn in (1, -1):
    res = circle_res(0, sgn, mp.mpf(1) / 2)
    # classical Res at z = +-i of 1/(1+z^2) is -+ i/2; GA residue = e1 * classical
    expect = E1 * (-sgn * Ic / 2)
    assert mp.mnorm(res - expect, 1) < mp.mpf(10) ** -12, sgn
total = circle_res(0, 0, 3)
assert mp.mnorm(total - (E1 * (-Ic / 2) + E1 * (Ic / 2)), 1) < mp.mpf(10) ** -12
print("7. 1/(1+z^2): Res at x = +-e2 is e1 (-+I/2), sum = big-circle integral = 0")

print("all checks passed")
