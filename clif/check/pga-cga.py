"""Checks for pga-cga.md.

Verifies the article's claims numerically with a small Clifford algebra
engine (basis blades as bitmasks, signature entries in {+1, -1, 0}).

PGA, Cl(3,0,1) with e0^2 = 0:

1. Plane embedding: p = a e1 + b e2 + c e3 + d e0 for ax+by+cz+d = 0
   satisfies p^2 = a^2+b^2+c^2 (= 1 for a unit normal).
2. Point as the meet of three planes:
   P(x,y,z) = e123 - x e023 + y e013 - z e012, with P(0,0,0) = e123.
   The e123 coefficient is the homogeneous weight (scaling is harmless).
3. Reflection: X -> p X p sends points, planes and lines to their
   mirror images (planes and lines pick up an overall sign, which is
   irrelevant for homogeneous coordinates).
4. Parallel planes q = n + d1 e0, r = n + d2 e0 give the motor
   M = rq = 1 + (d1-d2) n e0 with (n e0)^2 = 0, acting as the
   translation by 2(d1-d2) n.
5. Intersecting planes through the origin at angle th give
   M = rq = cos th - sin th e1e2, the rotation by 2th about their
   intersection line.
6. Motors: the even subalgebra is 8-dimensional (1, six bivectors,
   pseudoscalar) and closed, motors satisfy M ~M = 1, and a screw
   motion (rotation about a line off the origin, then a slide along
   it) is realized as a single motor.
7. Degeneracy: the pseudoscalar e0123 squares to 0 (not invertible).

CGA, Cl(4,1) with null basis no = (e- - e+)/2, ni = e- + e+:

8. no^2 = ni^2 = 0, no . ni = -1.
9. Point embedding P(x) = no + x + |x|^2/2 ni is null, and
   P(x) . P(y) = -|x-y|^2/2.
10. Sphere S = P(c) - r^2/2 ni and plane pi = n + d ni satisfy
    S^2 = r^2, pi^2 = 1, and P . S = 0 (resp. P . pi = 0) iff the point
    lies on it; the plane is the limit of spheres with r -> infinity.
11. Reflection X -> S X S in the unit sphere at the origin is the
    inversion x -> x/|x|^2; reflection in pi is the ordinary mirror.
12. Translator T = exp(-t ni/2) = 1 - t ni/2 (the square of t ni is 0)
    gives P(x) -> P(x+t); the dilator exp(lam/2 no^ni) scales by
    e^lam; both satisfy V ~V = 1.
13. Bivector counts: 6 = dim SE(3) in PGA, 10 = dim O(4,1) in CGA.
14. PGA sits inside CGA: e1, e2, e3, ni generate a subalgebra
    isomorphic to Cl(3,0,1) (16 linearly independent blades), and the
    CGA plane n + d ni matches the PGA plane n - d e0.
"""

import itertools
import math
import numpy as np

class GA:
    """Clifford algebra over R with signature sq[i] in {+1, -1, 0}."""

    def __init__(self, sq, names=None):
        self.sq = list(sq)
        self.n = len(self.sq)
        self.dim = 1 << self.n
        self.names = names or [f"e{i}" for i in range(self.n)]

    def blade_mul(self, a, b):  # basis blade product -> (sign, blade)
        sign = 1
        for i in range(self.n):
            if not (b >> i) & 1:
                continue
            if bin(a >> (i + 1)).count("1") % 2:  # swaps past higher bits of a
                sign = -sign
            if (a >> i) & 1:
                if self.sq[i] == 0:
                    return 0, 0
                sign *= self.sq[i]
                a &= ~(1 << i)
            else:
                a |= 1 << i
        return sign, a

    def zero(self):
        return np.zeros(self.dim)

    def scalar(self, c=1.0):
        z = self.zero()
        z[0] = c
        return z

    def e(self, *idx):
        z, blade, sign = self.zero(), 0, 1
        for i in idx:
            s, blade = self.blade_mul(blade, 1 << i)
            sign *= s
        z[blade] = sign
        return z

    def mul(self, x, y):
        z = self.zero()
        for a in np.nonzero(x)[0]:
            for b in np.nonzero(y)[0]:
                s, c = self.blade_mul(int(a), int(b))
                if s:
                    z[c] += s * x[a] * y[b]
        return z

    def wedge(self, x, y):
        z = self.zero()
        for a in np.nonzero(x)[0]:
            for b in np.nonzero(y)[0]:
                if int(a) & int(b):
                    continue
                s, c = self.blade_mul(int(a), int(b))
                if s:
                    z[c] += s * x[a] * y[b]
        return z

    def dot(self, x, y):  # symmetric part (scalar product)
        return (self.mul(x, y)[0] + self.mul(y, x)[0]) / 2

    def rev(self, x):  # reverse
        z = self.zero()
        for a in range(self.dim):
            k = bin(a).count("1")
            z[a] = x[a] * (-1) ** (k * (k - 1) // 2)
        return z

    def grade(self, x, k):
        z = self.zero()
        for a in range(self.dim):
            if bin(a).count("1") == k:
                z[a] = x[a]
        return z

    def sandwich(self, v, x):
        return self.mul(self.mul(v, x), self.rev(v))

    def exp(self, x, terms=60):
        z, t = self.scalar(1.0), self.scalar(1.0)
        for k in range(1, terms):
            t = self.mul(t, x) / k
            z = z + t
        return z

    def str(self, x, tol=1e-10):
        parts = []
        for a in range(self.dim):
            if abs(x[a]) < tol:
                continue
            if a == 0:
                parts.append(f"{x[a]:+.4g}")
            else:
                nm = "".join(self.names[i] for i in range(self.n) if (a >> i) & 1)
                parts.append(f"{x[a]:+.4g}{nm}")
        return " ".join(parts) or "0"


def check(label, cond):
    print(f"  [{'OK' if cond else 'NG'}] {label}")
    assert cond, label

def close(a, b, tol=1e-9):
    return np.allclose(a, b, atol=tol)

def proportional(x, y, tol=1e-9):  # equal as homogeneous coordinates
    i = int(np.argmax(np.abs(y)))
    return abs(y[i]) > tol and close(x * y[i] / x[i], y, tol)


# ============================================================
# PGA: Cl(3,0,1)
# ============================================================
pga = GA([0, 1, 1, 1], ["e0", "e1", "e2", "e3"])
E0, E1, E2, E3 = (pga.e(i) for i in range(4))

def plane(n, d):  # n.x + d = 0
    return n[0] * E1 + n[1] * E2 + n[2] * E3 + d * E0

def point(x):  # meet of the three planes x_a = const
    p = plane((1, 0, 0), -x[0])
    for a, i in ((1, 1), (2, 2)):
        p = pga.mul(p, plane(tuple(int(j == a) for j in range(3)), -x[i]))
    return p

print("=== 1-2. PGA: planes and points ===")
n = np.array([1.0, 2.0, -2.0]) / 3  # unit normal
check("p^2 = |n|^2 = 1", close(pga.mul(plane(n, 7), plane(n, 7)), pga.scalar(1)))
check("p^2 = a^2+b^2+c^2", close(pga.mul(plane((1, 2, -2), 7), plane((1, 2, -2), 7))[0], 9))

x = np.array([2.0, 3.0, 5.0])
expect = (pga.e(1, 2, 3) - x[0] * pga.e(0, 2, 3)
          + x[1] * pga.e(0, 1, 3) - x[2] * pga.e(0, 1, 2))
check("P(x) = e123 - x e023 + y e013 - z e012", close(point(x), expect))
check("P(0) = e123", close(point((0, 0, 0)), pga.e(1, 2, 3)))
print("   P(2,3,5) =", pga.str(point(x)))

print("=== 3. reflection X -> p X p ===")
p = plane((0, 0, 1), -1.0)  # x3 = 1
check("point mirrored", close(pga.mul(pga.mul(p, point(x)), p), point((2, 3, -3))))
q = plane((0, 1, 0), -4.0)  # x2 = 4
check("plane mirrored (up to sign)",
      proportional(pga.mul(pga.mul(p, q), p), q))          # parallel to the mirror
r = plane((0, 0, 1), -5.0)  # x3 = 5 -> x3 = -3
check("plane mirrored", proportional(pga.mul(pga.mul(p, r), p), plane((0, 0, 1), 3.0)))
line = pga.wedge(plane((1, 0, 0), 0.0), plane((0, 1, 0), 0.0))  # the x3 axis
check("line = e1e2", close(line, pga.e(1, 2)))
check("axis line fixed by the mirror x3=1", proportional(pga.mul(pga.mul(p, line), p), line))

print("=== 4. parallel planes -> translation ===")
d1, d2 = -1.0, -3.0
nz = (0, 0, 1)
q, r = plane(nz, d1), plane(nz, d2)
M = pga.mul(r, q)
check("M = 1 + (d1-d2) n e0", close(M, pga.scalar(1) + (d1 - d2) * pga.mul(E3, E0)))
check("(n e0)^2 = 0", close(pga.mul(pga.mul(E3, E0), pga.mul(E3, E0)), pga.zero()))
check("M ~M = 1", close(pga.mul(M, pga.rev(M)), pga.scalar(1)))
check("translation by 2(d1-d2) n", close(pga.sandwich(M, point(x)), point(x + 2 * (d1 - d2) * np.array(nz))))
print("   M =", pga.str(M), "-> P(2,3,5) becomes P(2,3,9)")

print("=== 5. intersecting planes -> rotation ===")
th = math.radians(37)
q2 = plane((1, 0, 0), 0.0)
r2 = plane((math.cos(th), math.sin(th), 0), 0.0)
M2 = pga.mul(r2, q2)
check("M = cos th - sin th e1e2", close(M2, math.cos(th) * pga.scalar(1) - math.sin(th) * pga.e(1, 2)))
c, s = math.cos(2 * th), math.sin(2 * th)
check("rotation by 2th about the x3 axis",
      close(pga.sandwich(M2, point((1, 0, 7))), point((c, s, 7))))

print("=== 6. motors ===")
even = [b for b in range(pga.dim) if bin(b).count("1") % 2 == 0]
check("even subalgebra is 8-dimensional", len(even) == 8)
prods = []
for a, b in itertools.product(even, repeat=2):
    z = pga.mul(pga.e(*[i for i in range(4) if (a >> i) & 1]),
                pga.e(*[i for i in range(4) if (b >> i) & 1]))
    prods.append(z)
check("even subalgebra is closed", all(close(z, sum(z[k] * pga.e(*[i for i in range(4) if (k >> i) & 1])
                                                   for k in even)) for z in prods))
# screw motion about the vertical line through x0, by angle th and slide dist
x0 = np.array([1.0, 0.0, 0.0])
th, dist = math.radians(90), 2.0
axis = pga.wedge(plane((1, 0, 0), -x0[0]), plane((0, 1, 0), -x0[1]))
check("axis^2 = -1", close(pga.mul(axis, axis), -pga.scalar(1)))
rotor = pga.exp(-th / 2 * axis)
slide = pga.scalar(1) + dist / 2 * pga.mul(E3, E0)
Mscrew = pga.mul(slide, rotor)
check("screw motor is a motor", close(pga.mul(Mscrew, pga.rev(Mscrew)), pga.scalar(1)))
check("rotation about the shifted axis",
      close(pga.sandwich(rotor, point((2.0, 0.0, 0.0))), point((1.0, 1.0, 0.0))))
check("screw = rotation + slide along the axis",
      close(pga.sandwich(Mscrew, point((2.0, 0.0, 0.0))), point((1.0, 1.0, dist))))

print("=== 7. degenerate pseudoscalar ===")
I4 = pga.e(0, 1, 2, 3)
check("e0123^2 = 0", close(pga.mul(I4, I4), pga.zero()))
check("bivectors: 6 = dim SE(3)", len([b for b in range(pga.dim) if bin(b).count("1") == 2]) == 6)


# ============================================================
# CGA: Cl(4,1)
# ============================================================
cga = GA([1, 1, 1, 1, -1], ["e1", "e2", "e3", "e+", "e-"])
F1, F2, F3, Fp, Fm = (cga.e(i) for i in range(5))
no = (Fm - Fp) / 2
ni = Fm + Fp

def vec(v):
    return v[0] * F1 + v[1] * F2 + v[2] * F3

def cpoint(v):
    v = np.asarray(v, dtype=float)
    return no + vec(v) + 0.5 * float(v @ v) * ni

def unpack(X):  # homogeneous null vector -> (weight, x)
    w = -cga.dot(X, ni)
    return w, np.array([X[1 << i] for i in range(3)]) / w

print("=== 8-9. CGA: null basis and points ===")
check("no^2 = 0", close(cga.mul(no, no), cga.zero()))
check("ni^2 = 0", close(cga.mul(ni, ni), cga.zero()))
check("no.ni = -1", close(cga.dot(no, ni), -1.0))
u, v = np.array([1.0, 2.0, -1.0]), np.array([0.0, 3.0, 4.0])
check("P(0) = no", close(cpoint((0, 0, 0)), no))
check("P(x)^2 = 0", close(cga.mul(cpoint(u), cpoint(u)), cga.zero()))
check("P(x).P(y) = -|x-y|^2/2", close(cga.dot(cpoint(u), cpoint(v)), -0.5 * (u - v) @ (u - v)))

print("=== 10. spheres and planes ===")
cen, rad = np.array([1.0, 0.0, 0.0]), 2.0
S = cpoint(cen) - 0.5 * rad ** 2 * ni
check("S^2 = r^2", close(cga.mul(S, S), rad ** 2 * cga.scalar(1)))
check("P.S = 0 on the sphere", close(cga.dot(cpoint(cen + np.array([0, rad, 0])), S), 0.0))
check("P.S != 0 off the sphere", not close(cga.dot(cpoint(u), S), 0.0))
nrm, dist0 = np.array([0.0, 0.0, 1.0]), 1.0
pi = vec(nrm) + dist0 * ni
check("pi^2 = 1", close(cga.mul(pi, pi), cga.scalar(1)))
check("P.pi = 0 on the plane", close(cga.dot(cpoint((5.0, -2.0, dist0)), pi), 0.0))
big = 1e6  # plane as the limit of spheres: center far away, radius equally large
Sbig = cpoint(nrm * (big + dist0)) - 0.5 * big ** 2 * ni
check("plane = limit of spheres", close(unpack(cga.mul(cga.mul(Sbig, cpoint(u)), Sbig))[1],
                                        unpack(cga.mul(cga.mul(pi, cpoint(u)), pi))[1], tol=1e-4))

print("=== 11. reflection = mirror and inversion ===")
S0 = no - 0.5 * ni  # unit sphere at the origin
check("inversion x -> x/|x|^2",
      close(unpack(cga.mul(cga.mul(S0, cpoint(u)), S0))[1], u / (u @ u)))
mirrored = u - 2 * (u @ nrm - dist0) * nrm
check("mirror in the plane", close(unpack(cga.mul(cga.mul(pi, cpoint(u)), pi))[1], mirrored))

print("=== 12. translator and dilator ===")
t = np.array([2.0, -1.0, 3.0])
T = cga.exp(-0.5 * cga.mul(vec(t), ni))
check("(t ni)^2 = 0", close(cga.mul(cga.mul(vec(t), ni), cga.mul(vec(t), ni)), cga.zero()))
check("T = 1 - t ni/2", close(T, cga.scalar(1) - 0.5 * cga.mul(vec(t), ni)))
check("T P(x) ~T = P(x+t)", close(cga.sandwich(T, cpoint(u)), cpoint(u + t)))
check("T ~T = 1", close(cga.mul(T, cga.rev(T)), cga.scalar(1)))
lam = math.log(3)
D = cga.exp(0.5 * lam * cga.wedge(no, ni))
check("dilator scales by e^lam", close(unpack(cga.sandwich(D, cpoint(u)))[1], math.exp(lam) * u))
check("D ~D = 1", close(cga.mul(D, cga.rev(D)), cga.scalar(1)))

print("=== 13-14. dimensions and PGA inside CGA ===")
check("bivectors: 10 = dim O(4,1)", len([b for b in range(cga.dim) if bin(b).count("1") == 2]) == 10)
gens = [F1, F2, F3, ni]
blades = []
for k in range(5):
    for idx in itertools.combinations(range(4), k):
        z = cga.scalar(1.0)
        for i in idx:
            z = cga.mul(z, gens[i])
        blades.append(z)
check("e1,e2,e3,ni generate 16 independent blades",
      np.linalg.matrix_rank(np.array(blades)) == 16)
check("ni^2 = 0 like e0", close(cga.mul(ni, ni), cga.zero()))
check("ni anticommutes with e_a",
      all(close(cga.mul(g, ni) + cga.mul(ni, g), cga.zero()) for g in (F1, F2, F3)))
# the CGA plane n + d ni and the PGA plane n - d e0 give the same mirror
cga_image = unpack(cga.mul(cga.mul(pi, cpoint(u)), pi))[1]
pga_plane = plane(tuple(nrm), -dist0)
pga_image = pga.mul(pga.mul(pga_plane, point(u)), pga_plane)
check("CGA plane n + d ni = PGA plane n - d e0", close(pga_image, point(cga_image)))

print("\nAll checks passed.")
