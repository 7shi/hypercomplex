"""Checks for article 02 (directed integration and the fundamental theorem).

Sign and orientation conventions of the fundamental theorem of geometric
calculus, checked by exact integration over concrete regions.

Cl_{2,0}(R) is realized as M_2(R) (every real 2x2 matrix is an element), and
Cl_{3,0}(R) as M_2(C) with the Pauli matrices as e1, e2, e3. A generic
algebra-valued function is therefore a matrix of polynomials.

2D (boundary counterclockwise, n outward unit normal, ds arc length,
dx = e1 dx + e2 dy the directed line element):

1. int_M D F dA = oint n F ds on the unit square and the unit disk.
2. n ds = I dx, so oint I dx F = int D F dA and oint dx F = -I int D F dA.
3. Bivector form: with dX = I dA (I = e1 e2, the orientation of the
   counterclockwise boundary), oint dx F = int (D . dX) F = -int (dX . D) F,
   where a . B = (a B - B a) / 2 = -(B . a) and D differentiates F only.
   Since D . I = -I D this is item 2 again.
4. For a vector F = a e1 + b e2 the scalar part of item 1 is the flux form
   of Green's theorem, int (a_x + b_y) = oint (a dy - b dx), and the
   bivector part is the circulation form, int (b_x - a_y) = oint (a dx + b dy).
5. Two-sided form: int [(G D) F + G (D F)] dA = oint G n F ds.
6. Cauchy: for monogenic even F, oint dz F = 0 with dz = e1 dx.

3D:

7. int_V D F dV = oint n F dS on the unit cube and the unit ball.
8. For a vector F the scalar part is the divergence theorem and the
   bivector part is int curl F dV = oint n x F dS.
9. Surface patch r(u, v) = (u, v, u v), 0 <= u, v <= 1, with
   dX = r_u ^ r_v du dv and the boundary counterclockwise in (u, v):
   oint dx F = int (D . dX) F = -int (dX . D) F, the same sign as item 3.
   Its scalar part for a vector F is Stokes' theorem.
"""

import sympy as sp

x, y, z, t, r, th, ph = sp.symbols("x y z t r theta phi", real=True)


def simp(M):
    return M.applyfunc(lambda e: sp.expand(sp.simplify(e)))


def integ(M, *lims):
    return simp(M.applyfunc(lambda e: sp.integrate(sp.expand(e), *lims)))


def eq(A, B):
    return simp(A - B) == sp.zeros(*A.shape)


def wedge_dot(B, a):
    """Inner product of a bivector B with a vector a: (B a - a B) / 2."""
    return (B * a - a * B) / 2


# ---------------------------------------------------------------- 2D
one2 = sp.eye(2)
e1 = sp.Matrix([[1, 0], [0, -1]])
e2 = sp.Matrix([[0, 1], [1, 0]])
I2 = e1 * e2
assert I2**2 == -one2


def D2(F):
    return e1 * F.diff(x) + e2 * F.diff(y)


def D2r(F):
    return F.diff(x) * e1 + F.diff(y) * e2


def poly2(seed):
    """A generic real 2x2 matrix whose entries are cubic polynomials."""
    import random

    rnd = random.Random(seed)
    mons = [x**i * y**j for i in range(4) for j in range(4 - i)]
    return sp.Matrix(2, 2, lambda i, j: sum(rnd.randint(-3, 3) * m for m in mons))


F2 = poly2(1)
G2 = poly2(2)


def square_boundary(integrand):
    """oint over the counterclockwise unit square boundary.

    integrand(pt, dxv, nds) returns the matrix to integrate in t in [0, 1],
    where pt = (x, y) on the edge, dxv = directed line element per dt and
    nds = outward normal times ds per dt."""
    edges = [
        ((t, 0), (1, 0)),  # bottom, left to right
        ((1, t), (0, 1)),  # right, upward
        ((1 - t, 1), (-1, 0)),  # top, right to left
        ((0, 1 - t), (0, -1)),  # left, downward
    ]
    total = sp.zeros(2, 2)
    for (px, py), (tx, ty) in edges:
        dxv = tx * e1 + ty * e2
        nds = ty * e1 - tx * e2  # (dy, -dx)
        total += integ(integrand((px, py), dxv, nds), (t, 0, 1))
    return simp(total)


def disk_boundary(integrand):
    px, py = sp.cos(t), sp.sin(t)
    dxv = -sp.sin(t) * e1 + sp.cos(t) * e2
    nds = sp.cos(t) * e1 + sp.sin(t) * e2
    return integ(integrand((px, py), dxv, nds), (t, 0, 2 * sp.pi))


def at(F, pt):
    return F.subs({x: pt[0], y: pt[1]}, simultaneous=True)


def over_square(M):
    return integ(M, (x, 0, 1), (y, 0, 1))


def over_disk(M):
    Mp = M.subs({x: r * sp.cos(th), y: r * sp.sin(th)}, simultaneous=True) * r
    return integ(Mp, (r, 0, 1), (th, 0, 2 * sp.pi))


regions2 = [("square", square_boundary, over_square), ("disk", disk_boundary, over_disk)]

# 1. int D F dA = oint n F ds
for name, bd, area in regions2:
    lhs = area(D2(F2))
    rhs = bd(lambda p, dxv, nds: nds * at(F2, p))
    assert eq(lhs, rhs), name
print("1. int D F dA = oint n F ds (square, disk)")

# 2. n ds = I dx
for name, bd, area in regions2:
    rhs = bd(lambda p, dxv, nds: I2 * dxv * at(F2, p))
    assert eq(area(D2(F2)), rhs), name
    rhs = bd(lambda p, dxv, nds: dxv * at(F2, p))
    assert eq(-I2 * area(D2(F2)), rhs), name
print("2. n ds = I dx, oint dx F = -I int D F dA")

# 3. bivector-dot form: which sign
for name, bd, area in regions2:
    rhs = bd(lambda p, dxv, nds: dxv * at(F2, p))
    dX = I2  # e1 e2
    dotD = lambda F: wedge_dot(dX, e1) * F.diff(x) + wedge_dot(dX, e2) * F.diff(y)
    assert eq(wedge_dot(dX, e1), I2 * e1)  # I . e_k = I e_k
    assert eq(area(dotD(F2)), -rhs), name
print("3. oint dx F = int (D . dX) F = -int (dX . D) F for dX = e1 e2 dA")

# 4. Green's theorem, two forms
av = 1 + x**2 * y + 2 * y**3
bv = x * y - 3 * x**3 + y**2
Fv = av * e1 + bv * e2
DFv = simp(D2(Fv))
assert eq(DFv, (av.diff(x) + bv.diff(y)) * one2 + (bv.diff(x) - av.diff(y)) * I2)
for name, bd, area in regions2:
    lhs = area(DFv)
    rhs = bd(lambda p, dxv, nds: nds * at(Fv, p))
    assert eq(lhs, rhs)
    scal = lambda M: sp.simplify(M.trace() / 2)
    biv = lambda M: sp.simplify((M * I2.inv()).trace() / 2)
    ap = lambda p: av.subs({x: p[0], y: p[1]}, simultaneous=True)
    bp = lambda p: bv.subs({x: p[0], y: p[1]}, simultaneous=True)
    # boundary parts in the classical forms: n ds = e1 dy - e2 dx
    flux = bd(lambda p, dxv, nds: (ap(p) * nds[0, 0] + bp(p) * nds[0, 1]) * one2)[0, 0]
    circ = bd(lambda p, dxv, nds: (ap(p) * dxv[0, 0] + bp(p) * dxv[0, 1]) * one2)[0, 0]
    assert sp.simplify(scal(lhs) - flux) == 0, name
    assert sp.simplify(biv(lhs) - circ) == 0, name
print("4. scalar part = flux form, bivector part = circulation form of Green")

# 5. two-sided form
for name, bd, area in regions2:
    lhs = area(D2r(G2) * F2 + G2 * D2(F2))
    rhs = bd(lambda p, dxv, nds: at(G2, p) * nds * at(F2, p))
    assert eq(lhs, rhs), name
print("5. int [(G D) F + G (D F)] dA = oint G n F ds")

# 6. Cauchy's theorem via dz = e1 dx
zz = x * one2 + y * I2
F3 = simp(zz**3)
assert D2(F3) == sp.zeros(2, 2)
for name, bd, area in regions2:
    assert eq(bd(lambda p, dxv, nds: e1 * dxv * at(F3, p)), sp.zeros(2, 2)), name
assert eq(e1 * (e1 * 1 + e2 * 0), one2) and eq(e1 * e2, I2)  # e1 dx = dx + I dy
print("6. monogenic F: oint dz F = 0 with dz = e1 dx")

# ---------------------------------------------------------------- 3D
one3 = sp.eye(2)
s1 = sp.Matrix([[0, 1], [1, 0]])
s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
s3 = sp.Matrix([[1, 0], [0, -1]])
E = [s1, s2, s3]
I3 = s1 * s2 * s3
assert I3 == sp.I * one3
X = (x, y, z)


def D3(F):
    return sum((E[k] * F.diff(X[k]) for k in range(3)), sp.zeros(2, 2))


def poly3(seed):
    import random

    rnd = random.Random(seed)
    mons = [x**i * y**j * z**k for i in range(3) for j in range(3 - i) for k in range(3 - i - j)]
    return sp.Matrix(
        2,
        2,
        lambda i, j: sum((rnd.randint(-3, 3) + sp.I * rnd.randint(-3, 3)) * m for m in mons),
    )


F3g = poly3(3)


def at3(F, pt):
    return F.subs(dict(zip(X, pt)), simultaneous=True)


def vec3(v):
    return sum((v[k] * E[k] for k in range(3)), sp.zeros(2, 2))


u, v = sp.symbols("u v", real=True)


def cube_boundary(F):
    total = sp.zeros(2, 2)
    for k in range(3):
        others = [i for i in range(3) if i != k]
        for val, sgn in ((1, 1), (0, -1)):
            pt = [None] * 3
            pt[k] = val
            pt[others[0]], pt[others[1]] = u, v
            n = [0, 0, 0]
            n[k] = sgn
            total += integ(vec3(n) * at3(F, pt), (u, 0, 1), (v, 0, 1))
    return simp(total)


def cube_volume(M):
    return integ(M, (x, 0, 1), (y, 0, 1), (z, 0, 1))


def ball_boundary(F):
    n = [sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)]
    return integ(vec3(n) * at3(F, n) * sp.sin(th), (th, 0, sp.pi), (ph, 0, 2 * sp.pi))


def ball_volume(M):
    pt = [r * sp.sin(th) * sp.cos(ph), r * sp.sin(th) * sp.sin(ph), r * sp.cos(th)]
    Mp = at3(M, pt) * r**2 * sp.sin(th)
    return integ(Mp, (r, 0, 1), (th, 0, sp.pi), (ph, 0, 2 * sp.pi))


# 7. int D F dV = oint n F dS
assert eq(cube_volume(D3(F3g)), cube_boundary(F3g))
assert eq(ball_volume(D3(F3g)), ball_boundary(F3g))
print("7. int D F dV = oint n F dS (cube, ball)")

# 8. divergence and curl for a vector field
fv = [x * y + z**2, y * z - x**2, x * z + y**3]
Fv3 = vec3(fv)
div = sum(fv[k].diff(X[k]) for k in range(3))
curl = [
    fv[2].diff(y) - fv[1].diff(z),
    fv[0].diff(z) - fv[2].diff(x),
    fv[1].diff(x) - fv[0].diff(y),
]
# D F = div + I curl  (D ^ F = I (D x F))
assert eq(D3(Fv3), div * one3 + I3 * vec3(curl))
lhs = cube_volume(D3(Fv3))
rhs = cube_boundary(Fv3)
assert eq(lhs, rhs)
print("8. D F = div F + I curl F; scalar part = divergence theorem,")
print("   bivector part = int curl F dV = oint n x F dS")

# 9. surface patch and Stokes
rr = sp.Matrix([u, v, u * v])
ru, rv = rr.diff(u), rr.diff(v)
dX = (vec3(ru) * vec3(rv) - vec3(rv) * vec3(ru)) / 2  # r_u ^ r_v
Gs = poly3(4)


def dotD(F):
    return sum((wedge_dot(dX, E[k]) * F.diff(X[k]) for k in range(3)), sp.zeros(2, 2))


surf = integ(at3(dotD(Gs), list(rr)), (u, 0, 1), (v, 0, 1))
edges = [((t, 0), (1, 0)), ((1, t), (0, 1)), ((1 - t, 1), (-1, 0)), ((0, 1 - t), (0, -1))]
line = sp.zeros(2, 2)
for (pu, pv), _ in edges:
    pt = [e.subs({u: pu, v: pv}, simultaneous=True) for e in rr]
    dpt = [sp.diff(e, t) for e in pt]
    line += integ(vec3(dpt) * at3(Gs, pt), (t, 0, 1))
line = simp(line)
if eq(line, surf):
    sign = 1
else:
    assert eq(line, -surf)
    sign = -1
assert sign == -1  # same convention as item 3
print("9. surface: oint dx F = int (D . dX) F with dX = r_u ^ r_v du dv")

# Stokes for a vector field: scalar part
line_s = sp.zeros(2, 2)
for (pu, pv), _ in edges:
    pt = [e.subs({u: pu, v: pv}, simultaneous=True) for e in rr]
    dpt = [sp.diff(e, t) for e in pt]
    fp = [f.subs(dict(zip(X, pt)), simultaneous=True) for f in fv]
    line_s += integ(sum(fp[k] * dpt[k] for k in range(3)) * one3, (t, 0, 1))
nvec = ru.cross(rv)
cp = [c.subs(dict(zip(X, list(rr))), simultaneous=True) for c in curl]
flux_curl = sp.integrate(sp.expand(sum(cp[k] * nvec[k] for k in range(3))), (u, 0, 1), (v, 0, 1))
assert sp.simplify(line_s[0, 0] - flux_curl) == 0
surf_v = integ(at3(dotD(Fv3), list(rr)), (u, 0, 1), (v, 0, 1))
assert sp.simplify((sign * surf_v).trace() / 2 - flux_curl) == 0
print("   scalar part = Stokes' theorem oint F . dx = int curl F . n dA")

print("all checks passed")
