"""Checks for homogeneous.md.

Verifies symbolically (sympy) and numerically (random points, numpy):

1. z = (x1+ix2)/(x3+ix4) has real part u = (x1x3+x2x4)/(x3²+x4²) and
   imaginary part v = (x2x3-x1x4)/(x3²+x4²).
2. The line (ut, vt, 1-t) from the north pole meets the unit sphere at
   t = 0 and t = 2/(u²+v²+1).
3. On S³, u²+v²+1 = 1/(x3²+x4²), and the intersection point is
   ξ = (2(x1x3+x2x4), 2(x2x3-x1x4), x1²+x2²-x3²-x4²), which lies on S².
4. As x3+ix4 -> 0 on S³, ξ -> (0, 0, 1).
5. With α = x1+ix2, β = x3+ix4: ξ1+iξ2 = 2αβ*, ξ3 = αα*-ββ*, and the
   form (αβ*+βα*, -i(αβ*-βα*), αα*-ββ*) agrees.
6. The remark: w = β/α projected from the south pole gives
   (2w, 1-|w|²)/(1+|w|²) = (2α*β, αα*-ββ*), differing from ξ only in the
   sign of ξ2.
7. ξ depends only on the ratio: multiplying (α, β) by e^{iω} leaves ξ fixed.
8. Fiber: two points of S³ with the same image differ by a phase. The
   image determines [α:β] (z = (ξ1+iξ2)/(1-ξ3) for ξ3 < 1, and β = 0 at
   the north pole), and a ratio-preserving λ between points of S³ has |λ| = 1.
9. The south-pole stereographic coordinate of the same point is 1/z*,
   not w = 1/z.
10. The boundary cases: β = 0 maps to (0, 0, 1), and in the remark α = 0
    gives (2α*β, |α|²-|β|²) = (0, -1), the south pole.
"""

import numpy as np
import sympy as sp

x1, x2, x3, x4, t = sp.symbols("x1 x2 x3 x4 t", real=True)
z = (x1 + sp.I * x2) / (x3 + sp.I * x4)
d = x3**2 + x4**2
u = (x1 * x3 + x2 * x4) / d
v = (x2 * x3 - x1 * x4) / d

# 1. real and imaginary parts
assert sp.simplify(sp.re(z) - u) == 0
assert sp.simplify(sp.im(z) - v) == 0
print("1. u, v are the real and imaginary parts of z: OK")

# 2. intersection with the unit sphere
U, V = sp.symbols("U V", real=True)
sols = sp.solve(sp.Eq((U * t) ** 2 + (V * t) ** 2 + (1 - t) ** 2, 1), t)
assert set(sols) == {0, 2 / (U**2 + V**2 + 1)}
print("2. t = 0, 2/(u²+v²+1): OK")

# 3. on S³ (substitute x1² = 1 - x2² - x3² - x4²)
s3 = {x1**2: 1 - x2**2 - x3**2 - x4**2}
den = sp.simplify(sp.expand(sp.numer(sp.together(u**2 + v**2 + 1))).subs(s3)
      / sp.denom(sp.together(u**2 + v**2 + 1)))
assert sp.simplify(den - 1 / d) == 0
k = 1 / d  # u²+v²+1 on S³
xi = [2 * u / k, 2 * v / k, 1 - 2 / k]
target = [2 * (x1 * x3 + x2 * x4), 2 * (x2 * x3 - x1 * x4),
          x1**2 + x2**2 - x3**2 - x4**2]
assert all(sp.simplify(a - b) == 0 for a, b in zip(xi[:2], target[:2]))
assert sp.simplify(sp.expand(xi[2] - target[2]).subs(s3)) == 0
norm = sp.expand(sum(c**2 for c in target) - (x1**2 + x2**2 + x3**2 + x4**2) ** 2)
assert norm == 0
print("3. u²+v²+1 = 1/(x3²+x4²) and ξ formula; |ξ|² = |x|⁴ = 1: OK")


def xi_np(p):
    a1, a2, a3, a4 = p
    return np.array([2 * (a1 * a3 + a2 * a4), 2 * (a2 * a3 - a1 * a4),
                     a1**2 + a2**2 - a3**2 - a4**2])


# 4. limit to the north pole
rng = np.random.default_rng(42)
for eps in [1e-2, 1e-4, 1e-6]:
    p = rng.normal(size=4)
    p[2:] *= eps / np.linalg.norm(p[2:])
    p[:2] *= np.sqrt(1 - eps**2) / np.linalg.norm(p[:2])
    assert np.allclose(xi_np(p), [0, 0, 1], atol=3 * eps)
print("4. x3+ix4 -> 0 gives ξ -> (0, 0, 1): OK")

# 5. complex pair form
al = x1 + sp.I * x2
be = x3 + sp.I * x4
c = 2 * al * sp.conjugate(be)
assert sp.simplify(sp.expand(sp.re(c)) - target[0]) == 0
assert sp.simplify(sp.expand(sp.im(c)) - target[1]) == 0
assert sp.expand(al * sp.conjugate(al) - be * sp.conjugate(be) - target[2]) == 0
f2 = [al * sp.conjugate(be) + be * sp.conjugate(al),
      -sp.I * (al * sp.conjugate(be) - be * sp.conjugate(al))]
assert all(sp.expand(a - b) == 0 for a, b in zip(f2, target[:2]))
print("5. ξ1+iξ2 = 2αβ*, ξ3 = αα*-ββ*, second form: OK")

# 6. remark: w = β/α from the south pole
for _ in range(1000):
    p = rng.normal(size=4)
    p /= np.linalg.norm(p)
    a, b = p[0] + 1j * p[1], p[2] + 1j * p[3]
    w = b / a
    # line from (0,0,-1) to (Re w, Im w, 0): (Re w s, Im w s, -1+s)
    s = 2 / (1 + abs(w) ** 2)
    q = np.array([w.real * s, w.imag * s, -1 + s])
    assert np.isclose(np.linalg.norm(q), 1)
    lhs = np.array([2 * w, 1 - abs(w) ** 2]) / (1 + abs(w) ** 2)
    assert np.allclose(lhs, [2 * np.conj(a) * b, abs(a) ** 2 - abs(b) ** 2])
    assert np.allclose(q, [lhs[0].real, lhs[0].imag, lhs[1].real])
    assert np.isclose(1 + abs(w) ** 2, 1 / abs(a) ** 2)
    x = xi_np(p)
    assert np.allclose(q, [x[0], -x[1], x[2]])
print("6. β/α from the south pole gives (2α*β, |α|²-|β|²), ξ2 sign flipped: OK")

# 7. phase invariance
for _ in range(100):
    p = rng.normal(size=4)
    p /= np.linalg.norm(p)
    a, b = p[0] + 1j * p[1], p[2] + 1j * p[3]
    e = np.exp(1j * rng.uniform(0, 2 * np.pi))
    a2, b2 = e * a, e * b
    assert np.allclose(xi_np(p), xi_np([a2.real, a2.imag, b2.real, b2.imag]))
print("7. ξ is invariant under (α, β) -> e^{iω}(α, β): OK")

# 8. fiber
for _ in range(1000):
    p = rng.normal(size=4)
    p /= np.linalg.norm(p)
    a, b = p[0] + 1j * p[1], p[2] + 1j * p[3]
    x = xi_np(p)
    zr = (x[0] + 1j * x[1]) / (1 - x[2])
    assert np.isclose(zr, a / b)
    # another point with the same ratio, rescaled onto S³
    lam = rng.normal() + 1j * rng.normal()
    a2, b2 = lam * a, lam * b
    n = np.sqrt(abs(a2) ** 2 + abs(b2) ** 2)
    a2, b2 = a2 / n, b2 / n
    l2 = a2 / a
    assert np.isclose(abs(l2), 1) and np.isclose(b2, l2 * b)
    assert np.allclose(xi_np([a2.real, a2.imag, b2.real, b2.imag]), x)
print("8. same image <=> same ratio <=> differ by a phase on S³: OK")

# 9. south-pole coordinate of the same point is 1/z*
for _ in range(1000):
    p = rng.normal(size=4)
    p /= np.linalg.norm(p)
    a, b = p[0] + 1j * p[1], p[2] + 1j * p[3]
    x = xi_np(p)
    zn = (x[0] + 1j * x[1]) / (1 - x[2])
    zs = (x[0] + 1j * x[1]) / (1 + x[2])
    assert np.isclose(zs, 1 / np.conj(zn))
print("9. south-pole coordinate of the same point is 1/z*: OK")

# 10. boundary cases
for f in rng.uniform(0, 2 * np.pi, 10):
    assert np.allclose(xi_np([np.cos(f), np.sin(f), 0, 0]), [0, 0, 1])
    b = np.exp(1j * f)
    assert np.allclose([2 * np.conj(0) * b, 0 - abs(b) ** 2], [0, -1])
print("10. β = 0 -> north pole, α = 0 -> south pole in the remark: OK")
