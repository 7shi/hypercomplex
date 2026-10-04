"""Checks for c2-to-s2.md.

Verifies symbolically (sympy) and numerically (random points, numpy):

1. z0* z1 is invariant under the global phase (z0, z1) -> e^{iω}(z0, z1),
   hence so are Re, Im, |z0|², |z1|².
2. With x = c Re(z0*z1), y = c Im(z0*z1), z = |z0|²-|z1|²,
   x²+y²+z² = |z0|⁴ + (c²-2)|z0|²|z1|² + |z1|⁴, so c = ±2 puts the image
   on the unit sphere.
3. The representative (cos θ, e^{iφ} sin θ) maps to
   (cos φ sin 2θ, sin φ sin 2θ, cos 2θ), and the general
   (e^{iα} cos θ, e^{iβ} sin θ) gives the same with φ = β-α.
4. The second form of the Hopf map
   (z0*z1 + z1*z0, -i(z0*z1 - z1*z0), z0*z0 - z1*z1) equals the first.
5. Surjectivity: for random points on S² the representative built from
   θ = arccos(z)/2, φ = atan2(y, x) maps back to the same point, and
   the poles are hit by (1, 0) and (0, 1).
6. Injectivity up to phase: for z > -1 the representative with positive
   real first component is recovered from the image as
   w0 = sqrt((1+z)/2), w1 = (x+iy)/(2 w0), and equals e^{-i arg z0}(z0, z1);
   for z = -1 every pair (0, e^{iβ}) is equivalent to (0, 1).
7. U² + V² = |z0|²|z1|² for U = Re(z0*z1), V = Im(z0*z1).
"""

import numpy as np
import sympy as sp

# 1. phase invariance
a0, a1, b0, b1, w = sp.symbols("a0 a1 b0 b1 omega", real=True)
z0 = a0 + sp.I * a1
z1 = b0 + sp.I * b1
p = sp.exp(sp.I * w)
inv = sp.conjugate(z0) * z1
inv2 = sp.conjugate(p * z0) * (p * z1)
assert sp.simplify(sp.expand(inv2 - inv)) == 0
print("1. z0*z1 is invariant under the global phase: OK")

# 2. coefficient c
c = sp.symbols("c", real=True)
x = c * sp.re(inv)
y = c * sp.im(inv)
z = sp.expand(z0 * sp.conjugate(z0) - z1 * sp.conjugate(z1))
n0 = a0**2 + a1**2
n1 = b0**2 + b1**2
lhs = sp.expand(x**2 + y**2 + z**2)
rhs = sp.expand(n0**2 + (c**2 - 2) * n0 * n1 + n1**2)
assert sp.simplify(lhs - rhs) == 0
assert sp.expand(rhs.subs(c, 2) - (n0 + n1) ** 2) == 0
print("2. x²+y²+z² = |z0|⁴+(c²-2)|z0|²|z1|²+|z1|⁴, c=2 gives (|z0|²+|z1|²)²: OK")

# 3. angle form
th, ph, al, be = sp.symbols("theta phi alpha beta", real=True)


def hopf(u0, u1):
    q = sp.conjugate(u0) * u1
    return [2 * sp.re(q), 2 * sp.im(q), sp.Abs(u0) ** 2 - sp.Abs(u1) ** 2]


target = [sp.cos(ph) * sp.sin(2 * th), sp.sin(ph) * sp.sin(2 * th), sp.cos(2 * th)]
got = hopf(sp.cos(th), sp.exp(sp.I * ph) * sp.sin(th))
assert all(sp.simplify(sp.expand_trig(g - t)) == 0 for g, t in zip(got, target))
got = hopf(sp.exp(sp.I * al) * sp.cos(th), sp.exp(sp.I * be) * sp.sin(th))
target_ab = [t.subs(ph, be - al) for t in target]
assert all(sp.simplify(sp.expand_trig(g - t)) == 0 for g, t in zip(got, target_ab))
print("3. angle form (cos φ sin 2θ, sin φ sin 2θ, cos 2θ), φ = β-α: OK")

# 4. second form
cz0, cz1 = sp.conjugate(z0), sp.conjugate(z1)
form2 = [cz0 * z1 + cz1 * z0, -sp.I * (cz0 * z1 - cz1 * z0), cz0 * z0 - cz1 * z1]
form1 = [2 * sp.re(inv), 2 * sp.im(inv), n0 - n1]
assert all(sp.simplify(sp.expand(f2 - f1)) == 0 for f1, f2 in zip(form1, form2))
print("4. second form equals the first: OK")

# 5. surjectivity (numerical)
rng = np.random.default_rng(42)


def hopf_np(u0, u1):
    q = np.conj(u0) * u1
    return np.array([2 * q.real, 2 * q.imag, abs(u0) ** 2 - abs(u1) ** 2])


for _ in range(1000):
    v = rng.normal(size=3)
    v /= np.linalg.norm(v)
    t = np.arccos(np.clip(v[2], -1, 1)) / 2
    f = np.arctan2(v[1], v[0])
    assert 0 <= t <= np.pi / 2
    assert np.allclose(hopf_np(np.cos(t), np.exp(1j * f) * np.sin(t)), v)
assert np.allclose(hopf_np(1, 0), [0, 0, 1])
assert np.allclose(hopf_np(0, 1), [0, 0, -1])
print("5. every point of S² is hit with θ ∈ [0, π/2]; poles from (1,0), (0,1): OK")

# 6. injectivity up to the global phase
for _ in range(1000):
    v = rng.normal(size=4)
    v /= np.linalg.norm(v)
    u0, u1 = v[0] + 1j * v[1], v[2] + 1j * v[3]
    x, y, z = hopf_np(u0, u1)
    w0 = np.sqrt((1 + z) / 2)
    w1 = (x + 1j * y) / (2 * w0)
    ph = np.exp(-1j * np.angle(u0))
    assert np.allclose([w0, w1], [ph * u0, ph * u1])
for b in rng.uniform(0, 2 * np.pi, 10):
    assert np.allclose(hopf_np(0, np.exp(1j * b)), [0, 0, -1])
print("6. the image determines the representative (z > -1); z = -1 is one class: OK")

# 7. relation among the invariants
assert sp.simplify(sp.expand(sp.re(inv) ** 2 + sp.im(inv) ** 2 - n0 * n1)) == 0
print("7. Re(z0*z1)² + Im(z0*z1)² = |z0|²|z1|²: OK")
