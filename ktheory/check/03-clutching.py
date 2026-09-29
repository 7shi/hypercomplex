"""Checks for ktheory/03-clutching.md (clutching functions and Hopf bundles).

Conventions: S^k is split into the northern hemisphere D_+ and the southern hemisphere
D_-; the point (x, v) of the trivial bundle over D_- is identified with (x, g(x) v) over
D_+ on the equator. An isomorphism E_g -> E_g' is given by A_+ and A_- with
g' = A_+ g A_-^{-1} on the equator.

1. Hopf bundles for K = R, C, H, O: with H(a, b) = (2 a b^*, |a|^2 - |b|^2), the fiber
   E_x = {v in K^2 | H(v) = |v|^2 x} over x = (c, r) in S^k is the real-linear subspace
   (1 - r) a = c b (r != 1), of real dimension dim K. On the equator a = c b, so the
   coordinates mu = a (on D_+) and lambda = b (on D_-) satisfy mu = L_c lambda: the
   clutching function is left multiplication by the unit number c, as in 01. For C this
   agrees with 01 (z = a / b, lambda (z, 1) = mu (1, w), mu = z lambda), and |z| < 1 is the
   southern hemisphere (r < 0).
2. The convention of hopf/04 (2 a^* b, left scalars): the conjugation (a, b) -> (a^*, b^*)
   maps our fiber over x to the fiber of hopf/04's map over the same x (as real bundles
   they are isomorphic). Its clutching is right multiplication by c^*, and conjugating the
   fibers turns it into L_c again: conj(conj(x) c^*) = c x. For C, as complex line bundles
   the clutching c^* has winding -1 (the conjugation is antilinear).
3. Winding numbers: winding of z^n is n, winding is additive under products, the
   "dog-walking" condition |g1 - g0| < |g0| keeps the winding, and a map with winding 0
   has a continuous logarithm (so it is homotopic to 1).
4. Tangent bundle of S^2: with the chart z (stereographic from the north pole) on D_- and
   w = 1/z on D_+, the components of a tangent vector satisfy v^w = -z^{-2} v^z, so the
   clutching function has winding -2. The field with v^z = 1 has v^w = -w^2: its only
   zero is the north pole, of order 2.
5. Evaluation at 1: L_c 1 = c, so x -> g(x) 1 is the identity map of S^{k-1}; adding a
   trivial summand puts it into a larger sphere where it is contractible (not checked
   numerically, only the formula).
6. det: diag(z, 1) and the Hopf bundle have det winding 1; E_g + trivial keeps det winding.
"""

import numpy as np

from common.octonion import L as OL
from common.octonion import R as OR

rng = np.random.default_rng(3)
dims = {1: "R", 2: "C", 4: "H", 8: "O"}


def conj(x):
    y = -x.copy()
    y[0] = x[0]
    return y


def mul(x, y, n):
    return sum(x[i] * OL[i][:n, :n] for i in range(n)) @ y


def Lm(c, n):
    return sum(c[i] * OL[i][:n, :n] for i in range(n))


def Rm(c, n):
    return sum(c[i] * OR[i][:n, :n] for i in range(n))


def hopf(a, b, n):
    return np.concatenate([2 * mul(a, conj(b), n), [a @ a - b @ b]])


print("1. Hopf bundles E_x = {v | H(v) = |v|^2 x}")
for n, name in dims.items():
    for _ in range(200):
        # random point x = (c, r) of S^n with r != 1 and a vector in the fiber
        x = rng.normal(size=n + 1)
        x /= np.linalg.norm(x)
        c, r = x[:n], x[n]
        b = rng.normal(size=n)
        a = mul(c, b, n) / (1 - r)
        v2 = a @ a + b @ b
        assert np.allclose(hopf(a, b, n), v2 * x)
        # converse: any v determines x = H(v)/|v|^2, and then (1 - r) a = c b
        a, b = rng.normal(size=n), rng.normal(size=n)
        v2 = a @ a + b @ b
        x = hopf(a, b, n) / v2
        assert np.isclose(np.linalg.norm(x), 1)
        c, r = x[:n], x[n]
        assert np.allclose((1 - r) * a, mul(c, b, n))
        # equator: a = c b, i.e. mu = L_c lambda
        c = rng.normal(size=n)
        c /= np.linalg.norm(c)
        b = rng.normal(size=n)
        a = Lm(c, n) @ b
        v2 = a @ a + b @ b
        assert np.allclose(hopf(a, b, n), v2 * np.concatenate([c, [0]]))
    print(f"  {name}: fiber (1-r) a = c b, H(v) = |v|^2 x, equator mu = L_c lambda: OK")

# complex case against 01: z = a/b, lambda (z, 1) = mu (1, w) with mu = z lambda
for _ in range(100):
    a, b = rng.normal(size=2) @ [1, 1j], rng.normal(size=2) @ [1, 1j]
    z, w = a / b, b / a
    lam, mu = b, a
    assert np.allclose([lam * z, lam], [a, b]) and np.allclose([mu, mu * w], [a, b])
    assert np.isclose(mu, z * lam)
    r = (abs(a) ** 2 - abs(b) ** 2) / (abs(a) ** 2 + abs(b) ** 2)
    assert (abs(z) < 1) == (r < 0)
    c = 2 * a * np.conj(b) / (abs(a) ** 2 + abs(b) ** 2)
    if np.isclose(abs(z), 1):
        assert np.isclose(c, z)
# on the equator |a| = |b|, 2 a b^* / |v|^2 = a / b
b = np.exp(0.3j)
a = np.exp(1.1j)
assert np.isclose(2 * a * np.conj(b) / 2, a / b)
print("  C: agrees with 01 (mu = z lambda, |z| < 1 is the southern hemisphere): OK")

print("2. hopf/04 convention (2 a^* b)")
for n, name in dims.items():
    for _ in range(100):
        c = rng.normal(size=n)
        c /= np.linalg.norm(c)
        b = rng.normal(size=n)
        a = Lm(c, n) @ b
        a2, b2 = conj(a), conj(b)  # conjugated pair
        # hopf/04's map on the conjugated pair gives the same base point c
        assert np.allclose(2 * mul(conj(a2), b2, n) / (a2 @ a2 + b2 @ b2), c)
        # coordinates: mu' = a2 = b2 c^* (right multiplication by c^*)
        assert np.allclose(a2, Rm(conj(c), n) @ b2)
        # conj R_{c^*} conj = L_c
        Cm = np.diag([1] + [-1] * (n - 1))
        assert np.allclose(Cm @ Rm(conj(c), n) @ Cm, Lm(c, n))
    print(f"  {name}: 2 a'^* b' = c for (a', b') = (a^*, b^*), clutching R_(c^*) ~ L_c: OK")


def winding(f, N=4000):
    t = np.linspace(0, 2 * np.pi, N + 1)
    ph = np.unwrap(np.angle(f(np.exp(1j * t))))
    return (ph[-1] - ph[0]) / (2 * np.pi)


print("3. winding numbers")
for n in range(-4, 5):
    assert np.isclose(winding(lambda z: z ** n), n)
assert np.isclose(winding(lambda z: -(z ** -2)), -2)
for _ in range(20):
    p = rng.normal(size=5) + 1j * rng.normal(size=5)
    q = rng.normal(size=5) + 1j * rng.normal(size=5)
    f = lambda z: np.polyval(p, z) * z ** -2
    g = lambda z: np.polyval(q, z) + 3 * z
    if min(abs(f(np.exp(1j * np.linspace(0, 7, 999))))) < 1e-2:
        continue
    assert np.isclose(winding(lambda z: f(z) * g(z)), winding(f) + winding(g), atol=1e-6)
print("  wind z^n = n, wind(-z^-2) = -2, additivity: OK")
# dog walking: |g1 - g0| < |g0| => same winding
for _ in range(50):
    g0 = lambda z: z ** 3 * (2 + 0.5 * z)
    eps = rng.normal(size=4) + 1j * rng.normal(size=4)
    g1 = lambda z: g0(z) + 0.2 * np.polyval(eps, z) / (1 + np.abs(np.polyval(eps, z)))
    t = np.exp(1j * np.linspace(0, 2 * np.pi, 2001))
    if np.all(abs(g1(t) - g0(t)) < abs(g0(t))):
        assert np.isclose(winding(g1), winding(g0))
        assert np.all((g1(t) / g0(t)).real > 0)
print("  |g1 - g0| < |g0| => g1/g0 has positive real part, same winding: OK")
# winding 0 => continuous log => homotopy g_s = exp(s log g) to 1
g = lambda z: (2 + z) * (3 + z ** -1) / (1 + 0.3 * z)
assert np.isclose(winding(g), 0)
t = np.linspace(0, 2 * np.pi, 4001)
lg = np.log(abs(g(np.exp(1j * t)))) + 1j * np.unwrap(np.angle(g(np.exp(1j * t))))
assert np.isclose(lg[0], lg[-1])
print("  winding 0 => periodic continuous log: OK")

print("4. tangent bundle of S^2")


def stereo(z):  # inverse stereographic projection from the north pole (hopf/homog)
    d = 1 + abs(z) ** 2
    return np.array([2 * z.real / d, 2 * z.imag / d, (abs(z) ** 2 - 1) / d])


def push(chart_inv, u, v, h=1e-6):  # tangent vector of the curve chart_inv(u + t v)
    return (chart_inv(u + h * v) - chart_inv(u - h * v)) / (2 * h)


for _ in range(50):
    z = rng.normal() + 1j * rng.normal()
    vz = rng.normal() + 1j * rng.normal()
    w = 1 / z
    vw = -(z ** -2) * vz
    t1 = push(stereo, z, vz)
    t2 = push(lambda w: stereo(1 / w), w, vw)
    assert np.allclose(t1, t2, atol=1e-6)
    assert np.isclose(t1 @ stereo(z), 0, atol=1e-6)
print("  v^w = -z^-2 v^z gives the same tangent vector of S^2: OK")
assert stereo(0j)[2] == -1 and abs(stereo(1e6 + 0j)[2] - 1) < 1e-9
print("  z = 0 is the south pole, z -> infinity the north pole: OK")
for _ in range(20):
    w = rng.normal() + 1j * rng.normal()
    assert np.isclose(-(w ** 2), -((1 / w) ** -2) * 1)
print("  v^z = 1 gives v^w = -w^2 (double zero at the north pole): OK")

print("5. evaluation at 1")
for n, name in dims.items():
    c = rng.normal(size=n)
    e = np.zeros(n)
    e[0] = 1
    assert np.allclose(Lm(c, n) @ e, c)
print("  L_c 1 = c: OK")

print("6. determinant")
det = np.vectorize(lambda z, f: np.linalg.det(f(z)), excluded={1})
assert np.isclose(winding(lambda z: det(z, lambda z: np.diag([z, 1]))), 1)
M = lambda z: np.array([[2, 1j], [0, 1]]) @ np.diag([z, 1]) @ np.array([[1, 0], [3, 1]])
assert np.isclose(winding(lambda z: det(z, M)), 1)
print("  wind det diag(z, 1) = 1, unchanged by constant factors: OK")
print("All checks passed.")
