"""Checks for 04-extension.md (fiber structure of the octonionic Hopf map).

Verifies the concrete example in the article:

1. Left multiplication by a unit octonion q does NOT preserve the fiber:
   H(qα, qβ) != H(α, β) for α = e1/√2, β = e4/√2, q = e2.
2. The fiber over (p, 0) is parametrized by right multiplication:
   β = αp gives H(α, β) = (p, 0) for any α with |α| = 1/√2,
   for a general unit octonion p (real part allowed).
3. The fiber over a general point (c, r) is given by β = αc/(1+r).
"""

import numpy as np
from common.octonion import L

def mul(x, y):
    M = sum(x[i] * L[i] for i in range(8))
    return M @ y

def conj(x):
    y = -x.copy()
    y[0] = x[0]
    return y

def H_O(a, b):
    # returns (2*conj(a)*b, |a|^2 - |b|^2)
    return 2 * mul(conj(a), b), a @ a - b @ b

def e(i):
    v = np.zeros(8)
    v[i] = 1.0
    return v

# alpha = e1/sqrt2, beta = e4/sqrt2, q = e2 (unit octonion)
s2 = np.sqrt(2)
alpha = e(1) / s2
beta = e(4) / s2
q = e(2)

img, real = H_O(alpha, beta)
print("H_O(alpha, beta)      :", img, real)
assert np.allclose(img, -e(5)) and np.isclose(real, 0)

qa = mul(q, alpha)
qb = mul(q, beta)
print("q*alpha, q*beta       :", qa * s2, qb * s2)
assert np.allclose(qa, -e(3) / s2) and np.allclose(qb, e(6) / s2)
img2, real2 = H_O(qa, qb)
print("H_O(q*alpha, q*beta)  :", img2, real2)
assert np.allclose(img2, e(5)) and np.isclose(real2, 0)

print("norms preserved:", np.linalg.norm(qa) - np.linalg.norm(alpha), np.linalg.norm(qb) - np.linalg.norm(beta))
print("same fiber?", np.allclose(img, img2) and np.allclose(real, real2))

print()
print("--- right-mult parametrization of the same fiber ---")

def fiber_point(alpha, p):
    # beta = alpha * p  (p: general unit octonion = target point (p, 0) on S^8)
    return alpha, mul(alpha, p)

p = -e(5)  # target point on S^8 found above: H_O(alpha,beta) = (-e5, 0)

for label, a in [("alpha = e1/sqrt2", alpha), ("a2 = (e1+e2)/2", (e(1)+e(2))/2)]:
    al, be = fiber_point(a, p)
    img, real = H_O(al, be)
    print(f"{label:20s} -> H_O = {img}, {real}")
    assert np.allclose(img, p) and np.isclose(real, 0)

print()
print("--- same construction with p having a real part ---")

rng0 = np.random.default_rng(1)
p2 = rng0.normal(size=8)
p2 /= np.linalg.norm(p2)  # general unit octonion, p2[0] != 0
assert abs(p2[0]) > 1e-3

for label, a in [("alpha = e1/sqrt2", alpha), ("a2 = (e1+e2)/2", (e(1)+e(2))/2)]:
    al, be = fiber_point(a, p2)
    img, real = H_O(al, be)
    ok = np.allclose(img, p2) and np.isclose(real, 0)
    print(f"{label:20s} -> H_O = (p2, 0)? {ok}")
    assert ok

print()
print("--- uniqueness: alpha (alpha* beta) = |alpha|^2 beta recovers beta = alpha p ---")

rng1 = np.random.default_rng(2)
for _ in range(3):
    a = rng1.normal(size=8)
    b = rng1.normal(size=8)
    n = np.sqrt(a @ a + b @ b)
    a /= n; b /= n
    ok = np.allclose(mul(a, mul(conj(a), b)), (a @ a) * b)
    print("x(x* y) = |x|^2 y ?", ok)
    assert ok

print()
print("--- general point (c, r) on S^8: |alpha|^2 = (1+r)/2, beta = alpha c / (1+r) ---")

rng = np.random.default_rng(0)
for _ in range(3):
    c = rng.normal(size=8)
    r = rng.uniform(-0.9, 0.9)
    c *= np.sqrt(1 - r**2) / np.linalg.norm(c)  # |c|^2 + r^2 = 1
    a = rng.normal(size=8)
    a *= np.sqrt((1 + r) / 2) / np.linalg.norm(a)
    beta = mul(a, c) / (1 + r)
    img, real = H_O(a, beta)
    ok = np.allclose(img, c) and np.isclose(real, r)
    print(f"target r = {r:+.4f} -> H_O matches (c, r)? {ok}")
    assert ok
