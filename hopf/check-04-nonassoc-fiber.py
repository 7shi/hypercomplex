"""Checks for 04-extension.md (fiber structure of the octonionic Hopf map).

Verifies the concrete example in the article:

1. Right multiplication by a unit octonion q does NOT preserve the fiber:
   H(αq, βq) != H(α, β) for α = e1/√2, β = e4/√2, q = e2.
2. The fiber over (p, 0) is parametrized by left multiplication:
   β = -pα gives H(α, β) = (p, 0) for any α with |α| = 1/√2.
"""

import sys
import numpy as np
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent / 'oct'))
from octonion import L

def mul(x, y):
    M = sum(x[i] * L[i] for i in range(8))
    return M @ y

def conj(x):
    y = -x.copy()
    y[0] = x[0]
    return y

def H_O(a, b):
    # returns (2*a*conj(b), |a|^2 - |b|^2)
    return 2 * mul(a, conj(b)), a @ a - b @ b

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

aq = mul(alpha, q)
bq = mul(beta, q)
img2, real2 = H_O(aq, bq)
print("H_O(alpha*q, beta*q)  :", img2, real2)

print("norms preserved:", np.linalg.norm(aq) - np.linalg.norm(alpha), np.linalg.norm(bq) - np.linalg.norm(beta))
print("same fiber?", np.allclose(img, img2) and np.allclose(real, real2))

print()
print("--- external (left-mult) parametrization of the same fiber ---")

def fiber_point(alpha, p):
    # beta = -p * alpha  (p: unit imaginary octonion = target point on S^8, real part 0)
    beta = -mul(p, alpha)
    return alpha, beta

p = -e(5)  # target point on S^8 found above: H_O(alpha,beta) = (-e5, 0)

for label, a in [("alpha = e1/sqrt2", alpha), ("a2 = (e1+e2)/2", (e(1)+e(2))/2)]:
    al, be = fiber_point(a, p)
    img, real = H_O(al, be)
    print(f"{label:20s} -> H_O = {img}, {real}")

print()
print("--- general point (c, r) on S^8: |alpha|^2 = (1+r)/2, beta* = alpha* c / (1+r) ---")

rng = np.random.default_rng(0)
for _ in range(3):
    c = rng.normal(size=8)
    r = rng.uniform(-0.9, 0.9)
    c *= np.sqrt(1 - r**2) / np.linalg.norm(c)  # |c|^2 + r^2 = 1
    a = rng.normal(size=8)
    a *= np.sqrt((1 + r) / 2) / np.linalg.norm(a)
    beta = conj(mul(conj(a), c) / (1 + r))
    img, real = H_O(a, beta)
    ok = np.allclose(img, c) and np.isclose(real, r)
    print(f"target r = {r:+.4f} -> H_O matches (c, r)? {ok}")
    assert ok
