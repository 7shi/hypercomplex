"""Checks for 04-extension.md.

Verifies numerically (random variables, numpy) the extension of Hopf fibrations:

1. Quaternionic (S^7 -> S^4):
   - H(α, β) = (2αβ*, |α|² - |β|²) maps S^7 to S^4 (norm is preserved as 1).
   - Fiber is preserved: H(αq, βq) == H(α, β) due to associativity.

2. Octonionic (S^15 -> S^8):
   - H(α, β) = (2αβ*, |α|² - |β|²) maps S^15 to S^8 (norm is preserved as 1).
   - Fiber is generally NOT preserved: H(αq, βq) != H(α, β) due to non-associativity.
"""

import numpy as np
import sys
from pathlib import Path

# Add the parent directory's 'oct' folder to sys.path to import octonion
sys.path.append(str(Path(__file__).parent.parent / 'oct'))
import octonion as oc

rng = np.random.default_rng(42)

def omul(x, y):
    res = np.zeros(8, dtype=float)
    for i in range(8):
        if x[i] != 0:
            res += x[i] * oc.L[i] @ y
    return res

def oconj(x):
    res = -x
    res[0] = x[0]
    return res

def qmul(a, b):
    # a, b are 4-vectors
    a_ext = np.zeros(8); a_ext[:4] = a
    b_ext = np.zeros(8); b_ext[:4] = b
    return omul(a_ext, b_ext)[:4]

def qconj(a):
    return oconj(np.concatenate([a, np.zeros(4)]))[:4]

def check_quaternionic():
    # S^7 -> S^4
    for _ in range(100):
        a = rng.normal(size=4)
        b = rng.normal(size=4)
        norm = np.sqrt(np.sum(a**2) + np.sum(b**2))
        a /= norm
        b /= norm
        
        # H(a, b)
        H_vec = 2 * qmul(a, qconj(b))
        H_real = np.sum(a**2) - np.sum(b**2)
        H_norm = np.sqrt(np.sum(H_vec**2) + H_real**2)
        assert np.isclose(H_norm, 1.0), f"Norm is {H_norm}"
        
        # Fiber check
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        aq = qmul(a, q)
        bq = qmul(b, q)
        
        Hq_vec = 2 * qmul(aq, qconj(bq))
        Hq_real = np.sum(aq**2) - np.sum(bq**2)
        
        assert np.allclose(H_vec, Hq_vec)
        assert np.isclose(H_real, Hq_real)

def check_octonionic():
    # S^15 -> S^8
    non_associative_count = 0
    for _ in range(100):
        a = rng.normal(size=8)
        b = rng.normal(size=8)
        norm = np.sqrt(np.sum(a**2) + np.sum(b**2))
        a /= norm
        b /= norm
        
        # H(a, b)
        H_vec = 2 * omul(a, oconj(b))
        H_real = np.sum(a**2) - np.sum(b**2)
        H_norm = np.sqrt(np.sum(H_vec**2) + H_real**2)
        assert np.isclose(H_norm, 1.0), f"Norm is {H_norm}"
        
        # Fiber check
        q = rng.normal(size=8)
        q /= np.linalg.norm(q)
        aq = omul(a, q)
        bq = omul(b, q)
        
        Hq_vec = 2 * omul(aq, oconj(bq))
        Hq_real = np.sum(aq**2) - np.sum(bq**2)
        
        # The real part is preserved because |aq| = |a||q| = |a|
        assert np.isclose(H_real, Hq_real)
        
        if not np.allclose(H_vec, Hq_vec):
            non_associative_count += 1
            
    assert non_associative_count > 0, "Expected non-associativity to break the fiber"

if __name__ == "__main__":
    check_quaternionic()
    check_octonionic()
    print("check-04-extension: all checks passed")
