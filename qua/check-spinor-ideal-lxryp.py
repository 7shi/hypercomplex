"""Checks for spinor-ideal.md: L_u R_v P -> uv (product extraction via the
left-ideal projection), generalizing L_x P -> x to L_x R_y P -> xy.
"""

import random

import numpy as np

random.seed(0)

BASIS4 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
NAMES = ["1", "i", "j", "k"]


def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return (w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
            w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2)


def Lmat(u):
    return np.array([qmul(u, b) for b in BASIS4]).T


def Rmat(u):
    return np.array([qmul(b, u) for b in BASIS4]).T


P = np.diag([1.0, 0.0, 0.0, 0.0])

# --- check 0: R_y P = L_y P (both send the identity column to y) -----------

ok = True
for v in BASIS4:
    ok &= np.allclose(Rmat(v) @ P, Lmat(v) @ P)
print("R_y P == L_y P for all basis y:", ok)

# --- check 1: for all 16 basis pairs u,v, L_u R_v P -> uv, other columns 0 --

ok = True
for iu, u in enumerate(BASIS4):
    for iv, v in enumerate(BASIS4):
        S = Lmat(u) @ Rmat(v) @ P
        target_col = np.array(qmul(u, v))
        match = np.allclose(S[:, 0], target_col) and np.allclose(S[:, 1:], 0)
        ok &= match
        if not match:
            print(f"MISMATCH u={NAMES[iu]} v={NAMES[iv]}")

print("all 16 basis pairs L_u R_v P -> uv:", ok)

# example cited in the article: L_i R_j P -> ij = k
S = Lmat(BASIS4[1]) @ Rmat(BASIS4[2]) @ P
print("L_i R_j P -> ij:", np.allclose(S[:, 0], qmul(BASIS4[1], BASIS4[2])),
      "  ij =", qmul(BASIS4[1], BASIS4[2]))

# --- check 2: linearity, L_x R_y P -> xy for random quaternions x, y -------

rng = np.random.default_rng(0)
ok = True
for _ in range(20):
    x = rng.standard_normal(4)
    y = rng.standard_normal(4)
    Lx = sum(x[i] * Lmat(BASIS4[i]) for i in range(4))
    Ry = sum(y[i] * Rmat(BASIS4[i]) for i in range(4))
    S = Lx @ Ry @ P
    xy = qmul(tuple(x), tuple(y))
    ok &= np.allclose(S[:, 0], xy) and np.allclose(S[:, 1:], 0)

print("random x,y: L_x R_y P -> xy (20 trials):", ok)

# --- check 3: the reduction L_x R_y P = L_{xy} P (associativity) -----------

ok = True
for _ in range(20):
    x = rng.standard_normal(4)
    y = rng.standard_normal(4)
    Lx = sum(x[i] * Lmat(BASIS4[i]) for i in range(4))
    Ry = sum(y[i] * Rmat(BASIS4[i]) for i in range(4))
    xy = qmul(tuple(x), tuple(y))
    Lxy = sum(xy[i] * Lmat(BASIS4[i]) for i in range(4))
    ok &= np.allclose(Lx @ Ry @ P, Lxy @ P)

print("L_x R_y P == L_{xy} P (reduction via associativity):", ok)
