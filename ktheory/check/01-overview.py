"""Checks for ktheory/01-overview.md (overview of K-theory).

1. Moebius band: two copies are glued by -I_2, and R(t pi) (rotation) is a path in
   GL(2, R) from I_2 to -I_2, so M + M is trivial. A single copy glued by -1 cannot be
   deformed to +1 even after adding trivial bundles: diag(-1, 1, ..., 1) has det -1.
2. Complex line bundles on S^1: e^{i pi t} is a path in C^x from 1 to -1, so the
   complex Moebius bundle is trivial.
3. Gluing functions by unit numbers: for R, C, H, O, left multiplication by x is
   invertible for x != 0 since |x y| = |x| |y|; the linear combination sum x_i e_i of the
   Cl_{0,k-1} generators (imaginary units) plus x_0 has (x_0 + sum x_i e_i) times its
   conjugate equal to |x|^2, and (sum x_i e_i)^2 = -|x|^2 for the imaginary part.
"""

import numpy as np

from common.octonion import L as OL

rng = np.random.default_rng(2)


def rot(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s], [s, c]])


print("1. Moebius band")
ts = np.linspace(0, 1, 101)
dets = [np.linalg.det(rot(t * np.pi)) for t in ts]
assert np.allclose(rot(0), np.eye(2)) and np.allclose(rot(np.pi), -np.eye(2))
assert min(dets) > 0.999
print("  R(t pi): I_2 -> -I_2 within GL(2,R), det = 1 along the path: OK")
for n in range(1, 6):
    d = np.diag([-1] + [1] * (n - 1))
    assert np.linalg.det(d) < 0
print("  det diag(-1,1,...,1) = -1 for n = 1..5 (not in the identity component): OK")

print("2. complex line bundles on S^1")
zs = np.exp(1j * np.pi * ts)
assert np.isclose(zs[0], 1) and np.isclose(zs[-1], -1) and np.allclose(abs(zs), 1)
print("  e^{i pi t}: 1 -> -1 within C^x: OK")

print("3. gluing by unit numbers")
dims = {1: "R", 2: "C", 4: "H", 8: "O"}
for n, name in dims.items():
    Ls = [m[:n, :n] for m in OL[:n]]  # left multiplications by 1, e_1, ..., e_{n-1}
    for _ in range(50):
        x = rng.normal(size=n)
        Lx = sum(xi * m for xi, m in zip(x, Ls))
        r2 = x @ x
        assert np.allclose(Lx.T @ Lx, r2 * np.eye(n))  # |x y| = |x| |y|
        assert np.isclose(abs(np.linalg.det(Lx)), r2 ** (n / 2))
        Im = sum(xi * m for xi, m in zip(x[1:], Ls[1:])) if n > 1 else np.zeros((1, 1))
        assert np.allclose(Im @ Im, -(x[1:] @ x[1:]) * np.eye(n))
    print(f"  {name}: L_x^T L_x = |x|^2, (sum x_i e_i)^2 = -|x|^2 on the imaginary part: OK")
print("All checks passed.")
