"""Checks for 02-nonion.md.

1. 2x2: the example u, v (theta = sqrt(-1)) satisfies u^2 = v^2 = -I,
   uv = -vu, det(zI + yv + xu) = z^2 + y^2 + x^2, and the 4 matrices obey the
   quaternion table; the matrix units give i = -mu + nu etc. (sign vs v).
   Also det(zI+yu+xv) for the swapped order, vu in the remark.
2. 3x3: u, v from the article satisfy u^3 = v^3 = I, vu = rho uv (and not
   rho^2 uv), det(zI + yu + xv) = z^3 + y^3 + x^3 (random samples), and the
   9 basis matrices u^a v^b equal the matrices listed in the article.
3. Swapping u and v keeps the determinant but turns rho into rho^2.
4. rho^2 uv is the shift matrix; the 27 elements rho^k u^a v^b form a group.
"""
import numpy as np

rng = np.random.default_rng(0)
th = 1j
I2 = np.eye(2)
v2 = np.array([[0, 1], [-1, 0]], dtype=complex)
u2 = np.array([[0, th], [th, 0]])
w2 = np.array([[-th, 0], [0, th]])
assert np.allclose(u2 @ u2, -I2) and np.allclose(v2 @ v2, -I2)
assert np.allclose(u2 @ v2, -v2 @ u2)
for _ in range(10):
    x, y, z = rng.normal(size=3)
    assert np.isclose(np.linalg.det(z * I2 + y * v2 + x * u2), x*x + y*y + z*z)
    assert np.isclose(np.linalg.det(z * I2 + y * u2 + x * v2), x*x + y*y + z*z)
assert np.allclose(u2 @ v2, w2) or np.allclose(u2 @ v2, -w2)
print("uv =", "w" if np.allclose(u2 @ v2, w2) else "-w", "(w = diag(-th, th))")
assert np.allclose(v2 @ u2, np.diag([th, -th]))
# matrix units: -mu + nu, theta(mu+nu), -theta lam + theta tau
lam, mu, nu, tau = (np.array(m, dtype=complex) for m in
    ([[1, 0], [0, 0]], [[0, 1], [0, 0]], [[0, 0], [1, 0]], [[0, 0], [0, 1]]))
print("-mu+nu == v ?", np.allclose(-mu + nu, v2), "; == -v ?", np.allclose(-mu + nu, -v2))
print("theta(mu+nu) == u ?", np.allclose(th * (mu + nu), u2))
print("-theta lam + theta tau == w ?", np.allclose(-th * lam + th * tau, w2))

r = np.exp(2j * np.pi / 3)
r2 = r * r
u = np.array([[0, 0, 1], [r, 0, 0], [0, r2, 0]])
v = np.array([[0, 0, 1], [r2, 0, 0], [0, r, 0]])
I3 = np.eye(3)
mp = np.linalg.matrix_power
assert np.allclose(mp(u, 3), I3) and np.allclose(mp(v, 3), I3)
assert np.allclose(v @ u, r * u @ v) and not np.allclose(v @ u, r2 * u @ v)
for _ in range(10):
    x, y, z = rng.normal(size=3)
    assert np.isclose(np.linalg.det(z * I3 + y * u + x * v), x**3 + y**3 + z**3)
    # swapped order has the same determinant
    assert np.isclose(np.linalg.det(z * I3 + y * v + x * u), x**3 + y**3 + z**3)
# swapping u and v: u v = rho^2 v u? i.e. the relation with roles exchanged has rho^2
assert np.allclose(u @ v, r2 * v @ u)

basis = {
 "v":   [[0, 0, 1], [r2, 0, 0], [0, r, 0]],
 "v2":  [[0, r, 0], [0, 0, r2], [1, 0, 0]],
 "uv":  [[0, r, 0], [0, 0, r], [r, 0, 0]],
 "uv2": [[1, 0, 0], [0, r2, 0], [0, 0, r]],
 "u2":  [[0, r2, 0], [0, 0, r], [1, 0, 0]],
 "u2v": [[r, 0, 0], [0, r2, 0], [0, 0, 1]],
 "u2v2":[[0, 0, r], [r, 0, 0], [0, r, 0]],
}
calc = {"v": v, "v2": v @ v, "uv": u @ v, "uv2": u @ v @ v, "u2": u @ u,
        "u2v": u @ u @ v, "u2v2": u @ u @ v @ v}
for k in basis:
    print(k, "matches article:", np.allclose(basis[k], calc[k]))
assert np.allclose(r2 * u @ v, [[0, 1, 0], [0, 0, 1], [1, 0, 0]])

els = []
for k in range(3):
    for a in range(3):
        for b in range(3):
            els.append(r**k * mp(u, a) @ mp(v, b))
def has(m): return any(np.allclose(m, e) for e in els)
assert all(has(p @ q) for p in els for q in els)
uniq = []
for e in els:
    if not any(np.allclose(e, f) for f in uniq): uniq.append(e)
print("group order:", len(uniq))
assert len(uniq) == 27
# counterexample: det = z^3+y^3+x^3 but neither commutation relation
U = np.diag([1, r, r2])
V = 2 ** (-1 / 3) * np.array([[0, 1, r2], [1, 0, 1], [1, r, 0]])
for _ in range(10):
    x, y, z = rng.normal(size=3)
    assert np.isclose(np.linalg.det(z * I3 + y * U + x * V), x**3 + y**3 + z**3)
assert not np.allclose(V @ U, r * U @ V) and not np.allclose(V @ U, r2 * U @ V)
# uv = rho S, vu = rho^2 S
S = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
assert np.allclose(u @ v, r * S) and np.allclose(v @ u, r2 * S)
# i = -v, j = u, k = uv: ij = k
assert np.allclose((-v2) @ u2, u2 @ v2) and np.allclose(v2 @ u2, -(u2 @ v2))
# u = I, v = 0 anticommute although b = 1
print("ok")
