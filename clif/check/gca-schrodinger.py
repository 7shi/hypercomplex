"""Checks for gca-schrodinger.md.

The article avoids matrices, but the operations X (shift, between
components) and Z (clock, within components) on an n-level system are
checked here as matrices.

1. X^n = Z^n = 1, ZX = w XZ and Z^l X^k = w^(kl) X^k Z^l
   with w = exp(2 pi i / n).
2. n = 2: X swaps the components, Z flips the sign of the second,
   (XZ)^2 = (iZ)^2 = (iX)^2 = -1,
   cos t + sin t XZ is the rotation between components and
   cos t + i sin t Z multiplies the components by e^(it), e^(-it).
3. n = 3: the explicit actions of ZX and XZ in the article, and
   xxz + xzx + zxx = 0, xzz + zxz + zzx = 0.
4. (aX + bZ)^n = a^n + b^n.
5. X^a Z^b (0 <= a, b < n) span all operations (n = 3).
6. Discrete Fourier transform phi_k = n^(-1/2) sum_j w^(-jk) psi_j:
   unitary, X w_k = w^(-k) w_k, Z w_k = w_(k+1); on phi, Z shifts and
   X multiplies the k-th component by w^(-k).  Examples n = 2 and n = 4.
7. Toy model hbar = 2 pi / n, x_j = hbar j, points numbered on both
   sides of the origin; for a Gaussian of width sqrt(hbar / 2):
   X ~ 1 - iP, Z ~ 1 + iQ, (ZX - XZ) psi ~ i hbar psi,
   (2 - X - X^-1) psi ~ P^2 psi  (errors small compared with hbar).
8. P = -i hbar d/dx multiplies e^(ikx) by hbar k; wavelength 2 pi / k;
   (2 pi)^2 / n = 2 pi hbar.
9. Multiplication by an arbitrary function at each point is a linear
   combination of 1, Z, ..., Z^(n-1).
10. Uncertainty: sigma_x sigma_p = hbar / 2 for Gaussians (table in the
    article), larger for another shape.
11. H = (2 - X - X^-1)/2 + (2 - Z - Z^-1)/2 has lowest levels close to
    hbar (k + 1/2) (values in the article for n = 256).
"""

import numpy as np

rng = np.random.default_rng(0)


def ops(n):
    w = np.exp(2j * np.pi / n)
    Z = np.diag(w ** np.arange(n))
    X = np.roll(np.eye(n), 1, axis=0)  # (X psi)_j = psi_{j-1}
    return w, X, Z


def mpow(A, k):
    return np.linalg.matrix_power(A, k)


# 1, 4
for n in range(2, 9):
    w, X, Z = ops(n)
    I = np.eye(n)
    psi = np.arange(n) + 1.0
    assert np.allclose((X @ psi)[1:], psi[:-1]) and (X @ psi)[0] == psi[-1]
    assert np.allclose(mpow(X, n), I) and np.allclose(mpow(Z, n), I)
    assert np.allclose(Z @ X, w * X @ Z)
    for k in range(n):
        for l in range(n):
            assert np.allclose(
                mpow(Z, l) @ mpow(X, k), w ** (k * l) * mpow(X, k) @ mpow(Z, l)
            )
    a, b = rng.normal(size=2)
    assert np.allclose(mpow(a * X + b * Z, n), (a**n + b**n) * I)
print("1, 4: relations and (aX+bZ)^n = a^n+b^n  ok")

# 2
w, X, Z = ops(2)
I = np.eye(2)
p = np.array([3.0, 5.0])
assert np.allclose(X @ p, [5, 3]) and np.allclose(Z @ p, [3, -5])
assert np.allclose(Z @ X @ p, [5, -3]) and np.allclose(X @ Z @ p, [-5, 3])
for J in (X @ Z, 1j * Z, 1j * X):
    assert np.allclose(J @ J, -I)
t = 0.7
q = rng.normal(size=2) + 1j * rng.normal(size=2)
c, s = np.cos(t), np.sin(t)
assert np.allclose(
    (c * I + s * X @ Z) @ q, [q[0] * c - q[1] * s, q[0] * s + q[1] * c]
)
assert np.allclose(
    (c * I + 1j * s * Z) @ q, [np.exp(1j * t) * q[0], np.exp(-1j * t) * q[1]]
)
print("2: 2-level system  ok")

# 3, 5
w, X, Z = ops(3)
q = rng.normal(size=3) + 1j * rng.normal(size=3)
assert np.allclose(X @ q, [q[2], q[0], q[1]])
assert np.allclose(Z @ q, [q[0], w * q[1], w**2 * q[2]])
assert np.allclose(Z @ X @ q, [q[2], w * q[0], w**2 * q[1]])
assert np.allclose(X @ Z @ q, [w**2 * q[2], q[0], w * q[1]])
assert np.allclose(X @ X @ Z + X @ Z @ X + Z @ X @ X, 0)
assert np.allclose(X @ Z @ Z + Z @ X @ Z + Z @ Z @ X, 0)
basis = [(mpow(X, a) @ mpow(Z, b)).ravel() for a in range(3) for b in range(3)]
assert np.linalg.matrix_rank(np.array(basis)) == 9
print("3, 5: 3-level system  ok")

# 6: discrete Fourier transform phi_k = n^(-1/2) sum_j w^(-jk) psi_j
for n in (2, 3, 4, 7):
    w, X, Z = ops(n)
    F = np.array([[w ** (-j * k) for j in range(n)] for k in range(n)]) / np.sqrt(n)
    Fi = F.conj().T
    assert np.allclose(F @ Fi, np.eye(n))  # inverse formula, sums preserved
    for k in range(n):
        wk = w ** (np.arange(n) * k)
        assert np.allclose(X @ wk, w ** (-k) * wk)
        assert np.allclose(Z @ wk, w ** (np.arange(n) * (k + 1)))
    assert np.allclose(F @ Z @ Fi, X)  # Z shifts phi
    assert np.allclose(F @ X @ Fi, np.linalg.inv(Z))  # X rotates phi backwards
w, X, Z = ops(2)
F = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
q = rng.normal(size=2)
assert np.allclose(F @ q, [(q[0] + q[1]) / np.sqrt(2), (q[0] - q[1]) / np.sqrt(2)])
w, X, Z = ops(4)
assert np.isclose(w, 1j)
F = np.array([[w ** (-j * k) for j in range(4)] for k in range(4)]) / 2
assert np.allclose(F @ (np.array([1, 1j, -1, -1j]) / 2), [0, 1, 0, 0])
assert np.allclose(F @ np.array([1, 0, 0, 0]), [0.5] * 4)
print("6: discrete Fourier transform  ok")

# 7: toy model, hbar = 2 pi / n, x_j = hbar j in (-pi, pi]
print("7: n, |X-(1-iP)|, |Z-(1+iQ)|, |(ZX-XZ)-i hbar|, |(2-X-X^-1)-P^2|/hbar")
prev = None
for n in (64, 256, 1024, 4096):
    w, X, Z = ops(n)
    hbar = 2 * np.pi / n
    assert np.isclose(w, np.exp(1j * hbar))
    j = (np.arange(n) + n // 2) % n - n // 2
    x = hbar * j
    assert np.allclose(np.diag(Z), np.exp(1j * x))
    s = np.sqrt(hbar / 2)  # concentrated in both x and p
    psi = np.exp(-(x**2) / (4 * s * s))
    d1 = -x / (2 * s * s) * psi
    d2 = (x**2 / (4 * s**4) - 1 / (2 * s * s)) * psi
    nrm = np.linalg.norm(psi)
    e1 = np.linalg.norm(X @ psi - (psi - hbar * d1)) / nrm
    e2 = np.linalg.norm(Z @ psi - (psi + 1j * x * psi)) / nrm
    e3 = np.linalg.norm((Z @ X - X @ Z) @ psi - 1j * hbar * psi) / nrm / hbar
    Xi = X.T
    e4 = np.linalg.norm((2 * psi - X @ psi - Xi @ psi) + hbar**2 * d2) / nrm / hbar
    print(f"  {n:5d} {e1:.2e} {e2:.2e} {e3:.2e} {e4:.2e}")
    if prev is not None:
        assert all(a < b for a, b in zip((e1, e2, e3, e4), prev))
    prev = (e1, e2, e3, e4)

# 8: P e^(ikx) = hbar k e^(ikx) (finite difference), p = hbar k = h / lambda
hbar, k, x0, d = 2 * np.pi / 256, 5, 0.3, 1e-6
f = lambda x: np.exp(1j * k * x)
Pf = -1j * hbar * (f(x0 + d) - f(x0 - d)) / (2 * d)
assert np.isclose(Pf, hbar * k * f(x0))
assert np.isclose(f(x0 + 2 * np.pi / k), f(x0))
assert np.isclose(hbar * k, 2 * np.pi * hbar / (2 * np.pi / k))
n = 256
assert np.isclose((2 * np.pi) ** 2 / n, 2 * np.pi * (2 * np.pi / n))
print("8: wave number and momentum  ok")

# 10: uncertainty, sigma_x sigma_p = hbar / 2 for Gaussians (n = 256)
n = 256
w, X, Z = ops(n)
hbar = 2 * np.pi / n
j = (np.arange(n) + n // 2) % n - n // 2
x = p = hbar * j
F = np.exp(-2j * np.pi * np.outer(np.arange(n), np.arange(n)) / n) / np.sqrt(n)
print(f"10: hbar = {hbar:.4f}; sigma_x, sigma_p, product / hbar")
for s in (0.05, np.sqrt(hbar / 2), 0.6):
    psi = np.exp(-(x**2) / (4 * s * s))
    psi /= np.linalg.norm(psi)
    phi = F @ psi
    sx = np.sqrt(np.sum(x**2 * abs(psi) ** 2))
    sp = np.sqrt(np.sum(p**2 * abs(phi) ** 2))
    print(f"  {sx:.3f} {sp:.3f} {sx * sp / hbar:.3f}")
    assert np.isclose(sx * sp / hbar, 0.5, atol=1e-3)
# a non-Gaussian shape gives a larger product
psi = np.where(abs(x) < 0.3, 1.0, 0.0)
psi /= np.linalg.norm(psi)
phi = F @ psi
assert np.sqrt(np.sum(x**2 * abs(psi) ** 2) * np.sum(p**2 * abs(phi) ** 2)) > hbar / 2

# 11: H = (2-X-X^-1)/2 + (2-Z-Z^-1)/2 has levels close to hbar (k + 1/2)
print("11: lowest levels of H in units of hbar")
for n in (64, 256, 1024):
    w, X, Z = ops(n)
    hbar = 2 * np.pi / n
    H = (4 * np.eye(n) - X - X.T - Z - Z.conj().T) / 2
    E = np.linalg.eigvalsh(H)[:4] / hbar
    print(f"  {n:5d}", np.round(E, 4))
assert np.allclose(E, [0.5, 1.5, 2.5, 3.5], atol=0.02)

# 9
n = 7
w, X, Z = ops(n)
v = rng.normal(size=n)
cb = [sum(v[j] * w ** (-j * b) for j in range(n)) / n for b in range(n)]
assert np.allclose(sum(cb[b] * mpow(Z, b) for b in range(n)), np.diag(v))
print("9: pointwise multiplication is a polynomial in Z  ok")
