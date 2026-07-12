"""Checks for 03-bloch-density.md.

Verifies numerically (random parameters, numpy):

1. Pauli products (σxσy = iσz etc., anticommutation), quaternion
   correspondence (-iσx, -iσy, -iσz follow the i, j, k rules), traces
   tr σi = 0, tr(σiσj) = 2δij, and coefficient extraction: a random
   Hermitian A equals ½{tr(A) I + Σ tr(Aσi) σi}.
2. Pure state ρ = ωω†: tr ρ = 1, ρ² = ρ, Bloch vector r from
   x = α*β+β*α, y = -i(α*β-β*α), z = α*α-β*β matches tr(ρσi), |r| = 1,
   and ρ is invariant under ω → ω e^{ic}.
3. (r·σ)² = |r|² I for arbitrary real r, and
   ρ² = ¼{(1+|r|²) I + 2 r·σ}.
4. Mixed state ρ = Σ pi ωiωi†: Bloch vector r = Σ pi ri, |r| < 1,
   purity tr ρ² = ½(1+|r|²).
5. Eigendecomposition: eigenvalues ½(1±|r|), eigenvectors are pure
   states with Bloch vectors ±r/|r| (antipodal), and
   ρ = λ+ρ+ + λ-ρ- reconstructs ρ.
6. Example: equal mixture of ω0 = (1,0)ᵀ and ω1 = (0,1)ᵀ gives ½I
   (r = 0), also obtained from the ±x antipodal pair; the superposition
   ω+ = (1,1)ᵀ/√2 gives ½(I+σx) (r = (1,0,0)) — same diagonal,
   different off-diagonal.
"""

import numpy as np

rng = np.random.default_rng(42)

I2 = np.eye(2, dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
paulis = (sigma_x, sigma_y, sigma_z)


def bloch(rho):
    return np.array([np.trace(rho @ s).real for s in paulis])


# 1. Pauli products, quaternion correspondence, traces, extraction
for s in paulis:
    assert np.allclose(s @ s, I2)
    assert np.isclose(np.trace(s), 0)
assert np.allclose(sigma_x @ sigma_y, 1j * sigma_z)
assert np.allclose(sigma_y @ sigma_z, 1j * sigma_x)
assert np.allclose(sigma_z @ sigma_x, 1j * sigma_y)
assert np.allclose(sigma_y @ sigma_x, -1j * sigma_z)
assert np.allclose(sigma_z @ sigma_y, -1j * sigma_x)
assert np.allclose(sigma_x @ sigma_z, -1j * sigma_y)
qi, qj, qk = -1j * sigma_x, -1j * sigma_y, -1j * sigma_z
assert np.allclose(qi @ qi, -I2)
assert np.allclose(qj @ qj, -I2)
assert np.allclose(qk @ qk, -I2)
assert np.allclose(qi @ qj, qk)
assert np.allclose(qj @ qk, qi)
assert np.allclose(qk @ qi, qj)
for i, si in enumerate(paulis):
    for j, sj in enumerate(paulis):
        assert np.isclose(np.trace(si @ sj), 2 * (i == j))

for _ in range(100):
    m = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    a = m + m.conj().T
    rebuilt = (np.trace(a) * I2
               + sum(np.trace(a @ s) * s for s in paulis)) / 2
    assert np.allclose(a, rebuilt)

    # 2. pure state
    v = rng.normal(size=4)
    v /= np.linalg.norm(v)
    alpha, beta = v[0] + 1j * v[1], v[2] + 1j * v[3]
    w = np.array([alpha, beta])
    rho = np.outer(w, w.conj())
    assert np.isclose(np.trace(rho), 1)
    assert np.allclose(rho @ rho, rho)
    ab = np.conj(alpha) * beta
    r = np.array([2 * ab.real, 2 * ab.imag, abs(alpha)**2 - abs(beta)**2])
    assert np.allclose(bloch(rho), r)
    assert np.isclose(np.linalg.norm(r), 1)
    assert np.allclose(rho, (I2 + sum(r[i] * paulis[i] for i in range(3))) / 2)
    phase = np.exp(1j * rng.uniform(0, 2 * np.pi))
    assert np.allclose(np.outer(w * phase, (w * phase).conj()), rho)

    # 3. (r·σ)² = |r|² I and ρ² formula for arbitrary r
    r3 = rng.normal(size=3)
    rs = sum(r3[i] * paulis[i] for i in range(3))
    assert np.allclose(rs @ rs, (r3 @ r3) * I2)
    rho3 = (I2 + rs) / 2
    assert np.allclose(rho3 @ rho3, ((1 + r3 @ r3) * I2 + 2 * rs) / 4)

    # 4. mixed state
    n = 3
    p = rng.dirichlet(np.ones(n))
    vs = rng.normal(size=(n, 4))
    ws = (vs[:, 0] + 1j * vs[:, 1], vs[:, 2] + 1j * vs[:, 3])
    ws = np.stack(ws, axis=1)
    ws /= np.linalg.norm(ws, axis=1, keepdims=True)
    rhos = [np.outer(wi, wi.conj()) for wi in ws]
    rho_mix = sum(p[i] * rhos[i] for i in range(n))
    r_mix = bloch(rho_mix)
    assert np.allclose(r_mix, sum(p[i] * bloch(rhos[i]) for i in range(n)))
    assert np.linalg.norm(r_mix) < 1
    assert np.isclose(np.trace(rho_mix @ rho_mix).real,
                      (1 + r_mix @ r_mix) / 2)

    # 5. eigendecomposition
    lam, vec = np.linalg.eigh(rho_mix)
    rn = np.linalg.norm(r_mix)
    assert np.allclose(lam, [(1 - rn) / 2, (1 + rn) / 2])
    assert np.allclose(bloch(np.outer(vec[:, 1], vec[:, 1].conj())), r_mix / rn)
    assert np.allclose(bloch(np.outer(vec[:, 0], vec[:, 0].conj())), -r_mix / rn)
    assert np.allclose(sum(lam[i] * np.outer(vec[:, i], vec[:, i].conj())
                           for i in range(2)), rho_mix)

# 6. mixture vs superposition example
w0 = np.array([1, 0], dtype=complex)
w1 = np.array([0, 1], dtype=complex)
rho_mix = (np.outer(w0, w0.conj()) + np.outer(w1, w1.conj())) / 2
assert np.allclose(rho_mix, I2 / 2)
assert np.allclose(bloch(rho_mix), [0, 0, 0])
wp = np.array([1, 1], dtype=complex) / np.sqrt(2)
wm = np.array([1, -1], dtype=complex) / np.sqrt(2)
assert np.allclose((np.outer(wp, wp.conj()) + np.outer(wm, wm.conj())) / 2,
                   rho_mix)
rho_sup = np.outer(wp, wp.conj())
assert np.allclose(rho_sup, (I2 + sigma_x) / 2)
assert np.allclose(bloch(rho_sup), [1, 0, 0])
assert np.allclose(np.diag(rho_mix), np.diag(rho_sup))
assert not np.allclose(rho_mix, rho_sup)

print("check-03-bloch-density: all checks passed")
