"""Checks for 07-real-qubit.md (two-qubit part).

Verifies numerically (random variables, numpy):

1. H, X, Z, CNOT, CZ, CCX are real matrices, the four Bell states are real
   vectors, and the determinant criterion αδ - βγ = 0 for separability works
   verbatim on real amplitudes.
2. CHSH with real states and real observables reaches the Tsirelson bound
   |S| = 2√2; the optimal Bloch directions 0°, 45°, 90°, 135° are the
   doubles of the state-space angles 0°, 22.5°, 45°, 67.5°.
3. Real amplitudes lose no computational power: the embedding
   ψ_R + iψ_I ↦ (ψ_R, ψ_I) turns a complex unitary A + iB into the real
   orthogonal matrix I⊗A + J⊗B (J = [[0,-1],[1,0]]) on one extra rebit,
   and the measurement probabilities agree.
4. Local tomography fails: real symmetric 4x4 matrices form a 10-dimensional
   space but {I, X, Z}⊗{I, X, Z} spans only 9; σ_y⊗σ_y is real symmetric and
   fills the gap. The states ρ± = (I⊗I ± σ_y⊗σ_y)/4 are distinct, positive
   semidefinite, unit trace, and give identical ⟨A⊗B⟩ for every pair of real
   observables.
"""

import numpy as np

rng = np.random.default_rng(42)

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=float)
Z = np.array([[1, 0], [0, -1]], dtype=float)
sy = np.array([[0, -1j], [1j, 0]])
H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
J = np.array([[0, -1], [1, 0]], dtype=float)

CNOT = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=float)
CZ = np.diag([1, 1, 1, -1]).astype(float)
CCX = np.eye(8)
CCX[[6, 7]] = CCX[[7, 6]]

s2 = np.sqrt(2)
bell = {
    "Phi+": np.array([1, 0, 0, 1]) / s2,
    "Phi-": np.array([1, 0, 0, -1]) / s2,
    "Psi+": np.array([0, 1, 1, 0]) / s2,
    "Psi-": np.array([0, 1, -1, 0]) / s2,
}

def real_state(n):
    v = rng.normal(size=n)
    return v / np.linalg.norm(v)

def sym_obs():
    # a random real symmetric 2x2 observable
    a = rng.normal(size=(2, 2))
    return a + a.T

def check_real_gates():
    for g in (H, X, Z, CNOT, CZ, CCX):
        assert np.allclose(g.imag if np.iscomplexobj(g) else 0, 0)
        assert np.allclose(g @ g.T, np.eye(g.shape[0]))
    # the Bell circuit stays inside the real amplitudes
    ket00 = np.array([1, 0, 0, 0], dtype=float)
    assert np.allclose(CNOT @ np.kron(H, I2) @ ket00, bell["Phi+"])
    for v in bell.values():
        assert np.allclose(v.imag if np.iscomplexobj(v) else 0, 0)
        a, b, c, d = v
        assert not np.isclose(a * d - b * c, 0)     # all four are entangled
    # determinant criterion on real amplitudes
    for _ in range(100):
        u, w = real_state(2), real_state(2)
        a, b, c, d = np.kron(u, w)
        assert np.isclose(a * d - b * c, 0)         # product states: det = 0
    for _ in range(100):
        a, b, c, d = real_state(4)
        sep = np.isclose(a * d - b * c, 0)
        # rank of the 2x2 coefficient matrix decides separability
        rank = np.linalg.matrix_rank(np.array([[a, b], [c, d]]), tol=1e-9)
        assert sep == (rank == 1)

def check_chsh():
    def obs(angle):
        # a real observable with Bloch axis at `angle` in the (z, x) plane
        return np.cos(angle) * Z + np.sin(angle) * X

    psi = bell["Psi-"]                       # the singlet, a real vector
    state_angles = np.radians([0, 45, 22.5, 67.5])
    a, ap, b, bp = [obs(2 * t) for t in state_angles]   # Bloch angles double

    def E(u, v):
        return psi @ np.kron(u, v) @ psi

    S = E(a, b) - E(a, bp) + E(ap, b) + E(ap, bp)
    assert np.isclose(abs(S), 2 * s2)
    # E is the negative cosine of the Bloch angle difference
    for _ in range(50):
        t1, t2 = rng.uniform(0, 2 * np.pi, 2)
        assert np.isclose(E(obs(t1), obs(t2)), -np.cos(t1 - t2))
    # every observable used lies in the y = 0 plane (real symmetric)
    for g in (a, ap, b, bp):
        assert np.allclose(g, g.T) and np.allclose(g.imag if np.iscomplexobj(g) else 0, 0)

def check_real_encoding():
    for _ in range(50):
        n = 4
        psi = rng.normal(size=n) + 1j * rng.normal(size=n)
        psi /= np.linalg.norm(psi)
        q, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
        U = q
        A, B = U.real, U.imag
        Ur = np.kron(I2, A) + np.kron(J, B)
        assert np.allclose(Ur @ Ur.T, np.eye(2 * n))          # real orthogonal
        enc = np.concatenate([psi.real, psi.imag])
        assert np.isclose(np.linalg.norm(enc), 1)
        # the encoding intertwines the two actions
        assert np.allclose(Ur @ enc, np.concatenate([(U @ psi).real, (U @ psi).imag]))
        # measuring the data register reproduces the complex probabilities
        out = Ur @ enc
        probs = out[:n] ** 2 + out[n:] ** 2
        assert np.allclose(probs, np.abs(U @ psi) ** 2)

def check_local_tomography():
    local = [np.kron(a, b) for a in (I2, X, Z) for b in (I2, X, Z)]
    yy = np.real(np.kron(sy, sy))
    assert np.allclose(np.kron(sy, sy).imag, 0)         # σ_y⊗σ_y is real
    assert np.allclose(yy, yy.T)                        # and symmetric
    basis = np.array([m.flatten() for m in local])
    assert np.linalg.matrix_rank(basis, tol=1e-9) == 9
    with_yy = np.vstack([basis, yy.flatten()])
    assert np.linalg.matrix_rank(with_yy, tol=1e-9) == 10
    assert 4 * 5 // 2 == 10                             # dim of real symmetric 4x4

    rho_p = (np.eye(4) + yy) / 4
    rho_m = (np.eye(4) - yy) / 4
    for rho in (rho_p, rho_m):
        assert np.isclose(np.trace(rho), 1)
        assert np.allclose(rho, rho.T)
        assert np.all(np.linalg.eigvalsh(rho) > -1e-12)
    assert not np.allclose(rho_p, rho_m)
    # indistinguishable by any pair of real observables
    for _ in range(200):
        A, B = sym_obs(), sym_obs()
        M = np.kron(A, B)
        assert np.isclose(np.trace(rho_p @ M), np.trace(rho_m @ M))
    # the reason: tr(Aσ_y) = 0 for every real symmetric A
    for _ in range(100):
        assert np.isclose(np.trace(sym_obs() @ sy), 0)
    # σ_y⊗σ_y itself separates them, but it is not a product of local
    # real observables, so it is unavailable to local measurement
    assert not np.isclose(np.trace(rho_p @ yy), np.trace(rho_m @ yy))
    assert not any(np.allclose(yy, c * m) for m in local for c in (1, -1))

if __name__ == "__main__":
    check_real_gates()
    check_chsh()
    check_real_encoding()
    check_local_tomography()
    print("check-07-two-rebit: all checks passed")
