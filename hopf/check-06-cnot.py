"""Checks for 06-clifford-gates.md (CNOT / two-qubit part).

Verifies numerically (random variables, numpy):

1. CNOT is unitary and self-inverse, and equals |0><0| ⊗ I + |1><1| ⊗ X.
2. Conjugation table: X⊗I ↦ X⊗X, Y⊗I ↦ Y⊗X, Z⊗I ↦ Z⊗I,
   I⊗X ↦ I⊗X, I⊗Y ↦ Z⊗Y, I⊗Z ↦ Z⊗Z; all 15 nontrivial two-qubit
   Pauli operators map to ± another one (CNOT is Clifford).
3. Local gates U⊗V transform M = [[α,β],[γ,δ]] as U M V^T, so
   C = 2|αδ - βγ| is invariant; CNOT creates entanglement: on a
   product state (a,b)⊗(c,d) the determinant becomes
   αδ - βγ = ab(c² - d²) ≠ 0 generically, so CNOT is not local.
4. CNOT(H⊗I)|00> = Φ+, whose quaternion Hopf image is -j with C = 1
   (helpers from check-05-two-qubit.py).
5. Axis tracking for U = CNOT(H⊗I): U(Z⊗I)U† = X⊗X, U(I⊗Z)U† = Z⊗Z;
   expectation values on Φ+: <X⊗X> = <Z⊗Z> = 1, <σ⊗I> = 0.
6. Pauli expansion ρ = (1/4)Σ c_ab σ_a⊗σ_b with
   c_ab = tr{ρ(σ_a⊗σ_b)} reconstructs ρ; random products of
   H⊗I, I⊗H, S⊗I, I⊗S, CNOT induce signed permutations of the
   15 nontrivial axes, while T⊗I does not.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent / 'oct'))
from octonion import L

rng = np.random.default_rng(42)

def omul(x, y):
    M = sum(x[i] * L[i] for i in range(8))
    return M @ y

def qmul(a, b):
    a_ext = np.zeros(8); a_ext[:4] = a
    b_ext = np.zeros(8); b_ext[:4] = b
    return omul(a_ext, b_ext)[:4]

def qconj(a):
    b = -a.copy()
    b[0] = a[0]
    return b

def c2q(alpha, beta):
    # q = alpha + beta j  (basis 1, i, j, k)
    return np.array([alpha.real, alpha.imag, beta.real, beta.imag])

def rand_state(n):
    v = rng.normal(size=n) + 1j * rng.normal(size=n)
    return v / np.linalg.norm(v)

def hopf_pair(psi):
    # psi = (α, β, γ, δ): q1 = α + βj, q2 = γ + δj
    return c2q(psi[0], psi[1]), c2q(psi[2], psi[3])

def H_H(q1, q2):
    return 2 * qmul(q1, qconj(q2)), q1 @ q1 - q2 @ q2

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], dtype=complex)
pauli4 = [I2, sx, sy, sz]

Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
S = np.diag([1, 1j])
T = np.diag([1, np.exp(1j * np.pi / 4)])
CNOT = np.array([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 1],
                 [0, 0, 1, 0]], dtype=complex)

pauli2q = [np.kron(a, b) for a in pauli4 for b in pauli4]  # 16, [0] = I⊗I

def check_cnot_basics():
    assert np.allclose(CNOT @ CNOT.conj().T, np.eye(4))
    assert np.allclose(CNOT @ CNOT, np.eye(4))
    p0 = np.array([[1, 0], [0, 0]], dtype=complex)
    p1 = np.array([[0, 0], [0, 1]], dtype=complex)
    assert np.allclose(CNOT, np.kron(p0, I2) + np.kron(p1, sx))

def check_conjugation_table():
    table = [
        (np.kron(sx, I2), np.kron(sx, sx)),
        (np.kron(sy, I2), np.kron(sy, sx)),
        (np.kron(sz, I2), np.kron(sz, I2)),
        (np.kron(I2, sx), np.kron(I2, sx)),
        (np.kron(I2, sy), np.kron(sz, sy)),
        (np.kron(I2, sz), np.kron(sz, sz)),
    ]
    for P, img in table:
        assert np.allclose(CNOT @ P @ CNOT.conj().T, img)
    # all 15 nontrivial Paulis map to ± another one
    for P in pauli2q[1:]:
        img = CNOT @ P @ CNOT.conj().T
        hits = [c for Q in pauli2q[1:] for c in (1, -1) if np.allclose(img, c * Q)]
        assert len(hits) == 1

def rand_unitary(n):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(A)
    return Q * (np.diag(R) / np.abs(np.diag(R)))

def check_entanglement_creation():
    # local gates: (U⊗V)Ψ has M' = U M V^T, det M' = det U det V det M
    for _ in range(50):
        psi = rand_state(4)
        M = psi.reshape(2, 2)
        U, V = rand_unitary(2), rand_unitary(2)
        M2 = (np.kron(U, V) @ psi).reshape(2, 2)
        assert np.allclose(M2, U @ M @ V.T)
        det, det2 = np.linalg.det(M), np.linalg.det(M2)
        assert np.isclose(det2, np.linalg.det(U) * np.linalg.det(V) * det)
        assert np.isclose(abs(det2), abs(det))
    for _ in range(50):
        a, b = rand_state(2)
        c, d = rand_state(2)
        psi = CNOT @ np.kron(np.array([a, b]), np.array([c, d]))
        det = psi[0] * psi[3] - psi[1] * psi[2]
        assert np.isclose(det, a * b * (c**2 - d**2))
    # concrete: |+0> becomes Φ+ (maximally entangled)
    plus = np.array([1, 1]) / np.sqrt(2)
    psi = CNOT @ np.kron(plus, np.array([1, 0], dtype=complex))
    assert np.isclose(abs(psi[0] * psi[3] - psi[1] * psi[2]), 0.5)

def check_bell_circuit():
    psi00 = np.zeros(4, dtype=complex); psi00[0] = 1
    psi = CNOT @ np.kron(Hd, I2) @ psi00
    bell = np.array([1, 0, 0, 1]) / np.sqrt(2)
    assert np.allclose(psi, bell)
    # Hopf image -j, C = 1 (05-entanglement.md)
    img, r = H_H(*hopf_pair(psi))
    assert np.allclose(img, [0, 0, -1, 0]) and np.isclose(r, 0)
    assert np.isclose(2 * abs(psi[0] * psi[3] - psi[1] * psi[2]), 1)

def check_axis_tracking():
    U = CNOT @ np.kron(Hd, I2)
    assert np.allclose(U @ np.kron(sz, I2) @ U.conj().T, np.kron(sx, sx))
    assert np.allclose(U @ np.kron(I2, sz) @ U.conj().T, np.kron(sz, sz))
    bell = np.array([1, 0, 0, 1]) / np.sqrt(2)
    assert np.isclose(bell.conj() @ np.kron(sx, sx) @ bell, 1)
    assert np.isclose(bell.conj() @ np.kron(sz, sz) @ bell, 1)
    for s in (sx, sy, sz):
        assert np.isclose(bell.conj() @ np.kron(s, I2) @ bell, 0)

def induced_15(U):
    # M[p,q] = (1/4) tr{P_p U P_q U†} over the 15 nontrivial Paulis
    M = np.empty((15, 15))
    for q, P in enumerate(pauli2q[1:]):
        img = U @ P @ U.conj().T
        for p, Q in enumerate(pauli2q[1:]):
            M[p, q] = np.trace(Q @ img).real / 4
    return M

def is_signed_permutation(M):
    A = np.abs(M) > 0.5
    return (np.allclose(np.abs(M) * A, np.abs(M), atol=1e-8)
            and np.allclose(np.abs(M)[A], 1)
            and (A.sum(axis=0) == 1).all() and (A.sum(axis=1) == 1).all())

def check_pauli_expansion():
    for _ in range(20):
        psi = rand_state(4)
        rho = np.outer(psi, psi.conj())
        c = np.array([np.trace(rho @ P) for P in pauli2q])
        assert np.allclose(c.imag, 0)
        recon = sum(c[i].real * pauli2q[i] for i in range(16)) / 4
        assert np.allclose(recon, rho)
    gens = [np.kron(Hd, I2), np.kron(I2, Hd), np.kron(S, I2), np.kron(I2, S), CNOT]
    for _ in range(20):
        word = [gens[i] for i in rng.integers(0, len(gens), size=6)]
        U = np.linalg.multi_dot(word) if len(word) > 1 else word[0]
        assert is_signed_permutation(induced_15(U))
    assert not is_signed_permutation(induced_15(np.kron(T, I2)))

if __name__ == "__main__":
    check_cnot_basics()
    check_conjugation_table()
    check_entanglement_creation()
    check_bell_circuit()
    check_axis_tracking()
    check_pauli_expansion()
    print("check-06-cnot: all checks passed")
