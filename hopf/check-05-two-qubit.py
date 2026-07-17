"""Checks for 05-entanglement.md (two-qubit part).

Verifies numerically (random variables, numpy):

1. Embedding rules: j z = z* j and (w j)* = -w j for complex z, w.
2. Component formula: q1 q2* = (αγ* + βδ*) + (βγ - αδ) j
   for q1 = α + βj, q2 = γ + δj.
3. Reversed embedding q = α + jβ gives j-coefficient β*γ - αδ*,
   which does NOT vanish for a separable counterexample.
4. Product states (α,β,γ,δ) = Φ_A ⊗ Φ_B have vanishing j,k components,
   and the image is the complex Hopf image of Φ_A = (a, b):
   H = (2ab*, |a|² - |b|²).
5. Converse: αδ - βγ = 0 implies vanishing j,k components, and
   Φ_A, Φ_B can be reconstructed so that Ψ = Φ_A ⊗ Φ_B.
6. Bloch vector of qubit A: x_A - i y_A = 2(αγ* + βδ*), z_A = |q1|² - |q2|²,
   and x_A² + y_A² + z_A² + C² = 1 with C = 2|αδ - βγ|.
7. Concrete examples (Bell states, meridian cos(t)|00> + sin(t)|11>),
   and the global phase e^{iφ} rotating the (j,k) components by 2φ.
8. Fiber: right multiplication by a unit quaternion preserves the image;
   on a product state it replaces the state of qubit B; Φ+ · j = Ψ-.
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

def q2c(q):
    return q[0] + 1j * q[1], q[2] + 1j * q[3]

def rand_state(n):
    v = rng.normal(size=n) + 1j * rng.normal(size=n)
    return v / np.linalg.norm(v)

def hopf_pair(psi):
    # psi = (α, β, γ, δ): q1 = α + βj, q2 = γ + δj
    return c2q(psi[0], psi[1]), c2q(psi[2], psi[3])

def H_H(q1, q2):
    return 2 * qmul(q1, qconj(q2)), q1 @ q1 - q2 @ q2

def check_embedding_rules():
    j = c2q(0, 1)
    for _ in range(20):
        z, w = rand_state(2)
        # j z = z* j
        assert np.allclose(qmul(j, c2q(z, 0)), c2q(0, z.conjugate()))
        # (w j)* = -w j
        assert np.allclose(qconj(c2q(0, w)), -c2q(0, w))

def check_component_formula():
    for _ in range(100):
        a, b, c, d = rand_state(4)
        lhs = qmul(c2q(a, b), qconj(c2q(c, d)))
        rhs = c2q(a * c.conjugate() + b * d.conjugate(), b * c - a * d)
        assert np.allclose(lhs, rhs)

def check_reversed_embedding():
    def c2q_rev(alpha, beta):
        # q = alpha + j beta = alpha + beta* j
        return c2q(alpha, beta.conjugate())
    # separable counterexample: (1, i)/√2 ⊗ (1, 1)/√2
    psi = np.kron(np.array([1, 1j]) / np.sqrt(2), np.array([1, 1]) / np.sqrt(2))
    a, b, c, d = psi
    prod = qmul(c2q_rev(a, b), qconj(c2q_rev(c, d)))
    w = prod[2] + 1j * prod[3]  # coefficient of j
    assert np.isclose(w, b.conjugate() * c - a * d.conjugate())
    assert np.isclose(w, 0.5j)  # does not vanish although separable
    # the correct embedding does vanish
    q1, q2 = hopf_pair(psi)
    assert np.allclose(qmul(q1, qconj(q2))[2:], 0)
    # general formula for the reversed embedding
    for _ in range(50):
        a, b, c, d = rand_state(4)
        prod = qmul(c2q_rev(a, b), qconj(c2q_rev(c, d)))
        w = prod[2] + 1j * prod[3]
        assert np.isclose(w, b.conjugate() * c - a * d.conjugate())

def check_product_states():
    for _ in range(100):
        phi_a = rand_state(2)
        phi_b = rand_state(2)
        q1, q2 = hopf_pair(np.kron(phi_a, phi_b))
        img, r = H_H(q1, q2)
        a, b = phi_a
        assert np.allclose(img[2:], 0)  # j,k components vanish
        assert np.isclose(img[0] + 1j * img[1], 2 * a * b.conjugate())
        assert np.isclose(r, abs(a)**2 - abs(b)**2)

def check_determinant_converse():
    for _ in range(100):
        # random state with αδ - βγ = 0
        a, b, c = rand_state(3)
        psi = np.array([a, b, c, b * c / a])
        psi /= np.linalg.norm(psi)
        q1, q2 = hopf_pair(psi)
        img, _ = H_H(q1, q2)
        assert np.allclose(img[2:], 0)
        # reconstruct Φ_A ⊗ Φ_B from the proportional rows
        lam = psi[2] / psi[0]
        phi_b = psi[:2]
        phi_a = np.array([1, lam])
        assert np.allclose(np.kron(phi_a, phi_b), psi)

def check_bloch_vector():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2)
    for _ in range(100):
        psi = rand_state(4)
        a, b, c, d = psi
        xA, yA, zA = (np.real(psi.conj() @ np.kron(s, I2) @ psi) for s in (sx, sy, sz))
        assert np.isclose(2 * (a * c.conjugate() + b * d.conjugate()), xA - 1j * yA)
        assert np.isclose(zA, abs(a)**2 + abs(b)**2 - abs(c)**2 - abs(d)**2)
        C = 2 * abs(a * d - b * c)
        assert np.isclose(xA**2 + yA**2 + zA**2 + C**2, 1)

def check_examples():
    s2 = np.sqrt(2)
    cases = [
        (np.array([1, 1, 1, 1]) / 2,      np.array([1, 0, 0, 0.]), 0),    # product
        (np.array([1, 0, 0, 1]) / s2,     np.array([0, 0, -1, 0.]), 0),   # Φ+ -> -j
        (np.array([1, 0, 0, -1]) / s2,    np.array([0, 0, 1, 0.]), 0),    # Φ- -> +j
        (np.array([0, 1, 1, 0]) / s2,     np.array([0, 0, 1, 0.]), 0),    # Ψ+ -> +j
        (np.array([0, 1, -1, 0]) / s2,    np.array([0, 0, -1, 0.]), 0),   # Ψ- -> -j
        (np.array([1, 0, 0, 1j]) / s2,    np.array([0, 0, 0, -1.]), 0),   # -> -k
    ]
    for psi, img_exp, r_exp in cases:
        img, r = H_H(*hopf_pair(psi))
        assert np.allclose(img, img_exp) and np.isclose(r, r_exp)
    # meridian cos(t)|00> + sin(t)|11> -> (0, 0, -sin 2t, 0; cos 2t), C = sin 2t
    for t in rng.uniform(0, np.pi, 20):
        psi = np.array([np.cos(t), 0, 0, np.sin(t)])
        img, r = H_H(*hopf_pair(psi))
        assert np.allclose(img, [0, 0, -np.sin(2 * t), 0])
        assert np.isclose(r, np.cos(2 * t))
        assert np.isclose(2 * abs(psi[0] * psi[3]), abs(np.sin(2 * t)))
    # global phase rotates the (j,k) components by 2φ
    for _ in range(20):
        psi = rand_state(4)
        phi = rng.uniform(0, 2 * np.pi)
        img, r = H_H(*hopf_pair(psi))
        img2, r2 = H_H(*hopf_pair(np.exp(1j * phi) * psi))
        w = img[2] + 1j * img[3]
        w2 = img2[2] + 1j * img2[3]
        assert np.allclose(img[:2], img2[:2]) and np.isclose(r, r2)
        assert np.isclose(w2, np.exp(2j * phi) * w)

def check_fiber():
    for _ in range(100):
        psi = rand_state(4)
        q1, q2 = hopf_pair(psi)
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        img, r = H_H(q1, q2)
        img2, r2 = H_H(qmul(q1, q), qmul(q2, q))
        assert np.allclose(img, img2) and np.isclose(r, r2)
    # on a product state, right multiplication replaces the state of qubit B
    for _ in range(50):
        phi_a = rand_state(2)
        phi_b = rand_state(2)
        q1, q2 = hopf_pair(np.kron(phi_a, phi_b))
        qc, qd = rand_state(2)
        q = c2q(qc, qd)
        phi_b2 = np.array(q2c(qmul(c2q(*phi_b), q)))
        q1n, q2n = hopf_pair(np.kron(phi_a, phi_b2))
        assert np.allclose(qmul(q1, q), q1n)
        assert np.allclose(qmul(q2, q), q2n)
    # Φ+ · j = Ψ-
    s2 = np.sqrt(2)
    q1, q2 = hopf_pair(np.array([1, 0, 0, 1]) / s2)
    j = c2q(0, 1)
    q1j, q2j = qmul(q1, j), qmul(q2, j)
    p1, p2 = hopf_pair(np.array([0, 1, -1, 0]) / s2)
    assert np.allclose(q1j, p1) and np.allclose(q2j, p2)

if __name__ == "__main__":
    check_embedding_rules()
    check_component_formula()
    check_reversed_embedding()
    check_product_states()
    check_determinant_converse()
    check_bloch_vector()
    check_examples()
    check_fiber()
    print("check-05-two-qubit: all checks passed")
