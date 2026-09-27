"""Checks for 05-entanglement.md (three-qubit part).

Verifies numerically (random variables, numpy):

1. Embedding C^4 ≅ O with complex structure = right multiplication by e1
   and C-basis {1, e2, e4, e6}: x = z1 + e2 z2 + e4 z3 + e6 z4;
   emb/unemb round trip, including the signs e2 e1 = -e3, e4 e1 = -e5,
   e6 e1 = +e7.
2. Right C-linearity: emb(z a) = emb(z) · a (right multiplication by a ∈ C).
3. Nested Cayley-Dickson trap: e4 (e2 z) = -e6 z* ≠ ±e6 z.
4. Identities: x*(xy) = |x|² y, x(x*y) = |x|² y, <xy, z> = <y, x*z>.
5. A|BC separable states Φ_A ⊗ X have vanishing e2..e7 components,
   o1 = x a, o2 = x b, and the image is (2a*b, |a|² - |b|²) for Φ_A = (a, b).
6. Converse: o1* o2 ∈ C implies the blocks are complex-proportional
   (all 2x2 minors vanish); generic states have nonzero e2..e7 components.
7. Concrete examples: |000>, |0>⊗Bell_BC (same image), GHZ -> e6,
   W -> -(2/3)(e2 + e4) with r = 1/3, Bell_AC⊗|0>_B -> e2, Bell_AB⊗|0>_C -> e4.
8. Bloch vector of qubit A via σ ⊗ I_4, the complex component formula
   π_C(o1* o2) = Σ u_k* v_k, and |r_A|² + C² = 1
   with C = |e2..e7 components of the image|.
8b. The image direction is not an invariant: e^{iπ/4} GHZ and the product
   state (|0>+i|1>)_C ⊗ Bell_AB share the image (e7; 0).
9. Fiber over a separable image point (hopf/04 construction): every point
   is A|BC separable with the same r_A; left multiplication does not
   preserve the image.
10. SU(3) as automorphisms fixing e1: a random U ∈ SU(3) acting on the
    C-basis (e2, e4, e6) preserves multiplication; det ≠ 1 breaks it;
    diagonal phases work iff θ1 + θ2 + θ3 = 0; the automorphism commutes
    with the Hopf map, acts on states as I_2 ⊗ U_BC (unitary), and a cyclic
    permutation maps the GHZ image e6 to the Bell_AC image e2, sending
    the GHZ state itself to Bell_AC ⊗ |0>_B and then to Bell_AB ⊗ |0>_C.
"""

import numpy as np
from common.octonion import L

rng = np.random.default_rng(0)

def mul(x, y):
    M = sum(x[i] * L[i] for i in range(8))
    return M @ y

def conj(x):
    y = -x.copy()
    y[0] = x[0]
    return y

def e(i):
    v = np.zeros(8)
    v[i] = 1.0
    return v

def c2o(z):
    # complex number as an octonion in the (1, e1) plane
    return z.real * e(0) + z.imag * e(1)

def emb(z):
    # C^4 -> O with C-basis {1, e2, e4, e6}, complex structure = right mult by e1
    v = np.zeros(8)
    for zi, b in zip(z, [0, 2, 4, 6]):
        v += mul(e(b), c2o(zi))
    return v

def unemb(v):
    # inverse of emb; e2 e1 = -e3, e4 e1 = -e5, e6 e1 = e7
    return np.array([v[0] + 1j * v[1], v[2] - 1j * v[3],
                     v[4] - 1j * v[5], v[6] + 1j * v[7]])

def H_O(o1, o2):
    return 2 * mul(conj(o1), o2), o1 @ o1 - o2 @ o2

def rand_state(n):
    v = rng.normal(size=n) + 1j * rng.normal(size=n)
    return v / np.linalg.norm(v)

def blocks(psi):
    # psi in C^8, basis index 4a + 2b + c for the qubits (A, B, C)
    return emb(psi[:4]), emb(psi[4:])

def minors(psi):
    # 2x2 minors of the 2x4 matrix (psi[:4]; psi[4:])
    z0, z1 = psi[:4], psi[4:]
    return np.array([z0[i] * z1[j] - z0[j] * z1[i]
                     for i in range(4) for j in range(i + 1, 4)])

def check_embedding():
    # sign convention on the 4th coordinate
    assert np.allclose(emb(np.array([0, 1j, 0, 0])), -e(3))
    assert np.allclose(emb(np.array([0, 0, 1j, 0])), -e(5))
    assert np.allclose(emb(np.array([0, 0, 0, 1j])), e(7))
    for _ in range(50):
        z = rng.normal(size=4) + 1j * rng.normal(size=4)
        assert np.allclose(unemb(emb(z)), z)
        v = rng.normal(size=8)
        assert np.allclose(emb(unemb(v)), v)

def check_c_linearity():
    for _ in range(50):
        a = complex(*rng.normal(size=2))
        z = rng.normal(size=4) + 1j * rng.normal(size=4)
        assert np.allclose(emb(z * a), mul(emb(z), c2o(a)))

def check_nested_trap():
    for _ in range(20):
        z = complex(*rng.normal(size=2))
        nested = mul(e(4), mul(e(2), c2o(z)))
        assert np.allclose(nested, -mul(e(6), c2o(z.conjugate())))
        if abs(z.imag) > 1e-6:
            assert not np.allclose(nested, mul(e(6), c2o(z)))
            assert not np.allclose(nested, -mul(e(6), c2o(z)))

def check_identities():
    for _ in range(50):
        x = rng.normal(size=8)
        y = rng.normal(size=8)
        z = rng.normal(size=8)
        n2 = x @ x
        assert np.allclose(mul(conj(x), mul(x, y)), n2 * y)
        assert np.allclose(mul(x, mul(conj(x), y)), n2 * y)
        assert np.isclose(mul(x, y) @ z, y @ mul(conj(x), z))

def check_separable():
    for _ in range(100):
        a, b = rand_state(2)
        X = rand_state(4)
        psi = np.kron([a, b], X)
        o1, o2 = blocks(psi)
        x = emb(X)
        assert np.allclose(o1, mul(x, c2o(a))) and np.allclose(o2, mul(x, c2o(b)))
        img, r = H_O(o1, o2)
        assert np.allclose(img[2:], 0)  # e2..e7 components vanish
        assert np.isclose(img[0] + 1j * img[1], 2 * a.conjugate() * b)
        assert np.isclose(r, abs(a)**2 - abs(b)**2)

def check_converse():
    for _ in range(100):
        o1 = rng.normal(size=8)
        c = complex(*rng.normal(size=2))
        o2 = mul(o1, c2o(c)) / (o1 @ o1)
        img, _ = H_O(o1, o2)
        assert np.allclose(img[2:], 0)  # o1* o2 ∈ C
        assert np.allclose(img, 2 * c2o(c))
        psi = np.concatenate([unemb(o1), unemb(o2)])
        assert np.allclose(minors(psi), 0)  # blocks are proportional
    # generic states are not separable
    count = 0
    for _ in range(20):
        img, _ = H_O(*blocks(rand_state(8)))
        if not np.allclose(img[2:], 0):
            count += 1
    assert count > 0

def state(indices, coeffs):
    psi = np.zeros(8, dtype=complex)
    for i, c in zip(indices, coeffs):
        psi[i] = c
    return psi / np.linalg.norm(psi)

def check_examples():
    def C_val(psi):
        img, _ = H_O(*blocks(psi))
        return np.linalg.norm(img[2:])
    # |000> and |0>⊗Bell_BC: both at the north pole (same fiber)
    for psi in [state([0], [1]), state([0, 3], [1, 1])]:
        img, r = H_O(*blocks(psi))
        assert np.allclose(img, 0) and np.isclose(r, 1)
        assert np.isclose(C_val(psi), 0)
    # GHZ = (|000> + |111>)/√2 -> (e6; 0), C = 1
    ghz = state([0, 7], [1, 1])
    img, r = H_O(*blocks(ghz))
    assert np.allclose(img, e(6)) and np.isclose(r, 0)
    assert np.isclose(C_val(ghz), 1)
    # W = (|001> + |010> + |100>)/√3 -> (-(2/3)(e2 + e4); 1/3), C = 2√2/3
    w = state([1, 2, 4], [1, 1, 1])
    img, r = H_O(*blocks(w))
    assert np.allclose(img, -(e(2) + e(4)) * 2 / 3) and np.isclose(r, 1 / 3)
    assert np.isclose(C_val(w), 2 * np.sqrt(2) / 3)
    # Bell_AC ⊗ |0>_B = (|000> + |101>)/√2 -> (e2; 0)
    img, r = H_O(*blocks(state([0, 5], [1, 1])))
    assert np.allclose(img, e(2)) and np.isclose(r, 0)
    # Bell_AB ⊗ |0>_C = (|000> + |110>)/√2 -> (e4; 0)
    img, r = H_O(*blocks(state([0, 6], [1, 1])))
    assert np.allclose(img, e(4)) and np.isclose(r, 0)

def check_image_direction_not_invariant():
    # e^{iπ/4} GHZ and the product state (|0>+i|1>)_C ⊗ Bell_AB
    # share the image (e7; 0) although only the former is genuinely 3-way
    ghz_phase = np.exp(1j * np.pi / 4) * state([0, 7], [1, 1])
    o1, o2 = blocks(ghz_phase)
    assert np.allclose(o1, (e(0) + e(1)) / 2)
    assert np.allclose(o2, (e(6) + e(7)) / 2)
    img, r = H_O(o1, o2)
    assert np.allclose(img, e(7)) and np.isclose(r, 0)
    # (|000> + i|001> + |110> + i|111>)/2
    prod = state([0, 1, 6, 7], [1, 1j, 1, 1j])
    C = np.array([1, 1j]) / np.sqrt(2)
    AB = np.array([1, 0, 0, 1]) / np.sqrt(2)
    assert np.allclose(prod, np.kron(AB, C))
    o1, o2 = blocks(prod)
    assert np.allclose(o1, (e(0) - e(3)) / 2)
    assert np.allclose(o2, (e(4) + e(7)) / 2)
    img, r = H_O(o1, o2)
    assert np.allclose(img, e(7)) and np.isclose(r, 0)

def bloch_A(psi):
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    I4 = np.eye(4)
    return np.array([np.real(psi.conj() @ np.kron(s, I4) @ psi)
                     for s in (sx, sy, sz)])

def check_bloch_vector():
    for _ in range(100):
        psi = rand_state(8)
        img, r = H_O(*blocks(psi))
        xA, yA, zA = bloch_A(psi)
        assert np.isclose(img[0] + 1j * img[1], xA + 1j * yA)
        assert np.isclose(r, zA)
        # π_C(o1* o2) = Σ u_k* v_k for the blocks Ψ0 = (u_k), Ψ1 = (v_k)
        o1, o2 = blocks(psi)
        prod = mul(conj(o1), o2)
        assert np.isclose(prod[0] + 1j * prod[1], psi[:4].conjugate() @ psi[4:])
        C = np.linalg.norm(img[2:])
        assert np.isclose(xA**2 + yA**2 + zA**2 + C**2, 1)

def check_fiber():
    # fiber over a separable image point (hopf/04 construction)
    for _ in range(50):
        a, b = rand_state(2)
        c_o = c2o(2 * a.conjugate() * b)
        r = abs(a)**2 - abs(b)**2
        o1 = rng.normal(size=8)
        o1 *= np.sqrt((1 + r) / 2) / np.linalg.norm(o1)
        o2 = mul(o1, c_o) / (1 + r)
        img, r2 = H_O(o1, o2)
        assert np.allclose(img, c_o) and np.isclose(r2, r)
        psi = np.concatenate([unemb(o1), unemb(o2)])
        assert np.allclose(minors(psi), 0)  # every fiber point is A|BC separable
        assert np.allclose(bloch_A(psi), [img[0], img[1], r])
    # left multiplication does not preserve the image
    count = 0
    for _ in range(20):
        o1, o2 = blocks(rand_state(8))
        q = rng.normal(size=8)
        q /= np.linalg.norm(q)
        img, _ = H_O(o1, o2)
        img2, _ = H_O(mul(q, o1), mul(q, o2))
        if not np.allclose(img, img2):
            count += 1
    assert count > 0

def phi_of(U):
    # C-linear map fixing 1 and acting by U on the C-basis (e2, e4, e6)
    def phi(v):
        z = unemb(v)
        return emb(np.concatenate([[z[0]], U @ z[1:]]))
    return phi

def is_automorphism(phi, trials=20):
    for _ in range(trials):
        x = rng.normal(size=8)
        y = rng.normal(size=8)
        if not np.allclose(phi(mul(x, y)), mul(phi(x), phi(y))):
            return False
    return True

def check_su3():
    # random SU(3)
    A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    U, _ = np.linalg.qr(A)
    U /= np.linalg.det(U) ** (1 / 3)
    assert np.isclose(np.linalg.det(U), 1)
    phi = phi_of(U)
    assert is_automorphism(phi)
    assert np.allclose(phi(e(1)), e(1))
    # det ≠ 1 (an element of U(3) \ SU(3)) breaks multiplication
    assert not is_automorphism(phi_of(U * np.exp(0.7j)))
    # diagonal phases: automorphism iff θ1 + θ2 + θ3 = 0
    # (e2 e^{iθ1})(e4 e^{iθ2}) = e6 e^{-i(θ1+θ2)}
    t1, t2 = rng.uniform(-np.pi, np.pi, 2)
    lhs = mul(mul(e(2), c2o(np.exp(1j * t1))), mul(e(4), c2o(np.exp(1j * t2))))
    assert np.allclose(lhs, mul(e(6), c2o(np.exp(-1j * (t1 + t2)))))
    t1, t2 = rng.uniform(-np.pi, np.pi, 2)
    D = np.diag(np.exp(1j * np.array([t1, t2, -(t1 + t2)])))
    assert is_automorphism(phi_of(D))
    assert not is_automorphism(phi_of(np.diag(np.exp(1j * np.array([t1, t2, t2])))))
    # commutes with the Hopf map, acts on states as I_2 ⊗ U_BC
    U4 = np.zeros((4, 4), dtype=complex)
    for k in range(4):
        zk = np.zeros(4, dtype=complex); zk[k] = 1
        U4[:, k] = unemb(phi(emb(zk)))
    assert np.allclose(U4 @ U4.conj().T, np.eye(4))  # U_BC is unitary
    for _ in range(20):
        psi = rand_state(8)
        o1, o2 = blocks(psi)
        img, r = H_O(o1, o2)
        img2, r2 = H_O(phi(o1), phi(o2))
        assert np.allclose(img2, phi(img)) and np.isclose(r2, r)
        assert np.allclose(img2[:2], img[:2])  # the (1, e1) part is fixed
        psi2 = np.kron(np.eye(2), U4) @ psi
        assert np.allclose(blocks(psi2)[0], phi(o1))
        assert np.allclose(blocks(psi2)[1], phi(o2))
        assert np.allclose(bloch_A(psi2), bloch_A(psi))  # A is untouched
    # cyclic permutation (z1, z2, z3) -> (z3, z1, z2) ∈ SU(3)
    P = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)
    assert np.isclose(np.linalg.det(P), 1)
    phiP = phi_of(P)
    assert is_automorphism(phiP)
    assert np.allclose(phiP(e(6)), e(2))  # GHZ image -> Bell_AC image
    assert np.allclose(phiP(e(2)), e(4))  # Bell_AC image -> Bell_AB image
    # on states, I_2 ⊗ U_P maps GHZ -> Bell_AC ⊗ |0>_B -> Bell_AB ⊗ |0>_C
    U4P = np.zeros((4, 4), dtype=complex)
    for k in range(4):
        zk = np.zeros(4, dtype=complex); zk[k] = 1
        U4P[:, k] = unemb(phiP(emb(zk)))
    step = np.kron(np.eye(2), U4P)
    ghz = state([0, 7], [1, 1])
    assert np.allclose(step @ ghz, state([0, 5], [1, 1]))
    assert np.allclose(step @ step @ ghz, state([0, 6], [1, 1]))

if __name__ == "__main__":
    check_embedding()
    check_c_linearity()
    check_nested_trap()
    check_identities()
    check_separable()
    check_converse()
    check_examples()
    check_image_direction_not_invariant()
    check_bloch_vector()
    check_fiber()
    check_su3()
    print("check-05-three-qubit: all checks passed")
