"""Checks for 01-quaternion.md.

Verifies numerically (random unit quaternions, numpy):

1. The expansion of ω k ω* into components
   2(ω0ω2+ω1ω3) i + 2(ω2ω3-ω0ω1) j + {(ω0²+ω3²)-(ω1²+ω2²)} k
   and that its real part vanishes and its norm is 1.
2. The complex-pair form: with α = ω3+ω0 i, β = ω1+ω2 i,
   ω k ω* = 2Re(α*β) i + 2Im(α*β) j + (|α|²-|β|²) k.
3. The embedding ω = (α0 k + α1) + j (β0 k + β1) reproduces ω, and the
   summary derivation ω k ω* = 2 α_ω β_ω* i + (|α|²-|β|²) k.
4. The Pauli-matrix form x = Ψ†σx Ψ, y = Ψ†σy Ψ, z = Ψ†σz Ψ with
   Ψ = (α, β)ᵀ gives the same coordinates.
5. Fibers: q = u+v k (u²+v²=1) fixes k, and (ωq) k (ωq)* = ω k ω*.
"""

import numpy as np

rng = np.random.default_rng(42)


def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.array([
        a0*b0 - a1*b1 - a2*b2 - a3*b3,
        a0*b1 + a1*b0 + a2*b3 - a3*b2,
        a0*b2 - a1*b3 + a2*b0 + a3*b1,
        a0*b3 + a1*b2 - a2*b1 + a3*b0,
    ])


def qconj(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])


K = np.array([0.0, 0.0, 0.0, 1.0])
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

for _ in range(100):
    w = rng.normal(size=4)
    w /= np.linalg.norm(w)
    w0, w1, w2, w3 = w

    wkw = qmul(qmul(w, K), qconj(w))

    # 1. component expansion, purity, unit norm
    expected = np.array([
        0.0,
        2*(w0*w2 + w1*w3),
        2*(w2*w3 - w0*w1),
        (w0**2 + w3**2) - (w1**2 + w2**2),
    ])
    assert np.allclose(wkw, expected)
    assert np.isclose(np.linalg.norm(wkw), 1)

    # 2. complex-pair form
    alpha = w3 + w0*1j
    beta = w1 + w2*1j
    ab = np.conj(alpha) * beta
    assert np.allclose(wkw, [0, 2*ab.real, 2*ab.imag, abs(alpha)**2 - abs(beta)**2])

    # 3. embedding: α_ω = α0 k + α1, β_ω = β0 k + β1, ω = α_ω + j β_ω
    a0, a1 = alpha.real, alpha.imag
    b0, b1 = beta.real, beta.imag
    aw = np.array([a1, 0, 0, a0])
    bw = np.array([b1, 0, 0, b0])
    J = np.array([0.0, 0.0, 1.0, 0.0])
    assert np.allclose(w, aw + qmul(J, bw))
    # ω k ω* = 2 α_ω β_ω* i + (|α|²-|β|²) k
    I = np.array([0.0, 1.0, 0.0, 0.0])
    summary = 2*qmul(qmul(aw, qconj(bw)), I) + (abs(alpha)**2 - abs(beta)**2)*K
    assert np.allclose(wkw, summary)

    # 4. Pauli-matrix coordinates
    psi = np.array([alpha, beta])
    xyz = [np.conj(psi) @ s @ psi for s in (sigma_x, sigma_y, sigma_z)]
    assert np.allclose(np.imag(xyz), 0)
    assert np.allclose(wkw[1:], np.real(xyz))

    # 5. fiber: q = u+v k fixes k and does not move the image
    t = rng.uniform(0, 2*np.pi)
    q = np.array([np.cos(t), 0, 0, np.sin(t)])
    assert np.allclose(qmul(qmul(q, K), qconj(q)), K)
    wq = qmul(w, q)
    assert np.allclose(qmul(qmul(wq, K), qconj(wq)), wkw)

print("check-01-quaternion: all checks passed")
