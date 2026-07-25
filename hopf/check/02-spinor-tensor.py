"""Checks for 02-spinor-tensor.md.

Verifies numerically (random parameters, numpy):

1. The rotor ω = (α0-α1 k) + j(β0-β1 k) gives the same Hopf image as
   article 01: ω k ω* = 2Re(α*β) i + 2Im(α*β) j + (|α|²-|β|²) k,
   with α = α0+iα1, β = β0+iβ1.  Also: this ω equals the ω of article
   01's summary (α0 k + α1 + β0 i + β1 j) right-multiplied by -k, i.e.
   both lie on the same fiber.
2. ωq with q = α0'+α1' k (α0-α1 k = γ(α0'-α1' k)) removes the k term:
   ωq = γ + (β0α0'+β1α1') j + (β0α1'-β1α0') i.
3. Trig form: ω = cosθ e^{-ak} + sinθ j e^{-bk}, q = e^{ak} gives
   ω' = cosθ + sinθ j e^{(a-b)k} and
   ω' k ω'* = sin2θ{cos(b-a) i + sin(b-a) j} + cos2θ k.
4. State vector ω' = (cos(θ/2), sin(θ/2) e^{i(b-a)})ᵀ maps by Pauli
   matrices to the Bloch vector (sinθ cos(b-a), sinθ sin(b-a), cosθ),
   and the full ω = (cos(θ/2) e^{ia}, sin(θ/2) e^{ib})ᵀ to the same.
5. Density matrix: ωω† = ½{I + sinθcos(b-a) σx + sinθsin(b-a) σy
   + cosθ σz}.
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


def rot_k(w):
    return qmul(qmul(w, [0, 0, 0, 1]), qconj(w))


K = np.array([0.0, 0.0, 0.0, 1.0])
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

for _ in range(100):
    v = rng.normal(size=4)
    v /= np.linalg.norm(v)
    a0, a1, b0, b1 = v
    alpha = a0 + a1*1j
    beta = b0 + b1*1j

    # ω = (α0-α1 k) + (β0 j - β1 i)
    w = np.array([a0, -b1, b0, -a1])
    wkw = rot_k(w)

    # 1. same Hopf image as article 01's formula, and same fiber as 01's ω
    ab = np.conj(alpha) * beta
    assert np.allclose(wkw, [0, 2*ab.real, 2*ab.imag, abs(alpha)**2 - abs(beta)**2])
    w01 = np.array([a1, b0, b1, a0])  # article 01: α0 k + α1 + β0 i + β1 j
    assert np.allclose(w, qmul(w01, [0, 0, 0, -1]))
    assert np.allclose(rot_k(w01), wkw)

    # 2. ωq removes the k term
    gamma = abs(alpha)
    ap0, ap1 = a0/gamma, a1/gamma
    q = np.array([ap0, 0, 0, ap1])
    wq = qmul(w, q)
    assert np.allclose(wq, [gamma, b0*ap1 - b1*ap0, b0*ap0 + b1*ap1, 0])

    # 3. trig form
    th = rng.uniform(0, np.pi)
    pa, pb = rng.uniform(0, 2*np.pi, 2)
    # ω = cosθ e^{-ak} + sinθ j e^{-bk} = cosθ(cos a - sin a k) + sinθ(cos b j - sin b i)
    wt = np.array([np.cos(th)*np.cos(pa), -np.sin(th)*np.sin(pb),
                   np.sin(th)*np.cos(pb), -np.cos(th)*np.sin(pa)])
    qt = np.array([np.cos(pa), 0, 0, np.sin(pa)])  # e^{ak}
    wp = qmul(wt, qt)
    # ω' = cosθ + sinθ j e^{(a-b)k} = cosθ + sinθ{cos(a-b) j + sin(a-b) i}
    assert np.allclose(wp, [np.cos(th), np.sin(th)*np.sin(pa-pb),
                            np.sin(th)*np.cos(pa-pb), 0])
    bloch = np.array([np.sin(2*th)*np.cos(pb-pa), np.sin(2*th)*np.sin(pb-pa),
                      np.cos(2*th)])
    assert np.allclose(rot_k(wp)[1:], bloch)
    assert np.allclose(rot_k(wt)[1:], bloch)

    # 4. state vector with half angle → Bloch vector
    psi = np.array([np.cos(th/2)*np.exp(1j*pa), np.sin(th/2)*np.exp(1j*pb)])
    psi_p = np.array([np.cos(th/2), np.sin(th/2)*np.exp(1j*(pb-pa))])
    bloch_half = np.array([np.sin(th)*np.cos(pb-pa), np.sin(th)*np.sin(pb-pa),
                           np.cos(th)])
    for p in (psi, psi_p):
        xyz = [np.conj(p) @ s @ p for s in (sigma_x, sigma_y, sigma_z)]
        assert np.allclose(np.imag(xyz), 0)
        assert np.allclose(np.real(xyz), bloch_half)

    # 5. density matrix decomposition
    dens = np.outer(psi, np.conj(psi))
    expected = (np.eye(2) + bloch_half[0]*sigma_x + bloch_half[1]*sigma_y
                + bloch_half[2]*sigma_z) / 2
    assert np.allclose(dens, expected)

print("check-02-spinor-tensor: all checks passed")
