"""Checks for PLAN.md section 4 (the boundary of real amplitudes).

Verifies numerically (numpy):

1. The Fourier boundary. The DFT matrix over Z_N is real iff N <= 2, because
   its entries are characters. Over (Z_2)^n every character takes values in
   {+1, -1}, and the character table divided by sqrt(2^n) is exactly H^(x)n.
   Among all finite abelian groups of the form Z_{d1} x ... x Z_{dk}, the
   Fourier matrix is real iff every d is at most 2.
2. Deutsch-Jozsa and Bernstein-Vazirani run entirely inside real amplitudes.
3. Simon's algorithm (period finding over (Z_2)^n) is real: every gate and
   every intermediate state vector is real, and the outcomes are orthogonal
   to the hidden period s.
4. Grover is a real rotation. The oracle and the diffusion operator are real
   reflections (det = -1 on the 2-dimensional invariant plane), their product
   is a rotation by 2*theta with sin(theta) = sqrt(M/N), and the amplitude on
   the marked subspace after k iterations is sin((2k+1)*theta).
5. Teleportation and superdense coding are real; a Pauli tensor product is a
   real matrix iff it contains an even number of Y factors, so the real Pauli
   group covers the stabilizers of the Steane, Shor and 5-qubit codes.
6. Shor is not natively real. QFT_N is complex for N >= 3, and order finding
   for a = 7 mod 15 puts its answer in the phases e^(2*pi*i*s/r). The rebit
   encoding U = A + iB -> I(x)A + J(x)B still runs it on one extra qubit with
   identical outcome probabilities: the power is preserved, the real form is
   not native.
"""

from functools import reduce
from itertools import product

import numpy as np

rng = np.random.default_rng(42)

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=float)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]], dtype=float)
H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
J = np.array([[0, -1], [1, 0]], dtype=float)

PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}


def kron(*ms):
    return reduce(np.kron, ms)


def is_real(m):
    return np.allclose(np.imag(m), 0)


def pauli(s):
    return kron(*[PAULI[c] for c in s])


def dft(n):
    """Fourier matrix of the cyclic group Z_n."""
    j, k = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    return np.exp(2j * np.pi * j * k / n) / np.sqrt(n)


def dft_abelian(dims):
    """Fourier matrix of Z_{d1} x ... x Z_{dk}."""
    return kron(*[dft(d) for d in dims])


# ---------------------------------------------------------------- 1. boundary

def check_fourier_boundary():
    # over Z_N the transform is real only for N = 1, 2
    for n in range(1, 13):
        assert is_real(dft(n)) == (n <= 2), n
    assert np.allclose(dft(2), H)

    # characters of (Z_2)^n are (-1)^(a.x), so they take values in {+1, -1};
    # the character table normalised is the Hadamard transform
    for n in range(1, 6):
        table = np.array(
            [[(-1) ** (bin(a & x).count("1")) for x in range(2**n)] for a in range(2**n)],
            dtype=float,
        )
        assert set(np.unique(table)) <= {-1.0, 1.0}
        Hn = kron(*[H] * n)
        assert np.allclose(table / np.sqrt(2**n), Hn)
        assert is_real(Hn) and np.allclose(Hn @ Hn.T, np.eye(2**n))

    # general finite abelian group: real iff every element has order <= 2
    for dims in [(2,), (3,), (4,), (2, 2), (2, 3), (2, 2, 2), (2, 4), (5,), (2, 2, 3)]:
        exponent_two = all(d <= 2 for d in dims)
        assert is_real(dft_abelian(dims)) == exponent_two, dims


# ------------------------------------------------------- 2-3. oracle algorithms

def phase_oracle(f, n):
    """diag((-1)^f(x)) on n qubits: real and diagonal."""
    return np.diag([(-1.0) ** f(x) for x in range(2**n)])


def xor_oracle(f, n, m):
    """|x>|y> -> |x>|y xor f(x)> on n + m qubits: a real permutation."""
    dim = 2 ** (n + m)
    P = np.zeros((dim, dim))
    for x in range(2**n):
        for y in range(2**m):
            P[(x << m) | (y ^ f(x)), (x << m) | y] = 1.0
    return P


def check_deutsch_jozsa():
    n = 3
    Hn = kron(*[H] * n)
    ket0 = np.zeros(2**n)
    ket0[0] = 1.0

    def run(f):
        Uf = phase_oracle(f, n)
        assert is_real(Uf) and is_real(Hn)
        psi = Hn @ ket0
        assert is_real(psi)
        psi = Uf @ psi
        assert is_real(psi)
        psi = Hn @ psi
        assert is_real(psi)
        return psi

    # constant: all the amplitude returns to |0...0>
    psi = run(lambda x: 0)
    assert np.isclose(abs(psi[0]) ** 2, 1.0)
    # balanced: |0...0> never appears
    half = {x for x in rng.permutation(2**n)[: 2 ** (n - 1)]}
    psi = run(lambda x: int(x in half))
    assert np.isclose(psi[0], 0.0)

    # Bernstein-Vazirani: f(x) = a.x recovers a in one query
    for a in range(2**n):
        psi = run(lambda x, a=a: bin(a & x).count("1") % 2)
        assert np.isclose(abs(psi[a]) ** 2, 1.0)
        assert is_real(psi)


def check_simon():
    n = 3
    s = 0b101                                  # the hidden period
    # f(x) = f(x xor s), two-to-one
    f = lambda x: min(x, x ^ s)

    Hn = kron(*[H] * n)
    Uf = xor_oracle(f, n, n)
    assert is_real(Uf) and np.allclose(Uf @ Uf.T, np.eye(2 ** (2 * n)))

    psi = np.zeros(2 ** (2 * n))
    psi[0] = 1.0
    for stage in (kron(Hn, np.eye(2**n)), Uf, kron(Hn, np.eye(2**n))):
        assert is_real(stage)
        psi = stage @ psi
        assert is_real(psi)                    # never leaves the real amplitudes

    # marginal distribution of the first register
    probs = (psi.reshape(2**n, 2**n) ** 2).sum(axis=1)
    assert np.isclose(probs.sum(), 1.0)
    for y in range(2**n):
        orthogonal = bin(y & s).count("1") % 2 == 0
        if orthogonal:
            assert np.isclose(probs[y], 1 / 2 ** (n - 1))
        else:
            assert np.isclose(probs[y], 0.0)   # y.s = 1 never occurs

    # the transform used is the Fourier transform of (Z_2)^n, hence real
    assert is_real(Hn)


# ------------------------------------------------------------------ 4. Grover

def check_grover():
    n, marked = 6, [11, 37]
    N, M = 2**n, len(marked)
    theta = np.arcsin(np.sqrt(M / N))

    Uf = phase_oracle(lambda x: int(x in marked), n)
    s = np.ones(N) / np.sqrt(N)
    D = 2 * np.outer(s, s) - np.eye(N)         # diffusion
    G = D @ Uf

    for m in (Uf, D, G):
        assert is_real(m) and np.allclose(m @ m.T, np.eye(N))

    good = np.zeros(N)
    good[marked] = 1 / np.sqrt(M)
    bad = np.ones(N)
    bad[marked] = 0
    bad /= np.linalg.norm(bad)
    B = np.column_stack([good, bad])           # orthonormal invariant plane

    assert np.allclose(s, np.sin(theta) * good + np.cos(theta) * bad)

    # both operators are reflections when restricted to the plane
    for m in (Uf, D):
        r = B.T @ m @ B
        assert np.allclose(m @ B, B @ r)       # the plane is invariant
        assert np.isclose(np.linalg.det(r), -1.0)

    # their product is a rotation by twice the angle between the mirrors
    g = B.T @ G @ B
    assert np.isclose(np.linalg.det(g), 1.0)
    rot = np.array([[np.cos(2 * theta), np.sin(2 * theta)],
                    [-np.sin(2 * theta), np.cos(2 * theta)]])
    assert np.allclose(g, rot)

    # amplitude on the marked subspace after k iterations is sin((2k+1)theta)
    psi = s.copy()
    for k in range(12):
        assert np.isclose(good @ psi, np.sin((2 * k + 1) * theta))
        assert is_real(psi)
        psi = G @ psi

    # the optimal number of iterations really does find a marked item
    k_opt = int(round((np.pi / 2 - theta) / (2 * theta)))
    psi = np.linalg.matrix_power(G, k_opt) @ s
    assert sum(psi[w] ** 2 for w in marked) > 0.99

    # the rotation angle is not a multiple of pi/4, so Grover is not Clifford
    assert not np.isclose((2 * theta) % (np.pi / 4), 0)


# ---------------------------------------------------- 5. Paulis, codes, protocols

def check_real_pauli_group():
    # a Pauli tensor product is real iff the number of Y factors is even
    for n in (1, 2, 3):
        for s in product("IXYZ", repeat=n):
            assert is_real(pauli("".join(s))) == (s.count("Y") % 2 == 0), s

    # single-qubit real Pauli group is D_4 (order 8), not Q_8: hopf/07
    elems = [c * m for c in (1, -1) for m in (I2, X, Z, X @ Z)]
    assert all(is_real(m) for m in elems)
    assert len({m.tobytes() for m in elems}) == 8

    codes = {
        # Steane [[7,1,3]]: CSS, X-type and Z-type from the Hamming code
        "steane": ["IIIXXXX", "IXXIIXX", "XIXIXIX",
                   "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"],
        # Shor [[9,1,3]]: CSS
        "shor": ["ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII",
                 "IIIIIIZZI", "IIIIIIIZZ",
                 "XXXXXXIII", "IIIXXXXXX"],
        # 5-qubit [[5,1,3]]: not CSS, but still free of single Y factors
        "five": ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"],
    }
    for name, gens in codes.items():
        ms = [pauli(g) for g in gens]
        for m in ms:
            assert is_real(m), (name, m)
            assert np.allclose(m @ m, np.eye(m.shape[0]))
        for a, b in product(ms, repeat=2):
            assert np.allclose(a @ b, b @ a), name        # abelian stabilizer
        # the code space has the expected dimension 2^(n-k)
        n = len(gens[0])
        proj = reduce(lambda p, m: p @ (np.eye(2**n) + m) / 2, ms, np.eye(2**n))
        assert np.isclose(np.trace(proj), 2 ** (n - len(gens)))

    # Y appears in the 5-qubit stabilizer only in even numbers, e.g. XZZXI*IXZZX
    prod5 = pauli("XZZXI") @ pauli("IXZZX")
    assert is_real(prod5)


def check_protocols():
    CNOT = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=float)
    bell = CNOT @ kron(H, I2) @ np.eye(4)[:, 0]
    assert is_real(bell)

    # superdense coding: the four encodings are real and give the Bell basis
    encodings = [I2, X, Z, X @ Z]
    states = [kron(e, I2) @ bell for e in encodings]
    for v in states:
        assert is_real(v)
    G = np.array([[u @ v for v in states] for u in states])
    assert np.allclose(np.abs(G), np.eye(4))          # mutually orthogonal

    # teleportation of a real qubit, using only real gates and X, Z corrections
    decode = kron(H, I2, I2) @ kron(CNOT, I2)
    for _ in range(50):
        a = rng.normal(size=2)
        a /= np.linalg.norm(a)
        psi = decode @ kron(a, bell)
        assert is_real(psi)
        for m in range(4):                            # the two measured bits
            branch = psi.reshape(4, 2)[m]
            p = branch @ branch
            correction = [I2, X, Z, Z @ X][m]
            assert np.allclose(correction @ branch / np.sqrt(p), a)
            assert np.isclose(p, 0.25)


# ---------------------------------------------------------- 6. Shor is complex

def rebit_encode_matrix(U):
    return np.kron(I2, U.real) + np.kron(J, U.imag)


def rebit_encode_state(v):
    return np.concatenate([v.real, v.imag])


def check_qft_is_complex():
    for n in range(3, 9):
        assert not is_real(dft(n))
    # the imaginary part is not a removable phase convention: for N >= 3 no
    # diagonal rephasing of rows and columns makes the DFT real
    F = dft(4)
    assert not is_real(F)
    # every entry of row 1 has a distinct argument, so a single global phase
    # cannot flatten it
    args = np.angle(F[1])
    assert len({round(float(a % (2 * np.pi)), 9) for a in args}) == 4


def check_shor_order_finding():
    t, w = 4, 4                      # counting qubits, work qubits
    a, mod, r = 7, 15, 4             # 7^4 = 2401 = 160*15 + 1
    assert pow(a, r, mod) == 1 and all(pow(a, k, mod) != 1 for k in range(1, r))

    dim_c, dim_w = 2**t, 2**w
    dim = dim_c * dim_w

    # modular exponentiation: |c>|y> -> |c>|a^c y mod 15>, a real permutation
    P = np.zeros((dim, dim))
    for c in range(dim_c):
        for y in range(dim_w):
            ny = (pow(a, c, mod) * y) % mod if y < mod else y
            P[c * dim_w + ny, c * dim_w + y] = 1.0
    assert is_real(P) and np.allclose(P @ P.T, np.eye(dim))

    Ht = kron(*[H] * t)
    iqft = dft(dim_c).conj().T
    assert is_real(Ht) and not is_real(iqft)          # the only complex stage

    psi = np.zeros(dim, dtype=complex)
    psi[0 * dim_w + 1] = 1.0                          # |0...0>|1>
    psi = kron(Ht, np.eye(dim_w)) @ psi
    assert is_real(psi)
    psi = P @ psi
    assert is_real(psi)                               # still real up to here
    psi = kron(iqft, np.eye(dim_w)) @ psi
    assert not is_real(psi)                           # complex from now on

    probs = (np.abs(psi.reshape(dim_c, dim_w)) ** 2).sum(axis=1)
    peaks = [c for c in range(dim_c) if probs[c] > 0.1]
    assert peaks == [k * dim_c // r for k in range(r)] == [0, 4, 8, 12]
    assert np.isclose(sum(probs[c] for c in peaks), 1.0)

    # the answer sits in the phases: U has eigenvalues exp(2*pi*i*s/r)
    U = np.zeros((dim_w, dim_w))
    for y in range(dim_w):
        U[(a * y) % mod if y < mod else y, y] = 1.0
    orbit = [1, 7, 4, 13]
    for s in range(r):
        v = sum(np.exp(-2j * np.pi * s * k / r) * np.eye(dim_w)[orbit[k]] for k in range(r))
        v /= np.linalg.norm(v)
        eig = np.exp(2j * np.pi * s / r)
        assert np.allclose(U @ v, eig * v)
        # only s = 0 and s = r/2 give real eigenvalues; the phases that carry
        # the order are genuinely complex
        assert is_real(eig) == (2 * s % r == 0), s
    # U is a real matrix, so its spectrum is closed under conjugation: the
    # phase s/r always comes paired with -s/r
    spectrum = np.linalg.eigvals(U)
    for lam in spectrum:
        assert np.min(np.abs(spectrum - lam.conj())) < 1e-9

    # the rebit encoding runs the same circuit with real amplitudes only
    full = kron(iqft, np.eye(dim_w)) @ P @ kron(Ht, np.eye(dim_w))
    Ur = rebit_encode_matrix(full)
    assert is_real(Ur) and np.allclose(Ur @ Ur.T, np.eye(2 * dim))

    start = np.zeros(dim, dtype=complex)
    start[0 * dim_w + 1] = 1.0
    out = Ur @ rebit_encode_state(start)
    probs_r = out[:dim] ** 2 + out[dim:] ** 2
    assert np.allclose(probs_r, np.abs(full @ start) ** 2)
    assert np.allclose(
        (probs_r.reshape(dim_c, dim_w)).sum(axis=1), probs
    )                                                 # same answer, one extra qubit


if __name__ == "__main__":
    check_fourier_boundary()
    check_deutsch_jozsa()
    check_simon()
    check_grover()
    check_real_pauli_group()
    check_protocols()
    check_qft_is_complex()
    check_shor_order_finding()
    print("check-real-algorithms: all checks passed")
