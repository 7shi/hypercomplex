"""Checks for 05-schrodinger.md.

Verifies the article's claims numerically with numpy:

1. Finite-dimensional Stone's theorem: for random Hermitian H,
   U(t) = exp(-itH) is unitary, satisfies the group law
   U(s+t) = U(s)U(t), conserves the norm, and the derivative
   (U(h) - U(-h))/(2h) at t = 0 recovers -iH.
2. Three forms of the kinetic term agree exactly:
   H0 = (2I - X - X^dag)/(2 eps^2) = F diag(E_k) F^dag with
   E_k = (1 - cos(eps p_k))/eps^2 = D^dag D/2 with the forward
   difference D = (X^dag - I)/eps; E_k does not depend on the
   momentum wrap, and 0 <= p_k^2/2 - E_k <= eps^2 p_k^4/24.
3. Second difference converges to the second derivative: for a
   unit-width Gaussian psi centered mid-circle, the relative error
   ||H0 psi - (-psi''/2)||/||psi|| decreases like eps^2 = 2 pi/n
   (halving per doubling of n) over n = 16, 32, ..., 256.
4. Diagonal matrices are polynomials in Z: V(Q) = sum_b c_b Z^b with
   inverse-DFT coefficients c_b reconstructs a random diagonal matrix.
5. Plane waves are stationary states: exp(-itH0) f_k = exp(-iE_k t) f_k
   for every k and several t.
6. Ehrenfest's theorem on the lattice: harmonic potential
   V = (x - c)^2/2 centered on the circle, displaced unit-width
   Gaussian packet (a coherent state).
   - Momentum side: |i<[H, P]> + <V'(Q)>| is tiny (limited only by the
     seam, like [Q, P]v = iv in check-04) at every sampled time.
   - Position side: |i<[H, Q]> - <P>| = O(eps^2), the lattice
     dispersion correction i[H0, Q] = sin(eps P)/eps = P - eps^2 P^3/6
     + ..., halving per doubling of n over n = 32..256.
   - The packet center <Q> follows the classical trajectory
     c + d cos t over a full period with O(eps^2) deviation.
7. n = 2 (Rabi oscillation): H0 = (I - sigma_x)/pi starting at e_0
   gives occupation |psi_1(t)|^2 = sin^2(t/pi).
"""

import numpy as np

def clock(n):
    w = np.exp(2j * np.pi / n)
    return np.diag(w ** np.arange(n))

def shift(n):  # X e_j = e_{j+1 mod n}
    return np.roll(np.eye(n), 1, axis=0).astype(complex)

def fourier(n):  # F[j, k] = w^(jk)/sqrt(n)
    w = np.exp(2j * np.pi / n)
    jk = np.outer(np.arange(n), np.arange(n))
    return w ** jk / np.sqrt(n)

def lattice_qp(n):  # Q = eps diag(j); P = F diag(p_k) F^dag, p_k wrapped
    eps = np.sqrt(2 * np.pi / n)
    j = np.arange(n)
    q = eps * j
    p = eps * np.where(j <= n // 2, j, j - n)  # wrapped to (-n eps/2, n eps/2]
    F = fourier(n)
    return eps, q, p, F, np.diag(q).astype(complex), F @ np.diag(p) @ F.conj().T

def kinetic(n):  # H0 = (2I - X - X^dag)/(2 eps^2)
    eps = np.sqrt(2 * np.pi / n)
    X = shift(n)
    return (2 * np.eye(n) - X - X.conj().T) / (2 * eps ** 2)

def expmi(H, t):  # exp(-itH) for Hermitian H, via diagonalization
    w, V = np.linalg.eigh(H)
    return (V * np.exp(-1j * t * w)) @ V.conj().T

def close(a, b):
    return np.allclose(a, b, atol=1e-12)

def mpow(m, k):
    return np.linalg.matrix_power(m, k)

# --- finite-dimensional Stone's theorem ------------------------------------------
def stone_checks():
    rng = np.random.default_rng(0)
    for n in [2, 5, 8]:
        A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        H = (A + A.conj().T) / 2
        U = expmi(H, 0.7)
        if not close(U @ U.conj().T, np.eye(n)):
            return False
        for s, t in [(0.3, 0.8), (-1.2, 0.5)]:
            if not close(expmi(H, s + t), expmi(H, s) @ expmi(H, t)):
                return False
        v = rng.normal(size=n) + 1j * rng.normal(size=n)
        if not np.isclose(np.linalg.norm(U @ v), np.linalg.norm(v)):
            return False
        h = 1e-6
        dU = (expmi(H, h) - expmi(H, -h)) / (2 * h)
        if not np.allclose(dU, -1j * H, atol=1e-6):
            return False
    return True

print("Stone (finite dim): exp(-itH) unitary, group law, norm-conserving,",
      "U'(0) = -iH:", stone_checks())

# --- three forms of the kinetic term ----------------------------------------------
def kinetic_form_checks(n):
    eps, q, p, F, Q, P = lattice_qp(n)
    X = shift(n)
    H0 = kinetic(n)
    Ek = (1 - np.cos(eps * p)) / eps ** 2
    Ek_unwrapped = (1 - np.cos(eps * q)) / eps ** 2  # q = eps j = unwrapped p
    if not close(Ek, Ek_unwrapped):
        return False
    if not close(F @ np.diag(Ek) @ F.conj().T, H0):
        return False
    D = (X.conj().T - np.eye(n)) / eps  # forward difference
    if not close(D.conj().T @ D / 2, H0):
        return False
    gap = p ** 2 / 2 - Ek  # 0 <= gap <= eps^2 p^4/24
    return np.all(gap >= -1e-12) and np.all(gap <= eps ** 2 * p ** 4 / 24 + 1e-12)

print("kinetic term: (2I - X - X^dag)/(2 eps^2) = F diag(E_k) F^dag = D^dag D/2,",
      "wrap-free, E_k = p_k^2/2 - O(eps^2 p^4) (n = 2..32):",
      all(kinetic_form_checks(n) for n in range(2, 33)))

# --- second difference -> second derivative ---------------------------------------
def second_difference_convergence():
    errs = []
    for n in [16, 32, 64, 128, 256]:
        eps, q, p, F, Q, P = lattice_qp(n)
        c = n * eps / 2
        psi = np.exp(-(q - c) ** 2 / 2)
        d2 = ((q - c) ** 2 - 1) * psi  # psi''
        err = np.linalg.norm(kinetic(n) @ psi + d2 / 2) / np.linalg.norm(psi)
        errs.append(err)
    ok = all(1.8 < x / y < 2.2 for x, y in zip(errs, errs[1:]))
    return ok and errs[-1] < 1e-2, errs

ok, errs = second_difference_convergence()
print("second difference -> -psi''/2: relative error ~ eps^2 = 2 pi/n",
      "(n = 16..256):", ok)
print("  errors:", ", ".join(f"{e:.2e}" for e in errs))

# --- diagonal matrices are polynomials in Z ---------------------------------------
def diagonal_z_checks(n):
    rng = np.random.default_rng(2)
    w = np.exp(2j * np.pi / n)
    v = rng.normal(size=n)
    b = np.arange(n)
    c = np.array([np.sum(v * w ** (-j * b)) / n for j in b])  # inverse DFT
    Z = clock(n)
    M = sum(c[k] * mpow(Z, k) for k in range(n))
    return close(M, np.diag(v).astype(complex))

print("diagonal V(Q) = sum_b c_b Z^b with inverse-DFT coefficients (n = 2..12):",
      all(diagonal_z_checks(n) for n in range(2, 13)))

# --- plane waves are stationary states --------------------------------------------
def plane_wave_checks(n):
    eps, q, p, F, Q, P = lattice_qp(n)
    H0 = kinetic(n)
    Ek = (1 - np.cos(eps * p)) / eps ** 2
    for t in [0.5, 2.0]:
        U = expmi(H0, t)
        for k in range(n):
            if not close(U @ F[:, k], np.exp(-1j * Ek[k] * t) * F[:, k]):
                return False
    return True

print("plane waves: exp(-itH0) f_k = exp(-iE_k t) f_k (n = 2..8):",
      all(plane_wave_checks(n) for n in range(2, 9)))

# --- Ehrenfest's theorem on the lattice -------------------------------------------
def ehrenfest_run(n, d=1.0):
    eps, q, p, F, Q, P = lattice_qp(n)
    c = n * eps / 2
    Vq = (q - c) ** 2 / 2
    H = kinetic(n) + np.diag(Vq)
    psi0 = np.exp(-(q - c - d) ** 2 / 2).astype(complex)
    psi0 /= np.linalg.norm(psi0)
    w, Vec = np.linalg.eigh(H)
    coef = Vec.conj().T @ psi0
    HQ, HP = H @ Q - Q @ H, H @ P - P @ H
    dV = np.diag(q - c).astype(complex)  # V'(Q)
    errQ = errP = errT = 0.0
    for t in np.linspace(0, 2 * np.pi, 101):
        psi = Vec @ (np.exp(-1j * t * w) * coef)
        expect = lambda A: (psi.conj() @ (A @ psi)).real
        errQ = max(errQ, abs((1j * psi.conj() @ (HQ @ psi)).real - expect(P)))
        errP = max(errP, abs((1j * psi.conj() @ (HP @ psi)).real + expect(dV)))
        errT = max(errT, abs(expect(Q) - (c + d * np.cos(t))))
    return errQ, errP, errT

def ehrenfest_checks():
    runs = [ehrenfest_run(n) for n in [32, 64, 128, 256]]
    errQ, errP, errT = zip(*runs)
    okP = max(errP) < 1e-6  # limited only by the seam
    okQ = all(1.8 < x / y < 2.2 for x, y in zip(errQ, errQ[1:]))  # O(eps^2)
    okT = all(1.8 < x / y < 2.2 for x, y in zip(errT, errT[1:])) and errT[-1] < 0.05
    return okP and okQ and okT, runs

ok, runs = ehrenfest_checks()
print("Ehrenfest (harmonic V, coherent packet): d<P>/dt = -<V'> to seam precision,",
      "d<Q>/dt = <P> + O(eps^2), <Q> tracks c + d cos t (n = 32..256):", ok)
for n, (eQ, eP, eT) in zip([32, 64, 128, 256], runs):
    print(f"  n = {n:3d}: |i<[H,Q]> - <P>| <= {eQ:.2e},",
          f"|i<[H,P]> + <V'>| <= {eP:.2e}, |<Q> - classical| <= {eT:.2e}")

# --- n = 2: Rabi oscillation -------------------------------------------------------
def rabi_checks():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    H = (np.eye(2) - sx) / np.pi  # eps^2 = pi
    e0 = np.array([1, 0], dtype=complex)
    for t in np.linspace(0, 10, 41):
        psi = expmi(H, t) @ e0
        if not np.isclose(abs(psi[1]) ** 2, np.sin(t / np.pi) ** 2, atol=1e-12):
            return False
    return True

print("n = 2: H0 = (I - sigma_x)/pi gives |psi_1(t)|^2 = sin^2(t/pi):",
      rabi_checks())
