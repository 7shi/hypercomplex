"""Checks for 04-weyl-algebra.md.

Verifies the article's claims numerically with numpy:

1. Discrete Weyl relation (n = 2..6): Z^b X^a = w^(ab) X^a Z^b with
   w = exp(2 pi i/n), and in group-commutator form
   (Z^b)(X^a)(Z^b)^(-1)(X^a)^(-1) = w^(ab) I.
2. Trace obstruction: tr(AB - BA) = 0 for random A, B, while
   tr(iI) = in != 0, so [A, B] = iI has no solution in M_n(C).
3. Discrete Fourier transform F = (w^(jk)/sqrt(n)): F is unitary,
   F X F^dag = Z, F Z F^dag = X^(-1), F^2 is the parity permutation
   e_j -> e_(-j), and F^4 = I; the k-th column f_k of F is an
   eigenvector of X with eigenvalue w^(-k) = exp(-i eps p_k).
4. Position and momentum: with eps = sqrt(2 pi/n), Q = eps diag(j),
   and P = F diag(p_k) F^dag where p_k = eps k is wrapped to the
   representative in (-n eps/2, n eps/2]; then Z = exp(i eps Q) and
   X = exp(-i eps P) exactly, and the wrap does not change
   exp(-i eps P) (the unwrapped choice gives P = F Q F^dag).
5. Scaled phases: for fixed real s, t and a = round(t/eps),
   b = round(s/eps), the phase w^(ab) = exp(i eps^2 ab) converges to
   exp(ist) as n grows: the error stays within the rounding bound
   (|s| + |t|) eps/2 + eps^2/4 at every n = 16, 64, ..., 4096, which
   itself tends to 0.
6. Exact Weyl relation on the lattice: for representable s = b eps,
   t = a eps, exp(isQ) = Z^b, exp(-itP) = X^a, and the group commutator
   exp(isQ) exp(-itP) exp(-isQ) exp(itP) = exp(ist) I holds to machine
   precision; for generic (non-lattice) s, t it fails by O(1) in
   operator norm, consistent with the trace obstruction (differentiating
   would need arbitrarily small s, t).
7. Canonical commutation relation in the limit: for a unit-width
   Gaussian v centered mid-lattice (away from the seam), the error
   ||[Q, P]v - iv|| decreases rapidly over n = 8, 16, 32, 64, even
   though tr[Q, P] = 0 forces [Q, P] != iI at every finite n.
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

def close(a, b):
    return np.allclose(a, b, atol=1e-12)

def mpow(m, k):
    return np.linalg.matrix_power(m, k)

# --- discrete Weyl relation ------------------------------------------------------
def weyl_discrete_checks(n):
    w = np.exp(2j * np.pi / n)
    X, Z = shift(n), clock(n)
    E = np.eye(n)
    for a in range(n):
        for b in range(n):
            Xa, Zb = mpow(X, a), mpow(Z, b)
            if not close(Zb @ Xa, w ** (a * b) * Xa @ Zb):
                return False
            comm = Zb @ Xa @ np.linalg.inv(Zb) @ np.linalg.inv(Xa)
            if not close(comm, w ** (a * b) * E):
                return False
    return True

print("discrete Weyl relation Z^b X^a = w^(ab) X^a Z^b (n = 2..6):",
      all(weyl_discrete_checks(n) for n in range(2, 7)))

# --- trace obstruction -----------------------------------------------------------
def trace_obstruction_checks():
    rng = np.random.default_rng(0)
    for n in range(2, 7):
        for _ in range(5):
            A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            B = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            if abs(np.trace(A @ B - B @ A)) > 1e-12:
                return False
            if np.isclose(np.trace(1j * np.eye(n)), 0):
                return False
    return True

print("trace obstruction: tr[A, B] = 0 but tr(iI) = in != 0:",
      trace_obstruction_checks())

# --- discrete Fourier transform --------------------------------------------------
def fourier_checks(n):
    X, Z, F = shift(n), clock(n), fourier(n)
    E = np.eye(n)
    w = np.exp(2j * np.pi / n)
    if not close(F @ F.conj().T, E):
        return False
    if not (close(F @ X @ F.conj().T, Z) and
            close(F @ Z @ F.conj().T, np.linalg.inv(X))):
        return False
    for k in range(n):  # columns of F are eigenvectors of X
        if not close(X @ F[:, k], w ** (-k) * F[:, k]):
            return False
    parity = np.eye(n)[:, (-np.arange(n)) % n]  # e_j -> e_{-j}
    return close(mpow(F, 2), parity) and close(mpow(F, 4), E)

print("DFT: F X F^dag = Z, F Z F^dag = X^(-1), X f_k = w^(-k) f_k,",
      "F^2 = parity, F^4 = I:",
      all(fourier_checks(n) for n in range(2, 7)))

# --- position and momentum on the lattice ----------------------------------------
def qp_exponential_checks(n):
    eps, q, p, F, Q, P = lattice_qp(n)
    X, Z = shift(n), clock(n)
    if not close(np.diag(np.exp(1j * eps * q)), Z):
        return False
    expP = F @ np.diag(np.exp(-1j * eps * p)) @ F.conj().T
    if not close(expP, X):
        return False
    # the unwrapped choice P = F Q F^dag gives the same exp(-i eps P)
    expP2 = F @ np.diag(np.exp(-1j * eps * q)) @ F.conj().T
    return close(expP2, X)

print("Z = exp(i eps Q), X = exp(-i eps P), wrap-independent (n = 2..12):",
      all(qp_exponential_checks(n) for n in range(2, 13)))

# --- scaled phases converge to exp(ist) ------------------------------------------
def phase_convergence_checks():
    rng = np.random.default_rng(1)
    for _ in range(5):
        s, t = rng.uniform(-2, 2, size=2)
        for n in [16, 64, 256, 1024, 4096]:
            eps = np.sqrt(2 * np.pi / n)
            a, b = round(t / eps), round(s / eps)
            w = np.exp(2j * np.pi / n)
            err = abs(w ** (a * b) - np.exp(1j * s * t))
            bound = (abs(s) + abs(t)) * eps / 2 + eps ** 2 / 4
            if err > bound:
                return False
        if err > 0.1:  # the bound at n = 4096 is already small
            return False
    return True

print("scaled phases: w^(ab) -> exp(ist) within the rounding bound:",
      phase_convergence_checks())

# --- exact Weyl relation for representable s, t ----------------------------------
def group_commutator(n, s, t):
    eps, q, p, F, Q, P = lattice_qp(n)
    expQ = np.diag(np.exp(1j * s * q))
    expP = F @ np.diag(np.exp(-1j * t * p)) @ F.conj().T
    return expQ @ expP @ expQ.conj().T @ expP.conj().T

def weyl_exact_checks(n):
    eps, q, p, F, Q, P = lattice_qp(n)
    X, Z = shift(n), clock(n)
    for a in range(n):
        for b in range(n):
            s, t = b * eps, a * eps
            if not close(np.diag(np.exp(1j * s * q)), mpow(Z, b)):
                return False
            expP = F @ np.diag(np.exp(-1j * t * p)) @ F.conj().T
            if not close(expP, mpow(X, a)):
                return False
            if not close(group_commutator(n, s, t),
                         np.exp(1j * s * t) * np.eye(n)):
                return False
    return True

def weyl_off_lattice_checks():
    for n in [8, 16, 32]:
        for s, t in [(0.9, 1.1), (np.sqrt(2), 0.5)]:
            dev = np.linalg.norm(group_commutator(n, s, t)
                                 - np.exp(1j * s * t) * np.eye(n), ord=2)
            if dev < 0.5:  # fails by O(1), not by a vanishing error
                return False
    return True

print("exact Weyl relation exp(isQ)exp(-itP) = exp(ist)exp(-itP)exp(isQ)",
      "for s, t on the lattice (n = 2..10):",
      all(weyl_exact_checks(n) for n in range(2, 11)))
print("off-lattice s, t: the relation fails by O(1) in operator norm:",
      weyl_off_lattice_checks())

# --- canonical commutation relation in the limit ---------------------------------
def ccr_convergence_checks():
    errs = []
    for n in [8, 16, 32, 64]:
        eps, q, p, F, Q, P = lattice_qp(n)
        if abs(np.trace(Q @ P - P @ Q)) > 1e-9:
            return False
        v = np.exp(-(q - n * eps / 2) ** 2 / 2)  # away from the seam
        v = v / np.linalg.norm(v)
        errs.append(np.linalg.norm((Q @ P - P @ Q) @ v - 1j * v))
    return all(x > y for x, y in zip(errs, errs[1:])) and errs[-1] < 1e-8

print("CCR in the limit: ||[Q, P]v - iv|| -> 0 for a centered Gaussian,",
      "tr[Q, P] = 0 at every n:", ccr_convergence_checks())
