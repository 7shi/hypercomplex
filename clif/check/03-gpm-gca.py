"""Checks for 03-gpm-gca.md.

Verifies the article's claims numerically with numpy:

1. Pauli matrices: the Hermitian condition on a general (a, b; c, d)
   forces a, d real and c = conj(b); the trace split t = (a+d)/2,
   z = (a-d)/2 leaves a traceless Hermitian part, and with b = x - iy
   the matrix decomposes over the basis {I, sx, sy, sz}; the Pauli
   matrices square to I, anticommute, and satisfy sx sy sz = iI;
   the opposite convention b = x + iy replaces sy by -sy and flips the
   cyclic products to -i; -i sx, -i sy, -i sz satisfy the quaternion
   relations i^2 = j^2 = k^2 = ijk = -1.
2. Biquaternion: with h = iI and i, j, k = -i sx, -i sy, -i sz, h commutes
   with i, j, k and h^2 = -1; the products hi, hj, hk square to +I and
   reproduce sx, sy, sz; sx sy sz = h; multiplying by -h strips h
   (-h(hi) = i etc.); and the 8 elements
   {1, h, i, hi, j, hj, k, hk} are linearly independent over R, hence a
   real basis of M_2(C) = C (x) H.
3. Clock and shift matrices (n = 2..6): X^n = Z^n = I, ZX = w XZ with
   w = exp(2 pi i/n); the n^2 monomials X^a Z^b are trace-orthogonal and
   traceless except for I, hence a basis of M_n(C); the power formula
   (X^a Z^b)^n = w^(ab n(n-1)/2) I holds.
4. n = 2 reduction: X = sx, Z = sz, i XZ = sy.
5. Nonion correspondence (n = 3, rho = w): u = ZX and v = Z^2 X reproduce
   Sylvester's generators, satisfy vu = rho uv, u^3 = v^3 = 1 and
   det(zI + yu + xv) = z^3 + y^3 + x^3 (random samples); rho^2 uv is
   the shift matrix quoted in the nonion article (= X^2 = X^(-1)); and
   every monomial u^a v^b is a phase multiple of X^(a+b) Z^(a+2b).
6. Generalized Clifford algebra: the extension rule -- attach (x) Z^(n-1)
   to the old generators and append 1 (x) X and 1 (x) c XZ with
   c = exp(-pi i (n-1)/n), c^n (XZ)^n = I -- yields generators satisfying
   e_i^n = 1 and e_j e_i = w e_i e_j (i < j) for n = 2..5, m = 2, 4, 6;
   the n^m monomials e_1^(a_1) ... e_m^(a_m) are linearly independent,
   so the algebra is the full matrix ring M_(n^(m/2))(C).
7. n = 2 case of the extension is the Jordan-Wigner construction (up to
   slot ordering): sigma_z strings with heads (sx, sz) for the first pair
   and (sx, -sy) for the appended pairs (c = -i).
8. Odd numbers of generators: the product e_1 ... e_m exchanges with e_i
   as w^(m+1-2i) e_i (e_1 ... e_m), so it is not central for n > 2; the
   alternating product zeta = e_1 e_2^(n-1) e_3 ... e_m (m odd) is
   central, zeta^n is a root of unity, and e_m is recovered as
   (e_1 e_2^(n-1) ... e_(m-1)^(n-1))^(-1) zeta.
9. Odd classification Cl_m = n M_(n^((m-1)/2))(C): the idempotents
   eps_k = (1/n) sum_j w^(-jk) (zeta/mu)^j (mu^n = zeta^n) are central,
   orthogonal, and sum to 1; the n^m monomials are linearly independent;
   and on each eigenspace of zeta the monomials span the full matrix
   ring of size n^((m-1)/2).
10. m = 1 and Fleury's multicomplex numbers: F = exp(i pi/n) Z satisfies
   F^n = -I (the complexified relation e^n = -1), generates the same
   commutative algebra as Z, and the idempotents eps_k split it into nC.
"""

import numpy as np

def clock(n):
    w = np.exp(2j * np.pi / n)
    return np.diag(w ** np.arange(n))

def shift(n):  # X e_j = e_{j+1 mod n}
    return np.roll(np.eye(n), 1, axis=0).astype(complex)

def close(a, b):
    return np.allclose(a, b, atol=1e-12)

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

# --- Pauli matrices from the general Hermitian matrix --------------------------
def pauli_checks():
    rng = np.random.default_rng(0)
    for _ in range(10):
        # general form with a, d real, c = conj(b) is Hermitian
        a, d = rng.normal(size=2)
        b = rng.normal() + 1j * rng.normal()
        h = np.array([[a, b], [b.conjugate(), d]])
        if not close(h, h.conj().T):
            return False
        # trace split and coefficients from b = x - iy
        t, z = (a + d) / 2, (a - d) / 2
        x, y = b.real, -b.imag
        h0 = h - t * I2
        if not (close(h0, h0.conj().T) and abs(np.trace(h0)) < 1e-12 and
                close(h, t * I2 + x * sx + y * sy + z * sz)):
            return False
    ps = [sx, sy, sz]
    if not all(close(p @ p, I2) for p in ps):
        return False
    if not all(close(ps[i] @ ps[j], -ps[j] @ ps[i])
               for i in range(3) for j in range(3) if i != j):
        return False
    if not close(sx @ sy @ sz, 1j * I2):
        return False
    # opposite convention b = x + iy: sy -> -sy, cyclic products flip to -i
    sy2 = -sy
    if not (close(sx @ sy2, -1j * sz) and close(sy2 @ sz, -1j * sx)
            and close(sz @ sx, -1j * sy2)):
        return False
    qi, qj, qk = -1j * sx, -1j * sy, -1j * sz
    return (all(close(q @ q, -I2) for q in (qi, qj, qk))
            and close(qi @ qj @ qk, -I2) and close(qi @ qj, qk))

print("Pauli basis of Hermitian matrices, relations, quaternions:",
      pauli_checks())

# --- biquaternion C (x) H --------------------------------------------------------
def biquaternion_checks():
    h = 1j * I2
    qi, qj, qk = -1j * sx, -1j * sy, -1j * sz
    # h commutes with i, j, k and h^2 = -1
    if not (close(h @ h, -I2)
            and all(close(h @ q, q @ h) for q in (qi, qj, qk))):
        return False
    # hi, hj, hk square to +I and reproduce the Pauli matrices
    pairs = [(h @ qi, sx), (h @ qj, sy), (h @ qk, sz)]
    if not all(close(a, b) and close(a @ a, I2) for a, b in pairs):
        return False
    # sx sy sz = h: the complex coefficient i is h inside the biquaternion
    if not close(sx @ sy @ sz, h):
        return False
    # multiplying by -h strips h: -h(hi) = i etc.
    if not all(close(-h @ s, q) for s, q in [(sx, qi), (sy, qj), (sz, qk)]):
        return False
    # {1, h, i, hi, j, hj, k, hk} is a real basis of M_2(C)
    basis = [I2, h, qi, h @ qi, qj, h @ qj, qk, h @ qk]
    stack = np.array([np.concatenate([m.flatten().real, m.flatten().imag])
                      for m in basis])
    return np.linalg.matrix_rank(stack) == 8

print("biquaternion: h central, hi/hj/hk = Pauli, real basis of M_2(C):",
      biquaternion_checks())

# --- clock and shift matrices ---------------------------------------------------
def clock_shift_checks(n):
    w = np.exp(2j * np.pi / n)
    X, Z = shift(n), clock(n)
    E = np.eye(n)
    if not (close(np.linalg.matrix_power(X, n), E) and
            close(np.linalg.matrix_power(Z, n), E) and
            close(Z @ X, w * X @ Z)):
        return False
    # monomials X^a Z^b: trace orthogonality and power formula
    mono = {}
    for a in range(n):
        for b in range(n):
            mono[a, b] = (np.linalg.matrix_power(X, a)
                          @ np.linalg.matrix_power(Z, b))
    for (a, b), m1 in mono.items():
        if (a, b) != (0, 0) and abs(np.trace(m1)) > 1e-12:
            return False
        p = np.linalg.matrix_power(m1, n)
        if not close(p, w ** (a * b * n * (n - 1) / 2) * E):
            return False
        for (a2, b2), m2 in mono.items():
            ip = np.trace(m1.conj().T @ m2)
            if not np.isclose(ip, n if (a, b) == (a2, b2) else 0, atol=1e-12):
                return False
    return True

print("clock/shift relations, X^a Z^b basis of M_n(C) (n = 2..6):",
      all(clock_shift_checks(n) for n in range(2, 7)))

print("n = 2: X = sx, Z = sz, i XZ = sy:",
      close(shift(2), sx) and close(clock(2), sz)
      and close(1j * shift(2) @ clock(2), sy))

# --- nonion correspondence (n = 3) ----------------------------------------------
def nonion_checks():
    n = 3
    rho = np.exp(2j * np.pi / n)
    X, Z = shift(n), clock(n)
    u, v = Z @ X, Z @ Z @ X
    u_syl = np.array([[0, 0, 1], [rho, 0, 0], [0, rho ** 2, 0]])
    v_syl = np.array([[0, 0, 1], [rho ** 2, 0, 0], [0, rho, 0]])
    if not (close(u, u_syl) and close(v, v_syl)):
        return False
    if not (close(v @ u, rho * u @ v)
            and close(np.linalg.matrix_power(u, 3), np.eye(3))
            and close(np.linalg.matrix_power(v, 3), np.eye(3))):
        return False
    # det(zI + yu + xv) = z^3 + y^3 + x^3
    rng = np.random.default_rng(1)
    for _ in range(10):
        x, y, z = rng.normal(size=3) + 1j * rng.normal(size=3)
        d = np.linalg.det(z * np.eye(3) + y * u + x * v)
        if not np.isclose(d, z ** 3 + y ** 3 + x ** 3, atol=1e-9):
            return False
    # rho^2 uv = shift matrix of the nonion article = X^2 = X^(-1)
    s = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=complex)
    if not (close(rho ** 2 * u @ v, s) and close(s, X @ X)):
        return False
    # u^a v^b is a phase multiple of X^(a+b) Z^(a+2b)
    for a in range(3):
        for b in range(3):
            m1 = (np.linalg.matrix_power(u, a) @ np.linalg.matrix_power(v, b))
            m2 = (np.linalg.matrix_power(X, (a + b) % 3)
                  @ np.linalg.matrix_power(Z, (a + 2 * b) % 3))
            phase = np.trace(m2.conj().T @ m1) / 3
            if not (np.isclose(abs(phase), 1) and close(m1, phase * m2)):
                return False
    return True

print("nonion: u = ZX, v = Z^2 X, vu = rho uv, det, rho^2 uv = shift:",
      nonion_checks())

# --- generalized Clifford algebra by tensor extension ---------------------------
def gca_generators(n, m):  # m even; returns generators as n^(m/2)-dim matrices
    X, Z = shift(n), clock(n)
    Zt = np.linalg.matrix_power(Z, n - 1)  # = Z^(-1)
    c = np.exp(-1j * np.pi * (n - 1) / n)  # c^n (XZ)^n = I
    gens = [X, Z]
    for _ in range(m // 2 - 1):
        d = len(gens[0])
        gens = [np.kron(e, Zt) for e in gens]
        gens += [np.kron(np.eye(d), X), np.kron(np.eye(d), c * X @ Z)]
    return gens

def gca_relations(n, m):
    w = np.exp(2j * np.pi / n)
    gens = gca_generators(n, m)
    E = np.eye(len(gens[0]))
    if not all(close(np.linalg.matrix_power(e, n), E) for e in gens):
        return False
    return all(close(gens[j] @ gens[i], w * gens[i] @ gens[j])
               for i in range(m) for j in range(i + 1, m))

print("extension generators satisfy e_i^n = 1, e_j e_i = w e_i e_j:",
      all(gca_relations(n, m) for n in range(2, 6) for m in (2, 4, 6)))

def gca_spans(n, m):  # n^m monomials linearly independent => M_(n^(m/2))(C)
    gens = gca_generators(n, m)
    mats = [np.eye(len(gens[0]), dtype=complex)]
    for e in gens:
        mats = [m1 @ np.linalg.matrix_power(e, a) for m1 in mats
                for a in range(n)]
    stack = np.array([m1.flatten() for m1 in mats])
    return np.linalg.matrix_rank(stack) == n ** m

print("monomials span the full matrix ring M_(n^(m/2))(C):",
      all(gca_spans(n, m) for n, m in
          [(2, 2), (2, 4), (2, 6), (3, 2), (3, 4), (4, 4), (5, 4), (3, 6)]))

def jordan_wigner_check():  # n = 2: sz strings; heads (sx, sz), then (sx, -sy)
    gens = gca_generators(2, 6)
    K = 3
    exp = []
    for k in range(K):
        for head in (sx, sz) if k == 0 else (sx, -sy):
            f = np.eye(1, dtype=complex)
            for slot in range(K):
                f = np.kron(f, head if slot == k else (sz if slot > k else I2))
            exp.append(f)
    return all(close(g, e) for g, e in zip(gens, exp))

print("n = 2 extension = Jordan-Wigner strings (sz tails, sx/sz/-sy heads):",
      jordan_wigner_check())

def alternating_product(gens, n):  # zeta = e_1 e_2^(n-1) e_3 ... e_m
    z = np.eye(len(gens[0]), dtype=complex)
    for i, e in enumerate(gens, start=1):
        z = z @ np.linalg.matrix_power(e, 1 if i % 2 else n - 1)
    return z

def odd_generator_checks(n, m):  # m odd
    w = np.exp(2j * np.pi / n)
    gens = gca_generators(n, m + 1)[:m]
    # exchange factor of the plain product e_1 ... e_m with e_i
    pi = np.eye(len(gens[0]), dtype=complex)
    for e in gens:
        pi = pi @ e
    for i, e in enumerate(gens, start=1):
        if not close(pi @ e, w ** (m + 1 - 2 * i) * e @ pi):
            return False
    if n > 2 and m >= 3 and all(close(pi @ e, e @ pi) for e in gens):
        return False  # the plain product must not be central for n > 2
    # the alternating product zeta is central, zeta^n is a root of unity
    z = alternating_product(gens, n)
    if not all(close(z @ e, e @ z) for e in gens):
        return False
    zn = np.linalg.matrix_power(z, n)
    if not (close(zn, zn[0, 0] * np.eye(len(zn)))
            and np.isclose(abs(zn[0, 0]), 1)):
        return False
    # e_m = (e_1 e_2^(n-1) ... e_(m-1)^(n-1))^(-1) zeta
    pre = alternating_product(gens[:-1], n) if m > 1 else np.eye(len(z))
    return close(gens[-1], np.linalg.inv(pre) @ z)

print("odd m: exchange factor w^(m+1-2i), alternating product central:",
      all(odd_generator_checks(n, m) for n in range(2, 6) for m in (1, 3, 5)))

def odd_decomposition_checks(n, m):  # m odd: Cl_m = n M_(n^((m-1)/2))(C)
    w = np.exp(2j * np.pi / n)
    gens = gca_generators(n, m + 1)[:m]
    d = len(gens[0])
    z = alternating_product(gens, n)
    lam = np.linalg.matrix_power(z, n)[0, 0]
    mu = lam ** (1 / n)  # any n-th root of zeta^n
    zh = z / mu  # zh^n = 1: same idempotents as the m = 1 case
    eps = [sum(w ** (-j * k) * np.linalg.matrix_power(zh, j)
               for j in range(n)) / n for k in range(n)]
    if not close(sum(eps), np.eye(d)):
        return False
    for k, p in enumerate(eps):
        for l, q in enumerate(eps):
            if not close(p @ q, p if k == l else np.zeros((d, d))):
                return False
        if not all(close(p @ e, e @ p) for e in gens):
            return False
    # monomials: linearly independent, and the full matrix ring per block
    mats = [np.eye(d, dtype=complex)]
    for e in gens:
        mats = [m1 @ np.linalg.matrix_power(e, a) for m1 in mats
                for a in range(n)]
    stack = np.array([m1.flatten() for m1 in mats])
    if np.linalg.matrix_rank(stack) != n ** m:
        return False
    for p in eps:
        vals, vecs = np.linalg.eigh(p)
        basis = vecs[:, np.isclose(vals, 1)]  # eigenspace, dim n^((m-1)/2)
        if basis.shape[1] != n ** ((m - 1) // 2):
            return False
        blocks = np.array([(basis.conj().T @ m1 @ basis).flatten()
                           for m1 in mats])
        if np.linalg.matrix_rank(blocks) != n ** (m - 1):
            return False
    return True

print("odd m: central idempotents split Cl_m into n M_(n^((m-1)/2))(C):",
      all(odd_decomposition_checks(n, m) for n, m in
          [(2, 1), (2, 3), (2, 5), (3, 1), (3, 3), (4, 3), (5, 3)]))

def fleury_checks(n):  # m = 1: complexified Fleury multicomplex numbers
    w = np.exp(2j * np.pi / n)
    Z = clock(n)
    F = np.exp(1j * np.pi / n) * Z  # e^n = -1
    if not close(np.linalg.matrix_power(F, n), -np.eye(n)):
        return False
    # powers of F span the same commutative algebra as powers of Z
    pf = np.array([np.linalg.matrix_power(F, j).flatten() for j in range(n)])
    pz = np.array([np.linalg.matrix_power(Z, j).flatten() for j in range(n)])
    if not (np.linalg.matrix_rank(pf) == np.linalg.matrix_rank(pz) == n
            and np.linalg.matrix_rank(np.vstack([pf, pz])) == n):
        return False
    # idempotents split the algebra into nC
    eps = [sum(w ** (-j * k) * np.linalg.matrix_power(Z, j)
               for j in range(n)) / n for k in range(n)]
    if not close(sum(eps), np.eye(n)):
        return False
    return all(close(eps[k] @ eps[l],
                     eps[k] if k == l else np.zeros((n, n)))
               for k in range(n) for l in range(n))

print("m = 1: F^n = -I (Fleury), same algebra as C[Z], idempotents -> nC:",
      all(fleury_checks(n) for n in range(2, 7)))
