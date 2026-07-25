"""Numerical checks for 04-sedenion.md.

Verifies: the Cayley-Dickson sedenions (basis products, 35 triples,
conjugate/norm, two-sided inverses, power associativity, flexibility,
failure of alternativity); the norm defect formula
|xy|^2 = |x|^2|y|^2 + 2 <(d,a,c), b> for x=(a,b), y=(c,d) and its
corollary that doubling an associative algebra preserves norm
multiplicativity; the standard zero divisor (e1+e10)(e5+e14) = 0 and the
failure of cancellation despite the existence of inverses; Moreno's
characterization of zero divisors (components imaginary, equal norm,
orthogonal) and the 4-dimensional annihilator; the correspondence
between basic triples (u1, u2, u3) and zero-divisor pairs
x = (u1, u2)/sqrt(2), y = (u1 u3, u2 u3)/sqrt(2); the diagonal G2 action
(derivations of O extend componentwise to derivations of S), the trivial
stabilizer of a pair, and the 14-dimensional tangent space of the pair
variety (Moreno's theorem: the pair space is homeomorphic to G2);
dim Der(S) = 14 with block-diagonal solutions; and the discrete
automorphisms eps: (a,b) -> (a,-b) and the order-3 psi rotating the
planes (q, q e8) by 2pi/3 (Aut(S) = G2 x S3), neither of which is a
diagonal lift from G2.
"""

from itertools import combinations

import numpy as np

rng = np.random.default_rng(16)


def conj(a):
    c = -a.copy()
    c[0] = a[0]
    return c


def mul(a, b):
    """Cayley-Dickson product: (a1,a2)(b1,b2) = (a1 b1 - b2* a2, b2 a1 + a2 b1*)."""
    n = len(a)
    if n == 1:
        return a * b
    h = n // 2
    a1, a2, b1, b2 = a[:h], a[h:], b[:h], b[h:]
    return np.concatenate([
        mul(a1, b1) - mul(conj(b2), a2),
        mul(b2, a1) + mul(a2, conj(b1)),
    ])


def basis(n, i):
    v = np.zeros(n)
    v[i] = 1.0
    return v


def norm(a):
    return np.sqrt(a @ a)


def im(a):
    b = a.copy()
    b[0] = 0.0
    return b


def rand(n):
    return rng.standard_normal(n)


def rand_unit_im(n=8):
    v = im(rand(n))
    return v / norm(v)


def assoc(a, b, c):
    return mul(mul(a, b), c) - mul(a, mul(b, c))


def pair(a, b):
    return np.concatenate([a, b])


def Lmat(x):
    n = len(x)
    return np.column_stack([mul(x, basis(n, j)) for j in range(n)])


def nullity(M, tol=1e-9):
    return M.shape[1] - np.linalg.matrix_rank(M, tol=tol)


def expm(A, k=20):
    B = A / 2**k
    S, T = np.eye(len(A)), np.eye(len(A))
    for m in range(1, 20):
        T = T @ B / m
        S = S + T
    for _ in range(k):
        S = S @ S
    return S


E8 = [basis(8, i) for i in range(8)]
E16 = [basis(16, i) for i in range(16)]
ONE16 = E16[0]

# ---------------------------------------------------------------- basics
print("== sedenion basics ==")
ok_sq = all(np.allclose(mul(E16[i], E16[i]), -ONE16) for i in range(1, 16))
ok_ac = all(np.allclose(mul(E16[i], E16[j]), -mul(E16[j], E16[i]))
            for i in range(1, 16) for j in range(1, 16) if i != j)
triples = set()
ok_unit = True
for i, j in combinations(range(1, 16), 2):
    p = mul(E16[i], E16[j])
    k = int(np.argmax(np.abs(p)))
    ok_unit &= bool(np.isclose(abs(p[k]), 1) and np.isclose(norm(p), 1) and k not in (0, i, j))
    triples.add(frozenset((i, j, k)))
ok_xor = all(np.allclose(np.abs(mul(E16[i], E16[j])), basis(16, i ^ j))
             for i in range(1, 16) for j in range(1, 16) if i != j)
print("e_i^2 = -1:", ok_sq, "/ anticommute:", ok_ac,
      "/ e_i e_j = ±e_k:", ok_unit, "/ #triples:", len(triples),
      "/ index rule e_i e_j = ±e_(i XOR j):", ok_xor)

x = rand(16)
y = rand(16)
print("(xy)* = y* x*:", np.allclose(conj(mul(x, y)), mul(conj(y), conj(x))))
print("x x* = |x|^2:", np.allclose(mul(x, conj(x)), (x @ x) * ONE16),
      "/ x^-1 x = x x^-1 = 1:",
      np.allclose(mul(conj(x) / (x @ x), x), ONE16)
      and np.allclose(mul(x, conj(x) / (x @ x)), ONE16))

# power associativity (spot check): all bracketings of x^4 agree
x2 = mul(x, x)
p4 = [mul(mul(x2, x), x), mul(x, mul(x2, x)), mul(x2, x2),
      mul(mul(x, x2), x), mul(x, mul(x, x2))]
print("power assoc (x^4 all bracketings):",
      all(np.allclose(p, p4[0]) for p in p4[1:]))

y = rand(16)
flex = np.allclose(mul(mul(x, y), x), mul(x, mul(y, x)))
alt = np.allclose(mul(mul(x, x), y), mul(x, mul(x, y)))
print("flexible (xy)x = x(yx):", flex, "/ alternative (xx)y = x(xy):", alt,
      "(expected False)")

# ------------------------------------------------- norm defect formula
print("\n== norm defect |xy|^2 - |x|^2|y|^2 ==")
# auxiliary identities in the octonions: <u, d*v> = <du, v>, <u, vc*> = <uc, v>
u, v, c, d = rand(8), rand(8), rand(8), rand(8)
print("octonion identities <u,d*v>=<du,v>, <u,vc*>=<uc,v>:",
      np.isclose(u @ mul(conj(d), v), mul(d, u) @ v)
      and np.isclose(u @ mul(v, conj(c)), mul(u, c) @ v))

ok_defect = True
for _ in range(20):
    a, b, c, d = rand(8), rand(8), rand(8), rand(8)
    x, y = pair(a, b), pair(c, d)
    lhs = mul(x, y) @ mul(x, y) - (x @ x) * (y @ y)
    rhs = 2 * (assoc(d, a, c) @ b)
    ok_defect &= bool(np.isclose(lhs, rhs))
print("|xy|^2 = |x|^2|y|^2 + 2<(d,a,c),b> (sedenions):", ok_defect)

# doubling the (associative) quaternions: defect vanishes, norm multiplicative
ok_q = True
for _ in range(20):
    x, y = rand(8), rand(8)
    ok_q &= bool(np.isclose(norm(mul(x, y)), norm(x) * norm(y)))
print("doubling H -> O keeps |xy| = |x||y|:", ok_q)

# ----------------------------------------------- standard zero divisor
print("\n== zero divisors ==")
x0 = E16[1] + E16[10]   # (e1, e2)
y0 = E16[5] + E16[14]   # (e5, e6)
print("(e1+e10)(e5+e14) = 0:", np.allclose(mul(x0, y0), 0),
      "/ (e5+e14)(e1+e10) = 0:", np.allclose(mul(y0, x0), 0))
x0inv = conj(x0) / (x0 @ x0)
print("x0 has inverse:", np.allclose(mul(x0inv, x0), ONE16),
      "/ cancellation fails, x0^-1 (x0 y0) = 0 != y0:",
      np.allclose(mul(x0inv, mul(x0, y0)), 0))

# Moreno's characterization: x = (a,b) is a zero divisor iff
# a, b imaginary, |a| = |b| != 0, <a,b> = 0
ok_zd = True
for _ in range(20):
    u1 = rand_unit_im()
    u2 = im(rand(8))
    u2 -= (u2 @ u1) * u1
    u2 /= norm(u2)
    ok_zd &= nullity(Lmat(pair(u1, u2))) == 4
print("orthonormal imaginary (a,b): dim ker L_x = 4:", ok_zd)

u1, u2 = rand_unit_im(), rand_unit_im()
u2 -= (u2 @ u1) * u1
u2 /= norm(u2)
viol = [
    pair(u1, 1.1 * u2),                    # |a| != |b|
    pair(u1, np.sqrt(0.5) * (u1 + u2)),    # not orthogonal
    pair(u1 + 0.1 * E8[0], u2),            # a not imaginary
    pair(u1, u2 + 0.1 * E8[0]),            # b not imaginary
    rand(16),                              # generic
]
print("violating each condition -> invertible:",
      all(nullity(Lmat(z)) == 0 for z in viol))

# ---------------------------------- annihilator = basic-triple partners
print("\n== annihilator and basic triples ==")
# kernel of L_x is {(u1 u3, u2 u3) : u3 in W}, W = (span 1,u1,u2,u1u2)^perp
u12 = mul(u1, u2)
X = pair(u1, u2)
Wbasis = list(np.linalg.svd(np.array([E8[0], u1, u2, u12]))[2][4:])
K = np.column_stack([pair(mul(u1, w), mul(u2, w)) for w in Wbasis])
print("(u1 u3, u2 u3) with u3 ⊥ 1,u1,u2,u1u2 annihilates x:",
      np.allclose(Lmat(X) @ K, 0),
      "/ spans the 4-dim kernel:",
      np.linalg.matrix_rank(K) == 4 == nullity(Lmat(X)))

# each annihilator element is again a zero divisor of the same kind
u3 = sum(rng.standard_normal() * w for w in Wbasis)
u3 /= norm(u3)
c, d = mul(u1, u3), mul(u2, u3)
print("partner components (u1u3, u2u3) again orthonormal imaginary:",
      np.isclose(c[0], 0) and np.isclose(d[0], 0)
      and np.isclose(norm(c), 1) and np.isclose(norm(d), 1)
      and np.isclose(c @ d, 0))

# basic triple (u1, u2, u3) -> zero-divisor pair, and recovery u3 = -u1 c
Y = pair(c, d)
print("basic triple gives x y = 0 and y x = 0:",
      np.allclose(mul(X, Y), 0) and np.allclose(mul(Y, X), 0),
      "/ u3 recovered as -u1 c:", np.allclose(-mul(u1, c), u3))

# ---------------------------------------------------- G2 acts diagonally
print("\n== G2 action ==")


def leibniz_solutions(n):
    """Basis of the derivation algebra of the dimension-n Cayley-Dickson algebra."""
    En = [basis(n, i) for i in range(n)]
    prod = [[mul(En[i], En[j]) for j in range(n)] for i in range(n)]
    rows = []
    for i in range(1, n):
        for j in range(1, n):
            eq = np.zeros((n, n, n))  # (component, D row k, D col m)
            eij = prod[i][j]
            for r in range(n):
                eq[r, r, :] += eij            # D(e_i e_j)
            for k in range(n):
                eq[:, k, i] -= prod[k][j]     # (D e_i) e_j
                eq[:, k, j] -= prod[i][k]     # e_i (D e_j)
            rows.append(eq.reshape(n, n * n))
    M = np.vstack(rows)
    _, s, vt = np.linalg.svd(M)
    rank = int(np.sum(s > 1e-9 * s[0]))
    return [v.reshape(n, n) for v in vt[rank:]]


DerO = leibniz_solutions(8)
print("dim Der(O) =", len(DerO))

# a derivation of O, extended diagonally, is a derivation of S
D = sum(rng.standard_normal() * B for B in DerO)
Dhat = np.block([[D, np.zeros((8, 8))], [np.zeros((8, 8)), D]])
ok_leib = all(np.allclose(Dhat @ mul(E16[i], E16[j]),
                          mul(Dhat @ E16[i], E16[j]) + mul(E16[i], Dhat @ E16[j]))
              for i in range(16) for j in range(16))
print("diagonal extension satisfies Leibniz on S:", ok_leib)

phi = expm(0.7 * Dhat)
ok_aut = all(np.allclose(phi @ mul(E16[i], E16[j]),
                         mul(phi @ E16[i], phi @ E16[j]))
             for i in range(16) for j in range(16))
print("exp(D^) is an automorphism of S:", ok_aut,
      "/ maps the pair (x0, y0) to a zero-divisor pair:",
      np.allclose(mul(phi @ x0, phi @ y0), 0))

# stabilizer of the pair (x0, y0) is trivial: D e1 = D e2 = D e5 = D e6 = 0 => D = 0
C = np.concatenate([np.column_stack([B @ E8[i] for B in DerO]) for i in (1, 2, 5, 6)])
print("stabilizer dim (D kills e1,e2,e5,e6):", nullity(C))

# tangent space of {(x,y) in S^15 x S^15 : xy = 0} at (x0, y0) has dim 14
X0, Y0 = x0 / np.sqrt(2), y0 / np.sqrt(2)
blockL = np.column_stack([mul(E16[j], Y0) for j in range(16)])
blockR = np.column_stack([mul(X0, E16[j]) for j in range(16)])
Mt = np.vstack([np.hstack([blockL, blockR]),
                np.hstack([X0, np.zeros(16)]),
                np.hstack([np.zeros(16), Y0])])
print("tangent dim of the pair variety:", nullity(Mt), "(= dim G2 = 14)")

# the infinitesimal G2 action D -> (D X0, D Y0) fills that tangent space
orbit = np.column_stack([np.concatenate([B @ X0[:8], B @ X0[8:], B @ Y0[:8], B @ Y0[8:]])
                         for B in DerO])
print("orbit map rank:", np.linalg.matrix_rank(orbit, tol=1e-9),
      "(= 14: injective, so the G2-orbit fills the 14-dim variety)")

# ------------------------------------------- Der(S) and Aut(S) = G2 x S3
print("\n== Der(S) and discrete automorphisms ==")
DerS = leibniz_solutions(16)
diag_ok = all(np.allclose(B[:8, 8:], 0) and np.allclose(B[8:, :8], 0)
              and np.allclose(B[:8, :8], B[8:, 8:]) for B in DerS)
print("dim Der(S) =", len(DerS), "/ all block-diagonal (D, D):", diag_ok)

eps = np.block([[np.eye(8), np.zeros((8, 8))], [np.zeros((8, 8)), -np.eye(8)]])
ok_eps = all(np.allclose(eps @ mul(E16[i], E16[j]),
                         mul(eps @ E16[i], eps @ E16[j]))
             for i in range(16) for j in range(16))
print("eps: (a,b) -> (a,-b) is an automorphism:", ok_eps,
      "/ moves e8 (so not a diagonal lift):", not np.allclose(eps @ E16[8], E16[8]))

# psi: rotate each plane (q, q e8), q in Im O, by 2pi/3; fix 1 and e8
cth, sth = np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3)
psi = np.zeros((16, 16))
psi[0, 0] = 1.0
psi[8, 8] = 1.0
for i in range(1, 8):
    qe = mul(E16[i], E16[8])       # = ±e_{i+8}
    k = int(np.argmax(np.abs(qe)))
    s = qe[k]
    psi[i, i] = cth
    psi[k, i] = sth * s
    psi[i, k] = -sth * s
    psi[k, k] = cth
ok_psi = all(np.allclose(psi @ mul(E16[i], E16[j]),
                         mul(psi @ E16[i], psi @ E16[j]))
             for i in range(16) for j in range(16))
print("psi (rotate (q, q e8) planes by 2pi/3) is an automorphism:", ok_psi,
      "/ psi^3 = 1:", np.allclose(np.linalg.matrix_power(psi, 3), np.eye(16)),
      "/ eps psi eps = psi^2:",
      np.allclose(eps @ psi @ eps, np.linalg.matrix_power(psi, 2)))
print("psi does not preserve the octonion subalgebra (not diagonal):",
      not np.allclose((psi @ E16[1])[8:], 0))
