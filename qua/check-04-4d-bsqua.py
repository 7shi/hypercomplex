"""Checks for 04-4d-bsqua.md.

Verifies the split-biquaternion correspondence for 4D quaternion rotation:

1. Basis pairing alpha = e1e2e3e4 beta for i:(e1e2,e4e3), j:(e1e3,e2e4),
   k:(e1e4,e3e2).
2. Properties of omega = e1e2e3e4: omega^2 = 1, commutes with all bivectors;
   the images i = e4e3, j = e2e4, k = e3e2 satisfy the quaternion relations.
3. The idempotents (1 +- omega)/2 are orthogonal and central in the even
   subalgebra, and each summand is a quaternion algebra (so the
   split-biquaternion is isomorphic to 2H).
4. T is a homomorphism on the even subalgebra: T(xy) = T(x)T(y).
5. For random unit vectors m, n (r = mn, r_L = nm*, r_R = m*n):
   - the quaternion rotation r_L q r_R matches the Clifford rotation
     (nm)v(mn) under 1 <-> e1, i <-> e2, j <-> e3, k <-> e4;
   - r_R = T(r), r_L* = T(r^dagger) (the two projections of the even
     subalgebra 2H give the right and left rotors up to conjugation), and,
     writing r = c + B (scalar + bivector), r_R = c + T(B) and
     r_L = c - T(B^dagger);
   - proof lemmas: T(e1 v) = q pairs vectors with quaternions, and
     e1 x e1 = x^dagger on the even subalgebra (conjugation by e1 is
     the omega-conjugation).
6. The L/R multiplication matrices: L_u R_v (16 products) span M4(R);
   e1 = L_iR_i, e2 = L_jR_i, e3 = L_kR_i, e4 = R_j generate Cl_{3,1}
   (squares +,+,+,-, pairwise anticommuting); omega = e1e2e3e4 = -R_k;
   grade samples e1e2 = -L_k, e1e2e3 = R_i.
7. The even subalgebra of Cl_{4,0} as Cl_{0,3}: u_a = e_1 e_{a+1} square
   to -1 and anticommute; u2u3 = e4e3 (=i), u3u1 = e2e4 (=j),
   u1u2 = e3e2 (=k), u1u2u3 = -omega; and omega*i = e1e2 = u1 etc., so
   grade 1 of Cl_{0,3} carries omega*i, omega*j, omega*k.
8. Double rotations: omega swaps paired bivectors, so the sums/differences
   e1e2 +- e4e3 are invariant/sign-flipped (self-dual/anti-self-dual); the
   normal-form rotation with angles (t1, t2) in the e1e2/e4e3 planes equals
   exp((t1-t2)/2 i) q exp((t1+t2)/2 i); and for a general bivector B,
   exp(-B) v exp(B) matches r_L q r_R with r_R = exp(T(B)),
   r_L = exp(-T(B^dagger)).
9. L_p R_{q*} for unit quaternions p, q is an SO(4) matrix (orthogonal,
   determinant 1), and (p, q), (-p, -q) give the same matrix: the tensor
   product SU(2) (x) SU(2) = SO(4) collapses the sign of the pair.
   L_p R_{p*} = I only for p = +-1, so the kernel of the 2-to-1 map
   (r_L, r_R) -> L_{r_L} R_{r_R} is {+-(1,1)}.
10. The even-subalgebra ladder: in Cl_{3,0} and Cl_{0,3}, e1e3 and e2e3
    square to -1 and anticommute, generating Cl_{0,2} = H; in Cl_{0,4},
    u_a = e1 e_{a+1} square to -1 and anticommute (Cl_{0,3}) and
    omega^2 = +1, matching Cl_{4,0}.
"""

import math
import random

random.seed(0)

# --- Clifford algebra Cl_{4,0}: multivector = {sorted index tuple: coeff} ---

def bmul(a, b, sq=None):
    """Product of two basis elements (sorted tuples) -> (sign, basis).

    sq maps index -> square of that generator; None means all +1 (Cl_{4,0}).
    """
    s, t = 1, list(a)
    for x in b:
        i = len(t)
        while i > 0 and t[i - 1] > x:
            i -= 1
            s = -s
        if i > 0 and t[i - 1] == x:
            t.pop(i - 1)
            if sq is not None:
                s *= sq[x]
        else:
            t.insert(i, x)
    return s, tuple(t)

def mul(x, y, sq=None):
    z = {}
    for a, ca in x.items():
        for b, cb in y.items():
            s, t = bmul(a, b, sq)
            z[t] = z.get(t, 0) + s * ca * cb
    return {k: v for k, v in z.items() if abs(v) > 1e-12}

def add(x, y):
    z = dict(x)
    for k, v in y.items():
        z[k] = z.get(k, 0) + v
    return {k: v for k, v in z.items() if abs(v) > 1e-12}

def smul(c, x):
    return {k: c * v for k, v in x.items()}

def sub(x, y):
    return add(x, smul(-1, y))

def eq(x, y, tol=1e-9):
    return all(abs(v) < tol for v in sub(x, y).values())

E = lambda *idx: {tuple(idx): 1}
ONE = E()
OMEGA = E(1, 2, 3, 4)

# --- quaternions: 4-tuples (w, x, y, z) -------------------------------------

def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return (w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
            w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2)

def qconj(a):
    return (a[0], -a[1], -a[2], -a[3])

def qeq(a, b, tol=1e-9):
    return all(abs(s - t) < tol for s, t in zip(a, b))

# --- check 1: basis pairing alpha = omega beta ------------------------------

pairs = [(E(1, 2), {(3, 4): -1}),   # i: (e1e2, e4e3)
         (E(1, 3), E(2, 4)),        # j: (e1e3, e2e4)
         (E(1, 4), {(2, 3): -1})]   # k: (e1e4, e3e2)
for alpha, beta in pairs:
    assert eq(alpha, mul(OMEGA, beta))

# --- check 2: omega properties and quaternion relations ---------------------

I, J, K = ({(3, 4): -1}, E(2, 4), {(2, 3): -1})  # e4e3, e2e4, e3e2

assert eq(mul(OMEGA, OMEGA), ONE)
bivectors = [E(i, j) for i in range(1, 5) for j in range(i + 1, 5)]
for b in bivectors:
    assert eq(mul(OMEGA, b), mul(b, OMEGA))

for u in (I, J, K):
    assert eq(mul(u, u), smul(-1, ONE))
assert eq(mul(I, J), K) and eq(mul(J, K), I) and eq(mul(K, I), J)

# --- check 3: idempotents and the 2H decomposition --------------------------

ep = smul(0.5, add(ONE, OMEGA))
em = smul(0.5, sub(ONE, OMEGA))
assert eq(mul(ep, ep), ep) and eq(mul(em, em), em)
assert eq(mul(ep, em), {}) and eq(add(ep, em), ONE)
for e in (ep, em):  # each summand e*H is a quaternion algebra with unit e
    ei, ej, ek = mul(e, I), mul(e, J), mul(e, K)
    assert eq(mul(ei, ei), smul(-1, e))
    assert eq(mul(ei, ej), ek) and eq(mul(ej, ek), ei) and eq(mul(ek, ei), ej)

# --- T: even subalgebra -> H, and the omega-conjugation dagger --------------
# In the article T acts on split-biquaternions; the script's T is the
# composite T(phi(x)), with the basis isomorphism phi baked into T_TABLE.

T_TABLE = {(): (1, 0, 0, 0), (1, 2, 3, 4): (1, 0, 0, 0),      # 1, omega -> 1
           (1, 2): (0, 1, 0, 0), (3, 4): (0, -1, 0, 0),       # omega*i, i
           (1, 3): (0, 0, 1, 0), (2, 4): (0, 0, 1, 0),        # omega*j, j
           (1, 4): (0, 0, 0, 1), (2, 3): (0, 0, 0, -1)}       # omega*k, k

def T(x):
    q = (0, 0, 0, 0)
    for b, c in x.items():
        q = tuple(s + c * t for s, t in zip(q, T_TABLE[b]))
    return q

def dagger(x):  # negate the omega-part: bases containing e1, and omega itself
    return {b: (-c if 1 in b else c) for b, c in x.items()}

# --- check 4: T is a homomorphism --------------------------------------------

even_basis = list(T_TABLE)
for _ in range(100):
    x = {b: random.uniform(-1, 1) for b in even_basis}
    y = {b: random.uniform(-1, 1) for b in even_basis}
    assert qeq(T(mul(x, y)), qmul(T(x), T(y)))

# --- check 5: rotation r_L q r_R matches (nm)v(mn), and rotor formulas ------

def rand_unit():
    v = [random.gauss(0, 1) for _ in range(4)]
    nrm = sum(c * c for c in v) ** 0.5
    return [c / nrm for c in v]

def vec(c):  # coefficients -> Clifford vector
    return {(i + 1,): c[i] for i in range(4)}

for _ in range(100):
    mc, nc = rand_unit(), rand_unit()
    m, n = vec(mc), vec(nc)
    mq, nq = tuple(mc), tuple(nc)

    # Clifford rotation (nm)v(mn) vs quaternion rotation (nm*)q(m*n)
    r = mul(m, n)
    rl, rr = qmul(nq, qconj(mq)), qmul(qconj(mq), nq)
    v = rand_unit()
    got = qmul(qmul(rl, tuple(v)), rr)
    want = mul(mul(mul(n, m), vec(v)), r)
    assert all(abs(want.get((i + 1,), 0) - got[i]) < 1e-9 for i in range(4))

    # r_R = T(r); r_L* = T(r^dagger); with r = c + B: r_R = c + T(B),
    # r_L = c - T(B^dagger)
    assert qeq(rr, T(r))
    assert qeq(qconj(rl), T(dagger(r)))
    c = r.get((), 0)
    B = {b: v2 for b, v2 in r.items() if len(b) == 2}
    assert qeq(rr, tuple(s + t for s, t in zip((c, 0, 0, 0), T(B))))
    assert qeq(rl, tuple(s - t for s, t in zip((c, 0, 0, 0), T(dagger(B)))))

# proof lemmas: T(e1 v) = q pairs vectors with quaternions, and conjugation
# by e1 realizes the omega-conjugation dagger on the even subalgebra
for _ in range(20):
    v = rand_unit()
    assert qeq(T(mul(E(1), vec(v))), tuple(v))
    x = {b: random.uniform(-1, 1) for b in even_basis}
    assert eq(mul(mul(E(1), x), E(1)), dagger(x))

# --- check 6: L/R matrices, M4(R), and the Cl_{3,1} generators --------------

import numpy as np

BASIS4 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]

def Lmat(u):
    return np.array([qmul(u, b) for b in BASIS4]).T

def Rmat(u):
    return np.array([qmul(b, u) for b in BASIS4]).T

I4 = np.eye(4)
Li, Lj, Lk = (Lmat(b) for b in BASIS4[1:])
Ri, Rj, Rk = (Rmat(b) for b in BASIS4[1:])

# product rules and commutativity of left and right actions
assert np.allclose(Li @ Lj, Lk) and np.allclose(Rj @ Ri, Rk)
for A in (Li, Lj, Lk):
    for B in (Ri, Rj, Rk):
        assert np.allclose(A @ B, B @ A)

# the 16 products L_u R_v are linearly independent and span M4(R)
prods = [Lmat(u) @ Rmat(v) for u in BASIS4 for v in BASIS4]
assert np.linalg.matrix_rank(np.array([P.flatten() for P in prods])) == 16

# generators of Cl_{3,1}: squares (+,+,+,-), pairwise anticommuting
gens = [Li @ Ri, Lj @ Ri, Lk @ Ri, Rj]
for g, s in zip(gens, (1, 1, 1, -1)):
    assert np.allclose(g @ g, s * I4)
for a in range(4):
    for b in range(a + 1, 4):
        assert np.allclose(gens[a] @ gens[b], -gens[b] @ gens[a])

# omega = e1e2e3e4 = -R_k, omega^2 = -1; grade samples
omega31 = gens[0] @ gens[1] @ gens[2] @ gens[3]
assert np.allclose(omega31, -Rk)
assert np.allclose(omega31 @ omega31, -I4)
assert np.allclose(gens[0] @ gens[1], -Lk)          # e1e2 = -L_k (grade 2)
assert np.allclose(gens[0] @ gens[1] @ gens[2], Ri)  # e1e2e3 = R_i (grade 3)

# --- check 7: Cl^0_{4,0} = Cl_{0,3}, grade assignment ------------------------

u123 = [E(1, 2), E(1, 3), E(1, 4)]
for a in range(3):
    assert eq(mul(u123[a], u123[a]), smul(-1, ONE))
    for b in range(a + 1, 3):
        ab, ba = mul(u123[a], u123[b]), mul(u123[b], u123[a])
        assert eq(ab, smul(-1, ba))
assert eq(mul(u123[1], u123[2]), I)  # u2u3 = e4e3
assert eq(mul(u123[2], u123[0]), J)  # u3u1 = e2e4
assert eq(mul(u123[0], u123[1]), K)  # u1u2 = e3e2
assert eq(mul(mul(u123[0], u123[1]), u123[2]), smul(-1, OMEGA))
for u, q in zip(u123, (I, J, K)):    # omega*i = e1e2 = u1 etc. (grade 1)
    assert eq(mul(OMEGA, q), u)

# --- check 8: double rotations and the self-dual decomposition ---------------

# omega swaps paired bivectors: e1e2 +- e4e3 are invariant / sign-flipped
E43 = {(3, 4): -1}  # e4e3
assert eq(mul(OMEGA, E(1, 2)), E43)
assert eq(mul(OMEGA, add(E(1, 2), E43)), add(E(1, 2), E43))
assert eq(mul(OMEGA, sub(E(1, 2), E43)), smul(-1, sub(E(1, 2), E43)))

def cl_exp(x, terms=40):
    z, term = dict(ONE), dict(ONE)
    for n in range(1, terms):
        term = smul(1 / n, mul(term, x))
        z = add(z, term)
    return z

def qexp_im(v):  # exp of a pure-imaginary quaternion
    t = (v[1] * v[1] + v[2] * v[2] + v[3] * v[3]) ** 0.5
    if t < 1e-15:
        return (1.0, 0.0, 0.0, 0.0)
    s = math.sin(t) / t
    return (math.cos(t), s * v[1], s * v[2], s * v[3])

BIVECTORS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
for _ in range(50):
    # normal form: angles t1 in the e1e2 plane, t2 in the e4e3 plane
    t1, t2 = random.uniform(-3, 3), random.uniform(-3, 3)
    B = add(smul(t1 / 2, E(1, 2)), smul(t2 / 2, E43))
    r, rinv = cl_exp(B), cl_exp(smul(-1, B))
    rl = qexp_im((0, (t1 - t2) / 2, 0, 0))
    rr = qexp_im((0, (t1 + t2) / 2, 0, 0))
    assert qeq(T(r), rr) and qeq(qconj(rl), T(dagger(r)))
    v = rand_unit()
    want = mul(mul(rinv, vec(v)), r)
    got = qmul(qmul(rl, tuple(v)), rr)
    assert all(abs(want.get((a + 1,), 0) - got[a]) < 1e-7 for a in range(4))

    # general rotor exp(B): r_R = exp(T(B)), r_L = exp(-T(B^dagger))
    B = {b: random.uniform(-1, 1) for b in BIVECTORS}
    r, rinv = cl_exp(B), cl_exp(smul(-1, B))
    rr = qexp_im(T(B))
    rl = qexp_im(tuple(-c for c in T(dagger(B))))
    assert qeq(T(r), rr) and qeq(qconj(rl), T(dagger(r)))
    v = rand_unit()
    want = mul(mul(rinv, vec(v)), r)
    got = qmul(qmul(rl, tuple(v)), rr)
    assert all(abs(want.get((a + 1,), 0) - got[a]) < 1e-7 for a in range(4))

# --- check 9: L_p R_{q*} is an SO(4) matrix, and the sign collapses ----------

for _ in range(20):
    p, q = tuple(rand_unit()), tuple(rand_unit())
    M = Lmat(p) @ Rmat(qconj(q))
    assert np.allclose(M.T @ M, I4) and abs(np.linalg.det(M) - 1) < 1e-9
    pn, qn = tuple(-c for c in p), tuple(-c for c in q)
    assert np.allclose(Lmat(pn) @ Rmat(qconj(qn)), M)

# kernel: L_p R_{p*} = I only for p = +-1 (p* = p^-1 for a unit quaternion)
for _ in range(20):
    p = tuple(rand_unit())
    if max(abs(c) for c in p[1:]) > 1e-3:
        assert not np.allclose(Lmat(p) @ Rmat(qconj(p)), I4)

# --- check 10: the even-subalgebra ladder in other signatures ----------------

for sq in ({a: 1 for a in (1, 2, 3)}, {a: -1 for a in (1, 2, 3)}):
    v1, v2 = E(1, 3), E(2, 3)  # Cl^0_{3,0} = Cl^0_{0,3} = Cl_{0,2} = H
    assert eq(mul(v1, v1, sq), smul(-1, ONE))
    assert eq(mul(v2, v2, sq), smul(-1, ONE))
    assert eq(mul(v1, v2, sq), smul(-1, mul(v2, v1, sq)))
    v3 = mul(v1, v2, sq)  # the third quaternion unit
    assert eq(mul(v3, v3, sq), smul(-1, ONE))

sq04 = {a: -1 for a in (1, 2, 3, 4)}  # Cl_{0,4}: same ladder as Cl_{4,0}
us = [E(1, a) for a in (2, 3, 4)]
for a in range(3):
    assert eq(mul(us[a], us[a], sq04), smul(-1, ONE))
    for b in range(a + 1, 3):
        assert eq(mul(us[a], us[b], sq04), smul(-1, mul(us[b], us[a], sq04)))
assert eq(mul(OMEGA, OMEGA, sq04), ONE)

print("All checks passed.")
