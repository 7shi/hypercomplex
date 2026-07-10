"""Numerical checks for 03-spin.md.

Spinor sign flip after one full rotation (theta=2pi) and recovery after
two; the outer product w*w^dagger transforming under the conjugation
action with global phase cancelling; w*w^dagger = (1/2)(I + n1*s1 +
n2*s2 + n3*s3) with (n1,n2,n3) a unit vector matching spherical
coordinates; the conjugation action giving the 3D rotation of
(n1,n2,n3), agreeing for +-U; the general form of su(2) as i times a
Hermitian matrix decomposed into Pauli matrices; the correspondence
i<->i*s3, j<->i*s2, k<->i*s1; Pauli matrix norm/orthogonality
(sa^2=I, anticommuting); the grade-2 factorization i=s1*s2, j=s3*s1,
k=s2*s3 and (s1*s2)^2=-I; the even subalgebra closing under the
product and excluding odd grades; the two-sided action p*x*q^-1 giving
per-plane rotation angles alpha -+ beta; general p, q giving a 4x4
orthogonal matrix with det=1; (-p,-q) giving the same rotation (double
cover); the left/right generators commuting and totaling 6 dimensions
(so(4) = sp(1) + sp(1)); the 4x4 matrix representations of left/right
action (Lu*Lv=Luv, Ru*Rv=Rvu, left-right commutativity) and the 16
products Lu*Rv spanning M4(R); p*x*q^-1 = Lp*R(q^-1); the doubled
gamma(v) satisfying the Cl(4,0) generating relations with 16
independent basis elements; the conjugation action
g*gamma(x)*g^-1=gamma(p*x*q^-1) and the split into half-spinors
(psi1,psi2) -> (p*psi1, q*psi2); the signs of the squares of the basis
Lu*Rv (6 one-sided give -I, 9 two-sided give +I); that mutually
anticommuting sets have size at most 5 (by exhaustive search); the
generators (LiRi, LjRi, LkRi, Rj) forming Cl(3,1) = M4(R) with volume
element -Rk; the classification of grade-2/3 products (including the
Cl(2,2)-type presentation); the projection P=(1/4)(I-LiRi-LjRj-LkRk)
extracting the real part (average of the conjugation action) so that
Lx*P has x in its first column; and the reconstruction formula
Lx = sum_u Ru*(Lx*P)*Ru^-1.
"""

import numpy as np

# quaternions as 2x2 complex matrices (same as the article): q = z + wj <-> [[z, w], [-conj(w), conj(z)]]
I2 = np.eye(2, dtype=complex)
qi = np.array([[1j, 0], [0, -1j]])
qj = np.array([[0, 1], [-1, 0]], dtype=complex)
qk = np.array([[0, 1j], [1j, 0]])

# Pauli matrices
s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]], dtype=complex)

def Q(a, b, c, d):
    return a*I2 + b*qi + c*qj + d*qk

def coords(x):  # matrix -> (a, b, c, d)
    return np.array([x[0,0].real, x[0,0].imag, x[0,1].real, x[0,1].imag])

def expm(A, terms=40):
    S, T = np.eye(*A.shape, dtype=complex), np.eye(*A.shape, dtype=complex)
    for n in range(1, terms):
        T = T @ A / n
        S = S + T
    return S

def comm(A, B):
    return A @ B - B @ A

rng = np.random.default_rng(0)

# multiplication table of the representation
print("i^2 = j^2 = k^2 = -1:", all(np.allclose(m@m, -I2) for m in (qi, qj, qk)))
print("ij=k, jk=i, ki=j:", np.allclose(qi@qj, qk) and np.allclose(qj@qk, qi) and np.allclose(qk@qi, qj))

# spinor: one-sided action with half angle
theta = 0.7
U = expm(qi*theta/2)
print("exp(i theta/2) = diag(e^{i theta/2}, e^{-i theta/2}):",
      np.allclose(U, np.diag([np.exp(1j*theta/2), np.exp(-1j*theta/2)])))
w = rng.standard_normal(2) + 1j*rng.standard_normal(2)
w = w / np.linalg.norm(w)
print("one full rotation (theta = 2 pi) flips the spinor sign:",
      np.allclose(expm(qi*2*np.pi/2) @ w, -w))
print("two full rotations (theta = 4 pi) restore the spinor:",
      np.allclose(expm(qi*4*np.pi/2) @ w, w))

# outer product: transforms by conjugation, phase cancels
W = np.outer(w, w.conj())
print("w w^dag is Hermitian with trace 1:",
      np.allclose(W, W.conj().T) and np.isclose(np.trace(W), 1))
print("(Uw)(Uw)^dag = U (w w^dag) U^-1:",
      np.allclose(np.outer(U@w, (U@w).conj()), U @ W @ U.conj().T))
w_ph = np.exp(1j*0.3) * w
print("overall phase e^{i phi} w cancels in the outer product:",
      np.allclose(np.outer(w_ph, w_ph.conj()), W))

# Pauli decomposition: w w^dag = (I + n1 s1 + n2 s2 + n3 s3)/2 with (n1, n2, n3) a unit vector
al, be = w
n3 = abs(al)**2 - abs(be)**2
n12 = 2*al.conjugate()*be  # n1 + i n2
n1, n2 = n12.real, n12.imag
print("w w^dag = (I + n1 s1 + n2 s2 + n3 s3)/2:",
      np.allclose(W, (I2 + n1*s1 + n2*s2 + n3*s3)/2))
print("n1^2 + n2^2 + n3^2 = 1:", np.isclose(n1**2 + n2**2 + n3**2, 1))
t0, p0 = 1.1, 2.3  # spherical parametrization a = cos(t/2), b = sin(t/2) e^{i p}
ws = np.array([np.cos(t0/2), np.sin(t0/2)*np.exp(1j*p0)])
Ws = np.outer(ws, ws.conj())
ns = (np.sin(t0)*np.cos(p0), np.sin(t0)*np.sin(p0), np.cos(t0))
print("spherical spinor -> (sin t cos p, sin t sin p, cos t):",
      np.allclose(Ws, (I2 + ns[0]*s1 + ns[1]*s2 + ns[2]*s3)/2))
print("det(n1 s1 + n2 s2 + n3 s3) = -(n1^2 + n2^2 + n3^2):",
      np.isclose(np.linalg.det(n1*s1 + n2*s2 + n3*s3), -(n1**2 + n2**2 + n3**2)))

# conjugation rotates (n1, n2, n3): full angle theta about the sigma_3 axis
al2, be2 = U @ w
n12_new = 2*al2.conjugate()*be2
print("U action: n3 fixed, n1 + i n2 -> e^{-i theta}(n1 + i n2):",
      np.isclose(abs(al2)**2 - abs(be2)**2, n3) and np.allclose(n12_new, np.exp(-1j*theta)*n12))
print("U s2 U^-1 = cos(theta) s2 + sin(theta) s1 (matches R_q(j) = j cos + k sin):",
      np.allclose(U @ s2 @ U.conj().T, np.cos(theta)*s2 + np.sin(theta)*s1))

def rot3(Um):  # 3x3 matrix of (n1, n2, n3) under conjugation by Um
    cols = [[(np.trace((Um @ s @ Um.conj().T) @ t) / 2).real for t in (s1, s2, s3)]
            for s in (s1, s2, s3)]
    return np.array(cols).T
uv = rng.standard_normal(4)
Uq = Q(*uv) / np.sqrt(sum(uv**2))  # random unit quaternion
R3 = rot3(Uq)
print("conjugation by a unit quaternion: orthogonal, det 1 (3D rotation):",
      np.allclose(R3.T @ R3, np.eye(3)) and np.isclose(np.linalg.det(R3), 1))
print("U and -U give the same rotation (double cover):", np.allclose(rot3(-Uq), R3))

# rem: quantum-mechanics-style derivation from su(2): A = [[ix, y+iz], [-y+iz, -ix]] = i(z s1 + y s2 + x s3)
x, y, z = rng.standard_normal(3)
A = np.array([[1j*x, y + 1j*z], [-y + 1j*z, -1j*x]])
H = np.array([[x, z - 1j*y], [z + 1j*y, -x]])
print("A = xi + yj + zk:", np.allclose(A, x*qi + y*qj + z*qk))
print("A = iH, H Hermitian traceless:",
      np.allclose(A, 1j*H) and np.allclose(H.conj().T, H) and np.isclose(np.trace(H), 0))
print("H = z s1 + y s2 + x s3:", np.allclose(H, z*s1 + y*s2 + x*s3))
print("i = i s3, j = i s2, k = i s1:",
      np.allclose(qi, 1j*s3) and np.allclose(qj, 1j*s2) and np.allclose(qk, 1j*s1))

# grade 1: unit length and mutual orthogonality (anticommutation)
print("s_a^2 = I:", all(np.allclose(s@s, I2) for s in (s1, s2, s3)))
print("s_a s_b + s_b s_a = 0 (a != b):",
      all(np.allclose(a@b, -b@a) for a, b in ((s1, s2), (s2, s3), (s3, s1))))

# grade 2: products of two Pauli matrices are the quaternion units
print("s1 s2 = i s3, s2 s3 = i s1, s3 s1 = i s2:",
      np.allclose(s1@s2, 1j*s3) and np.allclose(s2@s3, 1j*s1) and np.allclose(s3@s1, 1j*s2))
print("i = s1 s2, j = s3 s1, k = s2 s3:",
      np.allclose(qi, s1@s2) and np.allclose(qj, s3@s1) and np.allclose(qk, s2@s3))
print("(s1 s2)^2 = -I:", np.allclose((s1@s2) @ (s1@s2), -I2))

# even subalgebra {I, s1s2, s3s1, s2s3} is closed and isomorphic to H
even = [I2, s1@s2, s3@s1, s2@s3]
def in_span(m, basis):  # real-linear span (grading is over the reals)
    B = np.column_stack([np.concatenate([b.flatten().real, b.flatten().imag]) for b in basis])
    v = np.concatenate([m.flatten().real, m.flatten().imag])
    c, res, *_ = np.linalg.lstsq(B, v, rcond=None)
    return np.allclose(B @ c, v)
print("even subalgebra closed under products:",
      all(in_span(a@b, even) for a in even for b in even))
print("odd grades not in even subalgebra:", not any(in_span(s, even) for s in (s1, s2, s3)))

# conjugation on grade 1 gives the same rotation as on grade 2 (i passes through)
q = expm(qi*theta/2)
print("q s_a q^-1 = -i q (i s_a) q^-1 stays grade 1 (Hermitian traceless):",
      all(np.allclose((q@s@q.conj().T).conj().T, q@s@q.conj().T) and
          np.isclose(np.trace(q@s@q.conj().T), 0) for s in (s1, s2, s3)))

# two-sided action p x q^-1 with p = exp(i alpha), q = exp(i beta):
# (1,i)-plane rotates by alpha - beta, (j,k)-plane by alpha + beta
alpha, beta = 0.9, 0.4
p, q = expm(qi*alpha), expm(qi*beta)
a, b, c, d = rng.standard_normal(4)
X = Q(a, b, c, d)
Y = p @ X @ q.conj().T  # q^-1 = q^dagger for unit q
zc, wc = a + 1j*b, c + 1j*d
z2, w2 = np.exp(1j*(alpha-beta))*zc, np.exp(1j*(alpha+beta))*wc
print("p x q^-1 = e^{i(a-b)} z + e^{i(a+b)} w j:",
      np.allclose(Y, Q(z2.real, z2.imag, w2.real, w2.imag)))
print("beta = alpha fixes (1,i)-plane (conjugation):",
      np.allclose(expm(qi*alpha) @ I2 @ expm(qi*alpha).conj().T, I2))

# general p, q: x -> p x q^-1 is a 4D rotation (SO(4))
pv, qv = rng.standard_normal(4), rng.standard_normal(4)
P = Q(*pv) / np.sqrt(sum(pv**2))
R = Q(*qv) / np.sqrt(sum(qv**2))
basis = [I2, qi, qj, qk]
M4 = np.column_stack([coords(P @ e @ R.conj().T) for e in basis])
print("4x4 matrix of x -> p x q^-1: orthogonal, det 1:",
      np.allclose(M4.T @ M4, np.eye(4)) and np.isclose(np.linalg.det(M4), 1))
N4 = np.column_stack([coords((-P) @ e @ (-R).conj().T) for e in basis])
print("(-p, -q) gives the same rotation:", np.allclose(M4, N4))

# so(4) = sp(1) + sp(1): left and right multiplications commute and span a 6-dim space
def L(u):  # 4x4 matrix of x -> u x
    return np.column_stack([coords(u @ e) for e in basis])
def Rm(u):  # 4x4 matrix of x -> x u^dagger = -x u for pure unit u
    return np.column_stack([coords(e @ u.conj().T) for e in basis])
lefts = [L(u) for u in (qi, qj, qk)]
rights = [Rm(u) for u in (qi, qj, qk)]
print("left/right generators are antisymmetric (in so(4)):",
      all(np.allclose(m.T, -m) for m in lefts + rights))
print("[L_u, R_v] = 0 (associativity):",
      all(np.allclose(comm(a, b), 0) for a in lefts for b in rights))
print("L_i, L_j, L_k close like sp(1) ([L_i, L_j] = 2 L_k etc.):",
      np.allclose(comm(lefts[0], lefts[1]), 2*lefts[2]) and
      np.allclose(comm(lefts[1], lefts[2]), 2*lefts[0]) and
      np.allclose(comm(lefts[2], lefts[0]), 2*lefts[1]))
S = np.column_stack([m.flatten() for m in lefts + rights])
print("dim(span(L) + span(R)) = 6 = dim so(4):", np.linalg.matrix_rank(S) == 6)

# 4x4 real matrices of left/right multiplication: L_u x = ux, R_u x = xu (no inverse)
def Rr(u):  # 4x4 matrix of x -> x u
    return np.column_stack([coords(e @ u) for e in basis])
I4, Z4 = np.eye(4), np.zeros((4, 4))
Li, Lj, Lk = L(qi), L(qj), L(qk)
Ri, Rj, Rk = Rr(qi), Rr(qj), Rr(qk)
Li_art = np.array([[0,-1,0,0], [1,0,0,0], [0,0,0,-1], [0,0,1,0]], dtype=float)
Ri_art = np.array([[0,-1,0,0], [1,0,0,0], [0,0,0,1], [0,0,-1,0]], dtype=float)
print("L_i and R_i match the signed permutation matrices in the article:",
      np.allclose(Li, Li_art) and np.allclose(Ri, Ri_art))
print("L_u^2 = R_u^2 = -I:",
      all(np.allclose(m@m, -I4) for m in (Li, Lj, Lk, Ri, Rj, Rk)))
print("anticommutation within the L family and within the R family:",
      all(np.allclose(a@b, -b@a) for a, b in
          ((Li,Lj), (Lj,Lk), (Lk,Li), (Ri,Rj), (Rj,Rk), (Rk,Ri))))
print("L_i L_j = L_k (composition follows the product):", np.allclose(Li@Lj, Lk))
print("R_j R_i = R_k (right action reverses the order):", np.allclose(Rj@Ri, Rk))
print("L_u R_v = R_v L_u (associativity):",
      all(np.allclose(a@b, b@a) for a in (Li,Lj,Lk) for b in (Ri,Rj,Rk)))
print("left action alone closes at 4 dims (products stay in {I, L_i, L_j, L_k}):",
      np.allclose(Li@Lj, Lk) and
      np.linalg.matrix_rank(np.column_stack([m.flatten() for m in (I4, Li, Lj, Lk)])) == 4)
prods = [a@b for a in (I4, Li, Lj, Lk) for b in (I4, Ri, Rj, Rk)]
print("16 products L_u R_v are linearly independent (span M_4(R)):",
      np.linalg.matrix_rank(np.column_stack([m.flatten() for m in prods])) == 16)
print("p x q^{-1} = L_p R_{q^-1} as matrices:", np.allclose(M4, L(P) @ Rr(R.conj().T)))

# neutral-signature Clifford algebra: (L_i, L_j, L_kR_i, L_kR_j) generate Cl_{2,2} = M_4(R)
def subset_products(gens):
    mats = []
    for bits in range(2**len(gens)):
        m = np.eye(*gens[0].shape)
        for n in range(len(gens)):
            if bits >> n & 1:
                m = m @ gens[n]
        mats.append(m.flatten())
    return np.column_stack(mats)
# enumerate squares and anticommutation among the basis L_u R_v -> Cl_{3,1} = M_4(R)
from itertools import combinations
singles = [Li, Lj, Lk, Ri, Rj, Rk]
doubles = [a@b for a in (Li, Lj, Lk) for b in (Ri, Rj, Rk)]
print("squares: one-sided 6 give -I, two-sided 9 give +I:",
      all(np.allclose(m@m, -I4) for m in singles) and
      all(np.allclose(m@m, I4) for m in doubles))
def anti(a, b):
    return np.allclose(a@b, -b@a)
elems = singles + doubles  # the 15 non-identity basis elements
adj = [[anti(a, b) for b in elems] for a in elems]
def clique(idxs):
    return all(adj[a][b] for a, b in combinations(idxs, 2))
print("mutually anticommuting sets: some of size 5, none of size 6:",
      any(clique(c) for c in combinations(range(15), 5)) and
      not any(clique(c) for c in combinations(range(15), 6)))
e1, e2, e3, e4 = Li@Ri, Lj@Ri, Lk@Ri, Rj
five = [e1, e2, e3, e4, Rk]
print("{L_iR_i, L_jR_i, L_kR_i, R_j, R_k}: anticommuting, squares (+,+,+,-,-):",
      all(anti(a, b) for a, b in combinations(five, 2)) and
      all(np.allclose(m@m, s*I4) for m, s in zip(five, (1, 1, 1, -1, -1))))
print("volume element e1 e2 e3 e4 = -R_k (the 5th element):",
      np.allclose(e1@e2@e3@e4, -Rk))
print("2^4 = 16 products are linearly independent (Cl_{3,1} = M_4(R)):",
      np.linalg.matrix_rank(subset_products([e1, e2, e3, e4])) == 16)
grades = [(e1@e2, -Lk), (e1@e3, Lj), (e2@e3, -Li),                # grade 2: L_a
          (e1@e4, -Li@Rk), (e2@e4, -Lj@Rk), (e3@e4, -Lk@Rk),     # grade 2: L_a R_k
          (e1@e2@e3, Ri),                                         # grade 3: R_i
          (e1@e2@e4, -Lk@Rj), (e1@e3@e4, Lj@Rj), (e2@e3@e4, -Li@Rj)]  # grade 3: L_a R_j
print("grades 2 and 3 recover {L_a}, {L_a R_k}, {L_a R_j}, R_i (up to sign):",
      all(np.allclose(a, b) for a, b in grades))
g22 = [e1, e2, Rj, Rk]
print("choosing (L_iR_i, L_jR_i, R_j, R_k) instead: signature (+,+,-,-), Cl_{2,2}:",
      all(anti(a, b) for a, b in combinations(g22, 2)) and
      all(np.allclose(m@m, s*I4) for m, s in zip(g22, (1, 1, -1, -1))) and
      np.linalg.matrix_rank(subset_products(g22)) == 16)

# projection P: spinor in matrix form (left ideal), and restoration of the operator
P4 = np.diag([1., 0, 0, 0])
print("P = (I - L_iR_i - L_jR_j - L_kR_k)/4:",
      np.allclose(P4, (I4 - Li@Ri - Lj@Rj - Lk@Rk)/4))
LX = L(X)
print("L_x P has x in the first column and 0 elsewhere:",
      np.allclose((LX @ P4)[:, 0], coords(X)) and np.allclose((LX @ P4)[:, 1:], 0))
print("P extracts the real part = average over conjugations by 1, i, j, k:",
      np.allclose(P4 @ coords(X), [coords(X)[0], 0, 0, 0]) and
      np.allclose(sum(coords(u @ X @ np.linalg.inv(u)) for u in basis)/4,
                  [coords(X)[0], 0, 0, 0]))
print("restoration L_x = sum_u R_u (L_x P) R_u^{-1}:",
      np.allclose(LX, sum(Ru @ (LX @ P4) @ np.linalg.inv(Ru) for Ru in (I4, Ri, Rj, Rk))))

# doubling to H + H: gamma(v) = [[0, L_v], [L_conj(v), 0]] generates Cl_{4,0}
def gamma(v):
    return np.block([[Z4, L(v)], [L(v.conj().T), Z4]])
gammas = [gamma(m) for m in basis]
I8 = np.eye(8)
print("gamma(u)^2 = +I and the four generators anticommute (Cl_{4,0}):",
      all(np.allclose(g@g, I8) for g in gammas) and
      all(np.allclose(a@b, -b@a) for n, a in enumerate(gammas) for b in gammas[n+1:]))
v4 = rng.standard_normal(4)
V = Q(*v4)
print("gamma(v)^2 = |v|^2 I:", np.allclose(gamma(V) @ gamma(V), sum(v4**2)*I8))
mats = []
for bits in range(16):  # products over all subsets of the 4 generators
    m = I8.copy()
    for n in range(4):
        if bits >> n & 1:
            m = m @ gammas[n]
    mats.append(m.flatten())
print("2^4 = 16 products of generators are linearly independent (faithful Cl_{4,0}):",
      np.linalg.matrix_rank(np.column_stack(mats)) == 16)
print("gamma(u)gamma(v) = diag(L_{u conj(v)}, L_{conj(u) v}):",
      np.allclose(gammas[1] @ gammas[2],
                  np.block([[L(qi @ qj.conj().T), Z4], [Z4, L(qi.conj().T @ qj)]])))
print("volume element gamma(1)gamma(i)gamma(j)gamma(k) = diag(-I, I):",
      np.allclose(gammas[0] @ gammas[1] @ gammas[2] @ gammas[3],
                  np.block([[-I4, Z4], [Z4, I4]])))

# Spin(4) = diag(L_p, L_q): conjugation on grade 1 gives p x q^{-1}, spinors split in half
g8 = np.block([[L(P), Z4], [Z4, L(R)]])
g8inv = np.block([[L(P.conj().T), Z4], [Z4, L(R.conj().T)]])
print("g gamma(x) g^{-1} = gamma(p x q^{-1}):",
      np.allclose(g8 @ g8inv, I8) and
      np.allclose(g8 @ gamma(X) @ g8inv, gamma(P @ X @ R.conj().T)))
psi1, psi2 = Q(*rng.standard_normal(4)), Q(*rng.standard_normal(4))
psi = np.concatenate([coords(psi1), coords(psi2)])
print("spinor action splits into half-spinors: (psi1, psi2) -> (p psi1, q psi2):",
      np.allclose(g8 @ psi, np.concatenate([coords(P @ psi1), coords(R @ psi2)])))
g8m = np.block([[L(-P), Z4], [Z4, L(-R)]])
print("(-p, -q): same action on vectors, flips the spinor:",
      np.allclose(g8m @ gamma(X) @ np.linalg.inv(g8m), gamma(P @ X @ R.conj().T)) and
      np.allclose(g8m @ psi, -(g8 @ psi)))
