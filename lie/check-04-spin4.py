"""Numerical checks for 04-spin4.md.

The two-sided action p*x*q^-1 giving per-plane rotation angles alpha -+
beta; general p, q giving a 4x4 orthogonal matrix with det=1; (-p,-q)
giving the same rotation (double cover); the left/right generators
commuting and totaling 6 dimensions (so(4) = sp(1) + sp(1)); the 4x4
matrix representations of left/right action (Lu*Lv=Luv, Ru*Rv=Rvu,
left-right commutativity) and the 16 products Lu*Rv spanning M4(R);
p*x*q^-1 = Lp*R(q^-1); the signs of the squares of the basis Lu*Rv (6
one-sided give -I, 9 two-sided give +I); that mutually anticommuting
sets have size at most 5 (by exhaustive search); the generators (LiRi,
LjRi, LkRi, Rj) forming Cl(3,1) = M4(R) with volume element -Rk; the
classification of grade-2/3 products (including the Cl(2,2)-type
presentation); the projection P=(1/4)(I-LiRi-LjRj-LkRk) extracting the
real part (average of the conjugation action) so that Lx*P has x in
its first column; the reconstruction formula Lx = sum_u Ru*(Lx*P)*Ru^-1;
the doubled gamma(v) satisfying the Cl(4,0) generating relations with
16 independent basis elements; the conjugation action
g*gamma(x)*g^-1=gamma(p*x*q^-1) and the split into half-spinors
(psi1,psi2) -> (p*psi1, q*psi2).
"""

import numpy as np
from itertools import combinations

# quaternions as 2x2 complex matrices (same as the article): q = z + wj <-> [[z, w], [-conj(w), conj(z)]]
I2 = np.eye(2, dtype=complex)
qi = np.array([[1j, 0], [0, -1j]])
qj = np.array([[0, 1], [-1, 0]], dtype=complex)
qk = np.array([[0, 1j], [1j, 0]])

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
