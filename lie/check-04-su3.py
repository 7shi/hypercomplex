"""Numerical checks for 04-su3.md.

Properties of the Gell-Mann matrices (Hermitian, traceless,
normalization tr(la*lb)=2*delta_ab); Jx, Jy, Jz as block placements of
the so(2) generator; the upper-left block embedding of the Pauli
matrices; i*lambda_{2,5,7} matching +-J (so(3) subset su(3)); closure
of the bracket product; exp landing in SU(3)/SO(3); (i*l1)^2 != -I and
exp(i*theta*l1) not expressible as a real linear combination of
{I, i*la} (contrast with the SU(2) degeneracy); the Pauli
representation of quaternions q = aI + b(i*s3) + c(i*s2) + d(i*s1);
the failure of the Clifford relations for the Gell-Mann matrices
(l1^2 != I, {l1,l4} = l6); the C^3 spinor outer product
ww^dagger = I/3 + n_a*la/2 with constant |n|^2 = 4/3; the structure
constants f_abc and their total antisymmetry; [l3,l8]=0; the three
su(2) triples (I-, V-, U-spin); the block embeddings SO(2) subset
SO(3) and U(1) subset SU(2) subset SU(3); the conjugation action
preserving su(3) and giving an 8-dimensional rotation; the clock and
shift matrices generating M_3(C) as a generalized Clifford algebra
(U^3 = V^3 = I, VU = w UV), with Hermitian combinations returning to
Gell-Mann matrices; and SU(3) inside Spin(6): the Cl(6,0) half-spinor
representation realizes spin(6) = su(4), and the stabilizer of a
spinor is an 8-dimensional subalgebra acting on psi^perp as su(3).
"""

from functools import reduce
from itertools import permutations

import numpy as np

# Gell-Mann matrices
l1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
l2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]])
l3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
l4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
l5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]])
l6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
l7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]])
l8 = np.diag([1, 1, -2]).astype(complex) / np.sqrt(3)
lam = [l1, l2, l3, l4, l5, l6, l7, l8]
I3 = np.eye(3, dtype=complex)

def expm(A, terms=40):
    S, T = np.eye(*A.shape, dtype=complex), np.eye(*A.shape, dtype=complex)
    for n in range(1, terms):
        T = T @ A / n
        S = S + T
    return S

def comm(A, B):
    return A @ B - B @ A

# Hermitian, traceless, normalization tr(la lb) = 2 delta_ab
print("lambda_a Hermitian traceless:",
      all(np.allclose(m.conj().T, m) and np.isclose(np.trace(m), 0) for m in lam))
G = np.array([[np.trace(a @ b) for b in lam] for a in lam])
print("tr(la lb) = 2 delta_ab:", np.allclose(G, 2*np.eye(8)))

# so(3) basis: so(2) generator J placed in the three 2x2 blocks (planes 12, 23, 13)
J = np.array([[0, -1], [1, 0]], dtype=float)
Jx = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=float)
Jy = np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]], dtype=float)
Jz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
print("Jz, Jx = J in blocks (1,2), (2,3); Jy = -J in block (1,3):",
      np.allclose(Jz[:2, :2], J) and np.allclose(Jx[1:, 1:], J)
      and np.allclose(Jy[np.ix_([0, 2], [0, 2])], -J))

# SO(2) embeds in SO(3) as a block: exp(theta Jz) = diag(R(theta), 1)
theta = 0.7
c, s = np.cos(theta), np.sin(theta)
Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
print("exp(theta Jz) = block diag(SO(2), 1):", np.allclose(expm(Jz*theta).real, Rz))

# lambda_1..3 embed the Pauli matrices in the upper-left block
s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]], dtype=complex)
print("l1,l2,l3 = Pauli in upper-left block:",
      all(np.allclose(l[:2, :2], s) and np.allclose(l[2, :], 0) and np.allclose(l[:, 2], 0)
          for l, s in zip((l1, l2, l3), (s1, s2, s3))))

# so(3) inside su(3): i*lambda_{2,5,7} are the real rotation generators (up to sign)
print("i l2 = -Jz, i l5 = Jy, i l7 = -Jx:",
      np.allclose(1j*l2, -Jz) and np.allclose(1j*l5, Jy) and np.allclose(1j*l7, -Jx))

# su(3) basis: i*lambda_a anti-Hermitian traceless
print("i la anti-Hermitian traceless:",
      all(np.allclose((1j*m).conj().T, -1j*m) and np.isclose(np.trace(1j*m), 0) for m in lam))

# bracket closes in su(3): [x,y] anti-Hermitian traceless, and lies in span of i*lambda_a
rng = np.random.default_rng(0)
def su3(theta):
    return sum(t * 1j*m for t, m in zip(theta, lam))
x, y = su3(rng.standard_normal(8)), su3(rng.standard_normal(8))
b = comm(x, y)
print("[x,y] anti-Hermitian traceless:",
      np.allclose(b.conj().T, -b) and np.isclose(np.trace(b), 0))
coef = np.array([np.trace(m @ b).imag / 2 for m in lam])  # tr((i la)^dagger b)/2
print("[x,y] in span of i la:", np.allclose(su3(coef), b))

# exp of su(3) element is in SU(3): unitary with det 1
U = expm(x)
print("exp(x): U^dagger U = I, det = 1:",
      np.allclose(U.conj().T @ U, I3) and np.isclose(np.linalg.det(U), 1))

# SU(2) is degenerate (group lives in span{I, i*basis}); SU(3) separates
def span_residual(target, mats):
    A = np.stack([np.concatenate([m.real.ravel(), m.imag.ravel()]) for m in mats], axis=1)
    v = np.concatenate([target.real.ravel(), target.imag.ravel()])
    coef, *_ = np.linalg.lstsq(A, v, rcond=None)
    return np.linalg.norm(A @ coef - v)

print("(i l1)^2 = -(2/3)I - l8/sqrt3, not -I:",
      np.allclose((1j*l1) @ (1j*l1), -2/3*I3 - l8/np.sqrt(3)))
I2 = np.eye(2, dtype=complex)
x2 = sum(t * 1j*m for t, m in zip(rng.standard_normal(3), (s1, s2, s3)))
print("SU(2): exp(x) in real span of {I, i sigma_a}:",
      np.isclose(span_residual(expm(x2), [I2, 1j*s1, 1j*s2, 1j*s3]), 0))
print("SU(3): exp(i theta l1) NOT in real span of {I, i la}:",
      span_residual(expm(1j*theta*l1), [I3] + [1j*m for m in lam]) > 0.1)

# quaternion representation via Pauli matrices: q = aI + b(i s3) + c(i s2) + d(i s1)
def qrep(a, b, c, d):
    return a*I2 + b*1j*s3 + c*1j*s2 + d*1j*s1
def qmul(p, q):
    (a, b, c, d), (e, f, g, h) = p, q
    return (a*e - b*f - c*g - d*h, a*f + b*e + c*h - d*g,
            a*g - b*h + c*e + d*f, a*h + b*g - c*f + d*e)
p4, q4 = rng.standard_normal(4), rng.standard_normal(4)
print("aI + b(i s3) + c(i s2) + d(i s1) reproduces the quaternion product:",
      np.allclose(qrep(*p4) @ qrep(*q4), qrep(*qmul(p4, q4))))

# Gell-Mann matrices are not Clifford generators (unlike the Pauli matrices)
print("l1^2 = diag(1,1,0) != I:",
      np.allclose(l1 @ l1, np.diag([1, 1, 0])) and not np.allclose(l1 @ l1, I3))
print("{l1, l4} = l6 != 0:", np.allclose(l1 @ l4 + l4 @ l1, l6))

# spinor outer product in C^3: ww^dagger = I/3 + (1/2) sum n_a la, |n|^2 = 4/3
w3 = rng.standard_normal(3) + 1j*rng.standard_normal(3)
w3 /= np.linalg.norm(w3)
n = np.array([(w3.conj() @ m @ w3).real for m in lam])
print("ww^dagger = I/3 + n_a la / 2:",
      np.allclose(np.outer(w3, w3.conj()), I3/3 + sum(v*m for v, m in zip(n, lam))/2))
print("|n|^2 = 4/3 (constant):", np.isclose(n @ n, 4/3))

# structure constants: [la, lb] = 2i f_abc lc, f_abc = tr([la,lb] lc)/(4i)
f = np.zeros((8, 8, 8))
for a in range(8):
    for b_ in range(8):
        for c in range(8):
            f[a, b_, c] = (np.trace(comm(lam[a], lam[b_]) @ lam[c]) / 4j).real
known = {(1,2,3): 1, (1,4,7): .5, (1,5,6): -.5, (2,4,6): .5, (2,5,7): .5,
         (3,4,5): .5, (3,6,7): -.5, (4,5,8): np.sqrt(3)/2, (6,7,8): np.sqrt(3)/2}
print("f_abc matches known values:",
      all(np.isclose(f[a-1, b-1, c-1], v) for (a, b, c), v in known.items()))
print("f_abc totally antisymmetric:",
      np.allclose(f, -f.transpose(1, 0, 2)) and np.allclose(f, -f.transpose(0, 2, 1)))
recon = all(np.allclose(comm(lam[a], lam[b_]), 2j*sum(f[a, b_, c]*lam[c] for c in range(8)))
            for a in range(8) for b_ in range(8))
print("[la,lb] = 2i sum f_abc lc:", recon)
# no other independent nonzero components
mask = np.zeros((8, 8, 8), dtype=bool)
for idx in known:
    for p in permutations(idx):
        mask[p[0]-1, p[1]-1, p[2]-1] = True
print("f_abc = 0 outside listed (anti)symmetrizations:", np.allclose(f[~mask], 0))

# rank 2: lambda_3 and lambda_8 commute and are independent
print("[l3, l8] = 0:", np.allclose(comm(l3, l8), 0))

# [l4, l5] = i (l3 + sqrt3 l8)
print("[l4,l5] = i(l3 + sqrt3 l8):", np.allclose(comm(l4, l5), 1j*(l3 + np.sqrt(3)*l8)))

# three su(2) triples with Pauli relations [a,b]=2ic, [b,c]=2ia, [c,a]=2ib
V3 = (l3 + np.sqrt(3)*l8) / 2
U3 = (-l3 + np.sqrt(3)*l8) / 2
triples = [(l1, l2, l3), (l4, l5, V3), (l6, l7, U3)]
ok = all(np.allclose(comm(a, b), 2j*c) and np.allclose(comm(b, c), 2j*a)
         and np.allclose(comm(c, a), 2j*b) for a, b, c in triples)
print("su(2) triples (I, V, U spin):", ok)
print("V3 = diag(1,0,-1):", np.allclose(V3, np.diag([1, 0, -1])))

# SO(3) inside SU(3): exp of a real antisymmetric combo is real orthogonal with det 1
w = rng.standard_normal(3)
R3 = expm(w[0]*Jx + w[1]*Jy + w[2]*Jz + 0j)
print("exp(so(3) combo) real, orthogonal, det 1 (SO(3) in SU(3)):",
      np.allclose(R3.imag, 0) and np.allclose(R3.T @ R3, I3) and np.isclose(np.linalg.det(R3), 1))

# embeddings: U(1) in SU(2) in SU(3)
print("exp(i theta l3) = diag(e^it, e^-it, 1):",
      np.allclose(expm(1j*theta*l3), np.diag([np.exp(1j*theta), np.exp(-1j*theta), 1])))
A = expm(su3([*rng.standard_normal(3), 0, 0, 0, 0, 0]))
print("exp of l1..l3 combo = block diag(SU(2), 1):",
      np.isclose(A[2, 2], 1) and np.allclose(A[2, :2], 0) and np.allclose(A[:2, 2], 0)
      and np.allclose(A[:2, :2].conj().T @ A[:2, :2], np.eye(2))
      and np.isclose(np.linalg.det(A[:2, :2]), 1))

# conjugation x -> U x U^-1 preserves su(3) and acts as an 8-dim rotation
conj_ok = all(np.allclose((U @ (1j*m) @ U.conj().T).conj().T, -(U @ (1j*m) @ U.conj().T))
              and np.isclose(np.trace(U @ (1j*m) @ U.conj().T), 0) for m in lam)
print("U (i la) U^-1 stays in su(3):", conj_ok)
Ad = np.array([[np.trace(lam[a] @ U @ lam[b] @ U.conj().T).real / 2
                for b in range(8)] for a in range(8)])
print("Ad(U) as 8x8: orthogonal, det 1:",
      np.allclose(Ad.T @ Ad, np.eye(8)) and np.isclose(np.linalg.det(Ad), 1))

# clock and shift matrices: generalized Clifford algebra generating M_3(C)
omega = np.exp(2j*np.pi/3)
Vc = np.diag([1, omega, omega**2])
Us = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)
print("U^3 = V^3 = I, VU = w UV:",
      np.allclose(np.linalg.matrix_power(Us, 3), I3)
      and np.allclose(np.linalg.matrix_power(Vc, 3), I3)
      and np.allclose(Vc @ Us, omega * Us @ Vc))
print("U, V in SU(3), unitary but not Hermitian:",
      np.isclose(np.linalg.det(Us), 1) and np.isclose(np.linalg.det(Vc), 1)
      and np.allclose(Us.conj().T @ Us, I3) and np.allclose(Vc.conj().T @ Vc, I3)
      and not np.allclose(Us, Us.conj().T) and not np.allclose(Vc, Vc.conj().T))
mono = [np.linalg.matrix_power(Us, a) @ np.linalg.matrix_power(Vc, b)
        for a in range(3) for b in range(3)]
Gm = np.array([[np.trace(A_.conj().T @ B_) for B_ in mono] for A_ in mono])
print("tr((U^a V^b)^dagger U^c V^d) = 3 delta (basis of M_3(C)):",
      np.allclose(Gm, 3*np.eye(9)))
print("U^a V^b traceless except I:", all(np.isclose(np.trace(M), 0) for M in mono[1:]))
# 2x2 analogue: clock = s3, shift = s1, anticommuting (ordinary Clifford)
print("2x2: clock = s3, shift = s1, s3 s1 = -s1 s3:", np.allclose(s3 @ s1, -s1 @ s3))
# Hermitian combinations return to Gell-Mann matrices
print("U + U^2 = l1 + l4 + l6, i(U - U^2) = l2 - l5 + l7:",
      np.allclose(Us + Us @ Us, l1 + l4 + l6)
      and np.allclose(1j*(Us - Us @ Us), l2 - l5 + l7))
print("V + V^2, i(V - V^2) in span of l3, l8 (Cartan):",
      np.allclose(Vc + Vc @ Vc, 1.5*l3 + np.sqrt(3)/2*l8)
      and np.allclose(1j*(Vc - Vc @ Vc), np.sqrt(3)/2*l3 - 1.5*l8))

# SU(3) inside Spin(6): Cl(6,0) via 8x8 gamma matrices
def kron(*ms):
    return reduce(np.kron, ms)

I2c = np.eye(2, dtype=complex)
gam = [kron(s1, I2c, I2c), kron(s2, I2c, I2c), kron(s3, s1, I2c),
       kron(s3, s2, I2c), kron(s3, s3, s1), kron(s3, s3, s2)]
print("Cl(6,0) relations (6 anticommuting, square +I):",
      all(np.allclose(gam[i] @ gam[j] + gam[j] @ gam[i], 2*np.eye(8)*(i == j))
          for i in range(6) for j in range(6)))
g7 = -1j * reduce(lambda a_, b_: a_ @ b_, gam)
ev, evec = np.linalg.eigh(g7)
Wp = evec[:, np.isclose(ev, 1)]
print("chirality g7 Hermitian, g7^2 = I, half-spinor space = C^4:",
      np.allclose(g7, g7.conj().T) and np.allclose(g7 @ g7, np.eye(8))
      and Wp.shape[1] == 4)
biv = [gam[i] @ gam[j] / 2 for i in range(6) for j in range(i+1, 6)]
blk = [Wp.conj().T @ B_ @ Wp for B_ in biv]
Ab = np.stack([np.concatenate([B_.real.ravel(), B_.imag.ravel()]) for B_ in blk])
print("spin(6) on half-spinors: anti-Hermitian traceless, rank 15 = dim su(4):",
      all(np.allclose(B_.conj().T, -B_) and np.isclose(np.trace(B_), 0) for B_ in blk)
      and np.linalg.matrix_rank(Ab) == 15)
# stabilizer of a spinor psi: {X in spin(6) | X psi = 0} is su(3)
psi = rng.standard_normal(4) + 1j*rng.standard_normal(4)
psi /= np.linalg.norm(psi)
Mc = np.stack([B_ @ psi for B_ in blk]).T
_, sv, Vt = np.linalg.svd(np.concatenate([Mc.real, Mc.imag]))
stab = [sum(c_*B_ for c_, B_ in zip(v, blk)) for v in Vt[np.sum(sv > 1e-10):]]
print("stabilizer of psi: dimension 8 = dim su(3):", len(stab) == 8)
print("stabilizer closed under bracket:",
      all(np.isclose(span_residual(comm(X, Y), stab), 0) for X in stab for Y in stab))
Qm, _ = np.linalg.qr(np.column_stack([psi, np.eye(4, dtype=complex)[:, :3]]))
perp = Qm[:, 1:]
restr = [perp.conj().T @ X @ perp for X in stab]
Ar = np.stack([np.concatenate([B_.real.ravel(), B_.imag.ravel()]) for B_ in restr])
print("restriction to psi^perp: anti-Hermitian traceless on C^3, rank 8 -> su(3):",
      all(np.allclose(B_.conj().T, -B_) and np.isclose(np.trace(B_), 0) for B_ in restr)
      and np.linalg.matrix_rank(Ar) == 8)
