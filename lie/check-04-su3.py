"""Numerical checks for 04-su3.md.

Properties of the Gell-Mann matrices (Hermitian, traceless,
normalization tr(la*lb)=2*delta_ab); Lx, Ly, Lz as block placements of
the so(2) generator; the upper-left block embedding of the Pauli
matrices; i*lambda_{2,5,7} matching +-L (so(3) subset su(3)); closure
of the bracket product; exp landing in SU(3)/SO(3); (i*l1)^2 != -I and
exp(i*theta*l1) not expressible as a real linear combination of
{I, i*la} (contrast with the SU(2) degeneracy); the structure
constants f_abc and their total antisymmetry; [l3,l8]=0; the three
su(2) triples (I-, V-, U-spin); the block embeddings SO(2) subset
SO(3) and U(1) subset SU(2) subset SU(3); and the conjugation action
preserving su(3) and giving an 8-dimensional rotation.
"""

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
Lx = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=float)
Ly = np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]], dtype=float)
Lz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
print("Lz, Lx = J in blocks (1,2), (2,3); Ly = -J in block (1,3):",
      np.allclose(Lz[:2, :2], J) and np.allclose(Lx[1:, 1:], J)
      and np.allclose(Ly[np.ix_([0, 2], [0, 2])], -J))

# SO(2) embeds in SO(3) as a block: exp(theta Lz) = diag(R(theta), 1)
theta = 0.7
c, s = np.cos(theta), np.sin(theta)
Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
print("exp(theta Lz) = block diag(SO(2), 1):", np.allclose(expm(Lz*theta).real, Rz))

# lambda_1..3 embed the Pauli matrices in the upper-left block
s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]], dtype=complex)
print("l1,l2,l3 = Pauli in upper-left block:",
      all(np.allclose(l[:2, :2], s) and np.allclose(l[2, :], 0) and np.allclose(l[:, 2], 0)
          for l, s in zip((l1, l2, l3), (s1, s2, s3))))

# so(3) inside su(3): i*lambda_{2,5,7} are the real rotation generators (up to sign)
print("i l2 = -Lz, i l5 = Ly, i l7 = -Lx:",
      np.allclose(1j*l2, -Lz) and np.allclose(1j*l5, Ly) and np.allclose(1j*l7, -Lx))

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
R3 = expm(w[0]*Lx + w[1]*Ly + w[2]*Lz + 0j)
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
