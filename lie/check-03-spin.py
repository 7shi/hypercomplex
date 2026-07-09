import numpy as np

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

# multiplication table of the representation
print("i^2 = j^2 = k^2 = -1:", all(np.allclose(m@m, -I2) for m in (qi, qj, qk)))
print("ij=k, jk=i, ki=j:", np.allclose(qi@qj, qk) and np.allclose(qj@qk, qi) and np.allclose(qk@qi, qj))

# Pauli matrices from su(2): A = [[ix, y+iz], [-y+iz, -ix]] = i(z s1 + y s2 + x s3)
s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]], dtype=complex)
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
theta = 0.7
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
