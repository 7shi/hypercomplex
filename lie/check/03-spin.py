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
product and excluding odd grades.
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
