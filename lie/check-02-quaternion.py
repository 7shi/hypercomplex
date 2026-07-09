import numpy as np

# quaternions as 2x2 complex matrices: q = z + wj <-> [[z, -w], [conj(w), conj(z)]]
I2 = np.eye(2, dtype=complex)
qi = np.array([[1j, 0], [0, -1j]])
qj = np.array([[0, -1], [1, 0]], dtype=complex)
qk = np.array([[0, -1j], [-1j, 0]])

def Q(a, b, c, d):
    return a*I2 + b*qi + c*qj + d*qk

def expm(A, terms=40):
    S, T = np.eye(*A.shape, dtype=complex), np.eye(*A.shape, dtype=complex)
    for n in range(1, terms):
        T = T @ A / n
        S = S + T
    return S

def comm(A, B):
    return A @ B - B @ A

# multiplication table
print("i^2 = j^2 = k^2 = -1:", all(np.allclose(m@m, -I2) for m in (qi, qj, qk)))
print("ij=k, jk=i, ki=j:", np.allclose(qi@qj, qk) and np.allclose(qj@qk, qi) and np.allclose(qk@qi, qj))
print("ji=-k, kj=-i, ik=-j:", np.allclose(qj@qi, -qk) and np.allclose(qk@qj, -qi) and np.allclose(qi@qk, -qj))

# Sylvester conditions: traceless, det 1, mutually anticommuting (see derivation rem)
print("i,j,k traceless, det 1:",
      all(np.isclose(np.trace(m), 0) and np.isclose(np.linalg.det(m), 1) for m in (qi, qj, qk)))
print("i,j,k mutually anticommute:",
      all(np.allclose(a@b, -b@a) for a, b in ((qi, qj), (qj, qk), (qk, qi))))

# det = |q|^2 and multiplicativity of the norm
rng = np.random.default_rng(0)
pa, qa = rng.standard_normal(4), rng.standard_normal(4)
P, M = Q(*pa), Q(*qa)
print("det Q = |q|^2:", np.isclose(np.linalg.det(M).real, sum(qa**2)))
print("|pq| = |p||q|:", np.isclose(np.linalg.det(P@M).real, sum(pa**2)*sum(qa**2)))

# Euler formula: exp(u*theta) = cos(theta) + u sin(theta) for unit pure u
v = rng.standard_normal(3)
u = Q(0, *v) / np.linalg.norm(v)
theta = 0.7
print("u^2 = -1:", np.allclose(u@u, -I2))
print("exp(u t) = cos t + u sin t:", np.allclose(expm(u*theta), np.cos(theta)*I2 + np.sin(theta)*u))

# failure of exp(x)exp(y) = exp(x+y): exp(i pi/2) exp(j pi/2) = k != exp((i+j) pi/2)
lhs = expm(qi*np.pi/2) @ expm(qj*np.pi/2)
rhs = expm((qi+qj)*np.pi/2)
print("exp(i pi/2) exp(j pi/2) = k:", np.allclose(lhs, qk))
print("exp((i+j) pi/2) != k:", not np.allclose(rhs, qk))
c = np.cos(np.sqrt(2)*np.pi/2); s = np.sin(np.sqrt(2)*np.pi/2)
print("exp((i+j) pi/2) = cos + (i+j)/sqrt2 sin:", np.allclose(rhs, c*I2 + s*(qi+qj)/np.sqrt(2)))

# brackets: [i,j]=2k cyclic, and [p,q] = 2 (p x q) for pure p, q
print("[i,j]=2k, [j,k]=2i, [k,i]=2j:",
      np.allclose(comm(qi,qj), 2*qk) and np.allclose(comm(qj,qk), 2*qi) and np.allclose(comm(qk,qi), 2*qj))
pv, qv = rng.standard_normal(3), rng.standard_normal(3)
print("pq = -p.q + p x q (pure):",
      np.allclose(Q(0,*pv) @ Q(0,*qv), Q(-np.dot(pv,qv), *np.cross(pv,qv))))
print("[p,q] = 2 p x q (pure):",
      np.allclose(comm(Q(0,*pv), Q(0,*qv)), 2*Q(0,*np.cross(pv,qv))))

# SU(2): unit quaternion -> unitary with det 1; basis is anti-Hermitian traceless
U = M / np.sqrt(np.linalg.det(M).real)
print("unit q: A^dagger A = I, det = 1:",
      np.allclose(U.conj().T @ U, I2) and np.isclose(np.linalg.det(U).real, 1))
print("i,j,k anti-Hermitian traceless:",
      all(np.allclose(m.conj().T, -m) and np.isclose(np.trace(m), 0) for m in (qi, qj, qk)))

# Pauli matrices from the general form of su(2):
# traceless anti-Hermitian A = [[iz, y+ix], [-y+ix, -iz]] = i(x s1 + y s2 + z s3)
s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]], dtype=complex)
x, y, z = rng.standard_normal(3)
A = np.array([[1j*z, y + 1j*x], [-y + 1j*x, -1j*z]])
print("general form is anti-Hermitian traceless:",
      np.allclose(A.conj().T, -A) and np.isclose(np.trace(A), 0))
print("A = i(x s1 + y s2 + z s3):", np.allclose(A, 1j*(x*s1 + y*s2 + z*s3)))
UA = expm(A)
print("exp(A) in SU(2):",
      np.allclose(UA.conj().T @ UA, I2) and np.isclose(np.linalg.det(UA).real, 1))

# quaternion basis as another basis of su(2)
print("i = i s3, j = -i s2, k = -i s1:",
      np.allclose(qi, 1j*s3) and np.allclose(qj, -1j*s2) and np.allclose(qk, -1j*s1))

# grade 2: quaternion units factor into products of two Pauli matrices
print("s1 s2 = i s3, s2 s3 = i s1, s3 s1 = i s2:",
      np.allclose(s1@s2, 1j*s3) and np.allclose(s2@s3, 1j*s1) and np.allclose(s3@s1, 1j*s2))
print("i = s1 s2, j = s1 s3, k = s3 s2:",
      np.allclose(qi, s1@s2) and np.allclose(qj, s1@s3) and np.allclose(qk, s3@s2))
print("anticommute, (s1 s2)^2 = -I:",
      np.allclose(s1@s2, -s2@s1) and np.allclose((s1@s2) @ (s1@s2), -I2))

# conjugation R_q(x) = q x q^-1 with q = exp(i theta/2): rotation by theta about i-axis
q = expm(qi*theta/2)
def R(x): return q @ x @ q.conj().T  # q^-1 = q^dagger for unit q
print("R_q(i) = i:", np.allclose(R(qi), qi))
print("R_q(j) = j cos + k sin:", np.allclose(R(qj), np.cos(theta)*qj + np.sin(theta)*qk))
print("R_q(k) = -j sin + k cos:", np.allclose(R(qk), -np.sin(theta)*qj + np.cos(theta)*qk))
print("R_{-q} = R_q:", np.allclose((-q) @ qj @ (-q).conj().T, R(qj)))

# rotation matrix from R_q is in SO(3)
def coords(x):  # pure quaternion -> (x,y,z)
    return np.array([(x[0,0]/1j).real, -x[0,1].real, -(x[0,1]/1j).real])
Rot = np.column_stack([coords(R(m)) for m in (qi, qj, qk)])
print("R_q as 3x3: orthogonal, det 1:",
      np.allclose(Rot.T @ Rot, np.eye(3)) and np.isclose(np.linalg.det(Rot), 1))

# so(3): [Lx,Ly]=Lz cyclic; i/2 -> Lx etc. matches structure constants
Lx = np.array([[0,0,0],[0,0,-1],[0,1,0]], dtype=float)
Ly = np.array([[0,0,1],[0,0,0],[-1,0,0]], dtype=float)
Lz = np.array([[0,-1,0],[1,0,0],[0,0,0]], dtype=float)
print("[Lx,Ly]=Lz, [Ly,Lz]=Lx, [Lz,Lx]=Ly:",
      np.allclose(comm(Lx,Ly), Lz) and np.allclose(comm(Ly,Lz), Lx) and np.allclose(comm(Lz,Lx), Ly))
print("[i/2,j/2]=k/2 etc.:",
      np.allclose(comm(qi/2,qj/2), qk/2) and np.allclose(comm(qj/2,qk/2), qi/2) and np.allclose(comm(qk/2,qi/2), qj/2))
print("exp(i theta/2) gives rotation exp(theta Lx):", np.allclose(Rot, expm(theta*Lx).real))

# det exp(A) = exp(tr A)
A = Q(0, *rng.standard_normal(3))
print("det exp(A) = exp(tr A) = 1 for traceless A:", np.isclose(np.linalg.det(expm(A)).real, 1))
