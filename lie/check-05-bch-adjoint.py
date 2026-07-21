"""Numerical checks for 05-bch-adjoint.md.

BCH expansion: the 2nd-order truncation z = x + y + [x,y]/2 leaves an
O(t^3) residual and the 3rd-order truncation (with the 1/12 terms) an
O(t^4) residual; exp(x)exp(y) = exp(x+y) exactly for commuting x, y;
the group commutator exp(x)exp(y)exp(-x)exp(-y) = exp([x,y] + O(t^3)),
with the so(3) instance e^{eJx}e^{eJy}e^{-eJx}e^{-eJy} ~ e^{e^2 Jz};
the exact quaternion composition exp(ua)exp(vb) = exp(wc) with
cos c = cos a cos b - (u.v) sin a sin b, reproducing exp(i pi/2)
exp(j pi/2) = k, and its 2nd-order term ab(u x v) = [ua, vb]/2; the
Jacobi identity and its Leibniz form for random matrices; the adjoint
action Ad_g x = gxg^{-1} staying in su(3) and preserving brackets;
Ad_{gh} = Ad_g Ad_h; g exp(y) g^{-1} = exp(Ad_g y) exactly;
Ad_{exp z} = exp(ad_z) (Hadamard) both as a series on su(3) and as
8x8 matrices; ad_{i/2} = Jx, ad_{j/2} = Jy, ad_{k/2} = Jz on pure
quaternions, and Ad_{exp(theta i/2)} = exp(theta Jx); ker Ad = {+-1};
the matrix entries of ad_{e_a} being the structure constants
(su(2): Levi-Civita; su(3): -2 f_abc for e_a = i lambda_a); and
ad_{[x,y]} = [ad_x, ad_y] (ad is a representation).
"""

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

def expm(A, terms=60):
    S, T = np.eye(*A.shape, dtype=complex), np.eye(*A.shape, dtype=complex)
    for n in range(1, terms):
        T = T @ A / n
        S = S + T
    return S

def comm(A, B):
    return A @ B - B @ A

rng = np.random.default_rng(0)
def su3(theta):
    return sum(t * 1j*m for t, m in zip(theta, lam))

x, y = su3(rng.standard_normal(8)), su3(rng.standard_normal(8))

# BCH: 2nd-order truncation has O(t^3) residual, 3rd-order has O(t^4)
def bch2(t):
    return t*(x + y) + t**2/2*comm(x, y)
def bch3(t):
    return bch2(t) + t**3/12*(comm(x, comm(x, y)) + comm(y, comm(y, x)))
def res(z, t):
    return np.linalg.norm(expm(t*x) @ expm(t*y) - expm(z(t)))
r2 = res(bch2, 0.1) / res(bch2, 0.05)
r3 = res(bch3, 0.1) / res(bch3, 0.05)
print("2nd-order BCH residual is O(t^3) (halving t gives ~1/8):", 7 < r2 < 9)
print("3rd-order BCH residual is O(t^4) (halving t gives ~1/16):", 14 < r3 < 18)

# commuting case: exact exponential law (diagonal su(3) elements l3, l8)
d1, d2 = 0.8j*l3, -0.5j*l8
print("[x,y] = 0 restores exp(x)exp(y) = exp(x+y):",
      np.allclose(comm(d1, d2), 0) and np.allclose(expm(d1) @ expm(d2), expm(d1 + d2)))

# group commutator: exp(x)exp(y)exp(-x)exp(-y) = exp([x,y] + O(t^3))
def resc(t):
    P = expm(t*x) @ expm(t*y) @ expm(-t*x) @ expm(-t*y)
    return np.linalg.norm(P - expm(t**2*comm(x, y)))
rc = resc(0.1) / resc(0.05)
print("group commutator residual is O(t^3):", 7 < rc < 9)

# so(3) instance: J-rotations around x and y leave an eps^2 rotation around z
Jx = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=float)
Jy = np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]], dtype=float)
Jz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
print("[Jx, Jy] = Jz:", np.allclose(comm(Jx, Jy), Jz))
def resj(e):
    P = expm(e*Jx + 0j) @ expm(e*Jy + 0j) @ expm(-e*Jx + 0j) @ expm(-e*Jy + 0j)
    return np.linalg.norm(P - expm(e**2*Jz + 0j))
print("e^{eJx}e^{eJy}e^{-eJx}e^{-eJy} = e^{e^2 Jz} + O(e^3):",
      7 < resj(0.1) / resj(0.05) < 9)

# quaternions as 4-arrays (w, x, y, z)
def qmul(p, q):
    (a, b, c, d), (e, f, g, h) = p, q
    return np.array([a*e - b*f - c*g - d*h, a*f + b*e + c*h - d*g,
                     a*g - b*h + c*e + d*f, a*h + b*g - c*f + d*e])
def qexp(u, a):
    return np.concatenate([[np.cos(a)], np.sin(a)*u])
qi, qj, qk = np.eye(3)

# exact composition: exp(ua)exp(vb) = exp(wc) via cos c = cos a cos b - (u.v) sin a sin b
u = rng.standard_normal(3); u /= np.linalg.norm(u)
v = rng.standard_normal(3); v /= np.linalg.norm(v)
a, b = 0.9, 1.7
prod = qmul(qexp(u, a), qexp(v, b))
cosc = np.cos(a)*np.cos(b) - (u @ v)*np.sin(a)*np.sin(b)
imag = np.sin(a)*np.cos(b)*u + np.cos(a)*np.sin(b)*v + np.sin(a)*np.sin(b)*np.cross(u, v)
print("exp(ua)exp(vb): real part = cos c, imaginary part = w sin c:",
      np.isclose(prod[0], cosc) and np.allclose(prod[1:], imag)
      and np.isclose(np.linalg.norm(imag), np.sqrt(1 - cosc**2)))
c = np.arccos(cosc)
print("product = exp(wc) with w = imag/|imag|:",
      np.allclose(prod, qexp(imag/np.linalg.norm(imag), c)))

# the example: exp(i pi/2) exp(j pi/2) = k = exp(k pi/2)
print("exp(i pi/2)exp(j pi/2) = k:",
      np.allclose(qmul(qexp(qi, np.pi/2), qexp(qj, np.pi/2)), [0, 0, 0, 1]))

# 2nd-order BCH term for pure quaternions: [ua, vb]/2 = ab (u x v)
pu, pv = np.concatenate([[0], a*u]), np.concatenate([[0], b*v])
half = (qmul(pu, pv) - qmul(pv, pu)) / 2
print("[ua, vb]/2 = ab (u x v):",
      np.isclose(half[0], 0) and np.allclose(half[1:], a*b*np.cross(u, v)))

# Jacobi identity and Leibniz form for random complex matrices
X, Y, Z = (rng.standard_normal((3, 3)) + 1j*rng.standard_normal((3, 3)) for _ in range(3))
print("Jacobi: [X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0:",
      np.allclose(comm(X, comm(Y, Z)) + comm(Y, comm(Z, X)) + comm(Z, comm(X, Y)), 0))
print("Leibniz: [X,[Y,Z]] = [[X,Y],Z] + [Y,[X,Z]]:",
      np.allclose(comm(X, comm(Y, Z)), comm(comm(X, Y), Z) + comm(Y, comm(X, Z))))

# adjoint action on su(3): Ad_g x = g x g^{-1} stays in su(3), preserves brackets
g, h = expm(x), expm(y)
def Ad(g, m):
    return g @ m @ g.conj().T
w3 = su3(rng.standard_normal(8))
print("Ad_g x anti-Hermitian traceless (stays in su(3)):",
      np.allclose(Ad(g, x).conj().T, -Ad(g, x)) and np.isclose(np.trace(Ad(g, x)), 0))
print("Ad_g [x,y] = [Ad_g x, Ad_g y]:",
      np.allclose(Ad(g, comm(x, w3)), comm(Ad(g, x), Ad(g, w3))))
print("Ad_{gh} = Ad_g Ad_h:", np.allclose(Ad(g @ h, w3), Ad(g, Ad(h, w3))))
print("g exp(y) g^{-1} = exp(Ad_g y) exactly:",
      np.allclose(g @ expm(y) @ g.conj().T, expm(Ad(g, y))))

# Hadamard: exp(z) x exp(-z) = sum ad_z^n(x)/n!
z8 = su3(rng.standard_normal(8))
S, T = np.array(w3), np.array(w3)
for n in range(1, 60):
    T = comm(z8, T) / n
    S = S + T
print("Ad_{exp z} x = exp(ad_z) x (series):", np.allclose(Ad(expm(z8), w3), S))

# matrix version: 8x8 Ad(exp z) equals expm of the 8x8 ad_z
def coords(m):
    return np.array([np.trace(l @ m).imag / 2 for l in lam])
admat = np.stack([coords(comm(z8, 1j*l)) for l in lam], axis=1)
U = expm(z8)
Admat = np.stack([coords(Ad(U, 1j*l)) for l in lam], axis=1)
print("Ad(exp z) = expm(ad_z) as 8x8 matrices:", np.allclose(Admat, expm(admat + 0j).real))

# ad on pure quaternions: ad_{i/2} = Jx, ad_{j/2} = Jy, ad_{k/2} = Jz
def qad(u):
    pu = np.concatenate([[0], u]) / 2
    cols = []
    for e in (qi, qj, qk):
        pe = np.concatenate([[0], e])
        cols.append((qmul(pu, pe) - qmul(pe, pu))[1:])
    return np.stack(cols, axis=1)
print("ad_{i/2} = Jx, ad_{j/2} = Jy, ad_{k/2} = Jz:",
      np.allclose(qad(qi), Jx) and np.allclose(qad(qj), Jy) and np.allclose(qad(qk), Jz))

# Ad_{exp(theta i/2)} on (i,j,k) = exp(theta Jx) = rotation about the x-axis
theta = 0.7
q = qexp(qi, theta/2)
qc = q * [1, -1, -1, -1]
Adq = np.stack([qmul(qmul(q, np.concatenate([[0], e])), qc)[1:] for e in (qi, qj, qk)],
               axis=1)
Rx = np.array([[1, 0, 0],
               [0, np.cos(theta), -np.sin(theta)],
               [0, np.sin(theta), np.cos(theta)]])
print("Ad_{exp(theta i/2)} = exp(theta Jx) (rotation by theta):",
      np.allclose(Adq, expm(theta*Jx + 0j).real) and np.allclose(Adq, Rx))

# kernel of Ad: Ad_{-1} = id; a generic q != +-1 acts nontrivially
qm1 = np.array([-1, 0, 0, 0])
Adm1 = np.stack([qmul(qmul(qm1, np.concatenate([[0], e])), qm1)[1:] for e in (qi, qj, qk)],
                axis=1)
print("Ad_{-1} = identity:", np.allclose(Adm1, np.eye(3)))
print("Ad_q != identity for generic q:", not np.allclose(Adq, np.eye(3)))

# structure constants as matrix entries of ad
# su(2) with basis i/2, j/2, k/2: entries of ad are the Levi-Civita symbols
eps = np.zeros((3, 3, 3))
eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1
eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1
print("(ad_{e_a})_{cb} = eps_{abc} for su(2):",
      all(np.allclose(qad(u)[c, b_], eps[a_, b_, c])
          for a_, u in enumerate((qi, qj, qk)) for b_ in range(3) for c in range(3)))

# su(3) with basis e_a = i lambda_a: [e_a, e_b] = sum_c (-2 f_abc) e_c
f = np.zeros((8, 8, 8))
for a_ in range(8):
    for b_ in range(8):
        for c in range(8):
            f[a_, b_, c] = (np.trace(comm(lam[a_], lam[b_]) @ lam[c]) / 4j).real
admats = [np.stack([coords(comm(1j*la, 1j*lb)) for lb in lam], axis=1) for la in lam]
print("(ad_{e_a})_{cb} = -2 f_abc for su(3):",
      all(np.allclose(admats[a_], -2*f[a_].T) for a_ in range(8)))

# ad is a representation: ad_{[x,y]} = [ad_x, ad_y]
adx = np.stack([coords(comm(x, 1j*l)) for l in lam], axis=1)
ady = np.stack([coords(comm(y, 1j*l)) for l in lam], axis=1)
adxy = np.stack([coords(comm(comm(x, y), 1j*l)) for l in lam], axis=1)
print("ad_{[x,y]} = [ad_x, ad_y]:", np.allclose(adxy, comm(adx, ady)))
