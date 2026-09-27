"""Checks for dirac/01-pauli.md (Pauli spinors as elements of Cl_{3,0}^0 = H).

1. Matrix representation of Cl_{3,0}: sigma_k -> Pauli matrices, omega = s1 s2 s3 -> i I;
   reversion -> Hermitian conjugate.
2. The map psi = a0 + a_k omega sigma_k -> Psi = (a0 + i a3, -a2 + i a1)^T is the first
   column of the matrix of psi; it is an R-linear bijection H -> C^2 and
   sigma_k psi sigma_3 <-> sigma_k Psi,  psi omega sigma_3 <-> i Psi.
3. Left ideal: P = (1 + sigma_3)/2, the matrix of psi P is [Psi, 0]; sigma_3 P = P,
   so omega sigma_3 P = omega P (omega central).
4. Observables: psi psi~ = Psi^dagger Psi, psi sigma_3 psi~ = sum_k s_k sigma_k with
   s_k = Psi^dagger sigma_k Psi = <psi sigma_3 psi~ sigma_k>_0; psi e^{omega sigma_3 alpha}
   corresponds to e^{i alpha} Psi and gives the same psi sigma_3 psi~.
5. Quaternions: i = -omega sigma_1, j = -omega sigma_2, k = -omega sigma_3 satisfy
   i^2 = j^2 = k^2 = ijk = -1, sigma_3 = omega k, so psi sigma_3 psi~ = omega (psi k psi~).
   Rotor R = e^{-omega sigma_3 theta/2} rotates sigma_1 to cos sigma_1 + sin sigma_2.
6. Pauli equation: momentum p_k Psi = -i hbar d_k Psi <-> -hbar d_k psi omega sigma_3;
   (sigma.B) Psi <-> B psi sigma_3. In a uniform field B = B sigma_3 (no spatial dependence),
   hbar d_t psi omega sigma_3 = -(g q hbar / 4m) B psi sigma_3 is solved by
   psi = e^{omega B g q t / 4m} psi_0, and the spin direction rotates with
   Omega = -g q B / 2m (theta = Omega t in e^{-omega sigma_3 theta/2}).
   Cross-check against the matrix equation i hbar dPsi/dt = -(g q hbar/4m) sigma.B Psi.
7. Supplements: the second column of the matrix of psi is (-Psi_2^*, Psi_1^*); for B along sigma_3
   Psi(t) = diag(e^{-i Omega t/2}, e^{i Omega t/2}) Psi(0); the phase rotates e_1 to
   cos(2 alpha) e_1 - sin(2 alpha) e_2; for g = 2 the spin rotor is e^{q B omega sigma_3 t/2m}.
"""

import sympy as sp

from common.clifford import Alg, MV, mv

A3 = Alg(3, 0)
s = A3.e
w = A3.I
I2 = sp.eye(2)
pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]


def eq(A, B):
    return all(sp.simplify(sp.expand(v)) == 0 for v in (mv(A) - mv(B)).d.values())


def mat(X):
    M = sp.zeros(2)
    for m, v in X.d.items():
        P = sp.eye(2)
        for a in range(3):
            if m >> a & 1:
                P = P * pauli[a]
        M += v * P
    return M


def meq(M, N):
    return sp.simplify(sp.expand(M - N)) == sp.zeros(*M.shape)


a = sp.symbols("a0:4", real=True)
psi = a[0] + sum((a[k + 1] * w * s[k] for k in range(3)), A3.zero())


def Psi_of(X):
    """Spinor column from an even element via the explicit formula."""
    b0 = X.scalar()
    bk = [-(X * w * s[k]).scalar() for k in range(3)]
    return sp.Matrix([b0 + sp.I * bk[2], -bk[1] + sp.I * bk[0]])


def even_of(Psi):
    """Inverse map C^2 -> H."""
    u, v = [sp.expand(x) for x in Psi]
    a0, a3 = sp.re(u), sp.im(u)
    a2, a1 = -sp.re(v), sp.im(v)
    return a0 + a1 * w * s[0] + a2 * w * s[1] + a3 * w * s[2]


print("1. representation")
for k in range(3):
    assert meq(mat(s[k]), pauli[k])
assert meq(mat(w), sp.I * I2)
Rnd = A3.rnd(1, deg=0)
assert meq(mat(Rnd.rev()), mat(Rnd).H)
print("   sigma_k -> Pauli, omega -> iI, reversion -> Hermitian conjugate: ok")

print("2. psi <-> Psi")
Psi = Psi_of(psi)
assert meq(Psi, sp.Matrix([a[0] + sp.I * a[3], -a[2] + sp.I * a[1]]))
assert meq(mat(psi)[:, 0], Psi)
assert eq(even_of(Psi), psi)
for k in range(3):
    assert meq(Psi_of(s[k] * psi * s[2]), pauli[k] * Psi)
assert meq(Psi_of(psi * w * s[2]), sp.I * Psi)
assert eq(w * s[2], s[0] * s[1])
print("   first column, bijection, sigma_k psi sigma_3, psi omega sigma_3 = i Psi: ok")

print("3. left ideal")
P = (1 + s[2]) / 2
assert eq(P * P, P)
M = mat(psi * P)
assert meq(M[:, 0], Psi) and meq(M[:, 1], sp.zeros(2, 1))
assert eq(s[2] * P, P) and eq(w * s[2] * P, w * P)
assert eq(w * psi, psi * w)
print("   psi P = [Psi, 0], omega sigma_3 P = omega P: ok")

print("4. observables")
rho = (psi * psi.rev())
assert eq(rho, sum(x**2 for x in a))
assert sp.simplify((Psi.H * Psi)[0] - rho.scalar()) == 0
S = psi * s[2] * psi.rev()
assert eq(S, S.grade(1))
for k in range(3):
    sk = sp.expand((Psi.H * pauli[k] * Psi)[0])
    assert sp.simplify(sk - (S * s[k]).scalar()) == 0
    assert sp.simplify(sp.im(sk)) == 0
al = sp.Symbol("alpha", real=True)
ph = sp.cos(al) + w * s[2] * sp.sin(al)
psi2 = psi * ph
assert eq(psi2 * s[2] * psi2.rev(), S)
assert meq(Psi_of(psi2), (sp.cos(al) + sp.I * sp.sin(al)) * Psi)
print("   rho = Psi^+Psi, s_k = Psi^+ sigma_k Psi, phase invariance: ok")

print("5. quaternions and rotation")
qi, qj, qk = -w * s[0], -w * s[1], -w * s[2]
for q in (qi, qj, qk):
    assert eq(q * q, -1)
assert eq(qi * qj * qk, -1) and eq(qi * qj, qk)
assert eq(w * qk, s[2])
assert eq(S, w * (psi * qk * psi.rev()))
assert eq(psi, a[0] - a[1] * qi - a[2] * qj - a[3] * qk)
th = sp.Symbol("theta", real=True)
R = sp.cos(th / 2) - w * s[2] * sp.sin(th / 2)
assert eq(R * s[0] * R.rev(), sp.cos(th) * s[0] + sp.sin(th) * s[1])
assert eq(R * s[1] * R.rev(), -sp.sin(th) * s[0] + sp.cos(th) * s[1])
print("   quaternion relations, sigma_3 = omega k, rotor direction: ok")

print("6. Pauli equation, precession")
hbar, m, B, t, g = sp.symbols("hbar m B t g", positive=True)
q = sp.Symbol("q", real=True)
Bv = [sp.Symbol(f"B{k}", real=True) for k in range(3)]
Bvec = sum((Bv[k] * s[k] for k in range(3)), A3.zero())
sigB = sum((Bv[k] * pauli[k] for k in range(3)), sp.zeros(2))
assert meq(Psi_of(Bvec * psi * s[2]), sigB * Psi)
# momentum: -i hbar d_k on components <-> -hbar d_k psi omega sigma_3
f = [sp.Function(f"f{j}")(*sp.symbols("x1:4", real=True)) for j in range(4)]
xs = sp.symbols("x1:4", real=True)
psif = f[0] + sum((f[j + 1] * w * s[j] for j in range(3)), A3.zero())
for k in range(3):
    lhs = Psi_of((-hbar) * psif.map(lambda v: sp.diff(v, xs[k])) * w * s[2])
    rhs = -sp.I * hbar * Psi_of(psif).diff(xs[k])
    assert meq(lhs, rhs)
# uniform field along sigma_3
Om = -g * q * B / (2 * m)
phi = g * q * B * t / (4 * m)
U = sp.cos(phi) + w * s[2] * sp.sin(phi)  # e^{omega sigma_3 g q B t/4m}
psit = U * psi
lhs = hbar * psit.map(lambda v: sp.diff(v, t)) * w * s[2]
rhs = -(g * q * hbar / (4 * m)) * B * s[2] * psit * s[2]
assert eq(lhs, rhs)
# equals e^{-omega sigma_3 theta/2} with theta = Omega t
assert eq(U, (sp.cos(Om * t / 2) - w * s[2] * sp.sin(Om * t / 2)))
# matrix check: i hbar dPsi/dt = -(g q hbar/4m) B sigma_3 Psi
Pt = Psi_of(psit)
assert meq(sp.I * hbar * Pt.diff(t), -(g * q * hbar / (4 * m)) * B * pauli[2] * Pt)
St = psit * s[2] * psit.rev()
S0 = S
x0 = (S0 * s[0]).scalar()
y0 = (S0 * s[1]).scalar()
assert sp.simplify((St * s[0]).scalar() - (x0 * sp.cos(Om * t) - y0 * sp.sin(Om * t))) == 0
assert sp.simplify((St * s[1]).scalar() - (x0 * sp.sin(Om * t) + y0 * sp.cos(Om * t))) == 0
print("   sigma.B, momentum, precession Omega = -gqB/2m: ok")

print("7. supplements")
M = mat(psi)
assert meq(M[:, 1], sp.Matrix([-sp.conjugate(Psi[1]), sp.conjugate(Psi[0])]))
# example: B along sigma_3, Psi(t) = diag(e^{-i Om t/2}, e^{i Om t/2}) Psi(0)
Dg = sp.diag(sp.cos(Om * t / 2) - sp.I * sp.sin(Om * t / 2), sp.cos(Om * t / 2) + sp.I * sp.sin(Om * t / 2))
assert meq(Psi_of(psit), Dg * Psi)
# global phase rotates the frame e_1, e_2 by -2 alpha about s
E1, E2 = psi * s[0] * psi.rev(), psi * s[1] * psi.rev()
assert eq(psi2 * s[0] * psi2.rev(), sp.cos(2 * al) * E1 - sp.sin(2 * al) * E2)
# g = 2 spin rotor equals the magnetic rotor of em/05: e^{q B omega sigma_3 t / 2m}
assert eq(U.map(lambda v: v.subs(g, 2)), sp.cos(q * B * t / (2 * m)) + w * s[2] * sp.sin(q * B * t / (2 * m)))
print("   second column, diagonal example, frame rotation, g = 2 vs em/05: ok")
print("all ok")
