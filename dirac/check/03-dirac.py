"""Checks for dirac/03-dirac.md (the Dirac equation in Hestenes form).

1. Dirac representation: g^0^ = diag(I, -I), g^k^ = [[0, sk], [-sk, 0]] (upper index); they satisfy the Cl_{1,3}
   relations. For an even psi = phi + eta sigma_3 (phi = (psi + g0 psi g0)/2 commutes with gamma_0,
   eta in span{1, omega sigma_k}), Psi = (Phi(phi); Phi(eta)) with the map of 01.
   Then g_mu^ Psi <-> gamma_mu psi gamma_0, i Psi <-> psi omega sigma_3,
   g5^ = i g0^ g1^ g2^ g3^ Psi <-> psi sigma_3. The map is an R-linear bijection R^8 -> C^4.
2. Equation: gamma^mu^ (i hbar d_mu - q A_mu) Psi - m c Psi corresponds to
   (hbar D psi omega sigma_3 - (q/c) A psi - m c psi gamma_0) gamma_0, with
   A = phi gamma_0 + c A_k gamma_k and A_mu = (phi/c, -A_k) (symbolic functions).
3. Squaring: hbar D psi omega sigma_3 = m c psi gamma_0 implies (box + m^2 c^2/hbar^2) psi = 0
   (algebra: D psi = -(mc/hbar) psi gamma_0 omega sigma_3; gamma_0 commutes with omega sigma_3).
4. Plane waves: psi = psi0 e^{-omega sigma_3 p.x/hbar} gives hbar D psi omega sigma_3 = p psi;
   rest solution p = m c gamma_0 with psi0 commuting with gamma_0; boosted solution L psi0 with
   p = L (m c gamma_0) L~; right multiplication by sigma_1 maps solutions to solutions with q -> -q
   and gives the negative-energy plane wave psi0 sigma_1 e^{+omega sigma_3 p.x/hbar}.
5. Covariance: for a constant rotor R, D_x [R psi(R~ x R)] = R (D psi)(R~ x R) (random polynomial psi).
6. Bilinears: Psi^dagger Psi = <psi gamma_0 psi~ gamma_0>_0 = J.gamma_0 (J = psi gamma_0 psi~),
   Psi-bar Psi = <psi psi~>_0, Psi-bar g^mu^ Psi = J . gamma^mu,
   Psi-bar i g5^ Psi = -(omega coefficient of psi psi~), i.e. rho cos(beta) and -rho sin(beta).
7. Supplements: g5^ = [[0, I], [I, 0]]; gamma_0 P_+ = P_- gamma_0 and D(psi P_+) = (D psi) P_+;
   the boosted positive-energy solution has beta = 0 and J = rho L gamma_0 L~ (parallel to p);
   psi0 sigma_1 has psi psi~ = -rho (beta = pi) and the same J.
"""

import random

import sympy as sp

from common.clifford import Alg, MV, mv

S = Alg(1, 3)
g = S.e
g0 = g[0]
gr = S.er  # reciprocal frame
sg = [g[k] * g0 for k in (1, 2, 3)]
w = S.I
ws3 = w * sg[2]
pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
Z2, I2 = sp.zeros(2), sp.eye(2)
# standard Dirac representation: upper index g^k^ = [[0, sk], [-sk, 0]]
ghr = [sp.Matrix(sp.BlockMatrix([[I2, Z2], [Z2, -I2]]))] + \
      [sp.Matrix(sp.BlockMatrix([[Z2, p], [-p, Z2]])) for p in pauli]
gh = [ghr[0]] + [-x for x in ghr[1:]]  # lower index
g5h = sp.I * ghr[0] * ghr[1] * ghr[2] * ghr[3]


def eq(A, B):
    return all(sp.simplify(sp.expand(v)) == 0 for v in (mv(A) - mv(B)).d.values())


def meq(M, N):
    return sp.simplify(sp.expand(M - N)) == sp.zeros(*M.shape)


def Phi(X):
    b0 = X.scalar()
    bk = [-(X * w * sg[k]).scalar() for k in range(3)]
    return sp.Matrix([b0 + sp.I * bk[2], -bk[1] + sp.I * bk[0]])


def Psi_of(psi):
    phi = (psi + g0 * psi * g0) / 2
    eta = (psi - g0 * psi * g0) / 2 * sg[2]
    return sp.Matrix.vstack(Phi(phi), Phi(eta))


ebasis = [mv(1, S.neg)] + sg + [w * x for x in sg] + [w]
cs = sp.symbols("c0:8", real=True)
psi = sum((cs[i] * ebasis[i] for i in range(8)), S.zero())

print("1. Dirac representation and the map")
for m in range(4):
    for n in range(4):
        ac = gh[m] * gh[n] + gh[n] * gh[m]
        assert meq(ac, 2 * (1 if m == n == 0 else -1 if m == n else 0) * sp.eye(4))
Ps = Psi_of(psi)
J = sp.Matrix([[sp.diff(v, c) for c in cs] for v in Ps])
# R-linear bijection: real 8x8 matrix of (Re, Im) parts
Mr = sp.Matrix([[sp.re(J[i, j]) for j in range(8)] for i in range(4)] +
               [[sp.im(J[i, j]) for j in range(8)] for i in range(4)])
assert Mr.det() != 0
for m in range(4):
    assert meq(Psi_of(g[m] * psi * g0), gh[m] * Ps)
assert meq(Psi_of(psi * ws3), sp.I * Ps)
assert meq(Psi_of(psi * sg[2]), g5h * Ps)
assert eq(ws3 * g0, g0 * ws3)
print("   gamma_mu psi gamma_0, psi omega sigma_3 = i, psi sigma_3 = gamma_5: ok")

print("2. the equation with the field")
hbar, c, m_ = sp.symbols("hbar c m", positive=True)
q = sp.Symbol("q", real=True)
X = S.X
fs = [sp.Function(f"f{i}")(*X) for i in range(8)]
psif = sum((fs[i] * ebasis[i] for i in range(8)), S.zero())
ph = sp.Function("varphi")(*X)
Ak = [sp.Function(f"A{k}")(*X) for k in (1, 2, 3)]
A = ph * g0 + c * sum((Ak[k] * g[k + 1] for k in range(3)), S.zero())
A_low = [ph / c] + [-a for a in Ak]
ga = hbar * S.D(psif) * ws3 - (q / c) * A * psif - m_ * c * psif * g0
Pf = Psi_of(psif)
mat = sum((ghr[mu] * (sp.I * hbar * Pf.diff(X[mu]) - q * A_low[mu] * Pf) for mu in range(4)),
          sp.zeros(4, 1)) - m_ * c * Pf
assert meq(Psi_of(ga * g0), mat)
print("   ok")

print("3. squaring gives Klein-Gordon")
# if D psi = -(mc/hbar) psi g0 ws3 then D^2 psi = -(mc/hbar) (D psi) g0 ws3
Dpsi = -(m_ * c / hbar) * psi * g0 * ws3
D2 = -(m_ * c / hbar) * Dpsi * g0 * ws3
assert eq(D2, -(m_ * c / hbar)**2 * psi)
assert eq(-(ws3 * ws3), 1)
print("   ok")

print("4. plane waves")
p = sp.symbols("p0:4", real=True)
pv = sum((p[i] * g[i] for i in range(4)), S.zero())
px = p[0] * X[0] - p[1] * X[1] - p[2] * X[2] - p[3] * X[3]
assert sp.simplify((pv * sum((X[i] * g[i] for i in range(4)), S.zero())).scalar() - px) == 0
E = sp.cos(px / hbar) - ws3 * sp.sin(px / hbar)
pw = psi * E
assert eq(hbar * S.D(pw) * ws3, pv * pw)
# rest solution
psi0 = cs[0] + cs[4] * w * sg[0] + cs[5] * w * sg[1] + cs[6] * w * sg[2]
assert eq(m_ * c * g0 * psi0, m_ * c * psi0 * g0)
# boosted
eta = sp.Symbol("eta", real=True)
L = sp.cosh(eta / 2) + sg[0] * sp.sinh(eta / 2)
pb = L * (m_ * c * g0) * L.rev()
assert eq(pb, m_ * c * (sp.cosh(eta) * g0 + sp.sinh(eta) * g[1]))
assert eq(pb * (L * psi0), m_ * c * (L * psi0) * g0)
# sigma_1 on the right: solutions with q -> -q
ga1 = hbar * S.D(psif * sg[0]) * ws3 - (-q / c) * A * (psif * sg[0]) - m_ * c * psif * sg[0] * g0
assert eq(ga1, -ga * sg[0])
neg = psi0 * sg[0] * (sp.cos(px / hbar) + ws3 * sp.sin(px / hbar))
assert eq(psi0 * E * sg[0], neg)
# with p = m c gamma_0: neg has time factor e^{+omega sigma_3 m c x0/hbar}
negr = neg.map(lambda v: v.subs({p[0]: m_ * c, p[1]: 0, p[2]: 0, p[3]: 0}))
assert eq(hbar * S.D(negr) * ws3, m_ * c * negr * g0)
assert eq((-m_ * c * g0) * (psi0 * sg[0]), m_ * c * psi0 * sg[0] * g0)
print("   ok")

print("5. covariance")
th = sp.Rational(1, 3)
Rr = (sp.cosh(sp.Rational(1, 4)) + sg[1] * sp.sinh(sp.Rational(1, 4))) * (sp.cos(th) - w * sg[0] * sp.sin(th))
assert eq(Rr * Rr.rev(), 1)
x = sum((X[i] * g[i] for i in range(4)), S.zero())
xp = Rr.rev() * x * Rr
xpc = [(xp * gr[i]).scalar() for i in range(4)]  # components x'^i
sub = dict(zip(X, xpc))
F = S.rnd(3, deg=2, grades=(0, 2, 4))
lhs = S.D(Rr * F.map(lambda v: v.subs(sub, simultaneous=True)))
rhs = Rr * S.D(F).map(lambda v: v.subs(sub, simultaneous=True))
assert all(abs(sp.N(v)) < 1e-9 for v in (lhs - rhs).map(lambda v: sp.expand(v).subs(
    dict(zip(X, (sp.Rational(1, 2), sp.Rational(-1, 3), 2, sp.Rational(3, 7)))))).d.values())
print("   ok")

print("6. bilinears")
Jv = psi * g0 * psi.rev()
bar = Ps.H * gh[0]
assert sp.simplify((Ps.H * Ps)[0] - (Jv * g0).scalar()) == 0
assert sp.simplify((bar * Ps)[0] - (psi * psi.rev()).scalar()) == 0
for mu in range(4):
    assert sp.simplify((bar * ghr[mu] * Ps)[0] - (Jv * gr[mu]).scalar()) == 0
v5 = sp.expand((bar * sp.I * g5h * Ps)[0])
pp4 = sp.expand(-(psi * psi.rev() * w).scalar())  # coefficient of omega
assert sp.simplify(v5 + pp4) == 0
print("   Psi^+Psi = J.gamma_0, Psi-bar Psi = <psi psi~>_0, Psi-bar g^mu Psi = J.gamma^mu, Psi-bar i g5 Psi: ok")

print("7. supplements")
assert meq(g5h, sp.Matrix(sp.BlockMatrix([[Z2, I2], [I2, Z2]])))
Pp, Pm = (1 + sg[2]) / 2, (1 - sg[2]) / 2
assert eq(g0 * Pp, Pm * g0)
F = S.rnd(8, deg=1, grades=(0, 2, 4))
assert eq(S.D(F * Pp), S.D(F) * Pp)
# boosted positive-energy solution: beta = 0 and J parallel to p
psi0L = L * psi0
pp_ = psi0L * psi0L.rev()
assert eq(pp_, pp_.grade(0)) and sp.simplify(pp_.scalar() - sum(cs[i]**2 for i in (0, 4, 5, 6))) == 0
Jb = psi0L * g0 * psi0L.rev()
assert eq(Jb, pp_.scalar() * (L * g0 * L.rev()))
# negative-energy: beta = pi, same current
pn = psi0L * sg[0]
assert eq(pn * pn.rev(), -pp_) and eq(pn * g0 * pn.rev(), Jb)
print("   gamma_5 block form, Weyl projectors, beta and J of plane waves: ok")
print("all ok")
