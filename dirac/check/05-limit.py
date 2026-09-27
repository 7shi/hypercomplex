"""Checks for dirac/05-limit.md (non-relativistic limit and g = 2).

1. Split: psi = (psi_+ + psi_-) e^{-omega sigma_3 m c^2 t/hbar} with gamma_0 psi_+- gamma_0 = +-psi_+-.
   With Pi X = -hbar D_3 X omega sigma_3 - q A_vec X (D_3 = sum sigma_k d_k, A_vec = sum A_k sigma_k),
   gamma_0 (hbar D psi omega sigma_3 - (q/c) A psi - m c psi gamma_0) e^{omega sigma_3 m c^2 t/hbar}
   = (1/c) [ (hbar d_t psi_+ omega sigma_3 - q phi psi_+ - c Pi psi_-)
           + (hbar d_t psi_- omega sigma_3 - q phi psi_- - c Pi psi_+ + 2 m c^2 psi_-) ]
   and the two brackets are the parts commuting / anticommuting with gamma_0 (symbolic functions).
2. Pi maps the commuting part to the anticommuting part and vice versa.
3. Pi^2 X = -hbar^2 Lap X + q hbar ((div A) X + 2 (A.grad) X) omega sigma_3 + q^2 |A|^2 X
            - q hbar B X sigma_3,  B = curl A,
   and the last term comes from D_3 ^ A_vec = omega B (the wedge part of D_3 A_vec).
4. Matrix check in Cl_{3,0}: the map of 01 sends (Pi X) sigma_3 to (sigma . pi) Psi with pi = p - qA, and
   Pi^2 X to (pi^2 - q hbar sigma.B) Psi.
5. Pauli equation obtained: hbar d_t psi_+ omega sigma_3 = (1/2m) Pi^2 psi_+ + q phi psi_+, i.e. the
   01 form with g = 2. Plane-wave size estimate: for a free plane wave with momentum p along x_1,
   the exact ratio |psi_-| / |psi_+| = |p| / (m c + E/c) -> |p| / 2 m c.
"""

import sympy as sp

from common.clifford import Alg, MV, mv

S = Alg(1, 3)
g = S.e
g0 = g[0]
sg = [g[k] * g0 for k in (1, 2, 3)]
w = S.I
ws3 = w * sg[2]
X = S.X


def eq(A, B):
    return all(sp.simplify(sp.expand(v)) == 0 for v in (mv(A) - mv(B)).d.values())


hbar, c, m_ = sp.symbols("hbar c m", positive=True)
q = sp.Symbol("q", real=True)
t = sp.Symbol("t", real=True)
pe = [mv(1, S.neg)] + [w * x for x in sg]  # commute with gamma_0
fp = [sp.Function(f"u{i}")(*X) for i in range(4)]
fm = [sp.Function(f"v{i}")(*X) for i in range(4)]
psp = sum((fp[i] * pe[i] for i in range(4)), S.zero())
psm = sum((fm[i] * pe[i] for i in range(4)), S.zero()) * sg[2]
assert eq(g0 * psp * g0, psp) and eq(g0 * psm * g0, -psm)
ph = sp.Function("varphi")(*X)
Ak = [sp.Function(f"A{k}")(*X) for k in (1, 2, 3)]
A = ph * g0 + c * sum((Ak[k] * g[k + 1] for k in range(3)), S.zero())
Avec = sum((Ak[k] * sg[k] for k in range(3)), S.zero())


def d(F, a):
    return F.map(lambda v: sp.diff(v, X[a]))


def D3(F):
    return sum((sg[k] * d(F, k + 1) for k in range(3)), S.zero())


def dt(F):  # x0 = c t
    return c * d(F, 0)


def Pi(F):
    return -hbar * D3(F) * ws3 - q * Avec * F


print("1. split of the Dirac equation")
th = m_ * c * X[0] / hbar  # m c^2 t / hbar with x0 = c t
Em = sp.cos(th) - ws3 * sp.sin(th)
Ep = sp.cos(th) + ws3 * sp.sin(th)
psi = (psp + psm) * Em
full = g0 * (hbar * S.D(psi) * ws3 - (q / c) * A * psi - m_ * c * psi * g0) * Ep
Lp = hbar * dt(psp) * ws3 - q * ph * psp - c * Pi(psm)
Lm = hbar * dt(psm) * ws3 - q * ph * psm - c * Pi(psp) + 2 * m_ * c**2 * psm
assert eq(full, (Lp + Lm) / c)
assert eq(g0 * Lp * g0, Lp) and eq(g0 * Lm * g0, -Lm)
print("   ok")

print("2. Pi swaps the parts")
assert eq(g0 * Pi(psp) * g0, -Pi(psp))
assert eq(g0 * Pi(psm) * g0, -(-Pi(psm)))
print("   ok")

print("3. Pi^2 and the magnetic term")
Bk = [sp.diff(Ak[2], X[2]) - sp.diff(Ak[1], X[3]),
      sp.diff(Ak[0], X[3]) - sp.diff(Ak[2], X[1]),
      sp.diff(Ak[1], X[1]) - sp.diff(Ak[0], X[2])]
Bvec = sum((Bk[k] * sg[k] for k in range(3)), S.zero())
D3A = D3(Avec)
divA = sum(sp.diff(Ak[k], X[k + 1]) for k in range(3))
assert eq(D3A, divA + w * Bvec)
lap = sum((d(d(psp, k), k) for k in (1, 2, 3)), S.zero())
Agrad = sum((Ak[k] * d(psp, k + 1) for k in range(3)), S.zero())
A2 = sum(a**2 for a in Ak)
rhs = -hbar**2 * lap + q * hbar * (divA * psp + 2 * Agrad) * ws3 + q**2 * A2 * psp - q * hbar * Bvec * psp * sg[2]
assert eq(Pi(Pi(psp)), rhs)
print("   ok")

print("4. matrix check in Cl_{3,0}")
C = Alg(3, 0)
s, wc = C.e, C.I
pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
xs = C.X  # x0, x1, x2 used as spatial coordinates 1..3


def Phi(Y):
    b0 = Y.scalar()
    bk = [-(Y * wc * s[k]).scalar() for k in range(3)]
    return sp.Matrix([b0 + sp.I * bk[2], -bk[1] + sp.I * bk[0]])


fc = [sp.Function(f"h{i}")(*xs) for i in range(4)]
Y = fc[0] + sum((fc[i + 1] * wc * s[i] for i in range(3)), C.zero())
Ac = [sp.Function(f"a{k}")(*xs) for k in range(3)]
Acv = sum((Ac[k] * s[k] for k in range(3)), C.zero())


def PiC(F):
    return -hbar * sum((s[k] * F.map(lambda v: sp.diff(v, xs[k])) for k in range(3)), C.zero()) * wc * s[2] - q * Acv * F


Psi = Phi(Y)


def pik(V, k):
    return -sp.I * hbar * V.diff(xs[k]) - q * Ac[k] * V


sigpi = sum((pauli[k] * pik(Psi, k) for k in range(3)), sp.zeros(2, 1))
assert sp.simplify(Phi(PiC(Y) * s[2]) - sigpi) == sp.zeros(2, 1)
Bc = [sp.diff(Ac[2], xs[1]) - sp.diff(Ac[1], xs[2]),
      sp.diff(Ac[0], xs[2]) - sp.diff(Ac[2], xs[0]),
      sp.diff(Ac[1], xs[0]) - sp.diff(Ac[0], xs[1])]
pi2 = sum((pik(pik(Psi, k), k) for k in range(3)), sp.zeros(2, 1))
sigB = sum((Bc[k] * pauli[k] for k in range(3)), sp.zeros(2))
assert sp.simplify(sp.expand(Phi(PiC(PiC(Y))) - (pi2 - q * hbar * sigB * Psi))) == sp.zeros(2, 1)
print("   (Pi X) sigma_3 <-> sigma.pi, Pi^2 <-> pi^2 - q hbar sigma.B: ok")

print("5. size of the small part for a free plane wave")
p1 = sp.Symbol("p", positive=True)
Ek = sp.sqrt(m_**2 * c**4 + p1**2 * c**2)
# rest spinor 1 boosted along x1: L = cosh(eta/2) + sigma_1 sinh(eta/2), tanh(eta) = p c / E
eta = sp.asinh(p1 / (m_ * c))
L = sp.cosh(eta / 2) + sg[0] * sp.sinh(eta / 2)
plus = (L + g0 * L * g0) / 2
minus = (L - g0 * L * g0) / 2
ratio = sp.sqrt((minus * minus.rev()).scalar() * -1) / sp.sqrt((plus * plus.rev()).scalar())
target = p1 / (m_ * c + Ek / c)
for vals in ({p1: 1, m_: 1, c: 1}, {p1: sp.Rational(3, 7), m_: 2, c: 3}, {p1: 5, m_: sp.Rational(1, 3), c: 2}):
    assert abs(sp.N((ratio - target).subs(vals), 30)) < 1e-25
assert sp.limit(target / (p1 / (2 * m_ * c)), p1, 0) == 1
print("   |psi_-|/|psi_+| = p/(mc + E/c): ok")
print("all ok")
