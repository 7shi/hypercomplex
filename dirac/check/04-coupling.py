"""Checks for dirac/04-coupling.md (coupling to the electromagnetic field, gauge, current).

1. Gauge: with A' = A + D chi and psi' = psi e^{omega sigma_3 alpha}, alpha = -q chi/(hbar c),
   hbar D psi' omega sigma_3 - (q/c) A' psi' - m c psi' gamma_0
   = (hbar D psi omega sigma_3 - (q/c) A psi - m c psi gamma_0) e^{omega sigma_3 alpha}
   (symbolic functions). In the matrix picture this is Psi' = e^{-i q chi/(hbar c)} Psi.
2. Observables under the phase: psi' gamma_0 psi'~ = J, psi' gamma_3 psi'~ unchanged,
   psi' gamma_1 psi'~ = rho (cos 2alpha e_1 - sin 2alpha e_2), psi' gamma_2 psi'~ = rho (sin 2alpha e_1 + cos 2alpha e_2).
3. Current: <D(psi gamma_0 psi~)>_0 = 2 <(D psi) gamma_0 psi~>_0 (random polynomial psi);
   <psi omega sigma_3 psi~>_0 = 0 and psi gamma_0 omega sigma_3 psi~ is a trivector, so
   <A psi gamma_0 omega sigma_3 psi~>_0 = 0; hence D.J = 0 for solutions.
   J = psi gamma_0 psi~ is future-pointing timelike (J^2 = rho^2 > 0, J.gamma_0 = Psi^dagger Psi > 0)
   whenever psi psi~ != 0 (numerical samples).
4. Charge conjugation: psi -> psi sigma_1 maps solutions for charge q to solutions for -q, and
   the current of psi sigma_1 is (psi sigma_1) gamma_0 (psi sigma_1)~ = J (unchanged).
"""

import random

import sympy as sp

from common.clifford import Alg, MV, mv

S = Alg(1, 3)
g = S.e
g0 = g[0]
sg = [g[k] * g0 for k in (1, 2, 3)]
w = S.I
ws3 = w * sg[2]
X = S.X
ebasis = [mv(1, S.neg)] + sg + [w * x for x in sg] + [w]


def eq(A, B):
    return all(sp.simplify(sp.expand(v)) == 0 for v in (mv(A) - mv(B)).d.values())


hbar, c, m_ = sp.symbols("hbar c m", positive=True)
q = sp.Symbol("q", real=True)

print("1. gauge transformation")
fs = [sp.Function(f"f{i}")(*X) for i in range(8)]
psi = sum((fs[i] * ebasis[i] for i in range(8)), S.zero())
ph = sp.Function("varphi")(*X)
Ak = [sp.Function(f"A{k}")(*X) for k in (1, 2, 3)]
A = ph * g0 + c * sum((Ak[k] * g[k + 1] for k in range(3)), S.zero())
chi = sp.Function("chi")(*X)
al = -q * chi / (hbar * c)
U = sp.cos(al) + ws3 * sp.sin(al)


def dirac(psi, A):
    return hbar * S.D(psi) * ws3 - (q / c) * A * psi - m_ * c * psi * g0


Dchi = S.D(mv(chi, S.neg))
lhs = dirac(psi * U, A + Dchi)
rhs = dirac(psi, A) * U
assert eq(lhs, rhs)
assert eq(U * g0, g0 * U)
print("   ok")

print("2. frame under the phase")
cs = sp.symbols("c0:8", real=True)
psic = sum((cs[i] * ebasis[i] for i in range(8)), S.zero())
a = sp.Symbol("alpha", real=True)
Ua = sp.cos(a) + ws3 * sp.sin(a)
p2 = psic * Ua
E = [psic * g[m] * psic.rev() for m in range(4)]
assert eq(p2 * g0 * p2.rev(), E[0])
assert eq(p2 * g[3] * p2.rev(), E[3])
assert eq(p2 * g[1] * p2.rev(), sp.cos(2 * a) * E[1] - sp.sin(2 * a) * E[2])
assert eq(p2 * g[2] * p2.rev(), sp.sin(2 * a) * E[1] + sp.cos(2 * a) * E[2])
print("   ok")

print("3. conservation of the current")
P = S.rnd(5, deg=2, grades=(0, 2, 4))
lhs = S.D(P * g0 * P.rev()).scalar()
rhs = 2 * (S.D(P) * g0 * P.rev()).scalar()
assert sp.expand(lhs - rhs) == 0
assert sp.expand((P * ws3 * P.rev()).scalar()) == 0
T = P * g0 * ws3 * P.rev()
assert eq(T, T.grade(3))
Ar = S.rnd(6, deg=1, grades=(1,))
assert sp.expand((Ar * T).scalar()) == 0
r = random.Random(1)
for _ in range(20):
    ps = sum((r.randint(-3, 3) * b for b in ebasis), S.zero())
    pp = ps * ps.rev()
    if pp.scalar() == 0 and (pp * w).scalar() == 0:
        continue
    Jv = ps * g0 * ps.rev()
    J2 = (Jv * Jv).scalar()
    rho2 = pp.scalar()**2 + (pp * w).scalar()**2
    assert J2 == rho2 and (Jv * g0).scalar() > 0
print("   D.J = 0 identities, J timelike and future-pointing: ok")

print("4. charge conjugation")
lhs = hbar * S.D(psi * sg[0]) * ws3 + (q / c) * A * psi * sg[0] - m_ * c * psi * sg[0] * g0
assert eq(lhs, -dirac(psi, A) * sg[0])
pc = psic * sg[0]
assert eq(pc * g0 * pc.rev(), E[0])
print("   psi sigma_1 solves the equation with -q; its current is J: ok")
print("all ok")
