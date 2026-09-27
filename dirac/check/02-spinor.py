"""Checks for dirac/02-spinor.md (spacetime spinors, Spin+(1,3) = SL(2,C)).

1. sigma_k = gamma_k gamma_0, omega = gamma_0 gamma_1 gamma_2 gamma_3 = sigma_1 sigma_2 sigma_3,
   omega sigma_3 = gamma_2 gamma_1. The even subalgebra maps to M_2(C) by sigma_k -> Pauli,
   omega -> iI (algebra homomorphism, bijective on the 8-dim even part).
2. On even elements: the spacetime reversion X~ corresponds to the adjugate (Clifford conjugation
   of Cl_{3,0}), so X X~ = <X>_0 + <X>_4 corresponds to det(X) I; gamma_0 X~ gamma_0 corresponds
   to the Hermitian conjugate.
3. For a vector x, x gamma_0 = x_0 + x_k sigma_k is Hermitian with det = x^2, and
   (R x R~) gamma_0 = R (x gamma_0) (gamma_0 R~ gamma_0) corresponds to R X R^dagger.
4. Kernel: an even element commuting with all gamma_mu is a + b omega with b = 0 (omega
   anticommutes with vectors); R R~ = 1 gives a = +-1. Rotation by 2 pi: e^{-omega sigma_3 pi} = -1.
5. Boost to a future unit timelike u: R = (1 + u gamma_0)/sqrt(2(1 + u.gamma_0)) is a rotor with
   R gamma_0 R~ = u; then R~ L(gamma_0) R = gamma_0 means the rest is a spatial rotation.
6. Even elements commuting with gamma_0 are spanned by 1, omega sigma_k (the Pauli-even part).
7. Decomposition: psi psi~ = rho e^{omega beta}; R = e^{-omega beta/2} psi / sqrt(rho) is a rotor;
   psi gamma_mu psi~ = rho e_mu with e_mu = R gamma_mu R~ (numerically for random psi).
   Under psi -> L psi, psi psi~ is invariant and psi gamma_mu psi~ -> L (.) L~.
8. Supplements: reversion of omega sigma_k and sigma_k in Cl_{1,3}; for psi commuting with gamma_0,
   psi gamma_0 psi~ = rho gamma_0 and psi sigma_3 psi~ = (psi gamma_3 psi~) gamma_0; the boost rotor
   is a positive definite Hermitian matrix; sigma_1 sigma_1~ = -1; e^{omega sigma_3 alpha} fixes
   gamma_0 and gamma_3; e_0 and e_3 are orthogonal.
"""

import math
import random

import sympy as sp

from common.clifford import Alg, MV, mv

S = Alg(1, 3)
g = S.e
g0 = g[0]
sg = [g[k] * g0 for k in (1, 2, 3)]
w = S.I
pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]


def eq(A, B):
    return all(sp.simplify(sp.expand(v)) == 0 for v in (mv(A) - mv(B)).d.values())


def meq(M, N):
    return sp.simplify(sp.expand(M - N)) == sp.zeros(*M.shape)


# basis of the even part in terms of sigma's: 1, sigma_k, omega sigma_k, omega
ebasis = [mv(1, S.neg)] + sg + [w * x for x in sg] + [w]
emats = [sp.eye(2)] + pauli + [sp.I * p for p in pauli] + [sp.I * sp.eye(2)]


def coords(X):
    """Coordinates of an even element in ebasis (solve by the scalar product)."""
    out = []
    for b in ebasis:
        n = (b * b.rev()).scalar()  # +-1
        out.append(sp.expand((X * b.rev()).scalar() / n))
    return out


def mat(X):
    return sum((c * M for c, M in zip(coords(X), emats)), sp.zeros(2))


def rnd_even(seed):
    r = random.Random(seed)
    return sum((r.randint(-3, 3) * b for b in ebasis), S.zero())


print("1. even subalgebra and M_2(C)")
assert eq(w, sg[0] * sg[1] * sg[2])
assert eq(w * sg[2], g[2] * g[1])
for seed in range(5):
    X, Y = rnd_even(seed), rnd_even(seed + 100)
    assert eq(X, sum((c * b for c, b in zip(coords(X), ebasis)), S.zero()))
    assert meq(mat(X * Y), mat(X) * mat(Y))
print("   homomorphism on random elements: ok")

print("2. reversion = adjugate, psi psi~ = det")
for seed in range(5):
    X = rnd_even(seed)
    M = mat(X)
    assert meq(mat(X.rev()), M.adjugate())
    XX = X * X.rev()
    assert eq(XX, XX.grade(0) + XX.grade(4))
    assert meq(mat(XX), M.det() * sp.eye(2))
    assert meq(mat(g0 * X.rev() * g0), M.H)
print("   ok")

print("3. action on vectors = R X R^dagger")
xs = sp.symbols("x0:4", real=True)
x = sum((xs[m] * g[m] for m in range(4)), S.zero())
Xh = mat(x * g0)
assert meq(Xh, Xh.H)
assert sp.simplify(Xh.det() - (xs[0]**2 - xs[1]**2 - xs[2]**2 - xs[3]**2)) == 0
eta, th = sp.symbols("eta theta", real=True)
R = (sp.cosh(eta / 2) + sg[0] * sp.sinh(eta / 2)) * (sp.cos(th / 2) - w * sg[2] * sp.sin(th / 2))
assert eq(R * R.rev(), 1)
assert sp.simplify(mat(R).det() - 1) == 0
assert meq(mat((R * x * R.rev()) * g0), mat(R) * Xh * mat(R).H)
assert eq(R * x * R.rev() * g0, R * (x * g0) * (g0 * R.rev() * g0))
print("   ok")

print("4. kernel and 2 pi")
cs = sp.symbols("c0:8", real=True)
Z = sum((cs[i] * ebasis[i] for i in range(8)), S.zero())
eqs = []
for m in range(4):
    eqs += list((Z * g[m] - g[m] * Z).d.values())
sol = sp.solve(eqs, cs, dict=True)[0]
Zs = Z.map(lambda v: v.subs(sol))
assert eq(Zs, cs[0])
assert eq(sp.cos(sp.pi) - w * sg[2] * sp.sin(sp.pi), -1)
print("   commutant of vectors = scalars, e^{-omega sigma_3 pi} = -1: ok")

print("5. boost from gamma_0 to u")
u1, u2, u3 = sp.Rational(1, 2), sp.Rational(-2, 3), sp.Rational(1, 5)
u0 = sp.sqrt(1 + u1**2 + u2**2 + u3**2)
u = u0 * g0 + u1 * g[1] + u2 * g[2] + u3 * g[3]
assert sp.simplify((u * u).scalar() - 1) == 0
Rb = (1 + u * g0) / sp.sqrt(2 * (1 + u0))
assert eq(Rb * Rb.rev(), 1)
assert eq(Rb * g0 * Rb.rev(), u)
print("   ok")

print("6. commutant of gamma_0 in the even part")
eqs = list((Z * g0 - g0 * Z).d.values())
sol = sp.solve(eqs, cs, dict=True)[0]
Zs = Z.map(lambda v: v.subs(sol))
free = sorted(set().union(*[v.free_symbols for v in Zs.d.values()]), key=str)
assert free == [cs[0], cs[4], cs[5], cs[6]]  # 1, omega sigma_k
print("   span{1, omega sigma_k}: ok")

print("7. decomposition psi = sqrt(rho) e^{omega beta/2} R")
for seed in range(3):
    psi = rnd_even(seed + 7)
    pp = psi * psi.rev()
    a, b = pp.scalar(), -(pp * w).scalar()  # pp = a + b omega, omega^2 = -1
    rho = sp.sqrt(a**2 + b**2)
    beta = sp.atan2(b, a)
    Rm = (sp.cos(beta / 2) - w * sp.sin(beta / 2)) * psi / sp.sqrt(rho)
    assert eq(Rm * Rm.rev(), 1)
    assert eq(pp, rho * (sp.cos(beta) + w * sp.sin(beta)))
    for m in range(4):
        e = Rm * g[m] * Rm.rev()
        assert eq(psi * g[m] * psi.rev(), rho * e)
        assert eq(e, e.grade(1))
    L = R.map(lambda v: v.subs({eta: sp.Rational(1, 3), th: sp.Rational(2, 7)}))
    psi2 = L * psi
    assert eq(psi2 * psi2.rev(), pp)
    for m in range(4):
        assert eq(psi2 * g[m] * psi2.rev(), L * (psi * g[m] * psi.rev()) * L.rev())
print("   ok")

print("8. supplements")
# reversion of 1, omega sigma_k agrees with Cl_{3,0}; of sigma_k it does not
for k in range(3):
    assert eq((w * sg[k]).rev(), -(w * sg[k])) and eq(sg[k].rev(), -sg[k])
# for psi commuting with gamma_0: psi gamma_0 psi~ = rho gamma_0, psi sigma_3 psi~ = (psi gamma_3 psi~) gamma_0
cp = sp.symbols("d0:4", real=True)
pp = cp[0] + sum((cp[k + 1] * w * sg[k] for k in range(3)), S.zero())
rho = (pp * pp.rev()).scalar()
assert eq(pp * g0 * pp.rev(), rho * g0)
assert eq(pp * sg[2] * pp.rev(), (pp * g[3] * pp.rev()) * g0)
# the boost rotor is (1 + u0 + u_vec)/N: Hermitian positive definite matrix
Mb = mat(Rb)
assert meq(Mb, Mb.H)
assert all(sp.N(ev) > 0 for ev in Mb.eigenvals())
# psi = sigma_1: psi psi~ = -1 (beta = pi)
assert eq(sg[0] * sg[0].rev(), -1)
# right multiplication by e^{omega sigma_3 alpha} fixes gamma_0, gamma_3
al = sp.Symbol("alpha", real=True)
Sa = sp.cos(al) + w * sg[2] * sp.sin(al)
assert eq(Sa * g0 * Sa.rev(), g0) and eq(Sa * g[3] * Sa.rev(), g[3])
# e_3 orthogonal to e_0
for seed in range(3):
    psi = rnd_even(seed + 20)
    assert ((psi * g0 * psi.rev()) * (psi * g[3] * psi.rev()) + (psi * g[3] * psi.rev()) * (psi * g0 * psi.rev())).scalar() == 0
print("   ok")
print("all ok")
