"""Checks for em/06-potential.md (potential, gauge, energy-momentum in Cl_{1,3}).

1. A = phi gamma_0 + c sum A_k gamma_k: gamma_0 A = phi - cA (the P of 01).
   D ^ A = E + IcB with E = -grad phi - d_t A, B = curl A; D . A = d_0 phi + c div A.
2. Gauge: D ^ (D chi) = 0, D . (D chi) = d_0^2 chi - Lap chi.
3. D ^ (D ^ A) = 0; D (D ^ A) = D^2 A - D(D . A); in the Lorenz gauge F = DA and
   D^2 A = mu_0 c J gives (d_0^2 - Lap) phi = rho/eps_0, (d_0^2 - Lap) A = mu_0 J.
4. Static fields with div A = 0: F = DA = -D_3 P (the formula of 01).
5. Energy-momentum: T(a) = -(eps_0/2) F a F is a vector, a . T(b) = b . T(a),
   T(gamma_0) gamma_0 = (eps_0/2) F F^dagger = u + S/c (03) with
   F^dagger = gamma_0 F~ gamma_0, and for F = D ^ A
   with J defined by D . F = mu_0 c J: sum_mu d_mu T(gamma^mu) = -(1/c) F . J,
   where (1/c) F . J = (1/c) J . E gamma_0 + (rho E + J x B)_k gamma_k.
6. Maxwell stress: T^{ij} = gamma^i . T(gamma^j) = -tau_ij, and the momentum balance
   d_t(S_i/c^2) - sum_j d_j tau_ij = -(rho E + J x B)_i.
"""

import sympy as sp

from common.clifford import MV, Alg, eq, mv

c, eps0 = sp.symbols("c epsilon_0", positive=True)
mu0 = 1 / (eps0 * c**2)
S = Alg(1, 3)
g = S.e
g0 = g[0]
X = S.X
x0 = X[0]
Xs = X[1:]
sg = [g[k] * g0 for k in (1, 2, 3)]
I = S.I


def svec(V):
    return sum((V[k] * sg[k] for k in range(3)), S.zero())


def curl(V):
    return [sp.diff(V[2], Xs[1]) - sp.diff(V[1], Xs[2]),
            sp.diff(V[0], Xs[2]) - sp.diff(V[2], Xs[0]),
            sp.diff(V[1], Xs[0]) - sp.diff(V[0], Xs[1])]


def div(V):
    return sum(sp.diff(V[k], Xs[k]) for k in range(3))


def box(f):
    return sp.diff(f, x0, 2) - sum(sp.diff(f, v, 2) for v in Xs)


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def comp(M, b):
    return M.d.get(b, 0)


phi = sp.Function("phi")(*X)
Af = [sp.Function(f"A{k}")(*X) for k in (1, 2, 3)]
A = phi * g0 + c * sum((Af[k] * g[k + 1] for k in range(3)), S.zero())

# ---------------------------------------------------------------- 1.
assert eq(g0 * A, phi - c * svec(Af))
DA = S.D(A)
E = [-sp.diff(phi, v) - c * sp.diff(Af[k], x0) for k, v in enumerate(Xs)]  # c d_0 = d_t
B = curl(Af)
F = svec(E) + I * c * svec(B)
assert eq(DA.grade(2), F)
assert eq(DA.grade(0), sp.diff(phi, x0) + c * div(Af))
assert eq(DA.grade(4), 0)
print("1. gamma_0 A = phi - cA; D^A = E + IcB (E = -grad phi - d_t A, B = curl A); D.A = d0 phi + c div A")

# ---------------------------------------------------------------- 2.
chi = sp.Function("chi")(*X)
Dchi = S.D(mv(chi, S.neg))
assert eq(S.D(Dchi).grade(2), 0)
assert eq(S.D(Dchi).grade(0), box(chi))
print("2. D^(D chi) = 0, D.(D chi) = box chi")

# ---------------------------------------------------------------- 3.
DF = S.D(F)
assert eq(DF.grade(3), 0)
assert eq(DF, S.D(S.D(A)) - S.D(DA.grade(0)))
rho = div(E) * eps0
J = [(curl(B)[k] - sp.diff(E[k], x0) / c) * eps0 * c**2 for k in range(3)]
J4 = c * rho * g0 + sum((J[k] * g[k + 1] for k in range(3)), S.zero())
assert eq(DF.grade(1), mu0 * c * J4)
# Lorenz gauge: impose d_0 phi = -c div A by construction of the check below
lor = sp.diff(phi, x0) + c * div(Af)
# D^2 A = box A componentwise; D F = box A - D(D.A)
assert eq(S.D(S.D(A)), box(phi) * g0 + c * sum((box(Af[k]) * g[k + 1] for k in range(3)), S.zero()))
# so with D.A = 0: box phi = gamma_0 component of mu0 c J = mu0 c^2 rho = rho/eps0
resid = DF.grade(1) - (S.D(S.D(A)) - S.D(mv(lor, S.neg)))
assert eq(resid, 0)
assert sp.simplify(mu0 * c * c * rho - rho / eps0) == 0
print("3. D^(D^A) = 0, D(D^A) = D^2 A - D(D.A) = mu0 c J; Lorenz gauge: box phi = rho/eps0, box A = mu0 J")

# ---------------------------------------------------------------- 4.
phis = sp.Function("phi")(*Xs)
Afs = [sp.Function(f"A{k}")(*Xs) for k in (1, 2, 3)]
As = phis * g0 + c * sum((Afs[k] * g[k + 1] for k in range(3)), S.zero())
Ps = phis - c * svec(Afs)
D3 = lambda H: sum((sg[k] * S.d(H, k + 1) for k in range(3)), S.zero())
assert eq(S.D(As), -D3(Ps))
assert eq(S.D(As).grade(0), c * sum(sp.diff(Afs[k], Xs[k]) for k in range(3)))
print("4. static: DA = -D_3 P, scalar part c div A")

# ---------------------------------------------------------------- 5.
Ec, Bc = sp.symbols("E1:4", real=True), sp.symbols("B1:4", real=True)
Fc = svec(Ec) + I * c * svec(Bc)
T = lambda a, FF=Fc: -(eps0 / 2) * FF * a * FF
av = sp.symbols("a0:4", real=True)
bv = sp.symbols("b0:4", real=True)
va = sum((av[m] * g[m] for m in range(4)), S.zero())
vb = sum((bv[m] * g[m] for m in range(4)), S.zero())
assert eq(T(va), T(va).grade(1))
ip = lambda p, q: (p * q + q * p) / 2
assert eq(ip(va, T(vb)), ip(vb, T(va)))
dot = lambda p, q: sum(p[k] * q[k] for k in range(3))
u = eps0 / 2 * (dot(Ec, Ec) + c**2 * dot(Bc, Bc))
Sp = [w / mu0 for w in cross(Ec, Bc)]
# the Cl_{3,0} reversion F^dagger = E - IcB is gamma_0 F~ gamma_0 in Cl_{1,3}
Fdag = g0 * Fc.rev() * g0
assert eq(Fdag, svec(Ec) - I * c * svec(Bc))
assert eq(T(g0) * g0, eps0 / 2 * Fc * Fdag)
assert eq(T(g0) * g0, u + svec(Sp) / c)
# divergence
divT = sum((S.d(T(S.er[m], F), m) for m in range(4)), S.zero())
FJ = (F * J4 - J4 * F) / 2
assert eq(divT, -FJ / c)
Er = [sp.Symbol(f"E{k}r") for k in range(3)]
Br = [sp.Symbol(f"B{k}r") for k in range(3)]
Jr = sp.symbols("J1:4")
rr = sp.Symbol("rho_r")
Fr = svec(Er) + I * c * svec(Br)
J4r = c * rr * g0 + sum((Jr[k] * g[k + 1] for k in range(3)), S.zero())
f = ((Fr * J4r - J4r * Fr) / 2) / c
assert sp.simplify(comp(f, 1) - dot(Jr, Er) / c) == 0
for k in range(3):
    assert sp.simplify(comp(f, 1 << (k + 1)) - (rr * Er[k] + cross(Jr, Br)[k])) == 0
print("5. T(a) = -(eps0/2) F a F: vector, symmetric, T(gamma_0)gamma_0 = (eps0/2)FF^dagger = u + S/c; "
      "sum d_mu T(gamma^mu) = -(1/c)F.J, (1/c)F.J = (J.E/c) gamma_0 + (rho E + J x B)_k gamma_k")

# ---------------------------------------------------------------- 6.
# Maxwell stress tau_ij = eps0 (E_i E_j + c^2 B_i B_j - delta_ij (|E|^2 + c^2|B|^2)/2):
# T^{ij} = gamma^i . T(gamma^j) = -tau_ij, and for fields from potentials
# d_t(S_i/c^2) - sum_j d_j tau_ij = -(rho E + J x B)_i.
E2c, B2c = dot(Ec, Ec), dot(Bc, Bc)
tau = lambda EE, BB, i, j: eps0 * (EE[i] * EE[j] + c**2 * BB[i] * BB[j]
                                   - sp.Rational(1, 2) * (1 if i == j else 0) * (dot(EE, EE) + c**2 * dot(BB, BB)))
ip0 = lambda p, q: ((p * q + q * p) / 2).scalar()
for i in range(3):
    for j in range(3):
        assert sp.simplify(ip0(S.er[i + 1], T(S.er[j + 1])) + tau(Ec, Bc, i, j)) == 0
Sf = [w / mu0 for w in cross(E, B)]
for i in range(3):
    lhs = sp.diff(Sf[i] / c**2, x0) * c - sum(sp.diff(tau(E, B, i, j), Xs[j]) for j in range(3))
    force = rho * E[i] + cross(J, B)[i]
    assert sp.simplify(lhs + force) == 0
print("6. T^{ij} = -tau_ij; d_t(S/c^2)_i - d_j tau_ij = -(rho E + J x B)_i")
