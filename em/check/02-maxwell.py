"""Checks for em/02-maxwell.md (Maxwell's equations as one equation in Cl_{3,0}).

1. Dq = d_0 + D (x_0 = ct): the grades of Dq(E + IcB) are div E,
   d_0 E - c curl B, I(curl E + c d_0 B), Ic div B; with the SI Maxwell equations
   the equation Dq F = (rho - J/c)/eps_0 is equivalent to the four equations.
2. Coefficient squares: for generators a_l with a_l a_m = -a_m a_l (l != m),
   (d_0 - sum a_l d_l)(d_0 + sum a_l d_l) = d_0^2 - sum a_l^2 d_l^2 on generic functions:
   Cl_{3,0} (a_l = e_l, square +1) gives d_0^2 - Laplacian, Cl_{4,0} with
   h_l = e_0 e_l (square -1) gives the 4-dim Laplacian. Dqbar Dq = Dq Dqbar.
3. Dqbar applied to the source: scalar part (d_t rho + div J)/(c eps_0)
   (continuity), vector part -(grad rho)/eps_0 - mu_0 d_t J, bivector part
   I mu_0 c curl J; so (d_0^2 - Lap)E = -(grad rho)/eps_0 - mu_0 d_t J,
   (d_0^2 - Lap)B = mu_0 curl J for fields satisfying Maxwell's equations.
4. Static case: d_0 = 0 reduces Dq F to D F.
"""

import sympy as sp

from common.clifford import MV, Alg, eq, mv

c, eps0 = sp.symbols("c epsilon_0", positive=True)
mu0 = 1 / (eps0 * c**2)
P = Alg(3)
e, I = P.e, P.I
x0 = sp.Symbol("x0", real=True)
Xs = sp.symbols("x1:4", real=True)
X = (x0,) + tuple(Xs)


def vec(V):
    return sum((V[k] * e[k] for k in range(3)), P.zero())


def d(F, x):
    return F.map(lambda v: sp.diff(v, x))


def Dsp(F):
    return sum((e[k] * d(F, Xs[k]) for k in range(3)), P.zero())


def Dq(F):
    return d(F, x0) + Dsp(F)


def Dqbar(F):
    return d(F, x0) - Dsp(F)


def curl(V):
    return [sp.diff(V[2], Xs[1]) - sp.diff(V[1], Xs[2]),
            sp.diff(V[0], Xs[2]) - sp.diff(V[2], Xs[0]),
            sp.diff(V[1], Xs[0]) - sp.diff(V[0], Xs[1])]


def div(V):
    return sum(sp.diff(V[k], Xs[k]) for k in range(3))


def grad(f):
    return [sp.diff(f, v) for v in Xs]


def lapw(f):  # d_0^2 - Laplacian
    return sp.diff(f, x0, 2) - sum(sp.diff(f, v, 2) for v in Xs)


# ---------------------------------------------------------------- 1.
Ef = [sp.Function(f"E{k}")(*X) for k in (1, 2, 3)]
Bf = [sp.Function(f"B{k}")(*X) for k in (1, 2, 3)]
F = vec(Ef) + I * c * vec(Bf)
G = Dq(F)
cE, cB = curl(Ef), curl(Bf)
assert eq(G.grade(0), div(Ef))
assert eq(G.grade(1), vec([sp.diff(Ef[k], x0) - c * cB[k] for k in range(3)]))
assert eq(G.grade(2), I * vec([cE[k] + c * sp.diff(Bf[k], x0) for k in range(3)]))
assert eq(G.grade(3), I * c * div(Bf))
# d_0 = d_t / c: vector part (1/c) d_t E - c curl B = -J/(c eps0) <=> curl B = mu0 J + (1/c^2) d_t E
Jv = sp.symbols("J1:4")
dtE, curlB = sp.symbols("dtE1:4"), sp.symbols("curlB1:4")
for k in range(3):
    lhs = dtE[k] / c - c * curlB[k] + Jv[k] / (c * eps0)
    rhs = -c * (curlB[k] - mu0 * Jv[k] - dtE[k] / c**2)
    assert sp.simplify(lhs - rhs) == 0
print("1. Dq(E + IcB) = div E + (d0 E - c curl B) + I(curl E + c d0 B) + Ic div B; "
      "vector part = -J/(c eps0) <=> curl B = mu0 J + (1/c^2) d_t E")

# ---------------------------------------------------------------- 2.
f = sp.Function("f")(*X)
assert eq(Dqbar(Dq(mv(f))), lapw(f)) and eq(Dq(Dqbar(mv(f))), lapw(f))
H = sum((sp.Function(f"h{b}")(*X) * MV({b: 1}) for b in range(8)), P.zero())
assert eq(Dqbar(Dq(H)), H.map(lapw))
C4 = Alg(4)
h = [C4.e[0] * C4.e[l] for l in (1, 2, 3)]
Y = C4.X
g = sp.Function("g")(*Y)
Dh = lambda F: C4.d(F, 0) + sum((h[l] * C4.d(F, l + 1) for l in range(3)), C4.zero())
Dhbar = lambda F: C4.d(F, 0) - sum((h[l] * C4.d(F, l + 1) for l in range(3)), C4.zero())
assert eq(Dhbar(Dh(mv(g))), sum(sp.diff(g, v, 2) for v in Y))
print("2. (d0 - a.d)(d0 + a.d) = d0^2 - sum a_l^2 d_l^2: e_l (+1) -> d0^2 - Laplacian, "
      "h_l = e0 e_l (-1) -> 4-dim Laplacian")

# ---------------------------------------------------------------- 3.
rho = div(Ef) * eps0
J = [(cB[k] - sp.diff(Ef[k], x0) / c) * eps0 * c**2 for k in range(3)]
src = (rho - vec(J) / c) / eps0
Ssrc = Dqbar(src)
assert eq(Ssrc.grade(0), (sp.diff(rho, x0) + div(J) / c) / eps0)
assert eq(Ssrc.grade(0), 0)
# general sources
rf = sp.Function("rho")(*X)
Jf = [sp.Function(f"J{k}")(*X) for k in (1, 2, 3)]
Sg = Dqbar((rf - vec(Jf) / c) / eps0)
# d_t = c d_0
assert eq(Sg.grade(0), (c * sp.diff(rf, x0) + div(Jf)) / (c * eps0))
assert eq(Sg.grade(1), vec([-sp.diff(rf, v) / eps0 - mu0 * c * sp.diff(Jf[k], x0)
                            for k, v in enumerate(Xs)]))
assert eq(Sg.grade(2), I * mu0 * c * vec(curl(Jf)))
assert eq(Sg.grade(3), 0)
# wave equations for fields satisfying the homogeneous equations (derived from 1. with
# Faraday and div B = 0 imposed through B = curl A, E = -grad phi - d_t A)
phi = sp.Function("phi")(*X)
A = [sp.Function(f"A{k}")(*X) for k in (1, 2, 3)]
E2 = [-sp.diff(phi, v) - c * sp.diff(A[k], x0) for k, v in enumerate(Xs)]
B2 = curl(A)
F2 = vec(E2) + I * c * vec(B2)
rho2 = div(E2) * eps0
J2 = [(curl(B2)[k] - sp.diff(E2[k], x0) / c) * eps0 * c**2 for k in range(3)]
assert eq(Dq(F2), (rho2 - vec(J2) / c) / eps0)
for k in range(3):
    assert sp.simplify(lapw(E2[k]) - (-sp.diff(rho2, Xs[k]) / eps0 - mu0 * c * sp.diff(J2[k], x0))) == 0
    assert sp.simplify(lapw(B2[k]) - mu0 * curl(J2)[k]) == 0
print("3. Dqbar(rho - J/c)/eps0: scalar (d_t rho + div J)/(c eps0), vector -(grad rho)/eps0 - mu0 d_t J, "
      "bivector I mu0 c curl J; wave equations for E and B")

# ---------------------------------------------------------------- 4.
Es = [sp.Function(f"E{k}")(*Xs) for k in (1, 2, 3)]
Bs = [sp.Function(f"B{k}")(*Xs) for k in (1, 2, 3)]
Fs = vec(Es) + I * c * vec(Bs)
assert eq(Dq(Fs), Dsp(Fs))
print("4. static fields: Dq F = D F")
