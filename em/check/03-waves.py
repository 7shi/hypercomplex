"""Checks for em/03-waves.md (plane waves and energy flow in Cl_{3,0}).

1. Plane waves: for F = F(xi), xi = x_0 - k.x (|k| = 1),
   Dq F = (1 - k) F'. For F = E + IcB, (1 - k)F = 0 iff k.E = 0 and cB = k x E,
   and then F = (1 + k)E. Projectors P_pm = (1 pm k)/2.
2. F^2 = |E|^2 - c^2|B|^2 + 2Ic E.B; for F = (1 + k)E (E perp k), F^2 = 0.
3. Pseudoscalar exponential: F = (1 + k)E_0 e^{I theta}, theta = kappa xi, satisfies
   Dq F = 0 and has E = E_0 cos(theta) - (k x E_0) sin(theta) (circular),
   cB = k x E; (1 + k)E_0 e^{I theta} = (1 + k)E_0 e^{-I k theta}; the average of
   e^{I theta} and e^{-I theta} gives linear polarization.
4. Energy: (eps_0/2) F F^dagger = u + S/c with u = (eps_0/2)(|E|^2 + c^2|B|^2),
   S = E x B / mu_0; <Dq(F F^dagger)>_0 = 2 <(Dq F) F^dagger>_0, hence
   d_t u + div S = -J.E. For the plane wave F F^dagger = 2|E|^2 (1 + k).
"""

import sympy as sp

from common.clifford import MV, Alg, mv


def eq(A, B):
    # the unit vector k is parametrized by angles, so trigsimp is needed
    return all(sp.simplify(sp.trigsimp(sp.expand(v))) == 0 for v in (mv(A) - mv(B)).d.values())

c, eps0 = sp.symbols("c epsilon_0", positive=True)
mu0 = 1 / (eps0 * c**2)
P = Alg(3)
e, I = P.e, P.I
x0 = sp.Symbol("x0", real=True)
Xs = sp.symbols("x1:4", real=True)
X = (x0,) + Xs


def vec(V):
    return sum((V[k] * e[k] for k in range(3)), P.zero())


def comps(M):
    return [M.d.get(1 << k, 0) for k in range(3)]


def d(F, x):
    return F.map(lambda v: sp.diff(v, x))


def Dq(F):
    return d(F, x0) + sum((e[k] * d(F, Xs[k]) for k in range(3)), P.zero())


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return sum(a[k] * b[k] for k in range(3))


def div(V):
    return sum(sp.diff(V[k], Xs[k]) for k in range(3))


def curl(V):
    return [sp.diff(V[2], Xs[1]) - sp.diff(V[1], Xs[2]),
            sp.diff(V[0], Xs[2]) - sp.diff(V[2], Xs[0]),
            sp.diff(V[1], Xs[0]) - sp.diff(V[0], Xs[1])]


def rev(F):
    return F.rev()


# unit vector k via angles
th, ph = sp.symbols("theta_k phi_k", real=True)
kc = [sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)]
k = vec(kc)
assert eq(k * k, 1)
# orthonormal frame perpendicular to k
u1 = [sp.cos(th) * sp.cos(ph), sp.cos(th) * sp.sin(ph), -sp.sin(th)]
u2 = cross(kc, u1)
assert sp.simplify(dot(kc, u1)) == 0 and sp.simplify(dot(u1, u1) - 1) == 0

# ---------------------------------------------------------------- 1.
xi = x0 - dot(kc, Xs)
s = sp.Symbol("s")
Fp = [sp.Function(f"f{b}") for b in range(8)]
Fgen = sum((Fp[b](xi) * MV({b: 1}) for b in range(8)), P.zero())
Fder = sum((sp.diff(Fp[b](s), s).subs(s, xi) * MV({b: 1}) for b in range(8)), P.zero())
assert eq(Dq(Fgen), (1 - k) * Fder)
Pp, Pm = (1 + k) / 2, (1 - k) / 2
assert eq(Pp * Pp, Pp) and eq(Pm * Pm, Pm) and eq(Pp * Pm, 0) and eq(Pp + Pm, 1)
# (1 - k)(E + IcB) = 0 with E perpendicular to k and cB = k x E
a1, a2 = sp.symbols("a1 a2", real=True)
Eperp = [a1 * u1[i] + a2 * u2[i] for i in range(3)]
Fw = vec(Eperp) + I * vec(cross(kc, Eperp))
assert eq((1 - k) * Fw, 0) and eq(Fw, (1 + k) * vec(Eperp))
# conversely: solve (1 - k)(E + IcB) = 0 for k = e_3
Ev, cBv = sp.symbols("E1:4"), sp.symbols("cB1:4")
sol = sp.solve(list(((1 - e[2]) * (vec(Ev) + I * vec(cBv))).d.values()), list(Ev) + list(cBv), dict=True)
assert len(sol) == 1
sub = lambda v: sp.simplify(sp.sympify(v).subs(sol[0]))
assert sub(Ev[2]) == 0 and sub(cBv[2]) == 0
assert sub(cBv[0] + Ev[1]) == 0 and sub(cBv[1] - Ev[0]) == 0  # cB = e_3 x E
print("1. Dq F(x0 - k.x) = (1 - k)F'; (1 - k)(E + IcB) = 0 iff k.E = 0, cB = k x E; F = (1 + k)E")

# ---------------------------------------------------------------- 2.
Ec, Bc = sp.symbols("E1:4", real=True), sp.symbols("B1:4", real=True)
Fc = vec(Ec) + I * c * vec(Bc)
assert eq(Fc * Fc, dot(Ec, Ec) - c**2 * dot(Bc, Bc) + 2 * I * c * dot(Ec, Bc))
assert eq(Fw * Fw, 0)
print("2. F^2 = |E|^2 - c^2|B|^2 + 2Ic E.B; ((1 + k)E)^2 = 0 for E perp k")

# ---------------------------------------------------------------- 3.
kap = sp.Symbol("kappa", positive=True)
theta = kap * xi
expI = lambda a: sp.cos(a) + I * sp.sin(a)
E0 = vec(Eperp)
Fe = (1 + k) * E0 * expI(theta)
assert eq(Dq(Fe), 0)
Ee = comps(Fe.grade(1))
kxE0 = cross(kc, Eperp)
assert all(sp.simplify(Ee[i] - (Eperp[i] * sp.cos(theta) - kxE0[i] * sp.sin(theta))) == 0 for i in range(3))
cBe = comps((Fe.grade(2) * I).map(lambda v: -v))
assert all(sp.simplify(cBe[i] - cross(kc, Ee)[i]) == 0 for i in range(3))
Ik = I * k
assert eq(Ik * Ik, -1)
expIk = lambda a: sp.cos(a) + Ik * sp.sin(a)
assert eq(Fe, (1 + k) * E0 * expIk(-theta))
assert eq(Fe, (1 + k) * expIk(theta) * E0)
Fl = ((1 + k) * E0 * expI(theta) + (1 + k) * E0 * expI(-theta)) / 2
assert eq(Fl, (1 + k) * E0 * sp.cos(theta))
# |E| constant along the wave (circular)
assert sp.simplify(dot(Ee, Ee) - dot(Eperp, Eperp)) == 0
print("3. (1 + k)E0 e^{I theta}: Dq F = 0, E = E0 cos - (k x E0) sin, |E| const, "
      "= (1 + k)E0 e^{-Ik theta}; average with e^{-I theta} is linear")

# ---------------------------------------------------------------- 4.
T = eps0 / 2 * Fc * rev(Fc)
u = eps0 / 2 * (dot(Ec, Ec) + c**2 * dot(Bc, Bc))
Sp = [v / mu0 for v in cross(Ec, Bc)]
assert eq(T, u + vec(Sp) / c)
Ef = [sp.Function(f"E{i}")(*X) for i in (1, 2, 3)]
Bf = [sp.Function(f"B{i}")(*X) for i in (1, 2, 3)]
F = vec(Ef) + I * c * vec(Bf)
lhs = Dq(F * rev(F)).grade(0)
assert eq(lhs, 2 * (Dq(F) * rev(F)).grade(0))
# with Dq F = (rho - J/c)/eps0: <(rho - J/c) F^dagger>_0 = -J.E/c
rho = sp.Symbol("rho")
Jv = sp.symbols("J1:4")
assert eq(((rho - vec(Jv) / c) * rev(Fc)).grade(0), -dot(Jv, Ec) / c)
# scalar part of Dq T = d_0 u + div(S)/c = (d_t u + div S)/c
uf = eps0 / 2 * (dot(Ef, Ef) + c**2 * dot(Bf, Bf))
Sf = [v / mu0 for v in cross(Ef, Bf)]
assert eq((eps0 / 2) * lhs, sp.diff(uf, x0) + div(Sf) / c)
# Poynting theorem from Maxwell's equations
rhoM = div(Ef) * eps0
JM = [(curl(Bf)[i] - sp.diff(Ef[i], x0) / c) * eps0 * c**2 for i in range(3)]
# impose the homogeneous equations: use E = -grad phi - d_t A, B = curl A
phi = sp.Function("phi")(*X)
A = [sp.Function(f"A{i}")(*X) for i in (1, 2, 3)]
E2 = [-sp.diff(phi, v) - c * sp.diff(A[i], x0) for i, v in enumerate(Xs)]
B2 = curl(A)
J2 = [(curl(B2)[i] - sp.diff(E2[i], x0) / c) * eps0 * c**2 for i in range(3)]
u2 = eps0 / 2 * (dot(E2, E2) + c**2 * dot(B2, B2))
S2 = [v / mu0 for v in cross(E2, B2)]
assert sp.simplify(c * sp.diff(u2, x0) + div(S2) + dot(J2, E2)) == 0
Tw = Fw * rev(Fw)
assert eq(Tw, 2 * dot(Eperp, Eperp) * (1 + k))
print("4. (eps0/2)FF^dagger = u + S/c; <Dq(FF^dagger)>_0 = 2<(Dq F)F^dagger>_0; "
      "d_t u + div S = -J.E; plane wave FF^dagger = 2|E|^2(1 + k)")
