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
5. Wave fronts xi = const move with speed c along k; theta = kappa xi =
   omega t - k.x with omega = c kappa, f lambda = c; E x (k x E) = |E|^2 k;
   u = eps_0|E|^2/2 + |B|^2/(2 mu_0).
6. Parallel-plate capacitor: charging work Q^2 d/(2 eps_0 S) = (eps_0/2)E^2 S d.
7. Sunlight 1.4e3 W/m^2: <cos^2> = 1/2 gives E_0 ~ 1.0e3 V/m, B_0 ~ 3.4e-6 T.
8. Duality: Dq(F e^{I a}) = (Dq F) e^{I a}, F I = -cB + I E, (F e^{Ia})(F e^{Ia})^dagger = F F^dagger,
   (F e^{Ia})^2 = F^2 e^{2Ia}; u^2 - |S|^2/c^2 = (eps_0/2)^2((|E|^2 - c^2|B|^2)^2 + 4c^2(E.B)^2);
   1+1 dimensions: P_+ f(x_0 - x_1) + P_- g(x_0 + x_1) solves (d_0 + e_1 d_1)F = 0.
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

# ---------------------------------------------------------------- 5.
# Wave fronts: F(xi) at time t + dt and position x + c dt k equals F(xi) at (t, x).
tt, dt, om = sp.symbols("t dt omega", real=True)
xi_t = c * tt - dot(kc, Xs)
shifted = xi_t.subs({tt: tt + dt, **{Xs[i]: Xs[i] + c * dt * kc[i] for i in range(3)}})
assert sp.simplify(sp.trigsimp(shifted - xi_t)) == 0
# theta = kappa xi = omega t - k.x with omega = c kappa; f lambda = c
assert sp.expand(kap * xi_t - (c * kap * tt - dot([kap * v for v in kc], Xs))) == 0
fr, lam = c * kap / (2 * sp.pi), 2 * sp.pi / kap
assert sp.simplify(fr * lam - c) == 0
# E, B, k right-handed: E x (k x E) = |E|^2 k for E perp k
EkB = cross(Eperp, cross(kc, Eperp))
for i in range(3):
    assert sp.simplify(sp.trigsimp(sp.expand(EkB[i] - dot(Eperp, Eperp) * kc[i]))) == 0
# u = (eps0/2)|E|^2 + |B|^2/(2 mu0)
assert sp.simplify(u - (eps0 / 2 * dot(Ec, Ec) + dot(Bc, Bc) / (2 * mu0))) == 0
print("5. xi-planes move with speed c along k; omega = c kappa, f lambda = c; E x cB parallel to k; "
      "u = eps0|E|^2/2 + |B|^2/(2 mu0)")

# ---------------------------------------------------------------- 6.
# Parallel-plate capacitor: charging work int_0^Q (q d/(eps0 S)) dq = (eps0/2) E^2 S d
qv, Q, Sa, dd = sp.symbols("q Q S d", positive=True)
W = sp.integrate(qv * dd / (eps0 * Sa), (qv, 0, Q))
Ecap = Q / (eps0 * Sa)
assert sp.simplify(W - Q**2 * dd / (2 * eps0 * Sa)) == 0
assert sp.simplify(W - eps0 / 2 * Ecap**2 * Sa * dd) == 0
print("6. capacitor: W = Q^2 d/(2 eps0 S) = (eps0/2) E^2 S d")

# ---------------------------------------------------------------- 7.
# Sunlight: <cos^2> = 1/2; c eps0 E0^2/2 = 1.4e3 W/m^2 gives E0 ~ 1.0e3 V/m, B0 = E0/c ~ 3.4e-6 T
ph = sp.Symbol("phi", real=True)
assert sp.integrate(sp.cos(ph)**2, (ph, 0, 2 * sp.pi)) / (2 * sp.pi) == sp.Rational(1, 2)
eps0_n, c_n = 8.8541878188e-12, 299792458.0
E0 = (2 * 1.4e3 / (c_n * eps0_n)) ** 0.5
assert round(E0, -2) == 1.0e3 and round(E0 / c_n * 1e6, 1) == 3.4
print(f"7. sunlight 1.4e3 W/m^2: E0 = {E0:.0f} V/m, B0 = {E0 / c_n:.2e} T")

# ---------------------------------------------------------------- 8.
alpha = sp.Symbol("alpha", real=True)
eIa = sp.cos(alpha) + I * sp.sin(alpha)
Fd = F * eIa
assert eq(Dq(Fd), Dq(F) * eIa)
assert eq(Fc * I, -c * vec(Bc) + I * vec(Ec))
assert eq((Fc * eIa) * (Fc * eIa).rev(), Fc * Fc.rev())
assert eq((Fc * eIa) * (Fc * eIa), Fc * Fc * (sp.cos(2 * alpha) + I * sp.sin(2 * alpha)))
uu = eps0 / 2 * (dot(Ec, Ec) + c**2 * dot(Bc, Bc))
SS = [w / mu0 for w in cross(Ec, Bc)]
assert sp.simplify(uu**2 - dot(SS, SS) / c**2
                   - (eps0 / 2)**2 * ((dot(Ec, Ec) - c**2 * dot(Bc, Bc))**2 + 4 * c**2 * dot(Ec, Bc)**2)) == 0
fa, ga = sp.Function("fa"), sp.Function("ga")
x1 = Xs[0]
G1 = (1 + e[0]) / 2 * fa(x0 - x1) + (1 - e[0]) / 2 * ga(x0 + x1)
assert eq(d(G1, x0) + e[0] * d(G1, x1), 0)
print("8. duality F e^{I a}; u^2 - |S|^2/c^2 = (eps0/2)^2 |F^2|^2 >= 0; 1+1 d'Alembert solution with P_+-")
