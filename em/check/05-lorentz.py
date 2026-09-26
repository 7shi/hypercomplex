"""Checks for em/05-lorentz.md (Lorentz transformations as rotors in Cl_{1,3}).

1. Exponentials: sigma_1^2 = +1 gives e^{sigma_1 a} = cosh a + sigma_1 sinh a,
   (I sigma_3)^2 = -1 gives e^{I sigma_3 a} = cos a + I sigma_3 sin a (power series
   truncated and compared via the recurrence of the square).
2. Boost: R = e^{sigma_1 eta/2}, R R~ = 1, R gamma_0 R~ = cosh(eta) gamma_0 + sinh(eta) gamma_1,
   R gamma_1 R~ = sinh(eta) gamma_0 + cosh(eta) gamma_1, gamma_2, gamma_3 fixed; x^2 invariant.
   Rapidities add: e^{sigma_1 a/2} e^{sigma_1 b/2} = e^{sigma_1 (a+b)/2}.
3. Rotation: R = e^{-I sigma_3 theta/2} rotates gamma_1 -> cos gamma_1 + sin gamma_2,
   fixes gamma_0 and gamma_3; it commutes with gamma_0 and acts on sigma_k as a rotation.
4. Fields seen by an observer with gamma_0' = R gamma_0 R~ (velocity v e_1,
   v/c = tanh eta): F' = R~ F R has E'_1 = E_1, E'_2 = gamma(E_2 - v B_3),
   E'_3 = gamma(E_3 + v B_2), B'_1 = B_1, B'_2 = gamma(B_2 + v E_3/c^2),
   B'_3 = gamma(B_3 - v E_2/c^2). Equivalently F = sum E'_k sigma'_k + Ic sum B'_k sigma'_k
   with sigma'_k = R sigma_k R~.
5. Invariants: F^2 = (|E|^2 - c^2|B|^2) + 2Ic E.B commutes with R, so F'^2 = F^2.
6. Lorentz force: v = gamma(c gamma_0 + u_k gamma_k), v^2 = c^2; (q/c) F . v
   (F . v = (Fv - vF)/2) has gamma_0 component gamma q E.u / c and gamma_k components
   gamma q (E + u x B)_k; v . (F . v) = 0.
7. Doppler: for the plane wave F = (1 + sigma_1) E_perp f(x_0 - x_1) and the boost
   along e_1, R~ F R = e^{-eta} F and x_0 - x_1 = e^{-eta}(x'_0 - x'_1).
   The observer coordinates are x'_0 = gamma(x_0 - beta x_1), x'_1 = gamma(x_1 - beta x_0)
   with gamma = cosh(eta), beta = tanh(eta) (the Lorentz transformation quoted in em/04).
8. dx_1/dt = c tanh(eta) along gamma_0', t = cosh(eta) tau (time dilation);
   velocity addition (v1 + v2)/(1 + v1 v2/c^2) = c tanh(a + b), c with c, below c.
9. Boosting a pure electric field: B' = u x E'/c^2 with u = -v; Biot-Savart with
   I dl -> q u equals u x E_Coulomb/c^2.
10. J = c rho gamma_0 seen by the boosted observer: rho' = gamma rho, J'_1 = -gamma rho v.
11. gamma m c^2 = m c^2 + m u^2/2 + O(u^4).
12. Intuition: gamma'_mu = R gamma_mu R~ satisfy the same relations; R(gamma_0 +- gamma_1)R~
    = e^{+-eta}(gamma_0 +- gamma_1); dR/dtau = (q/2mc)FR gives m dU/dtau = (q/c)F . U for
    U = c R gamma_0 R~; hyperbolic motion in uniform E and cyclotron motion in uniform B.
13. sigma_1 sigma_2 - sigma_2 sigma_1 = 2 I sigma_3; e^{sigma_1 a/2} e^{sigma_2 b/2} has the rotation
    component I sigma_3 sinh(a/2) sinh(b/2). Canonical form: boosting along S with tanh 2eta = |S|/(cu)
    makes E' x B' = 0 (random integer fields, numerically).
"""

import sympy as sp

from common.clifford import MV, Alg, mv

c = sp.Symbol("c", positive=True)
S = Alg(1, 3)
g = S.e
g0 = g[0]
sg = [g[k] * g0 for k in (1, 2, 3)]
I = S.I


def eq(A, B):
    return all(sp.simplify(sp.trigsimp(sp.expand((v).rewrite(sp.exp)))) == 0 or
               sp.simplify(sp.trigsimp(sp.expand(v))) == 0
               for v in (mv(A) - mv(B)).d.values())


def svec(V):
    return sum((V[k] * sg[k] for k in range(3)), S.zero())


def rev(M):
    return M.rev()


def comp(M, b):
    return M.d.get(b, 0)


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


# ---------------------------------------------------------------- 1.
a = sp.Symbol("a", real=True)
N = 14


def expser(B, t):
    term, tot = mv(1, S.neg), mv(1, S.neg)
    for n in range(1, N):
        term = term * B * t / n
        tot = tot + term
    return tot


def trunc(expr, t):
    return sp.series(expr, t, 0, N).removeO()


E1 = expser(sg[0], a)
for b, v in (sp.cosh(a) + sg[0] * sp.sinh(a)).d.items():
    assert sp.expand(E1.d.get(b, 0) - trunc(v, a)) == 0
Is3 = I * sg[2]
assert eq(Is3 * Is3, -1) and eq(sg[0] * sg[0], 1)
E2 = expser(Is3, a)
for b, v in (sp.cos(a) + Is3 * sp.sin(a)).d.items():
    assert sp.expand(E2.d.get(b, 0) - trunc(v, a)) == 0
print("1. e^{sigma_1 a} = cosh a + sigma_1 sinh a, e^{I sigma_3 a} = cos a + I sigma_3 sin a")

# ---------------------------------------------------------------- 2.
eta = sp.Symbol("eta", real=True)
R = sp.cosh(eta / 2) + sg[0] * sp.sinh(eta / 2)
Rt = rev(R)
assert eq(Rt, sp.cosh(eta / 2) - sg[0] * sp.sinh(eta / 2))
assert eq(R * Rt, 1)
assert eq(R * g0 * Rt, sp.cosh(eta) * g0 + sp.sinh(eta) * g[1])
assert eq(R * g[1] * Rt, sp.sinh(eta) * g0 + sp.cosh(eta) * g[1])
assert eq(R * g[2] * Rt, g[2]) and eq(R * g[3] * Rt, g[3])
xp = R * S.x * Rt
assert eq(xp * xp, S.x * S.x)
bb = sp.Symbol("b", real=True)
Ra = sp.cosh(a / 2) + sg[0] * sp.sinh(a / 2)
Rb = sp.cosh(bb / 2) + sg[0] * sp.sinh(bb / 2)
assert eq(Ra * Rb, sp.cosh((a + bb) / 2) + sg[0] * sp.sinh((a + bb) / 2))
print("2. R = e^{sigma_1 eta/2}: gamma_0 -> cosh gamma_0 + sinh gamma_1, gamma_1 -> sinh gamma_0 + cosh gamma_1, "
      "x^2 invariant, rapidities add")

# ---------------------------------------------------------------- 3.
th = sp.Symbol("theta", real=True)
Rr = sp.cos(th / 2) - Is3 * sp.sin(th / 2)
assert eq(Rr * rev(Rr), 1)
assert eq(Rr * g[1] * rev(Rr), sp.cos(th) * g[1] + sp.sin(th) * g[2])
assert eq(Rr * g[2] * rev(Rr), -sp.sin(th) * g[1] + sp.cos(th) * g[2])
assert eq(Rr * g0 * rev(Rr), g0) and eq(Rr * g[3] * rev(Rr), g[3])
assert eq(Rr * sg[0] * rev(Rr), sp.cos(th) * sg[0] + sp.sin(th) * sg[1])
print("3. e^{-I sigma_3 theta/2}: gamma_1 -> cos gamma_1 + sin gamma_2, fixes gamma_0, gamma_3; same on sigma_k")

# ---------------------------------------------------------------- 4.
Ec, Bc = sp.symbols("E1:4", real=True), sp.symbols("B1:4", real=True)
F = svec(Ec) + I * c * svec(Bc)
Fp = rev(R) * F * R


def split(M):
    """E and B components of a bivector M = sum E_k sigma_k + Ic sum B_k sigma_k."""
    Ek = [sp.Symbol(f"e{k}") for k in range(3)]
    Bk = [sp.Symbol(f"f{k}") for k in range(3)]
    diff = M - (svec(Ek) + I * c * svec(Bk))
    sol = sp.solve(list(diff.d.values()), Ek + Bk, dict=True)[0]
    return [sol[v] for v in Ek], [sol[v] for v in Bk]


Ep, Bp = split(Fp)
gam = sp.cosh(eta)
vv = c * sp.tanh(eta)
exp_E = [Ec[0], gam * (Ec[1] - vv * Bc[2]), gam * (Ec[2] + vv * Bc[1])]
exp_B = [Bc[0], gam * (Bc[1] + vv * Ec[2] / c**2), gam * (Bc[2] - vv * Ec[1] / c**2)]
for k in range(3):
    assert sp.simplify((Ep[k] - exp_E[k]).rewrite(sp.exp)) == 0
    assert sp.simplify((Bp[k] - exp_B[k]).rewrite(sp.exp)) == 0
sgp = [R * s * Rt for s in sg]
assert eq(sum((Ep[k] * sgp[k] + I * c * Bp[k] * sgp[k] for k in range(3)), S.zero()), F)
assert eq(sgp[0], (R * g[1] * Rt) * (R * g0 * Rt))
print("4. F' = R~ F R: E'_par = E_par, E'_perp = gamma(E + v x B)_perp, B'_perp = gamma(B - v x E/c^2)_perp")

# ---------------------------------------------------------------- 5.
F2 = F * F
dot = lambda p, q: sum(p[k] * q[k] for k in range(3))
assert eq(F2, dot(Ec, Ec) - c**2 * dot(Bc, Bc) + 2 * I * c * dot(Ec, Bc))
assert eq(F2 * R, R * F2)
assert eq(Fp * Fp, F2)
print("5. F^2 = (|E|^2 - c^2|B|^2) + 2Ic E.B is invariant")

# ---------------------------------------------------------------- 6.
q = sp.Symbol("q", real=True)
u = sp.symbols("u1:4", real=True)
gL = 1 / sp.sqrt(1 - dot(u, u) / c**2)
vel = gL * (c * g0 + sum((u[k] * g[k + 1] for k in range(3)), S.zero()))
assert eq(vel * vel, c**2)
Fv = (F * vel - vel * F) / 2
force = q / c * Fv
assert eq(force.grade(1), force)
assert sp.simplify(comp(force, 1) - gL * q * dot(Ec, u) / c) == 0
lor = [Ec[k] + cross(u, Bc)[k] for k in range(3)]
for k in range(3):
    assert sp.simplify(comp(force, 1 << (k + 1)) - gL * q * lor[k]) == 0
vdot = ((vel * Fv + Fv * vel) / 2).grade(0)
assert eq(vdot, 0)
print("6. v = gamma(c gamma_0 + u gamma), v^2 = c^2; (q/c)F.v = gamma(q E.u/c) gamma_0 + gamma q(E + u x B)_k gamma_k; "
      "v.(F.v) = 0")

# ---------------------------------------------------------------- 7.
# Doppler: a plane wave along sigma_1, F = (1 + sigma_1) E_perp f(x_0 - x_1), seen by the
# observer boosted along e_1: R~ F R = e^{-eta} F (amplitude), and x_0 - x_1 = e^{-eta}(x'_0 - x'_1)
# with observer coordinates x'_0 = x . gamma'_0, x'_1 = -x . gamma'_1.
Ey, Ez = sp.symbols("E_y E_z", real=True)
Fw = (1 + sg[0]) * (Ey * sg[1] + Ez * sg[2])
assert eq(rev(R) * Fw * R, sp.exp(-eta) * Fw)
ip = lambda p, q: ((p * q + q * p) / 2).scalar()
g0p, g1p = R * g0 * Rt, R * g[1] * Rt
X = S.X
x0p = ip(S.x, g0p)
x1p = -ip(S.x, g1p)
assert sp.simplify((X[0] - X[1] - sp.exp(-eta) * (x0p - x1p)).rewrite(sp.exp)) == 0
assert sp.simplify(x0p - sp.cosh(eta) * (X[0] - sp.tanh(eta) * X[1])) == 0
assert sp.simplify(x1p - sp.cosh(eta) * (X[1] - sp.tanh(eta) * X[0])) == 0
print("7. plane wave along the boost: R~ F R = e^{-eta} F, x0 - x1 = e^{-eta}(x0' - x1') (Doppler factor); "
      "x0' = gamma(x0 - beta x1), x1' = gamma(x1 - beta x0)")

# ---------------------------------------------------------------- 8.
# gamma_0' = R gamma_0 R~ moves with dx_1/dt = c tanh(eta); time dilation: c tau gamma_0'
# has gamma_0 component c tau cosh(eta), i.e. t = gamma tau.
tau_ = sp.Symbol("tau", positive=True)
g0p = R * g0 * Rt
assert sp.simplify(comp(g0p, 2) / comp(g0p, 1) - sp.tanh(eta)) == 0
assert sp.simplify(comp(c * tau_ * g0p, 1) - c * tau_ * sp.cosh(eta)) == 0
assert sp.simplify(sp.cosh(eta) - 1 / sp.sqrt(1 - sp.tanh(eta)**2)) == 0
# velocity addition: (v1 + v2)/(1 + v1 v2/c^2); v2 = c gives c; |result| < c for |v1|, |v2| < c
v1, v2 = sp.symbols("v1 v2", real=True)
add = (v1 + v2) / (1 + v1 * v2 / c**2)
assert sp.simplify(add.subs({v1: c * sp.tanh(a), v2: c * sp.tanh(bb)}) - c * sp.tanh(a + bb)) == 0
assert sp.simplify(add.subs(v2, c) - c) == 0
assert sp.factor(c**2 - add**2) == sp.factor((c**2 - v1**2) * (c**2 - v2**2) / (c**2 + v1 * v2)**2 * c**2)
print("8. dx1/dt = c tanh(eta), t = cosh(eta) tau; (v1 + v2)/(1 + v1 v2/c^2): c with c, < c for |v1|, |v2| < c")

# ---------------------------------------------------------------- 9.
# A pure electric field (B = 0) seen by the boosted observer: with the charge velocity
# u = -v e_1 relative to the observer, B' = u x E'/c^2 exactly.
F0 = svec(Ec)
Ep0, Bp0 = split(rev(R) * F0 * R)
uvel = [-c * sp.tanh(eta), 0, 0]
uxE = cross(uvel, Ep0)
for k in range(3):
    assert sp.simplify((Bp0[k] - uxE[k] / c**2).rewrite(sp.exp)) == 0
# low speed: Biot-Savart with I dl -> q u gives mu0 q u x r/(4 pi r^3) = u x E_Coulomb/c^2
eps0 = sp.Symbol("epsilon_0", positive=True)
mu0 = 1 / (eps0 * c**2)
rv = sp.symbols("r1:4", real=True)
ul = sp.symbols("w1:4", real=True)
rn = sp.sqrt(sum(v**2 for v in rv))
Ecoul = [q * v / (4 * sp.pi * eps0 * rn**3) for v in rv]
BS = [mu0 * q / (4 * sp.pi) * w / rn**3 for w in cross(ul, rv)]
for k in range(3):
    assert sp.simplify(BS[k] - cross(ul, Ecoul)[k] / c**2) == 0
print("9. pure E boosted: B' = u x E'/c^2 with u = -v; Biot-Savart for q u equals u x E/c^2")

# ---------------------------------------------------------------- 10.
# Current seen by the boosted observer: J = c rho gamma_0 (charges at rest) has components
# c rho' = J . gamma_0' = gamma c rho and J'_1 = -J . gamma_1' = -gamma rho v.
rho_ = sp.Symbol("rho", real=True)
Jst = c * rho_ * g0
g1p = R * g[1] * Rt
ipv = lambda p, r: ((p * r + r * p) / 2).scalar()
assert sp.simplify(ipv(Jst, g0p) - sp.cosh(eta) * c * rho_) == 0
assert sp.simplify((-ipv(Jst, g1p) + sp.cosh(eta) * rho_ * c * sp.tanh(eta)).rewrite(sp.exp)) == 0
print("10. static charge seen by the moving observer: rho' = gamma rho, J'_1 = -gamma rho v")

# ---------------------------------------------------------------- 11.
# gamma m c^2 = m c^2 + m u^2/2 + O(u^4); d(gamma m u)/dt -> m du/dt as u -> 0
m_, us = sp.symbols("m u", positive=True)
ser = sp.series(m_ * c**2 / sp.sqrt(1 - us**2 / c**2), us, 0, 4).removeO()
assert sp.expand(ser - (m_ * c**2 + m_ * us**2 / 2)) == 0
print("11. gamma m c^2 = m c^2 + m u^2/2 + ...")

# ---------------------------------------------------------------- 12.
# Intuitive statements: the rotated basis gamma'_mu = R gamma_mu R~ satisfies the same
# relations; light directions gamma_0 +- gamma_1 are only rescaled by e^{+-eta};
# the rotor equation dR/dtau = (q/2mc) F R gives m dU/dtau = (q/c) F . U for U = c R gamma_0 R~;
# uniform E (hyperbolic motion) and uniform B (cyclotron, omega = qB/m in proper time).
gp = [R * gm * Rt for gm in g]
for a_ in range(4):
    assert eq(gp[a_] * gp[a_], g[a_] * g[a_])
    for b_ in range(a_):
        assert eq(gp[a_] * gp[b_], -(gp[b_] * gp[a_]))
for sgn in (1, -1):
    assert eq(R * (g0 + sgn * g[1]) * Rt, sp.exp(sgn * eta) * (g0 + sgn * g[1]))
Rs = S.rnd(7, deg=0, grades={0, 2, 4})
Fs = S.rnd(8, deg=0, grades={2})
kk = sp.Symbol("k")
Rdot = kk * Fs * Rs
assert eq(Rdot * g0 * rev(Rs) + Rs * g0 * rev(Rdot), kk * (Fs * (Rs * g0 * rev(Rs)) - (Rs * g0 * rev(Rs)) * Fs))
m, E_, B_, tau, h = sp.symbols("m E_0 B_0 tau h", positive=True)
lor = lambda FF, UU: m * UU.map(lambda w: sp.diff(w, tau)) - q / c * ((FF * UU - UU * FF) / 2)
al = q * E_ / (m * c)
assert eq(lor(E_ * sg[0], c * (g0 * sp.cosh(al * tau) + g[1] * sp.sinh(al * tau))), 0)
om = q * B_ / m
UB = c * (g0 * sp.cosh(h) + sp.sinh(h) * (g[1] * sp.cos(om * tau) - g[2] * sp.sin(om * tau)))
assert eq(lor(I * c * B_ * sg[2], UB), 0)
assert eq((I * c * sg[0]) * g0 - g0 * (I * c * sg[0]), 0)  # magnetic part does not act on c gamma_0
print("12. gamma'_mu satisfy the same relations; R(g0 +- g1)R~ = e^{+-eta}(g0 +- g1); rotor equation; "
      "uniform E: rapidity qE tau/mc; uniform B: rotation qB tau/m (clockwise about B)")

# ---------------------------------------------------------------- 13.
import random
assert eq(sg[0] * sg[1] - sg[1] * sg[0], 2 * I * sg[2])
a_, b_ = sp.symbols("a_ b_", real=True)
Rab = (sp.cosh(a_ / 2) + sg[0] * sp.sinh(a_ / 2)) * (sp.cosh(b_ / 2) + sg[1] * sp.sinh(b_ / 2))
assert eq(Rab, sp.cosh(a_ / 2) * sp.cosh(b_ / 2) + sg[0] * sp.sinh(a_ / 2) * sp.cosh(b_ / 2)
          + sg[1] * sp.cosh(a_ / 2) * sp.sinh(b_ / 2) + I * sg[2] * sp.sinh(a_ / 2) * sp.sinh(b_ / 2))
rnd = random.Random(3)
for _ in range(3):
    Ev = [rnd.randint(-5, 5) for _ in range(3)]
    Bv = [rnd.randint(-5, 5) for _ in range(3)]
    Sv = cross(Ev, Bv)
    uv = (sum(x * x for x in Ev) + sum(x * x for x in Bv)) / 2  # units c = 1
    sn = sum(x * x for x in Sv) ** 0.5
    if sn == 0:
        continue
    et = 0.5 * float(sp.atanh(sn / uv))
    nS = sum((Sv[k] / sn * sg[k] for k in range(3)), S.zero())
    Rn = float(sp.cosh(et / 2)) + nS * float(sp.sinh(et / 2))
    F4 = sum((Ev[k] * sg[k] + I * Bv[k] * sg[k] for k in range(3)), S.zero())
    Fp_ = rev(Rn) * F4 * Rn
    Epn = [float(Fp_.d.get((1 << (k + 1)) | 1, 0)) for k in range(3)]
    Bpn = []
    for k in range(3):
        mask, sgn_ = list((I * sg[k]).d.items())[0]
        Bpn.append(float(Fp_.d.get(mask, 0)) * float(sgn_))
    assert max(abs(x) for x in cross(Epn, Bpn)) < 1e-9
print("13. [sigma_1, sigma_2] = 2 I sigma_3, boost composition has a rotation part; "
      "boost along S with tanh 2eta = |S|/cu gives E' parallel to B'")
