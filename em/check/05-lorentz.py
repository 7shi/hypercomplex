"""Checks for em/05-lorentz.md (Lorentz transformations as rotors in Cl_{1,3}).

1. Exponentials: sigma_1^2 = +1 gives e^{sigma_1 a} = cosh a + sigma_1 sinh a,
   (I sigma_3)^2 = -1 gives e^{I sigma_3 a} = cos a + I sigma_3 sin a (power series
   truncated and compared via the recurrence of the square).
2. Boost: R = e^{sigma_1 phi/2}, R R~ = 1, R gamma_0 R~ = cosh(phi) gamma_0 + sinh(phi) gamma_1,
   R gamma_1 R~ = sinh(phi) gamma_0 + cosh(phi) gamma_1, gamma_2, gamma_3 fixed; x^2 invariant.
   Rapidities add: e^{sigma_1 a/2} e^{sigma_1 b/2} = e^{sigma_1 (a+b)/2}.
3. Rotation: R = e^{-I sigma_3 theta/2} rotates gamma_1 -> cos gamma_1 + sin gamma_2,
   fixes gamma_0 and gamma_3; it commutes with gamma_0 and acts on sigma_k as a rotation.
4. Fields seen by an observer with gamma_0' = R gamma_0 R~ (velocity v e_1,
   v/c = tanh phi): F' = R~ F R has E'_1 = E_1, E'_2 = gamma(E_2 - v B_3),
   E'_3 = gamma(E_3 + v B_2), B'_1 = B_1, B'_2 = gamma(B_2 + v E_3/c^2),
   B'_3 = gamma(B_3 - v E_2/c^2). Equivalently F = sum E'_k sigma'_k + Ic sum B'_k sigma'_k
   with sigma'_k = R sigma_k R~.
5. Invariants: F^2 = (|E|^2 - c^2|B|^2) + 2Ic E.B commutes with R, so F'^2 = F^2.
6. Lorentz force: v = gamma(c gamma_0 + u_k gamma_k), v^2 = c^2; (q/c) F . v
   (F . v = (Fv - vF)/2) has gamma_0 component gamma q E.u / c and gamma_k components
   gamma q (E + u x B)_k; v . (F . v) = 0.
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
phi = sp.Symbol("phi", real=True)
R = sp.cosh(phi / 2) + sg[0] * sp.sinh(phi / 2)
Rt = rev(R)
assert eq(Rt, sp.cosh(phi / 2) - sg[0] * sp.sinh(phi / 2))
assert eq(R * Rt, 1)
assert eq(R * g0 * Rt, sp.cosh(phi) * g0 + sp.sinh(phi) * g[1])
assert eq(R * g[1] * Rt, sp.sinh(phi) * g0 + sp.cosh(phi) * g[1])
assert eq(R * g[2] * Rt, g[2]) and eq(R * g[3] * Rt, g[3])
xp = R * S.x * Rt
assert eq(xp * xp, S.x * S.x)
bb = sp.Symbol("b", real=True)
Ra = sp.cosh(a / 2) + sg[0] * sp.sinh(a / 2)
Rb = sp.cosh(bb / 2) + sg[0] * sp.sinh(bb / 2)
assert eq(Ra * Rb, sp.cosh((a + bb) / 2) + sg[0] * sp.sinh((a + bb) / 2))
print("2. R = e^{sigma_1 phi/2}: gamma_0 -> cosh gamma_0 + sinh gamma_1, gamma_1 -> sinh gamma_0 + cosh gamma_1, "
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
gam = sp.cosh(phi)
vv = c * sp.tanh(phi)
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
