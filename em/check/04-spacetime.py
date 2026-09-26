"""Checks for em/04-spacetime.md (spacetime algebra Cl_{1,3} and Maxwell's equation).

1. Cl_{1,3}: gamma_0^2 = 1, gamma_k^2 = -1; D = sum gamma^mu d_mu with
   gamma^0 = gamma_0, gamma^k = -gamma_k; D x = 4, D^2 = d_0^2 - Laplacian.
2. Space-time split: sigma_k = gamma_k gamma_0 square to +1 and anticommute,
   I = sigma_1 sigma_2 sigma_3 = gamma_0 gamma_1 gamma_2 gamma_3; x gamma_0 = x_0 + x,
   gamma_0 x = x_0 - x, x^2 = x_0^2 - |x|^2; I anticommutes with vectors and
   commutes with even elements; gamma_0 D = Dq, D(gamma_0 H) = Dqbar H.
   The even subalgebra spanned by 1, sigma_k, I sigma_k, I is closed.
3. Bivectors: sigma_k contain gamma_0, I sigma_1 = gamma_3 gamma_2 etc. are spatial;
   sigma_k gamma_0 = -gamma_0 sigma_k, (I sigma_k) gamma_0 = gamma_0 (I sigma_k).
4. Current: J = c rho gamma_0 + sum J_k gamma_k, gamma_0 J = c(rho - J/c),
   D . J = d_t rho + div J.
5. For F = E + IcB (a bivector), D F has grades 1 and 3; gamma_0 (D . F) is the scalar
   and relative-vector part of Dq F (div E, d_0 E - c curl B), gamma_0 (D ^ F) is the
   rest (I(curl E + c d_0 B), Ic div B). With Maxwell's equations D . F = mu_0 c J,
   D ^ F = 0. D . (D . F) = 0 for any bivector F.
6. Light rays x_0(gamma_0 + n) are null; c gamma_0 + u_k gamma_k has square c^2 - |u|^2.
7. Observer split: E = (F - gamma_0 F gamma_0)/2 (planes containing gamma_0, anticommute with
   gamma_0), IcB = (F + gamma_0 F gamma_0)/2 (spatial planes, commute with gamma_0);
   vectors orthogonal to gamma_0 are the observer's space: x gamma_0 = x_0 + x.
8. I sigma_k is the orthogonal complement of the plane sigma_k: every vector of gamma_1 gamma_0
   is orthogonal to every vector of gamma_3 gamma_2, and F -> IF swaps the electric and magnetic planes;
   I(IF) = -F.
9. Boost x_0' = gamma(x_0 - beta x_1), x_1' = gamma(x_1 - beta x_0) preserves x_0^2 - x_1^2;
   simultaneous events (dx_0 = 0) get dx_0' = -gamma beta dx_1; for small beta, x_1' -> x_1 - v t.
"""

import sympy as sp

from common.clifford import MV, Alg, eq, mv

c, eps0 = sp.symbols("c epsilon_0", positive=True)
mu0 = 1 / (eps0 * c**2)
S = Alg(1, 3)
g = S.e
g0 = g[0]
X = S.X
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


def Dq(H):
    return S.d(H, 0) + sum((sg[k] * S.d(H, k + 1) for k in range(3)), S.zero())


def Dqbar(H):
    return S.d(H, 0) - sum((sg[k] * S.d(H, k + 1) for k in range(3)), S.zero())


# ---------------------------------------------------------------- 1.
assert eq(g0 * g0, 1) and all(eq(g[k] * g[k], -1) for k in (1, 2, 3))
assert eq(S.er[0], g0) and all(eq(S.er[k], -g[k]) for k in (1, 2, 3))
assert eq(S.D(S.x), 4)
f = sp.Function("f")(*X)
assert eq(S.D(S.D(mv(f))), sp.diff(f, X[0], 2) - sum(sp.diff(f, v, 2) for v in Xs))
print("1. Cl_{1,3}: gamma^0 = gamma_0, gamma^k = -gamma_k, D x = 4, D^2 = d0^2 - Laplacian")

# ---------------------------------------------------------------- 2.
for a in range(3):
    assert eq(sg[a] * sg[a], 1)
    for b in range(a):
        assert eq(sg[a] * sg[b], -(sg[b] * sg[a]))
assert eq(sg[0] * sg[1] * sg[2], I) and eq(g[0] * g[1] * g[2] * g[3], I)
xs = svec(Xs)
assert eq(S.x * g0, X[0] + xs) and eq(g0 * S.x, X[0] - xs)
assert eq(S.x * S.x, X[0]**2 - sum(v**2 for v in Xs))
for v in g:
    assert eq(I * v, -(v * I))
H = S.rnd(3, deg=1, grades={0, 2, 4})
assert eq(I * H, H * I)
Hg = S.rnd(4, deg=2)
assert eq(g0 * S.D(Hg), Dq(Hg))
assert eq(S.D(g0 * Hg), Dqbar(Hg))
even = [S.zero() + 1] + sg + [I * s for s in sg] + [I]
evmasks = {m for b in even for m in b.d}
assert len(evmasks) == 8 and all(grade % 2 == 0 for grade in (bin(m).count("1") for m in evmasks))
for a in even:
    for b in even:
        assert all(m in evmasks for m in (a * b).d)
print("2. sigma_k = gamma_k gamma_0: sigma^2 = 1, sigma1 sigma2 sigma3 = I = gamma_0123; "
      "x gamma_0 = x0 + x, x^2 = x0^2 - |x|^2; gamma_0 D = Dq, D(gamma_0 H) = Dqbar H; even part closed")

# ---------------------------------------------------------------- 3.
assert eq(I * sg[0], g[3] * g[2]) and eq(I * sg[1], g[1] * g[3]) and eq(I * sg[2], g[2] * g[1])
for s in sg:
    assert eq(s * g0, -(g0 * s)) and eq((I * s) * g0, g0 * (I * s))
print("3. I sigma_1 = gamma_3 gamma_2, I sigma_2 = gamma_1 gamma_3, I sigma_3 = gamma_2 gamma_1; "
      "sigma_k anticommute and I sigma_k commute with gamma_0")

# ---------------------------------------------------------------- 4.
rho = sp.Function("rho")(*X)
Jf = [sp.Function(f"J{k}")(*X) for k in (1, 2, 3)]
J = c * rho * g0 + sum((Jf[k] * g[k + 1] for k in range(3)), S.zero())
assert eq(g0 * J, c * (rho - svec(Jf) / c))
DJ = S.D(J)
# D . J = (1/c)... with x_0 = ct: d_0 (c rho) = d_t rho
assert eq(DJ.grade(0), c * sp.diff(rho, X[0]) + div(Jf))
print("4. J = c rho gamma_0 + J_k gamma_k: gamma_0 J = c(rho - J/c), D.J = c d0 rho + div J = d_t rho + div J")

# ---------------------------------------------------------------- 5.
Ef = [sp.Function(f"E{k}")(*X) for k in (1, 2, 3)]
Bf = [sp.Function(f"B{k}")(*X) for k in (1, 2, 3)]
F = svec(Ef) + I * c * svec(Bf)
assert all(eq(F.grade(r), 0) for r in (0, 1, 3, 4))
DF = S.D(F)
assert all(eq(DF.grade(r), 0) for r in (0, 2, 4))
dot, wedge = DF.grade(1), DF.grade(3)
cE, cB = curl(Ef), curl(Bf)
x0 = X[0]
assert eq(g0 * dot, div(Ef) + svec([sp.diff(Ef[k], x0) - c * cB[k] for k in range(3)]))
assert eq(g0 * wedge, I * svec([cE[k] + c * sp.diff(Bf[k], x0) for k in range(3)]) + I * c * div(Bf))
rhoM = div(Ef) * eps0
JM = [(cB[k] - sp.diff(Ef[k], x0) / c) * eps0 * c**2 for k in range(3)]
JM4 = c * rhoM * g0 + sum((JM[k] * g[k + 1] for k in range(3)), S.zero())
assert eq(dot, mu0 * c * JM4)
# D . (D . F) = 0 for a generic bivector
Bg = sum((sp.Function(f"F{m}")(*X) * MV({m: 1}, S.neg) for m in range(16) if bin(m).count("1") == 2), S.zero())
assert eq(S.D(S.D(Bg).grade(1)).grade(0), 0)
print("5. DF = D.F + D^F; gamma_0(D.F) = div E + (d0 E - c curl B), gamma_0(D^F) = I(curl E + c d0 B) + Ic div B; "
      "D.F = mu0 c J; D.(D.F) = 0")

# ---------------------------------------------------------------- 6.
# Classification: a light signal from the origin, x = x0 (gamma_0 + n) with |n| = 1, is null;
# a world line with speed |u| < c has timelike tangent c gamma_0 + u_k gamma_k.
th, ph, x0s = sp.symbols("theta phi x0", real=True)
nvec = [sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)]
light = x0s * (g0 + sum((nvec[k] * g[k + 1] for k in range(3)), S.zero()))
assert sp.simplify(sp.trigsimp((light * light).scalar())) == 0
uu = sp.symbols("u1:4", real=True)
tang = c * g0 + sum((uu[k] * g[k + 1] for k in range(3)), S.zero())
assert sp.expand((tang * tang).scalar() - (c**2 - sum(v**2 for v in uu))) == 0
print("6. light rays are null (x^2 = 0); the tangent c gamma_0 + u gamma has square c^2 - |u|^2 > 0 for |u| < c")

# ---------------------------------------------------------------- 7.
Ec_, Bc_ = sp.symbols("E1:4", real=True), sp.symbols("B1:4", real=True)
Fc_ = svec(Ec_) + I * c * svec(Bc_)
g0Fg0 = g0 * Fc_ * g0
assert eq((Fc_ - g0Fg0) / 2, svec(Ec_)) and eq((Fc_ + g0Fg0) / 2, I * c * svec(Bc_))
for s_ in sg:
    assert all(m & 1 for m in s_.d)          # sigma_k = gamma_k gamma_0 contains gamma_0
for s_ in sg:
    assert all(not (m & 1) for m in (I * s_).d)  # I sigma_k is a spatial plane
print("7. E = (F - g0 F g0)/2 on planes containing gamma_0, IcB = (F + g0 F g0)/2 on spatial planes")

# ---------------------------------------------------------------- 8.
ipv = lambda p_, q_: ((p_ * q_ + q_ * p_) / 2).scalar()
for (a1, a2), (b1, b2) in (((1, 0), (3, 2)), ((2, 0), (1, 3)), ((3, 0), (2, 1))):
    for u_ in (g[a1], g[a2]):
        for v_ in (g[b1], g[b2]):
            assert ipv(u_, v_) == 0
IF = I * Fc_
assert eq((IF - g0 * IF * g0) / 2, -c * svec(Bc_)) and eq((IF + g0 * IF * g0) / 2, I * svec(Ec_))
print("8. electric plane sigma_k and magnetic plane I sigma_k are orthogonal complements; IF swaps them")
assert eq(I * IF, -Fc_)

# ---------------------------------------------------------------- 9.
beta, t_, v_ = sp.symbols("beta t v", real=True)
x0_, x1_ = sp.symbols("x0 x1", real=True)
gam = 1 / sp.sqrt(1 - beta**2)
x0p, x1p = gam * (x0_ - beta * x1_), gam * (x1_ - beta * x0_)
assert sp.simplify(x0p**2 - x1p**2 - (x0_**2 - x1_**2)) == 0
dx1 = sp.Symbol("dx1", positive=True)
assert sp.simplify(gam * (0 - beta * dx1) + gam * beta * dx1) == 0
# first order in beta = v/c with x0 = ct: x1' ~ x1 - v t
lin = sp.series(x1p.subs({beta: v_ / c, x0_: c * t_}), v_, 0, 2).removeO()
assert sp.simplify(lin - (x1_ - v_ * t_)) == 0
print("9. boost preserves x0^2 - x1^2; simultaneous events: dx0' = -gamma beta dx1; x1' ~ x1 - v t for small v")
