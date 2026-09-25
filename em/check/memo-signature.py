"""Checks for the planning memo (em/MEMO.md): signature choice and key formulas.

1. Cl_{3,0}: with the paravector operator Dq = d_0 + sum e_k d_k (x_0 = ct) and
   F = E + I c B, the grades of Dq F are
   scalar: div E, vector: d_0 E - c curl B, bivector: I (curl E + c d_0 B),
   pseudoscalar: I c div B. With Maxwell's equations in SI units
   (d_0 = d_t / c, mu_0 = 1/(eps_0 c^2)) this is Dq F = (rho - J/c)/eps_0.
   Dqbar Dq = d_0^2 - Laplacian (Dqbar = d_0 - sum e_k d_k), and the scalar part
   of Dqbar applied to the source gives the continuity equation.
2. Cl_{3,0}: (eps_0/2) F F^dagger = (eps_0/2)(|E|^2 + c^2 |B|^2) + eps_0 c E x B.
3. Signature: paravector norm p pbar = p_0^2 - |p|^2; the block matrices
   gamma(p) = [[0, p], [pbar, 0]] with sigma^mu = (1, sigma_k),
   sigmabar^mu = (1, -sigma_k) satisfy {gamma^mu, gamma^nu} = 2 eta^{mu nu},
   eta = diag(+, -, -, -) (Weyl representation, Cl_{1,3}); with the lower block
   -pbar the signature is (-, +, +, +) (Cl_{3,1}).
4. Cl_{1,3} (gamma_0^2 = 1, gamma_k^2 = -1): sigma_k = gamma_k gamma_0 square to
   +1 and anticommute, sigma_1 sigma_2 sigma_3 = gamma_0 gamma_1 gamma_2 gamma_3,
   gamma_0 D = d_0 + sum sigma_k d_k and D gamma_0 = d_0 - sum sigma_k d_k on
   generic functions, D^2 = d_0^2 - Laplacian. In Cl_{3,1} (time generator
   squaring -1) D^2 = Laplacian - d_0^2. In Cl_{4,0}, h_l = e_0 e_l square to -1.
5. Cl_{1,3}: for the bivector F = sum E_k sigma_k + I c sum B_k sigma_k,
   gamma_0 D F reproduces item 1 under sigma_k -> e_k, and D F has only grades
   1 and 3; with Maxwell's equations D F = mu_0 c J, J = c rho gamma_0 + sum J_k gamma_k.
"""

import sympy as sp

from common.clifford import MV, Alg, eq, mv

c, eps0 = sp.symbols("c epsilon_0", positive=True)
x0, x1, x2, x3 = X = sp.symbols("x0:4", real=True)
Xs = X[1:]


def fields(name):
    return [sp.Function(f"{name}{k}")(*X) for k in (1, 2, 3)]


def curl(V):
    return [sp.diff(V[2], x2) - sp.diff(V[1], x3),
            sp.diff(V[0], x3) - sp.diff(V[2], x1),
            sp.diff(V[1], x1) - sp.diff(V[0], x2)]


def div(V):
    return sum(sp.diff(V[k], Xs[k]) for k in range(3))


Ef, Bf = fields("E"), fields("B")

# ---------------------------------------------------------------- 1.
P = Alg(3)
e, I = P.e, P.I


def vec(V):
    return sum((V[k] * e[k] for k in range(3)), MV())


def d(F, x):
    return F.map(lambda v: sp.diff(v, x))


def Dq(F):
    return d(F, x0) + sum((e[k] * d(F, Xs[k]) for k in range(3)), MV())


def Dqbar(F):
    return d(F, x0) - sum((e[k] * d(F, Xs[k]) for k in range(3)), MV())


F = vec(Ef) + I * c * vec(Bf)
G = Dq(F)
cE, cB = curl(Ef), curl(Bf)
assert eq(G.grade(0), div(Ef))
assert eq(G.grade(1), vec([sp.diff(Ef[k], x0) - c * cB[k] for k in range(3)]))
assert eq(G.grade(2), I * vec([cE[k] + c * sp.diff(Bf[k], x0) for k in range(3)]))
assert eq(G.grade(3), I * c * div(Bf))
# Maxwell's equations in SI: div E = rho/eps0, curl B = mu0 J + (1/c^2) d_t E
# with d_t = c d_0, mu0 = 1/(eps0 c^2): the vector part becomes -J/(c eps0).
rho = div(Ef) * eps0
J = [(cB[k] - sp.diff(Ef[k], x0) / c) * eps0 * c**2 for k in range(3)]
# the bivector and pseudoscalar parts vanish by the homogeneous equations
assert eq(G.grade(0) + G.grade(1), (rho - vec(J) / c) / eps0)
f = sp.Function("f")(*X)
assert eq(Dqbar(Dq(mv(f))), sp.diff(f, x0, 2) - sum(sp.diff(f, v, 2) for v in Xs))
# continuity: Dqbar Dq F has no scalar part since F has none, and the scalar part of
# Dqbar((rho - J/c)/eps0) is (d_0 rho + div J / c)/eps0 = (d_t rho + div J)/(c eps0)
src = (rho - vec(J) / c) / eps0
assert eq(Dqbar(src).grade(0), (sp.diff(rho, x0) + div(J) / c) / eps0)
assert eq(Dqbar(src).grade(0), 0)
print("1. Cl_{3,0}: Dq(E + IcB) = div E + (d0 E - c curl B) + I(curl E + c d0 B) + Ic div B "
      "= (rho - J/c)/eps0; Dqbar Dq = d0^2 - Laplacian; continuity")

# ---------------------------------------------------------------- 2.
E0 = sp.symbols("E1:4", real=True)
B0 = sp.symbols("B1:4", real=True)
Fc = vec(E0) + I * c * vec(B0)
T = eps0 / 2 * Fc * Fc.rev()
cross = [E0[1] * B0[2] - E0[2] * B0[1], E0[2] * B0[0] - E0[0] * B0[2], E0[0] * B0[1] - E0[1] * B0[0]]
assert eq(T, eps0 / 2 * (sum(v**2 for v in E0) + c**2 * sum(v**2 for v in B0)) + eps0 * c * vec(cross))
print("2. (eps0/2) F F^dagger = (eps0/2)(|E|^2 + c^2|B|^2) + eps0 c E x B")

# ---------------------------------------------------------------- 3.
p = sp.symbols("p0:4", real=True)
pv = p[0] + vec(p[1:])
pb = p[0] - vec(p[1:])
assert eq(pv * pb, p[0]**2 - p[1]**2 - p[2]**2 - p[3]**2)
s0 = sp.eye(2)
sk = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
sig = [s0] + sk
sigb = [s0] + [-m for m in sk]
Z = sp.zeros(2)
for lower, eta in ((1, sp.diag(1, -1, -1, -1)), (-1, sp.diag(-1, 1, 1, 1))):
    g = [sp.Matrix(sp.BlockMatrix([[Z, sig[m]], [lower * sigb[m], Z]])) for m in range(4)]
    for m in range(4):
        for n in range(4):
            assert g[m] * g[n] + g[n] * g[m] == 2 * eta[m, n] * sp.eye(4)
print("3. p pbar = p0^2 - |p|^2; [[0, sigma^mu], [sigmabar^mu, 0]] gives (+,-,-,-) (Cl_{1,3}), "
      "lower block -sigmabar gives (-,+,+,+) (Cl_{3,1})")

# ---------------------------------------------------------------- 4.
S = Alg(1, 3)
g0, g1, g2, g3 = S.e
sg = [S.e[k] * g0 for k in (1, 2, 3)]
for a in range(3):
    assert eq(sg[a] * sg[a], 1)
    for b in range(a):
        assert eq(sg[a] * sg[b], -(sg[b] * sg[a]))
assert eq(sg[0] * sg[1] * sg[2], S.I)
H = S.rnd(1, deg=2)
assert eq(g0 * S.D(H), S.d(H, 0) + sum((sg[k] * S.d(H, k + 1) for k in range(3)), S.zero()))
# D gamma_0 as an operator acting from the left: (sum gamma^mu gamma_0 d_mu) H
Dg0 = sum((S.er[a] * g0 * S.d(H, a) for a in range(4)), S.zero())
assert eq(Dg0, S.d(H, 0) - sum((sg[k] * S.d(H, k + 1) for k in range(3)), S.zero()))
assert eq(S.D(S.D(H)), S.lap(H))
assert eq(S.lap(mv(S.X[0] ** 2 + S.X[1] ** 2)), 0)  # d_0^2 - d_1^2: hyperbolic
M = Alg(3, 1)  # e_3 squares to -1 (time)
assert eq(M.D(M.D(mv(M.X[3] ** 2))), -2) and eq(M.D(M.D(mv(M.X[0] ** 2))), 2)
C4 = Alg(4)
for l in (1, 2, 3):
    assert eq((C4.e[0] * C4.e[l]) * (C4.e[0] * C4.e[l]), -1)
print("4. Cl_{1,3}: sigma_k = gamma_k gamma_0 square to +1, sigma1 sigma2 sigma3 = I, "
      "gamma_0 D = d0 + sigma.d, D gamma_0 = d0 - sigma.d, D^2 = d0^2 - Laplacian; "
      "Cl_{3,1}: D^2 = Laplacian - d0^2; Cl_{4,0}: (e0 e_l)^2 = -1")

# ---------------------------------------------------------------- 5.
EfS = [sp.Function(f"E{k}")(*S.X) for k in (1, 2, 3)]
BfS = [sp.Function(f"B{k}")(*S.X) for k in (1, 2, 3)]
vS = lambda V: sum((V[k] * sg[k] for k in range(3)), S.zero())
FS = vS(EfS) + S.I * c * vS(BfS)
assert eq(FS.grade(0), 0) and eq(FS.grade(1), 0) and eq(FS.grade(3), 0) and eq(FS.grade(4), 0)
GS = g0 * S.D(FS)
Y = S.X[1:]
curlS = lambda V: [sp.diff(V[2], Y[1]) - sp.diff(V[1], Y[2]),
                   sp.diff(V[0], Y[2]) - sp.diff(V[2], Y[0]),
                   sp.diff(V[1], Y[0]) - sp.diff(V[0], Y[1])]
divS = lambda V: sum(sp.diff(V[k], Y[k]) for k in range(3))
cES, cBS = curlS(EfS), curlS(BfS)
x0S = S.X[0]
assert eq(GS, divS(EfS) + vS([sp.diff(EfS[k], x0S) - c * cBS[k] for k in range(3)])
          + S.I * vS([cES[k] + c * sp.diff(BfS[k], x0S) for k in range(3)]) + S.I * c * divS(BfS))
DF = S.D(FS)
assert eq(DF.grade(0), 0) and eq(DF.grade(2), 0) and eq(DF.grade(4), 0)
mu0 = 1 / (eps0 * c**2)
rhoS = divS(EfS) * eps0
JS = [(cBS[k] - sp.diff(EfS[k], x0S) / c) * eps0 * c**2 for k in range(3)]
J4 = c * rhoS * g0 + sum((JS[k] * S.e[k + 1] for k in range(3)), S.zero())
assert eq(DF.grade(1), mu0 * c * J4)
print("5. Cl_{1,3}: F = E + IcB is a bivector, gamma_0 D F reproduces 1., "
      "D F has grades 1 and 3, vector part = mu0 c J with J = c rho gamma_0 + J_k gamma_k")
