"""Checks for em/01-static.md (static electric and magnetic fields in Cl_{3,0}).

1. Integrand: with r = x - y, K = r/|r|^3 and the source S = rho - J/c, the
   grades of K S/eps_0 are  vector: rho K/eps_0 (Coulomb), bivector:
   -(K ^ J)/(c eps_0) = I c mu_0 (J x K) (Biot-Savart, mu_0 = 1/(eps_0 c^2)),
   scalar: -(K . J)/(c eps_0).
2. Scalar part: K . J = div_y(J/|r|) - (div J)/|r| (K = grad_y |r|^{-1}), so the
   scalar part of the integral vanishes for div J = 0 and decaying J.
3. Potential: D|x|^{-1} = -x/|x|^3 and K is monogenic away from 0. With
   P = phi - cA, F = -DP; the scalar part of F is c div A.
4. Example with Gaussian sources rho = rho0 e^{-r^2}, J = curl(j0 e^{-r^2} e_3):
   the Newton potential psi = (sqrt(pi)/4) erf(r)/r solves Lap psi = -e^{-r^2};
   F = -D(phi - cA) has no scalar part, DF = (rho - J/c)/eps0, and the four
   static equations hold. Numerical quadrature of the integral formula at an
   exterior point agrees with -DP.
5. Grade parts in DF for F = E + IcB: div E, -c curl B, I curl E, Ic div B.
6. Integral forms: grades of nF (volume theorem) and of dx F and
   (n x grad)F (surface theorem) give Gauss, div B, Ampere and curl E laws.
"""

import numpy as np
import sympy as sp

from common.clifford import MV, Alg, eq, mv

c, eps0 = sp.symbols("c epsilon_0", positive=True)
mu0 = 1 / (eps0 * c**2)
P3 = Alg(3)
e, I = P3.e, P3.I
X = P3.X


def vec(V):
    return sum((V[k] * e[k] for k in range(3)), P3.zero())


def comps(M):
    return [M.d.get(1 << k, 0) for k in range(3)]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def curl(V):
    return [sp.diff(V[2], X[1]) - sp.diff(V[1], X[2]),
            sp.diff(V[0], X[2]) - sp.diff(V[2], X[0]),
            sp.diff(V[1], X[0]) - sp.diff(V[0], X[1])]


def div(V):
    return sum(sp.diff(V[k], X[k]) for k in range(3))


def grad(f):
    return [sp.diff(f, X[k]) for k in range(3)]


# ---------------------------------------------------------------- 1.
rv = sp.symbols("r1:4", real=True)
Jv = sp.symbols("J1:4", real=True)
rho = sp.Symbol("rho", real=True)
rr = sp.sqrt(sum(v**2 for v in rv))
K = vec(rv) / rr**3
S = rho - vec(Jv) / c
KS = K * S
assert eq(KS.grade(1), rho * K)
wedge = (K * vec(Jv)).grade(2)
assert eq(KS.grade(2), -wedge / c)
assert eq(KS.grade(2) / eps0, I * c * mu0 * vec(cross(Jv, comps(K))))
assert eq(KS.grade(0), -sum(rv[k] * Jv[k] for k in range(3)) / rr**3 / c)
assert eq(KS.grade(3), 0)
# a ^ b = I (a x b)
av, bv = sp.symbols("a1:4", real=True), sp.symbols("b1:4", real=True)
assert eq((vec(av) * vec(bv)).grade(2), I * vec(cross(av, bv)))
print("1. K(rho - J/c)/eps0: vector rho K/eps0, bivector -(K^J)/(c eps0) = I c mu0 (J x K), scalar -(K.J)/(c eps0)")

# ---------------------------------------------------------------- 2.
Y = sp.symbols("y1:4", real=True)
Jf = [sp.Function(f"J{k}")(*Y) for k in (1, 2, 3)]
dist = sp.sqrt(sum((X[k] - Y[k])**2 for k in range(3)))
Kc = [(X[k] - Y[k]) / dist**3 for k in range(3)]
for k in range(3):
    assert sp.simplify(sp.diff(1 / dist, Y[k]) - Kc[k]) == 0
lhs = sum(Kc[k] * Jf[k] for k in range(3))
divy = lambda V: sum(sp.diff(V[k], Y[k]) for k in range(3))
rhs = divy([Jf[k] / dist for k in range(3)]) - divy(Jf) / dist
assert sp.simplify(lhs - rhs) == 0
print("2. K.J = div_y(J/|x-y|) - (div J)/|x-y|, with K = grad_y |x-y|^{-1}")

# ---------------------------------------------------------------- 3.
r = sp.sqrt(sum(v**2 for v in X))
assert eq(P3.D(mv(1 / r)), -P3.x / r**3)
assert eq(P3.D(P3.x / r**3), 0) and eq(P3.Dr(P3.x / r**3), 0)
phif = sp.Function("phi")(*X)
Af = [sp.Function(f"A{k}")(*X) for k in (1, 2, 3)]
Fp = -P3.D(phif - c * vec(Af))
assert eq(Fp.grade(0), c * div(Af))
assert eq(Fp.grade(1), vec([-v for v in grad(phif)]))
assert eq(Fp.grade(2), I * c * vec(curl(Af)))
assert eq(Fp.grade(3), 0)
print("3. D|x|^{-1} = -x/|x|^3, x/|x|^3 monogenic; -D(phi - cA) = c div A + (-grad phi) + Ic curl A")

# ---------------------------------------------------------------- 4.
rho0, j0 = sp.symbols("rho_0 j_0", positive=True)
g = sp.exp(-r**2)
psi = sp.sqrt(sp.pi) / 4 * sp.erf(r) / r
lap = lambda f: sum(sp.diff(f, v, 2) for v in X)
assert sp.simplify(lap(psi) + g) == 0
rhoG = rho0 * g
JG = curl([0, 0, j0 * g])
phiG = rho0 / eps0 * psi
AG = [mu0 * v for v in curl([0, 0, j0 * psi])]
PG = phiG - c * vec(AG)
FG = -P3.D(PG)
assert eq(FG.grade(0), 0) and eq(FG.grade(3), 0)
assert eq(P3.D(FG), (rhoG - vec(JG) / c) / eps0)
EG = comps(FG.grade(1))
BG = [v / c for v in comps((FG.grade(2) * I).map(lambda v: -v))]  # IcB -> cB: -I(IcB)
assert eq(FG.grade(2), I * c * vec(BG))
assert sp.simplify(div(EG) - rhoG / eps0) == 0
assert all(sp.simplify(v) == 0 for v in curl(EG))
assert all(sp.simplify(curl(BG)[k] - mu0 * JG[k]) == 0 for k in range(3))
assert sp.simplify(div(BG)) == 0
print("4. Gaussian sources: Lap psi = -e^{-r^2}, F = -D(phi - cA) has no scalar part, "
      "DF = (rho - J/c)/eps0, div E = rho/eps0, curl E = 0, curl B = mu0 J, div B = 0")

# numerical quadrature of (1/4 pi eps0) int (x - y)/|x - y|^3 (rho - J/c) dV at an
# exterior point (units rho0 = j0 = eps0 = c = 1)
subs = {rho0: 1, j0: 1, eps0: 1, c: 1}
x_pt = np.array([3.1, 1.7, -2.2])
n = 121
L = 4.0
ax = np.linspace(-L, L, n)
h = ax[1] - ax[0]
y1, y2, y3 = np.meshgrid(ax, ax, ax, indexing="ij")
gy = np.exp(-(y1**2 + y2**2 + y3**2))
rho_n = gy
J_n = [-2 * y2 * gy, 2 * y1 * gy, 0 * gy]  # curl(e^{-r^2} e_3)
d = [x_pt[0] - y1, x_pt[1] - y2, x_pt[2] - y3]
d3 = (d[0]**2 + d[1]**2 + d[2]**2) ** 1.5
w = h**3 / (4 * np.pi)
Evec = [np.sum(d[k] * rho_n / d3) * w for k in range(3)]
scal = -np.sum(sum(d[k] * J_n[k] for k in range(3)) / d3) * w
biv = [-np.sum((d[a] * J_n[b] - d[b] * J_n[a]) / d3) * w for a, b in ((1, 2), (2, 0), (0, 1))]
pt = dict(zip(X, x_pt))
Eex = [float(v.subs(subs).subs(pt)) for v in comps(FG.grade(1))]
Bex = [float(v.subs(subs).subs(pt)) for v in BG]
assert abs(scal) < 1e-8
assert np.allclose(Evec, Eex, atol=1e-8)
# bivector part -(1/c) sum (d_a J_b - d_b J_a) e_a e_b equals I c B, i.e. (e2e3, e3e1, e1e2) -> cB
assert np.allclose(biv, Bex, atol=1e-8)
print("   quadrature at", x_pt, ": scalar", f"{scal:.1e}", "E", np.round(Evec, 6), "cB", np.round(biv, 6))

# ---------------------------------------------------------------- 5.
Ef = [sp.Function(f"E{k}")(*X) for k in (1, 2, 3)]
Bf = [sp.Function(f"B{k}")(*X) for k in (1, 2, 3)]
F = vec(Ef) + I * c * vec(Bf)
G = P3.D(F)
assert eq(G.grade(0), div(Ef))
assert eq(G.grade(1), -c * vec(curl(Bf)))
assert eq(G.grade(2), I * vec(curl(Ef)))
assert eq(G.grade(3), I * c * div(Bf))
print("5. D(E + IcB) = div E - c curl B + I curl E + Ic div B")

# ---------------------------------------------------------------- 6.
nv = sp.symbols("n1:4", real=True)
Ec, Bc = sp.symbols("E1:4", real=True), sp.symbols("B1:4", real=True)
Fc = vec(Ec) + I * c * vec(Bc)
dot = lambda a, b: sum(a[k] * b[k] for k in range(3))
nF = vec(nv) * Fc
assert eq(nF.grade(0), dot(nv, Ec))
assert eq(nF.grade(1), -c * vec(cross(nv, Bc)))
assert eq(nF.grade(2), I * vec(cross(nv, Ec)))
assert eq(nF.grade(3), I * c * dot(nv, Bc))
xF = vec(av) * Fc
assert eq(xF.grade(0), dot(av, Ec)) and eq(xF.grade(3), I * c * dot(av, Bc))
# (D . dX) F with dX = I n: e_k . (I n) = n x e_k
for k in range(3):
    ekIn = (e[k] * (I * vec(nv)) - (I * vec(nv)) * e[k]) / 2
    assert eq(ekIn, vec(cross(nv, comps(e[k]))))
nxD = lambda H: sum((vec(cross(nv, comps(e[k]))) * P3.d(H, k) for k in range(3)), P3.zero())
Hs = nxD(F)
assert eq(Hs.grade(0), dot(nv, curl(Ef)))
assert eq(Hs.grade(3), I * c * dot(nv, curl(Bf)))
print("6. nF = n.E - c n x B + I n x E + Ic n.B; <dx F>_0 = dx.E, <dx F>_3 = Ic dx.B; "
      "(n x grad)F: scalar n.curl E, pseudoscalar Ic n.curl B")
