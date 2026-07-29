"""Checks for MEMO.md (geometric calculus in Cl(2,0) and complex analysis).

Cl(2,0) is realized as M_2(R) with e1 = diag(1,-1), e2 = [[0,1],[1,0]],
I = e1 e2 (I^2 = -1). Algebra elements are decomposed on the basis
(1, e1, e2, I) by the trace, and the Dirac operator is

    D = e1 d/dx + e2 d/dy.

1. D^2 = Laplacian on an arbitrary algebra-valued function.
2. e1 D F = (dx + I dy) F = 2 dbar F and D (e1 F) = (dx - I dy) F = 2 d F,
   so multiplying by e1 turns D into the Wirtinger operators, on whichever
   side the e1 is placed.
3. e1 acts as complex conjugation on the even part: e1 I e1^{-1} = -I.
4. Even-valued F = u + I v: D F = 0 is the Cauchy-Riemann equations
   u_x = v_y, u_y = -v_x (left monogenic = holomorphic), while
   F D = 0 gives the conjugate equations (right monogenic = antiholomorphic).
5. Odd-valued G = e1 H (H even): D G = 0 iff H is antiholomorphic, so the
   two grades carry opposite conditions. The grading obstruction is item 3.
6. Components of a monogenic function are harmonic (immediate from item 1).
7. The vector field x = x e1 + y e2 has x^{-1} = x / |x|^2, which is both
   left and right monogenic away from the origin, and e1 x^{-1} = 1/zbar.
8. For a vector-valued F, D F = D.F + D^F splits into the divergence
   (scalar) and the curl (bivector), so monogenic means both vanish.
"""

import sympy as sp

x, y = sp.symbols("x y", real=True)

one = sp.eye(2)
e1 = sp.Matrix([[1, 0], [0, -1]])
e2 = sp.Matrix([[0, 1], [1, 0]])
I2 = e1 * e2  # pseudoscalar, I2^2 = -1

assert e1**2 == one and e2**2 == one and I2**2 == -one
assert e1 * e2 == -e2 * e1

basis = {"1": one, "e1": e1, "e2": e2, "I": I2}


def comps(M):
    """Coefficients of M on (1, e1, e2, I)."""
    M = sp.simplify(M)
    return {k: sp.simplify(sp.trace(M * B.inv()) / 2) for k, B in basis.items()}


def elem(s=0, a=0, b=0, p=0):
    """s + a e1 + b e2 + p I."""
    return s * one + a * e1 + b * e2 + p * I2


def D(F):
    """Dirac (vector) derivative acting from the left."""
    return e1 * sp.diff(F, x) + e2 * sp.diff(F, y)


def Dr(F):
    """Dirac derivative acting from the right."""
    return sp.diff(F, x) * e1 + sp.diff(F, y) * e2


# --- 1. D^2 = Laplacian ------------------------------------------------
f = [sp.Function(n)(x, y) for n in ("f0", "f1", "f2", "f3")]
F = elem(*f)
lap = sp.diff(F, x, 2) + sp.diff(F, y, 2)
assert sp.simplify(D(D(F)) - lap) == sp.zeros(2, 2)

# --- 2. e1 D and D e1 are the Wirtinger operators ----------------------
assert sp.simplify(e1 * D(F) - (sp.diff(F, x) + I2 * sp.diff(F, y))) == sp.zeros(2, 2)
assert sp.simplify(D(e1 * F) - (sp.diff(F, x) - I2 * sp.diff(F, y))) == sp.zeros(2, 2)

# --- 3. e1 conjugates I ------------------------------------------------
assert e1 * I2 * e1.inv() == -I2

# --- 4. even part: left monogenic = Cauchy-Riemann ---------------------
u, v = sp.Function("u")(x, y), sp.Function("v")(x, y)
Fe = u * one + v * I2  # even, identified with u + i v
c = comps(D(Fe))
assert c["1"] == 0 and c["I"] == 0  # D maps even to odd
assert sp.simplify(c["e1"] - (sp.diff(u, x) - sp.diff(v, y))) == 0
assert sp.simplify(c["e2"] - (sp.diff(u, y) + sp.diff(v, x))) == 0

cr = comps(Dr(Fe))  # right monogenic gives the conjugate system
assert sp.simplify(cr["e1"] - (sp.diff(u, x) + sp.diff(v, y))) == 0
assert sp.simplify(cr["e2"] - (sp.diff(u, y) - sp.diff(v, x))) == 0

# holomorphic example: (x + I y)^3
z3 = (x * one + y * I2) ** 3
assert sp.simplify(D(z3)) == sp.zeros(2, 2)
assert sp.simplify(Dr(z3)) != sp.zeros(2, 2)

# --- 5. odd part carries the conjugate condition -----------------------
H = u * one + v * I2
G = e1 * H  # odd
cg = comps(D(G))
assert cg["e1"] == 0 and cg["e2"] == 0  # D maps odd to even
assert sp.simplify(cg["1"] - (sp.diff(u, x) + sp.diff(v, y))) == 0
assert sp.simplify(cg["I"] - (sp.diff(v, x) - sp.diff(u, y))) == 0
# i.e. u_x = -v_y, u_y = v_x: H is antiholomorphic
zbar3 = (x * one - y * I2) ** 3
assert sp.simplify(D(e1 * zbar3)) == sp.zeros(2, 2)

# --- 6. harmonic components -------------------------------------------
for k, e in comps(z3).items():
    assert sp.simplify(sp.diff(e, x, 2) + sp.diff(e, y, 2)) == 0

# --- 7. the vector x^{-1} ---------------------------------------------
vec = x * e1 + y * e2
r2 = x**2 + y**2
assert sp.simplify(vec * vec - r2 * one) == sp.zeros(2, 2)  # x^2 = |x|^2
inv = vec / r2
assert sp.simplify(vec * inv - one) == sp.zeros(2, 2)
assert sp.simplify(D(inv)) == sp.zeros(2, 2)
assert sp.simplify(Dr(inv)) == sp.zeros(2, 2)
zbar_inv = (x * one - y * I2).inv()  # 1/zbar as an even element
assert sp.simplify(e1 * inv - zbar_inv) == sp.zeros(2, 2)

# --- 8. divergence and curl -------------------------------------------
a, b = sp.Function("a")(x, y), sp.Function("b")(x, y)
V = a * e1 + b * e2
cv = comps(D(V))
assert sp.simplify(cv["1"] - (sp.diff(a, x) + sp.diff(b, y))) == 0  # div
assert sp.simplify(cv["I"] - (sp.diff(b, x) - sp.diff(a, y))) == 0  # curl

print("all checks passed")
