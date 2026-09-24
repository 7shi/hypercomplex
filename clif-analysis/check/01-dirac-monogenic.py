"""Checks for 01-dirac-monogenic.md.

Geometric calculus in Cl_{2,0}(R) and complex analysis. Cl_{2,0}(R) is
realized as M_2(R) with e1 = diag(1,-1), e2 = [[0,1],[1,0]], I = e1 e2
(I^2 = -1). Algebra elements are decomposed on the basis (1, e1, e2, I) by
the trace, and the Dirac operator is

    D = e1 d/dx + e2 d/dy.

1. D^2 = Laplacian on an arbitrary algebra-valued function; with the
   signature Cl_{0,2}(R) the square is minus the Laplacian.
2. e1 x = z and x e1 = zbar for the vector x = x e1 + y e2 and z = x + I y;
   e1 D = dx + I dy = 2 dbar, so multiplying by e1 turns D into the
   Wirtinger operator.
3. Even-valued F = u + I v: D F = 0 is the Cauchy-Riemann equations
   u_x = v_y, u_y = -v_x (left monogenic = holomorphic), while F D = 0 gives
   the conjugate equations (right monogenic = antiholomorphic).
   Example: z^3 is left but not right monogenic.
4. e1 H = Hbar e1 for even H (e1 I e1^{-1} = -I). An odd function written as
   H e1 is monogenic iff H is holomorphic (D(H e1) = (D H) e1), and written
   as e1 H iff H is antiholomorphic: the grade alone does not decide it.
5. Components of a monogenic function are harmonic.
6. x^2 = |x|^2, x^{-1} = x / |x|^2 = (1/z) e1, and x^{-1} is both left and
   right monogenic away from the origin.
7. For a vector-valued F, D F = D.F + D^F splits into the divergence
   (scalar) and the curl (bivector); F D flips the sign of the curl, so
   for odd-valued functions left and right monogenicity are equivalent.
8. A scalar-valued monogenic function is constant. With the idempotent
   P = (1 + e1)/2, the minimal left ideal Cl P is 2-dimensional, F -> F P is
   injective on the even part (and on the odd part), and D(F P) = (D F) P,
   so even-valued monogenic functions correspond one-to-one to
   ideal-valued ones.
"""

import sympy as sp

x, y = sp.symbols("x y", real=True)

one = sp.eye(2)
zero = sp.zeros(2, 2)
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
    """Dirac operator acting from the left."""
    return e1 * sp.diff(F, x) + e2 * sp.diff(F, y)


def Dr(F):
    """Dirac operator acting from the right."""
    return sp.diff(F, x) * e1 + sp.diff(F, y) * e2


def is_zero(M):
    return sp.simplify(M) == zero


# --- 1. D^2 = Laplacian ------------------------------------------------
f = [sp.Function(n)(x, y) for n in ("f0", "f1", "f2", "f3")]
F = elem(*f)
lap = sp.diff(F, x, 2) + sp.diff(F, y, 2)
assert is_zero(D(D(F)) - lap)

# Cl_{0,2}(R) = H: generators i, j with i^2 = j^2 = -1 (as 2x2 complex)
qi = sp.Matrix([[sp.I, 0], [0, -sp.I]])
qj = sp.Matrix([[0, 1], [-1, 0]])
assert qi**2 == -one and qj**2 == -one and qi * qj == -qj * qi
g = sp.Function("g")(x, y)
Dq = lambda h: qi * sp.diff(h, x) + qj * sp.diff(h, y)
assert is_zero(Dq(Dq(g * one)) + (sp.diff(g, x, 2) + sp.diff(g, y, 2)) * one)

# --- 2. e1 x = z, x e1 = zbar, e1 D = 2 dbar ---------------------------
vec = x * e1 + y * e2
z = x * one + y * I2
zbar = x * one - y * I2
assert is_zero(e1 * vec - z)
assert is_zero(vec * e1 - zbar)
assert is_zero(e1 * D(F) - (sp.diff(F, x) + I2 * sp.diff(F, y)))

# --- 3. even part: left monogenic = Cauchy-Riemann ---------------------
u, v = sp.Function("u")(x, y), sp.Function("v")(x, y)
H = u * one + v * I2  # even, identified with u + i v
c = comps(D(H))
assert c["1"] == 0 and c["I"] == 0  # D maps even to odd
assert sp.simplify(c["e1"] - (sp.diff(u, x) - sp.diff(v, y))) == 0
assert sp.simplify(c["e2"] - (sp.diff(u, y) + sp.diff(v, x))) == 0

cr = comps(Dr(H))  # right monogenic gives the conjugate system
assert sp.simplify(cr["e1"] - (sp.diff(u, x) + sp.diff(v, y))) == 0
assert sp.simplify(cr["e2"] - (sp.diff(u, y) - sp.diff(v, x))) == 0

z3 = z**3
assert is_zero(D(z3))
assert not is_zero(Dr(z3))
assert is_zero(Dr(zbar**3))

# --- 4. e1 conjugates; odd part depends on the side of e1 --------------
assert e1 * I2 * e1.inv() == -I2
Hbar = u * one - v * I2
assert is_zero(e1 * H - Hbar * e1)
assert is_zero(D(H * e1) - D(H) * e1)  # H e1 monogenic iff H holomorphic
cg = comps(D(e1 * H))  # e1 H monogenic iff H antiholomorphic
assert sp.simplify(cg["1"] - (sp.diff(u, x) + sp.diff(v, y))) == 0
assert sp.simplify(cg["I"] - (sp.diff(v, x) - sp.diff(u, y))) == 0
assert is_zero(D(z3 * e1)) and not is_zero(D(e1 * z3))
assert is_zero(D(e1 * zbar**3))

# --- 5. harmonic components -------------------------------------------
for k, e in comps(z3).items():
    assert sp.simplify(sp.diff(e, x, 2) + sp.diff(e, y, 2)) == 0

# --- 6. the vector x^{-1} ---------------------------------------------
r2 = x**2 + y**2
assert is_zero(vec * vec - r2 * one)  # x^2 = |x|^2
inv = vec / r2
assert is_zero(vec * inv - one)
assert is_zero(inv - z.inv() * e1)  # x^{-1} = (1/z) e1
assert is_zero(D(inv)) and is_zero(Dr(inv))

# --- 7. divergence and curl -------------------------------------------
a, b = sp.Function("a")(x, y), sp.Function("b")(x, y)
V = a * e1 + b * e2
cv = comps(D(V))
assert sp.simplify(cv["1"] - (sp.diff(a, x) + sp.diff(b, y))) == 0  # div
assert sp.simplify(cv["I"] - (sp.diff(b, x) - sp.diff(a, y))) == 0  # curl
assert cv["e1"] == 0 and cv["e2"] == 0
cvr = comps(Dr(V))  # right action flips only the bivector part
assert sp.simplify(cvr["1"] - cv["1"]) == 0
assert sp.simplify(cvr["I"] + cv["I"]) == 0  # so D V = 0 iff V D = 0

# --- 8. scalar values do not suffice; minimal left ideal --------------
s = sp.Function("s")(x, y)
cs = comps(D(s * one))  # D s = e1 s_x + e2 s_y
assert cs["e1"] == sp.diff(s, x) and cs["e2"] == sp.diff(s, y)

P = (one + e1) / 2
assert P * P == P
# Cl P is spanned by (1+e1)/2 and e2 (1+e1)/2: first column of M_2(R)
ideal = [B * P for B in basis.values()]
assert sp.Matrix.hstack(*[M.reshape(4, 1) for M in ideal]).rank() == 2
# F -> F P is injective on the even part and on the odd part
p0, p1 = sp.symbols("p0 p1", real=True)
assert sp.solve(list(((p0 * one + p1 * I2) * P)), [p0, p1]) == {p0: 0, p1: 0}
assert sp.solve(list(((p0 * e1 + p1 * e2) * P)), [p0, p1]) == {p0: 0, p1: 0}
assert is_zero(D(H * P) - D(H) * P)

print("all checks passed")
