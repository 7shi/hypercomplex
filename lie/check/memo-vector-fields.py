"""Checks for lie/MEMO.md (球面上のベクトル場と角運動量).

Verifies symbolically (sympy):

1. In spherical coordinates the generators of the rotations are the
   Killing fields X1 = -sinφ ∂θ - cotθ cosφ ∂φ,
   X2 = cosφ ∂θ - cotθ sinφ ∂φ, X3 = ∂φ, reproducing dn/dt = e_a × n,
   with [X_a, X_b] = -ε_abc X_c (anti-homomorphism for a left action).
2. L_a = -i L_{X_a} satisfies [L_a, L_b] = i ε_abc L_c and
   L3 = -i∂φ, L± = e^{±iφ}(±∂θ + i cotθ ∂φ).
3. L² = L3² + (L+L- + L-L+)/2 equals -Δ_{S²} and gives
   L² Y_l^m = l(l+1) Y_l^m (l ≤ 2).

The complex coordinate w = β/α and the Möbius side are checked in
hopf/check/memo-mobius.py.
"""

import sympy as sp

I = sp.I
th, ph = sp.symbols("theta phi", real=True)

n = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])

# sample points for numerical comparison (symbolic simplify is unreliable here)
POINTS = [
    {th: sp.Rational(37, 50), ph: sp.Rational(11, 10)},
    {th: sp.Rational(9, 5), ph: sp.Rational(-3, 4)},
    {th: sp.Rational(1, 3), ph: sp.Rational(5, 2)},
]

def num_zero(expr, subs, tol=1e-20):
    return abs(complex(sp.N(expr.subs(subs)))) < tol

def is_zero(expr):
    """Zero test for expressions in (θ, φ).

    Half-angle identities defeat simplify(), so fall back to numerical
    evaluation at sample points.
    """
    if sp.simplify(expr) == 0:
        return True
    return all(num_zero(sp.expand(expr), pt, tol=1e-12) for pt in POINTS)

KILLING = {
    1: (-sp.sin(ph), -sp.cos(ph) / sp.tan(th)),
    2: (sp.cos(ph), -sp.sin(ph) / sp.tan(th)),
    3: (sp.Integer(0), sp.Integer(1)),
}

def apply_field(V, f):
    dth, dph = V
    return dth * sp.diff(f, th) + dph * sp.diff(f, ph)

def bracket(V1, V2):
    """[V1, V2] as a vector field, in components."""
    return tuple(
        sp.simplify(apply_field(V1, c2) - apply_field(V2, c1))
        for c1, c2 in zip(V1, V2)
    )

def check_killing_fields():
    axes = {1: (1, 0, 0), 2: (0, 1, 0), 3: (0, 0, 1)}
    e_th = sp.Matrix([sp.cos(th) * sp.cos(ph), sp.cos(th) * sp.sin(ph), -sp.sin(th)])
    e_ph = sp.Matrix([-sp.sin(ph), sp.cos(ph), 0])
    for a in (1, 2, 3):
        dn = sp.Matrix(axes[a]).cross(n)
        dth, dph = KILLING[a]
        assert is_zero(dn.dot(e_th) - dth)
        assert is_zero(dn.dot(e_ph) / sp.sin(th) - dph)
    # [X_a, X_b] = -ε_abc X_c
    for a, b, c in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
        got = bracket(KILLING[a], KILLING[b])
        for x, y in zip(got, KILLING[c]):
            assert is_zero(x + y)

# --- angular momentum --------------------------------------------------------

def L(a, f):
    return sp.simplify(-I * apply_field(KILLING[a], f))

def Lp(f):
    return sp.simplify(sp.exp(I * ph) * (sp.diff(f, th) + I * sp.cot(th) * sp.diff(f, ph)))

def Lm(f):
    return sp.simplify(sp.exp(-I * ph) * (-sp.diff(f, th) + I * sp.cot(th) * sp.diff(f, ph)))

def L2(f):
    return sp.simplify(L(3, L(3, f)) + (Lp(Lm(f)) + Lm(Lp(f))) / 2)

def laplacian(f):
    return sp.simplify(
        sp.diff(sp.sin(th) * sp.diff(f, th), th) / sp.sin(th)
        + sp.diff(f, ph, 2) / sp.sin(th) ** 2
    )

TEST = [
    sp.cos(th),
    sp.sin(th) * sp.exp(I * ph),
    sp.sin(th) ** 2 * sp.exp(2 * I * ph),
    sp.sin(th) * sp.cos(th) * sp.exp(-I * ph),
]

def check_angular_momentum():
    # [L_a, L_b] = i ε_abc L_c
    for a, b, c in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
        for f in TEST:
            assert is_zero(L(a, L(b, f)) - L(b, L(a, f)) - I * L(c, f))
    # explicit forms
    for f in TEST:
        assert is_zero(L(3, f) + I * sp.diff(f, ph))
        assert is_zero(L(1, f) + I * L(2, f) - Lp(f))
        assert is_zero(L(1, f) - I * L(2, f) - Lm(f))

def check_casimir():
    for f in TEST:
        assert is_zero(L2(f) + laplacian(f))
    # L² Y_l^m = l(l+1) Y_l^m
    for l in range(3):
        for m in range(-l, l + 1):
            Y = sp.simplify(sp.Ynm(l, m, th, ph).expand(func=True))
            assert is_zero(L2(Y) - l * (l + 1) * Y)
            assert is_zero(L(3, Y) - m * Y)

if __name__ == "__main__":
    check_killing_fields()
    check_angular_momentum()
    check_casimir()
    print("check-memo-vector-fields: all checks passed")
