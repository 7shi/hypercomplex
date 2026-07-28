"""Checks for hopf/MEMO.md (比を取ると行列が微分になる部分).

Verifies symbolically (sympy):

1. w = β/α = tan(θ/2)e^{iφ} = (n1+in2)/(1+n3) is the stereographic
   projection of the Bloch vector from the south pole (w = 0 at the north
   pole, w = ∞ at the south pole).
2. A generator X = [[a,b],[c,d]] acting linearly on (α,β) moves the ratio
   w by dw/dt = c + (d-a)w - b w², checked against direct exponentiation
   of exp(tX) and the Möbius form (c+dw)/(a+bw).
3. The generators -iσ_a/2 give the coefficients
   X1 = i(w²-1)/2, X2 = (w²+1)/2, X3 = i w, and these reproduce the
   rotation dn/dt = e_a × n of the Bloch vector.
4. Ladder operators map to σ+ ↦ -w²∂w, σ- ↦ ∂w, σ3 ↦ -2w ∂w, and σ±
   move w in a way no rotation of the Bloch sphere reproduces.
"""

import sympy as sp

I = sp.I
w = sp.Symbol("w")
t = sp.Symbol("t", real=True)
th, ph = sp.symbols("theta phi", real=True)

s1 = sp.Matrix([[0, 1], [1, 0]])
s2 = sp.Matrix([[0, -I], [I, 0]])
s3 = sp.Matrix([[1, 0], [0, -1]])
sp_ = (s1 + I * s2) / 2
sm_ = (s1 - I * s2) / 2

# Bloch vector and the stereographic coordinate
n = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
W = sp.tan(th / 2) * sp.exp(I * ph)

# coefficients of ∂w for the generators -i σ_a / 2
HOL = {1: I * (w**2 - 1) / 2, 2: (w**2 + 1) / 2, 3: I * w}

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

def dw_dt(X):
    """Motion of w = β/α induced by the generator X."""
    a, b, c, d = X[0, 0], X[0, 1], X[1, 0], X[1, 1]
    return sp.expand(c + (d - a) * w - b * w**2)

def field_from_bloch(axis):
    """dw/dt obtained from dn/dt = e_axis × n through w = (n1+in2)/(1+n3)."""
    a, b, c = sp.symbols("a b c")
    dn = sp.Matrix(axis).cross(n)
    Wabc = (a + I * b) / (1 + c)
    d = sum(sp.diff(Wabc, s) * v for s, v in zip((a, b, c), dn))
    return d.subs({a: n[0], b: n[1], c: n[2]})

def check_stereographic():
    assert is_zero(W - (n[0] + I * n[1]) / (1 + n[2]))
    # north pole (θ=0) ↦ 0, south pole (θ=π) ↦ ∞
    assert sp.limit(W.subs(ph, 0), th, 0) == 0
    assert sp.limit(W.subs(ph, 0), th, sp.pi, "-") is sp.oo
    # the |α|² + |β|² = 1 parametrization gives the same ratio
    a = sp.Symbol("a", positive=True)  # θ/2 ∈ (0, π/2)
    alpha, beta = sp.cos(a), sp.sin(a) * sp.exp(I * ph)
    assert sp.simplify(alpha**2 + sp.Abs(beta) ** 2 - 1) == 0
    assert is_zero((beta / alpha).subs(a, th / 2) - W)

def check_mobius():
    a, b, c, d = sp.symbols("a b c d")
    X = sp.Matrix([[a, b], [c, d]])
    # finite action: w ↦ (c + d w)/(a + b w)
    U = sp.eye(2) + t * X
    alpha, beta = 1, w  # (α, β) with β/α = w
    v = U * sp.Matrix([alpha, beta])
    moved = sp.simplify(sp.diff(v[1] / v[0], t).subs(t, 0))
    assert sp.simplify(moved - dw_dt(X)) == 0
    # same result from the exponentiated one-parameter group
    for M in (s1, s2, s3, sp_, sm_):
        G = sp.simplify(sp.exp(t * M))
        zt = (G[1, 0] + G[1, 1] * w) / (G[0, 0] + G[0, 1] * w)
        assert sp.simplify(sp.diff(zt, t).subs(t, 0) - dw_dt(M)) == 0

def check_pauli_fields():
    for a, M in ((1, s1), (2, s2), (3, s3)):
        assert sp.simplify(dw_dt(-I * M / 2) - HOL[a]) == 0
    # each generator moves the Bloch vector by dn/dt = e_a × n
    axes = {1: (1, 0, 0), 2: (0, 1, 0), 3: (0, 0, 1)}
    for a in (1, 2, 3):
        diff = field_from_bloch(axes[a]) - HOL[a].subs(w, W)
        for pt in POINTS:
            assert num_zero(diff, pt)

def check_ladder():
    assert sp.expand(dw_dt(sp_)) == -(w**2)
    assert sp.expand(dw_dt(sm_)) == 1
    assert sp.expand(dw_dt(s3)) == -2 * w  # -iσ3/2 ↦ i w ∂w
    # σ± do not move the Bloch sphere: no rotation reproduces them
    for M in (sp_, sm_):
        field = dw_dt(M).subs(w, W)
        for axis in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            diff = field_from_bloch(axis) - field
            assert not all(num_zero(diff, pt) for pt in POINTS)

if __name__ == "__main__":
    check_stereographic()
    check_mobius()
    check_pauli_fields()
    check_ladder()
    print("check-memo-mobius: all checks passed")
