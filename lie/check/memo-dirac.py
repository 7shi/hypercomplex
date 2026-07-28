"""Checks for lie/MEMO.md (ラプラシアンの平方根の節).

Verifies symbolically (sympy):

1. Flat case: D = σ·∇ on R³ satisfies D² = Δ on 2-component functions,
   which is the Pauli identity (σ·a)(σ·b) = a·b + iσ·(a×b) with a = b = ∇.
2. Curved case: on the unit S² with zweibein (∂θ, ∂φ/sinθ) and spin
   connection ∇θ = ∂θ, ∇φ = ∂φ - (i/2)cosθ σ3, the Dirac operator
   𝒟 = -i(σ1 ∇θ + σ2 ∇φ/sinθ) = -i(σ1(∂θ + cotθ/2) + σ2 ∂φ/sinθ)
   satisfies the Lichnerowicz formula 𝒟² = -Δ_spinor + R/4 with R = 2.
3. The spectrum ±(l+1) of 𝒟 (quoted, not derived here) corresponds to
   the half-integer j = l + 1/2 via (l+1)² = j(j+1) + 1/4.
"""

import sympy as sp

I = sp.I
s1 = sp.Matrix([[0, 1], [1, 0]])
s2 = sp.Matrix([[0, -I], [I, 0]])
s3 = sp.Matrix([[1, 0], [0, -1]])

# --- flat R³ -----------------------------------------------------------------

x, y, z = sp.symbols("x y z", real=True)

def check_flat():
    psi = sp.Matrix([sp.Function("f")(x, y, z), sp.Function("g")(x, y, z)])

    def D(p):
        return sum(
            (s * sp.Matrix([sp.diff(c, v) for c in p]) for s, v in
             ((s1, x), (s2, y), (s3, z))),
            sp.zeros(2, 1),
        )

    lap = sp.Matrix([sum(sp.diff(c, v, 2) for v in (x, y, z)) for c in psi])
    assert sp.simplify(D(D(psi)) - lap) == sp.zeros(2, 1)

# --- unit S² -----------------------------------------------------------------

th, ph = sp.symbols("theta phi", real=True)

def nabla_th(p):
    return sp.Matrix([sp.diff(c, th) for c in p])

def nabla_ph(p):
    return sp.Matrix([sp.diff(c, ph) for c in p]) - I * sp.cos(th) / 2 * (s3 * p)

def dirac(p):
    return -I * (s1 * nabla_th(p) + s2 * nabla_ph(p) / sp.sin(th))

def lap_spinor(p):
    return (
        nabla_th(nabla_th(p))
        + sp.cot(th) * nabla_th(p)
        + nabla_ph(nabla_ph(p)) / sp.sin(th) ** 2
    )

def check_curved():
    psi = sp.Matrix([sp.Function("f")(th, ph), sp.Function("g")(th, ph)])
    # the spin connection is already absorbed into the cotθ/2 term
    naive = s1 * sp.Matrix(
        [sp.diff(c, th) + sp.cot(th) * c / 2 for c in psi]
    ) + s2 * sp.Matrix([sp.diff(c, ph) for c in psi]) / sp.sin(th)
    assert sp.simplify(dirac(psi) + I * naive) == sp.zeros(2, 1)
    # Lichnerowicz: 𝒟² = -Δ_spinor + R/4, R = 2 for the unit sphere
    scal = sp.Integer(2)
    lhs = dirac(dirac(psi))
    rhs = -lap_spinor(psi) + scal / 4 * psi
    assert sp.simplify(sp.expand(lhs - rhs)) == sp.zeros(2, 1)

# --- half-integer bookkeeping ------------------------------------------------

def check_half_integer():
    l = sp.Symbol("l", nonnegative=True, integer=True)
    j = l + sp.Rational(1, 2)
    assert sp.simplify((l + 1) ** 2 - (j * (j + 1) + sp.Rational(1, 4))) == 0

if __name__ == "__main__":
    check_flat()
    check_curved()
    check_half_integer()
    print("check-memo-dirac: all checks passed")
