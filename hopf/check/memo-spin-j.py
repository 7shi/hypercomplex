"""Checks for hopf/MEMO.md (比で書いた場合 / ハイゼンベルク代数への縮約).

The memo's main line lives on C² (see memo-c2.py); this file checks the
ratio presentation F(α,β) = α^{2j} g(w), w = β/α, which the memo keeps as
a single section because it is where the Möbius view and the j → ∞
contraction land. memo-c2.py verifies that the two presentations agree.

Verifies symbolically (sympy):

1. J+ = ∂w, J3 = j - w∂w, J- = 2jw - w²∂w satisfy [J3,J±] = ±J±,
   [J+,J-] = 2J3 and the Casimir J² = j(j+1) on the polynomials of
   degree ≤ 2j (j = 1/2, 1, 3/2, 2).
2. The weight basis: J3 w^k = (j-k) w^k, J+ kills the top (k=0) and J-
   kills the bottom (k=2j), so the space is (2j+1)-dimensional and closed.
3. j = 1/2 reproduces the qubit: in the basis {1, w} ↔ {α, β} the matrices
   of J3, J+, J- are σ3/2, σ+, σ-.
4. Contraction j → ∞ with u = √(2j) w: J+/√(2j) = ∂u, J-/√(2j) → u and
   [J+,J-]/(2j) → 1, i.e. the Heisenberg algebra [a, a†] = 1 realized as
   differentiation and multiplication (Bargmann realization).
"""

import sympy as sp

w = sp.Symbol("w")
u = sp.Symbol("u")
j = sp.Symbol("j", positive=True)

def Jp(f, j=j):
    return sp.expand(sp.diff(f, w))

def J3(f, j=j):
    return sp.expand(j * f - w * sp.diff(f, w))

def Jm(f, j=j):
    return sp.expand(2 * j * w * f - w**2 * sp.diff(f, w))

def comm(A, B, f, j=j):
    return sp.expand(A(B(f, j), j) - B(A(f, j), j))

def casimir(f, j=j):
    return sp.expand(J3(J3(f, j), j) + (Jp(Jm(f, j), j) + Jm(Jp(f, j), j)) / 2)

SPINS = [sp.Rational(1, 2), sp.Integer(1), sp.Rational(3, 2), sp.Integer(2)]

def basis(jval):
    return [w**k for k in range(int(2 * jval) + 1)]

def check_algebra():
    g = sp.Function("g")(w)
    for f in (g, sp.Integer(1), w, w**2, w**3):
        assert sp.simplify(comm(J3, Jp, f) - Jp(f)) == 0
        assert sp.simplify(comm(J3, Jm, f) + Jm(f)) == 0
        assert sp.simplify(comm(Jp, Jm, f) - 2 * J3(f)) == 0
    for jval in SPINS:
        for f in basis(jval):
            assert sp.expand(casimir(f, jval) - jval * (jval + 1) * f) == 0

def check_weight_basis():
    for jval in SPINS:
        b = basis(jval)
        for k, f in enumerate(b):
            assert sp.expand(J3(f, jval) - (jval - k) * f) == 0
            # closed within the space (degrees stay within 0..2j)
            for image in (Jp(f, jval), Jm(f, jval)):
                assert image == 0 or sp.Poly(image, w).degree() <= 2 * jval
        assert Jp(b[0], jval) == 0      # top: J+ |j, j> = 0
        assert Jm(b[-1], jval) == 0     # bottom: J- |j, -j> = 0

def matrix_of(op, jval):
    b = basis(jval)
    M = sp.zeros(len(b), len(b))
    for col, f in enumerate(b):
        poly = sp.Poly(op(f, jval), w)
        for row in range(len(b)):
            M[row, col] = poly.coeff_monomial(w**row)
    return M

def check_qubit():
    half = sp.Rational(1, 2)
    # basis {1, w} corresponds to {α, β} through F = α^{2j} g(w)
    assert matrix_of(J3, half) == sp.diag(half, -half)
    assert matrix_of(Jp, half) == sp.Matrix([[0, 1], [0, 0]])   # σ+
    assert matrix_of(Jm, half) == sp.Matrix([[0, 0], [1, 0]])   # σ-
    # for j = 1 the ladder gives the usual √2 factors
    M = matrix_of(Jp, sp.Integer(1))
    assert M == sp.Matrix([[0, 1, 0], [0, 0, 2], [0, 0, 0]])

def check_contraction():
    s = sp.sqrt(2 * j)
    sub = {w: u / s}

    def rescale(f):
        """Rewrite an operator image in the variable u = √(2j) w."""
        return sp.expand(sp.simplify(f.subs(sub)))

    for f in (sp.Integer(1), u, u**2, u**3):
        g = f.subs(u, s * w)  # same function written in w
        a_f = rescale(Jp(g)) / s          # J+/√(2j)
        ad_f = rescale(Jm(g)) / s         # J-/√(2j)
        assert sp.simplify(a_f - sp.diff(f, u)) == 0            # exactly ∂u
        assert sp.limit(sp.simplify(ad_f), j, sp.oo) == sp.expand(u * f)
        # [J+, J-]/(2j) = 1 - (w ∂w)/j → 1
        c = sp.simplify(rescale(comm(Jp, Jm, g)) / (2 * j))
        assert sp.limit(c / f, j, sp.oo) == 1

if __name__ == "__main__":
    check_algebra()
    check_weight_basis()
    check_qubit()
    check_contraction()
    print("check-memo-spin-j: all checks passed")
