"""Checks for hopf/MEMO.md (C² 上の線形ベクトル場としての実現).

The main line of the memo: the generators of SU(2) act on the spinor space
C² linearly, so they are first-order differential operators in (α, β)
directly, with no ratio, no stereographic projection and no Riemann sphere.

Verifies symbolically (sympy):

1. A matrix X acting linearly on v = (α, β) generates the differential
   operator D_X = Σ X_ik v_k ∂_i, checked against exp(tX).
   D is an ANTI-homomorphism: [D_X, D_Y] = -D_[X,Y]. Transposing fixes the
   direction: E_X := D_{Xᵀ} satisfies [E_X, E_Y] = E_[X,Y].
2. The Pauli matrices give J+ = α∂β, J- = β∂α, J3 = (α∂α - β∂β)/2 with
   [J3, J±] = ±J±, [J+, J-] = 2J3.
3. The fiber action (α,β) ↦ (e^{iψ}α, e^{iψ}β) has generator i·N with
   N = α∂α + β∂β the Euler operator. N commutes with J±, J3, so the
   fiber charge is exactly the polynomial degree.
4. The homogeneous polynomials of degree 2j = N are closed under J±, J3,
   have dimension 2j+1, carry the weight basis (J3 eigenvalue j-k), and
   the Casimir J² = j(j+1) (j = 1/2, 1, 3/2, 2). Half-integer j needs no
   special treatment: only 2j appears.
5. j = 1/2 reproduces the qubit: in the basis {α, β} the matrices of
   J3, J+, J- are σ3/2, σ+, σ-.
6. Bridge to the ratio section: writing F(α,β) = α^{2j} g(w), w = β/α,
   pushes J± , J3 onto g as ∂w, 2jw - w²∂w, j - w∂w. The ratio form is a
   presentation of the same operators, not a separate construction.
"""

import sympy as sp

I = sp.I
al, be = sp.symbols("alpha beta")
w = sp.Symbol("w")
t = sp.Symbol("t")
psi = sp.Symbol("psi", real=True)
V = (al, be)

s1 = sp.Matrix([[0, 1], [1, 0]])
s2 = sp.Matrix([[0, -I], [I, 0]])
s3 = sp.Matrix([[1, 0], [0, -1]])
sp_ = sp.expand((s1 + I * s2) / 2)
sm_ = sp.expand((s1 - I * s2) / 2)

SPINS = [sp.Rational(1, 2), sp.Integer(1), sp.Rational(3, 2), sp.Integer(2)]


def D(X, f):
    """Generator of the linear action v ↦ exp(tX) v, acting on functions."""
    return sp.expand(sum(X[i, k] * V[k] * sp.diff(f, V[i])
                         for i in range(2) for k in range(2)))


def E(X, f):
    """Jordan-Schwinger map: D with X transposed, a genuine homomorphism."""
    return D(X.T, f)


def comm(A, B, f):
    return sp.expand(A(B(f)) - B(A(f)))


# --- the su(2) generators as differential operators on C² -----------------

def Jp(f):
    return sp.expand(al * sp.diff(f, be))


def Jm(f):
    return sp.expand(be * sp.diff(f, al))


def J3(f):
    return sp.expand((al * sp.diff(f, al) - be * sp.diff(f, be)) / 2)


def N(f):
    """Euler operator: counts the degree."""
    return sp.expand(al * sp.diff(f, al) + be * sp.diff(f, be))


def casimir(f):
    return sp.expand(J3(J3(f)) + (Jp(Jm(f)) + Jm(Jp(f))) / 2)


SAMPLES = [sp.Integer(1), al, be, al**2, al * be, be**3, al**2 * be,
           sp.Function("f")(al, be)]


def check_linear_action():
    """D_X is the generator of v ↦ exp(tX) v; it is an anti-homomorphism."""
    a, b, c, d = sp.symbols("a b c d")
    X = sp.Matrix([[a, b], [c, d]])
    f = sp.Function("f")(al, be)
    # moving the argument by exp(tX) and differentiating at t = 0
    for M in (s1, s2, s3, sp_, sm_, X):
        vt = (sp.eye(2) + t * M) * sp.Matrix([al, be])
        moved = f.subs({al: sp.Symbol("A"), be: sp.Symbol("B")}).subs(
            {sp.Symbol("A"): vt[0], sp.Symbol("B"): vt[1]})
        got = sp.expand(sp.diff(moved, t).subs(t, 0))
        assert sp.simplify(got - D(M, f)) == 0
    # [D_X, D_Y] = -D_[X,Y]  (anti-homomorphism)
    for X0, Y0 in ((s1, s2), (s3, sp_), (sp_, sm_)):
        XY = sp.expand(X0 * Y0 - Y0 * X0)
        lhs = sp.expand(D(X0, D(Y0, f)) - D(Y0, D(X0, f)))
        assert sp.simplify(lhs + D(XY, f)) == 0
    # transposing restores the homomorphism
    for X0, Y0 in ((s1, s2), (s3, sp_), (sp_, sm_)):
        XY = sp.expand(X0 * Y0 - Y0 * X0)
        lhs = sp.expand(E(X0, E(Y0, f)) - E(Y0, E(X0, f)))
        assert sp.simplify(lhs - E(XY, f)) == 0


def check_pauli_operators():
    """E of the Pauli matrices gives J+, J-, J3 with the right brackets."""
    f = sp.Function("f")(al, be)
    assert sp.simplify(E(sp_, f) - Jp(f)) == 0
    assert sp.simplify(E(sm_, f) - Jm(f)) == 0
    assert sp.simplify(E(s3, f) / 2 - J3(f)) == 0
    for g in SAMPLES:
        assert sp.simplify(comm(J3, Jp, g) - Jp(g)) == 0
        assert sp.simplify(comm(J3, Jm, g) + Jm(g)) == 0
        assert sp.simplify(comm(Jp, Jm, g) - 2 * J3(g)) == 0


def check_fiber_is_degree():
    """The fiber generator is i·N, and N is the degree."""
    f = sp.Function("f")(al, be)
    A, B = sp.symbols("A B")
    moved = f.subs({al: A, be: B}).subs(
        {A: sp.exp(I * psi) * al, B: sp.exp(I * psi) * be})
    got = sp.expand(sp.diff(moved, psi).subs(psi, 0))
    assert sp.simplify(got - I * N(f)) == 0
    # N counts the degree, and commutes with the su(2) generators
    for jval in SPINS:
        for f0 in basis(jval):
            assert sp.expand(N(f0) - 2 * jval * f0) == 0
    for g in SAMPLES:
        for A_ in (Jp, Jm, J3):
            assert sp.simplify(comm(N, A_, g)) == 0


def basis(jval):
    """Monomials of degree 2j, ordered by descending J3 weight."""
    n = int(2 * jval)
    return [al ** (n - k) * be**k for k in range(n + 1)]


def check_spin_j():
    """Degree 2j is closed, has dimension 2j+1, weights j-k, Casimir j(j+1)."""
    for jval in SPINS:
        b = basis(jval)
        assert len(b) == 2 * jval + 1
        for k, f in enumerate(b):
            assert sp.expand(J3(f) - (jval - k) * f) == 0
            assert sp.expand(casimir(f) - jval * (jval + 1) * f) == 0
            for image in (Jp(f), Jm(f)):
                assert image == 0 or sp.expand(N(image) - 2 * jval * image) == 0
        assert Jp(b[0]) == 0    # J+ |j, j> = 0
        assert Jm(b[-1]) == 0   # J- |j, -j> = 0


def matrix_of(op, jval):
    b = basis(jval)
    M = sp.zeros(len(b), len(b))
    for col, f in enumerate(b):
        img = sp.expand(op(f))
        for row, g in enumerate(b):
            M[row, col] = img.coeff(g) if img != 0 else 0
    return M


def check_qubit():
    """j = 1/2 in the basis {α, β} is literally the qubit."""
    half = sp.Rational(1, 2)
    assert matrix_of(J3, half) == sp.diag(half, -half)
    assert matrix_of(Jp, half) == sp.Matrix([[0, 1], [0, 0]])   # σ+
    assert matrix_of(Jm, half) == sp.Matrix([[0, 0], [1, 0]])   # σ-


def check_ratio_bridge():
    """F = α^{2j} g(w) pushes the same operators onto g in the ratio form."""
    j = sp.Symbol("j", positive=True)
    g = sp.Function("g")(w)

    def push(op):
        """Apply op to F = α^{2j} g(β/α), then strip α^{2j}."""
        F = al ** (2 * j) * g.subs(w, be / al)
        img = sp.simplify(op(F) / al ** (2 * j))
        return sp.simplify(img.subs(be, w * al).doit())

    assert sp.simplify(push(Jp) - sp.diff(g, w)) == 0
    assert sp.simplify(push(J3) - (j * g - w * sp.diff(g, w))) == 0
    assert sp.simplify(push(Jm) - (2 * j * w * g - w**2 * sp.diff(g, w))) == 0
    assert sp.simplify(push(N) - 2 * j * g) == 0


if __name__ == "__main__":
    check_linear_action()
    check_pauli_operators()
    check_fiber_is_degree()
    check_spin_j()
    check_qubit()
    check_ratio_bridge()
    print("check-memo-c2: all checks passed")
