"""Checks for tensor-from-complex.md.

C is represented by real 2x2 matrices (i -> [[0,-1],[1,0]]), and the
real tensor product C (x) C by the Kronecker product.

1. i(x)1 and 1(x)i commute, square to -1, and differ from +-each other;
   {1, i, j, ij} = {1(x)1, i(x)1, 1(x)i, i(x)i} is linearly independent
   (so ij is not of the form a+bi+cj).
2. The product formula (a+bi+cj)(d+ei+fj) of the article, both in the
   bicomplex basis and via the tensor-product computation.
3. (ij)^2 = 1, the zero divisors (1+ij)(1-ij) = 0, and that 1+ij is not
   invertible (determinant 0); the step "ij=1 -> j=-i".
4. The multiplication table (i(x)1)^2, (1(x)i)^2, (i(x)1)(1(x)i),
   (i(x)i)^2, and the rule (a(x)b)(c(x)d) = ac(x)bd for complex a..d.
5. Bilinearity: real coefficients move freely across (x).
6. (a+bi)(c+dj) = ac + ad j + bc i + bd ij matches
   (a+bi)(x)(c+di) under 1->1(x)1, i->i(x)1, j->1(x)i, ij->i(x)i.
"""
import itertools
import random
import sympy as sp

I2 = sp.eye(2)
Ic = sp.Matrix([[0, -1], [1, 0]])  # complex i


def cx(z):
    """complex number (sympy expr a+bI) -> 2x2 real matrix"""
    a, b = sp.re(z), sp.im(z)
    return a * I2 + b * Ic


def kron(A, B):
    return sp.kronecker_product(A, B)


ONE = kron(I2, I2)
Ei = kron(Ic, I2)   # i  <-> i (x) 1
Ej = kron(I2, Ic)   # j  <-> 1 (x) i
Eij = kron(Ic, Ic)  # ij <-> i (x) i
Z = sp.zeros(4, 4)

ok = True


def check(name, cond):
    global ok
    print(("OK  " if cond else "NG  ") + name)
    ok &= bool(cond)


# 1
check("i^2 = -1", Ei * Ei == -ONE)
check("j^2 = -1", Ej * Ej == -ONE)
check("ij = ji", Ei * Ej == Ej * Ei)
check("i != j, i != -j", Ei != Ej and Ei != -Ej)
check("ij = i(x)i", Ei * Ej == Eij)
M = sp.Matrix.hstack(*[X.reshape(16, 1) for X in (ONE, Ei, Ej, Eij)])
check("{1,i,j,ij} linearly independent", M.rank() == 4)

# 2
a, b, c, d, e, f = sp.symbols("a b c d e f", real=True)
lhs = (a * ONE + b * Ei + c * Ej) * (d * ONE + e * Ei + f * Ej)
rhs = (a*d - b*e - c*f) * ONE + (a*e + b*d) * Ei + (a*f + c*d) * Ej \
    + (b*f + c*e) * Eij
check("product formula (bicomplex)", sp.expand(lhs - rhs) == Z)
# tensor-product computation term by term: (x (x) y) with complex x,y
L = [(a, 1, 1), (b * sp.I, 1, 1), (c, 1, sp.I)]   # a(x)1 + bi(x)1 + c(x)i
R = [(d, 1, 1), (e * sp.I, 1, 1), (f, 1, sp.I)]
tot = Z
for (p, _, q), (r, _, s) in itertools.product(L, R):
    tot += kron(cx(p * r), cx(q * s))   # (p(x)q)(r(x)s) = pr (x) qs
check("product formula via tensor product", sp.expand(tot - rhs) == Z)

# 3
check("(ij)^2 = 1", Eij * Eij == ONE)
check("(1+ij)(1-ij) = 0", (ONE + Eij) * (ONE - Eij) == Z)
check("1+ij, 1-ij nonzero", ONE + Eij != Z and ONE - Eij != Z)
check("det(1+ij) = 0 (not invertible)", (ONE + Eij).det() == 0)
# ij = 1  ->  multiply by i:  i*ij = -j  and  i*1 = i  ->  j = -i
check("i*(ij) = -j", Ei * Eij == -Ej)

# 4
check("(i(x)1)^2 = -(1(x)1)", Ei * Ei == -ONE)
check("(1(x)i)^2 = -(1(x)1)", Ej * Ej == -ONE)
check("(i(x)1)(1(x)i) = i(x)i", Ei * Ej == Eij)
check("(i(x)i)^2 = 1(x)1", Eij * Eij == ONE)
random.seed(1)
good = True
for _ in range(20):
    al, be, ga, de = [random.randint(-5, 5) + random.randint(-5, 5) * sp.I
                      for _ in range(4)]
    good &= kron(cx(al), cx(be)) * kron(cx(ga), cx(de)) \
        == kron(cx(sp.expand(al * ga)), cx(sp.expand(be * de)))
check("(al(x)be)(ga(x)de) = al ga (x) be de", good)

# 5
x, y = sp.symbols("x y", real=True)
X, Y = cx(x + y * sp.I), cx(y - x * sp.I)
k = sp.Symbol("k", real=True)
check("kx(x)y = x(x)ky = k(x(x)y)",
      kron(k * X, Y) == kron(X, k * Y) == k * kron(X, Y))

# 6
lhs6 = (a * ONE + b * Ei) * (c * ONE + d * Ej)
rhs6 = kron(cx(a + b * sp.I), cx(c + d * sp.I))
check("(a+bi)(c+dj) <-> (a+bi)(x)(c+di)", sp.expand(lhs6 - rhs6) == Z)

print("ALL OK" if ok else "SOME CHECKS FAILED")
