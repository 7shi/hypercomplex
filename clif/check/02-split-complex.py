"""Checks for 02-split-complex.md.

Verifies the split-complex algebra R[j] (j^2 = +1) and the split-quaternion
section at the end of the article:

- the diagonal representation j = diag(1,-1), z = diag(x+y, x-y), the
  idempotents e, e* as the projections diag(1,0), diag(0,1), and their
  idempotency / orthogonality / completeness;
- the only idempotents of R[j] are 0, 1, e, e*, the zero divisors are
  exactly the x(1 +- j), and the isomorphism R[j] = R (+) R given by
  x + jy -> (x+y, x-y);
- the non-diagonal representation B = (0 1; 1 0) is the regular representation
  of "multiply by j" on the basis {1, j}, has eigenvalues +-1 and is similar to
  the diagonal one; it acts as (x, y) -> (y, x), the reflection in y = x;
- the same construction on C gives i = (0 -1; 1 0), acting as (x, y) -> (-y, x),
  the +90 degree rotation;
- ij = diag(-1,1) and ji = diag(1,-1), so ij = -ji, and k = ij satisfies
  k^2 = +1; together i^2 = -1, j^2 = k^2 = +1, ijk = +1, jk = -i and kj = i
  (so j and k anticommute, as Clifford generators must);
- a + bi + cj + dk = (a-d, -b+c; b+c, a+d) runs over all of M_2(R), so
  H' = M_2(R) = R(2) in the notation of 01;
- Pauli matrices: sigma_1 = j, sigma_3 = -k, -i sigma_2 = i, so (j, i, -k) are
  the real Pauli matrices (tau_1, tau_2, tau_3).
"""

import sympy as sp

# --- split-complex numbers as diagonal matrices ---------------------------------
x, y = sp.symbols("x y", real=True)
J = sp.Matrix([[1, 0], [0, -1]])
I2 = sp.eye(2)

print("j^2 = 1:", sp.simplify(J**2 - I2) == sp.zeros(2))

z = x * I2 + y * J
print("z = diag(x+y, x-y):", z == sp.diag(x + y, x - y))

e = (I2 + J) / 2
es = (I2 - J) / 2
print("e = diag(1,0), e* = diag(0,1):", e == sp.diag(1, 0) and es == sp.diag(0, 1))
print("idempotent:", e**2 == e and es**2 == es)
print("orthogonal:", e * es == sp.zeros(2) and es * e == sp.zeros(2))
print("complete:", e + es == I2)

# eigenvectors: je = e, je* = -e*
print("je = e, je* = -e*:", J * e == e and J * es == -es)

# --- idempotents of R[j], symbolically ------------------------------------------
sols = sp.solve([x**2 + y**2 - x, 2 * x * y - y], [x, y], dict=True)
got = {(s[x], s[y]) for s in sols}
print("idempotents are 0, 1, e, e*:",
      got == {(0, 0), (1, 0), (sp.Rational(1, 2), sp.Rational(1, 2)),
              (sp.Rational(1, 2), -sp.Rational(1, 2))})

# zero divisors: zz* = x^2 - y^2 = 0 <=> y = +-x, i.e. z = x(1 +- j)
conj = x * I2 - y * J
print("zz* = x^2 - y^2:", sp.expand(z * conj) == sp.expand((x**2 - y**2) * I2))
print("zz* = 0 <=> y = +-x:",
      sp.solve(sp.Eq(x**2 - y**2, 0), y) == [-x, x])
print("(1+j)(1-j) = 0 with both factors nonzero:",
      (I2 + J) * (I2 - J) == sp.zeros(2) and I2 + J != sp.zeros(2))

# --- ring isomorphism R[j] = R (+) R ---------------------------------------------
a, b, c, d = sp.symbols("a b c d", real=True)
z1, z2 = x * I2 + y * J, c * I2 + d * J
Phi = lambda m: (m[0, 0], m[1, 1])  # diag entries = (x+y, x-y)
p1, p2, pp = Phi(z1), Phi(z2), Phi(sp.expand(z1 * z2))
print("Phi is multiplicative:",
      sp.expand(pp[0] - p1[0] * p2[0]) == 0 and sp.expand(pp[1] - p1[1] * p2[1]) == 0)
print("Phi is additive:", Phi(z1 + z2) == (p1[0] + p2[0], p1[1] + p2[1]))

# --- the regular representation ---------------------------------------------------
def rep(images):
    """representation matrix: columns are the images of the basis vectors."""
    return sp.Matrix(images).T

# basis {1, j}: j*1 = j = (0,1), j*j = 1 = (1,0)
B = rep([(0, 1), (1, 0)])
print("j = (0 1; 1 0) on the basis {1, j}:", B == sp.Matrix([[0, 1], [1, 0]]))
print("B^2 = I, eigenvalues +-1:", B**2 == I2 and sorted(B.eigenvals()) == [-1, 1])
print("B is similar to the diagonal representation:",
      B.diagonalize()[1] in (sp.diag(1, -1), sp.diag(-1, 1)))
print("B acts as (x, y) -> (y, x), reflection in y = x:",
      B * sp.Matrix([x, y]) == sp.Matrix([y, x]) and B.det() == -1)

# basis {1, i}: i*1 = i = (0,1), i*i = -1 = (-1,0)
i_ = rep([(0, 1), (-1, 0)])
rot = lambda t: sp.Matrix([[sp.cos(t), -sp.sin(t)], [sp.sin(t), sp.cos(t)]])
print("i = (0 -1; 1 0) on the basis {1, i}:", i_ == sp.Matrix([[0, -1], [1, 0]]))
print("i acts as (x, y) -> (-y, x), the +90 degree rotation:",
      i_ * sp.Matrix([x, y]) == sp.Matrix([-y, x]) and i_ == rot(sp.pi / 2))

# --- split quaternions ----------------------------------------------------------
j_ = B
k_ = i_ * j_
print("ij = diag(-1, 1), ji = diag(1, -1):",
      i_ * j_ == sp.diag(-1, 1) and j_ * i_ == sp.diag(1, -1))
print("ij = -ji (their sum vanishes):", i_ * j_ + j_ * i_ == sp.zeros(2))
print("i^2 = -1, j^2 = k^2 = 1:", i_**2 == -I2 and j_**2 == I2 and k_**2 == I2)
print("ijk = +1 (while ijk = -1 in H):", i_ * j_ * k_ == I2)
print("jk = (0 1; -1 0) = -i:", j_ * k_ == sp.Matrix([[0, 1], [-1, 0]]) == -i_)
print("kj = i, so jk = -kj (j and k anticommute):", k_ * j_ == i_ and j_ * k_ == -(k_ * j_))

a, b, c, d = sp.symbols("a b c d", real=True)
q = a * I2 + b * i_ + c * j_ + d * k_
print("a + bi + cj + dk = (a-d, -b+c; b+c, a+d):",
      q == sp.Matrix([[a - d, -b + c], [b + c, a + d]]))

# the four matrices are linearly independent over R, so H' = M_2(R)
M = sp.Matrix([list(m) for m in (I2, i_, j_, k_)])  # rows = flattened matrices
print("1, i, j, k span M_2(R), hence H' = R(2):", M.rank() == 4)

# --- Pauli matrices --------------------------------------------------------------
s1 = sp.Matrix([[0, 1], [1, 0]])
s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
s3 = sp.Matrix([[1, 0], [0, -1]])
print("sigma_1 = j, sigma_3 = -k:", s1 == j_ and s3 == -k_)
print("-i sigma_2 = i:", -sp.I * s2 == i_)

# real Pauli matrices tau_1, tau_2, tau_3 = (j, i, -k)
t1, t2, t3 = s1, -sp.I * s2, s3
print("(tau_1, tau_2, tau_3) = (j, i, -k):", (t1, t2, t3) == (j_, i_, -k_))
print("tau_1^2 = tau_3^2 = I, tau_2^2 = -I:",
      t1**2 == I2 and t3**2 == I2 and t2**2 == -I2)
print("tau_1 tau_2 = tau_3:", t1 * t2 == t3)
