"""Checks for article 06 (summary: sorting the properties of complex analysis).

Cl_{n,0}(R) with generators e_0..e_{n-1} (e_a^2 = 1) on the bitmask blades of
common.clifford, as in 04 and 05. Paravector variable q = e_0 x, units h_0 = 1, h_l = e_0 e_l,
Fueter-type operator Dq = e_0 D = sum h_a d_a, right action F Dq = sum (d_a F) h_a.

1. sum_a h_a h_a = 1 - (n-1) = 2 - n, hence Dq q = 2 - n; q^{-1} = x^{-1} e_0 and
   Dq q^{-1} = (n-2)/|q|^2; zeta_l = x_l - x_0 h_l is left and right regular
   (n = 2..6).
2. Powers: the product rule Dq(fg) = (Dq f) g + f (Dq g) + sum_l [h_l, f] d_l g
   (random even f, g, n = 3, 4); for f = g = q:
   (Dq q) q + q (Dq q) = 2(2-n) q, sum_l h_l v h_l = (n-3) v,
   sum_l [h_l, q] h_l = 2(n-2) v, Dq q^2 = 2(2-n) x_0 (n = 2..6);
   n = 4: -4q + 4v = -4 x_0.
3. The even subalgebra Cl^0_{n,0} is commutative only for n = 2.
4. Left and right: for n = 2 and even F, F Dq = Dq F; with D itself,
   z = e_0 x has D z = 0 but z D = 2 e_0. F Dq = (F e_0) D (n = 3, 4).
5. Composition: for n = 2, a left regular even F has d_1 F = h_1 d_0 F, so
   dF = (dx_0 + dx_1 h_1) F' = dq F' without reordering, and conversely;
   the chain rule gives d(F o G) = dq G' F'(G) (checked for polynomials in q).
   For n = 3..5, h_m h_l - h_l h_m = 2 h_m h_l is invertible (so functions with
   d_l F = h_l d_0 F for all l are affine, F = qc + b) and Dq(qc + b) = (2-n) c; for n = 4, zeta_1 = x_1 - x_0 i has
   d_2 zeta_1 = 0 != h_2 d_0 zeta_1 = k; zeta_1(zeta_2) = -x_2 i is not regular.
"""

from common.clifford import MV, eq, grade, mv
from common.clifford import Alg as ClAlg

# ---------------------------------------------------------------- Cl_{n,0}


class Alg(ClAlg):
    """Cl_{n,0} with coordinates x_0..x_{n-1}, D and the paravector operator."""

    def __init__(self, n):
        super().__init__(n)
        self.h = [mv(1)] + [self.e[0] * self.e[l] for l in range(1, n)]
        self.q = self.e[0] * self.x
        self.v = self.q - self.X[0]

    def Dq(self, F):
        return sum((self.h[a] * self.d(F, a) for a in range(self.n)), MV())

    def Dqr(self, F):
        return sum((self.d(F, a) * self.h[a] for a in range(self.n)), MV())

    def rnd_even(self, seed, deg=2):
        return self.rnd(seed, deg, grades=range(0, self.n + 1, 2))


def comm(a, b):
    return a * b - b * a


# ---------------------------------------------------------------- 1.
for n in range(2, 7):
    C = Alg(n)
    assert eq(sum((h * h for h in C.h), MV()), 2 - n)
    assert eq(C.Dq(C.q), 2 - n)
    r2 = sum(t**2 for t in C.X)
    qinv = C.x * (1 / r2) * C.e[0]
    assert eq(qinv * C.q, 1)
    assert eq(C.Dq(qinv), mv((n - 2) / r2))
    for l in range(1, n):
        zeta = C.X[l] - C.X[0] * C.h[l]
        assert eq(C.Dq(zeta), 0) and eq(C.Dqr(zeta), 0)
print("1. sum h_a h_a = 1 - (n-1) = 2 - n, Dq q = 2 - n, Dq q^{-1} = (n-2)/|q|^2, "
      "Fueter variables x_l - x_0 h_l left and right regular (n = 2..6)")

# ---------------------------------------------------------------- 2.
for n in (3, 4):
    C = Alg(n)
    f, g = C.rnd_even(10 + n), C.rnd_even(20 + n)
    rhs = C.Dq(f) * g + f * C.Dq(g) + sum((comm(C.h[l], f) * C.d(g, l) for l in range(1, n)), MV())
    assert eq(C.Dq(f * g), rhs)
for n in range(2, 7):
    C = Alg(n)
    q, v, x0 = C.q, C.v, C.X[0]
    dim_part = C.Dq(q) * q + q * C.Dq(q)
    com_part = sum((comm(C.h[l], q) * C.d(q, l) for l in range(1, n)), MV())
    assert eq(dim_part, 2 * (2 - n) * q)
    assert eq(sum((C.h[l] * v * C.h[l] for l in range(1, n)), MV()), (n - 3) * v)
    assert eq(com_part, 2 * (n - 2) * v)
    assert eq(C.Dq(q * q), dim_part + com_part)
    assert eq(C.Dq(q * q), 2 * (2 - n) * x0)
    if n == 4:
        assert eq(dim_part, -4 * q) and eq(com_part, 4 * v) and eq(C.Dq(q * q), -4 * x0)
print("2. Dq(fg) product rule (n = 3, 4); Dq q^2 = 2(2-n) q + 2(n-2) v = 2(2-n) x0 (n = 2..6)")

# ---------------------------------------------------------------- 3.
for n in range(2, 7):
    even = [MV({b: 1}) for b in range(1 << n) if grade(b) % 2 == 0]
    commutative = all(eq(comm(a, b), 0) for a in even for b in even)
    assert commutative == (n == 2)
print("3. Cl^0_{n,0} is commutative only for n = 2 (n = 2..6)")

# ---------------------------------------------------------------- 4.
C = Alg(2)
F = C.rnd_even(1, deg=3)
assert eq(C.Dqr(F), C.Dq(F))
assert eq(C.D(C.q), 0) and eq(C.Dr(C.q), 2 * C.e[0])
for n in (3, 4):
    C = Alg(n)
    F = C.rnd_even(30 + n)
    assert eq(C.Dqr(F), C.Dr(F * C.e[0]))
print("4. n = 2: F Dq = Dq F for even F; D z = 0, z D = 2 e0; F Dq = (F e0) D (n = 3, 4)")

# ---------------------------------------------------------------- 5.
C = Alg(2)
q, h1 = C.q, C.h[1]
a, b = 2 + 3 * h1, -1 + h1
F = q * q * q * a + q * b + 5  # polynomial in q: left regular for n = 2
assert eq(C.Dq(F), 0)
assert eq(C.d(F, 1), h1 * C.d(F, 0))
# dF = dq F' without reordering: d_0 F = F', d_1 F = h_1 F'
Fp = C.d(F, 0)
assert eq(C.d(F, 1), h1 * Fp)
# converse: dF = dq c gives Dq F = c + h_1 h_1 c = 0 for any c (no commutativity used)
cc = C.rnd_even(7)
assert eq(cc + h1 * h1 * cc, 0)
# chain rule: G regular with paravector values, F o G regular, d(F o G) = dq G' F'(G)
G = q * q * (1 - 2 * h1) + q * 3
G0, G1 = G.d.get(0, 0), G.d.get(0b11, 0)  # coefficients of 1 and h_1 = e_0 e_1
assert eq(G, G0 + G1 * h1)
Fsub = lambda H: H * H * H * a + H * b + 5
FG = Fsub(G)
assert eq(C.Dq(FG), 0)
Fp_at_G = 3 * G * G * a + b  # F'(w) = 3 w^2 a + b (commutative in 2D)
assert eq(C.d(F, 0), 3 * q * q * a + b)
assert eq(C.d(FG, 0), C.d(G, 0) * Fp_at_G) and eq(C.d(FG, 1), h1 * C.d(G, 0) * Fp_at_G)
for n in (3, 4):
    C = Alg(n)
    i, j = C.h[1], C.h[2]
    zeta1 = C.X[1] - C.X[0] * i
    assert eq(C.Dq(zeta1), 0)
    # d_2 zeta_1 = 0, while h_2 d_0 zeta_1 = -ji (= k after the P_+ projection for n = 4)
    assert eq(C.d(zeta1, 2), 0) and eq(C.h[2] * C.d(zeta1, 0), -j * i) and not eq(-j * i, 0)
    if n == 4:
        # zeta_1 composed with zeta_2 = x_2 - x_0 j: real part x_2, i-part 0
        comp = -C.X[2] * i
        assert not eq(C.Dq(comp), 0) and eq(C.Dq(comp), -j * i)
# functions with d_l F = h_l d_0 F for all l: second derivatives force (h_m h_l - h_l h_m) F'' = 0
for n in (3, 4, 5):
    C = Alg(n)
    for l in range(1, n):
        for m in range(1, n):
            if l != m:
                c = C.h[m] * C.h[l] - C.h[l] * C.h[m]
                assert eq(c, 2 * C.h[m] * C.h[l]) and eq(c * c, -4)  # invertible
    cc = C.rnd_even(40 + n, deg=0)
    assert eq(C.Dq(C.q * cc + 7), (2 - n) * cc)
    for l in range(1, n):
        assert eq(C.d(C.q * cc, l), C.h[l] * C.d(C.q * cc, 0))
print("5a. n = 2: dF = dq F', chain rule d(F o G) = dq G' F'(G); n = 3..5: h_m h_l - h_l h_m invertible, "
      "Dq(qc + b) = (2-n) c")
print("5b. n = 2: d_1 F = h_1 d_0 F for regular F; n = 3, 4: d_2 zeta_1 = 0 != h_2 d_0 zeta_1; "
      "Dq(zeta_1 o zeta_2) = -ji != 0 (n = 4)")
