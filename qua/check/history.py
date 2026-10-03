"""Checks for history.md.

Verifies the formulas quoted in the historical overview:

1. Hamilton's relations i^2 = j^2 = k^2 = ijk = -1 and the derived
   products ij = k, jk = i, ki = j, ji = -k, kj = -i, ik = -j.
2. The product of pure quaternions p q = -p.q + p x q, the square
   q^2 = -(x^2 + y^2 + z^2) and q q* = x^2 + y^2 + z^2 > 0.
3. Hamilton's nabla (commuting partial symbols): nabla^2 is minus the
   Laplacian.
4. Gibbs/Heaviside separation: p.q = -(pq + qp)/2, p x q = (pq - qp)/2.
5. Rotation by a unit quaternion: q* = q^{-1}, q v q* is pure, has the
   same length as v and agrees with the Rodrigues formula for the angle
   theta when q = cos(theta/2) + n sin(theta/2); a non-unit q scales by
   |q|^2 (negative check).
6. Pauli matrices: sigma_a^2 = I (so they do not satisfy Hamilton's
   relations themselves), while i -> -i sigma_1, j -> -i sigma_2,
   k -> -i sigma_3 does, and the map is multiplicative on the basis.
7. Cl_{3,0}(R) through Pauli matrices: vectors square to |v|^2 > 0,
   p q = p.q + p^q (scalar + bivector), the bivectors square to -1 and
   -e2e3, -e3e1, -e1e2 satisfy Hamilton's relations; the rotor sandwich
   R v R~ with R R~ = 1 equals the quaternion sandwich q v q*.
8. Cl_{0,2}(R): two anticommuting grade-1 generators with square -1 and
   their product give the quaternion units.
9. 4D: x -> a x b^{-1} for unit a, b preserves the norm and has
   determinant +1; a = b fixes 1.
10. The number of coordinate planes n(n-1)/2 equals n only for n = 3
    (n >= 1).
11. Pseudovector law: (Qp) x (Qq) = det(Q) Q (p x q) for orthogonal Q,
    with an extra minus sign for a reflection (det Q = -1).
"""

import random

import sympy as sp

# --- quaternions: 4-tuples (w, x, y, z) --------------------------------------

def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return (w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
            w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2)

def qconj(a):
    return (a[0], -a[1], -a[2], -a[3])

def qadd(a, b):
    return tuple(s + t for s, t in zip(a, b))

def qsmul(c, a):
    return tuple(c * s for s in a)

def qeq(a, b):
    return all(sp.simplify(s - t) == 0 for s, t in zip(a, b))

ONE = (1, 0, 0, 0)
I, J, K = (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)
neg = lambda a: qsmul(-1, a)

def check(name, cond):
    print(("OK   " if cond else "FAIL ") + name)
    if not cond:
        raise SystemExit(1)

# 1. Hamilton's relations
check("i^2 = j^2 = k^2 = -1",
      all(qmul(u, u) == neg(ONE) for u in (I, J, K)))
check("ijk = -1", qmul(qmul(I, J), K) == neg(ONE))
check("ij = k, jk = i, ki = j",
      qmul(I, J) == K and qmul(J, K) == I and qmul(K, I) == J)
check("ji = -k, kj = -i, ik = -j",
      qmul(J, I) == neg(K) and qmul(K, J) == neg(I) and qmul(I, K) == neg(J))

# 2. product of pure quaternions
p1, p2, p3, q1, q2, q3, x, y, z = sp.symbols("p1 p2 p3 q1 q2 q3 x y z", real=True)
P = (0, p1, p2, p3)
Q = (0, q1, q2, q3)
pv, qv = sp.Matrix([p1, p2, p3]), sp.Matrix([q1, q2, q3])
dot = pv.dot(qv)
cross = pv.cross(qv)
check("pq = -p.q + p x q", qeq(qmul(P, Q), (-dot, *cross)))
V = (0, x, y, z)
check("q^2 = -(x^2+y^2+z^2)", qeq(qmul(V, V), (-(x**2 + y**2 + z**2), 0, 0, 0)))
check("q q* = x^2+y^2+z^2", qeq(qmul(V, qconj(V)), (x**2 + y**2 + z**2, 0, 0, 0)))

# 3. nabla^2 (partials commute, so treat them as commuting symbols)
dx, dy, dz = sp.symbols("dx dy dz")
N = (0, dx, dy, dz)
check("nabla^2 = -(dx^2+dy^2+dz^2)",
      qeq(qmul(N, N), (-(dx**2 + dy**2 + dz**2), 0, 0, 0)))

# 4. Gibbs/Heaviside separation
PQ, QP = qmul(P, Q), qmul(Q, P)
sym = qsmul(sp.Rational(-1, 2), qadd(PQ, QP))
anti = qsmul(sp.Rational(1, 2), qadd(PQ, neg(QP)))
check("-(pq+qp)/2 = p.q (scalar)", qeq(sym, (dot, 0, 0, 0)))
check("(pq-qp)/2 = p x q (pure)", qeq(anti, (0, *cross)))

# 5. rotation by a unit quaternion
th = sp.symbols("theta", real=True)
n1, n2 = sp.symbols("n1 n2", real=True)
n3 = sp.sqrt(1 - n1**2 - n2**2)
nv = sp.Matrix([n1, n2, n3])
q = (sp.cos(th / 2), *(sp.sin(th / 2) * nv))
check("|q| = 1", sp.simplify(sum(c**2 for c in q) - 1) == 0)
check("q* = q^{-1}", qeq(qmul(q, qconj(q)), ONE))
vv = sp.Matrix([x, y, z])
rot = qmul(qmul(q, V), qconj(q))
rod = vv * sp.cos(th) + nv.cross(vv) * sp.sin(th) + nv * nv.dot(vv) * (1 - sp.cos(th))
check("q v q* is pure", sp.simplify(rot[0]) == 0)
check("q v q* = Rodrigues rotation by theta",
      all(sp.simplify(sp.expand_trig(rot[a + 1] - rod[a])) == 0 for a in range(3)))
q2 = qsmul(2, (1, 1, 0, 0))
r2 = qmul(qmul(q2, (0, 1, 1, 0)), qconj(q2))
check("non-unit q scales by |q|^2 (negative check)",
      sum(c**2 for c in r2) == 8**2 * 2)

# 6. Pauli matrices
iC = sp.I
s1 = sp.Matrix([[0, 1], [1, 0]])
s2 = sp.Matrix([[0, -iC], [iC, 0]])
s3 = sp.Matrix([[1, 0], [0, -1]])
E = sp.eye(2)
check("sigma_a^2 = I", all(s * s == E for s in (s1, s2, s3)))
mi, mj, mk = -iC * s1, -iC * s2, -iC * s3
check("(-i sigma_a)^2 = -I", all(m * m == -E for m in (mi, mj, mk)))
check("ijk = -1 for -i sigma_a", mi * mj * mk == -E)
check("ij = k, jk = i, ki = j for -i sigma_a",
      mi * mj == mk and mj * mk == mi and mk * mi == mj)

def rho(a):
    return a[0] * E + a[1] * mi + a[2] * mj + a[3] * mk

check("rho(pq) = rho(p) rho(q)",
      sp.simplify(rho(qmul(P, Q)) - rho(P) * rho(Q)) == sp.zeros(2))

# 7. Cl_{3,0}(R) via Pauli matrices
e1, e2, e3 = s1, s2, s3

def vec(a):
    return a[0] * e1 + a[1] * e2 + a[2] * e3

pm, qm = vec(pv), vec(qv)
wedge = (pm * qm - qm * pm) / 2
check("v^2 = |v|^2 > 0", sp.simplify(vec(vv) ** 2 - (x**2 + y**2 + z**2) * E) == sp.zeros(2))
check("(pq+qp)/2 = p.q (positive scalar)",
      sp.simplify((pm * qm + qm * pm) / 2 - dot * E) == sp.zeros(2))
check("p^q is the bivector with components of p x q",
      sp.simplify(wedge - (cross[0] * e2 * e3 + cross[1] * e3 * e1 + cross[2] * e1 * e2))
      == sp.zeros(2))
B1, B2, B3 = e2 * e3, e3 * e1, e1 * e2
check("bivectors square to -1", all(b * b == -E for b in (B1, B2, B3)))
check("(-e2e3)(-e3e1)(-e1e2) = -1", (-B1) * (-B2) * (-B3) == -E)
check("-e2e3 = -i sigma_1 etc.", -B1 == mi and -B2 == mj and -B3 == mk)
R = sp.cos(th / 2) * E + sp.sin(th / 2) * (n1 * (-B1) + n2 * (-B2) + n3 * (-B3))
Rrev = sp.cos(th / 2) * E - sp.sin(th / 2) * (n1 * (-B1) + n2 * (-B2) + n3 * (-B3))
check("R R~ = 1", sp.simplify(R * Rrev - E) == sp.zeros(2))
check("R v R~ = vector of q v q*",
      sp.simplify(sp.expand_trig(R * vec(vv) * Rrev - vec(rod))) == sp.zeros(2))

# 8. Cl_{0,2}(R) = H with grade-1 generators i, j
check("Cl_{0,2}: i^2 = j^2 = -1, ij = -ji, (ij)^2 = -1, ij = k",
      mi * mi == -E and mj * mj == -E and mi * mj == -mj * mi
      and (mi * mj) ** 2 == -E and mi * mj == mk)

# 9. 4D: x -> a x b^{-1}
random.seed(0)

def unit():
    v = [random.gauss(0, 1) for _ in range(4)]
    r = sum(c * c for c in v) ** 0.5
    return tuple(c / r for c in v)

a, b = unit(), unit()
M = sp.Matrix([qmul(qmul(a, u), qconj(b)) for u in (ONE, I, J, K)]).T
check("x -> a x b^{-1} is orthogonal", (M.T * M - sp.eye(4)).norm() < 1e-12)
check("x -> a x b^{-1} has det +1", abs(M.det() - 1) < 1e-12)
check("a = b fixes 1", all(abs(c - d) < 1e-12
                            for c, d in zip(qmul(qmul(a, ONE), qconj(a)), ONE)))

# 10. planes vs vectors
check("n(n-1)/2 = n only for n = 3 (1 <= n <= 20)",
      [n for n in range(1, 21) if n * (n - 1) // 2 == n] == [3])

# 11. pseudovector law
def rand_orth(det):
    A = sp.Matrix(3, 3, lambda r, c: random.gauss(0, 1))
    Qm, _ = A.QRdecomposition()
    if (Qm.det() > 0) != (det > 0):
        Qm[:, 0] = -Qm[:, 0]
    return Qm

pn = sp.Matrix([random.gauss(0, 1) for _ in range(3)])
qn = sp.Matrix([random.gauss(0, 1) for _ in range(3)])
for d in (1, -1):
    Qm = rand_orth(d)
    lhs = (Qm * pn).cross(Qm * qn)
    rhs = Qm.det() * Qm * pn.cross(qn)
    check(f"(Qp) x (Qq) = det(Q) Q (p x q), det Q = {d:+d}",
          (lhs - rhs).norm() < 1e-9 and abs(Qm.det() - d) < 1e-9)
