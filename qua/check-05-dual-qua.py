"""Checks for 05-dual-qua.md.

Verifies dual quaternions H (x) D and their rigid-motion representation:

1. Dual numbers: eps^2 = 0, zero divisors eps a * eps b = 0, the inverse
   (x + eps y)^-1 = 1/x - eps y/x^2, exp(eps y) = 1 + eps y exactly, and
   f(x + eps y) = f(x) + eps y f'(x) for a polynomial and for sin/exp
   (evaluated by running the Taylor series in dual arithmetic).
2. Dual quaternion algebra: eps is central; the three conjugations satisfy
   bar(st) = bar(t) bar(s) (anti), (st)* = s* t* (homo), and
   bar((st)*) = bar(t*) bar(s*) (anti).
3. The norm sigma bar(sigma) is dual-number valued (no i, j, k parts) and
   multiplicative; sigma bar(sigma) = 1 iff |p| = 1 and <p, q> = 0.
4. Rigid motion: sigma = r + (eps/2) t r sends 1 + eps x to
   1 + eps (r x bar(r) + t) under the sandwich with bar(sigma*); the
   sandwich with plain bar(sigma) drops the translation entirely
   (negative check); composition sigma2 sigma1 realizes the SE(3) product
   (r2 r1, t2 + r2 t1 bar(r2)); +-sigma give the same action; and
   sigma bar(sigma) = 1.
5. Master equation with formal eps: for a planar triangle (angle sum pi,
   sides from the law of sines) the six-factor product
   e^{i(pi-alpha)/2} (1 + k eps b/2) e^{i(pi-gamma)/2} (1 + k eps a/2)
   e^{i(pi-beta)/2} (1 + k eps c/2) equals -1 exactly; perturbing one side
   or one angle breaks it (negative checks); and tracing a point 1 + eps x
   through the six motions returns it to the start.
6. Screw motion: exp((theta + eps d)/2 (l + eps m)) matches both the
   closed form cos(theta^/2) + sin(theta^/2) L (expanded by the f(x+eps y)
   rule) and the series exp; L^2 = -1; it equals the composition
   T(d l) T(x0) R T(-x0) for m = x0 x l (rotation about the displaced
   axis, then slide); points on the axis translate by d l; extraction of
   (theta, d, l, m) from a random unit sigma reconstructs it; theta = 0
   degenerates to a pure translation.
7. Cl+_{3,0,1}: with e0^2 = 0, the correspondence 1 <-> 1, i <-> e3e2,
   j <-> e1e3, k <-> e2e1, eps <-> e0e1e2e3, eps i <-> e0e1,
   eps j <-> e0e2, eps k <-> e0e3 is an algebra isomorphism (all 64 basis
   products plus random elements); e0123 commutes with all six bivectors,
   anticommutes with all four vectors, and squares to 0.
8. Pseudoscalar trichotomy: omega^2 = +1 in Cl_{4,0}, -1 in Cl_{3,1},
   0 in Cl_{3,0,1}.
"""

import math
import random

random.seed(0)

# --- dual numbers: pairs (x, y) = x + eps y ----------------------------------

def nmul(a, b):
    return (a[0] * b[0], a[0] * b[1] + a[1] * b[0])

def neq(a, b, tol=1e-9):
    return all(abs(s - t) < tol for s, t in zip(a, b))

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

def qdot(a, b):
    return sum(s * t for s, t in zip(a, b))

def qeq(a, b, tol=1e-9):
    return all(abs(s - t) < tol for s, t in zip(a, b))

Q1 = (1, 0, 0, 0)
Q0 = (0, 0, 0, 0)

def qexp_pure(v):  # exp of a pure quaternion
    n = math.sqrt(qdot(v, v))
    if n < 1e-300:
        return Q1
    return qadd((math.cos(n), 0, 0, 0), qsmul(math.sin(n) / n, v))

def rand_vec():  # random pure quaternion
    return (0,) + tuple(random.gauss(0, 1) for _ in range(3))

def rand_rotor():  # random unit quaternion
    v = tuple(random.gauss(0, 1) for _ in range(4))
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v)

def cross(a, b):  # cross product of pure quaternions
    return (0,
            a[2] * b[3] - a[3] * b[2],
            a[3] * b[1] - a[1] * b[3],
            a[1] * b[2] - a[2] * b[1])

# --- dual quaternions: pairs (p, q) = p + eps q ------------------------------

def dmul(a, b):
    return (qmul(a[0], b[0]), qadd(qmul(a[0], b[1]), qmul(a[1], b[0])))

def dadd(a, b):
    return (qadd(a[0], b[0]), qadd(a[1], b[1]))

def dsmul(c, a):  # scalar c: dual number pair or real
    if not isinstance(c, tuple):
        c = (c, 0)
    return (qsmul(c[0], a[0]), qadd(qsmul(c[0], a[1]), qsmul(c[1], a[0])))

def dbar(a):  # quaternion conjugate
    return (qconj(a[0]), qconj(a[1]))

def dstar(a):  # dual-number conjugate
    return (a[0], qsmul(-1, a[1]))

def deq(a, b, tol=1e-9):
    return qeq(a[0], b[0], tol) and qeq(a[1], b[1], tol)

D1 = (Q1, Q0)

def rand_dq():
    return (tuple(random.gauss(0, 1) for _ in range(4)),
            tuple(random.gauss(0, 1) for _ in range(4)))

def motion(r, t):  # rotate by rotor r, then translate by t (world frame)
    return (r, qsmul(0.5, qmul(t, r)))

def point(x):  # embed a point as 1 + eps x
    return (Q1, x)

def act(s, x):  # sandwich sigma (1 + eps x) bar(sigma*)
    return dmul(dmul(s, point(x)), dbar(dstar(s)))

# --- check 1: dual numbers ---------------------------------------------------

EPS = (0, 1)
assert neq(nmul(EPS, EPS), (0, 0))
for _ in range(20):
    a, b = random.gauss(0, 1), random.gauss(0, 1)
    assert neq(nmul((0, a), (0, b)), (0, 0))  # zero divisors
    x, y = random.gauss(0, 1) + 2, random.gauss(0, 1)
    inv = (1 / x, -y / x**2)
    assert neq(nmul((x, y), inv), (1, 0))

def nexp(z, terms=30):  # exp by its Taylor series in dual arithmetic
    s, t = (1, 0), (1, 0)
    for n in range(1, terms):
        t = nmul(t, (z[0] / n, z[1] / n))
        s = (s[0] + t[0], s[1] + t[1])
    return s

def nsin(z, terms=20):  # sin by its Taylor series in dual arithmetic
    s, t = z, z
    z2 = nmul(z, z)
    for n in range(1, terms):
        t = dsmul_n(-1 / ((2 * n) * (2 * n + 1)), nmul(t, z2))
        s = (s[0] + t[0], s[1] + t[1])
    return s

def dsmul_n(c, z):
    return (c * z[0], c * z[1])

for _ in range(20):
    y = random.gauss(0, 1)
    assert neq(nexp((0, y)), (1, y))  # exp(eps y) = 1 + eps y exactly
    x = random.gauss(0, 1)
    assert neq(nexp((x, y)), (math.exp(x), y * math.exp(x)))
    assert neq(nsin((x, y)), (math.sin(x), y * math.cos(x)))
    # polynomial p(t) = t^3 - 2t + 5, p'(t) = 3t^2 - 2
    z = (x, y)
    p = nmul(nmul(z, z), z)
    p = (p[0] - 2 * x + 5, p[1] - 2 * y)
    assert neq(p, (x**3 - 2 * x + 5, y * (3 * x**2 - 2)))

# --- check 2: centrality of eps and the three conjugations -------------------

DEPS = (Q0, Q1)
assert deq(dmul(DEPS, DEPS), (Q0, Q0))
for _ in range(50):
    s, t = rand_dq(), rand_dq()
    assert deq(dmul(DEPS, s), dmul(s, DEPS))
    st = dmul(s, t)
    assert deq(dbar(st), dmul(dbar(t), dbar(s)))
    assert deq(dstar(st), dmul(dstar(s), dstar(t)))
    assert deq(dbar(dstar(st)), dmul(dbar(dstar(t)), dbar(dstar(s))))
    assert deq(dbar(dstar(s)), dstar(dbar(s)))

# --- check 3: dual-valued norm and the unit condition ------------------------

def is_dual_scalar(s, tol=1e-9):
    return all(abs(c) < tol for c in s[0][1:]) and \
           all(abs(c) < tol for c in s[1][1:])

for _ in range(50):
    s, t = rand_dq(), rand_dq()
    ns, nt = dmul(s, dbar(s)), dmul(t, dbar(t))
    assert is_dual_scalar(ns)
    p, q = s
    assert abs(ns[0][0] - qdot(p, p)) < 1e-9        # real part |p|^2
    assert abs(ns[1][0] - 2 * qdot(p, q)) < 1e-9    # dual part 2<p,q>
    nst = dmul(dmul(s, t), dbar(dmul(s, t)))
    assert neq((nst[0][0], nst[1][0]),
               nmul((ns[0][0], ns[1][0]), (nt[0][0], nt[1][0])))
    # unit <=> |p| = 1 and <p,q> = 0
    r = rand_rotor()
    q0 = tuple(random.gauss(0, 1) for _ in range(4))
    q = qadd(q0, qsmul(-qdot(r, q0), r))
    assert deq(dmul((r, q), dbar((r, q))), D1)

# --- check 4: rigid motions and the double cover of SE(3) --------------------

for _ in range(50):
    r = rand_rotor()
    t, x = rand_vec(), rand_vec()
    s = motion(r, t)
    assert deq(dmul(s, dbar(s)), D1)  # unit
    moved = qadd(qmul(qmul(r, x), qconj(r)), t)
    assert deq(act(s, x), point(moved))
    # plain-bar sandwich drops the translation
    plain = dmul(dmul(s, point(x)), dbar(s))
    assert deq(plain, point(qmul(qmul(r, x), qconj(r))))
    # -sigma gives the same action
    assert deq(act(dsmul(-1, s), x), point(moved))
    # composition = SE(3) semidirect product
    r2, t2 = rand_rotor(), rand_vec()
    s2 = motion(r2, t2)
    comp = motion(qmul(r2, r), qadd(t2, qmul(qmul(r2, t), qconj(r2))))
    assert deq(dmul(s2, s), comp)

# --- check 5: the master equation holds exactly in the plane -----------------

def turn(phi):  # e^{i phi/2} as a dual quaternion
    return (qexp_pure((0, phi / 2, 0, 0)), Q0)

def edge(x):  # 1 + k eps x/2: translation by the vector k x
    return ((1, 0, 0, 0), (0, 0, 0, x / 2))

def master(alpha, beta, gamma, a, b, c):
    m = D1
    for f in (turn(math.pi - alpha), edge(b), turn(math.pi - gamma),
              edge(a), turn(math.pi - beta), edge(c)):
        m = dmul(m, f)
    return m

MINUS1 = (qsmul(-1, Q1), Q0)
for _ in range(50):
    alpha = random.uniform(0.2, 2.0)
    beta = random.uniform(0.2, min(2.0, math.pi - alpha - 0.2))
    gamma = math.pi - alpha - beta
    scale = random.uniform(0.5, 3.0)  # law of sines: a/sin(alpha) = const
    a, b, c = (scale * math.sin(x) for x in (alpha, beta, gamma))
    assert deq(master(alpha, beta, gamma, a, b, c), MINUS1)
    # negative checks: perturb one side, or one angle
    assert not deq(master(alpha, beta, gamma, a, b, c + 0.1), MINUS1)
    assert not deq(master(alpha + 0.1, beta, gamma, a, b, c), MINUS1)
    # tracing a point through the six motions returns it to the start
    x = rand_vec()
    y = x
    for f in (edge(c), turn(math.pi - beta), edge(a), turn(math.pi - gamma),
              edge(b), turn(math.pi - alpha)):  # applied right to left
        y = act(f, y)[1]
    assert qeq(y, x)

# --- check 6: screw motions --------------------------------------------------

def dexp(s, terms=40):  # series exp of a dual quaternion
    out, t = D1, D1
    for n in range(1, terms):
        t = dsmul(1 / n, dmul(t, s))
        out = dadd(out, t)
    return out

def screw(theta, d, l, m):  # cos(theta^/2) + sin(theta^/2) (l + eps m)
    ch = (math.cos(theta / 2), -d / 2 * math.sin(theta / 2))
    sh = (math.sin(theta / 2), d / 2 * math.cos(theta / 2))
    return dadd(dsmul(ch, D1), dsmul(sh, (l, m)))

for _ in range(50):
    l0 = rand_vec()
    l = qsmul(1 / math.sqrt(qdot(l0, l0)), l0)
    x0 = rand_vec()
    m = cross(x0, l)
    L = (l, m)
    assert deq(dmul(L, L), dsmul(-1, D1))  # L^2 = -1
    theta = random.uniform(0.3, 3.0)
    d = random.gauss(0, 1)
    # exp argument: (theta + eps d)/2 * L
    arg = dsmul((theta / 2, d / 2), L)
    s = screw(theta, d, l, m)
    assert deq(dexp(arg), s)
    assert deq(dmul(s, dbar(s)), D1)
    # equals rotation about the displaced axis, then slide: T(dl) T(x0) R T(-x0)
    R = (qexp_pure(qsmul(theta / 2, l)), Q0)
    comp = dmul(motion(Q1, qsmul(d, l)),
                dmul(motion(Q1, x0), dmul(R, motion(Q1, qsmul(-1, x0)))))
    assert deq(s, comp)
    # points on the axis translate by d l
    axis_pt = qadd(x0, qsmul(random.gauss(0, 1), l))
    assert deq(act(s, axis_pt), point(qadd(axis_pt, qsmul(d, l))))
    # extraction from a random unit sigma and reconstruction
    r = rand_rotor()
    if r[0] < 0:
        r = qsmul(-1, r)  # fix the sign so theta lies in (0, pi)
    q0 = tuple(random.gauss(0, 1) for _ in range(4))
    q = qadd(q0, qsmul(-qdot(r, q0), r))
    vn = math.sqrt(qdot(r, r) - r[0] ** 2)
    theta2 = 2 * math.atan2(vn, r[0])
    l2 = qsmul(1 / vn, (0,) + r[1:])
    sh, ch = math.sin(theta2 / 2), math.cos(theta2 / 2)
    d2 = -2 * q[0] / sh
    m2 = qsmul(1 / sh, qadd((0,) + q[1:], qsmul(-d2 / 2 * ch, l2)))
    assert abs(qdot(l2, m2)) < 1e-9  # the axis is a genuine line
    assert deq(screw(theta2, d2, l2, m2), (r, q))
    # theta = 0: pure translation
    assert deq(screw(0, d, l, Q0), motion(Q1, qsmul(d, l)))

# --- Clifford algebra: multivector = {sorted index tuple: coeff} -------------

def bmul(a, b, sq=None):
    """Product of two basis elements (sorted tuples) -> (sign, basis).

    sq maps index -> square of that generator; None means all +1.
    """
    s, t = 1, list(a)
    for x in b:
        i = len(t)
        while i > 0 and t[i - 1] > x:
            i -= 1
            s = -s
        if i > 0 and t[i - 1] == x:
            t.pop(i - 1)
            if sq is not None:
                s *= sq[x]
        else:
            t.insert(i, x)
    return s, tuple(t)

def mul(x, y, sq=None):
    z = {}
    for a, ca in x.items():
        for b, cb in y.items():
            s, t = bmul(a, b, sq)
            z[t] = z.get(t, 0) + s * ca * cb
    return {k: v for k, v in z.items() if abs(v) > 1e-12}

def add(x, y):
    z = dict(x)
    for k, v in y.items():
        z[k] = z.get(k, 0) + v
    return {k: v for k, v in z.items() if abs(v) > 1e-12}

def smul(c, x):
    return {k: c * v for k, v in x.items()}

def eq(x, y, tol=1e-9):
    z = dict(x)
    for k, v in y.items():
        z[k] = z.get(k, 0) - v
    return all(abs(v) < tol for v in z.values())

# --- check 7: dual quaternions = Cl+_{3,0,1} ---------------------------------

SQ301 = {0: 0, 1: 1, 2: 1, 3: 1}

# basis map: (real part 1,i,j,k | dual part 1,i,j,k) -> even basis of Cl_{3,0,1}
CL_BASIS = [{(): 1}, {(2, 3): -1}, {(1, 3): 1}, {(1, 2): -1},          # 1,i,j,k
            {(0, 1, 2, 3): 1}, {(0, 1): 1}, {(0, 2): 1}, {(0, 3): 1}]  # eps*...

def dq_to_cl(s):
    out = {}
    for c, b in zip(s[0] + s[1], CL_BASIS):
        out = add(out, smul(c, b))
    return out

# all 64 products of the 8 basis elements
DQ_BASIS = [(tuple(1 if i == n else 0 for i in range(4)), Q0) for n in range(4)] + \
           [(Q0, tuple(1 if i == n else 0 for i in range(4))) for n in range(4)]
for u in DQ_BASIS:
    for v in DQ_BASIS:
        assert eq(dq_to_cl(dmul(u, v)),
                  mul(dq_to_cl(u), dq_to_cl(v), SQ301))
for _ in range(50):  # and random elements
    s, t = rand_dq(), rand_dq()
    assert eq(dq_to_cl(dmul(s, t)), mul(dq_to_cl(s), dq_to_cl(t), SQ301))

E = lambda *idx: {tuple(idx): 1}
E0123 = E(0, 1, 2, 3)
assert eq(mul(E(0), E(0), SQ301), {})  # e0^2 = 0
assert eq(mul(E0123, E0123, SQ301), {})  # nilpotent pseudoscalar
for i in range(4):
    for j in range(i + 1, 4):  # commutes with all six bivectors
        assert eq(mul(E0123, E(i, j), SQ301), mul(E(i, j), E0123, SQ301))
    # anticommutes with all four vectors (centrality is even-only)
    assert eq(mul(E0123, E(i), SQ301), smul(-1, mul(E(i), E0123, SQ301)))

# --- check 8: pseudoscalar trichotomy ----------------------------------------

W4 = E(1, 2, 3, 4)
assert eq(mul(W4, W4, {1: 1, 2: 1, 3: 1, 4: 1}), {(): 1})    # Cl_{4,0}: +1
assert eq(mul(W4, W4, {1: 1, 2: 1, 3: 1, 4: -1}), {(): -1})  # Cl_{3,1}: -1
assert eq(mul(E0123, E0123, SQ301), {})                      # Cl_{3,0,1}: 0

print("All checks passed.")
