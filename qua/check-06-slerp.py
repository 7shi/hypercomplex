"""Checks for 06-slerp.md.

Verifies SLERP (spherical linear interpolation) of rotors:

1. Log and power: log(exp v) = v for |v| < pi, exp(log q) = q, the
   principal branch |log q| in [0, pi], and the one-direction exponent
   laws q^s q^t = q^{s+t}, q^0 = 1, q^1 = q, (q^t)^{-1} = q^{-t}.
2. Angle: cos Omega = <q0, q1> = Re(q0^{-1} q1); invariance under common
   left and right multiplication; Omega = theta/2 (the sandwich by the
   relative rotor turns a vector orthogonal to the axis by 2 Omega).
3. Geodesic: the slerp path has constant speed and total length Omega;
   a path through a random detour point is at least as long, and strictly
   longer when the detour leaves the great circle (minimality).
4. The power form q0 (q0^{-1} q1)^t and the sin-weight form agree for t
   in [-0.5, 1.5] (including extrapolation); endpoints match; for tiny
   Omega both agree with the plain lerp.
5. Invariance slerp(r q0, r q1; t) = r slerp(q0, q1; t), the right-hand
   version, and the symmetry slerp(q1, q0; 1-t) = slerp(q0, q1; t).
6. Uniform motion: q0^{-1} q(t) = exp(t theta n / 2); the world-frame
   axis q(t) n q(t)^{-1} is constant; finite differences give constant
   angular velocities 2 q' q^{-1} = theta q0 n q0^{-1} and
   2 q^{-1} q' = theta n.
7. Sign selection: Omega(q0, -q1) = pi - Omega(q0, q1); if <q0, q1> < 0
   the -q1 side is shorter; the midpoint rotations of the two paths
   differ (negative check); slerp(-q0, -q1; t) = -slerp(q0, q1; t).
8. NLERP: it lies on the same great arc (a reparameterized slerp); its
   angle satisfies tan psi(t) = t sin Omega / (1 - t + t cos Omega); for
   Omega = pi/2 the midpoint speed is twice the endpoint speed (not
   uniform, negative check).
9. General rotors: slerp in Cl+_{3,0} matches the quaternion slerp under
   the standard isomorphism; in Cl_{4,0} the power of a double-rotation
   rotor multiplies both plane angles by t (commuting factors, rotor
   condition u u~ = 1, sandwich turns e1 by t alpha); the linear
   combination of 1 and a genuine double-rotation rotor leaves a nonzero
   e1234 part in u u~, so no real rescaling makes it a rotor (negative
   check), while for a simple rotation (beta = 0) normalization does
   give a rotor (positive control).
10. ScLERP: sigma0 (sigma0^{-1} sigma1)^t = sigma0 screw(t theta, t d,
    l, m) interpolates rigid motions: endpoints match up to sign, the
    real part is the slerp of the real parts, and points on the screw
    axis translate linearly by t d l.
11. Degeneracies: q1 near q0 (tiny Omega) stays stable and matches the
    lerp; <q0, q1> = 0 gives two equal-length shortest paths.
"""

import math
import random

random.seed(0)

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

def qsub(a, b):
    return tuple(s - t for s, t in zip(a, b))

def qsmul(c, a):
    return tuple(c * s for s in a)

def qdot(a, b):
    return sum(s * t for s, t in zip(a, b))

def qnorm(a):
    return math.sqrt(qdot(a, a))

def qeq(a, b, tol=1e-9):
    return all(abs(s - t) < tol for s, t in zip(a, b))

Q1 = (1, 0, 0, 0)
Q0 = (0, 0, 0, 0)

def qexp_pure(v):  # exp of a pure quaternion
    n = qnorm(v)
    if n < 1e-300:
        return Q1
    return qadd((math.cos(n), 0, 0, 0), qsmul(math.sin(n) / n, v))

def qlog(q):  # principal log of a unit quaternion (axis undefined at q = -1)
    im = (0,) + q[1:]
    n = qnorm(im)
    if n < 1e-300:
        return Q0
    return qsmul(math.atan2(n, q[0]) / n, im)

def qpow(q, t):
    return qexp_pure(qsmul(t, qlog(q)))

def qang(a, b):  # angle Omega between unit quaternions
    return math.acos(max(-1.0, min(1.0, qdot(a, b))))

def slerp_pow(q0, q1, t):  # q0 (q0^{-1} q1)^t
    return qmul(q0, qpow(qmul(qconj(q0), q1), t))

def slerp_sin(q0, q1, t):  # sin-weight form (0 < Omega < pi)
    o = qang(q0, q1)
    return qsmul(1 / math.sin(o),
                 qadd(qsmul(math.sin((1 - t) * o), q0),
                      qsmul(math.sin(t * o), q1)))

def lerp(q0, q1, t):
    return qadd(qsmul(1 - t, q0), qsmul(t, q1))

def nlerp(q0, q1, t):
    v = lerp(q0, q1, t)
    return qsmul(1 / qnorm(v), v)

def rand_vec():  # random pure quaternion
    return (0,) + tuple(random.gauss(0, 1) for _ in range(3))

def rand_unit_pure():
    v = rand_vec()
    return qsmul(1 / qnorm(v), v)

def rand_rotor():  # random unit quaternion
    v = tuple(random.gauss(0, 1) for _ in range(4))
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v)

def cross(a, b):  # cross product of pure quaternions
    return (0,
            a[2] * b[3] - a[3] * b[2],
            a[3] * b[1] - a[1] * b[3],
            a[1] * b[2] - a[2] * b[1])

# --- check 1: log, exp, and powers -------------------------------------------

for _ in range(50):
    v = qsmul(random.uniform(0.01, math.pi - 0.01), rand_unit_pure())
    assert qeq(qlog(qexp_pure(v)), v)
    q = rand_rotor()
    assert qeq(qexp_pure(qlog(q)), q)
    assert qnorm(qlog(q)) <= math.pi + 1e-12  # principal branch
    s, t = random.uniform(-2, 2), random.uniform(-2, 2)
    assert qeq(qmul(qpow(q, s), qpow(q, t)), qpow(q, s + t))
    assert qeq(qpow(q, 0), Q1)
    assert qeq(qpow(q, 1), q)
    assert qeq(qpow(q, -t), qconj(qpow(q, t)))

# --- check 2: the angle Omega and the half-angle relation --------------------

for _ in range(50):
    q0, q1, r = rand_rotor(), rand_rotor(), rand_rotor()
    o = qang(q0, q1)
    rel = qmul(qconj(q0), q1)
    assert abs(math.cos(o) - rel[0]) < 1e-9  # cos Omega = Re(q0^{-1} q1)
    assert abs(qang(qmul(r, q0), qmul(r, q1)) - o) < 1e-9
    assert abs(qang(qmul(q0, r), qmul(q1, r)) - o) < 1e-9
    # theta = 2 Omega: the sandwich turns a vector orthogonal to the axis
    lg = qlog(rel)
    phi = qnorm(lg)
    if 0.05 < phi < math.pi / 2 - 0.05:  # keep 2 Omega below pi
        n = qsmul(1 / phi, lg)
        v0 = rand_vec()
        v = qsub(v0, qsmul(qdot(v0, n), n))
        v = qsmul(1 / qnorm(v), v)
        w = qmul(qmul(rel, v), qconj(rel))
        assert abs(qang(v, w) - 2 * o) < 1e-9

# --- check 3: constant speed and minimality of the great arc -----------------

for _ in range(10):
    q0, q1 = rand_rotor(), rand_rotor()
    o = qang(q0, q1)
    if o < 0.1 or o > math.pi - 0.1:
        continue
    N = 1000
    pts = [slerp_pow(q0, q1, i / N) for i in range(N + 1)]
    segs = [qnorm(qsub(pts[i + 1], pts[i])) for i in range(N)]
    assert abs(sum(segs) - o) < 1e-4          # total length = Omega
    assert max(segs) - min(segs) < 1e-9       # constant speed
    # detour through a random point is never shorter
    w = rand_rotor()
    detour = qang(q0, w) + qang(w, q1)
    assert detour >= o - 1e-9
    # strictly longer when w leaves the great circle
    e1 = q0
    e2 = qsub(q1, qsmul(qdot(q1, e1), e1))
    e2 = qsmul(1 / qnorm(e2), e2)
    perp = qsub(w, qadd(qsmul(qdot(w, e1), e1), qsmul(qdot(w, e2), e2)))
    if qnorm(perp) > 0.3:
        assert detour > o + 1e-6

# --- check 4: the two forms agree; endpoints; the lerp limit -----------------

for _ in range(50):
    q0, q1 = rand_rotor(), rand_rotor()
    o = qang(q0, q1)
    if o < 0.05 or o > math.pi - 0.05:
        continue
    for t in (-0.5, -0.1, 0, 0.3, 0.5, 0.7, 1, 1.5):
        assert qeq(slerp_pow(q0, q1, t), slerp_sin(q0, q1, t), 1e-8)
    assert qeq(slerp_pow(q0, q1, 0), q0)
    assert qeq(slerp_pow(q0, q1, 1), q1)
for _ in range(20):
    q0 = rand_rotor()
    q1 = qmul(q0, qexp_pure(qsmul(1e-4, rand_unit_pure())))
    for t in (0.25, 0.5, 0.75):
        assert qeq(slerp_pow(q0, q1, t), lerp(q0, q1, t), 1e-6)

# --- check 5: invariance and symmetry ----------------------------------------

for _ in range(50):
    q0, q1, r = rand_rotor(), rand_rotor(), rand_rotor()
    if qdot(q0, q1) < -0.99:
        continue
    t = random.uniform(0, 1)
    assert qeq(slerp_pow(qmul(r, q0), qmul(r, q1), t),
               qmul(r, slerp_pow(q0, q1, t)))
    assert qeq(slerp_pow(qmul(q0, r), qmul(q1, r), t),
               qmul(slerp_pow(q0, q1, t), r))
    assert qeq(slerp_pow(q1, q0, 1 - t), slerp_pow(q0, q1, t))

# --- check 6: fixed world axis and constant angular velocity -----------------

for _ in range(20):
    q0, q1 = rand_rotor(), rand_rotor()
    rel = qmul(qconj(q0), q1)
    lg = qlog(rel)
    phi = qnorm(lg)
    if phi < 0.1 or phi > math.pi - 0.1:
        continue
    n = qsmul(1 / phi, lg)
    theta = 2 * phi
    world_axis = qmul(qmul(q0, n), qconj(q0))
    h = 1e-6
    for t in (0.2, 0.5, 0.8):
        qt = slerp_pow(q0, q1, t)
        assert qeq(qmul(qconj(q0), qt), qexp_pure(qsmul(t * phi, n)))
        assert qeq(qmul(qmul(qt, n), qconj(qt)), world_axis)
        qd = qsmul(1 / (2 * h), qsub(slerp_pow(q0, q1, t + h),
                                     slerp_pow(q0, q1, t - h)))
        assert qeq(qsmul(2, qmul(qd, qconj(qt))),
                   qsmul(theta, world_axis), 1e-4)   # 2 q' q^{-1}
        assert qeq(qsmul(2, qmul(qconj(qt), qd)),
                   qsmul(theta, n), 1e-4)            # 2 q^{-1} q'

# --- check 7: the two lifts and the sign selection ---------------------------

for _ in range(50):
    q0, q1 = rand_rotor(), rand_rotor()
    o = qang(q0, q1)
    assert abs(qang(q0, qsmul(-1, q1)) - (math.pi - o)) < 1e-9
    if qdot(q0, q1) < -0.05:
        assert qang(q0, qsmul(-1, q1)) < o  # the flipped side is shorter
    if abs(qdot(q0, q1)) < 0.999:
        # the midpoints of the two paths are different rotations
        m1 = slerp_pow(q0, q1, 0.5)
        m2 = slerp_pow(q0, qsmul(-1, q1), 0.5)
        assert not (qeq(m1, m2, 1e-6) or qeq(m1, qsmul(-1, m2), 1e-6))
    t = random.uniform(0, 1)
    assert qeq(slerp_pow(qsmul(-1, q0), qsmul(-1, q1), t),
               qsmul(-1, slerp_pow(q0, q1, t)))

# --- check 8: NLERP shares the arc but not the speed -------------------------

for _ in range(50):
    q0, q1 = rand_rotor(), rand_rotor()
    o = qang(q0, q1)
    if o < 0.1 or o > math.pi - 0.1:
        continue
    t = random.uniform(0.05, 0.95)
    w = nlerp(q0, q1, t)
    psi = qang(q0, w)
    assert qeq(slerp_pow(q0, q1, psi / o), w, 1e-8)  # same great arc
    assert abs(psi - math.atan2(t * math.sin(o),
                                1 - t + t * math.cos(o))) < 1e-9

q0 = rand_rotor()
q1 = qmul(q0, qexp_pure(qsmul(math.pi / 2, rand_unit_pure())))  # Omega = pi/2
h = 1e-6
psi = lambda t: qang(q0, nlerp(q0, q1, t))
d0 = (psi(2 * h) - psi(0)) / (2 * h)
dm = (psi(0.5 + h) - psi(0.5 - h)) / (2 * h)
assert abs(dm / d0 - 2) < 1e-3  # twice as fast in the middle
assert abs(dm - d0) > 0.5       # negative check: not uniform

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

def rev(x):  # reversion: grade k gets (-1)^{k(k-1)/2}
    return {k: (v if len(k) % 4 in (0, 1) else -v) for k, v in x.items()}

def cl_exp(x, sq, terms=40):
    out, t = {(): 1}, {(): 1}
    for i in range(1, terms):
        t = smul(1 / i, mul(t, x, sq))
        out = add(out, t)
    return out

E = lambda *idx: {tuple(idx): 1}
CL1 = E()

# --- check 9: slerp for general rotors ---------------------------------------

SQ3 = {1: 1, 2: 1, 3: 1}
CL_BASIS3 = [{(): 1}, {(2, 3): -1}, {(1, 3): 1}, {(1, 2): -1}]  # 1, i, j, k

def q_to_cl(q):
    out = {}
    for c, b in zip(q, CL_BASIS3):
        out = add(out, smul(c, b))
    return out

def cl_log3(x):  # log of a rotor in Cl+_{3,0}: scalar + bivector
    biv = {k: v for k, v in x.items() if len(k) == 2}
    bn = math.sqrt(sum(v * v for v in biv.values()))
    if bn < 1e-300:
        return {}
    return smul(math.atan2(bn, x.get((), 0)) / bn, biv)

def cl_slerp3(r0, r1, t):
    u = mul(rev(r0), r1, SQ3)
    return mul(r0, cl_exp(smul(t, cl_log3(u)), SQ3), SQ3)

for _ in range(20):
    a, b = rand_rotor(), rand_rotor()
    assert eq(q_to_cl(qmul(a, b)), mul(q_to_cl(a), q_to_cl(b), SQ3))  # iso
    if qdot(a, b) < -0.99:
        continue
    t = random.uniform(0, 1)
    assert eq(q_to_cl(slerp_pow(a, b, t)), cl_slerp3(q_to_cl(a), q_to_cl(b), t))

SQ4 = {1: 1, 2: 1, 3: 1, 4: 1}

for _ in range(20):
    al, be = random.uniform(0.3, 2.5), random.uniform(0.3, 2.5)
    t = random.uniform(0.1, 0.9)
    B = add(smul(al / 2, E(1, 2)), smul(be / 2, E(3, 4)))
    u = cl_exp(B, SQ4)
    ut = cl_exp(smul(t, B), SQ4)
    # commuting planes: exp splits into the product of the two factors
    assert eq(ut, mul(cl_exp(smul(t * al / 2, E(1, 2)), SQ4),
                      cl_exp(smul(t * be / 2, E(3, 4)), SQ4), SQ4))
    assert eq(mul(ut, rev(ut), SQ4), CL1)  # rotor condition
    assert eq(cl_exp(smul(1, B), SQ4), u)  # endpoint
    # the sandwich turns e1 by t*alpha in the 12-plane
    v = mul(mul(rev(ut), E(1), SQ4), ut, SQ4)
    assert eq(v, {(1,): math.cos(t * al), (2,): math.sin(t * al)})
    # negative check: 1 + u cannot be rescaled into a rotor
    w = add(CL1, u)
    ww = mul(w, rev(w), SQ4)
    assert abs(ww.get((1, 2, 3, 4), 0)) > 1e-3  # nonzero e1234 part
    wn = smul(1 / math.sqrt(ww[()]), w)
    assert not eq(mul(wn, rev(wn), SQ4), CL1)
    # positive control: a simple rotation (beta = 0) does normalize
    u1 = cl_exp(smul(al / 2, E(1, 2)), SQ4)
    w1 = add(CL1, u1)
    wn1 = smul(1 / math.sqrt(mul(w1, rev(w1), SQ4)[()]), w1)
    assert eq(mul(wn1, rev(wn1), SQ4), CL1)

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

def motion(r, t):  # rotate by rotor r, then translate by t (world frame)
    return (r, qsmul(0.5, qmul(t, r)))

def point(x):  # embed a point as 1 + eps x
    return (Q1, x)

def act(s, x):  # sandwich sigma (1 + eps x) bar(sigma*)
    return dmul(dmul(s, point(x)), dbar(dstar(s)))

def screw(theta, d, l, m):  # cos(theta^/2) + sin(theta^/2) (l + eps m)
    ch = (math.cos(theta / 2), -d / 2 * math.sin(theta / 2))
    sh = (math.sin(theta / 2), d / 2 * math.cos(theta / 2))
    return dadd(dsmul(ch, D1), dsmul(sh, (l, m)))

# --- check 10: ScLERP --------------------------------------------------------

for _ in range(20):
    s0 = motion(rand_rotor(), rand_vec())
    s1 = motion(rand_rotor(), rand_vec())
    relq = dmul(dbar(s0), s1)  # unit, so the inverse is dbar
    if relq[0][0] < 0:
        relq = dsmul(-1, relq)  # fix the sign so theta lies in (0, pi)
    p, q = relq
    vn = math.sqrt(max(qdot(p, p) - p[0] ** 2, 0))
    if vn < 0.1:
        continue  # skip near-degenerate rotations
    theta = 2 * math.atan2(vn, p[0])
    l = qsmul(1 / vn, (0,) + p[1:])
    sh, ch = math.sin(theta / 2), math.cos(theta / 2)
    d = -2 * q[0] / sh
    m = qsmul(1 / sh, qadd((0,) + q[1:], qsmul(-d / 2 * ch, l)))
    assert deq(screw(theta, d, l, m), relq)
    # endpoints: t = 0 gives sigma0, t = 1 gives sigma1 up to sign
    assert deq(dmul(s0, screw(0, 0, l, m)), s0)
    end = dmul(s0, screw(theta, d, l, m))
    assert deq(end, s1) or deq(end, dsmul(-1, s1))
    q1r = qmul(s0[0], qexp_pure(qsmul(theta / 2, l)))  # sign-fixed rotor of s1
    x0 = cross(l, m)  # foot of the screw axis
    axis_pt = qadd(x0, qsmul(random.gauss(0, 1), l))
    for t in (0.3, 0.7):
        sc = dmul(s0, screw(t * theta, t * d, l, m))
        assert deq(dmul(sc, dbar(sc)), D1)  # unit at every t
        # the real part is the slerp of the rotor parts
        assert qeq(sc[0], slerp_pow(s0[0], q1r, t))
        # points on the axis translate linearly by t d l
        assert deq(act(screw(t * theta, t * d, l, m), axis_pt),
                   point(qadd(axis_pt, qsmul(t * d, l))))

# --- check 11: degeneracies --------------------------------------------------

for _ in range(20):
    q0 = rand_rotor()
    q1 = qmul(q0, qexp_pure(qsmul(1e-8, rand_unit_pure())))
    for t in (0.3, 0.7):
        s = slerp_pow(q0, q1, t)
        assert abs(qdot(s, s) - 1) < 1e-9  # still unit
        assert qeq(s, q0, 1e-6)
    q1 = qmul(q0, qexp_pure(qsmul(math.pi / 2, rand_unit_pure())))
    assert abs(qdot(q0, q1)) < 1e-9  # relative rotation angle pi
    assert abs(qang(q0, q1) - qang(q0, qsmul(-1, q1))) < 1e-9

print("All checks passed.")
