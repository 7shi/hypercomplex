"""Checks for spherical-coords.md.

Verifies symbolically (sympy):

1. The standard n-dimensional spherical coordinates (n = 2..6) give a unit
   vector.
2. 3D nesting (cos φ, sin φ cos φ', sin φ sin φ') is a unit vector.
3. 4D sequential nesting 1:(1:2) with angles φ1, φ2, φ3 is a unit vector,
   and the variant with sin φ4 in the last component is not.
4. 4D equal split 2:2 is a unit vector, its complex form
   (cos φ1 e^{iφ2}, sin φ1 e^{iφ3}) has norm 1, and the first-component
   realization and the Hadamard product decomposition hold.
5. 6D split 2:(2:2) is a unit vector, its complex form (3-level state)
   has norm 1, and the realization and Hadamard decomposition hold.
6. The unit quaternion rotor cos(θ/2) + sin(θ/2)(x i + y j + z k) has
   norm 1 when x²+y²+z² = 1 (the 1:(1:2) pattern with φ1 = θ/2).
7. The Hopf map of (cos φ1 e^{iφ2}, sin φ1 e^{iφ3}) is
   (sin 2φ1 cos δ, sin 2φ1 sin δ, cos 2φ1) with δ = φ3 - φ2: the zenith
   angle is 2φ1 and the azimuth is δ. Global phase leaves ΨΨ† unchanged.
"""

import sympy as sp

p = sp.symbols("phi1:7", real=True)
a, b, th = sp.symbols("alpha beta theta", real=True)
c, s = sp.cos, sp.sin


def is_unit(v):
    return sp.simplify(sp.trigsimp(sum(sp.expand(x * sp.conjugate(x)) for x in v)) - 1) == 0


# 1. standard spherical coordinates
for n in range(2, 7):
    v = []
    for k in range(n):
        x = sp.Integer(1)
        for j in range(k):
            x *= s(p[j])
        if k < n - 1:
            x *= c(p[k])
        v.append(x)
    assert is_unit(v), n
print("1. standard spherical coordinates (n = 2..6) give unit vectors: OK")

# 2. 3D nesting
assert is_unit([c(p[0]), s(p[0]) * c(p[1]), s(p[0]) * s(p[1])])
print("2. 3D nesting is a unit vector: OK")

# 3. 4D sequential nesting
v = [c(p[0]), s(p[0]) * c(p[1]), s(p[0]) * s(p[1]) * c(p[2]), s(p[0]) * s(p[1]) * s(p[2])]
assert is_unit(v)
v_bad = v[:3] + [s(p[0]) * s(p[1]) * s(p[3])]
assert not is_unit(v_bad)
print("3. 4D sequential nesting uses φ3 twice; sin φ4 would break the norm: OK")

# 4. 4D equal split
v = [c(p[0]) * c(p[1]), c(p[0]) * s(p[1]), s(p[0]) * c(p[2]), s(p[0]) * s(p[2])]
assert is_unit(v)
psi = [c(p[0]) * sp.exp(sp.I * p[1]), s(p[0]) * sp.exp(sp.I * p[2])]
assert all(sp.simplify(sp.expand_complex(psi[k] - (v[2 * k] + sp.I * v[2 * k + 1]))) == 0
           for k in range(2))
assert is_unit(psi)
real1 = [c(p[0]), s(p[0]) * sp.exp(sp.I * (p[2] - p[1]))]
assert all(sp.simplify(psi[k] - sp.exp(sp.I * p[1]) * real1[k]) == 0 for k in range(2))
had = [c(p[0]) * sp.exp(sp.I * p[1]), s(p[0]) * sp.exp(sp.I * p[2])]
assert all(sp.simplify(psi[k] - had[k]) == 0 for k in range(2))
print("4. 2:2 split, 2-level state, first-component realization, Hadamard: OK")

# 5. 6D split 2:(2:2)
v = [c(p[0]) * c(p[1]), c(p[0]) * s(p[1]),
     s(p[0]) * c(p[2]) * c(p[3]), s(p[0]) * c(p[2]) * s(p[3]),
     s(p[0]) * s(p[2]) * c(p[4]), s(p[0]) * s(p[2]) * s(p[4])]
assert is_unit(v)
psi = [c(p[0]) * sp.exp(sp.I * p[1]),
       s(p[0]) * c(p[2]) * sp.exp(sp.I * p[3]),
       s(p[0]) * s(p[2]) * sp.exp(sp.I * p[4])]
assert all(sp.simplify(sp.expand_complex(psi[k] - (v[2 * k] + sp.I * v[2 * k + 1]))) == 0
           for k in range(3))
assert is_unit(psi)
real1 = [c(p[0]), s(p[0]) * c(p[2]) * sp.exp(sp.I * (p[3] - p[1])),
         s(p[0]) * s(p[2]) * sp.exp(sp.I * (p[4] - p[1]))]
assert all(sp.simplify(psi[k] - sp.exp(sp.I * p[1]) * real1[k]) == 0 for k in range(3))
amp = [c(p[0]), s(p[0]) * c(p[2]), s(p[0]) * s(p[2])]
ph = [sp.exp(sp.I * p[1]), sp.exp(sp.I * p[3]), sp.exp(sp.I * p[4])]
assert all(sp.simplify(psi[k] - amp[k] * ph[k]) == 0 for k in range(3))
assert is_unit(amp)
print("5. 2:(2:2) split, 3-level state, realization, Hadamard: OK")

# 6. quaternion rotor
x, y, z = sp.symbols("x y z", real=True)
nrm = c(th / 2) ** 2 + s(th / 2) ** 2 * (x**2 + y**2 + z**2)
assert sp.simplify(nrm.subs(z**2, 1 - x**2 - y**2)) == 1
print("6. rotor cos(θ/2) + sin(θ/2)(xi+yj+zk) has norm 1: OK")

# 7. Hopf image of the 2-level state
psi = sp.Matrix([c(p[0]) * sp.exp(sp.I * p[1]), s(p[0]) * sp.exp(sp.I * p[2])])
q = sp.conjugate(psi[0]) * psi[1]
hop = [2 * sp.re(q), 2 * sp.im(q), sp.Abs(psi[0]) ** 2 - sp.Abs(psi[1]) ** 2]
d = p[2] - p[1]
target = [s(2 * p[0]) * c(d), s(2 * p[0]) * s(d), c(2 * p[0])]
assert all(sp.simplify(sp.expand_trig(sp.expand_complex(h - t))) == 0 for h, t in zip(hop, target))
w = sp.symbols("omega", real=True)
rho = psi * psi.H
rho2 = (sp.exp(sp.I * w) * psi) * (sp.exp(sp.I * w) * psi).H
assert sp.simplify(rho2 - rho) == sp.zeros(2, 2)
print("7. Hopf image has zenith 2φ1, azimuth δ = φ3-φ2; ΨΨ† is phase invariant: OK")
