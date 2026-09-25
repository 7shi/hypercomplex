"""Real Clifford algebras Cl_{p,q}(R) on bitmask blades, with sympy coefficients.

Generators e_0..e_{n-1} (n = p + q) square to +1 for a < p and to -1 for
a >= p. A basis blade is a bitmask; e.g. 0b101 is e_0 e_2. A multivector keeps
the bitmask `neg` of the generators squaring to -1, so Cl_{n,0} has neg = 0 and
Cl_{1,3} (e_0^2 = 1, e_1^2 = e_2^2 = e_3^2 = -1) has neg = 0b1110.
"""

import itertools
import random

import sympy as sp


def blade_mul(a, b, neg=0):
    """Product of basis blades (bitmasks): returns (sign, mask)."""
    s = 0
    t = a >> 1
    while t:
        s += bin(t & b).count("1")
        t >>= 1
    s += bin(a & b & neg).count("1")
    return (-1 if s & 1 else 1), a ^ b


def grade(m):
    return bin(m).count("1")


def _neg(a, b):
    if a and b and a != b:
        raise ValueError("multivectors from different signatures")
    return a | b


class MV:
    def __init__(self, d=None, neg=0):
        self.d = {k: v for k, v in (d or {}).items() if v != 0}
        self.neg = neg

    def __add__(self, o):
        o = mv(o, self.neg)
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, 0) + v
        return MV(d, _neg(self.neg, o.neg))

    __radd__ = __add__

    def __neg__(self):
        return MV({k: -v for k, v in self.d.items()}, self.neg)

    def __sub__(self, o):
        return self + (-mv(o, self.neg))

    def __rsub__(self, o):
        return mv(o, self.neg) - self

    def __mul__(self, o):
        o = mv(o, self.neg)
        neg = _neg(self.neg, o.neg)
        d = {}
        for a, x in self.d.items():
            for b, y in o.d.items():
                s, m = blade_mul(a, b, neg)
                d[m] = d.get(m, 0) + s * x * y
        return MV(d, neg)

    def __rmul__(self, o):
        return mv(o, self.neg) * self

    def __truediv__(self, c):
        return MV({k: v / c for k, v in self.d.items()}, self.neg)

    def map(self, f):
        return MV({k: f(v) for k, v in self.d.items()}, self.neg)

    def grade(self, k):
        return MV({m: v for m, v in self.d.items() if grade(m) == k}, self.neg)

    def rev(self):
        return MV({m: (-v if (grade(m) * (grade(m) - 1) // 2) % 2 else v) for m, v in self.d.items()}, self.neg)

    def scalar(self):
        return self.d.get(0, 0)

    def simp(self):
        return self.map(lambda v: sp.simplify(sp.expand(v)))


def mv(o, neg=0):
    return o if isinstance(o, MV) else MV({0: sp.sympify(o)}, neg)


def zero_expr(v):
    v = sp.cancel(sp.together(sp.expand(v)))
    return v == 0 or sp.simplify(v) == 0


def eq(A, B):
    return all(zero_expr(v) for v in (mv(A) - mv(B)).d.values())


class Alg:
    """Cl_{p,q} with coordinates x_0..x_{n-1} and the Dirac operator.

    D = sum e^a d_a uses the reciprocal frame e^a = e_a^{-1} = (e_a^2) e_a,
    so that D x = n and D^2 = sum (e_a^2) d_a^2 (`lap`).
    """

    def __init__(self, p, q=0, name="x"):
        self.p, self.q = p, q
        self.n = n = p + q
        self.neg = ((1 << n) - 1) ^ ((1 << p) - 1)
        self.sq = [1] * p + [-1] * q
        self.e = [MV({1 << a: 1}, self.neg) for a in range(n)]
        self.er = [self.sq[a] * self.e[a] for a in range(n)]
        self.X = sp.symbols(f"{name}0:{n}", real=True)
        self.x = sum((self.X[a] * self.e[a] for a in range(n)), MV({}, self.neg))
        self.r2 = sum(self.sq[a] * self.X[a] ** 2 for a in range(n))
        self.I = MV({(1 << n) - 1: 1}, self.neg)

    def zero(self):
        return MV({}, self.neg)

    def d(self, F, a):
        return F.map(lambda v: sp.diff(v, self.X[a]))

    def D(self, F):
        return sum((self.er[a] * self.d(F, a) for a in range(self.n)), self.zero())

    def Dr(self, F):
        return sum((self.d(F, a) * self.er[a] for a in range(self.n)), self.zero())

    def lap(self, F):
        return sum((self.sq[a] * self.d(self.d(F, a), a) for a in range(self.n)), self.zero())

    def rnd(self, seed, deg=2, grades=None):
        rnd = random.Random(seed)
        mons = [m for m in itertools.product(range(deg + 1), repeat=self.n) if sum(m) <= deg]
        blades = [b for b in range(1 << self.n) if grades is None or grade(b) in grades]
        return MV({b: sum(rnd.randint(-3, 3) * sp.prod([self.X[a] ** p for a, p in enumerate(m)])
                          for m in mons) for b in blades}, self.neg)
