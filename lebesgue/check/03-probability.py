"""Numerical checks for 03-probability.md.

Binary digits of a uniform point on [0,1] behave as independent fair
coin flips (each digit has probability 1/2, joint probabilities
factorize) and their running averages converge to 1/2 (Borel's normal
number theorem via the strong law of large numbers); E[XY] = E[X] E[Y]
for independent variables; Lebesgue-Stieltjes integration E[X] = int x
dF unifies discrete, continuous and mixed distributions; Radon-Nikodym
nu(A) = int_A h dmu for an exponential density; conditional expectation
on a die as the blockwise mean with orthogonality, best-prediction and
tower properties; exact path enumeration showing S_n, S_n^2 - n and a
martingale transform are martingales; the gambler's ruin solved exactly:
P(hit a before -b) = b/(a+b) and E[T] = ab; and the doubling strategy
with a finite horizon has expectation exactly 0.
"""

from fractions import Fraction
from itertools import product
import math

import numpy as np

# binary digits as functions on [0,1]: d_n(w) = floor(2^n w) mod 2.
# On a midpoint grid the measure of {d_n = 1} is 1/2 for each n, and
# P(d_1 = e_1, ..., d_k = e_k) = 2^-k factorizes as a product.
N = 2**16
w = (np.arange(N) + 0.5) / N

def digit(n):
    return (np.floor(2.0**n * w) % 2).astype(int)

digits = [digit(n) for n in range(1, 11)]
print("P(d_n = 1) = 1/2 for each digit:",
      all(np.isclose(d.mean(), 0.5) for d in digits))
joint = np.mean((digits[0] == 1) & (digits[1] == 0) & (digits[2] == 1))
print("P(d_1=1, d_2=0, d_3=1) = 1/8 (independence):",
      np.isclose(joint, 1 / 8))
print("E[d_1 d_2] = E[d_1] E[d_2] = 1/4:",
      np.isclose(np.mean(digits[0] * digits[1]), 0.25))

# strong law of large numbers / Borel normal numbers: for random
# points, the running average of the first n binary digits approaches
# 1/2, and the fraction of points straying from 1/2 shrinks with n
rng = np.random.default_rng(0)
bits = rng.integers(0, 2, size=(20_000, 2000))  # digits of 20000 points
avgs = {n: bits[:, :n].mean(axis=1) for n in (20, 200, 2000)}
stray = [np.mean(np.abs(avgs[n] - 0.5) > 0.05) for n in (20, 200, 2000)]
print("P(|avg of n digits - 1/2| > 0.05) decreases (SLLN):",
      stray[0] > stray[1] > stray[2] and stray[2] < 1e-3)

# independence for dice: E[XY] = E[X] E[Y] = 3.5^2 by exact enumeration
pairs = [(i, j) for i in range(1, 7) for j in range(1, 7)]
exy = Fraction(sum(i * j for i, j in pairs), 36)
print("two dice: E[XY] = E[X]E[Y] = 12.25:",
      exy == Fraction(49, 4))

# Lebesgue-Stieltjes expectation E[X] = int x dF as a Stieltjes sum
# sum x_i (F(x_{i+1}) - F(x_i)) on a fine grid, for three distributions:
# die (steps only), uniform on [0,1] (density only), and the rain
# example (jump 0.6 at 0 plus density 0.2 on [0,2]) with E[X] = 0.4
def stieltjes(F, a, b, n=400_001):
    x = np.linspace(a, b, n)
    return np.sum(x[:-1] * np.diff(F(x)))

F_die = lambda x: np.floor(np.clip(x, 0, 6)) / 6
F_unif = lambda x: np.clip(x, 0, 1)
F_rain = lambda x: np.where(x < 0, 0.0, np.minimum(0.6 + 0.2 * x, 1.0))
print("die: int x dF = sum x p(x) = 3.5:",
      np.isclose(stieltjes(F_die, -1, 7), 3.5, atol=1e-4))
print("uniform: int x dF = int x dx = 1/2:",
      np.isclose(stieltjes(F_unif, -1, 2), 0.5, atol=1e-4))
print("rain (mixed): int x dF = 0*0.6 + int_0^2 0.2 x dx = 0.4:",
      np.isclose(stieltjes(F_rain, -1, 3), 0.4, atol=1e-4))

# Radon-Nikodym: for nu = exponential(1) distribution and mu = Lebesgue,
# nu((a,b]) = F(b) - F(a) = int_a^b e^-x dx for several intervals
F_exp = lambda x: 1 - math.exp(-x)
ok = True
for a, b in ((0.0, 1.0), (0.5, 2.5), (1.0, 4.0)):
    x = np.linspace(a, b, 100_001)
    ok = ok and np.isclose(np.trapezoid(np.exp(-x), x),
                           F_exp(b) - F_exp(a), atol=1e-8)
print("Radon-Nikodym: nu((a,b]) = int_a^b (dnu/dmu) dx:", ok)

# conditional expectation on a fair die X(w) = w, w in {1,...,6}:
# given parity, E[X|G] = 4 on evens and 3 on odds
X = np.arange(1, 7, dtype=float)
even = np.array([False, True, False, True, False, True])
Y = np.where(even, X[even].mean(), X[~even].mean())
print("E[X|parity] = 4 on evens, 3 on odds:",
      X[even].mean() == 4 and X[~even].mean() == 3)
print("orthogonality: sum_B (X - E[X|G]) = 0 on each block:",
      np.isclose((X - Y)[even].sum(), 0)
      and np.isclose((X - Y)[~even].sum(), 0))

# best prediction: perturbing the value on a block increases the mean
# squared error by exactly P(B) * delta^2
mse = np.mean((X - Y) ** 2)
ok = True
for delta, block in ((0.5, even), (-0.8, ~even)):
    Yp = Y + delta * block
    ok = ok and np.isclose(np.mean((X - Yp) ** 2) - mse,
                           block.mean() * delta**2)
print("perturbing E[X|G] adds P(B) delta^2 to the error (projection):",
      ok)

# tower rule: G from {1,2},{3,4},{5,6} and coarser H from
# {1,2,3,4},{5,6} satisfy E[E[X|G]|H] = E[X|H]
G_blocks = [np.array([1, 1, 0, 0, 0, 0], bool),
            np.array([0, 0, 1, 1, 0, 0], bool),
            np.array([0, 0, 0, 0, 1, 1], bool)]
H_blocks = [np.array([1, 1, 1, 1, 0, 0], bool),
            np.array([0, 0, 0, 0, 1, 1], bool)]

def condexp(Z, blocks):
    out = np.zeros_like(Z)
    for B in blocks:
        out[B] = Z[B].mean()
    return out

lhs = condexp(condexp(X, G_blocks), H_blocks)
rhs = condexp(X, H_blocks)
print("tower rule E[E[X|G]|H] = E[X|H]:", np.allclose(lhs, rhs))

# martingales by exact enumeration of all sign paths of length 10:
# grouping paths by their first n steps, the conditional means satisfy
# E[S_{n+1}|F_n] = S_n and E[S_{n+1}^2 - (n+1)|F_n] = S_n^2 - n, and the
# transform (H.S)_n with predictable H_k = 2 if S_{k-1} < 0 else 1 has
# E[(H.S)_n] = 0 for all n
paths = np.array(list(product((-1, 1), repeat=10)))
S = np.cumsum(paths, axis=1)
ok_mart, ok_sq = True, True
for n in (0, 3, 6, 9):
    S_n = S[:, n - 1] if n > 0 else np.zeros(len(paths))
    for prefix in set(map(tuple, paths[:, :n])):
        idx = np.all(paths[:, :n] == prefix, axis=1)
        ok_mart = ok_mart and np.isclose(S[idx, n].mean(), S_n[idx][0])
        ok_sq = ok_sq and np.isclose((S[idx, n] ** 2).mean() - (n + 1),
                                     S_n[idx][0] ** 2 - n)
print("E[S_{n+1}|F_n] = S_n (martingale):", ok_mart)
print("E[S_{n+1}^2-(n+1)|F_n] = S_n^2-n (martingale):", ok_sq)

S_prev = np.hstack([np.zeros((len(paths), 1)), S[:, :-1]])
H = np.where(S_prev < 0, 2.0, 1.0)  # bets depend on the past only
HS = np.cumsum(H * paths, axis=1)
print("martingale transform: E[(H.S)_n] = 0 for all n:",
      np.allclose(HS.mean(axis=0), 0))

# gambler's ruin with a = 3, b = 2: exact linear systems for the
# hitting probability p_k of +a and the expected duration t_k, from
# p_k = (p_{k-1} + p_{k+1})/2 and t_k = 1 + (t_{k-1} + t_{k+1})/2;
# optional stopping gives p_0 = b/(a+b) = 2/5 and t_0 = ab = 6
a, b = 3, 2
inner = list(range(-b + 1, a))  # transient states
n_st = len(inner)
A = np.zeros((n_st, n_st))
rhs_p = np.zeros(n_st)
rhs_t = np.full(n_st, -1.0)
for i, k in enumerate(inner):
    A[i, i] = -1
    for step in (-1, 1):
        if k + step == a:
            rhs_p[i] -= 0.5  # absorbed at +a with value 1
        elif k + step != -b:
            A[i, inner.index(k + step)] = 0.5
p = np.linalg.solve(A, rhs_p)
t = np.linalg.solve(A, rhs_t)
p0, t0 = p[inner.index(0)], t[inner.index(0)]
print("ruin: P(hit +a before -b) = b/(a+b) = 2/5:",
      np.isclose(p0, b / (a + b)))
print("ruin: E[T] = ab = 6:", np.isclose(t0, a * b))
print("optional stopping: E[S_T] = pa - (1-p)b = 0 and E[S_T^2] = E[T]:",
      np.isclose(p0 * a - (1 - p0) * b, 0)
      and np.isclose(p0 * a**2 + (1 - p0) * b**2, t0))

# doubling strategy with horizon N: win at round k (probability 2^-k)
# nets +1, never winning (probability 2^-N) loses 2^N - 1, and the
# expectation is exactly 0 while P(win) = 1 - 2^-N -> 1
Nh = 20
e = sum(Fraction(1, 2**k) * 1 for k in range(1, Nh + 1)) \
    - Fraction(1, 2**Nh) * (2**Nh - 1)
print("doubling strategy: P(win) = 1 - 2^-N but E = 0 exactly:",
      e == 0 and 1 - Fraction(1, 2**Nh) > Fraction(999, 1000))
