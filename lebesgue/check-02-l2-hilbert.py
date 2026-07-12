"""Numerical checks for 02-l2-hilbert.md.

The ramp functions form an L^2 Cauchy sequence of continuous functions
whose limit is the discontinuous step 1_[1/2,1]; fattening the rationals
by intervals of widths eps/2^n gives a monotone Cauchy sequence of
Riemann-integrable indicators whose limit 1_U is dense yet of measure
<= eps; Cauchy-Schwarz; the parallelogram law holds in L^2 and fails in
L^1; orthonormality of the trigonometric system; Fourier coefficients of
f(x) = x with Bessel monotonicity and Parseval convergence; the Fourier
partial sum as best approximation; and Parseval yielding the Basel sum
1/1 + 1/4 + 1/9 + ... = pi^2/6.
"""

from fractions import Fraction
import math

import numpy as np

# ramp functions on [0,1]: f_n rises linearly from 0 to 1 over
# [1/2, 1/2 + 1/n]. Exact integrals: for m > n the functions differ only
# on an interval of width 1/n, so ||f_m - f_n||_2^2 <= 1/n (Cauchy),
# while ||f_n - step||_2^2 = 1/(3n) -> 0 with step = 1_[1/2,1].
x = np.linspace(0, 1, 200_001)
dx = x[1] - x[0]

def ramp(n):
    return np.clip(n * (x - 0.5), 0, 1)

def l2sq(u):
    return np.trapezoid(u * u, x)

diffs = [math.sqrt(l2sq(ramp(2 * n) - ramp(n))) for n in (10, 100, 1000)]
print("||f_2n - f_n||_2 decreases toward 0 (Cauchy in L^2):",
      all(a > b for a, b in zip(diffs, diffs[1:])) and diffs[-1] < 0.02)
step = (x >= 0.5).astype(float)
errs = [l2sq(ramp(n) - step) for n in (10, 100, 1000)]
print("||f_n - 1_[1/2,1]||_2^2 = 1/(3n) -> 0 (limit is the step):",
      all(np.isclose(e, 1 / (3 * n), atol=1e-4)
          for e, n in zip(errs, (10, 100, 1000))))

# fattened rationals: enumerate the rationals q_k in [0,1] and cover q_k
# by an open interval I_k of width eps/2^k (the covering from article 01).
# g_n = indicator of I_1 u ... u I_n is Riemann integrable;
# mu(union of all I_k) <= eps although the union is dense, and
# ||g_m - g_n||_2^2 <= sum_{k>n} eps/2^k = eps/2^n (Cauchy).
def rationals(count):
    seen, out = set(), []
    q = 1
    while len(out) < count:
        for p in range(q + 1):
            r = Fraction(p, q)
            if r not in seen:
                seen.add(r)
                out.append(r)
                if len(out) == count:
                    break
        q += 1
    return out

def union_measure(intervals):
    total, cur_a, cur_b = 0.0, None, None
    for a, b in sorted(intervals):
        if cur_b is None or a > cur_b:
            if cur_b is not None:
                total += cur_b - cur_a
            cur_a, cur_b = a, b
        else:
            cur_b = max(cur_b, b)
    return total + (cur_b - cur_a)

eps = 0.1
qs = rationals(2000)
ivals = [(float(q) - eps * 2.0**-k / 2, float(q) + eps * 2.0**-k / 2)
         for k, q in enumerate(qs, start=1)]
print("mu(I_1 u ... u I_n) <= eps (dense open set of small measure):",
      union_measure(ivals) <= eps)
gaps = np.diff(sorted(float(q) for q in qs))
print("the covered rationals fill [0,1] (max gap ~ 1/q_max):",
      gaps.max() < 0.02)
tails = [union_measure(ivals[n:]) for n in (5, 10, 20)]
print("||g_m - g_n||_2^2 <= eps/2^n (Cauchy):",
      all(t <= eps / 2**n for t, n in zip(tails, (5, 10, 20))))

# Cauchy-Schwarz on [0,1]: |<f,g>| <= ||f||_2 ||g||_2 for f = x^2,
# g = cos(3x), with strict inequality (f, g not proportional)
f, g = x**2, np.cos(3 * x)
inner = np.trapezoid(f * g, x)
print("Cauchy-Schwarz |<f,g>| <= ||f|| ||g||:",
      abs(inner) < math.sqrt(l2sq(f)) * math.sqrt(l2sq(g)))

# parallelogram law: ||f+g||^2 + ||f-g||^2 = 2(||f||^2 + ||g||^2)
# holds in L^2 but fails in L^1 for f = 1_[0,1/2], g = 1_[1/2,1]
# (exact values: 1^2 + 1^2 = 2 vs 2(1/4 + 1/4) = 1)
lhs = l2sq(f + g) + l2sq(f - g)
rhs = 2 * (l2sq(f) + l2sq(g))
print("parallelogram law holds in L^2:", np.isclose(lhs, rhs))
lhs1 = 1.0**2 + 1.0**2  # ||f+g||_1 = 1, ||f-g||_1 = 1
rhs1 = 2 * (0.5**2 + 0.5**2)  # ||f||_1 = ||g||_1 = 1/2
print("parallelogram law fails in L^1 (2 != 1):", lhs1 == 2 and rhs1 == 1)

# orthonormality of the trigonometric system on [-pi, pi]:
# {1/sqrt(2pi), cos(nx)/sqrt(pi), sin(nx)/sqrt(pi)} has Gram matrix I
t = np.linspace(-math.pi, math.pi, 400_001)
basis = [np.full_like(t, 1 / math.sqrt(2 * math.pi))]
for n in range(1, 6):
    basis.append(np.cos(n * t) / math.sqrt(math.pi))
    basis.append(np.sin(n * t) / math.sqrt(math.pi))
gram = np.array([[np.trapezoid(u * v, t) for v in basis] for u in basis])
print("trigonometric system is orthonormal (Gram = identity):",
      np.allclose(gram, np.eye(len(basis)), atol=1e-6))

# Fourier coefficients of f(x) = x on [-pi, pi]:
# <f, e_0> = <f, cos(nx)/sqrt(pi)> = 0 (odd function), and
# <f, sin(nx)/sqrt(pi)> = 2 sqrt(pi) (-1)^(n+1) / n
N = 60
coefs = np.array([np.trapezoid(t * np.sin(n * t) / math.sqrt(math.pi), t)
                  for n in range(1, N + 1)])
exact = np.array([2 * math.sqrt(math.pi) * (-1)**(n + 1) / n
                  for n in range(1, N + 1)])
print("<x, sin(nx)/sqrt(pi)> = 2 sqrt(pi) (-1)^(n+1)/n:",
      np.allclose(coefs, exact, atol=1e-5))
print("even components vanish (x is odd):",
      np.isclose(np.trapezoid(t / math.sqrt(2 * math.pi), t), 0, atol=1e-9)
      and all(np.isclose(np.trapezoid(t * np.cos(n * t), t), 0, atol=1e-9)
              for n in range(1, 6)))

# Bessel: partial sums of c_n^2 are increasing and bounded by
# ||f||^2 = int x^2 dx = 2 pi^3 / 3; the first values are
# 4pi ~ 12.57, +4pi/4 ~ 15.71, +4pi/9 ~ 17.10 as quoted in the article
norm_sq = 2 * math.pi**3 / 3
partials = np.cumsum(exact**2)
print("Bessel partial sums 12.57, 15.71, 17.10, ...:",
      np.allclose(partials[:3], [12.57, 15.71, 17.10], atol=5e-3))
print("Bessel: increasing and <= ||f||^2 = 2 pi^3/3:",
      all(a < b for a, b in zip(partials, partials[1:]))
      and partials[-1] < norm_sq)

# Parseval: ||f - S_N||^2 = ||f||^2 - sum c_n^2 -> 0, checked both by
# the identity and by direct numerical integration of (f - S_N)^2
def partial_sum(N):
    s = np.zeros_like(t)
    for n in range(1, N + 1):
        s += exact[n - 1] * np.sin(n * t) / math.sqrt(math.pi)
    return s

direct = np.trapezoid((t - partial_sum(20))**2, t)
print("||f - S_20||^2 = ||f||^2 - sum_{n<=20} c_n^2:",
      np.isclose(direct, norm_sq - partials[19], atol=1e-6))
tail_errs = [norm_sq - partials[n - 1] for n in (10, 30, 60)]
print("||f - S_N||^2 -> 0 (Parseval):",
      all(a > b > 0 for a, b in zip(tail_errs, tail_errs[1:]))
      and tail_errs[-1] < 4 * math.pi / 59)

# best approximation: S_N is the orthogonal projection onto trig
# polynomials, so perturbing any coefficient increases the L^2 error by
# exactly the squared perturbation
base_err = np.trapezoid((t - partial_sum(5))**2, t)
for delta, n in ((0.5, 1), (0.3, 3), (-0.7, 5)):
    pert = partial_sum(5) + delta * np.sin(n * t) / math.sqrt(math.pi)
    err = np.trapezoid((t - pert)**2, t)
    print(f"perturbing c_{n} by {delta} adds {delta}^2 to the error:",
          np.isclose(err - base_err, delta**2, atol=1e-6))

# Basel problem: Parseval for f(x) = x reads
# 2 pi^3/3 = sum 4 pi / n^2, i.e. sum 1/n^2 = pi^2/6
basel = sum(1 / n**2 for n in range(1, 200_001))
print("sum 1/n^2 = pi^2/6 (Basel problem):",
      np.isclose(basel, math.pi**2 / 6, atol=1e-5))
print("Parseval: sum 4 pi/n^2 = 2 pi^3/3 = int x^2 dx:",
      np.isclose(4 * math.pi * basel, norm_sq, atol=1e-4)
      and np.isclose(np.trapezoid(t * t, t), norm_sq, atol=1e-6))
