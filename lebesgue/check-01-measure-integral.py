"""Numerical checks for 01-measure-integral.md.

Riemann sums of the Dirichlet function depending on sample points
(rational vs irrational); covering a countable set by intervals of
widths eps/2^n whose total length is the geometric series eps; the
horizontal-slicing simple functions s_n for f(x)=x^2 and the monotone
convergence of their integrals to 1/3; the escaping-mass sequence
f_n = n*1_(0,1/n) where integral and pointwise limit disagree; and the
divergence of the lower bound 1/x - 1 for any dominating function g.
"""

from fractions import Fraction
import math

import numpy as np

# Riemann sums of the Dirichlet function 1_Q on [0,1]:
# with rational sample points every sum is 1, with irrational sample
# points every sum is 0, so refining the partition never converges.
# Sample points are exact Fractions (rational) or carry an offset
# 1/(sqrt(2)*N) known analytically to be irrational.
for N in (10, 100, 1000):
    rational_sum = sum(Fraction(1) * Fraction(1, N) for _ in range(N))  # f(k/N) = 1
    irrational_sum = sum(0 * (1 / N) for _ in range(N))  # f(k/N + 1/(sqrt(2)N)) = 0
    print(f"Dirichlet Riemann sums (N={N}): upper=1, lower=0:",
          rational_sum == 1 and irrational_sum == 0)

# covering a countable set: widths eps/2^n sum to eps (geometric series)
eps = 0.001
total = sum(eps / 2**n for n in range(1, 60))
print("sum eps/2^n = eps (countable set has measure zero):",
      np.isclose(total, eps))

# horizontal slicing of f(x) = x^2 on [0,1): the simple function
# s_n = floor(2^n f)/2^n has level sets A_k = [sqrt(k/2^n), sqrt((k+1)/2^n))
# with exact measures, so int s_n = sum_k (k/2^n) * mu(A_k).
def simple_integral(n):
    h = 2**n
    return sum((k / h) * (math.sqrt((k + 1) / h) - math.sqrt(k / h))
               for k in range(h))

vals = [simple_integral(n) for n in range(1, 21)]
print("int s_1 ~ 0.146, int s_2 ~ 0.232:",
      np.isclose(vals[0], 0.5 * (1 - 1 / math.sqrt(2)))
      and np.isclose(vals[0], 0.146, atol=5e-4)
      and np.isclose(vals[1], 0.232, atol=5e-4))
print("int s_n is monotone increasing:",
      all(a < b for a, b in zip(vals, vals[1:])))
print("int s_n -> 1/3 (monotone convergence to the Riemann value):",
      np.isclose(vals[-1], 1 / 3, atol=1e-5))

# error bound: s_n <= f with f - s_n <= 1/2^n at sample points
x = np.linspace(0, 1, 1001)[:-1]
for n in (1, 4, 8):
    s = np.floor(2**n * x**2) / 2**n
    ok = np.all(s <= x**2) and np.all(x**2 - s <= 1 / 2**n)
    print(f"0 <= f - s_{n} <= 1/2^{n} pointwise:", ok)

# escaping mass: f_n = n*1_(0,1/n) has integral n * (1/n) = 1 for all n,
# but f_n(x) -> 0 pointwise, so lim int f_n = 1 != 0 = int lim f_n
ns = np.arange(1, 101)
print("int f_n = 1 for all n:", np.allclose(ns * (1 / ns), 1))
print("f_n(x) -> 0 pointwise (f_n(x) = 0 for n > 1/x):",
      all(n * (x < 1 / n) == 0
          for x in (0.9, 0.5, 0.1, 0.001) for n in range(int(1 / x) + 1, int(1 / x) + 100)))

# no dominating function: any g >= f_n for all n satisfies
# g(x) >= sup_n f_n(x) >= 1/x - 1 on (0,1] (equality only at x = 1/n),
# and int_0^1 (1/x - 1) dx diverges like log(1/a) as a -> 0
def sup_fn(x):
    return math.ceil(1 / x) - 1  # largest n with x < 1/n

xs = np.linspace(1e-4, 1, 10000)
print("sup_n f_n(x) >= 1/x - 1 on (0,1]:",
      all(sup_fn(x) >= 1 / x - 1 for x in xs))
print("strict inequality between the endpoints 1/n:",
      all(sup_fn(x) > 1 / x - 1 for x in (0.9, 0.4, 0.3, 0.15, 0.0011)))
partials = [math.log(1 / a) + a - 1 for a in (1e-2, 1e-4, 1e-8)]  # int_a^1 (1/x-1) dx
print("int_a^1 (1/x - 1) dx = log(1/a) + a - 1 -> infinity:",
      all(a < b for a, b in zip(partials, partials[1:])) and partials[-1] > 17)

# sophomore's dream: int_0^1 x^-x dx = sum 1/n^n by term-by-term
# integration of exp(-x ln x), justified by monotone convergence
# (all terms (-x ln x)^n / n! are nonnegative on (0,1))
xs = np.linspace(1e-12, 1, 2_000_001)
for n in (1, 2, 3):
    term = np.trapezoid(xs**n * (-np.log(xs))**n, xs)
    print(f"int x^{n} (-ln x)^{n} dx = {n}!/{n+1}^{n+1}:",
          np.isclose(term, math.factorial(n) / (n + 1)**(n + 1), atol=1e-6))
series = sum(1 / n**n for n in range(1, 30))
integral = np.trapezoid(xs**(-xs), xs)
print("int_0^1 x^-x dx = sum 1/n^n (sophomore's dream):",
      np.isclose(integral, series, atol=1e-5))
terms = [(-0.5 * math.log(0.5))**n / math.factorial(n) for n in range(30)]
print("partial sums are monotone increasing (terms nonnegative):",
      all(t >= 0 for t in terms)
      and np.isclose(sum(terms), 0.5**-0.5))
