"""Numerical checks for 04-brownian-ito.md.

Scaled random walks sqrt(dt) S_{t/dt} approach the normal distribution
(CLT); Brownian increments over disjoint intervals have variance equal
to the elapsed time and are uncorrelated; the quadratic variation sum
(dB)^2 converges to t in L^2 with E[(Q_n - t)^2] = 2 t^2 / n, while
the total variation sum |dB| = sqrt(2 n t / pi) diverges; the Ito
isometry E[(int H dB)^2] = E[int H^2 dt] holds for the left-endpoint
approximation of H = B; the evaluation point matters: left-endpoint
sums of int B dB converge to (B_T^2 - T)/2 while midpoint sums
converge to B_T^2 / 2 (Stratonovich), the two limits differing by the
correction T/2; the Ito integral has mean zero; Ito's lemma holds
pathwise for f(x) = x^4 with error shrinking as the partition refines;
the Euler-Maruyama scheme for geometric Brownian motion converges to
the exact solution S_0 exp((mu - sigma^2/2) t + sigma B_t); and the
moments produced by Ito's lemma, E[B_t^4] = 3 t^2 and E[B_t^6] =
15 t^3, hold.
"""

import math

import numpy as np

rng = np.random.default_rng(0)
T = 1.0

# scaling limit: the walk sqrt(dt) S_{t/dt} at t = 1 takes n steps of
# size +-sqrt(1/n); its CDF approaches the standard normal CDF (CLT)
n_clt, paths = 10_000, 200_000
S = (2.0 * rng.binomial(n_clt, 0.5, size=paths) - n_clt) / math.sqrt(n_clt)
Phi = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))
print("scaled walk at t=1: P(S <= x) matches the normal CDF (CLT):",
      all(np.isclose(np.mean(S <= x), Phi(x), atol=0.01)
          for x in (-1.5, -0.5, 0.0, 0.5, 1.5)))

# Brownian increments: Var(B_t - B_s) = t - s, and increments over
# disjoint intervals are uncorrelated (independent increments)
n, paths = 400, 40_000
dt = T / n
B = np.cumsum(rng.normal(0.0, math.sqrt(dt), size=(paths, n)), axis=1)
inc1 = B[:, 99]                # B_{1/4}
inc2 = B[:, 299] - B[:, 99]    # B_{3/4} - B_{1/4}
print("Var(B_t - B_s) = t - s:",
      np.isclose(inc1.var(), 0.25, rtol=0.05)
      and np.isclose(inc2.var(), 0.50, rtol=0.05))
print("disjoint increments are uncorrelated:",
      abs(np.corrcoef(inc1, inc2)[0, 1]) < 0.02)
print("E[B_t^4] = 3t^2 at t = 1/4, 1:",
      np.isclose(np.mean(inc1**4), 3 * 0.25**2, rtol=0.05)
      and np.isclose(np.mean(B[:, -1] ** 4), 3.0, rtol=0.05))
print("E[B_1^6] = 15t^3 = 15:",
      np.isclose(np.mean(B[:, -1] ** 6), 15.0, rtol=0.1))

# quadratic variation: Q_n = sum (dB)^2 has mean T and L^2 error
# E[(Q_n - T)^2] = 2 T^2 / n -> 0, while the total variation
# sum |dB| has mean sqrt(2 n T / pi) -> infinity
ns = (100, 1000, 10_000)
qv_ok, tv_ok = True, True
for n_k in ns:
    m = 200_000 // n_k * 10
    dB = rng.normal(0.0, math.sqrt(T / n_k), size=(m, n_k))
    Q = (dB**2).sum(axis=1)
    TV = np.abs(dB).sum(axis=1)
    qv_ok = qv_ok and np.isclose(np.mean((Q - T) ** 2), 2 * T**2 / n_k,
                                 rtol=0.1)
    tv_ok = tv_ok and np.isclose(TV.mean(),
                                 math.sqrt(2 * n_k * T / math.pi), rtol=0.02)
print("quadratic variation: E[(Q_n - T)^2] = 2T^2/n -> 0:", qv_ok)
print("total variation: E[sum |dB|] = sqrt(2nT/pi) diverges:", tv_ok)

# one fine simulation reused below: B on a grid of 800 steps, with the
# coarse grid (400 steps) taking every other point so that midpoints
# are available
n2, n, paths = 800, 400, 20_000
dt = T / n
B0 = np.hstack([np.zeros((paths, 1)),
                np.cumsum(rng.normal(0.0, math.sqrt(T / n2),
                                     size=(paths, n2)), axis=1)])
Bc = B0[:, ::2]           # B at the coarse grid points t_0, ..., t_n
dBc = np.diff(Bc, axis=1)
left = Bc[:, :-1]         # B at left endpoints t_{k-1}
mid = B0[:, 1::2]         # B at the midpoints of the coarse intervals
BT = B0[:, -1]

# Ito isometry for the simple process H = B (left-endpoint steps):
# E[(sum B dB)^2] = E[sum B^2 dt] (= T^2/2 in the limit)
I_left = (left * dBc).sum(axis=1)
rhs = np.mean((left**2).sum(axis=1) * dt)
print("Ito isometry: E[(int B dB)^2] = E[int B^2 dt]:",
      np.isclose(np.mean(I_left**2), rhs, rtol=0.05)
      and np.isclose(np.mean(I_left**2), T**2 / 2, rtol=0.05))
print("Ito integral has mean zero: E[int B dB] = 0:",
      np.isclose(I_left.mean(), 0.0, atol=0.02))

# evaluation point: left sums converge pathwise to (B_T^2 - T)/2 (Ito)
# and midpoint sums to B_T^2/2 (Stratonovich); the limits differ by T/2
I_mid = (mid * dBc).sum(axis=1)
rms_left = math.sqrt(np.mean((I_left - (BT**2 - T) / 2) ** 2))
rms_mid = math.sqrt(np.mean((I_mid - BT**2 / 2) ** 2))
print("left sums -> (B_T^2 - T)/2, midpoint sums -> B_T^2/2:",
      rms_left < 0.1 and rms_mid < 0.1)
print("the two evaluations differ by the correction T/2:",
      np.isclose(np.mean(I_mid - I_left), T / 2, atol=0.02))

# Ito's lemma pathwise for f(x) = x^4:
# B_T^4 = 4 sum B^3 dB + 6 sum B^2 dt, with error -> 0 as the
# partition refines (compare a coarse and a fine partition)
def ito_lemma_rms(step):
    Bg = B0[:, ::step]
    dBg = np.diff(Bg, axis=1)
    lft = Bg[:, :-1]
    rhs = 4 * (lft**3 * dBg).sum(axis=1) \
        + 6 * (lft**2).sum(axis=1) * (step * T / n2)
    return math.sqrt(np.mean((BT**4 - rhs) ** 2))

err_coarse, err_fine = ito_lemma_rms(16), ito_lemma_rms(2)
print("Ito's lemma for x^4: pathwise error shrinks with refinement:",
      err_fine < err_coarse / 2 and err_fine < 0.05 * math.sqrt(105))

# geometric Brownian motion: Euler-Maruyama S += mu S dt + sigma S dB
# converges to the exact solution S_0 exp((mu - sigma^2/2) t + sigma B_t)
mu, sigma, S0 = 0.1, 0.4, 1.0
S_exact = S0 * np.exp((mu - sigma**2 / 2) * T + sigma * BT)

def euler_maruyama_rms(step):
    Bg = B0[:, ::step]
    dBg = np.diff(Bg, axis=1)
    h = step * T / n2
    S = np.full(paths, S0)
    for k in range(dBg.shape[1]):
        S = S * (1 + mu * h + sigma * dBg[:, k])
    return math.sqrt(np.mean((S - S_exact) ** 2))

em_coarse, em_fine = euler_maruyama_rms(16), euler_maruyama_rms(1)
print("GBM: Euler-Maruyama converges to S_0 exp((mu-s^2/2)t + s B_t):",
      em_fine < em_coarse / 2 and em_fine < 0.02)
