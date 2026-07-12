"""Numerical checks for 05-black-scholes.md.

Geometric Brownian motion S_T = S_0 exp((mu - sigma^2/2)T + sigma B_T)
has mean S_0 e^{mu T}; under the risk-neutral measure (drift r) the
discounted call payoff averages to the Black-Scholes formula
S_0 N(d_1) - K e^{-rT} N(d_2); the same price is recovered by Girsanov
reweighting of real-world paths with the Doleans exponential
Z = exp(-theta B_T - theta^2 T/2); real-measure Monte Carlo depends on
mu and misses the BS price, while risk-neutral pricing does not; discrete
delta hedging replication error shrinks as the rebalancing grid is
refined; and the call delta equals N(d_1).
"""

import math

import numpy as np

rng = np.random.default_rng(0)

S0, K, r, sigma, T = 100.0, 100.0, 0.05, 0.2, 1.0
mu = 0.12
paths = 200_000


def norm_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def bs_call(s0=S0, k=K, rate=r, vol=sigma, t=T):
    sqrt_t = math.sqrt(t)
    d1 = (math.log(s0 / k) + (rate + 0.5 * vol**2) * t) / (vol * sqrt_t)
    d2 = d1 - vol * sqrt_t
    return s0 * norm_cdf(d1) - k * math.exp(-rate * t) * norm_cdf(d2), d1, d2


bs_price, d1, d2 = bs_call()

# GBM exact solution: E[S_T] = S_0 e^{mu T}
B_T = rng.normal(0.0, math.sqrt(T), size=paths)
S_T_mu = S0 * np.exp((mu - 0.5 * sigma**2) * T + sigma * B_T)
print("GBM: E[S_T] = S_0 e^{mu T}:",
      np.isclose(S_T_mu.mean(), S0 * math.exp(mu * T), rtol=0.02))

# risk-neutral Monte Carlo vs BS formula
W_T = rng.normal(0.0, math.sqrt(T), size=paths)
S_T_Q = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * W_T)
mc_Q = math.exp(-r * T) * np.mean(np.maximum(S_T_Q - K, 0.0))
print("risk-neutral MC matches BS call price:",
      np.isclose(mc_Q, bs_price, rtol=0.02))

# real-measure MC depends on mu and does not match BS
S_T_P = S0 * np.exp((mu - 0.5 * sigma**2) * T + sigma * W_T)
mc_P = math.exp(-r * T) * np.mean(np.maximum(S_T_P - K, 0.0))
print("real-measure MC (drift mu) differs from BS:",
      abs(mc_P - bs_price) > 0.5)

# Girsanov reweighting: paths under P with weight Z_T recover BS
theta = (mu - r) / sigma
Z = np.exp(-theta * W_T - 0.5 * theta**2 * T)
# under P, S_T uses B_T = W_T here (same normals); weight converts to Q
S_T_P2 = S0 * np.exp((mu - 0.5 * sigma**2) * T + sigma * W_T)
mc_girsanov = math.exp(-r * T) * np.mean(Z * np.maximum(S_T_P2 - K, 0.0))
print("Girsanov-weighted real paths recover BS price:",
      np.isclose(mc_girsanov, bs_price, rtol=0.02))
print("E[Z_T] = 1:",
      np.isclose(Z.mean(), 1.0, atol=0.02))

# mu-independence of risk-neutral pricing: two real drifts, same Q price
mu2 = -0.05
theta2 = (mu2 - r) / sigma
Z2 = np.exp(-theta2 * W_T - 0.5 * theta2**2 * T)
S_T_P3 = S0 * np.exp((mu2 - 0.5 * sigma**2) * T + sigma * W_T)
mc_g2 = math.exp(-r * T) * np.mean(Z2 * np.maximum(S_T_P3 - K, 0.0))
print("Girsanov price independent of real drift mu:",
      np.isclose(mc_girsanov, mc_g2, rtol=0.03)
      and np.isclose(mc_g2, bs_price, rtol=0.03))

# call delta = N(d_1): finite difference of the BS formula
eps = 1e-4
c_up, _, _ = bs_call(s0=S0 + eps)
c_dn, _, _ = bs_call(s0=S0 - eps)
delta_fd = (c_up - c_dn) / (2 * eps)
print("call delta = N(d_1) (finite difference of BS formula):",
      np.isclose(delta_fd, norm_cdf(d1), rtol=1e-4))

# discrete delta-hedge replication error shrinks with finer rebalancing.
# Self-financing step under Q: dV = delta dS + r (V - delta S) dt;
# terminal V should approach the call payoff as the grid is refined.
def hedge_rms(n_steps, n_paths=15_000):
    dt = T / n_steps
    dW = rng.normal(0.0, math.sqrt(dt), size=(n_paths, n_steps))
    S = np.full(n_paths, S0)
    V = np.full(n_paths, bs_price)
    for k in range(n_steps):
        tau = T - k * dt
        if tau < 1e-12:
            delta = np.where(S > K, 1.0, 0.0)
        else:
            d1_k = (np.log(S / K) + (r + 0.5 * sigma**2) * tau) / (
                sigma * math.sqrt(tau))
            delta = 0.5 * (1.0 + np.vectorize(math.erf)(
                d1_k / math.sqrt(2.0)))
        dS = r * S * dt + sigma * S * dW[:, k]
        V = V + delta * dS + r * (V - delta * S) * dt
        S = S + dS
    payoff = np.maximum(S - K, 0.0)
    return math.sqrt(np.mean((V - payoff) ** 2))


err_coarse = hedge_rms(10)
err_fine = hedge_rms(200)
print("delta-hedge replication error shrinks with refinement:",
      err_fine < err_coarse / 2 and err_fine < 2.0)
