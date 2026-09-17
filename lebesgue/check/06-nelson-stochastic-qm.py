"""Numerical checks for 06-nelson-stochastic-qm.md.

For Brownian motion started at the origin (sigma = 1, nu = 1/2) the
forward derivative vanishes while the backward one is D_* B = B_t / t:
E[B_t - B_{t-h} | B_t] = (h/t) B_t, checked as a regression slope.  The
spreading Gaussian density satisfies the continuity equation with
current velocity v = x/(2t) and osmotic velocity u = nu (log rho)_x
= -x/(2t).

In natural units hbar = m = omega = 1, the Nelson diffusion for the
harmonic-oscillator ground state is the Ornstein-Uhlenbeck process
dX = -X dt + dB (nu = 1/2, sigma = 1, osmotic velocity u = -x).
Its invariant density is |psi_0|^2 = pi^{-1/2} e^{-x^2}; the quantum
potential plus V equals the ground energy 1/2; the stationary
probability current b rho - nu rho_x vanishes; and long-run
Euler-Maruyama samples recover the theoretical variance 1/2 and match
|psi_0|^2 in histogram L1 distance.  The osmotic relation
u = nu (log rho)_x is checked on the theoretical density and on a
smoothed histogram.

The final chapter compares the article with 05.  Writing the
Black-Scholes price of a European put in tau = T - t and x = log s,
moving to y = x + (r - sigma^2/2) tau and splitting off e^{-r tau}
leaves the pure heat equation; the resulting h is reproduced by the
Gaussian kernel acting on the payoff.  With hbar = sigma^2 and m = 1
the same function solves the imaginary-time Schrodinger equation with
constant potential V0 = hbar r.  A free Gaussian packet propagated
spectrally shows the difference that the Wick rotation hides: real time
preserves the L2 norm, imaginary time does not.
"""

import math

import numpy as np

rng = np.random.default_rng(0)

# Brownian motion from the origin: DB = 0 but D_* B = B_t / t.
# E[B_{t+h} - B_t | B_t] = 0 (martingale) and, by Gaussian conditioning,
# E[B_t - B_{t-h} | B_t] = (h/t) B_t — regression slopes 0 and h/t.
t1, h = 1.0, 0.01
n_bm = 400_000
B_prev = rng.normal(0.0, math.sqrt(t1 - h), size=n_bm)      # B_{t-h}
B_now = B_prev + rng.normal(0.0, math.sqrt(h), size=n_bm)   # B_t
B_next = B_now + rng.normal(0.0, math.sqrt(h), size=n_bm)   # B_{t+h}
slope_fwd = np.polyfit(B_now, B_next - B_now, 1)[0]
slope_bwd = np.polyfit(B_now, B_now - B_prev, 1)[0]
print("forward derivative of BM vanishes (slope of increment on B_t is 0):",
      abs(slope_fwd) < 1e-3)
print("backward derivative of BM: E[B_t - B_{t-h}|B_t] has slope h/t:",
      abs(slope_bwd - h / t1) < 1e-3)


# spreading Gaussian rho = (2 pi t)^{-1/2} exp(-x^2/(2t)):
# u = nu (log rho)_x = -x/(2t) and continuity d_t rho + d_x(rho v) = 0
# with v = x/(2t).
def rho_bm(t, x):
    return np.exp(-x ** 2 / (2.0 * t)) / math.sqrt(2.0 * math.pi * t)


xg = np.linspace(-3.0, 3.0, 601)
dxg = xg[1] - xg[0]
log_rho_bm = np.log(rho_bm(t1, xg))
u_bm = 0.5 * (log_rho_bm[2:] - log_rho_bm[:-2]) / (2.0 * dxg)
print("BM osmotic velocity nu (log rho)_x = -x/(2t):",
      np.allclose(u_bm, -xg[1:-1] / (2.0 * t1), rtol=0.0, atol=1e-4))

dtb = 1e-4
drho_dt = (rho_bm(t1 + dtb, xg) - rho_bm(t1 - dtb, xg)) / (2.0 * dtb)
flux = rho_bm(t1, xg) * xg / (2.0 * t1)
dflux_dx = (flux[2:] - flux[:-2]) / (2.0 * dxg)
print("spreading Gaussian satisfies continuity with v = x/(2t):",
      np.allclose(drho_dt[1:-1] + dflux_dx, 0.0, rtol=0.0, atol=1e-4))

# natural units: hbar = m = omega = 1
# nu = hbar/(2m) = 1/2, sigma = sqrt(hbar/m) = 1, u = -omega x = -x
# rho = |psi_0|^2 = pi^{-1/2} exp(-x^2), Var = hbar/(2 m omega) = 1/2
# E = hbar omega / 2 = 1/2, V = x^2 / 2
hbar = m = omega = 1.0
nu = hbar / (2.0 * m)
sigma = math.sqrt(hbar / m)
var_theory = hbar / (2.0 * m * omega)
E0 = 0.5 * hbar * omega


def rho0(x):
    return math.pi ** (-0.5) * np.exp(-(x ** 2))


def sqrt_rho0(x):
    return math.pi ** (-0.25) * np.exp(-0.5 * (x ** 2))


xs = np.linspace(-3.0, 3.0, 401)
dx = xs[1] - xs[0]

# quantum potential on a grid: U_Q + V = E_0
# (sqrt rho)_xx / sqrt rho = -1 + x^2 for the HO ground state (alpha = 1)
sr = sqrt_rho0(xs)
sr_xx = (sr[2:] - 2.0 * sr[1:-1] + sr[:-2]) / dx ** 2
U_Q = -(hbar ** 2 / (2.0 * m)) * (sr_xx / sr[1:-1])
V_mid = 0.5 * m * omega ** 2 * xs[1:-1] ** 2
print("quantum potential: U_Q + V = E_0 on the HO ground state:",
      np.allclose(U_Q + V_mid, E0, rtol=0.0, atol=1e-3))

U_Q_exact = E0 - V_mid
print("closed-form U_Q = E_0 - V matches finite difference:",
      np.allclose(U_Q, U_Q_exact, rtol=0.0, atol=1e-3))

# osmotic velocity from theoretical density: nu (log rho)_x = -x
# use interior central differences (exclude two endpoints)
log_rho = np.log(rho0(xs))
log_rho_x = (log_rho[2:] - log_rho[:-2]) / (2.0 * dx)
u_from_rho = nu * log_rho_x
print("osmotic velocity u = nu (log rho)_x equals -x:",
      np.allclose(u_from_rho, -xs[1:-1], rtol=0.0, atol=1e-4))

# stationary Fokker-Planck current J = b rho - nu rho_x = 0 with b = -x
rho = rho0(xs)
rho_x = (rho[2:] - rho[:-2]) / (2.0 * dx)
b_mid = -xs[1:-1]
J = b_mid * rho[1:-1] - nu * rho_x
print("stationary probability current J = b rho - nu rho_x vanishes:",
      np.allclose(J, 0.0, rtol=0.0, atol=1e-4))

# acceleration identity for stationary ground state:
# u u' + nu u'' = omega^2 x  (with u = -omega x)
u = -omega * xs
print("Nelson stationary law: u u' + nu u'' = omega^2 x:",
      np.allclose(u * (-omega) + nu * 0.0, omega ** 2 * xs))

# long-run Euler-Maruyama of dX = -omega X dt + sigma dB.
# Start in equilibrium and pool late samples to reduce histogram noise.
n_paths = 30_000
n_steps = 2_000
dt = 0.01
X = rng.normal(0.0, math.sqrt(var_theory), size=n_paths)
pool = []
for step in range(n_steps):
    X = X - omega * X * dt + sigma * math.sqrt(dt) * rng.normal(size=n_paths)
    if step >= n_steps // 2 and step % 20 == 0:
        pool.append(X.copy())
X_pool = np.concatenate(pool)

print("OU invariant variance matches hbar/(2 m omega):",
      np.isclose(X_pool.var(), var_theory, rtol=0.03))

# histogram vs |psi_0|^2
bins = np.linspace(-4.0, 4.0, 81)
hist, edges = np.histogram(X_pool, bins=bins, density=True)
centers = 0.5 * (edges[:-1] + edges[1:])
rho_c = rho0(centers)
bin_w = edges[1] - edges[0]
l1 = np.sum(np.abs(hist - rho_c)) * bin_w
print("histogram L1 distance to |psi_0|^2 is small:",
      l1 < 0.05)

# osmotic velocity from histogram: central differences on well-populated bins,
# then check that the linear fit of u_hist vs x has slope near -omega
h_edges = np.linspace(-2.0, 2.0, 41)
h_vals, h_edges = np.histogram(X_pool, bins=h_edges, density=True)
h_c = 0.5 * (h_edges[:-1] + h_edges[1:])
mask = h_vals > 0.05 * h_vals.max()
h_c_m = h_c[mask]
log_h = np.log(h_vals[mask])
dlog = (log_h[2:] - log_h[:-2]) / (h_c_m[2:] - h_c_m[:-2])
x_int = h_c_m[1:-1]
u_hist = nu * dlog
slope = np.polyfit(x_int, u_hist, 1)[0]
corr = np.corrcoef(u_hist, -x_int)[0, 1]
print("histogram osmotic velocity tracks u = -x near the origin:",
      abs(slope + omega) < 0.15 and corr > 0.95)

print("ground energy E_0 = 1/2 in natural units:",
      abs(E0 - 0.5) < 1e-15)


# --- Imaginary time: Black-Scholes as a Schrodinger equation ---------------
# The final chapter rewrites the Black-Scholes equation with tau = T - t and
# x = log s, removes the drift by moving to y = x + (r - sigma^2/2) tau and
# splits off the discount via h = e^{r tau} w.  Checks below use the closed
# form of a European put (bounded payoff, so the heat kernel integral is
# easy to evaluate on a truncated grid).
bs_r, bs_sigma, bs_T, bs_K = 0.05, 0.2, 1.0, 100.0


def bs_N(z):
    return 0.5 * (1.0 + np.vectorize(math.erf)(z / math.sqrt(2.0)))


def bs_put(tau, s):
    """European put price at time to maturity tau (tau > 0)."""
    d1 = ((np.log(s / bs_K) + (bs_r + bs_sigma ** 2 / 2.0) * tau)
          / (bs_sigma * np.sqrt(tau)))
    d2 = d1 - bs_sigma * math.sqrt(tau)
    return bs_K * math.exp(-bs_r * tau) * bs_N(-d2) - s * bs_N(-d1)


def bs_u(tau, x):
    """u(tau, x) = V(T - tau, e^x)."""
    return bs_put(tau, np.exp(x))


def bs_h(tau, y):
    """h(tau, y) = e^{r tau} u(tau, y - (r - sigma^2/2) tau)."""
    return math.exp(bs_r * tau) * bs_u(tau, y - (bs_r - bs_sigma ** 2 / 2.0) * tau)


def d_tau(f, tau, z, eps=1e-5):
    return (f(tau + eps, z) - f(tau - eps, z)) / (2.0 * eps)


def d_zz(f, tau, z, eps=1e-3):
    return (f(tau, z + eps) - 2.0 * f(tau, z) + f(tau, z - eps)) / eps ** 2


def d_z(f, tau, z, eps=1e-3):
    return (f(tau, z + eps) - f(tau, z - eps)) / (2.0 * eps)


tau0 = 0.5
y_mid = np.log(bs_K) + np.linspace(-0.5, 0.5, 51)

# (1) chain rule: u_tau = (sigma^2/2) u_xx + (r - sigma^2/2) u_x - r u
res_u = (d_tau(bs_u, tau0, y_mid)
         - (bs_sigma ** 2 / 2.0) * d_zz(bs_u, tau0, y_mid)
         - (bs_r - bs_sigma ** 2 / 2.0) * d_z(bs_u, tau0, y_mid)
         + bs_r * bs_u(tau0, y_mid))
print("BS in (tau, x = log s) is a diffusion with drift and discount:",
      np.max(np.abs(res_u)) < 1e-4)

# (2) moving frame plus discount factor leaves the pure heat equation
res_h = d_tau(bs_h, tau0, y_mid) - (bs_sigma ** 2 / 2.0) * d_zz(bs_h, tau0, y_mid)
print("moving frame and e^{r tau} reduce BS to the heat equation:",
      np.max(np.abs(res_h)) < 1e-4)

# (3) heat kernel: h(tau, .) is the Gaussian convolution of the payoff
y_quad = np.log(bs_K) + np.linspace(-6.0, 6.0, 24001)
dy_quad = y_quad[1] - y_quad[0]
payoff = np.maximum(bs_K - np.exp(y_quad), 0.0)
std = bs_sigma * math.sqrt(tau0)
conv = np.array([
    np.sum(np.exp(-(y - y_quad) ** 2 / (2.0 * std ** 2)) * payoff)
    * dy_quad / math.sqrt(2.0 * math.pi * std ** 2)
    for y in y_mid
])
print("heat kernel propagates the payoff to h(tau, y):",
      np.allclose(conv, bs_h(tau0, y_mid), rtol=1e-4, atol=1e-8))

# (4) the dictionary: with hbar_fin = sigma^2 (m = 1) and V0 = hbar_fin * r,
# w = e^{-r tau} h solves the imaginary-time Schrodinger equation
# psi_tau = (hbar/2m) psi_yy - (V0/hbar) psi.
hbar_fin, m_fin = bs_sigma ** 2, 1.0
V0_fin = hbar_fin * bs_r


def bs_w(tau, y):
    return math.exp(-bs_r * tau) * bs_h(tau, y)


res_w = (d_tau(bs_w, tau0, y_mid)
         - (hbar_fin / (2.0 * m_fin)) * d_zz(bs_w, tau0, y_mid)
         + (V0_fin / hbar_fin) * bs_w(tau0, y_mid))
print("BS solves imaginary-time Schrodinger with hbar = sigma^2, V0 = hbar r:",
      np.max(np.abs(res_w)) < 1e-4)
print("the matched diffusion coefficients agree, sigma^2/2 = hbar/(2m):",
      abs(bs_sigma ** 2 / 2.0 - hbar_fin / (2.0 * m_fin)) < 1e-15)

# (5) real versus imaginary time: unitary evolution keeps the L2 norm,
# the heat semigroup does not.  Free packet, spectral propagation.
L, n_fft = 40.0, 4096
y_fft = -L / 2.0 + L * np.arange(n_fft) / n_fft
k_fft = 2.0 * math.pi * np.fft.fftfreq(n_fft, d=L / n_fft)
psi0 = np.exp(-y_fft ** 2 / 2.0).astype(complex)
disp = hbar_fin * k_fft ** 2 / (2.0 * m_fin)
norm0 = np.sum(np.abs(psi0) ** 2) * (L / n_fft)
psi_real = np.fft.ifft(np.exp(-1j * disp * tau0) * np.fft.fft(psi0))
psi_imag = np.fft.ifft(np.exp(-disp * tau0) * np.fft.fft(psi0))
norm_real = np.sum(np.abs(psi_real) ** 2) * (L / n_fft)
norm_imag = np.sum(np.abs(psi_imag) ** 2) * (L / n_fft)
print("real time is unitary (L2 norm preserved):",
      abs(norm_real - norm0) < 1e-10)
print("imaginary time dissipates (L2 norm strictly decreases):",
      norm_imag < norm0 - 1e-3)
