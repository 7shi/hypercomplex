"""Numerical checks for double-cover.md.

Axis-angle redundancies behind the ball model of SO(3): R(-n,phi)=R(n,-phi),
R(n,theta)=R(-n,2pi-theta), and R(e_z,pi)=R(-e_z,pi)=diag(-1,-1,1); the unit
quaternion q = cos(theta/2) + u sin(theta/2) acting by conjugation as the
axis-angle rotation R(u,theta), with rho_{-q} = rho_q; the belt-trick
homotopy H_s(theta) = R(n(s,theta), theta) staying a loop for all s
(H_s(0) = H_s(2pi) = I), connecting the forward 2pi loop (s=0) to the
reversed one (s=pi) through small continuous steps; and lift endpoints:
the 2pi rotation loop lifts to an open path 1 -> -1, the 4pi loop lifts to
a closed loop ending at +1, and the lift of H_s ends at -1 for every s.
"""

import numpy as np

# quaternions as 2x2 complex matrices (matches 02-su2-so3.md)
I2 = np.eye(2, dtype=complex)
qi = np.array([[1j, 0], [0, -1j]])
qj = np.array([[0, 1], [-1, 0]], dtype=complex)
qk = np.array([[0, 1j], [1j, 0]])

def Q(a, b, c, d):
    return a*I2 + b*qi + c*qj + d*qk

def coords(x):  # pure quaternion matrix -> (b, c, d)
    return np.array([x[0,0].imag, x[0,1].real, x[0,1].imag])

def rot(n, theta):  # axis-angle rotation matrix R(n, theta) (Rodrigues)
    x, y, z = n
    K = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=float)
    return np.eye(3) + np.sin(theta)*K + (1 - np.cos(theta))*(K @ K)

def quat(n, theta):  # lift q = cos(theta/2) + n sin(theta/2)
    return np.cos(theta/2)*I2 + np.sin(theta/2)*Q(0, *n)

def rot_from_q(q):  # conjugation rho_q as a 3x3 matrix (unit q: q^-1 = q^dagger)
    return np.column_stack([coords(q @ m @ q.conj().T) for m in (qi, qj, qk)])

rng = np.random.default_rng(0)
v = rng.standard_normal(3)
n = v / np.linalg.norm(v)
theta = rng.uniform(0, np.pi)
ez = np.array([0.0, 0.0, 1.0])

# axis-angle redundancies: the warp on the boundary of the ball model
print("R(-n,phi) = R(n,-phi):", np.allclose(rot(-n, theta), rot(n, -theta)))
t2 = rng.uniform(np.pi, 2*np.pi)
print("R(n,theta) = R(-n,2pi-theta):", np.allclose(rot(n, t2), rot(-n, 2*np.pi - t2)))
print("R(e_z,pi) = R(-e_z,pi) = diag(-1,-1,1):",
      np.allclose(rot(ez, np.pi), np.diag([-1, -1, 1])) and
      np.allclose(rot(ez, np.pi), rot(-ez, np.pi)))

# unit quaternion <-> axis-angle: rho_q is R(n, theta), and rho_{-q} = rho_q
q = quat(n, theta)
print("q unit: q^dagger q = I:", np.allclose(q.conj().T @ q, I2))
print("rho_q = R(n,theta):", np.allclose(rot_from_q(q), rot(n, theta)))
print("rho_{-q} = rho_q:", np.allclose(rot_from_q(-q), rot_from_q(q)))

# belt-trick homotopy H_s(theta) = R(n(s,theta), theta)
def axis(s, th):
    return np.array([np.sin(s)*np.cos(th), np.sin(s)*np.sin(th), np.cos(s)])

ss = np.linspace(0, np.pi, 61)
ths = np.linspace(0, 2*np.pi, 121)
print("H_s(0) = H_s(2pi) = I for all s:",
      all(np.allclose(rot(axis(s, th), th), np.eye(3)) for s in ss for th in (0, 2*np.pi)))
print("H_0(theta) = R(e_z,theta):",
      all(np.allclose(rot(axis(0, th), th), rot(ez, th)) for th in ths))
print("H_pi(theta) = R(e_z,2pi-theta):",
      all(np.allclose(rot(axis(np.pi, th), th), rot(ez, 2*np.pi - th)) for th in ths))
grid = np.array([[rot(axis(s, th), th) for th in ths] for s in ss])
steps = max(np.abs(grid[1:] - grid[:-1]).max(), np.abs(grid[:, 1:] - grid[:, :-1]).max())
print("homotopy moves in small continuous steps:", steps < 0.2)

# lift endpoints: 2pi loop is open (ends at -1), 4pi loop closes (ends at +1)
print("lift of 2pi loop ends at -1:", np.allclose(quat(ez, 2*np.pi), -I2))
print("lift of 4pi loop ends at +1:", np.allclose(quat(ez, 4*np.pi), I2))
lifts = np.array([[quat(axis(s, th), th) for th in ths] for s in ss])
lsteps = max(np.abs(lifts[1:] - lifts[:-1]).max(), np.abs(lifts[:, 1:] - lifts[:, :-1]).max())
print("lift of H_s is continuous and projects to H_s:",
      lsteps < 0.2 and
      all(np.allclose(rot_from_q(lifts[i, j]), grid[i, j])
          for i in range(0, len(ss), 10) for j in range(0, len(ths), 10)))
print("lift endpoint stays -1 for all s:",
      all(np.allclose(quat(axis(s, 2*np.pi), 2*np.pi), -I2) for s in ss))
