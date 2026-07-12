"""Can the octonion projection P (03-oct-left-mul.md) be written in the
L_iR_i form used for the quaternion projection in lie/03-spin.md
(P = (I - LiRi - LjRj - LkRk)/4)?

Finding: L_iR_i is diagonal for every i, with -1 at position 0 and
position i and +1 elsewhere (true for both the quaternion and octonion
cases - it does not depend on associativity). Solving for coefficients
a0*I + a*sum(LiRi) = P with n generators gives a = 1/(2-2n), a0 = (2-n)*a.
n=3 (quaternion) makes a0 = -a = 1/4, matching the article's formula.
n=7 (octonion) gives a0 = 5/12, a = -1/12, an unequal split - so the
L_iR_i form still works, just with different (non-equal) coefficients.
"""

import numpy as np
from octonion import L, R, P

n = 7
LiRi = [L[i] @ R[i] for i in range(1, n + 1)]

print("each L_iR_i is diagonal:", all(np.allclose(M, np.diag(np.diag(M))) for M in LiRi))
print("diagonal pattern (-1 at 0 and i, +1 elsewhere):",
      all(np.array_equal(np.diag(M), [-1 if k in (0, i) else 1 for k in range(8)])
          for i, M in enumerate(LiRi, start=1)))

a = 1 / (2 - 2 * n)
a0 = (2 - n) * a
cand = a0 * np.eye(8) + a * sum(LiRi)
print(f"P == {a0:.4g}*I + {a:.4g}*sum(L_iR_i) for i=1..7:", np.allclose(cand, P))
print("as a fraction: P == (5*I - sum(L_iR_i))/12:", np.allclose((5 * np.eye(8) - sum(LiRi)) / 12, P))

# sanity check against the quaternion case (n=3): (I - LiRi - LjRj - LkRk)/4
n3 = 3
a3 = 1 / (2 - 2 * n3)
a03 = (2 - n3) * a3
print(f"formula reproduces quaternion case (n=3 -> a0={a03}, a={a3}):",
      np.isclose(a03, 1/4) and np.isclose(a3, -1/4))
