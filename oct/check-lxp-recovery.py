import numpy as np
from octonion import L, R, P

rng = np.random.default_rng(0)
x = rng.standard_normal(8)
Lx = sum(x[i]*L[i] for i in range(8))
S = Lx @ P

# reconstruction: sum_j R_j S R_j^{-1}  (R_0 = I, R_j^{-1} = -R_j for j>=1)
rec = S.copy()
for j in range(1,8):
    rec += R[j] @ S @ (-R[j])

print("max |rec - Lx| =", np.abs(rec - Lx).max())

# also check each term is exactly column j of Lx
ok = True
for j in range(1,8):
    T = R[j] @ S @ (-R[j])
    col = np.zeros((8,8)); col[:,j] = Lx[:,j]
    ok &= np.allclose(T, col)
print("each conjugate = single column of Lx:", ok)

# check R_j expressed via L's (so everything stays in the algebra <L_i>)
R1 = 0.5*(-L[1] + L[2]@L[3] + L[4]@L[5] + L[7]@L[6])
print("R1 formula matches:", np.allclose(R1, R[1]))
