import numpy as np
from octonion import L, R, P, basis_mul

# For every pair x,y in 1..7, check:
#   L_{xy} = sum_j R_j (L_x L_y P) R_j^{-1}   (R_0 = I, R_j^{-1} = -R_j for j>=1)
# where the target e_x e_y = sign * e_k comes from the Fano-plane algorithm
# (basis_mul), independent of the matrix construction.

ok = True
for x in range(1, 8):
    for y in range(1, 8):
        S = L[x] @ L[y] @ P
        rec = S.copy()
        for j in range(1, 8):
            rec += R[j] @ S @ (-R[j])

        sign, k = basis_mul(x, y)
        target = sign * L[k]

        match = np.array_equal(rec, target)
        ok &= match
        if not match:
            print(f"MISMATCH x={x} y={y}: e_{x}e_{y} = {'+' if sign>0 else '-'}e_{k}")

print("all 49 pairs (x,y in 1..7) match Fano-plane product:", ok)
