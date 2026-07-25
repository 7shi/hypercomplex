import numpy as np

# Fano triples (positive): e_a e_b = e_c cyclically
triples = [(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]
mul = {}  # (i,j) -> (sign, k) for imaginary units
for a,b,c in triples:
    for x,y,z in [(a,b,c),(b,c,a),(c,a,b)]:
        mul[(x,y)] = (1,z); mul[(y,x)] = (-1,z)

def basis_mul(i,j):  # returns (sign, k) for e_i e_j, indices 0..7 (0 = unit)
    if i==0: return (1,j)
    if j==0: return (1,i)
    if i==j: return (-1,0)
    s,k = mul[(i,j)]; return (s,k)

def Lmat(i):
    M = np.zeros((8,8), dtype=int)
    for j in range(8):
        s,k = basis_mul(i,j); M[k,j] = s
    return M

def Rmat(i):
    M = np.zeros((8,8), dtype=int)
    for j in range(8):
        s,k = basis_mul(j,i); M[k,j] = s
    return M

L = [Lmat(i) for i in range(8)]   # L[0] = I
R = [Rmat(i) for i in range(8)]   # R[0] = I
P = np.zeros((8,8), dtype=int); P[0,0] = 1
