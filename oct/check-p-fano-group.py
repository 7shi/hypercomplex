import numpy as np
from itertools import product
from octonion import triples, L, P

# the 8 terms of P: identity and -L_a L_b L_c for each Fano line
G = [np.eye(8)] + [-(L[a] @ L[b] @ L[c]) for a,b,c in triples]
names = ["I"] + [f"T{a}{b}{c}" for a,b,c in triples]

print("P = (1/8) sum(G):", np.allclose(sum(G)/8, P))
print("each g^2 = I:", all(np.allclose(g@g, np.eye(8)) for g in G))
print("each g fixes unit (g e0 = e0):", all(np.allclose(g[:,0], np.eye(8)[:,0]) for g in G))
print("gP = P for all g:", all(np.allclose(g@P, P) for g in G))

def find(M):
    for n,g in zip(names,G):
        if np.allclose(M, g): return n
        if np.allclose(M, -g): return "-"+n
    return None

print("closed & abelian (exact signs):",
      all(find(g@h) == find(h@g) and find(g@h) is not None and "-" not in find(g@h)
          for g,h in product(G,G)))
print("example: T123 * T145 =", find(G[1]@G[2]), "(lines 123,145,176 all pass through point 1)")

# absorption: L_i * T_line = L_j L_k  =>  L_i L_j P = L_k P  (Fano product)
print("L1 * T123 = L2 L3:", np.allclose(L[1]@G[1], L[2]@L[3]))
print("L2 L3 P = L1 P:", np.allclose(L[2]@L[3]@P, L[1]@P))

# cosets: {L_k g : g in G} are 8 distinct matrices, all projecting to L_k P
ok = True
for k in range(8):
    coset = [L[k]@g for g in G]
    ok &= all(np.allclose(m@P, L[k]@P) for m in coset)
    ok &= all(not np.allclose(coset[i], coset[j]) for i in range(8) for j in range(i))
print("each coset L_k G: 8 distinct elements, all == L_k P after projection:", ok)
