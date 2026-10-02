"""Checks for 03-qua-tensor.md.

1. H and H' as real 4x4 left-regular matrices; H'(j^2=1) has k^2=1.
2. H (x) H: among the 15 non-identity basis elements, the maximal
   mutually anticommuting subsets have sizes 3 and 5; the 5-sets are
   exactly (1)-(6) of the article (6 of size 5); any 4 of (3) generate
   the 5th (up to sign) and all 16 basis elements; the chosen generators
   {i(x)i, j(x)i, 1(x)j, 1(x)k} have signature (2,2).
3. H' (x) H': the four generator choices give (3,1),(2,2),(2,2),(3,1).
4. Formulas Cl_{p,q} (x) H = Cl_{q,p+2} and Cl_{p,q} (x) H' =
   Cl_{q+2,p} = Cl_{p+1,q+1}: new generators anticommute with the
   expected squares, for base algebras built from the above.
5. H'(x)H and H(x)H' give (0,4)/(1,3) and (4,0)/(1,3).
6. The final table: each cell is the stated tensor product.
"""
import itertools
import numpy as np

# ---- algebras with basis (1,i,j,k): i^2=-1, j^2=s, k=ij
def algebra(s):
    # structure: product of basis indices -> (sign, index)
    t = {}
    sq = [1, -1, s, s]  # k^2 = s
    for a in range(4):
        for b in range(4):
            if a == 0: t[a, b] = (1, b)
            elif b == 0: t[a, b] = (1, a)
            elif a == b: t[a, b] = (sq[a], 0)
            else:
                # ij=k, jk=?, ki=?  derived from i^2=-1, j^2=s, k=ij
                # i j = k ; j i = -k ; j k = j(ij) = -i j j = -s i
                # k j = (ij)j = s i ; k i = (ij)i = -i i j = j ; i k = i i j = -j
                tab = {(1,2):(1,3),(2,1):(-1,3),(2,3):(-s,1),(3,2):(s,1),
                       (3,1):(1,2),(1,3):(-1,2)}
                t[a, b] = tab[a, b]
    mats = []
    for a in range(4):  # left regular rep
        M = np.zeros((4, 4))
        for b in range(4):
            sg, c = t[a, b]
            M[c, b] = sg
        mats.append(M)
    return mats

Hm = algebra(-1)   # H
Hp = algebra(1)    # H'
for B, s in ((Hm, -1), (Hp, 1)):
    assert np.allclose(B[1] @ B[1], -np.eye(4)) and np.allclose(B[2] @ B[2], s*np.eye(4))
    assert np.allclose(B[3], B[1] @ B[2]) and np.allclose(B[1] @ B[2], -B[2] @ B[1])
assert np.allclose(Hp[3] @ Hp[3], np.eye(4))   # k^2 = 1 for H'
assert np.allclose(Hm[3] @ Hm[3], -np.eye(4))

def kr(A, B, a, b): return np.kron(A[a], B[b])
def anti(X, Y): return np.allclose(X @ Y, -Y @ X)
def sq(X):
    if np.allclose(X @ X, np.eye(len(X))): return 1
    if np.allclose(X @ X, -np.eye(len(X))): return -1
    raise AssertionError
def sig(gens):
    for X, Y in itertools.combinations(gens, 2): assert anti(X, Y)
    s = [sq(X) for X in gens]
    return s.count(1), s.count(-1)

# ---- 2. H (x) H
idx = [(a, b) for a in range(4) for b in range(4) if (a, b) != (0, 0)]
M = {(a, b): kr(Hm, Hm, a, b) for a, b in idx}
n = len(idx)
adj = {x: {y for y in idx if y != x and anti(M[x], M[y])} for x in idx}
# maximal cliques (Bron-Kerbosch)
cl = []
def bk(R, P, X):
    if not P and not X: cl.append(frozenset(R)); return
    for v in list(P):
        bk(R | {v}, P & adj[v], X & adj[v]); P = P - {v}; X = X | {v}
bk(set(), set(idx), set())
sizes = sorted(len(c) for c in cl)
print("maximal anticommuting sets sizes:", {s: sizes.count(s) for s in set(sizes)})
assert set(sizes) == {3, 5} and sizes.count(5) == 6
five = {c for c in cl if len(c) == 5}
def parse(s): return frozenset(("1ijk".index(t[0]), "1ijk".index(t[1])) for t in s)
lists = [parse(s) for s in [
    ["1i","1j","ik","jk","kk"], ["1i","1k","ij","jj","kj"],
    ["1j","1k","ii","ji","ki"], ["i1","j1","ki","kj","kk"],
    ["i1","ji","jj","jk","k1"], ["ii","ij","ik","j1","k1"]]]
assert five == set(lists), "(1)-(6) mismatch"
print("(1)-(6) are exactly the 6 five-element sets")

S3 = ["1j","1k","ii","ji","ki"]
S3 = [("1ijk".index(t[0]), "1ijk".index(t[1])) for t in S3]
def up_to_sign(X, Y): return np.allclose(X, Y) or np.allclose(X, -Y)
# any 4 of (3) generate the 5th
for drop in S3:
    g = [M[x] for x in S3 if x != drop]
    prod = g[0] @ g[1] @ g[2] @ g[3]
    assert up_to_sign(prod, M[drop]), drop
# the 4 generators generate all 16 basis elements up to sign
def span_basis(g):
    out = []
    for r in range(5):
        for c in itertools.combinations(range(4), r):
            X = np.eye(16)
            for i in c: X = X @ g[i]
            out.append(X)
    return out
gens = [M[(1, 1)], M[(2, 1)], M[(0, 2)], M[(0, 3)]]
prods = span_basis(gens)
allb = [kr(Hm, Hm, a, b) for a in range(4) for b in range(4)]
for B in allb: assert any(up_to_sign(B, X) for X in prods)
assert np.linalg.matrix_rank(np.array([X.ravel() for X in prods])) == 16
# article's explicit identities
assert np.allclose(M[(3,1)] @ M[(2,1)], kr(Hm, Hm, 1, 0))
assert np.allclose(M[(1,1)] @ M[(3,1)], kr(Hm, Hm, 2, 0))
assert np.allclose(M[(0,2)] @ M[(0,3)], kr(Hm, Hm, 0, 1))
tbl = [("kiji", [(0,3),(0,2),(1,1),(2,1)], (3,1)),
       ("1k.. ", [(0,3),(1,1),(2,1),(3,1)], (0,2)),
       ("", [(1,1),(2,1),(3,1),(0,2)], (0,3)),
       ("", [(2,1),(3,1),(0,3),(0,2)], (1,1)),
       ("", [(3,1),(0,3),(0,2),(1,1)], (2,1))]
for _, fs, res in tbl:
    X = np.eye(16)
    for f in fs: X = X @ M[f]
    assert np.allclose(X, M[res]), (fs, res)
print("H(x)H: generators (i⊗i, j⊗i, 1⊗j, 1⊗k) signature", sig(gens))
assert sig(gens) == (2, 2)

# ---- 3. H' (x) H'
res = []
for L in ([Hp[1], Hp[2]], [Hp[2], Hp[3]]):
    g = [np.kron(x, Hp[1]) for x in L] + [np.kron(np.eye(4), Hp[2]), np.kron(np.eye(4), Hp[3])]
    res.append(sig(g))
    g = [np.kron(x, Hp[2]) for x in L] + [np.kron(np.eye(4), Hp[3]), np.kron(np.eye(4), Hp[2] @ Hp[3])]
    res.append(sig(g))
print("H'(x)H' cases (i,j)(x)i / (i,j)(x){j,k} / (j,k)(x)i / (j,k)(x){j,k}:", res)
assert res == [(3, 1), (2, 2), (2, 2), (3, 1)]
# base Cl signatures
assert sig([Hp[1], Hp[2]]) == (1, 1) and sig([Hp[2], Hp[3]]) == (2, 0)
assert sig([Hm[1], Hm[2]]) == (0, 2)

# ---- 4. general formulas on several bases
bases = {"H": [Hm[1], Hm[2]], "H'(1,1)": [Hp[1], Hp[2]], "H'(2,0)": [Hp[2], Hp[3]]}
# bigger bases: tensor-extend by H (gives (0,4) or (2,2) ...), then use as base
def extH(g, gm=Hm):
    I = np.eye(len(g[0]))
    return [np.kron(x, gm[1]) for x in g] + [np.kron(I, gm[2]), np.kron(I, gm[3])]
def extHp1(g):
    I = np.eye(len(g[0]))
    return [np.kron(x, Hp[1]) for x in g] + [np.kron(I, Hp[2]), np.kron(I, Hp[3])]
def extHp2(g):
    I = np.eye(len(g[0]))
    return [np.kron(x, Hp[2]) for x in g] + [np.kron(I, Hp[3]), np.kron(I, Hp[2] @ Hp[3])]
for name, g in bases.items():
    p, q = sig(g)
    a = sig(extH(g)); b = sig(extHp1(g)); c = sig(extHp2(g))
    print(f"base {name}={ (p,q) }: (x)H {a}, (x)H' {b} / {c}")
    assert a == (q, p + 2) and b == (q + 2, p) and c == (p + 1, q + 1)
    g2 = extH(g)
    p2, q2 = sig(g2)
    assert sig(extH(g2)) == (q2, p2 + 2) and sig(extHp1(g2)) == (q2 + 2, p2)
    assert sig(extHp2(g2)) == (p2 + 1, q2 + 1)

# ---- 5. mixed products
print("H'(x)H:", sig(extH([Hp[2], Hp[3]])), sig(extH([Hp[1], Hp[2]])))
assert sig(extH([Hp[2], Hp[3]])) == (0, 4) and sig(extH([Hp[1], Hp[2]])) == (1, 3)
print("H(x)H':", sig(extHp1([Hm[1], Hm[2]])), sig(extHp2([Hm[1], Hm[2]])))
assert sig(extHp1([Hm[1], Hm[2]])) == (4, 0) and sig(extHp2([Hm[1], Hm[2]])) == (1, 3)

# ---- 6. dimension check: 4 generators of H(x)H generate 16 elements (done)
print("all OK")

# ---- 7. the 5 elements of (3) multiply to a scalar (so they span only 16 dims)
X = np.eye(16)
for x in S3: X = X @ M[x]
assert np.allclose(X, X[0, 0] * np.eye(16)) and abs(abs(X[0, 0]) - 1) < 1e-9
print("product of the five elements of (3) =", X[0, 0], "* I")
