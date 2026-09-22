"""Checks for spinor-ideal.md: the Cl_{3,1}(R) structure of M_4(R) built from
L_u R_v (squares, anticommuting sets, volume element, grade table), the
projection P = 1/4 (I - L_i R_i - L_j R_j - L_k R_k) and the recovery formula
L_x = sum_u R_u (L_x P) R_u^{-1}.
"""

import itertools

import numpy as np

BASIS4 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
NAMES = ["1", "i", "j", "k"]


def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return (w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
            w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2)


def Lmat(u):
    return np.array([qmul(u, b) for b in BASIS4], dtype=int).T


def Rmat(u):
    return np.array([qmul(b, u) for b in BASIS4], dtype=int).T


L = [Lmat(u) for u in BASIS4]   # L[0] = I
R = [Rmat(u) for u in BASIS4]   # R[0] = I
I4 = np.eye(4, dtype=int)
P = np.zeros((4, 4), dtype=int)
P[0, 0] = 1

# --- check 1: product rules and commutativity of L and R -------------------

ok = True
for a in range(4):
    for b in range(4):
        uv = qmul(BASIS4[a], BASIS4[b])
        Luv = sum(uv[m] * L[m] for m in range(4))
        vu = qmul(BASIS4[b], BASIS4[a])
        Rvu = sum(vu[m] * R[m] for m in range(4))
        ok &= np.array_equal(L[a] @ L[b], Luv)
        ok &= np.array_equal(R[a] @ R[b], Rvu)
        ok &= np.array_equal(L[a] @ R[b], R[b] @ L[a])
print("L_uL_v=L_{uv}, R_uR_v=R_{vu}, L_uR_v=R_vL_u:", ok)

# --- check 2: the 16 products L_u R_v span M_4(R) --------------------------

names = {}
mats = []
for a in range(4):
    for b in range(4):
        nm = ("I" if a == 0 and b == 0 else
              ("" if a == 0 else f"L_{NAMES[a]}") +
              ("" if b == 0 else f"R_{NAMES[b]}"))
        names[(a, b)] = nm
        mats.append((L[a] @ R[b]).flatten())
rank = np.linalg.matrix_rank(np.array(mats, dtype=float))
print("rank of the 16 L_u R_v in M_4(R) (expect 16):", rank)

# --- check 3: squares, 6 with -I (one-sided) and 9 with +I (two-sided) -----

minus, plus = [], []
for a in range(4):
    for b in range(4):
        if a == 0 and b == 0:
            continue
        S = L[a] @ R[b]
        S2 = S @ S
        if np.array_equal(S2, -I4):
            minus.append(names[(a, b)])
        elif np.array_equal(S2, I4):
            plus.append(names[(a, b)])
        else:
            print("unexpected square:", names[(a, b)])
print("square = -I (expect 6, one-sided):", len(minus), minus)
print("square = +I (expect 9, two-sided):", len(plus), plus)

# --- check 4: maximal mutually anticommuting sets (expect max size 5) ------


def anticommute(A, B):
    return np.array_equal(A @ B, -(B @ A))


elems = [(a, b) for a in range(4) for b in range(4) if not (a == 0 and b == 0)]
best = 0
best_sets = []
for size in range(2, 7):
    found = []
    for combo in itertools.combinations(elems, size):
        ms = [L[a] @ R[b] for a, b in combo]
        if all(anticommute(ms[p], ms[q])
               for p in range(size) for q in range(p + 1, size)):
            found.append(combo)
    if not found:
        break
    best, best_sets = size, found
print("largest mutually anticommuting set size (expect 5):", best)
print("number of such maximal sets:", len(best_sets))
cited = [(1, 1), (2, 1), (3, 1), (0, 2), (0, 3)]  # L_iR_i,L_jR_i,L_kR_i,R_j,R_k
print("cited set L_iR_i,L_jR_i,L_kR_i,R_j,R_k is one of them:",
      any(set(s) == set(cited) for s in best_sets))
print("its squares (expect +,+,+,-,-):",
      ["+" if np.array_equal((L[a] @ R[b]) @ (L[a] @ R[b]), I4) else "-"
       for a, b in cited])

# --- check 5: generators e_1..e_4, metric (+,+,+,-), volume element --------

e = [L[1] @ R[1], L[2] @ R[1], L[3] @ R[1], R[2]]
sq = [int(np.trace(x @ x) // 4) for x in e]
print("e_1..e_4 squares (expect [1,1,1,-1]):", sq)
print("e_m e_n = -e_n e_m (m != n):",
      all(anticommute(e[m], e[n]) for m in range(4) for n in range(m + 1, 4)))
omega = e[0] @ e[1] @ e[2] @ e[3]
print("omega = e_1e_2e_3e_4 == -R_k:", np.array_equal(omega, -R[3]))

# --- check 6: grade table (16 basis elements, up to sign) ------------------

table = {0: ["I"], 1: ["L_iR_i", "L_jR_i", "L_kR_i", "R_j"],
         2: ["L_i", "L_j", "L_k", "L_iR_k", "L_jR_k", "L_kR_k"],
         3: ["L_iR_j", "L_jR_j", "L_kR_j", "R_i"], 4: ["R_k"]}
by_name = {names[(a, b)]: L[a] @ R[b] for a, b in
           [(a, b) for a in range(4) for b in range(4)]}
ok = True
for g in range(5):
    got = []
    for combo in itertools.combinations(range(4), g):
        M = I4
        for idx in combo:
            M = M @ e[idx]
        hit = [nm for nm, B in by_name.items()
               if np.array_equal(M, B) or np.array_equal(M, -B)]
        got.append(hit[0] if len(hit) == 1 else f"?{hit}")
    ok &= sorted(got) == sorted(table[g])
    print(f"grade {g}: {sorted(got)}  (article: {sorted(table[g])})")
print("grade table matches the article:", ok)

# --- check 7: the projection P --------------------------------------------

Pform = (I4 - L[1] @ R[1] - L[2] @ R[2] - L[3] @ R[3])
print("P == 1/4 (I - L_iR_i - L_jR_j - L_kR_k):", np.array_equal(Pform, 4 * P))
print("P^2 == P:", np.array_equal(P @ P, P), " L_1 P == P:",
      np.array_equal(L[0] @ P, P))

rng = np.random.default_rng(0)
ok = True
for _ in range(20):
    x = tuple(rng.standard_normal(4))
    s = np.array(x)
    for u in BASIS4[1:]:
        s = s - np.array(qmul(qmul(u, x), u))
    ok &= np.allclose(s / 4, [x[0], 0, 0, 0])
print("1/4 (x - ixi - jxj - kxk) == x_0 (20 trials):", ok)

# --- check 8: recovery L_x = sum_u R_u (L_x P) R_u^{-1} --------------------

ok_sum = ok_col = True
for _ in range(20):
    x = rng.standard_normal(4)
    Lx = sum(x[m] * L[m] for m in range(4))
    S = Lx @ P
    rec = S.copy()
    for m in range(1, 4):
        T = R[m] @ S @ (-R[m])
        rec = rec + T
        col = np.zeros((4, 4))
        col[:, m] = Lx[:, m]
        ok_col &= np.allclose(T, col)
    ok_sum &= np.allclose(rec, Lx)
print("sum_u R_u (L_x P) R_u^{-1} == L_x (20 trials):", ok_sum)
print("each term is exactly one column of L_x in place:", ok_col)

# --- check 9: the recovery always lands on L_{A(1)}, for arbitrary A -------

ok = True
for _ in range(20):
    A = rng.standard_normal((4, 4))
    S = A @ P
    rec = S + sum(R[m] @ S @ (-R[m]) for m in (1, 2, 3))
    a1 = A[:, 0]
    ok &= np.allclose(S, sum(a1[m] * L[m] for m in range(4)) @ P)
    ok &= np.allclose(rec, sum(a1[m] * L[m] for m in range(4)))
print("AP == L_{A(1)}P and the sum returns L_{A(1)} for arbitrary A:", ok)

# --- check 10: f_1...f_5 = I, and the 5 generators are linearly independent -

f = [L[1] @ R[1], L[2] @ R[1], L[3] @ R[1], R[2], R[3]]
prod = I4
for M in f:
    prod = prod @ M
print("f_1...f_5 == I (so Cl_{3,2} is not faithfully represented):",
      np.array_equal(prod, I4))
print("rank of f_1..f_5 (linearly independent, expect 5):",
      np.linalg.matrix_rank(np.array([x.flatten() for x in f], dtype=float)))

# --- check 11: the octonion side of the comparison -------------------------

try:
    from common.octonion import L as LO, P as PO, basis_mul
except ImportError:  # pragma: no cover
    LO = None

if LO is not None:
    def omul(x, y):
        z = np.zeros(8)
        for m in range(8):
            for n in range(8):
                s, kk = basis_mul(m, n)
                z[kk] += s * x[m] * y[n]
        return z

    ok_proj, n_ne = True, 0
    for _ in range(20):
        x = rng.standard_normal(8)
        y = rng.standard_normal(8)
        Lx = sum(x[m] * LO[m] for m in range(8))
        Ly = sum(y[m] * LO[m] for m in range(8))
        xy = omul(x, y)
        Lxy = sum(xy[m] * LO[m] for m in range(8))
        ok_proj &= np.allclose(Lx @ Ly @ PO, Lxy @ PO)
        if not np.allclose(Lx @ Ly, Lxy):
            n_ne += 1
    print("octonions: L_x L_y P == L_{xy} P (20 trials):", ok_proj)
    print("octonions: L_x L_y != L_{xy} as operators (count/20):", n_ne)

    mono = []
    for r in range(7):
        for combo in itertools.combinations(range(1, 7), r):
            M = np.eye(8)
            for idx in combo:
                M = M @ LO[idx]
            mono.append(M.flatten())
    print("octonions: rank of the 64 monomials in L_1..L_6 (expect 64):",
          np.linalg.matrix_rank(np.array(mono, dtype=float)))

# --- check 12: the column-wise proof (L_x P = x c_1^T, c_1^T R_u^{-1} = c_u^T)

c = [np.eye(4, dtype=int)[:, m] for m in range(4)]
Rinv = [R[0]] + [-R[m] for m in (1, 2, 3)]
ok_row = all(np.array_equal(c[0] @ Rinv[m], c[m]) for m in range(4))
print("c_1^T R_u^{-1} == c_u^T:", ok_row)
ok_out = True
for _ in range(10):
    x = rng.standard_normal(4)
    Lx = sum(x[m] * L[m] for m in range(4))
    ok_out &= np.allclose(Lx @ P, np.outer(x, c[0]))
    for m in range(4):
        ok_out &= np.allclose(R[m] @ Lx @ P @ Rinv[m],
                              np.outer(R[m] @ x, c[m]))
print("L_x P == x c_1^T and R_u(L_xP)R_u^{-1} == (xu) c_u^T:", ok_out)
