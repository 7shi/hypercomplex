import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "oct"))
from octonion import L, R  # 8x8 matrices: L[i]x = e_i x, R[i]x = x e_i (L[0]=R[0]=I)

rng = np.random.default_rng(0)


def omul(x, y):
    """Octonion product of 8-component vectors."""
    return sum(x[i] * L[i] for i in range(8)) @ y


def conj(x):
    c = x.copy()
    c[1:] *= -1
    return c


def unit_imag(v7=None):
    """Random (or given) unit pure-imaginary octonion as 8-vector."""
    if v7 is None:
        v7 = rng.normal(size=7)
    x = np.zeros(8)
    x[1:] = v7 / np.linalg.norm(v7)
    return x


def unit_oct():
    x = rng.normal(size=8)
    return x / np.linalg.norm(x)


def exp_imag(theta, n):
    """exp(theta*n) for unit pure-imaginary n."""
    return np.cos(theta) * np.eye(8)[0] + np.sin(theta) * n


def map_matrix(f):
    """8x8 matrix of a linear map on octonions."""
    return np.array([f(np.eye(8)[j]) for j in range(8)]).T


def ok(label, cond):
    print(f"  {'OK ' if cond else 'NG '} {label}")
    assert cond, label


print("=== 7D: 純虚八元数 ===")
m = unit_imag()
# n at angle phi from m
phi = 0.7
t = unit_imag()
t -= (t @ m) * m
t /= np.linalg.norm(t)
n = np.cos(phi) * m + np.sin(phi) * t
v = np.zeros(8)
v[1:] = rng.normal(size=7)

# 鏡映 nvn = v - 2(v.n)n、線対称な鏡映は -nvn
nvn1 = omul(omul(n, v), n)
nvn2 = omul(n, omul(v, n))
ok("(nv)n = n(vn)  [交代性]", np.allclose(nvn1, nvn2))
ok("nvn = v - 2(v.n)n  [鏡映]", np.allclose(nvn1, v - 2 * (v @ n) * n))

# 2回の鏡映 n(mvm)n
def double_refl(v):
    return omul(omul(n, omul(omul(m, v), m)), n)

# span(m,n) に直交する成分は不変
w = np.zeros(8)
w[1:] = rng.normal(size=7)
w -= (w @ m) * m + (w @ t) * t
ok("n(mvm)n は span(m,n) の直交補空間で恒等", np.allclose(double_refl(w), w))

# span(m,n) 上では 2*phi の回転（m から n の向き）
Dm, Dt = double_refl(m), double_refl(t)
block = np.array([[m @ Dm, m @ Dt], [t @ Dm, t @ Dt]])
expected = np.array([[np.cos(2 * phi), -np.sin(2 * phi)], [np.sin(2 * phi), np.cos(2 * phi)]])
ok("n(mvm)n は span(m,n) 上で角度 2φ の回転", np.allclose(block, expected))

# 括弧を外すと別の変換になる
r = omul(n, m)  # nm
rc = conj(r)
ok("(nm)* = mn", np.allclose(rc, omul(m, n)))
sandwich = omul(omul(r, v), rc)  # (nm)v(mn)
ok("n(mvm)n ≠ (nm)v(mn)  [非結合性]", not np.allclose(double_refl(v), sandwich))

# (nm)v(mn) = rvr* は三重回転: 軸 u = Im(nm)/|Im(nm)| を固定し、
# 直交6次元は3平面に分かれ、すべて同じ角度で回転（等傾）
u = r.copy()
u[0] = 0
u /= np.linalg.norm(u)
S = map_matrix(lambda x: omul(omul(r, x), rc))
ok("rvr* は実部と軸 u を固定", np.allclose(S @ np.eye(8)[0], np.eye(8)[0]) and np.allclose(S @ u, u))
# u と実軸の直交補空間（6次元）に制限して M + M^T = 2cosψ I を確認
basis = []
for j in range(1, 8):
    b = np.eye(8)[j] - (np.eye(8)[j] @ u) * u
    for c in basis:
        b -= (b @ c) * c
    if np.linalg.norm(b) > 1e-9:
        basis.append(b / np.linalg.norm(b))
B = np.array(basis).T  # 8x6
M6 = B.T @ S @ B
sym = M6 + M6.T
cospsi = sym[0, 0] / 2
ok("rvr* は直交6次元で等傾（3平面が同角度）", np.allclose(sym, 2 * cospsi * np.eye(6)))
print(f"     回転角 ψ = {np.degrees(np.arccos(cospsi)):.4f}°, 2φ = {np.degrees(2*phi):.4f}°")

# 四元数部分代数 span(m, t, u) 上では両者は一致
for vq in (m, t, u):
    assert np.allclose(double_refl(vq), omul(omul(r, vq), rc))
print("  OK  span(m,n,m×n) 上では n(mvm)n = (nm)v(mn)  [四元数部分代数]")

print()
print("=== 7D: 三重回転と自己同型（G2） ===")
for theta, expect in [(0.7, False), (np.pi / 3, False), (np.pi, False),
                      (2 * np.pi / 3, True), (4 * np.pi / 3, True)]:
    nn = unit_imag()
    rr = exp_imag(theta / 2, nn)  # 回転角 theta の三重回転
    x, y = rng.normal(size=8), rng.normal(size=8)
    fx = omul(omul(rr, x), conj(rr))
    fy = omul(omul(rr, y), conj(rr))
    fxy = omul(omul(rr, omul(x, y)), conj(rr))
    is_auto = np.allclose(fxy, omul(fx, fy))
    r3 = omul(omul(rr, rr), rr)
    r3_real = np.allclose(r3[1:], 0)
    ok(f"θ={theta:.4f}: 自己同型={is_auto} (期待={expect}), r³∈ℝ: {r3_real}",
       is_auto == expect and r3_real == expect)

# Der(O) の次元 = 14 (= dim G2)
rows = []
for i in range(8):
    for j in range(8):
        w = L[i][:, j].astype(float)  # e_i e_j
        # D(e_i e_j) - (D e_i)e_j - e_i(D e_j) = 0 を D の成分について整理
        A = np.zeros((8, 64))
        for k in range(8):
            for l in range(8):
                A[k, k * 8 + l] += w[l]
        for k in range(8):
            for l in range(8):
                A[k, l * 8 + i] -= R[j][k, l]
                A[k, l * 8 + j] -= L[i][k, l]
        rows.append(A)
A = np.vstack(rows)
nullity = 64 - np.linalg.matrix_rank(A)
ok(f"dim Der(O) = {nullity} (= 14, G2)", nullity == 14)

print()
print("=== リー代数の閉包次元 ===")

def lie_closure(gens):
    elems = [g.astype(float) for g in gens]
    rank = np.linalg.matrix_rank(np.array([e.flatten() for e in elems]))
    while True:
        brackets = [elems[i] @ elems[j] - elems[j] @ elems[i]
                    for i in range(len(elems)) for j in range(i + 1, len(elems))]
        allmat = np.array([e.flatten() for e in elems + brackets])
        r2 = np.linalg.matrix_rank(allmat)
        if r2 == rank:
            return rank, np.array([e.flatten() for e in elems])
        _, _, Vt = np.linalg.svd(allmat)
        elems = [Vt[k].reshape(8, 8) for k in range(r2)]
        rank = r2

ad = [L[i] - R[i] for i in range(1, 8)]
d_ad, _ = lie_closure(ad)
d_r, _ = lie_closure([R[i] for i in range(1, 8)])
d_l, _ = lie_closure([L[i] for i in range(1, 8)])
ok(f"閉包 dim: 挟み込み(ad)={d_ad} (=21: SO(7))", d_ad == 21)
ok(f"閉包 dim: 右乗算={d_r} (=28: SO(8))", d_r == 28)
ok(f"閉包 dim: 左乗算={d_l} (=28: SO(8))", d_l == 28)

# 三重回転の合成は一般に三重回転ではない
# （三重回転の固有値は 1 と e^{±iψ}（3重）だが、合成では回転角が分裂する）
def conj_map(r):
    return map_matrix(lambda x: omul(omul(r, x), conj(r)))

M12 = conj_map(exp_imag(0.35, unit_imag())) @ conj_map(exp_imag(0.6, unit_imag()))
angles = np.sort(np.abs(np.angle(np.linalg.eigvals(M12))))
distinct = np.unique(np.round(angles[2:], 8))  # 実部・軸の固有値1を除く6個 → 3組
ok(f"三重回転の合成は三重回転ではない（回転角が {len(distinct)} 種に分裂）", len(distinct) > 1)
print(f"     分裂した回転角: {[f'{np.degrees(a):.4f}°' for a in distinct]}")

print()
print("=== 8D: 一般の八元数 ===")
m8 = unit_oct()
t8 = rng.normal(size=8)
t8 -= (t8 @ m8) * m8
t8 /= np.linalg.norm(t8)
n8 = np.cos(phi) * m8 + np.sin(phi) * t8
o = rng.normal(size=8)

# 鏡映 -n o* n
refl = lambda o, n: -omul(omul(n, conj(o)), n)
ok("-no*n = o - 2(o.n)n  [8次元の鏡映]", np.allclose(refl(o, n8), o - 2 * (o @ n8) * n8))

# 2回の鏡映（n → m の順）
E = lambda o: refl(refl(o, n8), m8)
ok("2回の鏡映 = m(n*on*)m",
   np.allclose(E(o), omul(omul(m8, omul(omul(conj(n8), o), conj(n8))), m8)))
w8 = rng.normal(size=8)
w8 -= (w8 @ m8) * m8 + (w8 @ t8) * t8
ok("span(m,n) の直交補空間（6次元）で恒等", np.allclose(E(w8), w8))
Em, Et = E(m8), E(t8)
block = np.array([[m8 @ Em, m8 @ Et], [t8 @ Em, t8 @ Et]])
ok("span(m,n) 上で角度 2φ の回転", np.allclose(block @ block.T, np.eye(2))
   and np.isclose(block[0, 0], np.cos(2 * phi)))
r8 = omul(m8, conj(n8))
ok("m(n*on*)m ≠ (mn*)o(n*m)  [非結合性]",
   not np.allclose(E(o), omul(omul(r8, o), omul(conj(n8), m8))))

print()
print("=== 8D: 片側からの積は等傾四重回転 ===")
theta = 0.7
nn = unit_imag()
rr = exp_imag(theta, nn)
Mright = map_matrix(lambda x: omul(x, rr))
Mleft = map_matrix(lambda x: omul(rr, x))
for name, M in [("右乗算 or", Mright), ("左乗算 ro", Mleft)]:
    ok(f"{name}: 直交行列", np.allclose(M.T @ M, np.eye(8)))
    ok(f"{name}: M+M^T = 2cosθ I  [全8次元が角度θで回転=等傾四重回転]",
       np.allclose(M + M.T, 2 * np.cos(theta) * np.eye(8)))

# n = e_1 の場合の回転面を確認
e = np.eye(8)
rr1 = exp_imag(theta, e[1])
c, s = np.cos(theta), np.sin(theta)
Mr = map_matrix(lambda x: omul(x, rr1))
Ml = map_matrix(lambda x: omul(rr1, x))
expect_r = np.eye(8)
expect_l = np.eye(8)
for (i, j), sign_r, sign_l in [((0, 1), 1, 1), ((2, 3), -1, 1), ((4, 5), -1, 1), ((7, 6), -1, 1)]:
    for M, sg in [(expect_r, sign_r), (expect_l, sign_l)]:
        M[i, i] = M[j, j] = c
        M[j, i] = sg * s
        M[i, j] = -sg * s
ok("xe^(θe1): (1,e1)面 +θ、(e2,e3),(e4,e5),(e7,e6)面 -θ", np.allclose(Mr, expect_r))
ok("e^(θe1)x: (1,e1),(e2,e3),(e4,e5),(e7,e6)面 すべて +θ", np.allclose(Ml, expect_l))

# 両側から挟む
Mrr = map_matrix(lambda x: omul(omul(rr, x), rr))
ok("rxr: (rx)r = r(xr)  [交代性]",
   np.allclose(Mrr, map_matrix(lambda x: omul(rr, omul(x, rr)))))
P1n = np.outer(e[0], e[0]) + np.outer(nn, nn)
ok("rxr は (1,n)面のみ 2θ 回転（他6次元は恒等）",
   np.allclose(Mrr @ (np.eye(8) - P1n), np.eye(8) - P1n)
   and np.isclose(e[0] @ Mrr @ e[0], np.cos(2 * theta)))
Mconj = map_matrix(lambda x: omul(omul(rr, x), conj(rr)))
sym6 = Mconj + Mconj.T
ok("rxr* は 1 と n を固定し、直交6次元の3面が 2θ 回転（三重回転）",
   np.allclose(Mconj @ e[0], e[0]) and np.allclose(Mconj @ nn, nn)
   and np.allclose((np.eye(8) - P1n) @ sym6 @ (np.eye(8) - P1n),
                   2 * np.cos(2 * theta) * (np.eye(8) - P1n)))

print()
print("=== 8D: 合成は畳めない ===")
a, b = unit_oct(), unit_oct()
ok("(oa)b ≠ o(ab)  [右乗算の合成は1回の乗算に戻らない]",
   not np.allclose(omul(omul(o, a), b), omul(o, omul(a, b))))
# 補足: 右乗算2つの合成は等傾性を保つ（<(xa)b, x> = <a, b*>|x|^2）が、
# 単一の右乗算 R_c なら 1 の行き先から c = ab となるため、上の不一致と矛盾。
# 3つ以上の合成では等傾性も失われる。
Mab = map_matrix(lambda x: omul(omul(x, a), b))
sym = Mab + Mab.T
ok("右乗算2つの合成は等傾回転のまま", np.allclose(sym, sym[0, 0] * np.eye(8)))
c3 = unit_oct()
M3 = map_matrix(lambda x: omul(omul(omul(x, a), b), c3))
sym3 = M3 + M3.T
ok("右乗算3つの合成は等傾回転ではない", not np.allclose(sym3, sym3[0, 0] * np.eye(8)))
A7, B7 = 0.6 * unit_imag(), 0.8 * unit_imag()
expA = exp_imag(0.6, A7 / 0.6)
expB = exp_imag(0.8, B7 / 0.8)
AB = A7 + B7
nAB = np.linalg.norm(AB)
expAB = exp_imag(nAB, AB / nAB)
ok("e^A e^B ≠ e^(A+B)", not np.allclose(omul(expA, expB), expAB))
rL, rR = unit_oct(), unit_oct()
ok("(r_L o)r_R ≠ r_L(o r_R)  [左右から掛ける順序で結果が変わる]",
   not np.allclose(omul(omul(rL, o), rR), omul(rL, omul(o, rR))))

print()
print("all checks passed")
