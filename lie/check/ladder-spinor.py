"""ladder-spinor.md の数式の数値検証。

- 昇降演算子 σ± の行列表示と交換関係
- σ± の射影とベクトルへの分解（σ+ = Pσ1 = σ1Q など）と σ±² = 0
- ヌルベクトル（等方ベクトル）と極大等方部分空間
- 純粋スピノルとしての P（消滅条件 σ+P = 0）
- Cl_{2m}(C) のフォック構成（生成消滅演算子、真空、次元 2^m）
- 消滅する等方部分空間の次元による純粋性の判定（m ≤ 3 は全スピノルが純粋、m = 4 で破れる）
"""

import numpy as np

I2 = np.eye(2, dtype=complex)
s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
s3 = np.array([[1, 0], [0, -1]], dtype=complex)


def eq(a, b):
    return np.allclose(a, b)


def comm(a, b):
    return a @ b - b @ a


def anticomm(a, b):
    return a @ b + b @ a


print("=== 昇降演算子の基本 ===")
sp = (s1 + 1j * s2) / 2
sm = (s1 - 1j * s2) / 2
print("σ+ =\n", sp.real.astype(int))
print("σ- =\n", sm.real.astype(int))
assert eq(sp, [[0, 1], [0, 0]])
assert eq(sm, [[0, 0], [1, 0]])
assert eq(sp.conj().T, sm), "σ+† = σ-"
print("σ+† = σ-:", eq(sp.conj().T, sm))

print("[σ3, σ±] = ±2σ±:", eq(comm(s3, sp), 2 * sp), eq(comm(s3, sm), -2 * sm))
print("[σ+, σ-] = σ3:", eq(comm(sp, sm), s3))
assert eq(comm(s3, sp), 2 * sp) and eq(comm(s3, sm), -2 * sm)
assert eq(comm(sp, sm), s3)

up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
print("σ+|↑> = 0:", eq(sp @ up, 0), " σ+|↓> = |↑>:", eq(sp @ dn, up))
print("σ-|↓> = 0:", eq(sm @ dn, 0), " σ-|↑> = |↓>:", eq(sm @ up, dn))
assert eq(sp @ up, 0) and eq(sp @ dn, up) and eq(sm @ dn, 0) and eq(sm @ up, up * 0 + dn)

print()
print("=== 射影とベクトルへの分解 ===")
P = (I2 + s3) / 2
Q = (I2 - s3) / 2
print("σ+ = Pσ1 = σ1Q:", eq(sp, P @ s1), eq(sp, s1 @ Q))
print("σ- = Qσ1 = σ1P:", eq(sm, Q @ s1), eq(sm, s1 @ P))
assert eq(sp, P @ s1) and eq(sp, s1 @ Q)
assert eq(sm, Q @ s1) and eq(sm, s1 @ P)

print("σ+σ- = P, σ-σ+ = Q:", eq(sp @ sm, P), eq(sm @ sp, Q))
print("σ±² = 0:", eq(sp @ sp, 0), eq(sm @ sm, 0))
print("P² = P, Q² = Q, P + Q = I, PQ = 0:",
      eq(P @ P, P), eq(Q @ Q, Q), eq(P + Q, I2), eq(P @ Q, 0))
assert eq(sp @ sm, P) and eq(sm @ sp, Q)
assert eq(sp @ sp, 0) and eq(sm @ sm, 0)
assert eq(P @ P, P) and eq(Q @ Q, Q) and eq(P + Q, I2) and eq(P @ Q, 0)

print("反交換 {σ+, σ-} = I（σ1, σ2 のクリフォード関係の書き換え）:",
      eq(anticomm(sp, sm), I2))
assert eq(anticomm(sp, sm), I2)

print()
print("=== ヌルベクトル（等方ベクトル） ===")
# 複素化したベクトル v = a σ1 + b σ2 + c σ3 に対し、v² = (a² + b² + c²) I
n_plus = np.array([1, 1j, 0], dtype=complex)   # σ1 + iσ2 = 2σ+
n_minus = np.array([1, -1j, 0], dtype=complex)  # σ1 - iσ2 = 2σ-


def gamma3(v):
    return v[0] * s1 + v[1] * s2 + v[2] * s3


print("n±·n± = 0（等方）:", eq(n_plus @ n_plus, 0), eq(n_minus @ n_minus, 0))
print("γ(n+)² = 0, γ(n-)² = 0:",
      eq(gamma3(n_plus) @ gamma3(n_plus), 0), eq(gamma3(n_minus) @ gamma3(n_minus), 0))
print("γ(n+) = 2σ+, γ(n-) = 2σ-:", eq(gamma3(n_plus), 2 * sp), eq(gamma3(n_minus), 2 * sm))
assert eq(n_plus @ n_plus, 0) and eq(n_minus @ n_minus, 0)
assert eq(gamma3(n_plus), 2 * sp) and eq(gamma3(n_minus), 2 * sm)
assert eq(gamma3(n_plus) @ gamma3(n_plus), 0)

# 3次元では極大等方部分空間の次元は1（n+ と n- は直交しないので同時には取れない）
print("n+·n- = 2 ≠ 0（両者は同じ等方部分空間に入れない）:", n_plus @ n_minus)
assert not np.isclose(n_plus @ n_minus, 0)

print()
print("=== 純粋スピノルとしての P ===")
print("σ+ P = 0（P の第1列 |↑> が n+ に消される）:", eq(sp @ P, 0))
print("σ- P ≠ 0:", not eq(sm @ P, 0))
assert eq(sp @ P, 0) and not eq(sm @ P, 0)
print("P = |↑><↑|（固有ベクトルの外積）:", eq(P, np.outer(up, up.conj())))
assert eq(P, np.outer(up, up.conj()))

# 左イデアル M2(C)P は第2列が 0 の行列全体（2次元 = スピノル空間）
basis = [np.zeros((2, 2), dtype=complex) for _ in range(4)]
for k, (i, j) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
    basis[k][i, j] = 1
ideal = [b @ P for b in basis]
rank = np.linalg.matrix_rank(np.array([m.flatten() for m in ideal]))
print("dim M2(C)P =", rank, "（スピノル空間 C² の次元と一致）")
assert rank == 2

print()
print("=== Cl_{2m}(C) のフォック構成 ===")


def kron(mats):
    out = np.array([[1]], dtype=complex)
    for m in mats:
        out = np.kron(out, m)
    return out


def fock(m):
    """m 個の生成消滅演算子とガンマ行列を返す（ジョルダン＝ウィグナー変換）。"""
    a = []
    for j in range(m):
        mats = [s3] * j + [sp] + [I2] * (m - j - 1)
        a.append(kron(mats))
    gam = []
    for j in range(m):
        gam.append(a[j] + a[j].conj().T)
        gam.append(-1j * (a[j] - a[j].conj().T))
    return a, gam


for m in range(1, 5):
    a, gam = fock(m)
    dim = 2 ** m
    # クリフォード関係 {γ_i, γ_j} = 2δ_ij
    ok_cl = all(eq(anticomm(gam[i], gam[j]), 2 * (i == j) * np.eye(dim))
                for i in range(2 * m) for j in range(2 * m))
    # 正準反交換関係 {a_i, a_j†} = δ_ij, {a_i, a_j} = 0
    ok_car = all(eq(anticomm(a[i], a[j].conj().T), (i == j) * np.eye(dim))
                 and eq(anticomm(a[i], a[j]), 0)
                 for i in range(m) for j in range(m))
    # w_j = (γ_{2j-1} + iγ_{2j})/2 = a_j
    ok_w = all(eq((gam[2 * j] + 1j * gam[2 * j + 1]) / 2, a[j]) for j in range(m))
    # 真空はすべての a_j に消される
    vac = np.zeros(dim, dtype=complex)
    vac[0] = 1
    ok_vac = all(eq(x @ vac, 0) for x in a)
    # 真空の射影 Π_j (I - iγ_{2j-1}γ_{2j})/2 は P = (I+σ3)/2 の一般化
    Pvac = np.eye(dim, dtype=complex)
    for j in range(m):
        Pvac = Pvac @ (np.eye(dim) - 1j * gam[2 * j] @ gam[2 * j + 1]) / 2
    ok_proj = eq(Pvac, np.outer(vac, vac.conj()))
    print(f"m={m}: dim={dim}, クリフォード関係={ok_cl}, CAR={ok_car}, "
          f"w_j=a_j={ok_w}, 真空の消滅={ok_vac}, 真空の射影={ok_proj}")
    assert ok_cl and ok_car and ok_w and ok_vac and ok_proj


def annihilator_dim(gam, psi):
    """{v ∈ C^{2m} : γ(v)ψ = 0} の次元。"""
    n = len(gam)
    M = np.array([g @ psi for g in gam]).T  # (dim × n)
    return n - np.linalg.matrix_rank(M, tol=1e-8)


def chirality(gam):
    m = len(gam) // 2
    G = np.eye(2 ** m, dtype=complex)
    for g in gam:
        G = G @ g
    return ((-1j) ** m) * G


print()
print("=== 純粋性の判定（消滅する等方部分空間の次元） ===")
rng = np.random.default_rng(0)
for m in range(1, 5):
    a, gam = fock(m)
    dim = 2 ** m
    G = chirality(gam)
    assert eq(G @ G, np.eye(dim))
    # 正カイラリティ（半スピノル）空間への射影
    Pplus = (np.eye(dim) + G) / 2
    half = np.linalg.matrix_rank(Pplus)
    vac = np.zeros(dim, dtype=complex)
    vac[0] = 1
    d_vac = annihilator_dim(gam, vac)
    dims = []
    for _ in range(20):
        psi = Pplus @ (rng.normal(size=dim) + 1j * rng.normal(size=dim))
        dims.append(annihilator_dim(gam, psi))
    print(f"m={m}: 半スピノル次元={half}, 真空の消滅次元={d_vac}, "
          f"一般の半スピノルの消滅次元={sorted(set(dims))}")
    assert d_vac == m
    assert half == dim // 2
    if m <= 3:
        assert set(dims) == {m}, "m ≤ 3 では半スピノルはすべて純粋"
    else:
        assert max(dims) < m, "m = 4 では一般の半スピノルは純粋でない"

print()
print("=== Cl_6 の純粋スピノルの安定化群（Spin(6) ≅ SU(4)） ===")
m = 3
a, gam = fock(m)
dim = 2 ** m
# spin(6) = span{γ_i γ_j / 2 (i<j)}（実係数で反エルミート）
gens = []
for i in range(6):
    for j in range(i + 1, 6):
        gens.append(gam[i] @ gam[j] / 2)
print("dim spin(6) =", len(gens))
assert len(gens) == 15
vac = np.zeros(dim, dtype=complex)
vac[0] = 1
# 真空（純粋スピノル）を消す実線形結合の次元
M = np.array([(g @ vac) for g in gens])  # (15 × 8) complex
Mr = np.hstack([M.real, M.imag])  # 実係数で解く
stab = 15 - np.linalg.matrix_rank(Mr, tol=1e-8)
print("純粋スピノルを固定する部分代数の次元 =", stab, "（su(3) の次元 8）")
assert stab == 8

print()
print("すべての検証に成功しました。")
