"""spherical-trig.md の数式の数値検証。"""

import numpy as np

rng = np.random.default_rng(0)

ok = True


def check(name, cond):
    global ok
    print(("OK  " if cond else "NG  ") + name)
    ok &= bool(cond)


def mul(p, q):
    a0, a1, a2, a3 = p
    b0, b1, b2, b3 = q
    return np.array([
        a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
        a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
        a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
        a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0,
    ])


def prod(*qs):
    r = ONE
    for q in qs:
        r = mul(r, q)
    return r


def conj(q):
    return q * np.array([1, -1, -1, -1])


def inner(q, r):
    return mul(q, conj(r))[0]


def ex(u, t):
    """e^{u t}（u は単位純虚四元数）"""
    return np.cos(t) * ONE + np.sin(t) * u


def close(x, y):
    return np.allclose(x, y, atol=1e-10)


ONE = np.array([1.0, 0, 0, 0])
I = np.array([0, 1.0, 0, 0])
J = np.array([0, 0, 1.0, 0])
K = np.array([0, 0, 0, 1.0])

# 四元数の基礎
check("i^2=j^2=k^2=-1", all(close(mul(u, u), -ONE) for u in (I, J, K)))
check("ij=k, jk=i, ki=j", close(mul(I, J), K) and close(mul(J, K), I) and close(mul(K, I), J))
q, r, s = rng.normal(size=(3, 4))
check("共役は反準同型", close(conj(mul(q, r)), mul(conj(r), conj(q))))
check("|qr|=|q||r|", np.isclose(np.linalg.norm(mul(q, r)), np.linalg.norm(q) * np.linalg.norm(r)))
check("<q,r> は成分の内積", np.isclose(inner(q, r), q @ r))
check("随伴性 <qr,s>=<q,s r̄>", np.isclose(inner(mul(q, r), s), inner(q, mul(s, conj(r)))))
check("随伴性 <rq,s>=<q,r̄ s>", np.isclose(inner(mul(r, q), s), inner(q, mul(conj(r), s))))
check("Re(qr)=Re(rq)（随伴性の証明で使用）", np.isclose(mul(q, r)[0], mul(r, q)[0]))

# i 軸周りの回転
th = 0.7
h = ex(I, th / 2)
check("e^{iθ/2} i e^{-iθ/2} = i", close(prod(h, I, conj(h)), I))
check("e^{iθ/2} j e^{-iθ/2} = cosθ j + sinθ k", close(prod(h, J, conj(h)), np.cos(th) * J + np.sin(th) * K))
check("e^{iθ/2} k e^{-iθ/2} = cosθ k - sinθ j", close(prod(h, K, conj(h)), np.cos(th) * K - np.sin(th) * J))
check("[i,j]=2k", close(mul(I, J) - mul(J, I), 2 * K))

# スピノルの周期
n = rng.normal(size=3)
n = np.r_[0, n / np.linalg.norm(n)]
check("n^2=-1", close(mul(n, n), -ONE))
check("e^{2πn/2}=-1, e^{4πn/2}=1", close(ex(n, np.pi), -ONE) and close(ex(n, 2 * np.pi), ONE))
check("±q は同じ回転", close(prod(-h, J, conj(-h)), prod(h, J, conj(h))))

# 回転の合成
q1, q2 = ex(n, 0.4), ex(I, 1.1)
v = np.r_[0, rng.normal(size=3)]
check("q2(q1 v q̄1)q̄2 = (q2q1)v(q2q1)‾", close(prod(q2, q1, v, conj(q1), conj(q2)), prod(mul(q2, q1), v, conj(mul(q2, q1)))))


# ランダムな球面三角形で公式を確認
def triangle(sign=-1):
    while True:
        P = rng.normal(size=(3, 3))
        A, B, C = (p / np.linalg.norm(p) for p in P)
        if np.sign(np.linalg.det(np.array([A, B, C]))) == sign:  # -1：外側から見て時計回り
            break
    side = lambda x, y: np.arccos(np.clip(x @ y, -1, 1))
    a, b, c = side(B, C), side(C, A), side(A, B)

    def angle(p, x, y):
        tx, ty = x - (x @ p) * p, y - (y @ p) * p
        return np.arccos(np.clip(tx @ ty / np.linalg.norm(tx) / np.linalg.norm(ty), -1, 1))

    return a, b, c, angle(A, B, C), angle(B, C, A), angle(C, A, B)


def master(a, b, c, al, be, ga):
    return prod(ex(I, (np.pi - al) / 2), ex(K, b / 2), ex(I, (np.pi - ga) / 2),
                ex(K, a / 2), ex(I, (np.pi - be) / 2), ex(K, c / 2))


tris = [triangle() for _ in range(200)]
check("マスター方程式 (4) の左辺 = -1（200個）", all(close(master(*t), -ONE) for t in tris))
check("反時計回りの三角形でも辺と角だけで同じ式が成り立つ",
      all(close(master(*triangle(+1)), -ONE) for _ in range(50)))

# 球の1/8
h8 = 0.5 * (ONE + I - J + K)
check("e^{iπ/4}e^{kπ/4} = (1+i-j+k)/2", close(mul(ex(I, np.pi / 4), ex(K, np.pi / 4)), h8))
u = np.r_[0, 1, -1, 1] / np.sqrt(3)
check("(1+i-j+k)/2 = e^{uπ/3}", close(ex(u, np.pi / 3), h8))
check("球の1/8で (4) = -1", close(master(*[np.pi / 2] * 6), -ONE))

# 無限小三角形：1次の係数
al, be = 0.9, 1.3
ga = np.pi - al - be
c = 1.0
a, b = c * np.sin(al) / np.sin(ga), c * np.sin(be) / np.sin(ga)
ei = lambda t: np.exp(1j * t)
check("平面三角形の閉合 b e^{i(π-γ)}e^{i(π-β)} + a e^{i(π-β)} + c = 0",
      abs(b * ei(np.pi - ga) * ei(np.pi - be) + a * ei(np.pi - be) + c) < 1e-12)
eps = 1e-5
lhs = prod(ex(I, (np.pi - al) / 2), ex(K, eps * b / 2), ex(I, (np.pi - ga) / 2),
           ex(K, eps * a / 2), ex(I, (np.pi - be) / 2), ex(K, eps * c / 2))
check("平面の角を使った無限小三角形で (4) = -1 + O(ε^2)", np.linalg.norm(lhs + ONE) < 10 * eps**2)

# (5) と (6)
ok5 = ok6 = True
for a, b, c, al, be, ga in tris:
    L5 = prod(ex(I, (np.pi - ga) / 2), ex(K, a / 2), ex(I, (np.pi - be) / 2))
    R5 = -prod(ex(K, -b / 2), ex(I, -(np.pi - al) / 2), ex(K, -c / 2))
    L6 = prod(ex(I, (np.pi - ga) / 2), ex(K, -a / 2), ex(I, (np.pi - be) / 2))
    R6 = -prod(ex(K, b / 2), ex(I, -(np.pi - al) / 2), ex(K, c / 2))
    ok5 &= close(L5, R5)
    ok6 &= close(L6, R6) and close(prod(I, L5, conj(I)), L6)
check("式 (5)", ok5)
check("式 (6) は (5) の i による共役", ok6)

# 内積の計算過程
a, b, c, al, be, ga = tris[0]
X = ex(I, -(np.pi - al) / 2)
check("左辺の内積 = cos a", np.isclose(inner(ex(K, a / 2), ex(K, -a / 2)), np.cos(a)))
r1 = inner(prod(ex(K, -b / 2), X, ex(K, -c / 2)), prod(ex(K, b / 2), X, ex(K, c / 2)))
r2 = inner(mul(ex(K, -b), X), mul(X, ex(K, c)))
r3 = inner(prod(conj(X), ex(K, -b), X), ex(K, c))
check("右辺の内積の変形（3段）", np.isclose(r1, r2) and np.isclose(r2, r3))
rot = prod(conj(X), ex(K, -b), X)
check("e^{-kb} の回転", close(rot, np.cos(b) * ONE + np.sin(b) * np.cos(al) * K + np.sin(b) * np.sin(al) * J))
check("右辺の内積 = cos b cos c + sin b sin c cos α", np.isclose(r3, np.cos(b) * np.cos(c) + np.sin(b) * np.sin(c) * np.cos(al)))

# 正弦定理・5要素の規則（Ā i A と B̄ i B）
okc = oks = okf = okstep = True
for a, b, c, al, be, ga in tris:
    okc &= np.isclose(np.cos(a), np.cos(b) * np.cos(c) + np.sin(b) * np.sin(c) * np.cos(al))
    A = prod(ex(I, (np.pi - ga) / 2), ex(K, a / 2), ex(I, (np.pi - be) / 2))
    B = -prod(ex(K, -b / 2), ex(I, -(np.pi - al) / 2), ex(K, -c / 2))
    AiA, BiB = prod(conj(A), I, A), prod(conj(B), I, B)
    okstep &= close(AiA, BiB)
    okstep &= np.isclose(inner(AiA, K), np.sin(a) * np.sin(be))
    okstep &= np.isclose(inner(BiB, K), np.sin(b) * np.sin(al))
    okstep &= np.isclose(inner(AiA, J), np.sin(a) * np.cos(be))
    okstep &= np.isclose(inner(BiB, J), np.cos(b) * np.sin(c) - np.sin(b) * np.cos(al) * np.cos(c))
    eka = ex(K, a / 2)
    okstep &= close(prod(conj(eka), I, eka), np.cos(a) * I - np.sin(a) * J)
    ekb = ex(K, b / 2)
    okstep &= close(prod(ekb, I, conj(ekb)), np.cos(b) * I + np.sin(b) * J)
    ekc = ex(K, c / 2)
    okstep &= close(prod(conj(ekc), J, ekc), np.cos(c) * J + np.sin(c) * I)
    oks &= np.isclose(np.sin(al) / np.sin(a), np.sin(be) / np.sin(b))
    okf &= np.isclose(np.sin(a) * np.cos(be), np.cos(b) * np.sin(c) - np.sin(b) * np.cos(c) * np.cos(al))
check("Ā i A = B̄ i B と各内積・途中の回転", okstep)
check("球面余弦定理（幾何から計算した辺と角）", okc)
check("球面正弦定理", oks)
check("5要素の規則", okf)

# 直角三角形とネイピアの法則・日の出方程式
okn = True
for _ in range(100):
    b, c = rng.uniform(0.1, np.pi / 2 - 0.1, 2)
    if c >= b:
        b, c = c, b
    a = np.arccos(np.cos(b) / np.cos(c))  # β=π/2 の余弦定理 cos b = cos a cos c
    al = np.arccos((np.cos(a) - np.cos(b) * np.cos(c)) / (np.sin(b) * np.sin(c)))
    okn &= np.isclose(np.cos(al), np.tan(c) / np.tan(b))
check("ネイピアの法則 cos α = tan c / tan b", okn)

phi, dl = np.radians(35), np.radians(23.44)
om = np.arccos(-np.tan(phi) * np.tan(dl))
check("cos(π-ω) = tanφ / tan(π/2-δ) と同値", np.isclose(np.cos(np.pi - om), np.tan(phi) / np.tan(np.pi / 2 - dl)))
# 天球上で直接計算：赤道座標で太陽の高度が0になる時角
H = np.arccos(-np.tan(phi) * np.tan(dl))
alt = np.arcsin(np.sin(phi) * np.sin(dl) + np.cos(phi) * np.cos(dl) * np.cos(H))
check("得られた時角で太陽の高度が0", abs(alt) < 1e-12)
check("24ω/π = 24 - (24/π)arccos(tanφ tanδ)",
      np.isclose(24 / np.pi * om, 24 - 24 / np.pi * np.arccos(np.tan(phi) * np.tan(dl))))
check("白夜の境界で ω = π", np.isclose(-np.tan(np.pi / 2 - dl) * np.tan(dl), -1))
check("極夜（φ > π/2+δ, δ<0）でも |tanφ tanδ| > 1",
      abs(np.tan(np.radians(80)) * np.tan(-dl)) > 1)

# トレース：A=i, B=(cos c)i-(sin c)j から始め、時計回りの三角形で各頂点が順に針の位置 i に来ること
def rot(q, v):
    return prod(q, np.r_[0, v], conj(q))[1:]


def side(x, y):
    return np.arccos(np.clip(x @ y, -1, 1))


def angle(p, x, y):
    tx, ty = x - (x @ p) * p, y - (y @ p) * p
    return np.arccos(np.clip(tx @ ty / np.linalg.norm(tx) / np.linalg.norm(ty), -1, 1))


okcw, okccw = True, True
for _ in range(100):
    c, al, b = rng.uniform(0.2, 2.5), rng.uniform(0.2, 2.8), rng.uniform(0.2, 2.5)
    A3, B3, tAB = np.array([1.0, 0, 0]), np.array([np.cos(c), -np.sin(c), 0]), np.array([0, -1.0, 0])
    for sgn in (1, -1):
        C3 = np.cos(b) * A3 + np.sin(b) * (np.cos(al) * tAB + sgn * np.sin(al) * np.cross(A3, tAB))
        cw = np.linalg.det(np.array([A3, B3, C3])) < 0
        a, be, ga = side(B3, C3), angle(B3, C3, A3), angle(C3, A3, B3)
        ops = [ex(K, c / 2), ex(I, (np.pi - be) / 2), ex(K, a / 2), ex(I, (np.pi - ga) / 2), ex(K, b / 2), ex(I, (np.pi - al) / 2)]
        q2, q4, q6 = ops[0], prod(*ops[2::-1]), prod(*ops[4::-1])
        Q = prod(*ops[::-1])
        r = (np.allclose(rot(q2, B3), A3) and np.allclose(rot(q4, C3), A3)
             and np.allclose(rot(q6, A3), A3) and np.allclose(rot(Q, B3), B3))
        if cw:
            okcw &= r
        else:
            okccw &= not r
check("トレース：時計回りでは B→C→A が順に i に来て最後に B が元の位置に戻る", okcw)
check("トレース：反時計回りでは同じ操作で頂点が i に来ない（向きの規約）", okccw)
t = tris[0]
check("(4) を k で共役した式（i→-i、反時計回り版）も右辺 -1",
      close(prod(K, master(*t), conj(K)), -ONE)
      and close(prod(K, ex(I, 0.3), conj(K)), ex(I, -0.3)) and close(prod(K, ex(K, 0.3), conj(K)), ex(K, 0.3)))
# 半角：微分 [i,v] と左右の乗算の効果
check("d/dt e^{ti}ve^{-ti}|0 = [i,v]：[i,j]=2k, [i,k]=-2j",
      close(mul(I, J) - mul(J, I), 2 * K) and close(mul(I, K) - mul(K, I), -2 * J))
th = 0.8
L, R = ex(I, th / 2), ex(I, -th / 2)
check("左右の乗算：j では同じ向きに θ/2 ずつ回り、i では打ち消し合う",
      close(mul(L, J), mul(J, R)) and close(prod(L, I, R), I))
check("左乗算だけでは純虚四元数が純虚にならない", abs(mul(L, I)[0]) > 1e-3)
b0, c0 = np.pi / 2, phi  # δ=0
check("δ=0：変形前の式から cos α = cos b sin c / (sin b cos c) = 0",
      np.isclose(np.cos(b0) * np.sin(c0) / (np.sin(b0) * np.cos(c0)), 0))

print("all OK" if ok else "FAILED")
