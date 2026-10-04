import numpy as np

from common.octonion import L  # 8x8 matrices: L[i]x = e_i x

rng = np.random.default_rng(0)


def ok(label, cond):
    print(f"  {'OK ' if cond else 'NG '} {label}")
    assert cond, label


def unit(n):
    return n / np.linalg.norm(n)


def refl(v, n):
    """v - 2(v.n)n"""
    return v - 2 * (v @ n) * n


# --- 八元数 ---
def omul(x, y):
    return sum(x[i] * L[i] for i in range(8)) @ y


def oconj(x):
    c = x.copy()
    c[1:] *= -1
    return c


def odot(v, w):
    return 0.5 * (omul(v, w) - omul(w, v) + omul(oconj(v), w) + omul(v, oconj(w)))[0]


print("=== 八元数（7次元）: v' = nvn, 内積 = -Re(vw) ===")
for _ in range(5):
    v = np.r_[0, rng.normal(size=7)]
    w = np.r_[0, rng.normal(size=7)]
    n = np.r_[0, unit(rng.normal(size=7))]
    ok("内積 -Re(vw)", np.isclose(-omul(v, w)[0], v @ w))
    ok("v' = nvn（括弧は交代性で不問）",
       np.allclose(omul(omul(n, v), n), refl(v, n)[:8] if False else v - 2 * (v @ n) * n)
       and np.allclose(omul(n, omul(v, n)), omul(omul(n, v), n)))
    ok("v' = -n v* n", np.allclose(-omul(omul(n, oconj(v)), n), omul(omul(n, v), n)))

print("=== 八元数（8次元）: o' = -n o* n, 内積の式 ===")
for _ in range(5):
    o = rng.normal(size=8)
    p = rng.normal(size=8)
    n = unit(rng.normal(size=8))
    ok("内積の式", np.isclose(odot(o, p), o @ p))
    ok("o' = -n o* n", np.allclose(-omul(omul(n, oconj(o)), n), refl(o, n)))
    ok("-n o* n の括弧は不問", np.allclose(omul(n, omul(oconj(o), n)), omul(omul(n, oconj(o)), n)))

print("=== 内積 Re(vw*) = (vw*+wv*)/2（四元数4次元・八元数8次元）、導出の再結合 ===")
for _ in range(5):
    o = rng.normal(size=8)
    p = rng.normal(size=8)
    n = unit(rng.normal(size=8))
    ok("Re(vw*) = v.w", np.isclose(omul(o, oconj(p))[0], o @ p))
    ok("(vw*+wv*)/2 = v.w", np.allclose(0.5 * (omul(o, oconj(p)) + omul(p, oconj(o))), o @ p * np.eye(8)[0]))
    ok("2(q.n) = qn*+nq*", np.allclose(omul(o, oconj(n)) + omul(n, oconj(o)), 2 * (o @ n) * np.eye(8)[0]))
    ok("(on*)n = o(n*n) = o", np.allclose(omul(omul(o, oconj(n)), n), o))

print("=== 四元数 = 八元数の部分代数 {1,e1,e2,e3}（i,j,k = e1,e2,e3） ===")
for _ in range(5):
    q = np.r_[rng.normal(size=4), np.zeros(4)]
    n = np.r_[unit(rng.normal(size=4)), np.zeros(4)]
    ok("4次元鏡映 -nq*n", np.allclose(-omul(omul(n, oconj(q)), n), refl(q, n)))

print("=== 7次元鏡映が4次元鏡映を含む（実部 -> e4、i,j,k -> e1,e2,e3） ===")
for _ in range(5):
    qw = rng.normal(size=4)  # (w,x,y,z)
    nw = unit(rng.normal(size=4))
    # 四元数側（全成分を八元数の{1,e1,e2,e3}に置いて計算）
    q = np.r_[qw, np.zeros(4)]
    n = np.r_[nw, np.zeros(4)]
    qp = -omul(omul(n, oconj(q)), n)[:4]
    # 八元数側（実部をe4へ）
    def emb(a):
        x = np.zeros(8)
        x[4] = a[0]
        x[1:4] = a[1:]
        return x
    v7, n7 = emb(qw), emb(nw)
    vp = omul(omul(n7, v7), n7)
    ok("e4 に埋め込んだ nvn が 4次元鏡映と一致", np.allclose(vp, emb(qp)))
    ok("e4 は i,j,k と反交換", all(np.allclose(omul(L[4][:, 0] * 0 + np.eye(8)[4], np.eye(8)[a]),
                                           -omul(np.eye(8)[a], np.eye(8)[4])) for a in (1, 2, 3)))

print("=== 鏡映の合成（回転）: 7次元・8次元は括弧が外せないことの確認 ===")
for _ in range(5):
    v = np.r_[0, rng.normal(size=7)]
    m = np.r_[0, unit(rng.normal(size=7))]
    n = np.r_[0, unit(rng.normal(size=7))]
    two = -omul(omul(m, oconj(omul(omul(n, oconj(v)), n))), m)  # -m(-nv*n)* m の形
    two = -omul(omul(m, oconj(-omul(omul(n, oconj(v)), n))), m)
    R = refl(refl(v, n), m)
    ok("7D: m(nvn)m = 2回の鏡映", np.allclose(omul(omul(m, omul(omul(n, v), n)), m), R))
    ok("7D: -m(-nv*n)*m = m(nvn)m", np.allclose(two, R))
    flat = omul(omul(omul(omul(m, n), v), n), m)  # ((mn)v)n)m: 括弧を外した左結合
    print("     mnvnm（左結合）と一致:", np.allclose(flat, R))

for _ in range(5):
    o = rng.normal(size=8)
    m = unit(rng.normal(size=8))
    n = unit(rng.normal(size=8))
    first = -omul(omul(n, oconj(o)), n)
    second = -omul(omul(m, oconj(first)), m)
    ok("8D: -m(-no*n)*m = m(n*on*)m", np.allclose(second, omul(omul(m, omul(omul(oconj(n), o), oconj(n))), m)))
    ok("8D: 2回の鏡映と一致", np.allclose(second, refl(refl(o, n), m)))

# --- 複素数 ---
print("=== 複素数 ===")
for _ in range(5):
    z = complex(*rng.normal(size=2))
    n = np.exp(1j * rng.uniform(0, 2 * np.pi))
    m = np.exp(1j * rng.uniform(0, 2 * np.pi))
    zv = np.array([z.real, z.imag])
    nv = np.array([n.real, n.imag])
    ok("内積 Re(a* b)", np.isclose((z.conjugate() * n).real, zv @ nv))
    zp = -z.conjugate() * n * n
    ok("z' = -z* n^2", np.allclose([zp.real, zp.imag], refl(zv, nv)))
    zpp = -m * (-n * z.conjugate() * n).conjugate() * m
    ok("合成 = z (n* m)^2", np.isclose(zpp, z * (n.conjugate() * m) ** 2))

# --- クリフォード代数 Cl_{n,0}(R) ---
print("=== クリフォード代数 Cl_{n,0}(R) ===")


def blade_mul(a, b):
    """基底ブレード（ビットマスク）の積: (符号, ブレード)。e_i^2 = +1"""
    s = 1
    x = a >> 1
    while x:
        if bin(x & b).count("1") % 2:
            s = -s
        x >>= 1
    return s, a ^ b


def cmul(X, Y):
    Z = {}
    for a, x in X.items():
        for b, y in Y.items():
            s, c = blade_mul(a, b)
            Z[c] = Z.get(c, 0) + s * x * y
    return Z


def cadd(X, Y, k=1.0):
    Z = dict(X)
    for b, y in Y.items():
        Z[b] = Z.get(b, 0) + k * y
    return Z


def vec(a):
    return {1 << i: float(x) for i, x in enumerate(a)}


def close(X, Y):
    return all(np.isclose(X.get(k, 0), Y.get(k, 0)) for k in set(X) | set(Y))


for dim in (2, 3, 4, 8):
    for _ in range(3):
        v = rng.normal(size=dim)
        w = rng.normal(size=dim)
        n = unit(rng.normal(size=dim))
        m = unit(rng.normal(size=dim))
        V, W, N, M = vec(v), vec(w), vec(n), vec(m)
        ok(f"Cl_{dim},0: 内積 (vw+wv)/2", close({0: v @ w}, {0: 0.5 * (cmul(V, W)[0] + cmul(W, V)[0])})
           and close(cadd(cmul(V, W), cmul(W, V)), {0: 2 * (v @ w)}))
        ok(f"Cl_{dim},0: v' = -nvn", close({k: -x for k, x in cmul(cmul(N, V), N).items()}, vec(refl(v, n))))
        two = {k: -x for k, x in cmul(cmul(M, {k: -x for k, x in cmul(cmul(N, V), N).items()}), M).items()}
        ok(f"Cl_{dim},0: 合成 = mnvnm", close(two, cmul(cmul(cmul(cmul(M, N), V), N), M)))
        ok(f"Cl_{dim},0: 合成 = 2回の鏡映", close(two, vec(refl(refl(v, n), m))))

print("=== 射影行列 ===")
n = unit(rng.normal(size=5))
P = np.outer(n, n)
ok("冪等", np.allclose(P @ P, P))
H = np.eye(5) - 2 * P
ok("I-2nn^T は対合かつ直交", np.allclose(H @ H, np.eye(5)) and np.allclose(H.T @ H, np.eye(5)))
ok("det = -1", np.isclose(np.linalg.det(H), -1))

print("すべて確認しました。")
