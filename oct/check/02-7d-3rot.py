"""Symbolic and numerical checks for oct/02-7d-3rot.md.

Verifies:
1. Product of pure imaginary octonions: vw = -v·w + v×w, and the 7D cross product
   defined as combinations of 3D cross products across the 7 Fano triples.
2. Associativity: (rx)r* = r(xr*) for r = exp(theta/2 * n) and pure imaginary x.
3. Component expansion: Symbolic verification that the expanded components of (rx)r*
   match the 7-component formula in the article and the 7D Rodrigues formula:
   exp(theta/2 * n) x exp(-theta/2 * n) = cos(theta) x + sin(theta) (n×x) + (1 - cos(theta)) n (n·x).
4. Geometric interpretation: For n = e_1, the component x_1 is invariant and
   rotations of angle theta occur synchronously in the 3 orthogonal planes
   (e_2, e_3), (e_4, e_5), and (e_7, e_6).
5. Numerical consistency check with random octonions.
"""

import numpy as np
import sympy as sp

from common.octonion import L, R, basis_mul, triples


def check_product_and_cross():
    print("=== 1. 純虚八元数の積と7次元外積の検証 ===")
    # Check that basis_mul gives v*w = -v·w + v×w
    # For basis imaginary units e_a, e_b (1..7):
    for a in range(1, 8):
        for b in range(1, 8):
            sign, c = basis_mul(a, b)
            if a == b:
                assert sign == -1 and c == 0, f"e_{a}^2 must be -1"
            else:
                assert c != 0, f"e_{a}*e_{b} must be imaginary for a != b"

    # Cross product from Fano triples:
    # 7 triples: 123, 145, 176, 246, 257, 347, 365
    v = [sp.Symbol(f"v{i}", real=True) for i in range(1, 8)]
    w = [sp.Symbol(f"w{i}", real=True) for i in range(1, 8)]

    # Product vw via L matrix
    L_v = sum((v[i - 1] * sp.Matrix(L[i]) for i in range(1, 8)), sp.zeros(8, 8))
    w_vec = sp.Matrix([0] + w)
    vw = L_v @ w_vec

    # Scalar part = - v · w
    v_dot_w = sum(v[i] * w[i] for i in range(7))
    assert sp.simplify(vw[0] - (-v_dot_w)) == 0, "Scalar part must be -v·w"

    # Imaginary part matches the 7D cross product definition
    cross_expected = [
        (v[1] * w[2] - v[2] * w[1]) + (v[3] * w[4] - v[4] * w[3]) + (v[6] * w[5] - v[5] * w[6]),
        (v[2] * w[0] - v[0] * w[2]) + (v[3] * w[5] - v[5] * w[3]) + (v[4] * w[6] - v[6] * w[4]),
        (v[0] * w[1] - v[1] * w[0]) + (v[3] * w[6] - v[6] * w[3]) + (v[5] * w[4] - v[4] * w[5]),
        (v[4] * w[0] - v[0] * w[4]) + (v[5] * w[1] - v[1] * w[5]) + (v[6] * w[2] - v[2] * w[6]),
        (v[0] * w[3] - v[3] * w[0]) + (v[6] * w[1] - v[1] * w[6]) + (v[2] * w[5] - v[5] * w[2]),
        (v[0] * w[6] - v[6] * w[0]) + (v[1] * w[3] - v[3] * w[1]) + (v[4] * w[2] - v[2] * w[4]),
        (v[5] * w[0] - v[0] * w[5]) + (v[1] * w[4] - v[4] * w[1]) + (v[2] * w[3] - v[3] * w[2]),
    ]
    for i in range(7):
        diff = sp.simplify(vw[i + 1] - cross_expected[i])
        assert diff == 0, f"Cross product component {i+1} mismatch"
    print("  OK: 純虚八元数の積の実部が-内積、虚部が7次元外積に一致")


def check_symbolic_components():
    print("=== 2. 単位八元数による回転の成分計算とロドリゲス公式の検証 ===")
    n = [sp.Symbol(f"n{i}", real=True) for i in range(1, 8)]
    x = [sp.Symbol(f"x{i}", real=True) for i in range(1, 8)]
    theta = sp.Symbol("theta", real=True)

    # Unit pure imaginary n: sum(n_i^2) = 1
    L_n = sum((n[i - 1] * sp.Matrix(L[i]) for i in range(1, 8)), sp.zeros(8, 8))
    R_n = sum((n[i - 1] * sp.Matrix(R[i]) for i in range(1, 8)), sp.zeros(8, 8))

    c = sp.cos(theta / 2)
    s = sp.sin(theta / 2)

    L_r = c * sp.eye(8) + s * L_n
    R_rc = c * sp.eye(8) - s * R_n

    x_vec = sp.Matrix([0] + x)

    # (rx)r* and r(xr*)
    rx_r = R_rc @ L_r @ x_vec
    r_xr = L_r @ R_rc @ x_vec

    # 1. Check alternativity: (rx)r* = r(xr*)
    assert sp.simplify(rx_r - r_xr) == sp.zeros(8, 1), "Alternativity failed: (rx)r* != r(xr*)"
    print("  OK: 結合性 (rx)r* = r(xr*) [交代代数の性質]")

    # 2. Compare against the explicit components from oct/02-7d-3rot.md
    cth = sp.cos(theta)
    sth = sp.sin(theta)
    n_dot_x = sum(n[i] * x[i] for i in range(7))

    expected_components = [
        cth * x[0] + sth * ((n[1]*x[2] - n[2]*x[1]) + (n[3]*x[4] - n[4]*x[3]) + (n[6]*x[5] - n[5]*x[6])) + (1 - cth) * n[0] * n_dot_x,
        cth * x[1] + sth * ((n[2]*x[0] - n[0]*x[2]) + (n[3]*x[5] - n[5]*x[3]) + (n[4]*x[6] - n[6]*x[4])) + (1 - cth) * n[1] * n_dot_x,
        cth * x[2] + sth * ((n[0]*x[1] - n[1]*x[0]) + (n[3]*x[6] - n[6]*x[3]) + (n[5]*x[4] - n[4]*x[5])) + (1 - cth) * n[2] * n_dot_x,
        cth * x[3] + sth * ((n[4]*x[0] - n[0]*x[4]) + (n[5]*x[1] - n[1]*x[5]) + (n[6]*x[2] - n[2]*x[6])) + (1 - cth) * n[3] * n_dot_x,
        cth * x[4] + sth * ((n[0]*x[3] - n[3]*x[0]) + (n[6]*x[1] - n[1]*x[6]) + (n[2]*x[5] - n[5]*x[2])) + (1 - cth) * n[4] * n_dot_x,
        cth * x[5] + sth * ((n[0]*x[6] - n[6]*x[0]) + (n[1]*x[3] - n[3]*x[1]) + (n[4]*x[2] - n[2]*x[4])) + (1 - cth) * n[5] * n_dot_x,
        cth * x[6] + sth * ((n[5]*x[0] - n[0]*x[5]) + (n[1]*x[4] - n[4]*x[1]) + (n[2]*x[3] - n[3]*x[2])) + (1 - cth) * n[6] * n_dot_x,
    ]

    # Verify each component using half-angle formulas and sum(n_i^2) = 1
    for i in range(7):
        diff = rx_r[i + 1] - expected_components[i]
        diff = diff.rewrite(sp.cos).expand()
        diff = diff.subs({
            sp.cos(theta / 2)**2: (1 + cth) / 2,
            sp.sin(theta / 2)**2: (1 - cth) / 2,
            sp.sin(theta / 2) * sp.cos(theta / 2): sth / 2,
        }).expand()
        diff = sp.simplify(diff.subs(n[6]**2, 1 - sum(n[k]**2 for k in range(6))))
        assert diff == 0, f"Component e_{i+1} does not match article formula"

    print("  OK: 7成分の展開式が記事中の明示式と完全に一致")


def check_example_e1():
    print("=== 3. 回転軸 n = e_1 の場合の幾何学的構造の検証 ===")
    theta = sp.Symbol("theta", real=True)
    cth = sp.cos(theta)
    sth = sp.sin(theta)

    # When n = e_1: n = (1, 0, 0, 0, 0, 0, 0)
    L_e1 = sp.Matrix(L[1])
    R_e1 = sp.Matrix(R[1])

    c = sp.cos(theta / 2)
    s = sp.sin(theta / 2)

    L_r = c * sp.eye(8) + s * L_e1
    R_rc = c * sp.eye(8) - s * R_e1
    M_rot = R_rc @ L_r

    # Restrict to 7D imaginary subspace (indices 1..7)
    M7 = M_rot[1:, 1:]
    # Substitute half-angle identities
    M7_simp = sp.zeros(7, 7)
    for r in range(7):
        for col in range(7):
            val = M7[r, col].rewrite(sp.cos).expand()
            val = val.subs({
                sp.cos(theta / 2)**2: (1 + cth) / 2,
                sp.sin(theta / 2)**2: (1 - cth) / 2,
                sp.sin(theta / 2) * sp.cos(theta / 2): sth / 2,
            }).expand()
            M7_simp[r, col] = sp.simplify(val)

    # Check axis preservation: e1 -> e1
    assert M7_simp[0, 0] == 1 and all(M7_simp[0, j] == 0 for j in range(1, 7))
    assert all(M7_simp[i, 0] == 0 for i in range(1, 7))
    print("  OK: 回転軸 e_1 方向は不変")

    # Check 3 rotation planes:
    # 1. Plane (e_2, e_3): [ [cos, -sin], [sin, cos] ]
    block_23 = sp.Matrix([[M7_simp[1, 1], M7_simp[1, 2]], [M7_simp[2, 1], M7_simp[2, 2]]])
    assert block_23 == sp.Matrix([[cth, -sth], [sth, cth]]), "(e_2, e_3) plane must rotate by theta"

    # 2. Plane (e_4, e_5): [ [cos, -sin], [sin, cos] ]
    block_45 = sp.Matrix([[M7_simp[3, 3], M7_simp[3, 4]], [M7_simp[4, 3], M7_simp[4, 4]]])
    assert block_45 == sp.Matrix([[cth, -sth], [sth, cth]]), "(e_4, e_5) plane must rotate by theta"

    # 3. Plane (e_7, e_6): note basis ordering in article
    # x_6' = cos theta x_6 + sin theta x_7
    # x_7' = -sin theta x_6 + cos theta x_7
    # in (e_7, e_6) coordinates:
    # [x_7']   [ cos theta  -sin theta ] [x_7]
    # [x_6'] = [ sin theta   cos theta ] [x_6]
    block_76 = sp.Matrix([[M7_simp[6, 6], M7_simp[6, 5]], [M7_simp[5, 6], M7_simp[5, 5]]])
    assert block_76 == sp.Matrix([[cth, -sth], [sth, cth]]), "(e_7, e_6) plane must rotate by theta"

    print("  OK: 3つの直交平面 (e_2,e_3), (e_4,e_5), (e_7,e_6) で角度 theta の等傾回転")


def check_numerical_random():
    print("=== 4. ランダムベクトルによる数値検証 ===")
    rng = np.random.default_rng(42)

    def oct_mul(a, b):
        return sum(a[i] * L[i] for i in range(8)) @ b

    for _ in range(50):
        # Random unit pure imaginary n
        n7 = rng.normal(size=7)
        n7 /= np.linalg.norm(n7)
        n = np.zeros(8)
        n[1:] = n7

        # Random pure imaginary x
        x7 = rng.normal(size=7)
        x = np.zeros(8)
        x[1:] = x7

        theta = rng.uniform(-np.pi, np.pi)

        # r = cos(theta/2) + sin(theta/2) n
        r = np.zeros(8)
        r[0] = np.cos(theta / 2)
        r[1:] = np.sin(theta / 2) * n7

        r_conj = r.copy()
        r_conj[1:] *= -1

        # Direct octonion multiplication: (rx)r* and r(xr*)
        rx_r = oct_mul(oct_mul(r, x), r_conj)
        r_xr = oct_mul(r, oct_mul(x, r_conj))

        assert np.allclose(rx_r, r_xr), "Associativity failed numerically"
        assert np.isclose(rx_r[0], 0.0, atol=1e-12), "Scalar component must be 0"

        # Rodrigues formula:
        # 7D cross product n x x
        n_cross_x = np.zeros(8)
        n_cross_x[1] = (n7[1]*x7[2] - n7[2]*x7[1]) + (n7[3]*x7[4] - n7[4]*x7[3]) + (n7[6]*x7[5] - n7[5]*x7[6])
        n_cross_x[2] = (n7[2]*x7[0] - n7[0]*x7[2]) + (n7[3]*x7[5] - n7[5]*x7[3]) + (n7[4]*x7[6] - n7[6]*x7[4])
        n_cross_x[3] = (n7[0]*x7[1] - n7[1]*x7[0]) + (n7[3]*x7[6] - n7[6]*x7[3]) + (n7[5]*x7[4] - n7[4]*x7[5])
        n_cross_x[4] = (n7[4]*x7[0] - n7[0]*x7[4]) + (n7[5]*x7[1] - n7[1]*x7[5]) + (n7[6]*x7[2] - n7[2]*x7[6])
        n_cross_x[5] = (n7[0]*x7[3] - n7[3]*x7[0]) + (n7[6]*x7[1] - n7[1]*x7[6]) + (n7[2]*x7[5] - n7[5]*x7[2])
        n_cross_x[6] = (n7[0]*x7[6] - n7[6]*x7[0]) + (n7[1]*x7[3] - n7[3]*x7[1]) + (n7[4]*x7[2] - n7[2]*x7[4])
        n_cross_x[7] = (n7[5]*x7[0] - n7[0]*x7[5]) + (n7[1]*x7[4] - n7[4]*x7[1]) + (n7[2]*x7[3] - n7[3]*x7[2])

        n_dot_x = np.dot(n7, x7)
        rodrigues = np.cos(theta) * x + np.sin(theta) * n_cross_x + (1 - np.cos(theta)) * n_dot_x * n

        assert np.allclose(rx_r, rodrigues, atol=1e-12), "Rodrigues formula mismatch"

    print("  OK: 50回のランダム試行で直接計算とロドリゲス公式が完全一致")


if __name__ == "__main__":
    check_product_and_cross()
    check_symbolic_components()
    check_example_e1()
    check_numerical_random()
    print("\nAll checks passed successfully!")
