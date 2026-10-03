"""exp-maclaurin-integral.md の数式の検証"""

from sympy import E, Integral, exp, factorial, simplify, symbols

x = symbols("x")


def check(name, expr):
    ok = simplify(expr) == 0
    print(("OK  " if ok else "NG  ") + name)
    assert ok


# 出発点：e^x = 1 + ∫_0^x e^t dt
t = symbols("t")
check("e^x = 1 + ∫_0^x e^t dt", exp(x) - (1 + Integral(exp(t), (t, 0, x)).doit()))

# 1 の n 重積分 ∫_0^x ∫_0^{x'} … 1 = x^n/n!（内側から順に積分する）
N = 8
xs = symbols(f"x0:{N + 1}")  # xs[0] = x, xs[k] = k 番目の積分変数


def nested(f, n):
    """∫_0^{x} ∫_0^{x'} … ∫_0^{x^{(n-1)}} f(x^{(n)}) dx^{(n)} … dx'"""
    g = f(xs[n])
    for k in range(n, 0, -1):
        g = Integral(g, (xs[k], 0, xs[k - 1])).doit()
    return g


for n in range(N + 1):
    check(f"1 の {n} 重積分 = x^{n}/{n}!",
          nested(lambda s: 1, n) - xs[0]**n / factorial(n))

# n 回代入後の等式：e^x = Σ_{k<n} x^k/k! + (e^x の n 重積分)
for n in range(1, 6):
    check(f"{n} 回代入した等式",
          exp(xs[0]) - (sum(xs[0]**k / factorial(k) for k in range(n))
                        + nested(exp, n)))

# 剰余項 R_n（e^x の n 重積分）の評価と収束
# x≧0 で 0 ≦ R_n ≦ e^x x^n/n!、x<0 で |R_n| ≦ |x|^n/n!
# n≦5 は上で n 重積分との一致を確かめたので、R_n = e^x - Σ_{k<n} x^k/k! で計算する
for x0 in [3, -3]:
    for n in [5, 10, 20]:
        R = float(E**x0 - sum(x0**k / factorial(k) for k in range(n)))
        bound = (float(E**x0) if x0 > 0 else 1.0) * abs(x0)**n / float(factorial(n))
        ok = abs(R) <= bound * (1 + 1e-12)
        print(("OK  " if ok else "NG  ") + f"|R_{n}({x0})| = {abs(R):.3e} ≦ {bound:.3e}")
        assert ok
