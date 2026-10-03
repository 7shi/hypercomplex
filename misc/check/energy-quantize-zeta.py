"""energy-quantize-zeta.md の数式の検証"""

import mpmath
from sympy import (
    Rational, Sum, cos, diff, exp, integrate, limit, oo, pi, simplify,
    symbols, zeta,
)

nu, T, k, c, h, L = symbols("nu T k c h L", positive=True)
E, x, kk = symbols("E x kk", positive=True)
n = symbols("n", integer=True, positive=True)


def check(name, expr):
    ok = simplify(expr) == 0
    print(("OK  " if ok else "NG  ") + name)
    assert ok


planck = 8 * pi * h * nu**3 / c**3 / (exp(h * nu / (k * T)) - 1)
rj = 8 * pi * nu**2 / c**3 * k * T
wien = 8 * pi * h * nu**3 / c**3 / exp(h * nu / (k * T))

# 低周波極限・高周波極限
check("低周波極限でレイリー・ジーンズ", limit(planck / rj, nu, 0) - 1)
check("高周波極限でウィーン", limit(planck / wien, nu, oo) - 1)

# モード密度
modes = 4 * pi * kk**2 / (pi / L) ** 3 / 8 * 2
check("モード数 L^3 k^2/π^2", modes - L**3 * kk**2 / pi**2)
dens = (kk**2 / pi**2).subs(kk, 2 * pi * nu / c) * diff(2 * pi * nu / c, nu)
check("モード密度 8πν^2/c^3", dens - 8 * pi * nu**2 / c**3)

# 連続的なエネルギーの平均
num = integrate(E * exp(-E / (k * T)), (E, 0, oo))
den = integrate(exp(-E / (k * T)), (E, 0, oo))
check("分母 = kT", den - k * T)
check("分子 = (kT)^2", num - (k * T) ** 2)
check("<E> = kT", num / den - k * T)

# 離散的なエネルギーの平均
q = symbols("q", positive=True)
check("Σ n q^n = q/(1-q)^2", q * diff(1 / (1 - q), q) - q / (1 - q) ** 2)
q0 = exp(-h * nu / (k * T))
avg = h * nu * (q / (1 - q) ** 2) / (1 / (1 - q))
check("<E> = hν/(e^{hν/kT}-1)",
      avg.subs(q, q0) - h * nu / (exp(h * nu / (k * T)) - 1))
# 数値でも級数を確認（hν/kT = 0.7）
r = Rational(7, 10)
s_num = sum(m * r * exp(-m * r) for m in range(400)).evalf(30)
s_den = sum(exp(-m * r) for m in range(400)).evalf(30)
err = abs(s_num / s_den - (r / (exp(r) - 1)).evalf(30))
print(("OK  " if err < 1e-20 else "NG  ") + f"<E> 級数の数値（誤差 {float(err):.1e}）")
assert err < 1e-20

# 全エネルギー密度
check("∫x^3 e^{-nx} = 6/n^4", integrate(x**3 * exp(-n * x), (x, 0, oo)) - 6 / n**4)
check("∫x^2 e^{-nx} = 2/n^3", integrate(x**2 * exp(-n * x), (x, 0, oo)) - 2 / n**3)
check("∫x e^{-nx} = 1/n^2", integrate(x * exp(-n * x), (x, 0, oo)) - 1 / n**2)
val = mpmath.quad(lambda t: t**3 / mpmath.expm1(t), [0, mpmath.inf])
err = abs(val - 6 * mpmath.zeta(4))
print(("OK  " if err < 1e-12 else "NG  ") + f"∫x^3/(e^x-1) = 6ζ(4)（数値積分、誤差 {float(err):.1e}）")
assert err < 1e-12
check("ζ(4) = π^4/90", zeta(4) - pi**4 / 90)
uT = 8 * pi * k**4 * T**4 / (c**3 * h**3) * 6 * zeta(4)
check("u(T) = 48πk^4T^4ζ(4)/(c^3h^3)", uT - 48 * pi * k**4 * T**4 / (c**3 * h**3) * zeta(4))
check("u(T) = 8π^5k^4T^4/(15c^3h^3)", uT - 8 * pi**5 * k**4 * T**4 / (15 * c**3 * h**3))
sigma = 2 * pi**5 * k**4 / (15 * c**2 * h**3)
check("u(T) = 4σT^4/c", uT - 4 * sigma / c * T**4)

# x^2 のフーリエ係数とパーセヴァル
a0 = integrate(x**2, (x, -pi, pi)) / pi
an = integrate(x**2 * cos(n * x), (x, -pi, pi)) / pi
check("a0 = 2π^2/3（定数項 a0/2 = π^2/3）", a0 - 2 * pi**2 / 3)
check("an = 4(-1)^n/n^2", simplify(an - 4 * (-1) ** n / n**2))
lhs = integrate(x**4, (x, -pi, pi)) / (2 * pi)
check("(1/2π)∫x^4 = π^4/5", lhs - pi**4 / 5)
rhs = a0**2 / 4 + Rational(1, 2) * Sum(16 / n**4, (n, 1, oo)).doit()
check("パーセヴァルの両辺", lhs - rhs)
