"""integration-by-parts.md の数式の検証"""

from sympy import Rational, cos, diff, integrate, pi, simplify, sin, symbols

x = symbols("x")


def check(name, expr):
    ok = simplify(expr) == 0
    print(("OK  " if ok else "NG  ") + name)
    assert ok


# 不定積分の結果は微分して被積分関数に戻ることで確かめる
r1 = x * sin(x) + cos(x)
check("∫x cos x dx", diff(r1, x) - x * cos(x))

r2 = -x**2 * cos(x) + 2 * x * sin(x) + 2 * cos(x)
check("∫x^2 sin x dx", diff(r2, x) - x**2 * sin(x))

# 途中式：∫x^2 d(cos x) = x^2 cos x - ∫2x cos x dx
check("∫x^2(-sin x)dx の途中式",
      integrate(-x**2 * sin(x), x)
      - (x**2 * cos(x) - integrate(2 * x * cos(x), x)))

# 次数が上がる変形：∫x^2 sin x dx = (x^3/3) sin x - (1/3)∫x^3 cos x dx
check("次数が上がる変形",
      integrate(x**2 * sin(x), x)
      - (sin(x) * x**3 / 3
         - Rational(1, 3) * integrate(x**3 * cos(x), x)))

# 定積分の部分積分：∫_0^π x cos x dx = [x sin x]_0^π - ∫_0^π sin x dx
check("定積分の部分積分",
      integrate(x * cos(x), (x, 0, pi))
      - ((x * sin(x)).subs(x, pi) - (x * sin(x)).subs(x, 0)
         - integrate(sin(x), (x, 0, pi))))
