"""variable-dependence.md の数式の検証"""

from sympy import Function, cos, diff, exp, simplify, sin, symbols

x, Y = symbols("x Y")


def check(name, expr):
    ok = simplify(expr) == 0
    print(("OK  " if ok else "NG  ") + name)
    assert ok


def both(name, f, y, result):
    """先に代入する手順と、全微分後に代入する手順を比べる"""
    first = diff(f.subs(Y, y), x)
    fx, fy = diff(f, x), diff(f, Y)
    later = (fx + fy * diff(y, x)).subs(Y, y)
    check(name + "：先に代入", first - result)
    check(name + "：全微分後に代入", later - result)


both("xy, y=x^2", x * Y, x**2, 3 * x**2)
both("sin(xy), y=e^x", sin(x * Y), exp(x),
     cos(x * exp(x)) * (exp(x) + x * exp(x)))
both("x^2+y^2, y=sin x", x**2 + Y**2, sin(x), 2 * x + 2 * sin(x) * cos(x))

# d/dt (∂L/∂q̇) の展開（L の例として多項式と三角関数を含む式）
t, Q, V = symbols("t Q V")
q = Function("q")(t)
L = Q**2 * V**3 + sin(Q * V) * t**2 + V * exp(t)
on_path = {Q: q, V: diff(q, t)}
lhs = diff(diff(L, V).subs(on_path), t)
rhs = (diff(L, V, Q) * diff(q, t) + diff(L, V, V) * diff(q, t, 2)
       + diff(L, V, t)).subs(on_path)
check("d/dt ∂L/∂q̇ の展開", lhs - rhs)
