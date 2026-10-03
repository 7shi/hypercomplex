"""epsilon-euler-lagrange.md の数式の検証"""

from sympy import Function, Rational, diff, exp, simplify, sin, symbols

t, eps, Q, V, m = symbols("t epsilon Q V m")
q = Function("q")(t)
eta = Function("eta")(t)


def check(name, expr):
    ok = simplify(expr) == 0
    print(("OK  " if ok else "NG  ") + name)
    assert ok


# ラグランジアンの例（Q, V は位置・速度の引数）
L = Q**2 * V**3 + sin(Q * V) * t + V * exp(Q)
LQ, LV = diff(L, Q), diff(L, V)


def on(expr, x, v):
    return expr.subs({Q: x, V: v})


# dL/dε|_{ε=0} = (∂L/∂q)η + (∂L/∂q̇)η̇
qe, ve = q + eps * eta, diff(q, t) + eps * diff(eta, t)
dL = diff(on(L, qe, ve), eps).subs(eps, 0)
check("dL/dε|_{ε=0}",
      dL - (on(LQ, q, diff(q, t)) * eta
            + on(LV, q, diff(q, t)) * diff(eta, t)))

# 部分積分の被積分関数の恒等式：η̇ L_q̇ = d/dt(η L_q̇) - η d/dt L_q̇
P = on(LV, q, diff(q, t))
check("部分積分の恒等式",
      diff(eta, t) * P - (diff(eta * P, t) - eta * diff(P, t)))

# L = T - V（T = m q̇^2/2）のとき ∂L/∂q = -V'(q)、d/dt ∂L/∂q̇ = m q̈
U = Function("U")
LN = Rational(1, 2) * m * V**2 - U(Q)
check("∂L/∂q = -∂V/∂q",
      on(diff(LN, Q), q, diff(q, t)) + diff(U(q), q))
check("d/dt ∂L/∂q̇ = m q̈",
      diff(on(diff(LN, V), q, diff(q, t)), t) - m * diff(q, t, 2))

# 基本補題の証明の η = (t-a)^2(b-t)^2 は端点で値も導関数も 0（0 との接続が C^1）
a, b = symbols("a b")
bump = (t - a)**2 * (b - t)**2
check("η(a)=η(b)=0", bump.subs(t, a) + bump.subs(t, b))
check("η'(a)=0", diff(bump, t).subs(t, a))
check("η'(b)=0", diff(bump, t).subs(t, b))
