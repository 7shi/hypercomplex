"""Draws the figures for gca-schrodinger.md as SVG (standard library only).

- gca-schrodinger-rotations.svg: rotation between components (a, b) and
  rotation within a component (a complex number psi).
- gca-schrodinger-two-level.svg: a 2-level system (psi_0, psi_1) under a
  rotation between components and a rotation within components.
- gca-schrodinger-three-level.svg: a 3-level system under X, Z, XZ and ZX;
  the results of XZ and ZX differ by the phase w (a third of a turn).
- gca-schrodinger-dial.svg: the n-level system (n = 8) on a dial; each small
  circle is the complex plane of one component.  X shifts the components,
  Z rotates the j-th component by the angle of w^j.
- gca-schrodinger-waves.svg: real parts of the waves w_1, w_2 for n = 8.
  The hollow point at the right end is the 0th component again.
- gca-schrodinger-spread.svg: bell-shaped |psi|^2 in position and in wave
  number for the widths in the table of the article (n = 256).
- gca-schrodinger-levels.svg: the lowest values of E for n = 256 (taken
  from clif/check/gca-schrodinger.py) against the continuous classical
  energy.

Usage: uv run clif/images/gca-schrodinger.py
PNG (for upload, not tracked):
    for f in clif/images/gca-schrodinger-*.svg; do
        rsvg-convert "$f" -o "${f%.svg}.png"
    done
"""

from cmath import exp, phase
from math import cos, exp as rexp, pi, radians, sin
from pathlib import Path

DIR = Path(__file__).parent
LIGHT, MAIN, OLD, AXIS, TEXT = "#9db7d5", "#1f5fa8", "#aaaaaa", "#888888", "#333333"


def save(name, w, h, body):
    head = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
        f'viewBox="0 0 {w} {h}" font-family="Noto Serif, BIZ UDPMincho, serif" font-size="15" fill="{TEXT}">'
    )
    bg = f'<rect width="{w}" height="{h}" fill="#ffffff"/>'
    path = DIR / f"gca-schrodinger-{name}.svg"
    path.write_text("\n".join([head, bg, *body, "</svg>"]) + "\n")
    print(path)


def line(x1, y1, x2, y2, color=AXIS, width=1, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
        f'stroke="{color}" stroke-width="{width}"{d}/>'
    )


def circle(cx, cy, r, stroke=LIGHT, fill="none", width=1.5):
    return (
        f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{width}"/>'
    )


def text(x, y, s, anchor="middle", size=15, color=TEXT, italic=False):
    st = ' font-style="italic"' if italic else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
        f'font-size="{size}" fill="{color}"{st}>{s}</text>'
    )


def arrow(cx, cy, z, scale, color=MAIN, width=2.5):
    """Arrow from (cx, cy) representing the complex number z (y axis up)."""
    x, y = cx + scale * z.real, cy - scale * z.imag
    a = phase(z)
    hl, hw = 9, 4
    bx, by = x - hl * cos(a), y + hl * sin(a)
    p1 = (bx - hw * sin(a), by - hw * cos(a))
    p2 = (bx + hw * sin(a), by + hw * cos(a))
    return [
        line(cx, cy, bx, by, color, width),
        f'<polygon points="{x:.2f},{y:.2f} {p1[0]:.2f},{p1[1]:.2f} '
        f'{p2[0]:.2f},{p2[1]:.2f}" fill="{color}"/>',
    ]


def arc(cx, cy, r, a0, a1, color=MAIN):
    """Counterclockwise arc from angle a0 to a1 (radians, y axis up)."""
    x0, y0 = cx + r * cos(a0), cy - r * sin(a0)
    x1, y1 = cx + r * cos(a1), cy - r * sin(a1)
    return (
        f'<path d="M{x0:.2f},{y0:.2f} A{r},{r} 0 0 0 {x1:.2f},{y1:.2f}" '
        f'fill="none" stroke="{color}" stroke-width="1.2"/>'
    )


def axes(cx, cy, r, xlabel, ylabel, italic=True):
    e = r + 16
    return [
        line(cx - e, cy, cx + e, cy),
        line(cx, cy + e, cx, cy - e),
        text(cx + e + 4, cy + 5, xlabel, "start", 14, AXIS, italic),
        text(cx, cy - e - 6, ylabel, "middle", 14, AXIS, italic),
    ]


def rotations():
    w, h, r = 520, 250, 70
    a0, th = radians(20), radians(50)
    z0, z1 = exp(1j * a0), exp(1j * (a0 + th))
    body = []
    for cx, title, xl, yl, it, before, after in (
        (135, "成分間の回転", "a", "b", True, "(a, b)", "回転後"),
        (395, "成分内の回転", "実部", "虚部", False, "ψ",
         'e<tspan dy="-6" font-size="11">iθ</tspan><tspan dy="6">ψ</tspan>'),
    ):
        cy = 140
        body.append(text(cx, 24, title, size=15))
        body += axes(cx, cy, r, xl, yl, it)
        body.append(circle(cx, cy, r))
        body += arrow(cx, cy, z0, r, OLD)
        body += arrow(cx, cy, z1, r)
        body.append(arc(cx, cy, 26, a0, a0 + th))
        m = a0 + th / 2
        body.append(text(cx + 38 * cos(m), cy - 38 * sin(m) + 5, "θ", size=14,
                         color=MAIN, italic=True))
        body.append(text(cx + (r + 8) * z0.real, cy - (r + 8) * z0.imag + 5,
                         before, "start", 14, AXIS, it or before == "ψ"))
        body.append(text(cx + (r + 8) * z1.real, cy - (r + 8) * z1.imag - 2,
                         after, "start", 14, MAIN, before == "ψ"))
    save("rotations", w, h, body)


def two_level():
    r, cw, rh, left, top = 46, 150, 156, 150, 34
    w, h = left + cw * 2 + 10, top + rh * 3
    p0 = (0.8 * exp(1j * radians(30)), 0.6 * exp(1j * radians(30)))
    t = radians(30)
    between = (p0[0] * cos(t) - p0[1] * sin(t), p0[0] * sin(t) + p0[1] * cos(t))
    within = (p0[0] * exp(1j * t), p0[1] * exp(-1j * t))
    body = [
        text(left + cw * (c + 0.5), 22, f'ψ<tspan dy="4" font-size="11" font-style="normal">{s}</tspan>', italic=True)
        for c, s in enumerate("01")
    ]
    rows = (("元の状態", p0, None), ("成分間の回転", between, p0),
            ("成分内の回転", within, p0))
    for i, (label, psi, old) in enumerate(rows):
        cy = top + rh * i + r + 14
        body.append(text(left - 22, cy + 5, label, "end"))
        if i:
            body.append(line(10, cy - r - 24, w - 10, cy - r - 24, "#dddddd"))
        for c in range(2):
            cx = left + cw * (c + 0.5)
            body.append(line(cx - r - 8, cy, cx + r + 8, cy))
            body.append(line(cx, cy + r + 8, cx, cy - r - 8))
            body.append(circle(cx, cy, r, width=1))
            arrows = [(psi[c], MAIN)] + ([(old[c], OLD)] if old else [])
            for z, color in sorted(arrows, key=lambda t: -abs(t[0])):
                body += arrow(cx, cy, z, r, color)
            body.append(text(cx, cy + r + 28,
                             f"割合 {abs(psi[c]) ** 2:.2f}", size=13, color=AXIS))
    save("two-level", w, h, body)


def waves():
    n, w, ph, left, right, amp = 8, 480, 130, 70, 24, 40
    step = (w - left - right) / n
    body = []
    for i, k in enumerate((1, 2)):
        mid = ph * i + 20 + amp
        px = lambda j: left + step * j
        py = lambda v: mid - amp * v
        body.append(
            f'<text x="{left - 16}" y="{mid + 5}" text-anchor="end">'
            f'<tspan font-style="italic">k</tspan> = {k}</text>'
        )
        body.append(line(px(0), mid, px(n), mid))
        pts = " ".join(
            f"{px(t / 20):.2f},{py(cos(2 * pi * k * t / 20 / n)):.2f}"
            for t in range(20 * n + 1)
        )
        body.append(
            f'<polyline points="{pts}" fill="none" stroke="{LIGHT}" stroke-width="2"/>'
        )
        for j in range(n + 1):
            x, y = px(j), py(cos(2 * pi * k * j / n))
            body.append(line(x, mid, x, y, MAIN, 1, "2 2"))
            body.append(circle(x, y, 4.5, MAIN, MAIN if j < n else "#ffffff"))
            body.append(text(x, mid + amp + 22, j % n, size=12, color=AXIS))
    save("waves", w, ph * 2 + 10, body)


def ring(cx, cy, R, vals, amax, label_r):
    """Components (z, color) as arrows in small circles placed on a ring."""
    n, out = len(vals), [circle(cx, cy, R, "#dddddd", width=1)]
    for j, (z, color) in enumerate(vals):
        a = 2 * pi * j / n
        sx, sy = cx + R * cos(a), cy - R * sin(a)
        out.append(
            text(cx + label_r * cos(a), cy - label_r * sin(a) + 4, j,
                 size=12, color=AXIS)
        )
        out.append(circle(sx, sy, amax, LIGHT, "#ffffff", 1))
        out += arrow(sx, sy, z, amax, color, 2.2)
    return out


def three_level():
    R, cw, rh, top, amax = 46, 210, 216, 48, 24
    w1 = exp(2j * pi / 3)
    orig = [(1.0 + 0j, "#d9731a"), (0.8 + 0j, MAIN), (0.6 + 0j, "#3a8f4a")]
    X = lambda v: [v[(j - 1) % 3] for j in range(3)]
    Z = lambda v: [(z * w1**j, c) for j, (z, c) in enumerate(v)]
    panels = (
        (0, 0, "元の並び", orig),
        (1, 0, "X（ずらす）", X(orig)),
        (2, 0, "Z（回す）", Z(orig)),
        (1, 1, "XZ（回してからずらす）", X(Z(orig))),
        (2, 1, "ZX（ずらしてから回す）", Z(X(orig))),
    )
    body = []
    for col, row, title, vals in panels:
        cx, cy = cw * (col + 0.5), top + rh * row + R + 40
        body.append(text(cx, top + rh * row - 22, title))
        body += ring(cx, cy, R, vals, amax, R + amax + 12)
    save("three-level", cw * 3, top + rh * 2 - 30, body)


def dial():
    n, R, cw, top, amax = 8, 82, 250, 44, 21
    ACCENT = "#d9731a"
    base = [1.0, 0.85, 0.7, 0.58, 0.58, 0.58, 0.7, 0.85]
    w1 = exp(2j * pi / n)
    orig = [(v + 0j, ACCENT if j == 0 else MAIN) for j, v in enumerate(base)]
    shifted = [orig[(j - 1) % n] for j in range(n)]
    turned = [(z * w1**j, c) for j, (z, c) in enumerate(orig)]
    body = []
    for i, (title, vals) in enumerate(
        (("元の並び", orig), ("X（ずらす）", shifted), ("Z（回す）", turned))
    ):
        cx, cy = cw * (i + 0.5), top + R + 30
        body.append(text(cx, 24, title))
        body += ring(cx, cy, R, vals, amax, R - 36)
    save("dial", cw * 3, top + 2 * R + 60, body)


def spread():
    pw, ph, left, gap, top, rng = 200, 62, 30, 40, 44, 2.0
    rows = ((0.050, 0.245), (0.111, 0.111), (0.600, 0.020))
    w, rh = left + pw * 2 + gap + 20, ph + 46
    body = []
    for c, title in enumerate(("場所", "波数")):
        body.append(text(left + (pw + gap) * c + pw / 2, 22, title))
    for i, sig in enumerate(rows):
        base = top + rh * i + ph
        for c, (s, name) in enumerate(zip(sig, "xp")):
            x0 = left + (pw + gap) * c
            px = lambda v: x0 + pw * (v + rng) / (2 * rng)
            pts = " ".join(
                f"{px(-rng + 2 * rng * t / 400):.2f},"
                f"{base - (ph - 6) * rexp(-(-rng + 2 * rng * t / 400) ** 2 / (2 * s * s)):.2f}"
                for t in range(401)
            )
            body.append(
                f'<polygon points="{px(-rng):.2f},{base} {pts} {px(rng):.2f},{base}" '
                f'fill="{LIGHT}" fill-opacity="0.5" stroke="{MAIN}" stroke-width="1.5"/>'
            )
            body.append(line(x0, base, x0 + pw, base))
            body.append(line(px(0), base, px(0), base + 4))
            body.append(text(px(0), base + 17, 0, size=11, color=AXIS))
            body.append(text(x0 + pw + 6, base + 4, name, "start", 14, AXIS, True))
            body.append(
                f'<text x="{x0 + 4:.2f}" y="{base - ph + 14:.2f}" font-size="13">'
                f'<tspan font-style="italic">σ</tspan>'
                f'<tspan dy="3" font-size="10" font-style="italic">{name}</tspan>'
                f'<tspan dy="-3"> = {s:.3f}</tspan></text>'
            )
    save("spread", w, top + rh * 3 - 14, body)


def levels():
    E = (0.4985, 1.4923, 2.4801, 3.4618)
    w, h, top, unit, ax = 400, 270, 40, 52, 70
    y = lambda e: h - 26 - unit * e
    body = [line(ax, y(0), ax, y(4.15))]
    body.append(text(ax, y(4.15) - 8, "E", size=14, color=AXIS, italic=True))
    for k in range(5):
        body.append(line(ax - 4, y(k), ax, y(k)))
        lab = "0" if k == 0 else ("ℏ" if k == 1 else f"{k}ℏ")
        body.append(text(ax - 9, y(k) + 5, lab, "end", 13, AXIS))
    c1, c2, bw = 135, 285, 90
    body.append(text(c1, 24, "古典力学"))
    body.append(text(c2, 24, "モデル"))
    body.append(
        f'<rect x="{c1 - bw / 2}" y="{y(4):.2f}" width="{bw}" '
        f'height="{y(0) - y(4):.2f}" fill="{LIGHT}" fill-opacity="0.6"/>'
    )
    body.append(line(c1 - bw / 2, y(0), c1 + bw / 2, y(0), MAIN, 1.5))
    for e in E:
        body.append(line(c2 - bw / 2, y(e), c2 + bw / 2, y(e), MAIN, 2.5))
        body.append(text(c2 + bw / 2 + 8, y(e) + 5, f"{e:.4f} ℏ", "start", 13))
    body.append(line(c2 - bw / 2, y(0), c2 + bw / 2, y(0), AXIS, 1, "2 2"))
    save("levels", w, h, body)


if __name__ == "__main__":
    rotations()
    two_level()
    three_level()
    dial()
    waves()
    spread()
    levels()
