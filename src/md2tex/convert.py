"""記事（Markdown）をLaTeXに変換する。"""

from __future__ import annotations

import re
import subprocess
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent

BOX = {
    "axm": "公理", "def": "定義", "thm": "定理", "cor": "系", "lem": "補題",
    "conj": "予想", "prop": "命題", "fml": "公式", "prf": "証明",
    "ex": "例", "exc": "問題", "rem": "注意",
}
COLOR = {
    "def": "green!45!black", "fml": "blue!55!black", "thm": "red!55!black",
    "ex": "orange!60!black", "rem": "gray!70!black",
}
DEFAULT_COLOR = "black!60"
# KaTeXに合わせ、ϕは\phi、φは\varphiとする（ϵ・εも同様）
GREEK = {
    "α": r"\alpha", "β": r"\beta", "γ": r"\gamma", "δ": r"\delta",
    "ϵ": r"\epsilon", "ε": r"\varepsilon", "ζ": r"\zeta", "η": r"\eta",
    "θ": r"\theta", "ι": r"\iota", "κ": r"\kappa", "λ": r"\lambda",
    "μ": r"\mu", "ν": r"\nu", "ξ": r"\xi", "π": r"\pi", "ρ": r"\rho",
    "σ": r"\sigma", "τ": r"\tau", "υ": r"\upsilon", "ϕ": r"\phi",
    "φ": r"\varphi", "χ": r"\chi", "ψ": r"\psi", "ω": r"\omega",
    "Γ": r"\Gamma", "Δ": r"\Delta", "Θ": r"\Theta", "Λ": r"\Lambda",
    "Ξ": r"\Xi", "Π": r"\Pi", "Σ": r"\Sigma", "Φ": r"\Phi", "Ψ": r"\Psi",
    "Ω": r"\Omega",
}
GREEK_RE = re.compile("[" + "".join(GREEK) + "]")

PREAMBLE = r"""\documentclass[a4paper,11pt]{ltjsarticle}
\usepackage[margin=25mm]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage[most]{tcolorbox}
\usepackage[hidelinks]{hyperref}
\tcbset{breakable,fonttitle=\bfseries\sffamily,boxrule=0.6pt,arc=1.5mm,
  left=3mm,right=3mm,top=2mm,bottom=2mm,before skip=4mm,after skip=4mm}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.6\baselineskip}
"""


def _tsv(name: str) -> list[list[str]]:
    return [ln.split("\t") for ln in (ROOT / name).read_text(encoding="utf-8").splitlines()]


class References:
    """[[slug]]を初出順の番号に置き換え、末尾の一覧を作る。"""

    def __init__(self) -> None:
        self.refs = tomllib.loads((ROOT / "refs.toml").read_text(encoding="utf-8"))
        self.master = tomllib.loads((ROOT / "refs-master.toml").read_text(encoding="utf-8"))
        self.slug2md = {p[1]: p[0] for p in _tsv("slugs.tsv") if len(p) == 2}
        self.md2title = {p[0]: p[1] for p in _tsv("md.tsv") if len(p) >= 2}
        self.md2url = {
            p[2]: "https://mathlog.info" + p[1]
            for p in _tsv("articles.tsv")
            if len(p) >= 4 and p[1].startswith("/")
        }
        self.cited: list[str] = []

    def cite(self, slug: str) -> str:
        if slug not in self.cited:
            self.cited.append(slug)
        return rf"\textsuperscript{{[{self.cited.index(slug) + 1}]}}"

    def bibliography(self) -> list[str]:
        if not self.cited:
            return []
        res = [
            r"\section*{参考文献}",
            r"\begin{enumerate}\renewcommand{\labelenumi}{[\arabic{enumi}]}",
        ]
        for slug in self.cited:
            r = {**self.master.get(slug, {}), **self.refs.get(slug, {})}
            md = self.slug2md.get(slug, "")
            text = r.get("citation") or r.get("title") or self.md2title.get(md, slug)
            url = r.get("url") or self.md2url.get(md)
            item = escape(re.sub(r",?\s*https?://\S+", "", text))
            if url:
                item += rf"\\ \url{{{url}}}"
            res.append(r"\item " + item)
        res.append(r"\end{enumerate}")
        return res


def math(s: str) -> str:
    return GREEK_RE.sub(lambda m: GREEK[m[0]] + " ", s)


def escape(s: str) -> str:
    s = s.replace("\\", r"\textbackslash{}")
    s = re.sub(r"([&%#_{}])", r"\\\1", s)
    return s.replace("~", r"\textasciitilde{}").replace("^", r"\textasciicircum{}")


def inline(s: str, refs: References) -> str:
    res = []
    for i, part in enumerate(re.split(r"(\$[^$]+\$)", s)):
        if i % 2:
            res.append("$" + math(part[1:-1]) + "$")
            continue
        for j, q in enumerate(re.split(r"(\*\*|\[\[[^\]]+\]\])", part)):
            if q == "**":
                res.append("\0")
            elif j % 2:
                res.append(refs.cite(q[2:-2]))
            else:
                res.append(escape(q))
    t = "".join(res)
    if t.count("\0") % 2:
        raise ValueError(f"**が対になっていません: {s}")
    return re.sub("\0(.*?)\0", r"\\textbf{\1}", t)


def image(src: Path, path: str, out: Path) -> str:
    """SVGはPDFに変換し、等倍（本文幅を超えるときは本文幅）で取り込む。"""
    img = src.parent / path
    if img.suffix != ".svg":
        raise ValueError(f"SVG以外の画像には対応していません: {path}")
    pdf = out / (img.stem + ".pdf")
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", str(pdf), str(img)], check=True)
    m = re.search(r'<svg[^>]*\bwidth="([\d.]+)"', img.read_text(encoding="utf-8"))
    if not m:
        raise ValueError(f"SVGの幅が読み取れません: {path}")
    w = float(m[1]) * 0.75  # px → bp
    return (
        "\\begin{center}\n"
        rf"\includegraphics[width=\dimexpr\ifdim {w}bp<\linewidth {w}bp\else\linewidth\fi\relax]{{{pdf.name}}}"
        "\n\\end{center}"
    )


def convert(src: Path, out: Path) -> Path:
    """記事srcを変換し、outに書き出したtexのパスを返す。"""
    out.mkdir(parents=True, exist_ok=True)
    refs = References()
    lines = src.read_text(encoding="utf-8").splitlines()
    body: list[str] = []
    counters: dict[str, int] = {}
    lists: list[int] = []  # 開いている箇条書きのインデント
    in_box = False

    def close_lists(to: int = 0) -> None:
        while len(lists) > to:
            lists.pop()
            body.append(r"\end{itemize}")

    i = 0
    while i < len(lines):
        ln = lines[i]
        i += 1
        if ln.strip() == "$$":
            close_lists()
            buf = []
            while lines[i].strip() != "$$":
                buf.append(math(lines[i]))
                i += 1
            i += 1
            body.append("\\[\n" + "\n".join(buf) + "\n\\]")
        elif ln.startswith("&&&"):
            close_lists()
            if in_box:
                if ln.strip() != "&&&":
                    raise ValueError(f"囲みが閉じていません: {ln}")
                body.append(r"\end{tcolorbox}")
                in_box = False
                continue
            m = re.match(r"&&&(\w*)\s*(.*?)(?:\s*\[([\w-]+)\])?\s*$", ln)
            typ, title = m[1], m[2]
            if typ and typ not in BOX:
                raise ValueError(f"未対応の囲みです: {ln}")
            head = inline(title, refs)
            if typ:
                counters[typ] = counters.get(typ, 0) + 1
                head = f"{BOX[typ]} {counters[typ]}" + (f"（{head}）" if head else "")
            opt = f"colframe={COLOR.get(typ, DEFAULT_COLOR)},colback=white"
            if head:
                opt += f",title={{{head}}}"
            body.append(rf"\begin{{tcolorbox}}[{opt}]")
            in_box = True
        elif m := re.match(r"(#{1,3}) (.*)", ln):
            close_lists()
            cmd = ("section", "subsection", "subsubsection")[len(m[1]) - 1]
            body.append(rf"\{cmd}{{{inline(m[2], refs)}}}")
        elif m := re.match(r"!\[(.*?)\]\((.*?)\)\s*$", ln):
            close_lists()
            body.append(image(src, m[2].split(" ")[0], out))
        elif m := re.match(r"( *)- (.*)", ln):
            ind = len(m[1])
            while lists and lists[-1] > ind:
                close_lists(len(lists) - 1)
            if not lists or lists[-1] < ind:
                lists.append(ind)
                body.append(r"\begin{itemize}")
            body.append(r"\item " + inline(m[2], refs))
        elif ln.strip() == "":
            close_lists()
            body.append("")
        elif re.match(r"\s*(\||>|```|\d+\. |#)", ln):
            raise ValueError(f"未対応の記法です: {ln}")
        else:
            close_lists()
            body.append(inline(ln, refs))
    close_lists()
    if in_box:
        raise ValueError("囲みが閉じていません")

    try:
        title = refs.md2title.get(str(src.relative_to(ROOT)), src.stem)
    except ValueError:
        title = src.stem
    tex = "\n".join(
        [
            PREAMBLE + rf"\title{{{inline(title, refs)}}}",
            r"\date{}",
            r"\begin{document}",
            r"\maketitle",
            *body,
            *refs.bibliography(),
            r"\end{document}",
            "",
        ]
    )
    path = out / (src.stem + ".tex")
    path.write_text(tex, encoding="utf-8")
    return path
