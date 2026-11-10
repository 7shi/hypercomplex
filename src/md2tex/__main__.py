"""Mathlog記法のMarkdown記事をLaTeXに変換するツールのエントリポイント。詳細はREADME.mdを参照。"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from md2tex.convert import ROOT, convert


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip())
    parser.add_argument("article", type=Path, help="変換する記事（.md）")
    parser.add_argument(
        "-o", "--out", type=Path, help="出力先ディレクトリ（既定：build/<記事名>/）"
    )
    parser.add_argument("--pdf", action="store_true", help="lualatexでPDFまで作る")
    args = parser.parse_args(argv)

    src = args.article.resolve()
    out = args.out or ROOT / "build" / src.stem
    tex = convert(src, out)
    print(tex)
    if args.pdf:
        for _ in range(2):  # 2回目でPDFのしおりが入る
            subprocess.run(
                ["lualatex", "-interaction=nonstopmode", "-halt-on-error", tex.name],
                cwd=out,
                check=True,
                stdout=subprocess.DEVNULL,
            )
        print(tex.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
