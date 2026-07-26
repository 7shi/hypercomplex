"""Mathlog記事とREADMEリンクを突き合わせるツールのエントリポイント。詳細はREADME.mdを参照。"""

from __future__ import annotations

import argparse

from articles import md, mathlog, merge


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip())
    subparsers = parser.add_subparsers(dest="command", required=True)

    mathlog.add_subparser(subparsers)
    md.add_subparser(subparsers)
    merge.add_subparser(subparsers)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
