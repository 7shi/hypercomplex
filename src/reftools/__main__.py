"""[[slug]]参照マーカーをrefs.tomlに集約するツールのエントリポイント。詳細はREADME.mdを参照。"""

from __future__ import annotations

import argparse

from reftools import build, check, html_format, refs_to_toml, show


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip())
    subparsers = parser.add_subparsers(dest="command", required=True)

    build.add_subparser(subparsers)
    check.add_subparser(subparsers)
    show.add_subparser(subparsers)
    html_format.add_subparser(subparsers)
    refs_to_toml.add_subparser(subparsers)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
