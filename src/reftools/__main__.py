"""[[slug]]参照マーカーをrefs.tomlに集約するツールのエントリポイント。詳細はREADME.mdを参照。"""

from __future__ import annotations

import argparse
from pathlib import Path

from reftools.build import build_command
from reftools.check import check_command, sync_command
from reftools.html_format import html_format_command
from reftools.paths import DEFAULT_OUTPUT
from reftools.refs_to_toml import refs_to_toml_command
from reftools.show import show_command


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip())
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser("build", help="write refs.toml (and optionally refs-url.txt/refs-file.txt)")
    build_parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT)
    build_parser.add_argument("--url-output", type=Path, default=None)
    build_parser.add_argument("--file-output", type=Path, default=None)
    build_parser.set_defaults(func=build_command)

    check_parser = subparsers.add_parser(
        "check",
        help="check body [[slug]] markers against refs/{ID}.toml keys, "
        "unpublished articles against project-wide known slugs, "
        "and refs.toml against refs-master.toml",
    )
    check_parser.set_defaults(func=check_command)

    sync_parser = subparsers.add_parser(
        "sync",
        help="same checks as \"check\", plus sync each slug's \"files\" list in "
        "refs-master.toml to match refs.toml",
    )
    sync_parser.set_defaults(func=sync_command)

    show_parser = subparsers.add_parser(
        "show",
        help="show [[slug]] markers in one article, in order, with their refs.toml/refs-master.toml info",
    )
    show_parser.add_argument("target", help="md path or canonical slug (slugs.tsv)")
    show_parser.set_defaults(func=show_command)

    html_format_parser = subparsers.add_parser(
        "format",
        help="pretty-print a Mathlog refs-panel HTML export with one tag per line",
    )
    html_format_parser.add_argument("input", type=Path, help="HTML file to format")
    html_format_group = html_format_parser.add_mutually_exclusive_group()
    html_format_group.add_argument("-o", "--output", type=Path, help="write result to this file")
    html_format_group.add_argument("--in-place", action="store_true", help="overwrite the input file")
    html_format_parser.set_defaults(func=html_format_command)

    refs_to_toml_parser = subparsers.add_parser(
        "toml",
        help="convert a formatted refs-panel HTML export (refs/*.html) to refs/*.toml",
    )
    refs_to_toml_parser.add_argument("input", type=Path, help="formatted refs HTML file")
    refs_to_toml_parser.add_argument("-o", "--output", type=Path, help="write result to this file")
    refs_to_toml_parser.set_defaults(func=refs_to_toml_command)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
