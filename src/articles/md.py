from __future__ import annotations

import argparse
import re
from pathlib import Path

from articles.paths import MD_TSV, ROOT, write_tsv

LINK_RE = re.compile(r"\[([^\]]+)\]\(([^)]+\.md)\)")


def extract_md() -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    seen: set[str] = set()

    for readme in sorted(ROOT.rglob("README.md")):
        text = readme.read_text(encoding="utf-8")
        base = readme.parent
        for title, href in LINK_RE.findall(text):
            href_path = href.split("#", 1)[0]
            if not href_path.endswith(".md"):
                continue
            resolved = (base / href_path).resolve()
            try:
                rel = resolved.relative_to(ROOT)
            except ValueError:
                continue
            key = rel.as_posix()
            stem = Path(key).stem
            if stem.isupper() or key in seen:
                continue
            seen.add(key)
            rows.append((key, title.strip()))

    rows.sort(key=lambda r: r[0])
    return rows


def add_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("md", help="README.md files → md.tsv")
    parser.set_defaults(func=md_command)


def md_command(args: argparse.Namespace) -> None:
    rows = extract_md()
    if not rows:
        raise SystemExit("no .md article links found in README.md files")

    lines = [f"{path}\t{title}" for path, title in rows]
    write_tsv(MD_TSV, "md\ttitle", lines)
    print(f"wrote {len(rows)} rows to {MD_TSV.relative_to(ROOT)}")
