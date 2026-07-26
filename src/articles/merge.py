from __future__ import annotations

import argparse
from pathlib import Path

from articles.paths import ARTICLES_TSV, MATHLOG_TSV, MD_TSV, ROOT, write_tsv


def iter_tsv_data_lines(path: Path):
    """Yield non-empty TSV lines, skipping a header row if present."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return
    start = 0
    first_col = lines[0].split("\t", 1)[0]
    if first_col in {"date", "md", "url", "title", "path"}:
        start = 1
    for line in lines[start:]:
        if line.strip():
            yield line


def load_mathlog_tsv(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    for line in iter_tsv_data_lines(path):
        date_, url, title = line.split("\t", 2)
        rows.append((date_, url, title))
    return rows


def load_md_tsv(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for line in iter_tsv_data_lines(path):
        p, title = line.split("\t", 1)
        rows.append((p, title))
    return rows


def match_md_index(mathlog_title: str, md_rows: list[tuple[str, str]], used: set[int]) -> int | None:
    """Return index into md_rows matching mathlog_title, or None.

    Prefer exact title match; fall back to mathlog_title.endswith(md_title),
    choosing the longest md title among unused rows.
    """
    exact: int | None = None
    for i, (_, md_title) in enumerate(md_rows):
        if i in used:
            continue
        if mathlog_title == md_title:
            exact = i
            break
    if exact is not None:
        return exact

    best_i: int | None = None
    best_len = -1
    for i, (_, md_title) in enumerate(md_rows):
        if i in used or not md_title:
            continue
        if mathlog_title.endswith(md_title) and len(md_title) > best_len:
            best_i = i
            best_len = len(md_title)
    return best_i


def merge_rows(
    mathlog_rows: list[tuple[str, str, str]],
    md_rows: list[tuple[str, str]],
) -> list[tuple[str, str, str, str]]:
    """Join on title. Output columns: date, url, md, title.

    title prefers the Mathlog title; falls back to the md title when
    Mathlog has no row. Unmatched sides leave the missing fields empty.
    """
    used: set[int] = set()
    out: list[tuple[str, str, str, str]] = []

    for date_, url, title in mathlog_rows:
        j = match_md_index(title, md_rows, used)
        if j is None:
            out.append((date_, url, "", title))
        else:
            used.add(j)
            md_path, _ = md_rows[j]
            out.append((date_, url, md_path, title))

    for i, (md_path, title) in enumerate(md_rows):
        if i not in used:
            out.append(("", "", md_path, title))

    return out


def add_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("merge", help="mathlog.tsv + md.tsv → articles.tsv (by title)")
    parser.set_defaults(func=merge_command)


def merge_command(args: argparse.Namespace) -> None:
    if not MATHLOG_TSV.is_file():
        raise SystemExit(f"missing {MATHLOG_TSV}; run: articles mathlog")
    if not MD_TSV.is_file():
        raise SystemExit(f"missing {MD_TSV}; run: articles md")

    mathlog_rows = load_mathlog_tsv(MATHLOG_TSV)
    md_rows = load_md_tsv(MD_TSV)
    rows = merge_rows(mathlog_rows, md_rows)

    lines = [f"{d}\t{url}\t{md}\t{title}" for d, url, md, title in rows]
    write_tsv(ARTICLES_TSV, "date\turl\tmd\ttitle", lines)

    matched = sum(1 for d, url, md, title in rows if d and md)
    only_ml = sum(1 for d, url, md, title in rows if d and not md)
    only_md = sum(1 for d, url, md, title in rows if md and not d)
    print(
        f"wrote {len(rows)} rows to {ARTICLES_TSV.relative_to(ROOT)} "
        f"(matched={matched}, mathlog_only={only_ml}, md_only={only_md})"
    )
