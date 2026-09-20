from __future__ import annotations

import argparse
from pathlib import Path

from articles.merge import iter_tsv_data_lines
from articles.paths import ARTICLES_TSV, ROOT


def load_articles_tsv(path: Path) -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []
    for line in iter_tsv_data_lines(path):
        date_, url, md, title = line.split("\t", 3)
        rows.append((date_, url, md, title))
    return rows


def find_pending(rows: list[tuple[str, str, str, str]]) -> list[tuple[str, str]]:
    """Return (md, title) for rows that are unpublished (no url) and
    unreviewed (no matching .txt review trace)."""
    pending: list[tuple[str, str]] = []
    for date_, url, md, title in rows:
        if not md or url:
            continue
        if (ROOT / md).with_suffix(".txt").exists():
            continue
        pending.append((md, title))
    return pending


def add_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "pending", help="articles.tsv → 未公開かつ未レビューの記事一覧"
    )
    parser.set_defaults(func=pending_command)


def pending_command(args: argparse.Namespace) -> None:
    if not ARTICLES_TSV.is_file():
        raise SystemExit(f"missing {ARTICLES_TSV}; run: articles merge")

    rows = load_articles_tsv(ARTICLES_TSV)
    pending = find_pending(rows)

    if not pending:
        print("未公開かつ未レビューの記事はありません")
        return

    for md, title in pending:
        print(f"{md}\t{title}")
