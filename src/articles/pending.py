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


def find_pending(
    rows: list[tuple[str, str, str, str]],
) -> tuple[
    list[tuple[bool, str, str]],
    list[tuple[bool, str, str]],
    list[tuple[bool, str, str]],
]:
    """Return (unreviewed, no_prompt, stale) lists of (published, md, title).

    unreviewed: no matching .txt review trace (takes priority).
    no_prompt: reviewed (.txt exists) but no matching -prompt.md.
    stale: the .txt exists but the md is older than it (review not reflected).
    """
    no_prompt: list[tuple[bool, str, str]] = []
    unreviewed: list[tuple[bool, str, str]] = []
    stale: list[tuple[bool, str, str]] = []
    for date_, url, md, title in rows:
        if not md:
            continue
        md_path = ROOT / md
        txt_path = md_path.with_suffix(".txt")
        prompt_path = md_path.with_name(f"{md_path.stem}-prompt.md")
        entry = (bool(url), md, title)
        if not txt_path.exists():
            unreviewed.append(entry)
        elif not prompt_path.exists():
            no_prompt.append(entry)
        elif md_path.exists() and md_path.stat().st_mtime < txt_path.stat().st_mtime:
            stale.append(entry)
    no_prompt.sort(key=lambda e: e[1])
    unreviewed.sort(key=lambda e: e[1])
    stale.sort(key=lambda e: e[1])
    return unreviewed, no_prompt, stale


def add_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "pending", help="articles.tsv → プロンプト未作成・未レビュー・レビュー未反映の記事一覧（[済]公開済・[未]未公開）"
    )
    parser.set_defaults(func=pending_command)


def pending_command(args: argparse.Namespace) -> None:
    if not ARTICLES_TSV.is_file():
        raise SystemExit(f"missing {ARTICLES_TSV}; run: articles merge")

    rows = load_articles_tsv(ARTICLES_TSV)
    unreviewed, no_prompt, stale = find_pending(rows)

    sections = [
        ("レビューされていない記事", unreviewed),
        ("レビューのプロンプトがない記事", no_prompt),
        ("レビューが反映されていない記事", stale),
    ]
    for i, (heading, entries) in enumerate(sections):
        if i:
            print()
        print(f"# {heading} ({len(entries)})")
        width = len(str(len(entries)))
        for n, (published, md, title) in enumerate(entries, 1):
            mark = "[済]" if published else "[未]"
            print(f"{n:>{width}}. {mark} {md} {title}")
