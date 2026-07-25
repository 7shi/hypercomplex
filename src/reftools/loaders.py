from __future__ import annotations

import tomllib
from pathlib import Path

from reftools.paths import MATHLOG_BASE, REFS_DIR, ROOT


def load_articles(path: Path) -> list[tuple[str, str, str]]:
    """Parse articles.tsv into (url, md, title) rows, in file order."""
    lines = path.read_text(encoding="utf-8").splitlines()
    rows: list[tuple[str, str, str]] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        url = cols[1] if len(cols) > 1 else ""
        md = cols[2] if len(cols) > 2 else ""
        title = cols[3] if len(cols) > 3 else ""
        if md:
            rows.append((url, md, title))
    return rows


def load_md_entries(path: Path) -> list[tuple[Path, str]]:
    """Return (md_path, article_id) pairs; article_id is '' if none."""
    return [
        (ROOT / md, url.rsplit("/", 1)[-1] if url else "")
        for url, md, _ in load_articles(path)
    ]


def load_md_list(path: Path) -> list[str]:
    """Return every md path listed in md.tsv (the full article list, published
    or not), in file order."""
    lines = path.read_text(encoding="utf-8").splitlines()
    result: list[str] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        md = line.split("\t")[0]
        if md:
            result.append(md)
    return result


def load_slugs_tsv(path: Path) -> dict[str, str]:
    """Return {md_path: slug} from slugs.tsv, the hand-maintained ledger of
    each article's canonical slug (used by other articles to cite it via
    [[slug]] even before it's published). "NONE" is kept as-is — it marks a
    file with no publication planned."""
    lines = path.read_text(encoding="utf-8").splitlines()
    result: dict[str, str] = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 2:
            continue
        result[cols[0]] = cols[1]
    return result


def load_md_urls(path: Path) -> dict[str, str]:
    """Map md path (relative, as in articles.tsv) -> full Mathlog url, for
    every published article in articles.tsv."""
    return {
        md: MATHLOG_BASE + url
        for url, md, _ in load_articles(path)
        if url
    }


def load_md_titles(path: Path) -> dict[str, str]:
    """Map md path (relative, as in articles.tsv) -> article title."""
    return {
        md: title
        for _, md, title in load_articles(path)
        if title
    }


def load_title_urls(path: Path) -> dict[str, str]:
    """Map article title -> full Mathlog url, from articles.tsv."""
    return {
        title: MATHLOG_BASE + url
        for url, _, title in load_articles(path)
        if url and title
    }


def load_refs_table(article_id: str) -> dict[str, dict]:
    refs_path = REFS_DIR / f"{article_id}.toml"
    if not refs_path.is_file():
        return {}
    with refs_path.open("rb") as f:
        return tomllib.load(f)
