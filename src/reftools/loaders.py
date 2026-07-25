from __future__ import annotations

import tomllib
from pathlib import Path

from reftools.paths import MATHLOG_BASE, REFS_DIR, ROOT


def load_md_entries(path: Path) -> list[tuple[Path, str]]:
    """Return (md_path, article_id) pairs; article_id is '' if none."""
    lines = path.read_text(encoding="utf-8").splitlines()
    entries: list[tuple[Path, str]] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        url = cols[1] if len(cols) > 1 else ""
        md = cols[2] if len(cols) > 2 else ""
        if md:
            article_id = url.rsplit("/", 1)[-1] if url else ""
            entries.append((ROOT / md, article_id))
    return entries


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
    lines = path.read_text(encoding="utf-8").splitlines()
    result: dict[str, str] = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 3:
            continue
        url, md = cols[1], cols[2]
        if url and md:
            result[md] = MATHLOG_BASE + url
    return result


def load_title_urls(path: Path) -> dict[str, str]:
    """Map article title -> full Mathlog url, from articles.tsv."""
    lines = path.read_text(encoding="utf-8").splitlines()
    result: dict[str, str] = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        _, url, _, title = cols[0], cols[1], cols[2], cols[3]
        if url and title:
            result[title] = MATHLOG_BASE + url
    return result


def load_refs_table(article_id: str) -> dict[str, dict]:
    refs_path = REFS_DIR / f"{article_id}.toml"
    if not refs_path.is_file():
        return {}
    with refs_path.open("rb") as f:
        return tomllib.load(f)
