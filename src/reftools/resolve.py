from __future__ import annotations

from pathlib import Path

from reftools.loaders import load_refs_table
from reftools.paths import MATHLOG_BASE, MATHLOG_CITATION_RE, ROOT, SLUG_RE


def resolve_slug(
    slug: str, refs_table: dict[str, dict], title_urls: dict[str, str]
) -> tuple[bool, dict]:
    """Return (defined, entry) for slug, with all attributes from refs/{ID}.toml.

    For type "mathlog", the citation's title is looked up in articles.tsv and
    the result is added as an explicit "url" key (refs/{ID}.toml itself has no
    url for mathlog entries).
    """
    entry = refs_table.get(slug)
    if entry is None:
        return False, {}
    entry = dict(entry)
    if entry.get("type") == "mathlog":
        m = MATHLOG_CITATION_RE.match(entry.get("citation", ""))
        if m:
            title = m.group(2)
            url = title_urls.get(title, "")
            if url:
                entry["url"] = url
    return True, entry


def collect_slugs_by_file(
    md_entries: list[tuple[Path, str]], title_urls: dict[str, str]
) -> list[tuple[Path, str, list[tuple[str, bool, dict]]]]:
    """Return (md_path, article_url, [(slug, defined, entry), ...]) per file."""
    result: list[tuple[Path, str, list[tuple[str, bool, dict]]]] = []
    for md_path, article_id in md_entries:
        text = md_path.read_text(encoding="utf-8")
        seen: set[str] = set()
        slugs: list[str] = []
        for slug in SLUG_RE.findall(text):
            if slug not in seen:
                seen.add(slug)
                slugs.append(slug)
        if not slugs:
            continue
        refs_table = load_refs_table(article_id) if article_id else {}
        resolved = [(slug, *resolve_slug(slug, refs_table, title_urls)) for slug in slugs]
        article_url = f"{MATHLOG_BASE}/articles/{article_id}" if article_id else ""
        result.append((md_path, article_url, resolved))
    result.sort(key=lambda e: e[0].relative_to(ROOT).as_posix())
    return result
