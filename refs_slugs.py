"""Collect [[slug]] reference markers from articles.tsv's md files into one TOML file.

For each md file (sorted by path) that has a Mathlog url in articles.tsv, its
per-article reference list refs/{ID}.toml is read and each [[slug]] marker
found in the md file is resolved to a type and url:
  - non-mathlog types: url is taken directly from refs/{ID}.toml
  - type "mathlog": the citation's title is looked up in articles.tsv to
    find the referenced article's own url
  - slugs with no matching key in refs/{ID}.toml are listed separately
    under that file's "undefined" field instead of "refs"

Md files without a Mathlog url (and thus no refs/{ID}.toml) fall back to a
plain slugs list. Files with no [[slug]] markers are omitted.

Also writes refs-url.txt: urls referenced under more than one distinct slug
(the same source cited under different slug names across articles), sorted
by url, one line each as "url = slug, slug, ...".

Also writes refs-file.txt: for each md file that is itself a Mathlog article,
the distinct slugs other articles use to cite it via type = "mathlog"
(useful for spotting slug drift and picking one canonical slug per target
article), sorted by md path, one line each as "md_path = slug, slug, ..." or
"md_path = NONE" when no article cites it.

Usage:
  refs_slugs.py                      # write refs.toml, refs-url.txt, refs-file.txt
  refs_slugs.py -o out.toml
"""

from __future__ import annotations

import argparse
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ARTICLES_TSV = ROOT / "articles.tsv"
REFS_DIR = ROOT / "refs"
DEFAULT_OUTPUT = ROOT / "refs.toml"
DEFAULT_URL_OUTPUT = ROOT / "refs-url.txt"
DEFAULT_FILE_OUTPUT = ROOT / "refs-file.txt"
MATHLOG_BASE = "https://mathlog.info"

SLUG_RE = re.compile(r"\[\[([^\]]+)\]\]")
MATHLOG_CITATION_RE = re.compile(r"^([^,]+), (.*), Mathlog, ")


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


def resolve_slug(
    slug: str, refs_table: dict[str, dict], title_urls: dict[str, str]
) -> tuple[bool, str, str]:
    """Return (defined, type, url) for slug."""
    entry = refs_table.get(slug)
    if entry is None:
        return False, "", ""
    ref_type = entry.get("type", "")
    if ref_type == "mathlog":
        m = MATHLOG_CITATION_RE.match(entry.get("citation", ""))
        if m:
            title = m.group(2)
            return True, ref_type, title_urls.get(title, "")
        return True, ref_type, ""
    return True, ref_type, entry.get("url", "")


def collect_slugs_by_file(
    md_entries: list[tuple[Path, str]], title_urls: dict[str, str]
) -> list[tuple[Path, str, list[tuple[str, bool, str, str]]]]:
    """Return (md_path, article_url, [(slug, defined, type, url), ...]) per file."""
    result: list[tuple[Path, str, list[tuple[str, bool, str, str]]]] = []
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


def toml_string(s: str) -> str:
    escaped = s.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def render_toml(entries: list[tuple[Path, str, list[tuple[str, bool, str, str]]]]) -> str:
    blocks: list[str] = []
    for md_path, article_url, refs in entries:
        key = toml_string(md_path.relative_to(ROOT).as_posix())
        lines = [f"[{key}]"]
        if article_url:
            lines.append(f"url = {toml_string(article_url)}")
            undefined = [slug for slug, defined, _, _ in refs if not defined]
            if undefined:
                items = ", ".join(toml_string(s) for s in undefined)
                lines.append(f"undefined = [{items}]")
        blocks.append("\n".join(lines) + "\n")

        if article_url:
            for slug, defined, ref_type, url in refs:
                if not defined:
                    continue
                r_lines = [f"[[{key}.refs]]", f"slug = {toml_string(slug)}"]
                if ref_type:
                    r_lines.append(f"type = {toml_string(ref_type)}")
                if url:
                    r_lines.append(f"url = {toml_string(url)}")
                blocks.append("\n".join(r_lines) + "\n")
        else:
            items = ",\n".join(f"  {toml_string(slug)}" for slug, _, _, _ in refs)
            blocks.append(f"slugs = [\n{items},\n]\n")
    return "\n".join(blocks)


def find_shared_urls(entries: list[tuple[Path, str, list[tuple[str, bool, str, str]]]]) -> list[tuple[str, list[str]]]:
    """Return (url, slugs) for urls cited under more than one distinct slug."""
    slugs_by_url: dict[str, list[str]] = {}
    for _, _, refs in entries:
        for slug, _, _, url in refs:
            if not url:
                continue
            slugs = slugs_by_url.setdefault(url, [])
            if slug not in slugs:
                slugs.append(slug)
    shared = [(url, slugs) for url, slugs in slugs_by_url.items() if len(slugs) > 1]
    shared.sort(key=lambda e: e[0])
    return shared


def render_shared_urls(shared: list[tuple[str, list[str]]]) -> str:
    return "".join(f"{url} = {', '.join(slugs)}\n" for url, slugs in shared)


def find_citing_slugs_by_file(
    md_entries: list[tuple[Path, str]],
    entries: list[tuple[Path, str, list[tuple[str, bool, str, str]]]],
) -> list[tuple[Path, list[str]]]:
    """For each md file with an article id, the distinct slugs other
    articles use to cite it via type = "mathlog", sorted by md path."""
    id_to_path = {article_id: md_path for md_path, article_id in md_entries if article_id}
    slugs_by_id: dict[str, list[str]] = {aid: [] for aid in id_to_path}

    for _, _, refs in entries:
        for slug, _, ref_type, url in refs:
            if ref_type != "mathlog" or not url:
                continue
            article_id = url.rsplit("/", 1)[-1]
            if article_id not in slugs_by_id:
                continue
            slugs = slugs_by_id[article_id]
            if slug not in slugs:
                slugs.append(slug)

    result = [(md_path, slugs_by_id[aid]) for aid, md_path in id_to_path.items()]
    result.sort(key=lambda e: e[0].relative_to(ROOT).as_posix())
    return result


def render_citing_slugs(citing: list[tuple[Path, list[str]]]) -> str:
    lines = []
    for md_path, slugs in citing:
        key = md_path.relative_to(ROOT).as_posix()
        value = ", ".join(slugs) if slugs else "NONE"
        lines.append(f"{key} = {value}\n")
    return "".join(lines)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--url-output", type=Path, default=DEFAULT_URL_OUTPUT)
    parser.add_argument("--file-output", type=Path, default=DEFAULT_FILE_OUTPUT)
    args = parser.parse_args(argv)

    md_entries = load_md_entries(ARTICLES_TSV)
    title_urls = load_title_urls(ARTICLES_TSV)
    entries = collect_slugs_by_file(md_entries, title_urls)
    total = sum(len(refs) for _, _, refs in entries)

    args.output.write_text(render_toml(entries), encoding="utf-8")
    print(f"wrote {total} slugs across {len(entries)} files to {args.output.relative_to(ROOT)}")

    shared = find_shared_urls(entries)
    args.url_output.write_text(render_shared_urls(shared), encoding="utf-8")
    print(f"wrote {len(shared)} shared urls to {args.url_output.relative_to(ROOT)}")

    citing = find_citing_slugs_by_file(md_entries, entries)
    args.file_output.write_text(render_citing_slugs(citing), encoding="utf-8")
    print(f"wrote {len(citing)} files to {args.file_output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
