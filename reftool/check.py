from __future__ import annotations

import argparse
import re
import tomllib
from pathlib import Path

from reftool.loaders import load_md_entries, load_md_list, load_refs_table, load_slugs_tsv
from reftool.paths import ARTICLES_TSV, DEFAULT_OUTPUT, MASTER_PATH, MD_TSV, REFS_DIR, ROOT, SLUGS_TSV, SLUG_RE
from reftool.toml_io import toml_string, toml_value


def find_slug_mismatches(
    md_entries: list[tuple[Path, str]],
) -> list[tuple[Path, str, list[str], list[str]]]:
    """For each published article (has refs/{ID}.toml), return (md_path,
    article_id, missing, unused) where missing = slugs used in the body but
    not defined in the toml, unused = toml keys not used in the body.
    Only entries with at least one mismatch are returned, sorted by md path.
    """
    result: list[tuple[Path, str, list[str], list[str]]] = []
    for md_path, article_id in md_entries:
        if not article_id:
            continue
        refs_path = REFS_DIR / f"{article_id}.toml"
        if not refs_path.is_file():
            continue
        refs_table = load_refs_table(article_id)
        text = md_path.read_text(encoding="utf-8")
        seen: set[str] = set()
        used: list[str] = []
        for slug in SLUG_RE.findall(text):
            if slug not in seen:
                seen.add(slug)
                used.append(slug)
        used_set = set(used)
        defined_set = set(refs_table.keys())
        missing = sorted(used_set - defined_set)
        unused = sorted(defined_set - used_set)
        if missing or unused:
            result.append((md_path, article_id, missing, unused))
    result.sort(key=lambda e: e[0].relative_to(ROOT).as_posix())
    return result


def render_mismatches(mismatches: list[tuple[Path, str, list[str], list[str]]]) -> str:
    lines: list[str] = []
    for md_path, article_id, missing, unused in mismatches:
        lines.append(f"== {md_path.relative_to(ROOT).as_posix()} ({article_id}) ==")
        if missing:
            lines.append(f"  未定義(本文にあるがtomlにない): {missing}")
        if unused:
            lines.append(f"  未使用(tomlにあるが本文にない): {unused}")
    return "\n".join(lines)


def find_missing_from_slugs_tsv(md_list: list[str], slugs_map: dict[str, str]) -> list[str]:
    """Return md.tsv paths that have no entry in slugs.tsv at all (not even
    "NONE"). slugs.tsv is meant to cover every file in md.tsv, so a missing
    entry means a new article hasn't been registered there yet."""
    return [md for md in md_list if md not in slugs_map]


def find_undefined_in_unpublished(
    md_entries: list[tuple[Path, str]], slugs_map: dict[str, str]
) -> list[tuple[Path, list[str]]]:
    """For each unpublished article (no Mathlog url, so no refs/{ID}.toml of
    its own), find [[slug]] markers in its body that aren't defined anywhere
    in the project — neither in any published article's refs/{ID}.toml, nor
    in refs-master.toml, nor as another article's reserved canonical slug in
    slugs.tsv. Unpublished articles have no refs/{ID}.toml of their own to
    check against, so this checks against every slug already known
    project-wide instead; it only catches typos and not-yet-registered
    references, not "unused" slugs (there is no per-file toml to compare
    against for that direction). slugs.tsv lets a [[slug]] citing a still-
    unpublished sibling article resolve early, via that sibling's reserved
    slug, rather than waiting for it to be published.
    """
    known: set[str] = set()
    for _, article_id in md_entries:
        if article_id:
            known |= load_refs_table(article_id).keys()
    if MASTER_PATH.is_file():
        with MASTER_PATH.open("rb") as f:
            known |= tomllib.load(f).keys()
    known |= {slug for slug in slugs_map.values() if slug != "NONE"}

    result: list[tuple[Path, list[str]]] = []
    for md_path, article_id in md_entries:
        if article_id:
            continue
        text = md_path.read_text(encoding="utf-8")
        seen: set[str] = set()
        used: list[str] = []
        for slug in SLUG_RE.findall(text):
            if slug not in seen:
                seen.add(slug)
                used.append(slug)
        undefined = sorted(set(used) - known)
        if undefined:
            result.append((md_path, undefined))
    result.sort(key=lambda e: e[0].relative_to(ROOT).as_posix())
    return result


def render_undefined_in_unpublished(entries: list[tuple[Path, list[str]]]) -> str:
    lines: list[str] = []
    for md_path, undefined in entries:
        lines.append(f"== {md_path.relative_to(ROOT).as_posix()} (未公開) ==")
        lines.append(f"  未定義(本文にあるがどこにも登録されていない): {undefined}")
    return "\n".join(lines)


def extract_title(entry: dict) -> str | None:
    """Derive the bibliographic title from a refs.toml entry's citation, by
    stripping the author prefix and the type-specific trailing fields
    (site/accessed/url for "website", journal/year/pages for "paper",
    publisher/year/pages for "book"). Returns None if the citation doesn't
    match the expected shape for its type. Used to check that refs.toml's
    citation agrees with refs-master.toml's hand-written "title" field.
    """
    citation = entry.get("citation")
    if citation is None:
        return None
    rem = citation
    author = entry.get("author")
    if author:
        prefix = ", ".join(author) + ", "
        if not rem.startswith(prefix):
            return None
        rem = rem[len(prefix) :]

    t = entry.get("type")
    suffix = ""
    if t == "website":
        tail_parts = []
        site = entry.get("site")
        if site:
            tail_parts.append(site)
        accessed = entry.get("accessed")
        if accessed is not None:
            tail_parts.append(f"閲覧日 {accessed.year}年{accessed.month}月{accessed.day}日")
        url = entry.get("url")
        if url:
            tail_parts.append(url)
        if tail_parts:
            suffix = ", " + ", ".join(tail_parts)
    elif t == "paper":
        journal = entry.get("journal")
        year = entry.get("year")
        tail_parts = [journal, str(year)]
        pages = entry.get("pages")
        if pages:
            tail_parts.append(str(pages))
        suffix = ", " + ", ".join(str(p) for p in tail_parts)
    elif t == "book":
        tail_parts = [entry.get("publisher"), str(entry.get("year")), str(entry.get("pages"))]
        suffix = ", " + ", ".join(str(p) for p in tail_parts)

    if suffix:
        if not rem.endswith(suffix):
            return None
        rem = rem[: -len(suffix)]
    return rem


# Fields present in refs.toml but not compared against refs-master.toml:
# "citation"/"accessed" have no counterpart there (replaced by "title"), and
# "files" is not required to stay in sync since refs-master.toml is
# hand-maintained rather than rebuilt from the articles.
MASTER_SKIP_KEYS = {"citation", "accessed", "files"}


def find_master_mismatches(refs: dict[str, dict], master: dict[str, dict]) -> list[str]:
    """Compare refs.toml against the hand-maintained refs-master.toml (the
    non-mathlog subset, keyed by the same slugs). Reports slugs missing from
    refs.toml, titles implied by refs.toml's citation that disagree with
    refs-master.toml's "title", and any other shared field that disagrees.

    A slug that refs.toml only knows as a body reference (an entry with just
    "files", no "citation") isn't resolved from any refs/{ID}.toml yet — e.g.
    it's cited only from an unpublished article. There's nothing to compare
    yet in that case, so it's skipped rather than reported as a mismatch on
    every field.
    """
    lines: list[str] = []
    for slug, m_entry in sorted(master.items()):
        r_entry = refs.get(slug)
        if r_entry is None:
            lines.append(f"master mismatch: {slug!r} defined in refs-master.toml but missing from refs.toml")
            continue
        if "citation" not in r_entry:
            continue

        title = extract_title(r_entry)
        m_title = m_entry.get("title")
        if title != m_title:
            lines.append(
                f"master mismatch: {slug!r} title differs: "
                f"refs.toml citation implies {title!r}, refs-master.toml has {m_title!r}"
            )

        for key, m_value in m_entry.items():
            if key in MASTER_SKIP_KEYS or key == "title":
                continue
            r_value = r_entry.get(key)
            if r_value != m_value:
                lines.append(
                    f"master mismatch: {slug!r} field {key!r} differs: "
                    f"refs.toml={r_value!r}, refs-master.toml={m_value!r}"
                )
    return lines


def sync_master_files(
    text: str, refs: dict[str, dict], master: dict[str, dict]
) -> tuple[str, list[str]]:
    """Rewrite each slug's "files = [...]" block in refs-master.toml's raw
    text to match refs.toml's "files" list for that slug, leaving every other
    field and the file's original slug order/formatting untouched. Slugs
    missing from refs.toml, or with no "files" key there, are left as-is
    (reported separately by find_master_mismatches). Returns the updated text
    and the list of slugs whose files list actually changed.
    """
    changed: list[str] = []
    for slug, m_entry in master.items():
        r_entry = refs.get(slug)
        if r_entry is None or "files" not in r_entry:
            continue
        new_files = r_entry["files"]
        if new_files == m_entry.get("files"):
            continue
        changed.append(slug)
        items = ",\n".join(f"  {toml_string(f)}" for f in new_files)
        new_block = f"files = [\n{items},\n]"
        pattern = re.compile(r'(\["' + re.escape(slug) + r'"\].*?)files = \[.*?\]', re.DOTALL)
        text = pattern.sub(lambda m: m.group(1) + new_block, text, count=1)
    return text, changed


def check_command(args: argparse.Namespace, sync: bool = False) -> None:
    md_entries = load_md_entries(ARTICLES_TSV)
    mismatches = find_slug_mismatches(md_entries)
    if not mismatches:
        print("mismatch nothing: all published articles match their refs/{ID}.toml")
    else:
        print(render_mismatches(mismatches))

    md_list = load_md_list(MD_TSV)
    slugs_map = load_slugs_tsv(SLUGS_TSV)
    missing_registrations = find_missing_from_slugs_tsv(md_list, slugs_map)
    if missing_registrations:
        print(f"未登録(md.tsvにあるがslugs.tsvにない): {missing_registrations}")
    else:
        print("registration nothing: every md.tsv file is registered in slugs.tsv")

    undefined = find_undefined_in_unpublished(md_entries, slugs_map)
    if undefined:
        print(render_undefined_in_unpublished(undefined))
    else:
        print("undefined nothing: all [[slug]] markers in unpublished articles are known project-wide")

    if not MASTER_PATH.is_file():
        return
    with DEFAULT_OUTPUT.open("rb") as f:
        refs = tomllib.load(f)
    with MASTER_PATH.open("rb") as f:
        master = tomllib.load(f)

    if sync:
        text = MASTER_PATH.read_text(encoding="utf-8")
        new_text, changed = sync_master_files(text, refs, master)
        if changed:
            MASTER_PATH.write_text(new_text, encoding="utf-8")
            print(f"synced files for {len(changed)} slug(s) in {MASTER_PATH.relative_to(ROOT)}: {changed}")
            with MASTER_PATH.open("rb") as f:
                master = tomllib.load(f)
        else:
            print("sync nothing: refs-master.toml files already match refs.toml")

    master_mismatches = find_master_mismatches(refs, master)
    if master_mismatches:
        print("\n".join(master_mismatches))
    else:
        print("master mismatch nothing: refs.toml agrees with refs-master.toml")


def sync_command(args: argparse.Namespace) -> None:
    check_command(args, sync=True)
