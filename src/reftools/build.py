from __future__ import annotations

import argparse
from pathlib import Path

from reftools.loaders import load_md_entries, load_title_urls
from reftools.paths import ARTICLES_TSV, ROOT
from reftools.resolve import collect_slugs_by_file
from reftools.toml_io import render_toml_by_slug

# Attributes that legitimately differ across refs/{ID}.toml exports of the
# same slug (the access date, and the citation string that embeds it) and so
# are excluded from the cross-file conflict check.
DATE_KEYS = {"citation", "accessed"}


def collect_by_slug(
    entries: list[tuple[Path, str, list[tuple[str, bool, dict]]]],
) -> tuple[dict[str, dict], list[str]]:
    """Group refs by slug into {slug: {**attrs, "files": [...]}}.

    All attributes other than DATE_KEYS must match exactly across every file
    that defines the slug; a mismatch is a naming collision. Every slug is
    checked before reporting, so a conflict in one slug never stops the scan
    of the rest — conflicting slugs are collected into the returned message
    list instead of being added to the result dict. Among the (otherwise
    identical) definitions of a slug with no conflict, the one with the
    latest "accessed" date is kept.
    """
    files_by_slug: dict[str, list[str]] = {}
    defs_by_slug: dict[str, list[tuple[str, dict]]] = {}
    for md_path, _, refs in entries:
        rel = md_path.relative_to(ROOT).as_posix()
        for slug, defined, entry in refs:
            files = files_by_slug.setdefault(slug, [])
            if rel not in files:
                files.append(rel)
            if defined:
                defs_by_slug.setdefault(slug, []).append((rel, entry))

    result: dict[str, dict] = {}
    conflicts: list[str] = []
    for slug, files in files_by_slug.items():
        files.sort()
        defs = defs_by_slug.get(slug, [])
        if not defs:
            result[slug] = {"files": files}
            continue

        base_rel, base_entry = defs[0]
        base_key = {k: v for k, v in base_entry.items() if k not in DATE_KEYS}
        mismatches = [
            (rel, entry)
            for rel, entry in defs[1:]
            if {k: v for k, v in entry.items() if k not in DATE_KEYS} != base_key
        ]
        if mismatches:
            lines = [f"slug conflict: {slug!r} differs across files (ignoring {sorted(DATE_KEYS)}):"]
            lines.append(f"  {base_rel}: {base_entry}")
            for rel, entry in mismatches:
                lines.append(f"  {rel}: {entry}")
            conflicts.append("\n".join(lines))
            # still merge below (latest "accessed" wins) so the build can
            # proceed; the conflict is reported, not silently ignored.

        chosen_rel, chosen_entry = base_rel, base_entry
        chosen_date = chosen_entry.get("accessed")
        for rel, entry in defs[1:]:
            date = entry.get("accessed")
            if date is not None and (chosen_date is None or date > chosen_date):
                chosen_rel, chosen_entry, chosen_date = rel, entry, date

        merged = dict(chosen_entry)
        merged["files"] = files
        result[slug] = merged
    return result, conflicts


def find_shared_urls(entries: list[tuple[Path, str, list[tuple[str, bool, dict]]]]) -> list[tuple[str, list[str]]]:
    """Return (url, slugs) for urls cited under more than one distinct slug."""
    slugs_by_url: dict[str, list[str]] = {}
    for _, _, refs in entries:
        for slug, _, entry in refs:
            url = entry.get("url", "")
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
    entries: list[tuple[Path, str, list[tuple[str, bool, dict]]]],
) -> list[tuple[Path, list[str]]]:
    """For each md file with an article id, the distinct slugs other
    articles use to cite it via type = "mathlog", sorted by md path."""
    id_to_path = {article_id: md_path for md_path, article_id in md_entries if article_id}
    slugs_by_id: dict[str, list[str]] = {aid: [] for aid in id_to_path}

    for _, _, refs in entries:
        for slug, _, entry in refs:
            url = entry.get("url", "")
            if entry.get("type") != "mathlog" or not url:
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


def build_command(args: argparse.Namespace) -> None:
    args.output = args.output.resolve()
    if args.url_output is not None:
        args.url_output = args.url_output.resolve()
    if args.file_output is not None:
        args.file_output = args.file_output.resolve()

    md_entries = load_md_entries(ARTICLES_TSV)
    title_urls = load_title_urls(ARTICLES_TSV)
    entries = collect_slugs_by_file(md_entries, title_urls)
    slug_map, conflicts = collect_by_slug(entries)

    args.output.write_text(render_toml_by_slug(slug_map), encoding="utf-8")
    print(f"wrote {len(slug_map)} slugs to {args.output.relative_to(ROOT)}")

    if args.url_output is not None:
        shared = find_shared_urls(entries)
        args.url_output.write_text(render_shared_urls(shared), encoding="utf-8")
        print(f"wrote {len(shared)} shared urls to {args.url_output.relative_to(ROOT)}")

    if args.file_output is not None:
        citing = find_citing_slugs_by_file(md_entries, entries)
        args.file_output.write_text(render_citing_slugs(citing), encoding="utf-8")
        print(f"wrote {len(citing)} files to {args.file_output.relative_to(ROOT)}")

    if conflicts:
        print(f"\n{len(conflicts)} slug conflict(s) found (output written using latest \"accessed\" per slug anyway):\n")
        print("\n\n".join(conflicts))
        raise SystemExit(1)
