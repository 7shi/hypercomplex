from __future__ import annotations

import argparse
import tomllib
from pathlib import Path

from reftools.loaders import load_md_urls, load_slugs_tsv
from reftools.paths import ARTICLES_TSV, DEFAULT_OUTPUT, MASTER_PATH, ROOT, SLUGS_TSV, SLUG_RE


def resolve_show_target(target: str, slugs_map: dict[str, str]) -> Path:
    """Resolve the "show" argument to an md path: if it names an existing
    file, use it directly; otherwise treat it as a canonical slug and look it
    up in slugs.tsv (md_path -> slug)."""
    candidate = ROOT / target
    if candidate.is_file():
        return candidate
    for md, slug in slugs_map.items():
        if slug == target:
            return ROOT / md
    raise SystemExit(f"error: {target!r} はファイルにもslugs.tsvのslugにも該当しません")


def extract_slugs_in_order(md_path: Path) -> list[str]:
    text = md_path.read_text(encoding="utf-8")
    seen: set[str] = set()
    result: list[str] = []
    for slug in SLUG_RE.findall(text):
        if slug not in seen:
            seen.add(slug)
            result.append(slug)
    return result


def has_resolved_info(entry: dict) -> bool:
    """True if a refs.toml entry carries anything beyond "files" (i.e. it was
    actually resolved from some refs/{ID}.toml, not just referenced)."""
    return any(k != "files" for k in entry)


def display_value(v) -> str:
    if isinstance(v, list):
        return ", ".join(display_value(x) for x in v)
    return str(v)


def render_entry(entry: dict) -> str:
    lines = [f"  {key}: {display_value(value)}" for key, value in entry.items() if key != "files"]
    return "\n".join(lines)


def build_slugs_reverse_map(slugs_map: dict[str, str]) -> dict[str, str]:
    """Return {slug: md_path} from slugs.tsv's {md_path: slug}, excluding
    "NONE" entries (files with no publication planned, so no canonical
    slug to reverse-lookup)."""
    return {slug: md for md, slug in slugs_map.items() if slug != "NONE"}


def add_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "show",
        help="show [[slug]] markers in one article, in order, with their refs.toml/refs-master.toml info",
    )
    parser.add_argument("target", help="md path or canonical slug (slugs.tsv)")
    parser.set_defaults(func=show_command)


def show_command(args: argparse.Namespace) -> None:
    slugs_map = load_slugs_tsv(SLUGS_TSV)
    slugs_reverse = build_slugs_reverse_map(slugs_map)
    md_urls = load_md_urls(ARTICLES_TSV)
    md_path = resolve_show_target(args.target, slugs_map)
    md_rel = md_path.relative_to(ROOT).as_posix()
    slugs = extract_slugs_in_order(md_path)

    master: dict[str, dict] = {}
    if MASTER_PATH.is_file():
        with MASTER_PATH.open("rb") as f:
            master = tomllib.load(f)
    refs: dict[str, dict] = {}
    if DEFAULT_OUTPUT.is_file():
        with DEFAULT_OUTPUT.open("rb") as f:
            refs = tomllib.load(f)

    own_slug = slugs_map.get(md_rel)
    if own_slug is None:
        print(f"{md_rel} (slugs.tsvに未登録、登録してください)")
    else:
        print(f"{md_rel} ({own_slug})")
    print(f"({len(slugs)} slugs)")
    for i, slug in enumerate(slugs, 1):
        print()
        if slug in master:
            print(f"{i}. [{slug}] (refs-master.toml)")
            print(render_entry(master[slug]))
        elif slug in refs and has_resolved_info(refs[slug]):
            print(f"{i}. [{slug}] (refs.toml)")
            print(render_entry(refs[slug]))
        elif slug in slugs_reverse:
            target_md = slugs_reverse[slug]
            print(f"{i}. [{slug}] (slugs.tsv)")
            print(f"  {target_md}")
            url = md_urls.get(target_md)
            if url:
                print(f"  url: {url}")
        else:
            print(f"{i}. [{slug}] 情報なし")
