"""Generate slug-fixes.md: a checklist of old slugs to fix by hand.

Old slugs (see SLUG.md for the decisions) appear in two places:
  - article bodies: [[old-slug]] markers in the local md files, which must
    also be fixed in the article body on Mathlog
  - reference registrations: slug labels in refs/{ID}.toml, which must be
    fixed in the reference panel on Mathlog (then re-exported or edited
    locally to match)

The rename map is embedded below. A few slug strings are ambiguous (the
same string points at different targets in different articles), so those
renames are scoped to a single md file via OVERRIDES; in md files without
an override entry such slugs are reported for manual review.

Usage:
  slug_fixes.py               # write slug-fixes.md
  slug_fixes.py -o out.md
"""

from __future__ import annotations

import argparse
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ARTICLES_TSV = ROOT / "articles.tsv"
REFS_DIR = ROOT / "refs"
DEFAULT_OUTPUT = ROOT / "slug-fixes.md"
MATHLOG_BASE = "https://mathlog.info"

SLUG_RE = re.compile(r"\[\[([^\]]+)\]\]")

# old slug -> canonical slug, applied in every article
RENAMES = {
    "7shi-hopf1": "7shi-h",
    "7shi-qhopf1": "7shi-h",
    "7shi-hopf": "7shi-h",
    "7shi-2": "7shi-s",
    "7shi-c": "7shi-coord",
    "7shi-pcoord": "7shi-coord",
    "7shi-pbq": "7shi-bq",
    "7shi-pdq": "7shi-bq",
    "7shi-qmat": "7shi-qcm",
    "oct": "7shi-oct1",
    "oct-7rot": "7shi-7rot",
    "7shi-ort": "7shi-7rot",
    "oct-cl6": "7shi-cl6",
    "oct-nonassoc": "7shi-nonassoc",
    "7shi-vo": "7shi-mir",
    "7shi-rot": "7shi-clrt",
    "7shi-oct480": "7shi-480",
    "wikipedia-s3": "wiki-s3",
    "wikipedia-dp": "wiki-dp",
    "ref1": "hasebe",
    "ref2": "tian",
    "gpm": "wiki-gpm",
    "gca": "wiki-gca",
    # 未公開記事の本文だけに現れる自著参照の揺れ
    "pauli-qua": "7shi-bq",
    "qua-tensor": "7shi-qt",
    "qua-nonion": "7shi-nonion",
    "clif-rep": "7shi-clif1",
    "clif-sc": "7shi-clif2",
    "lie-u1": "7shi-lie1",
    "lie-su3": "7shi-lie5",
    "7shi-oct-left-mul": "7shi-cl6",
    "7shi-bd": "7shi-bloch",
    "7shi-ext": "7shi-hopfext",
    "7shi-ent": "7shi-entangle",
    "7shi-spin": "7shi-lie3",
    "7shi-su2": "7shi-lie2",
    "7shi-tensor": "7shi-tp",
    "weyl-alg": "7shi-clif4",
    "nelson-sqm": "7shi-leb6",
}

# md path -> {old slug -> canonical slug}, for slug strings whose meaning
# depends on the article (7shi, wikipedia, 7shi-hc)
OVERRIDES = {
    "hopf/01-quaternion.md": {
        "7shi": "7shi-qrot",
        "wikipedia": "wiki-hopf",
        "7shi-hc": "7shi-homog",
    },
    "misc/spherical-coords.md": {
        "7shi": "7shi-h",
        "wikipedia": "wiki-nsphere",
    },
    # 7shi-oct はoct/02の参考文献登録でははてな「八元数と7次元の外積」の
    # 正準slug。oct/01（7shi-oct1）を指すのは以下の本文のみ
    "hopf/04-extension.md": {"7shi-oct": "7shi-oct1"},
    "hopf/05-entanglement.md": {"7shi-oct": "7shi-oct1"},
}

# ambiguous slug strings: flagged for review when found outside OVERRIDES
AMBIGUOUS = {"7shi", "wikipedia"}


def load_md_entries(path: Path) -> list[tuple[str, str]]:
    """Return (md, article_id) pairs from articles.tsv; article_id is '' if none."""
    lines = path.read_text(encoding="utf-8").splitlines()
    entries: list[tuple[str, str]] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        url = cols[1] if len(cols) > 1 else ""
        md = cols[2] if len(cols) > 2 else ""
        if md:
            entries.append((md, url.rsplit("/", 1)[-1] if url else ""))
    return entries


def collect_fixes(md: str, article_id: str) -> tuple[list[tuple[str, str, int, bool]], list[str]]:
    """Return ([(old, new, body_count, in_panel)], [ambiguous slugs])."""
    rename = RENAMES | OVERRIDES.get(md, {})
    body = (ROOT / md).read_text(encoding="utf-8")
    body_counts: dict[str, int] = {}
    for slug in SLUG_RE.findall(body):
        body_counts[slug] = body_counts.get(slug, 0) + 1

    panel: set[str] = set()
    refs_path = REFS_DIR / f"{article_id}.toml"
    if article_id and refs_path.is_file():
        with refs_path.open("rb") as f:
            panel = set(tomllib.load(f))

    fixes = [
        (old, new, body_counts.get(old, 0), old in panel)
        for old, new in rename.items()
        if body_counts.get(old, 0) or old in panel
    ]
    fixes.sort(key=lambda e: e[0])
    ambiguous = sorted(
        s for s in AMBIGUOUS - set(rename) if body_counts.get(s, 0) or s in panel
    )
    return fixes, ambiguous


def render(entries: list[tuple[str, str]]) -> tuple[str, int]:
    blocks: list[str] = [
        "# 旧slug修正チェックリスト",
        "",
        "`uv run slug_fixes.py` で生成。リネーム表は `slug_fixes.py`、決定の経緯は `SLUG.md` を参照。",
        "",
        "記事ごとに、Mathlog側の**本文**（`[[slug]]` の置換）と**参考文献登録**（ラベルの付け直し）、"
        "およびローカルの md と `refs/{ID}.toml` を修正する。",
        "",
    ]
    total = 0
    for md, article_id in sorted(entries):
        fixes, ambiguous = collect_fixes(md, article_id)
        if not fixes and not ambiguous:
            continue
        total += len(fixes)
        title = f"## {md}"
        if article_id:
            title += f" — {MATHLOG_BASE}/articles/{article_id}"
        else:
            title += " — Mathlog未公開（ローカルのみ）"
        blocks.append(title)
        blocks.append("")
        for old, new, body_count, in_panel in fixes:
            places = []
            if body_count:
                places.append(f"本文{body_count}箇所")
            if in_panel:
                places.append("参考文献登録")
            blocks.append(f"- [ ] `{old}` → `{new}`（{'、'.join(places)}）")
        for slug in ambiguous:
            blocks.append(f"- [ ] `{slug}`: 対象が文脈依存のslugです。要確認（リネーム表に未登録の出現）")
        blocks.append("")
    return "\n".join(blocks) + "\n", total


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    entries = load_md_entries(ARTICLES_TSV)
    text, total = render(entries)
    args.output.write_text(text, encoding="utf-8")
    print(f"wrote {total} fixes to {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
