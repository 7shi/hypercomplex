"""Apply the slug renames to local files: md bodies and refs/*.toml keys.

Uses the rename map in slug_fixes.py (RENAMES / OVERRIDES). For each md file
in articles.tsv, [[old-slug]] markers in the body are replaced with the
canonical slug, and section keys [old-slug] in the article's refs/{ID}.toml
are renamed. The Mathlog side is untouched; fix it by hand following
slug-fixes.md (generated from the pre-fix state — do not regenerate it
until the manual fixes are done).

Usage:
  slug_apply.py            # apply in place
  slug_apply.py -n         # dry run: report what would change
"""

from __future__ import annotations

import argparse

from slug_fixes import ARTICLES_TSV, OVERRIDES, REFS_DIR, RENAMES, ROOT, load_md_entries


def apply_body(text: str, rename: dict[str, str]) -> tuple[str, list[tuple[str, str, int]]]:
    changes = []
    for old, new in rename.items():
        marker = f"[[{old}]]"
        count = text.count(marker)
        if count:
            text = text.replace(marker, f"[[{new}]]")
            changes.append((old, new, count))
    return text, changes


def apply_toml(text: str, rename: dict[str, str]) -> tuple[str, list[tuple[str, str]]]:
    changes = []
    lines = text.splitlines(keepends=True)
    for old, new in rename.items():
        header = f"[{old}]"
        for i, line in enumerate(lines):
            if line.rstrip("\n") == header:
                lines[i] = f"[{new}]" + line[len(header):]
                changes.append((old, new))
    return "".join(lines), changes


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("-n", "--dry-run", action="store_true")
    args = parser.parse_args(argv)

    total = 0
    for md, article_id in load_md_entries(ARTICLES_TSV):
        rename = RENAMES | OVERRIDES.get(md, {})

        md_path = ROOT / md
        text = md_path.read_text(encoding="utf-8")
        new_text, changes = apply_body(text, rename)
        for old, new, count in changes:
            print(f"{md}: [[{old}]] -> [[{new}]] x{count}")
            total += count
        if changes and not args.dry_run:
            md_path.write_text(new_text, encoding="utf-8")

        toml_path = REFS_DIR / f"{article_id}.toml"
        if article_id and toml_path.is_file():
            text = toml_path.read_text(encoding="utf-8")
            new_text, key_changes = apply_toml(text, rename)
            for old, new in key_changes:
                print(f"{toml_path.relative_to(ROOT)}: [{old}] -> [{new}]")
                total += 1
            if key_changes and not args.dry_run:
                toml_path.write_text(new_text, encoding="utf-8")

    mode = "would change" if args.dry_run else "changed"
    print(f"{mode} {total} places")


if __name__ == "__main__":
    main()
