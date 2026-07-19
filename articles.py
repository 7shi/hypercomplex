"""Article index tools: extract Mathlog / README links and merge them.

Subcommands:
  mathlog  — mathlog.html → mathlog.tsv  (date, url, title)
  md       — README.md files → md.tsv    (path, title)
  merge    — join mathlog.tsv and md.tsv by title → articles.tsv

Merge matching (title key):
  1. exact equality
  2. fallback: mathlog_title.endswith(md_title)  (longest md title wins)
"""

from __future__ import annotations

import argparse
import re
from datetime import date, datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parent
HTML_PATH = ROOT / "mathlog.html"
MATHLOG_TSV = ROOT / "mathlog.tsv"
MD_TSV = ROOT / "md.tsv"
ARTICLES_TSV = ROOT / "articles.tsv"

# --- mathlog ---------------------------------------------------------------

ARTICLE_RE = re.compile(
    r'<div class="flex-shrink-0 text-muted">\s*(?P<date>[^<]+?)\s*</div>\s*'
    r"</div>\s*"
    r'<a class="my-1 text-break text-black text-truncate-3 lh-sm fw-bold" '
    r'href="(?P<url>[^"]+)">\s*(?P<title>[^<]+?)\s*</a>',
    re.DOTALL,
)

DAYS_AGO_RE = re.compile(r"^(?P<days>\d+)日前$")
MD_DATE_RE = re.compile(r"^(?P<month>\d+)月(?P<day>\d+)日$")
YMD_RE = re.compile(r"^(?P<year>\d+)年(?P<month>\d+)月(?P<day>\d+)日$")


def file_date(path: Path) -> date:
    """Local calendar date of path's last modification."""
    return datetime.fromtimestamp(path.stat().st_mtime).date()


def parse_date(raw: str, ref: date) -> str:
    """Convert Mathlog date text to yyyy/mm/dd using ref as the base day."""
    s = " ".join(raw.split())

    if m := DAYS_AGO_RE.fullmatch(s):
        d = ref - timedelta(days=int(m.group("days")))
        return d.strftime("%Y/%m/%d")

    if m := YMD_RE.fullmatch(s):
        d = date(int(m.group("year")), int(m.group("month")), int(m.group("day")))
        return d.strftime("%Y/%m/%d")

    if m := MD_DATE_RE.fullmatch(s):
        d = date(ref.year, int(m.group("month")), int(m.group("day")))
        return d.strftime("%Y/%m/%d")

    raise ValueError(f"unrecognized date format: {raw!r}")


def extract_mathlog(html: str, ref: date) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    seen: set[str] = set()
    for m in ARTICLE_RE.finditer(html):
        raw_date = " ".join(m.group("date").split())
        url = m.group("url").strip()
        title = " ".join(m.group("title").split())
        if url in seen:
            continue
        seen.add(url)
        rows.append((parse_date(raw_date, ref=ref), url, title))
    return rows


def write_tsv(path: Path, header: str, lines: list[str]) -> None:
    path.write_text(header + "\n" + "\n".join(lines) + "\n", encoding="utf-8")


def cmd_mathlog() -> None:
    ref = file_date(HTML_PATH)
    html = HTML_PATH.read_text(encoding="utf-8")
    rows = extract_mathlog(html, ref=ref)
    if not rows:
        raise SystemExit(f"no articles found in {HTML_PATH}")

    lines = [f"{d}\t{url}\t{title}" for d, url, title in rows]
    write_tsv(MATHLOG_TSV, "date\turl\ttitle", lines)
    print(f"wrote {len(rows)} rows to {MATHLOG_TSV.relative_to(ROOT)} (ref={ref})")


# --- md --------------------------------------------------------------------

LINK_RE = re.compile(r"\[([^\]]+)\]\(([^)]+\.md)\)")


def extract_md() -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    seen: set[str] = set()

    for readme in sorted(ROOT.rglob("README.md")):
        text = readme.read_text(encoding="utf-8")
        base = readme.parent
        for title, href in LINK_RE.findall(text):
            href_path = href.split("#", 1)[0]
            if not href_path.endswith(".md"):
                continue
            resolved = (base / href_path).resolve()
            try:
                rel = resolved.relative_to(ROOT)
            except ValueError:
                continue
            key = rel.as_posix()
            stem = Path(key).stem
            if stem.isupper() or key in seen:
                continue
            seen.add(key)
            rows.append((key, title.strip()))

    rows.sort(key=lambda r: r[0])
    return rows


def cmd_md() -> None:
    rows = extract_md()
    if not rows:
        raise SystemExit("no .md article links found in README.md files")

    lines = [f"{path}\t{title}" for path, title in rows]
    write_tsv(MD_TSV, "md\ttitle", lines)
    print(f"wrote {len(rows)} rows to {MD_TSV.relative_to(ROOT)}")


# --- merge -----------------------------------------------------------------

def iter_tsv_data_lines(path: Path):
    """Yield non-empty TSV lines, skipping a header row if present."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return
    start = 0
    first_col = lines[0].split("\t", 1)[0]
    if first_col in {"date", "md", "url", "title", "path"}:
        start = 1
    for line in lines[start:]:
        if line.strip():
            yield line


def load_mathlog_tsv(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    for line in iter_tsv_data_lines(path):
        date_, url, title = line.split("\t", 2)
        rows.append((date_, url, title))
    return rows


def load_md_tsv(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for line in iter_tsv_data_lines(path):
        p, title = line.split("\t", 1)
        rows.append((p, title))
    return rows


def match_md_index(mathlog_title: str, md_rows: list[tuple[str, str]], used: set[int]) -> int | None:
    """Return index into md_rows matching mathlog_title, or None.

    Prefer exact title match; fall back to mathlog_title.endswith(md_title),
    choosing the longest md title among unused rows.
    """
    exact: int | None = None
    for i, (_, md_title) in enumerate(md_rows):
        if i in used:
            continue
        if mathlog_title == md_title:
            exact = i
            break
    if exact is not None:
        return exact

    best_i: int | None = None
    best_len = -1
    for i, (_, md_title) in enumerate(md_rows):
        if i in used or not md_title:
            continue
        if mathlog_title.endswith(md_title) and len(md_title) > best_len:
            best_i = i
            best_len = len(md_title)
    return best_i


def merge_rows(
    mathlog_rows: list[tuple[str, str, str]],
    md_rows: list[tuple[str, str]],
) -> list[tuple[str, str, str, str]]:
    """Join on title. Output columns: date, url, md, title.

    title prefers the Mathlog title; falls back to the md title when
    Mathlog has no row. Unmatched sides leave the missing fields empty.
    """
    used: set[int] = set()
    out: list[tuple[str, str, str, str]] = []

    for date_, url, title in mathlog_rows:
        j = match_md_index(title, md_rows, used)
        if j is None:
            out.append((date_, url, "", title))
        else:
            used.add(j)
            md_path, _ = md_rows[j]
            out.append((date_, url, md_path, title))

    for i, (md_path, title) in enumerate(md_rows):
        if i not in used:
            out.append(("", "", md_path, title))

    return out


def cmd_merge() -> None:
    if not MATHLOG_TSV.is_file():
        raise SystemExit(f"missing {MATHLOG_TSV}; run: {Path(__file__).name} mathlog")
    if not MD_TSV.is_file():
        raise SystemExit(f"missing {MD_TSV}; run: {Path(__file__).name} md")

    mathlog_rows = load_mathlog_tsv(MATHLOG_TSV)
    md_rows = load_md_tsv(MD_TSV)
    rows = merge_rows(mathlog_rows, md_rows)

    lines = [f"{d}\t{url}\t{md}\t{title}" for d, url, md, title in rows]
    write_tsv(ARTICLES_TSV, "date\turl\tmd\ttitle", lines)

    matched = sum(1 for d, url, md, title in rows if d and md)
    only_ml = sum(1 for d, url, md, title in rows if d and not md)
    only_md = sum(1 for d, url, md, title in rows if md and not d)
    print(
        f"wrote {len(rows)} rows to {ARTICLES_TSV.relative_to(ROOT)} "
        f"(matched={matched}, mathlog_only={only_ml}, md_only={only_md})"
    )


# --- CLI -------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("mathlog", help="mathlog.html → mathlog.tsv")
    sub.add_parser("md", help="README.md files → md.tsv")
    sub.add_parser("merge", help="mathlog.tsv + md.tsv → articles.tsv (by title)")

    args = parser.parse_args(argv)
    if args.command == "mathlog":
        cmd_mathlog()
    elif args.command == "md":
        cmd_md()
    elif args.command == "merge":
        cmd_merge()
    else:
        parser.error(f"unknown command: {args.command}")


if __name__ == "__main__":
    main()
