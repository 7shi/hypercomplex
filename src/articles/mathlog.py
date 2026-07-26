from __future__ import annotations

import argparse
import re
from datetime import date, datetime, timedelta
from pathlib import Path

from articles.paths import HTML_PATH, MATHLOG_TSV, ROOT, write_tsv

ARTICLE_BLOCK_OPEN = '<div class="flex-grow-1 overflow-hidden">'
DIV_TAG_RE = re.compile(r"<div\b|</div>")

HTML_COMMENT_RE = re.compile(r"<!--.*?-->")

ARTICLE_RE = re.compile(
    r'<div class="flex-shrink-0 text-muted">\s*(?P<date>.+?)\s*</div>\s*'
    r"</div>\s*"
    r'<a class="my-1 text-break text-black text-truncate-3 lh-sm fw-bold" '
    r'href="(?P<url>[^"]+)">\s*(?P<title>[^<]+?)\s*</a>',
    re.DOTALL,
)


def iter_div_blocks(html: str, open_tag: str):
    """Yield the inner content of each balanced <div ...>...</div> block
    starting with open_tag (exact string match on the opening tag)."""
    start = 0
    while True:
        idx = html.find(open_tag, start)
        if idx == -1:
            return
        content_start = idx + len(open_tag)
        depth = 1
        for m in DIV_TAG_RE.finditer(html, content_start):
            depth += 1 if m.group() == "<div" else -1
            if depth == 0:
                yield html[content_start : m.start()]
                start = m.end()
                break
        else:
            return

SECONDS_AGO_RE = re.compile(r"^(?P<seconds>\d+)秒前$")
MINUTES_AGO_RE = re.compile(r"^(?P<minutes>\d+)分前$")
HOURS_AGO_RE = re.compile(r"^(?P<hours>\d+)時間前$")
DAYS_AGO_RE = re.compile(r"^(?P<days>\d+)日前$")
MD_DATE_RE = re.compile(r"^(?P<month>\d+)月(?P<day>\d+)日$")
YMD_RE = re.compile(r"^(?P<year>\d+)年(?P<month>\d+)月(?P<day>\d+)日$")


def file_datetime(path: Path) -> datetime:
    """Local timestamp of path's last modification."""
    return datetime.fromtimestamp(path.stat().st_mtime)


def parse_date(raw: str, ref: datetime) -> str:
    """Convert Mathlog date text to yyyy/mm/dd using ref as the base moment."""
    s = " ".join(raw.split())

    if m := SECONDS_AGO_RE.fullmatch(s):
        d = ref - timedelta(seconds=int(m.group("seconds")))
        return d.strftime("%Y/%m/%d")

    if m := MINUTES_AGO_RE.fullmatch(s):
        d = ref - timedelta(minutes=int(m.group("minutes")))
        return d.strftime("%Y/%m/%d")

    if m := HOURS_AGO_RE.fullmatch(s):
        d = ref - timedelta(hours=int(m.group("hours")))
        return d.strftime("%Y/%m/%d")

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


def extract_mathlog(html: str, ref: datetime) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    seen: set[str] = set()
    for block in iter_div_blocks(html, ARTICLE_BLOCK_OPEN):
        m = ARTICLE_RE.search(block)
        if m is None:
            continue
        raw_date = "".join(HTML_COMMENT_RE.sub("", m.group("date")).split())
        url = m.group("url").strip()
        title = " ".join(m.group("title").split())
        if url in seen:
            continue
        seen.add(url)
        rows.append((parse_date(raw_date, ref=ref), url, title))
    return rows


def add_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("mathlog", help="mathlog.html → mathlog.tsv")
    parser.set_defaults(func=mathlog_command)


def mathlog_command(args: argparse.Namespace) -> None:
    ref = file_datetime(HTML_PATH)
    html = HTML_PATH.read_text(encoding="utf-8")
    rows = extract_mathlog(html, ref=ref)
    if not rows:
        raise SystemExit(f"no articles found in {HTML_PATH}")

    lines = [f"{d}\t{url}\t{title}" for d, url, title in rows]
    write_tsv(MATHLOG_TSV, "date\turl\ttitle", lines)
    print(f"wrote {len(rows)} rows to {MATHLOG_TSV.relative_to(ROOT)} (ref={ref})")
