from __future__ import annotations

import argparse
import json
import re
import urllib.request
from datetime import datetime, timedelta, timezone

from articles.merge import load_mathlog_tsv
from articles.paths import MATHLOG_TSV, MATHLOG_URL, ROOT, write_tsv

NEXT_DATA_RE = re.compile(r'<script id="__NEXT_DATA__"[^>]*>(.*?)</script>', re.DOTALL)

JST = timezone(timedelta(hours=9))


def fetch_html(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as res:
        return res.read().decode("utf-8")


def extract_articles(html: str) -> list[tuple[str, str, str]]:
    """Extract (date, url, title) rows from the __NEXT_DATA__ JSON embedded in
    the server-rendered page. It holds only the latest articles (20 so far)."""
    m = NEXT_DATA_RE.search(html)
    if m is None:
        raise SystemExit("__NEXT_DATA__ not found in the fetched page")
    articles = json.loads(m.group(1))["props"]["pageProps"]["articles"]
    rows: list[tuple[str, str, str]] = []
    for a in articles:
        d = datetime.fromtimestamp(a["created_at"]["seconds"], JST).strftime("%Y/%m/%d %H:%M:%S")
        title = " ".join(a["title"].split())
        rows.append((d, f"/articles/{a['id']}", title))
    return rows


def add_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("fetch", help="mathlog.url → mathlog.tsv (latest articles only)")
    parser.set_defaults(func=fetch_command)


def fetch_command(args: argparse.Namespace) -> None:
    if not MATHLOG_TSV.is_file():
        raise SystemExit(f"missing {MATHLOG_TSV}; copy it with src/bookmarklets/mathlog_articles.url")
    url = MATHLOG_URL.read_text(encoding="utf-8").strip()
    fetched = extract_articles(fetch_html(url))
    if not fetched:
        raise SystemExit(f"no articles found in {url}")

    rows = {u: (d, t) for d, u, t in load_mathlog_tsv(MATHLOG_TSV)}
    old_urls = set(rows)
    for d, u, t in fetched:
        old = rows.get(u)
        if old is None:
            print(f"added: {d} {u} {t}")
        elif old != (d, t):
            print(f"updated: {u} {old[0]} {old[1]} → {d} {t}")
        rows[u] = (d, t)

    # The page holds only the latest articles. Existing rows newer than the
    # oldest fetched one should have been fetched too; report any missing.
    fetched_urls = {u for _, u, _ in fetched}
    oldest = min(d for d, _, _ in fetched)
    for u, (d, t) in rows.items():
        if u not in fetched_urls and d >= oldest:
            print(f"not found (deleted or unlisted?): {d} {u} {t}")
    # If none of the fetched articles was known, older ones may have been
    # missed between them and the existing rows.
    if fetched_urls.isdisjoint(old_urls):
        print("warning: no overlap with mathlog.tsv; articles may be missing. "
              "Copy the full list with src/bookmarklets/mathlog_articles.url")

    lines = [f"{d}\t{u}\t{t}" for u, (d, t) in sorted(rows.items(), key=lambda r: r[1][0], reverse=True)]
    write_tsv(MATHLOG_TSV, "date\turl\ttitle", lines)
    print(f"wrote {len(lines)} rows to {MATHLOG_TSV.relative_to(ROOT)} (fetched {len(fetched)})")
