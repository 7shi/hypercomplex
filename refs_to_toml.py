"""Convert a Mathlog references-panel HTML export (refs/*.html) to TOML.

The input is expected to be pre-formatted by html_format.py (one tag per
line). Each accordion item becomes one TOML table keyed by its label
(the same label used as a [[label]] citation marker in the articles).

Usage:
  refs_to_toml.py input.html            # write input.toml (extension swapped)
  refs_to_toml.py input.html -o out.toml
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

TYPE_MAP = {
    "Webサイト": "website",
    "書籍": "book",
    "論文": "paper",
    "Mathlogの記事": "mathlog",
}

# li key (Japanese) -> (toml key, is_list)
FIELD_MAP = {
    "著者": ("author", True),
    "サイト名": ("site", False),
    "サイトのURL": ("url", False),
    "雑誌名": ("journal", False),
    "ページ": ("pages", False),
    "出版社": ("publisher", False),
    "出版年": ("year", False),
    "出版日": ("published", False),
    "アクセス日": ("accessed", False),
}

DATE_RE = re.compile(r"(\d+)年(\d+)月(\d+)日")


def to_iso_date(s: str) -> str | None:
    m = DATE_RE.fullmatch(s.strip())
    if not m:
        return None
    y, mo, d = (int(g) for g in m.groups())
    return f"{y:04d}-{mo:02d}-{d:02d}"


def toml_string(s: str) -> str:
    escaped = s.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def toml_value(key: str, value: str) -> str:
    if key in ("accessed", "published"):
        iso = to_iso_date(value)
        if iso:
            return iso
    if key == "year":
        stripped = value.strip()
        if stripped.isdigit():
            return stripped
    return toml_string(value)


ITEM_START = re.compile(r'<div class="flex-shrink-0 me-2">\[(\d+)\]</div>')
TEXT_BREAK = re.compile(r'<div class="text-break">(.*)</div>')
LI_RE = re.compile(r"<li>([^：]+)：(.*)</li>")


def parse_entries(html: str) -> list[dict]:
    lines = html.splitlines()
    entries: list[dict] = []
    i = 0
    n = len(lines)
    while i < n:
        m = ITEM_START.search(lines[i])
        if not m:
            i += 1
            continue
        i += 1
        # citation text: either wrapped in a text-break div, or a bare
        # text line, immediately following the number marker.
        citation = ""
        tm = TEXT_BREAK.search(lines[i])
        if tm:
            citation = tm.group(1)
        else:
            citation = lines[i].strip()
        i += 1
        # advance to the <ul ...> that holds the structured fields
        while i < n and "<ul" not in lines[i]:
            i += 1
        i += 1
        fields: dict[str, str] = {}
        while i < n and "</ul>" not in lines[i]:
            lm = LI_RE.search(lines[i])
            if lm:
                fields[lm.group(1)] = lm.group(2)
            i += 1
        entries.append({"citation": citation, "fields": fields})
    return entries


def render_toml(entries: list[dict]) -> str:
    out: list[str] = []
    for entry in entries:
        fields = dict(entry["fields"])
        label = fields.pop("ラベル", None)
        if label is None:
            continue
        ref_type = TYPE_MAP.get(fields.pop("参考文献の種類", ""), "")

        out.append(f"[{label}]")
        if ref_type:
            out.append(f"type = {toml_string(ref_type)}")
        out.append(f"citation = {toml_string(entry['citation'])}")

        for ja_key, value in fields.items():
            if ja_key not in FIELD_MAP:
                continue
            toml_key, is_list = FIELD_MAP[ja_key]
            if is_list:
                items = [toml_string(v.strip()) for v in value.split(",")]
                out.append(f"{toml_key} = [{', '.join(items)}]")
            else:
                out.append(f"{toml_key} = {toml_value(toml_key, value)}")
        out.append("")
    return "\n".join(out).rstrip() + "\n"


def convert(html: str) -> str:
    return render_toml(parse_entries(html))


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("input", type=Path, help="formatted refs HTML file")
    parser.add_argument("-o", "--output", type=Path, help="write result to this file")
    args = parser.parse_args(argv)

    html = args.input.read_text(encoding="utf-8")
    result = convert(html)

    output = args.output or args.input.with_suffix(".toml")
    output.write_text(result, encoding="utf-8")


if __name__ == "__main__":
    main()
