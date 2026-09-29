from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
ARTICLES_TSV = ROOT / "articles.tsv"
MD_TSV = ROOT / "md.tsv"
SLUGS_TSV = ROOT / "slugs.tsv"
REFS_DIR = ROOT / "refs"
DEFAULT_OUTPUT = ROOT / "refs.toml"
MASTER_PATH = ROOT / "refs-master.toml"
MATHLOG_BASE = "https://mathlog.info"

SLUG_RE = re.compile(r"\[\[([^\]]+)\]\]")
BOX_LABEL_RE = re.compile(r"^&&&\w+\b.*\[([^\]]+)\]\s*$", re.MULTILINE)
MATHLOG_CITATION_RE = re.compile(r"^([^,]+), (.*), Mathlog, ")


def extract_slugs(text: str) -> list[str]:
    """Return the [[slug]] citation markers in text, deduplicated in order of
    first appearance. Labels of the article's own boxes (`&&&type title
    [label]`) are excluded: Mathlog uses the same [[label]] syntax for
    numbered references to them, but they are not bibliography slugs."""
    box_labels = set(BOX_LABEL_RE.findall(text))
    seen: set[str] = set()
    result: list[str] = []
    for slug in SLUG_RE.findall(text):
        if slug not in seen and slug not in box_labels:
            seen.add(slug)
            result.append(slug)
    return result
