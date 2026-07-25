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
MATHLOG_CITATION_RE = re.compile(r"^([^,]+), (.*), Mathlog, ")
