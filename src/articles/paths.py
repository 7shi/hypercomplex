from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
MATHLOG_TSV = ROOT / "mathlog.tsv"
MATHLOG_URL = ROOT / "mathlog.url"
MD_TSV = ROOT / "md.tsv"
ARTICLES_TSV = ROOT / "articles.tsv"


def write_tsv(path: Path, header: str, lines: list[str]) -> None:
    path.write_text(header + "\n" + "\n".join(lines) + "\n", encoding="utf-8")
