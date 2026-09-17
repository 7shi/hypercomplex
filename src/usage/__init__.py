"""usage.txtにトークン使用量を記録・集計するツール。詳細はREADME.mdを参照。"""

from __future__ import annotations

import ast
import re
from datetime import datetime
from pathlib import Path

from llm7shi.usage import Usage

USAGE_PATH = Path(__file__).resolve().parents[2] / "usage.txt"

DATE_RE = re.compile(r"^#\s*(\d{4}/\d{2}/\d{2})\s*$")
USAGE_RE = re.compile(r"^Usage\((\{.*\})\)\s*$")


def today() -> str:
    """今日の日付を見出しで使う`YYYY/MM/DD`形式で返す。"""
    return datetime.now().strftime("%Y/%m/%d")


def parse_usage_file(path: Path = USAGE_PATH) -> dict[str, Usage]:
    """usage.txtをパースし、日付をkey、その日の合計Usageをvalueとするdictを返す。

    日付はファイル内での出現順を保ちます。ファイルが存在しない場合は空のdictを返します。
    """
    if not path.exists():
        return {}

    totals: dict[str, Usage] = {}
    date = None
    for line in path.read_text().splitlines():
        if m := DATE_RE.match(line):
            date = m.group(1)
            continue
        if not (m := USAGE_RE.match(line)):
            continue
        if date is None:
            raise ValueError(f"{path}: 日付見出しの前にUsage行があります: {line}")
        # 手書きのJSON風（"）とreprの出力（'）が混在しうるのでliteral_evalで読む
        usage = Usage(raw=ast.literal_eval(m.group(1)))
        totals[date] = usage if date not in totals else totals[date] + usage
    return totals


def append_usage(usage: Usage, path: Path = USAGE_PATH, date: str | None = None) -> None:
    """Usageを指定日（省略時は今日）のセクションに追記する。

    日付の見出しがなければ末尾に作成し、あればそのセクションの末尾に追記します。
    """
    date = date or today()
    line = repr(usage)

    lines = path.read_text().splitlines() if path.exists() else []

    start = next((i for i, l in enumerate(lines) if (m := DATE_RE.match(l)) and m.group(1) == date), None)
    if start is None:
        if lines and lines[-1].strip():
            lines.append("")
        lines += [f"# {date}", line]
    else:
        end = next((i for i in range(start + 1, len(lines)) if DATE_RE.match(lines[i])), len(lines))
        while end > start + 1 and not lines[end - 1].strip():
            end -= 1
        lines.insert(end, line)

    path.write_text("\n".join(lines) + "\n")
