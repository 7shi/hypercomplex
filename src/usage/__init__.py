"""usage.jsonlにトークン使用量を記録・集計するツール。詳細はREADME.mdを参照。"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from llm7shi.usage import Usage

USAGE_PATH = Path(__file__).resolve().parents[2] / "usage.jsonl"


def today() -> str:
    """今日の日付をUTCで`YYYY/MM/DD`形式で返す。"""
    return datetime.now(timezone.utc).strftime("%Y/%m/%d")


def parse_usage_file(path: Path = USAGE_PATH) -> dict[str, dict[str, Usage]]:
    """usage.jsonlをパースし、日付をkey、モデル名をkeyとする合計Usageのdictをvalueとするdictを返す。

    日付は各レコードの`timestamp`をUTCに変換して求める。日付・モデル名はファイル内での出現順を保つ。
    ファイルが存在しない場合は空のdictを返す。
    """
    if not path.exists():
        return {}

    totals: dict[str, dict[str, Usage]] = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        timestamp = datetime.fromisoformat(record["timestamp"]).astimezone(timezone.utc)
        date = timestamp.strftime("%Y/%m/%d")
        model = record["model"]
        usage = Usage(raw={k: v for k, v in record.items() if k not in ("timestamp", "model")})
        by_model = totals.setdefault(date, {})
        by_model[model] = usage if model not in by_model else by_model[model] + usage
    return totals


def append_usage(usage: Usage, model: str, path: Path = USAGE_PATH, timestamp: datetime | None = None) -> None:
    """Usageをモデル名・タイムゾーン付きの生成日時とともにusage.jsonlに1行追記する。

    timestampを省略した場合は現在時刻（ローカルのタイムゾーン付き）を使う。
    """
    timestamp = timestamp or datetime.now().astimezone()
    record = {"timestamp": timestamp.isoformat(), "model": model, **usage.to_dict()}
    with path.open("a") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")
