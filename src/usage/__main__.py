"""usage.jsonlに記録されたトークン使用量を集計して表示するエントリポイント。詳細はREADME.mdを参照。"""

from __future__ import annotations

import argparse

from llm7shi.usage import Usage

from usage import USAGE_PATH, parse_usage_file, today


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.strip())
    parser.add_argument("-a", "--all", action="store_true",
                        help="日付ごとの合計をすべて表示する（デフォルト: 今日の分のみ）")
    args = parser.parse_args(argv)

    totals = parse_usage_file()
    if not totals:
        print(f"{USAGE_PATH}: 記録がありません")
        return 1

    if not args.all:
        date = today()
        if date not in totals:
            print(f"{date}の記録がありません")
            return 1
        print(date)
        for model, usage in totals[date].items():
            print(f"  {model} {usage.to_dict()}")
        return 0

    model_totals: dict[str, Usage] = {}
    for date, by_model in totals.items():
        print(date)
        for model, usage in by_model.items():
            print(f"  {model} {usage.to_dict()}")
            model_totals[model] = usage if model not in model_totals else model_totals[model] + usage
    print("=" * 10)
    for model, usage in model_totals.items():
        print(f"{model} {usage.to_dict()}")
    return 0


if __name__ == "__main__":
    exit(main())
