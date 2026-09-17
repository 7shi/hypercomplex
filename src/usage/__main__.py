"""usage.txtに記録されたトークン使用量を集計して表示するエントリポイント。詳細はREADME.mdを参照。"""

from __future__ import annotations

import argparse
from functools import reduce

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
        print(f"# {date}")
        print(totals[date])
        return 0

    for date, usage in totals.items():
        print(f"# {date}")
        print(usage)
    print()
    print("# Total")
    print(reduce(lambda a, b: a + b, totals.values()))
    return 0


if __name__ == "__main__":
    exit(main())
