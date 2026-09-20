"""Markdown記事をLLMにレビューさせるツールのエントリポイント。詳細はREADME.mdを参照。

指定された.mdファイルの全文がレビュープロンプトとともにモデルへ送られ、
結果は同じstemで拡張子を.txtに変えたファイルに書き出されます
（例: hopf/01.md -> hopf/01.txt）。
"""

from __future__ import annotations

import argparse
from pathlib import Path
from llm7shi import Client
from llm7shi.usage import append_usage, find_usage_file, format_usage_line, parse_usage_file, today

_DIR = Path(__file__).parent
COMMON = (_DIR / "COMMON.txt").read_text().strip()
PROMPT = (_DIR / "PROMPT.txt").read_text().strip()


REF_HEADER = (
    "以下は参照用の関連記事です。"
    "用語・記法・構成上の位置づけを確認するための文脈として使ってください。"
    "ただし、整合性を取るために参照記事の側を修正すべきだと判断した場合は、"
    "どのファイルのどの箇所をどう直すべきかを明示してレビューに含めてください。"
)
TARGET_HEADER = (
    "ここからが今回のレビュー対象の記事です。"
    "上記の参照用記事ではなく、この記事についてレビューしてください。"
)


def wrap_file(path: Path) -> str:
    text = path.read_text().strip()
    return f'<file name="{path.name}">\n{text}\n</file>'


def review_file(client: Client, path: Path, prompt: str, refs: list[Path]):
    contents = []
    if refs:
        contents += [REF_HEADER, *(wrap_file(ref_path) for ref_path in refs)]
        contents += [TARGET_HEADER]
    contents += [wrap_file(path), COMMON, prompt]
    response = client(contents)
    if response.usage:
        print(f"\n{response.usage}")
    return response.text.strip(), response.usage


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.strip())
    parser.add_argument("file", type=Path, help="レビュー対象の.mdファイル")
    parser.add_argument("-m", "--model", required=True,
                        help="ベンダープレフィックス付きのモデル名（例: openai:gpt-4.1-mini）")
    parser.add_argument("-p", "--prompt", type=Path,
                        help="レビュー観点を記述したプロンプトファイルのパス"
                             "（デフォルト: 組み込みの汎用プロンプト）")
    parser.add_argument("-r", "--ref", type=Path, action="append", default=[],
                        help="参照文脈として使う.md/.txtファイルのパス"
                             "（レビュー対象には含めない、複数指定可）")
    args = parser.parse_args()

    if args.file.suffix != ".md":
        parser.error(f"{args.file}: .mdファイルではありません")
    if not args.file.exists():
        parser.error(f"{args.file}: 見つかりません")

    for path in args.ref:
        if path.suffix not in (".md", ".txt"):
            parser.error(f"{path}: .md/.txtファイルではありません")
        if not path.exists():
            parser.error(f"{path}: 見つかりません")

    prompt = PROMPT
    if args.prompt:
        if not args.prompt.exists():
            parser.error(f"{args.prompt}: 見つかりません")
        prompt = args.prompt.read_text().strip()

    client = Client(model=args.model, show_params=False, keep_history=False)

    result, usage = review_file(client, args.file, prompt, args.ref)

    out_path = args.file.with_suffix(".txt")
    out_path.write_text(result + "\n")
    print(f"-> {out_path}")

    if usage:
        usage_path = find_usage_file()
        append_usage(usage, args.model, usage_path)
        print(f"-> {usage_path}")

        totals = parse_usage_file(usage_path)
        date = today()
        print(f"\n# {date}")
        for model, model_usage in totals[date].items():
            print(format_usage_line(model, model_usage))
    return 0


if __name__ == "__main__":
    exit(main())
