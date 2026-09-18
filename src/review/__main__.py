"""Markdown記事をLLMにレビューさせるツールのエントリポイント。詳細はREADME.mdを参照。

指定された.mdファイルはそれぞれ全文がレビュープロンプトとともにモデルへ送られ、
結果は同じstemで拡張子を.txtに変えたファイルに書き出されます
（例: hopf/01.md -> hopf/01.txt）。
"""

from __future__ import annotations

import argparse
from pathlib import Path
from llm7shi import Client
from llm7shi.usage import append_usage, find_usage_file

_DIR = Path(__file__).parent
COMMON = (_DIR / "COMMON.txt").read_text().strip()
PROMPT = (_DIR / "PROMPT.txt").read_text().strip()


REF_HEADER = "以下は参照用の関連記事です。レビュー対象ではなく、用語・記法・構成上の位置づけを確認するための文脈として使ってください。"
TARGET_HEADER = "ここからが今回のレビュー対象の記事です。上記の参照用記事ではなく、この記事についてレビューしてください。"


def review_file(client: Client, text: str, prompt: str, refs: list[str]):
    contents = []
    if refs:
        contents += [REF_HEADER, *refs]
        contents += [TARGET_HEADER]
    contents += [text, COMMON, prompt]
    response = client(contents)
    if response.usage:
        print(f"\n{response.usage}")
    return response.text.strip(), response.usage


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.strip())
    parser.add_argument("files", nargs="+", type=Path, help="レビュー対象の.mdファイル（複数指定可）")
    parser.add_argument("-m", "--model", required=True,
                        help="ベンダープレフィックス付きのモデル名（例: openai:gpt-4.1-mini）")
    parser.add_argument("-p", "--prompt", type=Path,
                        help="レビュー観点を記述したプロンプトファイルのパス"
                             "（デフォルト: 組み込みの汎用プロンプト）")
    parser.add_argument("-r", "--ref", type=Path, action="append", default=[],
                        help="参照文脈として使う.mdファイルのパス"
                             "（レビュー対象には含めない、複数指定可）")
    args = parser.parse_args()

    for path in args.files:
        if path.suffix != ".md":
            parser.error(f"{path}: .mdファイルではありません")
        if not path.exists():
            parser.error(f"{path}: 見つかりません")

    for path in args.ref:
        if path.suffix != ".md":
            parser.error(f"{path}: .mdファイルではありません")
        if not path.exists():
            parser.error(f"{path}: 見つかりません")

    prompt = PROMPT
    if args.prompt:
        if not args.prompt.exists():
            parser.error(f"{args.prompt}: 見つかりません")
        prompt = args.prompt.read_text().strip()

    refs = [path.read_text() for path in args.ref]

    client = Client(model=args.model, show_params=False, keep_history=False)

    usages = []
    for path in args.files:
        print()
        print("=" * 40)
        print(f"{path}: レビュー中")
        print("=" * 40)
        print()
        result, usage = review_file(client, path.read_text(), prompt, refs)
        if usage:
            usages.append(usage)

        out_path = path.with_suffix(".txt")
        out_path.write_text(result + "\n")
        print(f"-> {out_path}")

    if usages:
        total_usage = sum(usages)
        print(f"\n--- Total Usage ---\n{total_usage}")
        usage_path = find_usage_file()
        append_usage(total_usage, args.model, usage_path)
        print(f"-> {usage_path}")
    return 0


if __name__ == "__main__":
    exit(main())
