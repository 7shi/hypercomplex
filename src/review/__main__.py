"""Markdown記事をLLMにレビューさせるツールのエントリポイント。詳細はREADME.mdを参照。

指定された.mdファイルはそれぞれ全文がレビュープロンプトとともにモデルへ送られ、
結果は同じstemで拡張子を.txtに変えたファイルに書き出されます
（例: hopf/01.md -> hopf/01.txt）。
"""

from __future__ import annotations

import argparse
from pathlib import Path
from llm7shi import Client

COMMON = """
記事はMathlogの記法で書かれています。通常のMarkdownに加え、定義・定理・証明等を表す`&&&type title`〜`&&&`のブロックが使えます。

```
&&&def 三角関数
三角関数は...
&&&
```

typeは次の通りです。

- axm：公理
- def：定義
- thm：定理
- cor：定理の系
- lem：補題
- conj：予想
- prop：命題
- fml：公式
- prf：証明
- ex：具体例
- exc：問題
- rem：注意
- 未指定：単純な囲み枠

titleは省略可能で、省略すると見出しは表示されません。ブロックの中ではMarkdown（引用・箇条書き・文字装飾等）が使えますが、ブロックの中にブロックは基本的に入れられません。

特定の定義・定理を番号付きで参照する場合は、`&&&type title [label]`のようにラベルを貼り、`[[label]]`で参照します。アンカーテキストを変える場合は`[アンカーテキスト](#label)`とします。

これらはMathlogで正しく解釈される記法であり、Markdownの誤りではありません。ただし、内容の性質に対してtypeの選択が適切かどうかはレビュー対象です。単なる`&&&rem`（注釈）に丸めず、`&&&def`・`&&&fml`・`&&&ex`・`&&&prf`を能動的に使い分けることが望まれます。
""".strip()

PROMPT = """
添付は数式を含む日本語の技術記事です。以下の観点でレビューしてください。

- 数式・数学的主張の正しさ
- 説明の分かりやすさ、論理の飛躍の有無
- 用語・記法の一貫性
- 日本語表現の自然さ

指摘事項を箇条書きで挙げてください。問題がなければその旨を書いてください。
""".strip()


REF_HEADER = "以下は参照用の関連記事です。レビュー対象ではなく、用語・記法・構成上の位置づけを確認するための文脈として使ってください。"


def review_file(client: Client, text: str, prompt: str, refs: list[str]):
    contents = []
    if refs:
        contents += [REF_HEADER, *refs]
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

    total_usage = None
    for path in args.files:
        print()
        print("=" * 40)
        print(f"{path}: レビュー中")
        print("=" * 40)
        print()
        result, usage = review_file(client, path.read_text(), prompt, refs)
        if usage:
            total_usage = usage if total_usage is None else total_usage + usage

        out_path = path.with_suffix(".txt")
        out_path.write_text(result + "\n")
        print(f"-> {out_path}")

    if total_usage:
        print(f"\n--- Total Usage ---\n{total_usage}")
    return 0


if __name__ == "__main__":
    exit(main())
