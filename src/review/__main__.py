"""Markdown記事をLLMにレビューさせるツールのエントリポイント。詳細はREADME.mdを参照。

指定された.mdファイルの全文がレビュープロンプトとともにモデルへ送られ、
結果は同じstemで拡張子を.txtに変えたファイルに書き出されます
（例: hopf/01.md -> hopf/01.txt）。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml
from llm7shi import Client
from llm7shi.usage import append_usage, find_usage_file, print_today_totals

_DIR = Path(__file__).parent
ROOT = _DIR.resolve().parent.parent
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


def parse_refs(front_matter: object) -> list[str]:
    """front matterの`refs`をパスの一覧として取り出す。他のキーは無視する。"""
    if front_matter is None:
        return []
    if not isinstance(front_matter, dict):
        raise SystemExit("front matterがマッピングではありません")
    refs = front_matter.get("refs", [])
    if isinstance(refs, str):
        refs = [refs]
    if not isinstance(refs, list) or not all(isinstance(ref, str) for ref in refs):
        raise SystemExit("front matterの`refs`が文字列またはそのリストではありません")
    return refs


def split_front_matter(text: str) -> tuple[list[str], str]:
    """`---`で囲まれたYAML front matterを取り除き、(refs, 本文)を返す。

    front matterがなければ`([], text)`を返す。
    """
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        return [], text
    for i, line in enumerate(lines[1:], start=1):
        if line.strip() in ("---", "..."):
            try:
                front_matter = yaml.safe_load("\n".join(lines[1:i]))
            except yaml.YAMLError as e:
                raise SystemExit(f"front matterを解釈できません: {e}")
            return parse_refs(front_matter), "\n".join(lines[i + 1:])
    raise SystemExit("front matterが閉じられていません（終端の --- がありません）")


def wrap_file(path: Path) -> str:
    text = path.read_text().strip()
    return f'<file name="{path.name}">\n{text}\n</file>'


def review_file(client: Client, path: Path, prompt: str, refs: list[Path]) -> str:
    contents = []
    if refs:
        contents += [REF_HEADER, *(wrap_file(ref_path) for ref_path in refs)]
        contents += [TARGET_HEADER]
    contents += [wrap_file(path), COMMON, prompt]
    response = client(contents)
    return response.text.strip()


# 出力する場合はパスを入れる
USAGE_PATH = None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.strip())
    parser.add_argument("file", type=Path, help="レビュー対象の.mdファイル")
    parser.add_argument("-m", "--model", required=True,
                        help="ベンダープレフィックス付きのモデル名（例: openai:gpt-4.1-mini）")
    parser.add_argument("-p", "--prompt", type=Path,
                        help="レビュー観点を記述したプロンプトファイルのパス"
                             "（front matterの`refs`を参照記事として読む、"
                             "デフォルト: 対象の`<stem>-prompt.md`があればそれ、"
                             "なければ組み込みの汎用プロンプト）")
    parser.add_argument("-r", "--ref", type=Path, action="append", default=[],
                        help="参照文脈として使う.md/.txtファイルのパス"
                             "（レビュー対象には含めない、複数指定可）")
    parser.add_argument("--save-usage", action="store_true",
                        help="モデル名によらず使用量を記録する")
    args = parser.parse_args()

    global USAGE_PATH
    if args.model.startswith(("openai:", "gpt-")) or args.save_usage:
        USAGE_PATH = find_usage_file()

    if args.file.suffix != ".md":
        parser.error(f"{args.file}: .mdファイルではありません")
    if not args.file.exists():
        parser.error(f"{args.file}: 見つかりません")

    prompt_path = args.prompt
    if prompt_path:
        if prompt_path.suffix not in (".md", ".txt"):
            parser.error(f"{prompt_path}: .md/.txtファイルではありません")
        if not prompt_path.exists():
            parser.error(f"{prompt_path}: 見つかりません")
    else:
        # 対象記事と対になる<stem>-prompt.mdがあれば自動で使う
        default_path = args.file.with_name(f"{args.file.stem}-prompt.md")
        if default_path.is_file():
            prompt_path = default_path

    prompt = PROMPT
    refs: list[Path] = []
    if prompt_path:
        print(f"prompt: {prompt_path}")
        front_matter_refs, body = split_front_matter(prompt_path.read_text())
        # front matterのパスはリポジトリルートからの相対で書く
        refs += [ROOT / ref for ref in front_matter_refs]
        # front matterだけで本文が空なら汎用プロンプトを使う
        if body.strip():
            prompt = body.strip()

    # -rは後ろに追加し、front matterと重複するものは取り除く
    seen = {path.resolve() for path in refs}
    for path in args.ref:
        if path.resolve() not in seen:
            seen.add(path.resolve())
            refs.append(path)

    for path in refs:
        if path.suffix not in (".md", ".txt"):
            parser.error(f"{path}: .md/.txtファイルではありません")
        if not path.exists():
            parser.error(f"{path}: 見つかりません")

    client = Client(model=args.model, show_params=False, keep_history=False,
                    show_usage=True)

    for path in refs:
        print(f"ref: {path}")

    try:
        result = review_file(client, args.file, prompt, refs)
        out_path = args.file.with_suffix(".txt")
        out_path.write_text(result + "\n")
        print(f"-> {out_path}")
    finally:
        # 中断時も消費分を記録する（表示は正常終了時のみ）
        if client.usages and USAGE_PATH is not None:
            append_usage(sum(client.usages), args.model, USAGE_PATH)

    if client.usages and USAGE_PATH is not None:
        print()
        print_today_totals(USAGE_PATH, models=[args.model])
    return 0


if __name__ == "__main__":
    exit(main())
