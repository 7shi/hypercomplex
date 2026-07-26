# articles

Mathlog記事とREADMEリンクを突き合わせるツール。

## mathlog

`mathlog.html`（Mathlogの記事一覧ページのHTML）から記事一覧を抽出し、`mathlog.tsv`（日付・URL・タイトル）に書き出す。

## md

各`README.md`が持つ記事へのリンクを集め、`md.tsv`（パス・タイトル）に書き出す。

## merge

`mathlog.tsv`と`md.tsv`をタイトルで突き合わせ、`articles.tsv`（日付・URL・パス・タイトル）に書き出す。

マッチングはタイトルをキーに行う。

1. 完全一致
2. フォールバック: `mathlog_title.endswith(md_title)`（該当する`md_title`のうち最長のものを採用）

## 使い方

```
uv run articles mathlog
uv run articles md
uv run articles merge
```
