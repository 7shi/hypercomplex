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

## pending

`articles.tsv`から、次の3種類に分類して記事一覧を出力する。各記事は上から順に判定して最初に該当した1つにだけ入る。各行の先頭に公開済みは`[公開済]`、未公開（URLなし）は`[未公開]`を付ける。

1. レビューされていない記事：対になる`.txt`ファイルが存在しない
2. レビューのプロンプトがない記事：`.txt`は存在するが、対になる`-prompt.md`ファイルが存在しない
3. レビューが反映されていない記事：`.txt`と`-prompt.md`は存在するが、`.md`の更新時刻が`.txt`より古い

## 使い方

```
uv run articles mathlog
uv run articles md
uv run articles merge
uv run articles pending
```
