# articles

Mathlog記事とREADMEリンクを突き合わせるツール。

## fetch

`mathlog.url`に書かれたMathlogの記事一覧ページを取得し、ページに埋め込まれたJSON（`__NEXT_DATA__`）の記事を`mathlog.tsv`（日時・URL・タイトル）に反映する。日時は各記事の`created_at`を日本時間に直したもの（`yyyy/mm/dd hh:mm:ss`）。

取得できるのは最新の記事（現状20件）だけなので、既存の`mathlog.tsv`に無い記事を追加し、日時・タイトルが変わった記事を更新する差分方式とする。取得した記事が既存の行と1件も重ならないときは、取りこぼしの可能性を警告する。

`mathlog.tsv`を最初に作るときや全件を取り直すときは、記事一覧ページで全記事を表示させてからブックマークレット`src/bookmarklets/mathlog_articles.url`を実行し、`winclip -o mathlog.tsv`で保存する。出力の形式は`fetch`と同じ。

## md

各`README.md`が持つ記事へのリンクを集め、`md.tsv`（パス・タイトル）に書き出す。

## merge

`mathlog.tsv`と`md.tsv`をタイトルで突き合わせ、`articles.tsv`（日時・URL・パス・タイトル）に書き出す。

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
uv run articles fetch
uv run articles md
uv run articles merge
uv run articles pending
```
