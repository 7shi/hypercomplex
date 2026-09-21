# レビューと修正反映

記事をLLMにレビューさせる手順と、公開済み記事を修正してMathlogへ反映する手順をまとめます。ツールの仕様は [src/review/README.md](src/review/README.md) を参照してください。

## レビューの実行

### 1. 対象を確認する

```bash
uv run articles pending
```

未公開（Mathlog未掲載）かつ未レビュー（対になる`.txt`がない）の記事が一覧されます。レビュー待ちの記事はここに出ます。

### 2. プロンプトファイルを用意する

レビュー対象と同じディレクトリに `<stem>-prompt.md` を置くと、`review` が自動で読み込みます（例: `clif/04-weyl-algebra.md` → `clif/04-weyl-algebra-prompt.md`）。YAML front matterの`refs`に書いたファイルが参照記事として渡されます。

```markdown
---
refs:
  - clif/03-gpm-gca.md
  - lie/05-bch-adjoint.md
---

添付は「クリフォード代数」シリーズ第4回の記事です。以下の観点でレビューしてください。

（以下、レビュー観点）
```

- **参照記事（`refs`）** — 用語・記法・構成上の位置づけの整合性を照合させるための文脈です。パスはリポジトリルートからの相対で書きます。全文がコンテキストに入るため、**2件程度**に絞ります。増やすとトークン消費が膨らみ、レビューの焦点もぼやけます。
- **参照記事の選び方** — 対象記事が定義なしに使う記号・定理の出所（多くは直前記事）を最優先し、次に同じ概念の定義本体を持つ別シリーズの記事を選びます。`uv run articles pending` に載っている記事は内容が確定していないため、参照記事には選びません。
- **本文** — デフォルトプロンプト [src/review/PROMPT.txt](src/review/PROMPT.txt) を置き換えるので、その4観点（数式・分かりやすさ・用語記法の一貫性・日本語表現）を土台に、記事固有の重点項目を足す形で書きます。重点項目がなければ本文は空でよく、その場合はデフォルトプロンプトが使われます。

ボックス記法（`&&&`）とラベル参照（`[[label]]`）の説明は [src/review/COMMON.txt](src/review/COMMON.txt) が自動で渡すため、プロンプトに書く必要はありません。

Claude Code で作業する場合は、参照記事の選定とこのファイルの生成をスキル `review-refs` が行います。

### 3. 実行する

```bash
uv run review clif/04-weyl-algebra.md -m MODEL
```

読み込んだプロンプトと参照記事が`prompt:`・`ref:`として表示されます。結果は対象と同じ場所に`<stem>.txt`（例: `clif/04-weyl-algebra.txt`）として保存され、トークン使用量が`~/.local/state/llm7shi/usage.jsonl`に追記されます。

複数のモデルでレビューする場合は、`.txt`が上書きされるため、必要に応じて前の結果を退避します。

### 4. 結果を反映する

指摘は取捨選択して記事に反映します。誤指摘も混じるため、数式の指摘は検算して確認します（検証スクリプトは各ディレクトリの`check/`に置き、`uv run`で実行します）。

`<stem>.txt` は `uv run articles pending` のレビュー済み判定に使われるので、削除せず残します。

Claude Code で作業する場合は、`<stem>.txt` の検討と反映をスキル `review-apply` が行います（修正前に対応方針を提示し、判断を求めます）。

## 公開済み記事の修正手順

公開済みの記事を修正して Mathlog と同期する手順です。Mathlog には記事本文を外部から更新する API がないため、Mathlog 側への反映はブラウザ上で行います。

### 1. 本文のみの修正（誤字脱字・数式・説明など）

1. **ローカルの編集**: 対象の Markdown ファイル（例: `clif/01-representation.md`）を修正します。
2. **Mathlog URL の確認**: `uv run reftools show` を実行して対応する URL を確認します。
   ```bash
   uv run reftools show clif/01-representation.md
   ```
   出力の先頭に `url: https://mathlog.info/articles/...` が表示されます（[articles.tsv](articles.tsv) を直接参照しても構いません）。
3. **Mathlog 側の更新**: ブラウザで記事の編集画面を開き、修正内容を手動で反映・保存します。
   - ※ Mathlog のエディタへのコピペは、Visual Studio Code から行うと改行が崩れません。

### 2. 参考文献（slug）の変更を伴う修正

参考文献ラベル（slug）の改名・重複統一や書誌情報の変更を行う場合は、`refs/*.html` / `refs/*.toml` の同期や整合性検証が必要になります。

- まとめて修正する場合は、作業指示ファイル `mathlog_fix.md` を作成して [src/mathlog_fix.sh](src/mathlog_fix.sh) を利用します。
- 詳細な運用ルールや各ツールの仕様は [SLUG.md](SLUG.md) を参照してください。
