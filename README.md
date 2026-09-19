# hypercomplex

超複素数（八元数など）に関する記事と検証コードをまとめたリポジトリです。

## 記事

- [四元数](qua/README.md) — 四元数に関する記事。
- [八元数](oct/README.md) — 八元数の左作用に関する記事と検証コード。
- [クリフォード代数](clif/README.md) — クリフォード代数に関する記事。
- [クリフォード解析](clif-analysis/README.md) — 複素解析・四元数解析・クリフォード解析に関する記事。
- [リー群・リー代数の初歩](lie/README.md) — リー群とリー代数の記事シリーズ。
- [ホップファイブレーション](hopf/README.md) — ホップファイブレーションに関する記事。
- [ベクトルから八元数まで](vec-oct/README.md) — ベクトル・複素数・四元数・八元数・クリフォード代数に関する記事。
- [ルベーグ積分](lebesgue/README.md) — 測度とルベーグ積分に関する記事。
- [その他](misc/README.md) — 主要テーマに含まれない単発記事。

## 整理

- [全体構想](PLAN.md) — 幾何代数ベースで関係分野を再構築するという狙いと、今後の展開。
- [テーマ別記事リスト](THEMES.md) — ディレクトリ構成とは独立にテーマで記事を整理し、重複・統合候補を注記したリスト。将来の本への統合に向けた準備段階です。
- [slug管理](SLUG.md) — 参考文献ラベル（slug）の命名規約と運用ルール。

## ツール

- [src/](src/README.md) — 記事執筆・検証・参考文献管理に使うスクリプト類。
- [src/articles/](src/articles/README.md) — Mathlog記事とREADMEリンクを突き合わせるツール。
  - mathlog.tsv — mathlog.html から抽出した記事一覧（日付・URL・タイトル）。
  - md.tsv — README.md から抽出した記事一覧（パス・タイトル）。
  - articles.tsv — 両者をタイトルで突き合わせた結果。
- usage.jsonl — LLMのトークン使用量（[llm7shi](https://github.com/7shi/llm7shi)の`usage`コマンドで記録・集計）の記録。

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
