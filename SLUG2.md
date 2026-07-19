# slug管理 引継ぎ事項（セッション再開用）

作成: 2026-07-19。規約・台帳・決定の経緯は `SLUG.md` にすべて記録済み。このファイルは「今どこまで進んでいて、次に何をするか」だけをまとめる。

## 現在地

- **命名規約と正準slugは全記事分確定済み**（SLUG.md参照）。`slugs.tsv`（md → slug）が唯一の台帳。md.tsvの全54ファイルが対象で、公開予定のない qua/cd/grok45.md と qua/cd/memo.md は `NONE`。
- **ローカル一括修正は完了**。`slug_apply.py` で全md本文と `refs/*.toml` のキーを162箇所置換し、検証済み:
  - `refs-url.txt` が空（同一URLの複数slug引用が消滅）
  - `refs-file.txt` が `slugs.tsv` と全一致
  - `slug_apply.py -n` の再実行で0件（冪等・取り残しなし）
- **未完了**: Mathlog側の手動修正（下記）とその後続。

## 次にやること（SLUG.md「修正手順」の続き）

1. **手順2（ユーザーの手動作業・未着手）**: `slug-fixes.md` をチェックリストに、Mathlog側18記事を修正する。記事ごとに修正済みローカルmdを本文に貼り直し、参考文献登録のラベルを付け直す。セッション側の仕事は進捗確認とサポート。
2. **手順3**: 修正済み記事の参考文献パネルを再エクスポート（Mathlogで `mathlog_ref.url` ブックマークレット → `refs/*.html` → `html_format.py` + `refs_to_toml.py`）し、ローカルの `refs/*.toml` と一致することを確認。
3. **手順4**: マスターレジストリ統合。`refs/*.toml` を `refs/master.toml` に統合（正準slugがキー、自著Mathlog記事は `md = "..."` を持たせ `slugs.tsv` と整合、同一slugに異なる対象が混ざったらエラー）。未公開記事本文に残る**未定義の外部参照23箇所**（`wiki-ccr`, `wiki-schrodinger`, `bernevig`, `nielsen-chuang` 等、命名は規約適合）はここで定義を登録する。
4. **check_slugs.py の作成**（SLUG.md「未着手」参照）: mdファイルを渡すと `[[slug]]` を抽出し、マスターレジストリと突き合わせて定義済みは定義を表示、未定義は一覧して追加を促す。

## 重要な注意

- **`slug_fixes.py` を再実行しないこと**（手順2完了まで）。`slug-fixes.md` は修正前の状態から生成したチェックリストで、ローカルは修正済みのため、再実行すると空リストで上書きされる。
- `slug_apply.py` は適用済み。再実行しても0件で無害。
- 紛らわしいslug: `7shi-oct` ははてな「八元数と7次元の外積」、`7shi-oct1` は oct/01-octonion.md。`7shi-hc` ははてな「NotebookLMで多元数の記事を読む」（hopf/homogeneous.md は `7shi-homog`）。
- 今回の変更は**未コミット**。slug関連の変更一式と、着手前からあった qua/04-4d-bsqua.md の変更が混在しているので、コミット時は分けること。

## 関連ファイル

- `SLUG.md` — 規約・正準slug台帳・修正手順の本体
- `slugs.tsv` — 正準slug台帳（手動管理）
- `slug-fixes.md` — Mathlog手動修正チェックリスト（凍結中）
- `slug_fixes.py` — チェックリスト生成（リネーム表を内蔵、`slug_apply.py` からも参照）
- `slug_apply.py` — ローカル一括適用スクリプト（適用済み）
- `refs_slugs.py` → `refs.toml`, `refs-url.txt`, `refs-file.txt` — 検出・検証（再生成可）
