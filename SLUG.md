# slug管理の要件定義

記事本文中の `[[slug]]` 引用ラベルを一元管理するための作業ドキュメント。次回セッションはこれを読んで作業を再開する。

## 背景

各記事はMathlogで公開されており、Mathlogの参考文献パネル（記事ごとの引用ラベル一覧）を `refs/{記事ID}.toml` としてエクスポートしてある（`html_format.py` + `refs_to_toml.py` で HTML → TOML 変換）。この参照ラベル（slug）はMathlog側で記事ごとに独立して割り当てられているため、同じ実体（同じURL、あるいは同じ他記事）を指すのに記事ごとに異なるslugが使われている「揺れ」がある。

今後の記事執筆では、slugを一元管理して次を実現したい。

1. 記事に対するslug（他記事から引用される際の正準slug）を決める
2. 記事で参照されているslugをリストアップする
3. 既に定義されているものは定義を示し、未定義のものは定義を求める

これを実現するには、`refs/*.toml` を1つに統合する必要がある。**ただし統合の前に、既存の揺れ（重複）を解決するのが先。**

## 現状でできているもの（`refs_slugs.py`）

`articles.tsv` に載っている全mdファイルを対象に、以下を生成するスクリプト。実行: `uv run refs_slugs.py`

- **`refs.toml`** — mdファイルパスごとの見出し（`["path/to/file.md"]`）に、そのファイルが持つMathlog URL（`url`）と、本文中の `[[slug]]` を出現順に解決した一覧。
  - Mathlog URLを持つファイルは `[[key.refs]]` の配列テーブルで `slug`, `type`, `url` を持つ。`type = "mathlog"` の場合は citation中のタイトルを `articles.tsv` で検索し、参照先記事自身のURLを解決している。
  - `refs/{ID}.toml` にキーが存在しないslugは `refs` ではなく、その見出し直下の `undefined = [...]` フィールドに列挙される（未定義slugの検出）。
  - Mathlog URLを持たないファイル（Mathlog未公開など）は比較対象の `refs/{ID}.toml` がないため、単純な `slugs = [...]` のフォールバック形式。
- **`refs-url.txt`** — 同一URLが複数の異なるslugで引用されている（URL基準の揺れ）ものを、URLでソートして `url = slug, slug, ...` 形式で列挙。**重複解決の材料。**
- **`refs-file.txt`** — 自分たちのMathlog記事（mdファイル）ごとに、他記事から `type = "mathlog"` で引用される際に使われているslugの一覧を、mdパスでソートして `md_path = slug, slug, ...` （未引用は `NONE`）で列挙。**同じ記事に対するslugの揺れの検出。**

これらは重複（揺れ）を「検出」するところまでで、まだ「解決（正準slugを決めて統一する）」は行っていない。

## 命名規約（確定 2026-07-19）

既存slugから帰納的に導出した規約。既存slugの修正は最低限に留め、規約から大きく外れるもののみ修正する。

1. **自著記事**（Mathlog・はてなブログ・Colab/Gist）は `7shi-<略号>`。略号は内容由来（`mir`, `coord`）または連載番号（`lie1`, `lie2`）。例外は認めない（旧 `oct-*` 系も `7shi-*` に統一）。
2. **Wikipedia** は `wiki-<略号>`（既存約20件の多数派形式）。裸の `wikipedia` および `wikipedia-<略号>` は不可 → `wiki-*` に正規化。
3. **その他外部文献**は著者名（`nakajima`, `sylvester`, `tao-h`, `furey1`）または内容・サイト略号（`dim7-8`, `half`, `proj`）。`ref1` のような無意味な連番は不可。
4. 対象を特定しない汎用slug（`7shi` 単独、`wikipedia` 単独）は不可。
5. 新規割当の略号は自己記述的にする（原則3文字以上）。単文字略号（`7shi-h`, `7shi-s` 等）は既存のみ許容。
6. 新規slugはファイル名の内容部分を略したものを基本とする（連番・ディレクトリは含めない。例: `bloch-density`→`bloch`, `cd-to-matrix`→`cdm`）。確定後はファイル名が変わってもslugは変更しない（slugは公開済み記事に埋め込まれるため、ファイル名より強い安定性が必要）。

`slugs.tsv` は `md.tsv` の全ファイルを対象とする（ローカル管理のため非公開記事も含む）。公開予定のないファイル（qua/cd/grok45.md, qua/cd/memo.md）は slug に `NONE` を残置する。

## 正準slug台帳（確定 2026-07-19）

- 記事ファイルごとの正準slugはローカルの `slugs.tsv`（md → slug）で管理する。Mathlogには記事にslugを割り振る機能がないため、これが唯一の台帳となる。
- 旧slugは**Mathlog側を手動修正する**。修正対象は各記事の本文（`[[slug]]`）と参考文献登録のラベル。どの記事のどこを直すかのリストは `slug-fixes.md` を参照（`uv run slug_fixes.py` で生成）。ローカルのmd本文と `refs/*.toml` も同じ内容に合わせて修正する。

`refs/*.toml` 全体を機械的に照合した結果、同一slugが異なる対象を指す衝突は `7shi` / `wikipedia` / `7shi-hc` の3件のみで、以下の通り解決済み。それ以外の重複定義はすべて同一対象を指しており、統合可能。

### 自著Mathlog記事（公開済み30件）

★=今回新規割当（従来は被引用なし=NONE）。旧slugは他記事の本文・参考文献登録に残っているもので、`slug-fixes.md` の修正対象。

| 記事 | 正準slug | 旧slug（要修正） |
|---|---|---|
| hopf/01-quaternion.md | 7shi-h | 7shi-hopf1, 7shi-qhopf1, 7shi, 7shi-hopf（qua/01の参考文献登録のみ・本文未使用） |
| hopf/02-spinor-tensor.md | 7shi-s | 7shi-2 |
| hopf/03-bloch-density.md | 7shi-bloch ★ | |
| hopf/c2-to-s2.md | 7shi-c2s2 ★ | |
| hopf/homogeneous.md | 7shi-homog ★ | 7shi-hc（はてな記事に譲る） |
| lie/01-u1-so2.md | 7shi-lie1 | |
| lie/02-su2-so3.md | 7shi-lie2 ★ | |
| lie/magma-to-group.md | 7shi-group | |
| misc/energy-quantize-zeta.md | 7shi-zeta ★ | |
| misc/epsilon-euler-lagrange.md | 7shi-eleq | |
| misc/exp-maclaurin-integral.md | 7shi-macl ★ | |
| misc/integration-by-parts.md | 7shi-ibp | |
| misc/spherical-coords.md | 7shi-coord | 7shi-c, 7shi-pcoord |
| misc/variable-dependence.md | 7shi-diff | |
| oct/01-octonion.md | 7shi-oct1 | oct ※`7shi-oct` は既存（はてな「八元数と7次元の外積」）のため連載番号式 |
| oct/02-7d-3rot.md | 7shi-7rot | oct-7rot, 7shi-ort |
| oct/03-oct-left-mul.md | 7shi-cl6 | oct-cl6 |
| oct/nonassociativity.md | 7shi-nonassoc | oct-nonassoc |
| qua/01-pauli-qua.md | 7shi-bq | 7shi-pbq, 7shi-pdq |
| qua/02-nonion.md | 7shi-nonion ★ | ※外部slug `nonion`（Sylvester論文）とは別 |
| qua/03-qua-tensor.md | 7shi-qt ★ | |
| qua/cd/cd-to-matrix.md | 7shi-cdm ★ | |
| qua/cd/heuristic-cmatrix.md | 7shi-qcm | 7shi-qmat |
| qua/cd/matrix-to-pauli.md | 7shi-qp | |
| qua/cd/real-to-cmatrix.md | 7shi-rcm ★ | |
| qua/real-pauli-trace.md | 7shi-rpt ★ | |
| qua/spherical-trig.md | 7shi-strig ★ | |
| qua/tensor-from-complex.md | 7shi-tp | |
| vec-oct/01-reflection.md | 7shi-mir | 7shi-vo |
| vec-oct/geometric-product-exp.md | 7shi-vge | |

### 未公開記事（24件）

`slugs.tsv` に割当済み。連載は番号式（`7shi-clif1`〜`clif5`, `7shi-leb1`〜`leb6`, `7shi-lie3`〜`lie5`）、単発は内容略号。

| 記事 | 正準slug | 備考 |
|---|---|---|
| hopf/04-extension.md | 7shi-hopfext | |
| hopf/05-entanglement.md | 7shi-entangle | |
| hopf/06-clifford-gates.md | 7shi-qgate | 量子ゲート |
| qua/04-4d-bsqua.md | 7shi-4drot | 4次元回転 |
| qua/05-dual-qua.md | 7shi-dualq | 二重四元数 |
| qua/cd/bcmp-qua.md | 7shi-bcq | 双複素数→四元数 |
| qua/spinor-ideal.md | 7shi-ideal | |
| vec-oct/02-rotation.md | 7shi-vrot | ※`7shi-rot` は廃止済みslug（→7shi-clrt）のため再利用しない |
| qua/cd/grok45.md | NONE | 公開予定なし |
| qua/cd/memo.md | NONE | 公開予定なし |

### はてなブログ記事

| 記事（URL） | 正準slug | 廃止slug（alias） |
|---|---|---|
| クリフォード代数による回転（2020/03/07） | 7shi-clrt | 7shi-rot |
| 四元数による回転を実行列表現で考える（2020/09/10） | 7shi-qrot | 7shi |
| 八元数の定義が480通りあることの導出（2024/03/11） | 7shi-480 | 7shi-oct480 |
| NotebookLMで多元数の記事を読む（2024/06/13） | 7shi-hc | ※2記事で使用の多数派につき維持 |

### 規約逸脱の外部slug修正

| 旧slug | 新slug | 対象 |
|---|---|---|
| wikipedia（hopf/01のrefs） | wiki-hopf | en.wikipedia: Hopf fibration |
| wikipedia（misc/spherical-coordsのrefs） | wiki-nsphere | ja.wikipedia: 超球面 |
| wikipedia-s3（misc/spherical-coordsのrefs） | wiki-s3 | ja.wikipedia: 三次元球面 |
| wikipedia-dp（hopf/02のrefs） | wiki-dp | ja.wikipedia: 直積 (ベクトル) |
| ref1（oct/03のrefs） | hasebe | Hasebe, Hopf Maps, Lowest Landau Level, and Fuzzy Spheres, arXiv, 2010 |
| ref2（oct/03のrefs） | tian | Tian, Matrix Representations of Octonions and Their Applications, arXiv, 2000 |
| gpm（qua/02のrefs） | wiki-gpm | en.wikipedia: Generalizations of Pauli matrices |
| gca（qua/02のrefs） | wiki-gca | en.wikipedia: Generalized Clifford algebra |

### 未公開記事本文の揺れ（確定 2026-07-19）

未公開記事のローカル本文だけに現れる自著参照の揺れ。Mathlog修正は不要で、ローカル置換のみ（`slug-fixes.md` に含まれる）。

| 旧slug | 正準slug | 備考 |
|---|---|---|
| pauli-qua | 7shi-bq | |
| qua-tensor | 7shi-qt | |
| qua-nonion | 7shi-nonion | |
| clif-rep | 7shi-clif1 | |
| clif-sc | 7shi-clif2 | |
| lie-u1 | 7shi-lie1 | |
| lie-su3 | 7shi-lie5 | |
| 7shi-oct-left-mul | 7shi-cl6 | |
| 7shi-bd | 7shi-bloch | |
| 7shi-ext | 7shi-hopfext | |
| 7shi-ent | 7shi-entangle | |
| 7shi-spin | 7shi-lie3 | |
| 7shi-su2 | 7shi-lie2 | |
| 7shi-tensor | 7shi-tp | クロネッカー積の入門的説明の参照のため tensor-from-complex |
| weyl-alg | 7shi-clif4 | |
| nelson-sqm | 7shi-leb6 | |
| 7shi-oct | 7shi-oct1 | **hopf/04, hopf/05の本文のみ**。oct/02の参考文献登録の `7shi-oct` ははてな「八元数と7次元の外積」の正準slugなのでそのまま |

未公開記事の本文には、このほか未定義の外部参照（`wiki-ccr`, `bernevig`, `nielsen-chuang` 等、いずれも規約適合の命名）が残っているが、これらは揺れではなく定義待ち。マスターレジストリ整備時に定義を登録する。

## 修正手順（進行中）

1. ~~**ローカル一括修正**~~ — 完了（2026-07-19）。`uv run slug_apply.py` で全md本文162箇所と `refs/*.toml` のキーを置換。`refs_slugs.py` 再実行で揺れゼロ（refs-url.txt空、refs-file.txtがslugs.tsvと全一致）を確認済み。
2. **Mathlog側の手動修正（18記事）** — `slug-fixes.md` をチェックリストに、記事ごとに修正済みローカルmdを本文に貼り直し、参考文献登録のラベルを付け直す。
   - **注意**: `slug-fixes.md` は修正前の状態から生成したチェックリスト。ローカルは修正済みのため、手動作業が終わるまで `slug_fixes.py` を再実行しないこと（再実行すると空のリストで上書きされる）。
3. **整合確認** — Mathlog側を直した記事から参考文献パネルを再エクスポート（`html_format.py` + `refs_to_toml.py`）し、ローカルの `refs/*.toml` と一致することを確認。
4. **マスターレジストリ統合** — 下記参照。

## 未着手（次にやること）

### 1. 統合マスターレジストリの作成

`refs/*.toml`（記事ごとにバラバラ）を1つの `slug → 定義（type, citation, url等）` のマスターテーブルに統合する。namespaceは1つ（外部サイト／論文／書籍と、自分たちの他記事へのmathlog参照が同じslug空間を共有する）。

- 出力先ファイル名は未定（例: `refs/master.toml`）。
- 正準slugをキーとする。旧slugの手動修正（`slug-fixes.md`）が完了すれば旧slugは残らないため、aliasは持たない。
- 自著Mathlog記事のエントリには対象mdパス（`md = "..."`）も持たせる（`slugs.tsv` と整合させる）。
- 統合時、同一slugに異なる対象が混ざったらエラーにする（衝突検出の再確認）。

### 2. 新記事執筆時のワークフロー支援

新しいチェックスクリプト（例: `check_slugs.py`）を作り、対象mdファイルを渡すと以下を行うようにする。

1. 記事自身の正準slugをマスターレジストリに登録する（Mathlog URL確定前は仮登録も許容）。
2. 本文中の `[[slug]]` を抽出する（`refs_slugs.py` の `SLUG_RE` を流用可能）。
3. マスターレジストリと突き合わせ、定義済みslugは定義を表示、未定義slugは一覧表示して定義（type/citation/url）の追加を促す。

## 関連ファイル

- `slugs.tsv` — **正準slug台帳**（md → slug）。手動管理
- `slug_fixes.py` — 旧slugの修正箇所リストを生成するスクリプト（旧→新のリネーム表を内蔵）
- `slug-fixes.md` — 生成物。Mathlog側の手動修正チェックリスト（修正前状態から生成、作業完了まで再生成しないこと）
- `slug_apply.py` — リネーム表をローカル（md本文と `refs/*.toml`）に一括適用するスクリプト（適用済み）
- `refs_slugs.py` — 検出スクリプト本体
- `refs.toml`, `refs-url.txt`, `refs-file.txt` — 生成物（`uv run refs_slugs.py` で再生成可能）
- `refs/*.toml` — Mathlog記事ごとの参考文献エクスポート（`refs_to_toml.py` で `refs/*.html` から生成）
- `articles.tsv` — 記事一覧（date, url, md, title）。`articles.py` で生成・更新
