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

## 未着手（次にやること）

### 1. 重複解決（最優先）

`refs-url.txt` と `refs-file.txt` の各エントリについて、正準slugを1つ決める。

- 既存のSLUGの修正は最低限にする方針。よほど雑なもの以外は現状のslugを認め、そこからルールを決める。
- 決定結果をどこかに記録する必要がある（例: `slug → 正準slug` のマッピングファイル）。過去記事のMathlog上の表示自体は変更できない（公開済みのため）ので、これは「今後の管理台帳」として運用する前提。

#### 分析結果（未確定・要ユーザー確認）

`refs-file.txt`が拾う「対象記事1つに複数slug」以外に、**同じslug文字列が別々の対象を指す衝突**も見つかった。各`refs/{ID}.toml`はMathlog記事ごとにスコープが閉じているため今は問題にならないが、1つに統合すると上書き事故になる。

衝突3件:

| slug | 意味A | 意味B |
|---|---|---|
| `7shi` | hopf/01-quaternion.md（misc/spherical-coordsが引用時に使用） | 7shiのはてなブログ記事（複数記事の自身のrefsで使用） |
| `wikipedia` | Hopf fibration英語版（hopf/01自身のrefs） | 超球面 日本語版（misc/spherical-coords自身のrefs） |
| `7shi-hc` | hopf/homogeneous.md（hopf/01が引用時に使用） | 7shiのはてなブログ記事（vec-oct/01-reflection自身のrefs） |

`refs-file.txt`の多対1（揺れ）ケースについて、出現数の多数決や記事内容との一致度から検討した正準slug案（未確定）:

| 対象記事 | 候補（出現数） | 正準案 | 根拠（案） |
|---|---|---|---|
| hopf/01-quaternion.md | 7shi-h(2), 7shi-hopf1(1), 7shi-qhopf1(1), 7shi(1)※衝突 | 7shi-h | 多数決 |
| hopf/02-spinor-tensor.md | 7shi-s(1), 7shi-2(1) | 7shi-s | s=spinorで内容に即す |
| misc/spherical-coords.md | 7shi-c(1), 7shi-coord(1), 7shi-pcoord(1) | 7shi-coord | 最も自己記述的 |
| qua/01-pauli-qua.md | 7shi-bq(2), 7shi-pbq(1), 7shi-pdq(2) | 7shi-bq | bq=双四元数、記事タイトルと一致 |
| qua/cd/heuristic-cmatrix.md | 7shi-qmat(1), 7shi-qcm(1) | 7shi-qcm | ファイル名cmatrixと一致 |
| oct/02-7d-3rot.md | oct-7rot(1), 7shi-ort(1) | oct-7rot | oct/系は`oct-`接頭辞の家系ルールに合わせる |
| vec-oct/01-reflection.md | 7shi-mir(1), 7shi-vo(1) | 7shi-mir | mir=鏡映、ファイル名reflectionと一致 |

観察: `oct/`系の記事は既に`oct`, `oct-7rot`, `oct-cl6`, `oct-nonassoc`と`7shi-`を付けない家系ルールがあるように見える。これをそのまま採用ルール化できそう。

上表の案を採用すると衝突`7shi`は自動解消（misc/spherical-coordsが`7shi-h`に変わるため）。残り`wikipedia`と`7shi-hc`は別途どちらをリネームするか要判断（未確認）。

### 2. 統合マスターレジストリの作成

重複解決後、`refs/*.toml`（記事ごとにバラバラ）を1つの `slug → 定義（type, citation, url等）` のマスターテーブルに統合する。namespaceは1つ（外部サイト／論文／書籍と、自分たちの他記事へのmathlog参照が同じslug空間を共有する）。

- 出力先ファイル名は未定（例: `refs/master.toml`）。
- 重複解決で決めた正準slugをキーとする。

### 3. 新記事執筆時のワークフロー支援

新しいチェックスクリプト（例: `check_slugs.py`）を作り、対象mdファイルを渡すと以下を行うようにする。

1. 記事自身の正準slugをマスターレジストリに登録する（Mathlog URL確定前は仮登録も許容）。
2. 本文中の `[[slug]]` を抽出する（`refs_slugs.py` の `SLUG_RE` を流用可能）。
3. マスターレジストリと突き合わせ、定義済みslugは定義を表示、未定義slugは一覧表示して定義（type/citation/url）の追加を促す。

## 関連ファイル

- `refs_slugs.py` — 上記の検出スクリプト本体
- `refs.toml`, `refs-url.txt`, `refs-file.txt` — 生成物（`uv run refs_slugs.py` で再生成可能）
- `refs/*.toml` — Mathlog記事ごとの参考文献エクスポート（`refs_to_toml.py` で `refs/*.html` から生成）
- `articles.tsv` — 記事一覧（date, url, md, title）。`articles.py` で生成・更新
