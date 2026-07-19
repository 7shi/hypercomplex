# slug管理: 全般ルールと残作業

`SLUG.md` / `SLUG2.md` から、今後も有効な一般ルールと未完了の作業のみを抜き出したもの。完了済みの一括修正（旧slug→正準slugの対応表など）とその調査結果は含まない。

## 背景

各記事はMathlogで公開されており、Mathlogの参考文献パネル（記事ごとの引用ラベル一覧）を `refs/{記事ID}.toml` としてエクスポートしてある（`html_format.py` + `refs_to_toml.py` で HTML → TOML 変換）。この参照ラベル（slug）はMathlog側で記事ごとに独立して割り当てられているため、同じ実体を指すのに記事ごとに異なるslugが使われる「揺れ」が起こりうる。

## 命名規約

1. **自著記事**（Mathlog・はてなブログ・Colab/Gist）は `7shi-<略号>`。略号は内容由来（`mir`, `coord`）または連載番号（`lie1`, `lie2`）。例外は認めない。
2. **Wikipedia** は `wiki-<略号>`。裸の `wikipedia` および `wikipedia-<略号>` は不可。
3. **その他外部文献**は著者名（`nakajima`, `sylvester`, `tao-h`）または内容・サイト略号（`dim7-8`, `half`, `proj`）。`ref1` のような無意味な連番は不可。
4. 対象を特定しない汎用slug（`7shi` 単独、`wikipedia` 単独）は不可。
5. 新規割当の略号は自己記述的にする（原則3文字以上）。単文字略号（`7shi-h`, `7shi-s` 等）は既存のみ許容。
6. 新規slugはファイル名の内容部分を略したものを基本とする（連番・ディレクトリは含めない。例: `bloch-density`→`bloch`）。確定後はファイル名が変わってもslugは変更しない（slugは公開済み記事に埋め込まれるため、ファイル名より強い安定性が必要）。

## 正準slug管理の運用ルール

- 記事ファイルごとの正準slugはローカルの `slugs.tsv`（md → slug）で管理する。Mathlogには記事にslugを割り振る機能がないため、これが唯一の台帳となる。`slugs.tsv` は `md.tsv` の全ファイルを対象とし、公開予定のないファイルは `NONE` を残置する。
- 旧slugが見つかった場合は**Mathlog側を手動修正する**（本文の `[[slug]]` と参考文献登録のラベル）。ローカルのmd本文と `refs/*.toml` も同じ内容に合わせて修正する。この手動修正作業は `mathlog_fix.md` / `mathlog_fix.sh` で支援する（詳細は関連ファイル・ツール参照）。
- 同一slugが異なる対象を指す「衝突」は個別に解決が必要。同一対象を指す重複定義（表記揺れ）は統合可能。

## 残作業

### 1. 統合マスターレジストリの作成

`refs/*.toml`（記事ごとにバラバラ）を1つの `slug → 定義（type, citation, url等）` のマスターテーブルに統合する。namespaceは1つ（外部サイト／論文／書籍と、自分たちの他記事へのmathlog参照が同じslug空間を共有する）。

- 出力先ファイル名は未定（例: `refs/master.toml`）。
- 正準slugをキーとする。aliasは持たない。
- 自著Mathlog記事のエントリには対象mdパス（`md = "..."`）も持たせる（`slugs.tsv` と整合させる）。
- 統合時、同一slugに異なる対象が混ざったらエラーにする（衝突検出）。
- 未公開記事本文に残る未定義の外部参照（`wiki-ccr`, `bernevig`, `nielsen-chuang` 等、命名は規約適合済み）は、この統合時に定義を登録する。

### 2. 新記事執筆時のワークフロー支援

新しいチェックスクリプト（例: `check_slugs.py`）を作り、対象mdファイルを渡すと以下を行うようにする。

1. 記事自身の正準slugをマスターレジストリに登録する（Mathlog URL確定前は仮登録も許容）。
2. 本文中の `[[slug]]` を抽出する（`refs_slugs.py` の `SLUG_RE` を流用可能）。
3. マスターレジストリと突き合わせ、定義済みslugは定義を表示、未定義slugは一覧表示して定義（type/citation/url）の追加を促す。

## 関連ファイル・ツール

- `slugs.tsv` — 正準slug台帳（md → slug）。手動管理。
- `refs_slugs.py` — `build` サブコマンドで `refs.toml` を生成する。slugごとに見出しを立て、`type`/`url`（定義がなければ省略）と使用元mdファイル一覧（`files`）を持つ。同一slugが異なる`(type, url)`に解決される場合はビルド時にエラーとする。`--url-output`/`--file-output` で `refs-url.txt`（複数slugから引用される同一URLの検出）/ `refs-file.txt`（自著記事が他記事から引用されているslugの一覧）を追加生成できるが、これらは通常のbuildには含めず、必要なチェック時にのみ都度生成する使い捨てファイル。`check` サブコマンドは各記事本文の `[[slug]]` と対応する `refs/{ID}.toml` の過不足を検査する（「未使用」＝tomlにあるが本文にないslugは意図的な保持もあり、必ずしも修正対象ではない）。
- `refs/*.toml` — Mathlog記事ごとの参考文献エクスポート（`refs_to_toml.py` で `refs/*.html` から生成）。
- `mathlog_fix.md` — Mathlog側の参考文献パネルを手動修正する際の作業リスト。何を直すか（表記ゆれの統一先、slug衝突の解消、ラベル改名など）は文脈依存のヒューリスティックな判断が必要で自動生成できないため、都度手で書く使い捨てファイル（処理後に削除する）。記事ごとに `## <mdファイル> — <Mathlog記事URL>` を見出しとし、その下に修正内容を箇条書きする（例: `- <slug>: <field> = <新値>`、改名は `- <旧slug> → <新slug>`）。
- `mathlog_fix.sh` — `mathlog_fix.md` を見出しごとのブロックに分割し、各ブロックで対象mdファイルをエディタで開き、Mathlog記事URLをクリップボードにコピーし、ブロック本文（修正内容）を表示して手動修正の完了を待つ。完了後 `refs/{ID}.html` をクリップボードから取得し、`html_format.py --in-place` で整形する。
- `articles.tsv` — 記事一覧（date, url, md, title）。`articles.py` で生成・更新。
