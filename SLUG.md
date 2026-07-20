# slug管理: 全般ルールと運用

`SLUG.md` / `SLUG2.md` から、今後も有効な一般ルールのみを抜き出したもの。完了済みの一括修正（旧slug→正準slugの対応表など）とその調査結果は含まない。

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

- 記事ファイルごとの正準slugはローカルの `slugs.tsv`（md → slug）で管理する。Mathlogには記事にslugを割り振る機能がないため、これが唯一の台帳となる。`slugs.tsv` は `md.tsv` の全ファイルを対象とし、公開予定のないファイルは `NONE` を残置する。未公開の他記事を `[[slug]]` で参照する場合は、この台帳に予約された正準slugを使う（公開を待たずに参照できる）。
- 外部文献（Mathlog以外）の正準定義は `refs-master.toml` で手動管理する（後述）。
- 旧slugが見つかった場合は**Mathlog側を手動修正する**（本文の `[[slug]]` と参考文献登録のラベル）。ローカルのmd本文と `refs/*.toml` も同じ内容に合わせて修正する。この手動修正作業は `mathlog_fix.md` / `mathlog_fix.sh` で支援する（詳細は関連ファイル・ツール参照）。
- 同一slugが異なる対象を指す「衝突」は個別に解決が必要。同一対象を指す重複定義（表記揺れ）は統合可能。

## 外部文献の統合管理（refs-master.toml）

`refs/*.toml`（記事ごとにバラバラ）は `refs_slugs.py build` で `refs.toml` に集約されるが、これはあくまで既存記事からの機械的な集約結果であり、正準定義そのものではない。外部文献（type != "mathlog"）の正準定義は `refs-master.toml` に集約し、手動で維持する。

- slugをキーとする。aliasは持たない。
- `type`/`url`/`author`/`site`/`journal`/`year`/`pages`/`publisher`等に加え、`citation`/`accessed`の代わりに `title` を持つ（`citation` は閲覧日入りの文字列でありマスター管理に不向きなため）。
- `type` ごとに持つべきフィールドの組み合わせは以下の通り（記事公開後、`refs.toml`の`citation`からこの組み合わせで`title`を機械的に復元できるかを `refs_slugs.py check` が検証するため、実体と一致させる必要がある）。
  - `website`: `url`, `author`（任意）, `site`（任意）
  - `paper`: `author`, `journal`, `year`, `pages`（任意）
  - `book`: `author`, `publisher`, `year`, `pages`
- 自著Mathlog記事へのslug（type = "mathlog"）はここには含めない。定義はMathlog自身が正なので、二重管理しない。
- 新規の外部参照が未公開記事本文に現れたら、`refs_slugs.py check` の「未定義」検出をきっかけに、この`refs-master.toml`へ手で追記する。
- `refs.toml`（機械生成）と矛盾していないかは `refs_slugs.py check` が自動検証する（後述）。

## 関連ファイル・ツール

- `slugs.tsv` — 正準slug台帳（md → slug）。手動管理。
- `refs-master.toml` — 外部文献（type != "mathlog"）の正準定義。手動管理、自動生成しない。
- `refs_slugs.py` — 以下のサブコマンドを持つ。
  - `build`: `refs.toml` を生成する。slugごとに見出しを立て、`type`/`url`（定義がなければ省略）と使用元mdファイル一覧（`files`）を持つ。同一slugが異なる`(type, url)`に解決される場合はビルド時にエラーとする。`--url-output`/`--file-output` で `refs-url.txt`（複数slugから引用される同一URLの検出）/ `refs-file.txt`（自著記事が他記事から引用されているslugの一覧）を追加生成できるが、これらは通常のbuildには含めず、必要なチェック時にのみ都度生成する使い捨てファイル。
  - `check`: 以下をすべて実行する。
    1. 各公開済み記事本文の `[[slug]]` と対応する `refs/{ID}.toml` の過不足を検査する（「未使用」＝tomlにあるが本文にないslugは意図的な保持もあり、必ずしも修正対象ではない）。
    2. `md.tsv` の全ファイルが `slugs.tsv` に登録されているか確認する。
    3. 未公開記事（`refs/{ID}.toml` を持たない）本文の `[[slug]]` を、プロジェクト全体で既知のslug（全公開記事の `refs/{ID}.toml` のキー ∪ `refs-master.toml` のキー ∪ `slugs.tsv` に予約された正準slug）と突き合わせ、どこにも定義のないslugを報告する。
    4. `refs.toml`（機械生成）と `refs-master.toml`（手動管理）を突き合わせ、titleおよびその他共有フィールドの矛盾を報告する。
  - `show <mdパス|slug>`: 引数（mdパスまたは`slugs.tsv`の正準slug）で記事を1つ特定し、ファイル名と正準slug（`slugs.tsv`未登録なら未登録である旨）を表示した上で、本文の `[[slug]]` を出現順に列挙する。各slugは `refs-master.toml`（優先）または `refs.toml` の内容があれば表示し、なければ「情報なし」と表示する。
- `refs/*.toml` — Mathlog記事ごとの参考文献エクスポート（`refs_to_toml.py` で `refs/*.html` から生成）。
- `mathlog_fix.md` — Mathlog側の参考文献パネルを手動修正する際の作業リスト。何を直すか（表記ゆれの統一先、slug衝突の解消、ラベル改名など）は文脈依存のヒューリスティックな判断が必要で自動生成できないため、都度手で書く使い捨てファイル（処理後に削除する）。記事ごとに `## <mdファイル> — <Mathlog記事URL>` を見出しとし、その下に修正内容を箇条書きする（例: `- <slug>: <field> = <新値>`、改名は `- <旧slug> → <新slug>`）。
- `mathlog_fix.sh` — `mathlog_fix.md` を見出しごとのブロックに分割し、各ブロックで対象mdファイルをエディタで開き、Mathlog記事URLをクリップボードにコピーし、ブロック本文（修正内容）を表示して手動修正の完了を待つ。完了後 `refs/{ID}.html` をクリップボードから取得し、`html_format.py --in-place` で整形する。
- `articles.tsv` — 記事一覧（date, url, md, title）。`articles.py` で生成・更新。
- `md.tsv` — 全記事ファイル一覧（md, title）。公開・未公開を問わず全ファイルを含む。
- `Makefile` — `make md`（`articles.py md`）、`make merge`（`articles.py merge`）、`make check`（`refs_slugs.py check`）のショートカット。`make help` で一覧表示。
