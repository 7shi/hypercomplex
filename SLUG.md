# slug管理: 全般ルールと運用

`SLUG.md` / `SLUG2.md` から、今後も有効な一般ルールのみを抜き出したもの。完了済みの一括修正（旧slug→正準slugの対応表など）とその調査結果は含まない。

## 背景

各記事はMathlogで公開されており、Mathlogの参考文献パネル（記事ごとの引用ラベル一覧）を `refs/{記事ID}.toml` としてエクスポートしてある（`reftools format` + `reftools toml` で HTML → TOML 変換）。この参照ラベル（slug）はMathlog側で記事ごとに独立して割り当てられているため、同じ実体を指すのに記事ごとに異なるslugが使われる「揺れ」が起こりうる。

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
- 旧slugが見つかった場合は**Mathlog側を手動修正する**（本文の `[[slug]]` と参考文献登録のラベル）。ローカルのmd本文と `refs/*.toml` も同じ内容に合わせて修正する。この手動修正作業は `mathlog_fix-refs.md`・`mathlog_fix.md` / `src/mathlog_fix.sh` で支援する（詳細は関連ファイル・ツール参照）。
- 同一slugが異なる対象を指す「衝突」は個別に解決が必要。同一対象を指す重複定義（表記揺れ）は統合可能。

## 外部文献の統合管理（refs-master.toml）

`refs/*.toml`（記事ごとにバラバラ）は `reftools build` で `refs.toml` に集約されるが、これはあくまで既存記事からの機械的な集約結果であり、正準定義そのものではない。外部文献（type != "mathlog"）の正準定義は `refs-master.toml` に集約し、手動で維持する。

- slugをキーとする。aliasは持たない。
- `type`/`url`/`author`/`site`/`journal`/`year`/`pages`/`publisher`等に加え、`citation`/`accessed`の代わりに `title` を持つ（`citation` は閲覧日入りの文字列でありマスター管理に不向きなため）。
- `type` ごとに持つべきフィールドの組み合わせは以下の通り（記事公開後、`refs.toml`の`citation`からこの組み合わせで`title`を機械的に復元できるかを `reftools check` が検証するため、実体と一致させる必要がある）。
  - `website`: `url`, `author`（任意）, `site`（任意）
  - `paper`: `author`, `journal`, `year`, `pages`（任意）
  - `book`: `author`, `publisher`, `year`, `pages`
- 自著Mathlog記事へのslug（type = "mathlog"）はここには含めない。定義はMathlog自身が正なので、二重管理しない。
- 新規の外部参照が未公開記事本文に現れたら、`reftools check` の「未定義」検出をきっかけに、この`refs-master.toml`へ手で追記する。
- `refs.toml`（機械生成）と矛盾していないかは `reftools check` が自動検証する（後述）。

## 更新手順

記事の公開・修正や参照の追加に伴ってrefsを更新するときは、次の順で行う。

1. **Mathlogの参考文献パネルを取り込む**（公開・修正した記事がある場合）
   - 新規公開は `bash src/mathlog_new.sh <md>` で行う（操作の手順は [README.md](README.md) の「新規公開」）。投稿後にクリップボードの参考文献パネルを一時ファイル `mathlog_new-refs.html` に保存し、記事一覧（手順2）を取り込んでURLが決まったら `refs/{ID}.html` に移し、整形・TOML変換から `reftools build` / `sync` / `check` まで行う。
   - 参考文献パネルのHTMLを `refs/{記事ID}.html` として保存し、`uv run reftools format --in-place` で整形してから `uv run reftools toml` で `refs/{ID}.toml` を生成する。
   - 本文や参考文献パネルの修正を伴う場合は、リポジトリ直下に `mathlog_fix-refs.md` を書いて `bash src/mathlog_fix.sh --refs` を実行する（HTMLの取り込み・整形・TOML変換に続けて `reftools build` / `check` まで行う）。処理後に `mathlog_fix-refs.md` を削除する。参考文献パネルに変更がない本文のみの修正は `mathlog_fix.md` に書いて `--no-refs` で実行する。
2. **記事一覧を更新する**
   - Mathlogの記事一覧が変わった場合は `make fetch` で `mathlog.tsv` に差分を取り込む（`mathlog_new.sh` はこれを含む）。全件を取り直すときは、記事一覧ページでブックマークレット `src/bookmarklets/mathlog_articles.url` を実行し、`winclip -o mathlog.tsv` で保存する。
   - `make md` → `make merge` で `md.tsv` / `articles.tsv` を更新する。
3. **集約・検証する**
   - `make build` で `refs.toml` を再生成する。
   - `make check` で検証する。
4. **指摘に対応する**
   - 未定義の外部文献は `refs-master.toml` に手で追記する。
   - 未登録の記事ファイルは `slugs.tsv` に正準slug（公開予定がなければ `NONE`）を登録する。
   - 旧slug・衝突はMathlog側を手動修正し、手順1からやり直す。
5. **同期する**
   - `make sync` で `refs-master.toml` の各slugの `files` を `refs.toml` に合わせる。

手順2（`make fetch` を除く）〜5は `make all`（`md merge build sync`）でまとめて実行できる。

機械生成されるファイル（`mathlog.tsv`・`md.tsv`・`articles.tsv`・`refs.toml`、`refs-master.toml` の `files`）は直接編集しない。内容を変えたいときは生成元（Mathlogの記事一覧・各 `README.md`・`refs/*.toml` 等）を直してから再生成する。

## 関連ファイル・ツール

- `slugs.tsv` — 正準slug台帳（md → slug）。手動管理。
- `refs-master.toml` — 外部文献（type != "mathlog"）の正準定義。手動管理、自動生成しない。
- `reftools`（`src/reftools/` パッケージ、`uv run reftools`で呼び出す） — 以下のサブコマンドを持つ。
  - `build`: `refs.toml` を生成する。slugごとに見出しを立て、`type`/`url`（定義がなければ省略）と使用元mdファイル一覧（`files`）を持つ。同一slugが異なる`(type, url)`に解決される場合はビルド時にエラーとする。`--url-output`/`--file-output` で `refs-url.txt`（複数slugから引用される同一URLの検出）/ `refs-file.txt`（自著記事が他記事から引用されているslugの一覧）を追加生成できるが、これらは通常のbuildには含めず、必要なチェック時にのみ都度生成する使い捨てファイル。
  - `check`: 以下をすべて実行する。
    1. 各公開済み記事本文の `[[slug]]` と対応する `refs/{ID}.toml` の過不足を検査する（「未使用」＝tomlにあるが本文にないslugは意図的な保持もあり、必ずしも修正対象ではない）。
    2. `md.tsv` の全ファイルが `slugs.tsv` に登録されているか確認する。
    3. 未公開記事（`refs/{ID}.toml` を持たない）本文の `[[slug]]` を、プロジェクト全体で既知のslug（全公開記事の `refs/{ID}.toml` のキー ∪ `refs-master.toml` のキー ∪ `slugs.tsv` に予約された正準slug）と突き合わせ、どこにも定義のないslugを報告する。
    4. `refs.toml`（機械生成）と `refs-master.toml`（手動管理）を突き合わせ、titleおよびその他共有フィールドの矛盾を報告する。
  - `sync`: `check` と同じ検査を行った上で、`refs-master.toml` の各slugの `files` を `refs.toml` に合わせて書き換える。
  - `show <mdパス|slug>`: 引数（mdパスまたは`slugs.tsv`の正準slug）で記事を1つ特定し、ファイル名と正準slug（`slugs.tsv`未登録なら未登録である旨）を表示した上で、本文の `[[slug]]` を出現順に列挙する。各slugは `refs-master.toml`（優先）または `refs.toml` の内容があれば表示し、どちらもなければ `slugs.tsv` の逆引き（他記事の正準slugであれば、未公開でもそのmdパス。公開済みなら `articles.tsv` から取得したURLも併記）を試し、それもなければ「情報なし」と表示する。
- `refs/*.toml` — Mathlog記事ごとの参考文献エクスポート（`reftools toml` で `refs/*.html` から生成）。
- `mathlog_fix-refs.md` / `mathlog_fix.md` — Mathlog側の記事を手動修正する際の作業リスト。参考文献パネルの修正を伴う記事は `mathlog_fix-refs.md`（`--refs` 用）、本文のみの修正は `mathlog_fix.md`（`--no-refs` 用）に書く。2つを分けることで、参照の変更が確定せず反映を保留している記事と、すぐに反映できる記事を同時に管理できる。何を直すか（表記ゆれの統一先、slug衝突の解消、ラベル改名など）は文脈依存のヒューリスティックな判断が必要で自動生成できないため、都度手で書く使い捨てファイル（処理後に削除する）。書式は両者で共通で、記事ごとに `## <mdファイル> — <Mathlog記事URL>` を見出しとし、その下に修正内容を箇条書きする（例: `- <slug>: <field> = <新値>`、改名は `- <旧slug> → <新slug>`）。
- `src/mathlog_fix.sh` — リポジトリ直下で実行し、作業リスト（`--refs` なら `mathlog_fix-refs.md`、`--no-refs` なら `mathlog_fix.md`）を見出しごとのブロックに分割し、各ブロックで対象mdファイルをエディタで開き、Mathlog記事URLをクリップボードにコピーし、ブロック本文（修正内容）を表示して手動修正の完了を待つ。完了後 `refs/{ID}.html` をクリップボードから取得し、`reftools format --in-place` で整形して `reftools toml` で `refs/{ID}.toml` に変換する。全ブロックの処理後に `reftools build` と `reftools check` を実行する。`--refs` と `--no-refs` のどちらかの指定が必須で、省略するとusageを表示して終了する。`--no-refs` では `mathlog_fix.md` を読み、クリップボードからの取り込み・整形・TOML変換を省く（`reftools build` / `check` は実行する）。
- `articles.tsv` — 記事一覧（date, url, md, title）。`src/articles/`（`articles`コマンド）で生成・更新。
- `md.tsv` — 全記事ファイル一覧（md, title）。公開・未公開を問わず全ファイルを含む。
- `src/mathlog_new.sh` — 未公開の記事1本を新規投稿する手順を順に案内する（手順は上記「更新手順」の1）。事前に `md.tsv`・`slugs.tsv` への登録と未公開であることを確かめ、`reftools show` で登録すべき参考文献を表示する。
- `mathlog.url` — Mathlogの記事一覧ページのURL。`articles fetch` が取得する。
- `src/bookmarklets/mathlog_articles.url` — Mathlogの記事一覧ページで実行し、表示中の記事を `mathlog.tsv` の形式（日時・URL・タイトル）でクリップボードにコピーする。日時は各記事の `created_at` を日本時間に直したもの（`yyyy/mm/dd hh:mm:ss`）。全件を取り直すときに、全記事を表示させてから実行する。
- `Makefile` — `make fetch`/`md`/`merge`（`articles` の同名サブコマンド）、`make build`/`sync`/`check`（`reftools` の同名サブコマンド）、`make all`（`md merge build sync`）のショートカット。`make help` で一覧表示。
