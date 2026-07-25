# reftool

記事中の`[[slug]]`参照マーカーを、`articles.tsv`が指すmdファイルから集めて1つのTOMLファイルに集約するツール。

## build

Mathlogのurlを持つ各mdファイル（`articles.tsv`でパスがソートされる）について、記事ごとの参照リスト`refs/{ID}.toml`を読み込み、mdファイル中に見つかった各`[[slug]]`マーカーをそのtomlの全属性に解決する。

- mathlog以外のtype: `type`, `citation`, `url`, `accessed`, `author`, `site`, `journal`, `pages`, `year`, `publisher`など、属性をそのまま採用する。
- type = "mathlog": 同様の属性に加え、`citation`のタイトルを`articles.tsv`で検索し、参照先記事自身のurlを明示的な`url`として追加する。

出力はファイル単位ではなくslug単位でグループ化する。各slugは`[slug]`見出しの下に統合済みの属性と、そのslugを使う全mdファイルの`files`リストを持つ。slugはプロジェクト全体で正準であるべきなので、`citation`/`accessed`以外の全属性は同じslugを定義する全ファイルで完全に一致していなければならない。不一致は命名衝突とみなす。全slugをチェックしてから報告する（1つのslugでの衝突が残りのスキャンを止めることはない）。`refs.toml`は衝突があっても書き出される（衝突しているslugは、各ファイルの中で最新の`accessed`日付を持つ定義にフォールバックする）が、その後で見つかった衝突がすべてまとめて表示され、コマンドは非0で終了する。`refs/{ID}.toml`を手で直せるように、衝突箇所を確認すること。`citation`/`accessed`フィールドは、エクスポート日によって正当に異なりうるため、この衝突チェック自体からは除外される。

Mathlogのurlを持たないmdファイル（したがって`refs/{ID}.toml`も持たない）でも、その`[[slug]]`マーカーは未解決のまま`files`リストに寄与する。`[[slug]]`マーカーを持たないファイルは省略される。

`refs-url.txt`も書き出す。2つ以上の異なるslugで参照されているurl（同じ出典が別々の記事で違うslug名で引用されている状態）を、urlでソートして「url = slug, slug, ...」の形式で1行ずつ列挙する。

`refs-file.txt`も書き出す。それ自身がMathlog記事であるmdファイルごとに、他の記事がtype = "mathlog"でそれを引用する際に使っている異なるslugを列挙する（slugのブレを見つけ、対象記事1件につき正準slugを1つ選ぶのに便利）。mdパスでソートし、「md_path = slug, slug, ...」または引用が1つもない場合は「md_path = NONE」の形式で1行ずつ出力する。

## check

公開済みの各記事（`refs/{ID}.toml`を持つもの）について、md本文で使われている`[[slug]]`マーカーと`refs/{ID}.toml`で定義されているキーを比較し、不一致（本文にあるがtomlにない、tomlにあるが本文にない）を報告する。

`md.tsv`に列挙された全ファイルが`slugs.tsv`（手動管理のmdパス→正準slug台帳、`SLUG.md`参照）に登録されているかどうかも確認し、未登録のファイルを報告する。

未公開の記事（Mathlogのurlがまだなく、自分の`refs/{ID}.toml`を持たない）については、代わりにその`[[slug]]`マーカーを、プロジェクト全体で既知の全slug（全公開記事の`refs/{ID}.toml`のキー、`refs-master.toml`のキー、`slugs.tsv`の予約済みslugの和集合）と比較し、どこにも未定義のものを報告する。これはタイプミスや未登録の参照を捕捉するだけである（比較対象のファイル単位のtomlが存在しないため、「未使用」側の検出はできない）。`slugs.tsv`があることで、まだ未公開の姉妹記事を引用する`[[slug]]`が、公開を待たずにその記事の予約済み正準slugを通じて早期に解決できる。

`refs.toml`と`refs-master.toml`の相互チェックも行う。`refs-master.toml`は非mathlog系slugすべてを網羅する手動管理のサブセットである（mathlogエントリはその出典データがMathlog自体にあるため除外される）。`refs-master.toml`が定義する各slugについて、`refs.toml`が同じslugを定義していること、その`citation`が同じ`title`を暗示していること（website型なら`author`/`site`/`accessed`/`url`から、paper型なら`journal`/`year`/`pages`から、book型なら`publisher`/`year`/`pages`から`title`を復元する同じ方法で）、その他の共有フィールドがすべて完全に一致することを確認する。不一致は表示される。`refs-master.toml`自体はこのスクリプトによって再生成されることはない。

## sync

`check`と同じチェックを行った上で、`refs-master.toml`内の各slugの`files`リストを`refs.toml`のものに合わせて書き換える。他のフィールドや`refs-master.toml`自身のslug順序・整形はそのまま保つ。`refs.toml`に存在しないslugはそのまま放置する（マスター不一致として`check`側ですでに報告されている）。

## show

1つの記事（mdパスまたは`slugs.tsv`の正準slugで指定）について、本文中に現れる`[[slug]]`マーカーを最初に登場した順に表示し、それぞれについて次の優先順位で情報源を探して表示する。

1. `refs-master.toml`に定義があればそのエントリ（優先）
2. `refs.toml`に何らかの解決済み属性を持つエントリがあればそれ
3. そのslugが別の（未公開の可能性もある）記事自身の正準slugであれば、`slugs.tsv`の逆引きでmdパスを表示（公開済みなら`articles.tsv`のMathlog urlも添える）
4. どこにも情報が見つからなければ、その旨を表示する

## 使い方

```
uv run reftool build                      # refs.tomlのみを書き出す
uv run reftool build -o out.toml
uv run reftool build --url-output refs-url.txt --file-output refs-file.txt
uv run reftool check
uv run reftool sync
uv run reftool show oct/01-octonion.md
uv run reftool show 7shi-oct1
```
