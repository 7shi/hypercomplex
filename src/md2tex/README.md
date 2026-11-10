# md2tex

Mathlog記法のMarkdown記事をLaTeXに変換するツール。紙面での確認用で、Mathlogの表示を再現するものではない。

```
uv run md2tex clif/gca-schrodinger.md --pdf
```

`build/<記事名>/`（追跡対象外）に`<記事名>.tex`と図のPDFを書き出す。出力先は`-o`で変えられる。`--pdf`を付けると`lualatex`を2回実行してPDFまで作る。

## 必要なもの

- `rsvg-convert`：図（SVG）をPDFに変換する
- LuaLaTeXと`luatexja`（`ltjsarticle`）・`tcolorbox`：`--pdf`を付けるときだけ使う

## 変換の内容

- 題名：`md.tsv`のタイトルを使う。載っていなければファイル名とする
- 見出し：`#`・`##`・`###`を`\section`・`\subsection`・`\subsubsection`にする
- 囲み：`&&&type title`〜`&&&`を色付きの枠にし、「定義 1（title）」のように種類ごとの通し番号を付ける。種類のない囲みはtitleだけの枠にする。ラベル（`[label]`）は読み捨てる
- 数式：`$…$`と`$$`の行で挟んだ別行立ての数式をそのまま渡す。Unicodeのギリシャ文字はコマンドに置き換える（`ϕ`は`\phi`、`φ`は`\varphi`）
- 図：`![alt](images/xxx.svg)`をPDFに変換し、等倍（96 dpi換算）で取り込む。本文幅を超えるものは本文幅に縮める
- 参照：`[[slug]]`を初出順の番号にし、末尾に参考文献の一覧を付ける。引用文は`refs.toml`の`citation`、無ければ`refs-master.toml`の`title`、それも無ければ`slugs.tsv`と`md.tsv`から引いた記事のタイトルを使う
- その他：箇条書き（入れ子可）と`**太字**`

表・引用・コードブロック・番号付きリスト・`[text](url)`形式のリンクなどには対応していない。表・引用・コードブロック・番号付きリストは、その行を示して止まる。
