#!/bin/bash
set -eu

usage() {
    echo "usage: $0 <md>" >&2
    echo "  未公開の記事をMathlogに新規投稿し、記事一覧と参考文献をローカルに同期する" >&2
    exit 1
}

[ $# -eq 1 ] || usage
md="${1#./}"
[ -f "$md" ] || { echo "$md がありません" >&2; exit 1; }

# 投稿時に取った参考文献パネルを、URL（記事ID）が決まるまで置いておく一時ファイル。
# 1行目に対象のmdパスを書き、再実行時に照合する。
tmp=mathlog_new-refs.html

lookup_url() {
    awk -F'\t' -v md="$md" '$3 == md { print $2 }' articles.tsv
}

if [ -e "$tmp" ]; then
    tmp_md=$(head -n 1 "$tmp")
    if [ "$tmp_md" != "$md" ]; then
        echo "$tmp は $tmp_md のものです。処理を終えるか削除してから実行してください。" >&2
        exit 1
    fi
    echo "$tmp が残っているため、参考文献の取得を省いて再開します。"
else
    grep -q "^$md	" md.tsv || { echo "$md が md.tsv にありません（README.md にリンクを追加して make md）" >&2; exit 1; }
    grep -q "^$md	" slugs.tsv || { echo "$md が slugs.tsv に登録されていません" >&2; exit 1; }
    [ -z "$(lookup_url)" ] || { echo "$md は公開済みです: $(lookup_url)" >&2; exit 1; }

    echo "## $md - $(awk -F'\t' -v md="$md" '$1 == md { print $2 }' md.tsv)"
    code "$md"
    echo
    uv run reftools show "$md"
    echo
    echo "本文を貼り付けて参考文献を登録したら、参考文献パネルをコピーしてから投稿してください。"
    read -p "投稿が完了したら[Enter]を押してください。"
    refs=$(winclip -o)
    if [[ "$refs" != *accordion* ]]; then
        echo "クリップボードの内容が参考文献パネルではありません。参考文献パネルをコピーしてから再実行してください。" >&2
        exit 1
    fi
    printf '%s\n%s\n' "$md" "$refs" > "$tmp"
fi

url=$(lookup_url)
if [ -z "$url" ]; then
    uv run articles fetch
    uv run articles md
    uv run articles merge
    url=$(lookup_url)
    if [ -z "$url" ]; then
        echo "articles.tsv に $md のURLがありません。Mathlogのタイトルと README.md のリンク文字列を確認してください。" >&2
        head -n 4 mathlog.tsv >&2
        exit 1
    fi
fi
echo "URL: https://mathlog.info$url"

ref="refs/$(basename "$url").html"
tail -n +2 "$tmp" > "$ref"
uv run reftools format "$ref" --in-place
uv run reftools toml "$ref"
rm "$tmp"

uv run reftools build
uv run reftools sync
uv run reftools check
