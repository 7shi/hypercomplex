#!/bin/bash
set -eu

usage() {
    echo "usage: $0 --refs|--no-refs" >&2
    echo "  --refs     mathlog_fix-refs.md を処理し、修正後に参考文献をクリップボードから取り込む" >&2
    echo "  --no-refs  mathlog_fix.md を処理し、参考文献の取り込みを省く（本文のみの修正）" >&2
    exit 1
}

[ $# -eq 1 ] || usage
case "$1" in
    --refs) no_refs=false; list=mathlog_fix-refs.md ;;
    --no-refs) no_refs=true; list=mathlog_fix.md ;;
    *) usage ;;
esac

mapfile -t lines < "$list"

total=0
for line in "${lines[@]}"; do
    [[ "$line" == "## "* ]] && total=$((total + 1))
done

index=0
file=""
url=""
content=()

flush() {
    if [ -n "$file" ]; then
        echo
        index=$((index + 1))
        echo "## ($index/$total) $file - $url"
        code "$file"
        winclip "$url"
        printf '%s\n' "${content[@]}"
        if $no_refs; then
            read -p "修正が完了したら[Enter]を押してください。"
            return
        fi
        echo
        uv run reftools show "$file"
        echo
        read -p "修正が完了したら、参考文献をコピーして[Enter]を押してください。"
        ref="refs/$(basename "$url").html"
        winclip -o "$ref"
        uv run reftools format "$ref" --in-place
        uv run reftools toml "$ref"
    fi
}

for line in "${lines[@]}"; do
    if [[ "$line" == "## "* ]]; then
        flush
        rest="${line#\#\# }"
        file="${rest%% — *}"
        url="${rest#* — }"
        content=()
    elif [ -n "$line" ]; then
        content+=("$line")
    fi
done
flush
uv run reftools build
uv run reftools check
