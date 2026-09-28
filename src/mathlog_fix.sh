#!/bin/bash
set -eu

no_refs=false
for arg in "$@"; do
    case "$arg" in
        --no-refs) no_refs=true ;;
        *) echo "usage: $0 [--no-refs]" >&2; exit 1 ;;
    esac
done

mapfile -t lines < mathlog_fix.md

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
        read -p "修正が完了したら[Enter]を押してください。"
        if $no_refs; then
            return
        fi
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
