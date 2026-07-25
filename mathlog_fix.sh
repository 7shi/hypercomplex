#!/bin/bash
set -eu

mapfile -t lines < mathlog_fix.md

file=""
url=""
content=()

flush() {
    if [ -n "$file" ]; then
        code "$file"
        winclip "$url"
        printf '%s\n' "${content[@]}"
        read -p "修正が完了したら[Enter]を押してください。"
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
