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
        uv run html_format.py "$ref" --in-place
        uv run refs_to_toml.py "$ref"
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
uv run refs_slugs.py build
uv run refs_slugs.py check
