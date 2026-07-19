#!/bin/bash
set -eu

mapfile -t lines < mathlog_fix.md

for line in "${lines[@]}"; do
    rest="${line#\#\# }"
    file="${rest%% — *}"
    url="${rest#* — }"
    ref="refs/$(basename "$url").html"

    code "$file"
    winclip "$url"
    read -p "$line"
    winclip -o "$ref"
    uv run html_format.py "$ref" --in-place
done
