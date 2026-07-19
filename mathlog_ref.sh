#!/bin/bash
set -eu

mapfile -t lines < mathlog.tsv

for line in "${lines[@]:1}"; do
    IFS=$'\t' read -r date url title <<< "$line"
    ref="refs/$(basename "$url").html"
    [ -e "$ref" ] && continue

    read -p "https://mathlog.info$url $title"
    winclip -o "$ref"
    uv run html_format.py "$ref" --in-place
done
