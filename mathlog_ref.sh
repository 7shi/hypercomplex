#!/bin/bash
set -eu

mapfile -t lines < mathlog.tsv

for line in "${lines[@]:1}"; do
    IFS=$'\t' read -r date url title <<< "$line"
    ref="refs/$(basename "$url").html"
    [ -e "$ref" ] && continue

    winclip "https://mathlog.info$url"
    read -p "$url $title"
    winclip -o "$ref"
    uv run html_format.py "$ref" --in-place
    uv run refs_to_toml.py "$ref"
done
uv run refs_slugs.py build
uv run refs_slugs.py check
