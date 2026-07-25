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
    uv run reftools format "$ref" --in-place
    uv run reftools toml "$ref"
done
uv run reftools build
uv run reftools check
