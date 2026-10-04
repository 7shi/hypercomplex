#!/bin/bash
set -eu

mapfile -t lines < mathlog.tsv

for line in "${lines[@]:1}"; do
    IFS=$'\t' read -r date url title <<< "$line"
    ref="refs/$(basename "$url").html"
    [ -e "$ref" ] && continue

    winclip "https://mathlog.info$url"
    prompt="$url $title"
    while :; do
        read -p "$prompt"
        refs=$(winclip -o)
        [[ "$refs" == *accordion* ]] && break
        prompt="クリップボードの内容が参考文献パネルではありません。コピーし直して[Enter]を押してください。"
    done
    printf '%s\n' "$refs" > "$ref"
    uv run reftools format "$ref" --in-place
    uv run reftools toml "$ref"
done
uv run reftools build
uv run reftools check
