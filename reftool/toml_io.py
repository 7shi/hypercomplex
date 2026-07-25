from __future__ import annotations


def toml_string(s: str) -> str:
    escaped = s.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def toml_value(v) -> str:
    if isinstance(v, str):
        return toml_string(v)
    if isinstance(v, list):
        return "[" + ", ".join(toml_value(x) for x in v) + "]"
    return str(v)  # int, date: valid as-is in TOML


def render_toml_by_slug(slug_map: dict[str, dict]) -> str:
    blocks: list[str] = []
    for slug in sorted(slug_map):
        entry = slug_map[slug]
        lines = [f"[{toml_string(slug)}]"]
        for key, value in entry.items():
            if key == "files":
                continue
            lines.append(f"{key} = {toml_value(value)}")
        items = ",\n".join(f"  {toml_string(f)}" for f in entry["files"])
        lines.append(f"files = [\n{items},\n]")
        blocks.append("\n".join(lines) + "\n")
    return "\n".join(blocks)
