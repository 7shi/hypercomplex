"""Pretty-print an HTML fragment/file with one tag per line.

Elements whose content is text-only (no element children) are kept on a
single line together with their closing tag; elements with element
children are expanded, with each child indented one level deeper.

Usage:
  html_format.py input.html            # print formatted result to stdout
  html_format.py input.html -o out.html
  html_format.py input.html --in-place
"""

from __future__ import annotations

import argparse
from html.parser import HTMLParser
from pathlib import Path

VOID_TAGS = {
    "area", "base", "br", "col", "embed", "hr", "img", "input",
    "link", "meta", "param", "source", "track", "wbr",
}

INDENT = "  "


class Element:
    def __init__(self, tag: str, attrs: list[tuple[str, str | None]]):
        self.tag = tag
        self.attrs = attrs
        self.children: list["Element | Text"] = []

    def open_tag(self) -> str:
        parts = [self.tag]
        for name, value in self.attrs:
            if value is None:
                parts.append(name)
            else:
                parts.append(f'{name}="{value}"')
        return f"<{' '.join(parts)}>"

    def close_tag(self) -> str:
        return f"</{self.tag}>"


class Text:
    def __init__(self, data: str):
        self.data = data


class TreeBuilder(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.root = Element("#root", [])
        self.stack: list[Element] = [self.root]

    def handle_starttag(self, tag, attrs):
        elem = Element(tag, attrs)
        self.stack[-1].children.append(elem)
        if tag not in VOID_TAGS:
            self.stack.append(elem)

    def handle_startendtag(self, tag, attrs):
        self.stack[-1].children.append(Element(tag, attrs))

    def handle_endtag(self, tag):
        for i in range(len(self.stack) - 1, 0, -1):
            if self.stack[i].tag == tag:
                del self.stack[i:]
                return

    def handle_data(self, data):
        self.stack[-1].children.append(Text(data))

    def handle_entityref(self, name):
        self.stack[-1].children.append(Text(f"&{name};"))

    def handle_charref(self, name):
        self.stack[-1].children.append(Text(f"&#{name};"))

    def handle_comment(self, data):
        self.stack[-1].children.append(Text(f"<!--{data}-->"))

    def handle_decl(self, decl):
        self.root.children.append(Text(f"<!{decl}>"))


def is_text_only(children: list["Element | Text"]) -> bool:
    return all(isinstance(c, Text) for c in children)


def flat_text(children: list["Element | Text"]) -> str:
    return " ".join(" ".join(c.data.split()) for c in children if isinstance(c, Text)).strip()


def render(node: Element, depth: int, out: list[str]) -> None:
    prefix = INDENT * depth
    if not node.children:
        if node.tag in VOID_TAGS:
            out.append(f"{prefix}{node.open_tag()}")
        else:
            out.append(f"{prefix}{node.open_tag()}{node.close_tag()}")
        return

    if is_text_only(node.children):
        text = flat_text(node.children)
        out.append(f"{prefix}{node.open_tag()}{text}{node.close_tag()}")
        return

    out.append(f"{prefix}{node.open_tag()}")
    for child in node.children:
        if isinstance(child, Text):
            text = " ".join(child.data.split())
            if text:
                out.append(f"{INDENT * (depth + 1)}{text}")
        else:
            render(child, depth + 1, out)
    out.append(f"{prefix}{node.close_tag()}")


def format_html(html: str) -> str:
    builder = TreeBuilder()
    builder.feed(html)
    builder.close()

    out: list[str] = []
    for child in builder.root.children:
        if isinstance(child, Text):
            text = " ".join(child.data.split())
            if text:
                out.append(text)
        else:
            render(child, 0, out)
    return "\n".join(out) + "\n"


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("input", type=Path, help="HTML file to format")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-o", "--output", type=Path, help="write result to this file")
    group.add_argument("--in-place", action="store_true", help="overwrite the input file")
    args = parser.parse_args(argv)

    html = args.input.read_text(encoding="utf-8")
    formatted = format_html(html)

    if args.in_place:
        args.input.write_text(formatted, encoding="utf-8")
    elif args.output:
        args.output.write_text(formatted, encoding="utf-8")
    else:
        print(formatted, end="")


if __name__ == "__main__":
    main()
