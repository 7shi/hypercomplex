.PHONY: help all md merge build check sync

help:
	@echo "使用可能なターゲット:"
	@echo "  all    - md merge build check"
	@echo "  md     - uv run articles.py md"
	@echo "  merge  - uv run articles.py merge"
	@echo "  build  - uv run refs_slugs.py build"
	@echo "  check  - uv run refs_slugs.py check"
	@echo "  sync   - uv run refs_slugs.py sync"

all: md merge build check

md:
	uv run articles.py md

merge:
	uv run articles.py merge

build:
	uv run refs_slugs.py build

check:
	uv run refs_slugs.py check

sync:
	uv run refs_slugs.py sync
