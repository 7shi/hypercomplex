.PHONY: help all mathlog md merge build check sync

help:
	@echo "使用可能なターゲット:"
	@echo "  all    - md merge build check"
	@echo "  md     - uv run articles.py md"
	@echo "  merge  - uv run articles.py merge"
	@echo "  build  - uv run reftools build"
	@echo "  check  - uv run reftools check"
	@echo "  sync   - uv run reftools sync"

all: md merge build check

mathlog md merge:
	uv run articles.py $@

build check sync:
	uv run reftools $@
