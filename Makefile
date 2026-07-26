.PHONY: help all mathlog md merge build sync check

help:
	@echo "使用可能なターゲット:"
	@echo "  all    - md merge build check"
	@echo "  md     - uv run articles md"
	@echo "  merge  - uv run articles merge"
	@echo "  build  - uv run reftools build"
	@echo "  sync   - uv run reftools sync"
	@echo "  check  - uv run reftools check"

all: md merge build sync

mathlog md merge:
	uv run articles $@

build sync check:
	uv run reftools $@
