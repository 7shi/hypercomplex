.PHONY: help md merge check

help:
	@echo "使用可能なターゲット:"
	@echo "  md     - uv run articles.py md"
	@echo "  merge  - uv run articles.py merge"
	@echo "  check  - uv run refs_slugs.py check"

md:
	uv run articles.py md

merge:
	uv run articles.py merge

check:
	uv run refs_slugs.py check
