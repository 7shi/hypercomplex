# src

記事執筆・検証・参考文献管理に使うスクリプト類をまとめたディレクトリです。

- [articles/](articles/) — Mathlog記事とREADMEリンクを突き合わせるツール。詳細は[articles/README.md](articles/README.md)参照。
- [common/](common/) — 各記事の検証コードから共有されるPythonモジュール（`octonion.py`など）。
- [reftools/](reftools/) — 参考文献slugの集約・整合性チェックツール。詳細は[reftools/README.md](reftools/README.md)参照。
- [bookmarklets/](bookmarklets/) — Mathlog関連のブックマークレット集。
- [mathlog_ref.sh](mathlog_ref.sh) — `mathlog.tsv`の未取得記事について、参考文献パネルのHTMLをクリップボード経由で取り込み、整形・TOML化する。
- [mathlog_fix.sh](mathlog_fix.sh) — `mathlog_fix.md`に列挙した既存記事を開き直し、本文修正と参考文献の再取り込みを行う。
