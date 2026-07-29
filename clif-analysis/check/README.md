# 検証コード

クリフォード解析の記事の数式を検証するPythonスクリプトです。プロジェクトルートで `uv sync` により環境を構築し、`uv run clif-analysis/check/<ファイル名>` で実行します。

- [memo-geometric-calculus.py](memo-geometric-calculus.py) — [MEMO.md](../MEMO.md)の第1回（$C\ell_{2,0}(\mathbb R)$のディラック作用素とモノジェニック関数）の数式の検証。$C\ell_{2,0}(\mathbb R)$を$M_2(\mathbb R)$で実現し、ディラック作用素$D=e_1\partial_x+e_2\partial_y$が$D^2=\Delta$を満たすこと、$e_1$を掛けるとヴィルティンガー微分になること、偶部分に値を取る関数のモノジェニック性がコーシー＝リーマンの方程式と同値であること、奇部分では反正則性に入れ替わること、成分が調和であること、$\boldsymbol x^{-1}$が左右ともモノジェニックで$e_1\boldsymbol x^{-1}=1/\bar z$となること、ベクトル値関数に対して$DF$が発散と回転に分解することを確認します。
