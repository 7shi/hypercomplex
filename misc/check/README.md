# 検証コード

プロジェクトルートで `uv run misc/check/<ファイル名>` により実行します。

- [energy-quantize-zeta.py](energy-quantize-zeta.py) — energy-quantize-zeta.md の数式の検証。プランクの公式の低周波・高周波極限、モード数とモード密度$8\pi\nu^2/c^3$、連続分布の平均エネルギー$kT$、量子化した平均エネルギー$h\nu/(e^{h\nu/kT}-1)$（記号計算と級数の数値）、$\int_0^\infty x^3e^{-nx}dx=6/n^4$と$\int_0^\infty x^3/(e^x-1)dx=6\zeta(4)$、全エネルギー密度とシュテファン＝ボルツマン定数、$x^2$のフーリエ係数（$a_0/2$規約）とパーセヴァルの等式による$\zeta(4)=\pi^4/90$を確認します。
