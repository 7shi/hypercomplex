# 検証コード

プロジェクトルートで `uv run misc/check/<ファイル名>` により実行します。

- [energy-quantize-zeta.py](energy-quantize-zeta.py) — energy-quantize-zeta.md の数式の検証。プランクの公式の低周波・高周波極限、モード数とモード密度$8\pi\nu^2/c^3$、連続分布の平均エネルギー$kT$、量子化した平均エネルギー$h\nu/(e^{h\nu/kT}-1)$（記号計算と級数の数値）、$\int_0^\infty x^3e^{-nx}dx=6/n^4$と$\int_0^\infty x^3/(e^x-1)dx=6\zeta(4)$、全エネルギー密度とシュテファン＝ボルツマン定数、$x^2$のフーリエ係数（$a_0/2$規約）とパーセヴァルの等式による$\zeta(4)=\pi^4/90$を確認します。
- [variable-dependence.py](variable-dependence.py) — variable-dependence.md の数式の検証。3つの具体例（$xy,\ y=x^2$／$\sin(xy),\ y=e^x$／$x^2+y^2,\ y=\sin x$）について、先に代入して微分した結果と、偏微分に$dy/dx$を掛けて合成した結果が記事の値に一致することと、オイラー＝ラグランジュ方程式の$\frac{d}{dt}\frac{\partial L}{\partial\dot q}$の展開を確認します。
- [integration-by-parts.py](integration-by-parts.py) — integration-by-parts.md の数式の検証。$\int x\cos x\,dx$ と $\int x^2\sin x\,dx$ の結果を微分して被積分関数に戻ること、途中式 $\int x^2(-\sin x)\,dx=x^2\cos x-\int2x\cos x\,dx$、次数が上がる変形の等式、定積分の部分積分の具体例を確認します。
