# 検証コード

プロジェクトルートで `uv run hopf/check/<ファイル名>` により実行します。

- [01-quaternion.py](01-quaternion.py) — $ω\mathbf kω^*$の成分展開とノルム、複素数ペア・パウリ行列による表式の一致、ファイバー$q=u+v\mathbf k$の検証。
- [02-spinor-tensor.py](02-spinor-tensor.py) — $ωq$による$\mathbf k$項の消去、三角関数表示での回転結果、ブロッホベクトルと密度行列の検証。01-quaternion.md の$ω$との対応（右から$-\mathbf k$を掛けたもの）も確認します。
- [03-bloch-density.py](03-bloch-density.py) — パウリ行列の積とエルミート行列の展開、純粋状態の性質（べき等性・位相不変性）、混合状態のブロッホベクトルと純粋度、固有分解、混合と重ね合わせの例の検証。
- [04-extension.py](04-extension.py) — $H(\alpha, \beta)=(2\alpha\beta^*, |\alpha|^2-|\beta|^2)$ の四元数・八元数への拡張によるノルム1の保存の検証、および右乗算によるファイバーが四元数では結合性により保存され、八元数では非結合性により破壊されることの確認。
- [04-nonassoc-fiber.py](04-nonassoc-fiber.py) — 八元数で右乗算がファイバーを保存しない具体例（$\alpha=e_1/\sqrt2$等）と、左乗算$\beta=-p\alpha$によるファイバーの構成、一般の点$(c, r)$への拡張の検証。
- [05-two-qubit.py](05-two-qubit.py) — $\mathbb{C}^4\cong\mathbb{H}^2$の詰め込みと$q_1q_2^*$の成分公式、分離可能⟺$\alpha\delta-\beta\gamma=0$⟺$\mathbf j,\mathbf k$成分の消失、関係式$x_A^2+y_A^2+z_A^2+C^2=1$、ベル状態の像と全体位相・ファイバーの検証。
- [05-three-qubit.py](05-three-qubit.py) — $\mathbb{C}^4\cong\mathbb{O}$の複素基底$\{1,e_2,e_4,e_6\}$による詰め込みと左$\mathbb{C}$-線形性、A|BC分離可能⟺$e_2,\dots,e_7$成分の消失、GHZ・W状態の像、ファイバー$S^7$の分離可能性、$e_1$を固定する自己同型群$\operatorname{SU}(3)$の検証。
- [06-single-qubit.py](06-single-qubit.py) — 回転子$R_{\boldsymbol n}(\theta)$の共役作用とロドリゲスの回転公式の一致、6ゲートの回転子分解、パウリ群と四元数群$Q_8$の対応、$\langle H,S\rangle$が位数24の正八面体の回転群を生成すること、6状態の軌道、$T$が非クリフォードであることの検証。
- [06-cnot.py](06-cnot.py) — CNOTの共役作用表と15軸の符号付き置換、局所ゲートでの$C$の不変性ともつれの生成、ベル回路$\mathrm{CNOT}(H\otimes\mathrm I)|00\rangle=\Phi_+$とホップ像$-\mathbf j$、軸の追跡$Z\otimes\mathrm I\mapsto X\otimes X$、密度行列の15成分展開の検証。
