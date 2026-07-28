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
- [07-real-qubit.py](07-real-qubit.py) — 実振幅のブロッホベクトル$(2\alpha\beta,0,\alpha^2-\beta^2)$と実数ホップ写像の一致、$\zeta^2$による二重被覆、密度行列の第1行が描く半径$\frac12$の円、実振幅のボルンの規則、$O(2)$が底空間の$O(2)$へ2対1に落ちること、実パウリ群$\cong D_4$と$Q_8$の区別、実クリフォード群16元と正八面体群の指数3の部分群、$S,T$が$y=0$平面を保たないことの検証。
- [07-two-rebit.py](07-two-rebit.py) — 主要ゲートとベル状態が実であること、分離可能性$\alpha\delta-\beta\gamma=0$、実観測量のみでCHSHがツィレルソン限界$2\sqrt2$に到達すること、補助レビットによる複素振幅の符号化$I\otimes A+J\otimes B$、局所観測量が9次元しか張らず$\sigma_y\otimes\sigma_y$が残り1次元を埋めること、$\rho_\pm$が局所測定で区別できないことの検証。
- [memo-c2.py](memo-c2.py) — MEMO.md の主線（$\mathbb{C}^2$上の線形ベクトル場としての実現）の検証。線形作用$v\mapsto e^{tX}v$の生成子$D_X=\sum X_{ik}v_k\partial_{v_i}$が指数関数から再導出されることと反準同型（$[D_X,D_Y]=-D_{[X,Y]}$）、転置で準同型に揃うこと、$J_+=\alpha\partial_\beta,\ J_-=\beta\partial_\alpha,\ J_3=(\alpha\partial_\alpha-\beta\partial_\beta)/2$の交換関係、ファイバーの生成子がオイラー作用素$N=\alpha\partial_\alpha+\beta\partial_\beta$であり$N$の固有値$2j$が同次式の次数に一致すること、次数$2j$の同次多項式が$2j+1$次元でカシミール$j(j+1)$を持つこと（$j\le2$）、$j=1/2$で$\{\alpha,\beta\}$の行列が$\sigma_3/2,\sigma_\pm$に一致すること、$F=\alpha^{2j}g(\beta/\alpha)$への分解が比の表示$\partial_w,\ j-w\partial_w,\ 2jw-w^2\partial_w$を再現することの検証。
- [memo-mobius.py](memo-mobius.py) — MEMO.md の「比で書いた場合」の検証。$w=\beta/\alpha=\tan(\theta/2)e^{i\varphi}=(n_1+in_2)/(1+n_3)$が南極からの立体射影であること、生成子$X$による$w$の動き$c+(d-a)w-bw^2$（メビウス変換と指数関数からの再導出）、$-i\sigma_a/2$に対応する係数が$d\boldsymbol n/dt=\boldsymbol e_a\times\boldsymbol n$を再現すること、昇降演算子の像$\partial_w,\ -w^2\partial_w$がブロッホ球の回転としては存在しないことの検証。
- [memo-spin-j.py](memo-spin-j.py) — MEMO.md の比の表示とハイゼンベルク代数への縮約の検証。$J_+=\partial_w,\ J_3=j-w\partial_w,\ J_-=2jw-w^2\partial_w$の交換関係・カシミール$j(j+1)$・ウェイト基底と梯子の端（$j\le2$）、$j=1/2$で$\{1,w\}\leftrightarrow\{\alpha,\beta\}$の行列が$\sigma_3/2,\sigma_\pm$に一致すること、$u=\sqrt{2j}\,w$による$j\to\infty$の縮約で$[a,a^\dagger]=1$（掛け算と微分）に至ることの検証。
