# 検証コード

プロジェクトルートで `uv sync` により環境を構築し、`uv run qua/check/<ファイル名>` で実行します。

- [04-4d-bsqua.py](04-4d-bsqua.py) — 04-4d-bsqua.md の数式の数値検証。$\omega$による2ベクトルの対と冪等元$\dfrac{1\pm\omega}2$による$\mathbb H\oplus\mathbb H$分解、対角表現の成分ごとの積と$T$の準同型性、等傾回転の合成による二重回転（$e^{\alpha i}qe^{\beta i}$の2平面の回転角$\alpha\pm\beta$）、正規形・一般の回転面での$r_L\,q\,r_R$とクリフォード代数の回転の一致、$L_uR_v$が張る$M_4(\mathbb R)\cong\operatorname{Cl}_{3,1}(\mathbb R)$、テンソル積での符号の相殺（2対1対応の核）、偶部分代数の階梯を確認します。
- [05-dual-qua.py](05-dual-qua.py) — 05-dual-qua.md の数式の数値検証。二重数の演算と$f(x+\varepsilon y)=f(x)+\varepsilon y f'(x)$、3種の共役と二重数値ノルム、剛体変換のサンドイッチ公式と$\mathrm{SE}(3)$の合成則（通常の共役では並進が消える負のチェックを含む）、接平面モデル（$\varepsilon$倍角の回転子による並進と基点を通る軸の回転）、ねじ運動の指数表示とねじ軸の抽出、$C\ell^+_{3,0,1}(\mathbb R)$との同型（全$64$積）、擬スカラーの三分法を確認します。
- [06-slerp.py](06-slerp.py) — 06-slerp.md の数式の数値検証。対数写像と冪の指数法則、なす角と半角関係、大円弧の等速性と最短性、冪形式とsin重み形式の一致（外挿を含む）、左右不変性と対称性、世界座標軸と角速度の一定性、符号選択と2本の経路（中間の姿勢が異なる負のチェック）、NLERPの不等速性、$C\ell^+_{3,0}(\mathbb R)$・$C\ell_{4,0}(\mathbb R)$での回転子の補間（線形結合が回転子にならない負のチェック）、ScLERPの端点・実部・軸上の並進、退化の数値安定性を確認します。
- [spinor-ideal-lxryp.py](spinor-ideal-lxryp.py)
