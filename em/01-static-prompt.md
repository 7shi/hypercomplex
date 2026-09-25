---
refs:
  - clif-analysis/05-clifford.md
  - clif-analysis/02-fundamental-theorem.md
---

添付は「幾何代数による電磁気学」シリーズ第1回の記事です。「クリフォード解析」シリーズ（全6回）の物理側への展開として、$\operatorname{Cl}_{3,0}(\mathbb R)$のディラック作用素$D$、コーシー核と基本解、2つの基本定理を電磁場に使い、クーロンの法則とビオ＝サバールの法則を1本の積分公式$F=\frac1{4\pi\varepsilon_0}\int\frac{\boldsymbol x-\boldsymbol y}{|\boldsymbol x-\boldsymbol y|^3}(\rho-\boldsymbol J/c)\,dV$にまとめます。参照記事として、核・$|S^{n-1}|$・基本解・$D|\boldsymbol x|^{2-n}$の出所である第5回と、領域と曲面の基本定理・有向面素の出所である第2回を添付します。以下の観点でレビューしてください。

- 数式・数学的主張の正しさ
    - `&&&prop スカラー部の消滅`の証明（有界な台の仮定、$\boldsymbol y=\boldsymbol x$の特異点を小球で除く極限の扱い）が十分か
    - `&&&thm 静的なマクスウェル方程式`の証明で、「$D$は核にだけ作用し、源は核の右に残るので非可換性は問題にならない」とする説明と、基本解の性質（第5回の`&&&rem 基本解`）への依拠の仕方が正確か。超関数に立ち入らないという線引きが明確か
    - `&&&rem スカラー部が残る場合`（$\nabla\cdot\boldsymbol J\ne0$で$\Delta\Phi=-\nabla\cdot\boldsymbol J/\varepsilon_0c$）の議論が妥当か
- 説明の分かりやすさ、論理の飛躍の有無
    - 実験則から積分公式、微分形、積分形、ポテンシャルへ進む節の順序が自然か。「磁場は2ベクトル」の節と`&&&rem 擬ベクトルの解消`（空間反転での符号）の説明が、擬ベクトルの説明として不正確になっていないか
- 用語・記法の一貫性
    - 第5回との照合：核$E=\boldsymbol x/|\boldsymbol x|^n$を電場との衝突を避けて分数のまま書いた扱い、$4\pi=|S^2|$、基本解$\boldsymbol x/(|S^{n-1}||\boldsymbol x|^n)$、`&&&fml 核とポテンシャル`の$n=3$の場合（$D|\boldsymbol x|^{-1}=-\boldsymbol x/|\boldsymbol x|^3$）の引き方
    - 第2回との照合：領域の基本定理$\int_VDF\,dV=\oint\boldsymbol nF\,dS$と曲面の基本定理$\oint d\boldsymbol x\,F=\int_S(D\cdot d\boldsymbol X)F$、$d\boldsymbol X=I\boldsymbol n\,dA$、$(D\cdot d\boldsymbol X)F=((\boldsymbol n\times\nabla)F)dA$、「回転の体積分」の引用が規約（向き・内積の順序）と一致しているか。積分形の4成分（とくに$\oint\boldsymbol n\times\boldsymbol B\,dS=\mu_0\int\boldsymbol J\,dV$）の符号
    - SI単位系の係数（$\varepsilon_0\mu_0c^2=1$、$F=\boldsymbol E+Ic\boldsymbol B$の$c$）が全体で一貫しているか
- 日本語表現の自然さ
- プロジェクト方針
    - 後続記事の内容を予告する前方参照がないか（ゲージの自由度を「本記事では扱いません」と先送りした箇所の書き方を含む）
    - `&&&`ブロックの種類の使い分け（def・thm・prop・fml・prf・rem）が内容に合っているか。本文と定理の主張はですます調、`&&&prf`内だけである調になっているか

次の計算は別途、記号計算と数値積分で検証済みです：被積分関数のグレード分けとビオ＝サバールの係数、スカラー部についての恒等式、$-D(\varphi-c\boldsymbol A)$の成分、ガウス型の源での$F$のスカラー部の消滅・$DF=(\rho-\boldsymbol J/c)/\varepsilon_0$・4本の式・源の外の点での積分公式、$\boldsymbol nF$・$d\boldsymbol x\,F$・$(\boldsymbol n\times\nabla)F$のグレード成分。検算よりも、論証の構成と述べ方の正確さを優先してください。

指摘事項を箇条書きで挙げてください。問題がなければその旨を書いてください。
