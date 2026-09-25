静電場を与えるクーロンの法則と、静磁場を与えるビオ＝サバールの法則は、電磁気学の出発点となる2つの実験則です。ベクトル解析では、前者は電荷密度と変位ベクトルの積、後者は電流密度と変位ベクトルのベクトル積を使う、別々の式として書かれます。幾何代数では、両者は電荷密度と電流密度をまとめた1つの量に同じ核を掛けて積分する1本の公式にまとまり、その結果のベクトル部が電場、2ベクトル部が磁場になります。核は[[7shi-cla5]]のコーシー核$\boldsymbol x/|\boldsymbol x|^n$の$n=3$の場合であり、クーロン定数$1/4\pi\varepsilon_0$の$4\pi$は単位球面$S^2$の面積です。この公式にディラック作用素$D$を掛けると、静的な場のマクスウェル方程式4本が1本の式のグレード成分として得られます。磁場がもともと2ベクトルであることから、ベクトル積が擬ベクトルとして例外的に扱われる事情（[[7shi-hist]]）も解消されます。

# 設定

$\operatorname{Cl}_{3,0}(\mathbb R)$の生成元を$e_1,e_2,e_3$（$e_k^2=1$、$k\ne l$なら$e_ke_l=-e_le_k$）とし、擬スカラーを$I=e_1e_2e_3$とします。$I^2=-1$で、$I$はすべての元と可換です。空間の点を$\boldsymbol x=\sum_kx_ke_k$と書き、ディラック作用素を

$$
D=\sum_{k=1}^3e_k\partial_k
$$

とします。$D^2=\Delta$です（[[7shi-cla1]]）。本記事の量はすべて時間によらない静的なものとします。

ベクトル$\boldsymbol a,\boldsymbol b$の幾何積は$\boldsymbol a\boldsymbol b=\boldsymbol a\cdot\boldsymbol b+\boldsymbol a\wedge\boldsymbol b$と分かれ、3次元では外積とベクトル積が擬スカラーで結ばれます。

$$
\boldsymbol a\wedge\boldsymbol b=I(\boldsymbol a\times\boldsymbol b)
$$

たとえば$e_1\wedge e_2=e_1e_2=Ie_3$です。ベクトル値の関数$\boldsymbol V$に$D$を作用させると、同じ分解から発散と回転が現れます。$I$は$D$とも可換なので、2ベクトル値の関数$I\boldsymbol V$についても続けて書けます。

&&&fml ベクトルと2ベクトルの微分 [fml-dv]
$$
D\boldsymbol V=\nabla\cdot\boldsymbol V+I\,\nabla\times\boldsymbol V,\qquad
D(I\boldsymbol V)=I\,\nabla\cdot\boldsymbol V-\nabla\times\boldsymbol V
$$
&&&

$D$はベクトルを掛ける作用素なので、グレードを1つ上げる部分と1つ下げる部分を持ちます。ベクトルからはスカラーと2ベクトルが、2ベクトルからはベクトルと擬スカラーが出ます。

物理定数は$\mathrm{SI}$単位系で書き、真空の誘電率$\varepsilon_0$、透磁率$\mu_0$、光速$c$の間には$\varepsilon_0\mu_0c^2=1$が成り立ちます。

# 2つの実験則

電荷密度$\rho$と電流密度$\boldsymbol J$は$C^1$級で、ある有界な領域の外で$0$になるとします。積分は$\boldsymbol y$について全空間で取り、$dV$は$\boldsymbol y$の体積要素です。

&&&fml クーロンの法則とビオ＝サバールの法則
$$
\boldsymbol E(\boldsymbol x)=\frac1{4\pi\varepsilon_0}\int\rho(\boldsymbol y)\frac{\boldsymbol x-\boldsymbol y}{|\boldsymbol x-\boldsymbol y|^3}\,dV,\qquad
\boldsymbol B(\boldsymbol x)=\frac{\mu_0}{4\pi}\int\boldsymbol J(\boldsymbol y)\times\frac{\boldsymbol x-\boldsymbol y}{|\boldsymbol x-\boldsymbol y|^3}\,dV
$$
&&&

どちらの式にも同じベクトル値の核$(\boldsymbol x-\boldsymbol y)/|\boldsymbol x-\boldsymbol y|^3$が現れます。違うのは、核に掛ける源がスカラーの$\rho$かベクトルの$\boldsymbol J$か、そして掛け方が数の積かベクトル積かという点です。ビオ＝サバールの法則は定常電流、すなわち$\nabla\cdot\boldsymbol J=0$を満たす電流に対して使います。電荷が保存されるので、$\nabla\cdot\boldsymbol J\ne0$なら電荷密度が時間とともに変わり、静的という前提が崩れるためです。

# 磁場は2ベクトル

ベクトル積を外積に戻します。$\boldsymbol r=\boldsymbol x-\boldsymbol y$と書くと、$\boldsymbol J\times\boldsymbol r=-I(\boldsymbol J\wedge\boldsymbol r)=I(\boldsymbol r\wedge\boldsymbol J)$です。ビオ＝サバールの法則の両辺に$Ic$を掛け、$I^2=-1$と$\mu_0c=1/\varepsilon_0c$を使うと

$$
Ic\boldsymbol B(\boldsymbol x)=-\frac1{4\pi\varepsilon_0c}\int\frac{\boldsymbol r\wedge\boldsymbol J(\boldsymbol y)}{|\boldsymbol r|^3}\,dV
$$

となります。右辺にはベクトル積も擬スカラーも現れず、2つのベクトル$\boldsymbol r$と$\boldsymbol J$が張る向き付きの面の重ね合わせだけが残ります。磁場は、源点$\boldsymbol y$での電流と、源点から観測点$\boldsymbol x$への変位が張る面として、2ベクトルの量で書けます。

&&&rem 擬ベクトルの解消
[[7shi-hist]]で見たとおり、ギブスとヘヴィサイドのベクトル解析は四元数の積を内積とベクトル積に分け、ベクトル積は鏡映に対してベクトルと異なる符号の変化を示す擬ベクトルとして扱われました。磁場$\boldsymbol B$はその代表例です。空間反転$\boldsymbol x\mapsto-\boldsymbol x$で$\boldsymbol r$と$\boldsymbol J$はともに符号を変えるので、$\boldsymbol r\wedge\boldsymbol J$は符号を変えません。これは2ベクトルとして当然の変換で、例外的な扱いは要りません。向きを固定した擬スカラー$I$を使って2ベクトル$Ic\boldsymbol B$をベクトル$\boldsymbol B$で表示すると、その表示に擬ベクトルの変換則が現れます。擬ベクトルという扱いは、2ベクトルの変換則をベクトル表示で読んだものです。この表示が可能なのは、3次元で2ベクトルとベクトルがともに3成分を持つためです。
&&&

係数の$c$は、$\boldsymbol E$と$c\boldsymbol B$の単位を揃えるために入れています。

# 積分公式

電場の式は$\rho$に核を掛け、磁場の式は核と$\boldsymbol J$の外積を取っています。幾何積は内積と外積をまとめて含むので、源を$\rho-\boldsymbol J/c$というスカラーとベクトルの和（パラベクトル）にまとめて核に掛ければ、2つの式を同時に扱えます。

&&&def 静的な場の積分公式 [def-F]
$$
F(\boldsymbol x)=\frac1{4\pi\varepsilon_0}\int\frac{\boldsymbol x-\boldsymbol y}{|\boldsymbol x-\boldsymbol y|^3}\Bigl(\rho(\boldsymbol y)-\frac{\boldsymbol J(\boldsymbol y)}c\Bigr)dV
$$
&&&

被積分関数を$\boldsymbol r=\boldsymbol x-\boldsymbol y$でグレードに分けます。

$$
\boldsymbol r\Bigl(\rho-\frac{\boldsymbol J}c\Bigr)=\rho\,\boldsymbol r-\frac1c\,\boldsymbol r\cdot\boldsymbol J-\frac1c\,\boldsymbol r\wedge\boldsymbol J
$$

ベクトル部$\rho\boldsymbol r$はクーロンの法則の被積分関数、2ベクトル部$-\boldsymbol r\wedge\boldsymbol J/c$は前節の$Ic\boldsymbol B$の被積分関数です。残るのはスカラー部$-\boldsymbol r\cdot\boldsymbol J/c$で、これが積分すると消えることを確かめます。

&&&prop スカラー部の消滅 [prop-scalar]
$\nabla\cdot\boldsymbol J=0$なら、すべての$\boldsymbol x$で

$$
\int\frac{(\boldsymbol x-\boldsymbol y)\cdot\boldsymbol J(\boldsymbol y)}{|\boldsymbol x-\boldsymbol y|^3}\,dV=0
$$

が成り立ちます。
&&&

&&&prf
$\boldsymbol y$についての勾配を$\nabla_{\boldsymbol y}$と書くと、$\nabla_{\boldsymbol y}|\boldsymbol x-\boldsymbol y|^{-1}=(\boldsymbol x-\boldsymbol y)/|\boldsymbol x-\boldsymbol y|^3$である。積の微分から

$$
\frac{(\boldsymbol x-\boldsymbol y)\cdot\boldsymbol J}{|\boldsymbol x-\boldsymbol y|^3}
=\nabla_{\boldsymbol y}\cdot\frac{\boldsymbol J}{|\boldsymbol x-\boldsymbol y|}-\frac{\nabla\cdot\boldsymbol J}{|\boldsymbol x-\boldsymbol y|}
$$

となり、仮定から右辺の第2項は$0$である。第1項を、$\boldsymbol J$の台と点$\boldsymbol x$を内部に含む十分大きな球から$\boldsymbol x$を中心とする半径$\varepsilon$の小球を除いた領域で積分し、発散定理（[[7shi-cla2]]）を使う。外側の球面上では$\boldsymbol J=0$である。小球面上では$|\boldsymbol J/|\boldsymbol x-\boldsymbol y||\le\max|\boldsymbol J|/\varepsilon$で面積は$4\pi\varepsilon^2$なので、境界積分は$\varepsilon\to0$で$0$に近づく。左辺の被積分関数は大きさが$|\boldsymbol x-\boldsymbol y|^{-2}$程度で積分可能なので、小球を除いた積分は$\varepsilon\to0$で全空間の積分に収束する。
&&&

$\nabla\cdot\boldsymbol J=0$は、ビオ＝サバールの法則を使う前提そのものでした。したがって定常電流では$F$にスカラー部が残らず、$F$はベクトル部と2ベクトル部だけを持ちます。

&&&thm 静的な場の統一
$\nabla\cdot\boldsymbol J=0$のとき、[積分公式](#def-F)の$F$は

$$
F=\boldsymbol E+Ic\boldsymbol B
$$

で、ベクトル部がクーロンの法則の電場、2ベクトル部がビオ＝サバールの法則の磁場の$Ic$倍です。
&&&

核は[[7shi-cla5]]のコーシー核$E(\boldsymbol x)=\boldsymbol x/|\boldsymbol x|^n$の$n=3$の場合です。記号$E$は電場と紛らわしいので、本記事では核を分数のまま書きます。クーロン定数の$4\pi$は、[[7shi-cla5]]の正規化定数$|S^{n-1}|$の$n=3$の値$|S^2|=4\pi$です。電磁気学では$4\pi$を「全立体角」と説明しますが、それは単位球面の面積のことであり、積分公式の定数として自然に現れます。

# 静的なマクスウェル方程式

[[7shi-cla5]]の注意「基本解」で見たとおり、$\boldsymbol x/(|S^{n-1}||\boldsymbol x|^n)$は$D$の基本解です。$n=3$では

$$
D\,\frac{\boldsymbol x}{4\pi|\boldsymbol x|^3}=\delta(\boldsymbol x)
$$

を表します。[積分公式](#def-F)では、$D$は$\boldsymbol x$についての微分なので核にだけ作用し、源$\rho-\boldsymbol J/c$は$\boldsymbol y$の関数として核の右に残ります。したがって$F$は、源の$1/\varepsilon_0$倍を$D$の基本解で積分したものです。

&&&thm 静的なマクスウェル方程式 [thm-static]
$$
DF=\frac1{\varepsilon_0}\Bigl(\rho-\frac{\boldsymbol J}c\Bigr)
$$
&&&

&&&prf
基本解の性質により、台がコンパクトな$C^1$級の関数$G$に対して$D\int\frac{\boldsymbol x-\boldsymbol y}{4\pi|\boldsymbol x-\boldsymbol y|^3}G(\boldsymbol y)\,dV=G(\boldsymbol x)$である。ここで$G$は左から核が掛かる形で積分されており、$D$は核を通り越さずに核へ直接作用するので、非可換性は問題にならない。$G=(\rho-\boldsymbol J/c)/\varepsilon_0$とすればよい。核の微分は特異点の外では$0$であり、$G$を再現する寄与は特異点から来る。したがってこの等式は、積分記号の下で通常の微分を取って得られるものではなく、基本解の性質として既知の結果を引用したものである。$F$を最後の節のポテンシャルの$D$微分として書き、$D^2=\Delta$とラプラシアンの基本解に帰着させても同じ結論を得る。本記事では超関数の枠組みには立ち入らない。
&&&

左辺を[ベクトルと2ベクトルの微分](#fml-dv)で成分に分けると

$$
D(\boldsymbol E+Ic\boldsymbol B)=\nabla\cdot\boldsymbol E-c\,\nabla\times\boldsymbol B+I\,\nabla\times\boldsymbol E+Ic\,\nabla\cdot\boldsymbol B
$$

です。右辺はスカラーとベクトルだけなので、グレードごとに比べて4本の式が得られます。ベクトル部では$-c\,\nabla\times\boldsymbol B=-\boldsymbol J/\varepsilon_0c$の両辺を$-c$で割り、$1/\varepsilon_0c^2=\mu_0$を使います。

&&&fml 静的な場の4本の式
| グレード | 式 | 名前 |
|---|---|---|
| スカラー | $\nabla\cdot\boldsymbol E=\rho/\varepsilon_0$ | ガウスの法則 |
| ベクトル | $\nabla\times\boldsymbol B=\mu_0\boldsymbol J$ | アンペールの法則 |
| 2ベクトル | $\nabla\times\boldsymbol E=0$ | 静電場の渦なし |
| 擬スカラー | $\nabla\cdot\boldsymbol B=0$ | 磁気単極子の不在 |
&&&

電場の発散と回転、磁場の発散と回転は、ベクトル解析では別々の4本の式ですが、ここでは$D$を1回作用させた結果の4つのグレードです。源を持つのはスカラー部とベクトル部だけで、2ベクトル部と擬スカラー部が$0$になるのは源$\rho-\boldsymbol J/c$に2ベクトルと擬スカラーの成分がないためです。磁気単極子があればその密度は擬スカラーの位置に入ります。

&&&rem スカラー部が残る場合
この注意では、二階微分を通常の意味で扱うため$\boldsymbol J$を$C^2$級とします。$\nabla\cdot\boldsymbol J\ne0$の電流で[積分公式](#def-F)を作ると、$F$にスカラー部$\Phi$が残ります。[静的なマクスウェル方程式](#thm-static)の証明はそのまま成り立つので$DF=(\rho-\boldsymbol J/c)/\varepsilon_0$は変わりませんが、ベクトル部が$\nabla\Phi-c\,\nabla\times\boldsymbol B=-\boldsymbol J/\varepsilon_0c$となり、アンペールの法則が$\nabla\Phi$の分だけずれます。両辺の発散を取ると$\nabla\cdot(\nabla\times\boldsymbol B)=0$から$\Delta\Phi=-\nabla\cdot\boldsymbol J/\varepsilon_0c$です。静的な枠組みでは、この項は定常電流の条件で消すしかありません。
&&&

# 積分形

ガウスの法則とアンペールの法則には積分形があります。[[7shi-cla2]]の2つの基本定理

$$
\int_VDF\,dV=\oint_{\partial V}\boldsymbol nF\,dS,\qquad
\oint_{\partial S}d\boldsymbol x\,F=\int_S(D\cdot d\boldsymbol X)F
$$

に$F=\boldsymbol E+Ic\boldsymbol B$を入れ、グレード成分を取ると得られます。領域・曲面・境界の条件と向きは[[7shi-cla2]]と同じとします。

## 領域の基本定理

外向きの単位法線$\boldsymbol n$との積を成分に分けます。

$$
\boldsymbol nF=\boldsymbol n\cdot\boldsymbol E-c\,\boldsymbol n\times\boldsymbol B+I\,\boldsymbol n\times\boldsymbol E+Ic\,\boldsymbol n\cdot\boldsymbol B
$$

左辺には[静的なマクスウェル方程式](#thm-static)を代入します。$Q=\int_V\rho\,dV$は$V$の中の電荷です。

&&&fml 領域の基本定理のグレード成分
$$
\oint_{\partial V}\boldsymbol E\cdot\boldsymbol n\,dS=\frac Q{\varepsilon_0},\qquad
\oint_{\partial V}\boldsymbol B\cdot\boldsymbol n\,dS=0,
$$

$$
\oint_{\partial V}\boldsymbol n\times\boldsymbol E\,dS=0,\qquad
\oint_{\partial V}\boldsymbol n\times\boldsymbol B\,dS=\mu_0\int_V\boldsymbol J\,dV
$$
&&&

スカラー部はガウスの法則の積分形、擬スカラー部は閉曲面を通る正味の磁束が$0$であることです。残る2つは、[[7shi-cla2]]の「回転の体積分」$\int_V\nabla\times\boldsymbol V\,dV=\oint\boldsymbol n\times\boldsymbol V\,dS$を$\nabla\times\boldsymbol E=0$と$\nabla\times\boldsymbol B=\mu_0\boldsymbol J$に当てはめたものにあたります。

## 曲面の基本定理

アンペールの法則の積分形は、閉曲線$\partial S$に沿った線積分です。こちらは曲面の基本定理から出ます。[[7shi-cla2]]で見たとおり、有向面素を$d\boldsymbol X=I\boldsymbol n\,dA$と書くと$(D\cdot d\boldsymbol X)F=\bigl((\boldsymbol n\times\nabla)F\bigr)dA$です。左辺の$d\boldsymbol x\,F$のスカラー部は$d\boldsymbol x\cdot\boldsymbol E$、擬スカラー部は$Ic\,d\boldsymbol x\cdot\boldsymbol B$であり、右辺の対応する成分は$\boldsymbol n\cdot(\nabla\times\boldsymbol E)$と$Ic\,\boldsymbol n\cdot(\nabla\times\boldsymbol B)$です。

&&&fml 曲面の基本定理のグレード成分
$$
\oint_{\partial S}\boldsymbol E\cdot d\boldsymbol x=0,\qquad
\oint_{\partial S}\boldsymbol B\cdot d\boldsymbol x=\mu_0\int_S\boldsymbol J\cdot\boldsymbol n\,dA
$$
&&&

スカラー部は静電場の周回積分が消えること、擬スカラー部はアンペールの法則の積分形です。

微分形の4本の式が$DF$のグレード成分だったのと同じく、積分形の式は2つの基本定理のグレード成分として並びます。ガウスの法則が体積と閉曲面の関係、アンペールの法則が曲面と閉曲線の関係になるのは、$\boldsymbol E$の源（スカラー）と$\boldsymbol B$の源（ベクトル）のグレードの違いに対応しています。

# ポテンシャル

[[7shi-cla5]]の公式「核とポテンシャル」で、核はラプラシアンの基本解の$D$微分として得られました。$n=3$では

$$
D\frac1{|\boldsymbol x|}=-\frac{\boldsymbol x}{|\boldsymbol x|^3}
$$

です。$D$は$\boldsymbol x$についての微分なので、$\boldsymbol x$を$\boldsymbol x-\boldsymbol y$に替えても同じ式が成り立ちます。これを[積分公式](#def-F)に使うと、$F$を1つの関数の$D$微分として書けます。

&&&def ポテンシャル
$$
P(\boldsymbol x)=\frac1{4\pi\varepsilon_0}\int\frac{\rho(\boldsymbol y)-\boldsymbol J(\boldsymbol y)/c}{|\boldsymbol x-\boldsymbol y|}\,dV
=\varphi-c\boldsymbol A
$$

ここで

$$
\varphi=\frac1{4\pi\varepsilon_0}\int\frac{\rho(\boldsymbol y)}{|\boldsymbol x-\boldsymbol y|}\,dV,\qquad
\boldsymbol A=\frac{\mu_0}{4\pi}\int\frac{\boldsymbol J(\boldsymbol y)}{|\boldsymbol x-\boldsymbol y|}\,dV
$$

はスカラーポテンシャルとベクトルポテンシャルです（$1/4\pi\varepsilon_0c=\mu_0c/4\pi$）。
&&&

核$1/|\boldsymbol x-\boldsymbol y|$はスカラーなので、$D$を積分の中に入れても源との順序は問題になりません。

&&&fml 場とポテンシャル
$$
F=-DP=-\nabla\varphi+c\,\nabla\cdot\boldsymbol A+Ic\,\nabla\times\boldsymbol A
$$
&&&

グレードを比べると$\boldsymbol E=-\nabla\varphi$、$\boldsymbol B=\nabla\times\boldsymbol A$という通常の関係が得られます。スカラー部は$c\,\nabla\cdot\boldsymbol A$なので、[スカラー部の消滅](#prop-scalar)は$\nabla\cdot\boldsymbol A=0$と同じことです。ベクトルポテンシャルを上の積分で定めると、定常電流に対しては発散が自動的に$0$になります。

$D^2=\Delta$なので、[静的なマクスウェル方程式](#thm-static)は$\Delta P=-(\rho-\boldsymbol J/c)/\varepsilon_0$、すなわち成分ごとのポアソン方程式

$$
\Delta\varphi=-\frac\rho{\varepsilon_0},\qquad\Delta\boldsymbol A=-\mu_0\boldsymbol J
$$

と同じ内容になります。上の積分表示による$\boldsymbol A$は、$\nabla\cdot\boldsymbol A=0$（クーロンゲージ）を満たします。$\boldsymbol A$に勾配$\nabla\chi$を加えても$\boldsymbol B$は変わりませんが、$-D\bigl(\varphi-c(\boldsymbol A+\nabla\chi)\bigr)=F+c\Delta\chi$となり、$F=-DP$のスカラー部は変わります。本記事ではポテンシャルをこの積分表示に固定して扱います。

# まとめ

クーロンの法則とビオ＝サバールの法則は、同じ核$(\boldsymbol x-\boldsymbol y)/|\boldsymbol x-\boldsymbol y|^3$を使う2つの積分でした。源を$\rho-\boldsymbol J/c$にまとめて幾何積で核に掛けると、1本の積分公式$F=\frac1{4\pi\varepsilon_0}\int\frac{\boldsymbol x-\boldsymbol y}{|\boldsymbol x-\boldsymbol y|^3}(\rho-\boldsymbol J/c)\,dV$になり、ベクトル部が電場、2ベクトル部が磁場の$Ic$倍、スカラー部は定常電流の条件で$0$になります。

- **核と定数**：核は[[7shi-cla5]]のコーシー核の$n=3$の場合であり、$4\pi$は$|S^2|$です。
- **磁場**：$\boldsymbol r\wedge\boldsymbol J$の重ね合わせとして2ベクトルで書け、擬ベクトルという例外的な扱いは不要になります。
- **微分形**：核が$D$の基本解であることから$DF=(\rho-\boldsymbol J/c)/\varepsilon_0$が成り立ち、そのスカラー・ベクトル・2ベクトル・擬スカラー部が静的な場の4本の式です。
- **積分形**：ガウスの法則は領域の基本定理のスカラー部、アンペールの法則は曲面の基本定理の擬スカラー部です。
- **ポテンシャル**：$F=-DP$、$P=\varphi-c\boldsymbol A$と書け、スカラー部の消滅は$\nabla\cdot\boldsymbol A=0$と言い換えられます。
