幾何微積分の基本定理[[7shi-cla2]]にモノジェニック性を代入して、コーシーの積分定理・積分公式・留数定理を導きます。積分公式の核は、ディラック作用素について左右ともにモノジェニックなベクトルの逆元$\boldsymbol x^{-1}$です[[7shi-cla1]]。小さな円の上では核と外向きの法線の積が半径の逆数になり、周に沿った積分から$S^1$の長さ$2\pi$が現れます。複素解析の公式に現れる$2\pi i$は、この周長と、接線を法線に回す擬スカラー$\omega$に分かれます。

# 準備

[[7shi-cla1]]・[[7shi-cla2]]と同じく、$\operatorname{Cl}_{2,0}(\mathbb R)$の生成元を$e_1,e_2$（$e_1^2=e_2^2=1$、$e_1e_2=-e_2e_1$）、擬スカラーを$\omega=e_1e_2$（$\omega^2=-1$）とし、ディラック作用素を

$$
D=e_1\partial_x+e_2\partial_y
$$

とします。左からの作用を$DF=e_1\partial_xF+e_2\partial_yF$、右からの作用を$GD=(\partial_xG)e_1+(\partial_yG)e_2$と書き、$DF=0$を満たす$F$を左モノジェニック、$GD=0$を満たす$G$を右モノジェニックと呼びます。単にモノジェニックと言えば左モノジェニックを指します。位置ベクトルは$\boldsymbol x=xe_1+ye_2$、対応する複素数は偶部分の元$z=e_1\boldsymbol x=x+\omega y$です。

[[7shi-cla2]]の設定を引き継ぎます。平面の有界な領域$M$の境界$\partial M$は区分的に滑らかな閉曲線（穴があれば複数）とし、領域を左手に見る向き（外側の境界は反時計回り、穴の境界は時計回り）にたどります。$\boldsymbol n$は外向きの単位法線、$ds$は弧長の要素、$d\boldsymbol x$は有向線素です。関数$F,G$は、$M$とその境界を含む開集合上で定義され、$\operatorname{Cl}_{2,0}(\mathbb R)$に値を取る$C^1$級関数とします。特異点を扱う箇所では、その点を除いた開集合上でこの条件を課します。

[[7shi-cla2]]で示した次の式を使います。一般の領域については、[[7shi-cla2]]と同じく証明の概略にとどめます。

&&&fml 基本定理と両側形式
$$
\begin{aligned}
\int_M DF\,dA&=\oint_{\partial M}\boldsymbol nF\,ds \\
\int_M\bigl((GD)F+G(DF)\bigr)dA&=\oint_{\partial M}G\,\boldsymbol nF\,ds \\
\boldsymbol n\,ds&=\omega\,d\boldsymbol x
\end{aligned}
$$
&&&

# コーシーの積分定理

基本定理の左辺は$DF$の積分なので、$F$がモノジェニックなら消えます。

&&&thm コーシーの積分定理 [thm-cauchy]
$M$上で$DF=0$なら

$$
\oint_{\partial M}\boldsymbol nF\,ds=0
$$
&&&

$F$の値は偶部分に限らず、$\operatorname{Cl}_{2,0}(\mathbb R)$の任意の元でかまいません。法線を掛けた境界積分の形は、次元によらず同じ形で書けます。

## 有向線素による形

$\boldsymbol n\,ds=\omega\,d\boldsymbol x$を代入して$\omega^{-1}=-\omega$を左から掛けると$\oint d\boldsymbol x\,F=0$です。さらに$e_1$を左から掛けます。$e_1$は定数なので積分の中に入れることができ、有向線素に$e_1$を掛けたものは$z=e_1\boldsymbol x$の微分

$$
dz=e_1\,d\boldsymbol x=dx+\omega\,dy
$$

になります。

&&&fml 複素数による形
$$
\oint_{\partial M}dz\,F=0
$$
&&&

偶部分に値を取る$F$では$dz$と$F$が可換なので、$\oint F\,dz=0$とも書けます。これが複素解析で知られるコーシーの積分定理の形です。複素解析では、この定理は複素数値の線積分についての命題として述べられます。ここでは、基本定理の左辺がモノジェニック性で消えるという1行の帰結であり、$dz$は法線$\boldsymbol n\,ds$を$\omega$と$e_1$で書き換えたものにあたります。

## ベクトル値の関数

ベクトル値の関数$F=ae_1+be_2$では、[[7shi-cla2]]のグリーンの定理の2つの形により、境界積分のスカラー部と2ベクトル部がそれぞれ流束と循環です。

$$
\oint_{\partial M}\boldsymbol nF\,ds=\oint_{\partial M}(a\,dy-b\,dx)+\omega\oint_{\partial M}(a\,dx+b\,dy)
$$

モノジェニックなベクトル場は発散と回転がともに消えるので、$\partial M$を横切る流束と$\partial M$に沿った循環がともに$0$です。積分定理はこの2つを1つにまとめた形になっています。

## 両側の形

両側形式の左辺は、$G$が右モノジェニック、$F$が左モノジェニックなら消えます。

&&&thm 両側の積分定理 [thm-cauchy2]
$M$上で$GD=0$かつ$DF=0$なら

$$
\oint_{\partial M}G\,\boldsymbol nF\,ds=0
$$
&&&

$G=1$とすれば[積分定理](#thm-cauchy)に戻ります。積分公式を導くには、$G$として右モノジェニックな関数を選び、特異点のまわりを除いた領域にこの形を当てはめます。

# 円環上の計算

積分公式では、特異点を中心とする小さな円を除いた領域を使います。その境界は円なので、[[7shi-cla2]]の長方形による証明がそのままは使えません。ここでは極座標で計算し、円板から同心の円板を除いた円環であれば、1次元の基本定理だけで両側形式が示せることを見ます。

## 極座標のディラック作用素

点$\boldsymbol a$を中心とする極座標を、次のように取ります。

$$
\boldsymbol x-\boldsymbol a=r\boldsymbol n,\qquad
\boldsymbol n=\cos\theta\,e_1+\sin\theta\,e_2,\qquad
\boldsymbol t=-\sin\theta\,e_1+\cos\theta\,e_2
$$

$\boldsymbol n$は中心から外向きの単位ベクトル、$\boldsymbol t$は反時計回りの単位接線です。$\partial_\theta\boldsymbol n=\boldsymbol t$、$\partial_\theta\boldsymbol t=-\boldsymbol n$であり、$\boldsymbol n\omega=\cos\theta\,e_2-\sin\theta\,e_1=\boldsymbol t$です。

&&&fml 極座標のディラック作用素
$$
DF=\boldsymbol n\Bigl(\partial_r+\frac \omega r\,\partial_\theta\Bigr)F,\qquad
GD=(\partial_rG)\,\boldsymbol n+\frac1r(\partial_\theta G)\,\boldsymbol t
$$
&&&

&&&prf
$\partial_r=\cos\theta\,\partial_x+\sin\theta\,\partial_y$、$r^{-1}\partial_\theta=-\sin\theta\,\partial_x+\cos\theta\,\partial_y$を$\partial_x,\partial_y$について解くと

$$
\partial_x=\cos\theta\,\partial_r-\frac{\sin\theta}r\,\partial_\theta,\qquad
\partial_y=\sin\theta\,\partial_r+\frac{\cos\theta}r\,\partial_\theta
$$

であり、$D=e_1\partial_x+e_2\partial_y=\boldsymbol n\,\partial_r+r^{-1}\boldsymbol t\,\partial_\theta$となる。$\boldsymbol t=\boldsymbol n\omega$を代入すれば第1式を得る。右からの作用も係数が同じなので第2式となる。
&&&

微分が動径方向と角度方向に分かれ、それぞれの係数が法線$\boldsymbol n$と接線$\boldsymbol t$になっています。$\boldsymbol n$は可逆なので、$DF=0$は$\partial_rF=-(\omega/r)\partial_\theta F$と同値です。偶部分に値を取る関数では、これが極座標のコーシー＝リーマンの方程式にあたります。

## 円環での両側形式

&&&fml 極座標の恒等式
$$
\partial_r\bigl(r\,G\boldsymbol nF\bigr)+\partial_\theta\bigl(G\boldsymbol tF\bigr)
=r\bigl((GD)F+G(DF)\bigr)
$$
&&&

&&&prf
積の微分法則と$\partial_\theta\boldsymbol n=\boldsymbol t$、$\partial_\theta\boldsymbol t=-\boldsymbol n$より

$$
\begin{aligned}
\partial_r(rG\boldsymbol nF)&=G\boldsymbol nF+r(\partial_rG)\boldsymbol nF+rG\boldsymbol n\,\partial_rF \\
\partial_\theta(G\boldsymbol tF)&=(\partial_\theta G)\boldsymbol tF-G\boldsymbol nF+G\boldsymbol t\,\partial_\theta F
\end{aligned}
$$

である。足し合わせると$G\boldsymbol nF$が打ち消し、残りは

$$
r\Bigl((\partial_rG)\boldsymbol n+\frac1r(\partial_\theta G)\boldsymbol t\Bigr)F
+rG\Bigl(\boldsymbol n\,\partial_r+\frac1r\boldsymbol t\,\partial_\theta\Bigr)F
$$

となる。極座標のディラック作用素の式より、これは$r\bigl((GD)F+G(DF)\bigr)$である。
&&&

面積要素は$dA=r\,dr\,d\theta$です。恒等式を$\varepsilon\le r\le R$、$0\le\theta\le2\pi$で積分すると、$\theta$の項は周期性により1次元の基本定理で消え、$r$の項は1次元の基本定理により両端の値の差になります。

&&&thm 円環での両側形式 [thm-annulus]
中心$\boldsymbol a$、半径$\varepsilon<R$の円環$M$で

$$
\int_M\bigl((GD)F+G(DF)\bigr)dA
=\oint_{r=R}G\boldsymbol nF\,ds-\oint_{r=\varepsilon}G\boldsymbol nF\,ds
$$

ここで$\boldsymbol n$はどちらの円でも中心から外向きの法線、$ds=r\,d\theta$とします。
&&&

内側の円では、円環から見た外向きの法線は$-\boldsymbol n$なので、これは両側形式を円環に当てはめたものと一致します。$F,G$は閉円環を含む開集合上で$C^1$級であれば十分で、中心$\boldsymbol a$で定義されている必要はありません。$G=1$、$DF=0$とすると、$\oint_{r=\rho}\boldsymbol nF\,ds$が半径$\rho$によらないことがわかります。

# コーシーの積分公式

## 小円上の核

[[7shi-cla1]]で見たとおり、$\boldsymbol x^{-1}=\boldsymbol x/|\boldsymbol x|^2$は原点を除いて左右ともにモノジェニックです。$D$は平行移動と可換なので、$(\boldsymbol x-\boldsymbol a)^{-1}$も$\boldsymbol a$を除いて左右ともにモノジェニックです。

$\boldsymbol a$を中心とする半径$r$の円の上では、$\boldsymbol x-\boldsymbol a=r\boldsymbol n$なので

&&&fml 円の上の核と法線
$$
(\boldsymbol x-\boldsymbol a)^{-1}\boldsymbol n=\boldsymbol n(\boldsymbol x-\boldsymbol a)^{-1}=\frac{\boldsymbol n\boldsymbol n}r=\frac1r
$$
&&&

となります。核と法線は平行なので、積はベクトルではなくスカラーになります。この円の上で積分すると

$$
\oint_{r=\rho}(\boldsymbol x-\boldsymbol a)^{-1}\boldsymbol n\,ds=\frac1\rho\cdot2\pi\rho=2\pi
$$

であり、半径によりません。核が距離に反比例して減衰する分と、円周が半径に比例して伸びる分が打ち消し合い、単位円$S^1$の長さ$2\pi$だけが残ります。

## コーシー＝ポンペイウの公式

モノジェニックとは限らない$F$に対して、境界の値と$DF$の積分から内部の値を表す公式が得られます。

&&&thm コーシー＝ポンペイウの公式 [thm-pompeiu]
$\boldsymbol a$が$M$の内部の点なら

$$
F(\boldsymbol a)=\frac1{2\pi}\left(\oint_{\partial M}(\boldsymbol x-\boldsymbol a)^{-1}\boldsymbol nF\,ds
-\int_M(\boldsymbol x-\boldsymbol a)^{-1}DF\,dA\right)
$$
&&&

&&&prf
$G=(\boldsymbol x-\boldsymbol a)^{-1}$とし、$M$から中心$\boldsymbol a$、半径$\varepsilon$の閉円板を除いた領域$M_\varepsilon$に両側形式を当てはめる。$M_\varepsilon$上では$GD=0$なので

$$
\int_{M_\varepsilon}G(DF)\,dA=\oint_{\partial M}G\boldsymbol nF\,ds-\oint_{r=\varepsilon}G\boldsymbol nF\,ds
$$

である。小円の項では、$M_\varepsilon$から見た外向きの法線が$-\boldsymbol n$であることを符号に取り込んだ。小円の上では$G\boldsymbol n=1/\varepsilon$、$ds=\varepsilon\,d\theta$なので

$$
\oint_{r=\varepsilon}G\boldsymbol nF\,ds=\int_0^{2\pi}F(\boldsymbol a+\varepsilon\boldsymbol n)\,d\theta
\ \longrightarrow\ 2\pi F(\boldsymbol a)\qquad(\varepsilon\to0)
$$

となる。左辺では$|G|=1/r$、$dA=r\,dr\,d\theta$より、極座標で見た被積分関数$G(DF)\,r$は$\boldsymbol a$の近くで有界であり、$\varepsilon\to0$で$M$上の積分に収束する。以上を整理すれば公式を得る。

$M$が$\boldsymbol a$を中心とする円板なら、$M_\varepsilon$は円環であり、[円環での両側形式](#thm-annulus)により1次元の基本定理だけで示される。一般の$M$では、[[7shi-cla2]]の両側形式の証明と同じく概略である。
&&&

$DF=0$なら面積分の項が消え、境界の値だけで内部の値が決まります。

&&&cor コーシーの積分公式 [cor-cauchy]
$M$上で$DF=0$なら、$M$の内部の点$\boldsymbol a$で

$$
F(\boldsymbol a)=\frac1{2\pi}\oint_{\partial M}(\boldsymbol x-\boldsymbol a)^{-1}\boldsymbol nF\,ds
$$
&&&

モノジェニック関数は境界の値で決まり、核は$(\boldsymbol x-\boldsymbol a)^{-1}/2\pi$です。[[7shi-cla1]]で先送りにした、$\boldsymbol x^{-1}$による積分での特徴付けがこれにあたります。$2\pi$は小円の極限から来ており、その値は核と法線の積$1/\varepsilon$を周に沿って積分した$S^1$の長さです。

&&&rem 右モノジェニックな関数
[コーシー＝ポンペイウの公式](#thm-pompeiu)の証明で、両側形式の$G$の位置に$F$を、$F$の位置に$(\boldsymbol x-\boldsymbol a)^{-1}$を置くと、右からの作用に対する公式が得られます。

$$
F(\boldsymbol a)=\frac1{2\pi}\left(\oint_{\partial M}F\boldsymbol n(\boldsymbol x-\boldsymbol a)^{-1}ds
-\int_M(FD)(\boldsymbol x-\boldsymbol a)^{-1}dA\right)
$$

核が左右ともにモノジェニックであるため、どちらの側にも同じ核が使えます。
&&&

## 複素数による形

核と法線の積を複素数で書き直します。$\boldsymbol a$に対応する複素数を$w=e_1\boldsymbol a$とすると、$z-w=e_1(\boldsymbol x-\boldsymbol a)$より$(\boldsymbol x-\boldsymbol a)^{-1}=(z-w)^{-1}e_1$です。

&&&fml 核と複素数の対応
$$
(\boldsymbol x-\boldsymbol a)^{-1}\boldsymbol n\,ds=\omega^{-1}\frac{dz}{z-w}
$$
&&&

&&&prf
$\boldsymbol n\,ds=\omega\,d\boldsymbol x$と$d\boldsymbol x=e_1\,dz$より、$(\boldsymbol x-\boldsymbol a)^{-1}\boldsymbol n\,ds=(z-w)^{-1}e_1\omega e_1\,dz$である。$e_1\omega e_1=e_1e_1e_2e_1=-\omega=\omega^{-1}$であり、偶部分の元どうしは可換なので右辺を得る。
&&&

偶部分に値を取るモノジェニック関数$F$では、[積分公式](#cor-cauchy)は次の形になります。

$$
F(\boldsymbol a)=\frac1{2\pi \omega}\oint_{\partial M}\frac{F\,dz}{z-w}
$$

複素解析のコーシーの積分公式です。そこでは$2\pi i$が1つの定数として現れますが、ここでは出所の異なる2つの因子に分かれます。$2\pi$は核と法線の積から来る$S^1$の長さであり、$\omega$は法線付きの境界要素と有向線素を結ぶ関係$\boldsymbol n\,ds=\omega\,d\boldsymbol x$から来ています。これは、接線を法線へ時計回りに$90^\circ$回す左乗法です。

## 平均値の性質

$M$を$\boldsymbol a$を中心とする半径$R$の円板とすると、境界の上で$(\boldsymbol x-\boldsymbol a)^{-1}\boldsymbol n=1/R$なので、[積分公式](#cor-cauchy)の核と法線の積が定数$1/R$になります。

&&&cor 平均値の性質
$\boldsymbol a$を中心とする半径$R$の閉円板を含む開集合上で$DF=0$なら

$$
F(\boldsymbol a)=\frac1{2\pi R}\oint_{r=R}F\,ds=\frac1{2\pi}\int_0^{2\pi}F(\boldsymbol a+R\boldsymbol n)\,d\theta
$$
&&&

中心での値は、円周上の値の平均です。この場合は$M_\varepsilon$が円環なので、証明は1次元の基本定理だけで完結しています。

## 基本解

$F$がある閉円板の外で$0$になる$C^1$級関数なら、$M$をその円板を内部に含む領域に取ると境界の項が消え、[コーシー＝ポンペイウの公式](#thm-pompeiu)は次の形になります。$E=\boldsymbol x^{-1}/2\pi$と置くと、$\boldsymbol x^{-1}$は奇関数なので$-E(\boldsymbol x-\boldsymbol a)=E(\boldsymbol a-\boldsymbol x)$です。

&&&fml 台がコンパクトな関数の復元
$$
F(\boldsymbol a)=-\frac1{2\pi}\int(\boldsymbol x-\boldsymbol a)^{-1}DF\,dA=\int E(\boldsymbol a-\boldsymbol x)\,DF(\boldsymbol x)\,dA
$$
&&&

この積分作用素は、台がコンパクトな$C^1$級関数$F$を、その$DF$から復元します。

&&&rem 超関数による解釈
超関数の言葉では、この式は$E$が$DE=ED=\delta$を満たす$D$の基本解であることを表します（$ED=\delta$が上の式、$DE=\delta$が右からの作用に対する公式にあたります）。本記事では超関数の枠組みには立ち入りません。
&&&

# 留数

## 留数の定義

点$\boldsymbol a$の近く（$\boldsymbol a$自身を除く）でモノジェニックな$F$を考えます。[円環での両側形式](#thm-annulus)で$G=1$とすると、$\boldsymbol a$を中心とする小円上の$\oint\boldsymbol nF\,ds$は半径によりません。

&&&def 留数
$F$が$0<|\boldsymbol x-\boldsymbol a|<\rho$でモノジェニックなとき、$0<\varepsilon<\rho$に対して

$$
\operatorname{Res}_{\boldsymbol a}F=\frac1{2\pi}\oint_{|\boldsymbol x-\boldsymbol a|=\varepsilon}\boldsymbol nF\,ds
$$

を$\boldsymbol a$における$F$の**留数**と定めます。値は$\varepsilon$によりません。
&&&

複素解析では留数をローラン展開の係数として導入するのが通例ですが、ここでは小円上の積分そのものを定義とします。冪級数を使わないため、級数展開を前提としません。

この留数はクリフォード代数に値を取ります。偶部分に値を取る関数でも、法線を左から掛けるため留数は奇部分に入り、複素解析の留数とは$e_1$による対応を介して一致します（後述）。

&&&thm 留数定理
$M$とその境界を含む開集合$U$と、$M$の内部の相異なる有限個の点$\boldsymbol a_1,\dots,\boldsymbol a_m$を取ります。$F$が$U\setminus\{\boldsymbol a_1,\dots,\boldsymbol a_m\}$上の$C^1$級関数で、$M\setminus\{\boldsymbol a_1,\dots,\boldsymbol a_m\}$上でモノジェニックなら

$$
\oint_{\partial M}\boldsymbol nF\,ds=2\pi\sum_{k=1}^m\operatorname{Res}_{\boldsymbol a_k}F
$$
&&&

&&&prf
各$\boldsymbol a_k$を中心とする小さな閉円板を$M$から除いた領域に[積分定理](#thm-cauchy)を当てはめる。小円の上では、除いた領域から見た外向きの法線は中心から内向きなので、小円の項は$-\oint\boldsymbol nF\,ds=-2\pi\operatorname{Res}_{\boldsymbol a_k}F$となる。これらと$\partial M$の項の和が$0$であることから主張を得る。
&&&

## 核とその導関数の留数

核に右から定数$c\in\operatorname{Cl}_{2,0}(\mathbb R)$を掛けた$(\boldsymbol x-\boldsymbol a)^{-1}c$はモノジェニックです（$D$は右の定数に作用しません）。円の上の核と法線の関係$\boldsymbol n(\boldsymbol x-\boldsymbol a)^{-1}=1/\varepsilon$から

$$
\operatorname{Res}_{\boldsymbol a}\bigl((\boldsymbol x-\boldsymbol a)^{-1}c\bigr)=\frac1{2\pi}\cdot\frac1\varepsilon\cdot2\pi\varepsilon\,c=c
$$

です。核を偏微分したものも、$D$が偏微分と可換なのでモノジェニックですが、留数は$0$です。

&&&prop 核の導関数の留数
$\partial^\alpha$を$\partial_x,\partial_y$を1回以上合成した偏微分とすると

$$
\operatorname{Res}_{\boldsymbol a}\bigl(\partial^\alpha(\boldsymbol x-\boldsymbol a)^{-1}c\bigr)=0
$$
&&&

&&&prf
$\boldsymbol a=0$としてよい。$\partial^\alpha$が$m\ge1$回の微分なら、$K=\partial^\alpha\boldsymbol x^{-1}$は$-1-m$次の同次関数であり、$K(\varepsilon\boldsymbol n)=\varepsilon^{-1-m}K(\boldsymbol n)$である。$ds=\varepsilon\,d\theta$より

$$
\oint_{|\boldsymbol x|=\varepsilon}\boldsymbol nKc\,ds=\varepsilon^{-m}\oint_{|\boldsymbol x|=1}\boldsymbol nKc\,ds
$$

となる。左辺は$\varepsilon$によらないので、右辺の積分は$0$でなければならない。
&&&

したがって$F$が$\boldsymbol a$の近くで

$$
F=(\boldsymbol x-\boldsymbol a)^{-1}c+\sum_\alpha\partial^\alpha(\boldsymbol x-\boldsymbol a)^{-1}c_\alpha+(\boldsymbol a\text{でもモノジェニックな関数})
$$

と有限和で書けるなら、$\operatorname{Res}_{\boldsymbol a}F=c$です。最後の項の留数が$0$であることは[積分定理](#thm-cauchy)によります。留数は、核$(\boldsymbol x-\boldsymbol a)^{-1}$そのものの係数として読み取れます。

偶部分に値を取る$F$について、複素解析の留数を$\operatorname{res}_wF=\frac1{2\pi \omega}\oint F\,dz$（$w=e_1\boldsymbol a$、積分は$\boldsymbol a$を中心とする小円上）と書いて区別します。$\boldsymbol n\,ds=\omega\,d\boldsymbol x=\omega e_1\,dz=-e_1\omega\,dz$より

$$
\operatorname{Res}_{\boldsymbol a}F=-e_1\omega\cdot\frac1{2\pi}\oint dz\,F=e_1\operatorname{res}_wF
$$

であり、複素解析の留数に$e_1$を左から掛けたものになります。たとえば$1/z=\boldsymbol x^{-1}e_1$の留数は、核の係数を読み取って$e_1$です。複素解析での値$1$と$e_1$だけ違うのは、偶部分の関数を核（ベクトル）の係数として書いたためであり、[[7shi-cla1]]で見た$e_1$による偶部分と奇部分の対応と同じものです。

&&&ex 有理関数の留数
偶部分の関数$F=1/(1+z^2)$の特異点は$z=\pm \omega$、すなわち$\boldsymbol x=e_1(\pm \omega)=\pm e_2$です。$z\mp \omega=e_1(\boldsymbol x\mp e_2)$より$(z\mp \omega)^{-1}=(\boldsymbol x\mp e_2)^{-1}e_1$なので、部分分数分解は

$$
F=\frac1{2\omega}\left(\frac1{z-\omega}-\frac1{z+\omega}\right)
=(\boldsymbol x-e_2)^{-1}e_1\frac1{2\omega}-(\boldsymbol x+e_2)^{-1}e_1\frac1{2\omega}
$$

となります。核の係数を読み取ると

$$
\operatorname{Res}_{e_2}F=-\frac{e_1\omega}2,\qquad\operatorname{Res}_{-e_2}F=\frac{e_1\omega}2
$$

です。複素解析の値$\mp \omega/2$に$e_1$を左から掛けたものになっています。

原点を中心とする半径$R>1$の円の上の$\oint\boldsymbol nF\,ds$は、留数定理により$2\pi\times$（留数の和）$=0$です。一方、この積分は$R$によらず、$|F|$は$1/R^2$程度で減衰するので積分は$1/R$程度で$0$に近づきます。どちらの見方でも$0$になります。
&&&

# まとめ

基本定理$\int_MDF\,dA=\oint_{\partial M}\boldsymbol nF\,ds$にモノジェニック性を代入すると、積分定理$\oint\boldsymbol nF\,ds=0$が得られます。$\boldsymbol n\,ds=\omega\,d\boldsymbol x$と$dz=e_1\,d\boldsymbol x$で書き換えると、複素解析の$\oint F\,dz=0$になります。

積分公式の核はベクトルの逆元$(\boldsymbol x-\boldsymbol a)^{-1}$です。小円の上では核と法線の積が$1/\varepsilon$になり、両側形式で小円の極限を取るとコーシー＝ポンペイウの公式が得られます。モノジェニック関数ではコーシーの積分公式、円板では平均値の性質になり、境界の項を消せば核が$D$の基本解であることがわかります。

| 因子 | 出所 |
|---|---|
| $2\pi$ | 核と法線の積$1/\varepsilon$を小円に沿って積分した、$S^1$の長さ |
| $\omega$ | 接線を法線に回す左乗法$\boldsymbol n\,ds=\omega\,d\boldsymbol x$ |

留数は小円上の積分として定義し、留数定理は積分定理から直ちに従います。核の導関数の留数は同次性により$0$です。したがって、核とその導関数の有限和に、特異点でもモノジェニックな関数を加えた形で書ける場合、留数は核$(\boldsymbol x-\boldsymbol a)^{-1}$そのものの右側の係数として読み取れます。偶部分の関数では、複素解析の留数に$e_1$を左から掛けたものです。

対数の枝や多価性は、定義域の位相を含む追加の議論になるため、本記事では扱いません。
