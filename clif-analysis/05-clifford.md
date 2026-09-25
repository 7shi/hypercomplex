複素解析をディラック作用素の核として組み直した枠組み[[7shi-cla1]][[7shi-cla2]][[7shi-cla3]]と、その4次元版である四元数解析[[7shi-cla4]]を、一般の次元$n$で述べ直します。コーシーの積分公式の核は$\boldsymbol x/|\boldsymbol x|^n$、正規化定数は単位球面$S^{n-1}$の面積$2\pi^{n/2}/\Gamma(n/2)$になり、$n=2$で[[7shi-cla3]]の$\boldsymbol x^{-1}$と$2\pi$、$n=4$で[[7shi-cla4]]の核と$2\pi^2$に戻ります。複素解析の$2\pi i$は、$S^{n-1}$の面積と、有向境界要素を法線に移す擬スカラーの2つの因子に分かれます。恒等写像の正則性を壊す係数、逆元の微分に現れる係数、フューターの定理のずれの係数は、どれも同じ数$n-2$です。ラプラシアンを反復して$1/z$から核を作るときの回数$(n-2)/2$も、この数から次数の差として決まります。

# 一般次元のディラック作用素

## 設定

[[7shi-cla2]]の一般次元の設定をそのまま使います（$n\ge2$）。$\operatorname{Cl}_{n,0}(\mathbb R)$の生成元を$e_0,e_1,\dots,e_{n-1}$（$e_a^2=1$、$a\ne b$なら$e_ae_b=-e_be_a$）とし、位置ベクトルとディラック作用素を

$$
\boldsymbol x=\sum_{a=0}^{n-1}x_ae_a,\qquad D=\sum_{a=0}^{n-1}e_a\partial_a
$$

とします。添字を$0$から振るのは、後で[[7shi-cla4]]と同じく$e_0$を特別な方向として使うためです。$D^2=\Delta$（$n$次元のラプラシアン）は[[7shi-cla1]]と同じ計算で成り立ちます。

&&&def モノジェニック関数
$\mathbb R^n$の開集合上の$\operatorname{Cl}_{n,0}(\mathbb R)$値の$C^1$級関数$F$が$DF=0$を満たすとき、$F$を**左モノジェニック**と呼びます。右からの作用を$FD=\sum_a(\partial_aF)e_a$とし、$FD=0$を満たすとき**右モノジェニック**と呼びます。単にモノジェニックと言えば左モノジェニックを指します。
&&&

$D^2=\Delta$より、$C^2$級のモノジェニック関数の各成分は$n$次元の調和関数です。

&&&rem 外微分と余微分
[[7shi-cla1]]で見たとおり、$k$ベクトル値の$F$に対して$DF=D\cdot F+D\wedge F$と分かれます。基底の$k$ベクトル$e_{a_1}\wedge\dots\wedge e_{a_k}$を$k$次微分形式$dx_{a_1}\wedge\dots\wedge dx_{a_k}$と同一視すると、$D\wedge F$は外微分$dF$そのものです。$D\cdot F$は余微分$\delta F$の符号を変えたものになります。ここで$\delta$は$d$の形式的な随伴で、ホッジの星印を$\star F=\tilde FI$（$\tilde F$は積の順序を反転するリバージョン、$I=e_0e_1\cdots e_{n-1}$は擬スカラー）として$\delta=(-1)^{n(k+1)+1}\star d\star$と書けます。したがって

$$
D=d-\delta,\qquad D^2=-(d\delta+\delta d)=\Delta
$$

です（$d^2=\delta^2=0$）。$d\delta+\delta d$はホッジ・ラプラシアンで、ユークリッド空間では$-\Delta$に一致します。文献ではディラック作用素を$d+\delta$と書くことが多いですが、これはベクトルの作用の違いによります。本記事の$e_a^2=+1$では、ベクトルを左から掛ける作用が外積と内部積の和になり、$D=d-\delta$となります。$e_a^2=-1$の規約では外積から内部積を引いた作用になり、$d+\delta$が現れます。規約を一貫させれば、議論は対応する符号変更のもとで進められます（[[7shi-cla1]]）。[[7shi-cla2]]の微分形式との比較で述べた「次数の異なる形式をまとめて扱う」ことを、作用素の側で行ったものがこの$D$にあたります。
&&&

## 恒等写像と逆元

[[7shi-cla4]]では、生成元の挟み込み$\sum_ae_ae_0e_a=(2-n)e_0$を一般の$n$で示し、そこから$D(e_0\boldsymbol x)=(2-n)e_0$と

$$
D\boldsymbol x^{-1}=\frac{n-2}{|\boldsymbol x|^2}
$$

を得ました。$e_0\boldsymbol x$と$\boldsymbol x^{-1}$がモノジェニックになるのは$n=2$のときだけです。以下の議論では、この$n-2$という係数が繰り返し現れます。

# コーシー核

## 核の同次性

$\boldsymbol x^{-1}=\boldsymbol x/|\boldsymbol x|^2$の分母の冪を$n$に替えると、モノジェニックになります。

&&&prop コーシー核 [prop-kernel]
$$
E(\boldsymbol x)=\frac{\boldsymbol x}{|\boldsymbol x|^n}
$$

は$\mathbb R^n\setminus\{0\}$で左右ともにモノジェニックです。
&&&

&&&prf
$D\boldsymbol x=n$、$D|\boldsymbol x|^{-n}=-n\boldsymbol x/|\boldsymbol x|^{n+2}$より

$$
DE=\frac{D\boldsymbol x}{|\boldsymbol x|^n}+(D|\boldsymbol x|^{-n})\boldsymbol x=\frac n{|\boldsymbol x|^n}-\frac{n\boldsymbol x^2}{|\boldsymbol x|^{n+2}}=0
$$

である。$\boldsymbol x^2=|\boldsymbol x|^2$を使った。$|\boldsymbol x|^{-n}$はスカラーなので、右からの作用$ED$も同じ計算になる。
&&&

$n=2$では$E=\boldsymbol x^{-1}$であり、[[7shi-cla3]]の核に一致します。$n=4$では[[7shi-cla4]]のコーシー＝フューター核$\bar q/|q|^4=(\boldsymbol x/|\boldsymbol x|^4)e_0$の、$e_0$を掛ける前の形です。

$E$は$-(n-1)$次の同次関数です。この次数は、[[7shi-cla3]]・[[7shi-cla4]]と同じ理由で必然的に決まります。半径$\varepsilon$の球面の面積は$\varepsilon^{n-1}$に比例するので、同次な核で、球面上の積分が半径によらない$0$でない値になることを要求すると、次数は$-(n-1)$でなければなりません。

## ラプラシアンの基本解との関係

核は、スカラーの調和関数を$D$で微分しても得られます。

&&&fml 核とポテンシャル
$$
D|\boldsymbol x|^{2-n}=(2-n)\frac{\boldsymbol x}{|\boldsymbol x|^n}\quad(n\ge3),\qquad
D\log|\boldsymbol x|=\frac{\boldsymbol x}{|\boldsymbol x|^2}\quad(n=2)
$$

$|\boldsymbol x|^{2-n}$と$\log|\boldsymbol x|$は、それぞれ$\mathbb R^n\setminus\{0\}$、$\mathbb R^2\setminus\{0\}$の調和関数です。
&&&

&&&prf
スカラー関数$\phi$に対して$D\phi=\sum_ae_a\partial_a\phi$は勾配である。$\partial_a|\boldsymbol x|^p=p\,x_a|\boldsymbol x|^{p-2}$より$D|\boldsymbol x|^p=p\,\boldsymbol x|\boldsymbol x|^{p-2}$であり、$p=2-n$とすれば第1式を得る。$\log|\boldsymbol x|=\frac12\log|\boldsymbol x|^2$より$D\log|\boldsymbol x|=\boldsymbol x/|\boldsymbol x|^2$である。調和性は、$E$がモノジェニックであることから$\Delta\phi=D(D\phi)=0$として従う。
&&&

$\log|\boldsymbol x|$は1価の実関数であり、[[7shi-cla3]]で扱わないとした複素対数の枝や多価性とは別のものです。

$D^2=\Delta$なので、核は「ラプラシアンの基本解の$D$微分」として位置づけられます。

&&&rem 基本解
$n\ge3$では$\Phi=|\boldsymbol x|^{2-n}/\bigl((2-n)|S^{n-1}|\bigr)$、$n=2$では$\Phi=\log|\boldsymbol x|/2\pi$がラプラシアンの基本解（$\Delta\Phi=\delta$）です。ここで$|S^{n-1}|$は単位球面の面積です。上の式から$D\Phi=E/|S^{n-1}|$であり、$D(D\Phi)=\Delta\Phi=\delta$は、$E/|S^{n-1}|$が$D$の基本解であることを表します。[[7shi-cla3]]では$E=\boldsymbol x^{-1}/2\pi$を$D$の基本解として述べましたが、これは$\log|\boldsymbol x|/2\pi$の$D$微分です。本記事でも超関数の枠組みには立ち入りません。
&&&

# 球面の面積

積分公式の正規化定数を求めます。以下、$|S^{n-1}|$は$\mathbb R^n$の単位球面の$(n-1)$次元の面積を表します。

&&&fml 単位球面の面積
$$
|S^{n-1}|=\frac{2\pi^{n/2}}{\Gamma(n/2)}
$$
&&&

&&&prf
ガウス積分$\int_{-\infty}^\infty e^{-t^2}dt=\sqrt\pi$の$n$個の積として

$$
\int_{\mathbb R^n}e^{-|\boldsymbol x|^2}dV=\pi^{n/2}
$$

である。極座標$dV=\rho^{n-1}d\rho\,d\Omega$（$d\Omega$は$S^{n-1}$の面積要素）で書くと、左辺は

$$
|S^{n-1}|\int_0^\infty\rho^{n-1}e^{-\rho^2}d\rho=|S^{n-1}|\cdot\frac12\Gamma\!\left(\frac n2\right)
$$

となる（$s=\rho^2$と置換した）。両辺を比べれば公式を得る。
&&&

$\Gamma(1/2)=\sqrt\pi$、$\Gamma(1)=1$、$\Gamma(s+1)=s\Gamma(s)$から値が求まります。

[[7shi-cla2]]の基本定理$\int_MDF\,dV=\oint_{\partial M}\boldsymbol nF\,dS$を単位球$B^n$で$F=\boldsymbol x$として使うと、$D\boldsymbol x=n$、球面上で$\boldsymbol n\boldsymbol x=\boldsymbol n^2=1$より

$$
n\,|B^n|=|S^{n-1}|
$$

です（$|B^n|$は単位球の体積）。[[7shi-cla2]]では$n=2,3$の例として示したものが、一般の$n$で球の体積を与えます。

&&&ex 低次元の値
| $n$ | $\lvert S^{n-1}\rvert$ | $\lvert B^n\rvert$ |
|---|---|---|
| 2 | $2\pi$ | $\pi$ |
| 3 | $4\pi$ | $4\pi/3$ |
| 4 | $2\pi^2$ | $\pi^2/2$ |
| 5 | $8\pi^2/3$ | $8\pi^2/15$ |
| 6 | $\pi^3$ | $\pi^3/6$ |

$n=2$が[[7shi-cla3]]の$2\pi$、$n=4$が[[7shi-cla4]]の$2\pi^2$です。
&&&

# 積分公式

以下、$M$は区分的に滑らかな境界を持つ有界な領域で、$M$とその境界を含む開集合$U$を取り、$F$は$U$上の$C^1$級関数とします。

## 小球面上の核

点$\boldsymbol a$を中心とする半径$\varepsilon$の球面の上では、$\boldsymbol x-\boldsymbol a=\varepsilon\boldsymbol n$（$\boldsymbol n$は外向きの単位法線）です。

&&&fml 球面の上の核と法線
$$
E(\boldsymbol x-\boldsymbol a)\,\boldsymbol n=\boldsymbol n\,E(\boldsymbol x-\boldsymbol a)=\frac{\varepsilon\boldsymbol n^2}{\varepsilon^n}=\frac1{\varepsilon^{n-1}}
$$
&&&

面積要素は$dS=\varepsilon^{n-1}d\Omega$なので

$$
\oint_{|\boldsymbol x-\boldsymbol a|=\varepsilon}E(\boldsymbol x-\boldsymbol a)\,\boldsymbol n\,dS=|S^{n-1}|
$$

であり、半径によりません。核が距離の$n-1$乗に反比例して減衰する分と、球面の面積が半径の$n-1$乗に比例して伸びる分が打ち消し合い、単位球面の面積だけが残ります。

## コーシー＝ポンペイウの公式

証明には、[[7shi-cla2]]の両側形式を使います。$C^1$級の関数$G,F$に対して

$$
\int_M\bigl((GD)F+G(DF)\bigr)dV=\oint_{\partial M}G\,\boldsymbol nF\,dS
$$

です。

&&&thm コーシー＝ポンペイウの公式 [thm-pompeiu]
$\boldsymbol a$が$M$の内部の点なら

$$
F(\boldsymbol a)=\frac1{|S^{n-1}|}\left(\oint_{\partial M}E(\boldsymbol x-\boldsymbol a)\,\boldsymbol nF\,dS-\int_ME(\boldsymbol x-\boldsymbol a)\,DF\,dV\right)
$$
&&&

&&&prf
[[7shi-cla3]]・[[7shi-cla4]]と同じ手順による。$G=E(\boldsymbol x-\boldsymbol a)$とし、$M$から中心$\boldsymbol a$、半径$\varepsilon$の閉球を除いた領域$M_\varepsilon$に[[7shi-cla2]]の両側形式を当てはめる。$M_\varepsilon$上では$GD=0$なので

$$
\int_{M_\varepsilon}G(DF)\,dV=\oint_{\partial M}G\boldsymbol nF\,dS-\oint_{|\boldsymbol x-\boldsymbol a|=\varepsilon}G\boldsymbol nF\,dS
$$

である。小球面の項では、$M_\varepsilon$から見た外向きの法線が$-\boldsymbol n$であることを符号に取り込んだ。小球面の上では$G\boldsymbol n=1/\varepsilon^{n-1}$、$dS=\varepsilon^{n-1}d\Omega$なので

$$
\oint_{|\boldsymbol x-\boldsymbol a|=\varepsilon}G\boldsymbol nF\,dS=\oint_{S^{n-1}}F(\boldsymbol a+\varepsilon\boldsymbol n)\,d\Omega
\ \longrightarrow\ |S^{n-1}|\,F(\boldsymbol a)\qquad(\varepsilon\to0)
$$

となる。左辺では$|G|=1/\rho^{n-1}$（$\rho=|\boldsymbol x-\boldsymbol a|$）、$dV=\rho^{n-1}d\rho\,d\Omega$より、被積分関数$G(DF)\rho^{n-1}$は$\boldsymbol a$の近くで有界であり、$\varepsilon\to0$で$M$上の積分に収束する。以上を整理すれば公式を得る。$M_\varepsilon$への両側形式の適用は、[[7shi-cla2]]と同じく概略である。
&&&

&&&cor コーシーの積分公式 [cor-cauchy]
$M$上で$DF=0$なら、$M$の内部の点$\boldsymbol a$で

$$
F(\boldsymbol a)=\frac1{|S^{n-1}|}\oint_{\partial M}E(\boldsymbol x-\boldsymbol a)\,\boldsymbol nF\,dS
$$
&&&

[[7shi-cla3]]・[[7shi-cla4]]の公式と式の形は同じで、核と定数が次元に応じて替わるだけです。

&&&rem 右モノジェニックな関数と基本解
両側形式で$G$の位置に$F$を、$F$の位置に$E(\boldsymbol x-\boldsymbol a)$を置くと、右からの作用に対する公式

$$
F(\boldsymbol a)=\frac1{|S^{n-1}|}\left(\oint_{\partial M}F\boldsymbol n\,E(\boldsymbol x-\boldsymbol a)\,dS-\int_M(FD)\,E(\boldsymbol x-\boldsymbol a)\,dV\right)
$$

が得られます。核が左右ともにモノジェニックであるため、どちらの側にも同じ核が使えます。[[7shi-cla3]]と同じく、台がコンパクトな関数では境界の項が消え、$E/|S^{n-1}|$が$D$の基本解であることを表します。
&&&

## 平均値の性質

$M$を$\boldsymbol a$を中心とする半径$R$の球とすると、境界の上で$E(\boldsymbol x-\boldsymbol a)\boldsymbol n=1/R^{n-1}$が定数になります。

&&&cor 平均値の性質
$\boldsymbol a$を中心とする半径$R$の閉球を含む開集合上で$DF=0$なら

$$
F(\boldsymbol a)=\frac1{|S^{n-1}|R^{n-1}}\oint_{|\boldsymbol x-\boldsymbol a|=R}F\,dS=\frac1{|S^{n-1}|}\oint_{S^{n-1}}F(\boldsymbol a+R\boldsymbol n)\,d\Omega
$$
&&&

中心での値は、球面上の値の平均です。

# $2\pi i$の分解

[[7shi-cla3]]では、コーシーの積分公式を複素数で書き直したときに現れる$2\pi i$を、出所の異なる2つの因子$2\pi$と$I$に分けました。$2\pi$は小円の上で核と法線の積を積分した$S^1$の長さであり、$I$は法線付きの境界要素と有向線素を結ぶ関係$\boldsymbol n\,ds=I\,d\boldsymbol x$から来ていました。一般の次元では、前者は$|S^{n-1}|$に置き換わりました。後者の置き換えを調べます。

## 有向境界要素

2次元の境界は曲線なので、向きを持つ要素は線素$d\boldsymbol x$（ベクトル）でした。$n$次元の領域の境界は$n-1$次元の超曲面なので、向きを持つ要素は$(n-1)$ベクトルになります。[[7shi-cla2]]の3次元の有向面素$d\boldsymbol X=\boldsymbol r_u\wedge\boldsymbol r_v\,du\,dv$を一般化します。

&&&def 有向境界要素
境界を$\boldsymbol r(u_1,\dots,u_{n-1})$とパラメーター表示し

$$
d\boldsymbol X=\boldsymbol r_{u_1}\wedge\dots\wedge\boldsymbol r_{u_{n-1}}\,du_1\cdots du_{n-1}
$$

とします。パラメーターの順序は、外向きの法線を先頭に置いた$\boldsymbol n\wedge d\boldsymbol X$が擬スカラー$I=e_0e_1\cdots e_{n-1}$の正の倍になるように選びます。
&&&

2次元では、境界を領域を左手に見る向きにたどるという[[7shi-cla2]]の規約がこれにあたります。たとえば単位円の点$e_0$では$\boldsymbol n=e_0$、左手に領域を見る向きの接線は$e_1$で、$\boldsymbol n\wedge e_1=I$です。

&&&fml 法線と有向境界要素
$$
\boldsymbol n\,dS=d\boldsymbol X\,I^{-1}
$$
&&&

&&&prf
$d\boldsymbol X$の因子はすべて接ベクトルなので$\boldsymbol n$と直交し、$\boldsymbol nd\boldsymbol X=\boldsymbol n\wedge d\boldsymbol X$である。$|\boldsymbol n|=1$と直交性より$|\boldsymbol n\wedge d\boldsymbol X|=|d\boldsymbol X|=dS$であり、向きの約束から$\boldsymbol n\,d\boldsymbol X=I\,dS$である。左から$\boldsymbol n$を掛けて$d\boldsymbol X=\boldsymbol nI\,dS$、右から$I^{-1}$を掛けて公式を得る。
&&&

$n=2$では$I^{-1}=-I$で、ベクトルは$I$と反交換するので、$d\boldsymbol x\,I^{-1}=I\,d\boldsymbol x$となり、[[7shi-cla2]]の$\boldsymbol n\,ds=I\,d\boldsymbol x$に一致します。$n=3$では$I$がすべての元と可換なので、$d\boldsymbol X=I\boldsymbol n\,dS$となり、[[7shi-cla2]]の有向面素の表示に一致します。

## 積分公式の書き換え

[積分公式](#cor-cauchy)を有向境界要素で書くと

$$
F(\boldsymbol a)=\frac1{|S^{n-1}|}\oint_{\partial M}E(\boldsymbol x-\boldsymbol a)\,d\boldsymbol X\,I^{-1}F
$$

です。$2\pi i$の2つの因子は、次のように置き換わります。

| 因子 | $n=2$ | 一般の$n$ | 出所 |
|---|---|---|---|
| 定数 | $2\pi$ | $\lvert S^{n-1}\rvert$ | 核と法線の積を小球面上で積分した値 |
| 擬スカラー | $I$ | $I=e_0\cdots e_{n-1}$ | 有向境界要素を法線に移す$\boldsymbol n\,dS=d\boldsymbol X\,I^{-1}$ |

$I^{-1}$は、どの次元でも積分の外に出せます。$n$が奇数なら$I$はすべての元と可換です。$n$が偶数なら$I$は偶部分の元と可換であり、$E$（ベクトル）と$d\boldsymbol X$（$(n-1)$ベクトル）はどちらも奇数グレードなので、積$E\,d\boldsymbol X$は偶部分に入り$I^{-1}$と可換です。したがって$E\,d\boldsymbol X\,I^{-1}F=I^{-1}E\,d\boldsymbol X\,F$であり

$$
F(\boldsymbol a)=\frac1{|S^{n-1}|\,I}\oint_{\partial M}E(\boldsymbol x-\boldsymbol a)\,d\boldsymbol X\,F
$$

と書けます。$\frac1{|S^{n-1}|\,I}=I^{-1}/|S^{n-1}|$は積分に左から掛けます。$n=3$では$4\pi I$が$2\pi i$の位置に来ます。

2次元に特有なのは、この先の書き換えです。[[7shi-cla3]]では、$d\boldsymbol x$がベクトルであるため$e_1$を掛けて$dz$に移し、偶部分全体が可換な$\mathbb C$であることを使って核・$dz$・$F$を並べ替え、通常の複素線積分$\frac1{2\pi i}\oint\frac{F\,dz}{z-w}$の形に戻しました。境界要素がベクトルの線素になることと、偶部分が可換であることの2つがそろうのは2次元だけです。

# 留数

[[7shi-cla3]]の留数の定義は、そのまま一般の次元に持ち上がります。点$\boldsymbol a$の近く（$\boldsymbol a$自身を除く）でモノジェニックな$F$について、$\boldsymbol a$を中心とする2つの球面に挟まれた球殻に基本定理（$G=1$の両側形式）を当てはめると、小球面上の$\oint\boldsymbol nF\,dS$は半径によりません（球殻への適用は[[7shi-cla2]]と同じく概略です）。

&&&def 留数
$F$が$0<|\boldsymbol x-\boldsymbol a|<\rho$でモノジェニックなとき、$0<\varepsilon<\rho$に対して

$$
\operatorname{Res}_{\boldsymbol a}F=\frac1{|S^{n-1}|}\oint_{|\boldsymbol x-\boldsymbol a|=\varepsilon}\boldsymbol nF\,dS
$$

を$\boldsymbol a$における$F$の**留数**と定めます。
&&&

留数定理も[[7shi-cla3]]と同じ形で成り立ちます。$M$とその境界を含む開集合$U$と、$M$の内部の相異なる有限個の点$\boldsymbol a_1,\dots,\boldsymbol a_m$を取り、$F$が$U\setminus\{\boldsymbol a_1,\dots,\boldsymbol a_m\}$上の$C^1$級関数で、$M\setminus\{\boldsymbol a_1,\dots,\boldsymbol a_m\}$上でモノジェニックなら、$\oint_{\partial M}\boldsymbol nF\,dS=|S^{n-1}|\sum_k\operatorname{Res}_{\boldsymbol a_k}F$です。特異点の周りの小球を除いた領域に基本定理を当てはめれば得られます。

&&&prop 核とその導関数の留数
$c\in\operatorname{Cl}_{n,0}(\mathbb R)$を定数とし、$\partial^\alpha$を$\partial_0,\dots,\partial_{n-1}$を1回以上合成した偏微分とすると

$$
\operatorname{Res}_{\boldsymbol a}\bigl(E(\boldsymbol x-\boldsymbol a)\,c\bigr)=c,\qquad
\operatorname{Res}_{\boldsymbol a}\bigl(\partial^\alpha E(\boldsymbol x-\boldsymbol a)\,c\bigr)=0
$$
&&&

&&&prf
$D$は右の定数に作用せず、偏微分と可換なので、どちらもモノジェニックである。第1式は小球面の上で$\boldsymbol nE(\boldsymbol x-\boldsymbol a)=1/\varepsilon^{n-1}$、$dS=\varepsilon^{n-1}d\Omega$から従う。第2式では、$\boldsymbol a=0$として$K=\partial^\alpha E$が$m\ge1$回の微分なら$-(n-1)-m$次の同次関数であるから

$$
\oint_{|\boldsymbol x|=\varepsilon}\boldsymbol nKc\,dS=\varepsilon^{-m}\oint_{|\boldsymbol x|=1}\boldsymbol nKc\,dS
$$

となる。左辺は$\varepsilon$によらないので、右辺の積分は$0$でなければならない。
&&&

したがって[[7shi-cla3]]と同じく、$F$が$\boldsymbol a$の近くで核とその導関数の有限和に、$\boldsymbol a$でもモノジェニックな関数を加えた形で書ける場合には、留数は核$E(\boldsymbol x-\boldsymbol a)$そのものの右側の係数として読み取れます。変わるのは核の形と、定義の分母の$|S^{n-1}|$だけです。

# フューター＝ソーの定理

[[7shi-cla4]]のフューターの定理は、複素解析の正則関数$f$から4次元の正則関数$\Delta\tilde f$を作るものでした。一般の次元でも同じ構成を考えます。軸対称性を使うので、[[7shi-cla4]]と同じく$e_0$を掛けた変数に移ります。

## パラベクトル変数

&&&def パラベクトル変数
$h_l=e_0e_l$（$l=1,\dots,n-1$）とし

$$
q=e_0\boldsymbol x=x_0+\sum_{l=1}^{n-1}x_lh_l,\qquad\mathcal D=e_0D=\partial_0+\sum_{l=1}^{n-1}h_l\partial_l
$$

とし、右からの作用を$F\mathcal D=\partial_0F+\sum_l(\partial_lF)h_l$とします。[[7shi-cla4]]と同じく、$\mathcal DF=0$を満たす$F$を**左正則**、$F\mathcal D=0$を満たす$F$を**右正則**と呼びます。
&&&

$h_l^2=-1$で、$h_l$どうしは反交換します。[[7shi-cla4]]と同じく、これは[[7shi-clif1]]の同型$\operatorname{Cl}_{n,0}^0(\mathbb R)\cong\operatorname{Cl}_{0,n-1}(\mathbb R)$で$h_l$が生成元になる構成です。スカラーと$h_l$の1次結合$q$は**パラベクトル**と呼ばれます。$n=2$では$q$は複素変数$z$に対応し、$n=4$では[[7shi-cla4]]と同じく$P_+$の成分に射影すると四元数の変数になります。

$e_0$は可逆なので、左正則性$\mathcal DF=0$は左モノジェニック性$DF=0$と同値です。右からの作用では$(Fe_0)D=F\mathcal D$なので、右正則性$F\mathcal D=0$は$F$ではなく$Fe_0$の右モノジェニック性にあたり、$FD=0$とは異なります。$\mathcal D$の係数は偶部分に属するので、偶部分に値を取る関数は偶部分に移ります。パラベクトル変数の理論は、同じ次元のベクトル変数の理論を、左は同じ方程式、右は値の対応$F\mapsto Fe_0$で書き直したものです。

## 軸対称な関数

$q=x_0+\boldsymbol q$（$\boldsymbol q=\sum_lx_lh_l$）と分け、$r=|\boldsymbol q|$、$\omega=\boldsymbol q/r$とします。$\omega^2=-1$です。[[7shi-cla4]]と同じ条件（複素共役で不変な開集合上の正則関数で$f(\bar z)=\overline{f(z)}$）を満たす$f=u+iv$に対して、$\tilde f=u(x_0,r)+\omega\,v(x_0,r)$と置きます。定義域は$\Omega_U=\{x_0+\boldsymbol q:x_0+i|\boldsymbol q|\in U\}$（$U$は$f$の定義域）です。実軸の上（$r=0$）では、[[7shi-cla4]]と同じく$u$が第2変数について偶、$v$が奇であることから、$\tilde f$を$r^2$の関数として滑らかに延長します。$f(z)=z^k$なら$\tilde f=q^k$、$f(z)=1/z$なら$\tilde f=q^{-1}$です。

&&&fml 軸対称な関数の微分 [fml-axial]
$A,B$を$x_0,r$の実数値$C^1$級関数とすると、$r>0$で

$$
\mathcal D(A+\omega B)=(A+\omega B)\mathcal D
=\partial_0A-\partial_rB-\frac{n-2}rB+\omega\,(\partial_rA+\partial_0B)
$$
&&&

&&&prf
[[7shi-cla4]]の証明で、虚数単位の個数を$3$から$n-1$に替える。$\varphi=B/r$とおくと$\omega B=\boldsymbol q\varphi$であり

$$
\sum_{l=1}^{n-1}h_l\partial_l(\boldsymbol q\varphi)
=\sum_lh_l^2\,\varphi+\frac{\boldsymbol q^2}r\,\partial_r\varphi
=-(n-1)\varphi-r\,\partial_r\varphi
$$

となる。$r\,\partial_r\varphi=\partial_rB-B/r$を代入すると$-\partial_rB-(n-2)B/r$である。残りの項と右からの作用は[[7shi-cla4]]と同じである。
&&&

コーシー＝リーマンの方程式により、$\tilde f$に当てはめると

$$
\mathcal D\tilde f=-(n-2)\frac vr
$$

が残ります。係数$n-2$は、$\omega$に直交する虚数単位の個数$(n-1)-1$です。$f(z)=z$では$v=r$なので$\mathcal Dq=2-n$であり、冒頭で引いた$D(e_0\boldsymbol x)=(2-n)e_0$に一致します。恒等写像の正則性を壊す係数と、フューターの定理のずれの係数は、同じものです。

## 定理

&&&thm フューター＝ソーの定理
$n$が偶数なら、上の条件を満たす正則関数$f$に対して、$\Delta^{(n-2)/2}\tilde f$は左右ともに正則です。
&&&

$n=4$では[[7shi-cla4]]のフューターの定理、$n=2$では$\tilde f=f$そのものです。一般の偶数$n$での証明は本記事では扱いません。

$n=4$の証明が短く済んだのは、$v$の2変数での調和性だけから、ずれ$v/r$が4次元の調和関数になったためです。$x_0,r$だけの関数$g$に対して$n$次元のラプラシアンは$\Delta g=\partial_0^2g+\partial_r^2g+\frac{n-2}r\partial_rg$と書けるので

$$
\Delta\frac vr=\frac{\partial_0^2v+\partial_r^2v}r+(n-4)\left(\frac{\partial_rv}{r^2}-\frac v{r^3}\right)
$$

となり、$v$が調和でも第2項は一般には残ります。$v$の調和性だけから、一般の$f$に対して$v/r$の調和性が従うのは$n=4$の場合です。$n\ge6$の偶数では、$\Delta$を1回作用させてもずれは一般には消えません（$f(z)=z$の$v/r=1$のように、個々の関数では消えることもあります）。一般の$f$に対してずれを消すには、定理のとおり$\Delta^{(n-2)/2}$を用います。その証明には、軸対称な関数に$\Delta$を反復したときの形を追う別の議論が要ります。

&&&rem 奇数次元
$n$が奇数なら指数$(n-2)/2$は半整数になります。適切な関数・超関数のクラスと定義域の条件のもとで、分数ラプラシアン$(-\Delta)^{(n-2)/2}$を用いる拡張が知られています。これは非局所的な作用素であり、定義には追加の準備が必要なため、本記事では扱いません。
&&&

## $1/z$から核へ

$f(z)=1/z$の場合は、定理の結論を直接確かめられます。

&&&fml 逆元の反復ラプラシアン
$n=2k+2$のとき

$$
\Delta^kq^{-1}=(-4)^k(k!)^2\,\frac{\bar q}{|q|^n}
$$

ここで$\bar q=\boldsymbol xe_0$であり、$\bar q/|q|^n=E(\boldsymbol x)\,e_0$は左右ともに正則です。
&&&

&&&prf
$q^{-1}=\boldsymbol x^{-1}e_0=\bar q/|q|^2$であり、$\bar q$の成分は座標の1次式である。$\rho=|\boldsymbol x|^2$の関数$g(\rho)$に対して$\partial_ag=2x_ag'$、$\Delta g=4\rho g''+2ng'$なので、座標$x_j$について

$$
\Delta\bigl(x_jg\bigr)=x_j\Delta g+2\partial_jg=x_j\bigl(4\rho g''+(2n+4)g'\bigr)
$$

である。$g=\rho^p$なら$4\rho g''+(2n+4)g'=4p\,(p+n/2)\,\rho^{p-1}$となる。$p=-1$から始めて$k$回繰り返すと、$j=1,\dots,k$について$p=-j$で係数$4(-j)(n/2-j)=-4j(k+1-j)$が掛かり、その積は$(-4)^k(k!)^2$である。$\rho^{-1-k}=|q|^{-n}$なので公式を得る。$\bar q/|q|^n=(\boldsymbol x/|\boldsymbol x|^n)e_0$は[コーシー核](#prop-kernel)と右の$e_0$の積であり、$\mathcal D(Ee_0)=e_0(DE)e_0=0$、$(Ee_0)\mathcal D=ED=0$である。
&&&

$n=4$（$k=1$）では$\Delta q^{-1}=-4\bar q/|q|^4$となり、[[7shi-cla4]]の結果に一致します。$n=6$では$\Delta^2q^{-1}=64\,\bar q/|q|^6$です。次数で見ると、$q^{-1}$の$-1$次から$\Delta$を$k$回作用させて$-1-2k=-(n-1)$次になり、球面上の積分から要求された核の次数にちょうど届きます。ラプラシアンの反復で逆元から核を作る場合、その回数$(n-2)/2$は、$-1$次の逆元と$-(n-1)$次の核の次数の差の半分として決まります。一般の$\tilde f$が正則になることは次数だけからは従わず、上の定理によります。

# まとめ

一般の次元$n$では、コーシーの積分公式の核は$E=\boldsymbol x/|\boldsymbol x|^n$、定数は単位球面の面積$|S^{n-1}|=2\pi^{n/2}/\Gamma(n/2)$です。核は球面の面積の伸びと打ち消し合う$-(n-1)$次の同次関数であり、ラプラシアンの基本解を$D$で微分したものです。積分公式、平均値の性質、留数の定義と核の導関数の留数が$0$であることは、[[7shi-cla3]]・[[7shi-cla4]]と同じ形で成り立ちます。

| 項目 | $n=2$ | $n=4$ | 一般の$n$ |
|---|---|---|---|
| 核 | $\boldsymbol x^{-1}$ | $\boldsymbol x/\lvert\boldsymbol x\rvert^4$ | $\boldsymbol x/\lvert\boldsymbol x\rvert^n$ |
| 核の次数 | $-1$ | $-3$ | $-(n-1)$ |
| ポテンシャル | $\log\lvert\boldsymbol x\rvert$ | $\lvert\boldsymbol x\rvert^{-2}$ | $\lvert\boldsymbol x\rvert^{2-n}$（$n\ge3$） |
| 定数 | $2\pi$ | $2\pi^2$ | $2\pi^{n/2}/\Gamma(n/2)$ |
| 境界要素と法線 | $\boldsymbol n\,ds=I\,d\boldsymbol x$ | $\boldsymbol n\,dS=d\boldsymbol X\,I^{-1}$ | $\boldsymbol n\,dS=d\boldsymbol X\,I^{-1}$ |
| $\mathcal Dq$ | $0$ | $-2$ | $2-n$ |
| ずれを消す作用素 | 不要 | $\Delta$ | $\Delta^{(n-2)/2}$（$n$が偶数） |

複素解析の$2\pi i$は、$|S^{n-1}|$と擬スカラー$I$に分かれ、どの次元でも$\frac1{|S^{n-1}|\,I}$の形で積分の外に出せます。これが通常の複素線積分$\frac1{2\pi i}\oint\frac{F\,dz}{z-w}$に戻るのは、境界要素がベクトルの線素になり、偶部分が可換な$\mathbb C$になる2次元の場合です。

恒等写像と逆元が正則でないこと（$D(e_0\boldsymbol x)=(2-n)e_0$、$D\boldsymbol x^{-1}=(n-2)/|\boldsymbol x|^2$）と、フューターの定理のずれ$-(n-2)v/r$は、同じ係数$n-2$で表されます。いずれも$n=2$でだけ消え、2次元で恒等写像$z$と逆数$1/z$がそのまま正則だったことに対応します。ラプラシアンの反復で$1/z$から核を作るには、次数を$n-2$だけ下げる必要があり、その回数は$(n-2)/2$です。

球面モノジェニックス（球面調和関数のモノジェニック版）による展開、$n$が奇数の場合のフューター＝ソーの定理は、本記事では扱いません。
