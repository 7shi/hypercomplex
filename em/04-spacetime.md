[[7shi-em2]]と[[7shi-em3]]では、3次元の$\operatorname{Cl}_{3,0}(\mathbb R)$に時間の作用素$\partial_0$を加えて、マクスウェル方程式を$\mathcal DF=(\rho-\boldsymbol J/c)/\varepsilon_0$と書きました。本記事では時間と空間を1つのベクトル空間にまとめた時空代数$\operatorname{Cl}_{1,3}(\mathbb R)$に移ります。ディラック作用素$D$は$D^2=\partial_0^2-\Delta$を満たし、電磁場$F$は時空の2ベクトル、電流$J$は時空のベクトルになって、マクスウェル方程式は$DF=\mu_0cJ$と書けます。時間方向の生成元$\gamma_0$を掛けると、[[7shi-em2]]の式に戻ります。この操作は、[[7shi-cla4]]で$e_0$を掛けて四元数の変数と作用素を作った操作と、基準方向の生成元を掛けて偶部分に移すという構造を共有しています。違いは空間の生成元の2乗の符号です。空間の生成元の2乗が$+1$なら、偶部分の生成元の2乗は$-1$となり、楕円型の四元数解析が現れます。空間の生成元の2乗が$-1$なら、偶部分の生成元の2乗は$+1$となり、双曲型の電磁気学が現れます。

# 時空代数

## 生成元と計量

$\operatorname{Cl}_{1,3}(\mathbb R)$の生成元を$\gamma_0,\gamma_1,\gamma_2,\gamma_3$とします。

&&&def 時空代数
$$
\gamma_0^2=1,\qquad\gamma_k^2=-1\ (k=1,2,3),\qquad\mu\ne\nu\text{ なら }\gamma_\mu\gamma_\nu=-\gamma_\nu\gamma_\mu
$$

で生成される実クリフォード代数$\operatorname{Cl}_{1,3}(\mathbb R)$を**時空代数**と呼びます。
&&&

時空の点を$x_0=ct$と空間座標$x_1,x_2,x_3$で

$$
x=x_0\gamma_0+\sum_{k=1}^3x_k\gamma_k
$$

と書くと、$x^2=x_0^2-x_1^2-x_2^2-x_3^2$です。ベクトルの2乗が符号数$(+,-,-,-)$のミンコフスキー計量を与えます。

## ディラック作用素

[[7shi-cla1]]の$\operatorname{Cl}_{n,0}(\mathbb R)$では、作用素の係数に生成元をそのまま使いました。一般には、微分作用素の係数には座標基底$\gamma_\mu$に対する相反基底$\gamma^\mu$（$\gamma^\mu\cdot\gamma_\nu=\delta^\mu{}_\nu$）を使います。本記事の正規直交基底では、相反基底は各生成元の逆元$\gamma^\mu=\gamma_\mu^{-1}$に一致します。

&&&def 時空のディラック作用素
$$
\gamma^0=\gamma_0,\qquad\gamma^k=-\gamma_k,\qquad
D=\sum_{\mu=0}^3\gamma^\mu\partial_\mu=\gamma_0\partial_0-\sum_{k=1}^3\gamma_k\partial_k
$$
&&&

このため$\gamma^\mu\gamma_\mu=1$から$Dx=\sum_\mu\gamma^\mu\gamma_\mu=4$となり、[[7shi-cla5]]の$D\boldsymbol x=n$もそのまま保たれます。$\operatorname{Cl}_{n,0}(\mathbb R)$では生成元がそれ自身の逆元なので、相反基底は生成元に一致していました。$D$の2乗は、[[7shi-cla1]]と同じ計算で係数の2乗の和になります。

&&&fml 波動作用素
$$
D^2=\sum_\mu(\gamma^\mu)^2\partial_\mu^2=\partial_0^2-\Delta
$$
&&&

[[7shi-em2]]では共役との積$\bar{\mathcal D}\mathcal D$として得た波動作用素が、ここでは$D$そのものの2乗として現れます。$\operatorname{Cl}_{n,0}(\mathbb R)$の$D^2=\Delta$と同じく、$D$は$\partial_0^2-\Delta$の平方根です。

# 時空分割

## 相対ベクトル

時間方向の$\gamma_0$を基準にして、時空の量を$\operatorname{Cl}_{3,0}(\mathbb R)$の量に分けます。

&&&def 相対ベクトル
$$
\sigma_k=\gamma_k\gamma_0\qquad(k=1,2,3)
$$
&&&

&&&prop 偶部分代数と$\operatorname{Cl}_{3,0}(\mathbb R)$
$\sigma_k^2=1$、$k\ne l$なら$\sigma_k\sigma_l=-\sigma_l\sigma_k$で、$\sigma_1\sigma_2\sigma_3=\gamma_0\gamma_1\gamma_2\gamma_3$です。$1,\sigma_k,\sigma_k\sigma_l\ (k<l),\sigma_1\sigma_2\sigma_3$は偶部分代数$\operatorname{Cl}_{1,3}^0(\mathbb R)$の基底をなし、$\operatorname{Cl}_{1,3}^0(\mathbb R)\cong\operatorname{Cl}_{3,0}(\mathbb R)$です。
&&&

&&&prf
$\sigma_k^2=\gamma_k\gamma_0\gamma_k\gamma_0=-\gamma_k^2\gamma_0^2=1$である。$k\ne l$なら$\gamma_k\gamma_0\gamma_l\gamma_0=-\gamma_k\gamma_l$であり、これは$k,l$について反対称だから$\sigma_k\sigma_l=-\sigma_l\sigma_k$となる。$\sigma_1\sigma_2\sigma_3=(-\gamma_1\gamma_2)\gamma_3\gamma_0=-\gamma_1\gamma_2\gamma_3\gamma_0=\gamma_0\gamma_1\gamma_2\gamma_3$である（$\gamma_0$を3つの生成元の前に移すと符号が3回変わる）。$\sigma_k$は$\operatorname{Cl}_{3,0}(\mathbb R)$の生成元の関係を満たし、それらの積の8個は偶部分代数の8次元の基底$1,\gamma_\mu\gamma_\nu\ (\mu<\nu),\gamma_0\gamma_1\gamma_2\gamma_3$と、符号を除いて一致する。したがって、$e_k\mapsto\sigma_k$で定まる代数準同型$\operatorname{Cl}_{3,0}(\mathbb R)\to\operatorname{Cl}_{1,3}^0(\mathbb R)$は基底を基底に移すので、同型である。
&&&

$\operatorname{Cl}_{3,0}(\mathbb R)$の擬スカラー$I=\sigma_1\sigma_2\sigma_3$は、そのまま時空代数の擬スカラー$\gamma_0\gamma_1\gamma_2\gamma_3$です。$I^2=-1$で、$I$は時空のベクトルとは反可換、偶部分代数の元とは可換です。[[7shi-em1]]から[[7shi-em3]]で使った$\operatorname{Cl}_{3,0}(\mathbb R)$の計算は、すべて偶部分代数の中の計算として読み直せます。以下、$\operatorname{Cl}_{3,0}(\mathbb R)$のベクトルは$\sigma_k$で書き、$\boldsymbol x=\sum_kx_k\sigma_k$などとします。

## 位置と作用素の分割

時空のベクトルに$\gamma_0$を掛けると、偶部分代数のパラベクトルになります。

&&&fml 位置の分割
$$
x\gamma_0=x_0+\boldsymbol x,\qquad\gamma_0x=x_0-\boldsymbol x,\qquad
x^2=x\gamma_0\,\gamma_0x=x_0^2-|\boldsymbol x|^2
$$
&&&

ミンコフスキー計量は、パラベクトル$p=x_0+\boldsymbol x$とその共役$\bar p=x_0-\boldsymbol x$の積$p\bar p$として現れます。

&&&fml 作用素の分割 [fml-split]
$$
\gamma_0D=\partial_0+\sum_k\sigma_k\partial_k=\mathcal D,\qquad
D(\gamma_0H)=\partial_0H-\sum_k\sigma_k\partial_kH=\bar{\mathcal D}H
$$
&&&

&&&prf
$\gamma_0\gamma^k=-\gamma_0\gamma_k=\gamma_k\gamma_0=\sigma_k$より第1式を得る。第2式では$\gamma^k\gamma_0=-\gamma_k\gamma_0=-\sigma_k$である。
&&&

[[7shi-em2]]の$\mathcal D$と$\bar{\mathcal D}$は、時空のディラック作用素に$\gamma_0$を左から掛けたものと、右から掛けたものです。したがって$\bar{\mathcal D}\mathcal DH=D(\gamma_0\gamma_0DH)=D^2H$であり、[[7shi-em2]]の$\bar{\mathcal D}\mathcal D=\partial_0^2-\Delta$は$D^2=\partial_0^2-\Delta$の言い換えです。

# 係数の2乗の符号

[[7shi-cla4]]では、$\operatorname{Cl}_{4,0}(\mathbb R)$で位置ベクトルに$e_0$を掛けて四元数の変数$q=e_0\boldsymbol x$を作り、フューター作用素を$\mathcal D=e_0D=\partial_0+\sum_lh_l\partial_l$（$h_l=e_0e_l$）としました。上の[作用素の分割](#fml-split)は、基準方向の生成元を作用素に左から掛けて偶部分に移すという同じ構造を、$\operatorname{Cl}_{1,3}(\mathbb R)$で実現したものです。積の順序が$h_l=e_0e_l$と$\sigma_k=\gamma_k\gamma_0$で逆に見えるのは、時空では相反基底が$\gamma^k=-\gamma_k$になるためです。四元数側の係数は$e_0e_l=h_l$、時空側の係数は$\gamma_0\gamma^k=-\gamma_0\gamma_k=\sigma_k$で、どちらも基準方向の生成元を左から掛けています。ただし四元数側では偶部分$\operatorname{Cl}_{4,0}^0(\mathbb R)\cong\mathbb H\oplus\mathbb H$の一方を選んでいたのに対し、時空側では偶部分の全体が$\operatorname{Cl}_{3,0}(\mathbb R)$です。

偶部分に移った生成元の2乗は、空間の生成元の2乗の符号で決まります。$e_0^2=\gamma_0^2=1$なので

$$
h_l^2=(e_0e_l)^2=-e_0^2e_l^2=-1,\qquad
\sigma_k^2=(\gamma_k\gamma_0)^2=-\gamma_k^2\gamma_0^2=+1
$$

です。次の表の$\Delta$は3次元の空間のラプラシアンです。

| | 代数 | 空間の生成元の2乗 | 偶部分の生成元 | その2乗 | $\bar{\mathcal D}\mathcal D$ | 型 |
|---|---|---|---|---|---|---|
| 四元数解析（[[7shi-cla4]]） | $\operatorname{Cl}_{4,0}(\mathbb R)$ | $e_l^2=+1$ | $h_l=e_0e_l$ | $-1$ | $\partial_0^2+\Delta$ | 楕円型 |
| 電磁気学 | $\operatorname{Cl}_{1,3}(\mathbb R)$ | $\gamma_k^2=-1$ | $\sigma_k=\gamma_k\gamma_0$ | $+1$ | $\partial_0^2-\Delta$ | 双曲型 |

どちらも「2乗が$+1$の基準方向の生成元を掛けて偶部分に移す」操作で、得られる作用素は$\partial_0+\sum(\text{偶部分の生成元})\partial$という同じ形をしています。この形の中での違いは空間の生成元の2乗の符号だけで、それが偶部分の生成元の2乗の符号を反転させ、ラプラシアンと波動作用素を分けます。[[7shi-em2]]では$\operatorname{Cl}_{3,0}(\mathbb R)$のまま、共役との積として波動作用素を得ました。時空全体のディラック作用素$D$そのものの2乗として波動作用素を得るには、ユークリッド型の符号数からミンコフスキー型の符号数に替える必要があります。

# マクスウェル方程式

## 場と電流

$\operatorname{Cl}_{3,0}(\mathbb R)$の$F=\boldsymbol E+Ic\boldsymbol B$を偶部分代数の元として読み直します。[[7shi-em2]]と同じく、場と源は必要な回数だけ連続微分可能とします。

&&&def 時空の電磁場と電流
$$
F=\sum_kE_k\sigma_k+Ic\sum_kB_k\sigma_k,\qquad
J=c\rho\,\gamma_0+\sum_kJ_k\gamma_k
$$
&&&

$\sigma_k=\gamma_k\gamma_0$も$I\sigma_k$も時空の2ベクトルなので、$F$は時空の2ベクトルです。$I\sigma_1=\gamma_3\gamma_2$、$I\sigma_2=\gamma_1\gamma_3$、$I\sigma_3=\gamma_2\gamma_1$であり、時空の2ベクトルの$\binom42=6$成分のうち、$\gamma_0$を含む面の3成分が電場、空間の面の3成分が磁場です。電場と磁場は1つの2ベクトルの成分で、その分け方は時間方向$\gamma_0$の選び方に依存します。

電流$J$は時空のベクトルです。$\gamma_0$を掛けると

$$
\gamma_0J=c\rho+\sum_kJ_k\gamma_0\gamma_k=c\Bigl(\rho-\frac{\boldsymbol J}c\Bigr)
$$

となり、[[7shi-em2]]の源が現れます。

## 1本の式

&&&thm 時空のマクスウェル方程式
$$
DF=\mu_0cJ
$$

は[[7shi-em2]]の$\mathcal DF=(\rho-\boldsymbol J/c)/\varepsilon_0$と同値です。
&&&

&&&prf
$\gamma_0$は可逆なので、$DF=\mu_0cJ$は$\gamma_0DF=\mu_0c\gamma_0J$と同値である。[作用素の分割](#fml-split)より左辺は$\mathcal DF$、右辺は$\mu_0c^2(\rho-\boldsymbol J/c)=(\rho-\boldsymbol J/c)/\varepsilon_0$である。
&&&

係数$\mu_0c$は$1/\varepsilon_0c$とも書けます。$J$の時間成分を$c\rho$としたのは、$J$の4成分の単位を揃えるためです。

## 内積と外積への分解

$D$はベクトルを係数とする作用素なので、2ベクトル$F$に作用させるとグレードが1つ下がる部分と1つ上がる部分に分かれます。[[7shi-cla1]]で見た分解$DF=D\cdot F+D\wedge F$です。

&&&fml 2本の式
$$
D\cdot F=\mu_0cJ,\qquad D\wedge F=0
$$
&&&

右辺$\mu_0cJ$はベクトルだけなので、$DF$のベクトル部$D\cdot F$とトリベクトル部$D\wedge F$を別々に比べられます。$\gamma_0$を掛けて$\operatorname{Cl}_{3,0}(\mathbb R)$に戻すと

$$
\gamma_0(D\cdot F)=\nabla\cdot\boldsymbol E+\bigl(\partial_0\boldsymbol E-c\,\nabla\times\boldsymbol B\bigr),\qquad
\gamma_0(D\wedge F)=I\bigl(\nabla\times\boldsymbol E+c\,\partial_0\boldsymbol B\bigr)+Ic\,\nabla\cdot\boldsymbol B
$$

であり、[[7shi-em2]]のスカラー部とベクトル部が$D\cdot F$から、2ベクトル部と擬スカラー部が$D\wedge F$から来ています。次の表の「スカラー」「ベクトル」などは、時空ではなく偶部分代数を$\operatorname{Cl}_{3,0}(\mathbb R)$と見たときのグレードです。

| 時空代数 | $\operatorname{Cl}_{3,0}(\mathbb R)$のグレード | 式 |
|---|---|---|
| $D\cdot F=\mu_0cJ$（ベクトル） | スカラー | ガウスの法則 |
| | ベクトル | アンペール＝マクスウェルの法則 |
| $D\wedge F=0$（トリベクトル） | 2ベクトル | ファラデーの法則 |
| | 擬スカラー | $\nabla\cdot\boldsymbol B=0$ |

時空代数では、源を持つ2本と源を持たない2本が、それぞれ1本の式にまとまります。

## 連続の式と波動方程式

電流の発散は$D\cdot J=\partial_0(c\rho)+\sum_k\partial_kJ_k=\partial_t\rho+\nabla\cdot\boldsymbol J$です。

&&&prop 連続の式
$D\cdot F=\mu_0cJ$の解について$D\cdot J=0$、すなわち$\partial_t\rho+\nabla\cdot\boldsymbol J=0$が成り立ちます。
&&&

&&&prf
2ベクトル$F$について$C_{\mu\nu}=\gamma^\mu\cdot(\gamma^\nu\cdot F)$と置くと、$D\cdot(D\cdot F)=\sum_{\mu,\nu}\partial_\mu\partial_\nu C_{\mu\nu}$である。ベクトル$a,b,u,v$について$a\cdot\bigl(b\cdot(u\wedge v)\bigr)=(b\cdot u)(a\cdot v)-(b\cdot v)(a\cdot u)$は$a,b$について反対称なので、線形性により$C_{\mu\nu}=-C_{\nu\mu}$である。偏微分は交換するので

$$
\sum_{\mu,\nu}\partial_\mu\partial_\nu C_{\mu\nu}=\frac12\sum_{\mu,\nu}\partial_\mu\partial_\nu\bigl(C_{\mu\nu}+C_{\nu\mu}\bigr)=0
$$

となる。
&&&

$D$をもう一度掛けると$D^2F=\mu_0c\,DJ=\mu_0c\,(D\cdot J+D\wedge J)=\mu_0c\,D\wedge J$です。左辺は$(\partial_0^2-\Delta)F$で、[[7shi-em2]]の電場と磁場の波動方程式がこの1本の2ベクトルの式にまとまります。

# 符号数の選択

時空代数には、本記事の$\operatorname{Cl}_{1,3}(\mathbb R)$（時間の2乗が$+1$）と、$\operatorname{Cl}_{3,1}(\mathbb R)$（空間の2乗が$+1$、時間の2乗が$-1$）の2つの流儀があります。

&&&rem $\operatorname{Cl}_{1,3}(\mathbb R)$を選ぶ理由
- **前半との接続**：$\operatorname{Cl}_{3,1}(\mathbb R)$では時間の生成元の2乗が$\gamma_0^2=-1$です。$\gamma_k\gamma_0$の2乗は$+1$のままで、偶部分代数はやはり$\operatorname{Cl}_{3,0}(\mathbb R)$と同型ですが、$\gamma_0\gamma_0=-1$のために$D(\gamma_0\gamma_0DH)=-D^2H$となり、$D^2=\Delta-\partial_0^2$と全体の符号が逆になります。相反基底が$\gamma^0=-\gamma_0$、$\gamma^k=\gamma_k$になるので、同じ$\sigma_k=\gamma_k\gamma_0$を使うと$\gamma_0D=\bar{\mathcal D}$、$D(\gamma_0H)=\mathcal DH$と対応も入れ替わります。[[7shi-cla4]]の$e_0$と同じく2乗が$+1$の生成元を掛けて偶部分に移すという操作をそのまま保てるのは、$\operatorname{Cl}_{1,3}(\mathbb R)$のほうです。
- **パラベクトルの計量**：[[7shi-em2]]の$\bar{\mathcal D}\mathcal D=\partial_0^2-\Delta$と、パラベクトルのノルム$p\bar p=x_0^2-|\boldsymbol x|^2$は、ともに$(+,-,-,-)$の符号です。
- **パウリ行列からの拡張**：下の例のように、パウリ行列からワイル表現のガンマ行列が符号の追加なしに得られます。
- **代数としての違い**：[[7shi-clif1]]の分類では$\operatorname{Cl}_{1,3}(\mathbb R)\cong M_2(\mathbb H)$、$\operatorname{Cl}_{3,1}(\mathbb R)\cong M_4(\mathbb R)$で、実代数としては同型ではありません。複素化すればともに$M_4(\mathbb C)$で、共通の複素化の中では、$\operatorname{Cl}_{1,3}(\mathbb R)$の生成元を複素スカラー$i$倍したものを$\operatorname{Cl}_{3,1}(\mathbb R)$の生成元として取れます。この$i$は擬スカラー$I$ではなく、すべての元と可換な複素数の単位です。[[7shi-lie4]]と[[7shi-ideal]]が$\operatorname{Cl}_{3,1}(\mathbb R)$を使うのは$\operatorname{SO}(4)$や実表現の文脈で、本シリーズとは目的が異なります。
&&&

&&&ex パウリ行列からの構成
パウリ行列による$\operatorname{Cl}_{3,0}(\mathbb R)$の表現（[[7shi-bq]]・[[7shi-lie3]]）で、パラベクトル$p$を

$$
\gamma(p)=\begin{pmatrix}0&p\\\bar p&0\end{pmatrix}
$$

と4次の行列に置きます。$\bar p$は複素共役やエルミート共役ではなく、パラベクトルの共役（ベクトル部の符号の反転）です。$p\bar p=\bar pp$はスカラーなので、$\gamma(p)^2=(p\bar p)I_4$です。$p=1,\sigma_k$の像はワイル表現のガンマ行列そのもので、$\gamma(1)^2=I_4$、$\gamma(\sigma_k)^2=-I_4$から符号数$(+,-,-,-)$が出ます。$\operatorname{Cl}_{3,1}(\mathbb R)$にするには、下のブロックを$-\bar p$に変える符号が1つ要ります。
&&&

$\operatorname{Cl}_{1,3}(\mathbb R)$の代償は、時空のベクトルとしての$\gamma_k$の2乗が$-1$で、[[7shi-em1]]から[[7shi-em3]]の$\boldsymbol x^2=|\boldsymbol x|^2$という規約から外れることです。空間のベクトルは相対ベクトル$\sigma_k$（2乗$+1$）として保たれるので、$\operatorname{Cl}_{3,0}(\mathbb R)$での計算はそのまま使えます。

# まとめ

時空代数$\operatorname{Cl}_{1,3}(\mathbb R)$で、マクスウェル方程式は$DF=\mu_0cJ$の1本の式になりました。

- **作用素**：逆元を係数にしたディラック作用素は$D^2=\partial_0^2-\Delta$を満たし、波動作用素の平方根です。
- **時空分割**：$\sigma_k=\gamma_k\gamma_0$が偶部分代数で$\operatorname{Cl}_{3,0}(\mathbb R)$の生成元になり、$\gamma_0D=\mathcal D$、$D(\gamma_0H)=\bar{\mathcal D}H$です。[[7shi-em2]]の式は$DF=\mu_0cJ$に$\gamma_0$を掛けたものです。
- **場と電流**：$F$は時空の2ベクトルで、$\gamma_0$を含む面の成分が電場、空間の面の成分が磁場です。$J=c\rho\gamma_0+\sum J_k\gamma_k$は時空のベクトルです。
- **2本の式**：$D\cdot F=\mu_0cJ$が源を持つ2本を、$D\wedge F=0$が源を持たない2本をまとめます。
- **係数の2乗の符号**：$\gamma_0$を掛ける操作は、[[7shi-cla4]]の$e_0$を掛ける操作と「基準方向の生成元を掛けて偶部分に移す」構造を共有します。空間の生成元の2乗の符号が偶部分の生成元の2乗の符号を反転させ、楕円型と双曲型を分けます。
