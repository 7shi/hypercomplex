時空代数の偶部分代数におけるディラック方程式の幾何学的定式化を扱い、共変性と平面波解を導出します。

# 概要

前回は、時空代数$\operatorname{Cl}_{1,3}(\mathbb R)$の偶部分代数の元としてディラックスピノル$\psi$を導入し、ローレンツ変換が左から回転子$\psi\mapsto R\psi$として作用することを見ました。[[7shi-dirac2]]

本記事では、$\psi$が従う方程式を扱います。相対論的なエネルギーと運動量の関係からクライン＝ゴルドン方程式を作り、その作用素$\square=\partial_0^2-\Delta$の平方根としてディラック作用素$D$を使います。[[7shi-em4]]

値を偶部分代数に取るため右からの因子が必要になり、ディラック方程式は$\hbar D\psi\,\omega\sigma_3=mc\,\psi\gamma_0$と書けます。行列形式$i\hbar\hat\gamma^\mu\partial_\mu\Psi=mc\Psi$との対応を示し、右から掛かる因子が左からのローレンツ変換と可換であることから、方程式の共変性が直ちに従うことを確かめます。最後に平面波解を求め、正と負のエネルギーの解を調べます。

# 1階の方程式

## クライン＝ゴルドン方程式

時空代数の記法に従い、$x_0=ct$、$D=\sum_\mu\gamma^\mu\partial_\mu$（$\gamma^0=\gamma_0$、$\gamma^k=-\gamma_k$）とし、$D^2=\partial_0^2-\Delta$を$\square$と書きます。擬スカラーは$\omega=\gamma_0\gamma_1\gamma_2\gamma_3$、相対ベクトルは$\sigma_k=\gamma_k\gamma_0$です。[[7shi-em4]]

質量$m$の自由粒子のエネルギー$E$と運動量$\boldsymbol p$は、特殊相対論では

$$
E^2=|\boldsymbol p|^2c^2+m^2c^4
$$

を満たします。固有速度$U$（$U^2=c^2$）を使えば、4元運動量$p=mU$について$p^2=m^2c^2$と書け、$p\cdot\gamma_0=E/c$がその観測者の見るエネルギー、残りの成分が運動量です。[[7shi-em5]]

シュレーディンガー方程式と同じく、$E\to i\hbar\,\partial_t$、$\boldsymbol p\to-i\hbar\nabla$と置き換えると

$$
\Bigl(\square+\frac{m^2c^2}{\hbar^2}\Bigr)\phi=0
$$

を得ます。これを**クライン＝ゴルドン方程式**と呼びます。[[7shi-clif5]]

この方程式は時間について2階なので、初期条件として$\phi$とその時間微分の両方を与える必要があります。また、$\phi^*\phi$のような正の量が保存される密度になりません。保存される量を作ると、$\phi$と時間微分の組み合わせになり、正にも負にもなり得ます。ディラックは、確率密度を正に保つため、時間について1階の方程式を求めました。

## 作用素の平方根

1階の作用素で、2乗すると$\square$になるものは、すでに手元にあります。ディラック作用素$D$です。ラプラシアンの平方根として$D$を組んだのと同じ発想で、係数の2乗の符号を替えた$\operatorname{Cl}_{1,3}(\mathbb R)$では$D^2=\square$となります。[[7shi-em4]][[7shi-cla1]][[7shi-em7]]

$D$をスピノルに作用させる方程式を作ります。$D$はベクトルを係数とする作用素なので、偶部分代数の元$\psi$に作用させると奇部分の元$D\psi$になります。質量の項も奇部分の元でなければならず、$\psi$の左からベクトルを掛けるとローレンツ変換と両立しないので（左側は$\psi\mapsto R\psi$の作用を受けます）、右から掛けます。複素数の虚数単位$i$の役割も、以前の記事と同じく右から掛かる$\omega\sigma_3$が担います。[[7shi-dirac1]]

&&&def ディラック方程式
$$
\hbar\,D\psi\,\omega\sigma_3=mc\,\psi\gamma_0
$$
&&&

右から掛かる$\omega\sigma_3=\gamma_2\gamma_1$と$\gamma_0$は可換です。この方程式を満たす$\psi$は、クライン＝ゴルドン方程式も満たします。

&&&prop クライン＝ゴルドン方程式との関係
ディラック方程式の解$\psi$は、$\bigl(\square+m^2c^2/\hbar^2\bigr)\psi=0$を満たします。
&&&

&&&prf
$(\omega\sigma_3)^{-1}=-\omega\sigma_3$より、方程式は$D\psi=-\frac{mc}\hbar\psi\gamma_0\,\omega\sigma_3$と同値である。$D$は左から作用するので

$$
D^2\psi=-\frac{mc}\hbar(D\psi)\gamma_0\,\omega\sigma_3=\frac{m^2c^2}{\hbar^2}\psi\gamma_0\,\omega\sigma_3\,\gamma_0\,\omega\sigma_3=\frac{m^2c^2}{\hbar^2}\psi\,(\omega\sigma_3)^2=-\frac{m^2c^2}{\hbar^2}\psi
$$

となる（$\gamma_0$と$\omega\sigma_3$は可換、$\gamma_0^2=1$）。
&&&

右からの因子$\gamma_0$と$\omega\sigma_3$が可換で、それぞれの2乗が$1$と$-1$であることが、2乗して符号が正しく合う理由です。

# 行列形式との対応

## ディラック表現

行列形式では、$4\times4$の複素行列$\hat\gamma^\mu$と4成分の列ベクトル$\Psi$で、ディラック方程式を

$$
i\hbar\sum_\mu\hat\gamma^\mu\partial_\mu\Psi=mc\Psi
$$

と書きます。ここでは標準的なディラック表現

$$
\hat\gamma^0=\begin{pmatrix}I&0\\0&-I\end{pmatrix},\qquad
\hat\gamma^k=\begin{pmatrix}0&\hat\sigma_k\\-\hat\sigma_k&0\end{pmatrix}
$$

を使います。$\hat\gamma_\mu$は添字を下げたもので、$\hat\gamma_0=\hat\gamma^0$、$\hat\gamma_k=-\hat\gamma^k$です。$\hat\gamma_\mu$は$\operatorname{Cl}_{1,3}(\mathbb R)$の生成元の関係を満たします。以前に紹介したパウリ行列からの構成はワイル表現と呼ばれる別の表現で、両者は$\mathbb C^4$の基底の取り替えで移り合います。また$\hat\gamma_5=i\hat\gamma^0\hat\gamma^1\hat\gamma^2\hat\gamma^3=\begin{pmatrix}0&I\\I&0\end{pmatrix}$と置きます。[[7shi-em4]]

## スピノルの対応

ディラックスピノル$\psi$を、$\gamma_0$と可換な部分と反可換な部分に分けます。$\gamma_0$で挟む操作$\psi\mapsto\gamma_0\psi\gamma_0$は、$\gamma_0$と可換な部分を変えず、反可換な部分の符号を変えます。前回見たとおり、$\gamma_0$と可換な偶部分代数の元は$1,\omega\sigma_k$が張るパウリスピノルの空間です。反可換な部分は、パウリスピノルの空間の元に右から$\sigma_3$を掛けて得られます（$\sigma_3$は$\gamma_0$と反可換）。[[7shi-dirac2]]

$$
\psi=\phi+\eta\,\sigma_3,\qquad\phi=\frac12(\psi+\gamma_0\psi\gamma_0),\qquad\eta\,\sigma_3=\frac12(\psi-\gamma_0\psi\gamma_0)
$$

$\phi$と$\eta$はパウリスピノルの空間の元です。パウリスピノルの対応をそれぞれに使って$\mathbb C^2$の元を作り、上下に並べます。以下、その対応を$\phi\mapsto|\phi\rangle$と書きます。[[7shi-dirac1]]

&&&def ディラックスピノルの対応
$$
\psi=\phi+\eta\,\sigma_3\quad\longmapsto\quad\Psi=\begin{pmatrix}|\phi\rangle\\|\eta\rangle\end{pmatrix}\in\mathbb C^4
$$
&&&

実8次元の偶部分代数と、実8次元の$\mathbb C^4$の間の実線形な全単射です。

&&&prop 作用の翻訳 [prop-dirac-action]
1. $\hat\gamma_\mu\Psi$には$\gamma_\mu\psi\gamma_0$が対応します（$\hat\gamma^\mu\Psi$には$\gamma^\mu\psi\gamma_0$）。
2. $i\Psi$には$\psi\,\omega\sigma_3$が対応します。
3. $\hat\gamma_5\Psi$には$\psi\sigma_3$が対応します。
&&&

&&&prf
パウリスピノルの空間における作用の翻訳を使う。すなわち$|\sigma_k\phi\sigma_3\rangle=\hat\sigma_k|\phi\rangle$、$|\phi\,\omega\sigma_3\rangle=i|\phi\rangle$である。[[7shi-dirac1]]

1：$\gamma_0\psi\gamma_0=\phi-\eta\sigma_3$は$(|\phi\rangle,-|\eta\rangle)^T=\hat\gamma_0\Psi$に対応する。$\gamma_k\psi\gamma_0=\gamma_k\gamma_0\,\gamma_0\psi\gamma_0=\sigma_k\phi-\sigma_k\eta\sigma_3$である。$\sigma_k\phi$は$\gamma_0$と反可換で、$\sigma_k\phi=(\sigma_k\phi\sigma_3)\sigma_3$と書けば下の成分$\hat\sigma_k|\phi\rangle$を与える。$-\sigma_k\eta\sigma_3$は$\gamma_0$と可換で、上の成分$-\hat\sigma_k|\eta\rangle$を与える。したがって$\gamma_k\psi\gamma_0$は$(-\hat\sigma_k|\eta\rangle,\hat\sigma_k|\phi\rangle)^T=\hat\gamma_k\Psi$に対応する。

2：$\omega\sigma_3$は$\sigma_3$と可換だから$\psi\,\omega\sigma_3=\phi\,\omega\sigma_3+(\eta\,\omega\sigma_3)\sigma_3$であり、$(i|\phi\rangle,i|\eta\rangle)^T$に対応する。

3：$\psi\sigma_3=\eta+\phi\sigma_3$は$(|\eta\rangle,|\phi\rangle)^T=\hat\gamma_5\Psi$に対応する。
&&&

&&&thm 行列形式との同値性
$\psi\mapsto\Psi$のもとで、$i\hbar\sum_\mu\hat\gamma^\mu\partial_\mu\Psi=mc\Psi$と$\hbar D\psi\,\omega\sigma_3=mc\,\psi\gamma_0$は同値です。
&&&

&&&prf
[作用の翻訳](#prop-dirac-action)の1・2より、$i\hbar\hat\gamma^\mu\partial_\mu\Psi$には$\gamma^\mu(\hbar\,\partial_\mu\psi\,\omega\sigma_3)\gamma_0$が対応し、$\mu$について和を取ると$\hbar D\psi\,\omega\sigma_3\gamma_0$である。したがって行列形式の方程式は$\hbar D\psi\,\omega\sigma_3\gamma_0=mc\psi$と同値で、右から$\gamma_0$を掛ければ$\hbar D\psi\,\omega\sigma_3=mc\,\psi\gamma_0$となる。
&&&

行列形式の4本の複素数の方程式は、ヘステネス形式では偶部分代数の1本の式です。$\hat\gamma^\mu$を掛けることは$\gamma^\mu$で挟むことに、$i$を掛けることは右から$\omega\sigma_3$を掛けることになり、行列の側の複素数と$\hat\gamma$行列の区別は、代数の側では右から掛けるか左から掛けるかの区別になります。

## 観測量との対応

分解$\psi=\sqrt\rho\,e^{\omega\beta/2}R$と枠$\psi\gamma_\mu\tilde\psi=\rho e_\mu$は、行列形式でよく使われる2次式に対応します。$\bar\Psi=\Psi^\dagger\hat\gamma^0$と書くと

$$
\Psi^\dagger\Psi=J\cdot\gamma_0,\qquad
\bar\Psi\hat\gamma^\mu\Psi=J\cdot\gamma^\mu,\qquad
\bar\Psi\Psi=\rho\cos\beta,\qquad
\bar\Psi\,i\hat\gamma_5\Psi=-\rho\sin\beta
$$

です（$J=\psi\gamma_0\tilde\psi=\rho e_0$）。行列形式の確率密度$\Psi^\dagger\Psi$は、流れ$J$の$\gamma_0$成分です。$J$は未来向きの時間的ベクトルなので、確率密度はどの観測者から見ても正で、クライン＝ゴルドン方程式にあった困難はここで解消されています。行列形式で別々の2次式として現れる量が、ヘステネス形式では密度$\rho$、角$\beta$、枠$e_\mu$という幾何的な量にまとまります。[[7shi-dirac2]]

# 共変性

ディラック方程式がローレンツ変換で形を保つことを確かめます。電磁気学での議論と同じく、場を能動的に変換します。[[7shi-em5]]

&&&thm 共変性 [thm-cov]
一定の回転子$R$について、$\psi$がディラック方程式の解なら

$$
\psi'(x)=R\,\psi(\tilde RxR)
$$

も解です。
&&&

&&&prf
$x'=\tilde RxR$の成分は$x'^\nu=\gamma^\nu\cdot(\tilde RxR)=(R\gamma^\nu\tilde R)\cdot x$である。ベクトル$a$について$D(a\cdot x)=a$だから、連鎖律より、偶部分代数に値を取る関数$f$に対して

$$
D\bigl(f(x')\bigr)=\sum_\nu D(x'^\nu)\,(\partial_\nu f)(x')=\sum_\nu R\gamma^\nu\tilde R\,(\partial_\nu f)(x')
$$

となる。$f=R\psi$と置くと$\tilde RR=1$より$D\psi'(x)=R\,(D\psi)(x')$である。したがって

$$
\hbar D\psi'\,\omega\sigma_3-mc\,\psi'\gamma_0=R\bigl(\hbar D\psi\,\omega\sigma_3-mc\,\psi\gamma_0\bigr)(x')=0
$$

となる。
&&&

証明の要点は、ローレンツ変換が左から、方程式の中の$\omega\sigma_3$と$\gamma_0$が右から掛かり、互いに干渉しないことです。行列形式では、ローレンツ変換に対して$S^{-1}\hat\gamma^\mu S=\sum_\nu\Lambda^\mu{}_\nu\hat\gamma^\nu$を満たす行列$S$を構成して共変性を示します。ヘステネス形式では、その$S$が回転子$R$そのもので、$\hat\gamma^\mu$を変換する式は$R$で挟む式$R\gamma^\nu\tilde R$として証明の中に現れています。右から掛かる$\gamma_0$は観測者の時間軸ではなく、スピノルの成分を読む基準です。どの観測者の座標で書いても、同じ$\gamma_0$と$\omega\sigma_3$が右に現れます。[[7shi-dirac2]]

# 平面波解

## 正のエネルギー

4元運動量$p$（時空のベクトル）と一定のスピノル$\psi_0$により

$$
\psi=\psi_0\,e^{-\omega\sigma_3\,p\cdot x/\hbar}
$$

と置きます。行列形式の$e^{-ip\cdot x/\hbar}$にあたり、時間には$e^{-\omega\sigma_3Et/\hbar}$（$E=c\,p\cdot\gamma_0$）として依存します。$\partial_\mu(p\cdot x)=p\cdot\gamma_\mu$と$\sum_\mu\gamma^\mu(p\cdot\gamma_\mu)=p$、および$\omega\sigma_3$がその指数関数と可換であることから

$$
\hbar D\psi\,\omega\sigma_3=-p\,\psi\,(\omega\sigma_3)^2=p\,\psi
$$

です。したがってディラック方程式は、一定のスピノルについての代数方程式

$$
p\,\psi_0=mc\,\psi_0\gamma_0
$$

になります。$\psi_0$が可逆なら、左から$p$を掛けて$p^2\psi_0=mc\,p\psi_0\gamma_0=m^2c^2\psi_0\gamma_0^2$より$p^2=m^2c^2$で、エネルギーと運動量の関係が得られます。

&&&ex 静止した解
$p=mc\gamma_0$（$E=mc^2$、$\boldsymbol p=0$）なら、条件は$\gamma_0\psi_0=\psi_0\gamma_0$で、$\psi_0$は$\gamma_0$と可換、すなわちパウリスピノルです。解は

$$
\psi=\psi_0\,e^{-\omega\sigma_3mc^2t/\hbar}
$$

で、行列形式では下の2成分が$0$の$\Psi=(|\psi_0\rangle,0)^T e^{-imc^2t/\hbar}$です。$\psi_0$の4つの実数の自由度が、静止した粒子のスピンの状態（2つの複素数）にあたります。[[7shi-dirac1]]
&&&

運動する解は、静止した解をブーストして得られます。ブーストの回転子$L$について$p=L(mc\gamma_0)\tilde L$と置き、$\psi_0=L\phi_0$（$\phi_0$は$\gamma_0$と可換）とすると

$$
p\,\psi_0=L\,mc\gamma_0\,\tilde LL\,\phi_0=mc\,L\phi_0\gamma_0=mc\,\psi_0\gamma_0
$$

です。これは[共変性](#thm-cov)の平面波の場合です。このとき$\psi_0\tilde\psi_0$は$\phi_0\tilde\phi_0>0$なので$\beta=0$で、流れは$J=\psi_0\gamma_0\tilde\psi_0=\rho L\gamma_0\tilde L$です。流れの向き$e_0$は4元運動量$p=mc\,L\gamma_0\tilde L$の向きと一致し、粒子の運動の向きを表します。

## 負のエネルギー

$p$の向きを逆にした$-p$に対しても、代数方程式$(-p)\psi_0=mc\,\psi_0\gamma_0$は解を持ちます。$\psi_0$を右から$\sigma_1$倍すれば得られます。$\sigma_1$は$\gamma_0$とも$\omega\sigma_3$とも反可換なので、$p\psi_0=mc\psi_0\gamma_0$の解$\psi_0$について

$$
(-p)(\psi_0\sigma_1)=-mc\,\psi_0\gamma_0\sigma_1=mc\,(\psi_0\sigma_1)\gamma_0
$$

です。対応する平面波は

$$
\psi_0\sigma_1\,e^{+\omega\sigma_3\,p\cdot x/\hbar}=\psi_0\,e^{-\omega\sigma_3\,p\cdot x/\hbar}\,\sigma_1
$$

で、時間には$e^{+\omega\sigma_3Et/\hbar}$として依存します。行列形式の$e^{+iEt/\hbar}$、すなわちエネルギー$-E$の解です。

一般に、自由粒子の方程式の解$\psi$に右から$\sigma_1$を掛けたものはまた解になります。$D$は左から作用し、$\sigma_1$は右からの因子$\omega\sigma_3$と$\gamma_0$の両方と反可換で、両辺に同じ符号の変化を与えるからです。右からの作用で、正のエネルギーの解と負のエネルギーの解が移り合います。

各$p$（$p^2=m^2c^2$、未来向き）について、静止系で見れば、正のエネルギーの解の$\psi_0$は$\gamma_0$と可換な実4次元、負のエネルギーの解の$\psi_0$は反可換な実4次元で、合わせて偶部分代数の実8次元を尽くします。行列形式でいえば、スピンの2状態が正と負のエネルギーのそれぞれにあり、4成分になります。

負のエネルギーの解は$(\psi_0\sigma_1)\widetilde{(\psi_0\sigma_1)}=-\psi_0\tilde\psi_0$より$\beta=\pi$を持ちます。流れは$(\psi_0\sigma_1)\gamma_0\widetilde{(\psi_0\sigma_1)}=\psi_0\gamma_0\tilde\psi_0$で、正のエネルギーの解と同じ未来向きの$J$です。負のエネルギーの解の物理的な解釈（空孔理論や反粒子）は、本記事では扱いません。

# 質量0の場合とワイルスピノル

&&&rem 質量0のディラック方程式
$m=0$ではディラック方程式は$D\psi=0$となり、右からの因子は消えます。これは真空のマクスウェル方程式$DF=0$と同じ作用素の核で、違いは値の空間だけです（$F$は2ベクトル、$\psi$は偶部分代数全体）。さらに時間によらない$\psi$では、$\gamma_0D=\partial_0+D_3$（$D_3=\sum_k\sigma_k\partial_k$）より$D_3\psi=0$となり、偶部分代数を$\operatorname{Cl}_{3,0}(\mathbb R)$と見れば、3次元のモノジェニック関数（楕円型）です。クリフォード解析から電磁場を経てディラック方程式へ、同じ形の作用素の核を、係数の符号と値の空間を替えながら扱ってきたことになります。[[7shi-em4]][[7shi-cla5]]
&&&

&&&rem ワイルスピノル
$\hat\gamma_5$は右からの$\sigma_3$に対応するので、行列形式のカイラリティの射影$(1\pm\hat\gamma_5)/2$は右からの$P_\pm=(1\pm\sigma_3)/2$です。$D$は左から作用するので$D(\psi P_\pm)=(D\psi)P_\pm$で、$m=0$の方程式$D\psi=0$は$\psi P_+$と$\psi P_-$に分かれます。これらを**ワイルスピノル**と呼びます。質量の項$\psi\gamma_0$は$\gamma_0P_\pm=P_\mp\gamma_0$より2つを混ぜます。
&&&

&&&rem マヨラナ形式
本シリーズは符号数$(+,-,-,-)$の$\operatorname{Cl}_{1,3}(\mathbb R)\cong M_2(\mathbb H)$を使っています。逆の符号数の$\operatorname{Cl}_{3,1}(\mathbb R)$は実行列環$M_4(\mathbb R)$と同型で、ガンマ行列をすべて実行列に取り、ディラック方程式を実4成分のスピノルで書く流儀があります（マヨラナ表現）。ヘステネス形式も複素数を使わずに実数だけで書く定式化ですが、スピノルを列ベクトルでなく偶部分代数の元として扱う点で異なります。本シリーズではマヨラナ形式には立ち入りません。[[7shi-em4]][[7shi-clif1]]
&&&

# まとめ

ディラック方程式をヘステネス形式で書き、行列形式と対応させました。

- **方程式**：$\square$の平方根$D$を使い、$\hbar D\psi\,\omega\sigma_3=mc\,\psi\gamma_0$と書きます。右から掛かる$\omega\sigma_3$と$\gamma_0$は可換で、2乗の符号が$-1$と$+1$なので、2回使うとクライン＝ゴルドン方程式になります。
- **行列形式との対応**：$\psi=\phi+\eta\sigma_3$（$\phi,\eta$は$\gamma_0$と可換）を$\Psi=(|\phi\rangle,|\eta\rangle)^T$に移すと、$\hat\gamma_\mu\Psi$は$\gamma_\mu\psi\gamma_0$、$i\Psi$は$\psi\,\omega\sigma_3$、$\hat\gamma_5\Psi$は$\psi\sigma_3$に対応し、$i\hbar\hat\gamma^\mu\partial_\mu\Psi=mc\Psi$と同値です。確率密度$\Psi^\dagger\Psi$は未来向きの流れ$J=\psi\gamma_0\tilde\psi$の時間成分です。
- **共変性**：$\psi'(x)=R\psi(\tilde RxR)$はまた解です。ローレンツ変換は左から、方程式の中の因子は右から掛かるので干渉しません。
- **平面波**：$\psi=\psi_0e^{-\omega\sigma_3p\cdot x/\hbar}$は$p\psi_0=mc\psi_0\gamma_0$に帰着し、$p^2=m^2c^2$を与えます。静止した解はパウリスピノルで、運動する解はそのブーストです。右から$\sigma_1$を掛けると負のエネルギーの解が得られ、$\beta=\pi$を持ちます。
- **質量0**：$D\psi=0$は真空のマクスウェル方程式と同じ作用素の核で、右からの射影$(1\pm\sigma_3)/2$でワイルスピノルに分かれます。

&&& ディラック方程式
$$
\hbar D\psi\,\omega\sigma_3=mc\,\psi\gamma_0
$$
&&&

&&& 行列形式との対応
$\psi=\phi+\eta\sigma_3$を$\Psi=(|\phi\rangle,|\eta\rangle)^T$に移すと、次のように対応します。
$$
\hat\gamma_\mu\Psi\leftrightarrow\gamma_\mu\psi\gamma_0,\qquad i\Psi\leftrightarrow\psi\,\omega\sigma_3,\qquad\hat\gamma_5\Psi\leftrightarrow\psi\sigma_3
$$
&&&

&&& 平面波解
$\psi=\psi_0e^{-\omega\sigma_3p\cdot x/\hbar}$は次の代数方程式に帰着し、$\psi_0$が可逆なら$p^2=m^2c^2$が従います。
$$
p\psi_0=mc\,\psi_0\gamma_0
$$
&&&
