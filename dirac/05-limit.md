電磁場中のディラック方程式の非相対論的極限からパウリ方程式を導出し、スピンのg因子が2となる代数的機構を解明します。

# 概要

パウリ方程式のスピンと磁場の結合の係数$g$を、以前の記事では現象論的に置きました。[[7shi-dirac1]]

本記事では、電磁場の中のディラック方程式から、速さが光速より十分小さい極限でパウリ方程式を導き、$g=2$を得ます。[[7shi-dirac4]]

ディラックスピノルを時間軸$\gamma_0$と可換な部分と反可換な部分に分け、静止エネルギーの振動を分離すると、反可換な部分は$v/c$の程度に小さくなります。それを消去すると、空間のディラック作用素$D_3$とベクトルポテンシャルから作った1階の作用素$\Pi$の2乗が運動項として残ります。幾何積$D_3\boldsymbol A=D_3\cdot\boldsymbol A+D_3\wedge\boldsymbol A$の外積の部分が磁場の2ベクトル$\omega\boldsymbol B$で、これがスピンと磁場の結合を係数$g=2$で与えます。

# 時間軸による分解

## 可換な部分と反可換な部分

電磁場の中のディラック方程式

$$
\hbar\,D\psi\,\omega\sigma_3-\frac qcA\psi=mc\,\psi\gamma_0,\qquad A=\varphi\gamma_0+c\sum_kA_k\gamma_k
$$

から出発します。[[7shi-dirac4]]

観測者の時間軸$\gamma_0$を選び、ディラックスピノルを

$$
\psi=\psi_++\psi_-,\qquad\gamma_0\psi_\pm\gamma_0=\pm\psi_\pm
$$

と分けます。$\psi_+$は$\gamma_0$と可換、$\psi_-$は反可換な部分です。ディラックスピノルの対応$\psi=\phi+\eta\sigma_3\mapsto\Psi=(|\phi\rangle,|\eta\rangle)^T$では、$\psi_+=\phi$が$\Psi$の上の2成分、$\psi_-=\eta\sigma_3$が下の2成分にあたります。[[7shi-dirac3]]

右から掛かる$\omega\sigma_3$は$\gamma_0$と可換なので、分解を保ちます。左から掛かる相対ベクトル$\sigma_k$は$\gamma_0$と反可換なので、2つの部分を入れ替えます。

## 分割した方程式

作用素の分割$\gamma_0D=\partial_0+D_3$（$D_3=\sum_k\sigma_k\partial_k$）と、ポテンシャルの分割$\gamma_0A=\varphi-c\boldsymbol A$（$\boldsymbol A=\sum_kA_k\sigma_k$）を使うため、方程式に左から$\gamma_0$を掛けます。[[7shi-em4]][[7shi-em6]]

$$
\hbar(\partial_0+D_3)\psi\,\omega\sigma_3-\frac qc(\varphi-c\boldsymbol A)\psi=mc\,\gamma_0\psi\gamma_0
$$

右辺は$mc(\psi_+-\psi_-)$です。ここで、次の作用素を置きます。

&&&def 作用素$\Pi$
$$
\Pi X=-\hbar\,D_3X\,\omega\sigma_3-q\boldsymbol AX
$$
&&&

$\Pi$は、$\pi_k\psi=-\hbar\,\partial_k\psi\,\omega\sigma_3-qA_k\psi$に$\sigma_k$を掛けて足したもの$\sum_k\sigma_k\pi_k$で、行列形式の$\sum_k\hat\sigma_k\hat\pi_k$（$\hat\pi_k=-i\hbar\partial_k-qA_k$）にあたります。$\sigma_k$を含むので、$\Pi$は$\psi_+$と$\psi_-$を入れ替えます。[[7shi-dirac1]]

質量$m$の粒子の静止エネルギー$mc^2$に伴う速い振動を分離するため

$$
\psi=(\psi_++\psi_-)\,e^{-\omega\sigma_3mc^2t/\hbar}
$$

と置き直します（以下の$\psi_\pm$は、この因子を除いた部分を表します）。静止した解$\psi_0e^{-\omega\sigma_3mc^2t/\hbar}$と同じ因子で、$e^{-\omega\sigma_3mc^2t/\hbar}$は$\gamma_0$と可換なので、分解を保ちます。[[7shi-dirac3]]

&&&fml 分割した方程式 [fml-split-eq]
$$
\hbar\,\partial_t\psi_+\,\omega\sigma_3=q\varphi\,\psi_++c\,\Pi\psi_-
$$

$$
\hbar\,\partial_t\psi_-\,\omega\sigma_3=q\varphi\,\psi_-+c\,\Pi\psi_+-2mc^2\psi_-
$$
&&&

&&&prf
$x_0=ct$より$\hbar\,\partial_0=\frac\hbar c\partial_t$である。$\psi_\pm e^{-\omega\sigma_3mc^2t/\hbar}$の時間微分は$\bigl(\partial_t\psi_\pm-\psi_\pm\,\omega\sigma_3\frac{mc^2}\hbar\bigr)e^{-\omega\sigma_3mc^2t/\hbar}$であり、右から$\omega\sigma_3$を掛けると$(\omega\sigma_3)^2=-1$より

$$
\frac\hbar c\,\partial_t\bigl(\psi_\pm e^{\cdots}\bigr)\omega\sigma_3=\Bigl(\frac\hbar c\,\partial_t\psi_\pm\,\omega\sigma_3+mc\,\psi_\pm\Bigr)e^{\cdots}
$$

となる。ほかの項は右からの因子$e^{\cdots}$をそのまま保つので、方程式全体を$e^{\cdots}$で割ったものは

$$
\frac\hbar c\,\partial_t(\psi_++\psi_-)\omega\sigma_3+mc(\psi_++\psi_-)-\Pi(\psi_++\psi_-)-\frac qc\varphi(\psi_++\psi_-)=mc(\psi_+-\psi_-)
$$

である。$\gamma_0$と可換な部分と反可換な部分をそれぞれ比べ、$c$を掛ければよい。$\Pi$は2つの部分を入れ替えるので、$\Pi\psi_-$は可換な部分、$\Pi\psi_+$は反可換な部分に入る。
&&&

ここまでは近似を含まず、ディラック方程式と同値です。$\psi_+$の式から静止エネルギー$mc^2$が消え、$\psi_-$の式には$2mc^2$が残りました。

# 非相対論極限

## 小さい成分の消去

粒子の運動エネルギーとポテンシャルエネルギーが静止エネルギー$mc^2$より十分小さい場合を考えます。$\psi_\pm$の時間変化は運動エネルギー程度の振動数を持つので、$\psi_-$の式の左辺$\hbar\,\partial_t\psi_-\,\omega\sigma_3$と$q\varphi\psi_-$は、右辺の$2mc^2\psi_-$に比べて小さくなります。これらを無視すると

$$
\psi_-\approx\frac1{2mc}\Pi\psi_+
$$

です。$\Pi$は運動量程度の大きさ（$mv$）を持つので、$\psi_-$は$\psi_+$の$v/2c$倍程度です。$\psi_+$を**大きい成分**、$\psi_-$を**小さい成分**と呼びます。

&&&ex 自由粒子の平面波
正のエネルギーの平面波$\psi_0=L\phi_0$（$L=e^{\sigma_1\eta/2}$、$\phi_0$は$\gamma_0$と可換）では、$L=\cosh\frac\eta2+\sigma_1\sinh\frac\eta2$の$\cosh\frac\eta2$の項が大きい成分、$\sigma_1\sinh\frac\eta2$の項が小さい成分を与えます。大きさの比は$\tanh\frac\eta2$で、運動量$p=mc\sinh\eta$、エネルギー$E=mc^2\cosh\eta$で書けば

$$
\tanh\frac\eta2=\frac{\sinh\eta}{1+\cosh\eta}=\frac p{mc+E/c}
$$

です。$p\ll mc$では$E\approx mc^2$で、比は$p/2mc\approx v/2c$になります。[[7shi-dirac3]]
&&&

## パウリ方程式の回収

$\psi_-\approx\Pi\psi_+/2mc$を$\psi_+$の式に代入すると

$$
\hbar\,\partial_t\psi_+\,\omega\sigma_3=\frac1{2m}\Pi^2\psi_++q\varphi\,\psi_+
$$

を得ます。$\psi_+$は$\gamma_0$と可換で、パウリスピノルの空間の元です。残るのは$\Pi^2$の計算です。[[7shi-dirac1]]

&&&prop $\Pi$の2乗 [prop-pi2]
$\gamma_0$と可換な$X$について

$$
\Pi^2X=\sum_k\pi_k^2X-q\hbar\,\boldsymbol BX\sigma_3,\qquad\boldsymbol B=\nabla\times\boldsymbol A
$$

が成り立ちます。ここで$\pi_kX=-\hbar\,\partial_kX\,\omega\sigma_3-qA_kX$です。
&&&

&&&prf
$\Pi^2X=-\hbar D_3(\Pi X)\omega\sigma_3-q\boldsymbol A\,\Pi X$に$\Pi X=-\hbar D_3X\omega\sigma_3-q\boldsymbol AX$を代入すると、$(\omega\sigma_3)^2=-1$と$D_3^2=\Delta$より

$$
\Pi^2X=-\hbar^2\Delta X+q\hbar\bigl(D_3(\boldsymbol AX)+\boldsymbol AD_3X\bigr)\omega\sigma_3+q^2|\boldsymbol A|^2X
$$

である。積の微分により$D_3(\boldsymbol AX)=\sum_{k,l}\sigma_k\sigma_l(\partial_kA_l)X+\sum_{k,l}\sigma_k\sigma_lA_l\partial_kX$であり、第1項は$\boldsymbol A$だけを微分した$(D_3\boldsymbol A)X$である。第2項と$\boldsymbol AD_3X=\sum_{k,l}\sigma_l\sigma_kA_l\partial_kX$を合わせると、$\sigma_k\sigma_l+\sigma_l\sigma_k=2\delta_{kl}$より$2\sum_kA_k\partial_kX$となる。幾何積を内積と外積に分けると$D_3\boldsymbol A=\nabla\cdot\boldsymbol A+D_3\wedge\boldsymbol A=\nabla\cdot\boldsymbol A+\omega\boldsymbol B$だから

$$
\Pi^2X=-\hbar^2\Delta X+q\hbar\Bigl((\nabla\cdot\boldsymbol A)X+2\sum_kA_k\partial_kX\Bigr)\omega\sigma_3+q^2|\boldsymbol A|^2X+q\hbar\,\omega\boldsymbol BX\,\omega\sigma_3
$$

となる。最後の項は、$\omega$が偶部分代数の元と可換で$\omega^2=-1$だから$-q\hbar\boldsymbol BX\sigma_3$である。残りの項は、$\pi_k^2X=-\hbar^2\partial_k^2X+q\hbar\bigl((\partial_kA_k)X+2A_k\partial_kX\bigr)\omega\sigma_3+q^2A_k^2X$の$k$についての和に等しい。
&&&

&&&thm 非相対論極限
小さい成分を消去すると、大きい成分$\psi_+$は

$$
\hbar\,\partial_t\psi_+\,\omega\sigma_3=\frac1{2m}\sum_k\pi_k^2\psi_++q\varphi\,\psi_+-\frac{q\hbar}{2m}\boldsymbol B\psi_+\sigma_3
$$

に従います。これはパウリ方程式で$g=2$と置いたものです。[[7shi-dirac1]]
&&&

パウリ方程式の磁場の項は$-\frac{gq\hbar}{4m}\boldsymbol B\psi\sigma_3$でした。$g=2$でこれが$-\frac{q\hbar}{2m}\boldsymbol B\psi\sigma_3$になります。現象論的に置いた係数が、ディラック方程式から決まりました。[[7shi-dirac1]]

## g=2の出所

磁場の項は、[$\Pi$の2乗](#prop-pi2)の計算のうち、幾何積$D_3\boldsymbol A$の外積の部分$D_3\wedge\boldsymbol A=\omega\boldsymbol B$から出ています。時空のポテンシャル$A$から$F=D\wedge A$として電磁場を得たように、その空間部分が、ここで$D_3\wedge\boldsymbol A=\omega\nabla\times\boldsymbol A$として現れた磁場の2ベクトルです。[[7shi-em6]]

スカラーの運動項$\sum_k\pi_k^2$だけから出発したのでは、磁場とスピンの結合は現れません。ディラック方程式は1階の作用素から出発するので、その2乗を計算すると、ベクトル$\sigma_k$どうしの積が内積（$\sum_k\pi_k^2$）と外積に分かれ、外積の部分が磁場の2ベクトルを与えます。ベクトル解析では、パウリ行列の積の公式$\hat\sigma_k\hat\sigma_l=\delta_{kl}I+i\sum_m\varepsilon_{klm}\hat\sigma_m$を使って$(\hat{\boldsymbol\sigma}\cdot\hat{\boldsymbol\pi})^2=\hat{\boldsymbol\pi}^2-q\hbar\hat{\boldsymbol\sigma}\cdot\boldsymbol B$と計算する部分です。幾何代数では、これは幾何積を内積と外積に分けることそのものです。

係数も同じ見方で決まります。$\Pi^2$の中で外積の部分は、$D_3(\boldsymbol AX)$から1回だけ現れ、その係数$q\hbar$がそのまま$\frac1{2m}$倍されて磁場の項になります。パウリ方程式の記号では$\frac{gq\hbar}{4m}=\frac{q\hbar}{2m}$、すなわち$g=2$です。[[7shi-dirac1]]

以前に見たとおり、$g=2$のとき、一様な磁場の中でスピンの向きを回す方程式$\partial_t\psi=\frac q{2m}\omega\boldsymbol B\psi$は、固有速度の向きを回す方程式と同じ形になります。ディラック方程式に従う粒子は、速さが光速より十分小さい範囲では、磁場の中でスピンの向きと運動の向きが同じ角速度で回ります。[[7shi-dirac1]][[7shi-em5]]

&&&rem 近似の範囲と実際の値
本記事の近似は、$\psi_-$の式で$\hbar\,\partial_t\psi_-$と$q\varphi\psi_-$を落としたもので、$(v/c)^2$の程度の補正（運動エネルギーの相対論的な補正、スピン軌道相互作用、ダーウィン項と呼ばれる項）を含みません。これらは同じ分解を次の次数まで進めると得られますが、本記事では扱いません。

また、電子の$g$の実測値は$2.0023\ldots$で、$2$からわずかにずれています。このずれは電磁場そのものを量子化した理論（量子電磁力学）で説明されるもので、本シリーズでは扱いません。
&&&

# 行列形式との照合

ディラックスピノルの対応では、$\psi_+=\phi$が上の成分$|\phi\rangle$、$\psi_-=\eta\sigma_3$が下の成分$|\eta\rangle=|\psi_-\sigma_3\rangle$を与えます。作用の翻訳より、$\gamma_0$と可換な$X$について$(\Pi X)\sigma_3=\sum_k\sigma_k(\pi_kX)\sigma_3$は$\sum_k\hat\sigma_k\hat\pi_k|X\rangle$に対応します。[[7shi-dirac3]][[7shi-dirac1]]

したがって

$$
\psi_-\approx\frac1{2mc}\Pi\psi_+\quad\longleftrightarrow\quad|\eta\rangle\approx\frac1{2mc}\bigl(\hat{\boldsymbol\sigma}\cdot\hat{\boldsymbol\pi}\bigr)|\phi\rangle
$$

であり、行列形式でよく知られた小さい成分の近似式です。$\Pi^2$は$(\hat{\boldsymbol\sigma}\cdot\hat{\boldsymbol\pi})^2=\hat{\boldsymbol\pi}^2-q\hbar\,\hat{\boldsymbol\sigma}\cdot\boldsymbol B$に対応します。行列形式では上下の2成分への分け方として現れる分解が、ヘステネス形式では時間軸$\gamma_0$との可換・反可換として現れています。この分け方は観測者の時間軸の選択に依存し、電場と磁場を$\frac12(F\mp\gamma_0F\gamma_0)$に分けたのと同じ形です。非相対論極限は、特定の観測者に対して粒子がゆっくり動く場合の近似です。[[7shi-em4]]

# まとめ

電磁場の中のディラック方程式から、非相対論極限でパウリ方程式を導きました。

- **分解**：$\psi=\psi_++\psi_-$（$\gamma_0\psi_\pm\gamma_0=\pm\psi_\pm$）と分け、静止エネルギーの因子$e^{-\omega\sigma_3mc^2t/\hbar}$を分離すると、方程式は[分割した方程式](#fml-split-eq)の2本と同値になります。$\Pi X=-\hbar D_3X\omega\sigma_3-q\boldsymbol AX$は2つの部分を入れ替えます。
- **小さい成分**：非相対論極限では$\psi_-\approx\Pi\psi_+/2mc$で、$\psi_-$は$\psi_+$の$v/2c$倍程度です。
- **パウリ方程式**：消去すると$\hbar\,\partial_t\psi_+\omega\sigma_3=\frac1{2m}\Pi^2\psi_++q\varphi\psi_+$で、$\Pi^2=\sum_k\pi_k^2-q\hbar\boldsymbol B(\cdot)\sigma_3$より、パウリ方程式が$g=2$で得られます。
- **$g=2$の出所**：幾何積$D_3\boldsymbol A$の外積の部分$D_3\wedge\boldsymbol A=\omega\boldsymbol B$が、スピンと磁場の結合を与えます。1階の作用素を2乗すると、内積の部分が運動項、外積の部分が磁場の項になります。
- **歳差**：$g=2$のとき、一様な磁場の中でスピンの向きと運動の向きは同じ角速度で回ります。

&&& 時間軸による分解
$$
\psi=\psi_++\psi_-,\qquad\gamma_0\psi_\pm\gamma_0=\pm\psi_\pm
$$
&&&

&&& 小さい成分
$\Pi X=-\hbar\,D_3X\,\omega\sigma_3-q\boldsymbol AX$とすると、非相対論極限では次のようになります。
$$
\psi_-\approx\frac1{2mc}\Pi\psi_+
$$
&&&

&&& パウリ方程式（$g=2$）
小さい成分を消去すると、$\psi_+$について次の式が得られます。
$$
\hbar\,\partial_t\psi_+\,\omega\sigma_3=\frac1{2m}\sum_k\pi_k^2\psi_++q\varphi\,\psi_+-\frac{q\hbar}{2m}\boldsymbol B\psi_+\sigma_3
$$
&&&
