パウリスピノルを実クリフォード代数$\operatorname{Cl}_{3,0}(\mathbb R)$の偶部分代数の元（四元数）として再定式化し、パウリ方程式とスピンの歳差運動を幾何学的に解説します。

# 概要

量子力学の標準的な定式化では、電子のスピンは複素数を2つ並べた列ベクトル$\Psi\in\mathbb C^2$で表され、パウリ行列がそれに作用します。本記事では、この列ベクトルを$\operatorname{Cl}_{3,0}(\mathbb R)$の偶部分代数の元$\psi$に置き換えます。偶部分代数は四元数と同型なので、パウリスピノルは四元数そのものです。対応の要点は虚数単位$i$の読み替えにあり、$i$倍は$\psi$に右から2ベクトル$\omega\sigma_3=\sigma_1\sigma_2$を掛ける操作になります。この表示では、確率密度が$\psi\tilde\psi$、スピンの向きが$\psi\sigma_3\tilde\psi$として読み取れ、大域位相はスピン軸まわりの回転として見えます。最後に、パウリ方程式をこの表示に翻訳し、一様な磁場の中のスピンの歳差を回転子の時間発展として解きます。

このような表示は、スピノルを列ベクトルでなく代数の元として扱う定式化で、ヘステネス（D. Hestenes）によって整備されました。本シリーズでは、これを**ヘステネス形式**と呼びます。

# 準備

## 記号

電磁気学のこれまでの記事と同じく3次元の実クリフォード代数$\operatorname{Cl}_{3,0}(\mathbb R)$を使いますが、本シリーズでは生成元を$\sigma_1,\sigma_2,\sigma_3$と書きます（$\sigma_k^2=1$、$k\ne l$なら$\sigma_k\sigma_l=-\sigma_l\sigma_k$）。擬スカラーは

$$
\omega=\sigma_1\sigma_2\sigma_3
$$

と書き、$\omega^2=-1$で、$\omega$はすべての元と可換です。emシリーズでは擬スカラーを$i$、複素数の虚数単位を$j$と書きましたが、本シリーズでは量子力学の慣用に合わせて、複素数の虚数単位を$i$、擬スカラーを$\omega$と書きます。emシリーズの式を引くときは$F=\boldsymbol E+\omega c\boldsymbol B$のように書き換えます。[[7shi-em1]][[7shi-em3]]

ベクトル$\boldsymbol a=\sum_ka_k\sigma_k$と2ベクトルの対応は、以前の記事と同じく$\omega\boldsymbol a$です。[[7shi-em1]]

回転子の作用は、emシリーズや[[7shi-lie3]]と同じく$x\mapsto Rx\tilde R$とします。

たとえば

$$
\omega\sigma_1=\sigma_2\sigma_3,\qquad\omega\sigma_2=\sigma_3\sigma_1,\qquad\omega\sigma_3=\sigma_1\sigma_2
$$

で、$\omega\sigma_3$は$\sigma_3$に垂直な$x_1x_2$平面を表します。反転$\tilde X$は基底の積の順序を逆にする操作で、ベクトルを変えず、2ベクトルと擬スカラーの符号を変えます。空間のディラック作用素を$D_3=\sum_k\sigma_k\partial_k$と書きます。

代数の元と、それを表す行列を区別するため、行列表示と、列スピノルに作用する演算子にはハットを付けます。パウリ行列は

$$
\hat\sigma_1=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
\hat\sigma_2=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad
\hat\sigma_3=\begin{pmatrix}1&0\\0&-1\end{pmatrix}
$$

です。$\sigma_k\mapsto\hat\sigma_k$は$\operatorname{Cl}_{3,0}(\mathbb R)\cong M_2(\mathbb C)$を与え、擬スカラーは$\hat\omega=\hat\sigma_1\hat\sigma_2\hat\sigma_3=iI$に移ります。反転はエルミート共役に移ります（$\hat\sigma_k$はエルミートで、積の順序が逆になるため）。[[7shi-bq]][[7shi-lie3]]

## スピン

電子は、空間の中の位置のほかに、向きを持つ内部の自由度を持っています。銀原子のビームを不均一な磁場に通すと、ビームが2本に分かれます（シュテルン＝ゲルラッハの実験）。磁場から受ける力が原子の持つ小さな磁石の向きで決まり、その成分が測定すると2通りの値しか取らないためです。基底状態の銀原子では、主に不対電子のスピンに由来する磁気モーメントによって、この2本への分裂が生じます。この内部の自由度を**スピン**と呼びます。選んだ測定軸に沿うスピンの成分は、測定すると$+\hbar/2$または$-\hbar/2$のどちらかになります。これは、状態の向きが2通りしかないという意味ではありません。本記事で扱う電子のようなスピン1/2粒子の状態は、各点で2つの複素数の組$\Psi=(\Psi_1,\Psi_2)^T$で表され、$\Psi$を**パウリスピノル**と呼びます。

本シリーズでは、シュレーディンガー方程式$i\hbar\,\partial_t\Psi=\hat H\Psi$と、$\Psi^\dagger\Psi$を確率密度とする解釈を前提とします。[[7shi-clif5]][[7shi-born]]

スピンについて使うのは次の事実だけです。

- $\hat\sigma_k$は、$x_k$軸方向のスピンを測る観測量（の$2/\hbar$倍）です。
- 空間の回転は、$\operatorname{SU}(2)$の行列を左から掛けることで$\Psi$に作用します。[[7shi-lie3]]

# 偶部分代数による表示

## 対応の定義

$\operatorname{Cl}_{3,0}(\mathbb R)$の偶部分代数$\operatorname{Cl}_{3,0}^0(\mathbb R)$は$1,\omega\sigma_1,\omega\sigma_2,\omega\sigma_3$が張る4次元の空間で、四元数と同型です。双四元数の記事における対応$\mathbf i\cong-i\hat\sigma_1$、$\mathbf j\cong-i\hat\sigma_2$、$\mathbf k\cong-i\hat\sigma_3$は、本シリーズの記号では次のように表されます。[[7shi-bq]][[7shi-lie3]]の基底の対応とは異なりますが、同じ四元数代数の別の表現です。

$$
\mathbf i=-\omega\sigma_1,\qquad\mathbf j=-\omega\sigma_2,\qquad\mathbf k=-\omega\sigma_3
$$

実際$(\omega\sigma_k)^2=\omega^2\sigma_k^2=-1$、$(-\omega\sigma_1)(-\omega\sigma_2)=\omega^2\sigma_1\sigma_2=-\omega\sigma_3$で、$\mathbf i\mathbf j=\mathbf k$が成り立ちます。

&&&def パウリスピノルの対応 [def-map]
偶部分代数の元

$$
\psi=a_0+a_1\,\omega\sigma_1+a_2\,\omega\sigma_2+a_3\,\omega\sigma_3\qquad(a_0,\dots,a_3\in\mathbb R)
$$

に、列ベクトル

$$
\Psi=\begin{pmatrix}a_0+ia_3\\-a_2+ia_1\end{pmatrix}\in\mathbb C^2
$$

を対応させます。
&&&

この対応は実線形で、実4次元の空間$\operatorname{Cl}_{3,0}^0(\mathbb R)$と$\mathbb C^2$の間の全単射です。成分の並べ方は恣意的に見えますが、次の命題が示すとおり、$\Psi$は$\psi$を表す行列の第1列です。

&&&prop 作用の翻訳 [prop-action]
[対応](#def-map)$\psi\mapsto\Psi$について、次が成り立ちます。

1. $\Psi$は行列$\hat\psi$の第1列です。
2. $\hat\sigma_k\Psi$には$\sigma_k\psi\sigma_3$が対応します（$k=1,2,3$）。
3. $i\Psi$には$\psi\,\omega\sigma_3$が対応します。
&&&

&&&prf
$\hat\omega=iI$より

$$
\hat\psi=a_0I+i\sum_ka_k\hat\sigma_k=\begin{pmatrix}a_0+ia_3&a_2+ia_1\\-a_2+ia_1&a_0-ia_3\end{pmatrix}
$$

であり、第1列は$\Psi$である。第2列は第1列の成分から$(-\Psi_2^*,\Psi_1^*)^T$として決まるので、偶部分代数の元は第1列で決まる。$\boldsymbol e=(1,0)^T$と置くと$\Psi=\hat\psi\boldsymbol e$であり、$\hat\sigma_3\boldsymbol e=\boldsymbol e$、$\widehat{\omega\sigma_3}\boldsymbol e=i\hat\sigma_3\boldsymbol e=i\boldsymbol e$である。$\sigma_k\psi\sigma_3$と$\psi\,\omega\sigma_3$はともに偶部分代数の元で、その第1列は

$$
\hat\sigma_k\hat\psi\hat\sigma_3\boldsymbol e=\hat\sigma_k\Psi,\qquad
\hat\psi\,\widehat{\omega\sigma_3}\,\boldsymbol e=i\Psi
$$

となる。
&&&

パウリ行列$\hat\sigma_k$を掛けることは左から$\sigma_k$、右から$\sigma_3$を掛けることに、虚数単位$i$を掛けることは右から$\omega\sigma_3$を掛けることになります。$\sigma_k\psi$は奇部分の元になるので、右から$\sigma_3$を掛けて偶部分代数に戻しています。四元数の記号では$\psi=a_0-a_1\mathbf i-a_2\mathbf j-a_3\mathbf k$であり、パウリスピノルは四元数そのものです。ここでの同一視は実ベクトル空間としてのものです。

## 左イデアルとの関係

以前の記事では、射影$P=(1+\sigma_3)/2$を右から掛けて、行列の第1列だけを残した左イデアル$\operatorname{Cl}_{3,0}(\mathbb R)P$の元としてスピノルを取りました。偶部分代数による表示は、この表示の言い直しです。[[7shi-lie3]][[7shi-ideal]]

&&&prop 偶部分と左イデアル
$\psi\mapsto\psi P$は$\operatorname{Cl}_{3,0}^0(\mathbb R)$から$\operatorname{Cl}_{3,0}(\mathbb R)P$への実線形な全単射で、$\widehat{\psi P}=\begin{pmatrix}\Psi&\boldsymbol 0\end{pmatrix}$です。さらに

$$
(\psi\,\omega\sigma_3)P=\omega\,(\psi P)
$$

が成り立ちます。
&&&

&&&prf
$\hat P=\operatorname{diag}(1,0)$なので、$\hat\psi\hat P$は$\hat\psi$の第1列$\Psi$だけを残した行列である。偶部分代数の元は第1列で決まるから単射であり、左イデアルの元は第1列$\mathbb C^2$で決まる実4次元の空間なので全単射である。$\sigma_3P=P$と、$\omega$が中心の元であることから、$\psi\,\omega\sigma_3P=\psi\,\omega P=\omega\,\psi P$となる。
&&&

左イデアルの上では、右からの$\omega\sigma_3$が中心の元$\omega$の作用になり、$\hat\omega=iI$として行列の虚数単位に一致します。以前の記事の命題「偶部分と極小左イデアルの対応」は、$\operatorname{Cl}_{2,0}(\mathbb R)$で同じことを$P=(1+e_1)/2$について示したものです。ただし、そちらでは偶部分の元を$F$、左イデアルの元を$\psi=FP$と書いており、本シリーズとは文字の役割が逆です。本シリーズでは、偶部分代数の元を$\psi$、列ベクトルを$\Psi$と書きます（[[7shi-lie3]]では列スピノルを$\omega$と書いています）。[[7shi-cla1]]

## 右から掛かる虚数単位

$i$倍が右からの$\omega\sigma_3$になることは、本記事の中心となる論点です。

$(\omega\sigma_3)^2=-1$なので、$\omega\sigma_3$は複素数の虚数単位と同じ代数的な役割を果たし、位相因子は

$$
e^{i\alpha}\Psi\quad\longleftrightarrow\quad\psi\,e^{\omega\sigma_3\alpha}=\psi\,(\cos\alpha+\omega\sigma_3\sin\alpha)
$$

と対応します。ただし$\omega\sigma_3$は$\sigma_1$や$\sigma_2$とは反可換で、偶部分代数の中心の元ではありません。それでも複素数のスカラー倍として振る舞えるのは、右から掛かるからです。回転による左乗算は、右からの$\omega\sigma_3$の乗算と可換です。パウリ行列の作用$\psi\mapsto\sigma_k\psi\sigma_3$も、$\sigma_3$と$\omega\sigma_3$が可換なので、右からの$\omega\sigma_3$の乗算と可換です。この意味で、右乗算$\omega\sigma_3$は複素構造を与え、複素線形な作用素は、この右乗算と可換な実線形作用素として表されます。

電磁気学の記事では、複素指数関数の虚数単位の役割を擬スカラー$\omega$が担いました。$\omega$はすべての元と可換なので、どちらから掛けても同じです。ここで$i$の役割を担う$\omega\sigma_3$は、特定の平面$\sigma_1\sigma_2$を表す2ベクトルであり、掛ける側が意味を持ちます。量子力学の虚数単位は、スピノルの側では、1つの平面を選んだ2ベクトルとして現れます。どの平面を選ぶかは$\hat\sigma_3$を対角にする表現の選択に対応していて、[対応](#def-map)の成分の並べ方もこの選択に由来します。[[7shi-em3]]

# 観測量

## 確率密度と回転子

$\psi$の反転$\tilde\psi$は、行列表示ではエルミート共役$\hat\psi^\dagger$に対応します。$\psi\tilde\psi$はスカラーで、次のように列スピノルのノルムの2乗に一致します。

&&&fml 確率密度
$$
\rho=\psi\tilde\psi=a_0^2+a_1^2+a_2^2+a_3^2=\Psi^\dagger\Psi
$$
&&&

四元数のノルムの2乗が、$\Psi$の確率密度です。$\rho\ne0$なら$R=\psi/\sqrt\rho$は$R\tilde R=1$を満たす偶部分代数の元、すなわち回転子です。[[7shi-em5]]

$$
\psi=\sqrt\rho\,R
$$

パウリスピノルは、密度の平方根と回転子の積に分解されます。$\psi$の4つの実数の自由度のうち、1つが密度、残りの3つが回転子です。

## スピンの向き

回転子$R$は$\sigma_3$を回して新しい向き$R\sigma_3\tilde R$を作ります。これがスピンの向きです。

&&&fml スピンの向き
$$
\psi\sigma_3\tilde\psi=\rho\,\boldsymbol s,\qquad\boldsymbol s=R\sigma_3\tilde R=\sum_ks_k\sigma_k,\qquad
\rho\,s_k=\Psi^\dagger\hat\sigma_k\Psi
$$
&&&

&&&prf
$\psi\sigma_3\tilde\psi$は反転で変わらない奇数グレードの元なので、ベクトルである。成分（$\langle X\rangle_0$はスカラー部分）$\langle\psi\sigma_3\tilde\psi\,\sigma_k\rangle_0$を行列で計算すると$\frac12\operatorname{tr}(\hat\psi\hat\sigma_3\hat\psi^\dagger\hat\sigma_k)$である。$\hat\sigma_3=2\hat P-I$と書くと、$\hat\psi\hat P\hat\psi^\dagger=\Psi\Psi^\dagger$より

$$
\tfrac12\operatorname{tr}\bigl((2\Psi\Psi^\dagger-\hat\psi\hat\psi^\dagger)\hat\sigma_k\bigr)=\Psi^\dagger\hat\sigma_k\Psi-\tfrac12\rho\operatorname{tr}\hat\sigma_k=\Psi^\dagger\hat\sigma_k\Psi
$$

となる（$\hat\psi\hat\psi^\dagger=\rho I$）。$\rho>0$の点では$\boldsymbol s=R\sigma_3\tilde R$が定義され、$\boldsymbol s^2=R\sigma_3\tilde RR\sigma_3\tilde R=1$で単位ベクトルである。
&&&

$\rho s_k=\Psi^\dagger\hat\sigma_k\Psi$は、$\hat\sigma_k$の期待値密度です。物理的なスピン成分の密度は$(\hbar/2)\rho s_k$です。単位ベクトル$\boldsymbol s$はブロッホ球の記事のブロッホベクトルにあたります。行列形式では期待値を$\Psi^\dagger\hat\sigma_k\Psi$として計算しますが、ヘステネス形式では、スピノルが基準の向き$\sigma_3$を回した結果$\boldsymbol s$そのものが観測量です。[[7shi-bloch]]

ホップファイブレーションの記事では、単位四元数$q$から$q\mathbf kq^*$を作る写像としてホップ写像を導入しました。[[7shi-h]]

本記事の記号では$\sigma_3=\omega\mathbf k$で、四元数の共役は反転にあたるので

$$
\psi\sigma_3\tilde\psi=\omega\,(\psi\mathbf k\tilde\psi)
$$

です。ホップ写像に対応するのは単位回転子$R$についての$\boldsymbol s=\omega(R\mathbf k\tilde R)$で、一般の$\psi$の式はそれに$\rho$を掛けたものです。スピンの向きはホップ写像と同じ形の式で与えられ、両者は$\omega$を掛ける分だけ違います。$R\mathbf k\tilde R$は$\operatorname{Cl}_{3,0}(\mathbb R)$の2ベクトルで、スピンの向き$\boldsymbol s$はそれに垂直なベクトルです。通常の四元数表示では純虚四元数を3次元ベクトルと同一視しますが、ここではそれを2ベクトルとして埋め込み、対応するベクトルを$\omega$倍で区別しています。回転の面（2ベクトル）と回転の軸（ベクトル）は、異なるグレードの量として区別されます。

## 大域位相とスピン軸まわりの回転

$\Psi$に位相$e^{i\alpha}$を掛けても、確率密度と期待値は変わりません。ヘステネス形式では、これは$\psi\mapsto\psi e^{\omega\sigma_3\alpha}$であり

$$
\psi e^{\omega\sigma_3\alpha}\,\sigma_3\,e^{-\omega\sigma_3\alpha}\tilde\psi=\psi\sigma_3\tilde\psi
$$

です。$e^{\omega\sigma_3\alpha}$は$\sigma_3$を軸とする回転子で、$\sigma_3$を動かさないからです。回転子$R$が$\sigma_3$を$\boldsymbol s$に移すとき、$Re^{\omega\sigma_3\alpha}$も$\sigma_3$を$\boldsymbol s$に移します。ホップファイバーは、同じ$\boldsymbol s$を与える回転子の円周です。[[7shi-h]]

位相が変えるものもあります。回転子は$\sigma_3$だけでなく$\sigma_1,\sigma_2$も回し、正規直交基底$e_k=R\sigma_k\tilde R$を作ります。$e_3=\boldsymbol s$で、$e_1,e_2$は$\boldsymbol s$に垂直な平面の中の2本です。位相を掛けると

$$
Re^{\omega\sigma_3\alpha}\,\sigma_1\,e^{-\omega\sigma_3\alpha}\tilde R=e_1\cos2\alpha-e_2\sin2\alpha
$$

となり、$e_1,e_2$はスピン軸$\boldsymbol s$のまわりに角$-2\alpha$だけ回ります（$e^{\omega\sigma_3\alpha}=e^{-\omega\sigma_3(-2\alpha)/2}$で、回転の角は指数の係数の2倍です）。確率密度とスピンの向きはこの回転で変わりません。$\alpha$が位置・時間に依存しない大域位相なら、通常の量子力学と同じくすべての期待値が不変なので、大域位相は観測にかかりません。大域位相は、スピン軸まわりの枠の回転として見えています。$e_1,e_2$は独立に観測される物理的な軸ではなく、位相を枠の回転として見ているだけです。

## 回転と半角

空間の回転は、スピノルに左から回転子$U$を掛ける操作$\psi\mapsto U\psi$です。行列では$\operatorname{SU}(2)$の行列を掛けることにあたります。このとき$\rho$は変わらず、スピンの向きは

$$
U\psi\,\sigma_3\,\widetilde{U\psi}=U(\psi\sigma_3\tilde\psi)\tilde U
$$

より、ベクトルとしての回転$\boldsymbol s\mapsto U\boldsymbol s\tilde U$を受けます。スピノルは片側から、ベクトルは両側から回転子を受けます。角$\theta$の回転の回転子$U=e^{-\omega\boldsymbol n\theta/2}$は半角で書かれ、$\theta=2\pi$では$U=-1$です。$2\pi$回転でスピノルは$-\psi$になり、スピンの向きは元に戻ります。[[7shi-lie3]][[7shi-cover]]

左からの作用（回転）と右からの作用（位相）は、どちらも回転子を掛ける操作ですが、意味が違います。左からの作用は空間の中の回転で、観測量$\boldsymbol s$を動かします。右からの作用は、基準として選んだ$\sigma_3$のまわりの回転で、観測量を動かしません。

# パウリ方程式

## 行列形式

電荷$q$、質量$m$のスピン1/2粒子が、スカラーポテンシャル$\varphi$とベクトルポテンシャル$\boldsymbol A$の電磁場の中にあるとします。スピンを持つ粒子は小さな磁石としても振る舞い、磁場$\boldsymbol B$の中でスピンの向きに応じたエネルギーを持ちます。この項を含めたシュレーディンガー方程式を**パウリ方程式**と呼びます。

&&&def パウリ方程式（行列形式）
$$
i\hbar\,\partial_t\Psi=\frac1{2m}\sum_k\hat\pi_k^2\Psi+q\varphi\Psi-\frac{gq\hbar}{4m}\sum_kB_k\hat\sigma_k\Psi,\qquad
\hat\pi_k=-i\hbar\,\partial_k-qA_k
$$
&&&

$\hat\pi_k$は運動量$-i\hbar\,\partial_k$からベクトルポテンシャルの寄与を引いたもので、自由粒子の運動項に電磁場を加えた形です。最後の項がスピンと磁場の結合で、$\frac{\hbar}2\hat\sigma_k$をスピン、$\frac{gq}{2m}\cdot\frac\hbar2\hat\sigma_k$を磁気モーメントとして、磁気モーメントと磁場の内積の符号を変えたものがエネルギーになります。係数$g$は、スピンの大きさに対する磁気モーメントの比を表す無次元の数で、電子では実験的に$g\approx2$です。本記事では$g$を現象論的に置いた定数として扱います。[[7shi-clif5]]

## 翻訳

[作用の翻訳](#prop-action)を各項に使います。$i$倍は右からの$\omega\sigma_3$、$\hat\sigma_k$は左から$\sigma_k$、右から$\sigma_3$を掛ける操作なので

$$
\hat\pi_k\Psi\ \longleftrightarrow\ \pi_k\psi=-\hbar\,\partial_k\psi\,\omega\sigma_3-qA_k\psi,\qquad
\sum_kB_k\hat\sigma_k\Psi\ \longleftrightarrow\ \boldsymbol B\psi\sigma_3
$$

です（$\boldsymbol B=\sum_kB_k\sigma_k$）。

&&&fml パウリ方程式（ヘステネス形式） [fml-pauli]
$$
\hbar\,\partial_t\psi\,\omega\sigma_3=\frac1{2m}\sum_k\pi_k^2\psi+q\varphi\psi-\frac{gq\hbar}{4m}\boldsymbol B\psi\sigma_3
$$
&&&

実数の係数と、左右から掛かる代数の元だけで書かれています。右から掛かるのは$\omega\sigma_3$と$\sigma_3$だけで、どちらも基準の向き$\sigma_3$に結び付いています。

## 一様な磁場の中の歳差

一様な磁場$\boldsymbol B=B\boldsymbol b$（$\boldsymbol b$は単位ベクトル）の中で、スピンの向きの時間変化を調べます。時間に依存しない一様磁場を考え、状態が$\Psi(\boldsymbol x,t)=f(\boldsymbol x,t)\chi(t)$と軌道部分とスピン部分に分離されているとします。軌道ハミルトニアンは$\chi$に作用しないので、$\chi(t)$の時間発展は磁場との結合項だけで記述できます。以下、$\chi$に対応する偶部分代数の元を$\psi(t)$と書き、これについて

$$
\hbar\,\partial_t\psi\,\omega\sigma_3=-\frac{gq\hbar}{4m}\boldsymbol B\psi\sigma_3
$$

を解きます。右から$\sigma_3$を掛けると$\omega\sigma_3\sigma_3=\omega$で、$\omega$は中心の元なので、両辺に$-\omega$を掛けて

$$
\partial_t\psi=\frac{gq}{4m}\,\omega\boldsymbol B\,\psi
$$

を得ます。右から掛かる因子が消え、2ベクトル$\omega\boldsymbol B$が左から$\psi$を回す方程式になりました。

&&&fml スピンの歳差 [fml-precession]
$$
\psi(t)=e^{-\omega\boldsymbol b\,\Omega t/2}\,\psi(0),\qquad\Omega=-\frac{gqB}{2m}
$$

スピンの向きは$\boldsymbol s(t)=e^{-\omega\boldsymbol b\Omega t/2}\,\boldsymbol s(0)\,e^{\omega\boldsymbol b\Omega t/2}$で、$\boldsymbol b$を軸として角速度$\Omega$で回ります。
&&&

回転子$e^{-\omega\boldsymbol b\theta/2}$は、$\boldsymbol b$に垂直な面$\omega\boldsymbol b$の中で角$\theta$の回転を与えます（ローレンツ変換の記事における$e^{-\omega\sigma_3\theta/2}$が$x_1$軸を$x_2$軸へ回す向き）。スピンの向きが磁場のまわりを一定の角速度で回るこの運動を、**歳差**と呼びます。電子（$q=-e$）では$\Omega=geB/2m>0$で、$\boldsymbol B$の矢の先の側から見て反時計回りです。[[7shi-em5]]

&&&ex $\boldsymbol b=\sigma_3$の場合
$\boldsymbol s(0)=\sigma_1$から出発すると、$\boldsymbol s(t)=\sigma_1\cos\Omega t+\sigma_2\sin\Omega t$です。行列形式では$\Psi(t)=\operatorname{diag}(e^{-i\Omega t/2},e^{i\Omega t/2})\Psi(0)$となり、2つの成分の位相差が$\Omega t$の速さで進みます。$\Psi^\dagger\hat\sigma_1\Psi$と$\Psi^\dagger\hat\sigma_2\Psi$がこの位相差の余弦と正弦を与えます。ヘステネス形式では、この位相差の進行が、そのままスピンの向きの回転として読めます。
&&&

## 固有速度の回転との比較

ローレンツ変換の記事では、荷電粒子の固有速度$U=cR\gamma_0\tilde R$の回転子が$dR/d\tau=\frac q{2mc}FR$に従うことを見ました。電場がなく磁場だけがある場合、次のようになります。[[7shi-em5]]

$$
\frac{dR}{d\tau}=\frac q{2m}\,\omega\boldsymbol B\,R
$$

で、磁場は粒子の速度の向きを回します（サイクロトロン運動）。スピンの方程式は

$$
\partial_t\psi=\frac g2\cdot\frac q{2m}\,\omega\boldsymbol B\,\psi
$$

であり、同じ2ベクトル$\omega\boldsymbol B$が左から掛かる同じ形の方程式です。違いは係数$g/2$だけで、$g=2$なら、速さが$c$より十分小さい範囲（$\tau\approx t$）で、スピンの向きと速度の向きは同じ角速度で回ります。磁場はどちらに対しても、それに垂直な面の回転を生成する2ベクトルとして働いています。

# まとめ

パウリスピノルを$\operatorname{Cl}_{3,0}(\mathbb R)$の偶部分代数の元として扱いました。

- **対応**：列ベクトル$\Psi=(a_0+ia_3,-a_2+ia_1)^T$は偶部分代数の元$\psi=a_0+\sum_ka_k\,\omega\sigma_k$、すなわち四元数に対応し、$\Psi$は$\hat\psi$の第1列です。$\hat\sigma_k\Psi$は$\sigma_k\psi\sigma_3$に、$i\Psi$は$\psi\,\omega\sigma_3$に対応します。$\psi\mapsto\psi P$で、左イデアルによる表示と結ばれます。
- **虚数単位**：量子力学の$i$は、右から掛かる2ベクトル$\omega\sigma_3=\sigma_1\sigma_2$です。右乗算と可換な作用に対して、スカラーとして振る舞います。
- **観測量**：$\psi=\sqrt\rho\,R$と分解すると、$\rho=\psi\tilde\psi$が確率密度、$\psi\sigma_3\tilde\psi=\rho R\sigma_3\tilde R$がスピンの向きです。スピンの向きはホップ写像と同じ形の式で与えられ、ベクトルと2ベクトルの区別はグレードで付きます。
- **大域位相**：$\psi\mapsto\psi e^{\omega\sigma_3\alpha}$は、スピン軸まわりに枠$e_1,e_2$を回す操作で、確率密度とスピンの向きを変えません。
- **パウリ方程式**：$\hbar\,\partial_t\psi\,\omega\sigma_3=\frac1{2m}\sum_k\pi_k^2\psi+q\varphi\psi-\frac{gq\hbar}{4m}\boldsymbol B\psi\sigma_3$と書けます。スピン部分の時間発展は$\partial_t\psi=\frac{gq}{4m}\omega\boldsymbol B\psi$となり、スピンは角速度$-gqB/2m$で磁場のまわりを歳差運動します。$g=2$なら、固有速度の回転子の方程式と同じ形です。

&&& スピノルの対応
列ベクトル$\Psi=(a_0+ia_3,-a_2+ia_1)^T$に偶部分代数の元$\psi$を対応させると、作用は次のように対応します。
$$
\psi=a_0+\sum_ka_k\,\omega\sigma_k,\qquad\hat\sigma_k\Psi\leftrightarrow\sigma_k\psi\sigma_3,\qquad i\Psi\leftrightarrow\psi\,\omega\sigma_3
$$
&&&

&&& 確率密度とスピンの向き
$\psi=\sqrt\rho\,R$と分解すると、次のようになります。
$$
\rho=\psi\tilde\psi,\qquad\psi\sigma_3\tilde\psi=\rho\,\boldsymbol s,\qquad\boldsymbol s=R\sigma_3\tilde R
$$
&&&

&&& パウリ方程式
$$
\hbar\,\partial_t\psi\,\omega\sigma_3=\frac1{2m}\sum_k\pi_k^2\psi+q\varphi\psi-\frac{gq\hbar}{4m}\boldsymbol B\psi\sigma_3
$$
&&&

&&& 一様な磁場の中の歳差運動
一様な磁場の中で、スピンの向きの変化を与える磁場の項だけを残すと、次の式になります。
$$
\partial_t\psi=\frac{gq}{4m}\,\omega\boldsymbol B\,\psi
$$
&&&
