時空代数の偶部分代数におけるディラックスピノルの幾何学的構造を調べ、ローレンツ変換の二重被覆と観測量の枠分解を解説します。

# 概要

前回の記事では、パウリスピノルを$\operatorname{Cl}_{3,0}(\mathbb R)$の偶部分代数（四元数）の元$\psi$として扱い、空間の回転が左からの回転子、位相が右からの$e^{\omega\sigma_3\alpha}$として作用することを見ました。本記事では、これを時空代数$\operatorname{Cl}_{1,3}(\mathbb R)$に上げます。時空代数の偶部分代数は$\operatorname{Cl}_{3,0}(\mathbb R)\cong M_2(\mathbb C)$と同型で、双四元数にあたります。その回転子の群$\operatorname{Spin}^+(1,3)$は$\operatorname{SL}(2,\mathbb C)$と同型で、ローレンツ変換の群$\operatorname{SO}^+(1,3)$を二重に覆います。偶部分代数の一般の元をディラックスピノルとし、密度・角$\beta$・回転子の積$\psi=\sqrt\rho\,e^{\omega\beta/2}R$に分解します。回転子が基底を回した枠$e_\mu=R\gamma_\mu\tilde R$から、流れとスピンの向きが読み取れます。[[7shi-dirac1]]

# 時空代数の偶部分代数

## 記号

電磁気学の記事と同じく、時空代数$\operatorname{Cl}_{1,3}(\mathbb R)$の生成元を$\gamma_0,\gamma_1,\gamma_2,\gamma_3$（$\gamma_0^2=1$、$\gamma_k^2=-1$）、相対ベクトルを$\sigma_k=\gamma_k\gamma_0$とします。擬スカラーは前回の記事と同じく$\omega$と書きます。

$$
\omega=\gamma_0\gamma_1\gamma_2\gamma_3=\sigma_1\sigma_2\sigma_3
$$

$\omega^2=-1$で、$\omega$は偶部分代数の元とは可換、ベクトルとは反可換です。$\sigma_k$は偶部分代数の中で$\operatorname{Cl}_{3,0}(\mathbb R)$の生成元の関係を満たし、$\operatorname{Cl}_{1,3}^0(\mathbb R)\cong\operatorname{Cl}_{3,0}(\mathbb R)$です。空間の面を表す2ベクトルは、次のように書けます。[[7shi-em4]][[7shi-dirac1]]

$$
\omega\sigma_1=\gamma_3\gamma_2,\qquad\omega\sigma_2=\gamma_1\gamma_3,\qquad\omega\sigma_3=\gamma_2\gamma_1
$$

です。前回の記事で右から掛かる虚数単位だった$\omega\sigma_3$は、時空代数では$x_1x_2$平面の2ベクトル$\gamma_2\gamma_1$です。

反転$\tilde X$は、時空代数の基底の積の順序を逆にする操作です。偶部分代数の元に対しては、スカラーと擬スカラーを変えず、2ベクトルの符号を変えます。注意が要るのは$\sigma_k$で、$\operatorname{Cl}_{3,0}(\mathbb R)$ではベクトルなので反転で変わりませんが、時空代数では2ベクトル$\gamma_k\gamma_0$なので$\tilde\sigma_k=\gamma_0\gamma_k=-\sigma_k$です。$1$と$\omega\sigma_k$については、どちらの代数で反転しても同じ結果になります。前回の記事における$\tilde\psi$は$1,\omega\sigma_k$の実係数の結合にだけ使ったので、以下ではそのまま時空代数の反転として読めます。[[7shi-dirac1]]

## 双四元数と行列

$\sigma_k\mapsto\hat\sigma_k$（パウリ行列）、$\omega\mapsto iI$により、偶部分代数は$M_2(\mathbb C)$と同型です。双四元数の記事では、$M_2(\mathbb C)$を四元数の係数を複素数に広げた双四元数として扱いました。偶部分代数は実8次元で、その元は双四元数1つにあたります。[[7shi-bq]]

&&&prop 反転と行列 [prop-rev]
偶部分代数の元$X$について、次が成り立ちます。

1. $\tilde X$は余因子行列$\operatorname{adj}\hat X$に対応します。
2. $X\tilde X$はスカラーと擬スカラーの和で、$\det\hat X$に対応します（$\omega\leftrightarrow i$）。
3. $\gamma_0\tilde X\gamma_0$はエルミート共役$\hat X^\dagger$に対応します。
&&&

&&&prf
反転と$\operatorname{adj}$はともに積の順序を逆にする実線形写像である（$2\times2$行列では$\operatorname{adj}(AB)=\operatorname{adj}B\operatorname{adj}A$）。$\tilde\sigma_k=-\sigma_k$であり、跡が$0$の$2\times2$行列では$\operatorname{adj}M=(\operatorname{tr}M)I-M=-M$だから$\operatorname{adj}\hat\sigma_k=-\hat\sigma_k$である。偶部分代数は$\sigma_k$で生成されるので、2つの写像は一致する。$\hat X\operatorname{adj}\hat X=(\det\hat X)I$より2が従う。$X\tilde X$は反転で変わらないので、偶部分代数の中ではグレード0と4の成分だけを持つ。

3について、$\gamma_0\sigma_k\gamma_0=\gamma_0\gamma_k=-\sigma_k$なので、$X\mapsto\gamma_0\tilde X\gamma_0$は積の順序を逆にする実線形写像で、$\sigma_k$を$\sigma_k$に移す。エルミート共役も積の順序を逆にする実線形写像で、$\hat\sigma_k$を$\hat\sigma_k$に移すから、2つの写像は一致する。
&&&

前回の記事では$\operatorname{Cl}_{3,0}(\mathbb R)$の反転がエルミート共役に対応しました。時空代数の反転は、$\operatorname{Cl}_{3,0}(\mathbb R)$で見るとクリフォード共役（ベクトルと2ベクトルの符号を変える操作）で、行列では余因子行列です。エルミート共役を得るには、さらに$\gamma_0$で挟む必要があります。エルミート共役は時間軸$\gamma_0$を選んではじめて定まる操作であり、観測者によらない反転とは区別されます。[[7shi-dirac1]]

# 回転子の群

## ベクトルへの作用と行列

ローレンツ変換の記事では、偶部分代数の元$R$で$R\tilde R=1$を満たすものを回転子と呼び、$x\mapsto Rx\tilde R$がローレンツ変換になることを見ました。この作用を行列で書きます。時空のベクトル$x$に$\gamma_0$を掛けるとパラベクトル$x\gamma_0=x_0+\boldsymbol x$になり、その行列表示を考えます。[[7shi-em4]][[7shi-em5]]

$$
\hat X=x_0I+\sum_kx_k\hat\sigma_k=\begin{pmatrix}x_0+x_3&x_1-ix_2\\x_1+ix_2&x_0-x_3\end{pmatrix}
$$

はエルミート行列で、$\det\hat X=x_0^2-|\boldsymbol x|^2=x^2$です。

&&&fml 回転子の作用の行列表示
$$
(Rx\tilde R)\gamma_0\ \longleftrightarrow\ \hat R\hat X\hat R^\dagger
$$
&&&

&&&prf
$\gamma_0^2=1$より$(Rx\tilde R)\gamma_0=R(x\gamma_0)(\gamma_0\tilde R\gamma_0)$であり、[反転と行列](#prop-rev)の3より$\gamma_0\tilde R\gamma_0$は$\hat R^\dagger$に対応する。
&&&

エルミート行列$\hat X$を$\hat R\hat X\hat R^\dagger$に移す変換は、相対論的な量子力学の教科書でスピノルとローレンツ変換を結ぶ標準的な式です。$\det(\hat R\hat X\hat R^\dagger)=|\det\hat R|^2\det\hat X$なので、$|\det\hat R|=1$なら$x^2$が保たれます。

## Spin(1,3)とSL(2,C)

&&&def $\operatorname{Spin}^+(1,3)$
時空代数の回転子の全体

$$
\operatorname{Spin}^+(1,3)=\{R\in\operatorname{Cl}_{1,3}^0(\mathbb R)\mid R\tilde R=1\}
$$

を$\operatorname{Spin}^+(1,3)$と書きます。
&&&

$R\tilde R=1$なら$R^{-1}=\tilde R$で、$\widetilde{RS}=\tilde S\tilde R$より回転子の積は回転子なので、$\operatorname{Spin}^+(1,3)$は群です。一方、空間の向きと時間の向きを保つローレンツ変換（行列式$1$で、未来向きの時間的ベクトルを未来向きに移す線形変換）の群を$\operatorname{SO}^+(1,3)$と書きます。

&&&thm 二重被覆
1. $R\mapsto\hat R$は$\operatorname{Spin}^+(1,3)$から$\operatorname{SL}(2,\mathbb C)$への群同型です。
2. $R\mapsto(x\mapsto Rx\tilde R)$は$\operatorname{Spin}^+(1,3)$から$\operatorname{SO}^+(1,3)$への全射な準同型で、核は$\{1,-1\}$です。
&&&

&&&prf
1：[反転と行列](#prop-rev)の2より、$R\tilde R=1$は$\det\hat R=1$と同値である。偶部分代数から$M_2(\mathbb C)$への写像は代数の同型だから、群の同型を与える。

2（核）：すべてのベクトル$x$について$Rx\tilde R=x$なら、$R\tilde R=1$より$Rx=xR$であり、$R$はすべての生成元と可換である。偶部分代数で$\sigma_k=\gamma_k\gamma_0$のすべてと可換な元は、$M_2(\mathbb C)$の中心にあたる$a+b\omega$の形である。$\omega$はベクトルと反可換だから$b=0$で、$R\tilde R=a^2=1$より$R=\pm1$である。

2（像）：$Rx\tilde R$の像が$\operatorname{SO}^+(1,3)$に入ることは次のように分かる。$R\mapsto\hat R$は連結な群$\operatorname{SL}(2,\mathbb C)$と同型なので、$R$は$1$と連続につながり、行列式と時間の向きは連続的に変わらない。

2（全射）：$L\in\operatorname{SO}^+(1,3)$を取り、$u=L(\gamma_0)$と置く。$u$は未来向きの単位時間的ベクトルで、$u\cdot\gamma_0\ge1$である。

$$
B=\frac{1+u\gamma_0}{\sqrt{2(1+u\cdot\gamma_0)}}
$$

と置くと、$u\gamma_0+\gamma_0u=2u\cdot\gamma_0$より$B\tilde B=\frac{(1+u\gamma_0)(1+\gamma_0u)}{2(1+u\cdot\gamma_0)}=1$であり、$u\gamma_0u=2(u\cdot\gamma_0)u-\gamma_0$より

$$
(1+u\gamma_0)\gamma_0(1+\gamma_0u)=(\gamma_0+u)(1+\gamma_0u)=2u+u\gamma_0u+\gamma_0=2(1+u\cdot\gamma_0)u
$$

となるので$B\gamma_0\tilde B=u$である。$L'(x)=\tilde BL(x)B$は$\gamma_0$を動かさないローレンツ変換で、$\gamma_0$に直交する空間$\gamma_1,\gamma_2,\gamma_3$の向きを保つ回転を与える。空間の回転は$\gamma_0$と可換な回転子$V=e^{-\omega\hat{\boldsymbol n}\theta/2}$で$x\mapsto Vx\tilde V$と書ける。したがって$L(x)=BVx\tilde V\tilde B$であり、$R=BV$が$L$を与える。[[7shi-em5]][[7shi-lie3]]
&&&

証明の中の$B$は、$\gamma_0$を$u$に移すブーストの回転子です。ブースト$e^{\sigma_1\eta/2}$は、$u=\cosh\eta\,\gamma_0+\sinh\eta\,\gamma_1$の場合の$B$にあたります。ローレンツ変換はブーストと空間の回転の合成に分解され、行列ではエルミートで正定値の行列とユニタリ行列の積（極分解）にあたります。[[7shi-em5]]

## 空間の回転への制限

回転子のうち$\gamma_0$と可換なもの$V$は、時間軸を動かさない空間の回転を与えます。$\gamma_0\tilde V\gamma_0=\tilde V=V^{-1}$なので、[反転と行列](#prop-rev)の3より$\hat V^\dagger=\hat V^{-1}$で、$\hat V$はユニタリです。$\gamma_0$と可換な偶部分代数の元は$1,\omega\sigma_1,\omega\sigma_2,\omega\sigma_3$の実係数の結合で、パウリスピノルの空間$\operatorname{Cl}_{3,0}^0(\mathbb R)\cong\mathbb H$そのものです。回転子を$\gamma_0$と可換なものに制限すると、二重被覆は$\operatorname{SU}(2)\to\operatorname{SO}(3)$に戻ります。[[7shi-dirac1]][[7shi-lie3]][[7shi-cover]]

二重被覆の2対1の性質は、空間の回転で見るのが分かりやすいです。$x_3$軸まわりの角$\theta$の回転子$e^{-\omega\sigma_3\theta/2}$は、$\theta=2\pi$で$-1$になります。ベクトルは$2\pi$回転で元に戻りますが、回転子の経路は$1$から$-1$へつながります。$\operatorname{SO}^+(1,3)$の中の閉じた経路（$2\pi$回転）が、$\operatorname{Spin}^+(1,3)$では閉じない経路に持ち上がります。以前の記事で見た$\operatorname{SO}(3)$の縮められないループと同じ構造が、ブーストを含む非コンパクトな群でもそのまま現れます。極分解により、$\operatorname{SL}(2,\mathbb C)$は$\operatorname{SU}(2)$とエルミートで正定値な行列式$1$の行列の集まり（ブーストの全体で、$\mathbb R^3$と同じ形の空間）の直積の形をしており、縮められないループはブーストの側には現れません。[[7shi-cover]]

# ディラックスピノル

## スピノルの変換

前回の記事では、パウリスピノル$\psi$に左から空間の回転子を掛けました。時空代数では、ローレンツ変換の回転子を左から掛けます。[[7shi-dirac1]]

&&&def ディラックスピノル
偶部分代数$\operatorname{Cl}_{1,3}^0(\mathbb R)$の元$\psi$を**ディラックスピノル**と呼び、回転子$R$によるローレンツ変換は$\psi\mapsto R\psi$として作用するものとします。
&&&

偶部分代数は実8次元で、行列形式のディラックスピノル（複素4成分）と同じ自由度を持ちます。ベクトルが両側から$x\mapsto Rx\tilde R$の変換を受けるのに対し、スピノルは片側からだけ変換を受けます。$R$と$-R$はベクトルに同じ変換を与えますが、スピノルには符号の違う変換を与えます。スピノルはローレンツ変換の群$\operatorname{SO}^+(1,3)$でなく、それを二重に覆う$\operatorname{Spin}^+(1,3)$の作用を受ける量です。

右側には何も作用しないので、スピノルに右から掛ける操作はローレンツ変換と可換です。前回の記事で虚数単位を右からの$\omega\sigma_3$として扱えたのと同じく、右側は観測者によらない構造を置く場所として使えます。[[7shi-dirac1]]

## 密度・角・回転子への分解

&&&prop スピノルの分解 [prop-decomp]
ディラックスピノル$\psi$について$\psi\tilde\psi$はスカラーと擬スカラーの和です。$\psi\tilde\psi\ne0$なら、$\rho>0$と実数$\beta$により$\psi\tilde\psi=\rho e^{\omega\beta}$と書け

$$
\psi=\sqrt\rho\,e^{\omega\beta/2}R,\qquad R\tilde R=1
$$

と分解されます。
&&&

&&&prf
$\psi\tilde\psi$は偶部分代数の元で反転で変わらないから、グレード0と4の成分だけを持ち、$\psi\tilde\psi=a+b\omega$と書ける。$(a,b)\ne(0,0)$なら$\rho=\sqrt{a^2+b^2}$、$a=\rho\cos\beta$、$b=\rho\sin\beta$と置くと、$\omega^2=-1$より$\psi\tilde\psi=\rho e^{\omega\beta}$である。$R=e^{-\omega\beta/2}\psi/\sqrt\rho$と置くと、$\tilde\omega=\omega$で$\omega$は偶部分代数の元と可換だから

$$
R\tilde R=\frac1\rho e^{-\omega\beta/2}\psi\tilde\psi e^{-\omega\beta/2}=e^{-\omega\beta}e^{\omega\beta}=1
$$

となる。
&&&

行列では$\psi\tilde\psi$は$\det\hat\psi$にあたり、$\psi\tilde\psi\ne0$は$\hat\psi$が可逆であることを意味します。自由度は、密度$\rho$が1、角$\beta$が1、回転子$R$が6（ローレンツ変換の自由度）で、合わせて8です。パウリスピノル$\psi=\sqrt\rho\,R$と比べると、回転子が空間の回転からローレンツ変換に広がり、新たに角$\beta$が加わっています。[[7shi-dirac1]]

ローレンツ変換$\psi\mapsto L\psi$では$L\psi\,\widetilde{L\psi}=L\psi\tilde\psi\tilde L=\psi\tilde\psi$なので、$\rho$と$\beta$は変わらず、回転子だけが$R\mapsto LR$と変わります。

&&&rem 角$\beta$
$\beta$は観測者によらない量ですが、その物理的な意味は定まっていません。文献ではイヴォン＝高林角（Yvon–Takabayashi angle）と呼ばれます。たとえば$\psi=1$なら$\beta=0$、$\psi=\sigma_1$なら$\psi\tilde\psi=\sigma_1(-\sigma_1)=-1$で$\beta=\pi$です。本シリーズでは$\beta$の解釈には立ち入りません。
&&&

## 枠と観測量

回転子$R$は時空の正規直交基底$\gamma_\mu$を回して、新しい正規直交基底$e_\mu=R\gamma_\mu\tilde R$を作ります。スピノルからは、この枠が直接読み取れます。

&&&fml スピノルの作る枠
$$
\psi\gamma_\mu\tilde\psi=\rho\,e_\mu,\qquad e_\mu=R\gamma_\mu\tilde R\qquad(\mu=0,1,2,3)
$$
&&&

&&&prf
$\omega$はベクトルと反可換なので$e^{\omega\beta/2}e_\mu=e_\mu e^{-\omega\beta/2}$であり

$$
\psi\gamma_\mu\tilde\psi=\rho\,e^{\omega\beta/2}R\gamma_\mu\tilde Re^{\omega\beta/2}=\rho\,e^{\omega\beta/2}e_\mu e^{\omega\beta/2}=\rho\,e_\mu
$$

となる。
&&&

角$\beta$は枠に現れません。$e_0$は未来向きの単位時間的ベクトルで、$\rho e_0=\psi\gamma_0\tilde\psi$は時空のベクトルとして変換されます。ディラック方程式の理論では、この$J=\psi\gamma_0\tilde\psi$を確率の流れとして使い、$J\cdot\gamma_0$をその観測者の見る確率密度とします。$e_3$はスピンの向きで、$\rho e_3=\psi\gamma_3\tilde\psi$です。スピンの向き$e_3$は流れ$e_0$と直交します。

パウリスピノルの$\psi\sigma_3\tilde\psi$との関係を確かめます。$\psi$が$\gamma_0$と可換（パウリスピノル）なら、$\tilde\psi$も$\gamma_0$と可換なので

$$
\psi\gamma_0\tilde\psi=\psi\tilde\psi\,\gamma_0=\rho\,\gamma_0,\qquad
\psi\sigma_3\tilde\psi=\psi\gamma_3\gamma_0\tilde\psi=(\psi\gamma_3\tilde\psi)\,\gamma_0=\rho\,e_3\gamma_0
$$

です。パウリスピノルは流れが$\gamma_0$の向き、つまりこの観測者に対して静止したスピノルで、以前に求めたスピンの向きは、時空のベクトル$e_3$を相対ベクトルとして読んだものです。一般のディラックスピノルでは、流れ$e_0$が$\gamma_0$から傾き、スピンの向き$e_3$もそれに合わせて傾きます。[[7shi-dirac1]]

ローレンツ変換$\psi\mapsto L\psi$のもとで、$\psi\gamma_\mu\tilde\psi\mapsto L(\psi\gamma_\mu\tilde\psi)\tilde L$で、観測量はベクトルとして変換されます。スピノルの片側の作用から、両側の作用を受ける量が2次式として作られるという構造は、スピノルの外積と同じです。[[7shi-lie3]]

## 右からの作用

スピノルに右から回転子$S$を掛けると、$\psi\gamma_\mu\tilde\psi\mapsto\psi(S\gamma_\mu\tilde S)\tilde\psi$で、枠の基準にしている基底$\gamma_\mu$が取り替えられます。とくに$S=e^{\omega\sigma_3\alpha}$は$\gamma_0$と$\gamma_3$を動かさず、$\gamma_1,\gamma_2$を$x_1x_2$平面の中で回します。したがって流れ$e_0$とスピンの向き$e_3$は変わらず、$e_1,e_2$だけがスピンの向きのまわりに回ります。大域位相の役割が、時空でもそのまま保たれています。[[7shi-dirac1]]

左からの作用と右からの作用は次のように分担されます。

| | 作用 | 観測量への影響 |
|---|---|---|
| 左から$R$ | ローレンツ変換（観測者によらない時空の回転） | 枠$e_\mu$を回す |
| 右から$S$ | 基準の基底$\gamma_\mu$の取り替え | 基準の取り方を変える |

右から掛けるときに$\gamma_0$や$\sigma_3$といった特定の基底が現れますが、それは観測者の選択ではなく、スピノルの成分を読むための基準の選択です。左からのローレンツ変換と可換なので、どの観測者に対しても同じ基準が使えます。

# まとめ

時空代数の偶部分代数の元としてディラックスピノルを扱いました。

- **偶部分代数**：$\operatorname{Cl}_{1,3}^0(\mathbb R)\cong\operatorname{Cl}_{3,0}(\mathbb R)\cong M_2(\mathbb C)$で、ディラックスピノルは双四元数1つにあたります。時空の反転は余因子行列に、$\gamma_0$で挟んだ反転はエルミート共役に対応し、$\psi\tilde\psi$は$\det\hat\psi$です。
- **二重被覆**：回転子の群$\operatorname{Spin}^+(1,3)$は$\operatorname{SL}(2,\mathbb C)$と同型で、$x\mapsto Rx\tilde R$は$\operatorname{SO}^+(1,3)$への2対1の全射です。行列では$\hat X\mapsto\hat R\hat X\hat R^\dagger$です。$\gamma_0$と可換な回転子に制限すると$\operatorname{SU}(2)\to\operatorname{SO}(3)$に戻ります。
- **分解**：$\psi\tilde\psi=\rho e^{\omega\beta}$から$\psi=\sqrt\rho\,e^{\omega\beta/2}R$と分解され、ローレンツ変換$\psi\mapsto L\psi$は回転子だけを変えます。
- **枠**：$\psi\gamma_\mu\tilde\psi=\rho\,e_\mu$で、$e_0$が流れの向き、$e_3$がスピンの向きです。パウリスピノルは$e_0=\gamma_0$の場合にあたり、$\psi\sigma_3\tilde\psi$は$\rho e_3\gamma_0$です。[[7shi-dirac1]]
- **右からの作用**：右から掛ける回転子は基準の基底を取り替え、$e^{\omega\sigma_3\alpha}$は流れとスピンの向きを変えずに$e_1,e_2$を回します。

&&& 偶部分代数と回転子の群
$$
\operatorname{Cl}_{1,3}^0(\mathbb R)\cong\operatorname{Cl}_{3,0}(\mathbb R)\cong M_2(\mathbb C),\qquad\operatorname{Spin}^+(1,3)\cong\operatorname{SL}(2,\mathbb C)
$$
&&&

&&& スピノルの分解
$\psi\tilde\psi\ne0$なら、$\rho>0$と実数$\beta$により次のように書けます。
$$
\psi\tilde\psi=\rho e^{\omega\beta},\qquad\psi=\sqrt\rho\,e^{\omega\beta/2}R,\qquad R\tilde R=1
$$
&&&

&&& 枠
分解の回転子$R$により、次が成り立ちます。
$$
\psi\gamma_\mu\tilde\psi=\rho\,e_\mu,\qquad e_\mu=R\gamma_\mu\tilde R\qquad(\mu=0,1,2,3)
$$
&&&
