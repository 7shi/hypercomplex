電磁場と結合したディラック方程式の幾何学的定式化、局所ゲージ変換、および流れの保存則を扱います。

# 概要

前回は、自由粒子のディラック方程式を$\hbar D\psi\,\omega\sigma_3=mc\,\psi\gamma_0$と書きました。[[7shi-dirac3]]

本記事では、電荷$q$の粒子を電磁場の中に置きます。時空のポテンシャル$A$を使うと、行列形式の最小結合$i\hbar\partial_\mu\to i\hbar\partial_\mu-qa_\mu$は、左から$A$を掛ける項$-\frac qcA\psi$を加えることになります。[[7shi-em6]]

ポテンシャルのゲージ変換$A\mapsto A+D\chi$に対して、スピノルは右からの回転子$\psi\mapsto\psi e^{\omega\sigma_3\alpha}$で変換されます。これは大域位相を時空の点ごとに変えたもので、流れの向き$e_0$とスピンの向き$e_3$を変えずに、$e_1,e_2$をスピン軸のまわりに回します。[[7shi-dirac1]]

最後に、流れ$J=\psi\gamma_0\tilde\psi$が保存されること、右からの$\sigma_1$が電荷の符号を反転することを示します。

# 最小結合

## 行列形式

電磁気学の記法を使い、スカラーポテンシャル$\varphi$とベクトルポテンシャル$\boldsymbol A=\sum_kA_k\sigma_k$から、時空のベクトル

$$
A=\varphi\,\gamma_0+c\sum_kA_k\gamma_k
$$

を作ります。$F=D\wedge A$が電磁場で、$F=\boldsymbol E+\omega c\boldsymbol B$です（emシリーズの擬スカラー$i$を本シリーズでは$\omega$と書きます）。$A$の成分を$c$で割って添字を下げた$a_\mu=\gamma_\mu\cdot A/c$（空間成分の$A_k$と区別するため$a$と書きます）、すなわち

$$
a_0=\frac\varphi c,\qquad a_k=-A_k\quad(k=1,2,3)
$$

を使うと、$A=c\sum_\mu a_\mu\gamma^\mu$です。[[7shi-em6]]

古典力学では、正準運動量$\boldsymbol p$と、速度に結び付く運動学的運動量$\boldsymbol\pi$は、ベクトルポテンシャルの分だけずれた$\boldsymbol\pi=\boldsymbol p-q\boldsymbol A$で関係し、ハミルトニアンには$q\varphi$が加わります。パウリ方程式の$\hat\pi_k=-i\hbar\partial_k-qA_k$と$q\varphi$がこれにあたります。4元運動量で書けば、$i\hbar\partial_\mu$を$i\hbar\partial_\mu-qa_\mu$に替える操作で、**最小結合**と呼ばれます。行列形式のディラック方程式は

$$
\sum_\mu\hat\gamma^\mu\bigl(i\hbar\,\partial_\mu-qa_\mu\bigr)\Psi=mc\Psi
$$

となります。[[7shi-dirac1]]

## ヘステネス形式

&&&thm 電磁場の中のディラック方程式
スピノルの対応$\psi\mapsto\Psi$のもとで、上の行列形式の方程式は

$$
\hbar\,D\psi\,\omega\sigma_3-\frac qcA\psi=mc\,\psi\gamma_0
$$

と同値である。[[7shi-dirac3]]
&&&

&&&prf
作用の翻訳より、$\hat\gamma^\mu(-qa_\mu\Psi)$には$-q\gamma^\mu a_\mu\psi\gamma_0$が対応し、$\mu$について和を取ると$-\frac qcA\psi\gamma_0$である。自由粒子の場合と合わせて、行列形式の方程式は$\bigl(\hbar D\psi\,\omega\sigma_3-\frac qcA\psi\bigr)\gamma_0=mc\psi$と同値であり、右から$\gamma_0$を掛ければよい。[[7shi-dirac3]]
&&&

結合の項$-\frac qcA\psi$は、ベクトル$A$を左から掛けます。ローレンツ力$m\,dU/d\tau=\frac qcF\cdot U$と同じく係数は$q/c$で、ポテンシャル$A$の単位（$\varphi$と同じ）に合わせた形です。$A$は$D$と同じく左から作用するので、共変性はそのまま保たれます。一定の回転子$R$について、$\psi'(x)=R\psi(\tilde RxR)$、$A'(x)=RA(\tilde RxR)\tilde R$と置けば、$\frac qcA'\psi'=R\bigl(\frac qcA\psi\bigr)(\tilde RxR)$で、方程式の各項が同じ形で$R$を左に出します。[[7shi-em5]][[7shi-em6]][[7shi-dirac3]]

# ゲージ変換

## スピノルの変換

スカラー関数$\chi$によるゲージ変換$A\mapsto A+D\chi$が電磁場$F=D\wedge A$を変えないことを見ました。ディラック方程式は$F$でなく$A$を含むので、方程式を保つにはスピノルも変換する必要があります。[[7shi-em6]]

&&&thm ゲージ変換
$A\mapsto A+D\chi$と同時に

$$
\psi\mapsto\psi\,e^{\omega\sigma_3\alpha},\qquad\alpha=-\frac{q\chi}{\hbar c}
$$

と変換すると、電磁場の中のディラック方程式は保たれる。
&&&

&&&prf
$\psi'=\psi e^{\omega\sigma_3\alpha}$と置く。$\alpha$は位置の関数で、$e^{\omega\sigma_3\alpha}$の微分は$\partial_\mu\alpha\;\omega\sigma_3e^{\omega\sigma_3\alpha}$である。$\partial_\mu\alpha$はスカラーなので

$$
D\psi'=(D\psi)e^{\omega\sigma_3\alpha}+(D\alpha)\,\psi\,\omega\sigma_3e^{\omega\sigma_3\alpha}
$$

となる。右から$\omega\sigma_3$を掛けると、$\omega\sigma_3$は$e^{\omega\sigma_3\alpha}$と可換で$(\omega\sigma_3)^2=-1$だから

$$
\hbar D\psi'\,\omega\sigma_3=\bigl(\hbar D\psi\,\omega\sigma_3-\hbar(D\alpha)\psi\bigr)e^{\omega\sigma_3\alpha}
$$

である。結合の項は$-\frac qc(A+D\chi)\psi'=\bigl(-\frac qcA\psi-\frac qc(D\chi)\psi\bigr)e^{\omega\sigma_3\alpha}$で、$\hbar D\alpha=-\frac qcD\chi$により余分な項は打ち消し合う。質量の項は、$\gamma_0$が$\omega\sigma_3=\gamma_2\gamma_1$と可換なので$\psi'\gamma_0=\psi\gamma_0e^{\omega\sigma_3\alpha}$である。したがって方程式全体が右から$e^{\omega\sigma_3\alpha}$を掛けた形になり、解は解に移る。
&&&

行列形式では、$\psi e^{\omega\sigma_3\alpha}$は$e^{i\alpha}\Psi$なので、ゲージ変換は$\Psi\mapsto e^{-iq\chi/\hbar c}\Psi$です。位相を位置によって変えると、微分から余分な項$(D\alpha)\psi$が出ます。ポテンシャルのゲージ変換で加わる$D\chi$が、ちょうどその項を打ち消します。ポテンシャルを導入しない自由方程式は、定数位相の変換には不変ですが、一般の局所位相変換には不変ではありません。局所位相変換にも対応させるため、位相の微分から生じる項を補うポテンシャル$A$を導入します。この見方から電磁場との結合を導く考え方を、**ゲージ原理**と呼びます。

## 枠の回転としてのゲージ

この節では$\psi\tilde\psi\ne0$の領域を考え、前回の分解と枠を使います。ゲージ変換は右から掛かる回転子なので、枠$\psi\gamma_\mu\tilde\psi=\rho e_\mu$への影響が読み取れます。$e^{\omega\sigma_3\alpha}$は$\gamma_0$と$\gamma_3$を動かさず、$x_1x_2$平面の中で$\gamma_1,\gamma_2$を回すので

$$
e_0\mapsto e_0,\qquad e_3\mapsto e_3,\qquad
e_1\mapsto e_1\cos2\alpha-e_2\sin2\alpha,\qquad
e_2\mapsto e_1\sin2\alpha+e_2\cos2\alpha
$$

です。密度$\rho$と角$\beta$も変わりません（$\psi e^{\omega\sigma_3\alpha}\widetilde{\psi e^{\omega\sigma_3\alpha}}=\psi\tilde\psi$）。[[7shi-dirac2]]

流れの向き$e_0$とスピンの向き$e_3$はゲージによらない観測量です。残りの$e_1,e_2$は、$e_1,e_2$の張る面（スピン軸$e_3$に垂直な面）の中で点ごとに任意に回せる向きで、その向きそのものはゲージ依存であり、単独では観測量になりません。大域位相がスピン軸まわりの枠の回転として見えたのと同様に、ゲージ変換はその回転の角を点ごとに変えたもので、ポテンシャル$A$に現れたゲージの自由度が、スピノルの側では$e_1,e_2$の回転の自由度として現れます。[[7shi-dirac1]][[7shi-em6]]

emシリーズでは、ゲージ変換による$DA$の変化はスカラー部$D\cdot A$だけに現れ、$F$は変わらないことを見ました。本シリーズでは、同じ自由度がスピノルの右側、すなわち観測者によらない基準の側で、$\omega\sigma_3$の面の回転として働いています。スピノルの値に対する左からのローレンツ作用と右からの位相作用は、左右から掛かるので可換です。ただし、場のローレンツ変換は引数も変えるので、局所位相$\alpha$もスカラー場として引き戻して扱います。

# 流れの保存

前回は、$J$を確率の流れとして使うには保存が必要であることを述べました。ここで方程式から示します。行列形式の確率密度$\Psi^\dagger\Psi$は流れ$J=\psi\gamma_0\tilde\psi$の$\gamma_0$成分です。$J$は時空のベクトルで、ゲージ変換で変わりません。[[7shi-dirac3]]

&&&thm 流れの保存
電磁場の中のディラック方程式の解について

$$
D\cdot J=0,\qquad J=\psi\gamma_0\tilde\psi
$$

が成り立つ。
&&&

&&&prf
$D\cdot J$は$DJ$のスカラー部である。$D$を$\psi$と$\tilde\psi$の両方に作用させると

$$
\langle D(\psi\gamma_0\tilde\psi)\rangle_0=\langle(D\psi)\gamma_0\tilde\psi\rangle_0+\sum_\mu\langle\gamma^\mu\psi\gamma_0\partial_\mu\tilde\psi\rangle_0
$$

である。スカラー部は反転と巡回的な並べ替えで変わらないので、第2項は$\sum_\mu\langle\partial_\mu\psi\,\gamma_0\tilde\psi\gamma^\mu\rangle_0=\sum_\mu\langle\gamma^\mu\partial_\mu\psi\,\gamma_0\tilde\psi\rangle_0$となり、第1項に等しい。したがって$D\cdot J=2\langle(D\psi)\gamma_0\tilde\psi\rangle_0$である。

方程式より$D\psi=-\frac1\hbar\bigl(mc\,\psi\gamma_0+\frac qcA\psi\bigr)\omega\sigma_3$であり、$\gamma_0$と$\omega\sigma_3$は可換だから

$$
(D\psi)\gamma_0\tilde\psi=-\frac1\hbar\Bigl(mc\,\psi\,\omega\sigma_3\tilde\psi+\frac qcA\,\psi\gamma_0\omega\sigma_3\tilde\psi\Bigr)
$$

となる。$\psi\,\omega\sigma_3\tilde\psi$は反転で符号が変わるので、スカラー部を持たない。$\gamma_0\omega\sigma_3=\gamma_0\gamma_2\gamma_1$は3ベクトルで反転で符号が変わるから、$\psi\gamma_0\omega\sigma_3\tilde\psi$は反転で符号が変わる奇数グレードの元、すなわち3ベクトルである。ベクトル$A$と3ベクトルの積はグレード2と4の成分しか持たず、スカラー部は$0$である。したがって両項のスカラー部は$0$であり、$D\cdot J=0$となる。
&&&

証明で使ったのは、右から掛かる因子の形だけです。$\omega\sigma_3$が2ベクトルで、$\gamma_0\omega\sigma_3$が3ベクトルであることから、質量の項と結合の項がスカラー部に寄与しないことが従います。$J=J^0\gamma_0+\sum_kJ^k\gamma_k$と成分に分け、$J^0=\Psi^\dagger\Psi$、$\boldsymbol j=c\sum_kJ^k\sigma_k$（確率の流れ）と置くと、$D\cdot J=0$は$\partial_tJ^0+\nabla\cdot\boldsymbol j=0$の形の連続の式です。

電磁気学では、電流$J_{\mathrm{em}}=c\rho_{\mathrm{em}}\gamma_0+\sum_kJ_k\gamma_k$の保存$D\cdot J_{\mathrm{em}}=0$を、マクスウェル方程式から$D\cdot(D\cdot F)=0$として導きました。ここでは保存則がディラック方程式そのものから出ます。確率密度に電荷$q$を掛けたものが電荷密度（$\rho_{\mathrm{em}}=qJ^0$、$\boldsymbol J_{\mathrm{em}}=q\boldsymbol j$）なので、$cqJ$が$J_{\mathrm{em}}$と同じ形の電流を与え、同じ保存則を満たします。ディラック方程式の解が作る電磁場（マクスウェル方程式との連立）は、本記事では扱いません。[[7shi-em4]]

# 電荷の反転

自由粒子の方程式の解に右から$\sigma_1$を掛けると、正のエネルギーの解が負のエネルギーの解に移ることを見ました。電磁場がある場合、同じ操作は電荷の符号を変えます。[[7shi-dirac3]]

&&&prop 電荷の反転
$\psi$が電荷$q$の方程式$\hbar D\psi\,\omega\sigma_3-\frac qcA\psi=mc\,\psi\gamma_0$の解なら、$\psi\sigma_1$は電荷$-q$の方程式の解である。流れは$\psi\sigma_1\gamma_0\widetilde{\psi\sigma_1}=J$で変わらない。
&&&

&&&prf
$\sigma_1$は$\omega\sigma_3$とも$\gamma_0$とも反可換なので

$$
\hbar D(\psi\sigma_1)\,\omega\sigma_3+\frac qcA\psi\sigma_1-mc\,\psi\sigma_1\gamma_0
=-\Bigl(\hbar D\psi\,\omega\sigma_3-\frac qcA\psi-mc\,\psi\gamma_0\Bigr)\sigma_1=0
$$

である。流れについては、$\tilde\sigma_1=-\sigma_1$と$\sigma_1\gamma_0\sigma_1=-\gamma_0$より$\psi\sigma_1\gamma_0(-\sigma_1)\tilde\psi=\psi\gamma_0\tilde\psi$である。
&&&

$\psi\sigma_1$の流れ$J$は同じで電荷が$-q$なので、電流$cqJ$は符号が反転します。行列形式では、この操作は複素共役を含む変換（荷電共役）として書かれます。ヘステネス形式では複素共役という操作は現れず、右からの$\sigma_1$が、$i$の役割の$\omega\sigma_3$と反可換であることによって、$i\mapsto-i$の効果を担っています。その物理的な解釈（反粒子）には、本記事では立ち入りません。

# まとめ

ディラック方程式を電磁場と結合させ、ゲージ変換と保存則を扱いました。

- **最小結合**：$A=\varphi\gamma_0+c\sum_kA_k\gamma_k$により、電磁場の中のディラック方程式は$\hbar D\psi\,\omega\sigma_3-\frac qcA\psi=mc\,\psi\gamma_0$です。$A$は$D$と同じく左から作用し、共変性は保たれます。[[7shi-em6]]
- **ゲージ変換**：$A\mapsto A+D\chi$に対して$\psi\mapsto\psi e^{\omega\sigma_3\alpha}$、$\alpha=-q\chi/\hbar c$で方程式は保たれます。行列形式の$\Psi\mapsto e^{-iq\chi/\hbar c}\Psi$です。
- **枠の回転**：ゲージ変換は右からの回転子で、流れ$e_0$、スピン$e_3$、$\rho$、$\beta$を変えず、$e_1,e_2$をスピン軸のまわりに点ごとに回します。大域位相を局所化したもので、ポテンシャルのゲージの自由度にあたります。
- **流れの保存**：$J=\psi\gamma_0\tilde\psi$について$D\cdot J=0$が方程式から従います。右から掛かる因子が2ベクトル$\omega\sigma_3$と3ベクトル$\gamma_0\omega\sigma_3$であることが要点です。
- **電荷の反転**：右から$\sigma_1$を掛けると電荷$q$の解が電荷$-q$の解に移り、流れは変わりません。

&&& 電磁場の中のディラック方程式
$$
\hbar D\psi\,\omega\sigma_3-\frac qcA\psi=mc\,\psi\gamma_0
$$
&&&

&&& ゲージ変換
$A\mapsto A+D\chi$に対して、次の変換で方程式は保たれる。
$$
\psi\mapsto\psi\,e^{\omega\sigma_3\alpha},\qquad\alpha=-\frac{q\chi}{\hbar c}
$$
&&&

&&& 流れの保存
電磁場の中のディラック方程式の解について、次が成り立つ。
$$
D\cdot J=0,\qquad J=\psi\gamma_0\tilde\psi
$$
&&&
