平坦なディラック作用素を動径方向と球面方向に分けて球面上のディラック作用素を取り出し、その固有値と、核が消えることを示します。

# 概要

以前の記事では、クリフォード代数$\operatorname{Cl}_{n,0}(\mathbb R)$に値を取る関数に作用するディラック作用素$D=\sum_ae_a\partial_a$を扱い、$D^2=\Delta$（$n$次元のラプラシアン）となることを見ました。$D$は、スカラーの関数の上では書けないラプラシアンの平方根を、値の空間をクリフォード代数に広げることで実現したものです。同じ記事では、球面モノジェニックス（球面調和関数のモノジェニック版）を扱わずに残しました。[[7shi-cla5]]

以前の記事では、アティヤ＝シンガーの指数定理が閉多様体（コンパクトで境界のない多様体）上の作用素についての定理であることを紹介し、平坦な$D$を、その主表象$\boldsymbol\xi\mapsto\boldsymbol\xi\cdot$を理解する局所モデルとして見ました。平坦な空間$\mathbb R^n$は閉多様体ではありません。閉多様体の最も簡単な例は球面です。[[7shi-kth1]]

本記事で扱う問いは次のとおりです。

> 球面の上で、ラプラシアンの平方根にあたる1階の作用素はどう書けるか。その固有値と核はどうなるか。

曲がった空間のディラック作用素は、通常、各点に正規直交な接ベクトルの組（正規直交枠）を取り、その回転を補正する項（スピン接続）を加えて組み立てます。本記事ではこの構成から出発せず、平坦な$D$を球面の方向に分解して球面上の作用素を取り出します。球面は平坦な空間の中にあるので、周囲の空間から誘導される構成を使えば、スピノルは偶部分代数に値を取る球面上の1価の関数として書け、スピン接続は、分解で現れる定数のずれとして読めます。

前提は次のとおりです。

- **一般次元のディラック作用素**：$\operatorname{Cl}_{n,0}(\mathbb R)$（$e_a^2=1$、添字は$0$から$n-1$）、$\boldsymbol x=\sum_ax_ae_a$、$D=\sum_ae_a\partial_a$、$D^2=\Delta$、モノジェニック関数の定義、$h_l=e_0e_l$が$h_l^2=-1$を満たして互いに反交換すること[[7shi-cla5]]
- **偶部分代数**：$\operatorname{Cl}_{n,0}^0(\mathbb R)\cong\operatorname{Cl}_{0,n-1}(\mathbb R)$[[7shi-clif1]]
- **主表象**：平坦な$D$の$\partial_a$を$\xi_a$に置き換えた$\boldsymbol\xi\mapsto\boldsymbol\xi\cdot$[[7shi-kth1]]
- **二重被覆**：回転を1周させると、それを表す単位四元数は$-1$までしか進まないこと[[7shi-cover]]

本記事は次の順に進みます。

1. $\boldsymbol xD$を動径方向の微分と球面方向の作用素$Γ$に分け、$Γ$と球面のラプラシアンの関係を求めます。
2. 次数$k$の球面モノジェニックスが$Γ$の固有関数になり、固有値が$-k$と$k+n-1$の2系列に分かれることを示します。
3. 接ベクトルを掛ける演算と、それと両立する共変微分を定め、球面上のディラック作用素が$\frac{n-1}2-Γ$となることを導きます。固有値は$\pm\bigl(k+\frac{n-1}2\bigr)$に揃います。
4. 円周と2次元の球面で、正規直交枠に移した形と比べ、定数のずれがスピン接続の項にあたることと、半整数の固有値が現れる仕組みを確かめます。
5. リヒネロビッチの公式を示し、$n\ge3$で核が$0$になることを導きます。

証明するのは、$Γ$の性質、球面モノジェニックスの固有値と空間の次元、共変微分が接ベクトルの掛け算と両立すること、ディラック作用素の形、リヒネロビッチの公式と核の消失です。フィッシャー分解と固有関数が完全系をなすこと、一般の$n$で本記事の共変微分が標準的なスピン接続に一致すること、球面のスピン構造についての事実、単位球面のスカラー曲率の値、フリードリヒの評価は主張として使います。球面上のディラック作用素は$L^2$上では非有界な作用素なので、前回の記事の有界なフレドホルム作用素の枠組みに当てはめるには、定義域などの準備が要ります。本記事では核だけを調べます。[[7shi-kth5]]

# 平坦なディラック作用素の分解

## 動径方向と球面方向

以下、$n\ge2$とし、以前の記事の設定を使います。[[7shi-cla5]]

$\mathbb R^n$の原点の周りの微分作用素を、原点からの距離を変える方向と、距離を保つ方向に分けます。前者は**オイラー作用素**$E=\sum_ax_a\partial_a$で、$r=|\boldsymbol x|$として$E=r\partial_r$です。次数$k$の同次多項式$P$に対して$EP=kP$となります。後者は、座標平面$e_ae_b$の中の回転の方向の微分です。

$$
L_{ab}=x_a\partial_b-x_b\partial_a
$$

$L_{ab}$は、ベクトル場$\boldsymbol v_{ab}=x_ae_b-x_be_a$の方向の微分です。$\boldsymbol v_{ab}$は$\boldsymbol x$と直交するので、$L_{ab}$は原点を中心とする各球面に接する方向の微分であり、球面上の関数に作用します。

$D$に左から$\boldsymbol x$を掛けると、この2種類の微分に分かれます。

&&&prop 平坦なディラック作用素の分解 [prop-decomp]
$$
Γ=\sum_{a<b}e_ae_b\,L_{ab}
$$

と置くと、次が成り立ちます。

$$
\boldsymbol xD=E+Γ
$$
&&&

&&&prf
$\boldsymbol xD=\sum_{a,b}x_ae_ae_b\partial_b$である。$a=b$の項は$e_a^2=1$より$\sum_ax_a\partial_a=E$になる。$a\ne b$の項は$e_be_a=-e_ae_b$を使って$a<b$の対にまとめると、$x_ae_ae_b\partial_b+x_be_be_a\partial_a=e_ae_b(x_a\partial_b-x_b\partial_a)$となる。
&&&

$\boldsymbol x\ne0$では$\boldsymbol x^{-1}=\boldsymbol x/r^2$なので、$D=\boldsymbol x^{-1}(E+Γ)$です。$E$は同次多項式の次数を数え、$Γ$は同次多項式の次数を変えないので、$E$と$Γ$は可換です。

&&&rem ガンマ作用素
$Γ$は、クリフォード解析の文献で**ガンマ作用素**と呼ばれる作用素にあたり、$\boldsymbol x\wedge D$と書かれることもあります。生成元の2乗を$-1$とする規約をとる文献が多く、その場合は符号が逆の作用素を$Γ$と呼ぶことがあります。本記事では上の定義で通します。
&&&

## 球面のラプラシアン

$D^2=\Delta$の分解を調べます。スカラーの関数については、回転の方向の2階微分の和が、極座標のラプラシアンの角度の部分になります。

&&&fml 球面のラプラシアン
$$
Δ_S=\sum_{a<b}L_{ab}^2,\qquad r^2\Delta=E^2+(n-2)E+Δ_S
$$
&&&

&&&prf
$a\ne b$なら$\partial_b(x_b\,\cdot)=1+x_b\partial_b$、$(x_a\partial_b)^2=x_a^2\partial_b^2$などから次を得る。

$$
L_{ab}^2=x_a^2\partial_b^2+x_b^2\partial_a^2-x_a\partial_a-x_b\partial_b-2x_ax_b\partial_a\partial_b
$$

$a\ne b$の対の全体で和を取って半分にすると、$x_a^2\partial_a^2$の形の項は打ち消し合い、次のようになる。

$$
\sum_{a<b}L_{ab}^2=r^2\Delta-(n-1)E-\sum_{a,b}x_ax_b\partial_a\partial_b
$$

$E^2=\sum_{a,b}x_a\partial_a(x_b\partial_b)=E+\sum_{a,b}x_ax_b\partial_a\partial_b$を代入すると公式を得る。
&&&

$r=1$で$E=r\partial_r$の項を除いた$Δ_S$は、単位球面$S^{n-1}$上のラプラシアンです。$L_{ab}$は定数の基底を変えず、各成分の係数関数に作用するので、$Δ_S$はクリフォード代数に値を取る関数の各成分に作用します。

$Γ$は$Δ_S$と次の関係にあります。

&&&prop ガンマ作用素の2乗 [prop-gamma2]
$$
Γ^2=(n-2)Γ-Δ_S
$$
&&&

&&&prf
$\boldsymbol xF$の微分は$\partial_a(\boldsymbol xF)=e_aF+\boldsymbol x\partial_aF$であり、$e_a\boldsymbol x=-\boldsymbol xe_a+2x_a$より次が成り立つ。

$$
D(\boldsymbol xF)=nF-\boldsymbol xDF+2EF
$$

また$E\partial_aF=\partial_a(EF)-\partial_aF$より$ED=D(E-1)$である。この2式と[[prop-decomp]]から次を得る。

$$
(E+Γ)^2F=\boldsymbol xD(\boldsymbol xDF)=n\,\boldsymbol xDF-r^2D^2F+2\boldsymbol xD(E-1)F
=\bigl(n(E+Γ)-r^2\Delta+2(E+Γ)(E-1)\bigr)F
$$

$E$と$Γ$は可換なので、整理すると$r^2\Delta=(E+Γ)(E-Γ+n-2)=E^2+(n-2)E+(n-2)Γ-Γ^2$である。前の公式と比べれば命題を得る。
&&&

&&&rem 平方完成
[[prop-gamma2]]を平方完成すると、次のようになります。

$$
\left(Γ-\frac{n-2}2\right)^2=-Δ_S+\frac{(n-2)^2}4
$$

成分ごとのラプラシアン$-Δ_S$に定数を加えたものの平方根は、この平方完成から$Γ-\frac{n-2}2$として得られます。しかし、これはスピノルの共変微分から作る球面上のディラック作用素とは一致しません。後で見るように$Γ$の固有値は$-k$と$k+n-1$（$k\ge0$）なので、この作用素の固有値は$-\bigl(k+\frac{n-2}2\bigr)$と$k+\frac n2$です。$n\ge3$では正負で揃わず、$n=2$では整数全体になります。ディラック作用素の定数のずれは$\frac{n-1}2$で、その2乗は成分ごとの$-Δ_S$ではなく、スピノルに合わせたラプラシアン（接続ラプラシアン）と曲率の項に分かれます。
&&&

## 対称性

$\operatorname{Cl}_{n,0}(\mathbb R)$に、リバージョン$\tilde A$（積の順序の反転）を使って内積$\langle A,B\rangle=\langle\tilde AB\rangle_0$を入れます。$\langle\ \rangle_0$はスカラー部分です。基底$e_{a_1}\cdots e_{a_k}$（$a_1<\dots<a_k$）は、$e_a^2=1$より$\tilde e_{a_1\cdots a_k}e_{a_1\cdots a_k}=1$を満たすので、この内積について正規直交基底です。$S^{n-1}$上の$\operatorname{Cl}_{n,0}(\mathbb R)$値の関数には、$dΩ$を球面の面積要素として、内積を次のように入れます。

$$
(F,G)=\int_{S^{n-1}}\langle F,G\rangle\,dΩ
$$

&&&prop ガンマ作用素の対称性 [prop-symm]
球面上の滑らかな関数$F,G$について、$(F,ΓG)=(ΓF,G)$です。
&&&

&&&prf
2ベクトル$B=e_ae_b$（$a\ne b$）を左から掛ける演算は、$\tilde B=e_be_a=-B$より$\langle BA,C\rangle=\langle\tilde A\tilde BC\rangle_0=-\langle A,BC\rangle$を満たす。

$L_{ab}f$は、$e_ae_b$平面の中の角度$t$の回転$ρ_t$による$f\circρ_t$の$t=0$での微分である。面積要素は回転で変わらないので$\int f\circρ_t\,dΩ$は$t$によらず、$\int L_{ab}f\,dΩ=0$である。これを$f=\langle F,G\rangle$に当てはめると、ライプニッツ則から$(L_{ab}F,G)=-(F,L_{ab}G)$を得る。

$L_{ab}$と$e_ae_b$を左から掛ける演算は可換なので、$(F,e_ae_bL_{ab}G)=-(e_ae_bF,L_{ab}G)=(L_{ab}e_ae_bF,G)=(e_ae_bL_{ab}F,G)$である。$a<b$で和を取ればよい。
&&&

$Γ$の固有値は実数で、異なる固有値の固有関数は直交します。平坦な$D$の側では、$Γ$を対称にするために虚数単位を掛ける必要はありません。2乗が$-1$になる2ベクトルが、その役を担っています。

# 球面モノジェニックス

## 固有値の2系列

$Γ$は同次多項式の次数を変えないので、固有関数は同次多項式の中から探せます。[[prop-decomp]]によれば、$D=0$と$E=k$から$Γ=-k$が従います。以前の記事で扱わなかった次の対象が、ちょうどその条件を満たします。[[7shi-cla5]]

&&&def 球面モノジェニックス
$\operatorname{Cl}_{n,0}(\mathbb R)$に値を取る次数$k$の同次多項式$P$で、$DP=0$を満たすものの全体を$M_k$と書きます。$M_k$の元を単位球面$S^{n-1}$に制限したものを、次数$k$の**球面モノジェニックス**と呼びます。
&&&

&&&prop 球面モノジェニックスの固有値 [prop-eigen]
$P\in M_k$なら、次が成り立ちます。

$$
ΓP=-kP,\qquad Γ(\boldsymbol xP)=(k+n-1)\,\boldsymbol xP
$$
&&&

&&&prf
[[prop-decomp]]より$0=\boldsymbol xDP=(E+Γ)P=kP+ΓP$である。[[prop-gamma2]]の証明の式から$D(\boldsymbol xP)=nP-\boldsymbol xDP+2EP=(n+2k)P$なので、$(E+Γ)(\boldsymbol xP)=\boldsymbol xD(\boldsymbol xP)=(n+2k)\boldsymbol xP$である。$\boldsymbol xP$は次数$k+1$の同次多項式なので、$E(\boldsymbol xP)=(k+1)\boldsymbol xP$を引けばよい。
&&&

$P\ne0$なら、$\boldsymbol x$は$\boldsymbol x\ne0$で可逆なので$\boldsymbol xP\ne0$です。$\boldsymbol xP$はモノジェニックではなく、$D(\boldsymbol xP)=(n+2k)P$です。固有値は$0,-1,-2,\dots$の系列と、$n-1,n,n+1,\dots$の系列に分かれます。

## 空間の次元

$M_k$の次元は、$x_0=0$の超平面への制限から数えられます。

&&&prop 球面モノジェニックスの次元 [prop-dim]
$M_k$の実次元は、次のとおりです。

$$
\dim M_k=2^n\binom{k+n-2}{n-2}
$$

偶部分代数$\operatorname{Cl}_{n,0}^0(\mathbb R)$に値を取るものに限った$M_k^0$の実次元は、その半分の$2^{n-1}\binom{k+n-2}{n-2}$です。
&&&

&&&prf
$D'=\sum_{l=1}^{n-1}e_l\partial_l$とすると$D=e_0\partial_0+D'$であり、$e_0^2=1$より$DP=0$は$\partial_0P=-e_0D'P$と同値である。$P=\sum_jx_0^jp_j$（$p_j$は$x_1,\dots,x_{n-1}$の多項式）と展開すると、この条件は次の漸化式になる。

$$
p_{j+1}=\frac1{j+1}(-e_0D')p_j
$$

$D'$は$x_0$を含まないので$x_0$の掛け算と可換である。したがって$P$は$p_0=P|_{x_0=0}$から一意に決まる。逆に$x_1,\dots,x_{n-1}$の次数$k$の同次多項式$p$から$p_j=(-e_0D')^jp/j!$と置くと、$D'$は次数を1下げるので$j>k$で$p_j=0$となり、$P=\sum_jx_0^jp_j$は次数$k$の同次多項式で漸化式を満たす。よって$M_k$は、$\operatorname{Cl}_{n,0}(\mathbb R)$に値を取る$n-1$変数の次数$k$の同次多項式の空間と同型である。単項式の個数は$\binom{k+n-2}{n-2}$、値の空間の次元は$2^n$である。$-e_0D'$は偶部分と奇部分をそれぞれ保つので、偶部分に値を取るものに限れば$2^{n-1}$倍になる。
&&&

$n=2$では$\dim M_k=4$で、$k$によりません。$n=3$では$\dim M_k=8(k+1)$です。

&&&rem 固有関数の完全性
$D^2=\Delta$より、$M_k$の元は調和多項式です。$\operatorname{Cl}_{n,0}(\mathbb R)$に値を取る次数$k$の調和多項式の空間は$M_k\oplus\boldsymbol xM_{k-1}$（$M_{-1}=0$とする）と分かれ（フィッシャー分解）、球面調和関数は球面上の2乗可積分な関数の空間で完全系をなすことが知られています。したがって$Γ$の固有値は$-k$と$k+n-1$（$k\ge0$）で尽きます。本記事では、この完全性を主張として使います。
&&&

# 球面上のディラック作用素

## 接ベクトルの掛け算

球面上の作用素を、球面自身の言葉で書き直します。単位球面上の点$\boldsymbol x$で、$\boldsymbol x$と直交する正規直交なベクトル$\boldsymbol t_1,\dots,\boldsymbol t_{n-1}$（接空間の正規直交基底）を取ります。

&&&prop 接方向による表示 [prop-tangent]
単位球面上で、次が成り立ちます。

$$
Γ=-\sum_{i=1}^{n-1}(\boldsymbol t_i\boldsymbol x)\,\partial_{\boldsymbol t_i}
$$

$\partial_{\boldsymbol t}$は$\boldsymbol t$方向の微分です。
&&&

&&&prf
$\boldsymbol t_1,\dots,\boldsymbol t_{n-1},\boldsymbol x$は$\mathbb R^n$の正規直交基底であり、$D=\sum_ae_a\partial_a$は正規直交基底の取り方によらないので、$D=\sum_i\boldsymbol t_i\partial_{\boldsymbol t_i}+\boldsymbol x\partial_{\boldsymbol x}$である。$r=1$で$\partial_{\boldsymbol x}=E$、$\boldsymbol x^2=1$、$\boldsymbol x\boldsymbol t_i=-\boldsymbol t_i\boldsymbol x$より、$\boldsymbol xD=-\sum_i\boldsymbol t_i\boldsymbol x\,\partial_{\boldsymbol t_i}+E$となる。[[prop-decomp]]と比べればよい。
&&&

$Γ$の中で接方向の微分に掛かっているのは、接ベクトル$\boldsymbol t$そのものではなく、法線$\boldsymbol x$との積$\boldsymbol t\boldsymbol x$です。そこで、球面上の接ベクトル$X$を掛ける演算を次のように定めます。

$$
X\bullet ψ=X\boldsymbol x\,ψ
$$

$X$と$\boldsymbol x$は直交するので、$X\boldsymbol x$は2ベクトルです。直交する接ベクトル$X,Y$について$(X\boldsymbol x)^2=-X^2\boldsymbol x^2=-|X|^2$、$(X\boldsymbol x)(Y\boldsymbol x)=-XY=YX=-(Y\boldsymbol x)(X\boldsymbol x)$なので、$\boldsymbol t_i\boldsymbol x$は2乗が$-1$で互いに反交換します。球面の接空間の上のクリフォード代数は、生成元の2乗が$-1$の$\operatorname{Cl}_{0,n-1}(\mathbb R)$として働きます。

点$\boldsymbol x=e_0$では、接ベクトル$e_l$（$l\ge1$）を掛ける演算は$e_le_0=-h_l$です。以前の記事のパラベクトル変数の生成元$h_l=e_0e_l$と符号を除いて一致し、偶部分代数$\operatorname{Cl}_{n,0}^0(\mathbb R)\cong\operatorname{Cl}_{0,n-1}(\mathbb R)$を生成します。[[7shi-cla5]][[7shi-clif1]]

$X\bullet$も$Γ$も偶部分と奇部分をそれぞれ保つので、以下、球面上の作用素は$\operatorname{Cl}_{n,0}^0(\mathbb R)$に値を取る滑らかな関数に作用させます。偶部分代数は、各点で$X\mapsto X\boldsymbol x$により接空間のクリフォード代数と同一視でき、その上に左からの積で作用します。既約なスピノルを最初から選ぶ代わりに、偶部分代数全体を値の空間とするモデルを使い、その値を取る関数を本記事での**スピノル**と呼びます。右から定数を掛ける演算は、以下の作用素とすべて可換です。

## 共変微分

$-Γ$は、接ベクトルの掛け算と接方向の微分を組み合わせた形をしています。しかし、関数を接方向にそのまま微分する$\partial_X$は、接ベクトルの掛け算と両立しません。接ベクトル場$Y$と関数$ψ$について、$\partial_X(Y\bullet ψ)$を計算します。

$$
\partial_X(Y\boldsymbol xψ)=(\partial_XY)\boldsymbol xψ+YXψ+Y\boldsymbol x\,\partial_Xψ
$$

中央の項は、法線$\boldsymbol x$自身が$X$方向に動くことから来ます（$\partial_X\boldsymbol x=X$）。また$\partial_XY$は一般に接ベクトルではありません。$Y\cdot\boldsymbol x=0$を微分すると$(\partial_XY)\cdot\boldsymbol x=-X\cdot Y$なので、接方向への射影を$\nabla_XY=\partial_XY+(X\cdot Y)\boldsymbol x$と書くと、$(\partial_XY)\boldsymbol x=(\nabla_XY)\boldsymbol x-X\cdot Y$です。$\nabla_XY$は球面上のベクトル場の微分（レビ＝チビタ接続）です。まとめると次のようになります。

$$
\partial_X(Y\bulletψ)=(\nabla_XY)\bulletψ+Y\bullet\partial_Xψ+(YX-X\cdot Y)ψ
$$

最後の項$Y\wedge X$がライプニッツ則からのずれです。スピノルの微分には、ベクトル場の微分とライプニッツ則で結びつくこと、すなわち$\nabla_X(Y\bulletψ)=(\nabla_XY)\bulletψ+Y\bullet\nabla_Xψ$を要請します。$\partial_X$に補正$cX\boldsymbol x$を加えると、ずれは$Y\wedge X+c(X\boldsymbol xY\boldsymbol x-Y\boldsymbol xX\boldsymbol x)=Y\wedge X+c(YX-XY)=(1+2c)\,Y\wedge X$となり、$c=-\frac12$で消えます。ただし$n=2$では接空間が1次元で$Y\wedge X=0$となり、この要請では$c$が決まりません。$n=2$でも一般の$n$と同じ$c=-\frac12$を採ります。

&&&def スピノルの共変微分
球面上の点$\boldsymbol x$での接ベクトル$X$について、次のように定めます。

$$
\nabla_Xψ=\partial_Xψ-\frac12X\boldsymbol x\,ψ
$$
&&&

&&&prop 共変微分の性質 [prop-compat]
球面上の接ベクトル場$X,Y$とスピノル$ψ,χ$について、次が成り立ちます。

1. $\nabla_X(Y\bulletψ)=(\nabla_XY)\bulletψ+Y\bullet\nabla_Xψ$
2. $\partial_X\langle ψ,χ\rangle=\langle\nabla_Xψ,χ\rangle+\langle ψ,\nabla_Xχ\rangle$
&&&

&&&prf
1は上の計算で$c=-\frac12$とした場合である。2は、$\partial_X$についてライプニッツ則が成り立つことと、補正$-\frac12X\boldsymbol x$が2ベクトルで、[[prop-symm]]の証明で見たとおり左から掛ける演算が内積について反対称であることから従う。
&&&

&&&rem 標準的なスピン接続との関係
接ベクトルの掛け算と両立し内積を保つという要請だけでは、共変微分は一意に決まりません。接ベクトルの掛け算と可換で、内積について反対称な零階の項を加える余地が残ります。上の$\nabla$は、正規直交枠とレビ＝チビタ接続から組む標準的なスピン接続と一致することが知られています。平坦な空間の中の超曲面で、スピノルの微分を法線で補正する公式（スピノルのガウスの公式）の、球面の場合にあたります。一般の$n$での一致は本記事では主張に留め、$n=2,3$では後で枠に移して確かめます。

本記事の構成は、周囲のユークリッド空間の標準的なスピン構造から球面に誘導されるものを使っています。スピン構造を抽象的に定義してから出発する必要はありませんが、その選択が構成に組み込まれています。$n\ge3$では球面は単連結で、スピン構造は1つしかないことが知られています。円周には2つあり、本記事の構成（$n=2$で$c=-\frac12$とすること）は、円板から誘導されるほうにあたります。
&&&

## ディラック作用素

平坦な$D=\sum_ae_a\partial_a$にならい、接ベクトルの掛け算と共変微分を組み合わせます。

&&&def 球面上のディラック作用素
$$
D_Sψ=\sum_{i=1}^{n-1}\boldsymbol t_i\bullet\nabla_{\boldsymbol t_i}ψ
$$

右辺は正規直交基底$\boldsymbol t_i$の取り方によりません。
&&&

以前の記事ではパラベクトル変数の作用素を$\mathcal D$と書きましたが、本記事では混同を避けて$D_S$と書きます。[[7shi-cla5]]

&&&prop ディラック作用素の分解による表示 [prop-dirac]
$$
D_S=\frac{n-1}2-Γ
$$
&&&

&&&prf
$D_Sψ=\sum_i\boldsymbol t_i\boldsymbol x\,\partial_{\boldsymbol t_i}ψ-\frac12\sum_i(\boldsymbol t_i\boldsymbol x)^2ψ$である。第1項は[[prop-tangent]]より$-Γψ$、第2項は$(\boldsymbol t_i\boldsymbol x)^2=-1$より$\frac{n-1}2ψ$である。
&&&

$Γ$から$D_S$へのずれ$\frac{n-1}2$は、共変微分の補正$-\frac12X\boldsymbol x$を$n-1$個の接方向について足したものです。定数と対称な$Γ$の差なので、$D_S$も対称です。

$D_S$の1階の部分は$\sum_i\boldsymbol t_i\boldsymbol x\,\partial_{\boldsymbol t_i}$です。計量で余接ベクトルを接ベクトルと同一視し、$\partial_{\boldsymbol t_i}$を接ベクトル$\boldsymbol\xi$の成分$\boldsymbol\xi\cdot\boldsymbol t_i$に置き換えると、$\boldsymbol\xi\mapsto\boldsymbol\xi\boldsymbol x=\boldsymbol\xi\bullet$になります。以前の記事で平坦な$D$の主表象として見たベクトルを掛ける写像は、球面の上では接ベクトルを掛ける演算として現れます。[[7shi-kth1]]

&&&thm 球面上のディラック作用素の固有値 [thm-spectrum]
$P\in M_k^0$（偶部分代数に値を取る次数$k$のモノジェニックな同次多項式）について、次が成り立ちます。

$$
D_SP=\left(k+\frac{n-1}2\right)P,\qquad D_S(\boldsymbol xPe_0)=-\left(k+\frac{n-1}2\right)\boldsymbol xPe_0
$$

固有関数の完全性を認めると、$D_S$の固有値は$\pm\bigl(k+\frac{n-1}2\bigr)$（$k=0,1,2,\dots$）で尽き、各固有値の固有空間の実次元は$2^{n-1}\binom{k+n-2}{n-2}$です。
&&&

&&&prf
[[prop-eigen]]と[[prop-dirac]]から、$D_SP=\bigl(\frac{n-1}2+k\bigr)P$、$D_S(\boldsymbol xP)=\bigl(\frac{n-1}2-(k+n-1)\bigr)\boldsymbol xP$である。$\boldsymbol xP$は奇部分に値を取るので、右から$e_0$を掛けて偶部分に値を取る関数に移す。$D_S$は左から作用するので、右からの定数倍と可換である。

$D$は偶部分と奇部分を入れ替えるので、$M_k$の元の偶部分と奇部分は、それぞれモノジェニックである。奇部分に値を取る$Q\in M_k$は$Qe_0\in M_k^0$と対応し、$\boldsymbol xQ=\boldsymbol x(Qe_0)e_0$である。完全性を認めれば、$M_k$と$\boldsymbol xM_k$（$k\ge0$）の元は$\operatorname{Cl}_{n,0}(\mathbb R)$値の関数の完全系をなすので、偶部分を取り出すと、$M_k^0$と$\boldsymbol xM_k^0e_0$の元が偶部分代数に値を取る関数の完全系をなす。正の固有値$k+\frac{n-1}2$を与えるのは$M_k^0$だけで、$P\mapsto\boldsymbol xPe_0$は単射なので、負の固有値の固有空間も同じ次元を持つ。次元は[[prop-dim]]による。
&&&

$Γ$の2系列$-k$と$k+n-1$は、$\frac{n-1}2$だけずらすことで、原点に関して対称な$\pm\bigl(k+\frac{n-1}2\bigr)$に揃いました。$n$が偶数なら固有値は半整数、$n$が奇数なら整数です。固有値の絶対値の最小値は$\frac{n-1}2>0$なので、$0$は固有値ではありません。

偶部分代数は、$\operatorname{Cl}_{0,n-1}(\mathbb R)$の既約加群（最小左イデアル）の直和に分かれます。各成分への射影は右からの積で与えられるので$D_S$と可換であり、値をそれぞれの既約加群$W$に制限できます。$W$に値を取る関数についても同様に数えられ、固有値は同じ$\pm\bigl(k+\frac{n-1}2\bigr)$、各固有値の実重複度は$\dim W\binom{k+n-2}{n-2}$です。

# 低次元の球面

## 円周

$n=2$では$S^1$上の作用素になります。$\boldsymbol x=\cosφ\,e_0+\sinφ\,e_1$、$J=e_0e_1$（$J^2=-1$）とすると、$L_{01}=\partial_φ$なので$Γ=J\partial_φ$です。偶部分代数は$1$と$J$で張られ、複素数と同型です。

&&&ex 円周上のディラック作用素
[[prop-dirac]]より、次のようになります。

$$
D_S=\frac12-J\partial_φ
$$

$ψ=e^{Jmφ}c$（$m$は整数、$c$は偶部分代数の定数）について$D_Sψ=\bigl(m+\frac12\bigr)ψ$です。固有値$m+\frac12$の全体は、$m\ge0$の$k+\frac12$と、$m=-(k+1)$の$-\bigl(k+\frac12\bigr)$に分かれ、[[thm-spectrum]]の$n=2$の場合にあたります。
&&&

円周では、フーリエ級数の完全性から、固有関数はこれで尽きます。核についても直接確かめられます。$D_Sψ=0$に左から$J$を掛けると$\partial_φψ=-\frac12Jψ$となり、解は$ψ=e^{-Jφ/2}c$です。$φ$が$2\pi$進むと$e^{-Jπ}=-1$が掛かるので、$c\ne0$なら$ψ$は円周上の1価の関数になりません。したがって核は$0$です。[[7shi-leb2]]

円周上を動く点の接ベクトル$\boldsymbol t=-\sinφ\,e_0+\cosφ\,e_1$について$\boldsymbol t\boldsymbol x=-J$なので、共変微分は$\nabla_{\boldsymbol t}=\partial_φ+\frac12J$です。回転子$U=e^{-Jφ/2}$は$Ue_0\tilde U=\boldsymbol x$、$Ue_1\tilde U=\boldsymbol t$を満たし、定点での基底を動く点の基底に移します。$ψ=Uψ_f$と置くと、$\tilde U\partial_φU=-\frac12J$が補正$\frac12J$を打ち消して$\tilde U\nabla_{\boldsymbol t}(Uψ_f)=\partial_φψ_f$となり、$D_S$は$-J\partial_φ$に移ります。枠の側の成分$ψ_f$は、$U$の符号の反転を引き受けて、$φ$の1周で符号を変える関数になります。そこでは固有値$μ$の固有関数が$e^{Jμφ}$（$μ$は半整数）であり、半整数の固有値は、この反周期性から来ます。円周のもう一方のスピン構造では、枠の側の成分を周期的に取ることになり、固有値は整数で、定数関数が核に入ります。本記事の構成は、反周期的なほう（円板から誘導されるほう）を選んでいます。

## 2次元の球面

$n=3$では、$x_2$の軸を極に取って$\boldsymbol x=\sinθ\cosφ\,e_0+\sinθ\sinφ\,e_1+\cosθ\,e_2$とします。$\hatθ=\partial_θ\boldsymbol x$と$\hatφ=\partial_φ\boldsymbol x/\sinθ$は正規直交な接ベクトルです。回転子を次のように置きます。

$$
U=e^{-e_0e_1φ/2}\,e^{-e_2e_0θ/2}
$$

$U$は$e_2e_0$平面で$θ$、$e_0e_1$平面で$φ$の回転を表し、$Ue_0\tilde U=\hatθ$、$Ue_1\tilde U=\hatφ$、$Ue_2\tilde U=\boldsymbol x$を満たします。$ψ=Uψ_f$と置き、極を除いた範囲で、枠$(\hatθ,\hatφ)$の側の成分$ψ_f$で$D_S$を書き直します。

&&&prop 2次元の球面の枠による表示 [prop-frame]
$f_1=e_0e_2$、$f_2=e_1e_2$とすると、枠の側では共変微分とディラック作用素が次の形になります。

$$
\tilde U\nabla_{\hatθ}U=\partial_θ,\qquad
\tilde U\nabla_{\partial_φ}U=\partial_φ-\frac12\cosθ\,e_0e_1
$$

$$
\tilde UD_SU=f_1\left(\partial_θ+\frac{\cotθ}2\right)+\frac{f_2}{\sinθ}\partial_φ
$$

$\nabla_{\partial_φ}=\sinθ\,\nabla_{\hatφ}$は座標$φ$の方向の共変微分です。
&&&

&&&prf
$\tilde U\hatθ\boldsymbol xU=e_0e_2=f_1$、$\tilde U\hatφ\boldsymbol xU=e_1e_2=f_2$である。

$θ$方向。$e^{-e_2e_0θ/2}$は$e_2e_0$と可換なので$\tilde U\partial_θU=-\frac12e_2e_0=\frac12f_1$であり、共変微分の補正$-\frac12\tilde U\hatθ\boldsymbol xU=-\frac12f_1$と打ち消し合う。

$φ$方向。$\tilde U\partial_φU=-\frac12e^{e_2e_0θ/2}e_0e_1e^{-e_2e_0θ/2}$である。$e^{e_2e_0θ/2}$による回転は$e_0$を$\cosθ\,e_0+\sinθ\,e_2$に移し、$e_1$を動かさないので、これは$-\frac12(\cosθ\,e_0e_1+\sinθ\,e_2e_1)$である。補正は$-\frac12\sinθ\,f_2=-\frac12\sinθ\,e_1e_2$なので、$e_2e_1$の項が打ち消し合い、$-\frac12\cosθ\,e_0e_1$が残る。

ディラック作用素は$f_1\partial_θ+\frac{f_2}{\sinθ}\bigl(\partial_φ-\frac12\cosθ\,e_0e_1\bigr)$となる。$f_2e_0e_1=e_1e_2e_0e_1=-e_0e_2=-f_1$より、最後の項は$\frac{\cotθ}2f_1$である。
&&&

$f_1,f_2$は2乗が$-1$で反交換するので、$f_1\mapsto-iσ_1$、$f_2\mapsto-iσ_2$（$σ_a$はパウリ行列）で$\operatorname{Cl}_{0,2}(\mathbb R)$を表現できます。このとき$e_0e_1=-f_1f_2\mapsto iσ_3$であり、[[prop-frame]]は次の形になります。

$$
\nabla_φ=\partial_φ-\frac i2\cosθ\,σ_3,\qquad
D_S=-i\left(σ_1\left(\partial_θ+\frac{\cotθ}2\right)+\frac{σ_2}{\sinθ}\partial_φ\right)
$$

これは、正規直交枠$(\partial_θ,\partial_φ/\sinθ)$とスピン接続から組み立てる、$S^2$上のディラック作用素の標準的な形です。$-\frac i2\cosθ\,σ_3$がスピン接続の項で、$\frac{\cotθ}2$はそこから来ます。

同じ共変微分を、2通りの自明化で見たことになります。$\mathbb R^3$の定まった基底で成分を取ると、共変微分の補正は$-\frac12X\boldsymbol x$で、ディラック作用素に現れる零階の項は定数$\frac{n-1}2=1$にまとまります。球面上を動く枠で成分を取ると、補正はスピン接続の項$-\frac12\cosθ\,e_0e_1$になり、ディラック作用素の零階の項は$\frac{\cotθ}2$になります。平坦な空間の基底は各点で共通なので、周囲の空間から誘導される構成では、スピノルを球面上の1価の関数として書けます。

回転子$U$は$φ$の1周で$e^{-e_0e_1π}=-1$を掛けられます。単位四元数が回転の1周で$-1$までしか進まないのと同じ事情です。$ψ$が1価なら、枠の側の成分$ψ_f=\tilde Uψ$は$φ$の1周で符号を変えます。スピノル自体が2価なのではなく、枠を持ち上げた回転子の符号の変化を、成分の反周期性が補っています。正規直交枠で書いたスピノルに現れる符号の2価性は、平坦な空間の基底から動く枠に乗り換えるときに生じます。[[7shi-cover]]

&&&ex 2次元の球面の固有値
[[thm-spectrum]]で$n=3$とすると、固有値は$\pm(k+1)$（$k\ge0$）で、固有空間の実次元は$4(k+1)$です。偶部分代数$\operatorname{Cl}_{3,0}^0(\mathbb R)$は四元数と同型で、$\operatorname{Cl}_{0,2}(\mathbb R)$の既約加群そのものです。右から$e_0e_1$を掛ける演算を虚数単位とみなすと、これは$D_S$と可換で、偶部分代数は$\mathbb C^2$になります。複素次元で数えた重複度は$2(k+1)$です。$j=k+\frac12$と置くと、固有値は$\pm\bigl(j+\frac12\bigr)$、重複度は$2j+1$となり、$j$は半整数です。
&&&

# リヒネロビッチの公式

## 接続ラプラシアン

平坦な場合は$D^2=\Delta$でした。球面上では、$D_S^2$を共変微分の2階の和と比べます。接空間の正規直交基底は、一般には球面全体で連続に取れません（$S^2$では毛玉の定理による）。そこで代わりに、球面全体で定義された回転の場$\boldsymbol v_{ab}=x_ae_b-x_be_a$を使います。[[7shi-kth3]]

&&&lem 回転の場の和 [lem-rot]
単位球面上で、ベクトル$\boldsymbol u,\boldsymbol w$について次が成り立ちます。

$$
\sum_{a<b}(\boldsymbol v_{ab}\cdot\boldsymbol u)(\boldsymbol v_{ab}\cdot\boldsymbol w)=\boldsymbol u\cdot\boldsymbol w-(\boldsymbol x\cdot\boldsymbol u)(\boldsymbol x\cdot\boldsymbol w)
$$

したがって、接ベクトルについて双線形な式$B$に対して$\sum_{a<b}B(\boldsymbol v_{ab},\boldsymbol v_{ab})=\sum_iB(\boldsymbol t_i,\boldsymbol t_i)$です。
&&&

&&&prf
$\boldsymbol v_{ab}\cdot\boldsymbol u=x_au_b-x_bu_a$である。$a,b$の全体で和を取って半分にすると$\frac12\sum_{a,b}(x_au_b-x_bu_a)(x_aw_b-x_bw_a)=|\boldsymbol x|^2\,\boldsymbol u\cdot\boldsymbol w-(\boldsymbol x\cdot\boldsymbol u)(\boldsymbol x\cdot\boldsymbol w)$であり、$|\boldsymbol x|=1$とすればよい。右辺は接空間への射影であり、$\sum_i(\boldsymbol t_i\cdot\boldsymbol u)(\boldsymbol t_i\cdot\boldsymbol w)$に等しい。
&&&

とくに$\sum_{a<b}(\boldsymbol v_{ab}\boldsymbol x)L_{ab}=\sum_i(\boldsymbol t_i\boldsymbol x)\partial_{\boldsymbol t_i}=-Γ$、$\sum_{a<b}|\boldsymbol v_{ab}|^2=n-1$です。

共変微分の2階の和は、通常、正規直交な接ベクトル場$\boldsymbol t_i$を局所的に取って次のように定めます。これを**接続ラプラシアン**と呼びます。

$$
\nabla^*\nabla ψ=-\sum_{i=1}^{n-1}\bigl(\nabla_{\boldsymbol t_i}\nabla_{\boldsymbol t_i}ψ-\nabla_{\nabla_{\boldsymbol t_i}\boldsymbol t_i}ψ\bigr)
$$

第2項は、枠$\boldsymbol t_i$自身が曲がることの補正です。回転の場を使うと、この補正が消えます。

&&&prop 回転の場による接続ラプラシアン [prop-rough-rot]
$$
\nabla^*\nabla ψ=-\sum_{a<b}\nabla_{\boldsymbol v_{ab}}\nabla_{\boldsymbol v_{ab}}ψ
$$
&&&

&&&prf
接ベクトル場$X,Y$について$H(X,Y)=\nabla_X\nabla_Yψ-\nabla_{\nabla_XY}ψ$と置く。$H$は$X$について関数倍で線形であり、$Y$についても$\nabla_X(fY)=(Xf)Y+f\nabla_XY$と$\nabla_X(f\nabla_Yψ)=(Xf)\nabla_Yψ+f\nabla_X\nabla_Yψ$から$(Xf)\nabla_Yψ$の項が打ち消し合って関数倍で線形である。したがって各点での$H$の値は、その点での$X,Y$だけで決まる双線形な式であり、[[lem-rot]]より$\sum_iH(\boldsymbol t_i,\boldsymbol t_i)=\sum_{a<b}H(\boldsymbol v_{ab},\boldsymbol v_{ab})$である。

$\partial_{\boldsymbol v_{ab}}\boldsymbol v_{ab}=L_{ab}\boldsymbol v_{ab}=-(x_ae_a+x_be_b)$で、その接方向への射影は$\nabla_{\boldsymbol v_{ab}}\boldsymbol v_{ab}=-(x_ae_a+x_be_b)+(x_a^2+x_b^2)\boldsymbol x$である。$a<b$で和を取ると$-(n-1)\boldsymbol x+(n-1)\boldsymbol x=0$なので、$\sum_{a<b}\nabla_{\nabla_{\boldsymbol v_{ab}}\boldsymbol v_{ab}}ψ=0$であり、命題を得る。
&&&

&&&prop 接続ラプラシアンの表示 [prop-rough]
$$
\nabla^*\nabla=-Δ_S-Γ+\frac{n-1}4
$$
&&&

&&&prf
[[prop-rough-rot]]の表示を使う。$\boldsymbol v=\boldsymbol v_{ab}$、$L=L_{ab}$と書くと$\nabla_{\boldsymbol v}=L-\frac12\boldsymbol v\boldsymbol x$なので、次のようになる。

$$
\nabla_{\boldsymbol v}\nabla_{\boldsymbol v}ψ=L^2ψ-\frac12\bigl(L(\boldsymbol v\boldsymbol x)\bigr)ψ-\boldsymbol v\boldsymbol x\,Lψ+\frac14(\boldsymbol v\boldsymbol x)^2ψ
$$

$a<b$で和を取る。第1項は$Δ_S$、第3項は[[lem-rot]]の後の式から$Γψ$、第4項は$(\boldsymbol v\boldsymbol x)^2=-|\boldsymbol v|^2$より$-\frac{n-1}4ψ$である。第2項では、$L\boldsymbol x=\boldsymbol v$と$L\boldsymbol v_{ab}=-(x_ae_a+x_be_b)$から$\sum_{a<b}L(\boldsymbol v\boldsymbol x)=-(n-1)\boldsymbol x^2+\sum_{a<b}|\boldsymbol v_{ab}|^2=0$である。符号を反転すれば命題を得る。
&&&

## 公式

&&&thm リヒネロビッチの公式 [thm-lich]
単位球面$S^{n-1}$上で、次が成り立ちます。

$$
D_S^2=\nabla^*\nabla+\frac{(n-1)(n-2)}4
$$
&&&

&&&prf
[[prop-dirac]]と[[prop-gamma2]]から次を得る。

$$
D_S^2=Γ^2-(n-1)Γ+\frac{(n-1)^2}4=-Δ_S-Γ+\frac{(n-1)^2}4
$$

[[prop-rough]]を引くと$\frac{(n-1)^2}4-\frac{n-1}4=\frac{(n-1)(n-2)}4$が残る。
&&&

単位球面$S^m$のスカラー曲率は$R=m(m-1)$であることが知られています。$m=n-1$とすると$R=(n-1)(n-2)$なので、公式は$D_S^2=\nabla^*\nabla+\frac R4$という一般の形になります。平坦な場合は$e_a^2=1$から$D^2=\Delta$でしたが、球面側では接ベクトルの掛け算の2乗が負になるため、$D_S^2$の主要部は非負の接続ラプラシアンになります。その接続ラプラシアンと$D_S^2$との差が、曲率の項$\frac R4$です。

## 核の消失

&&&prop 核の消失 [prop-kernel]
$n\ge3$とします。$S^{n-1}$上の滑らかなスピノル$ψ$が$D_Sψ=0$を満たすなら、$ψ=0$です。
&&&

&&&prf
$\nabla_{\boldsymbol v_{ab}}=L_{ab}-\frac12\boldsymbol v_{ab}\boldsymbol x$は反対称である。$L_{ab}$は[[prop-symm]]の証明で見たとおり反対称であり、$\boldsymbol v_{ab}\boldsymbol x$は2ベクトルなので左から掛ける演算は反対称である。したがって$(ψ,\nabla^*\nabla ψ)=\sum_{a<b}\|\nabla_{\boldsymbol v_{ab}}ψ\|^2\ge0$である（$\|ψ\|^2=(ψ,ψ)$）。$D_S$は対称なので、[[thm-lich]]から次を得る。

$$
\|D_Sψ\|^2=(ψ,D_S^2ψ)=\sum_{a<b}\|\nabla_{\boldsymbol v_{ab}}ψ\|^2+\frac{(n-1)(n-2)}4\|ψ\|^2
$$

$D_Sψ=0$なら右辺の各項は$0$であり、$n\ge3$では$(n-1)(n-2)>0$なので$ψ=0$である。
&&&

証明の式から、固有値$λ$について$λ^2\ge\frac{(n-1)(n-2)}4$も従います。この評価は固有関数の完全性を使わずに得られます。

$n=2$ではスカラー曲率が$0$なので、この議論では$\nabla ψ=0$までしか言えません。円周では、共変微分が$0$になる$ψ=e^{-Jφ/2}c$が1価でないことから、核が$0$になりました。

[[thm-spectrum]]の固有値の最小値は$λ^2=\frac{(n-1)^2}4$で、リヒネロビッチの評価$\frac{(n-1)(n-2)}4$より真に大きい値です。評価は核が消えることを示すには十分ですが、最小の固有値を与えるものではありません。

&&&rem フリードリヒの評価
$m\ge2$次元の閉リーマンスピン多様体上のディラック作用素について、スカラー曲率が正の定数$R_0$以上なら、固有値は$λ^2\ge\frac m{4(m-1)}R_0$を満たすことが知られています（フリードリヒの評価）。$m=n-1$、$R_0=(n-1)(n-2)$とすると右辺は$\frac{(n-1)^2}4$となり、球面では絶対値が最小の固有値で等号が成り立ちます。証明は本記事では扱いません。
&&&

# まとめ

平坦なディラック作用素$D$を、$\boldsymbol xD=E+Γ$と動径方向と球面方向に分け、球面方向の作用素$Γ$から球面上のディラック作用素を取り出しました。接ベクトル$X$を掛ける演算を$X\boldsymbol x$とすると、スピノルは偶部分代数に値を取る球面上の関数になり、それと両立する共変微分の補正$-\frac12X\boldsymbol x$が、$Γ$の固有値を$\frac{n-1}2$だけずらします。球面モノジェニックスとそれに$\boldsymbol x$を掛けた関数から、固有値$\pm\bigl(k+\frac{n-1}2\bigr)$を求めました。円周と2次元の球面では、動く枠に移すとディラック作用素の零階の項が定数からスピン接続による項に変わり、枠の側の成分が1周で符号を変えることを確かめました。リヒネロビッチの公式を示し、$n\ge3$で核が$0$になることを導きました。

&&& 分解
$$
\boldsymbol xD=E+Γ,\qquad Γ=\sum_{a<b}e_ae_b(x_a\partial_b-x_b\partial_a),\qquad Γ^2=(n-2)Γ-Δ_S
$$
&&&

&&& 球面モノジェニックス
$P$が次数$k$のモノジェニックな同次多項式なら、次が成り立ちます。
$$
ΓP=-kP,\qquad Γ(\boldsymbol xP)=(k+n-1)\,\boldsymbol xP,\qquad\dim M_k=2^n\binom{k+n-2}{n-2}
$$
&&&

&&& 球面上のディラック作用素
接ベクトルの掛け算を$X\bulletψ=X\boldsymbol xψ$、共変微分を$\nabla_Xψ=\partial_Xψ-\frac12X\boldsymbol xψ$とすると、次のようになります。
$$
D_S=\sum_i\boldsymbol t_i\bullet\nabla_{\boldsymbol t_i}=\frac{n-1}2-Γ
$$
固有値は$\pm\bigl(k+\frac{n-1}2\bigr)$（$k\ge0$）です。
&&&

&&& 2次元の球面
枠$(\hatθ,\hatφ)$に移すと、次の形になります。
$$
D_S=f_1\left(\partial_θ+\frac{\cotθ}2\right)+\frac{f_2}{\sinθ}\partial_φ
$$
固有値は$\pm(k+1)$、複素の重複度は$2(k+1)$です。
&&&

&&& リヒネロビッチの公式
$$
D_S^2=\nabla^*\nabla+\frac R4,\qquad R=(n-1)(n-2)
$$
$n\ge3$では$R>0$から$D_S$の核は$0$です。
&&&

| 項目 | 平坦な空間から取り出した形 | 正規直交枠の形（$n=3$） |
|:---|:---|:---|
| スピノル | 偶部分代数に値を取る1価の関数 | $φ$の1周で符号を変える成分 |
| ディラック作用素の零階の項 | 定数$\frac{n-1}2$ | $\frac{\cotθ}2$（スピン接続$-\frac12\cosθ\,e_0e_1$から） |
| 固有値 | $\pm\bigl(k+\frac{n-1}2\bigr)$ | $\pm(k+1)$ |
