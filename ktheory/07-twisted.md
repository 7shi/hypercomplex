2次元の球面上のディラック作用素をホップ束の冪で捩り、核が同次多項式で書けることと、指数が貼り合わせ関数の回転数から決まることを示します。

# 概要

前回の記事では、平坦なディラック作用素を球面の方向に分解して、球面上のディラック作用素$D_S$を取り出しました。2次元の球面$S^2$では、固有値は$\pm(k+1)$（$k\ge0$）で、リヒネロビッチの公式から核は$0$でした。[[7shi-kth6]]

以前の記事では、$S^2$上の複素直線束が貼り合わせ関数$z^n$で作る$H^n$（$n\in\mathbb Z$）で尽くされ、ホップ束$H$の回転数が$+1$であることを示しました。また、テプリッツ作用素$T_f$の指数が記号の回転数の符号を反転したもの$-\operatorname{wind}(f)$に等しく、$\tilde K(S^2)$の生成元$[H]-1$に$-1$が対応することを見ました。[[7shi-kth3]][[7shi-kth5]]

$D_S$の核が$0$なので、$D_S$から指数を作っても$0$にしかなりません。本記事で扱う問いは次のとおりです。

> 球面上のディラック作用素に核が現れるのはどのような場合か。そのとき、核と余核の次元の差は、貼り合わせ関数の回転数とどう関係するか。

前回の記事のスピノルは、球面全体で定義された偶部分代数値の関数でした。本記事では、スピノルを北側と南側の2つの領域で別々に与え、貼り合わせ関数$z^n$で貼り合わせます。これを$H^n$で**捩る**と言います。貼り合わせと両立するように共変微分を補正すると、捩ったディラック作用素が得られます。補正の項は、前回の記事の注意「標準的なスピン接続との関係」で述べた、共変微分に加える余地のある零階の項にあたります。

前提は次のとおりです。

- **前回の記事から**：$\operatorname{Cl}_{3,0}(\mathbb R)$（$e_a^2=1$）、単位球面上の$\boldsymbol x$、接ベクトルの掛け算$X\bulletψ=X\boldsymbol xψ$、共変微分$\nabla_Xψ=\partial_Xψ-\frac12X\boldsymbol xψ$とその性質、$D_S=\sum_i\boldsymbol t_i\bullet\nabla_{\boldsymbol t_i}$、回転の場$\boldsymbol v_{ab}=x_ae_b-x_be_a$、回転子$U$による枠の形、接続ラプラシアン$\nabla^*\nabla$、リヒネロビッチの公式$D_S^2=\nabla^*\nabla+\frac12$[[7shi-kth6]]
- **貼り合わせ**：南半球$D_-$の座標$v$に貼り合わせ関数$g$を掛けたものが北半球$D_+$の座標になる向き、回転数、$H^n=E_{z^n}$、ホップ束のファイバー$E_x$とホップ写像$h(α,β)=(2αβ^*,\ |α|^2-|β|^2)$[[7shi-kth3]][[7shi-homog]]
- **テプリッツ作用素**：$\operatorname{ind}T_f=-\operatorname{wind}(f)$[[7shi-kth5]]
- **ストークスの定理**：曲面上の回転の積分が境界に沿った周回積分に等しいこと[[7shi-cla2]]

前回の記事では$\operatorname{Cl}_{n,0}(\mathbb R)$の$n$を次元に使いましたが、本記事では周囲のユークリッド空間の次元を$3$に固定し、記号$n$は捩りの冪$H^n$に使います。

本記事は次の順に進みます。

1. $D_S$と反交換する作用素（カイラリティ）を作り、スピノルを2つに分けて指数を定めます。正のカイラリティの部分がホップ束$H$そのものであることを示します。
2. スピノルを$z^n$で貼り合わせ、貼り合わせと両立する共変微分と、その曲率を調べます。曲率の積分が回転数で決まることを示し、曲率が一定の接続を求めます。
3. 捩ったディラック作用素のリヒネロビッチの公式を導き、片方のカイラリティで核が消えることを示します。
4. 曲率が一定の接続について、変数分離で核を求め、核の元が次数$|n|-1$の同次多項式で書けることを示します。指数は$-n$で、テプリッツ作用素の指数と一致します。

証明するのは、カイラリティの性質、正負のカイラリティとホップ束の対応、曲率の積分、捩ったディラック作用素のリヒネロビッチの公式、核の次元と具体形、指数の値です。フーリエ級数の完全性は主張として使います。前回の記事と同じく、ディラック作用素を以前の記事の有界なフレドホルム作用素の枠組みに当てはめる準備には立ち入らず、核の次元の差を指数と呼びます。指数が接続の選び方によらないことも、本記事では扱いません。

# カイラリティ

## 接平面の面積要素

有限次元の線形写像では、核と余核の次元の差は定義域と値域の次元の差で決まりました。ディラック作用素に意味のある指数を持たせるには、作用素を、ある空間から別の空間への写像に分ける必要があります。そのために、$D_S$と反交換し、2乗が$1$になる作用素を探します。$D_S$は接ベクトルの掛け算と共変微分からなるので、その両方と相性のよい作用素が候補です。[[7shi-kth1]]

単位球面上の点$\boldsymbol x$で、正規直交な接ベクトル$\boldsymbol t_1,\boldsymbol t_2$を$\boldsymbol t_1\boldsymbol t_2\boldsymbol x=ω$（$ω=e_0e_1e_2$）となる向きに取ります。このとき$\boldsymbol t_1\boldsymbol t_2=ω\boldsymbol x$であり、$ω\boldsymbol x$は接平面の単位2ベクトル（面積要素）です。$\boldsymbol x$は$\boldsymbol t_i$と反交換するので、2つの接ベクトルの掛け算を続けると、面積要素を左から掛ける演算になります。

$$
\boldsymbol t_1\bullet\boldsymbol t_2\bulletψ=\boldsymbol t_1\boldsymbol x\boldsymbol t_2\boldsymbol xψ=-\boldsymbol t_1\boldsymbol t_2ψ=-ω\boldsymbol xψ
$$

$ω$は$\operatorname{Cl}_{3,0}(\mathbb R)$のすべての元と可換なので、$(ω\boldsymbol x)^2=ω^2\boldsymbol x^2=-1$です。2乗が$-1$なので、そのままでは固有値$\pm1$で空間を分けられません。そこで虚数単位を掛けます。前回の記事の例「2次元の球面の固有値」で見たとおり、右から$J=e_0e_1$を掛ける演算は$D_S$と可換で、偶部分代数$\operatorname{Cl}_{3,0}^0(\mathbb R)$を複素2次元の空間$\mathbb C^2$にします。以下、この複素構造を使い、複素数$a+bi$を右から掛けることを$ψ\mapstoψ(a+bJ)$とします。$1$と$J$が張る部分を$\mathbb C_J$と書きます。[[7shi-kth6]]

&&&def カイラリティ
スピノル$ψ$に対して、次のように定めます。

$$
γψ=ω\boldsymbol x\,ψJ
$$

$γψ=ψ$となるスピノルを**正のカイラリティ**、$γψ=-ψ$となるスピノルを**負のカイラリティ**と呼び、各点でそれぞれの値の空間を$S^+_{\boldsymbol x}$、$S^-_{\boldsymbol x}$と書きます。
&&&

$γ$は前回の記事のガンマ作用素$Γ$とは別のものです。$S^\pm$は、以前の記事の片側のずらし$S$とも関係ありません。[[7shi-kth5]]

&&&prop カイラリティの性質 [prop-chirality]
球面上の接ベクトル$X$とスピノル$ψ,χ$について、次が成り立ちます。

1. $γ^2ψ=ψ$
2. $γ(X\bulletψ)=-X\bulletγψ$
3. $\nabla_X(γψ)=γ\nabla_Xψ$
4. $γ(ψJ)=(γψ)J$、$\langle γψ,χ\rangle=\langle ψ,γχ\rangle$

したがって$D_Sγ=-γD_S$です。
&&&

&&&prf
1. $γ^2ψ=(ω\boldsymbol x)^2ψJ^2=(-1)ψ(-1)=ψ$である。

2. $\boldsymbol xX=-X\boldsymbol x$と$\boldsymbol x^2=1$から、$X\boldsymbol x\,ω\boldsymbol x=ωX\boldsymbol x\boldsymbol x=ωX$、$ω\boldsymbol x\,X\boldsymbol x=-ωX\boldsymbol x\boldsymbol x=-ωX$である。よって$X\boldsymbol x$と$ω\boldsymbol x$は反交換する。

3. $\partial_X\boldsymbol x=X$より$\partial_X(ω\boldsymbol xψJ)=ωXψJ+ω\boldsymbol x(\partial_Xψ)J$である。2の計算から$-\frac12X\boldsymbol x\,ω\boldsymbol xψJ=-\frac12ωXψJ$なので、$\nabla_X(γψ)=\frac12ωXψJ+ω\boldsymbol x(\partial_Xψ)J$となる。一方$γ\nabla_Xψ=ω\boldsymbol x(\partial_Xψ)J-\frac12ω\boldsymbol xX\boldsymbol xψJ=ω\boldsymbol x(\partial_Xψ)J+\frac12ωXψJ$であり、両者は等しい。

4. 前半は左からの積と右からの積の可換性による。後半について、スカラー部分は巡回的に入れ替えられるので$\langle ψJ,χ\rangle=\langle\tilde J\tildeψχ\rangle_0=\langle\tildeψχ\tilde J\rangle_0=-\langle ψ,χJ\rangle$であり、右から$J$を掛ける演算は反対称である。2ベクトル$ω\boldsymbol x$を左から掛ける演算も反対称であり、両者は可換なので、積$γ$は対称である。

$D_S=\sum_i\boldsymbol t_i\bullet\nabla_{\boldsymbol t_i}$なので、2と3から$D_Sγ=-γD_S$を得る。
&&&

$γ$は対称で$γ^2=1$なので、$S^+_{\boldsymbol x}$と$S^-_{\boldsymbol x}$は直交し、$\operatorname{Cl}_{3,0}^0(\mathbb R)=S^+_{\boldsymbol x}\oplus S^-_{\boldsymbol x}$と分かれます。$γ$は複素構造と可換なので、どちらも複素部分空間です。

前回の記事の回転子$U=e^{-Jφ/2}e^{-e_2e_0θ/2}$で枠の側の成分$ψ=Uψ_f$に移すと、$\tilde U\boldsymbol xU=e_2$より$\tilde Uω\boldsymbol xU=ωe_2=J$なので、$γ$は$ψ_f\mapsto Jψ_fJ$になります。$f_1=e_0e_2$、$f_2=e_1e_2$とすると、$J1J=-1$、$JJJ=-J$、$Jf_1J=f_1$、$Jf_2J=f_2$です。$f_2=f_1J$なので、枠の側では$S^+$が$f_1\mathbb C_J$、$S^-$が$\mathbb C_J$で、どちらも複素1次元です。[[7shi-kth6]]

## 指数

[[prop-chirality]]の$D_Sγ=-γD_S$により、$D_S$は正のカイラリティのスピノルを負のカイラリティのスピノルに、負を正に写します。

&&&def ディラック作用素の指数
ディラック作用素$D$が$γ$と反交換するとき、正のカイラリティのスピノルに制限したものを$D^+$、負のカイラリティのスピノルに制限したものを$D^-$と書きます。核が有限次元なら、次の整数を$D$の**指数**と呼びます。

$$
\operatorname{ind}D=\dim_{\mathbb C}\ker D^+-\dim_{\mathbb C}\ker D^-
$$
&&&

$D$が対称なら、正のカイラリティの$ψ$と負のカイラリティの$χ$について$(D^+ψ,χ)=(ψ,D^-χ)$です。したがって、滑らかなスピノルの範囲では、$\ker D^-$は$D^+$の像に直交するものの全体で、$D^+$の余核の役割を果たします。以前の記事の指数$\dim\ker T-\dim\operatorname{coker}T$と同じ形です。ただし、余核を$\ker D^-$と同一視することの正当化には、前回の記事で断ったとおり定義域などの準備が要るので、本記事では上の差を指数の定義として使います。[[7shi-kth5]][[7shi-kth6]]

前回の記事で見たとおり$D_S$の核は$0$なので、$\operatorname{ind}D_S=0$です。

## ホップ束との関係

正のカイラリティの空間$S^+_{\boldsymbol x}$は、各点で$\mathbb C^2$の中の複素直線です。球面の各点に$\mathbb C^2$の直線を割り当てたものは、以前の記事のホップ束と同じ形をしています。実際、両者は一致します。[[7shi-kth3]]

$\mathbb C_J$の元$α,β$に対して、次のように置きます。

$$
Φ(α,β)=f_1α-β
$$

$1,J,f_1,f_2$は偶部分代数の基底なので、$Φ$は$\mathbb C_J^2$から偶部分代数への同型で、右からの積$Φ(α,β)c=Φ(αc,βc)$（$c\in\mathbb C_J$）について複素線形です。以前の記事のホップ写像は、$\mathbb C_J$を$\mathbb C$とみなし、球面の点$\boldsymbol x$を$(x_0+x_1J,\ x_2)$と表すと、$h(α,β)=(2αβ^*,\ |α|^2-|β|^2)$です。$|α|^2+|β|^2=1$で$h(α,β)=\boldsymbol x$となる$(α,β)$の全体は、$\boldsymbol x$の上のホップ束のファイバー$E_{\boldsymbol x}$の単位ベクトルです。[[7shi-kth3]][[7shi-homog]]

&&&prop 正のカイラリティとホップ束 [prop-hopf-chirality]
$E_{\boldsymbol x}$の元$(α,β)$について$γΦ(α,β)=Φ(α,β)$です。したがって$S^+_{\boldsymbol x}=Φ(E_{\boldsymbol x})$、$S^-_{\boldsymbol x}=Φ(E_{-\boldsymbol x})$です。
&&&

&&&prf
$Φ$は複素線形なので、$|α|^2+|β|^2=1$の場合に示せばよい。$c=2αβ^*$、$r=|α|^2-|β|^2$とすると$x_0e_0+x_1e_1=e_0(x_0+x_1J)$より$\boldsymbol x=e_0c+re_2$であり、$ωe_0=f_2=f_1J$、$ωe_2=J$から$ω\boldsymbol x=f_1Jc+rJ$である。$\mathbb C_J$の元$d$について$f_1d=d^*f_1$、また$f_1Jf_1=J$、$Jf_1=-f_1J$である。これらを使うと、次のようになる。

$$
ω\boldsymbol x(f_1α-β)J=-c^*α+f_1cβ+rf_1α+rβ=f_1(cβ+rα)+(rβ-c^*α)
$$

$cβ+rα=α(2|β|^2+|α|^2-|β|^2)=α$、$rβ-c^*α=β(|α|^2-|β|^2-2|α|^2)=-β$なので、右辺は$f_1α-β$である。

$E_{\boldsymbol x}$と$E_{-\boldsymbol x}$は$\mathbb C^2$の直交する2本の直線で、$γ$の定義の$\boldsymbol x$を$-\boldsymbol x$に替えると符号が変わるから、$Φ(E_{-\boldsymbol x})$は$γ$の固有値$-1$の空間に入る。$Φ$は同型なので、$Φ(E_{\boldsymbol x})$と$Φ(E_{-\boldsymbol x})$はどちらも複素1次元で、それぞれ$S^+_{\boldsymbol x}$と$S^-_{\boldsymbol x}$に一致する。
&&&

正のカイラリティのスピノルは、ホップ束の切断を定まった同型$Φ$で写したものです。負のカイラリティの側を貼り合わせ関数で書くため、北極と南極をそれぞれ除いた領域で、$E_{\boldsymbol x}$の単位ベクトルを1つずつ選びます。$θ,φ$を前回の記事の極座標（$x_0+x_1J=\sinθ\,e^{Jφ}$、$x_2=\cosθ$）とします。

$$
(α_+,β_+)=\left(\cos\fracθ2,\ \sin\fracθ2\,e^{-Jφ}\right),\qquad
(α_-,β_-)=\left(\cos\fracθ2\,e^{Jφ},\ \sin\fracθ2\right)
$$

どちらも$h(α_\pm,β_\pm)=\boldsymbol x$を満たし、$(α_+,β_+)=(α_-,β_-)e^{-Jφ}$です。$\boldsymbol x$で書くと次のようになり、$(α_+,β_+)$は南極を除く領域$U_+$（$x_2>-1$）で、$(α_-,β_-)$は北極を除く領域$U_-$（$x_2<1$）で滑らかです。

$$
(α_+,β_+)=\left(\sqrt{\frac{1+x_2}2},\ \frac{x_0-x_1J}{\sqrt{2(1+x_2)}}\right),\qquad
(α_-,β_-)=\left(\frac{x_0+x_1J}{\sqrt{2(1-x_2)}},\ \sqrt{\frac{1-x_2}2}\right)
$$

$β_-$は正の実数で、$(α_-,β_-)$は以前の記事の南半球の座標$λ=β$が正になる単位ベクトルです。同様に$(α_+,β_+)$は北半球の座標$μ=α$が正になる単位ベクトルです。[[7shi-kth3]]

&&&cor カイラリティの空間の貼り合わせ関数 [cor-spinor-hopf]
複素直線束として、正のカイラリティのスピノルの束は$H$、負のカイラリティのスピノルの束は$H^{-1}$と同型です。
&&&

&&&prf
正のカイラリティのスピノル$χ$を$U_\pm$で$χ=Φ(α_\pm,β_\pm)c_\pm$（$c_\pm\in\mathbb C_J$）と表す。$Φ(α_+,β_+)=Φ(α_-,β_-)e^{-Jφ}$より$c_+=e^{Jφ}c_-$であり、赤道$e^{Jφ}=z$の上で貼り合わせ関数は$z$である。

負のカイラリティでは、$(-β^*,α^*)$が$E_{-\boldsymbol x}$の元であることを使う（$h(-β^*,α^*)=(-2β^*α,\ |β|^2-|α|^2)=-h(α,β)$）。$χ=Φ(-β_\pm^*,α_\pm^*)c_\pm$と表すと、$(-β_+^*,α_+^*)=(-β_-^*,α_-^*)e^{Jφ}$より$c_+=e^{-Jφ}c_-$であり、貼り合わせ関数は$z^{-1}$である。
&&&

前回の記事のスピノルは、偶部分代数に値を取る関数として、全体では自明な$\mathbb C^2$の束の切断でした。それが$H\oplus H^{-1}$に分かれることになります。以前の記事で$H\oplus H^{-1}$が自明であることを回転で示したことと整合します。[[7shi-kth4]]

# ホップ束で捩る

## 貼り合わせ

前回の記事のスピノルは、球面全体で定義された関数です。共変微分に零階の項を加えれば作用素は変わりますが、関数の空間は自明な束の切断のままです。貼り合わせ関数のねじれを取り込むには、関数の空間そのものを変える必要があります。そこで、以前の記事の貼り合わせをスピノルに施します。[[7shi-kth3]]

&&&def 捩ったスピノル
$n\in\mathbb Z$とします。$U_+$と$U_-$の上で定義された、偶部分代数に値を取る滑らかな関数の組$ψ=(ψ_+,ψ_-)$で、両方の領域が重なる部分（両極を除いた部分）で次を満たすものを、$H^n$で**捩ったスピノル**と呼びます。

$$
ψ_+=ψ_-e^{nJφ}
$$
&&&

赤道の上では右から$e^{nJφ}=z^n$を掛けることなので、閉じた半球$D_\pm$に制限すると、以前の記事の向き（南半球の座標$v$に$g$を掛けたものが北半球の座標）での貼り合わせ関数$z^n$そのものです。$n=0$なら$ψ_+=ψ_-$で、前回の記事のスピノルに戻ります。

右から$e^{nJφ}$を掛ける演算は、接ベクトルを左から掛ける演算と可換です。また、$e^{nJφ}$は$J$と可換なので、$γ$とも可換です。さらに内積を保ちます（$\langle ψg,χg\rangle=\langle\tilde g\tildeψχg\rangle_0=\langle\tildeψχg\tilde g\rangle_0$）。したがって、接ベクトルの掛け算、カイラリティによる分解、各点の内積$\langle ψ,χ\rangle$と長さ$|ψ|$は、捩ったスピノルについても球面全体で意味を持ちます。

## 接続

微分は貼り合わせと両立しません。$ψ_+=ψ_-e^{nJφ}$を微分すると、$φ$の微分から余分な項が出ます。

$$
\partial_Xψ_+=(\partial_Xψ_-)e^{nJφ}+n(\partial_Xφ)\,ψ_+J
$$

$\partial_Xφ=\nablaφ\cdot X$と書くと、$φ$の勾配は$\nablaφ=\boldsymbol v_{01}/(x_0^2+x_1^2)$です（$\boldsymbol v_{01}=x_0e_1-x_1e_0$）。余分な項は$ψJ$の実数倍なので、各領域の共変微分に同じ形の項を加えて打ち消します。

右から$J$を掛ける演算は、接ベクトルの掛け算と可換で、内積について反対称です。前回の記事の注意「標準的なスピン接続との関係」で、共変微分に加える余地があると述べた零階の項は、ちょうどこの性質を持つものでした。[[7shi-kth6]]

&&&def 捩った共変微分とディラック作用素
$U_\pm$上の接ベクトル場$\boldsymbol a_\pm$で、両極を除いた部分で次を満たすものを取ります。

$$
\boldsymbol a_+-\boldsymbol a_-=-n\nablaφ
$$

$U_\pm$の上で、捩った共変微分とディラック作用素を次のように定めます。

$$
\nabla^{\boldsymbol a}_Xψ_\pm=\nabla_Xψ_\pm+(\boldsymbol a_\pm\cdot X)\,ψ_\pmJ,\qquad
D_{\boldsymbol a}ψ_\pm=\sum_{i=1}^2\boldsymbol t_i\bullet\nabla^{\boldsymbol a}_{\boldsymbol t_i}ψ_\pm
$$
&&&

$\nabla^{\boldsymbol a}_Xψ_+=(\nabla^{\boldsymbol a}_Xψ_-)e^{nJφ}$となるので、$\nabla^{\boldsymbol a}_Xψ$と$D_{\boldsymbol a}ψ$は捩ったスピノルです。$\boldsymbol a_\pm=\sum_i(\boldsymbol a_\pm\cdot\boldsymbol t_i)\boldsymbol t_i$より、ディラック作用素は次のようになります。

$$
D_{\boldsymbol a}ψ=D_Sψ+\boldsymbol a\boldsymbol x\,ψJ
$$

共変微分に加えた項$(\boldsymbol a\cdot X)R_J$（$R_J$は右から$J$を掛ける演算）は、接ベクトルの掛け算と可換で反対称なので、前回の記事の命題「共変微分の性質」（接ベクトルの掛け算とのライプニッツ則、内積の保存）は$\nabla^{\boldsymbol a}$についても成り立ちます。[[7shi-kth6]]

&&&prop 捩ったディラック作用素の性質 [prop-twisted]
1. $D_{\boldsymbol a}γ=-γD_{\boldsymbol a}$
2. 捩ったスピノル$ψ,χ$について$(ψ,D_{\boldsymbol a}χ)=(D_{\boldsymbol a}ψ,χ)$
&&&

&&&prf
1. [[prop-chirality]]より$D_S$は$γ$と反交換する。$\boldsymbol a\boldsymbol x$は接ベクトルの掛け算なので[[prop-chirality]]の2から$ω\boldsymbol x$と反交換し、右からの$J$どうしは可換なので、$ψ\mapsto\boldsymbol a\boldsymbol xψJ$も$γ$と反交換する。

2. 前回の記事の補題「回転の場の和」により、$D_{\boldsymbol a}=\sum_{a<b}\boldsymbol v_{ab}\bullet\nabla^{\boldsymbol a}_{\boldsymbol v_{ab}}$である。$\langle ψ,χ\rangle$は球面全体の関数なので、前回の記事の命題「ガンマ作用素の対称性」の証明と同じく$\int\partial_{\boldsymbol v_{ab}}\langle ψ,χ\rangle\,dΩ=0$であり、内積の保存から$(ψ,\nabla^{\boldsymbol a}_{\boldsymbol v}χ)=-(\nabla^{\boldsymbol a}_{\boldsymbol v}ψ,χ)$である（$\boldsymbol v=\boldsymbol v_{ab}$）。$\boldsymbol v\boldsymbol x$を左から掛ける演算は反対称なので、ライプニッツ則から次を得る。

$$
(ψ,\boldsymbol v\bullet\nabla^{\boldsymbol a}_{\boldsymbol v}χ)=-(\boldsymbol v\bulletψ,\nabla^{\boldsymbol a}_{\boldsymbol v}χ)=\bigl((\nabla_{\boldsymbol v}\boldsymbol v)\bulletψ+\boldsymbol v\bullet\nabla^{\boldsymbol a}_{\boldsymbol v}ψ,\ χ\bigr)
$$

$a<b$で和を取ると、前回の記事で見た$\sum_{a<b}\nabla_{\boldsymbol v_{ab}}\boldsymbol v_{ab}=0$により第1項が消える。
&&&

$\boldsymbol a$の符号を反転すると捩りの向きが逆になります。これを右からの積で表しておきます。$K=f_2$と置きます。$K$は長さ$1$の2ベクトルで、$K^2=-1$、$KJ=-JK$を満たします。

&&&prop 捩りの向きの反転 [prop-flip]
$\boldsymbol a_\pm$が$H^n$の捩りの条件を満たすとき、$ψ\mapstoψK$は$H^n$で捩ったスピノルを$H^{-n}$で捩ったスピノルに写し、正負のカイラリティを入れ替えます。$H^{-n}$の側の接続を$-\boldsymbol a_\pm$とすると、$D_{-\boldsymbol a}(ψK)=(D_{\boldsymbol a}ψ)K$が成り立ちます。
&&&

&&&prf
$KJ=-JK$より$e^{nJφ}K=Ke^{-nJφ}$なので、$ψ_+K=(ψ_-K)e^{-nJφ}$である。$(-\boldsymbol a_+)-(-\boldsymbol a_-)=n\nablaφ$なので、$-\boldsymbol a_\pm$は$H^{-n}$の捩りの条件を満たす。$γ(ψK)=ω\boldsymbol xψKJ=-(γψ)K$である。$D_S$は左から作用するので$D_S(ψK)=(D_Sψ)K$、また$-\boldsymbol a\boldsymbol xψKJ=(\boldsymbol a\boldsymbol xψJ)K$であり、2つを足せば$D_{-\boldsymbol a}(ψK)=(D_{\boldsymbol a}ψ)K$を得る。
&&&

$ψ\mapstoψK$は$(ψc)K=(ψK)c^*$（$c\in\mathbb C_J$）を満たすので、複素線形ではなく複素反線形です。逆写像は右から$-K$を掛ける写像なので、複素部分空間を同じ次元の複素部分空間に写します。

## 曲率

接続$\boldsymbol a_\pm$は、貼り合わせの条件だけでは決まりません。どの選び方にもよらない量を探します。

&&&def 曲率
$\boldsymbol t_1\boldsymbol t_2\boldsymbol x=ω$となる正規直交な接ベクトル$\boldsymbol t_1,\boldsymbol t_2$について、次の関数$F$を接続の**曲率**と呼びます。

$$
F=\boldsymbol t_2\cdot\partial_{\boldsymbol t_1}\boldsymbol a-\boldsymbol t_1\cdot\partial_{\boldsymbol t_2}\boldsymbol a
$$
&&&

$F$は$\sum_i\boldsymbol t_i\wedge\partial_{\boldsymbol t_i}\boldsymbol a$の面積要素$\boldsymbol t_1\wedge\boldsymbol t_2=ω\boldsymbol x$の係数なので、接平面の中での$\boldsymbol t_1,\boldsymbol t_2$の回転によりません。$\boldsymbol a$を球面の外に延長して3次元の回転$\nabla\times\boldsymbol a$を取ると、定ベクトル$\boldsymbol u,\boldsymbol w$について$(\nabla\times\boldsymbol a)\cdot(\boldsymbol u\times\boldsymbol w)=\boldsymbol w\cdot\partial_{\boldsymbol u}\boldsymbol a-\boldsymbol u\cdot\partial_{\boldsymbol w}\boldsymbol a$なので、$\boldsymbol t_1\times\boldsymbol t_2=\boldsymbol x$より$F=(\nabla\times\boldsymbol a)\cdot\boldsymbol x$です。勾配$\nablaφ$や回転$\nabla\times$の$\nabla$はベクトル解析の記号で、共変微分とは別のものです。右辺には球面に沿った方向の微分だけが現れるので、延長の仕方によりません。両極を除いた部分では$\boldsymbol a_+-\boldsymbol a_-$が関数$φ$の勾配で、勾配の回転は$0$なので、$\boldsymbol a_+$から計算した$F$と$\boldsymbol a_-$から計算した$F$は一致し、$F$は球面全体の関数です。

&&&prop 曲率の積分 [prop-curvature-integral]
$H^n$の捩りの条件を満たす任意の接続について、次が成り立ちます。

$$
\int_{S^2}F\,dΩ=-2\pi n=-2\pi\operatorname{wind}(z^n)
$$
&&&

&&&prf
以前の記事のストークスの定理を北半球$D_+$と南半球$D_-$に使う。法線は外向きの$\boldsymbol x$とする。$D_+$の境界の赤道は$φ$が増える向き、$D_-$の境界は$φ$が減る向きにたどられるので、次を得る。

$$
\int_{S^2}F\,dΩ=\oint\boldsymbol a_+\cdot d\boldsymbol x-\oint\boldsymbol a_-\cdot d\boldsymbol x=-n\oint\nablaφ\cdot d\boldsymbol x=-2\pi n
$$

周回積分は$φ$が増える向きの赤道に沿ったもので、$\oint\nablaφ\cdot d\boldsymbol x$は$φ$の1周の増加量$2\pi$である。
&&&

曲率の積分は接続の選び方によらず、貼り合わせ関数の回転数だけで決まります。以前の記事では回転数を偏角の変化量として定義し、積分による表示を注意で述べるに留め、曲率としての表示は扱いませんでした。貼り合わせの差$\boldsymbol a_+-\boldsymbol a_-$が偏角の勾配であることが、両者を結びつけています。[[7shi-kth3]]

符号が逆になるのは、本記事の接続の書き方と向きの規約のもとで、貼り合わせを$ψ_+=ψ_-z^n$としたために、接続の差が$-n\nablaφ$になることによります。

## 曲率が一定の接続

球面はどの向きにも対称なので、曲率が一定の接続が自然な候補です。[[prop-curvature-integral]]と球面の面積$4\pi$から、曲率が一定なら$F=-\frac n2$です。

$e_2$の軸の周りの回転で変わらない形$\boldsymbol a=f(x_2)\boldsymbol v_{01}$で探します。$f(x_2)(x_0e_1-x_1e_0)$をそのまま3次元のベクトル場と見て回転を計算すると、$\nabla\times\boldsymbol a=(-x_0f',\ -x_1f',\ 2f)$なので、単位球面上で次のようになります。

$$
F=-(1-x_2^2)f'+2x_2f=-\frac d{dx_2}\bigl((1-x_2^2)f\bigr)
$$

$F=c$（定数）とすると$(1-x_2^2)f=-cx_2+k$（$k$は定数）です。$U_+$で$f$が滑らかなら左辺は北極$x_2=1$で$0$なので$k=c$、同様に$U_-$では$k=-c$です。これから$f_+=\frac c{1+x_2}$、$f_-=-\frac c{1-x_2}$となります。$x_0^2+x_1^2=1-x_2^2$に注意すると、貼り合わせの条件は次のとおりです。

$$
\boldsymbol a_+-\boldsymbol a_-=\frac{2c}{1-x_2^2}\boldsymbol v_{01}=-n\frac{\boldsymbol v_{01}}{x_0^2+x_1^2}
$$

したがって$c=-\frac n2$で、[[prop-curvature-integral]]と一致します。

&&&fml 曲率が一定の接続 [fml-constant]
$$
\boldsymbol a_+=-\frac n2\frac{\boldsymbol v_{01}}{1+x_2},\qquad
\boldsymbol a_-=\frac n2\frac{\boldsymbol v_{01}}{1-x_2},\qquad
F=-\frac n2
$$
&&&

以下、この接続を使い、捩ったディラック作用素を$D_n$と書きます。$D_0=D_S$です。

$\boldsymbol v_{01}=\partial_φ\boldsymbol x$で、$|\boldsymbol v_{01}|^2=\sin^2θ$なので、$\boldsymbol a_+\cdot\partial_φ\boldsymbol x=-\frac n2(1-\cosθ)$、$\boldsymbol a_-\cdot\partial_φ\boldsymbol x=\frac n2(1+\cosθ)$です。前回の記事の命題「2次元の球面の枠による表示」の枠の側の成分$ψ_+=Uψ_f$に移します。[[7shi-kth6]]

右から$J$を掛ける演算は$U$と可換なので、$U_+$の側で次の形になります。

$$
\tilde UD_nU=f_1\left(\partial_θ+\frac{\cotθ}2\right)+\frac{f_2}{\sinθ}\bigl(\partial_φ+A(θ)R_J\bigr),\qquad A(θ)=-\frac n2(1-\cosθ)
$$

$R_J$は右から$J$を掛ける演算です。

&&&rem 射影から誘導される接続
$n=1$の接続は、ホップ束$H$の各ファイバーを$\mathbb C^2$の直線と見て、微分をその直線へ直交射影して得られる接続と一致します。$U_-$の単位ベクトル$\boldsymbol e=(α_-,β_-)$について、$\mathbb C^2$の内積で$\boldsymbol e$方向の成分を取ると、次のようになります。

$$
α_-^*\partial_φα_-+β_-^*\partial_φβ_-=\cos^2\fracθ2\,J,\qquad
α_-^*\partial_θα_-+β_-^*\partial_θβ_-=0
$$

したがって$\boldsymbol ec$の微分の射影は$\boldsymbol e\bigl(\partial_Xc+(\boldsymbol a_-\cdot X)cJ\bigr)$となり、$\cos^2\fracθ2=\frac12(1+\cosθ)$は$n=1$の$\boldsymbol a_-\cdot\partial_φ\boldsymbol x$に一致します。$H^n$の接続は、この接続を$n$倍したものです。
&&&

# リヒネロビッチの公式

## 曲率の項

前回の記事では$D_S^2=\nabla^*\nabla+\frac R4$（2次元の単位球面では$\frac R4=\frac12$）を示し、右辺が正であることから核が$0$になることを導きました。捩ると、ここに接続の曲率の項が加わります。[[7shi-kth6]]

捩った接続ラプラシアン$\nabla^{\boldsymbol a*}\nabla^{\boldsymbol a}$を、前回の記事と同じく局所的な正規直交枠で定めます。前回の記事の命題「回転の場による接続ラプラシアン」の証明は共変微分の線形性だけを使っているので、$\nabla^{\boldsymbol a*}\nabla^{\boldsymbol a}=-\sum_{a<b}\nabla^{\boldsymbol a}_{\boldsymbol v_{ab}}\nabla^{\boldsymbol a}_{\boldsymbol v_{ab}}$も成り立ちます。

&&&thm 捩ったディラック作用素のリヒネロビッチの公式 [thm-lich-twisted]
$H^n$の捩りの条件を満たす任意の接続について、次が成り立ちます。

$$
D_{\boldsymbol a}^2=\nabla^{\boldsymbol a*}\nabla^{\boldsymbol a}+\frac12-Fγ
$$

とくに曲率が一定の接続では$D_n^2=\nabla^{\boldsymbol a*}\nabla^{\boldsymbol a}+\frac12+\frac n2γ$です。
&&&

&&&prf
$\operatorname{div}\boldsymbol a=\sum_i\boldsymbol t_i\cdot\partial_{\boldsymbol t_i}\boldsymbol a$、$\nabla_{\boldsymbol a}=\sum_i(\boldsymbol a\cdot\boldsymbol t_i)\nabla_{\boldsymbol t_i}$と書く。

接続ラプラシアン。$\boldsymbol v=\boldsymbol v_{ab}$について、$R_J$は$\nabla_{\boldsymbol v}$と可換で$R_J^2=-1$なので、次のようになる。

$$
\nabla^{\boldsymbol a}_{\boldsymbol v}\nabla^{\boldsymbol a}_{\boldsymbol v}ψ=\nabla_{\boldsymbol v}\nabla_{\boldsymbol v}ψ+\partial_{\boldsymbol v}(\boldsymbol a\cdot\boldsymbol v)\,ψJ+2(\boldsymbol a\cdot\boldsymbol v)\nabla_{\boldsymbol v}ψJ-(\boldsymbol a\cdot\boldsymbol v)^2ψ
$$

$a<b$で和を取る。前回の記事の補題「回転の場の和」から$\sum(\boldsymbol a\cdot\boldsymbol v)^2=|\boldsymbol a|^2$、$\sum(\boldsymbol a\cdot\boldsymbol v)\nabla_{\boldsymbol v}=\nabla_{\boldsymbol a}$、$\sum(\partial_{\boldsymbol v}\boldsymbol a)\cdot\boldsymbol v=\operatorname{div}\boldsymbol a$である。また$\sum\partial_{\boldsymbol v}\boldsymbol v$の接方向の成分は$0$で、$\boldsymbol a$は接ベクトルなので$\sum\boldsymbol a\cdot\partial_{\boldsymbol v}\boldsymbol v=0$である。したがって次を得る。

$$
\nabla^{\boldsymbol a*}\nabla^{\boldsymbol a}=\nabla^*\nabla-(\operatorname{div}\boldsymbol a)R_J-2\nabla_{\boldsymbol a}R_J+|\boldsymbol a|^2
$$

ディラック作用素。$Tψ=\boldsymbol a\boldsymbol xψJ$とすると$D_{\boldsymbol a}=D_S+T$である。$(\boldsymbol a\boldsymbol x)^2=-|\boldsymbol a|^2$より$T^2=|\boldsymbol a|^2$である。接ベクトルの掛け算とのライプニッツ則から、次のようになる。

$$
D_STψ+TD_Sψ=\sum_i\boldsymbol t_i\boldsymbol x(\nabla_{\boldsymbol t_i}\boldsymbol a)\boldsymbol xψJ+\sum_i(\boldsymbol t_i\boldsymbol x\boldsymbol a\boldsymbol x+\boldsymbol a\boldsymbol x\boldsymbol t_i\boldsymbol x)\nabla_{\boldsymbol t_i}ψJ
$$

接ベクトル$\boldsymbol u,\boldsymbol w$について$\boldsymbol u\boldsymbol x\boldsymbol w\boldsymbol x=-\boldsymbol u\boldsymbol w$である。第2項は$-\sum_i(\boldsymbol t_i\boldsymbol a+\boldsymbol a\boldsymbol t_i)\nabla_{\boldsymbol t_i}ψJ=-2\nabla_{\boldsymbol a}ψJ$となる。第1項の係数は$-\sum_i\boldsymbol t_i\nabla_{\boldsymbol t_i}\boldsymbol a=-\operatorname{div}\boldsymbol a-\sum_i\boldsymbol t_i\wedge\nabla_{\boldsymbol t_i}\boldsymbol a=-\operatorname{div}\boldsymbol a-Fω\boldsymbol x$である（$\nabla_{\boldsymbol t_i}\boldsymbol a$は$\partial_{\boldsymbol t_i}\boldsymbol a$の接方向の成分で、面積要素の係数は変わらない）。$ω\boldsymbol xψJ=γψ$なので、次を得る。

$$
D_{\boldsymbol a}^2=D_S^2-(\operatorname{div}\boldsymbol a)R_J-Fγ-2\nabla_{\boldsymbol a}R_J+|\boldsymbol a|^2
$$

$D_S^2=\nabla^*\nabla+\frac12$と上の接続ラプラシアンの式を比べれば定理を得る。曲率が一定の接続では$F=-\frac n2$である。
&&&

## 片側での核の消失

$γ$は正のカイラリティで$1$、負のカイラリティで$-1$なので、曲率が一定の接続では、零階の項$\frac12+\frac n2γ$が次のように分かれます。

$$
D_n^2=\begin{cases}\nabla^{\boldsymbol a*}\nabla^{\boldsymbol a}+\dfrac{1+n}2&(\text{正のカイラリティ})\\[2ex]\nabla^{\boldsymbol a*}\nabla^{\boldsymbol a}+\dfrac{1-n}2&(\text{負のカイラリティ})\end{cases}
$$

片方のカイラリティでは曲率の項が$\frac R4=\frac12$に加わり、もう片方では$\frac12$から差し引かれます。

&&&prop 片側での核の消失 [prop-vanish]
$ψ$が$H^n$で捩ったスピノルで$D_nψ=0$とします。

1. $n\ge0$なら$ψ$は負のカイラリティで、$n\le0$なら正のカイラリティです。$n=0$なら$ψ=0$です。
2. $n=-1$なら$ψ$は平行（すべての接ベクトル$X$で$\nabla^{\boldsymbol a}_Xψ=0$）で、$\ker D_{-1}$は複素1次元以下です。$n=1$についても同様です。
&&&

&&&prf
$ψ^\pm=\frac12(ψ\pmγψ)$と置く。[[prop-twisted]]の1より$D_nψ^+$は負、$D_nψ^-$は正のカイラリティなので、$D_nψ=0$なら$D_nψ^\pm=0$である。前回の記事の命題「核の消失」の証明と同じく、$\nabla^{\boldsymbol a}_{\boldsymbol v_{ab}}$の反対称性から$(ψ^+,\nabla^{\boldsymbol a*}\nabla^{\boldsymbol a}ψ^+)=\sum_{a<b}\|\nabla^{\boldsymbol a}_{\boldsymbol v_{ab}}ψ^+\|^2$である。[[prop-twisted]]の2と[[thm-lich-twisted]]から次を得る。

$$
0=\|D_nψ^+\|^2=\sum_{a<b}\|\nabla^{\boldsymbol a}_{\boldsymbol v_{ab}}ψ^+\|^2+\frac{1+n}2\|ψ^+\|^2
$$

$n\ge0$なら$ψ^+=0$である。負のカイラリティについても同様で、$n\le0$なら$ψ^-=0$である。

$n=-1$では係数が$0$になり、$\nabla^{\boldsymbol a}_{\boldsymbol v_{ab}}ψ=0$がすべての$a<b$で成り立つ。前回の記事の補題「回転の場の和」から、各点で$\boldsymbol v_{ab}$は接平面を張るので、$ψ$は平行である。平行なスピノルは長さが一定なので（内積の保存）、ある点で$0$なら至るところ$0$である。したがって、点$\boldsymbol x_0$での値を取る写像は、平行なスピノルの空間から$\boldsymbol x_0$での正のカイラリティの空間への複素線形な単射である。後者は複素1次元なので、平行なスピノルの空間も複素1次元以下である。$n=1$は負のカイラリティで同様である。
&&&

捩らない場合は、どちらのカイラリティでも$\frac R4=\frac12>0$が核を消しました（前回の記事）。捩ると、$|n|=1$では片方のカイラリティで曲率の項が$\frac R4$をちょうど打ち消し、核は平行なスピノルになります。$|n|\ge2$では曲率の項が$\frac R4$を上回り、曲率の項を差し引く側については、この非負性の議論だけでは核の有無や次元を決められません。その側の核の次元は、次の節で直接計算します。

# 核の計算

## 変数分離

正のカイラリティの核を、枠の形で求めます。両極を除いた部分で$ψ_+=Uψ_f$と置くと、正のカイラリティでは$ψ_f=f_1v$（$v$は$\mathbb C_J$値）です。$U$は$φ$の1周で符号を変えるので、$ψ_+$が1価なら$v$は$φ$の1周で符号を変えます（前回の記事）。したがって$v$は、$μ$を半整数として$e^{Jμφ}$でフーリエ展開されます。[[7shi-kth6]]

$f_1Jf_1=J$と、$\mathbb C_J$の元が$J$と可換であることを使うと、枠の形の$D_n$は次のようになります。

$$
\tilde UD_nU(f_1v)=-\left(\partial_θ+\frac{\cotθ}2\right)v+\frac J{\sinθ}\bigl(\partial_φv+A(θ)vJ\bigr)
$$

$v=g(θ)e^{Jμφ}$とすると、$J(Jμ+AJ)=-(μ+A)$より、$D_n(Uf_1v)=0$は$g$についての1階の常微分方程式になります。

$$
g'+\frac{\cotθ}2g+\frac{μ+A(θ)}{\sinθ}g=0
$$

&&&prop モードの解 [prop-modes]
上の方程式の解は、定数倍を除いて次の関数です。

$$
g_μ(θ)=\sin^{-μ-\frac12}\fracθ2\,\cos^{μ-\frac12-n}\fracθ2
$$

$g_μ$が$0<θ<π$で有界であるのは、$n+\frac12\leμ\le-\frac12$のときに限ります。そのような半整数$μ$の個数は、$n\le-1$なら$-n$、$n\ge0$なら$0$です。
&&&

&&&prf
$s=\sin\fracθ2$、$c=\cos\fracθ2$とすると$\sinθ=2sc$、$\cotθ=\frac{c^2-s^2}{2sc}$、$A=-ns^2$である。$\frac{g_μ'}{g_μ}=-\left(μ+\frac12\right)\frac c{2s}-\left(μ-\frac12-n\right)\frac s{2c}$を通分すると$\frac{-2μ-c^2+s^2+2ns^2}{4sc}$となり、$-\frac{\cotθ}2-\frac{μ+A}{\sinθ}=\frac{-(c^2-s^2)-2μ+2ns^2}{4sc}$と一致する。1階の線形の方程式なので、解は定数倍を除いて一意である。$θ\to0$で$s\to0$、$θ\toπ$で$c\to0$なので、有界であるのは2つの指数がともに$0$以上のときである。
&&&

## 核の元

有界なモードが、実際に球面全体で滑らかな捩ったスピノルを与えることを確かめます。鍵は、[[prop-hopf-chirality]]の$Φ$で局所的な単位ベクトル$(α_+,β_+)$を写したものが、枠の形で単純に書けることです。$e_2e_0f_1=1$と$e^{-Jφ/2}f_1=f_1e^{Jφ/2}$から、次のようになります。

$$
Φ(α_+,β_+)=\cos\fracθ2\,f_1-\sin\fracθ2\,e^{-Jφ}=Uf_1e^{-Jφ/2}
$$

&&&thm 捩ったディラック作用素の核 [thm-kernel]
$n=-m$（$m\ge1$）とします。$P(α,β)$を$\mathbb C_J$係数の次数$m-1$の同次多項式として、次の組は$D_n$の核に入る正のカイラリティのスピノルです。

$$
ψ_\pm=Φ(α_\pm,β_\pm)\,P(α_\pm,β_\pm)
$$

$\ker D_n^+$の元はすべてこの形で、$\ker D_n^-=0$です。$H^n$で捩った場合のそれぞれの複素次元は次のとおりです。

$$
\dim_{\mathbb C}\ker D_n^+=\max(-n,0),\qquad\dim_{\mathbb C}\ker D_n^-=\max(n,0)
$$
&&&

&&&prf
$n=-m\le-1$とする。

捩ったスピノルであること。$(α_\pm,β_\pm)$は$U_\pm$で滑らかで、$Φ$は定まった線形写像なので、$ψ_\pm$は$U_\pm$で滑らかである。$(α_+,β_+)=(α_-,β_-)e^{-Jφ}$と、$Φ$の複素線形性、$P$の同次性から、$ψ_+=ψ_-e^{-Jφ}e^{-(m-1)Jφ}=ψ_-e^{nJφ}$である。[[prop-hopf-chirality]]により正のカイラリティである。

核に入ること。単項式$P=α^{m-1-p}β^p$（$0\le p\le m-1$）について、$α_+^{m-1-p}β_+^p=\cos^{m-1-p}\fracθ2\sin^p\fracθ2\,e^{-pJφ}$なので、$ψ_+=Uf_1g_μe^{Jμφ}$（$μ=-p-\frac12$）である。$-μ-\frac12=p$、$μ-\frac12-n=m-1-p$なので、[[prop-modes]]の$g_μ$そのものであり、両極を除いて$D_nψ_+=0$となる。$D_n$は貼り合わせと両立するので、重なりでは$D_nψ_-=(D_nψ_+)e^{-nJφ}=0$でもある。$D_nψ_\pm$は$U_\pm$で連続なので、北極でも南極でも$0$である。

すべてであること。$\ker D_n^+$の元$ψ$を$ψ_+=Uf_1v$と書き、$v$の$φ$についてのフーリエ係数を$v_μ(θ)$とする。部分積分により$v_μ$は[[prop-modes]]の方程式を満たすので$v_μ=c_μg_μ$である。$U$は単位の回転子で、貼り合わせは長さを保つから$|v|=|ψ|$であり、$|ψ|$は球面上の連続関数なので有界である。したがって$|v_μ|$も有界で、$c_μ\ne0$となるのは有界な$g_μ$に限る。これは上の単項式の$p=0,\dots,m-1$に対応する$m$個である。フーリエ級数の完全性により、$v$はこれらのモードの和で、$ψ$は上の形になる。

次元。単項式に対応するモードは互いに異なる$μ$を持つので、$m$個の単項式は独立であり、$\dim_{\mathbb C}\ker D_n^+=m$である。

$n\ge0$では、[[prop-vanish]]の1より$\ker D_n^+=0$である。負のカイラリティについては、[[prop-flip]]を$H^{-n}$で捩ったスピノル（接続$-\boldsymbol a_\pm$）に使うと、$ψ\mapstoψK$が$\ker D_{-n}^+$を$\ker D_n^-$に写し、複素次元を保つので、$\dim\ker D_n^-=\dim\ker D_{-n}^+=\max(n,0)$である。
&&&

$n=0$では両方の核が$0$で、前回の記事の結果に戻ります。$n=-1$では$P$は定数で、核は$Φ(α,β)$の定数倍です。$|Φ(α,β)|^2=|α|^2+|β|^2=1$で、[[prop-vanish]]の2のとおり平行です。どこでも$0$にならない$Φ(α_\pm,β_\pm)$が$H^{-1}$で捩った正のカイラリティのスピノルになることは、[[cor-spinor-hopf]]の$H\otimes H^{-1}$が自明であることの言い換えです。

## 同次多項式

[[thm-kernel]]の核の元は、ホップファイブレーション$S^3\to S^2$の言葉で読めます。$S^3$の点$(α,β)$が$\boldsymbol x$の上にあるとき、同じファイバーの点は$(α,β)e^{Jt}$で、$Φ(α,β)P(α,β)$は$e^{mJt}$倍になります。$U_+$と$U_-$で$e^{Jt}=e^{-Jφ}$だけ異なる点を選んだことが、貼り合わせ$ψ_+=ψ_-e^{-mJφ}$として現れています。核の元は、$S^3$上の関数$Φ(α,β)P(α,β)$のうち、ファイバーを1点ずつ選んで球面に降ろしたものです。

南半球の座標で書くと、$(α_-,β_-)=(z,1)/\sqrt{1+|z|^2}$です。$z=\cot\frac θ2e^{Jφ}$は以前の記事の比$α/β$です。[[7shi-homog]]

したがって次のようになります。

$$
P(α_-,β_-)=\frac{P(z,1)}{(1+|z|^2)^{(m-1)/2}}
$$

$P(z,1)$は$z$の$m-1$次以下の多項式で、その全体は$1,z,\dots,z^{m-1}$で張られます。以前の記事の言葉では、$H^{-m}$で捩ったディラック作用素の核は、ハーディ空間の$1,z,\dots,z^{m-1}$で張られる部分と同じ形をしています。後者は$T_{z^{-m}}=(S^*)^m$の核です。[[7shi-kth5]]

[[cor-spinor-hopf]]から、正のカイラリティで$H^n$に捩ったスピノルの束は$H\otimes H^n=H^{n+1}$です。$n=-m$なら$H^{-(m-1)}$で、その切断のうち核に入るものが、$\mathbb C^2$の座標$α,β$の次数$m-1$の同次多項式で与えられます。前回の記事では、球面上の関数を$\mathbb R^3$の多項式（球面モノジェニックス）から作りました。捩ったスピノルは球面全体の関数ではないので、多項式の構造は$\mathbb R^3$の側ではなく、ホップファイブレーションの上の$\mathbb C^2$の側に現れます。[[7shi-kth6]]

# 指数と回転数

[[thm-kernel]]から、指数が求まります。

&&&cor 捩ったディラック作用素の指数 [cor-index]
$$
\operatorname{ind}D_n=\max(-n,0)-\max(n,0)=-n
$$
&&&

以前の記事の量と並べます。$H^n$の貼り合わせ関数は$z^n$で、回転数は$n$でした。テプリッツ作用素の指数は$\operatorname{ind}T_{z^n}=-n$でした。[[7shi-kth3]][[7shi-kth5]]

曲率の積分は$-2\pi n$です（[[prop-curvature-integral]]）。したがって次が成り立ちます。

$$
\operatorname{ind}D_n=\frac1{2\pi}\int_{S^2}F\,dΩ=-\operatorname{wind}(z^n)=\operatorname{ind}T_{z^n}
$$

$D_n$と$T_{z^n}$は、球面上の微分作用素と円周上のハーディ空間の作用素で、まったく別のものです。それでも、同じ貼り合わせ関数$z^n$から同じ指数$-n$が得られます。以前の記事の準同型$φ$で書くと、どちらも$[H^n]-1=n([H]-1)$に$-φ([H^n]-1)=-n$を対応させます。[[7shi-kth4]][[7shi-kth5]]

指数は核の次元の差なので、核がどちらのカイラリティに現れるかは指数だけからは読めません。本記事の接続では、核は$n$の符号に応じて正と負のカイラリティのどちらか一方だけに現れ、指数の絶対値がそのまま核の次元になります。

&&&rem 符号の規約
指数の符号は、どちらのカイラリティを正と呼ぶかで変わります。$\boldsymbol t_1\boldsymbol t_2\boldsymbol x=ω$となる向きの接ベクトルについて、接空間の面積要素に虚数単位を掛けた$i\,\boldsymbol t_1\bullet\boldsymbol t_2\bullet$でカイラリティを定める規約もあります。これは本記事の$-γ$にあたり、その規約では指数は$+n$です。本記事では、以前の記事のテプリッツ作用素の指数と向きが揃うように$γ$を選びました。[[7shi-kth5]]

曲率の積分と回転数の符号が逆になるのは、本記事の規約（接続を$\nabla+(\boldsymbol a\cdot X)R_J$と書き、曲率を$\boldsymbol t_1\boldsymbol t_2\boldsymbol x=ω$の向きで定めること）のもとで、貼り合わせを$ψ_+=ψ_-z^n$と定めたためです。[[prop-curvature-integral]]の証明で見たとおり、接続の差が$-n\nablaφ$になります。
&&&

# まとめ

2次元の球面上のディラック作用素を、スピノルを貼り合わせ関数$z^n$で貼り合わせることで$H^n$に捩りました。接平面の面積要素に複素構造を組み合わせたカイラリティ$γ$でスピノルを分けると、正のカイラリティの部分は写像$Φ(α,β)=f_1α-β$でホップ束$H$と一致しました。貼り合わせと両立する接続の曲率の積分は、接続によらず回転数で決まり、曲率が一定の接続は$e_2$の軸の周りの対称性から求まりました。リヒネロビッチの公式に曲率の項が加わり、片方のカイラリティで核が消えます。変数分離で残りの核を求めると、核は$\mathbb C^2$の座標の同次多項式で書け、指数は$-n$となってテプリッツ作用素の指数と一致しました。

&&& カイラリティ
$$
γψ=ω\boldsymbol xψJ,\qquad γ^2=1,\qquad D_Sγ=-γD_S
$$
正のカイラリティの空間は$Φ(E_{\boldsymbol x})$で、束として$H$と同型です。負のカイラリティは$H^{-1}$と同型です。
&&&

&&& 捩ったスピノルと接続
$$
ψ_+=ψ_-e^{nJφ},\qquad\nabla^{\boldsymbol a}_Xψ=\nabla_Xψ+(\boldsymbol a\cdot X)ψJ,\qquad\boldsymbol a_+-\boldsymbol a_-=-n\nablaφ
$$
&&&

&&& 曲率
曲率の積分は接続によりません。曲率が一定の接続では$F=-\frac n2$です。
$$
\int_{S^2}F\,dΩ=-2\pi n
$$
&&&

&&& リヒネロビッチの公式
$$
D_n^2=\nabla^{\boldsymbol a*}\nabla^{\boldsymbol a}+\frac12+\frac n2γ
$$
&&&

&&& 核と指数
曲率が一定の接続で考えます。$n=-m\le-1$では、$P$を次数$m-1$の同次多項式として、核は次の形のスピノルです。
$$
ψ=Φ(α,β)P(α,β)
$$
$n\ge1$では核は負のカイラリティにあり、指数は次のとおりです。
$$
\operatorname{ind}D_n=-n=\operatorname{ind}T_{z^n}
$$
&&&

| 量 | $H^n$（貼り合わせ関数$z^n$）での値 |
|:---|:---|
| 回転数$\operatorname{wind}(z^n)$ | $n$ |
| 曲率の積分$\frac1{2\pi}\int F\,dΩ$ | $-n$ |
| テプリッツ作用素の指数$\operatorname{ind}T_{z^n}$ | $-n$ |
| 捩ったディラック作用素の核（正、負のカイラリティ） | $\max(-n,0)$、$\max(n,0)$ |
| 捩ったディラック作用素の指数$\operatorname{ind}D_n$ | $-n$ |
