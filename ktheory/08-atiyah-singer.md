アティヤ＝シンガーの指数定理を主張として述べ、これまで球面で計算した例がその特殊な場合にあたることを確かめます。

# 概要

前回の記事では、2次元の球面$S^2$上のディラック作用素をホップ束の冪$H^n$で捩り、曲率が一定の接続について核を求めました。[[7shi-kth7]]

指数は次のようになりました。

$$
\operatorname{ind}D_n=\frac1{2\pi}\int_{S^2}F\,dΩ=-\operatorname{wind}(z^n)=\operatorname{ind}T_{z^n}=-n
$$

捩ったディラック作用素$D_n$は、どの$n$でも1階の部分が同じで、異なるのは零階の項とスピノルを貼り合わせる関数だけです。それでも指数は$n$によって変わり、その値は貼り合わせ関数の回転数、すなわち束のねじれで決まりました。

以前の記事では、アティヤ＝シンガーの指数定理を、閉多様体上の作用素の指数（解析的指数）が作用素の最高階の部分（主表象）から位相的に定まる量（位相的指数）に等しいという主張として紹介しました。また、球面上のディラック作用素の1階の部分から、主表象として接ベクトルを掛ける演算$\boldsymbol\xi\mapsto\boldsymbol\xi\bullet$を取り出しました。[[7shi-kth1]][[7shi-kth6]]

本記事で扱う問いは次のとおりです。

> 閉多様体上のディラック作用素の指数は、何から決まるか。球面で計算してきた指数は、一般の主張の中でどう位置づけられるか。

指数定理は、解析の側で直接数えにくい量を、位相の側の計算に置き換えます。たとえば次のような使い方があります。

- **解の個数**：方程式を解かずに、核と余核の次元の差が分かります。余核が消えることも分かれば、解の空間の次元が求まります。前回の記事で求めた核が同次多項式の個数で数えられたことは、その最も簡単な例です（リーマン＝ロッホの定理）。
- **曲率の制約**：以前の記事では、球面ではスカラー曲率が正であることから、ディラック作用素の核が消えることを示しました。これを裏返すと、位相的な量$\int_M\hat A(M)$が$0$でない多様体には、スカラー曲率が至るところ正の計量が入らないことが分かります。
- **ベクトル場の障害**：オイラー標数が$0$でない閉多様体には、どこでも$0$にならない接ベクトル場がありません。以前の記事で主張に留めた毛玉の定理が、この形で得られます。

前提は次のとおりです。

- **前回の記事から**：カイラリティ$γψ=ω\boldsymbol xψJ$と、右から$J=e_0e_1$を掛ける複素構造、$Φ(α,β)=f_1α-β$、正負のカイラリティのスピノルの束が$H$と$H^{-1}$に同型であること、捩ったスピノル$ψ_+=ψ_-e^{nJφ}$と接続$\nabla^{\boldsymbol a}_X=\nabla_X+(\boldsymbol a\cdot X)R_J$、曲率$F$、枠の形の$D_n$、核と指数$\operatorname{ind}D_n=-n$、指数の定義$\dim_{\mathbb C}\ker D^+-\dim_{\mathbb C}\ker D^-$[[7shi-kth7]]
- **球面上のディラック作用素**：$\operatorname{Cl}_{n,0}(\mathbb R)$の偶部分代数に値を取るスピノル、接ベクトルの掛け算$X\bulletψ=X\boldsymbol xψ$、主表象$\boldsymbol\xi\bullet$、球面のラプラシアン$Δ_S$、回転子$U$、リヒネロビッチの公式と核の消失[[7shi-kth6]]
- **K群**：加群$W$の上のパラベクトル$g_W$を貼り合わせ関数とする束$E_W$、延長できる加群が自明な束を与えること、$\tilde K(S^2)\cong\mathbb Z$、アティヤ＝ボット＝シャピロの定理（主張）[[7shi-kth4]]
- **貼り合わせと指数**：貼り合わせの向きとホップ束$H$の回転数$+1$、複素直線束として$TS^2\cong H^{-2}$、テプリッツ作用素の指数$\operatorname{ind}T_f=-\operatorname{wind}(f)$[[7shi-kth3]][[7shi-kth5]]

本記事は次の順に進みます。

1. 主表象を定め、ディラック作用素の主表象を正のカイラリティから負のカイラリティへの写像と見ると、各点で貼り合わせ関数になり、既約なスピノルでは球面のK群の生成元を与えることを示します。これを踏まえて位相的指数を言葉で導入します。
2. 指数定理を主張として述べ、$S^2$では前回の記事の指数$-n$が、向きとカイラリティを揃えたうえで$H^n$の第1チャーン数に一致することを確かめます。四元数のホップ束の例を挙げ、以前の記事で主張に留めた安定な非自明性（複素ベクトル束として）が指数で検出できることを示します。
3. ガウス＝ボネの定理・符号数定理・リーマン＝ロッホの定理・$\hat A$種数への特殊化を表にまとめ、$S^2$で具体的に確かめられるものを確かめます。
4. 消えないベクトル場が主表象の変形を与えることから、偶数次元の球面の毛玉の定理を導きます。

証明するのは、主表象と貼り合わせ関数の対応、$S^2$での符号の照合、$S^2$の調和形式によるオイラー標数の計算、捩ったディラック作用素の核の方程式がコーシー＝リーマン方程式になること、主張を認めたうえで導ける帰結（四元数のホップ束の複素ベクトル束としての安定な非自明性の検出、正のスカラー曲率と$\hat A$種数、毛玉の定理）です。指数定理そのもの、位相的指数の定義とその性質、特性類（$\operatorname{ch}$・$\hat A$など）の一般論、ホッジ理論による$d+δ$の指数とオイラー標数の関係、$S^2$上で回転が$0$の接ベクトル場が勾配であること（閉1形式の完全性）は主張として使います。指数定理の証明（熱核による方法、同境による方法、K理論による方法）は、本記事では扱いません。

# 主表象

## 最高階の部分

前回の記事の$D_n$は、どの$n$でも$D_S$に零階の項$\boldsymbol a\boldsymbol xψJ$を加えたものでした。零階の項は、接続の選び方によっても変わります。指数がこれらの細部によらず束のねじれだけで決まるなら、作用素の側で指数に関わるのは、最高階の部分だけのはずです。そこで、1階の微分作用素から最高階の部分を取り出します。[[7shi-kth7]]

&&&def 主表象
多様体上の1階の微分作用素$P$が、各点の近くで、接空間の正規直交基底$\boldsymbol t_i$と線形写像$A_i$を使って次の形に書けるとする。

$$
Pψ=\sum_iA_i\,\partial_{\boldsymbol t_i}ψ+(\text{零階の項})
$$

計量で余接ベクトルを接ベクトルと同一視し、接ベクトル$\boldsymbol\xi$に対して次の線形写像を$P$の**主表象**と呼ぶ。

$$
σ_P(\boldsymbol\xi)=\sum_i(\boldsymbol\xi\cdot\boldsymbol t_i)\,A_i
$$

$\boldsymbol\xi\ne0$で$σ_P(\boldsymbol\xi)$が可逆なとき、$P$を**楕円型**と呼ぶ。
&&&

微分$\partial_{\boldsymbol t_i}$を成分$\boldsymbol\xi\cdot\boldsymbol t_i$に置き換えたもので、以前の記事の規約（$i\xi_a$ではなく$\xi_a$に置き換える）に合わせています。$i\xi_a$に置き換える通常の規約との違いは主表象を定数$i$倍することだけで、以下で扱うK群の元は変わりません。零階の項は主表象に現れません。[[7shi-kth1]]

&&&ex 球面上のディラック作用素
以前の記事で見たとおり、$D_S=\sum_i\boldsymbol t_i\bullet\nabla_{\boldsymbol t_i}$の主表象は$σ(\boldsymbol\xi)=\boldsymbol\xi\bullet$である。$(\boldsymbol\xi\bullet)^2=(\boldsymbol\xi\boldsymbol x)^2=-|\boldsymbol\xi|^2$なので、$\boldsymbol\xi\ne0$なら可逆で、$D_S$は楕円型である。前回の記事の$D_{\boldsymbol a}=D_S+\boldsymbol a\boldsymbol xR_J$は零階の項だけが異なるので、どの$n$、どの接続でも主表象は$\boldsymbol\xi\bullet$である。[[7shi-kth6]][[7shi-kth7]]
&&&

## カイラリティとの関係

前回の記事の命題「カイラリティの性質」により、接ベクトルの掛け算はカイラリティ$γ$と反交換します。したがって$\boldsymbol\xi\bullet$は、各点で正のカイラリティの空間$S^+_{\boldsymbol x}$を負のカイラリティの空間$S^-_{\boldsymbol x}$に写します。どちらも複素1次元なので、$\boldsymbol\xi$を接空間の単位円の上で動かすと、$S^+_{\boldsymbol x}$から$S^-_{\boldsymbol x}$への線形同型の族、すなわち単位円から$\mathbb C^\times$への写像が得られます。その回転数を調べます。[[7shi-kth7]]

&&&prop 主表象の回転数 [prop-symbol-wind]
$S^2$の点$\boldsymbol x$で、$\boldsymbol t_1\boldsymbol t_2\boldsymbol x=ω$となる正規直交な接ベクトルを取り、$\boldsymbol\xi=\cos t\,\boldsymbol t_1+\sin t\,\boldsymbol t_2$とする。$S^\pm_{\boldsymbol x}$の単位ベクトル$u_\pm$を1つずつ選ぶと、$c(t)\in\mathbb C_J$によって$\boldsymbol\xi\bullet u_+=u_-\,c(t)$と書け、$|c(t)|=1$で、$t\mapsto c(t)$の回転数は$-1$である。
&&&

&&&prf
$u_\pm$を別の単位ベクトルに取り替えると、$c(t)$には$t$によらない長さ$1$の複素数が掛かるだけなので、回転数は変わらない。$\boldsymbol t_1,\boldsymbol t_2$を向きを保って取り替えると、$t$が定数だけずれるだけなので、これも回転数を変えない。各点の近くでは、接ベクトルと$u_\pm$を$\boldsymbol x$について連続に選べるので、$c(t)$は$\boldsymbol x$について連続に変わり、回転数は局所的に一定である。球面は連結なので、両極を除いた点で示せばよい。

両極を除いた点では、前回の記事と同じく$\boldsymbol t_1=\hatθ$、$\boldsymbol t_2=\hatφ$と取れる（$\hatθ\hatφ\boldsymbol x=ω$）。回転子$U$で枠の側に移すと$\tilde U\hatθ\boldsymbol xU=f_1$、$\tilde U\hatφ\boldsymbol xU=f_2=f_1J$であり、$S^+_{\boldsymbol x}=Uf_1\mathbb C_J$、$S^-_{\boldsymbol x}=U\mathbb C_J$である。$u_+=Uf_1$、$u_-=U$とする。$f_1^2=-1$、$f_1Jf_1=J$より、次のようになる。

$$
\boldsymbol\xi\bullet Uf_1=U(\cos t\,f_1+\sin t\,f_1J)f_1=U(-\cos t+J\sin t)=-U\,e^{-Jt}
$$

したがって$c(t)=-e^{-Jt}$であり、回転数は$-1$である。[[7shi-kth6]][[7shi-kth7]]
&&&

接空間$T_{\boldsymbol x}S^2\cong\mathbb R^2$の単位円板を2つ用意し、単位円で貼り合わせると2次元の球面になります。$c(t)$はその赤道の上の貼り合わせ関数で、回転数$\pm1$なので、以前の記事の$\tilde K(S^2)\cong\mathbb Z$の生成元を与えます。主表象は、各点の接空間の上で、K群の生成元を与える貼り合わせ関数になっています。回転数の符号は、単位円をたどる向き（$\boldsymbol t_1$から$\boldsymbol t_2$へ）によります。[[7shi-kth4]]

## 貼り合わせ関数との一致

同じことが一般の次元で成り立ちます。点を1つ固定すると、接ベクトルの掛け算は接空間の上のクリフォード代数の作用なので、代数の言葉だけで述べられます。

以前の記事では、$S^k$上の束を、$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群$W$の上のパラベクトル$g_W(\boldsymbol x)=x_0I+x_1J_1+\dots+x_{k-1}J_{k-1}$を貼り合わせ関数として作りました。一方、ディラック作用素のスピノルは、カイラリティで2つに分かれ、接ベクトルの掛け算は一方を他方に写します。そこで、$\operatorname{Cl}_{0,k}(\mathbb R)$の生成元$e_1,\dots,e_k$（$e_a^2=-1$、互いに反交換する）が、ベクトル空間$S=S^+\oplus S^-$の上で、$S^+$と$S^-$を入れ替える線形写像として作用しているとします。これが、$k$次元の接空間の上の接ベクトルの掛け算とカイラリティの代数的なモデルです。[[7shi-kth4]]

&&&prop 主表象と貼り合わせ関数 [prop-symbol-clutching]
$h_l=e_1^{-1}e_l$（$l=2,\dots,k$）は$S^+$を保ち、$h_l^2=-1$を満たして互いに反交換する。したがって$S^+$は、$h_l$を生成元の作用とする$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群である。$\boldsymbol\xi=\sum_a\xi_ae_a$について、次が成り立つ。

$$
e_1^{-1}\boldsymbol\xi\big|_{S^+}=\xi_1+\xi_2h_2+\dots+\xi_kh_k=g_{S^+}(\xi_1,\dots,\xi_k)
$$
&&&

&&&prf
$e_1^{-1}$と$e_l$はどちらも$S^+$と$S^-$を入れ替えるので、積$h_l$は$S^+$を保つ。$e_1^{-1}=-e_1$と反交換性から、$h_l^2=e_1e_le_1e_l=-e_1^2e_l^2=-1$であり、$l\ne l'$なら$h_lh_{l'}=e_1e_le_1e_{l'}=-e_1^2e_le_{l'}=e_le_{l'}$、同様に$h_{l'}h_l=e_{l'}e_l$なので、$h_lh_{l'}=-h_{l'}h_l$である。最後の式は$e_1^{-1}e_1=1$と$h_l$の定義による。
&&&

$e_1^{-1}$は$S^-$から$S^+$への定まった同型なので、$\boldsymbol\xi\mapsto\boldsymbol\xi|_{S^+}$と$\boldsymbol\xi\mapsto g_{S^+}(\boldsymbol\xi)$は、値域の同一視を除いて同じ写像です。以前の記事では、ディラック作用素の主表象と球面の生成元を与える貼り合わせ関数を「同じ形の写像」と述べ、作用させる加群を適切に選べば対応させられると述べました。その加群は、正のカイラリティの空間$S^+$です。以前の記事の記号との対応は$x_0=\xi_1$、$x_j=\xi_{j+1}$、$J_j=h_{j+1}$です。[[7shi-kth1]][[7shi-kth4]]

貼り合わせ関数がK群の生成元を与えるかどうかは、加群$S$によります。$S$として偶数次元の既約な複素スピノルを取ると、この貼り合わせ関数は複素K群の生成元を与えます。$k=2$では$S^+$は複素1次元で、[[prop-symbol-wind]]の$c(t)$になります。$k=4$では$S^+$は$\operatorname{Cl}_{0,3}(\mathbb R)\cong2\mathbb H$の既約加群で、貼り合わせ関数は四元数のホップ束のもの（またはその共役）です。一般の次数付きの加群（既約なスピノルの直和や、後で扱う微分形式全体など）では、その加群に対応するK群の元を与え、それは生成元の倍数や$0$になることもあります。

## 位相的指数

主表象は、各点の接空間の上で球面のK群の元を定めました。これを多様体全体で集めると、余接束の上のK群（2つの束と、コンパクト集合の外での両者の同型を扱うもの）の元が得られることが知られています。主表象の定義域と値域の束（正負のカイラリティのスピノルの束）と、$\boldsymbol\xi\ne0$での同型$σ(\boldsymbol\xi)$を組にしたものです。多様体をユークリッド空間に埋め込み、ボット周期性を使ってこの元から整数を取り出す手続きがあり、その整数を**位相的指数**と呼びます。捩ったディラック作用素の主表象は$σ_{D_E}(\boldsymbol\xi)=\boldsymbol\xi\bullet\otimes\operatorname{id}_E$で、定義域と値域の束に$E$のテンソル積が入ることから、$E$のねじれがK群の元に入ります。余接束の上のK群と位相的指数の定義には、本記事では立ち入りません。

本記事で使うのは、次の3つの性質です。1つ目は、位相的指数が主表象の定めるK群の元だけで決まり、零階の項によらないことです。2つ目は、コンパクト集合の外での可逆性を保つ連続変形では、主表象の定めるK群の元が変わらないことです。以前の記事で、貼り合わせ関数をホモトピーで変形しても同型な束が得られたことにあたります。3つ目は、以前の記事で延長できる加群が自明な束を与えたのと同様に、$\boldsymbol\xi=0$を含めてどこでも可逆な写像が定めるK群の元は$0$で、位相的指数も$0$になることです。いずれも主張として使います。[[7shi-kth3]][[7shi-kth4]]

&&&rem テプリッツ作用素
以前の記事のテプリッツ作用素の指数$\operatorname{ind}T_f=-\operatorname{wind}(f)$も、記号$f$（貼り合わせ関数）の位相的な量から指数が決まるという同じ形の定理です。テプリッツ作用素は円周上の作用素なので、閉多様体上の微分作用素についての指数定理とは設定が異なりますが、その一種として扱えることが知られています。本記事では扱いません。[[7shi-kth5]]
&&&

# 指数定理

## 主張

指数定理を述べるための設定を整えます。$M$を向き付けられた$2m$次元の閉リーマン多様体で、スピン構造を持つものとします。スピン構造は、各点の既約なスピノルの空間を、接ベクトルの正規直交な組の取り替えと両立するように、多様体全体で貼り合わせるための構造です。存在の条件や構成には、本記事では立ち入りません。各点の向きに合った正規直交な接ベクトル$\boldsymbol t_1,\dots,\boldsymbol t_{2m}$について、スピノルの上の次の作用素は$2$乗が$1$で、ディラック作用素と反交換することが知られています。

$$
γ=i^m\,\boldsymbol t_1\bullet\cdots\boldsymbol t_{2m}\bullet
$$

これを標準的なカイラリティとし、固有値$\pm1$の空間を正負のカイラリティとします。$M$上のエルミート計量を持つ複素ベクトル束$E$と、その計量を保つ接続を取ると、前回の記事と同様に、スピノルの束と$E$のテンソル積の切断に作用する捩ったディラック作用素$D_E$が定まり、$γ$と反交換します。$D_E$は対称で、$D_E^-$は$D_E^+$の随伴になります。指数は前回の記事と同じく$\operatorname{ind}D_E=\dim_{\mathbb C}\ker D_E^+-\dim_{\mathbb C}\ker D_E^-$とします。[[7shi-kth7]]

一般の楕円型微分作用素$P$の解析的指数は$\dim\ker P-\dim\ker P^*$（$P^*$は随伴）で、閉多様体の上ではどちらの核も有限次元であることが知られています。$D_E^+$では、これが上の指数になります。

&&&thm アティヤ＝シンガーの指数定理
1. 閉多様体上の楕円型微分作用素について、解析的指数は位相的指数に等しくなる。
2. 上の設定で、次が成り立つ。

$$
\operatorname{ind}D_E=\int_M\hat A(M)\operatorname{ch}(E)
$$
&&&

証明は本記事では扱いません。2は、1の位相的指数を微分形式の積分で書き表したものです。

右辺の$\operatorname{ch}(E)$（**チャーン指標**の形式）と$\hat A(M)$（**$\hat A$形式**）は、それぞれ$E$の接続の曲率と$M$の曲率から作る微分形式（特性形式）で、いくつかの次数の部分の和です。$\hat A(M)$の最高次の部分の積分$\int_M\hat A(M)$を**$\hat A$種数**と呼びます。積分に効くのは、$2m$次の部分だけです。その積分は接続や計量の選び方によらず、束と多様体だけで決まることが知られています。低い次数の部分は次のとおりです。

$$
\operatorname{ch}(E)=\operatorname{rk}E+c_1(E)+\operatorname{ch}_2(E)+\cdots,\qquad
\hat A(M)=1-\frac{p_1(M)}{24}+\cdots
$$

$c_1(E)$は2次の**第1チャーン形式**、$\operatorname{ch}_2(E)=\frac12c_1(E)^2-c_2(E)$は4次の部分（$c_2$は第2チャーン形式）、$p_1(M)$は4次の**第1ポントリャーギン形式**です。$\hat A(M)$は$4$の倍数の次数の部分だけを持ちます。特性類の一般論には、本記事では立ち入りません。

とくに2次元（$m=1$）では、$\hat A(M)$の$0$次の部分$1$だけが効き、次のようになります。

$$
\operatorname{ind}D_E=\int_Mc_1(E)
$$

複素直線束$E$の接続が、局所的な切断$s$と実の1次微分形式$A$によって$\nabla s=i\,A\,s$と書けるとき、$c_1(E)=-\frac1{2\pi}dA$です。$dA$は向きを決めると関数と面積要素の積になり、積分の符号は向きによります。

## 2次元の球面

前回の記事の$S^2$で、両辺を照合します。そのためには、主張の中の標準的なカイラリティと向きを、前回の記事の$γ$と揃える必要があります。

前回の記事では、$\boldsymbol t_1\boldsymbol t_2\boldsymbol x=ω$となる向き（外向きの法線$\boldsymbol x$に対して反時計回り、以下**外向きの向き**）の接ベクトルについて、$\boldsymbol t_1\bullet\boldsymbol t_2\bulletψ=-ω\boldsymbol xψ$を示しました。[[7shi-kth7]]

虚数単位$i$は右から$J$を掛ける演算なので、順序を入れ替えた$(\boldsymbol t_2,\boldsymbol t_1)$について次のようになります。

$$
i\,\boldsymbol t_2\bullet\boldsymbol t_1\bulletψ=\boldsymbol t_2\boldsymbol x\boldsymbol t_1\boldsymbol xψJ=\boldsymbol t_1\boldsymbol t_2ψJ=ω\boldsymbol xψJ=γψ
$$

したがって、前回の記事の$γ$は、外向きと逆の向き$(\boldsymbol t_2,\boldsymbol t_1)$での標準的なカイラリティです。この向きは、南側の座標$z=α/β$を複素座標とする向きと一致します。実際、前回の記事の極座標で$z=\cot\fracθ2e^{Jφ}$なので、$θ$が増えると$|z|$は減り、$φ$が増えると$z$は反時計回りに回ります。外向きの向きの組$(\partial_θ,\partial_φ)$は、$z$の平面では内向きと反時計回りの組になり、$z$の平面の標準的な向きと逆です。[[7shi-kth7]][[7shi-homog]]

$H^n$の接続を$c_1$の式に当てはめます。前回の記事の捩った共変微分は、捩りの部分に$(\boldsymbol a\cdot X)R_J$を加えたもので、$R_J$は$i$を掛けることなので、$A$は$\boldsymbol a$を1次微分形式と見たもの$A(X)=\boldsymbol a\cdot X$です。$dA(\boldsymbol t_1,\boldsymbol t_2)=\boldsymbol t_2\cdot\partial_{\boldsymbol t_1}\boldsymbol a-\boldsymbol t_1\cdot\partial_{\boldsymbol t_2}\boldsymbol a$は前回の記事の曲率$F$の定義そのものなので、外向きの向きでは$dA=F\,dΩ$です。[[7shi-kth7]]

外向きと逆の向きで積分すると符号が反転し、次のようになります。

$$
\int_{S^2}c_1(H^n)=\frac1{2\pi}\int_{S^2}F\,dΩ=-n
$$

中辺の積分は、前回の記事と同じく外向きの$dΩ$で取ったものです。指数定理の右辺は$-n$で、前回の記事で求めた$\operatorname{ind}D_n=-n$と一致します。

&&&rem 符号の規約
向きを外向きに取ると、標準的なカイラリティは$i\,\boldsymbol t_1\bullet\boldsymbol t_2\bullet=-γ$になり、正負のカイラリティが入れ替わって指数は$+n$になります。同時に積分の向きも反転して$\int c_1(H^n)=+n$となるので、指数定理の両辺はそろって符号を変えます。前回の記事の注意「符号の規約」で述べた$+n$は、この規約での値です。本記事では、前回の記事の$γ$を保ち、$z$を複素座標とする向きで$c_1$を積分します。この向きでは$\int c_1(H)=-1$で、ホップ束$H$（$\mathbb C^2$の直線を割り当てる束）の第1チャーン数として知られている値と一致します。[[7shi-kth7]]
&&&

$z$の向きでは、貼り合わせ関数$f$の束$E_f$について$\int c_1(E_f)=-\operatorname{wind}(f)$です。以前の記事のテプリッツ作用素の指数と並べると、$\operatorname{ind}T_f=-\operatorname{wind}(f)=\int c_1(E_f)$となり、テプリッツ作用素、捩ったディラック作用素、第1チャーン数の3つが同じ符号で揃います。[[7shi-kth5]]

前回の記事では、指数を曲率が一定の接続についてだけ求め、指数が接続の選び方によらないことは扱いませんでした。指数定理の右辺は接続によらないので、指数定理を認めれば、$H^n$の捩りの条件を満たすどの接続についても指数は$-n$です。前回の記事では、曲率の積分が接続によらず回転数で決まることを示しましたが、指数定理はその積分が指数そのものであることを主張します。[[7shi-kth7]]

## 四元数のホップ束

4次元の球面$S^4$で、指数定理の2は次の形になります。$S^4$のポントリャーギン形式の積分は$0$であることが知られているので、$\hat A(S^4)$の4次の部分は積分に効かず、$\operatorname{ch}(E)$の4次の部分だけが残ります。$c_1(E)=0$となる束では、次のようになります。

$$
\operatorname{ind}D_E=\int_{S^4}\operatorname{ch}_2(E)=-\int_{S^4}c_2(E)
$$

&&&ex 四元数のホップ束で捩る [ex-quaternion]
以前の記事の$W=\mathbb H$（$\mathbf i,\mathbf j,\mathbf k$を左から掛ける作用）の束$E_W$は、四元数のホップ束である。右から$i$を掛ける演算は左からの積と可換なので、$E_W$は複素2次元の束と見なせる。この束は$c_1=0$、$\int c_2=\pm1$（符号は向きによる）を満たすことが知られている。したがって、指数定理を認めると、四元数のホップ束で捩ったディラック作用素の指数は$\pm1$である。$S^2$で$H$に捩った場合の指数$-1$と同じく、K群の生成元に絶対値$1$の指数が対応する。[[7shi-kth4]]

向きを固定した$S^4$上で、第2チャーン数が$k\ge1$の$\operatorname{SU}(2)$の束を考える。曲率が反自己双対という方程式を満たす接続（物理ではインスタントンと呼ばれる）をゲージ同値で割ったモジュライ空間の次元$8k-3$も、線形化した方程式の作用素に指数定理を当てはめて計算されることが知られている。本記事では扱わない。
&&&

以前の記事では、四元数のホップ束が自明束を足しても自明にならないことを、主張に留めました。複素ベクトル束として（$\tilde K(S^4)$の中で）のこの性質は、指数でも検出できます。[[7shi-kth4]]

&&&prop 四元数のホップ束の安定な非自明性 [prop-stable-quaternion]
指数定理と[[ex-quaternion]]の値$\int c_2(E_W)=\pm1$を認めると、どの$k\ge0$についても、$E_W\oplus\underline{\mathbb C}^k$は自明束ではない。
&&&

&&&prf
$\operatorname{ch}$は直和について和になり、自明束$\underline{\mathbb C}^N$の接続を曲率$0$に取れば$\operatorname{ch}(\underline{\mathbb C}^N)=N$である。$\int_{S^4}p_1=0$より$\hat A(S^4)$の4次の部分の積分は$0$なので、指数定理の右辺は次のようになる。

$$
\int_{S^4}\hat A(S^4)\operatorname{ch}(E_W\oplus\underline{\mathbb C}^k)=\int_{S^4}\operatorname{ch}_2(E_W)-\frac{k+2}{24}\int_{S^4}p_1=\mp1
$$

同様に、自明束$\underline{\mathbb C}^{k+2}$で捩った場合の右辺は$0$である。右辺は束の同型類だけで決まるので、$E_W\oplus\underline{\mathbb C}^k$が$\underline{\mathbb C}^{k+2}$と同型なら両者の指数は一致し、$\mp1=0$となって矛盾する。
&&&

K群の言葉では、$\tilde K(S^4)$の中で$[E_W]-2\ne0$ということです。証明で実際に使ったのは$\int\operatorname{ch}_2$が自明束を足しても変わらないことなので、第2チャーン数の積分が安定な不変量であることを認めれば、指数を経由しなくても同じ結論が出ます。指数定理が加えるのは、この位相的な量が作用素の核の次元の差として現れるという読み方です。アティヤ＝ボット＝シャピロの定理で主張に留めた部分のうち、$k=4$の複素K群の生成元が$0$でないことを、指数の側から確かめたことになります。実ベクトル束としての非自明性（$\widetilde{KO}(S^4)$）には、本記事では立ち入りません。[[7shi-kth4]]

# 特殊化

## 4つの作用素

指数定理の1は、どの楕円型作用素にも当てはまります。幾何に現れる主な作用素はいずれもディラック型の作用素として扱え、指数定理の1から次の公式が得られることが知られています。スピン多様体の上のディラック作用素を捩ったものとして2の形に書くには、追加の条件や設定が要ります。

| 作用素 | 解析的指数 | 位相的な側 | 名前 |
|:---|:---|:---|:---|
| $d+δ$（偶数次の形式から奇数次の形式へ） | オイラー標数$\chi(M)$ | $\int_Me(M)$（オイラー形式） | ガウス＝ボネの定理 |
| 符号数作用素 | 符号数$τ(M)$ | $\int_ML(M)$、$L=1+\frac{p_1}3+\cdots$ | ヒルツェブルフの符号数定理 |
| ドルボー作用素$\bar\partial_E+\bar\partial_E^*$（複素多様体、$E$は正則ベクトル束） | $\sum_q(-1)^q\dim H^q(M,E)$ | $\int_M\operatorname{Td}(M)\operatorname{ch}(E)$、$\operatorname{Td}=1+\frac{c_1}2+\cdots$ | ヒルツェブルフ＝リーマン＝ロッホの定理 |
| ディラック作用素$D$ | $\operatorname{ind}D$ | $\int_M\hat A(M)$ | $\hat A$種数の整数性 |

$d$は外微分、$δ$はその形式的な随伴（余微分）です。以前の記事では、ユークリッド空間のディラック作用素が、$e_a^2=+1$の規約で$D=d-δ$と書けることを見ました。符号数作用素は、$4l$次元の多様体の上で、$d+δ$を定義域の分け方だけ替えたものです。符号数$τ(M)$は、中間次元の形式の積分$\int α\wedgeβ$が定める2次形式の符号数です。前回の記事ではスピノルに$χ$を使いましたが、本記事では$\chi$をオイラー標数に使います。[[7shi-cla5]]

表の一般的な公式の証明は、本記事では扱いません。以下、$S^2$で確かめられるものを確かめます。

## ガウス＝ボネの定理

$d+δ$の指数がオイラー標数に等しいことは、閉多様体の上では調和形式（$d+δ$の核）の空間の次元が各次数のベッチ数に等しいこと（ホッジの定理）から従うことが知られています。$S^2$では、両辺を直接計算できます。

$S^2$上の微分形式を、ベクトル解析の言葉で書きます。0次の形式は関数$f$、1次の形式は接ベクトル場$\boldsymbol a$（$X\mapsto\boldsymbol a\cdot X$と同一視）、2次の形式は関数$g$と外向きの面積要素の積$g\,dΩ$です。$\nabla f$は球面に沿った勾配、$\operatorname{div}\boldsymbol a$は球面に沿った発散とします。外微分は次のようになります。

$$
df=\nabla f,\qquad d\boldsymbol a=F(\boldsymbol a)\,dΩ,\qquad F(\boldsymbol a)=(\nabla\times\boldsymbol a)\cdot\boldsymbol x
$$

$F(\boldsymbol a)$は前回の記事の曲率と同じ式で、$\boldsymbol a$を球面の外に延長した3次元の回転の法線成分です。勾配の回転は$0$なので、$F(\nabla f)=0$です。[[7shi-kth7]]

$δ$は$L^2$内積について$d$の随伴です。発散定理（球面には境界がないので境界の項は現れない）から$(\nabla f,\boldsymbol a)=-(f,\operatorname{div}\boldsymbol a)$、また$\nabla\times(g\boldsymbol a)=g\nabla\times\boldsymbol a+\nabla g\times\boldsymbol a$の法線成分を積分すると、ストークスの定理（境界のない閉曲面では左辺の積分が$0$）と$(\nabla g\times\boldsymbol a)\cdot\boldsymbol x=(\boldsymbol x\times\nabla g)\cdot\boldsymbol a$から$(g,F(\boldsymbol a))=(-\boldsymbol x\times\nabla g,\boldsymbol a)$となります。[[7shi-cla2]]

したがって、$δ$は次のようになります。

$$
δ\boldsymbol a=-\operatorname{div}\boldsymbol a,\qquad δ(g\,dΩ)=-\boldsymbol x\times\nabla g
$$

$\boldsymbol x\times\nabla g$は$\nabla g$を接平面の中で直角に回したもので、$|\boldsymbol x\times\nabla g|=|\nabla g|$です。以前の記事の球面のラプラシアン$Δ_S$は$\operatorname{div}\nabla$に等しく、部分積分により$(f,Δ_Sf)=-\|\nabla f\|^2$です。[[7shi-kth6]]

&&&prop 2次元の球面のオイラー標数 [prop-euler-s2]
$S^2$上の$d+δ$の、偶数次の形式から奇数次の形式への部分を$P^+$、その逆を$P^-$とする。$\ker P^+$は定数関数と$dΩ$の定数倍で張られる2次元の空間で、$\ker P^-=0$である。したがって$\operatorname{ind}(d+δ)=2$である。
&&&

&&&prf
$P^+(f,g\,dΩ)=\nabla f-\boldsymbol x\times\nabla g$、$P^-\boldsymbol a=(-\operatorname{div}\boldsymbol a,\ F(\boldsymbol a)\,dΩ)$である。

$\ker P^+$。$\nabla f=\boldsymbol x\times\nabla g$とする。上の随伴の関係と$F(\nabla f)=0$から、$\|\nabla f\|^2=(\nabla f,\boldsymbol x\times\nabla g)=-(g,F(\nabla f))=0$である。よって$\nabla f=0$であり、$|\boldsymbol x\times\nabla g|=|\nabla g|$から$\nabla g=0$でもある。球面は連結なので$f,g$は定数である。逆に定数の組は$P^+$で$0$に写る。

$\ker P^-$。$\operatorname{div}\boldsymbol a=0$、$F(\boldsymbol a)=0$とする。$S^2$上で回転が$0$の接ベクトル場は、ある関数の勾配であることが知られている（$S^2$が単連結であることによる）。$\boldsymbol a=\nabla f$とすると$Δ_Sf=\operatorname{div}\nabla f=0$なので、$\|\nabla f\|^2=-(f,Δ_Sf)=0$であり、$\boldsymbol a=0$である。
&&&

位相的な側では、$S^2$のベッチ数は$b_0=b_2=1$、$b_1=0$で$\chi(S^2)=2$です。オイラー形式は、2次元ではガウス曲率$K$を使って$\frac1{2\pi}K\,dA$で、単位球面では$K=1$、面積$4\pi$なので積分は$2$です。調和形式の数え上げ、ベッチ数、曲率の積分が、すべて$2$で一致します。

以前の記事では、$TS^2$の貼り合わせ関数から作った接ベクトル場の北極の零点が$+2$回まわることを見て、ポアンカレ＝ホップの定理の値$\chi(S^2)=2$と一致することを述べました。また、$z$を複素座標とする向きでは、複素直線束$TS^2\cong H^{-2}$の第1チャーン数は$\int c_1(H^{-2})=2$です。接束のねじれ、零点の回転数、調和形式の個数が、同じ数$2$で結びついています。[[7shi-kth3]]

## リーマン＝ロッホの定理

前回の記事では、$H^{-m}$で捩ったディラック作用素の核が、南側の座標で$P(z,1)/(1+|z|^2)^{(m-1)/2}$の形に書けることを見ました。$P(z,1)$が$z$の多項式であることは、核の方程式が正則関数の方程式であることを示唆します。実際、正のカイラリティの側で、核の方程式はコーシー＝リーマン方程式になります。[[7shi-kth7]]

南極を含む領域$U_-$の上で、正のカイラリティの捩ったスピノル$ψ$は、$\mathbb C_J$値の関数$q$を使って次のように書けます。

$$
ψ_-=Φ(α_-,β_-)\,(1+|z|^2)^{(n+1)/2}q
$$

$Φ(α_-,β_-)$は各点で$S^+_{\boldsymbol x}$を張り、係数は正の関数なので、この書き方は一意です。$z=α_-/β_-$は、前回の記事の南側の座標です。

&&&prop 核の方程式とコーシー＝リーマン方程式 [prop-dbar]
$U_-$の上で、$D_nψ=0$は$\dfrac{\partial q}{\partial\bar z}=0$と同値である。ここで$\dfrac{\partial}{\partial\bar z}=\dfrac12\left(\dfrac{\partial}{\partial(\operatorname{Re}z)}+J\dfrac{\partial}{\partial(\operatorname{Im}z)}\right)$である。
&&&

&&&prf
両極を除いた部分で示せば、連続性から$U_-$全体で成り立つ。前回の記事の枠の形を使うため、$U_+$の側に移す。$ψ_+=ψ_-e^{nJφ}$であり、前回の記事で見た$Φ(α_+,β_+)=Uf_1e^{-Jφ/2}$と$(α_+,β_+)=(α_-,β_-)e^{-Jφ}$から$Φ(α_-,β_-)=Uf_1e^{Jφ/2}$である。$1+|z|^2=\sin^{-2}\fracθ2$なので、$ψ_+=Uf_1v$、次のように置ける（$\mathbb C_J$の元は互いに可換）。

$$
v=e^{J(n+\frac12)φ}\,W\,q,\qquad W=\sin^{-(n+1)}\fracθ2
$$

前回の記事で見た枠の形は、$A=-\frac n2(1-\cosθ)$として次のとおりである。

$$
\tilde UD_nU(f_1v)=-\left(\partial_θ+\frac{\cotθ}2\right)v+\frac J{\sinθ}\bigl(\partial_φv+AvJ\bigr)
$$

$\partial_θW=-\frac{n+1}2\cot\fracθ2\,W$、$\partial_φe^{J(n+\frac12)φ}=J(n+\frac12)e^{J(n+\frac12)φ}$を代入すると、右辺は$e^{J(n+\frac12)φ}W$と次の量の積になる。

$$
-\partial_θq+\frac J{\sinθ}\partial_φq+\left(\frac{n+1}2\cot\fracθ2-\frac{\cotθ}2-\frac{n+\frac12+A}{\sinθ}\right)q
$$

括弧の中に$\sinθ$を掛けると、$\sinθ\cot\fracθ2=1+\cosθ$より$\frac{n+1}2(1+\cosθ)-\frac{\cosθ}2-n-\frac12+\frac n2(1-\cosθ)=0$である。一方、$z=re^{Jφ}$、$r=\cot\fracθ2$とすると、$\frac{\partial}{\partial\bar z}=\frac{e^{Jφ}}2\left(\partial_r+\frac Jr\partial_φ\right)$、$\partial_r=-2\sin^2\fracθ2\,\partial_θ$、$\tan\fracθ2/(2\sin^2\fracθ2)=1/\sinθ$から、次のようになる。

$$
\frac{\partial q}{\partial\bar z}=-e^{Jφ}\sin^2\fracθ2\left(\partial_θq-\frac J{\sinθ}\partial_φq\right)
$$

したがって$\tilde UD_nU(f_1v)=e^{J(n-\frac12)φ}W\sin^{-2}\fracθ2\,\dfrac{\partial q}{\partial\bar z}$であり、係数は$0$にならないので主張を得る。[[7shi-kth7]]
&&&

したがって、$\ker D_n^+$は、$z$の平面全体で正則な関数$q$のうち、$ψ$が北極（$z=\infty$）まで滑らかに延びるものに対応します。$q=z^j$なら$|ψ|=(1+|z|^2)^{(n+1)/2}|z|^j$で、$z\to\infty$で有界なのは$n\le-1$かつ$j\le-n-1$のときだけです。前回の記事で求めた核の元$Φ(α,β)P(α,β)$は、$q=P(z,1)$にあたります。

この結果を、正則な直線束の言葉で読みます。前回の記事で見たとおり、正のカイラリティで$H^n$に捩ったスピノルの束は$L=H^{n+1}$です。[[prop-dbar]]は、$D_n^+$が、$z$の正則関数を正則な切断とする複素構造での$L$の$\bar\partial$作用素と、0でない因子を除いて一致することを示しています。$\ker D_n^-$は$\bar\partial$の余核にあたります。$S^2$でのリーマン＝ロッホの定理は、次数$\deg L=\int c_1(L)$の正則な直線束について、$\bar\partial$の指数が$\deg L+1$であることを主張します。[[7shi-kth7]]

$z$の向きでは$\deg H^{n+1}=-(n+1)$なので、$\bar\partial$の指数は$-(n+1)+1=-n$となり、前回の記事の$\operatorname{ind}D_n=-n$と一致します。表の公式では、$\int\operatorname{Td}(S^2)\operatorname{ch}(L)=\int c_1(L)+\frac12\int c_1(TS^2)$で、$TS^2\cong H^{-2}$から$\frac12\int c_1(TS^2)=1$です。$\int\hat A(S^2)\operatorname{ch}(H^n)=\int c_1(H^n)=-n$と、$\bar\partial$として数えた$\int c_1(H^{n+1})+1=-n$は、$S^+\cong H$と$\int c_1(H)=-\frac12\int c_1(TS^2)$によって同じ値になります。[[7shi-kth3]]

$n\le-1$では、ここでは$m=-n$と置くと、正則な切断は次数$m-1$以下の多項式で、その個数は$m=\deg L+1$です。リーマン＝ロッホの定理だけから分かるのは核と余核の次元の差ですが、前回の記事の片側の核の消失（$n\le0$で$\ker D_n^-=0$）と合わせると、指数$\deg L+1$がそのまま核の次元になります。前回の記事で変数分離により求めた核の次元は、こうして次数と種数（$S^2$では$0$）から得られます。

## $\hat A$種数

以前の記事では、球面上のディラック作用素について、リヒネロビッチの公式$D_S^2=\nabla^*\nabla+\frac R4$とスカラー曲率$R>0$から、核が$0$になることを示しました。一般の閉スピン多様体でも、同じ形の公式$D^2=\nabla^*\nabla+\frac R4$が成り立つことが知られています。これを指数定理と組み合わせます。[[7shi-kth6]]

&&&prop 正のスカラー曲率と$\hat A$種数 [prop-ahat]
一般の閉スピン多様体でのリヒネロビッチの公式と指数定理を認める。$M$が偶数次元の閉スピン多様体で、スカラー曲率が至るところ正の計量を持つなら、$\int_M\hat A(M)=0$である。
&&&

&&&prf
以前の記事の命題「核の消失」の証明と同じく、$D$が対称であることと$(ψ,\nabla^*\nabla ψ)=\|\nabla ψ\|^2\ge0$から$\|Dψ\|^2\ge\frac14\min R\,\|ψ\|^2$であり、$Dψ=0$なら$ψ=0$である。よって$\ker D^+=\ker D^-=0$で$\operatorname{ind}D=0$である。指数定理の2を、自明な直線束（曲率$0$の接続、$\operatorname{ch}=1$）で捩った場合に使えば$\int\hat A(M)=0$を得る。[[7shi-kth6]]
&&&

球面ではスカラー曲率が正なので$\int\hat A(S^{2m})=0$で、これは球面のポントリャーギン形式の積分が$0$であることとも整合します。この命題の価値は対偶にあります。

&&&ex K3曲面
K3曲面と呼ばれる4次元の閉スピン多様体は、符号数$τ=-16$を持つことが知られている。4次元では、符号数定理から$τ=\frac13\int p_1$なので$\int p_1=-48$で、$\int\hat A=-\frac1{24}\int p_1=2$である。$0$でないので、[[prop-ahat]]により、K3曲面にはスカラー曲率が至るところ正の計量が入らない。
&&&

指数は整数なので、閉スピン多様体では$\int\hat A(M)$は整数でなければなりません。これも位相への制約になります。

&&&ex 複素射影平面
複素射影平面$\mathbb CP^2$は$\int p_1=3$を持ち、符号数定理から$τ=1$である。$\int\hat A=-\frac18$は整数でないので、$\mathbb CP^2$はスピン構造を持たない。4次元の閉スピン多様体では、同じ計算から$\int\hat A=-\frac τ8$が整数なので、符号数は$8$で割り切れる。さらに、4次元のスピノルの四元数構造から指数が偶数になることも使うと、$16$で割り切れること（ロホリンの定理）が導けることが知られている。
&&&

# 毛玉の定理

## 微分形式と接ベクトルの掛け算

以前の記事では、偶数次元の球面$S^{2m}$にどこでも$0$にならない連続な接ベクトル場がないこと（毛玉の定理）を主張に留め、$S^2$の場合だけを貼り合わせ関数の回転数から示しました。$m\ge2$の場合を、指数定理から導きます。鍵は、消えない接ベクトル場が、主表象にもう1つの反交換する生成元を与えることです。[[7shi-kth2]][[7shi-kth3]]

ホッジの定理から、$d+δ$の指数は$\chi(S^{2m})=2$です（$S^{2m}$のベッチ数は$b_0=b_{2m}=1$、ほかは$0$であることが知られています）。$d+δ$の主表象を、以前の記事の接ベクトルの掛け算で書きます。

以前の記事では、$S^{n-1}$上の接ベクトル$X$を$X\boldsymbol x$に対応させると、接空間のクリフォード代数が偶部分代数$\operatorname{Cl}_{n,0}^0(\mathbb R)$と同一視できることを見ました。$n=2m+1$とし、点$\boldsymbol x$で正規直交な接ベクトル$\boldsymbol t_1,\dots,\boldsymbol t_{2m}$を取ります。$p$次の形式$\boldsymbol t_{i_1}\wedge\cdots\wedge\boldsymbol t_{i_p}$（$i_1<\dots<i_p$）を、次の偶部分代数の元に対応させます。

$$
(\boldsymbol t_{i_1}\boldsymbol x)(\boldsymbol t_{i_2}\boldsymbol x)\cdots(\boldsymbol t_{i_p}\boldsymbol x)
$$

どちらの側も$2^{2m}$次元で、この対応は各点で線形同型です。以下、形式の係数を複素数に広げます。[[7shi-kth6]]

&&&lem 形式の上の演算 [lem-forms]
上の同一視のもとで、次が成り立つ。

1. 接ベクトル$\boldsymbol\xi$の掛け算$ψ\mapsto\boldsymbol\xi\boldsymbol xψ$は、外積から内部積を引いた$\boldsymbol\xi\wedge-ι_{\boldsymbol\xi}$である。
2. $εψ=\boldsymbol xψ\boldsymbol x$は、$p$次の形式に$(-1)^p$を掛ける演算である。
&&&

&&&prf
1. 線形性から$\boldsymbol\xi=\boldsymbol t_j$の場合を示せばよい。$\boldsymbol t_i\boldsymbol x$どうしは反交換し、$(\boldsymbol t_j\boldsymbol x)^2=-1$である。$j$が$i_1,\dots,i_p$に含まれないとき、$\boldsymbol t_j\boldsymbol x$を正しい位置まで動かすと、飛び越えた個数の符号が付き、これは$\boldsymbol t_j\wedge$の符号と一致する。$j=i_q$のとき、$\boldsymbol t_j\boldsymbol x$を$q$番目まで動かすと符号$(-1)^{q-1}$が付き、$(\boldsymbol t_j\boldsymbol x)^2=-1$で消えるので、全体で$-(-1)^{q-1}$倍の$\boldsymbol t_{i_q}$を除いた形式になる。これは$-ι_{\boldsymbol t_j}$である。

2. $\boldsymbol x$は$\boldsymbol t_i$と反交換するので$\boldsymbol x(\boldsymbol t_i\boldsymbol x)=-(\boldsymbol t_i\boldsymbol x)\boldsymbol x$であり、$\boldsymbol x^2=1$から$\boldsymbol xψ\boldsymbol x$は$p$個の因子を通過して$(-1)^pψ$になる。
&&&

$d$の主表象は$\boldsymbol\xi\wedge$、$δ$の主表象は$-ι_{\boldsymbol\xi}$であることが知られています（$d$は微分に外積を、$δ$は微分に内部積を組み合わせた作用素）。したがって[[lem-forms]]の1から、$d+δ$の主表象は$\boldsymbol\xi\bullet$そのもので、ディラック作用素と同じです。$ε$は偶数次と奇数次の形式を分ける演算で、$d+δ$の指数を定める分け方です。

## ベクトル場による変形

$\boldsymbol v$を$S^{2m}$上のどこでも$0$にならない接ベクトル場とし、形式に次の演算を施します。

$$
B_{\boldsymbol v}ψ=i\,\boldsymbol xψ\boldsymbol v
$$

$ψ$が偶部分代数に値を取るので、$\boldsymbol xψ\boldsymbol v$も偶部分代数に値を取ります。$B_{\boldsymbol v}=i\,εR_{\boldsymbol v\boldsymbol x}$（$R_{\boldsymbol v\boldsymbol x}$は右から$\boldsymbol v\boldsymbol x$を掛ける演算）であり、左からの積$\boldsymbol\xi\bullet$とは別の側から作用します。

&&&prop ベクトル場による主表象の変形 [prop-deform]
接ベクトル$\boldsymbol\xi$と$0\le t\le1$について、次が成り立つ。

1. $B_{\boldsymbol v}$は$ε$と反交換し、偶数次の形式と奇数次の形式を入れ替える。
2. $B_{\boldsymbol v}$は$\boldsymbol\xi\bullet$と反交換し、$B_{\boldsymbol v}^2=-|\boldsymbol v|^2$である。
3. $(\boldsymbol\xi\bullet+tB_{\boldsymbol v})^2=-(|\boldsymbol\xi|^2+t^2|\boldsymbol v|^2)$である。とくに$σ_t(\boldsymbol\xi)=\boldsymbol\xi\bullet+tB_{\boldsymbol v}$は、$\boldsymbol\xi\ne0$ならどの$t$でも可逆で、$t=1$では$\boldsymbol\xi=0$でも可逆である。
&&&

&&&prf
1. $\boldsymbol x^2=1$より$εB_{\boldsymbol v}ψ=i\,\boldsymbol x\boldsymbol xψ\boldsymbol v\boldsymbol x=iψ\boldsymbol v\boldsymbol x$、$B_{\boldsymbol v}εψ=i\,\boldsymbol x\boldsymbol xψ\boldsymbol x\boldsymbol v=iψ\boldsymbol x\boldsymbol v$であり、$\boldsymbol x\boldsymbol v=-\boldsymbol v\boldsymbol x$から両者は符号が逆である。

2. $B_{\boldsymbol v}(\boldsymbol\xi\boldsymbol xψ)=i\,\boldsymbol x\boldsymbol\xi\boldsymbol xψ\boldsymbol v=-i\,\boldsymbol\xiψ\boldsymbol v$、$\boldsymbol\xi\boldsymbol xB_{\boldsymbol v}ψ=i\,\boldsymbol\xi\boldsymbol x\boldsymbol xψ\boldsymbol v=i\,\boldsymbol\xiψ\boldsymbol v$なので反交換する。$B_{\boldsymbol v}^2ψ=i^2\boldsymbol x\boldsymbol xψ\boldsymbol v\boldsymbol v=-|\boldsymbol v|^2ψ$である。

3. 2と$(\boldsymbol\xi\bullet)^2=-|\boldsymbol\xi|^2$から、交差項が消えて主張の式を得る。右辺が$0$でなければ、$-(|\boldsymbol\xi|^2+t^2|\boldsymbol v|^2)^{-1}σ_t(\boldsymbol\xi)$が逆写像である。
&&&

$B_{\boldsymbol v}$は零階の項なので、$d+δ+tB_{\boldsymbol v}$は$d+δ$と同じ主表象を持ちます。$σ_t$は$d+δ+tB_{\boldsymbol v}$の主表象ではなく（$t>0$では$\boldsymbol\xi$について1次式でない）、主表象の定めるK群の元を調べるための写像の連続変形です。[[prop-deform]]の意味は、主表象$\boldsymbol\xi\bullet$を、$\boldsymbol\xi\ne0$での可逆性を保ったまま、$\boldsymbol\xi=0$まで可逆な$σ_1$に変形できることです。

これは、以前の記事の命題「延長できる加群は自明な束を与える」と同じ形をしています。そこでは、加群の上に1つ多い生成元$J_k$があれば、$\cos t\,g_W+\sin t\,J_k$によって赤道全体を1点に寄せられ、貼り合わせ関数から作った束は自明になりました。ここでは、各点で$\boldsymbol v/|\boldsymbol v|$から作った$B$が、2乗が$-1$で$\boldsymbol\xi\bullet$と反交換する1つ多い生成元にあたります。単位ベクトル$\boldsymbol\xi$について$\cos s\,\boldsymbol\xi\bullet+\sin s\,B_{\boldsymbol v/|\boldsymbol v|}$は2乗が$-1$で、$s=0$の主表象から$s=\frac\pi2$の$\boldsymbol\xi$によらない写像に変形されます。以前の記事では、クリフォード代数$\operatorname{Cl}_{0,k}(\mathbb R)$の加群から$k$本の正規直交な接ベクトル場を構成しましたが、ここでは逆に、接ベクトル場が生成元として現れています。[[7shi-kth4]][[7shi-kth2]]

&&&thm 毛玉の定理 [thm-hairy-ball]
$m\ge1$とする。$S^{2m}$上には、どこでも$0$にならない連続な接ベクトル場は存在しない。
&&&

&&&prf
次の主張を認める。$d+δ$の指数は$\chi(S^{2m})=2$である（ホッジの定理と球面のベッチ数）。指数定理の1により、解析的指数は位相的指数に等しい。位相的指数については、「位相的指数」の節で述べた3つの性質（主表象の定めるK群の元で決まること、コンパクト集合の外での可逆性を保つ連続変形で変わらないこと、どこでも可逆な写像では$0$になること）を認める。

消えない連続な接ベクトル場があれば、少し動かして消えない滑らかな接ベクトル場$\boldsymbol v$が取れる。[[prop-deform]]により、$d+δ$の主表象$\boldsymbol\xi\bullet$は$\boldsymbol\xi\ne0$で可逆なまま$σ_1$に変形できる。$\boldsymbol\xi\ne0$では常に可逆なので、この変形はコンパクト集合（零切断）の外での可逆性を保ち、K群の元は変わらない。$σ_1$は$\boldsymbol\xi=0$でも可逆なので、その元は$0$である。したがって$d+δ$の位相的指数は$0$であり、解析的指数$2$と矛盾する。
&&&

$B_{\boldsymbol v}$は形式の上では$B_{\boldsymbol v}=-i(\boldsymbol v^\flat\wedge+ι_{\boldsymbol v})$（$\boldsymbol v^\flat$は$\boldsymbol v$を1次形式と見たもの）と書け、球面の位置ベクトル$\boldsymbol x$を使わずに、一般のリーマン多様体の上で同じ反交換関係を満たします。証明に使ったのはこの構成と$\chi(S^{2m})\ne0$だけなので、同じ議論から、オイラー標数が$0$でない閉多様体には消えない接ベクトル場がないことが従います。奇数次元の球面では$\chi=0$で、以前の記事で見たとおり消えない接ベクトル場が実際に構成できます。$m=1$の場合は、以前の記事で回転数から直接示したものと同じ結論です。[[7shi-kth2]][[7shi-kth3]]

# まとめ

アティヤ＝シンガーの指数定理を主張として述べ、これまで球面で計算してきた結果がその特殊な場合にあたることを確かめました。ディラック作用素の主表象$\boldsymbol\xi\bullet$は、正のカイラリティから負のカイラリティへの写像と見ると、各点で正のカイラリティの空間を加群とする貼り合わせ関数そのもので、既約なスピノルでは球面のK群の生成元を与えます。$S^2$では、前回の記事の$γ$が$z$を複素座標とする向きでの標準的なカイラリティであることを確かめ、指数$-n$が$H^n$の第1チャーン数に一致することを見ました。指数定理を認めると、四元数のホップ束の複素ベクトル束としての安定な非自明性が指数で検出でき、正のスカラー曲率と$\hat A$種数の関係、偶数次元の球面の毛玉の定理が導けます。$S^2$では、オイラー標数$2$を調和形式から数え、捩ったディラック作用素の核の方程式がコーシー＝リーマン方程式になることを示して、前回の記事の結果をリーマン＝ロッホの定理として読みました。

&&& 主表象と貼り合わせ関数
主表象は各点で次の形になり、$S^+$を$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群とする貼り合わせ関数に一致する。
$$
σ(\boldsymbol\xi)=\boldsymbol\xi\bullet,\qquad e_1^{-1}\boldsymbol\xi\big|_{S^+}=g_{S^+}(\boldsymbol\xi)
$$
&&&

&&& 指数定理
偶数次元の閉スピン多様体$M$と、エルミート計量を保つ接続を持つ複素ベクトル束$E$について、次が成り立つ。
$$
\operatorname{ind}D_E=\int_M\hat A(M)\operatorname{ch}(E)
$$
&&&

&&& 2次元の球面
$c_1$を$z$を複素座標とする向きで積分し、中辺の$dΩ$は外向きとすると、次が成り立つ。
$$
\operatorname{ind}D_n=\int_{S^2}c_1(H^n)=\frac1{2\pi}\int_{S^2}F\,dΩ=-n
$$
正のカイラリティでは、核の方程式は$\partial q/\partial\bar z=0$で、指数はリーマン＝ロッホの定理の$\deg H^{n+1}+1=-n$である。
&&&

&&& 特殊化
$S^2$の調和形式は定数関数と面積要素で、$\operatorname{ind}(d+δ)=\chi(S^2)=2$である。スカラー曲率が正の閉スピン多様体では$\int\hat A(M)=0$である。
&&&

&&& 毛玉の定理
消えない接ベクトル場$\boldsymbol v$は、主表象と反交換する生成元$B_{\boldsymbol v}ψ=i\,\boldsymbol xψ\boldsymbol v$を与える。
$$
(\boldsymbol\xi\bullet+tB_{\boldsymbol v})^2=-(|\boldsymbol\xi|^2+t^2|\boldsymbol v|^2)
$$
$\chi(S^{2m})=2\ne0$なので、$S^{2m}$にはそのような$\boldsymbol v$がない。
&&&

K理論と指数定理を扱った一連の記事では、ベクトルをクリフォード積で掛ける写像$\boldsymbol x\mapsto\boldsymbol x\cdot$が、次の4つの場面に同じ形で現れました。

| 場面 | 写像 | 得られるもの |
|:---|:---|:---|
| 球面上のベクトル場 | $\boldsymbol x\mapsto e_i\boldsymbol x$ | $S^{n-1}$上の正規直交な接ベクトル場（ラドン＝フルヴィッツ数） |
| ベクトル束の貼り合わせ | $\boldsymbol x\mapsto x_0+\sum_lx_le_l$ | 球面のK群の生成元（複素数・四元数・八元数のホップ束） |
| ディラック作用素の主表象 | $\boldsymbol\xi\mapsto\boldsymbol\xi\bullet$ | 各点で、正のカイラリティを加群とする貼り合わせ関数（既約なスピノルでは生成元） |
| 指数 | 主表象から位相的指数へ | $\operatorname{ind}D_E=\int\hat A(M)\operatorname{ch}(E)$（$S^2$では$-n$） |

毛玉の定理の証明では、1つ目の場面の接ベクトル場が、3つ目の場面の主表象に加わる生成元として戻ってきます。
