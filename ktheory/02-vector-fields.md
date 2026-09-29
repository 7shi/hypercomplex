球面上に互いに直交する接ベクトル場が何本取れるかを、クリフォード代数の加群と分類表から求め、数体系の次元が$1,2,4,8$に限られることを導きます。

# 概要

2次元の球面$S^2$には、どこでも$0$にならない連続な接ベクトル場がありません（毛玉の定理）。一方、円周$S^1$には$\boldsymbol x=(x_1,x_2)\mapsto(-x_2,x_1)$という消えない接ベクトル場があり、3次元の球面$S^3$には単位四元数$q$に$\mathbf iq,\mathbf jq,\mathbf kq$を対応させる3本の接ベクトル場があって、各点で接空間の正規直交基底をなします。球面の次元によって、取れる接ベクトル場の本数は大きく異なります。本記事では、次の問いを扱います。

> $S^{n-1}\subset\mathbb R^n$の上に、各点で線形独立な連続な接ベクトル場は何本取れるか。

前回の記事では、ベクトル束を安定同値で見ると、$S^2$の接束$TS^2$は自明束と区別できなくなることを見ました。$TS^2\oplus\underline{\mathbb R}\cong\underline{\mathbb R}^3$だからです。同じ議論はどの次元の球面でも成り立ち、球面の接束はつねに自明束と安定同値です。接ベクトル場の本数は、安定同値では捨てられてしまう接束の違いを測る量です。[[7shi-kth1]]

この問いのうち、1次式で与えられる正規直交な接ベクトル場の最大本数は、クリフォード代数の分類表から読めます。アダムスの定理によれば、これは連続な接ベクトル場一般の最大本数でもあります。$S^3$の例で接ベクトル場を与えた$\mathbf i,\mathbf j,\mathbf k$は、2乗して$-1$になり互いに反交換する元で、四元数の左からの積はクリフォード代数$\operatorname{Cl}_{0,3}(\mathbb R)$の作用になっています。一般に、$\mathbb R^n$にクリフォード代数$\operatorname{Cl}_{0,k}(\mathbb R)$が作用すれば、基底を適切に選ぶことで、生成元を左から掛ける写像が$S^{n-1}$上に$k$本の正規直交な接ベクトル場を与えます。$\mathbb R^n$にこの作用が入るかどうかは、分類表から読める既約加群の次元で決まります。[[7shi-clif1]]

本記事は次の順に進みます。

1. 接ベクトル場と平行化可能性を定義し、数の積から作られる例を見る
2. クリフォード代数の加群から正規直交な接ベクトル場が得られること、逆に1次式で与えられる正規直交な接ベクトル場がクリフォード代数の作用を定めることを示す
3. 分類表から既約加群の次元を読み、接ベクトル場の本数（ラドン＝フルヴィッツ数）を求める
4. すべての接ベクトル場が取れる（平行化可能な）のは$S^0,S^1,S^3,S^7$の場合に限られることを見る。この議論から、ノルムの乗法性を持つ数体系の次元が$1,2,4,8$に限られること（フルヴィッツの定理の次元の部分）が従う

前提は、線形代数（直交行列・内積）とクリフォード代数の分類表です。位相については連続性しか使いません。構成で得られる本数が連続な接ベクトル場の本数の最大であること（アダムスの定理）と、毛玉の定理は主張として述べ、証明は本記事では扱いません。

フルヴィッツの定理は、以前の記事で証明なしに述べました。本記事で示すのは、その次元の部分です。[[7shi-lie7]][[7shi-sedenion]]

# 球面の接ベクトル場

## 接ベクトル場と平行化可能性

球面$S^{n-1}=\{\boldsymbol x\in\mathbb R^n\mid|\boldsymbol x|=1\}$の点$\boldsymbol x$での接空間は、$\boldsymbol x$に直交するベクトル全体$\{\boldsymbol v\in\mathbb R^n\mid\boldsymbol v\perp\boldsymbol x\}$で、$n-1$次元です。各点に接ベクトルを連続に選んだものが接ベクトル場です。本数を数えるには、各点で「別々の方向を向いている」ことを求める必要があります。

&&&def 接ベクトル場
連続な写像$\boldsymbol v:S^{n-1}\to\mathbb R^n$で、各点で$\boldsymbol v(\boldsymbol x)\perp\boldsymbol x$を満たすものを、$S^{n-1}$上の**接ベクトル場**と呼びます。$k$本の接ベクトル場$\boldsymbol v_1,\dots,\boldsymbol v_k$が、各点$\boldsymbol x$で$\boldsymbol v_1(\boldsymbol x),\dots,\boldsymbol v_k(\boldsymbol x)$が線形独立になるとき、これらを**線形独立な接ベクトル場**と呼びます。
&&&

線形独立な接ベクトル場には、各点でグラム＝シュミットの直交化を施せます。直交化の手続きは内積と割り算だけでできていて、線形独立なら分母が$0$にならないので、結果は$\boldsymbol x$について連続です。したがって、線形独立な$k$本が取れることと、各点で正規直交な$k$本が取れることは同じです。以下では、正規直交な接ベクトル場を扱います。

接空間は$n-1$次元なので、本数は$n-1$を超えません。上限の$n-1$本が取れる場合が特別です。

&&&def 平行化可能
$S^{n-1}$上に$n-1$本の線形独立な接ベクトル場が取れるとき、$S^{n-1}$は**平行化可能**であると言います。
&&&

$n-1$本の線形独立な接ベクトル場$\boldsymbol v_1,\dots,\boldsymbol v_{n-1}$があれば、$(\boldsymbol x,c_1,\dots,c_{n-1})\mapsto c_1\boldsymbol v_1(\boldsymbol x)+\dots+c_{n-1}\boldsymbol v_{n-1}(\boldsymbol x)$が各点で$\mathbb R^{n-1}$から接空間への線形同型を与え、接束$TS^{n-1}$は自明束$\underline{\mathbb R}^{n-1}$と同型になります。逆に、自明束との同型があれば、$\underline{\mathbb R}^{n-1}$の定数の切断$(1,0,\dots,0),\dots,(0,\dots,0,1)$の像が線形独立な接ベクトル場になります。平行化可能とは、接束が自明であることです。

一方、各点で$\mathbb R^n$は接空間と法線$\mathbb R\boldsymbol x$の直和に分かれるので、$S^2$の場合と同じ議論により、どの$n$でも

$$
TS^{n-1}\oplus\underline{\mathbb R}\cong\underline{\mathbb R}^n
$$

が成り立ちます。接束は、自明束を1つ足せばつねに自明になります。接ベクトル場の本数は、この足し算で消えてしまう違いを数えています。

## 数の積が与える例

低い次元では、数の積から接ベクトル場が得られます。

&&&ex 複素数・四元数・八元数
- **$S^1$**：$\mathbb R^2=\mathbb C$と見て、$z\mapsto iz$と置きます。$z=x_1+x_2i$なら$iz=-x_2+x_1i$で、$(x_1,x_2)\cdot(-x_2,x_1)=0$です。$|iz|=|z|=1$なので、1本の正規直交な接ベクトル場です。
- **$S^3$**：$\mathbb R^4=\mathbb H$と見て、$q\mapsto\mathbf iq,\ \mathbf jq,\ \mathbf kq$と置きます。3本の接ベクトル場が各点で正規直交で、$S^3$は平行化可能です。
- **$S^7$**：$\mathbb R^8=\mathbb O$と見て、7つの虚数単位$e_1,\dots,e_7$を左から掛ける$x\mapsto e_1x,\dots,e_7x$と置きます。7本の接ベクトル場が各点で正規直交で、$S^7$は平行化可能です。

直交性は、いずれもノルムの乗法性$|uq|=|u||q|$から従います。この計算は後の[[thm-hurwitz]]の証明と同じなので、ここでは結果だけを示します。
&&&

$S^0=\{\pm1\}$は0次元で、接空間は$\{0\}$です。$n-1=0$本の接ベクトル場で足りるので、$S^0$も平行化可能と見なします。

$n$が偶数なら、$\mathbb R^n=\mathbb C^{n/2}$と見て各成分に$i$を掛ける$\boldsymbol z\mapsto i\boldsymbol z$が、$S^{n-1}$上の消えない接ベクトル場を与えます。$n$が奇数の場合はそうはいきません。

&&&rem 毛玉の定理
偶数次元の球面$S^{2m}$（$n=2m+1$）には、どこでも$0$にならない連続な接ベクトル場が1本もないことが知られています。これを**毛玉の定理**と呼びます。証明は本記事では扱いません。$m=1$の場合は、前回の記事で$TS^2$が自明でないことの根拠として使いました。[[7shi-kth1]]
&&&

# クリフォード加群が与えるベクトル場

## 左からの積

数の例で接ベクトル場を与えた写像$\boldsymbol x\mapsto u\boldsymbol x$（$u$は虚数単位）は、$\mathbb R^n$の線形変換です。これらの線形変換には、共通する2つの性質があります。

- **長さを保つ**：$|u\boldsymbol x|=|\boldsymbol x|$で、接ベクトル場の各ベクトルが単位ベクトルになります。
- **2乗すると$-1$**：$u(u\boldsymbol x)=-\boldsymbol x$です。八元数でも、2つの元だけからなる積には結合法則が成り立つので（アルティンの定理）、$u(u\boldsymbol x)=(u^2)\boldsymbol x$です。

さらに、異なる虚数単位$u,v$を掛ける線形変換は反交換します（八元数を含めて、ノルムの乗法性から従います。後の[[thm-hurwitz]]の証明と同じ計算です）。この3つの条件を満たす線形変換の組は、クリフォード代数の直交的な作用を定めます。

代数の元が線形変換として作用するベクトル空間を、その代数の**加群**と呼びます。$\mathbb R^n$が$\operatorname{Cl}_{0,k}(\mathbb R)$の加群であるとは、$n$次正方行列$J_1,\dots,J_k$が

$$
J_i^2=-I,\qquad J_iJ_j=-J_jJ_i\quad(i\ne j)
$$

を満たし、生成元$e_i$が$J_i$として作用することです。$J_i$が直交行列であるとき、この作用は**直交的**であると言います。

&&&lem 直交で2乗が$-1$の行列 [lem-antisym]
$J$を直交行列とすると、$J^2=-I$と$J^T=-J$は同値です。このとき、すべての$\boldsymbol x$について$\langle J\boldsymbol x,\boldsymbol x\rangle=0$です。
&&&

&&&prf
直交行列なので$J^T=J^{-1}$である。$J^2=-I$は$J^{-1}=-J$と同値だから、$J^T=-J$と同値である。このとき$\langle J\boldsymbol x,\boldsymbol x\rangle=\langle\boldsymbol x,J^T\boldsymbol x\rangle=-\langle\boldsymbol x,J\boldsymbol x\rangle$より、$\langle J\boldsymbol x,\boldsymbol x\rangle=0$となる。
&&&

&&&thm クリフォード加群が与える接ベクトル場 [thm-module-fields]
$\mathbb R^n$が$\operatorname{Cl}_{0,k}(\mathbb R)$の直交的な加群で、生成元が$J_1,\dots,J_k$として作用するなら、

$$
\boldsymbol x\mapsto J_1\boldsymbol x,\ \dots,\ \boldsymbol x\mapsto J_k\boldsymbol x
$$

は$S^{n-1}$上の$k$本の正規直交な接ベクトル場です。
&&&

&&&prf
[[lem-antisym]]より$\langle J_i\boldsymbol x,\boldsymbol x\rangle=0$なので、$J_i\boldsymbol x$は接ベクトルである。内積は$J_i^T=-J_i$を使って

$$
\langle J_i\boldsymbol x,J_j\boldsymbol x\rangle=\langle\boldsymbol x,J_i^TJ_j\boldsymbol x\rangle=-\langle\boldsymbol x,J_iJ_j\boldsymbol x\rangle
$$

となる。$i=j$なら$J_i^2=-I$より$|\boldsymbol x|^2=1$である。$i\ne j$なら$(J_iJ_j)^T=J_j^TJ_i^T=J_jJ_i=-J_iJ_j$より$J_iJ_j$は反対称行列で、[[lem-antisym]]の証明と同じく$\langle\boldsymbol x,J_iJ_j\boldsymbol x\rangle=0$となる。
&&&

接ベクトル場であることには2乗して$-1$になることが、互いに直交することには反交換することが対応しています。

&&&rem 直交性の仮定
直交的でない加群でも、内積を取り替えれば作用を直交的にできます。$\pm J_{i_1}\cdots J_{i_r}$（$i_1<\dots<i_r$）の全体は有限個の行列からなる群をなすので、この群の元$g$全体について$\langle g\boldsymbol x,g\boldsymbol y\rangle$を平均したものを新しい内積にすれば、各$J_i$はその内積について直交変換になります。この内積について正規直交基底を選び直せば、各$J_i$は標準の内積について直交行列で表されるので、標準の単位球面$S^{n-1}$の上で接ベクトル場が得られます。以下では直交的な加群だけを考えます。
&&&

## 1次式の接ベクトル場からクリフォード関係へ

[[thm-module-fields]]の接ベクトル場は、$\boldsymbol x$の1次式です。逆に、1次式で与えられる正規直交な接ベクトル場は、必ずクリフォード代数の作用から来ます。

&&&thm 1次式の接ベクトル場 [thm-linear-fields]
$n$次正方行列$A_1,\dots,A_k$について、$\boldsymbol x\mapsto A_1\boldsymbol x,\dots,A_k\boldsymbol x$が$S^{n-1}$上の$k$本の正規直交な接ベクトル場なら、$A_1,\dots,A_k$は直交行列で

$$
A_i^2=-I,\qquad A_iA_j=-A_jA_i\quad(i\ne j)
$$

を満たします。すなわち、$\mathbb R^n$は$A_i$を生成元の作用とする$\operatorname{Cl}_{0,k}(\mathbb R)$の直交的な加群です。
&&&

&&&prf
仮定は、$|\boldsymbol x|=1$のとき$\langle A_i\boldsymbol x,\boldsymbol x\rangle=0$、$\langle A_i\boldsymbol x,A_j\boldsymbol x\rangle=\delta_{ij}$となることである。両辺とも$\boldsymbol x$について2次式なので、$\boldsymbol x$を定数倍して、すべての$\boldsymbol x\in\mathbb R^n$について

$$
\langle A_i\boldsymbol x,\boldsymbol x\rangle=0,\qquad\langle A_i\boldsymbol x,A_j\boldsymbol x\rangle=\delta_{ij}|\boldsymbol x|^2
$$

が成り立つ。行列$B$について$\langle B\boldsymbol x,\boldsymbol x\rangle=\frac12\boldsymbol x^T(B+B^T)\boldsymbol x$であり、対称行列$S$の2次形式$\boldsymbol x^TS\boldsymbol x$がすべての$\boldsymbol x$で$0$なら、$\boldsymbol x$を$\boldsymbol x+\boldsymbol y$に置き換えて展開すれば$\boldsymbol x^TS\boldsymbol y=0$となり$S=0$である（偏極）。これを使うと、1つ目の式から$A_i+A_i^T=0$、2つ目の式から

$$
A_i^TA_i=I,\qquad A_i^TA_j+A_j^TA_i=0\quad(i\ne j)
$$

が得られる。よって$A_i$は直交行列で、$A_i^T=-A_i$と[[lem-antisym]]から$A_i^2=-I$である。2つ目の式に$A_i^T=-A_i$を代入すれば$A_iA_j+A_jA_i=0$となる。
&&&

2つの定理を合わせると、$S^{n-1}$上に1次式で与えられる$k$本の正規直交な接ベクトル場が取れることと、$\mathbb R^n$が$\operatorname{Cl}_{0,k}(\mathbb R)$の直交的な加群になることは同じです。1次式で与えられる正規直交な接ベクトル場の本数を求める問題は、$\mathbb R^n$にどれだけ大きなクリフォード代数が作用するかという代数の問題になりました。

# 既約加群の次元

## 分類表から読む

$\mathbb R^n$が$\operatorname{Cl}_{0,k}(\mathbb R)$の加群になるかどうかは、$n$で決まります。作用で保たれる自明でない部分空間を持たない加群を**既約加群**と呼びます。$\operatorname{Cl}_{0,k}(\mathbb R)$の有限次元の加群は既約加群の直和に分解でき、既約加群は分類表の行列環から読めるからです。

分類表を導いた記事の「ピノルとスピノル」で見たように、行列環$M_m(\mathbb F)$の既約加群は列ベクトルの空間$\mathbb F^m$で、実次元は$m\dim\mathbb F$です。直和型$2M_m(\mathbb F)$には既約加群が2つありますが、どちらも$\mathbb F^m$で次元は同じです。[[7shi-clif1]]

$\operatorname{Cl}_{0,k}(\mathbb R)$の既約加群の実次元を$a_k$と書き、分類表の第1行から読むと次のとおりです。

$$
\begin{array}{c|ccccccccc}
k&0&1&2&3&4&5&6&7&8\\\hline
\operatorname{Cl}_{0,k}(\mathbb R)&\mathbb R&\mathbb C&\mathbb H&2\mathbb H&M_2(\mathbb H)&M_4(\mathbb C)&M_8(\mathbb R)&2M_8(\mathbb R)&M_{16}(\mathbb R)\\
a_k&1&2&4&4&8&8&8&8&16
\end{array}
$$

分類表の8周期性$\operatorname{Cl}_{0,k+8}(\mathbb R)\cong\operatorname{Cl}_{0,k}(\mathbb R)\otimes M_{16}(\mathbb R)$により、行列のサイズが16倍になります。

&&&fml 既約加群の次元の周期性
$$
a_{k+8}=16\,a_k
$$
&&&

$\operatorname{Cl}_{0,k}(\mathbb R)$の有限次元の加群は既約加群の直和になることが知られています（本記事では主張として使います）。したがって、$\mathbb R^n$が$\operatorname{Cl}_{0,k}(\mathbb R)$の加群になるなら、$a_k$は$n$を割り切ります。逆に$a_k$が$n$を割り切れば、既約加群$\mathbb R^{a_k}$を$n/a_k$個並べた直和として$\mathbb R^n$が加群になります。

&&&prop 加群になる条件 [prop-divides]
$\mathbb R^n$が$\operatorname{Cl}_{0,k}(\mathbb R)$の加群になるのは、$a_k$が$n$を割り切るときに限ります。
&&&

$a_k$は$k$について単調に増えるので、$k$を大きくしていくと、どこかで$a_k$が$n$を割り切らなくなります。

## 既約加群の具体形

既約加群は、数の左からの積で具体的に作れます。加群を既約加群の直和に分解すると各成分の次元は$a_k$なので、次元がちょうど$a_k$の加群は成分が1つしかなく、既約です。

&&&ex 既約加群の構成
- **$k=1$**：$\mathbb C$に$J_1=i$を掛けます。
- **$k=2,3$**：$\mathbb H$に$\mathbf i,\mathbf j$（$k=3$では$\mathbf k$も）を左から掛けます。$k=3$では擬スカラー$\omega=J_1J_2J_3$が$\mathbf i\mathbf j\mathbf k=-1$を掛ける写像、すなわちスカラー$-1$になります。符号を反転した$-\mathbf i,-\mathbf j,-\mathbf k$を掛ければ$\omega=+1$となり、これが直和型$2\mathbb H$のもう1つの既約加群です。
- **$k=4,\dots,7$**：$\mathbb O$に虚数単位$e_1,\dots,e_k$を左から掛けます。$k=7$では$k=3$と同じく、擬スカラーがスカラーになります。
- **$k=8$**：$\mathbb R^{16}=\mathbb O\oplus\mathbb O$に、$e_1,\dots,e_7$の左からの積（左作用）$L_1,\dots,L_7$を使って次のように作用させます。

$$
J_i=\begin{pmatrix}L_i&0\\0&-L_i\end{pmatrix}\ (i=1,\dots,7),\qquad J_8=\begin{pmatrix}0&-I\\I&0\end{pmatrix}
$$

$J_8$は$J_1,\dots,J_7$と反交換し、$J_8^2=-I$です。
&&&

次元の増え方が速いことが見て取れます。$k=4$から$k=7$までは8次元の$\mathbb O$で足りますが、$k=8$で次元を倍にしても、生成元は1つしか増えません。$k$が1増えるごとに$a_k$は$1$倍か$2$倍になり、8増えると16倍になります。生成元の数は1ずつしか増えないので、大きな$k$では$a_k$が$k$よりはるかに大きくなります。

# ラドン＝フルヴィッツ数

クリフォード加群から構成できる$S^{n-1}$上の正規直交な接ベクトル場の最大本数は、$\mathbb R^n$に作用するクリフォード代数$\operatorname{Cl}_{0,k}(\mathbb R)$の$k$の最大値です。

&&&def ラドン＝フルヴィッツ数
正の整数$n$に対して、$a_k$が$n$を割り切る最大の$k$に$1$を足したものを$\rho(n)$と書き、**ラドン＝フルヴィッツ数**と呼びます。

$$
\rho(n)-1=\max\{k\mid a_k\text{は}n\text{を割り切る}\}
$$
&&&

[[thm-module-fields]]と[[prop-divides]]により、$S^{n-1}$上には$\rho(n)-1$本の正規直交な接ベクトル場が取れます。[[thm-linear-fields]]により、1次式で与えられる正規直交な接ベクトル場は、これより多くは取れません。

$a_k$はすべて2の冪なので、$a_k$が$n$を割り切るかどうかは、$n$が2で何回割り切れるかだけで決まります。

&&&prop ラドン＝フルヴィッツ数の公式 [prop-rho]
$n=2^{4a+b}m$（$m$は奇数、$0\le b\le3$）と書くと

$$
\rho(n)=8a+2^b
$$

です。
&&&

&&&prf
$a_k=2^{d_k}$と書くと、表から$d_0,\dots,d_7=0,1,2,2,3,3,3,3$で、$a_{k+8}=16a_k$より$d_{k+8}=d_k+4$である。$d_k$は$k$について単調に増加し、$a_k$が$n$を割り切ることは$d_k\le4a+b$と同値である。$k=8a+j$（$0\le j\le7$）なら$d_k=4a+d_j$なので、条件は$d_j\le b$となる。$j=8$（$d_8=4>b$）や$k\ge8(a+1)$では条件を満たさないので、最大の$k$は$8a+j_{\max}$で、$j_{\max}$は$d_j\le b$を満たす最大の$j$である。$b=0,1,2,3$に対して$j_{\max}=0,1,3,7=2^b-1$だから、$\rho(n)=8a+2^b$である。
&&&

奇数の$n$では$\rho(n)=1$で、接ベクトル場は構成できません。これは毛玉の定理と合っています。$n$が2の冪のときの値は次のとおりです。

$$
\begin{array}{c|ccccccccc}
n&1&2&4&8&16&32&64&128&256\\\hline
\rho(n)&1&2&4&8&9&10&12&16&17
\end{array}
$$

$n=1,2,4,8$では$\rho(n)=n$で、上限の$n-1$本が取れます。$n=16$では接空間は15次元ですが、1次式で与えられる正規直交な接ベクトル場は8本までです。以後は$n$の増え方に比べて$\rho(n)$はゆっくりとしか増えません。

ここまでは、構成できる本数を求めました。1次式に限らない連続な接ベクトル場でも、これより多くは取れないことが知られています。

&&&thm アダムスの定理
$S^{n-1}$上の線形独立な連続な接ベクトル場の本数の最大値は、$\rho(n)-1$です。
&&&

上限の証明はK理論を用いるもので、本記事では扱いません。$n$が奇数の場合（$\rho(n)-1=0$）が毛玉の定理にあたります。

# 平行化可能な球面とフルヴィッツの定理

## 平行化可能な球面

$S^{n-1}$が平行化可能になるには、$n-1$本の接ベクトル場が必要です。

&&&prop 平行化可能な球面 [prop-parallel]
$\rho(n)=n$となるのは、$n=1,2,4,8$のときに限ります。
&&&

&&&prf
$n=2^{4a+b}m$（$m$は奇数、$0\le b\le3$）と書く。$a=0$なら$\rho(n)=2^b\le2^bm=n$で、等号は$m=1$のときに限るから、$n=1,2,4,8$である。$a\ge1$なら$n\ge16^a2^b$で、$16^a-1\ge15a>8a$より$n-\rho(n)\ge2^b(16^a-1)-8a>0$となり、$\rho(n)<n$である。
&&&

$n=1,2,4,8$の場合は、実数・複素数・四元数・八元数による例がちょうどこれにあたります。アダムスの定理（接ベクトル場）を認めれば、平行化可能な球面は$S^0,S^1,S^3,S^7$に限られます。以前の記事では、単位元つきの連続な積を持てる球面が同じ$S^0,S^1,S^3,S^7$に限られること（ホップ不変量1についてのアダムスの定理）を紹介しました。[[7shi-lie7]]

## フルヴィッツの定理の次元の部分

[[prop-parallel]]の$n=1,2,4,8$は、数体系の次元の並びと同じです。この一致は偶然ではなく、ノルムの乗法性を持つ数体系からは、必ず1次式の接ベクトル場が作られます。

&&&thm ノルムの乗法性を持つ代数の次元 [thm-hurwitz]
$\mathbb R^n$に単位元$1$を持つ双線形な積が定義され、通常の長さについて$|uv|=|u||v|$がすべての$u,v$で成り立つなら、$n=1,2,4,8$です。
&&&

&&&prf
$|u\boldsymbol x|^2=|u|^2|\boldsymbol x|^2$の$u$を$u+v$に置き換えて展開し、$|u\boldsymbol x|^2$と$|v\boldsymbol x|^2$の項を消すと

$$
\langle u\boldsymbol x,v\boldsymbol x\rangle=\langle u,v\rangle|\boldsymbol x|^2
$$

が得られる。$|1\cdot\boldsymbol x|=|1||\boldsymbol x|$より$|1|=1$である。$1$に直交する部分空間の正規直交基底を$u_1,\dots,u_{n-1}$とすると、上式で$v=1$と置けば$\langle u_i\boldsymbol x,\boldsymbol x\rangle=\langle u_i,1\rangle|\boldsymbol x|^2=0$、$u=u_i$、$v=u_j$と置けば$\langle u_i\boldsymbol x,u_j\boldsymbol x\rangle=\delta_{ij}|\boldsymbol x|^2$である。積は双線形なので$\boldsymbol x\mapsto u_i\boldsymbol x$は1次式であり、これらは$S^{n-1}$上の$n-1$本の正規直交な接ベクトル場になる。[[thm-linear-fields]]により$\mathbb R^n$は$\operatorname{Cl}_{0,n-1}(\mathbb R)$の加群であり、$\rho(n)-1\ge n-1$となる。$\rho(n)\le n$と合わせて$\rho(n)=n$だから、[[prop-parallel]]より$n=1,2,4,8$である。
&&&

証明では結合法則を使っていません。1次式の接ベクトル場を与えるのは、ノルムの乗法性だけです。また、構成した接ベクトル場が1次式なので、[[thm-linear-fields]]だけで結論が得られ、アダムスの定理には依存していません。以前の記事で述べたフルヴィッツの定理（ノルムの乗法性を持つ単位元付きの実代数は、同型を除いて実数・複素数・四元数・八元数に限られる）のうち、次元が$1,2,4,8$に限られる部分がこれで示せました。各次元で代数が同型を除いて1つに決まることには、本記事では立ち入りません。[[7shi-lie7]][[7shi-sedenion]]

十六元数でノルムの乗法性が崩れることも、この観点から見直せます。16次元では$\rho(16)=9$で、1次式の正規直交な接ベクトル場は8本しか取れません。虚数単位の15個の左からの積が正規直交な接ベクトル場になることは、クリフォード代数の既約加群の次元から不可能です。

&&&rem 八元数の7番目の虚数単位
正の3つ組を$123,145,176,246,257,347,365$とする規約（以前の記事と同じもの）では、$k=7$の既約加群を与える八元数の左からの積について、$L_1L_2\cdots L_6=L_7$が成り立ちます。以前の記事では、$L_1,\dots,L_6$で$\operatorname{Cl}_{0,6}(\mathbb R)$を生成すると、7番目の虚数単位$L_7$が擬スカラー$\omega$として現れることを見ました。$L_1,\dots,L_7$を$\operatorname{Cl}_{0,7}(\mathbb R)$の生成元の作用と見れば、その擬スカラーは$L_1\cdots L_6L_7=L_7^2=-I$で、スカラーとして作用します。この規約で定めた八元数への作用は、直和型$\operatorname{Cl}_{0,7}(\mathbb R)\cong2M_8(\mathbb R)$の2つの既約加群のうち、擬スカラーが$-I$として作用する方です。[[7shi-cl6]]
&&&

# まとめ

$S^{n-1}$上の正規直交な接ベクトル場を、クリフォード代数の加群から構成しました。

- **接ベクトル場とクリフォード加群**：$\mathbb R^n$が$\operatorname{Cl}_{0,k}(\mathbb R)$の直交的な加群なら、生成元の作用$\boldsymbol x\mapsto J_i\boldsymbol x$が$k$本の正規直交な接ベクトル場を与えます。逆に、1次式で与えられる正規直交な接ベクトル場は、クリフォード代数の作用を定めます。$J_i^2=-I$が接することに、反交換が直交することに対応します。
- **既約加群の次元**：分類表の第1行から、$\operatorname{Cl}_{0,k}(\mathbb R)$の既約加群の次元$a_k$を読みました。$\mathbb R^n$が加群になるのは$a_k$が$n$を割り切るときです。
- **ラドン＝フルヴィッツ数**：$\rho(n)-1$本の接ベクトル場が構成でき、$n=2^{4a+b}m$なら$\rho(n)=8a+2^b$です。これが連続な接ベクトル場の本数の最大であることがアダムスの定理です。
- **平行化可能性とフルヴィッツの定理**：$\rho(n)=n$となるのは$n=1,2,4,8$だけです。ノルムの乗法性を持つ代数は$n-1$本の1次式の接ベクトル場を与えるので、その次元も$1,2,4,8$に限られます。

既約加群と、その出所をまとめます。

| $k$ | $\operatorname{Cl}_{0,k}(\mathbb R)$ | $a_k$ | 既約加群 | 接ベクトル場を持つ球面 |
|:---:|:---:|:---:|:---|:---|
| 0 | $\mathbb R$ | 1 | $\mathbb R$ | $S^0$ |
| 1 | $\mathbb C$ | 2 | 複素数 | $S^1$（平行化可能） |
| 2 | $\mathbb H$ | 4 | 四元数 | $S^3$ |
| 3 | $2\mathbb H$ | 4 | 四元数 | $S^3$（平行化可能） |
| 4〜6 | $M_2(\mathbb H)$・$M_4(\mathbb C)$・$M_8(\mathbb R)$ | 8 | 八元数 | $S^7$ |
| 7 | $2M_8(\mathbb R)$ | 8 | 八元数 | $S^7$（平行化可能） |
| 8 | $M_{16}(\mathbb R)$ | 16 | 八元数の倍加 | $S^{15}$（8本） |

各行の球面は、既約加群そのものの単位球面$S^{a_k-1}$で、$k$本の正規直交な接ベクトル場を持ちます。
