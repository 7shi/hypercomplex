[[7shi-clif1]]では、実クリフォード代数の型が8周期で、複素クリフォード代数の型が2周期で繰り返すことを示し、この周期性が「K理論におけるボット周期性定理と対応することが知られています」と述べるに留めました。本記事では、K理論が何を数える理論なのかを見渡します。空間の各点にベクトル空間を連続に並べたものをベクトル束と呼び、ベクトル束を直和で足し、差を許して群にしたものがK群です。メビウスの帯・球面の接束・ホップ束を例にこの流れをたどり、円周の上ではメビウスの帯が位数2の元になることを確かめます。そのうえで、球面のK群の周期性（ボット周期性）とクリフォード代数の分類表の対応、作用素の指数がK群の元として現れることを、主張として紹介します。

本記事は概観が目的で、証明するのはメビウスの帯まわりの初等的な事実に限ります。ボット周期性定理、アティヤ＝ボット＝シャピロの定理、アティヤ＝シンガーの指数定理は主張として述べ、証明は本記事では扱いません。

# ベクトル束

## 定義

&&&def ベクトル束
空間$X$の各点$x$に$n$次元の実ベクトル空間$E_x$を連続に割り当てたものを、$X$上の階数$n$の**ベクトル束**と呼び、$E$と書きます。各点の近くでは、$E$は積$U\times\mathbb R^n$（$U$は$x$の近傍）の形をしているものとします。複素ベクトル空間$\mathbb C^n$を割り当てたものを複素ベクトル束と呼びます。

$X$全体で積$X\times\mathbb R^n$になっているものを**自明束**と呼び、$\underline{\mathbb R}^n$と書きます。2つのベクトル束$E,F$の間に、各点で線形同型$E_x\to F_x$を与え、$x$について連続に変わる写像があるとき、$E$と$F$は**同型**であると言い、$E\cong F$と書きます。
&&&

「連続に割り当てる」「局所的に積の形をしている」の形式的な定義には立ち入らず、以下では例を通して具体的に扱います。

ベクトル束$E$の**切断**とは、各点$x$に$E_x$の元$s(x)$を連続に対応させるものです。自明束$\underline{\mathbb R}^n$の切断は$\mathbb R^n$に値を取る連続関数にほかならず、たとえば定数関数$s(x)=(1,0,\dots,0)$は、どこでも$0$にならない切断です。したがって、どこでも$0$にならない切断を持たないベクトル束は自明束と同型になりません。これが、自明でないことを示す最も手軽な方法です。

## 円周上のベクトル束

円周$S^1$を区間$[0,2\pi]$の両端を同一視したものと見ます。区間の上の自明束$[0,2\pi]\times\mathbb R^n$を用意し、両端のファイバーを可逆行列$g\in\operatorname{GL}(n,\mathbb R)$で貼り合わせます。すなわち、$\theta=2\pi$のファイバーの元$v$を、$\theta=0$のファイバーの元$gv$と同一視します。こうして得られる$S^1$上のベクトル束を$E_g$と書き、$g$を**貼り合わせ行列**と呼びます。円周上のベクトル束はすべてこの形で得られることが知られています（区間上のベクトル束が自明であることによります）。

$E_g$の切断は、連続関数$s:[0,2\pi]\to\mathbb R^n$で両端の条件

$$
s(0)=g\,s(2\pi)
$$

を満たすものです。$g=I_n$なら$s(0)=s(2\pi)$となり、$E_{I_n}$は自明束$\underline{\mathbb R}^n$です。

&&&ex 円柱とメビウスの帯 [ex-moebius]
階数1の場合、貼り合わせ行列は$0$でない実数です。$g=1$で貼り合わせたものが円柱$\underline{\mathbb R}$、$g=-1$で貼り合わせたものが**メビウスの帯**$M$です。$M$の切断は$s(0)=-s(2\pi)$を満たす関数で、$[0,2\pi]$を周期的に延長すれば$s(\theta+2\pi)=-s(\theta)$と書けます。たとえば$s(\theta)=\cos(\theta/2)$は$M$の切断です。
&&&

&&&prop メビウスの帯は自明でない [prop-moebius]
$M$の切断は必ず零点を持つ。したがって$M\not\cong\underline{\mathbb R}$である。
&&&

&&&prf
$s$を$M$の切断とすると、$s(0)=-s(2\pi)$より$s(0)$と$s(2\pi)$は異符号か、ともに$0$である。どちらの場合も、中間値の定理により$[0,2\pi]$のどこかで$s=0$となる。自明束$\underline{\mathbb R}$は零点のない切断$s=1$を持つので、$M$とは同型でない。
&&&

$S^1$上のベクトル束の同型も、貼り合わせ行列の言葉で書けます。$E_g$から$E_{g'}$への同型は、各$\theta$での可逆行列$A(\theta)$で$\theta$について連続なものです。$\theta=2\pi$の元$v$と$\theta=0$の元$gv$は同じ点を表すので、その像$A(2\pi)v$と$A(0)gv$も$E_{g'}$で同じ点を表す必要があります。すなわち$g'A(2\pi)=A(0)g$で、

$$
g'=A(0)\,g\,A(2\pi)^{-1}
$$

が同型の条件です。

&&&prop 貼り合わせ行列の変形 [prop-homotopy]
$g$から$g'$へ$\operatorname{GL}(n,\mathbb R)$の中で連続に変形できるなら、$E_g\cong E_{g'}$である。複素ベクトル束でも同様である。
&&&

&&&prf
$g$から$g'$への変形を$g_t$（$0\le t\le1$、$g_0=g$、$g_1=g'$）とし、$A(\theta)=g_{\theta/2\pi}^{-1}\,g$と置く。$A(\theta)$は連続で可逆であり、$A(0)=I_n$、$A(2\pi)=g'^{-1}g$だから、$A(0)\,g\,A(2\pi)^{-1}=g\cdot g^{-1}g'=g'$となり、同型の条件を満たす。
&&&

## 球面上のベクトル束

球面$S^k$も同じ方法で扱えます。$S^k$を北半球$D_+$と南半球$D_-$に分け、それぞれの上の自明束を、赤道$S^{k-1}$の各点$\boldsymbol x$で可逆行列$g(\boldsymbol x)$により貼り合わせます。$g:S^{k-1}\to\operatorname{GL}(n,\mathbb R)$を**貼り合わせ関数**と呼びます。$k=1$では赤道が2点$S^0=\{\pm1\}$からなり、一方の点での貼り合わせを$1$に揃えれば、円周の場合の貼り合わせ行列に戻ります。

&&&ex 球面の接束 [ex-tangent]
$S^2$の各点$\boldsymbol x$に接平面$T_{\boldsymbol x}S^2=\{\boldsymbol v\in\mathbb R^3\mid\boldsymbol v\perp\boldsymbol x\}$を割り当てたものは、階数2の実ベクトル束です。これを**接束**と呼び、$TS^2$と書きます。$TS^2$の切断は$S^2$上の接ベクトル場です。

$S^2$上には、どこでも$0$にならない連続な接ベクトル場が存在しないことが知られています（毛玉の定理）。したがって$TS^2$は自明束$\underline{\mathbb R}^2$と同型ではありません。毛玉の定理の証明は本記事では扱いません。
&&&

&&&ex ホップ束 [ex-hopf]
$\mathbb C^2$の原点を通る複素直線$L$全体は、同次座標$[\alpha:\beta]$で表され、$z=\alpha/\beta$により球面$S^2$と同一視されます（[[7shi-homog]]・[[7shi-c2s2]]）。各点$L$に、その直線$L$自身を1次元の複素ベクトル空間として割り当てたものが、$S^2$上の複素ベクトル束になります。これを**ホップ束**と呼び、$H$と書きます。

$\beta\ne0$の範囲では$L$は$(z,1)$で張られ、ファイバーの元は$\lambda(z,1)$と書けます。$\alpha\ne0$の範囲では$w=\beta/\alpha$として$\mu(1,w)$と書けます。両方の範囲が重なるところで$\lambda(z,1)=\mu(1,w)$を比べると

$$
\mu=z\lambda
$$

です。赤道$|z|=1$で2つの座標$\lambda,\mu$が単位複素数$z$を掛けることで移り合うので、$H$は貼り合わせ関数$g(z)=z$で作ったベクトル束です。

ファイバーの中の長さ1のベクトルだけを集めると、$\mathbb C^2$の単位ベクトル全体、すなわち$S^3$になります。各ベクトルをそれが属する直線$L\in S^2$に送る写像が、ホップファイブレーション$S^3\to S^2$です。$H$が自明でないことは知られていますが、本記事では主張に留めます。
&&&

# 直和と安定同値

## 直和

2つのベクトル束$E,F$の**直和**$E\oplus F$は、各点のベクトル空間の直和$E_x\oplus F_x$を割り当てたベクトル束です。階数は足し算になります。$S^1$上では、貼り合わせ行列がブロック対角行列$\operatorname{diag}(g,g')$になり、$E_g\oplus E_{g'}=E_{\operatorname{diag}(g,g')}$です。

&&&ex 接束と法線 [ex-tangent-normal]
$S^2$の各点$\boldsymbol x$で、$\mathbb R^3$は接平面と法線の直和$\mathbb R^3=T_{\boldsymbol x}S^2\oplus\mathbb R\boldsymbol x$に分かれます。$(\boldsymbol v,t)\mapsto\boldsymbol v+t\boldsymbol x$は各点で線形同型で、$\boldsymbol x$について連続なので

$$
TS^2\oplus\underline{\mathbb R}\cong\underline{\mathbb R}^3
$$

です。自明でない$TS^2$に自明束を1つ足すと、自明束になります。
&&&

&&&prop メビウスの帯2枚 [prop-mm]
$M\oplus M\cong\underline{\mathbb R}^2$である。
&&&

&&&prf
$M\oplus M$の貼り合わせ行列は$-I_2$である。回転行列

$$
R(\varphi)=\begin{pmatrix}\cos\varphi&-\sin\varphi\\\sin\varphi&\cos\varphi\end{pmatrix}
$$

は$\varphi\in[0,\pi]$で$R(0)=I_2$から$R(\pi)=-I_2$まで$\operatorname{GL}(2,\mathbb R)$の中を連続に動く。[[prop-homotopy]]により$E_{-I_2}\cong E_{I_2}=\underline{\mathbb R}^2$である。
&&&

$M$は1枚では自明でないのに、2枚重ねると自明になります。2枚の帯のねじれを、平面の半回転でまとめてほどけるからです。

## 安定同値

&&&def 安定同値
ベクトル束$E,F$が、ある$k$について$E\oplus\underline{\mathbb R}^k\cong F\oplus\underline{\mathbb R}^k$を満たすとき、$E$と$F$は**安定同値**であると言います。
&&&

[[ex-tangent-normal]]から、$TS^2$は$\underline{\mathbb R}^2$と安定同値です。$TS^2$と$\underline{\mathbb R}^2$の違いは、自明束を足すと見えなくなります。安定同値は同型より粗い見方で、接ベクトル場の有無のような違いを捨てています。その代わりに、次に見るとおり、足し算と引き算ができる群の構造が得られます。

# K群

## 定義

ベクトル束の同型類の全体は、直和を演算として、可換なモノイドになります。単位元は階数0のベクトル束です（[[7shi-group]]）。しかし、階数が足し算になるので、階数0の束以外に逆元はありません。自然数から整数を作るときと同じく、形式的な差を導入して群にします。

&&&def K群
空間$X$上の実ベクトル束の組$(E,F)$を「$E-F$」を表すものと考え、

$$
(E,F)\sim(E',F')\iff\text{あるベクトル束}G\text{について}\ E\oplus F'\oplus G\cong E'\oplus F\oplus G
$$

で同一視したものの全体を$KO(X)$と書き、$(E,F)$の類を$[E]-[F]$と書きます。演算は$([E]-[F])+([E']-[F'])=[E\oplus E']-[F\oplus F']$で、$KO(X)$は可換群になります。複素ベクトル束から同じように作った群を$K(X)$と書きます。
&&&

自然数の場合は$(a,b)\sim(c,d)\iff a+d=b+c$で十分ですが、ベクトル束では$G$を足してから比べる必要があります。[[ex-tangent-normal]]のように$TS^2\oplus\underline{\mathbb R}\cong\underline{\mathbb R}^2\oplus\underline{\mathbb R}$でありながら$TS^2\not\cong\underline{\mathbb R}^2$となり、両辺から同じものを取り除く簡約律が成り立たないためです。$KO(S^2)$では$[TS^2]-[\underline{\mathbb R}^2]=0$で、K群は$TS^2$と自明束を区別しません。

球面や円周のように有界で閉じた空間では、$G$を自明束に取り替えてよいことが知られています。このとき$[E]-[F]=[E']-[F']$は$E\oplus F'$と$E'\oplus F$の安定同値にあたり、K群は「安定同値で見たベクトル束の差」の群です。以下、階数$n$の自明束の類$[\underline{\mathbb R}^n]$を単に$n$と書きます。

$X$が連結なら、$[E]-[F]$に階数の差$\operatorname{rank}E-\operatorname{rank}F\in\mathbb Z$を対応させる写像があります。その核を**簡約K群**と呼び、$\widetilde{KO}(X)$（複素では$\tilde K(X)$）と書きます。$KO(X)\cong\mathbb Z\oplus\widetilde{KO}(X)$で、$\widetilde{KO}(X)$の元は階数$n$の$E$を使って$[E]-n$と書けます。簡約K群は、階数を除いて、ベクトル束が自明からどれだけずれているかを記録します。

## 1点と円周

&&&ex 1点 [ex-point]
$X$が1点$\mathrm{pt}$なら、ベクトル束は1つのベクトル空間で、同型類は次元で決まります。$[\mathbb R^m]-[\mathbb R^n]$は整数$m-n$と同一視され

$$
KO(\mathrm{pt})\cong\mathbb Z,\qquad\widetilde{KO}(\mathrm{pt})=0
$$

です。複素でも$K(\mathrm{pt})\cong\mathbb Z$です。K群の元は、1点の上では「ベクトル空間の次元の差」です。一般の空間$X$の上のK群は、次元の差という考え方を、$X$の上に連続に広げたものと見ることができます。
&&&

&&&thm 円周のK群 [thm-circle]
$\widetilde{KO}(S^1)\cong\mathbb Z_2$で、生成元は$[M]-1$である。
&&&

&&&prf
まず$[M]-1$の位数が2であることを示す。[[prop-mm]]より$2([M]-1)=[M\oplus M]-2=0$である。

次に$[M]-1\ne0$を示す。$[M]-1=0$なら、あるベクトル束$G$について$M\oplus G\cong\underline{\mathbb R}\oplus G$となる。$G=E_h$とすると、両辺の貼り合わせ行列は$\operatorname{diag}(-1,h)$と$\operatorname{diag}(1,h)$で、行列式は$-\det h$と$\det h$であり符号が異なる。一方、同型の条件$g'=A(0)\,g\,A(2\pi)^{-1}$から

$$
\det g'=\frac{\det A(0)}{\det A(2\pi)}\det g
$$

であり、$\det A(\theta)$は$[0,2\pi]$上で$0$にならない連続関数だから符号を変えない。したがって$\det g'$と$\det g$は同符号で、同型は貼り合わせ行列の行列式の符号を保つ。これは矛盾だから$[M]-1\ne0$である。

最後に、$\widetilde{KO}(S^1)$の元がこの2つに限られることを示す。$\operatorname{GL}(n,\mathbb R)$のうち行列式が正のものの全体は、連続につながっていることが知られている。これを認めると、[[prop-homotopy]]により、$\det g>0$なら$E_g\cong\underline{\mathbb R}^n$、$\det g<0$なら$g$は$\operatorname{diag}(-1,1,\dots,1)$に連続に変形できて$E_g\cong M\oplus\underline{\mathbb R}^{n-1}$となる。よって$\widetilde{KO}(S^1)$の元は$0$か$[M]-1$である。
&&&

複素ベクトル束では事情が異なります。階数1の複素ベクトル束$E_{-1}$（$-1$で貼り合わせた複素直線束）は、$\mathbb C^\times$の中の道$e^{i\pi t}$（$0\le t\le1$）で$-1$を$1$に変形できるので、[[prop-homotopy]]により自明です。一般に$\operatorname{GL}(n,\mathbb C)$は連続につながっていることが知られており、$S^1$上の複素ベクトル束はすべて自明で、$\tilde K(S^1)=0$となります。実数では$-1$を$1$につなぐ道が$\mathbb R^\times$の中にないのに対し、複素数では$\mathbb C^\times$の中を回ってつなげます。実と複素の違いは、円周の上で既に現れています。

# 球面のK群とボット周期性

## 周期性

球面のK群は次のように知られています。

&&&thm 球面のK群（ボット周期性）
$n\ge1$について、複素K群は

$$
\tilde K(S^n)\cong\begin{cases}\mathbb Z&(n\text{が偶数})\\0&(n\text{が奇数})\end{cases}
$$

で、周期2で繰り返す。実K群は周期8で繰り返し、$n=1,\dots,8$では次のとおりである。

$$
\begin{array}{c|cccccccc}
n&1&2&3&4&5&6&7&8\\\hline
\widetilde{KO}(S^n)&\mathbb Z_2&\mathbb Z_2&0&\mathbb Z&0&0&0&\mathbb Z
\end{array}
$$
&&&

周期性の部分が**ボット周期性定理**で、証明は本記事では扱いません。$n=1$の$\mathbb Z_2$は[[thm-circle]]で確かめたもので、$\tilde K(S^1)=0$とあわせて、表の最初の列を与えます。

周期は[[7shi-clif1]]の分類表と同じです。実クリフォード代数の型は$p-q\pmod8$で決まり、複素クリフォード代数の型は$n\pmod2$で決まりました。

## 単位数による貼り合わせ

周期が一致するのは偶然ではありません。球面のK群の生成元は、クリフォード代数から作ることができます。

$S^k$上のベクトル束を作るには、赤道$S^{k-1}$上の貼り合わせ関数が必要です。$\boldsymbol x=(x_0,\dots,x_{k-1})\in S^{k-1}$に対して、$\operatorname{Cl}_{0,k-1}(\mathbb R)$の生成元$e_1,\dots,e_{k-1}$（$e_l^2=-1$、$l\ne m$なら$e_le_m=-e_me_l$）を使って

$$
g(\boldsymbol x)=x_0+x_1e_1+\dots+x_{k-1}e_{k-1}
$$

と置きます。これはスカラーとベクトルの和で、[[7shi-cla5]]でパラベクトルと呼んだ形です。$\boldsymbol v=\sum_lx_le_l$とすると、反交換する項が打ち消し合って

$$
\boldsymbol v^2=\sum_lx_l^2e_l^2=-\sum_lx_l^2
$$

となるので、

$$
(x_0+\boldsymbol v)(x_0-\boldsymbol v)=x_0^2-\boldsymbol v^2=\sum_{i=0}^{k-1}x_i^2=1
$$

です。したがって$g(\boldsymbol x)$は可逆で、$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群$W$に作用させれば$g:S^{k-1}\to\operatorname{GL}(W)$が得られます。可逆性は、生成元の1次結合の2乗が2次形式$-|\boldsymbol v|^2$になること、すなわちクリフォード代数が2次形式の平方根を与えることから来ています。

&&&ex 実数・複素数・四元数・八元数
$k=1,2,4,8$で$W$を次のように取ると、$g(\boldsymbol x)$は長さ1の数（単位数）を掛ける写像になります。

- **$k=1$**：$\operatorname{Cl}_{0,0}(\mathbb R)=\mathbb R$、$W=\mathbb R$で、$g(x_0)=x_0=\pm1$です。2点$S^0$での貼り合わせが$+1$と$-1$なので、メビウスの帯$M$が得られます。
- **$k=2$**：$\operatorname{Cl}_{0,1}(\mathbb R)\cong\mathbb C$、$W=\mathbb C$で、$e_1=i$とすれば$g=x_0+x_1i$は単位複素数を掛ける写像です。[[ex-hopf]]のホップ束が得られます。
- **$k=4$**：$\operatorname{Cl}_{0,3}(\mathbb R)$は四元数$\mathbb H$の$\mathbf i,\mathbf j,\mathbf k$の左からの積で$W=\mathbb H$に作用し、$g$は単位四元数を掛ける写像です。得られるのは四元数のホップ束です（[[7shi-hopfext]]）。
- **$k=8$**：$\operatorname{Cl}_{0,7}(\mathbb R)$は八元数$\mathbb O$の7つの虚数単位の左からの積で$W=\mathbb O$に作用し、$g$は単位八元数を左から掛ける写像です。得られるのは八元数のホップ束です（[[7shi-hopfext]]）。
&&&

これらはそれぞれ$\widetilde{KO}(S^1)$、$\tilde K(S^2)$、$\widetilde{KO}(S^4)$、$\widetilde{KO}(S^8)$の生成元を与えることが知られています（$\widetilde{KO}(S^2)$の生成元も、ホップ束を実ベクトル束と見たものから得られます）。一般の$k$についても、$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群から上の$g$で作ったベクトル束が$\widetilde{KO}(S^k)$を生成し、$\widetilde{KO}(S^k)$そのものがクリフォード加群の分類から計算できます。これを**アティヤ＝ボット＝シャピロの定理**と呼びます。証明は本記事では扱いません。

## 分類表との比較

$\widetilde{KO}(S^k)$を、生成元を作る代数$\operatorname{Cl}_{0,k-1}(\mathbb R)$の型と並べます。代数の型は[[7shi-clif1]]の分類表の第1行（$p=0$）から引きます。

$$
\begin{array}{c|cccccccc}
k&1&2&3&4&5&6&7&8\\\hline
\operatorname{Cl}_{0,k-1}(\mathbb R)&\mathbb R&\mathbb C&\mathbb H&2\mathbb H&\mathbb H(2)&\mathbb C(4)&\mathbb R(8)&2\mathbb R(8)\\
\widetilde{KO}(S^k)&\mathbb Z_2&\mathbb Z_2&0&\mathbb Z&0&0&0&\mathbb Z
\end{array}
$$

$\mathbb Z$が現れる$k=4,8$は、代数が直和型$2\mathbb H,2\mathbb R(8)$になる位置です。直和型の代数には既約加群が2つあり、擬スカラーの符号で区別されます（[[7shi-clif1]]の「擬スカラーと型の判別」）。複素でも同じ対応が成り立ち、$\operatorname{Cl}_{k-1}(\mathbb C)$は$k$が偶数のときに直和型となり、ちょうどそのとき$\tilde K(S^k)\cong\mathbb Z$です。$k=1,2$の$\mathbb Z_2$は実数に特有のもので、その現れ方には既約加群の次元が関わりますが、本記事では扱いません。

# 指数とK理論

## 有限次元の指数

有限次元のベクトル空間の間の線形写像$T:V\to W$に対して、核$\ker T$と余核$\operatorname{coker}T=W/\operatorname{im}T$の次元の差を**指数**と呼びます。

&&&fml 有限次元の指数
$$
\operatorname{ind}T=\dim\ker T-\dim\operatorname{coker}T=\dim V-\dim W
$$
&&&

$\dim\ker T+\dim\operatorname{im}T=\dim V$と$\dim\operatorname{coker}T=\dim W-\dim\operatorname{im}T$から、この等式が従います。核と余核の次元は$T$ごとに変わりますが、差は$T$によらず$V,W$の次元だけで決まり、正方行列なら常に$0$です。指数は$[V]-[W]$、すなわち$K(\mathrm{pt})\cong\mathbb Z$の元です（[[ex-point]]）。

## 無限次元の指数

無限次元では、指数は$T$によって変わる量になります。2乗和が有限な数列の空間で、右へのずらし

$$
S(a_0,a_1,a_2,\dots)=(0,a_0,a_1,\dots)
$$

を考えます。$S\boldsymbol a=0$なら$\boldsymbol a=0$なので$\ker S=0$です。$S$の像は先頭が$0$の数列全体なので、余核は先頭の成分だけを残す1次元の空間です。したがって$\operatorname{ind}S=0-1=-1$です。同じ空間から同じ空間への写像なのに、指数が$0$になりません。

核と余核が有限次元になる作用素を**フレドホルム作用素**と呼びます。その指数は、作用素を連続に変形しても変わらない整数であることが知られています。指数は、解析的に定義されながら位相的な性質を持つ量です。さらに、空間$X$の点で連続に変わるフレドホルム作用素の族を考えると、各点の核と余核の差が$K(X)$の元を定めることが知られています。1点の上の次元の差が$K(\mathrm{pt})$の元であったことを、$X$の上に広げたものです。

## アティヤ＝シンガーの指数定理

微分作用素の指数、すなわち微分方程式の解の空間の次元から余核の次元を引いたものは、解析の側で定義される量です。**アティヤ＝シンガーの指数定理**は、閉じた多様体上の楕円型微分作用素について、この**解析的指数**が、作用素の最高階の部分（主表象）から位相的に定まる**位相的指数**に等しいことを主張します。位相的指数は、主表象が定めるK群の元から計算されます。証明は本記事では扱いません。

ディラック作用素はその代表例です。[[7shi-cla5]]のディラック作用素$D=\sum_ae_a\partial_a$（$\operatorname{Cl}_{n,0}(\mathbb R)$、$e_a^2=1$）で$\partial_a$を変数$\xi_a$に置き換えると、ベクトル$\boldsymbol\xi=\sum_a\xi_ae_a$を掛ける写像

$$
\boldsymbol\xi\mapsto\boldsymbol\xi\cdot
$$

が得られます。これが$D$の主表象です（文献によっては虚数単位の定数倍を付けます）。$\boldsymbol\xi^2=|\boldsymbol\xi|^2$から、$\boldsymbol\xi\ne0$ならこの写像は可逆です。左から$e_0$を掛けると$e_0\boldsymbol\xi=\xi_0+\sum_l\xi_l\,e_0e_l$となり、$h_l=e_0e_l$は$h_l^2=-1$を満たして互いに反交換するので（[[7shi-cla5]]）、これは球面の生成元を作ったパラベクトル$g(\boldsymbol x)$と同じ形です。ディラック作用素の主表象と、球面のK群の生成元を与える貼り合わせ関数は、同じ写像です。

# まとめ

K理論の基本的な概念を、定義と最小限の例で見渡しました。

- **ベクトル束**：空間の各点にベクトル空間を連続に割り当てたものです。円周や球面の上では、自明束を貼り合わせ行列（貼り合わせ関数）で貼り合わせて作れます。メビウスの帯、球面の接束、ホップ束が自明でない例です。
- **直和と安定同値**：直和で階数が足されます。$TS^2\oplus\underline{\mathbb R}\cong\underline{\mathbb R}^3$、$M\oplus M\cong\underline{\mathbb R}^2$のように、自明束を足すと消える違いがあり、それを無視するのが安定同値です。
- **K群**：ベクトル束の同型類のモノイドに形式的な差を導入した群です。1点では次元の差$\mathbb Z$で、円周では$\widetilde{KO}(S^1)\cong\mathbb Z_2$（生成元$[M]-1$）、$\tilde K(S^1)=0$です。
- **ボット周期性**：球面のK群は、実では周期8、複素では周期2で繰り返します。
- **指数**：核と余核の次元の差です。有限次元では$\dim V-\dim W$で決まり、無限次元では作用素の位相的な性質を反映します。アティヤ＝シンガーの指数定理は、微分作用素の指数が主表象から位相的に計算できることを主張します。

K理論の側の概念と、クリフォード代数の側の対応物を並べます。

| K理論 | クリフォード代数 |
|:---|:---|
| 球面のK群の周期（実8・複素2） | 分類表の周期（$p-q\bmod8$・$n\bmod2$） |
| $S^k$の生成元の貼り合わせ関数 | $\operatorname{Cl}_{0,k-1}(\mathbb R)$加群上のパラベクトル$x_0+\sum_lx_le_l$ |
| 貼り合わせ関数の可逆性 | $\boldsymbol v^2=-\lvert\boldsymbol v\rvert^2$（2次形式の平方根） |
| $S^1,S^2,S^4,S^8$の生成元 | 実数・複素数・四元数・八元数の単位数を掛ける写像 |
| $\widetilde{KO}(S^k)$に$\mathbb Z$が現れる位置 | $\operatorname{Cl}_{0,k-1}(\mathbb R)$が直和型になる位置 |
| ディラック作用素の主表象 | ベクトル$\boldsymbol\xi$を掛ける写像$\boldsymbol\xi\mapsto\boldsymbol\xi\cdot$ |
