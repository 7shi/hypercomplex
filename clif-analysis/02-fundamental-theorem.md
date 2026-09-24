ディラック作用素$D$で微分した関数を領域上で積分すると、境界上の積分に書き換えられます。1変数の微積分学の基本定理を多次元へ持ち上げたこの式は、クリフォード代数に値を取る関数に対して1つの式として書けます。平面のグリーンの定理、空間の発散定理とストークスの定理は、通常は別々の定理として示されますが、ここでは領域の基本定理とその曲面版という2つの式のグレード成分として得られます。法線付きの境界要素$\boldsymbol n\,dS$や接線方向を持つ線素$d\boldsymbol x$のように、積分の要素に向きを代数の元として持たせるため、この種の積分を有向積分と呼びます。代数に値を取る関数の積分は、基底に関する成分ごとの積分として定めます。

# 1次元の基本定理

区間$[a,b]$上の関数$F$について、微積分学の基本定理は次の式です。

$$
\int_a^b F'(x)\,dx=F(b)-F(a)
$$

右辺は、区間の境界である2点$a,b$での値の差です。右端では外側が正の向き、左端では負の向きなので、境界の各点に**外向きの向き**$n=\pm1$を割り当てると、右辺は$n F$の境界上の和になります。

$$
\int_a^b F'(x)\,dx=\sum_{\text{境界}}nF=(+1)F(b)+(-1)F(a)
$$

左辺は領域全体での微分の積分、右辺は境界だけでの値です。この形を2次元・3次元へ持ち上げます。

# 2次元の基本定理

[[7shi-cla1]]と同じく、$\operatorname{Cl}_{2,0}(\mathbb R)$の生成元を$e_1,e_2$（$e_1^2=e_2^2=1$、$e_1e_2=-e_2e_1$）、擬スカラーを$I=e_1e_2$（$I^2=-1$）とし、ディラック作用素を

$$
D=e_1\partial_x+e_2\partial_y
$$

とします。$D$は左から$DF=e_1\partial_xF+e_2\partial_yF$と作用させます。

平面の有界な領域$M$の境界$\partial M$は、区分的に滑らかな閉曲線（穴があれば複数）とします。関数$F$は、$M$とその境界を含む開集合上で定義され、$\operatorname{Cl}_{2,0}(\mathbb R)$に値を取る$C^1$級関数とします。境界上の各点で外向きの単位法線ベクトルを$\boldsymbol n=n_1e_1+n_2e_2$、弧長の要素を$ds$と書きます。1次元の$n=\pm1$を、ベクトルとしての外向きに置き換えたものです。

&&&thm 2次元の基本定理 [thm-2d]
$$
\int_M DF\,dA=\oint_{\partial M}\boldsymbol nF\,ds
$$
&&&

左辺の$DF$と右辺の$\boldsymbol nF$はどちらも代数の中の積であり、積が非可換なので、$\boldsymbol n$は$D$と同じく$F$の左に置きます。

## 長方形での証明

まず$M$が長方形$[a,b]\times[c,d]$の場合を示します。

&&&prf
左辺を$e_1\partial_xF$と$e_2\partial_yF$の項に分ける。$e_1$は定数なので積分の外に出せ、$x$についての積分に1次元の基本定理を使うと

$$
\int_c^d\!\!\int_a^b e_1\partial_xF\,dx\,dy=\int_c^d\bigl(e_1F(b,y)-e_1F(a,y)\bigr)\,dy
$$

となる。右辺の第1項は右の辺（外向きの法線$\boldsymbol n=e_1$）、第2項は左の辺（$\boldsymbol n=-e_1$）での$\boldsymbol nF$の積分である。同様に$e_2\partial_yF$の項は、上の辺（$\boldsymbol n=e_2$）と下の辺（$\boldsymbol n=-e_2$）での$\boldsymbol nF$の積分を与える。4つの辺を合わせると$\oint\boldsymbol nF\,ds$になる。
&&&

$x$方向の微分は左右の辺に、$y$方向の微分は上下の辺に対応し、それぞれの係数$e_1,e_2$がそのまま辺の法線になっています。係数を生成元に取ったことで、1次元の基本定理を方向ごとに当てはめた結果が、法線ベクトルを掛ける形にまとまります。

## 一般の領域

領域を小さな長方形に分割して、それぞれで定理を当てはめて足し合わせます。隣り合う2つの長方形が共有する辺では、一方から見た外向きの法線と他方から見た外向きの法線が逆向きなので、その辺での$\boldsymbol nF$の積分は打ち消し合います。残るのは元の領域の境界だけです。

これで、有限個の長方形の和で表せる領域について定理が示されます。境界が曲線の場合は長方形に有限分割できないため、近似した領域をとって極限を取る必要があり、体積積分だけでなく、法線を掛けた境界積分の収束も示さなければなりません。この部分の詳細には立ち入らず、一般の領域については証明の概略にとどめます。

&&&ex 円板と位置ベクトル
単位円板$M$で$F=\boldsymbol x=xe_1+ye_2$とします。$D\boldsymbol x=e_1e_1+e_2e_2=2$なので、左辺は$2\times\pi=2\pi$です。境界の単位円上では外向きの法線が$\boldsymbol n=\boldsymbol x$なので、$\boldsymbol n\boldsymbol x=\boldsymbol x^2=|\boldsymbol x|^2=1$となり、右辺は円周の長さ$2\pi$です。
&&&

# 有向線素

境界は、領域を左手に見る向きにたどります。外側の境界では反時計回り、穴の境界では時計回りです。

&&&def 有向線素
境界をたどるときの接線方向の変位をベクトル

$$
d\boldsymbol x=e_1\,dx+e_2\,dy
$$

で表し、**有向線素**と呼びます。長さは$|d\boldsymbol x|=ds$で、向きが境界のたどる向きを表します。
&&&

外向きの法線は接線を時計回りに$90^\circ$回したものです。領域を左手に見てたどると、進行方向の右手が外側だからです。この回転は$I$を左から掛ける操作で表せます。

&&&fml 法線と有向線素
$$
\boldsymbol n\,ds=I\,d\boldsymbol x
$$
&&&

&&&prf
$Ie_1=e_1e_2e_1=-e_2$、$Ie_2=e_1e_2e_2=e_1$より

$$
I\,d\boldsymbol x=e_1\,dy-e_2\,dx
$$

である。接線方向$(dx,dy)$を時計回りに$90^\circ$回すと$(dy,-dx)$であり、これが$\boldsymbol n\,ds$である。穴の境界でも、領域を左手に見る向きにたどる限り同じである。
&&&

[基本定理](#thm-2d)に代入し、$I^{-1}=-I$を左から掛けると、境界積分を有向線素で書いた形が得られます。

&&&fml 有向線素による形
$$
\oint_{\partial M}d\boldsymbol x\,F=-I\int_M DF\,dA
$$
&&&

測度$d\boldsymbol x$がベクトルであり、被積分関数と代数の積で結ばれています。積が非可換なので、$d\boldsymbol x$を$F$のどちら側に置くかで積分の値が変わります。

# グリーンの定理

ベクトル値の関数$F=ae_1+be_2$に[基本定理](#thm-2d)を当てはめます。左辺は[[7shi-cla1]]で見たとおり、発散と回転に分かれます。

$$
DF=(a_x+b_y)+(b_x-a_y)I
$$

右辺では、$\boldsymbol n\,ds=e_1\,dy-e_2\,dx$を使って

$$
\begin{aligned}
\boldsymbol nF\,ds
&=(e_1\,dy-e_2\,dx)(ae_1+be_2) \\
&=(a\,dy-b\,dx)+(a\,dx+b\,dy)I
\end{aligned}
$$

となります。両辺のスカラー部と$I$の係数を比べると、2つの式が得られます。

&&&fml グリーンの定理の2つの形
$$
\begin{aligned}
\int_M(a_x+b_y)\,dA&=\oint_{\partial M}(a\,dy-b\,dx) \\
\int_M(b_x-a_y)\,dA&=\oint_{\partial M}(a\,dx+b\,dy)
\end{aligned}
$$
&&&

第1式は、発散の積分が境界を横切る流束に等しいことを表す平面の発散定理です。第2式は、回転の積分が境界に沿った循環に等しいことを表し、通常グリーンの定理と呼ばれる形です。ベクトル解析ではこの2つを別々の定理として扱いますが、ここでは1つの式のスカラー部と2ベクトル部です。微分の側で$DF$が発散と回転を1つにまとめていたこと、および法線を掛ける積$\boldsymbol nF$が内積と外積を1つにまとめていることが、積分の側にそのまま引き継がれています。

内積とベクトル積を別々の演算として扱うと[[7shi-hist]]、対応する積分公式も別々の式になります。幾何積を使うと、それらが1つの式にまとまります。

値をスカラー$\varphi$に取る場合、[基本定理](#thm-2d)は$\int_M D\varphi\,dA=\oint\boldsymbol n\varphi\,ds$となり、勾配の積分が境界上の法線の重み付き積分に等しいという式になります。

# 3次元の基本定理

$\operatorname{Cl}_{3,0}(\mathbb R)$の生成元を$e_1,e_2,e_3$（$e_k^2=1$、異なるもの同士は反交換）とし、ディラック作用素を

$$
D=e_1\partial_1+e_2\partial_2+e_3\partial_3
$$

とします。2次元と同じ理由で$D^2=\Delta$が成り立ちます。

## 立体

空間の有界な領域$V$の境界$\partial V$は区分的に滑らかな閉曲面とし、外向きの単位法線ベクトルを$\boldsymbol n$、面積の要素を$dS$と書きます。関数$F$は、$V$とその境界を含む開集合上で定義され、$\operatorname{Cl}_{3,0}(\mathbb R)$に値を取る$C^1$級関数とします。

&&&thm 3次元の基本定理 [thm-3d]
$$
\int_V DF\,dV=\oint_{\partial V}\boldsymbol nF\,dS
$$
&&&

証明は2次元と同じです。直方体では、$e_k\partial_kF$の項が$x_k$方向に向かい合う2つの面での$\pm e_kF$の積分を与えます。一般の領域へは、2次元と同じく分割と近似によって広げます（概略）。

&&&ex 球と位置ベクトル
単位球$V$で$F=\boldsymbol x$とすると、$D\boldsymbol x=e_1e_1+e_2e_2+e_3e_3=3$なので、左辺は$3\times\frac43\pi=4\pi$です。境界の単位球面上では$\boldsymbol n\boldsymbol x=\boldsymbol x^2=1$なので、右辺は球面の面積$4\pi$です。
&&&

ベクトル値の関数$F$に当てはめます。3次元では擬スカラー$I=e_1e_2e_3$がすべての元と可換で$I^2=-1$を満たし、2ベクトルとベクトル積が次のように対応します。

$$
\boldsymbol a\wedge\boldsymbol b=I(\boldsymbol a\times\boldsymbol b)
$$

基底で確かめると、$e_1\wedge e_2=e_1e_2=e_1e_2e_3e_3=Ie_3$で、$e_1\times e_2=e_3$に一致します。$D$についても同様に、ベクトル値の$F$に対して

$$
DF=D\cdot F+D\wedge F=\nabla\cdot F+I(\nabla\times F),\qquad
\boldsymbol nF=\boldsymbol n\cdot F+I(\boldsymbol n\times F)
$$

となります。[基本定理](#thm-3d)のスカラー部と2ベクトル部を比べると、2つの式が得られます。

&&&fml 発散定理と回転の体積分
$$
\begin{aligned}
\int_V\nabla\cdot F\,dV&=\oint_{\partial V}\boldsymbol n\cdot F\,dS \\
\int_V\nabla\times F\,dV&=\oint_{\partial V}\boldsymbol n\times F\,dS
\end{aligned}
$$
&&&

第1式はガウスの発散定理です。第2式は回転の体積分を境界の積分で表すもので、発散定理ほど知られてはいませんが、同じ式の2ベクトル部として同時に得られます。

## 曲面

空間内の曲面$S$とその境界の曲線$\partial S$を考えます。曲面では、立体の境界のように「立体の外側」を指す法線は自動的には決まりません。そこで曲面の向きを選び、それを2ベクトルで表します。

曲面は、1枚の正則な曲面片とします。すなわち、正方形$0\le u,v\le1$から重複なく表示する$C^2$級のパラメーター表示$\boldsymbol r(u,v)$を持ち、$\boldsymbol r_u=\partial\boldsymbol r/\partial u$、$\boldsymbol r_v=\partial\boldsymbol r/\partial v$が$\boldsymbol r_u\wedge\boldsymbol r_v\ne0$を満たすものとします。関数$F$は、曲面の近傍で定義された$C^1$級関数とします。

&&&def 有向面素
$\boldsymbol r_u,\boldsymbol r_v$が張る平行四辺形の向きを持つ面を

$$
d\boldsymbol X=\boldsymbol r_u\wedge\boldsymbol r_v\,du\,dv
$$

と書き、**有向面素**と呼びます。境界は、パラメーターの正方形を反時計回りにたどる向きに合わせます。
&&&

2ベクトルを積分に組み込むために、ベクトル$\boldsymbol a$と2ベクトル$B$の内積を次で定めます。

&&&def ベクトルと2ベクトルの内積
$$
\boldsymbol a\cdot B=\frac12(\boldsymbol aB-B\boldsymbol a),\qquad
B\cdot\boldsymbol a=\frac12(B\boldsymbol a-\boldsymbol aB)=-\boldsymbol a\cdot B
$$
&&&

$B=\boldsymbol b\wedge\boldsymbol c$に対しては、基底で確かめられるとおり次のベクトルになります。

$$
\boldsymbol a\cdot(\boldsymbol b\wedge\boldsymbol c)=(\boldsymbol a\cdot\boldsymbol b)\,\boldsymbol c-(\boldsymbol a\cdot\boldsymbol c)\,\boldsymbol b
$$

$D$との内積は、係数の側だけで内積を取り、偏微分は$F$に作用させるものとします。

$$
(D\cdot d\boldsymbol X)F=\sum_k(e_k\cdot d\boldsymbol X)\,\partial_kF
$$

&&&thm 曲面の基本定理 [thm-surf]
$$
\oint_{\partial S}d\boldsymbol x\,F=\int_S(D\cdot d\boldsymbol X)F
$$
&&&

&&&prf
内積の公式より$e_k\cdot(\boldsymbol r_u\wedge\boldsymbol r_v)=(e_k\cdot\boldsymbol r_u)\boldsymbol r_v-(e_k\cdot\boldsymbol r_v)\boldsymbol r_u$であり、連鎖律$\partial_uF=\sum_k(e_k\cdot\boldsymbol r_u)\partial_kF$を使うと

$$
\sum_k\bigl(e_k\cdot(\boldsymbol r_u\wedge\boldsymbol r_v)\bigr)\partial_kF
=\boldsymbol r_v\,\partial_uF-\boldsymbol r_u\,\partial_vF
$$

となる。$\partial_u\boldsymbol r_v=\partial_v\boldsymbol r_u$なので、これは$\partial_u(\boldsymbol r_vF)-\partial_v(\boldsymbol r_uF)$に等しい。正方形上で積分し、それぞれの項に1次元の基本定理を使うと

$$
\int_S(D\cdot d\boldsymbol X)F
=\int_0^1\bigl[\boldsymbol r_vF\bigr]_{u=0}^{u=1}dv-\int_0^1\bigl[\boldsymbol r_uF\bigr]_{v=0}^{v=1}du
$$

となる。一方、境界を反時計回りにたどると、有向線素は下の辺（$v=0$、$u$が増える向き）で$\boldsymbol r_u\,du$、右の辺（$u=1$、$v$が増える向き）で$\boldsymbol r_v\,dv$、上の辺と左の辺ではそれぞれ逆向きの$-\boldsymbol r_u\,du$、$-\boldsymbol r_v\,dv$である。4辺を合わせると上の右辺に一致する。
&&&

証明に使ったのは1変数の基本定理だけで、曲面が空間の中でどう曲がっているかは$\boldsymbol r_u,\boldsymbol r_v$を通じて自動的に処理されています。$(D\cdot d\boldsymbol X)$に現れるのは曲面に沿った方向の微分だけで、曲面から外れる方向の微分は内積によって落ちます。したがって結果は、$F$を曲面の外へどう延長したかにはよりません。

1枚の曲面片で表示できない向き付けられた曲面は、曲面片を貼り合わせて扱います。隣り合う曲面片が共有する辺は互いに逆向きにたどられるので、その辺での積分は打ち消し合います。

&&&rem 内積の順序
$\boldsymbol a\cdot B=-B\cdot\boldsymbol a$なので、$(d\boldsymbol X\cdot D)F$と書くと符号が反転します。本記事では$D\cdot d\boldsymbol X$の順に揃えます。
&&&

平面そのものを曲面と見ると、$\boldsymbol r(u,v)=ue_1+ve_2$で$d\boldsymbol X=e_1e_2\,dA$です。$e_3$の項は内積で消えて

$$
D\cdot(e_1e_2)=e_2\partial_1-e_1\partial_2=-e_1e_2(e_1\partial_1+e_2\partial_2)
$$

となり、右辺は2次元の擬スカラーと作用素による$-ID$です。したがって[曲面の基本定理](#thm-surf)は、2次元の有向線素による形$\oint d\boldsymbol x\,F=-I\int DF\,dA$に一致します。

## ストークスの定理

3次元では有向面素を単位法線$\boldsymbol n$で$d\boldsymbol X=I\boldsymbol n\,dA$と書けます（$\boldsymbol r_u\wedge\boldsymbol r_v=I(\boldsymbol r_u\times\boldsymbol r_v)$）。$I$はすべての元と可換なので

$$
\boldsymbol a\cdot(I\boldsymbol n)=\frac I2(\boldsymbol a\boldsymbol n-\boldsymbol n\boldsymbol a)=I(\boldsymbol a\wedge\boldsymbol n)=I^2(\boldsymbol a\times\boldsymbol n)=\boldsymbol n\times\boldsymbol a
$$

となり、$(D\cdot d\boldsymbol X)F=\bigl((\boldsymbol n\times\nabla)F\bigr)dA$です。ベクトル値の$F$に対して[曲面の基本定理](#thm-surf)のスカラー部を取ると、左辺は$d\boldsymbol x\cdot F$、右辺は$(\boldsymbol n\times\nabla)\cdot F=\boldsymbol n\cdot(\nabla\times F)$です。

&&&fml ストークスの定理
$$
\oint_{\partial S}F\cdot d\boldsymbol x=\int_S(\nabla\times F)\cdot\boldsymbol n\,dA
$$
&&&

# 一般次元と両側形式

直方体での証明は次元によらないため、$n$次元でも同じ形が成り立ちます。$\operatorname{Cl}_{n,0}(\mathbb R)$の生成元$e_1,\dots,e_n$で$D=\sum_ke_k\partial_k$とし、区分的に滑らかな境界を持つ$n$次元の有界な領域$M$の外向きの単位法線を$\boldsymbol n$、境界の面積要素を$dS$とします。関数$F$は、$M$とその境界を含む開集合上で定義され、$\operatorname{Cl}_{n,0}(\mathbb R)$に値を取る$C^1$級関数とします。

&&&thm $n$次元の基本定理
$$
\int_M DF\,dV=\oint_{\partial M}\boldsymbol nF\,dS
$$
&&&

積が非可換なので、関数を左右から挟む形も用意しておきます。$G$も$F$と同じ条件を満たす関数とし、右からの作用を$GD=\sum_k(\partial_kG)e_k$とします。

&&&thm 両側形式
$$
\int_M\bigl((GD)F+G(DF)\bigr)dV=\oint_{\partial M}G\,\boldsymbol nF\,dS
$$
&&&

&&&prf
積の微分法則より$\partial_k(Ge_kF)=(\partial_kG)e_kF+Ge_k\,\partial_kF$であり、$k$について和を取ると被積分関数に一致する。直方体では、$\partial_k(Ge_kF)$の積分が$x_k$方向に向かい合う2つの面での$\pm Ge_kF$の積分を与え、$\pm e_k$がその面の外向きの法線である。一般の領域へは、元の形と同じく分割と近似によって広げる（概略）。
&&&

$G=1$とすれば元の形に戻ります。法線$\boldsymbol n$は$G$と$F$の間に挟まれ、$D$が$G$に右から、$F$に左から作用した結果と対応しています。可換な場合は$\boldsymbol n$の位置を気にする必要はありませんが、非可換な代数に値を取る関数ではこの位置が意味を持ちます。

曲面の場合と同様に、$n$次元空間の中の一般の次元の曲面に対しても基本定理が成り立ちますが、本記事では扱いません。

&&&rem 微分形式との比較
微分形式を使うと、ストークスの定理$\int_Md\omega=\int_{\partial M}\omega$が各種の積分定理を1つにまとめます。通常の提示では、発散定理と回転の積分定理は次数の異なる形式$\omega$に対する別々の適用として現れます。微分形式でも次数の異なる形式をまとめて扱うことはできますが、本記事の基本定理では、幾何積によって、1つのベクトル値関数$F$から発散と回転、およびそれぞれに対応する境界項が一度に得られます。
&&&

# まとめ

1次元の基本定理は、境界の各点に外向きの符号$n=\pm1$を割り当てると$\int F'dx=\sum nF$と書けます。$\operatorname{Cl}_{n,0}(\mathbb R)$の生成元を係数とするディラック作用素$D$を使うと、この形がそのまま多次元に持ち上がり、$\int_M DF\,dV=\oint_{\partial M}\boldsymbol nF\,dS$となります。長方形（直方体）では、方向ごとに1次元の基本定理を当てはめるだけで証明でき、係数$e_k$がそのまま境界の法線になります。一般の領域へは分割と近似で広げます。

| 次元 | 形 | $F$ | 成分 | 定理 |
|---|---|---|---|---|
| 2 | $\int DF\,dA=\oint\boldsymbol nF\,ds$ | ベクトル | スカラー部 | 平面の発散定理 |
| 2 | 同上 | ベクトル | 2ベクトル部 | グリーンの定理 |
| 3 | $\int DF\,dV=\oint\boldsymbol nF\,dS$ | ベクトル | スカラー部 | ガウスの発散定理 |
| 3 | 同上 | ベクトル | 2ベクトル部 | 回転の体積分 |
| 3 | $\oint d\boldsymbol x\,F=\int(D\cdot d\boldsymbol X)F$ | ベクトル | スカラー部 | ストークスの定理 |

平面では$\boldsymbol n\,ds=I\,d\boldsymbol x$により、境界積分を有向線素$d\boldsymbol x$で書いた形$\oint d\boldsymbol x\,F=-I\int DF\,dA$にもなります。空間内の曲面では外向きの法線の代わりに有向面素$d\boldsymbol X$を使い、$\oint d\boldsymbol x\,F=\int(D\cdot d\boldsymbol X)F$となります。平面をそのまま曲面と見ると、この2つは一致します。

ベクトル解析で別々に扱われる積分定理は、領域の基本定理とその曲面版という2つの式のグレード成分として現れます。微分の側の$DF$が発散と回転をまとめていたことが、積分の側では法線との積$\boldsymbol nF$が内積と外積をまとめることに対応しています。非可換な代数に値を取る関数のために、$\boldsymbol n$を$G$と$F$で挟む両側形式も成り立ちます。
