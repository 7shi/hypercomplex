パウリ行列をエルミート行列の基底として導入し、四元数の行列表現との関係を確認します。[[7shi-bq]] パウリ行列をクロック行列とシフト行列の組と見なして一般の次数$n$に拡張したものが一般化パウリ行列で、$n=3$の場合はシルベスターの九元数に一致します。[[7shi-nonion]][[wiki-gpm]] 一般化パウリ行列を行列表現として持つ枠組みが一般化クリフォード代数です。[[wiki-gca]] クリフォード代数がテンソル積によって拡張できたのと同様に[[7shi-qt]]、一般化クリフォード代数でも生成元を2個増やすことが行列環とのテンソル積に対応することを示します。生成元が奇数個の場合は、擬スカラーに代わる中心元によって行列環の直和へ分解されます。生成元が1個の場合はフルーリーの多重複素数の複素化に相当します。

&&& 凡例
$\mathbb F$を成分とする$n$次の全行列環$M_n(\mathbb F)$を$\mathbb F(n)$と表記し（例：$\mathbb C(2) \cong M_2(\mathbb C)$）、代数$A$の$n$個の直和を$nA$と略記します。
&&&

# パウリ行列

一般の2×2複素行列と、その成分の複素共役を取って転置した行列（エルミート共役）を考えます。

$$
A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}, \quad
A^\dagger = \begin{pmatrix} a^* & c^* \\ b^* & d^* \end{pmatrix}
$$

エルミート共役で変わらない行列（$A^\dagger=A$）を**エルミート行列**と呼びます。成分を比較すると、条件は対角成分$a,d$が実数であること、非対角成分が互いに複素共役であること（$b=c^*$）です。実パラメーターで数えると計4つです。

このエルミート行列を、トレースを担う部分（単位行列の実数倍）と、トレースが$0$の部分に分解します。

$$
w=\dfrac{a+d}2,\quad z=\dfrac{a-d}2
$$
とおくと
$$
\begin{pmatrix} a & c^* \\ c & d \end{pmatrix}
= w\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}
+ \begin{pmatrix} z & c^* \\ c & -z \end{pmatrix}
$$

$w,z$は実数なので、どちらの部分もエルミート性を保ちます。$c=x+iy$とおいて成分ごとに行列を取り出します。

$$
\begin{pmatrix} z & x-iy \\ x+iy & -z \end{pmatrix}
= x\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
+ y\begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}
+ z\begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$

ここに現れる3つの行列が**パウリ行列**です。

&&&def パウリ行列
$$
\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
\sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad
\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$
&&&

つまりパウリ行列は、トレースが$0$のエルミート行列の実基底であり、単位行列$I$を加えるとエルミート行列全体の実基底となります。

パウリ行列は次の関係式を満たします（直接計算で確認できます）。

$$
\sigma_x^2 = \sigma_y^2 = \sigma_z^2 = I, \quad
\sigma_x\sigma_y\sigma_z = iI
$$
$$
\sigma_x\sigma_y = -\sigma_y\sigma_x = i\sigma_z, \quad
\sigma_y\sigma_z = -\sigma_z\sigma_y = i\sigma_x, \quad
\sigma_z\sigma_x = -\sigma_x\sigma_z = i\sigma_y
$$

## クリフォード代数の表現として

パウリ行列は互いに反交換して2乗が$+1$、つまり符号数$(3,0)$のクリフォード代数の生成元の関係式を満たします。実係数の積和で生成される代数は$1,\sigma_x,\sigma_y,\sigma_z,\sigma_y\sigma_z,\sigma_z\sigma_x,\sigma_x\sigma_y,\sigma_x\sigma_y\sigma_z$を基底とする実8次元で、前々回確認した分類の

$$
\operatorname{Cl}_{3,0}(\mathbb R) \cong \mathbb C(2)
$$

の行列表現そのものです。[[7shi-clif1]]

一方、係数を複素数に広げると、$\sigma_y = i\sigma_x\sigma_z$より$\sigma_y$は$\sigma_x,\sigma_z$の積から作れるため、生成元は$\sigma_x,\sigma_z$の2個で足ります。これは複素クリフォード代数として

$$
\operatorname{Cl}_2(\mathbb C) \cong \mathbb C(2)
$$

に対応します。

## 双四元数との関係

係数を複素数に広げると生成元2個で足りるという見方は、前々回の分類と対応しています。2個の生成元からなる実クリフォード代数$\operatorname{Cl}_{0,2}(\mathbb R)$は四元数体$\mathbb H$と同型であり、複素クリフォード代数は実クリフォード代数の複素化として得られました。[[7shi-clif1]] これらを合わせると、次の同型が得られます。

$$
\operatorname{Cl}_2(\mathbb C)
\cong \mathbb C \otimes \operatorname{Cl}_{0,2}(\mathbb R)
\cong \mathbb C \otimes \mathbb H
$$

この$\mathbb C \otimes \mathbb H$は**双四元数 (biquaternion)** と呼ばれます。[[7shi-bq]] 複素数と四元数の虚数単位を区別するため、複素数の虚数単位を$h$と書き、テンソル積の記号は省略します。$h$は係数側の元なので四元数の$i,j,k$と可換で、$h^2=-1$を満たします。一般の双四元数は実8次元の元

$$
a_0+a_1h+a_2i+a_3hi+a_4j+a_5hj+a_6k+a_7hk
\quad (a_0,\dots,a_7 \in \mathbb R)
$$

として表されます。

パウリ行列は2乗が$+I$でした。双四元数の側でこれに対応するのは、2乗が$+1$となる元です。

$$
(hi)^2 = h^2i^2 = (-1)(-1) = 1
$$

$hj,hk$も同様なので、次の対応が得られます。

&&&fml パウリ行列と双四元数の対応
$$
\sigma_x \cong hi, \quad \sigma_y \cong hj, \quad \sigma_z \cong hk
$$
&&&

冒頭で確認した関係式は、この対応によって成分計算なしに再現できます。

$$
\begin{aligned}
\sigma_x\sigma_y &\cong (hi)(hj) = h^2(ij) = -k = h(hk) \cong i\sigma_z \\
\sigma_x\sigma_y\sigma_z &\cong (hi)(hj)(hk) = h^3(ijk) = (-h)(-1) = h \cong iI
\end{aligned}
$$

2行目は、複素係数の$i$が双四元数では$h$に対応することを示しています。

双四元数の基底のうち$1,i,j,k$が張る実部分代数は、四元数$\mathbb H$そのものです。つまり双四元数の中には四元数が含まれています。これをパウリ行列の側から取り出します。

## 四元数との関係

四元数の虚数単位$i,j,k$は、対応$\sigma_x \cong hi,\ \sigma_y \cong hj,\ \sigma_z \cong hk$から$h$を取り除けば得られます。$h^2=-1$より、$-h$を掛ければ$h$が消えます。

$$
-h(hi) = -h^2i = i, \quad -h(hj) = j, \quad -h(hk) = k
$$

複素係数の$i$が双四元数の$h$に対応していたので、$-h$を掛けることはパウリ行列を$-i$倍することに対応します。

&&&fml 四元数の行列表現
$$
i \cong -i\sigma_x, \quad j \cong -i\sigma_y, \quad k \cong -i\sigma_z
$$
&&&

行列の側で直接確認すると、2乗は$-I$に反転し

$$
(-i\sigma_x)^2 = (-i\sigma_y)^2 = (-i\sigma_z)^2 = -I
$$

積は

$$
(-i\sigma_x)(-i\sigma_y) = -\sigma_x\sigma_y = -i\sigma_z, \quad
(-i\sigma_x)(-i\sigma_y)(-i\sigma_z) = (-i)^3(iI) = -I
$$

となり、四元数の関係式$i^2=j^2=k^2=ijk=-1,\ ij=k$がそのまま再現されます。

&&&rem
2乗して$+I$になるエルミートなパウリ行列と、2乗して$-1$になる四元数の虚数単位は、このように$\pm h$倍（$\pm i$倍）で行き来します。詳細は以前の記事を参照してください。[[7shi-bq]]
&&&

# クロック行列とシフト行列

複素係数では$\sigma_x,\sigma_z$の2個の生成元で足りるという見方を出発点に、一般化を進めます。まず$\sigma_z,\sigma_x$を、成分の形ではなく基底ベクトル$e_0,e_1$への作用で捉え直します。

- $\sigma_z$は$e_0,e_1$を$\pm1$倍します。$\pm1$は1の平方根です。
- $\sigma_x$は$e_0,e_1$を入れ替えます。すなわち巡回シフトです。

「1の平方根を対角に並べる」「基底を巡回シフトする」という操作は、そのまま一般の次数$n$に拡張できます。$\omega=e^{2\pi i/n}$を1の原始$n$乗根として、次の2つの行列を定義します。

&&&def クロック行列とシフト行列
$$
Z = \begin{pmatrix}
1 & & & \\
& \omega & & \\
& & \ddots & \\
& & & \omega^{n-1}
\end{pmatrix}, \quad
X = \begin{pmatrix}
0 & & & 1 \\
1 & 0 & & \\
& \ddots & \ddots & \\
& & 1 & 0
\end{pmatrix}
$$
&&&

基底ベクトル$e_0,\dots,e_{n-1}$への作用で書けば、添え字を$n$を法として

$$
Ze_j = \omega^j e_j, \quad Xe_j = e_{j+1}
$$

です。$Z$が「時計の文字盤」のように1の$n$乗根を順に並べることから**クロック行列**、$X$が基底を巡回させることから**シフト行列**と呼ばれます。[[wiki-gpm]]

&&&fml クロック行列とシフト行列の関係式
$$
X^n = Z^n = I, \quad ZX = \omega XZ
$$
&&&

&&&prf
$X^ne_j=e_{j+n}=e_j$、$Z^ne_j=\omega^{jn}e_j=e_j$より$X^n=Z^n=I$。交換関係は

$$
ZXe_j = Ze_{j+1} = \omega^{j+1}e_{j+1}, \quad
XZe_j = \omega^jXe_j = \omega^je_{j+1}
$$

より$ZX=\omega XZ$となる。
&&&

$n=2$では$X=\sigma_x,\ Z=\sigma_z,\ \omega=-1$で、関係式は$\sigma_z\sigma_x=-\sigma_x\sigma_z$、すなわち反交換関係に帰着します。パウリ行列の反交換関係は、$\omega$によるねじれた交換関係の$n=2$の場合です。

&&&rem
$n=2$だけを見ていると、反交換関係の$-1$は「掛ける順序を入れ替えると符号が反転する」という符号の現象に見えます。しかし$n\ge3$では$\omega$は複素数となり、入れ替えで付くのは符号ではなく単位円上の位相です。この視点から振り返ると、$n=2$の$-1$も1の原始2乗根$e^{2\pi i/2}$、すなわち位相の特別な場合であり、反交換関係は位相のねじれがたまたま実数に収まった姿だったと分かります。
&&&

## 一般化パウリ行列

クロック行列とシフト行列の単項式を**一般化パウリ行列**と呼びます。

&&&def 一般化パウリ行列
$$
X^aZ^b \quad (0 \le a,b < n)
$$
&&&

計算のために、単項式の入れ替えと冪の公式を用意します。$ZX=\omega XZ$を繰り返し使うと、$b$個の$Z$が$a$個の$X$を追い越すたびに$\omega$が1つ付くので

$$
Z^bX^a = \omega^{ab}X^aZ^b
$$

となります。これを使って$k$個のブロック$X^aZ^b$を並べ替えると、冪の公式が得られます。

&&&fml 単項式の冪
$$
(X^aZ^b)^k = \omega^{ab\,k(k-1)/2}X^{ka}Z^{kb}
$$
特に$k=n$では$(X^aZ^b)^n = \omega^{ab\,n(n-1)/2}I$となる。
&&&

$n^2$個の一般化パウリ行列は、行列環$\mathbb C(n)$の基底をなします。

&&&thm 一般化パウリ行列は$\mathbb C(n)$の基底
$X^aZ^b\ (0 \le a,b < n)$は互いにトレース内積で直交し、$\mathbb C(n)$の複素基底をなす。
&&&

&&&prf
$X,Z$はユニタリなので$(X^aZ^b)^\dagger = Z^{-b}X^{-a}$であり、トレース内積は

$$
\operatorname{Tr}\left[(X^aZ^b)^\dagger(X^{a'}Z^{b'})\right]
= \omega^{-b(a'-a)}\operatorname{Tr}\left(X^{a'-a}Z^{b'-b}\right)
$$

となる。ここで$X^cZ^d$のトレースは、$c\not\equiv0$なら対角成分がすべて$0$、$c\equiv0,\ d\not\equiv0$なら$\operatorname{Tr}Z^d=\sum_j\omega^{jd}=0$（1の冪根の和）となる。よって内積は$(a,b)=(a',b')$のとき$n$、それ以外は$0$となり、$n^2$個の単項式は直交する。$\dim M_n(\mathbb C)=n^2$と一致するため基底となる。
&&&

$n=2$の一般化パウリ行列は$I,\ \sigma_x,\ \sigma_z,\ \sigma_x\sigma_z=-i\sigma_y$で、位相（1の冪根や$\pm i$倍）を除けば単位行列とパウリ行列に一致します。

## シルベスターの九元数との対応

$n=3$の場合を、シルベスターの九元数と比較します。[[7shi-nonion]] 九元数の記事では$\rho=e^{2\pi i/3}$として、関係式

$$
vu = \rho uv, \quad u^3 = v^3 = 1
$$

を満たす生成元として、次の行列を確認しました。

$$
u = \begin{pmatrix} 0 & 0 & 1 \\ \rho & 0 & 0 \\ 0 & \rho^2 & 0 \end{pmatrix}, \quad
v = \begin{pmatrix} 0 & 0 & 1 \\ \rho^2 & 0 & 0 \\ 0 & \rho & 0 \end{pmatrix}
$$

これらはクロック行列とシフト行列の単項式です。成分計算で直接確認できます（$\rho=\omega$）。

$$
u = ZX = \rho XZ, \quad v = Z^2X = \rho^2XZ^2
$$

つまり$u,v$は位相倍を除いて一般化パウリ行列$XZ,\ XZ^2$であり、関係式もクロック行列とシフト行列の関係式から従います。

&&&prf
$vu$と$uv$を$Z^bX^a=\omega^{ab}X^aZ^b$で整理する。

$$
vu = Z^2(XZ)X = \omega^{-1}Z^3X^2 = \omega^{-1}X^2, \quad
uv = Z(XZ^2)X = \omega^{-2}Z^3X^2 = \omega^{-2}X^2
$$

よって$vu=\omega uv=\rho uv$となる。冪は単項式の冪の公式より

$$
u^3 = \rho^3(XZ)^3 = \omega^{3\cdot2/2}X^3Z^3 = I, \quad
v^3 = \rho^6(XZ^2)^3 = \omega^{2\cdot3\cdot2/2}X^3Z^6 = I
$$

となる。
&&&

さらに一般の単項式$u^av^b$は、位相倍を除いて$X^{a+b}Z^{a+2b}$に一致します。指数の対応$(a,b)\mapsto(a+b,\ a+2b)$は$3$を法として可逆なので、九元数の$9$個の基底は一般化パウリ行列$X^aZ^b$の$9$個と、位相倍を除いて1対1に対応します。すなわち九元数の実体は$\mathbb C(3)$です。

&&&rem
九元数の記事に現れた「$\rho^2uv$がシフト行列になる」という関係も$\rho^2uv=\rho^2\cdot\omega^{-2}X^2=X^2$と確認できます。$X^2=X^{-1}$は逆向きの巡回シフトで、九元数の記事のシフト行列はこの向きで書かれています。
&&&

四元数がパウリ行列（$\mathbb C(2)$の基底）と結びついたのと同じ形で、九元数は一般化パウリ行列（$\mathbb C(3)$の基底）と結びつきます。

# 一般化クリフォード代数

パウリ行列を複素係数で見たとき、生成元$\sigma_x,\sigma_z$は「2乗が$1$」「$-1$倍のねじれで交換」という関係式を満たしました。クロック行列とシフト行列はこれを「$n$乗が$1$」「$\omega$倍のねじれで交換」に置き換えたものです。この関係式を公理として抽象化します。

&&&def 一般化クリフォード代数
$\omega=e^{2\pi i/n}$とする。$m$個の生成元$e_1,\dots,e_m$が関係式

$$
e_i^n = 1, \quad e_je_i = \omega e_ie_j \quad (i<j)
$$

を満たすとき、生成される複素代数を**一般化クリフォード代数**と呼び、$\operatorname{Cl}^{(n)}_m(\mathbb C)$と表記する。[[wiki-gca]]
&&&

関係式によって生成元の積は順序を整理できるため、基底は単項式

$$
e_1^{a_1}e_2^{a_2}\cdots e_m^{a_m} \quad (0 \le a_i < n)
$$

の$n^m$個で、次元は$n^m$です。

&&&rem $n=2$は複素クリフォード代数
$n=2$では$\omega=-1$で、関係式は$e_i^2=1$と反交換関係となり、生成元$m$個の複素クリフォード代数に一致します。

$$
\operatorname{Cl}^{(2)}_m(\mathbb C) \cong \operatorname{Cl}_m(\mathbb C)
$$

複素数体上では生成元の2乗の符号が区別されない（$e$を$ie$に取り直せば反転できる）ため、$e_i^2=1$として一般性を失いません。[[7shi-clif1]]
&&&

&&&rem 生成元が1個の場合
$m=1$では関係式は$e^n=1$だけで、可換な代数$\mathbb C[e]/(e^n-1)$となります。$n=2$では、分解型複素数と同じ冪等元$\dfrac{1\pm e}2$による直和分解$2\mathbb C$が起こります。[[7shi-clif2]] 一般の$n$でも、冪等元

$$
\varepsilon_k = \frac1n\sum_{j=0}^{n-1}\omega^{-jk}e^j \quad (0 \le k < n)
$$

によって$n\mathbb C$へ直和分解されます。

この$m=1$の代数は、$e^n=-1$を満たす生成元による実$n$次元の可換代数**フルーリーの多重複素数**の複素化と同型です。[[wiki-fmc]] 複素係数では$e$を1の$2n$乗根倍で取り直せば$e^n=-1$を$e^n=1$に移せるためです。フルーリーの多重複素数は$n=2$では複素数$\mathbb C$そのもので、$n=3$では三元数（実3次元の多元数）の一例を与えます。複素数を重ねたテンソル積として構成されるセグレの多重複素数[[7shi-qt]]とは別の体系です。
&&&

## 行列表現

クロック行列とシフト行列は、$e_1=X,\ e_2=Z$として2生成元の関係式を満たします。

$$
X^n = Z^n = I, \quad ZX = \omega XZ
$$

単項式$X^aZ^b$が$\mathbb C(n)$の基底をなすことから、この表現は全行列環に一致します。

$$
\operatorname{Cl}^{(n)}_2(\mathbb C) \cong \mathbb C(n)
$$

$n=2$では$\operatorname{Cl}_2(\mathbb C)\cong\mathbb C(2)$（パウリ行列）、$n=3$では九元数$\cong\mathbb C(3)$です。クリフォード代数の行列表現がパウリ行列だったのと同じ意味で、一般化クリフォード代数の行列表現が一般化パウリ行列です。

# テンソル積による拡張

クリフォード代数では、$\mathbb H,\mathbb H'$とのテンソル積によって生成元を2個ずつ増やせました。[[7shi-qt]] 一般化クリフォード代数でも、生成元を2個増やすことは$\mathbb C(n)$とのテンソル積に対応します。

## 生成元の構成

$\operatorname{Cl}^{(n)}_m(\mathbb C)$の生成元$e_1,\dots,e_m$に対して、$\mathbb C(n)$とのテンソル積の中に$m+2$個の生成元を構成します。既存の生成元にはクロック行列の逆元$Z^{-1}=Z^{n-1}$を付け、新しい生成元はシフト行列$X$と、積$XZ$の位相倍から作ります。

$$
\begin{array}{ccc}

\begin{array}{c|ccc}
\operatorname{Cl}^{(n)}_m & e_1 & \cdots & e_m
\end{array} &

\xrightarrow{⊗\,\mathbb C(n)} &

\begin{array}{c|ccccc}
\operatorname{Cl}^{(n)}_m & e_1 & \cdots & e_m & 1 & 1 \\
& ⊗ & \cdots & ⊗ & ⊗ & ⊗ \\
\mathbb C(n) & Z^{-1} & \cdots & Z^{-1} & X & cXZ
\end{array} \\

\operatorname{Cl}^{(n)}_m(\mathbb C) & &
\operatorname{Cl}^{(n)}_{m+2}(\mathbb C)

\end{array}
$$

ここで$c$は$c^n(XZ)^n=I$、すなわち単項式の冪の公式から$c^n=\omega^{-n(n-1)/2}$を満たす位相で、$n$が奇数なら1の$n$乗根$c=\omega^{-(n-1)/2}$と取れます。

&&&rem
$n$が偶数のときは$c^n=-1$となるため、$c$は1の$2n$乗根から取ります。一様には$c=e^{-\pi i(n-1)/n}$と書けます。$n=2$では$c=\pm i$で、$\pm i\sigma_x\sigma_z=\pm\sigma_y$として3つ目のパウリ行列が現れます。
&&&

&&&prf
新しい生成元を$\tilde e_i = e_i \otimes Z^{-1}\ (i \le m)$、$\tilde e_{m+1} = 1 \otimes X$、$\tilde e_{m+2} = 1 \otimes cXZ$とおき、関係式を確認する。

**$n$乗**：$\tilde e_i^n = e_i^n \otimes Z^{-n} = 1 \otimes 1$、$\tilde e_{m+1}^n = 1 \otimes X^n = 1 \otimes 1$、$\tilde e_{m+2}^n = 1 \otimes c^n(XZ)^n = 1 \otimes 1$（最後は$c$の定義による）。

**既存どうし**（$i<j\le m$）：テンソル積の右因子は共通なので

$$
\tilde e_j\tilde e_i = e_je_i \otimes Z^{-2} = \omega\, e_ie_j \otimes Z^{-2} = \omega\,\tilde e_i\tilde e_j
$$

**新旧の間**（$i\le m$）：$ZX=\omega XZ$より$XZ^{-1}=\omega Z^{-1}X$、また$(XZ)Z^{-1}=X=\omega Z^{-1}(XZ)$となる。すなわち$N=X,\ cXZ$のどちらでも$NZ^{-1}=\omega Z^{-1}N$が成り立つので

$$
\tilde e_{m+a}\tilde e_i = e_i \otimes NZ^{-1} = \omega\, e_i \otimes Z^{-1}N = \omega\,\tilde e_i\tilde e_{m+a}
$$

**新しい生成元どうし**：$(cXZ)X = cX(ZX) = c\omega X^2Z = \omega X(cXZ)$より

$$
\tilde e_{m+2}\tilde e_{m+1} = \omega\,\tilde e_{m+1}\tilde e_{m+2}
$$

よって$\tilde e_1,\dots,\tilde e_{m+2}$は生成元$m+2$個の一般化クリフォード代数の関係式を満たす。
&&&

構成した生成元は、テンソル積の全体を生成します。

&&&prf
まず$\tilde e_{m+1}^a\tilde e_{m+2}^b$は位相倍を除いて$1 \otimes X^{a+b}Z^b$であり、指数の対応$(a,b)\mapsto(a+b,\ b)$は可逆なので、$1 \otimes X^aZ^b$がすべて（位相倍を除いて）得られる。次に

$$
\tilde e_1^{a_1}\cdots\tilde e_m^{a_m}
= e_1^{a_1}\cdots e_m^{a_m} \otimes Z^{-(a_1+\cdots+a_m)}
$$

に$1 \otimes X^aZ^b$を掛ければ、右因子は任意の$X^aZ^b$に調整できる。よって「$\operatorname{Cl}^{(n)}_m(\mathbb C)$の基底」と「$\mathbb C(n)$の基底」の任意の組み合わせが得られ、生成される代数はテンソル積全体と一致する。
&&&

## 拡張公式

生成元の構成と次元の一致$n^{m+2}=n^m \cdot n^2$から、拡張公式が得られます。

&&&fml 一般化クリフォード代数の$⊗\,\mathbb C(n)$による拡張
$$
\operatorname{Cl}^{(n)}_{m+2}(\mathbb C)
\cong \operatorname{Cl}^{(n)}_m(\mathbb C) \otimes \mathbb C(n)
$$
&&&

$\operatorname{Cl}^{(n)}_2(\mathbb C)\cong\mathbb C(n)$を出発点として公式を繰り返し適用すれば、行列環とのテンソル積の規則$\mathbb F(m) \otimes \mathbb F(n) \cong (\mathbb F \otimes \mathbb F)(mn)$より[[7shi-clif1]]、生成元が偶数個の場合の分類が完成します。

&&&thm 偶数個の生成元の行列表現
$$
\operatorname{Cl}^{(n)}_{2k}(\mathbb C)
\cong \mathbb C(n)^{\otimes k}
\cong \mathbb C(n^k)
$$
&&&

クリフォード代数では基底（生成元）を任意に増やせて、その行列表現は$\mathbb H,\mathbb H'$とのテンソル積で作れました。一般化クリフォード代数でも同じことが起きています。生成元を増やすことは、行列表現のレベルではクロネッカー積で実現され、代数のレベルでは$\mathbb C(n)$とのテンソル積による拡張になっています。

## 複素クリフォード代数との対応

$n=2$では拡張公式は

$$
\operatorname{Cl}_{m+2}(\mathbb C) \cong \operatorname{Cl}_m(\mathbb C) \otimes \mathbb C(2)
$$

となり、クリフォード代数の分類で確認した複素の2周期性に一致します。[[7shi-clif1]] このとき生成元の構成を$k=3$（生成元6個）まで繰り返すと、次の形になります（$c=-i$とした場合）。

$$
\begin{alignedat}{3}
\tilde e_1 &= \sigma_x \otimes \sigma_z \otimes \sigma_z, &\quad
\tilde e_3 &= 1 \otimes \sigma_x \otimes \sigma_z, &\quad
\tilde e_5 &= 1 \otimes 1 \otimes \sigma_x, \\
\tilde e_2 &= \sigma_z \otimes \sigma_z \otimes \sigma_z, &
\tilde e_4 &= 1 \otimes (-\sigma_y) \otimes \sigma_z, &
\tilde e_6 &= 1 \otimes 1 \otimes (-\sigma_y)
\end{alignedat}
$$

各生成元はパウリ行列1個と、それに続く$\sigma_z$の列（しっぽ）から成ります。

&&&rem
この形の生成元の構成は、ヨルダン＝ウィグナー変換として知られるものと（しっぽの向きを除いて）同じです。一般の$n$での構成は、その$\omega$版への一般化に相当します。
&&&

## 奇数個の生成元

クリフォード代数の分類では、擬スカラー（すべての生成元の積）が中心に入るかどうかが直和型の出現を決めました。[[7shi-clif1]] 一般化クリフォード代数では事情が異なります。単項式と$e_i$の交換を関係式で整理すると、$e_i$が$e_j^{c_j}$を追い越すたびに$j>i$なら$\omega^{c_j}$、$j<i$なら$\omega^{-c_j}$が付くので

$$
(e_1^{c_1}\cdots e_m^{c_m})\,e_i
= \omega^{d_i}\,e_i\,(e_1^{c_1}\cdots e_m^{c_m}), \quad
d_i = \sum_{j>i}c_j - \sum_{j<i}c_j
$$

となります。生成元をすべて掛けた積$e_1e_2\cdots e_m$では$d_i=m+1-2i$です。$n=2$なら$m$が奇数のとき$d_i$は常に偶数となって中心に入りますが、$n>2$では中心に入りません。代わりに中心に入るのは、冪を交互に反転させた積です。

&&&thm 中心元
$m$が奇数のとき、積

$$
\zeta = e_1e_2^{n-1}e_3e_4^{n-1}\cdots e_m
$$

は中心元であり、$\zeta^n$は$1$の冪根（スカラー）となる。
&&&

&&&prf
指数は$c_j=1\ (j\text{は奇数}),\ c_j=n-1\equiv-1\ (j\text{は偶数})$の交互列である。$\pm1$の交互列の和は、項数が偶数なら$0$、項数が奇数で端が$+1$なら$1$となる。$e_i$の左側$c_1,\dots,c_{i-1}$は$+1$で始まる交互列、右側$c_{i+1},\dots,c_m$は$+1$で終わる交互列（$m$が奇数のため）で、項数$i-1$と$m-i$の偶奇は一致する。よって両側の和は等しく、$d_i\equiv0\pmod n$、すなわち$\zeta$はすべての$e_i$と交換する。

また$\zeta^n$を関係式で並べ替えても生じるのは$\omega$の冪だけなので、ある整数$s$について

$$
\zeta^n = \omega^s\,e_1^{nc_1}\cdots e_m^{nc_m} = \omega^s
$$

となる。
&&&

$\mu^n=\zeta^n$（$=\omega^s$）となる1の冪根$\mu$を取れば$(\mu^{-1}\zeta)^n=1$なので、生成元が1個の場合とまったく同じ冪等元の構成によって、$\zeta$が生成する中心の部分代数$\mathbb C[\zeta]$は$n\mathbb C$へ直和分解されます。これを使えば、奇数個の場合の分類が完成します。

&&&thm 奇数個の生成元の行列表現
$$
\operatorname{Cl}^{(n)}_{2k+1}(\mathbb C)
\cong \operatorname{Cl}^{(n)}_{2k}(\mathbb C) \otimes n\mathbb C
\cong n\,\mathbb C(n^k)
$$
&&&

&&&prf
$m=2k+1$とおく。$e_1,\dots,e_{m-1}$は生成元$m-1$個の関係式を満たすので$\operatorname{Cl}^{(n)}_{m-1}(\mathbb C)$を生成し、また$\zeta^a$は位相倍を除いて単項式$e_1^ae_2^{-a}\cdots e_m^a$なので、$1,\zeta,\dots,\zeta^{n-1}$は線形独立で$\dim\mathbb C[\zeta]=n$となる。$\zeta$は中心元なので、掛け算$x \otimes \zeta^a \mapsto x\zeta^a$は代数準同型

$$
\operatorname{Cl}^{(n)}_{m-1}(\mathbb C) \otimes \mathbb C[\zeta]
\to \operatorname{Cl}^{(n)}_m(\mathbb C)
$$

を定める。$\zeta$の定義を逆に解くと$e_m = (e_1e_2^{n-1}\cdots e_{m-1}^{n-1})^{-1}\zeta$と表せるため、この準同型は全射である。次元は$n^{m-1} \cdot n = n^m$で一致するため同型となり、$\mathbb C[\zeta] \cong n\mathbb C$と偶数個の分類$\operatorname{Cl}^{(n)}_{m-1}(\mathbb C) \cong \mathbb C(n^{(m-1)/2})$を代入すれば主張が従う。
&&&

$k=0$（生成元1個）の場合は、$\zeta=e_1$で前述の$n\mathbb C$への直和分解そのものです。$n=2$では$e^{n-1}=e$より$\zeta$は擬スカラー$e_1e_2\cdots e_m$に戻り、複素クリフォード代数の奇数の場合の分類

$$
\operatorname{Cl}_{2k+1}(\mathbb C) \cong 2\,\mathbb C(2^k)
$$

に一致します。[[7shi-clif1]] $n>2$で擬スカラーの役割を引き継ぐのは、冪を交互に反転させた$\zeta$だというわけです。

# まとめ

本記事では、パウリ行列とその一般化を通じて、一般化クリフォード代数を導入しました。

&&&
パウリ行列：トレースが$0$の2×2エルミート行列の実基底
$$
i \cong -i\sigma_x, \quad j \cong -i\sigma_y, \quad k \cong -i\sigma_z
$$
$$
\operatorname{Cl}_{3,0}(\mathbb R) \cong \mathbb C(2) \cong \operatorname{Cl}_2(\mathbb C) \cong \mathbb C \otimes \mathbb H
$$
一般化パウリ行列：$\mathbb C(n)$の複素基底$X^aZ^b \ (0 \le a,b < n)$
$$
X^n = Z^n = I, \quad ZX = \omega XZ\ (\omega=e^{2\pi i/n})
$$
一般化クリフォード代数
$$
e_i^n = 1, \quad e_je_i = \omega e_ie_j \ (i<j), \quad
\operatorname{Cl}^{(n)}_2(\mathbb C) \cong \mathbb C(n)
$$
テンソル積による拡張
$$
\operatorname{Cl}^{(n)}_{m+2}(\mathbb C)
\cong \operatorname{Cl}^{(n)}_m(\mathbb C) \otimes \mathbb C(n), \quad
\operatorname{Cl}^{(n)}_{2k}(\mathbb C) \cong \mathbb C(n^k)
$$
奇数個の生成元：中心元による直和分解
$$
\zeta = e_1e_2^{n-1}e_3e_4^{n-1}\cdots e_m, \quad
\operatorname{Cl}^{(n)}_{2k+1}(\mathbb C) \cong n\,\mathbb C(n^k)
$$
&&&

クリフォード代数のテンソル積による拡張が、一般化クリフォード代数でも同じ形で機能することを確認しました。$n=2$の複素2周期性は、この拡張公式の特殊な場合です。奇数個の生成元では、擬スカラーに代わって冪を交互に反転させた中心元が直和分解を与え、生成元が1個の場合はフルーリーの多重複素数の複素化という可換な特殊例が現れます。
