シリーズ：[四元数の行列表現](https://mathlog.info/series/PXPuUuQLYZk6HHho9eP8)

四元数の虚数単位を複素2次正方行列で表現して、そこからパウリ行列を構成します。

# 四元数の基本的性質

四元数は実部と3つの虚数単位を持つ虚部からなる数体系で、一般的に以下のように表されます。

&&&def 四元数
$$
a + bi + cj + dk \quad (a, b, c, d \in \mathbb{R})
$$

ここで虚数単位$i,j,k$は以下の関係式を満たす。

$$
i^2 = j^2 = k^2 = -1
$$
$$
ij = -ji = k, \quad jk = -kj = i, \quad ki = -ik = j
$$
&&&

四元数は次のように2つの複素数の組として表現することもできます。

$$
a + bi + cj + dk = (a + bi) + (c + di)j
$$

この形式は外側が$X + Yj$の形になっています。この表現から、四元数の行列表現を考える際、まず外側の虚数単位$j$を複素数の虚数単位$i$に対応付けることから始めます。

&&&rem
行列表現には任意性があるため、扱いやすいようにまず$j$を固定するということです。
&&&

# 四元数の行列表現

四元数の虚数単位$i,j,k$を$2 \times 2$行列で表現することを考えます。これらの行列は四元数全体の集合を表す$\mathbb{H}$を添字にして$I_H,J_H,K_H$と表記します。

&&&rem 本記事での方針
単純に大文字にすると$I$が単位行列と衝突します。それを避けるため$I_H$などと表記します。
&&&

先ほどの動機付けから、まず$J_H$を複素数の虚数単位$i$の行列表現として以下のように定義します。

&&&def 四元数の虚数単位$j$の行列表現
$$
J_H = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
$$
&&&

&&&prf
$J_H$が複素数$i$の行列表現であることを確認します。実平面上のベクトルに$J_H$を作用させます。

$$
J_H \begin{pmatrix} x \\ y \end{pmatrix}
= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix}
= \begin{pmatrix} -y \\ x \end{pmatrix}
$$

これは複素数の演算$i(x+yi)=-y+xi$に対応します。
&&&

$J_H$を複素数の虚数単位$i$と同一視することで、四元数の行列表現の出発点とします。

次に、四元数の性質を満たすためには、これらの行列は以下の条件を満たす必要があります。（$I$は単位行列）

$$
{I_H}^2 = {J_H}^2 = {K_H}^2 = -I
$$
$$
\begin{alignedat}{2}
I_H J_H &= -J_H I_H &&= K_H \\
J_H K_H &= -K_H J_H &&= I_H \\
K_H I_H &= -I_H K_H &&= J_H \\
\end{alignedat}
$$

## $I_H$の導出

$I_H$を一般的な$2 \times 2$行列として

$$
I_H = \begin{pmatrix} a & b \\ c & d \end{pmatrix}
$$

とおきます。条件$I_H J_H = -J_H I_H$より

$$
\begin{aligned}
\begin{pmatrix} a & b \\ c & d \end{pmatrix} \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
&= -\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} a & b \\ c & d \end{pmatrix} \\
\begin{pmatrix} b & -a \\ d & -c \end{pmatrix}
&= \begin{pmatrix} c & d \\ -a & -b \end{pmatrix}
\end{aligned}
$$

よって$b=c,\ d=-a$より

$$
I_H = \begin{pmatrix} a & b \\ b & -a \end{pmatrix}
$$

## ${I_H}^2 = -I$の条件

$$
{I_H}^2
= \begin{pmatrix} a & b \\ b & -a \end{pmatrix}^2
= \begin{pmatrix} a^2 + b^2 & 0 \\ 0 & a^2 + b^2 \end{pmatrix}
= (a^2 + b^2)I
$$

${I_H}^2 = -I$より$a^2 + b^2 = -1$

この条件を満たす実数$a,b$は存在しないため、複素数の範囲で解を求めます。$\theta$を任意の実数とすれば

$$
-1 = i^2(\sin^2\theta+\cos^2\theta) = (i\sin\theta)^2+(i\cos\theta)^2
$$

よって

$$
a = i\sin\theta, \quad b = i\cos\theta
$$

とパラメーター表示できます。

&&&rem パラメーター表示の任意性
パラメーター表示には任意性があるため、例えば

$$
a = -i\cos\theta, \quad b = i\sin\theta
$$

のような組み合わせも解となります。これは$θ$の位相がずれるだけで本質的な違いではありません。
&&&

## $K_H$の導出

$$
K_H = I_H J_H = \begin{pmatrix} a & b \\ b & -a \end{pmatrix} \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} b & -a \\ -a & -b \end{pmatrix}
$$

パラメーター表示を用いると

$$
\begin{aligned}
I_H &= i \begin{pmatrix} \sin\theta &  \cos\theta \\  \cos\theta & -\sin\theta \end{pmatrix} \\
K_H &= i \begin{pmatrix} \cos\theta & -\sin\theta \\ -\sin\theta & -\cos\theta \end{pmatrix}
\end{aligned}
$$

$I_H,J_H$が四元数$i,j$と同じ性質を持つことから、その積によって得られた$K_H$は自動的に$k$と同じ性質を持ちます。

$$
{K_H}^2 = (I_H J_H)^2 = I_H J_H I_H J_H = -I_H I_H J_H J_H = -(-I)(-I) = -I
$$

&&&rem
四元数で$k$を含む演算は、$k=ij$により$i,j$のみの計算に還元できます。
&&&

## 行列表現の固定

$\theta$には任意性があります。取り扱いを容易にするため、$\sin\theta=0$となるように$\theta$を固定します。$\theta=\pi$で固定すれば、四元数の虚数単位の行列表現は以下のようになります。

&&&def 四元数の虚数単位の行列表現
\begin{alignedat}{3}
I_H &= i &&\begin{pmatrix} \sin\pi & \cos\pi \\ \cos\pi & -\sin\pi \end{pmatrix}
   &&= \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix} \\
J_H &&&
   &&= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \\
K_H &= i &&\begin{pmatrix} \cos\pi & -\sin\pi \\ -\sin\pi & -\cos\pi \end{pmatrix}
   &&= \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
\end{alignedat}
&&&

&&&rem
$θ=0$とする方が素直に思えますが、後で構成するパウリ行列と辻褄を合わせるため$θ=π$としています。
&&&

# 共役四元数とエルミート共役

複素数では虚部の符号を反転することで**共役複素数**を得ます。共役を$*$で表し、以下のように定義します。

&&&def 共役複素数
$$
(a+bi)^* = a - bi
$$
&&&

四元数の場合も、同様に虚部の符号を反転することで**共役四元数**を定義します。

&&&def 共役四元数
$$
(a+bi+cj+dk)^* = a - bi - cj - dk
$$
&&&

複素数の虚数単位$i$の行列表現でもある$J_H$は、転置によって符号が反転します。このことから、複素数の行列表現では転置が共役に対応付けられます。

&&&prf 転置で符号反転
$$
{J_H}^\mathrm{T}
= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}^\mathrm{T}
= \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}
= -J_H
$$
&&&

一方、$I_H,K_H$は転置で変化しません。このような性質を持つ行列を**対称行列**と呼びます。

&&&prf 対称行列
$$
\begin{alignedat}{3}
{I_H}^\mathrm{T}
 &= \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}^\mathrm{T}
&&= \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}
&&= I_H \\
{K_H}^\mathrm{T}
 &= \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}^\mathrm{T}
&&= \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
&&= K_H
\end{alignedat}
$$
&&&

成分が純虚数であることから、成分の共役を取ることで符号が反転します。

&&&prf 共役で符号反転
$$
\begin{alignedat}{3}
{I_H}^*
 &= \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}^*
&&= \begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix}
&&= -I_H \\
{K_H}^*
 &= \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}^*
&&= \begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}
&&= -K_H
\end{alignedat}
$$
&&&

$J_H$は実行列のため共役で変化しないことから、$I_H,J_H,K_H$をすべて符号反転させるには、共役と転置を同時に行う必要があります。このような操作を$\dagger$で表し、**エルミート共役**と呼びます。

$$
\begin{alignedat}{4}
{I_H}^\dagger
 &= ({I_H}^*)^\mathrm{T}
&&= \begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix}^\mathrm{T}
&&= \begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix}
&&= -I_H \\
{J_H}^\dagger
 &= ({J_H}^*)^\mathrm{T}
&&= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}^\mathrm{T}
&&= \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}
&&= -J_H \\
{K_H}^\dagger
 &= ({K_H}^*)^\mathrm{T}
&&= \begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}^\mathrm{T}
&&= \begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}
&&= -K_H
\end{alignedat}
$$

よって、四元数の行列表現のエルミート共役は、四元数の共役に対応します。

&&&rem 実行列表現での転置
$2×2$の複素行列は、成分の複素数を行列表現に置き換えることで、$4×4$の実行列に変換できます。同じ性質を持つ異なる表現は、同型を表す$\cong$で結びます。

$$
\begin{pmatrix}
  a+bi & c+di \\
  e+fi & g+hi
\end{pmatrix} \cong
\begin{pmatrix}
  a & -b & c & -d \\
  b &  a & d &  c \\
  e & -f & g & -h \\
  f &  e & h & g
\end{pmatrix}
$$

このように変換した実行列の転置を取ります。

$$
\begin{pmatrix}
  a & -b & c & -d \\
  b &  a & d &  c \\
  e & -f & g & -h \\
  f &  e & h & g
\end{pmatrix}^\mathrm{T}
=
\begin{pmatrix}
   a & b &  e & f \\
  -b & a & -f & e \\
   c & d &  g & h \\
  -d & c & -h & g
\end{pmatrix}
$$

これを複素行列に戻せば、エルミート共役と一致します。

$$
\begin{pmatrix}
   a & b &  e & f \\
  -b & a & -f & e \\
   c & d &  g & h \\
  -d & c & -h & g
\end{pmatrix} \\
\cong
\begin{pmatrix}
  a-bi & e-fi \\
  c-di & g-hi
\end{pmatrix} \\
=
\begin{pmatrix}
  a+bi & c+di \\
  e+fi & g+hi
\end{pmatrix}^\dagger
$$

つまり、複素行列は実行列の圧縮表現であると考えることで、転置に共役を伴うエルミート共役が自然な操作であると解釈できます。
&&&

# パウリ行列の構成

四元数の虚数単位の行列表現は2乗で$-I$になります。それらを$i$倍すれば、2乗で$I$になる行列が得られます。そうして得られた行列を$\sigma_1,\sigma_2,\sigma_3$と表記して**パウリ行列**と呼びます。

&&&def パウリ行列
$$
\begin{alignedat}{3}
\sigma_1 &:= i I_H
         &&= i \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}
         &&= \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \\
\sigma_2 &:= i J_H
         &&= i \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
         &&= \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix} \\
\sigma_3 &:= i K_H
         &&= i \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
         &&= \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
\end{alignedat}
$$
$$
{\sigma_1}^2={\sigma_2}^2={\sigma_3}^2=I
$$
&&&

&&&rem テンソル積
$i J_H$のような構造は、複素数と四元数のテンソル積$i⊗j$として一般化されます。つまり、パウリ行列は複素数と四元数のテンソル積として理解できます。[[7shi-bq]]
&&&

パウリ行列はエルミート共役で変化しません。このような性質を持つ行列を**エルミート行列**と呼びます。

&&&prf エルミート行列
$$
\begin{alignedat}{3}
{\sigma_1}^\dagger
 &= \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}^\dagger
&&= \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
&&= \sigma_1 \\
{\sigma_2}^\dagger
 &= \begin{pmatrix} 0 & i \\ -i & 0 \end{pmatrix}^\dagger
&&= \begin{pmatrix} 0 & i \\ -i & 0 \end{pmatrix}
&&= \sigma_2 \\
{\sigma_3}^\dagger
 &= \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}^\dagger
&&= \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
&&= \sigma_3
\end{alignedat}
$$
&&&

&&&rem エルミート行列と対称行列
「実行列表現での転置」で見たように、複素行列を実行列に変換すればエルミート共役は転置となるため、エルミート行列は対称行列となります。つまり、実行列における対称行列の概念を複素行列に拡張したのがエルミート行列です。

なお、$\sigma_1$と$\sigma_3$は成分に虚数を含まないため、そのままでも対称行列です。
&&&

&&&rem 反エルミート行列
$I_H,J_H,K_H$はエルミート共役によって符号が反転することから、**反エルミート行列**と呼ばれます。

$I_H,J_H,K_H$に掛けた$i$が共役によって符号反転することから、符号反転が相殺することでパウリ行列がエルミート行列になると解釈できます。
&&&

パウリ行列$\sigma_1,\sigma_2,\sigma_3$の性質を確認します。四元数の行列表現に還元することで、成分計算の手間が省けます。

$$
\begin{alignedat}{4}
\sigma_1\sigma_2 &= (iI_H)(iJ_H) &&= ii(I_H J_H) &&= i(iK_H) &&= i\sigma_3 = -K_H \\
\sigma_2\sigma_3 &= (iJ_H)(iK_H) &&= ii(J_H K_H) &&= i(iI_H) &&= i\sigma_1 = -I_H \\
\sigma_3\sigma_1 &= (iK_H)(iI_H) &&= ii(K_H I_H) &&= i(iJ_H) &&= i\sigma_2 = -J_H \\
\end{alignedat}
$$
$$
\begin{alignedat}{5}
\sigma_2\sigma_1 &= (iJ_H)(iI_H) &&= ii(J_H I_H) &&= -i(iK_H) &&= -i\sigma_3 = K_H &&= -\sigma_1\sigma_2 \\
\sigma_3\sigma_2 &= (iK_H)(iJ_H) &&= ii(K_H J_H) &&= -i(iI_H) &&= -i\sigma_1 = I_H &&= -\sigma_2\sigma_3 \\
\sigma_1\sigma_3 &= (iI_H)(iK_H) &&= ii(I_H K_H) &&= -i(iJ_H) &&= -i\sigma_2 = J_H &&= -\sigma_3\sigma_1 \\
\end{alignedat}
$$

&&&rem $i$の由来
$\sigma_1\sigma_2=i\sigma_3$という結果だけ見ると、どこからか$i$が湧き出したように見えます。（成分計算すればそうなるというのは当然として）

$(iI_H)(iJ_H)=i(iK_H)$という計算過程により、元々含まれていた$i$が片方だけ出てきたと解釈できます。
&&&

結果をまとめます。積によって$i$が現れ、四元数の虚数単位に対応します。

&&&thm パウリ行列の積の関係
\begin{aligned}
\sigma_1\sigma_2 &= -\sigma_2\sigma_1 = i\sigma_3 = -K_H \\
\sigma_2\sigma_3 &= -\sigma_3\sigma_2 = i\sigma_1 = -I_H \\
\sigma_3\sigma_1 &= -\sigma_1\sigma_3 = i\sigma_2 = -J_H \\
\end{aligned}
&&&

&&&rem 因数分解とクリフォード代数
$I_H=\sigma_3\sigma_2$などの関係から、四元数の虚数単位を因数分解したのがパウリ行列だと解釈できます。この考え方は四元数よりもパウリ行列が基本的な構造であるという見方を反映しており、クリフォード代数につながります。
&&&

&&&rem ホッジスター
$i$は外積代数におけるホッジスターに対応します。

$$
i\sigma_1
\cong \star e_1
= e_2 \wedge e_3
\cong \sigma_2\sigma_3
$$
&&&

パウリ行列は量子力学や量子計算において、電子のスピンや量子ビットの状態を表すための基本的な演算子として重要な役割を果たします。

# 双四元数

行列表現における係数の虚数単位を$h$と表現することで、行列によらず四元数の拡張としてパウリ行列と同型の代数を表現できます。これを**双四元数**と呼びます。[[7shi-bq]]

&&&def 双四元数
$$
h \cong iI = \sigma_1\sigma_2\sigma_3,\quad
hi \cong iI_H = \sigma_1,\quad
hj \cong iJ_H = \sigma_2,\quad
hk \cong iK_H = \sigma_3
$$

一般形：

$$
\begin{aligned}
&a_0 + a_1hi + a_2hj + a_3hk + a_4i + a_5j + a_6k + a_7h \\
&\cong a_0I + a_1\sigma_1 + a_2\sigma_2 + a_3\sigma_3 + a_4\sigma_3\sigma_2 + a_5\sigma_1\sigma_3 + a_6\sigma_2\sigma_1 + a_7\sigma_1\sigma_2\sigma_3 \\
&= a_0I + a_1\sigma_1 + a_2\sigma_2 + a_3\sigma_3 - a_4i\sigma_1 - a_5i\sigma_2 - a_6i\sigma_3 + a_7iI \\
\end{aligned}
$$
&&&

&&&rem テンソル積
双四元数は、テンソル積の構造を「四元数の係数を複素数に拡張」という形で簡略化した表記だと解釈できます。

$$
\sigma_2 = iJ_H \cong i⊗j \cong hj
$$
&&&

双四元数とパウリ行列の対応を把握すれば、パウリ行列と四元数の関係も理解が深まるでしょう。

&&&ex
$$
{\sigma_1}^2 \cong (hi)^2 = h^2i^2 = (-1)(-1) = 1
\quad \because (hi)^2 \cong (iI_H)^2=i^2{I_H}^2=I
$$
$$
\sigma_1\sigma_2 \cong (hi)(hj) = hh(ij) = h(hk) \cong i\sigma_3
$$
&&&
