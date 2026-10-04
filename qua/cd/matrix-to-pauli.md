四元数の虚数単位を複素2次正方行列で表現し、そこからパウリ行列および双四元数との代数的な対応関係を導出します。

シリーズ：[四元数の行列表現](https://mathlog.info/series/PXPuUuQLYZk6HHho9eP8)

&&& 改訂履歴
- 2026.10.04 一般の四元数の行列表現とエルミート共役の対応を明示し、パラメータ表示の選択理由を補い、双四元数の定義とホッジスターの符号の扱いを修正
&&&

# 概要

四元数の虚数単位$i, j, k$は、反交換関係と2乗が$-1$となる性質を持ちます。これらを$2 \times 2$複素行列として表現すると、パウリ行列との対応を具体的に記述できます。本記事では、虚数単位$j$を実数平面上の虚数単位（$90^\circ$回転行列）として固定する発想から出発し、$I_H, J_H, K_H$の行列表現をパラメータ表示を通じて導出します。さらに、複素数の虚数単位$i$をスカラー倍として掛けることでパウリ行列$σ_1, σ_2, σ_3$を構成し、双四元数やクリフォード代数との接続を整理します。[[7shi-bq]]

複素数と2次正方行列の基本演算を前提とします。回転や量子状態への応用は扱いません。

# 四元数の基本的性質

四元数は実部と3つの虚数単位を持つ虚部からなる数体系で、一般的に以下のように表されます。

&&&def 四元数
$$
a + bi + cj + dk \quad (a, b, c, d \in \mathbb{R})
$$

ここで虚数単位$i,j,k$は以下の関係式を満たします。

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

&&&rem $j$の固定と任意性
行列表現には任意性があるため、ここでは計算しやすい$j$の表現を先に選びます。
&&&

# 四元数の行列表現

四元数の虚数単位$i,j,k$を$2 \times 2$行列で表現することを考えます。これらの行列は四元数全体の集合を表す$\mathbb{H}$を添字にして$I_H,J_H,K_H$と表記します。

&&&rem 本記事での方針
単純に大文字にすると$I$が単位行列と衝突します。それを避けるため$I_H$などと表記します。
&&&

先ほどの動機付けから、複素数$i$の実行列表現と同じ行列を、四元数の$j$の像$J_H$として採用します。

&&&def 四元数の虚数単位$j$の行列表現
$$
J_H = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
$$
&&&

&&&ex 複素数の$i$倍との対応
実平面上のベクトルに$J_H$を作用させます。

$$
J_H \begin{pmatrix} x \\ y \end{pmatrix}
= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix}
= \begin{pmatrix} -y \\ x \end{pmatrix}
$$

これは複素数の演算$i(x+yi)=-y+xi$に対応します。
&&&

この$J_H$を四元数の行列表現の出発点とします。行列の成分やスカラー倍に現れる$i$は複素数の虚数単位で、四元数の$i$の像は$I_H$です。

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

$I_H$を、複素数$a,b,c,d$を成分とする一般的な$2 \times 2$行列としておきます。

$$
I_H = \begin{pmatrix} a & b \\ c & d \end{pmatrix}
$$

条件$I_H J_H = -J_H I_H$より、次の式が得られます。

$$
\begin{aligned}
\begin{pmatrix} a & b \\ c & d \end{pmatrix} \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
&= -\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} a & b \\ c & d \end{pmatrix} \\
\begin{pmatrix} b & -a \\ d & -c \end{pmatrix}
&= \begin{pmatrix} c & d \\ -a & -b \end{pmatrix}
\end{aligned}
$$

よって$b=c,\ d=-a$となり、次の形になります。

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

${I_H}^2 = -I$より$a^2 + b^2 = -1$となります。

この条件を満たす実数$a,b$は存在しないため、複素数の範囲で解を求めます。$θ$を任意の実数とすれば、次が成り立ちます。

$$
-1 = i^2(\sin^2θ+\cos^2θ) = (i\sinθ)^2+(i\cosθ)^2
$$

そこで、解の族として次を選びます。

$$
a = i\sinθ, \quad b = i\cosθ
$$

これは複素解の全体ではなく、$a,b$がともに純虚数となる解の族です。この選択により、$I_H,K_H$は反エルミート行列になります。

&&&rem パラメーター表示の任意性
パラメーター表示には任意性があるため、例えば次のような組み合わせも解となります。

$$
a = -i\cosθ, \quad b = i\sinθ
$$

これは$θ$の位相がずれるだけで本質的な違いではありません。
&&&

## $K_H$の導出

$$
K_H = I_H J_H = \begin{pmatrix} a & b \\ b & -a \end{pmatrix} \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} b & -a \\ -a & -b \end{pmatrix}
$$

パラメーター表示を用いると、次のようになります。

$$
\begin{aligned}
I_H &= i \begin{pmatrix} \sinθ &  \cosθ \\  \cosθ & -\sinθ \end{pmatrix} \\
K_H &= i \begin{pmatrix} \cosθ & -\sinθ \\ -\sinθ & -\cosθ \end{pmatrix}
\end{aligned}
$$

$I_H,J_H$が四元数$i,j$と同じ性質を持つことから、その積によって得られた$K_H$は自動的に$k$と同じ性質を持ちます。

$$
{K_H}^2 = (I_H J_H)^2 = I_H J_H I_H J_H = -I_H I_H J_H J_H = -(-I)(-I) = -I
$$

&&&rem $k$の還元
四元数で$k$を含む演算は、$k=ij$により$i,j$のみの計算に還元できます。
&&&

## 行列表現の固定

$θ$には任意性があります。取り扱いを容易にするため、$\sinθ=0$となるように$θ$を固定します。$θ=π$で固定すれば、四元数の虚数単位の行列表現は以下のようになります。

&&&def 四元数の虚数単位の行列表現
\begin{alignedat}{3}
I_H &= i &&\begin{pmatrix} \sinπ & \cosπ \\ \cosπ & -\sinπ \end{pmatrix}
   &&= \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix} \\
J_H &&&
   &&= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \\
K_H &= i &&\begin{pmatrix} \cosπ & -\sinπ \\ -\sinπ & -\cosπ \end{pmatrix}
   &&= \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
\end{alignedat}
&&&

&&&rem $θ=π$の選択理由
$θ=0$とする方が素直に思えますが、その場合は後で定義する標準的なパウリ行列に対して$iI_H=-σ_1,\ iJ_H=σ_2,\ iK_H=-σ_3$となります。3つとも符号をそろえるため$θ=π$としています。
&&&

虚数単位の行列を実係数で線形結合して、一般の四元数を行列で表します。

&&&def 四元数の行列表現
$$
ρ(a+bi+cj+dk)=aI+bI_H+cJ_H+dK_H\quad(a,b,c,d\in\mathbb R)
$$
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

係数$a,b,c,d$は実数なので、エルミート共役で係数は変わりません。よって、一般の四元数$q=a+bi+cj+dk$について次が成り立ち、四元数の行列表現のエルミート共役は四元数の共役に対応します。

$$
ρ(q)^\dagger=aI-bI_H-cJ_H-dK_H=ρ(q^*)
$$

&&&rem 実行列表現での転置
$2×2$の複素行列は、成分の複素数を行列表現に置き換えることで、$4×4$の実行列に変換できます。この変換は加法と乗法を保ち、複素行列の代数と、その像である実行列の部分代数との同型を与えます。この同型による対応を$\cong$で表します。

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

四元数の虚数単位の行列表現は2乗で$-I$になります。それらを$i$倍すれば、2乗で$I$になる行列が得られます。そうして得られた行列を$σ_1,σ_2,σ_3$と表記して**パウリ行列**と呼びます。

&&&def パウリ行列
$$
\begin{alignedat}{3}
σ_1 &:= i I_H
         &&= i \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}
         &&= \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \\
σ_2 &:= i J_H
         &&= i \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
         &&= \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix} \\
σ_3 &:= i K_H
         &&= i \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
         &&= \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
\end{alignedat}
$$
$$
{σ_1}^2={σ_2}^2={σ_3}^2=I
$$
&&&

$σ_1,σ_2,σ_3$は、それぞれ$σ_x,σ_y,σ_z$とも表記されます。

&&&rem テンソル積
実数体上のテンソル積$\mathbb C\otimes_{\mathbb R}\mathbb H$から$M_2(\mathbb C)$への同型$z\otimes q\mapsto z\,ρ(q)$のもとで、$σ_2=iJ_H$は$i⊗j$に対応します。左の$i$は複素数、右の$j$は四元数です。このように、パウリ行列は複素数と四元数のテンソル積の元として理解できます。[[7shi-bq]]
&&&

パウリ行列はエルミート共役で変化しません。このような性質を持つ行列を**エルミート行列**と呼びます。

&&&prf エルミート行列
$$
\begin{alignedat}{3}
{σ_1}^\dagger
 &= \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}^\dagger
&&= \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
&&= σ_1 \\
{σ_2}^\dagger
 &= \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}^\dagger
&&= \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}
&&= σ_2 \\
{σ_3}^\dagger
 &= \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}^\dagger
&&= \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
&&= σ_3
\end{alignedat}
$$
&&&

&&&rem エルミート行列と対称行列
「実行列表現での転置」で見たように、複素行列を実行列に変換すればエルミート共役は転置となるため、エルミート行列は対称行列となります。つまり、実行列における対称行列の概念を複素行列に拡張したのがエルミート行列です。

なお、$σ_1$と$σ_3$は成分に虚数を含まないため、そのままでも対称行列です。
&&&

&&&rem 反エルミート行列
$I_H,J_H,K_H$はエルミート共役によって符号が反転することから、**反エルミート行列**と呼ばれます。

$I_H,J_H,K_H$に掛けた$i$が共役によって符号反転することから、符号反転が相殺することでパウリ行列がエルミート行列になると解釈できます。
&&&

パウリ行列$σ_1,σ_2,σ_3$の性質を確認します。四元数の行列表現に還元することで、成分計算の手間が省けます。

$$
\begin{alignedat}{4}
σ_1σ_2 &= (iI_H)(iJ_H) &&= ii(I_H J_H) &&= i(iK_H) &&= iσ_3 = -K_H \\
σ_2σ_3 &= (iJ_H)(iK_H) &&= ii(J_H K_H) &&= i(iI_H) &&= iσ_1 = -I_H \\
σ_3σ_1 &= (iK_H)(iI_H) &&= ii(K_H I_H) &&= i(iJ_H) &&= iσ_2 = -J_H \\
\end{alignedat}
$$
$$
\begin{alignedat}{5}
σ_2σ_1 &= (iJ_H)(iI_H) &&= ii(J_H I_H) &&= -i(iK_H) &&= -iσ_3 = K_H &&= -σ_1σ_2 \\
σ_3σ_2 &= (iK_H)(iJ_H) &&= ii(K_H J_H) &&= -i(iI_H) &&= -iσ_1 = I_H &&= -σ_2σ_3 \\
σ_1σ_3 &= (iI_H)(iK_H) &&= ii(I_H K_H) &&= -i(iJ_H) &&= -iσ_2 = J_H &&= -σ_3σ_1 \\
\end{alignedat}
$$

&&&rem $i$の由来
$σ_1σ_2=iσ_3$の係数$i$は、四元数の行列表現をパウリ行列に置き換える際のスカラー倍に由来します。$(iI_H)(iJ_H)=i(iK_H)$と書けば、その対応を追えます。
&&&

結果をまとめます。積によって$i$が現れ、四元数の虚数単位に対応します。

&&&fml パウリ行列の積の関係
\begin{aligned}
σ_1σ_2 &= -σ_2σ_1 = iσ_3 = -K_H \\
σ_2σ_3 &= -σ_3σ_2 = iσ_1 = -I_H \\
σ_3σ_1 &= -σ_1σ_3 = iσ_2 = -J_H \\
\end{aligned}
&&&

&&&rem 因数分解とクリフォード代数
$I_H=σ_3σ_2$などの関係から、四元数の虚数単位を因数分解したのがパウリ行列だと解釈できます。クリフォード代数$\operatorname{Cl}_{3,0}(\mathbb R)$の見方では、パウリ行列をベクトルに対応する生成元とし、四元数を偶部分代数$\operatorname{Cl}_{3,0}^0(\mathbb R)$として捉えられます。この意味で、四元数の虚数単位をベクトルの積へ分解していると解釈できます。
&&&

&&&rem ホッジスター
基底ベクトル$e_r$をパウリ行列$σ_r$に対応させると、1次の元に対する外積代数のホッジスターは$i$倍に対応します。

$$
iσ_1
\cong \star e_1
= e_2 \wedge e_3
\cong σ_2σ_3
$$

ただし、2次の元では$\star(e_1\wedge e_2)=e_3$に対して$i(σ_1σ_2)=-σ_3$となり、符号が逆になります。
&&&

# 双四元数

行列表現における係数の虚数単位を$h$と表現することで、行列によらず四元数の拡張として、パウリ行列が生成する$M_2(\mathbb C)$と同型の代数を表現できます。これを**双四元数**と呼びます。[[7shi-bq]]

&&&def 双四元数
双四元数は、実数体上のテンソル積$\mathbb C\otimes_{\mathbb R}\mathbb H$です。複素数側の虚数単位を$h$と書き、次の規則に従います。

$$
h^2=-1,\quad hi=ih,\quad hj=jh,\quad hk=kh
$$

一般形は次の通りです。

$$
a_0 + a_1hi + a_2hj + a_3hk + a_4i + a_5j + a_6k + a_7h\quad(a_0,\ldots,a_7\in\mathbb R)
$$
&&&

&&&fml 双四元数とパウリ行列の対応
$$
h \cong iI = σ_1σ_2σ_3,\quad
hi \cong iI_H = σ_1,\quad
hj \cong iJ_H = σ_2,\quad
hk \cong iK_H = σ_3
$$

一般形は次のように対応します。

$$
\begin{aligned}
&a_0 + a_1hi + a_2hj + a_3hk + a_4i + a_5j + a_6k + a_7h \\
&\cong a_0I + a_1σ_1 + a_2σ_2 + a_3σ_3 + a_4σ_3σ_2 + a_5σ_1σ_3 + a_6σ_2σ_1 + a_7σ_1σ_2σ_3 \\
&= a_0I + a_1σ_1 + a_2σ_2 + a_3σ_3 - a_4iσ_1 - a_5iσ_2 - a_6iσ_3 + a_7iI \\
\end{aligned}
$$
&&&

&&&rem テンソル積
双四元数は、四元数の係数を複素数に拡張した代数です。複素数側の虚数単位を$h$と書き、テンソル積の記号を省略しています。

$$
σ_2 = iJ_H \cong i⊗j \cong hj
$$
&&&

&&&ex 双四元数における演算例
$$
{σ_1}^2 \cong (hi)^2 = h^2i^2 = (-1)(-1) = 1
\quad \because (hi)^2 \cong (iI_H)^2=i^2{I_H}^2=I
$$
$$
σ_1σ_2 \cong (hi)(hj) = hh(ij) = h(hk) \cong iσ_3
$$
&&&

# まとめ

本記事では、虚数単位$j$を$90^\circ$回転行列として固定する発想から出発し、四元数の虚数単位の複素2次正方行列による表現$I_H,J_H,K_H$をパラメーター表示を通じて導きました。パラメーターの任意性は、標準的なパウリ行列と符号がそろうように固定しました。また、複素行列を実行列の圧縮表現とみなすことで、転置に共役を伴うエルミート共役が自然な操作として現れることを見ました。

四元数の虚数単位の行列表現を$i$倍すると、2乗が$I$になるパウリ行列$σ_1,σ_2,σ_3$が得られます。逆に、パウリ行列の積は四元数の虚数単位を与えるため、四元数の虚数単位を因数分解したのがパウリ行列だと解釈できます。行列の係数に現れる虚数単位を$h$と表すと、行列によらず四元数の拡張として、パウリ行列が生成する$M_2(\mathbb C)$と同型の代数（双四元数）が得られます。

&&& 四元数の虚数単位の行列表現
$$
I_H = \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}, \quad
J_H = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \quad
K_H = \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
$$
&&&

&&& パウリ行列の積
$$
σ_1σ_2 = -σ_2σ_1 = iσ_3 = -K_H, \quad
σ_2σ_3 = -σ_3σ_2 = iσ_1 = -I_H, \quad
σ_3σ_1 = -σ_1σ_3 = iσ_2 = -J_H
$$
&&&

&&& 双四元数との対応
$$
h \cong iI, \quad hi \cong σ_1, \quad hj \cong σ_2, \quad hk \cong σ_3
$$
&&&
