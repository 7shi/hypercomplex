パウリ行列から導かれる実パウリ行列が2次正方実行列の基底をなす構造を整理し、トレースとフロベニウス内積による成分抽出を解説します。

&&& 改訂履歴
- 2026.10.04 基底であることの論証と交差項が消える理由を補い、フロベニウス内積を定義して全係数の抽出式と直交基底としての解釈を追加した
&&&

# 概要

パウリ行列は通常複素数を含みますが、$\sigma_2$を$-i\sigma_2$に置き換えることで、すべてが実行列となる組を構成できます。本記事ではこの組$\{\tau_1, \tau_2, \tau_3\}$を実パウリ行列と呼びます。単位行列と実パウリ行列の組$\{I, \tau_1, \tau_2, \tau_3\}$が2次正方実行列空間$M_2(\mathbb{R})$の基底をなすことを確認します。さらに、各基底のトレースが0になる性質を利用して、転置行列との積のトレース$\operatorname{tr}(A^\mathsf{T}B)$が係数ベクトルの標準内積に対応することと、実パウリ行列が直交基底をなすことを示します。

行列の積・転置・トレースと、ベクトル空間の基底を前提とします。量子力学の知識は必要ありません。

# パウリ行列

通常のパウリ行列$\sigma_1, \sigma_2, \sigma_3$は以下のように定義されます。[[7shi-bq]]

$$
\sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad \sigma_2 = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad \sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$

これらは$2 \times 2$の複素行列であり、以下の関係式を満たします。

$$
\begin{aligned}
\sigma_1^2 &= \sigma_2^2 = \sigma_3^2 = I \\
\sigma_1\sigma_2 &= -\sigma_2\sigma_1 = i\sigma_3 \\
\sigma_2\sigma_3 &= -\sigma_3\sigma_2 = i\sigma_1 \\
\sigma_3\sigma_1 &= -\sigma_1\sigma_3 = i\sigma_2
\end{aligned}
$$

# 実パウリ行列の構成

上記の関係式$\sigma_1\sigma_2 = i\sigma_3$を出発点とし、この式を実数のみの行列で表現することを考えます。

両辺に$-i$を掛けます。

$$
-i\sigma_1\sigma_2 = \sigma_3
$$

スカラー$-i$は行列と交換するので、次のように書けます。

$$
\sigma_1(-i\sigma_2) = \sigma_3
$$

$\sigma_1,\sigma_3$は実行列ですが、$-i\sigma_2$も実行列となります。

$$
-i\sigma_2 = -i\begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
$$

よってすべて実行列の計算となります。

$$
\underbrace{\pmatrix{0&1\\1&0}}_{\sigma_1}\underbrace{\pmatrix{0&-1\\1&0}}_{-i\sigma_2}=\underbrace{\pmatrix{1&0\\0&-1}}_{\sigma_3}
$$

この関係式に着目し、実パウリ行列$\tau_1, \tau_2, \tau_3$を以下のように定義します。

&&&def 実パウリ行列
$$
\begin{alignedat}{2}
\tau_1 &:= &\sigma_1 &= \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \\
\tau_2 &:= &-i\sigma_2 &= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \\
\tau_3 &:= &\sigma_3 &= \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
\end{alignedat}
$$
&&&

この定義では、$\tau_1, \tau_2, \tau_3$ は全て実行列であり、以下の関係式を満たします。

$$
\begin{aligned}
\tau_1^2 = I&,\ \tau_2^2=-I,\ \tau_3^2=I \\
\tau_1\tau_2 &= -\tau_2\tau_1 = \tau_3 \\
\tau_1\tau_3 &= -\tau_3\tau_1 = \tau_2 \\
\tau_2\tau_3 &= -\tau_3\tau_2 = \tau_1
\end{aligned}
$$

## $\tau_2$の持つ意味

この定義で得られた$\tau_2$は、複素数を実2次元ベクトルとみなしたときの、虚数単位$i$による乗算の表現行列と見なすことができます。実際、$\tau_2$の2乗は$i^2 = -1$に対応します。

$$
\tau_2^2
= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
= \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}
= -I
$$

また、複素数$x+yi$をベクトル$\begin{pmatrix} x \\ y \end{pmatrix}$と同一視すれば、$\tau_2$をこのベクトルに作用させることは、複素数に$i$を掛ける操作に対応します。

$$
\tau_2 \begin{pmatrix} x \\ y \end{pmatrix}
= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix}
= \begin{pmatrix} -y \\ x \end{pmatrix}
\leftrightarrow i(x+yi) = -y + ix
$$

このように、$\tau_2$は虚数単位$i$の性質を表現する行列であると言えます。

# 2次正方実行列の基底

単位行列$I$と実パウリ行列$\tau_1, \tau_2, \tau_3$を適切に線形結合することで、任意の2次正方実行列を表現できます。

&&&prop 実行列空間の基底
$\{I, \tau_1, \tau_2, \tau_3\}$は2次正方実行列の空間$M_2(\mathbb R)$の**基底**をなす。
&&&

&&&prf
$(i, j)$ 成分のみが 1 で、それ以外の成分が 0 であるような行列を$E_{ij}\ (i,j\in\{1,2\})$と表記して、実パウリ行列の線形結合で表す。

$$
\begin{aligned}
E_{11} &= \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} = \frac{1}{2}(I + \tau_3) \\
E_{12} &= \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} = \frac{1}{2}(\tau_1 - \tau_2) \\
E_{21} &= \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix} = \frac{1}{2}(\tau_1 + \tau_2) \\
E_{22} &= \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix} = \frac{1}{2}(I - \tau_3)
\end{aligned}
$$

任意の2次正方実行列は、これらの線形結合として表せる。$p,q,r,s$ を実数とすれば、次のようになる。

$$
\begin{aligned}
\begin{pmatrix} p & q \\ r & s \end{pmatrix}
&= pE_{11} + qE_{12} + rE_{21} + sE_{22} \\
&= \frac p2(I+\tau_3) + \frac q2(\tau_1-\tau_2) + \frac r2(\tau_1+\tau_2) + \frac s2(I-\tau_3) \\
&= \frac{p+s}2 I + \frac{q+r}2 \tau_1 + \frac{-q+r}2 \tau_2 + \frac{p-s}2 \tau_3
\end{aligned}
$$

したがって$\{I, \tau_1, \tau_2, \tau_3\}$は$M_2(\mathbb R)$を張る。$\dim_{\mathbb R}M_2(\mathbb R)=4$であり、この組は4個の行列からなるため、基底をなす。
&&&

&&&rem 生成元
ベクトル空間の基底には4個の行列が必要ですが、積も用いれば、$\tau_1\tau_2=\tau_3,\ \tau_1^2=I$より$\tau_1,\tau_2$の2個から$\tau_3,I$が生成され、$M_2(\mathbb R)$全体が得られます。
&&&

# トレース

実パウリ行列のトレースは$0$です。

$$
\begin{aligned}
\operatorname{tr}(\tau_1)&=\operatorname{tr}\pmatrix{0&1\\1&0}=0 \\
\operatorname{tr}(\tau_2)&=\operatorname{tr}\pmatrix{0&-1\\1&0}=0 \\
\operatorname{tr}(\tau_3)&=\operatorname{tr}\pmatrix{1&0\\0&-1}=0
\end{aligned}
$$

これにより、2次正方実行列$A$を単位行列と実パウリ行列の線形結合で表せば、トレースには単位行列の項しか寄与しないことが分かります。

$$
\begin{aligned}
\operatorname{tr}(A)
&=\operatorname{tr}(a_0I+a_1\tau_1+a_2\tau_2+a_3\tau_3) \\
&=\operatorname{tr}(a_0I) \\
&=\operatorname{tr}\pmatrix{a_0&0\\0&a_0} \\
&=2a_0
\end{aligned}
$$

つまり、トレースには単位行列の項だけが寄与し、その値を2で割ることで係数$a_0$を取り出せます。

## 内積

2次正方実行列の成分を並べたベクトルの内積は、トレースで書けます。

&&&def フロベニウス内積
実行列$A,B\in M_2(\mathbb R)$に対して、次の量を**フロベニウス内積**と呼ぶ。
$$
\langle A,B\rangle_F
:=\operatorname{tr}(A^\mathsf TB)
=\sum_{u,v=1}^2 A_{uv}B_{uv}
$$
&&&

2次正方実行列$A,B$を単位行列と実パウリ行列の線形結合で表します。係数は実数とします。

$$
\begin{aligned}
A&=a_0I+a_1\tau_1+a_2\tau_2+a_3\tau_3=\pmatrix{a_0+a_3&a_1-a_2\\a_1+a_2&a_0-a_3} \\
B&=b_0I+b_1\tau_1+b_2\tau_2+b_3\tau_3=\pmatrix{b_0+b_3&b_1-b_2\\b_1+b_2&b_0-b_3}
\end{aligned}
$$

$A$を転置すれば、係数$a_2$の符号を反転することになります。

$$
A^\mathsf{T}=\pmatrix{a_0+a_3&a_1+a_2\\a_1-a_2&a_0-a_3}=a_0I+a_1\tau_1-a_2\tau_2+a_3\tau_3
$$

$A^\mathsf{T} B$ のトレースを取ります。異なる実パウリ行列の積は別の実パウリ行列の符号倍になり、単位行列と実パウリ行列の積も実パウリ行列なので、トレースは$0$です。したがって、展開した項のうち、同じ基底行列どうしの積だけがトレースに寄与します。実パウリ行列の2乗により単位行列が現れるのに注意します。

$$
\begin{aligned}
\operatorname{tr}(A^\mathsf{T} B)
&=\operatorname{tr}\left[(a_0I+a_1\tau_1-a_2\tau_2+a_3\tau_3)(b_0I+b_1\tau_1+b_2\tau_2+b_3\tau_3)\right] \\
&=\operatorname{tr}(a_0b_0I^2+a_1b_1\tau_1^2-a_2b_2\tau_2^2+a_3b_3\tau_3^2) \\
&=\operatorname{tr}(a_0b_0I+a_1b_1I+a_2b_2I+a_3b_3I) \\
&=2(a_0b_0+a_1b_1+a_2b_2+a_3b_3)
\end{aligned}
$$

同じ添え字の係数を掛けて、それらを足した値の2倍が得られました。これは基底$\{I,\tau_1,\tau_2,\tau_3\}$に関する係数をベクトルと見なしたときの内積の2倍に相当します。

&&&fml トレースと内積
$$
\langle A,B\rangle_F=\operatorname{tr}(A^\mathsf TB)=2(a_0b_0+a_1b_1+a_2b_2+a_3b_3)
$$
&&&

この基底はフロベニウス内積について直交基底ですが、各基底行列のノルムの2乗は$2$なので、正規直交基底ではありません。各基底行列を$\sqrt2$で割ると、正規直交基底になります。係数の内積が2倍になるのはこのためです。

直交性から、各係数を内積で取り出せます。

&&&fml 係数の抽出
$A=a_0I+a_1\tau_1+a_2\tau_2+a_3\tau_3$に対して、次が成り立つ。
$$
a_0=\frac12\operatorname{tr}(A),\qquad
a_r=\frac12\operatorname{tr}(\tau_r^\mathsf{T}A)
\quad(r=1,2,3)
$$
特に$\tau_2^\mathsf T=-\tau_2$なので、$a_2=-\frac12\operatorname{tr}(\tau_2A)$である。転置を付けることで、3つの係数の抽出式を同じ形にできる。
&&&

&&&rem ベクトルの内積
縦ベクトル$\mathbf a, \mathbf b$の内積が転置との積で与えられることに類似しています。

$$
\mathbf a \cdot \mathbf b = {\mathbf a}^\mathsf{T} \mathbf b
$$
&&&

# まとめ

本記事では、パウリ行列の$\sigma_2$に$-i$を掛けて、すべて実行列となる実パウリ行列$\tau_1,\tau_2,\tau_3$を構成し、単位行列と合わせた$\{I,\tau_1,\tau_2,\tau_3\}$が2次正方実行列の基底をなすことを確認しました。実パウリ行列のトレースはすべて$0$なので、トレースを2で割ると単位行列の係数が得られます。転置との積のトレース$\operatorname{tr}(A^\mathsf TB)$（フロベニウス内積）は、基底$\{I,\tau_1,\tau_2,\tau_3\}$に関する係数ベクトルの標準内積の2倍になり、この基底は直交基底をなします。

&&& 実パウリ行列
$$
\tau_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
\tau_2 = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \quad
\tau_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$
これらは次の関係を満たす。
$$
\tau_1^2 = I, \quad \tau_2^2 = -I, \quad \tau_3^2 = I
$$
$$
\tau_1\tau_2 = -\tau_2\tau_1 = \tau_3, \quad
\tau_1\tau_3 = -\tau_3\tau_1 = \tau_2, \quad
\tau_2\tau_3 = -\tau_3\tau_2 = \tau_1
$$
&&&

&&& トレースと内積
$A = a_0 I + a_1 \tau_1 + a_2 \tau_2 + a_3 \tau_3,\ B = b_0 I + b_1 \tau_1 + b_2 \tau_2 + b_3 \tau_3$とすると、次のようになる。
$$
\operatorname{tr}(A) = 2a_0, \quad
\operatorname{tr}(A^\mathsf{T}B) = 2(a_0b_0 + a_1b_1 + a_2b_2 + a_3b_3)
$$
&&&

&&& 係数の抽出
$$
a_0=\frac12\operatorname{tr}(A),\quad
a_r=\frac12\operatorname{tr}(\tau_r^\mathsf{T}A)\quad(r=1,2,3)
$$
&&&
