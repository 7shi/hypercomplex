本記事では、量子力学でよく知られたパウリ行列から、成分が全て実数である**実パウリ行列**と呼ばれる行列の組を定義します。この実パウリ行列と単位行列は二次正方実行列の基底をなします。この基底で表したときの係数によって**トレース**を解釈します。

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

$-i$を複素行列である$\sigma_2$に作用させます。

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

$$
\begin{alignedat}{2}
\tau_1 &:= &\sigma_1 &= \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \\
\tau_2 &:= &-i\sigma_2 &= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \\
\tau_3 &:= &\sigma_3 &= \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
\end{alignedat}
$$

この定義では、$\tau_1, \tau_2, \tau_3$ は全て実行列であり、以下の関係式を満たします。

$$
\begin{aligned}
\tau_1^2 = I&,\ \tau_2^2=-I,\ \tau_3^2=I \\
\tau_1\tau_2 &= -\tau_2\tau_1 = \tau_3 \\
\tau_1\tau_3 &= -\tau_3\tau_1 = \tau_2 \\
\tau_2\tau_3 &= -\tau_3\tau_2 = \tau_1
\end{aligned}
$$

##$\tau_2$の持つ意味

この定義で得られた$\tau_2$は、虚数単位$i$の$2 \times 2$行列表現と見なすことができます。実際、$\tau_2$の二乗は$i^2 = -1$に対応します。

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

# 二次正方実行列の基底

単位行列$I$と実パウリ行列$\tau_1, \tau_2, \tau_3$を適切に線形結合することで、任意の二次正方実行列の各成分を表現できます。

まず、$(i, j)$ 成分のみが 1 で、それ以外の成分が 0 であるような行列を$E_{ij}$と表記して、実パウリ行列の線形結合で表します。

$$
\begin{aligned}
E_{11} &= \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} = \frac{1}{2}(I + \tau_3) \\
E_{12} &= \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} = \frac{1}{2}(\tau_1 - \tau_2) \\
E_{21} &= \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix} = \frac{1}{2}(\tau_1 + \tau_2) \\
E_{22} &= \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix} = \frac{1}{2}(I - \tau_3)
\end{aligned}
$$

任意の二次正方実行列は、これらの線形結合として表せます。$p,q,r,s$ を実数とすれば

$$
\begin{aligned}
\begin{pmatrix} p & q \\ r & s \end{pmatrix}
&= pE_{11} + qE_{12} + rE_{21} + sE_{22} \\
&= \frac p2(I+\tau_3) + \frac q2(\tau_1-\tau_2) + \frac r2(\tau_1+\tau_2) + \frac s2(I-\tau_3) \\
&= \frac{p+s}2 I + \frac{q+r}2 \tau_1 + \frac{-q+r}2 \tau_2 + \frac{p-s}2 \tau_3
\end{aligned}
$$

このことは$\{I, \tau_1, \tau_2, \tau_3\}$が二次正方実行列の**基底**をなすことを示しています。

&&&rem 生成元
$\tau_1\tau_2=\tau_3,\ \tau_1^2=I$より、$\tau_1,\tau_2$から$\tau_3,I$が生成されます。
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

これにより、二次正方実行列$A$を単位行列と実パウリ行列の線形結合で表せば、トレースには単位行列の項しか寄与しないことが分かります。

$$
\begin{aligned}
\operatorname{tr}(A)
&=\operatorname{tr}(a_0I+a_1\tau_1+a_2\tau_2+a_3\tau_3) \\
&=\operatorname{tr}(a_0I) \\
&=\operatorname{tr}\pmatrix{a_0&0\\0&a_0} \\
&=2a_0
\end{aligned}
$$

つまり、トレースを取ることで、二次正方実行列から単位行列の成分を抽出することになります。

## 内積

二次正方実行列$A,B$を単位行列と実パウリ行列の線形結合で表します。（係数は実数）

$$
\begin{aligned}
A&=a_0I+a_1\tau_1+a_2\tau_2+a_3\tau_3=\pmatrix{a_0+a_3&a_1-a_2\\a_1+a_2&a_0-a_3} \\
B&=b_0I+b_1\tau_1+b_2\tau_2+b_3\tau_3=\pmatrix{b_0+b_3&b_1-b_2\\b_1+b_2&b_0-b_3}
\end{aligned}
$$

$A$を転置すれば、$a_2$を符号反転することになります。

$$
A^\mathsf{T}=\pmatrix{a_0+a_3&a_1+a_2\\a_1-a_2&a_0-a_3}=a_0I+a_1\tau_1-a_2\tau_2+a_3\tau_3
$$

$A^\mathsf{T} B$ のトレースを取ります。実パウリ行列の2乗により単位行列が現れるのに注意します。

$$
\begin{aligned}
\operatorname{tr}(A^\mathsf{T} B)
&=\operatorname{tr}\left[(a_0I+a_1\tau_1-a_2\tau_2+a_3\tau_3)(b_0I+b_1\tau_1+b_2\tau_2+b_3\tau_3)\right] \\
&=\operatorname{tr}(a_0b_0I^2+a_1b_1\tau_1^2-a_2b_2\tau_2^2+a_3b_3\tau_3^2) \\
&=\operatorname{tr}(a_0b_0I+a_1b_1I+a_2b_2I+a_3b_3I) \\
&=2(a_0b_0+a_1b_1+a_2b_2+a_3b_3)
\end{aligned}
$$

同じ添え字の係数を掛けて、それらを足した値の2倍が得られました。これは実パウリ行列の係数をベクトルと見なしたときの内積の2倍に相当します。

&&&rem ベクトルの内積
縦ベクトル$\mathbf a, \mathbf b$の内積が転置との積で与えられることに類似しています。

$$
\mathbf a \cdot \mathbf b = {\mathbf a}^\mathsf{T} \mathbf b
$$
&&&

# まとめ

二次正方実行列を単位行列と実パウリ行列の線形結合で表すことで、積のトレースが内積（の2倍）に相当すると解釈できます。
