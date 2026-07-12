前回[[7shi-s]]、状態ベクトルとパウリ行列からブロッホベクトルと密度行列を導出しました。今回はパウリ行列によって何が表現されるかという視点から、純粋状態と混合状態を説明します。量子情報への接続は予備知識程度にとどめます。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

# エルミート行列の基底

パウリ行列を再掲します。

&&&def パウリ行列
$$
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix}
,\ \sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix}
,\ \sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}
$$
&&&

エルミート行列（エルミート共役を取っても変わらない行列 $A^\dagger=A$）の一般形を考えます。対角成分は実数で、非対角成分は互いに複素共役です。実数$a,d,x,y$を使って成分を書き下します。

$$
A=\begin{pmatrix}a&x-iy\\x+iy&d\end{pmatrix}
$$

これを単位行列とパウリ行列に分解します。

$$
A
=\frac{a+d}2\begin{pmatrix}1&0\\0&1\end{pmatrix}
+x\begin{pmatrix}0&1\\1&0\end{pmatrix}
+y\begin{pmatrix}0&-i\\i&0\end{pmatrix}
+\frac{a-d}2\begin{pmatrix}1&0\\0&-1\end{pmatrix}
$$

係数を$w=\frac{a+d}2,\ z=\frac{a-d}2$と置き直せば、任意のエルミート行列が単位行列とパウリ行列の実係数線形結合で一意に表せることが分かります。

&&&fml エルミート行列の展開
$$
A=w\,\mathrm I+x\,\sigma_x+y\,\sigma_y+z\,\sigma_z\quad(w,x,y,z\in\mathbb R)
$$
&&&

逆に、この形の行列は常にエルミートです。つまり単位行列とパウリ行列は、$2\times2$エルミート行列がなす実 4 次元空間の基底です。パウリ行列の係数$(x,y,z)$が 3 次元の座標を与えるというのが、前回までの計算の背景にある構造です。

## トレースによる係数の抽出

行列から係数を取り出すには、トレース（対角成分の和）を使います。準備としてパウリ行列の積を確認します。

&&&fml パウリ行列の積
$$
\sigma_x^2=\sigma_y^2=\sigma_z^2=\mathrm I
,\quad
\sigma_x\sigma_y=i\sigma_z
,\quad
\sigma_y\sigma_z=i\sigma_x
,\quad
\sigma_z\sigma_x=i\sigma_y
$$
順序を入れ替えると符号が反転します（反交換性）。
$$
\sigma_y\sigma_x=-i\sigma_z
,\quad
\sigma_z\sigma_y=-i\sigma_x
,\quad
\sigma_x\sigma_z=-i\sigma_y
$$
&&&

&&&rem 四元数との関係
$-i$を掛けた$-i\sigma_x,\ -i\sigma_y,\ -i\sigma_z$は$\mathbf i,\mathbf j,\mathbf k$と同じ積の規則に従います。例えば$(-i\sigma_x)(-i\sigma_y)=-\sigma_x\sigma_y=-i\sigma_z$は$\mathbf{ij}=\mathbf k$に、$(-i\sigma_x)^2=-\sigma_x^2=-\mathrm I$は$\mathbf i^2=-1$に対応します。
&&&

パウリ行列はトレースがゼロです。積のトレースは、同じもの同士では$\mathrm{tr}\,\mathrm I=2$、異なるもの同士ではパウリ行列の$\pm i$倍となるためゼロです。

$$
\mathrm{tr}\,\sigma_x=\mathrm{tr}\,\sigma_y=\mathrm{tr}\,\sigma_z=0
,\quad
\mathrm{tr}(\sigma_i\sigma_j)=2\delta_{ij}\quad(i,j=x,y,z)
$$

これを使えば、展開した行列の右からパウリ行列を掛けてトレースを取ることで、対応する係数だけが残ります。

$$
\mathrm{tr}(A\sigma_x)
=\mathrm{tr}(w\,\sigma_x+x\,\mathrm I-iy\,\sigma_z+iz\,\sigma_y)
=2x
$$

&&&fml 係数の抽出
$$
w=\frac12\mathrm{tr}\,A
,\quad
x=\frac12\mathrm{tr}(A\sigma_x)
,\quad
y=\frac12\mathrm{tr}(A\sigma_y)
,\quad
z=\frac12\mathrm{tr}(A\sigma_z)
$$
&&&

# 純粋状態

前回導出した密度行列を、この視点から見直します。状態ベクトル$ω$から作った密度行列はエルミート行列です。

$$
ρ=ωω^\dagger
=\begin{pmatrix}α\\β\end{pmatrix}\begin{pmatrix}α^*&β^*\end{pmatrix}
=\begin{pmatrix}αα^*&αβ^*\\βα^*&ββ^*\end{pmatrix}
\quad(|α|^2+|β|^2=1)
$$

係数を抽出します。

$$
\begin{alignedat}{2}
\frac12\mathrm{tr}\,ρ&=\frac12(αα^*+ββ^*)&&=\frac12 \\
\frac12\mathrm{tr}(ρ\sigma_x)&=\frac12(αβ^*+βα^*)&&=\frac12(α^*β+β^*α) \\
\frac12\mathrm{tr}(ρ\sigma_y)&=\frac i2(αβ^*-βα^*)&&=-\frac i2(α^*β-β^*α) \\
\frac12\mathrm{tr}(ρ\sigma_z)&=\frac12(αα^*-ββ^*)&&=\frac12(α^*α-β^*β)
\end{alignedat}
$$

パウリ行列の係数として、初回[[7shi-h]]に$ω\mathbf kω^*$の係数として求めた$x,y,z$の$\frac12$倍が現れました。

&&&fml 純粋状態の密度行列
$$
ρ=ωω^\dagger=\frac12(\mathrm I+x\,\sigma_x+y\,\sigma_y+z\,\sigma_z)
$$
$(x,y,z)$はブロッホベクトル
$$
x=α^*β+β^*α
,\quad
y=-i(α^*β-β^*α)
,\quad
z=α^*α-β^*β
$$
&&&

以降の計算のため、ベクトルの内積の形で略記します。

&&&def 内積表記
$$
\boldsymbol r=(x,y,z)
,\quad
\boldsymbol σ=(\sigma_x,\sigma_y,\sigma_z)
,\quad
\boldsymbol r\cdot\boldsymbol σ=x\,\sigma_x+y\,\sigma_y+z\,\sigma_z
$$
$$
ρ=\frac12(\mathrm I+\boldsymbol r\cdot\boldsymbol σ)
$$
&&&

このように 1 つの状態ベクトルから作られる密度行列を**純粋状態**と呼びます。

## べき等性

純粋状態の密度行列は 2 乗しても変わりません。$ω^\daggerω=|α|^2+|β|^2=1$がスカラーであることから直ちに従います。

$$
ρ^2=(ωω^\dagger)(ωω^\dagger)=ω(ω^\daggerω)ω^\dagger=ωω^\dagger=ρ
$$

同じことをブロッホベクトルの側から見ます。反交換性により交差項が消えるため、$\boldsymbol r\cdot\boldsymbol σ$の 2 乗はスカラー倍になります。

$$
\begin{aligned}
(\boldsymbol r\cdot\boldsymbol σ)^2
&=x^2\sigma_x^2+y^2\sigma_y^2+z^2\sigma_z^2 \\
&\quad+xy(\sigma_x\sigma_y+\sigma_y\sigma_x)+yz(\sigma_y\sigma_z+\sigma_z\sigma_y)+zx(\sigma_z\sigma_x+\sigma_x\sigma_z) \\
&=(x^2+y^2+z^2)\,\mathrm I \\
&=|\boldsymbol r|^2\,\mathrm I
\end{aligned}
$$

これを使って$ρ^2$を計算します。

$$
ρ^2
=\frac14(\mathrm I+\boldsymbol r\cdot\boldsymbol σ)^2
=\frac14\{(1+|\boldsymbol r|^2)\,\mathrm I+2\,\boldsymbol r\cdot\boldsymbol σ\}
$$

$ρ^2=ρ$となる条件は$\frac14(1+|\boldsymbol r|^2)=\frac12$、つまり$|\boldsymbol r|=1$です。べき等性は、ブロッホベクトルが単位球面上にあることと同値です。初回に$ω\mathbf kω^*$の係数の 2 乗和が$1$となることを確認しましたが、それをパウリ行列の言葉で言い直したことになります。

## 位相の消去

状態ベクトルに位相$e^{ic}$（絶対値$1$の複素数）を掛けても、密度行列は変わりません。

$$
(ωe^{ic})(ωe^{ic})^\dagger=ω\,e^{ic}e^{-ic}\,ω^\dagger=ωω^\dagger=ρ
$$

&&&rem ホップファイブレーションとの関係
状態ベクトル$ω$は$S^3$上の点で、位相$e^{ic}$の自由度がファイバーの円周に対応します。密度行列ではこの円周が最初から潰れており、$ρ$はファイバーを区別せず底空間$S^2$（ブロッホ球面）の点を直接表します。つまり純粋状態の密度行列とは、ホップファイブレーションの像を行列で表現したものです。
&&&

# 混合状態

複数の純粋状態を古典的な確率で混ぜた状態を考えます。重ね合わせ（状態ベクトルの和）ではなく、確率$p_i$で状態$ω_i$が用意されているという状況です。

&&&def 混合状態の密度行列
$$
ρ=\sum_i p_i\,ω_iω_i^\dagger\quad\Bigl(p_i>0,\ \sum_i p_i=1\Bigr)
$$
&&&

各項をブロッホベクトル$\boldsymbol r_i$で展開すれば、和が単位行列の部分とパウリ行列の部分に分かれます。

$$
ρ=\sum_i p_i\cdot\frac12(\mathrm I+\boldsymbol r_i\cdot\boldsymbol σ)
=\frac12\Bigl(\mathrm I+\sum_i p_i\boldsymbol r_i\cdot\boldsymbol σ\Bigr)
$$

つまり混合状態のブロッホベクトルは、各純粋状態のブロッホベクトルの重み付き平均です。

$$
\boldsymbol r=\sum_i p_i\boldsymbol r_i
$$

単位ベクトルの重み付き平均は、方向が異なれば長さが$1$未満になります（三角不等式）。純粋状態がブロッホ球の表面を指すのに対して、混合状態は内部を指します。[[quantumuniverse]]

&&&ex 混合と重ね合わせ
基底となる 2 つの状態を用意します（量子情報の記法では$|0\rangle,|1\rangle$と書かれます）。

$$
ω_0=\begin{pmatrix}1\\0\end{pmatrix}
,\quad
ω_1=\begin{pmatrix}0\\1\end{pmatrix}
$$

ブロッホベクトルはそれぞれ$(0,0,1),(0,0,-1)$、つまり北極と南極です。等確率$\frac12$ずつで混合します。

$$
ρ=\frac12ω_0ω_0^\dagger+\frac12ω_1ω_1^\dagger
=\frac12\begin{pmatrix}1&0\\0&1\end{pmatrix}
,\quad
\boldsymbol r=(0,0,0)
$$

ブロッホベクトルは原点（球の中心）を指します。これを**最大混合状態**と呼びます。

一方、等係数の重ね合わせは 1 本の状態ベクトルなので純粋状態です。

$$
ω_+=\frac1{\sqrt2}\begin{pmatrix}1\\1\end{pmatrix}
,\quad
ρ_+=ω_+ω_+^\dagger=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}=\frac12(\mathrm I+\sigma_x)
,\quad
\boldsymbol r=(1,0,0)
$$

対角成分（$z$成分の情報）はどちらも同じですが、非対角成分（$x,y$成分）が異なります。混合でも重ね合わせでも「半々」であることは共通していて、その違いはパウリ行列の係数の違いとして現れます。
&&&

## 純粋度

べき等性の計算で使った式は、混合状態でもそのまま成り立ちます。

$$
ρ^2=\frac14\{(1+|\boldsymbol r|^2)\,\mathrm I+2\,\boldsymbol r\cdot\boldsymbol σ\}
$$

トレースを取れば、パウリ行列の部分が消えます。

&&&fml 純粋度
$$
\mathrm{tr}\,ρ^2=\frac12(1+|\boldsymbol r|^2)
$$
&&&

球の表面$|\boldsymbol r|=1$で最大値$1$（純粋状態）、中心で最小値$\frac12$（最大混合状態）となります。ブロッホベクトルの長さは、状態がどれだけ純粋かを測る指標です。

## 固有分解

$(\boldsymbol r\cdot\boldsymbol σ)^2=|\boldsymbol r|^2\,\mathrm I$かつ$\mathrm{tr}(\boldsymbol r\cdot\boldsymbol σ)=0$より、$\boldsymbol r\cdot\boldsymbol σ$の固有値は$\pm|\boldsymbol r|$です。したがって$ρ=\frac12(\mathrm I+\boldsymbol r\cdot\boldsymbol σ)$の固有値は次のようになります。

$$
λ_\pm=\frac12(1\pm|\boldsymbol r|)
$$

対応する固有ベクトルは、$\pm\dfrac{\boldsymbol r}{|\boldsymbol r|}$方向のブロッホベクトルを持つ純粋状態、つまり球面上の対蹠点のペアです。これらを$ρ_\pm$とすれば、$\boldsymbol r\ne\boldsymbol 0$の混合状態は直交する 2 つの純粋状態の混合として書けます。

$$
ρ_\pm=\frac12\Bigl(\mathrm I\pm\frac{\boldsymbol r}{|\boldsymbol r|}\cdot\boldsymbol σ\Bigr)
,\quad
λ_+ρ_++λ_-ρ_-
=\frac12\{(λ_++λ_-)\,\mathrm I+(λ_+-λ_-)\frac{\boldsymbol r}{|\boldsymbol r|}\cdot\boldsymbol σ\}
=ρ
$$

&&&rem 分解の非一意性
同じ密度行列を与える混合の作り方は 1 通りではありません。例えば最大混合状態$\frac12\mathrm I$は、北極と南極の等確率混合としても、$x$軸上の対蹠点$ω_\pm=\frac1{\sqrt2}\begin{pmatrix}1\\\pm1\end{pmatrix}$の等確率混合としても得られます。密度行列は「どう混ぜたか」の情報を持たず、統計的に区別できない混合を同じ 1 点で表します。
&&&

# 測定との関係

量子情報では、エルミート行列$A$で表される観測量の期待値が密度行列から計算されます。

$$
\langle A\rangle=\mathrm{tr}(ρA)
$$

純粋状態では$\mathrm{tr}(ωω^\dagger A)=ω^\dagger Aω$となり、前回計算した$ω^\dagger\sigma_xω$などと同じ形です。特にパウリ行列を観測量とすれば、期待値はブロッホベクトルの成分そのものです。

$$
\langle\sigma_x\rangle=\mathrm{tr}(ρ\sigma_x)=x
,\quad
\langle\sigma_y\rangle=\mathrm{tr}(ρ\sigma_y)=y
,\quad
\langle\sigma_z\rangle=\mathrm{tr}(ρ\sigma_z)=z
$$

また、$ω_0$と$ω_1$のどちらであるかを判別する測定（$z$測定）の確率は、対角成分から得られます。

$$
p_0=\frac{1+z}2
,\quad
p_1=\frac{1-z}2
$$

ブロッホベクトルとは、パウリ行列で表される観測量の期待値を並べた座標だと言えます。

# まとめ

単位行列とパウリ行列は$2\times2$エルミート行列の基底で、パウリ行列の係数が 3 次元の座標を与えます。純粋状態の密度行列はブロッホ球の表面を指し、位相（ファイバー）が潰れることから、ホップファイブレーションの像の行列表現とみなせます。混合状態はブロッホベクトルの重み付き平均としてブロッホ球の内部に広がり、ファイブレーションの枠組みの外側の世界を表します。

&&&rem
パウリ行列の係数として取り出される情報（位相差、混合の度合い）が観測量の期待値と一致するという点に、この表現が量子情報で使われる理由があります。
&&&
