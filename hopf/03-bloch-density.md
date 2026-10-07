2×2エルミート行列の基底としてのパウリ行列の性質を整理し、ブロッホ球における純粋状態と混合状態の幾何学的構造および量子測定との関係を解き明かします。

シリーズ：[ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

&&& 改訂履歴
- 2026.10.04 混合状態の定義を異なる純粋状態の確率的な混合に限定、密度行列の条件と$|\boldsymbol r|\le1$の対応を追加、測定と混合の説明を修正
&&&

# 概要

前回は、状態ベクトル（スピノル）とそのエルミート共役のテンソル積からブロッホベクトルと密度行列を導出しました。[[7shi-s]]

この表現は量子情報において中心的な役割を果たしますが、なぜ任意の密度行列がパウリ行列で展開できるのか、またブロッホ球の表面にある純粋状態に対して、球の「内部」に位置する混合状態がどのような代数的・幾何学的意味を持つのかという全体像を整理する問いが残されています。

本記事では、まず$2\times2$エルミート行列全体が単位行列とパウリ行列の実係数線形結合で一意に表せること、およびパウリ行列の直交性（トレース）から各係数を抽出できることを確認します。次に、状態ベクトルから作られる純粋状態の密度行列が冪等性（$ρ^2=ρ$）と単位長（$|\boldsymbol r|=1$）を持つことを見直し、複数の純粋状態を確率的に混ぜた混合状態がブロッホベクトルの重み付き平均としてブロッホ球の内部（$|\boldsymbol r|<1$）に対応することを明らかにします。さらに、古典的な確率混合と量子的な重ね合わせの違いを具体例で比較し、最後に任意の方向への測定の期待値と観測確率がブロッホベクトルの射影として簡明に求まる構造を示します。

前回までの四元数によるホップファイブレーション、状態ベクトル、パウリ行列、および密度行列の定義を前提とします。本記事では1量子ビット系の密度行列とブロッホ球の代数的・幾何学的な側面に集中し、複数量子ビットのエンタングルメントやPOVMなどの一般化測定理論は扱いません。

# エルミート行列の基底

パウリ行列を再掲します。

&&&def パウリ行列
$$
σ_x=\begin{pmatrix}0&1\\1&0\end{pmatrix}, \quad
σ_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix}, \quad
σ_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}
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

係数を$w=\dfrac{a+d}2,\ z=\dfrac{a-d}2$と置き直せば、任意のエルミート行列が単位行列とパウリ行列の実係数線形結合で一意に表せることが分かります。

&&&fml エルミート行列のパウリ展開
$$
A=w\,\mathrm I+x\,σ_x+y\,σ_y+z\,σ_z\quad(w,x,y,z\in\mathbb R)
$$
&&&

逆に、この形の行列は常にエルミートです。つまり単位行列とパウリ行列は、$2\times2$エルミート行列がなす実4次元空間の基底です。パウリ行列の係数$(x,y,z)$が3次元の座標を与えるというのが、前回までの計算の背景にある構造です。

&&&rem ベクトルの行列表現
単位行列$\mathrm I$がスカラー$1$の行列表現であるのに対し、パウリ行列$σ_x,σ_y,σ_z$は3次元空間の基底ベクトル$(1,0,0),(0,1,0),(0,0,1)$の行列表現とみなせます。線形結合$x\,σ_x+y\,σ_y+z\,σ_z$は、この対応のもとでベクトル$(x,y,z)$を行列として表現したものです。

ブロッホベクトルは、この対応によりパウリ行列の線形結合と同一視されたベクトルです。
&&&

## パウリ行列の積

パウリ行列の積を確認します。

&&&fml パウリ行列の積
$$
σ_x^2=σ_y^2=σ_z^2=\mathrm I
,\quad
σ_xσ_y=iσ_z
,\quad
σ_yσ_z=iσ_x
,\quad
σ_zσ_x=iσ_y
$$
順序を入れ替えると符号が反転する（反交換性）。
$$
σ_yσ_x=-iσ_z
,\quad
σ_zσ_y=-iσ_x
,\quad
σ_xσ_z=-iσ_y
$$
&&&

&&&rem 四元数との関係
$-i$を掛けた$-iσ_x,\ -iσ_y,\ -iσ_z$は$\mathbf i,\mathbf j,\mathbf k$と同じ積の規則に従います。
$$
(-iσ_x)(-iσ_y)=-σ_xσ_y=-iσ_z
\quad\Leftrightarrow\quad
\mathbf{ij}=\mathbf k
$$

これは複素数$\mathbb{C}$と四元数$\mathbb{H}$のテンソル積である双四元数$\mathbb{C}\otimes\mathbb{H}$の構造を反映しています。$h$は$\mathbf i,\mathbf j,\mathbf k$と可換で$h^2=-1$を満たす、双四元数の複素数部分の虚数単位です。[[7shi-bq]]
$$
h \cong σ_xσ_yσ_z=i\mathrm I,\quad h^2=-1
$$
$$
h\mathbf i \cong i(-iσ_x) = σ_x, \quad
h\mathbf j \cong i(-iσ_y) = σ_y, \quad
h\mathbf k \cong i(-iσ_z) = σ_z
$$
$$
(h\mathbf i)(h\mathbf j)=h(h\mathbf k)
\quad\Leftrightarrow\quad
σ_xσ_y=iσ_z
$$
&&&

## トレースによる係数の抽出

パウリ行列はトレースがゼロです。積のトレースは、同じもの同士では$\mathrm{tr}\,\mathrm I=2$、異なるもの同士ではパウリ行列の$\pm i$倍となるためゼロです。

$$
\mathrm{tr}\,σ_x=\mathrm{tr}\,σ_y=\mathrm{tr}\,σ_z=0
,\quad
\mathrm{tr}(σ_iσ_j)=2\delta_{ij}\quad(i,j=x,y,z)
$$

行列から係数を取り出すには、トレース（対角成分の和）を使用します。上で確認した積とトレースの関係を使えば、展開した行列の右からパウリ行列を掛けてトレースを取ることで、対応する係数だけが残ります。

$$
\mathrm{tr}(Aσ_x)
=\mathrm{tr}(w\,σ_x+x\,\mathrm I-iy\,σ_z+iz\,σ_y)
=2x
$$

&&&fml 係数の抽出
$$
w=\frac12\mathrm{tr}\,A, \quad
x=\frac12\mathrm{tr}(Aσ_x), \quad
y=\frac12\mathrm{tr}(Aσ_y), \quad
z=\frac12\mathrm{tr}(Aσ_z)
$$
&&&

&&&rem 双四元数におけるトレースの対応
双四元数において、トレースは四元数としてのスカラー部$\operatorname{Sc}$（$1$の係数で、$h$を含みうる）の2倍に対応します。
$$
A = w + xh\mathbf i + yh\mathbf j + zh\mathbf k
$$
$$
2\operatorname{Sc}(Ah\mathbf i)
=2\operatorname{Sc}(wh\mathbf i + x + y\mathbf k - z\mathbf j)
=2x
$$
&&&

# 純粋状態

前回導出した密度行列を、この視点から見直します。状態ベクトル$\Psi$から作った密度行列はエルミート行列です。

$$
ρ=\Psi\Psi^\dagger
=\begin{pmatrix}α\\β\end{pmatrix}\begin{pmatrix}α^*&β^*\end{pmatrix}
=\begin{pmatrix}αα^*&αβ^*\\βα^*&ββ^*\end{pmatrix}
\quad(|α|^2+|β|^2=1)
$$

前回見たとおり$ρ=\frac12(\mathrm I+x\,σ_x+y\,σ_y+z\,σ_z)$なので、$ρ$のパウリ展開の係数はブロッホベクトルの各成分の半分です。$\mathrm{tr}(σ_iσ_j)=2\delta_{ij}$によってこの$\frac12$が相殺され、トレースを取るとブロッホベクトルの成分をそのまま抽出できます。

$$
\begin{alignedat}{3}
&\mathrm{tr}\,ρ   &&=αα^*+ββ^*&&=1 \\
&\mathrm{tr}(ρσ_x)&&=αβ^*+βα^*&&=α^*β+β^*α \\
&\mathrm{tr}(ρσ_y)&&=i(αβ^*-βα^*)&&=-i(α^*β-β^*α) \\
&\mathrm{tr}(ρσ_z)&&=αα^*-ββ^*&&=α^*α-β^*β
\end{alignedat}
$$

初回の記事で$ω\mathbf kω^*$の係数として求めた$x,y,z$が、$\mathrm{tr}(ρσ_i)$としてそのまま現れます。[[7shi-h]]

&&&fml 純粋状態の密度行列
$$
ρ=\Psi\Psi^\dagger=\frac12(\mathrm I+x\,σ_x+y\,σ_y+z\,σ_z)
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
\boldsymbol σ=(σ_x,σ_y,σ_z)
,\quad
\boldsymbol r\cdot\boldsymbol σ=x\,σ_x+y\,σ_y+z\,σ_z
$$
$$
ρ=\frac12(\mathrm I+\boldsymbol r\cdot\boldsymbol σ)
$$
&&&

&&&rem 記法の比較
$\boldsymbol r\cdot\boldsymbol σ$は、アインシュタインの縮約記法による$r^iσ_i$に対応します。前者は添字を隠した書き方、後者は明示した書き方です。

なお、後で$\boldsymbol r_i$という表記も登場しますが、そちらの添字$i$は成分ではなく、複数のブロッホベクトルを区別する番号です。
&&&

このように1つの状態ベクトルから作られる密度行列を**純粋状態**と呼びます。

## 冪等性

純粋状態の密度行列は2乗しても変わりません。$\Psi^\dagger\Psi=|α|^2+|β|^2=1$がスカラーであることから直ちに従います。

$$
ρ^2=(\Psi\Psi^\dagger)(\Psi\Psi^\dagger)=\Psi(\Psi^\dagger\Psi)\Psi^\dagger=\Psi\Psi^\dagger=ρ
$$

同じことをブロッホベクトルの側から見ます。反交換性により交差項が消えるため、$\boldsymbol r\cdot\boldsymbol σ$の2乗はスカラー倍になります。

$$
\begin{aligned}
(\boldsymbol r\cdot\boldsymbol σ)^2
&=x^2σ_x^2+y^2σ_y^2+z^2σ_z^2 \\
&\quad+xy(σ_xσ_y+σ_yσ_x)+yz(σ_yσ_z+σ_zσ_y)+zx(σ_zσ_x+σ_xσ_z) \\
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

$ρ^2=ρ$となる条件は$\frac14(1+|\boldsymbol r|^2)=\frac12$、つまり$|\boldsymbol r|=1$です。冪等性は、ブロッホベクトルが単位球面上にあることと同値です。初回に$ω\mathbf kω^*$の係数の2乗和が$1$となることを確認しましたが、それをパウリ行列の言葉で言い直したことになります。

&&&rem 単位行列の役割
$ρ$に含まれる単位行列$\mathrm I$の係数$\frac12$は、ブロッホベクトルに追加された自由度ではありません。パウリ行列は無トレースなので、$\mathrm{tr}\,ρ=1$（規格化条件$|α|^2+|β|^2=1$）から常に$\frac12$に固定される値です。状態の向きの情報はすべて$x,y,z$側が担っています。
&&&

## 位相の消去

状態ベクトルに位相$e^{ic}$（絶対値$1$の複素数）を掛けても、密度行列は変わりません。

$$
(\Psi e^{ic})(\Psi e^{ic})^\dagger=\Psi\,e^{ic}e^{-ic}\,\Psi^\dagger=\Psi\Psi^\dagger=ρ
$$

&&&rem ホップファイブレーションとの関係
状態ベクトル$\Psi$は$S^3$上の点で、位相$e^{ic}$の自由度がファイバーの円周に対応します。密度行列ではこの円周が最初から潰れており、$ρ$はファイバーを区別せず底空間$S^2$（ブロッホ球面）の点を直接表します。つまり純粋状態の密度行列とは、ホップファイブレーションの像を行列で表現したものです。
&&&

# 混合状態

複数の純粋状態$\Psi_i$を古典的な確率$p_i$で混ぜた状態を考えます。

&&&def 密度行列の確率混合
$$
ρ=\sum_i p_i\,\Psi_i\Psi_i^\dagger\quad\Bigl(p_i>0,\ \sum_i p_i=1\Bigr)
$$
各$\Psi_i$は規格化されているとする。$\Psi_i\Psi_i^\dagger$に異なる純粋状態が含まれるとき、$ρ$を**混合状態**と呼ぶ。
&&&

すべての$\Psi_i\Psi_i^\dagger$が同じ純粋状態なら、$ρ$はその純粋状態のままです。

各項をブロッホベクトル$\boldsymbol r_i$で展開すれば、和が単位行列の部分とパウリ行列の部分に分かれます。

$$
ρ=\sum_i p_i\cdot\frac12(\mathrm I+\boldsymbol r_i\cdot\boldsymbol σ)
=\frac12\Bigl(\mathrm I+\sum_i p_i\boldsymbol r_i\cdot\boldsymbol σ\Bigr)
$$

つまり混合状態のブロッホベクトルは、各純粋状態のブロッホベクトルの重み付き平均です。

$$
\boldsymbol r=\sum_i p_i\boldsymbol r_i
$$

単位ベクトルの重み付き平均は、方向が異なれば長さが$1$未満になります（三角不等式）。純粋状態がブロッホ球の表面を指すのに対して、混合状態は内部を指します。

&&&rem 重み付き平均としての解釈
最初に述べた「ベクトルの行列表現」に対応させると、混合は単位行列の係数$\frac12$を保ったまま、パウリ行列側の係数だけを$\sum_i p_i\boldsymbol r_i$に置き換える操作です。行列としては複雑に見えても、ベクトル表現では$\boldsymbol r_i$を重み$p_i$で足し合わせているだけであり、上で述べた「長さが$1$未満になる」という性質も、単位ベクトルの加重平均が一般に短くなるという素朴なベクトルの性質がそのまま現れたものです。
&&&

## 混合と重ね合わせ

「2つの状態を半々で組み合わせる」という操作には2通りの意味があります。密度行列を確率で混ぜる**混合**と、状態ベクトルを足し合わせる**重ね合わせ**です。混ぜる対象（行列かベクトルか）が異なるため、結果も異なります。具体例で確認します。

基底となる2つの状態を用意します。対応するブロッホベクトルも示します。

$$
\begin{aligned}
\Psi_\uparrow&=\begin{pmatrix}1\\0\end{pmatrix}, && \boldsymbol r_\uparrow=(0,0,1) \\
\\
\Psi_\downarrow&=\begin{pmatrix}0\\1\end{pmatrix}, && \boldsymbol r_\downarrow=(0,0,-1)
\end{aligned}
$$

&&&rem 量子情報のブラケット記法
量子情報の記法では$|\uparrow\rangle,|\downarrow\rangle$または$|0\rangle,|1\rangle$と書かれます。
&&&

これらのブロッホベクトルはそれぞれ北極と南極を指します。

&&&rem 北極が基準になる理由
$\Psi_\uparrow$が北極に対応するのは偶然ではありません。ホップファイブレーションは初回の記事において、$\mathbf k$を回転の始点として$ω\mathbf kω^*$と定義されていました。恒等回転$ω=1$（$ω_0=1$，他は$0$）は始点$\mathbf k$（北極）をそのまま返します。初回の対応で四元数から直接作る複素数の組は$(i,0)^T$で、前回と同様に共通の位相$i$を除けば$\Psi_\uparrow=(1,0)^T$になります。つまり$\Psi_\uparrow$が北極を指すのは、$\mathbf k$を起点に選んだホップファイブレーションの定義がそのまま反映された結果です。[[7shi-h]]
&&&

### 混合

密度行列そのものを確率$\frac12$ずつで混ぜます。

$$
ρ_{\mathrm{mix}}=\frac12\Psi_\uparrow\Psi_\uparrow^\dagger+\frac12\Psi_\downarrow\Psi_\downarrow^\dagger
=\frac12\begin{pmatrix}1&0\\0&1\end{pmatrix}
,\quad
\boldsymbol r=(0,0,0)
$$

ブロッホベクトルは原点（球の中心）を指します。これを**最大混合状態**と呼びます。この作り方では、$\Psi_\uparrow,\Psi_\downarrow$という2つの純粋状態のうちどちらか一方が確率$\frac12$で実現している、という古典的な不確実性を表します。

### 重ね合わせ

一方、状態ベクトル$\Psi_\uparrow$と$\Psi_\downarrow$自体を等係数$c$で足し合わせます。

$$
\Psi_{\mathrm{sup}}=c\Psi_\uparrow+c\Psi_\downarrow = \begin{pmatrix}c\\c\end{pmatrix}\quad(c>0)
$$

規格化条件$|c|^2+|c|^2=1$から$c=\frac1{\sqrt2}$です。

$$
\Psi_{\mathrm{sup}}=\frac1{\sqrt2}\begin{pmatrix}1\\1\end{pmatrix}
$$

これは1本の状態ベクトルなので純粋状態です。

$$
ρ_{\mathrm{sup}}=\Psi_{\mathrm{sup}}\Psi_{\mathrm{sup}}^\dagger=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}=\frac12(\mathrm I+σ_x)
,\quad
\boldsymbol r=(1,0,0)
$$

ブロッホベクトルは$x$軸上の球面を指します。混合状態の$\boldsymbol r=(0,0,0)$とは異なる、確定した1つの向きです。

### 比較

対角成分（$z$成分の情報）は$ρ_{\mathrm{mix}}$と$ρ_{\mathrm{sup}}$で同じですが、非対角成分（$x,y$成分）が異なります。この例では、混合によって非対角成分が消えて原点（球の内部）を指すのに対し、重ね合わせでは非対角成分が残り球面上の点を指します。「半々」という確率的な言葉は共通していても、混合は密度行列の線形結合、重ね合わせは状態ベクトルの線形結合であり、両者は演算として別物です。

&&&rem 測定における差異
混合と重ね合わせの物理的な差異は、測定において現れます。両者は$z$軸方向の測定では区別できませんが、$x$軸方向の測定では異なる確率を与えます。
&&&

## 純粋度

冪等性の計算で使った式は、混合状態でもそのまま成り立ちます。

$$
ρ^2=\frac14\{(1+|\boldsymbol r|^2)\,\mathrm I+2\,\boldsymbol r\cdot\boldsymbol σ\}
$$

トレースを取れば、パウリ行列の部分が消えます。

&&&def 純粋度
密度行列$ρ$に対して$\mathrm{tr}\,ρ^2$を**純粋度**と呼ぶ。2準位系では次のようになる。
$$
\mathrm{tr}\,ρ^2=\frac12(1+|\boldsymbol r|^2)
$$
&&&

球の表面$|\boldsymbol r|=1$で最大値$1$（純粋状態）、中心で最小値$\frac12$（最大混合状態）となります。ブロッホベクトルの長さは、状態がどれだけ純粋かを測る指標です。

## 固有分解

混合状態$ρ$は、互いに直交する2つの純粋状態$ρ_+,ρ_-$を重み$λ_\pm$で混合したものとして書けます。

$$
ρ=λ_+ρ_++λ_-ρ_-
$$

以下では、$ρ_\pm$と$λ_\pm$を具体的に構成します。

&&&rem 「固有分解」という名前の由来
$ρ_\pm$は純粋状態なので、規格化された状態ベクトル$\Psi_\pm$を使って$ρ_\pm=\Psi_\pm\Psi_\pm^\dagger$と書けます。このとき上の式は次のようになります。

$$
ρ=λ_+\Psi_+\Psi_+^\dagger+λ_-\Psi_-\Psi_-^\dagger
$$

これは$ρ$を固有値$λ_\pm$と固有ベクトル$\Psi_\pm$で分解した形となるため、**固有分解**と呼ばれます。

$\boldsymbol r\ne\boldsymbol 0$のとき、$ρ_\pm=\Psi_\pm\Psi_\pm^\dagger$から計算したブロッホベクトルは、$ρ$のブロッホベクトル$\boldsymbol r$と同じ向きか逆向きかのいずれかになります。そのため、行列$ρ$からブロッホベクトル$\boldsymbol r$が既に求まっていれば、固有ベクトル$\Psi_\pm$を改めて計算する必要はありません。
&&&

密度行列$ρ$のブロッホベクトルを$\boldsymbol r$とします。

$$
ρ=\frac12(\mathrm I+\boldsymbol r\cdot\boldsymbol σ)
$$

正規化したブロッホベクトル$\boldsymbol n$を使って、対蹠点方向の密度行列$ρ_\pm$を定義します。

$$
ρ_\pm=\frac12(\mathrm I\pm\boldsymbol n\cdot\boldsymbol σ) \quad \left(\boldsymbol n=\frac{\boldsymbol r}{|\boldsymbol r|},\ |\boldsymbol r|\neq 0\right)
$$

$(\boldsymbol n\cdot\boldsymbol σ)^2=\mathrm I$より、$ρ_\pm$は冪等（$ρ_\pm^2=ρ_\pm$）かつトレース$1$なので、ブロッホベクトル$\pm\boldsymbol n$を持つ純粋状態の密度行列です。

重み$λ_\pm$で$ρ_+,ρ_-$を混ぜたときのブロッホベクトルは、$ρ_\pm$のブロッホベクトルが$\pm\boldsymbol n$であることから、次の条件を満たす必要があります。

$$
λ_+\boldsymbol n-λ_-\boldsymbol n=\boldsymbol r
\quad\Leftrightarrow\quad
λ_+-λ_-=|\boldsymbol r|
$$

これと規格化条件$λ_++λ_-=1$より、重みが求まります。

$$
λ_\pm=\frac12(1\pm|\boldsymbol r|)
$$

トレース$1$のエルミート行列が密度行列となる条件は正半定値であること（固有値が$0$以上であること）で、ここでは$λ_\pm\ge0$、すなわち$|\boldsymbol r|\le1$と同値です。ブロッホ球の球体全体が、ちょうど密度行列全体に対応します。

&&&rem 固有値と退化
$λ_\pm$は$ρ$の固有値と一致します。

$\boldsymbol r\neq\boldsymbol 0$のとき$λ_+\neq λ_-$なので、固有値は非退化であり、対応する固有空間はそれぞれ1次元で、$ρ_\pm$がその固有ベクトル方向を表します。

$\boldsymbol r=\boldsymbol 0$（最大混合状態）の場合は$\boldsymbol n=\boldsymbol r/|\boldsymbol r|$が定義できず、$ρ$から一意に定まる軸はありません。このとき$ρ=\frac12\mathrm I$となって固有値は$λ_+=λ_-=\frac12$に縮退し、任意の単位ベクトル$\boldsymbol n$について$\frac12ρ_++\frac12ρ_-=\frac12\mathrm I$と分解できます。
&&&

実際に混合を計算してみると、元の密度行列$ρ$が得られます。

$$
λ_+ρ_++λ_-ρ_-
=\frac12\left\{(λ_++λ_-)\,\mathrm I+(λ_+-λ_-)\,\boldsymbol n\cdot\boldsymbol σ\right\}
=\frac12(\mathrm I+\boldsymbol r\cdot\boldsymbol σ)=ρ
$$

さらに$(\boldsymbol n\cdot\boldsymbol σ)^2=\mathrm I$より$ρ_+ρ_-=\frac14\{\mathrm I-(\boldsymbol n\cdot\boldsymbol σ)^2\}=0$なので、$ρρ_\pm=λ_\pmρ_\pm$が成り立ちます。$ρ_\pm$は互いに直交する固有ベクトルへの射影であり、単なる2つの純粋状態の混合とは異なります。

密度行列は「どう混ぜたか」の情報を持たないため、同じ密度行列を与える混合のさせ方は1通りではありません。

&&&ex 異なる混合から生じる同じ密度行列
$\boldsymbol r_1=(1,0,0)$（$x$軸）と$\boldsymbol r_2=(0,0,1)$（北極）を確率$\frac12$ずつで混合する。
$$
\frac12ρ_1+\frac12ρ_2=\frac14(\mathrm I+σ_x)+\frac14(\mathrm I+σ_z)=\frac12\mathrm I+\frac14(σ_x+σ_z)
$$
ブロッホベクトルは$\boldsymbol r=\bigl(\frac12,0,\frac12\bigr)$である。この$ρ$を固有分解すると、軸は$\boldsymbol n=\boldsymbol r/|\boldsymbol r|=\frac1{\sqrt2}(1,0,1)$となり、$\boldsymbol r_1,\boldsymbol r_2$のどちらとも異なる。$ρ$だけを見ても、それが$x$軸と$z$軸の混合から作られたことは分からず、混合に使った元の方向の情報は失われている。
&&&

&&&rem 重み付き平均による直接計算
この例のブロッホベクトル$\boldsymbol r=\bigl(\frac12,0,\frac12\bigr)$は、$\boldsymbol r=\frac12\boldsymbol r_1+\frac12\boldsymbol r_2$から直接求まります。行列を展開して計算するまでもなく、混合状態のブロッホベクトルが重み付き平均であるという一般則そのままです。
&&&

# 測定との関係

密度行列が表す系に対して、ある方向を選んで**測定**を行えば、$\boldsymbol n$方向・$-\boldsymbol n$方向に対応する2値の結果（$+1$または$-1$）が確率的に得られます。この節では、多数回測定したときの平均である**期待値**から始め、そこから1回ごとの結果の確率を取り出します。

## 期待値とトレース

単位ベクトル$\boldsymbol n=(n_x,n_y,n_z)$が指す方向への測定を考えます。ここでの$\boldsymbol n$は任意に選ぶ測定方向で、固有分解で用いた$\boldsymbol r/|\boldsymbol r|$とは限りません。

方向$\boldsymbol n$を行列$A=\boldsymbol n\cdot\boldsymbol σ$で表します。本記事では、観測量$A$の測定結果は固有値$\pm1$のいずれかで、その期待値は$\mathrm{tr}(ρA)$で与えられるという量子測定の規則を採用します。$\mathrm{tr}\,σ_i=0$と$\mathrm{tr}(σ_iσ_j)=2\delta_{ij}$を使えば、このトレースはブロッホベクトル$\boldsymbol r$の$\boldsymbol n$方向への符号付きの射影成分、つまり内積$\boldsymbol r\cdot\boldsymbol n$になります。

&&&fml 期待値のトレース表示
$$
\mathrm{tr}(ρA)
=\mathrm{tr}\left(\frac12(\mathrm I+\boldsymbol r\cdot\boldsymbol σ)(\boldsymbol n\cdot\boldsymbol σ)\right)
=\boldsymbol r\cdot\boldsymbol n
$$
&&&

特に、$\boldsymbol n$として$x$軸、$y$軸、$z$軸それぞれの単位ベクトルを取れば、$A$はパウリ行列$σ_x,σ_y,σ_z$そのものになり、期待値はブロッホベクトルの各成分と一致します。

$$
\mathrm{tr}(ρσ_x)=x
,\quad
\mathrm{tr}(ρσ_y)=y
,\quad
\mathrm{tr}(ρσ_z)=z
$$

## 期待値から確率へ

期待値は、1回の測定結果そのものではありません。方向$\boldsymbol n$への測定を行えば、1回ごとの結果は$+1$または$-1$のどちらかの値として確率的に観測されます。同じ状態を多数回測定したときの平均が期待値$\boldsymbol r\cdot\boldsymbol n$です。

&&&rem 状態の収縮
ここで扱う理想的な射影測定を考えます。例えば$\boldsymbol n$が$z$軸方向なら、結果$+1,-1$に対応する状態は、ブロッホベクトルが北極・南極を指す純粋状態の密度行列、つまり$ρ_\uparrow=\frac12(\mathrm I+σ_z),\ ρ_\downarrow=\frac12(\mathrm I-σ_z)$です。物理では、測定で得られた結果に応じて、密度行列は$ρ_\uparrow$または$ρ_\downarrow$へ収縮すると考えます。そのため、測定の間に状態の時間発展などがなければ、同じ測定を連続で行うと、測定直後の状態が対応する固有状態になるため、2回目以降も1回目と同じ結果が得られます。

したがって、上で述べた「同じ状態を多数回測定」とは、この意味での連続測定ではありません。期待値を得るには、同じ状態$ρ$を独立に何度も用意し、そのたびに1回だけ測定するという手続きが必要です。
&&&

$+1$が出る確率を$p_+$、$-1$が出る確率を$p_-$とすれば、期待値は「値×確率」の総和として定義されるので、次が成り立ちます。

$$
\boldsymbol r\cdot\boldsymbol n=(+1)\,p_+ +(-1)\,p_-=p_+-p_-
$$

これに確率の規格化条件$p_++p_-=1$を合わせた連立方程式を解けば、期待値から確率が求まります。

$$
p_+=\frac{1+\boldsymbol r\cdot\boldsymbol n}2
,\quad
p_-=\frac{1-\boldsymbol r\cdot\boldsymbol n}2
$$

この変換は、期待値の範囲$[-1,1]$を確率の範囲$[0,1]$へ写すスケール調整に相当します。$\boldsymbol r=\boldsymbol n$なら確実に$+1$、$\boldsymbol r=-\boldsymbol n$なら確実に$-1$、直交または原点（$\boldsymbol r\cdot\boldsymbol n=0$）なら半々です。

&&&rem 三角関数で見る2乗
純粋状態（$|\boldsymbol r|=1$）に限れば、$\boldsymbol r$と$\boldsymbol n$のなす角$θ$を使って$\boldsymbol r\cdot\boldsymbol n=\cosθ$と書けます。半角の公式$\cos^2\dfracθ2=\dfrac{1+\cosθ}2$より、次のようになります。
$$
p_+=\frac{1+\cosθ}2=\cos^2\fracθ2
$$

期待値を確率へ写すこの1次式の変換は、純粋状態では半角公式によって$\cos$の2乗の形に書けます。
&&&

## 混合と重ね合わせの違い

混合と重ね合わせの例に戻り、測定の期待値$\boldsymbol r\cdot\boldsymbol n$を使って両者を比較します。混合状態$ρ_{\mathrm{mix}}=\frac12\mathrm I$のブロッホベクトルは$\boldsymbol r=(0,0,0)$、重ね合わせ$ρ_{\mathrm{sup}}$のブロッホベクトルは$\boldsymbol r=(1,0,0)$でした。

$z$軸方向の期待値はどちらも$0$で一致します。

$$
\mathrm{tr}(ρ_{\mathrm{mix}}σ_z)=0
,\quad
\mathrm{tr}(ρ_{\mathrm{sup}}σ_z)=0
$$

ところが$x$軸方向の期待値は一致しません。

$$
\mathrm{tr}(ρ_{\mathrm{mix}}σ_x)=0
,\quad
\mathrm{tr}(ρ_{\mathrm{sup}}σ_x)=1
$$

確率で言い換えると、$z$軸方向の測定結果はどちらの状態でも$p_+=p_-=1/2$（半々）です。$x$軸方向では、重ね合わせは$p_+=1,\ p_-=0$、つまり確実に$+1$が観測されるのに対し、混合状態は$x$軸でも$p_+=p_-=1/2$のままです。

最大混合状態は、$\boldsymbol r=(0,0,0)$より$\boldsymbol r\cdot\boldsymbol n=0$がどの$\boldsymbol n$についても成り立ちます。つまり$z$軸で半々だったのは$x$軸でもたまたま同じ確率になったわけではなく、どの方向に測定しても半々になります。「北極と南極の混合」という材料からは$z$軸方向の不確実性しか見えませんが、$ρ_{\mathrm{mix}}$自体はあらゆる方向に対して不定です。

一方、重ね合わせ$ρ_{\mathrm{sup}}$は$x$軸方向の測定で確実に$+1$を与える純粋状態です。$z$軸方向（$x$軸に直交する方向）の測定で結果が半々になることは、混合状態であることを意味しません。

# まとめ

単位行列とパウリ行列は2×2エルミート行列の基底で、パウリ行列の係数が3次元の座標を与えます。純粋状態の密度行列はブロッホ球の表面を指し、位相（ファイバー）が潰れることから、ホップファイブレーションの像の行列表現とみなせます。混合状態はブロッホベクトルの重み付き平均としてブロッホ球の内部に広がります。

ブロッホベクトルの各成分は各軸方向の測定の期待値を与え、その向きと長さに相対位相や純粋度の情報が反映されています。この点に、この表現が量子情報で使われる理由があります。
