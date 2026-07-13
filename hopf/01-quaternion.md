ホップファイブレーションは射影の一種で、三次元球面（四次元空間内の超球面）上の点を二次元球面（三次元空間内の球面）上に射影します。四元数による計算方法を確認して、パウリ行列を求めます。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

改訂履歴：

- 2026.07.13 ファイバーの局所自明化を具体形で追加
- 2025.01.24 符号調整を行わずに内積・外積と対応付けるよう修正

# 概要

単位四元数$q$について、四元数共役を$q^*$とすれば、純虚四元数$p$に対して$qpq^*$は回転を表します。回転軸を$n$、回転角を$θ$とすれば、パラメーターは以下の通りです。[[7shi]]

&&&def 単位四元数による回転子
$$
\begin{aligned}
n&=n_x\,\mathbf{i}+n_y\,\mathbf{j}+n_z\,\mathbf{k}\quad(n_x^2+n_y^2+n_z^2=1,\ n^2=-1) \\
q&=\cos\frac\theta2+n\sin\frac\theta2=e^{nθ/2}
\end{aligned}
$$
&&&

Wikipedia からホップファイブレーションと四元数の関係を引用します。[[wikipedia]]

> In fact, this identifies the group of versors with the group of rotations of $\mathbb{R}^3$, modulo the fact that the versors $q$ and $-q$ determine the same rotation. As noted above, the rotations act transitively on $\mathbb{S}^2$, and the set of versors $q$ which fix a given right versor $p$ have the form $q = u + vp$, where $u$ and $v$ are real numbers with $u^2 + v^2 = 1$. This is a circle subgroup. For concreteness, one can take $p = \mathbf k$, and then the Hopf fibration can be defined as the map sending a versor $ω$ to $ω\mathbf kω^*$. All the quaternions $ωq$, where $q$ is one of the circle of versors that fix $\mathbf k$, get mapped to the same thing (which happens to be one of the two $180°$ rotations rotating $\mathbf k$ to the same place as $ω$ does).

日本語訳

> 実際、単位四元数$q$と$-q$が同じ回転を決定するという事実を除いて、これは単位四元数の群を$\mathbb{R}^3$の回転群と同一視している。上記のように、回転は$\mathbb{S}^2$に推移的に作用し、与えられた純虚単位四元数$p$を固定する単位四元数$q$の集合は$q = u + vp$の形をしている。ここで$u$と$v$は$u^2 + v^2 = 1$を満たす実数である。これは円周部分群である。具体的には$p = \mathbf k$とすることができ、そうするとホップファイブレーションは、単位四元数$ω$を$ω\mathbf kω^*$に送る写像として定義できる。$q$が$\mathbf k$を固定する単位四元数の円周の1つであれば、すべての四元数$ωq$は同じもの（$\mathbf k$を$ω$と同じ場所に回転させる 2 つの$180°$回転のうちの 1 つ）に写される。

2 点に絞って説明します。

1. 回転: 始点を$\mathbf k$に固定すれば、単位四元数$ω$での回転により、ホップファイブレーションが定義できます。
2. ファイバー: ホップファイブレーションは単射ではなく、ある一点の逆像は円周となります。

# 回転

始点を$\mathbf k$に固定すれば、単位四元数$ω$での回転により、ホップファイブレーションが定義できます。

> ホップファイブレーションは、単位四元数$ω$を$ω\mathbf kω^*$に送る写像として定義できる。[[wikipedia]]

&&&def 回転子 ω
$$
ω=ω_0+ω_1\mathbf{i}+ω_2\mathbf{j}+ω_3\mathbf{k}\quad(ω_0^2+ω_1^2+ω_2^2+ω_3^2=1)
$$
&&&

$ω$によって虚数単位$\mathbf{k}$を回転させます。

$$
\begin{aligned}
ω\mathbf{k}ω^*
&=(ω_0+ω_1\mathbf{i}+ω_2\mathbf{j}+ω_3\mathbf{k})\mathbf{k}(ω_0+ω_1\mathbf{i}+ω_2\mathbf{j}+ω_3\mathbf{k})^* \\
&=(ω_0\mathbf{k}-ω_1\mathbf{j}+ω_2\mathbf{i}-ω_3)(ω_0-ω_1\mathbf{i}-ω_2\mathbf{j}-ω_3\mathbf{k}) \\
&=  ω_0\mathbf{k}(ω_0-ω_1\mathbf{i}-ω_2\mathbf{j}-ω_3\mathbf{k}) \\
&\ -ω_1\mathbf{j}(ω_0-ω_1\mathbf{i}-ω_2\mathbf{j}-ω_3\mathbf{k}) \\
&\ +ω_2\mathbf{i}(ω_0-ω_1\mathbf{i}-ω_2\mathbf{j}-ω_3\mathbf{k}) \\
&\ -ω_3          (ω_0-ω_1\mathbf{i}-ω_2\mathbf{j}-ω_3\mathbf{k}) \\
&=  ω_0^2\mathbf{k} -ω_0ω_1\mathbf{j} +ω_0ω_2\mathbf{i} +ω_0ω_3 \\
&\ -ω_0ω_1\mathbf{j} -ω_1^2\mathbf{k} -ω_1ω_2 +ω_1ω_3\mathbf{i} \\
&\ +ω_0ω_2\mathbf{i} +ω_1ω_2 -ω_2^2\mathbf{k} +ω_2ω_3\mathbf{j} \\
&\ -ω_0ω_3 +ω_1ω_3\mathbf{i} +ω_2ω_3\mathbf{j} +ω_3^2\mathbf{k} \\
&=  ( ω_0ω_3 -ω_1ω_2 +ω_1ω_2 -ω_0ω_3) \\
&\ +( ω_0ω_2 +ω_1ω_3 +ω_0ω_2 +ω_1ω_3)\mathbf{i} \\
&\ +(-ω_0ω_1 -ω_0ω_1 +ω_2ω_3 +ω_2ω_3)\mathbf{j} \\
&\ +( ω_0^2  -ω_1^2  -ω_2^2  +ω_3^2 )\mathbf{k} \\
&= 2(ω_0ω_2+ω_1ω_3)\mathbf{i} +2(ω_2ω_3-ω_0ω_1)\mathbf{j} +\{(ω_0^2+ω_3^2)-(ω_1^2+ω_2^2)\}\mathbf{k}
\end{aligned}
$$

$ω\mathbf{k}ω^*$は実部を持たない純虚四元数です。係数の 2 乗和を確認します。

$$
\begin{aligned}
&\{2(ω_0ω_2+ω_1ω_3)\}^2+\{2(ω_2ω_3-ω_0ω_1)\}^2+\{(ω_0^2+ω_3^2)-(ω_1^2+ω_2^2)\}^2 \\
=&4ω_0^2ω_2^2+8ω_0ω_1ω_2ω_3+4ω_1^2ω_3^2+4ω_2^2ω_3^2-8ω_0ω_1ω_2ω_3+4ω_0^2ω_1^2 \\
&+ω_0^4+2ω_0^2ω_3^2+ω_3^4-2(ω_0^2ω_1^2+ω_0^2ω_2^2+ω_3^2ω_1^2+ω_3^2ω_2^2)+ω_1^4+2ω_1^2ω_2^2+ω_2^4 \\
=&ω_0^4+ω_1^4+ω_2^4+ω_3^4+2ω_0^2ω_1^2+2ω_0^2ω_2^2+2ω_0^2ω_3^2+2ω_1^2ω_2^2+2ω_1^2ω_3^2+2ω_2^2ω_3^2 \\
=&(ω_0^2+ω_1^2+ω_2^2+ω_3^2)^2 \\
=&1
\end{aligned}
$$

よって$ω\mathbf{k}ω^*$は純虚単位四元数であり、単位球面上の一点を指すことが確認できました。

&&&rem ノルムの乗法性
四元数のノルムは乗法的$|pq|=|p||q|$であるため、展開しなくても$|ω\mathbf{k}ω^*|=|ω||\mathbf{k}||ω^*|=1$が直ちに得られます。
&&&

## 内積・外積

ここで、2次元ベクトル間の内積と外積（ウェッジ積）を思い出してみます。

内積：$\begin{pmatrix}α_0\\α_1\end{pmatrix}\cdot\begin{pmatrix}β_0\\β_1\end{pmatrix}=α_0β_0+α_1β_1$

外積：$\begin{pmatrix}α_0\\α_1\end{pmatrix}\wedge\begin{pmatrix}β_0\\β_1\end{pmatrix}=\star(α_0β_1-α_1β_0)$

$ω\mathbf{k}ω^*$の$\mathbf{i},\mathbf{j}$の係数が内積・外積と対応付けられるように並べ替えます。

$$
\begin{aligned}
ω\mathbf{k}ω^*
&= 2(ω_0ω_2+ω_1ω_3)\mathbf{i} +2(ω_2ω_3-ω_0ω_1)\mathbf{j} +\{(ω_0^2+ω_3^2)-(ω_1^2+ω_2^2)\}\mathbf{k} \\
&= 2(ω_3ω_1+ω_0ω_2)\mathbf{i} +2(ω_3ω_2-ω_0ω_1)\mathbf{j} +\{(ω_3^2+ω_0^2)-(ω_1^2+ω_2^2)\}\mathbf{k} \\
\end{aligned}
$$

ここで以下のように対応付ければ、

$$
\begin{pmatrix}α_0\\α_1\end{pmatrix}=\begin{pmatrix}ω_3\\ω_0\end{pmatrix}
,\quad
\begin{pmatrix}β_0\\β_1\end{pmatrix}=\begin{pmatrix}ω_1\\ω_2\end{pmatrix}
$$

$ω\mathbf{k}ω^*$の係数に内積と外積の2倍が含まれることが分かります。

$$
ω\mathbf{k}ω^*
= 2(α_0β_0+α_1β_1)\mathbf{i} +2(α_0β_1-α_1β_0)\mathbf{j} +\{(α_0^2+α_1^2)-(β_0^2+β_1^2)\}\mathbf{k}
$$

## 複素数への分割

内積・外積を求めるのに使用したベクトルを複素数で表現することで、$ω$を2つの複素数へ分割します。これらを$α,β$とします。

&&&def 単位四元数の複素数への分割
$$
\begin{alignedat}{3}
α&=&α_0&+α_1i&&=ω_3+ω_0i \\
β&=&β_0&+β_1i&&=ω_1+ω_2i
\end{alignedat}
$$
&&&

左因子を複素共役にして積を求めれば、内積と外積が得られます。

$$
\begin{aligned}
α^*β
&=(α_0-α_1i)(β_0+β_1i) \\
&=(α_0β_0+α_1β_1)+(α_0β_1-α_1β_0)i
\end{aligned}
$$

$ω\mathbf{k}ω^*$を書き換えます。

$$
ω\mathbf{k}ω^*=2\mathrm{Re}(α^*β)\mathbf{i}+2\mathrm{Im}(α^*β)\mathbf{j}+(|α|^2-|β|^2)\mathbf{k}
$$

よく使われる形のホップファイブレーションの表式が四元数によって得られました。

## まとめ

ここまでの結果をまとめます。

&&&def 複素数ペアと四元数の対応関係
$$
\begin{aligned}
α&=α_0+α_1i \\
β&=β_0+β_1i \\
ω&=ω_0+ω_1\mathbf{i}+ω_2\mathbf{j}+ω_3\mathbf{k}
\end{aligned}
$$
$$
\begin{pmatrix}α_0\\α_1\end{pmatrix}=\begin{pmatrix}ω_3\\ω_0\end{pmatrix}
,\quad
\begin{pmatrix}β_0\\β_1\end{pmatrix}=\begin{pmatrix}ω_1\\ω_2\end{pmatrix}
$$
$$
\begin{aligned}
ω&=α_1+β_0\mathbf{i}+β_1\mathbf{j}+α_0\mathbf{k} \\
 &=α_0\mathbf{k}+α_1+β_0\mathbf{i}+β_1\mathbf{j} \\
 &=(α_0\mathbf{k}+α_1)+\mathbf{j}(β_0\mathbf{k}+β_1)
\end{aligned}
$$
&&&
&&&def 複素数ペアの四元数への埋め込み
$$
\begin{aligned}
α_ω&=α_0\mathbf{k}+α_1 \\
β_ω&=β_0\mathbf{k}+β_1 \\
ω&=α_ω+\mathbf{j}β_ω
\end{aligned}
$$
&&&

$ω$で$\mathbf{k}$を回転させます。

$$
\begin{aligned}
ω\mathbf{k}ω^*
&=(α_ω+\mathbf{j}β_ω)\mathbf{k}(α_ω+\mathbf{j}β_ω)^* \\
&=(α_ω+β_ω^*\mathbf{j})\mathbf{k}(α_ω^*-\mathbf{j}β_ω) \\
&=(α_ω\mathbf{k}+β_ω^*\mathbf{i})(α_ω^*-\mathbf{j}β_ω) \\
&=α_ω\mathbf{k}α_ω^*-α_ω\mathbf{k}\mathbf{j}β_ω+β_ω^*\mathbf{i}α_ω^*-β_ω^*\mathbf{i}\mathbf{j}β_ω \\
&=α_ω\mathbf{k}α_ω^*+α_ω\mathbf{i}β_ω+β_ω^*\mathbf{i}α_ω^*-β_ω^*\mathbf{k}β_ω \\
&=α_ωα_ω^*\mathbf{k}+α_ωβ_ω^*\mathbf{i}+β_ω^*α_ω\mathbf{i}-β_ω^*β_ω\mathbf{k} \\
&=2α_ωβ_ω^*\mathbf{i}+(|α|^2-|β|^2)\mathbf{k} \\
&=2(α_0\mathbf{k}+α_1)(-β_0\mathbf{k}+β_1)\mathbf{i}+(|α|^2-|β|^2)\mathbf{k} \\
&=2(α_0β_0+α_0β_1\mathbf{k}-α_1β_0\mathbf{k}+α_1β_1)\mathbf{i}+(|α|^2-|β|^2)\mathbf{k} \\
&=2(α_0β_0+α_1β_1)\mathbf{i}+2(α_0β_1-α_1β_0)\mathbf{j}+(|α|^2-|β|^2)\mathbf{k} \\
&=2\mathrm{Re}(α^*β)\mathbf{i}+2\mathrm{Im}(α^*β)\mathbf{j}+(|α|^2-|β|^2)\mathbf{k}
\end{aligned}
$$

# 実部と虚部の分離

複素共役での虚部の符号反転によって、実部と虚部が分離できます。分離した虚部に$-i$を掛けることで実数化します。

&&&fml 実部と虚部の分離
$$
\begin{alignedat}{4}
2\mathrm{Re}(α^*β)&=&&α^*β+(α^*β)^*&&=&&α^*β+β^*α \\
2\mathrm{Im}(α^*β)&=-i(&&α^*β-(α^*β)^*)&&=-i(&&α^*β-β^*α)
\end{alignedat}
$$
&&&

これを $ω\mathbf{k}ω^*$ に適用します。

$$
\begin{aligned}
ω\mathbf{k}ω^*
&=2\mathrm{Re}(α^*β)\mathbf{i}+2\mathrm{Im}(α^*β)\mathbf{j}+(|α|^2-|β|^2)\mathbf{k} \\
&=(α^*β+β^*α)\mathbf{i}-i(α^*β-β^*α)\mathbf{j}+(α^*α-β^*β)\mathbf{k}
\end{aligned}
$$

&&&rem 虚数単位 i の区別
複素数の$i$と四元数の$\mathbf i$を表記上区別しています。$-i(α^*β-β^*α)$は実数となるため、得られる結果に複素数と四元数が混在することはありません。
&&&

# パウリ行列

複素数$α,β$を成分とする複素ベクトル$\Psi$とそのエルミート共役（成分の複素共役を取る転置）を定義します。

&&&def 複素ベクトル Ψ とそのエルミート共役
$$
\Psi=\begin{pmatrix}α\\β\end{pmatrix},\ \Psi^\dagger=\begin{pmatrix}α^*&β^*\end{pmatrix}
$$
&&&

$ω\mathbf{k}ω^*$の$i,j,k$の係数を$x,y,z$とします。

$$
ω\mathbf{k}ω^*
=(\underbrace{α^*β+β^*α}_{x})\mathbf{i}\underbrace{-i(α^*β-β^*α)}_{y}\mathbf{j}+(\underbrace{α^*α-β^*β}_{z})\mathbf{k}
$$

各項に含まれる複素共役を$\Psi^\dagger$に対応させ、ベクトルの内積で表現します。内積の右因子は、$\Psi$から変換するための行列$\sigma_x,\sigma_y,\sigma_z$を作って挟みます。

$$
\begin{alignedat}{5}
x
 &=&&α^*β+β^*α &
 &=\begin{pmatrix}α^*&β^*\end{pmatrix}\begin{pmatrix}β\\α\end{pmatrix} &
 &=\begin{pmatrix}α^*&β^*\end{pmatrix}\begin{pmatrix}0&1\\1&0\end{pmatrix}\begin{pmatrix}α\\β\end{pmatrix} &
 &=\Psi^\dagger\sigma_x\Psi \\
y&=-i(&&α^*β-β^*α) &
 &=\begin{pmatrix}α^*&β^*\end{pmatrix}\begin{pmatrix}\begin{array}{r}-iβ\\iα\end{array}\end{pmatrix} &
 &=\begin{pmatrix}α^*&β^*\end{pmatrix}\begin{pmatrix}0&-i\\i&0\end{pmatrix}\begin{pmatrix}α\\β\end{pmatrix} &
 &=\Psi^\dagger\sigma_y\Psi \\
z&=&&α^*α-β^*β &
 &=\begin{pmatrix}α^*&β^*\end{pmatrix}\begin{pmatrix}\begin{array}{r}α\\-β\end{array}\end{pmatrix} &
 &=\begin{pmatrix}α^*&β^*\end{pmatrix}\begin{pmatrix}1&0\\0&-1\end{pmatrix}\begin{pmatrix}α\\β\end{pmatrix} &
 &=\Psi^\dagger\sigma_z\Psi \\
\end{alignedat}
$$

これらの変換行列が**パウリ行列**です。

&&&def パウリ行列
$$
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix}
,\ \sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix}
,\ \sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}
$$
&&&

パウリ行列による計算は、$ω\mathbf{k}ω^*$と同等の計算を形式的かつ簡潔にまとめます。

&&&fml パウリ行列による座標変換
$$\Psi\mapsto\begin{pmatrix}
  \Psi^\dagger\sigma_x\Psi \\
  \Psi^\dagger\sigma_y\Psi \\
  \Psi^\dagger\sigma_z\Psi\end{pmatrix}$$
&&&

量子情報では、$\Psi$を**状態ベクトル**、パウリ行列による座標変換によって構成される単位球を**ブロッホ球**と呼びます。

# ファイバー

ホップファイブレーションは全射ですが単射ではなく、ある一点の逆像は円周となります。この円周は**ファイバー**と呼ばれます。

> 与えられた純虚単位四元数$p$を固定する単位四元数$q$の集合は$q = u + vp$の形をしている。ここで$u$と$v$は$u^2 + v^2 = 1$を満たす実数である。これは円周部分群である。[[wikipedia]]

純虚単位四元数の 2 乗は -1 となります。

$$
p^2=-1
$$

$q$によって$p$を回転させます。

$$
\begin{aligned}
qpq^*
&=(u+vp)p(u+vp)^* \\
&=(u+vp)p(u-vp) \\
&=(u+vp)(u-vp)p \\
&=(u^2+v^2)p \\
&=p
\end{aligned}
$$

$p$ が固定されていることが確認できました。

> 具体的には$p = \mathbf k$とすることができ、そうするとホップファイブレーションは、単位四元数$ω$を$ω\mathbf kω^*$に送る写像として定義できる。$q$が$\mathbf k$を固定する単位四元数の円周の1つであれば、すべての四元数$ωq$は同じもの（$\mathbf k$を$ω$と同じ場所に回転させる 2 つの$180°$回転のうちの 1 つ）に写される。

$p=\mathbf{k}$とおけば、$q\mathbf{k}q^*=\mathbf{k}$より
$$
(ωq)\mathbf{k}(ωq)^*=ω(q\mathbf{k}q^*)ω^*=ω\mathbf{k}ω^*
$$

となるため、$ωq$は$ω$と同じ点$ω\mathbf{k}ω^*$に写されます。$ωq$の軌跡がこの点$ω\mathbf{k}ω^*$の逆像であり、円周となります。ホップファイブレーションによって、超球面$S^3$は球面$S^2$上に円周$S^1$を貼り付けたファイバー束として捉えられます。

$$
S^{1}\hookrightarrow S^{3}\xrightarrow {\ \text{Hopf} \,} S^{2}
$$

この束は、球面$S^2$の上に円周をただ並べただけの直積$S^2\times S^1$とは異なります。

## 積構造とのずれ

$S^3$上の点$a$を「球面上の位置$b$」と「ファイバー上の位置$c$」の組$(b,c)$に分解しようとすると、そのずれが現れます。$S^2$を北半球・南半球の2枚のパッチに分ければ、それぞれのパッチの上では$a\mapsto(b,c)$を連続に定義できます。しかし赤道での貼り合わせを見ると、北半球側の$c$座標と南半球側の$c$座標は単純には一致せず、赤道を1周する間に1回転ねじれています。このねじれのせいで、$(b,c)$の対応を$S^2$全体にわたって連続に選び直すことができず、$S^3$を$S^2\times S^1$という単純な直積とみなすことはできません。

実際に$a\mapsto(b,c)$を書き下してみます。$a=ω$に対応する複素数ペアを$(α,β)$とすれば、$b=(x,y,z)$は既出の対応関係から次のように書けます。

$$
x=α^*β+β^*α,\quad y=-i(α^*β-β^*α),\quad z=α^*α-β^*β
$$

ファイバー方向は、$q=u+v\mathbf{k}\ (u^2+v^2=1)$による右からの積$ω\mapsto ωq$で生成されます。$ω=α_ω+\mathbf{j}β_ω$（$α_ω=α_0\mathbf{k}+α_1,\ β_ω=β_0\mathbf{k}+β_1$）に対して、

$$
ωq=(α_ω+\mathbf{j}β_ω)(u+v\mathbf{k})=α_ω(u+v\mathbf{k})+\mathbf{j}β_ω(u+v\mathbf{k})
$$

$u+v\mathbf{k}$は、複素数$ζ=u+vi$を$\{1,\mathbf{k}\}$に埋め込んだもの（$1\mapsto1,\ i\mapsto\mathbf{k}$という環同型による像）に一致します。$α_ω,β_ω$も同じ埋め込みによる$α,β$の像なので、埋め込みが積を保つことから$α_ω(u+v\mathbf{k})$は$αζ$の像、$β_ω(u+v\mathbf{k})$は$βζ$の像になります。すなわち$ωq$に対応する複素数ペアは$(αζ,βζ)$です。

&&&def ファイバー方向の作用
$$
ω\mapsto ωq \iff (α,β)\mapsto(ζα,ζβ)\quad(ζ=u+vi,\ |ζ|=1)
$$
&&&

$b$は$(α,β)$の比だけで決まるため、この位相$ζ$の分だけ動いても$b$は変化しません。$ζ$こそがファイバー上の位置$c$です。

$α\neq0$となるパッチ（$z\neq-1$、南極を除く全域）では、$c_N=α/|α|$とおけます。$|α|^2=(1+z)/2$と$α^*β=(x+iy)/2$を使えば、

$$
α=\sqrt{\frac{1+z}2}\,c_N,\qquad β=\frac{x+iy}{\sqrt{2(1+z)}}\,c_N
$$

これが北側パッチ上の$a\mapsto(b,c_N)$の具体形です。同様に$β\neq0$となるパッチ（$z\neq1$、北極を除く全域）では$c_S=β/|β|$とおけて、

$$
β=\sqrt{\frac{1-z}2}\,c_S,\qquad α=\frac{x-iy}{\sqrt{2(1-z)}}\,c_S
$$

両方のパッチが重なる領域（$z\neq\pm1$）で$c_N$と$c_S$を比較します。$x+iy=re^{iφ}$（$r=\sqrt{1-z^2}$、$φ$は$S^2$上の方位角）とおいて$α$の2つの表式を等しく置けば、

$$
\sqrt{\frac{1+z}2}\,c_N=\frac{x-iy}{\sqrt{2(1-z)}}\,c_S=\frac{re^{-iφ}}{\sqrt{2(1-z)}}\,c_S=\sqrt{\frac{1+z}2}\,e^{-iφ}c_S
$$

より、貼り合わせの関係式が得られます。

$$
c_S=e^{iφ}c_N
$$

$φ$は$S^2$の方位角そのものなので、赤道を1周して$φ$が$0$から$2π$まで動くと$c_S/c_N=e^{iφ}$はちょうど1回転します。これが「北半球側の$c$座標と南半球側の$c$座標が赤道を1周する間に1回転ねじれる」ことの具体的な中身であり、このねじれ（巻き数1の遷移関数）のために$(b,c)$を$S^2$全体で連続に選び直すことができません。
