ホップファイブレーションは射影の一種で、三次元球面（四次元空間内の超球面）上の点を二次元球面（三次元空間内の球面）上に射影します。四元数による計算方法を確認して、パウリ行列を求めます。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

改訂履歴：

- 2026.07.15 複素数ホップ写像を定義、まとめを追加
- 2026.07.13 ノルムの乗法性の注記を追加、ファイバーの行間を補完
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

## 結果まとめ

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

計算の過程で、$α_ω,β_ω$を展開する直前で、以下の結果が得られています。

$$
ω\mathbf{k}ω^*=2α_ωβ_ω^*\mathbf{i}+(|α|^2-|β|^2)\mathbf{k}
$$

この結果より、複素数のペア$(α,β)$から$(2αβ^*,|α|^2-|β|^2)$への複素数ホップ写像を定義します。

&&&def 複素数ホップ写像
単位球面上の点 $|\alpha|^2 + |\beta|^2 = 1$ は、以下の写像によって2次元球面上の点に射影されます。
$$
H(\alpha, \beta) = \left(2\alpha \beta^*, |\alpha|^2 - |\beta|^2\right)
$$
&&&

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

> 具体的には$p = \mathbf k$とすることができ、そうするとホップファイブレーションは、単位四元数$ω$を$ω\mathbf kω^*$に送る写像として定義できる。$q$が$\mathbf k$を固定する単位四元数の円周の1つであれば、すべての四元数$ωq$は同じもの（$\mathbf k$を$ω$と同じ場所に回転させる 2 つの$180°$回転のうちの 1 つ）に写される。[[wikipedia]]

$p=\mathbf{k}$とおけば、$q\mathbf{k}q^*=\mathbf{k}$より
$$
(ωq)\mathbf{k}(ωq)^*=ω(q\mathbf{k}q^*)ω^*=ω\mathbf{k}ω^*
$$

となるため、$ωq$は$ω$と同じ点$ω\mathbf{k}ω^*$に写されます。$ωq$の軌跡がこの点$ω\mathbf{k}ω^*$の逆像であり、円周となります。ホップファイブレーションによって、超球面$S^3$は球面$S^2$上に円周$S^1$を貼り付けたファイバー束として捉えられます。

$$
\underbrace{S^{1}}_{q}
\hookrightarrow
\underbrace{S^{3}}_{ω}
\xrightarrow{H}
\underbrace{S^{2}}_{ω\mathbf{k}ω^*}
$$

&&&rem
この束は、「超球面$S^3$上の点」を「球面$S^2$上の点」と「円周$S^1$上の点」の組として一意に分解できるような、単純な直積$S^2\times S^1$ではありません。
&&&

# まとめ

単位四元数$ω$を2つの埋め込まれた複素数$α_ω,β_ω$で表現して、回転$ω\mathbf{k}ω^*$を成分展開します。

$$
ω\mathbf{k}ω^*
=2α_ωβ_ω^*\mathbf{i}+(|α|^2-|β|^2)\mathbf{k}
$$

この結果より、複素数のペアから複素数と実数のペアへの複素数ホップ写像を定義します。

$$
H(\alpha, \beta) = \left(2\alpha \beta^*, |\alpha|^2 - |\beta|^2\right)
$$

$α_ω,β_ω$を成分展開して計算を続ければ、結果が2つの複素数$α,β$の内積・外積として書けます。

$$
ω\mathbf{k}ω^*
=2\,\underbrace{\mathrm{Re}(α^*β)}_{\text{内積}}\,\mathbf{i}+2\,\underbrace{\mathrm{Im}(α^*β)}_{\text{外積}}\,\mathbf{j}+(|α|^2-|β|^2)\mathbf{k}
$$

$\mathrm{Re},\mathrm{Im}$を展開して、虚数単位$\mathbf{i},\mathbf{j},\mathbf{k}$の係数から$x,y,z$成分を求めます。

$$
\begin{aligned}
ω\mathbf{k}ω^*
&=(α^*β+β^*α)\mathbf{i}-i(α^*β-β^*α)\mathbf{j}+(α^*α-β^*β)\mathbf{k} \\
(x,y,z)
&=(α^*β+β^*α,-i(α^*β-β^*α),α^*α-β^*β)
\end{aligned}
$$

複素ベクトル$\Psi=(α,β)^T$とそのエルミート共役$\Psi^\dagger$による内積$\Psi^\dagger\sigma_i\Psi$として整理すると、変換行列としてパウリ行列$\sigma_x,\sigma_y,\sigma_z$が現れます。

$$
\begin{alignedat}{3}
x&=\Psi^\dagger\sigma_x\Psi &
 &=\begin{pmatrix}α^*&β^*\end{pmatrix}\begin{pmatrix}0&1\\1&0\end{pmatrix}\begin{pmatrix}α\\β\end{pmatrix} &
 &=α^*β+β^*α \\
y&=\Psi^\dagger\sigma_y\Psi &
 &=\begin{pmatrix}α^*&β^*\end{pmatrix}\begin{pmatrix}0&-i\\i&0\end{pmatrix}\begin{pmatrix}α\\β\end{pmatrix} &
 &=-i(α^*β-β^*α) \\
z&=\Psi^\dagger\sigma_z\Psi &
 &=\begin{pmatrix}α^*&β^*\end{pmatrix}\begin{pmatrix}1&0\\0&-1\end{pmatrix}\begin{pmatrix}α\\β\end{pmatrix} &
 &=α^*α-β^*β \\
\end{alignedat}
$$

単位純虚四元数$p$は、単位四元数$q=u+vp\ (u^2+v^2=1)$による回転で固定されます。

$$
qpq^* = p
$$

$ωq$による回転は、$ω$による回転と同じ点に移されます。

$$
(ωq)p(ωq)^* = ω(qpq^*)ω^* = ωpω^*
$$

$q$は円周をなす自由度であり、$ωq$は$ω$と同じ点$ωpω^*$に写されます。したがって$ωq$の軌跡がその点の逆像（ファイバー）である円周となり、超球面$S^3$は球面$S^2$上に円周$S^1$を貼り付けたファイバー束として捉えられます。

$$
\underbrace{S^{1}}_{q}
\hookrightarrow
\underbrace{S^{3}}_{ω}
\xrightarrow{H}
\underbrace{S^{2}}_{ω\mathbf{k}ω^*}
$$

&&&rem
ホップファイブレーションの定式化はいくつかあり、射影直線によるものが一般的です。[[dim7-8]][[7shi-hc]]
&&&
