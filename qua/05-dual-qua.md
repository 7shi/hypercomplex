球面三角法の記事では、球面三角形を一周する操作を四元数の積で表した「マスター方程式」から球面三角法の諸公式を導きました。その具体例として、辺の長さを無限小パラメーター$\varepsilon$で縮めた無限小三角形を計算したところ、$\varepsilon$の1次の項に平面三角形の閉合条件が現れ、無限小の領域では回転が並進として扱えることが示唆されました。

本稿では、この「無限小の領域では回転が並進として扱える」という観察を敷衍します。出発点は無限小の読み替えです。回転角を限りなく小さい実数とする代わりに、$\varepsilon^2=0$を満たす形式的な基底$\varepsilon$を使って$\varepsilon a$と表すと、「回転が並進と見なせる」はテイラー近似ではなく厳密な等式になります。実数にこの$\varepsilon$を添加した数が二重数、四元数の係数を二重数に取り替えた数が二重四元数 (dual quaternion) です。二重四元数を使うと並進は接空間の縛りから解放され、3次元の剛体変換（回転と並進の合成）が回転子と同じ「挟み込み」の形で統一的に表現できます。

最後に、二重四元数が退化した計量を持つクリフォード代数$C\ell_{3,0,1}(\mathbb R)$の偶部分代数と同型であることを示します。これは射影幾何代数 (PGA) や共形幾何代数 (CGA) と呼ばれる枠組みへの入り口となります。

四元数の基本的な演算（共役・指数関数・共役作用による回転）は既知とします。内積は$\langle q,r\rangle:=\mathrm{Re}(q\overline r)$で定め、$|q|=\sqrt{\langle q,q\rangle}$と書きます。

# 回転から並進へ

## 微小角の回転

単位球面上の点$i$の近くを考えます。回転子$e^{k\theta/2}$による共役作用は

$$
e^{k\theta/2} i e^{-k\theta/2} = \cos\theta\, i + \sin\theta\, j
$$

と、点$i$を球面に沿って$j$の方向へ動かします。角$\theta$が小さいとき、テイラー展開の1次で打ち切ると

$$
\cos\theta\, i + \sin\theta\, j \approx i + \theta j
$$

となります。球面上の移動が、点$i$で球面に接する平面（接平面）の中を$\theta j$だけずれる並進のように見えています。「無限小の領域では回転が並進として扱える」とはこの現象のことです。ただし$\approx$の背後には$O(\theta^2)$の誤差があり、この段階ではあくまで近似の話です。

## 無限小の読み替え

ここで「小さい角」の扱いを切り替えます。$\theta$を限りなく小さい実数と見なす代わりに、$\theta = \varepsilon a$とおきます。$a$は普通の実数、$\varepsilon$は$i,j,k$と可換で

$$
\varepsilon^2 = 0
$$

を満たす形式的な基底です。$0$でないのに2乗すると$0$になる、実数には存在しない元を新たに添加するわけです。すると指数関数の級数は2次以降の項が文字通り$0$になり、2項で厳密に閉じます。

$$
e^{k\varepsilon a/2} = 1 + \frac{k\varepsilon a}2
$$

近似記号$\approx$の出る幕はありません。1次で打ち切っても誤差が出ないため、以降の計算はすべて厳密な恒等式の変形になります。

点の側も$\varepsilon$を使って書きます。接平面上の点を、基点$i$からのずれ$x \in \operatorname{span}(j,k)$によって

$$
i + \varepsilon x
$$

と表します。四元数共役を係数ごとに施す共役を$\overline{\phantom{q}}$と書くと、$x$が$i$と直交すること（$ix + xi = -2\langle i,x \rangle = 0$）から

$$
(i + \varepsilon x)\overline{(i + \varepsilon x)} = (i + \varepsilon x)(-i - \varepsilon x) = 1 - \varepsilon(ix + xi) = 1
$$

となり、$i + \varepsilon x$は厳密にノルム$1$を保ちます。接平面が、近似としてではなく「球面の$\varepsilon$近傍」として球面上に乗っている形です。

この設定で、角を$\varepsilon$倍した回転子による共役作用を計算します。$ki - ik = 2j$より

$$
e^{k\varepsilon a/2} (i + \varepsilon x) e^{-k\varepsilon a/2}
= \left(1 + \frac{k\varepsilon a}2\right)(i + \varepsilon x)\left(1 - \frac{k\varepsilon a}2\right)
= i + \varepsilon x + \frac{\varepsilon a}2 (ki - ik)
= i + \varepsilon (x + aj)
$$

有限の角の回転が点を球面上の遠くへ運んだのに対し、角を$\varepsilon$倍した回転子は、接平面のすべての点を一斉に$aj$だけずらします。剰余項はありません。回転が、文字通りの**並進**に読み替えられました。同様に、$ji - ij = -2k$より$e^{j\varepsilon b/2}$は$-bk$方向の並進を与え、2方向を組み合わせれば接平面の任意の並進が得られます。

## 接平面の剛体変換

一方、基点$i$を通る軸の回転$e^{i\theta/2}$は、角を$\varepsilon$倍することなくそのまま使えます。

$$
e^{i\theta/2} (i + \varepsilon x) e^{-i\theta/2} = i + \varepsilon\, e^{i\theta/2} x e^{-i\theta/2}
$$

基点は動かず、ずれ$x$が接平面の中で角$\theta$だけ回転します。まとめると、通常の回転子から次の2種類の変換が得られました。

* 基点を通る軸の回転$e^{i\theta/2}$は、接平面の**回転**
* 基点と直交する軸の$\varepsilon$倍角の回転$e^{k\varepsilon a/2},\ e^{j\varepsilon b/2}$は、接平面の**並進**

これらの積の全体は、平面の合同変換のうち向きを保つもの、すなわち2次元剛体変換群$\mathrm{SE}(2)$を二重被覆する群をなします。「無限小の領域では回転が並進として扱える」という観察は、$\varepsilon^2=0$の下では近似ではなく厳密な代数法則です。

&&&rem
球面三角法の記事の無限小三角形は、まさにこの接平面の中で「進んで、曲がって」を繰り返して一周する計算でした。そこでは$O(\varepsilon^2)$の剰余項を持ち回りながら平面三角形の閉合条件を導きましたが、本稿の読み替えの下では、同じ計算が剰余項なしで成り立ちます。
&&&

ところで、ここまでの計算に現れた係数は、いずれも「実数＋$\varepsilon\times$実数」の形をしていました。次節でこの数の体系を切り出して整備し、そのうえで「基点が$i$に固定され、並進が接平面の2次元に限られる」という縛りを外していきます。

# 二重数

前節で使った$\varepsilon$の規則を、まず実数係数の範囲で数の体系として切り出します。

&&&def 二重数
実数$x,y$と、$\varepsilon^2=0$を満たす基底$\varepsilon$によって

$$
z = x + \varepsilon y
$$

と表される数を**二重数** (dual number) と呼び、その全体を$\mathbb D$と書きます。積は分配法則と$\varepsilon^2=0$から定まります。

$$
(x_1 + \varepsilon y_1)(x_2 + \varepsilon y_2) = x_1 x_2 + \varepsilon (x_1 y_2 + y_1 x_2)
$$
&&&

$\varepsilon$は$0$でないのに2乗すると$0$になる冪零元です。そのため二重数には零因子があります（$(\varepsilon a)(\varepsilon b)=0$）。複素数の共役と同様に共役$z^* := x - \varepsilon y$を定めると、$zz^* = x^2$となります。$x \neq 0$のとき、そしてそのときに限り$z$は可逆です。

$$
(x + \varepsilon y)^{-1} = \frac{z^*}{zz^*} = \frac1x - \varepsilon \frac{y}{x^2}
$$

## 関数の拡張

二重数の最も重要な性質は、関数の引数に$\varepsilon$の項を入れると微分が現れることです。

&&&fml
多項式または収束する冪級数で定義される関数$f$に対して

$$
f(x + \varepsilon y) = f(x) + \varepsilon y f'(x)
$$
&&&

&&&prf
単項式で確認すれば、線形性により多項式・冪級数に拡張されます。二項展開で$\varepsilon^2$以上の項が消えることから

$$
(x + \varepsilon y)^n = x^n + \varepsilon n x^{n-1} y
$$

となり、$(x^n)' = nx^{n-1}$と比較すれば主張が従います。
&&&

特に指数関数では、前節でも使った次の関係が厳密に成り立ちます。

$$
e^{\varepsilon y} = 1 + \varepsilon y
$$

前節の読み替えを支えていたのは、この「1次で打ち切っても誤差が出ない」性質です。

&&&rem
数値計算では、この性質を利用して関数値と微分係数を同時に求める手法が**自動微分**（前進モード）として実用されています。
&&&

## 2次元実代数の三分法

実数に基底を1つ添加して得られる2次元の代数は、添加する基底の平方によって3種類に分かれます。

| 基底の平方 | 代数 | $e^{\theta \cdot \text{基底}}$の作用 | 単位元の軌跡 $zz^*=1$ |
|:--:|:--:|:--:|:--:|
| $i^2 = -1$ | 複素数$\mathbb C$ | 回転 $\cos\theta + i\sin\theta$ | 円 |
| $j^2 = +1$ | 分解型複素数$\mathbb C'$ | 双曲回転 $\cosh\theta + j\sinh\theta$ | 双曲線 |
| $\varepsilon^2 = 0$ | 二重数$\mathbb D$ | せん断 $1 + \varepsilon\theta$ | 2直線 $x = \pm 1$ |

円と双曲線が退化して2本の平行な直線になったものが二重数の「単位円」です。$e^{\varepsilon\theta}$を掛ける操作は、実部を保ったまま$\varepsilon$部を$\theta x$だけずらす、せん断（シアー）として働きます。回転角が「進んでも戻ってこない」この構造が、前節で回転が並進として振る舞ったことの正体です。

# 二重四元数の定義

四元数の係数を二重数に取り替えます。前節までの計算で実際に使っていた体系を、ここで正式に定義します。

&&&def 二重四元数
四元数$p, q \in \mathbb H$と、$\varepsilon^2 = 0$かつ$i,j,k$と可換な基底$\varepsilon$によって

$$
\sigma = p + \varepsilon q
$$

と表される数を**二重四元数** (dual quaternion) と呼びます。実ベクトル空間としては8次元で、テンソル積として$\mathbb H \otimes \mathbb D$と書けます。$p$を**実部分**、$q$を**二重部分**と呼びます。
&&&

これは、四元数の係数を別の2次元代数に取り替える系列の3つ目にあたります。

| 係数 | テンソル積 | 名称 |
|:--:|:--:|:--:|
| $\mathbb C$ | $\mathbb H \otimes \mathbb C$ | 双四元数 (biquaternion) |
| $\mathbb C'$ | $\mathbb H \otimes \mathbb C'$ | 分解型双四元数 (split-biquaternion) |
| $\mathbb D$ | $\mathbb H \otimes \mathbb D$ | 二重四元数 (dual quaternion) |

&&&rem 用語について
「双四元数」はハミルトンによる複素係数の四元数$\mathbb H \otimes \mathbb C$の名称で、パウリ行列の記事や4次元回転の記事で扱いました。一方、クリフォードが1873年に "biquaternion" の名で導入した対象は、今日の二重四元数にあたります。[[clifford-biq]] 名称の衝突を避けるため、本稿では dual quaternion の訳語として「二重四元数」を使います。なお、二重数 (dual number) は「双対数」とも訳されます。
&&&

&&&rem 記号について
本稿の$\varepsilon$は$\varepsilon^2=0$を満たす形式的な基底であり、無限小の数ではありません。球面三角法の記事の$\varepsilon$は無限小の展開パラメーター、4次元回転の記事の$\varepsilon = \frac{1+\omega}2$は冪等元（$\varepsilon^2 = \varepsilon$）で、それぞれ別物です。
&&&

## 3つの共役

二重四元数には、四元数由来の共役と二重数由来の共役という独立な2つの共役があり、その合成を含めて3種類の共役が定まります。$\sigma = p + \varepsilon q$に対して

$$
\overline\sigma := \overline p + \varepsilon \overline q, \quad
\sigma^* := p - \varepsilon q, \quad
\overline{\sigma^*} = \overline p - \varepsilon \overline q
$$

$\overline{(\sigma^*)} = (\overline\sigma)^*$なので、合成の順序は問いません。積との関係は次のようになります。

$$
\overline{\sigma\tau} = \overline\tau\,\overline\sigma, \quad
(\sigma\tau)^* = \sigma^*\tau^*, \quad
\overline{(\sigma\tau)^*} = \overline{\tau^*}\,\overline{\sigma^*}
$$

四元数共役は積の順序を反転し（反準同型）、二重数共役は保ちます（準同型）。したがって合成は反準同型です。3種類も共役があるのは冗長に思えるかもしれませんが、剛体変換の表現では3つ目の$\overline{\sigma^*}$が本質的な役割を果たします。

## ノルムと単位二重四元数

四元数共役とのペアリングでノルムを定めます。

$$
\sigma\overline\sigma = (p + \varepsilon q)(\overline p + \varepsilon \overline q)
= p\overline p + \varepsilon (p\overline q + q\overline p)
= |p|^2 + 2\varepsilon \langle p, q \rangle
$$

結果は実数ではなく**二重数**になることに注意してください。$\tau\overline\tau$が中心的な二重数スカラーであることから、乗法性$(\sigma\tau)\overline{(\sigma\tau)} = \sigma(\tau\overline\tau)\overline\sigma = (\sigma\overline\sigma)(\tau\overline\tau)$が成り立ちます。

&&&def 単位二重四元数
$\sigma\overline\sigma = 1$を満たす二重四元数を**単位二重四元数**と呼びます。上の計算から、この条件は2つの実条件

$$
|p| = 1, \quad \langle p, q \rangle = 0
$$

と同値です。このとき$\sigma^{-1} = \overline\sigma$となります。
&&&

8次元から2つの条件で絞られるため、単位二重四元数のなす群は$8-2=6$次元です。これは3次元の剛体変換のなす群$\mathrm{SE}(3)$の次元（回転3＋並進3）と一致します。次節でこの一致が偶然でないことを見ます。

# 剛体変換の表現

接平面モデルでは、基点が$i$に固定され、並進は接平面の2次元に限られていました。本節では基点を取り替えて、3次元空間全体の並進と剛体変換を表します。

## 点の埋め込み

3次元空間の点$x$（純虚四元数）を、二重四元数

$$
1 + \varepsilon x
$$

として埋め込みます。接平面モデルの$i + \varepsilon x$と同じ形ですが、基点を$i$から$1$に取り替え、$x$の範囲を接平面の2次元から純虚四元数の3次元全体に広げています。先頭の$1$は一見不自然ですが、同次座標の「重み」に相当し、この形が保たれることが変換の整合性を担保します。後で見るPGAでは、この埋め込みが自然に導かれます。

## 並進と回転

まず並進です。接平面モデルで並進を担ったのは、角を$\varepsilon$倍した回転子$1 + \frac\varepsilon2(ak)$でした。そこで純虚四元数$t$に対して同じ形の元

$$
T = 1 + \frac\varepsilon2 t
$$

をとり、$t$方向の並進を担わせます。ただし、挟み方には修正が必要です。接平面モデルで共役が並進を生んだのは、基点$i$が$j,k$と反交換したからでした。基点$1$はすべての元と可換なので、同じ共役では

$$
T (1 + \varepsilon x) \overline{T} = 1 + \varepsilon x
$$

と、点がまったく動きません。ここで3つ目の共役の出番です。$\overline t = -t$より$\overline{T^*} = 1 + \frac\varepsilon2 t = T$となることに注意して、$\overline{T^*}$で挟むと

$$
T (1 + \varepsilon x) \overline{T^*}
= \left(1 + \frac\varepsilon2 t\right)(1 + \varepsilon x)\left(1 + \frac\varepsilon2 t\right)
= 1 + \varepsilon \left(x + \frac t2 + \frac t2\right)
= 1 + \varepsilon (x + t)
$$

となり、$t$だけの並進が実現します。係数の$\frac12$は、回転子の半角と同じく「両側から半分ずつ」寄与するために必要です。

回転は通常の回転子$r$（単位四元数）そのままです。$r$は二重部分を持たないので$\overline{r^*} = \overline r$であり

$$
r (1 + \varepsilon x) \overline r = 1 + \varepsilon\, r x \overline r
$$

と、実部分の$1$を保ったまま$x$だけが回転します。

これらを合成すると、一般の剛体変換が得られます。

&&&thm 剛体変換のサンドイッチ公式
単位回転子$r$と純虚四元数$t$に対して

$$
\sigma = r + \frac\varepsilon2 t r = \left(1 + \frac\varepsilon2 t\right) r
$$

とおくと、$\sigma$は単位二重四元数であり

$$
\sigma (1 + \varepsilon x) \overline{\sigma^*} = 1 + \varepsilon (r x \overline r + t) \tag{1}
$$

すなわち「$r$で回転してから$t$だけ並進する」剛体変換を表します。
&&&

&&&prf
$\overline{(tr)} = \overline r\, \overline t = -\overline r t$より$\overline\sigma = \overline r - \frac\varepsilon2 \overline r t$、二重部分の符号を反転して$\overline{\sigma^*} = \overline r + \frac\varepsilon2 \overline r t$です。まず$\sigma\overline\sigma = 1 + \frac\varepsilon2(t r \overline r - r \overline r t) = 1$なので単位です。挟み込みを直接展開します。

$$
\begin{aligned}
\sigma (1 + \varepsilon x) \overline{\sigma^*}
&= \left(r + \varepsilon \left(r x + \frac{t r}2\right)\right)\left(\overline r + \frac\varepsilon2 \overline r t\right) \\
&= r \overline r + \varepsilon \left(\frac12 r \overline r t + r x \overline r + \frac12 t r \overline r\right) \\
&= 1 + \varepsilon (r x \overline r + t)
\end{aligned}
$$
&&&

&&&rem なぜ第3の共役が必要か
一般の$\sigma$を通常の四元数共役$\overline\sigma$で挟むと、並進が両側から逆符号で寄与して打ち消されてしまいます。

$$
\sigma (1 + \varepsilon x) \overline{\sigma}
= \left(r + \varepsilon \left(r x + \frac{t r}2\right)\right)\left(\overline r - \frac\varepsilon2 \overline r t\right)
= 1 + \varepsilon\, r x \overline r
$$

回転だけが残り、並進の情報が完全に失われます。二重数共役$*$を組み合わせて二重部分の符号を反転させることで、両側の寄与が同符号になり、並進が生き残ります。3つの共役の存在意義がここにあります。
&&&

## 合成とSE(3)の二重被覆

2つの剛体変換$\sigma_1 = (1 + \frac\varepsilon2 t_1) r_1$、$\sigma_2 = (1 + \frac\varepsilon2 t_2) r_2$を合成します。操作の順序は回転子と同じく右から左です。$r_2 t_1 = (r_2 t_1 \overline{r_2}) r_2$を使って$r_2$を$t_1$の右へ通過させると

$$
\sigma_2 \sigma_1 = r_2 r_1 + \frac\varepsilon2 (t_2 + r_2 t_1 \overline{r_2})\, r_2 r_1 \tag{2}
$$

となります。これは「回転は合成$r_2 r_1$、並進は$t_1$を$r_2$で回してから$t_2$を加える」という$\mathrm{SE}(3)$の合成規則そのものです。実際、(1)を2回適用した結果$x \mapsto r_2(r_1 x \overline{r_1} + t_1)\overline{r_2} + t_2$と一致します。

また、$-\sigma$は$\sigma$と同じ変換を与えます（(1)で符号が両側から2回掛かるため）。逆に同じ変換を与える単位二重四元数は$\pm\sigma$に限られるので、単位二重四元数のなす群から$\mathrm{SE}(3)$への対応は2対1の全射準同型、すなわち単位二重四元数は$\mathrm{SE}(3)$の**二重被覆**をなします。単位四元数が$\mathrm{SO}(3)$を二重被覆する構図が、並進込みでそのまま持ち上がった形です。

球面三角法の記事では「無限小の領域において、回転が並進として扱える」と述べるに留まりました。二重四元数はこの観察を代数として実装したものです。$\varepsilon$を「無限小」と読めば接空間（曲面に接する平面）の中の近似計算ですが、$\varepsilon^2=0$という規則だけを残して形式的に扱えば、接空間の縛りは外れ、空間全体で厳密に成り立つ剛体変換の代数が得られます。

# ねじ運動

剛体変換の表示$\sigma = (1 + \frac\varepsilon2 t) r$は「回転してから並進」という分解でした。一方、通常の回転子が軸と角度によって$e^{l\theta/2}$と1本の指数関数で書けたように、剛体変換も1本の指数関数で書けます。その幾何学的な意味が**ねじ運動** (screw motion)、すなわち「1本の直線を軸として回転しながら、軸方向に滑る」運動です。

&&&def 二重角とねじ軸
実数$\theta, d$に対して$\hat\theta := \theta + \varepsilon d$を**二重角**と呼びます。また、単位純虚四元数$l$と、$\langle l, m \rangle = 0$を満たす純虚四元数$m$に対して

$$
L := l + \varepsilon m
$$

を**ねじ軸**と呼びます。$L^2 = l^2 + \varepsilon(lm + ml) = -1 - 2\varepsilon\langle l,m \rangle = -1$となり、$L$は虚数単位と同じ平方を持ちます。
&&&

&&&rem プリュッカー座標
組$(l, m)$は空間内の直線を表します。点$x_0$を通り方向$l$を向く直線に対して$m := x_0 \times l$とおくと、$m$は直線上の点の取り方によりません（$(x_0 + sl) \times l = x_0 \times l$）。$l$は直線の方向、$m$は原点から見た直線の「モーメント」で、この表し方は直線のプリュッカー座標と呼ばれます。逆に$\langle l,m \rangle = 0$なら$x_0 = l \times m$がこの直線上の点（原点からの垂線の足）を与えます。
&&&

&&&thm ねじ運動の指数表示
$L^2 = -1$なので、オイラーの公式と同じ級数計算により

$$
\exp\left(\frac{\hat\theta}2 L\right) = \cos\frac{\hat\theta}2 + L \sin\frac{\hat\theta}2 \tag{3}
$$

が成り立ちます。二重角の三角関数を$f(x + \varepsilon y) = f(x) + \varepsilon y f'(x)$で展開すると

$$
\exp\left(\frac{\hat\theta}2 L\right)
= \underbrace{\left(\cos\frac\theta2 + l \sin\frac\theta2\right)}_{\text{実部分}}
+ \varepsilon \underbrace{\left(-\frac d2 \sin\frac\theta2 + \frac d2 \cos\frac\theta2\, l + \sin\frac\theta2\, m\right)}_{\text{二重部分}}
$$

これは、直線$(l, m)$を軸とする角$\theta$の回転と、軸方向への距離$d$の滑りを合成した剛体変換（ねじ運動）を表します。
&&&

&&&prf
$m = x_0 \times l$とおき、ねじ運動を剛体変換の合成として組み立てて(3)と一致することを示します。軸を原点に平行移動して回転し、元に戻してから滑る、という合成です。並進$T(v) := 1 + \frac\varepsilon2 v$と回転子$r := e^{l\theta/2} = \cos\frac\theta2 + l\sin\frac\theta2$を使って

$$
T(dl)\, T(x_0)\, r\, T(-x_0)
$$

を計算します。純虚四元数$u, v$の積の反対称部分が外積を与えること（$uv - vu = 2\,u \times v$）に注意すると

$$
T(x_0)\, r\, T(-x_0) = r + \frac\varepsilon2 (x_0 r - r x_0)
= r + \frac\varepsilon2 \sin\frac\theta2 (x_0 l - l x_0)
= r + \varepsilon \sin\frac\theta2\, m
$$

さらに$T(dl)$を掛けると、$l r = \cos\frac\theta2\, l + \sin\frac\theta2\, l^2$より

$$
T(dl) \left(r + \varepsilon \sin\frac\theta2\, m\right)
= r + \varepsilon \left(-\frac d2 \sin\frac\theta2 + \frac d2 \cos\frac\theta2\, l + \sin\frac\theta2\, m\right)
$$

となり、(3)の展開と一致します。
&&&

特殊な場合を確認します。$d = 0,\ m = 0$（原点を通る軸、滑りなし）では通常の回転子$e^{l\theta/2}$に戻ります。$\theta = 0$では$\frac{\hat\theta}2 L = \frac{\varepsilon d}2 l$（$\varepsilon^2=0$で$m$の項は消えます）となり

$$
\exp\left(\frac{\varepsilon d}2 l\right) = 1 + \frac\varepsilon2 dl
$$

と、方向$l$への距離$d$の純並進が得られます。並進とは「無限遠に去った軸の回りのねじ運動」の極限であり、回転角$\theta$が消えて滑り$d$だけが残った形です。

&&&thm 剛体変換のねじ表示（シャールの定理）
任意の単位二重四元数$\sigma = p + \varepsilon q$は、符号$\pm$を除いて(3)の形に書けます。したがって、任意の3次元剛体変換はあるねじ運動として実現できます。
&&&

&&&prf
$\pm$の符号を選んで$\mathrm{Re}(p) \geq 0$とします。$p \neq \pm1$のとき、実部分$p$は単位四元数なので$\theta = 2\operatorname{atan2}(|\mathrm{Im}(p)|, \mathrm{Re}(p))$、$l = \mathrm{Im}(p)/|\mathrm{Im}(p)|$と一意に分解できます。二重部分を(3)の展開と比較して

$$
d = -\frac{2\,\mathrm{Re}(q)}{\sin(\theta/2)}, \quad
m = \frac1{\sin(\theta/2)}\left(\mathrm{Im}(q) - \frac d2 \cos\frac\theta2\, l\right)
$$

と$d, m$が求まります。残る要件$\langle l, m \rangle = 0$は、単位条件$\langle p, q \rangle = 0$から従います。実際、$\langle p,q \rangle = \cos\frac\theta2\,\mathrm{Re}(q) + \sin\frac\theta2 \langle l, \mathrm{Im}(q) \rangle = 0$に上の$d$を代入すると$\langle l, \mathrm{Im}(q) \rangle = \frac d2 \cos\frac\theta2$となり、$\langle l, m \rangle = 0$が得られます。$p = \pm 1$のときは前述の純並進（$\theta = 0$）の場合です。
&&&

&&&ex 回転と並進の合成のねじ軸
$k$軸回りの角$\theta$の回転と、それに垂直な並進$t = ai$の合成を考えます。

$$
\sigma = \left(1 + \frac\varepsilon2 ai\right) e^{k\theta/2}
= e^{k\theta/2} + \varepsilon\, \frac a2 \left(\cos\frac\theta2\, i - \sin\frac\theta2\, j\right)
$$

（$ik = -j$を使いました。）二重部分の実部（スカラー部）は$0$なので$d = 0$、すなわち滑りのない純粋な回転です。ねじ軸は

$$
m = \frac a2 \left(\cot\frac\theta2\, i - j\right), \quad
x_0 = l \times m = k \times m = \frac a2\, i + \frac a2 \cot\frac\theta2\, j
$$

つまりこの変換は、原点からずれた点$x_0$を通る鉛直軸回りの、角$\theta$の回転です。「平面上の回転と並進の合成は、中心をずらした回転になる」という初等幾何の事実が、ねじ軸の抽出として機械的に再現されました。$\theta \to 0$とすると$\cot\frac\theta2 \to \infty$で軸は無限遠に逃げ、純並進に退化します。
&&&

&&&rem
2つの姿勢の間の補間は、回転子では球面線形補間 (SLERP) として知られますが、二重四元数では(3)の二重角$\hat\theta$を等分するねじ運動の補間 (ScLERP) に拡張され、回転と並進が同一のねじ軸上で同時に補間されます。回転子の補間は次の記事で扱います。
&&&

# クリフォード代数との対応

四元数のテンソル積の記事では、非退化な計量（生成元の平方が$\pm1$）を持つクリフォード代数を構成しました。二重四元数の$\varepsilon^2 = 0$は、平方が$0$になる生成元、すなわち**退化した計量**を許すことに対応します。

&&&def 退化クリフォード代数 $C\ell_{3,0,1}(\mathbb R)$
生成元$e_1, e_2, e_3$（$e_a^2 = +1$）と$e_0$（$e_0^2 = 0$）が互いに反交換するとして生成される$2^4 = 16$次元の代数を$C\ell_{3,0,1}(\mathbb R)$と書きます。符号数の3つ組$(p,q,r) = (3,0,1)$は、平方が$+1, -1, 0$となる生成元の個数を表します。
&&&

偶数個の生成元の積が張る偶部分代数$C\ell^+_{3,0,1}(\mathbb R)$は、スカラー$1$、6つの2ベクトル$e_ae_b$、擬スカラー$e_0e_1e_2e_3$の8次元です。これが二重四元数と同型になります。

&&&thm
次の対応は代数の同型$\mathbb H \otimes \mathbb D \cong C\ell^+_{3,0,1}(\mathbb R)$を与えます。

$$
\begin{array}{c|cccc}
 & 1 & i & j & k \\ \hline
1 & 1 & e_3e_2 & e_1e_3 & e_2e_1 \\
\varepsilon & e_0e_1e_2e_3 & e_0e_1 & e_0e_2 & e_0e_3
\end{array}
$$
&&&

&&&prf
代表的な積を確認します（すべての基底の組み合わせは検証コードで確認します）。まず$e_1, e_2, e_3$の部分は通常の3次元クリフォード代数なので、その偶部分代数が四元数と同型になるのは既知の構図です。

$$
(e_3e_2)^2 = -e_3e_3e_2e_2 = -1, \quad
(e_3e_2)(e_1e_3) = e_3e_2e_1e_3 = e_3e_3e_2e_1 = e_2e_1
$$

は$i^2 = -1,\ ij = k$を再現します。$E := e_0e_1e_2e_3$とおくと、$e_0$の平方が$0$なので

$$
E^2 = e_0e_1e_2e_3\,e_0e_1e_2e_3 = -e_0e_0(e_1e_2e_3)^2 = 0
$$

と$\varepsilon^2 = 0$を再現します（$e_0$を左端に集める並べ替えで符号が付きますが、$e_0e_0 = 0$の因子を含むため結果は$0$です）。また

$$
E\,(e_3e_2) = e_0e_1e_2e_3e_3e_2 = e_0e_1e_2e_2 = e_0e_1
$$

は$\varepsilon \cdot i = \varepsilon i$に対応します。整合性の確認として$(e_0e_1)^2 = -e_0e_0e_1e_1 = 0$は$(\varepsilon i)^2 = \varepsilon^2 i^2 = 0$と一致します。
&&&

$\varepsilon$に対応する擬スカラー$E = e_0e_1e_2e_3$が二重四元数の$\varepsilon$として振る舞うには、冪零性に加えて「すべての元と可換」であることが必要です。4次元の擬スカラーは、1ベクトルとは反交換（$e_a E = -E e_a$）ですが、2ベクトルとは可換です。したがって$E$は代数全体の中心には入りませんが、**偶部分代数の中では中心的**であり、$\varepsilon$の役割を果たせます。

## 擬スカラーの三分法

4次元クリフォード代数の偶部分代数は$\mathbb H \otimes (\text{2次元代数})$の形になり、どの2次元代数が現れるかは擬スカラー$\omega$の平方で決まります。

| クリフォード代数 | $\omega^2$ | 偶部分代数 |
|:--:|:--:|:--:|
| $C\ell_{3,1}(\mathbb R)$ | $-1$ | 双四元数 $\mathbb H \otimes \mathbb C$ |
| $C\ell_{4,0}(\mathbb R)$ | $+1$ | 分解型双四元数 $\mathbb H \otimes \mathbb C' \cong \mathbb H \oplus \mathbb H$ |
| $C\ell_{3,0,1}(\mathbb R)$ | $0$ | 二重四元数 $\mathbb H \otimes \mathbb D$ |

二重数の節で見た2次元実代数の三分法（$-1$: 複素数、$+1$: 分解型複素数、$0$: 二重数）が、偶部分代数の分類としてそのまま再登場しています。$\omega^2 = +1$の場合に冪等元$\frac{1\pm\omega}2$による直和分解$\mathbb H \oplus \mathbb H$が生じることは4次元回転の記事で見た通りです。$\omega^2 = 0$の場合は直和に分解する代わりに、冪零なイデアル$\varepsilon\mathbb H$（並進部分）が現れます。回転群$\mathrm{Spin}(4)$を担った代数と剛体変換を担う代数が、擬スカラーの平方ひとつを付け替えた隣人どうしだったわけです。

## PGAへの入り口

本稿では$C\ell_{3,0,1}(\mathbb R)$の偶部分代数だけを使いましたが、16次元の代数全体を使う枠組みは**射影幾何代数** (PGA) と呼ばれます。そこでは1ベクトルが平面を、2ベクトルが直線を、3ベクトルが点を表し、平面による鏡映を挟み込みで表すところから、回転も並進もすべての剛体変換が組み立てられます。点の埋め込み$1 + \varepsilon x$や、ねじ軸$l + \varepsilon m$が直線を表したことは、この枠組みの断片が二重四元数に映り込んだものです。さらに計量の取り替えで球面や反転を扱う共形幾何代数 (CGA) へも道が続きます。これらは今後の記事で扱う予定です。

# まとめ

- 「無限小の領域では回転が並進として扱える」という観察は、回転角を$\varepsilon$倍（$\varepsilon^2=0$）に読み替えることで厳密な代数法則に昇格します。角を$\varepsilon$倍した回転子は接平面の点$i + \varepsilon x$を並進させ、基点を通る軸の回転と合わせて、接平面の剛体変換群$\mathrm{SE}(2)$が通常の回転子の共役作用として実現します。
- 二重数は$\varepsilon^2 = 0$を満たす基底を実数に添加した代数で、$f(x + \varepsilon y) = f(x) + \varepsilon y f'(x)$により「1次で打ち切っても誤差が出ない」無限小計算を厳密に実装します。
- 四元数の係数を二重数に取り替えた二重四元数$\mathbb H \otimes \mathbb D$は、双四元数$\mathbb H \otimes \mathbb C$・分解型双四元数$\mathbb H \otimes \mathbb C'$と並ぶテンソル積の族の一員です。
- 基点を$i$から$1$に取り替えた点$1 + \varepsilon x$を、剛体変換$\sigma = (1 + \frac\varepsilon2 t) r$の第3の共役$\overline{\sigma^*}$で挟むことで、3次元全体の「回転＋並進」が実現します。単位二重四元数は$\mathrm{SE}(3)$を二重被覆します。
- 任意の剛体変換は二重角とねじ軸による指数関数$\exp(\frac{\hat\theta}2 L)$、すなわちねじ運動として書けます（シャールの定理）。並進は軸が無限遠に去ったねじ運動の極限です。
- 二重四元数は退化クリフォード代数の偶部分代数$C\ell^+_{3,0,1}(\mathbb R)$と同型であり、擬スカラーの平方$-1, +1, 0$に応じて双四元数・分解型双四元数・二重四元数が現れるという三分法が成り立ちます。これはPGA・CGAへの入り口です。

本稿の計算はすべて検証コードで数値的に確認しています。
