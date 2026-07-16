テレンス・タオ教授の四元数と球面三角法に関するブログ記事[[tao-h]]の主な内容を、前提知識などを補って解説します。一部数式に変更がありますが、結論には影響しません。

## 概要

ハミルトンの四元数は、複素数の非可換な拡張であり、球面三角法の基本法則を導くための強力なツールとなります。本稿では、四元数の基本的な性質から出発し、球面上の回転作用を通じて、球面三角法の諸公式（余弦定理、正弦定理、5要素の公式）を導出します。また、その応用として日の出方程式についても解説します。

## 四元数の基礎

四元数は$t + xi + yj + zk$の形の数です。ここで$t,x,y,z$は実数であり、$i,j,k$は$-1$の平方根で、以下の関係を満たします。

$$
i^2 = j^2 = k^2 = -1, \quad ij=k, \quad jk=i, \quad ki=j
$$

これらは非可換（$ij \neq ji$）ですが、結合法則は成り立ちます。四元数全体は斜体（可除環）を形成し、すべての非ゼロ四元数は逆元を持ちます。

### 共役と分解

複素数と同様に、四元数には共役があります。

$$
\overline{t+xi+yj+zk} := t-xi-yj-zk
$$

これは反準同型となります（$\overline{qr} = \overline{r}\ \overline{q}$）。四元数$q$は実部と虚部に分解できます。

$$
\mathrm{Re} (q) := \frac{q + \overline{q}}{2}, \quad \mathrm{Im} (q) := \frac{q - \overline{q}}{2}
$$

&&&rem 虚部
複素数での$\mathrm{Im}$は虚数単位$i$の係数のみを表しますが、四元数には複数の虚数単位$i,j,k$があるため、虚数単位も含みます。

- 複素数:$\mathrm{Im}(1+2i) = 2$
- 四元数:$\mathrm{Im}(1+2i+3j+4k) = 2i+3j+4k$
&&&

### 内積とノルム

内積は次のように定義されます。

$$
\langle q, r \rangle := \mathrm{Re} (q \overline{r})
$$

これは対称かつ正定値であり、$1,i,j,k$は正規直交基底を形成します。ノルムは

$$
|q| = \sqrt{\langle q,q \rangle} = \sqrt{t^2 + x^2 + y^2 + z^2}
$$

で与えられ、乗法的な性質$|qr| = |q| |r|$を持ちます。

単位四元数の集合$\mathrm{U}(1,\mathbb{H}) = \{ q \in \mathbb{H}: |q|=1\}$は群をなし、$\mathrm{SU}(2)$と同型です。

### 随伴性

内積には以下の随伴性が成り立ちます。これは線形代数における随伴作用素の概念に対応します。

$$
\langle qr, s \rangle = \langle q, s\overline{r} \rangle, \quad
\langle rq, s \rangle = \langle q, \overline{r}s \rangle
$$

&&&prf
$$
\begin{aligned}
\langle qr, s \rangle
&= \mathrm{Re}(qr\overline{s})
 = \langle q, \overline{r\overline{s}} \rangle
 = \langle q, s\overline{r} \rangle \\

\langle rq, s \rangle
&= \mathrm{Re}(rq\overline{s})
 = \mathrm{Re}(q\overline{s}r)
 = \langle q, \overline{\overline{s}r} \rangle
 = \langle q, \overline{r}s \rangle
\end{aligned}
$$
&&&

### オイラーの公式

$i,j,k$は$-1$の平方根であるため、実数$\theta$に対して以下のオイラーの公式が成り立ちます。

$$
e^{i\theta} = \cos \theta + i \sin \theta, \quad e^{j\theta} = \cos \theta + j \sin \theta, \quad e^{k\theta} = \cos \theta + k \sin \theta
$$

## 球面上の回転作用

単位四元数$q \in \mathrm{U}(1,\mathbb{H})$は、純虚四元数（$\mathbb{R}^3$と同一視）に対して共役によって作用します。

$$
v \mapsto q v \overline{q}
$$

この作用は$\mathbb{R}^3$の内積と向きを保つため、回転変換となります。

&&&rem
群論における**共役**（conjugation）は通常$gxg^{-1}$の形を指します。単位四元数においては$\overline{q} = q^{-1}$が成り立つため、$qv\overline{q}$は$qvq^{-1}$と等価であり、共役作用と呼ばれます。なお、四元数$\overline{q}$自体も「共役（四元数）」と呼ばれるため、文脈による区別が必要です。
&&&

### 二重被覆と半角

この作用において、$q$と$-q$は同じ回転を引き起こします。

$$
(-q)v\overline{(-q)} = qv\overline{q}
$$

&&&rem
直感的には、$q$とその共役$\overline{q}$とで両側から挟むことで、符号の違いがキャンセルされるためだと解釈できます。
&&&

特に、$q=1$と$q=-1$はどちらも恒等変換（何もしない変換）となります。群論の言葉では、この作用の核は$\{1, -1\}$であり、作用は忠実ではありません。しかし、同じ変換を引き起こすのは常に$\pm q$のペアに限られるため、この対応は$\mathrm{SU}(2)$から$\mathrm{SO}(3)$への「2対1」の写像（二重被覆）となります。

&&&
群$G$の集合$X$への作用が**忠実** (faithful) であるとは、異なる群元が異なる変換を引き起こすこと、つまり作用が単射であることを意味します。また、恒等変換として作用する群元の集合を作用の**核** (kernel) と呼びます。作用が忠実であることは、核が単位元のみ（$\{1\}$）であることと同値です。本件の場合、核が非自明な元$-1$を含むため、忠実ではありません。
&&&

例えば、$i$軸周りの$\theta$回転は、$e^{i\theta/2}$による共役で表されます。

$$
\begin{aligned}
e^{i\theta/2} i e^{-i\theta/2} &= i \\
e^{i\theta/2} j e^{-i\theta/2} &= (\cos\theta) j + (\sin\theta) k \\
e^{i\theta/2} k e^{-i\theta/2} &= (\cos\theta) k - (\sin\theta) j
\end{aligned}
$$

&&&prf
$$
\begin{aligned}
e^{i\theta/2} i e^{-i\theta/2}
&= e^{i\theta/2} e^{-i\theta/2} i = i \\

e^{i\theta/2} j e^{-i\theta/2}
&= (\cos(\theta/2) + \sin(\theta/2) i) j (\cos(\theta/2) - \sin(\theta/2) i) \\
&= (\cos(\theta/2) j + \sin(\theta/2) k) (\cos(\theta/2) - \sin(\theta/2) i) \\
&= \cos^2(\theta/2) j - \sin^2(\theta/2) j + 2\sin(\theta/2) \cos(\theta/2) k \\
&= (\cos \theta) j + (\sin \theta) k \\

e^{i\theta/2} k e^{-i\theta/2}
&= (\cos(\theta/2) + \sin(\theta/2) i) k (\cos(\theta/2) - \sin(\theta/2) i) \\
&= (\cos(\theta/2) k - \sin(\theta/2) j) (\cos(\theta/2) - \sin(\theta/2) i) \\
&= \cos^2(\theta/2) k - \sin^2(\theta/2) k - 2 \sin(\theta/2) \cos(\theta/2) j \\
&= (\cos \theta) k - (\sin \theta) j
\end{aligned}
$$
&&&

ここで角度が$\theta/2$となっているのは、リー代数の構造に由来します。3次元回転群$\mathrm{SO}(3)$のリー代数$\mathfrak{so}(3)$の標準的な基底$L_x, L_y, L_z$は、交換関係$[L_x, L_y] = L_z$（およびその巡回置換）を満たします。一方、四元数の虚部がなすリー代数においては、交換関係は$[i, j] = ij-ji = 2k$となります。係数の$2$を相殺して$\mathfrak{so}(3)$の構造と一致させるには、基底を$1/2$倍して$i/2, j/2, k/2$とする必要があります（実際、$[i/2, j/2] = k/2$となります）。

正規化された基底$i/2$を採用すれば、角度$\theta$の回転は$e^{\theta(i/2)} = e^{i\theta/2}$となります。

&&&rem
$\theta/2$という「半分の角度」があるというより、回転の生成子が$i$ではなく$i/2$である（基底の方が半分である）と解釈するのが、リー代数の構造としては本質的です。このように、交換関係における係数$2$のズレが、指数写像における$1/2$の因子として現れています。

指数の肩は$i\frac{\theta}{2}$のように角度が半分になっているように表記されることが多いですが、直感的には、共役作用$q v \overline{q}$が左右両側から半分ずつ作用すると解釈できます。具体的には、左側からの作用$q=e^{i\theta/2}$と、右側からの作用$\overline{q} = e^{-i\theta/2}$の効果が合わさって、回転角度は$\theta$となります。
&&&

### スピノルの構造

リー群とリー代数の対応において、生成子$X$によって生成される回転変換は、指数写像$e^{\theta X}$で表されます。四元数の場合、これは$\mathrm{SU}(2)$の元であり、幾何学的には**スピノル**（spinor）と呼ばれます。スピノルは「1回転（$360^\circ$）しても元に戻らず、2回転（$720^\circ$）して初めて元に戻る」という特異な性質を持っています。（後で具体例を見ます）

生成子$X$は回転軸を表す純虚四元数です。例えば、$X=i/2$は$i$軸周りの回転を表します。$X$は任意の方向を指定することができます。スピノルは$e^{\theta X}$で生成されるため、回転軸（3次元ベクトル）と回転角度の4パラメーターを含みます。

回転軸$n$をノルム$1$の単位純虚四元数として表す場合、回転角度を$\theta$としてスピノルは$e^{\theta n/2}$となります。$n^2=-1$より、指数関数において$n$は虚数単位$i,j,k$と同様に振る舞います。

&&&prf 単位純虚四元数の2乗
正規化された回転軸$n$の成分を$n_x, n_y, n_z\ (n_x^2+n_y^2+n_z^2=1)$とおきます。

$$
\begin{aligned}
n^2
&= (n_x i + n_y j + n_z k)^2 \\
&= n_x^2 i^2 + n_y^2 j^2 + n_z^2 k^2 + n_x n_y (ij+ji) + n_y n_z (jk+kj) + n_z n_x (ki+ik) \\
&= -n_x^2 - n_y^2 - n_z^2 \\
&= -(n_x^2 + n_y^2 + n_z^2) \\
&= -1
\end{aligned}
$$
&&&

スピノルは$e^{\theta n/2} = \cos(\theta/2) + \sin(\theta/2) n$によって生成されることから、三角関数の引数に渡される位相が半分になります。そのため、1回転$\theta=2\pi$を表すスピノルは$-1$となります。

$$
e^{2\pi n/2} = e^{\pi n} = -1
$$

2回転$\theta=4\pi$によって$1$に戻ります。

$$
e^{4\pi n/2} = e^{2\pi n} = 1
$$

これらは、既に見たように$1$と$-1$が恒等変換であることを示しています。

&&&rem
2回転で元に戻る回転は、数学的な概念に留まらず、腕やベルトなどで物理的に再現できます。[[wiki-trick]]
&&&

### 回転の合成と順序

スピノルの積は、ベクトルに対する回転（共役作用）の合成に対応しています。あるベクトル$v$に対して、まず回転$q_1$を適用し、次に回転$q_2$を適用する場合を考えます。

1. $v \mapsto q_1 v \overline{q_1}$
2. $(q_1 v \overline{q_1}) \mapsto q_2 (q_1 v \overline{q_1}) \overline{q_2} = (q_2 q_1) v \overline{(q_2 q_1)}$

このように、回転の合成はスピノルの積$q_2 q_1$として表されます。$v$を包むように操作が重なるため、$q_2 q_1$のように右から左に並べる必要があります。

## 球面三角形とマスター方程式

単位球面$S^2 \subset \mathbb{R}^3$上の球面三角形を考えます。

空間に固定された直交座標系$(i,j,k)$と、原点を中心に自由に回転できる球体を考えます。球体とは独立して固定された針が設置されており、球体上の一点$(1,0,0)$（$i$軸上の点）を指します。

&&&rem
針は球体の回転には連動しません。球体が回転スタンドに置かれており、針は回転スタンドに取り付けられているイメージです。
&&&

スピノルによる回転操作は以下のように対応します。

*  $e^{i\theta/2}$：球体を$i$軸（針の軸）周りに反時計回りで$\theta$回転させます。針が指す点は変わらず、球面の向きだけが変わります。（方向転換に相当）
*  $e^{k\theta/2}$：球体を$k$軸周りに反時計回りで$\theta$回転させます。球面上の点を横方向に移動させるのに使用します。（自転に相当）

### 三角形のトレース

単位球面上には、頂点$A, B, C$を持ち、辺の長さ（中心角）が$c, a, b$、頂点の角が$\alpha, \beta, \gamma$である球面三角形が描かれています。頂点は時計回りに配置されているものとします。

&&&rem
単位球面では辺の長さと中心角（ラジアン）が一致します。また、頂点が反時計回りに配置されている場合、$i$軸周りの回転を逆向きにする必要があります。
&&&

球体を回転させて、各頂点を順番に針の位置まで持ってくる操作を考えます。

1. 針が点$A$を指し、点$B$から見て針が自転方向にある状態からスタートします。
2. $e^{kc/2}$：球体を$k$軸周りに回転させ、距離$c$だけ離れた点$B$を針の位置まで持ってきます。
3. $e^{i(\pi-\beta)/2}$：球体を$i$軸周りに$\pi-\beta$回転させ、次の点$C$から見て針を自転方向に合わせます。
4. $e^{ka/2}$：球体を$k$軸周りに回転させ、距離$a$だけ離れた点$C$を針の位置まで持ってきます。
5. $e^{i(\pi-\gamma)/2}$：球体を$i$軸周りに$\pi-\gamma$回転させ、次の点$A$から見て針を自転方向に合わせます。
6. $e^{kb/2}$：球体を$k$軸周りに回転させ、距離$b$だけ離れた点$A$を針の位置まで持ってきます。
7. $e^{i(\pi-\alpha)/2}$：球体を$i$軸周りに$\pi-\alpha$回転させ、最初の向き（点$B$から見て針が自転方向にある状態）に戻します。

三角形の内角を$\alpha, \beta, \gamma$としたとき、頂点で方向転換する角度は$\pi - \alpha$,$\pi - \beta$,$\pi - \gamma$となります。方向転換は三角形の外側で行われるため、回転角度は内角の**補角**（外角）となります。

### 球面三角法のマスター方程式

三角形を一周して元の位置と向きに戻る操作を、「球体を回転させる操作の積」としてスピノルで表現します。最終的に球面の向きが元通りになったとき、スピノルとしては一周分の回転（$2\pi$回転）が蓄積されています。既に見たように、スピノルは1回転（$360^\circ$）しても元に戻らず、2回転（$720^\circ$）して初めて元に戻るという特異な性質を持っています。そのため、右辺は符号が反転した$-1$となります。

&&&def 球面三角法のマスター方程式
$$
e^{i(\pi-\alpha)/2} e^{kb/2} e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2} e^{kc/2} = -1 \tag{4}
$$
&&&

これは球面三角形の辺$a,b,c$と角$\alpha,\beta,\gamma$が満たすべき制約を課しています。この式から各種公式が導出できるため、本稿では球面三角法の「マスター方程式」と呼びます。

&&&rem
前述のようにスピノルは操作順に右から左へと並べます（例：$q_1$の後に$q_2$を適用する場合は$q_2 q_1$）。元記事では接空間での操作を想定して並べ方が逆になっていますが、本稿ではベクトルと同じように共役作用によって球体を回転させるモデルに変更しました。また、右辺の符号も異なりますが、公式の導出には影響しません。
&&&

右辺が$-1$となることを、具体例で確認します。

&&&ex 無限小三角形（ユークリッド極限）
辺の長さ$a,b,c$を$\varepsilon a, \varepsilon b, \varepsilon c$に置き換え、$\varepsilon$について1次の項までテイラー展開します。$(e^{k\varepsilon x/2} \approx 1 + k\varepsilon x/2)$

$$
\begin{aligned}
&e^{i(\pi-\alpha)/2} \left(1+\frac{k\varepsilon b}{2}\right) e^{i(\pi-\gamma)/2} \left(1+\frac{k\varepsilon a}{2}\right) e^{i(\pi-\beta)/2} \left(1+\frac{k\varepsilon c}{2}\right) \\
= &e^{i(3\pi - (\alpha+\beta+\gamma))/2}
   \left(1 + \frac{k\varepsilon b}{2} e^{i(\pi-\gamma)} e^{i(\pi-\beta)} + \frac{k\varepsilon a}{2} e^{i(\pi-\beta)} + \frac{k\varepsilon c}{2} \right) + O(\varepsilon^2) \\
\end{aligned}
$$

平面三角形では内角の和が$\alpha+\beta+\gamma = \pi$となることを利用します。$\varepsilon$について1次の項まで残します。

$$
\begin{aligned}
&(-1) \left(1 + \frac{k\varepsilon b}{2} e^{i(\pi-\gamma)} e^{i(\pi-\beta)} + \frac{k\varepsilon a}{2} e^{i(\pi-\beta)} + \frac{k\varepsilon c}{2} \right) \\
= &-1 - \frac{k\varepsilon}{2}(b e^{i(\pi-\gamma)} e^{i(\pi-\beta)} + a e^{i(\pi-\beta)} + c)
\end{aligned}
$$

第2項の括弧の中は、複素平面で表現された3辺のベクトルの和であるため、$0$となります。（これは無限小の領域において、回転が並進として扱えることを示唆します）

よって、無限小三角形において(4)は$-1$に収束します。&&&

球面上の任意の三角形は無限小三角形からの連続変形で構成できるため、$-1$は不変に保たれます。

&&&ex 直角二等辺三角形（球の1/8）
球の8分の1を切り取った三角形を考えます。

$$
a = b = c = \frac{\pi}{2}, \quad \alpha = \beta = \gamma = \frac{\pi}{2}
$$

このとき各項は以下のようになります。

$$
e^{k\pi/4} = \frac{1}{\sqrt{2}}(1+k), \quad e^{i(\pi-\pi/2)/2} = e^{i\pi/4} = \frac{1}{\sqrt{2}}(1+i)
$$

これらを掛け合わせると、

$$
e^{i\pi/4} e^{k\pi/4} = \frac{1}{2}(1+i)(1+k) = \frac{1}{2}(1+i-j+k)
$$

これは正規化された純虚四元数$u = \frac{1}{\sqrt{3}}(i-j+k)$を用いて$e^{u\pi/3}$と書けます。

$$
e^{u\pi/3}
= \cos \frac{\pi}{3} + \sin \frac{\pi}{3}u
= \frac{1}{2} + \frac{1}{2}(i-j+k)
$$

式 (4) の左辺はこれの3乗となるため、

$$
(e^{u\pi/3})^3 = e^{u\pi} = -1
$$

となり、同様に$-1$となることが確認できます。
&&&

## 球面余弦定理の導出

式 (4) を変形して、球面余弦定理を導きます。両辺に左から$e^{-kb/2} e^{-i(\pi-\alpha)/2}$、右から$e^{-kc/2}$を掛けます。

$$
e^{i(\pi-\alpha)/2} e^{kb/2} e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2} e^{kc/2} = -1
$$
$$
e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2} = - e^{-kb/2} e^{-i(\pi-\alpha)/2} e^{-kc/2} \tag{5}
$$

次に、(5)の両辺を$i$で共役変換します（$iq\bar{i}$）。$i$と反交換する$k$の成分は符号が反転するため、$e^{k\phi}$は$e^{-k\phi}$になります。

$$
e^{i(\pi-\gamma)/2} e^{-ka/2} e^{i(\pi-\beta)/2} = - e^{kb/2} e^{-i(\pi-\alpha)/2} e^{kc/2} \tag{6}
$$

式 (5) と (6) の、左辺同士、右辺同士の内積を取ります。これらは等しくなります。

**左辺の計算**:
随伴性$\langle qr, s \rangle = \langle q, s\overline{r} \rangle,\ \langle rq, s \rangle = \langle q, \overline{r}s \rangle$を用いて両端の因子を消去します。

$$
\begin{aligned}
&\langle e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2}, e^{i(\pi-\gamma)/2} e^{-ka/2} e^{i(\pi-\beta)/2} \rangle \\
= &\langle e^{ka/2}, e^{-ka/2} \rangle
= \mathrm{Re}(e^{ka/2} e^{ka/2})
= \mathrm{Re}(e^{ka}) = \cos a
\end{aligned}
$$

**右辺の計算**:
同様に随伴性を用いて整理します。
$$
\begin{aligned}
&\langle e^{-kb/2} e^{-i(\pi-\alpha)/2} e^{-kc/2}, e^{kb/2} e^{-i(\pi-\alpha)/2} e^{kc/2} \rangle \\
= &\langle e^{-kb} e^{-i(\pi-\alpha)/2}, e^{-i(\pi-\alpha)/2} e^{kc} \rangle \\
= &\langle e^{i(\pi-\alpha)/2} e^{-kb} e^{-i(\pi-\alpha)/2}, e^{kc} \rangle
\end{aligned}
$$
ここで、$e^{kc} = \cos c + (\sin c) k$です。また、内積の左因子は$e^{-kb} = \cos b - (\sin b) k$を$i$軸周りに$(\pi-\alpha)$回転させたものです。
$$
\begin{aligned}
e^{i(\pi-\alpha)/2} e^{-kb} e^{-i(\pi-\alpha)/2}
&= \cos b - (\sin b) (\cos(\pi-\alpha) k - \sin(\pi-\alpha) j) \\
&= \cos b - (\sin b) (-\cos(\alpha) k - \sin(\alpha) j) \\
&= \cos b + (\sin b \cos \alpha) k + (\sin b \sin \alpha) j
\end{aligned}
$$
これらの内積を計算します。実部同士と$k$成分同士の積の和となります。
$$
\begin{aligned}
&\langle \cos b + (\sin b \cos \alpha) k + (\sin b \sin \alpha) j, \cos c + (\sin c) k \rangle \\
= &\cos b \cos c + \sin b \cos \alpha \sin c
\end{aligned}
$$

左辺と右辺が等しいことから、球面余弦定理が得られます。

&&&thm 球面余弦定理
$$
\cos a = \cos b \cos c + \sin b \sin c \cos \alpha
$$
&&&

微小な三角形の極限では、ユークリッド幾何学の余弦定理に収束します。$a,b,c$を$\varepsilon a, \varepsilon b, \varepsilon c$に置き換え、2次の項までテイラー展開（$\cos x \approx 1 - x^2/2, \sin x \approx x$）を行います。

$$
\begin{aligned}
\cos \varepsilon a &= \cos \varepsilon b \cos \varepsilon c + \sin \varepsilon b \sin \varepsilon c \cos \alpha \\
1 - \frac{(\varepsilon a)^2}{2} &= \left(1 - \frac{(\varepsilon b)^2}{2}\right)\left(1 - \frac{(\varepsilon c)^2}{2}\right) + (\varepsilon b)(\varepsilon c) \cos \alpha \\
1 - \frac{\varepsilon^2 a^2}{2} &= 1 - \frac{\varepsilon^2 b^2}{2} - \frac{\varepsilon^2 c^2}{2} + \varepsilon^2 bc \cos \alpha
\end{aligned}
$$

整理すると、おなじみの公式が得られます。

$$
a^2 = b^2 + c^2 - 2bc \cos \alpha
$$

## 球面正弦定理の導出

次に、式 (5) の両辺$A=B$から得られる関係$\overline{A}iA = \overline{B}iB$を利用します。この両辺と$k$との内積をとります。

$$
\langle \overline{A}iA, k \rangle = \langle \overline{B}iB, k \rangle
$$

これを計算することで、球面正弦定理が得られます。

&&&thm 球面正弦定理
$$
\frac{\sin \alpha}{\sin a} = \frac{\sin \beta}{\sin b}
$$
&&&

&&&prf
$$
\begin{aligned}
\langle \overline{A}iA, k \rangle
&= \langle e^{-i(\pi-\beta)/2} e^{-ka/2} e^{-i(\pi-\gamma)/2} i e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2}, k \rangle \\
&= \langle e^{-ka/2} i e^{ka/2}, e^{i(\pi-\beta)/2} k e^{-i(\pi-\beta)/2} \rangle \\
&= \langle (\cos a) i - (\sin a) j, -(\cos \beta) k - (\sin \beta) j \rangle \\
&= \sin a \sin \beta \\

\langle \overline{B}iB, k \rangle
&= \langle e^{kc/2} e^{i(\pi-\alpha)/2} e^{kb/2} i e^{-kb/2} e^{-i(\pi-\alpha)/2} e^{-kc/2}, k \rangle \\
&= \langle e^{kb/2} i e^{-kb/2}, e^{-i(\pi-\alpha)/2} k e^{i(\pi-\alpha)/2} \rangle \\
&= \langle (\cos b) i + (\sin b) j, -(\cos \alpha) k + (\sin \alpha) j \rangle \\
&= \sin b \sin \alpha
\end{aligned}
$$
$$
\begin{aligned}
\langle \overline{A}iA, k \rangle &= \langle \overline{B}iB, k \rangle \\
\sin a \sin \beta &= \sin b \sin \alpha \\
\frac{\sin \beta}{\sin b} &= \frac{\sin \alpha}{\sin a}
\end{aligned}
$$
&&&

こちらも辺$a,b,c$の無限小極限によって平面上の正弦定理に収束します。

$$
\frac{\sin \alpha}{a} = \frac{\sin \beta}{b}
$$

## 5要素の規則

同様に、$\overline{A}iA = \overline{B}iB$の両辺と$j$との内積をとります。

$$
\begin{aligned}
\langle \overline{A}iA, j \rangle &= \langle \overline{B}iB, j \rangle \\
\end{aligned}
$$

これを計算することで、5要素の規則が得られます。

&&&fml 5要素の規則 (five-part rules)
$$
\sin a \cos \beta = \cos b \sin c - \sin b \cos c \cos \alpha
$$
&&&

&&&prf
$$
\begin{aligned}
\langle \overline{A}iA, j \rangle
&= \langle e^{-i(\pi-\beta)/2} e^{-ka/2} e^{-i(\pi-\gamma)/2} i e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2}, j \rangle \\
&= \langle e^{-ka/2} i e^{ka/2}, e^{i(\pi-\beta)/2} j e^{-i(\pi-\beta)/2} \rangle \\
&= \langle (\cos a) i - (\sin a) j, -(\cos \beta) j + (\sin \beta) k \rangle \\
&= \sin a \cos \beta \\

\langle \overline{B}iB, j \rangle
&= \langle e^{kc/2} e^{i(\pi-\alpha)/2} e^{kb/2} i e^{-kb/2} e^{-i(\pi-\alpha)/2} e^{-kc/2}, j \rangle \\
&= \langle e^{i(\pi-\alpha)/2} e^{kb/2} i e^{-kb/2} e^{-i(\pi-\alpha)/2}, e^{-kc/2} j e^{kc/2} \rangle \\
&= \langle (\cos b) i - (\sin b)(\cos \alpha) j + (\sin b)(\sin \alpha) k, (\cos c) j + (\sin c) i \rangle \\
&= \cos b \sin c - \sin b \cos \alpha \cos c
\end{aligned}
$$
$$
\begin{aligned}
\langle \overline{A}iA, j \rangle &= \langle \overline{B}iB, j \rangle \\
\sin a \cos \beta &= \cos b \sin c - \sin b \cos \alpha \cos c
\end{aligned}
$$
&&&

## 応用：日の出方程式

5要素の規則において、$\beta = \pi/2$（直角三角形）の場合を考えると、左辺は$0$になります。
$$
0 = \cos b \sin c - \sin b \cos c \cos \alpha
$$
これを整理すると、ネイピアの法則の1つが得られます。

&&&fml ネイピアの法則の1つ
$$
\cos \alpha = \frac{\tan c}{\tan b} \tag{7}
$$
&&&

これを用いて、地球上の緯度$\phi$における日の出時刻を求めることができます。

&&& 用語解説

- **赤緯**: 天の赤道（地球の赤道を天球に投影した大円）から天体がどれだけ北（または南）に離れているかを示す角度です。太陽の赤緯は季節によって変化し、夏至には最大$+23.5^\circ$、冬至には最小$-23.5^\circ$となり、 春分・秋分には$0^\circ$となります。

- **時角**: 天体が子午線（真南）を通過してからどれだけ回転したかを示す角度です。太陽の場合、正午（南中時）に$0^\circ$となり、西（午後）へ行くとプラス、東（午前）へ行くとマイナスになります。

- **白夜**: 太陽の赤緯$\delta$が正のときは、自転軸の北側（北極点）が太陽の方向を向いているため、北極点$(\phi=\pi/2)$から$\delta$の範囲内$(\phi > \pi/2-\delta)$で、太陽は沈みません。

- **極夜**: 太陽の赤緯$\delta$が負のときは、自転軸の北側（北極点）が太陽と反対の方向を向いているため、北極点$(\phi=\pi/2)$から$-\delta$の範囲内$(\phi > \pi/2+\delta)$で、太陽は昇りません。
&&&

### 設定

天球上で以下の直角三角形$\triangle PSN$を考えます。

* 頂点 $P$: 北極星
* 頂点 $S$: 日の出時刻の地平線上の太陽
* 頂点 $N$: 北極の真下の地平線上の点
* 直角 $\angle N = \pi/2$（子午線と地平線の交角）

各辺と角は以下のようになります。

* 斜辺 $b$: 北極星から太陽までの角度。太陽の赤緯を$\delta$とすると、$b = PS = \pi/2 - \delta$。
* 垂直辺 $c$: 地平線から北極星までの角度。観測者の緯度$\phi$に等しい。$c = PN = \phi$。
* 角 $\alpha$: 北極星における角度。太陽の時角（正午からの角度）の絶対値を$\omega$とすると、真北基準の角は$\alpha = \angle P = \pi - \omega$。

&&&rem
日の出の時角は負の値ですが、ここでは三角形の幾何学的性質を扱うため、その大きさ$\omega$を用います。
&&&

### 導出

これらを式 (7) に代入します。

$$
\cos(\pi-\omega) = \frac{\tan \phi}{\tan(\pi/2-\delta)}
$$

$\cos(\pi-\omega) = -\cos \omega,\ \tan(\pi/2-\delta) = 1/\tan \delta$を用いて整理すると、日の出方程式が得られます。

&&&fml 日の出方程式
$$
\cos \omega = - \tan \phi \tan \delta
$$
&&&

左辺を$\omega$とした式も示します。

$$
\omega = \arccos(-\tan \phi \tan \delta)
$$

$\omega$は日の出から正午までに太陽が移動する角度を表すため、日照時間（日の出から日の入りまでの時間）は$2\omega$となります。24時間で一周するため、日照時間は$24 \times \frac{2\omega}{2\pi} = \frac{24}{\pi} \omega$となります。

$\arccos(-x) = \pi - \arccos x$の関係を用いると、日照時間は次のように表されます。（ただし、白夜の地域$\phi > \pi/2 - \delta$を除きます）

$$
\begin{aligned}
\frac{24}{\pi} \omega
&= \frac{24}{\pi}\arccos(-\tan \phi \tan \delta) \\
&= \frac{24}{\pi}(\pi - \arccos(\tan \phi \tan \delta)) \\
&= 24 - \frac{24}{\pi}\arccos(\tan \phi \tan \delta) \\
\end{aligned}
$$

&&&ex
北半球$(\phi > 0)$において、以下の関係が成り立ちます。

*   春分・秋分$(\delta = 0)$:$\cos \omega = 0$より$\omega = \pi/2$（12時間）となり、昼夜の長さが等しくなります。
*   夏$(\delta > 0)$:$\cos \omega < 0$より$\omega > \pi/2$となり、昼が長くなります。
*   冬$(\delta < 0)$:$\cos \omega > 0$より$\omega < \pi/2$となり、昼が短くなります。

夏$(\delta > 0)$に白夜の境界$(\phi = \pi/2 - \delta)$において

$$
\begin{aligned}
\cos \omega
&= -\tan \phi \tan \delta \\
&= -\tan(\pi/2-\delta) \tan \delta \\
&= -\frac{\tan \delta}{\tan \delta} \\
&= -1 \\
\omega &= \arccos(-1) = \pi \\
\frac{24}{\pi} \omega &= 24
\end{aligned}
$$

となるため、日照時間は$24$時間となります。
&&&
