これまでの記事で暗黙のうちに使ってきた「リー群とリー代数の対応」そのものを主題として、2つの道具を整備します。1つは、方向の異なる元で崩れる指数法則のずれを括弧積の言葉で完全に記述する**ベイカー・キャンベル・ハウスドルフの公式**（BCH公式）、もう1つは、たびたび登場した共役作用をリー代数の言葉で捉え直す**随伴作用**です。

シリーズ: [リー群・リー代数の初歩](https://mathlog.info/series/LENUzX64bZ63y462z5I9)

# 前提知識

前回までの記事の内容を前提とします。要点をまとめます。

- リー群は滑らかな構造を持つ群であり、単位元における接空間としてリー代数が得られる。リー代数から指数写像$\exp$を通してリー群の元が得られる
- 本シリーズで扱った$\mathfrak{u}(1)$や$\mathfrak{su}(2)$は通常の積については閉じていないが、括弧積$[x,y]=xy-yx$については閉じている。基底同士の括弧積の係数（構造定数）がリー代数を特徴付ける
- 単位四元数全体はリー群$\operatorname{Sp}(1) \cong \operatorname{SU}(2)$。リー代数$\mathfrak{sp}(1) \cong \mathfrak{su}(2)$は純虚四元数全体で、括弧積は$[i,j]=2k$（巡回的）、純虚四元数を3次元ベクトルと見れば$[p,q] = 2\,(p \times q)$
- 同じ方向の元同士では指数法則$\exp(u\theta_1)\exp(u\theta_2) = \exp(u(\theta_1+\theta_2))$が成り立つが、方向が異なると一般には$\exp(x)\exp(y) \ne \exp(x+y)$[[7shi-lie2]]
- 共役作用$\rho_q(x) = qxq^{-1}$は純虚四元数の3次元空間の回転を与え、$q \mapsto \rho_q$は2対1（二重被覆）。$\mathfrak{so}(3)$の基底$J_x, J_y, J_z$との対応$i/2 \Leftrightarrow J_x,\ j/2 \Leftrightarrow J_y,\ k/2 \Leftrightarrow J_z$で構造定数が一致し、$\mathfrak{su}(2) \cong \mathfrak{so}(3)$[[7shi-lie2]]

&&&rem 本記事の範囲
以下では有限次元の実リー群・実リー代数を扱います。具体的な計算はすべて行列群と四元数群で行い、一般論はそこから見える構造の説明に用います。
&&&

# 指数法則のずれを測る

オイラーの公式$\exp(u\theta) = \cos\theta + u\sin\theta$（$u$は単位純虚四元数）を使って、$\exp(i)\exp(j)$と$\exp(i+j)$を計算し比較します。まず積の方は、

$$
\begin{aligned}
\exp(i)\exp(j)
&= (\cos1 + i\sin1)(\cos1 + j\sin1) \\
&= \cos^2 1 + (i+j)\sin1\cos1 + ij\sin^2 1
\end{aligned}
$$

$ij = k$なので、

$$
\exp(i)\exp(j) = \cos^2 1 + (i+j)\sin1\cos1 + k\sin^2 1
$$

$k$成分$\sin^2 1$が現れます。一方$\exp(i+j)$を計算するには、引数を「単位純虚四元数 × 実数」の形に直す必要があります。$i+j$の絶対値は$\sqrt2$なので、単位純虚四元数$u = (i+j)/\sqrt2$を使って$i+j = u\sqrt2$と書け、

$$
\exp(i+j) = \cos\sqrt2 + \frac{i+j}{\sqrt2}\sin\sqrt2
$$

これは$1, i, j$だけの線形結合で、$k$成分を持ちません。$\exp(i)\exp(j)$には$k$成分$\sin^2 1 \ne 0$があるのに対し$\exp(i+j)$にはないので、

$$
\exp(i)\exp(j) \ne \exp(i+j)
$$

方向の異なる元$i, j$に対して、指数法則$\exp(x)\exp(y) = \exp(x+y)$は崩れます。

とはいえ、$\exp(i)\exp(j)$は単位四元数同士の積なので、それ自身も単位四元数です。つまり何らかの単位純虚四元数$w$と実数$\gamma$によって$\exp(w\gamma)$の形に書けているはずです。そこで問いを立て直します。積を

$$
\exp(x)\exp(y) = \exp(z)
$$

と書いたとき、$z$は$x + y$からどれだけずれるのでしょうか。

## 2次までの展開

$x, y$が小さいとして、$\exp$の級数を2次の項まで展開して掛け合わせます。3次以上の項は$\cdots$にまとめます。

$$
\begin{aligned}
\exp(x)\exp(y)
&= \left(1 + x + \frac{x^2}{2} + \cdots\right)\left(1 + y + \frac{y^2}{2} + \cdots\right) \\
&= 1 + x + y + \frac{x^2}{2} + xy + \frac{y^2}{2} + \cdots
\end{aligned}
$$

これを$\exp(z)$と等しく置きます。1次の項の比較から$z$の1次部分は$x + y$なので、2次の補正項を$w$として$z = x + y + w + \cdots$と書けば、

$$
\exp(z) = 1 + z + \frac{z^2}{2} + \cdots
= 1 + x + y + w + \frac{(x+y)^2}{2} + \cdots
$$

となります。積が非可換であることに注意して$(x+y)^2 = x^2 + xy + yx + y^2$と展開し、2次の項を比較すれば

$$
w = xy - \frac{xy + yx}{2} = \frac{xy - yx}{2} = \frac{1}{2}[x, y]
$$

が得られます。つまり2次までの範囲で次が成り立ちます。

$$
\exp(x)\exp(y) = \exp\left(x + y + \frac{1}{2}[x,y] + \cdots\right)
$$

指数法則のずれは2次から始まり、その最初の項は括弧積そのものです。$x, y$が可換なら$[x,y] = 0$となって（後で見るように高次の補正もすべて消えて）指数法則が復活します。$\operatorname{U}(1)$で指数法則が成り立ったのはこのためです。以前は括弧積を「積の交換によるずれを測る量」として導入しましたが、同じ量が「指数法則のずれ」の主要項としても現れました。

## ベイカー・キャンベル・ハウスドルフの公式

3次まで進めると次のようになります（計算は省略します）。

&&&fml ベイカー・キャンベル・ハウスドルフの公式（BCH公式）
$$
\exp(x)\exp(y) = \exp\left(x + y + \frac{1}{2}[x,y]
+ \frac{1}{12}[x,[x,y]] + \frac{1}{12}[y,[y,x]] + \cdots\right)
$$
&&&

3次の項には、括弧積の結果とさらに括弧積を取る**入れ子の括弧積**が現れました。重要なのは係数の具体的な値ではなく、高次の項はどこまで行っても括弧積の入れ子だけで書けるという事実です。

すべての項を入れ子の括弧積で明示的に書き下す公式がディンキンにより与えられています。

&&&fml ディンキンの公式
$$
z = \sum_{n=1}^{\infty} \frac{(-1)^{n-1}}{n} \sum_{\substack{r_i+s_i>0 \\ r_i,\,s_i \in \mathbb{Z}_{\ge 0} \\ 1\le i\le n}} \frac{[x^{r_1}y^{s_1}x^{r_2}y^{s_2}\cdots x^{r_n}y^{s_n}]}{\left(\sum_{i=1}^n(r_i+s_i)\right) r_1!\,s_1!\cdots r_n!\,s_n!}
$$
&&&

$[x^{r_1}y^{s_1}\cdots]$は「$x$との括弧積を$r_1$回、続けて$y$との括弧積を$s_1$回、…」という手順で作る入れ子の括弧積を表す記法で、括弧の向きは

$$
[a_1 a_2 \cdots a_m] = [a_1,[a_2,\ldots,[a_{m-1},a_m]\ldots]]
$$

と定めます（$[x^1] = x$、$[x^2y^1] = [x,[x,y]]$など）。各項の次数は$\sum_i(r_i+s_i)$で決まります。同じ次数の項を集めて整理すると、先ほどのBCH公式の係数が得られます。4次の項をまとめると$-\frac{1}{24}[y,[x,[x,y]]]$の1つだけになりますが、5次以降は入れ子の組み合わせ方が増え、項数が急速に増えます。

&&&rem
実際の計算で高次項まで使うことは稀で、この公式の価値はもっぱら「括弧積だけで書ける」という構造的な事実にあります。
&&&

BCH公式やディンキンの公式で$z$を与えるこの無限級数（**BCH級数**）は、$x, y$が十分小さい範囲でのみ収束が保証される、局所的なものです。リー代数が単位元近傍の局所構造を記述する道具であることと対応しています。

&&&rem
冒頭の$\exp(i)\exp(j)$のように$x, y$が小さくない場合でも、$\mathfrak{sp}(1) \cong \mathfrak{su}(2)$に限れば級数に頼らない厳密な計算方法が後で登場します。
&&&

## 括弧積が積を決める

BCH公式の意味を考えます。右辺に現れるのは$x, y$の線形結合と括弧積だけで、行列の積$xy$そのものは現れません。つまり、リー代数を「ベクトル空間と括弧積の組」としてだけ知っていれば（個々の元がどんな行列か、そもそも行列かどうかすら知らなくても）単位元の近くでの群の積が計算できてしまいます。

以前、$\operatorname{SU}(2)$と$\operatorname{SO}(3)$はリー群としては2対1で異なるのにリー代数としては同型であることから、「リー代数は局所的な構造しか見ない」と述べました。BCH公式はその裏返しです。**局所的な構造なら、括弧積だけで完全に決まります。**$\operatorname{SU}(2)$と$\operatorname{SO}(3)$の違いは、単位元の近くには一切現れない大域的な形（二重被覆）だけの違いだったことが、ここで裏付けられます。

&&&rem リー対応
「リー代数からリー群を局所的に復元できる」という事実は、リー理論の基本定理の1つです。この対応により、連結なリー群の局所的な分類はリー代数の分類に帰着します。ただし大域的な群を区別するには、被覆や中心による商といった情報がさらに必要です。
&&&

## 群の交換子

BCH公式の応用として、「行って、行って、戻って、戻る」という往復の組み合わせを計算します。

$$
\exp(x)\exp(y)\exp(-x)\exp(-y)
$$

前半と後半をそれぞれBCH公式でまとめます。$[-x,-y] = [x,y]$に注意すれば、

$$
\begin{aligned}
\exp(x)\exp(y) &= \exp\Bigl(\underbrace{x + y + \tfrac{1}{2}[x,y] + \cdots}_{A}\Bigr) \\
\exp(-x)\exp(-y) &= \exp\Bigl(\underbrace{-x - y + \tfrac{1}{2}[x,y] + \cdots}_{B}\Bigr)
\end{aligned}
$$

となります。$A, B$にもう一度BCH公式を使うと、1次と2次の項は

$$
A + B = [x,y] + \cdots, \quad
\frac{1}{2}[A, B] = \frac{1}{2}[x+y,\ -(x+y)] + \cdots = 0 + \cdots
$$

なので（$\cdots$は3次以上）、次が得られます。

&&&fml 群の交換子
$$
\exp(x)\exp(y)\exp(-x)\exp(-y) = \exp([x,y] + \cdots)
$$
&&&

左辺は群における「非可換性の測定装置」です。単位元から出発して右から元を掛け足していく経路として読むと、$x$方向に進み、$y$方向に進み、来た道を逆順ではなく同順で戻ることに当たります。可換なら出発点に戻りますが、非可換なら戻りきれません。その閉じ残りが、2次の精度でちょうど括弧積$[x,y]$の方向に現れます。

$\mathfrak{so}(3)$で具体的に確かめます。$x$軸・$y$軸周りの角度$\varepsilon$の回転行列は

$$
R_x(\varepsilon) = \begin{pmatrix}1&0&0\\0&\cos\varepsilon&-\sin\varepsilon\\0&\sin\varepsilon&\cos\varepsilon\end{pmatrix}, \quad
R_y(\varepsilon) = \begin{pmatrix}\cos\varepsilon&0&\sin\varepsilon\\0&1&0\\-\sin\varepsilon&0&\cos\varepsilon\end{pmatrix}
$$

です。$[J_x, J_y] = J_z$なので、上の一般論から$R_x(\varepsilon)R_y(\varepsilon)R_x(-\varepsilon)R_y(-\varepsilon) = I + \varepsilon^2 J_z + O(\varepsilon^3)$となるはずです。実際に4つの行列を掛け合わせて$\varepsilon^2$の項まで展開すると

$$
R_x(\varepsilon)R_y(\varepsilon)R_x(-\varepsilon)R_y(-\varepsilon) =
\begin{pmatrix}1&-\varepsilon^2&0\\ \varepsilon^2&1&0\\0&0&1\end{pmatrix} + O(\varepsilon^3)
$$

となり、$(1,2)$・$(2,1)$成分に$\mp\varepsilon^2$が現れます。これは2次まで見れば$z$軸周りの角度$\varepsilon^2$の回転$\exp(\varepsilon^2 J_z) = I + \varepsilon^2 J_z + \cdots$と一致します（有限の$\varepsilon$で厳密にこの回転になるわけではなく、3次以上のずれが残ります）。なお、この行列積はベクトルに対しては右から順に作用します。$x$軸周りの微小角$\varepsilon$の回転と$y$軸周りの微小角$\varepsilon$の回転を往復させると、2次の精度で$z$軸周りの角度$\varepsilon^2$の回転が残ります。$x$軸と$y$軸の回転だけを組み合わせて$z$軸の回転が作り出せるのは、この仕組みによります。

&&&rem 縦列駐車
「2つの操作を往復させると、どちらの操作でも直接はできない第3の動きが残る」という構図は、縦列駐車に喩えられます。「ハンドルを切る」と「前後に進む」の2操作しかない車が、操作の往復や、その往復をさらに組み合わせることで、直接には動けない方向への移動を作り出します。制御理論では、交換子で生まれる方向まで含めて車や機械の到達できる範囲が決まります。
&&&

## 四元数による厳密な合成

一般のリー代数ではBCH公式は無限級数ですが、$\mathfrak{sp}(1) \cong \mathfrak{su}(2)$では級数に頼らない厳密な合成則が得られます。数体系（四元数）の中では、積を直接計算できるからです。

単位純虚四元数$u, v$と実数$\alpha, \beta$に対して、オイラーの公式と純虚四元数の積の分解$uv = -u \cdot v + u \times v$を使って積を展開します。

$$
\begin{aligned}
&\exp(u\alpha)\exp(v\beta) \\
&= (\cos\alpha + u\sin\alpha)(\cos\beta + v\sin\beta) \\
&= \underbrace{\cos\alpha\cos\beta - (u \cdot v)\sin\alpha\sin\beta}_{\text{実部}}
+ \underbrace{u\sin\alpha\cos\beta + v\cos\alpha\sin\beta + (u \times v)\sin\alpha\sin\beta}_{\text{虚部}}
\end{aligned}
$$

&&&ex $\exp(i)\exp(j)$
$$
u = i,\ v = j,\ \alpha = \beta = 1
$$
に当てはめると、$u \cdot v = 0,\ u \times v = k$より実部は$\cos^2 1$、虚部は$(i+j)\sin1\cos1 + k\sin^2 1$となり、直接計算した結果と一致します。
&&&

積は再び単位四元数なので、実部を$a$、虚部を$V$として$a + V = \cos\gamma + w\sin\gamma$の形で書けます。ここから合成後の角度$\gamma$と方向$w$を取り出すには、次のように選べます。

&&&fml 合成後の角度と方向
$V \ne 0$のとき
$$
\gamma = \operatorname{atan2}(|V|, a) \in (0, \pi), \qquad w = \frac{V}{|V|}
$$
&&&

実部だけでは$\gamma$が一意に定まらないため、$|V|$と$a$の両方を使って枝を選びます。$V = 0$の場合は積が$\pm1$となり、方向$w$は定まりません（任意に取れます）。

このように、角度と方向が閉じた式で求まります。角度の式を取り出しておきます。

&&&fml 四元数の合成角の公式
$$
\cos\gamma = \cos\alpha\cos\beta - (u \cdot v)\sin\alpha\sin\beta
$$
&&&

これは任意の実数$\alpha, \beta$について成り立つ積の公式ですが、球面三角形の辺長としての解釈を与えると球面余弦定理になります。単位四元数全体は3次元球面$S^3$なので、$1,\ \exp(u\alpha),\ \exp(u\alpha)\exp(v\beta)$の3点を頂点とする球面三角形が描けます。単位四元数を掛ける操作は$S^3$上の等長変換なので、$\exp(u\alpha)$から積までの距離は$1$から$\exp(v\beta)$までの距離に等しく、3辺の長さはそれぞれ$\alpha, \beta, \gamma$です。

頂点$\exp(u\alpha)$における内角は、両隣の頂点への方向を左から$\exp(u\alpha)^{-1}$を掛けて$1$に揃えると、$\exp(-u\alpha)$の方向$-u$と$\exp(v\beta)$の方向$v$のなす角です。$u, v$のなす角を$\phi$（$u \cdot v = \cos\phi$）とすれば、内角は$\pi - \phi$となります。

球面余弦定理は、辺の長さ$a, b, c$と、辺$b, c$に挟まれた内角$A$について

$$
\cos a = \cos b \cos c + \sin b \sin c \cos A
$$

という形の公式です（記号の重複を避けるため、内角を$A$と書いています）[[7shi-strig]]。上の対応$a = \gamma,\ b = \alpha,\ c = \beta,\ A = \pi - \phi$を代入すると、$\cos A = -\cos\phi = -(u \cdot v)$より、合成角の公式に一致します。

&&&rem
四元数の指数関数の積からは、球面余弦定理の他にも正弦定理や5要素の規則が導出できます。[[7shi-strig]]
&&&

BCH公式との整合も確認できます。BCH公式が与えるのは積そのものの虚部ではなく指数の引数$z = w\gamma$ですが、上の表示から

$$
z = \frac{\gamma}{\sin\gamma}\,V
$$

であり、単位元近傍では$\gamma/\sin\gamma = 1 + O(\gamma^2)$なので、$z$と虚部$V$の違いは3次以上に現れます。したがって2次までは虚部をそのまま比較できます。虚部
$$
u\sin\alpha\cos\beta + v\cos\alpha\sin\beta + (u \times v)\sin\alpha\sin\beta
$$
で$\alpha, \beta$が小さいとすると、$\sin\alpha\cos\beta \approx \alpha,\ \cos\alpha\sin\beta \approx \beta,\ \sin\alpha\sin\beta \approx \alpha\beta$より、虚部は
$$
u\alpha + v\beta + (u \times v)\alpha\beta + \cdots
$$
と展開されます。3項目は、括弧積が外積の2倍であることから

$$
\frac{1}{2}[u\alpha, v\beta] = \frac{1}{2}\alpha\beta \cdot 2\,(u \times v) = \alpha\beta\,(u \times v)
$$

となり、BCH公式の2次項と一致します。BCH級数の先頭部分が、四元数では外積の項として顔を出しているわけです。

BCH級数には収束範囲の制約がありましたが、$\exp(u\alpha)\exp(v\beta) = \cos\gamma + w\sin\gamma$という合成則そのものは任意の$\alpha, \beta$で成り立ちます。四元数では、実部と3次元の虚部だけで積が記述でき、そこから合成後の方向と角度、さらに指数の引数までが簡潔に取り出せます。一般のリー群では、同じように扱いやすい閉じた式が得られるとは限りません。

# ヤコビ恒等式

BCH公式の3次の項には$[x,[x,y]]$のような入れ子の括弧積が現れました。入れ子の括弧積には、すべてのリー代数で成り立つ基本的な恒等式があります。

&&&fml ヤコビ恒等式
$$
[x,[y,z]] + [y,[z,x]] + [z,[x,y]] = 0
$$
&&&

&&&prf
括弧積を定義通り展開すると
$$
\begin{aligned}
[x,[y,z]] &= x(yz-zy) - (yz-zy)x = xyz - xzy - yzx + zyx \\
[y,[z,x]] &= yzx - yxz - zxy + xzy \\
[z,[x,y]] &= zxy - zyx - xyz + yxz
\end{aligned}
$$
となります。3つの和を項ごとに見ると、すべての項が符号違いの相手と対になって打ち消し合います。
&&&

この証明は、$xyz$のような3文字の積が括弧の付け方によらず1通りに定まる、つまり積の**結合法則**に依存しています。結合法則が崩れる代数では、この証明は足場ごと崩れます。

&&&rem
証明で使ったのは結合法則だけで、積$xy=yx$のような交換法則は仮定していません。ヤコビ恒等式は非可換な結合代数でも成り立ちます。
&&&

## ライプニッツ則としての読み替え

ヤコビ恒等式は、反対称性$[z,x] = -[x,z]$を使って移項すると、次の形に書き直せます。

$$
[x,[y,z]] = [[x,y],z] + [y,[x,z]]
$$

「$x$との括弧積を取る」という操作を$D_x$と書けば、この式は

$$
D_x[y,z] = [D_x y, z] + [y, D_x z]
$$

となります。これは積の微分におけるライプニッツ則（$(fg)' = f'g + fg'$）と同じ形です。つまりヤコビ恒等式は、「$x$との括弧積を取る操作が、括弧積に対して微分のように振る舞う」ことを主張しています。

この「$x$との括弧積を取る操作」こそ、次節の主役である随伴作用のリー代数版です。

&&&rem 抽象的なリー代数
ここまでリー代数は「行列や四元数の集合」として具体的に扱ってきましたが、公理的には「双線形・反対称でヤコビ恒等式を満たす括弧積を持つベクトル空間」として定義されます。行列の交換子$[x,y] = xy - yx$は、この公理を満たす代表例です。BCH公式が「括弧積だけで積を復元できる」ことを保証しているため、この抽象化によって失われる局所的な情報はありません。
&&&

# 随伴作用

単位四元数$q$の共役作用$\rho_q(x) = qxq^{-1}$は、純虚四元数の3次元空間の回転を与えました。純虚四元数全体はリー代数$\mathfrak{sp}(1)$そのものであり、これは**リー群$G$の元$g$が、自分自身のリー代数$\mathfrak{g}$に共役で作用する**という構図の例になっています。この作用を**随伴作用**と呼び、$\operatorname{Ad}_g$と書きます。

一般のリー群では$g \in G$と$x \in \mathfrak{g}$を直接掛けることはできないため、定義は共役写像の微分として与えます。

&&&def リー群の随伴作用
$$
C_g(h) = ghg^{-1} \quad (g, h \in G), \qquad
\operatorname{Ad}_g = (dC_g)_e
$$
&&&

$C_g$は単位元を単位元に写すので、その微分$(dC_g)_e$は単位元での接空間、すなわちリー代数$\mathfrak{g}$上の線形写像になります。行列群や四元数群では、この微分は次のように書けます。

$$
\operatorname{Ad}_g(x) = gxg^{-1} \quad (g \in G,\ x \in \mathfrak{g})
$$

結果がリー代数に留まることは、この形からも確かめられます。$x \in \mathfrak{g}$に対して$\exp(tx)$は単位元を通る$G$内の曲線で、$x$はその$t = 0$における速度ベクトルです。共役で写した$g\exp(tx)g^{-1}$も単位元を通る$G$内の曲線であり、その速度ベクトルは

$$
\left.\frac{d}{dt}\,g\exp(tx)g^{-1}\right|_{t=0} = gxg^{-1}
$$

となります。単位元における速度ベクトルは接空間の元、すなわちリー代数の元なので、$gxg^{-1} \in \mathfrak{g}$です。以前は純虚性を計算して個別に確認しましたが、実は「リー代数は単位元における接空間」という一般論だけから従うのです。

&&&rem 相似変換との関係
$\operatorname{Ad}_g(x) = gxg^{-1}$は、線形代数の相似変換$A \mapsto PAP^{-1}$の共役元$P$を$G$の元に限定したものです。相似変換が積を保つ性質がそのまま上式の括弧積保存の根拠になっており、$\mathfrak{g}$からはみ出さないのは$\mathfrak{g}$が単位元での接空間だからです。
&&&

随伴作用は2つの構造を保ちます。第一に、群の積と整合します。

$$
\operatorname{Ad}_{gh}(x) = (gh)x(gh)^{-1} = g(hxh^{-1})g^{-1} = (\operatorname{Ad}_g \circ \operatorname{Ad}_h)(x)
$$

第二に、括弧積を保ちます。積の間に$g^{-1}g = 1$を挟むだけで確かめられます。

$$
\operatorname{Ad}_g[x,y] = g(xy-yx)g^{-1} = (gxg^{-1})(gyg^{-1}) - (gyg^{-1})(gxg^{-1}) = [\operatorname{Ad}_g x, \operatorname{Ad}_g y]
$$

&&&rem 括弧の省略
$\operatorname{Ad}_g(x)$と$\operatorname{Ad}_g x$は同じ意味です。括弧を省略することもあります。
&&&

つまり$\operatorname{Ad}_g$は、リー代数の構造（線形性と括弧積）を保つ可逆な変換、すなわちリー代数の**自己同型**です。

もう1つ、群のレベルでの関係を挙げておきます。$\exp$の級数の各項$y^n$に対して$g y^n g^{-1} = (gyg^{-1})^n$が成り立つので、

$$
C_g(\exp y) = g\exp(y)g^{-1} = \exp(gyg^{-1}) = \exp(\operatorname{Ad}_g y)
$$

となります。群の元の共役は、指数写像の引数への随伴作用に翻訳されます。一般の3つの積をBCH公式で展開すると無限級数になりますが、共役という特別な組み合わせでは厳密に閉じた式になります。

## 微分としての$\operatorname{ad}$

リー群の議論をリー代数に落とすときの常套手段は、単位元近傍での微分でした。随伴作用にもこれを適用します。$g$を単位元を通る曲線$\exp(tz)$として、$\operatorname{Ad}_{\exp(tz)}(x)$を$t = 0$で微分します。積の微分法則より

$$
\begin{aligned}
\left.\frac{d}{dt}\operatorname{Ad}_{\exp(tz)}(x)\right|_{t=0}
&= \left.\frac{d}{dt}\exp(tz)\,x\exp(-tz)\right|_{t=0} \\
&= \left.z\exp(tz)\,x\exp(-tz) - \exp(tz)\,x\exp(-tz)\,z\right|_{t=0} \\
&= zx - xz \\
&= [z, x]
\end{aligned}
$$

となり、括弧積が現れました。この線形写像を$\operatorname{ad}_z$と書きます。

&&&def リー代数の随伴作用
$$
\operatorname{ad}_z(x) = [z, x]
$$
&&&

大文字の$\operatorname{Ad}$がリー群の随伴作用、小文字の$\operatorname{ad}$がそのリー代数版です。**共役作用の無限小版は、括弧積そのものだった**わけです。前節の$D_x$も、この記法では$\operatorname{ad}_x$と書けます。

$$
\operatorname{ad}_x[y,z] = [\operatorname{ad}_x y, z] + [y, \operatorname{ad}_x z]
$$

このライプニッツ則は、「$\operatorname{Ad}_g$が括弧積を保つ」という性質の無限小版として導くこともできます。$g = \exp(tx)$として$\operatorname{Ad}_g[y,z] = [\operatorname{Ad}_g y, \operatorname{Ad}_g z]$の両辺を$t = 0$で微分すると、右辺は積の微分法則によって2項に分かれ、ちょうど上の式になります。ヤコビ恒等式は、共役作用が括弧積を保つことの、リー代数への写し絵でもあるのです。

## $\operatorname{Ad}$と$\operatorname{ad}$をつなぐ指数写像

リー群とリー代数が指数写像で結ばれるように、$\operatorname{Ad}$と$\operatorname{ad}$も指数写像で結ばれます。
$$
F(t) = \operatorname{Ad}_{\exp(tz)}(x) = \exp(tz)\,x\exp(-tz)
$$
と置いて$t$で微分すると、
$$
\begin{aligned}
F'(t)
&= z\exp(tz)\,x\exp(-tz) - \exp(tz)\,x\exp(-tz)\,z \\
&= z\,F(t) - F(t)\,z \\
&= [z, F(t)] \\
&= \operatorname{ad}_z F(t)
\end{aligned}
$$
となります。これはスカラーの微分方程式$y'(t) = ay(t)$の$a$を線形写像$\operatorname{ad}_z$に置き換えた形になっています。スカラーの場合、初期値$y(0) = c$のもとで解は$y(t) = \exp(at)c$と指数関数で書けました。$F$についても同じ形の解が期待できます。

$$
F(t) = \exp(t\operatorname{ad}_z)\,x
$$

&&&prf 解の形が正しいことの確認
$\operatorname{ad}_z$はリー代数上の線形写像（リー代数からリー代数自身への線形変換）なので、行列と同じように和・スカラー倍・合成（べき乗）が定義できます。そこで$\exp(t\operatorname{ad}_z)$を、行列の指数関数と同じ級数
$$
\exp(t\operatorname{ad}_z) = \operatorname{id}_{\mathfrak{g}} + t\operatorname{ad}_z + \frac{t^2}{2}\operatorname{ad}_z^2 + \cdots
$$
で定義します（$\operatorname{id}_{\mathfrak{g}}$は$\mathfrak{g}$上の恒等写像、$\operatorname{ad}_z^2$は$\operatorname{ad}_z$を2回合成する意味）。この級数を項別に$t$で微分すると
$$
\frac{d}{dt}\exp(t\operatorname{ad}_z) = \operatorname{ad}_z + t\operatorname{ad}_z^2 + \cdots = \operatorname{ad}_z\exp(t\operatorname{ad}_z)
$$
となるので、$F(t) = \exp(t\operatorname{ad}_z)\,x$は確かに$F'(t) = \operatorname{ad}_z F(t)$と$F(0) = x$を満たします。線形常微分方程式の解の一意性により、これが$F(t)$に一致します。
&&&

$t = 1$と置けば、次の関係が得られます。

&&&fml $\operatorname{Ad}$と$\operatorname{ad}$の関係
$$
\operatorname{Ad}_{\exp z} = \exp(\operatorname{ad}_z)
$$
&&&

級数で書き下せば次のようになります。

&&&lem アダマールの補題
$$
\exp(z)\,x\exp(-z) = x + [z, x] + \frac{1}{2!}[z, [z, x]] + \frac{1}{3!}[z, [z, [z, x]]] + \cdots
$$
&&&

左辺は群の元による共役、右辺は括弧積の入れ子だけです。BCH公式と同じく「積の情報が括弧積で書ける」形ですが、こちらは任意の$z$で収束する厳密な級数です。

まとめると、リー代数の元$z$を$\exp$で持ち上げるとリー群の元$\exp z$になり、$z$の括弧積$\operatorname{ad}_z$を$\exp$で持ち上げると$\exp z$の共役作用$\operatorname{Ad}_{\exp z}$になります。元と作用のどちらのレベルでも、同じ指数写像がリー代数とリー群の橋になっています。

## $\mathfrak{su}(2)$での具体形

一般論はここまでにして、$\mathfrak{sp}(1) \cong \mathfrak{su}(2)$で$\operatorname{ad}$の具体的な形を求めます。基底$i, j, k$に対する$\operatorname{ad}_{i/2}$の作用を、括弧積$[i,j] = 2k$などから計算します。

$$
\operatorname{ad}_{i/2}(i) = 0, \quad
\operatorname{ad}_{i/2}(j) = \frac{1}{2}[i,j] = k, \quad
\operatorname{ad}_{i/2}(k) = \frac{1}{2}[i,k] = -j
$$

純虚四元数$xi + yj + zk$を3次元ベクトル$(x,y,z)$と同一視して、この線形写像を行列で書きます。$j \mapsto k,\ k \mapsto -j$は$y$成分を$z$成分へ、$z$成分を$-y$成分へ送るので、

$$
\operatorname{ad}_{i/2} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & -1 \\ 0 & 1 & 0 \end{pmatrix} = J_x
$$

$\mathfrak{so}(3)$の基底$J_x$そのものです。同様に$\operatorname{ad}_{j/2} = J_y$、$\operatorname{ad}_{k/2} = J_z$となります。

以前、対応$i/2 \Leftrightarrow J_x,\ j/2 \Leftrightarrow J_y,\ k/2 \Leftrightarrow J_z$で構造定数が一致することを計算で確認し、同型$\mathfrak{su}(2) \cong \mathfrak{so}(3)$を得ました。その対応の背後にある構造が、いま明らかになりました。**この同型対応は、写像$x \mapsto \operatorname{ad}_x$そのものです。**

さらに$\operatorname{Ad}_{\exp z} = \exp(\operatorname{ad}_z)$を使えば、先の中心的な結果（共役作用が回転を与えること）が1行で再導出されます。$q = \exp(\theta i/2)$に対して

$$
\rho_q = \operatorname{Ad}_q = \operatorname{Ad}_{\exp(\theta i/2)} = \exp(\theta\operatorname{ad}_{i/2}) = \exp(\theta J_x)
$$

となり、右辺は$x$軸周りの角度$\theta$の回転行列です。以前は$\rho_q(j) = j\cos\theta + k\sin\theta$のような成分計算で確認しましたが、その背後には$\operatorname{Ad}$と$\operatorname{ad}$をつなぐ一般的な仕組みが働いていたのです。謎めいて見えた半角も、$\exp$の引数を$i \cdot (\theta/2)$ではなく$(i/2) \cdot \theta$と括り直せば、基底の対応$\operatorname{ad}_{i/2} = J_x$の帰結として整理されます。

## 二重被覆の構造的な理由

$g \mapsto \operatorname{Ad}_g$は群の準同型でした（$\operatorname{Ad}_{gh} = \operatorname{Ad}_g \circ \operatorname{Ad}_h$）。以前見た二重被覆$\operatorname{SU}(2) \to \operatorname{SO}(3)$は、まさにこの準同型$q \mapsto \operatorname{Ad}_q$です。準同型の対応が「何対1」になるかは、**核**（恒等変換に写る元の集合）で決まります。

$\operatorname{Ad}_q$が恒等変換になるのは、$qxq^{-1} = x$、つまり$q$がすべての純虚四元数と可換になるときです。実数はすべての四元数と可換なので、このとき$q$は「実部＋純虚部」であるすべての四元数と可換になります。すべての四元数と可換な四元数は実数に限られ、単位四元数の中で実数は$q = \pm 1$の2つだけです。

$$
\ker \operatorname{Ad} = \{\pm 1\}
$$

核が2元なので、対応は2対1です。以前は「$q$と$-q$が同じ回転を与える」ことを式変形で見ましたが、随伴作用の言葉では「核が$\{\pm 1\}$である」という構造的な事実に整理されます。

&&&rem 中心
すべての元と可換な元の集合を群の**中心**と呼びます。$\operatorname{SU}(2)$の中心が$\{\pm 1\}$であることが、二重被覆の源です。連結なリー群では、随伴作用の核は中心に一致することが知られています。
&&&

## 随伴表現

$g \mapsto \operatorname{Ad}_g$は、リー群の元を「リー代数というベクトル空間上の可逆な線形変換」、つまり行列に対応させる準同型です。このように群の元を行列（線形変換）で実現する対応を**表現**と呼び、特にこの表現を**随伴表現**と呼びます。リー代数側の$x \mapsto \operatorname{ad}_x$も同様に、リー代数の随伴表現と呼びます。

表現として見ると、その次元（作用する空間の次元）はリー代数自身の次元です。$\mathfrak{su}(2)$の随伴表現は3次元で、リー群側の像はちょうど$\operatorname{SO}(3)$でした。

$\operatorname{ad}$を行列として書き下すと、その成分は構造定数です。基底$e_a$に対して構造定数を
$$
[e_a, e_b] = \sum_c C_{abc}\,e_c
$$
で定めれば、
$$
\operatorname{ad}_{e_a}(e_b) = \sum_c C_{abc}\,e_c
$$
なので、$\operatorname{ad}_{e_a}$を行列と見たときの$(c, b)$成分が$C_{abc}$です。以前の注釈で「リー代数の性質は構造定数によって完全に特徴付けられる」と述べましたが、随伴表現は、構造定数を並べて作った行列でリー代数自身を表現し直したものと言えます。$\operatorname{ad}_{i/2} = J_x$の行列成分が$\mathfrak{su}(2)$の構造定数だったのが、その実例です。基底$i, j, k$では構造定数は$2\epsilon_{abc}$ですが、基底を$e_1 = i/2,\ e_2 = j/2,\ e_3 = k/2$に取ると構造定数はレヴィ＝チヴィタ記号$\epsilon_{abc}$そのものとなり、係数$2$が消えます。

ヤコビ恒等式にも、随伴表現の言葉での3つ目の顔があります。

&&&fml 随伴表現の準同型性
$$
\operatorname{ad}_{[x,y]} = [\operatorname{ad}_x, \operatorname{ad}_y] = \operatorname{ad}_x\operatorname{ad}_y - \operatorname{ad}_y\operatorname{ad}_x
$$
&&&

左辺と右辺は、どちらも「リー代数の元を受け取って別の元を返す」線形変換なので、等しいかどうかは任意の元$z$に作用させて比べられます。左辺の$\operatorname{ad}_{[x,y]}$は「$[x,y]$との括弧積を取る」操作なので

$$
\operatorname{ad}_{[x,y]}(z) = [[x,y],z]
$$

です。右辺の$\operatorname{ad}_x\operatorname{ad}_y$は合成、つまり「先に$\operatorname{ad}_y$、次に$\operatorname{ad}_x$を作用させる」操作なので

$$
(\operatorname{ad}_x\operatorname{ad}_y - \operatorname{ad}_y\operatorname{ad}_x)(z) = [x,[y,z]] - [y,[x,z]]
$$

です。したがって両辺が等しいという主張は

$$
[[x,y],z] = [x,[y,z]] - [y,[x,z]]
$$

に帰着しますが、これはライプニッツ則$[x,[y,z]] = [[x,y],z] + [y,[x,z]]$の右辺第2項を左辺に移項しただけの式です。つまり、ヤコビ恒等式の書き換えにほかなりません。

随伴表現の準同型性が語っているのは、写像$x \mapsto \operatorname{ad}_x$の性質です。左辺は「リー代数側で括弧積$[x,y]$を作ってから、$\operatorname{ad}$で行列に写す」、右辺は「先に$x, y$をそれぞれ行列$\operatorname{ad}_x, \operatorname{ad}_y$に写してから、行列の交換子を取る」という手順で、どちらの順路でも同じ行列に到着します。つまり$\operatorname{ad}$は、リー代数の括弧積を行列の交換子に写す、構造を保つ写像（準同型）です。括弧積の計算を随伴表現に写して行ってよいことの保証になります。ただし、一般には忠実とは限りません。

$$
\ker \operatorname{ad} = \{x \mid [x,y] = 0 \text{ for all } y\} = Z(\mathfrak{g})
$$

であり、この核（リー代数の**中心**）に入る元の情報は失われます。可換なリー代数では$\operatorname{ad}$はすべて$0$になり、何も見えません。$\mathfrak{su}(2)$では中心が$0$なので、随伴表現に写しても情報は失われません。

まとめると、ヤコビ恒等式には3つの読み方があることになります。

1. 入れ子の括弧積の巡回和が消える恒等式（行列や四元数では結合法則の帰結、抽象的なリー代数では公理）
2. $\operatorname{ad}_x$がライプニッツ則を満たす（括弧積の「微分」として振る舞う）
3. $x \mapsto \operatorname{ad}_x$が括弧積を交換子に写す（随伴「表現」の名の由来）

&&&rem ルートへの布石
随伴表現は、リー代数の分類理論の主役でもあります。複素半単純リー代数では、カルタン部分代数$\mathfrak{h}$の元$h$による作用$\operatorname{ad}_h$を同時に調べます。共通固有ベクトル$v$に対して$[h, v] = \alpha(h)\,v$と書いたとき、この非零の線形汎関数$\alpha$が**ルート**と呼ばれる分類の基本データです。
&&&

# 結合法則という土台

本記事で行った行列・四元数による具体計算を見渡すと、すべてが積の**結合法則**の上に建っていることに気づきます。

- BCH公式の導出では、$\exp$の級数同士の積を、因子の順序を保って展開・整理しました。3文字以上の積を括弧の付け方を気にせず扱えるのは、結合法則あってこそです
- ヤコビ恒等式の証明は、$xyz$のような積が1通りに定まることに全面的に依存していました
- $\operatorname{Ad}_{gh} = \operatorname{Ad}_g \circ \operatorname{Ad}_h$や$g\exp(y)g^{-1} = \exp(\operatorname{Ad}_g y)$も、積の組み替えを自由に行う計算でした

そもそも群の定義自体が結合法則を含んでいるので、これは当然のことです。しかし当然だからこそ、普段は意識に上りません。なお、抽象的なリー代数では結合的な積そのものを仮定せず、ヤコビ恒等式を公理として課します。結合法則が支えているのは、あくまでこうした具体的な代数での計算です。

# まとめ

方向の異なる元では崩れる指数法則のずれは、BCH公式

$$
\exp(x)\exp(y) = \exp\left(x + y + \frac{1}{2}[x,y] + \cdots\right)
$$

によって括弧積の言葉で完全に記述されます。高次項はすべて括弧積の入れ子で書けるため、リー代数（ベクトル空間と括弧積の組）を知れば単位元近傍での群の積が復元できます。リー代数が見るのは局所構造だけですが（二重被覆の例で見た通り）、逆に局所構造は括弧積で完全に決まります。応用として、群の往復
$$
\exp(x)\exp(y)\exp(-x)\exp(-y) = \exp([x,y] + \cdots)
$$
は括弧積を群のレベルで取り出します。また$\mathfrak{sp}(1) \cong \mathfrak{su}(2)$では、四元数の積がBCH級数を経ない厳密な合成則を与え、$S^3$上の球面三角形を通して球面余弦定理につながります。

入れ子の括弧積を律するのがヤコビ恒等式$[x,[y,z]] + [y,[z,x]] + [z,[x,y]] = 0$です。証明は結合法則を使った展開だけで済み、ライプニッツ則の形にも、「$\operatorname{ad}$が括弧積を交換子に写す」という主張にも読み替えられます。

たびたび登場した共役作用は、リー群$G$が自分のリー代数$\mathfrak{g}$に作用する随伴作用$\operatorname{Ad}_g(x) = gxg^{-1}$として統一されます。その無限小版が$\operatorname{ad}_z = [z,\ \cdot\ ]$で、両者は$\operatorname{Ad}_{\exp z} = \exp(\operatorname{ad}_z)$（アダマールの補題）で結ばれます。$\mathfrak{su}(2)$では$\operatorname{ad}_{i/2} = J_x$となり、同型$\mathfrak{su}(2) \cong \mathfrak{so}(3)$の対応$i/2 \Leftrightarrow J_x$は$x \mapsto \operatorname{ad}_x$そのもの、二重被覆$q \mapsto \rho_q$は随伴表現$q \mapsto \operatorname{Ad}_q$そのものでした。対応が2対1になる理由も、核$\ker\operatorname{Ad} = \{\pm 1\}$（群の中心）として構造的に整理されます。$\operatorname{ad}$の行列成分は構造定数であり、随伴表現の次元はリー代数の次元です。

これらの道具はすべて結合法則の上に成り立っています。
