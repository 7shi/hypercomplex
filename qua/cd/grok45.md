リー群の記事では、四元数$q = z + wj$に対して複素行列表現

$$
\begin{pmatrix}
z & w \\
-w^* & z^*
\end{pmatrix}
$$

を天下りに導入しました[[7shi-lie2]]。本記事では、複素数の性質からこの種の表現がどのように動機付けられるかを、ヒューリスティックな導出としてまとめます。厳密な一意性の証明ではなく、「なぜこの形を取るのが自然か」を辿る道筋です。

出発点は次の2つです。

1. 複素数の積で片方を共役にすると、2次元ベクトルの**内積**と**外積**が同時に得られる
2. 複素数の行列表現の成分を複素数にすると、**双複素数**が得られる

この2つを組み合わせると、行列表現から四元数が現れます。

# 複素数の行列表現

復習から始めます。複素数$z = a + ib$（$a,b \in \mathbb{R}$）は、次の2×2実行列と対応します。

$$
z = a + ib
\quad \Leftrightarrow \quad
\begin{pmatrix}
a & -b \\
b & a
\end{pmatrix}
$$

この対応は積を保ち、絶対値の2乗は行列式に一致します。

$$
|z|^2 = a^2 + b^2 = \det\begin{pmatrix} a & -b \\ b & a \end{pmatrix}
$$

図形的には、複素数の積は平面$\mathbb{R}^2$上の回転と拡大縮小であり、その線形変換が上の行列です。

# ① 共役との積：内積と外積

2つの複素数$\alpha = a_1 + a_2 i$、$\beta = b_1 + b_2 i$を、平面上のベクトル$(a_1, a_2)$、$(b_1, b_2)$と同一視します。通常の積$\alpha\beta$は回転と拡大を混ぜた演算で、ベクトルの幾何そのものではありません。

一方、**片方を共役にしてから掛ける**と、様子が変わります。

$$
\alpha^*\beta
= (a_1 - a_2 i)(b_1 + b_2 i)
= (a_1 b_1 + a_2 b_2) + (a_1 b_2 - a_2 b_1)\,i
$$

実部は内積、虚部は（符号付き面積としての）外積です。

$$
\alpha^*\beta
= \alpha \cdot \beta + i\,(\alpha \wedge \beta)
$$

&&&def 複素数の共役積と内積・外積
$$
\begin{aligned}
\alpha \cdot \beta &= a_1 b_1 + a_2 b_2 = \operatorname{Re}(\alpha^*\beta) \\
\alpha \wedge \beta &= a_1 b_2 - a_2 b_1 = \operatorname{Im}(\alpha^*\beta)
\end{aligned}
$$
&&&

特に$\beta = \alpha$とすれば、$|\alpha|^2 = \alpha^*\alpha$となり、行列表現の行列式$a^2 + b^2$が共役との積として読み直せます。

$$
|\alpha|^2 = \alpha^*\alpha = \det\begin{pmatrix} a_1 & -a_2 \\ a_2 & a_1 \end{pmatrix}
$$

共役は、複素数の積のなかにユークリッド幾何（長さ・角度・符号付き面積）を取り出すためのスイッチです。この観察が、後で行列の成分を拡張するときの指針になります。

&&&rem 四元数での内積・外積
純虚四元数の積でも同様の分解が現れます。$p,q$を純虚四元数とすれば$pq = -\,p\cdot q + p\times q$であり、実部が内積（の符号反転）、虚部が外積です。2次元では共役積が同じ役割を果たしており、四元数はその3次元版と見ることができます。
&&&

# ② 成分の複素化：双複素数

複素数の行列表現で、実数成分$a,b$を複素数$z,w$に置き換えてみます。

$$
\begin{pmatrix}
a & -b \\
b & a
\end{pmatrix}
\quad\longrightarrow\quad
\begin{pmatrix}
z & -w \\
w & z
\end{pmatrix}
$$

こうして得られる4実次元の代数が**双複素数**です。一般の元は$z + w\,j'$の形に書け、$j'$を「もう一つの虚数単位」と見なせます。ただし、この$j'$は行列

$$
j' \Leftrightarrow \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
$$

として既存の複素構造の**外側**にあり、成分の$i$とは可換です。実際、

$$
\begin{pmatrix} i & 0 \\ 0 & i \end{pmatrix}
\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
=
\begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}
=
\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
\begin{pmatrix} i & 0 \\ 0 & i \end{pmatrix}
$$

となり、$i$と$j'$は交換します。双複素数は可換な代数です。

&&&rem 双複素数とテンソル積
双複素数は$\mathbb{C}\otimes\mathbb{C}$としても定義できます。2つの複素構造が互いに可換に並ぶため、可換性が残ります[[wiki-bc]]。
&&&

この形の行列式を計算すると、

$$
\det\begin{pmatrix} z & -w \\ w & z \end{pmatrix} = z^2 + w^2
$$

となります。実数のとき$a^2 + b^2$だったものが$z^2 + w^2$に化けており、これは一般に非負実数とは限りません。$z = 1$、$w = i$とすれば$z^2 + w^2 = 0$となり、零因子も現れます。ユークリッド的な絶対値$|z|^2 + |w|^2$とは別物です。

素朴な複素化は、代数としては双複素数を与えますが、複素数が持っていた「行列式＝絶対値の2乗」という幾何的な性質は壊れてしまいます。

## 双複素数と共役の折衷 (Terra)

ここまでの行列表現そのものから、共役をどこへ入れるかを探る見方もあります。複素数$a = a_1 + a_2 i$とその共役は、それぞれ

$$
a \longmapsto
\begin{pmatrix}
a_1 & -a_2 \\
a_2 & a_1
\end{pmatrix},
\qquad
a^* \longmapsto
\begin{pmatrix}
a_1 & a_2 \\
-a_2 & a_1
\end{pmatrix}
$$

と表せます。すなわち共役は、対角成分を保ったまま非対角成分の符号を反転する操作です。

この対応で実成分を複素数にした双複素数を

$$
a = a_1 + a_2 i + (a_3 + a_4 i)j'
= z + wj'
$$

と書きます。ここで$z = a_1 + a_2i$、$w = a_3 + a_4i$です。$i$と$j'$の両方を反転する共役を$*$とすれば、双複素数とその共役は

$$
a \longmapsto
\begin{pmatrix}
z & -w \\
w & z
\end{pmatrix},
\qquad
a^* \longmapsto
\begin{pmatrix}
z^* & w^* \\
-w^* & z^*
\end{pmatrix}
$$

となります。前者の非対角成分と後者の対角成分を採る、という折衷をすると、

$$
\begin{pmatrix}
z^* & -w \\
w^* & z
\end{pmatrix}
=
\begin{pmatrix}
(a_1+a_2i)^* & -(a_3+a_4i) \\
(a_3+a_4i)^* & a_1+a_2i
\end{pmatrix}
$$

が得られます。行列式は

$$
z^*z - (-w)w^* = |z|^2 + |w|^2
$$

であり、双複素数では失われたユークリッドノルムが回復します。また、この4実次元の行列全体は積について閉じ、四元数の行列表現の一種になります。これは「双複素数に共役を付け足す」と考える代わりに、元の行列表現と共役された行列表現の間から四元数の形を読み取る見方です。

&&&rem 規約の違い
この折衷で得た行列は、次節で用いる$Q$とは対角成分の$z,z^*$が入れ替わった規約です。$i$に対応する基底行列の符号、ひいては$ij=k$か$ij=-k$かという向きが反転するだけで、いずれも同型な四元数代数を表します。
&&&

# ①と②の組み合わせ：四元数

①では、ユークリッド的な長さが共役との積$\alpha^*\alpha$として現れることを見ました。②の素朴な複素化では、行列式が$z^2 + w^2$になってこの性質が失われました。

では、行列式が

$$
|z|^2 + |w|^2 = zz^* + ww^*
$$

となるよう、成分に共役を入れてみます。雛形は複素数の行列そのもの

$$
\begin{pmatrix} a & -b \\ b & a \end{pmatrix}
$$

です。上段の形$(z\ -w)$はそのままにし、下段だけを共役で閉じる（すなわち$w$を$w^*$、$z$を$z^*$に置き換える）と、次の形が得られます。

&&&def 共役を入れた複素行列
$$
Q =
\begin{pmatrix}
z & -w \\
w^* & z^*
\end{pmatrix}
$$
&&&

これは$\begin{pmatrix}a & -b \\ b & a\end{pmatrix}$の**素直な拡張**です。$w$が実数（したがって$w^* = w$、$z^* = z$）のときは、もとの実行列にそのまま戻ります。

行列式は狙い通り、

$$
\det Q = zz^* - (-w)w^* = zz^* + ww^* = |z|^2 + |w|^2
$$

となり、常に非負の実数です。実数のときの$a^2 + b^2$が、共役積を通じて自然に4実成分のユークリッドノルムへ拡張されています。①で見た「共役が長さを取り出す」という仕組みが、行列の成分配置に書き込まれた形です。

&&&rem リー群の記事との規約
リー群の記事[[7shi-lie2]]では$\begin{pmatrix}z & w \\ -w^* & z^*\end{pmatrix}$を採用しています。これは本文の$Q$で$w$を$-w$に置き換えたものにすぎず、現れる行列の集合は同一です。リー群側では、1行目が$(z\ w)$と読めること、および基底がパウリ行列に$i$を掛けるだけで揃うこと（$i \leftrightarrow i\sigma_3,\ j \leftrightarrow i\sigma_2,\ k \leftrightarrow i\sigma_1$）を優先した規約です。本記事では、複素数の行列からの拡張として自然な形を主に用います。
&&&

## 基底の取り出し

$z = a + bi$、$w = c + di$（$a,b,c,d \in \mathbb{R}$）とおけば、

$$
Q =
\begin{pmatrix}
a+bi & -c-di \\
c-di & a-bi
\end{pmatrix}
=
a\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}
+ b\begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}
+ c\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
+ d\begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}
$$

となります。4つの基底行列を

$$
1 \Leftrightarrow \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}, \quad
i \Leftrightarrow \begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}, \quad
j \Leftrightarrow \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \quad
k \Leftrightarrow \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}
$$

と書けば、$Q$は$a + bi + cj + dk$に対応します。

$j$の行列は、複素数の虚数単位$i$の実行列表現$\begin{pmatrix}0 & -1 \\ 1 & 0\end{pmatrix}$そのものです。素朴な複素化で現れた$j'$と同じ配置が、共役を入れたあとも「成分間の回転」として残っています。

## 積の規則

これらの行列が四元数の関係式を満たすことを確認します。

$$
\begin{aligned}
i^2
&= \begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}^2
= \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}
= -1
\\
j^2
&= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}^2
= \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}
= -1
\\
k^2
&= \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}^2
= \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}
= -1
\end{aligned}
$$

$$
\begin{aligned}
ij
&= \begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}
\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
= \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}
= k
\\
ji
&= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
\begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}
= \begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix}
= -k
\end{aligned}
$$

$ij = -ji = k$が成り立ち、同様に$jk = -kj = i$、$ki = -ik = j$も得られます。すなわち、この行列のなす代数は四元数$\mathbb{H}$です。

双複素数では$i$と$j'$が可換だったのに対し、ここでは対角の$i$と非対角の$j$が反交換します。共役を片方の対角成分だけに入れたことで、$i$の作用が上下の成分で符号反転し、非可換性が生まれたわけです。

&&&rem パウリ行列との対応
パウリ行列を
$$
\sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
\sigma_2 = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad
\sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$
とすれば、上の基底は
$$
i \leftrightarrow i\sigma_3, \quad
j \leftrightarrow -i\sigma_2, \quad
k \leftrightarrow -i\sigma_1
$$
となります。$j,k$にマイナスが付くのは、複素数の行列を素直に拡張した規約の帰結です。リー群の記事の規約（$w$の符号反転）では、このマイナスが消えて$i\sigma_3, i\sigma_2, i\sigma_1$と揃います。
&&&

&&&rem なぜ対角が$z$と$z^*$なのか
$i$に対応する行列$\operatorname{diag}(i,-i)$をベクトル$\begin{pmatrix}r\\s\end{pmatrix}$に掛けると、$\begin{pmatrix}ir\\-is\end{pmatrix}$となります。上成分と下成分で$i$の向きが逆です。これが$jz = z^*j$という四元数の関係（$j$を越えると複素係数が共役になる）の行列版であり、下段の$w^*$も同じ共役構造から来ています。
&&&

## $q = z + wj$との対応

四元数を$q = z + wj$と書いて積を展開すると、共役の必然性が代数側からも見えます。

$$
\begin{aligned}
(z_1 + w_1 j)(z_2 + w_2 j)
&= z_1 z_2 + z_1 w_2\,j + w_1 j\,z_2 + w_1 j\,w_2\,j \\
&= z_1 z_2 + z_1 w_2\,j + w_1 z_2^*\,j + w_1 w_2\,j^2 \\
&= (z_1 z_2 - w_1 w_2^*) + (z_1 w_2 + w_1 z_2^*)\,j
\end{aligned}
$$

ここで$jz = z^*j$と$j^2 = -1$を使いました。右辺に現れる$w_2^*$や$z_2^*$が、行列表現の下段の共役成分と正確に対応します。実際、行列の積

$$
\begin{pmatrix} z_1 & -w_1 \\ w_1^* & z_1^* \end{pmatrix}
\begin{pmatrix} z_2 & -w_2 \\ w_2^* & z_2^* \end{pmatrix}
=
\begin{pmatrix}
z_1 z_2 - w_1 w_2^* &
-(z_1 w_2 + w_1 z_2^*) \\
(z_1 w_2 + w_1 z_2^*)^* &
(z_1 z_2 - w_1 w_2^*)^*
\end{pmatrix}
$$

は、上で得た$(z_1 z_2 - w_1 w_2^*) + (z_1 w_2 + w_1 z_2^*)\,j$の行列表現そのものです。共役は天下りではなく、$j$をくぐるときの規則として積の中に強制されます。

また、四元数の共役$q^* = z^* - wj$は行列のエルミート共役に対応します。

$$
Q^\dagger
=
\begin{pmatrix}
z^* & w \\
-w^* & z
\end{pmatrix}
=
\begin{pmatrix}
z^* & -(-w) \\
(-w)^* & z
\end{pmatrix}
$$

右辺は$z^* + (-w)j = z^* - wj$の行列表現です。

# 双複素数と四元数の分岐

ここまでの道筋を整理します。

| | 素朴な複素化 | 共役を入れた複素化 |
|---|---|---|
| 行列 | $\begin{pmatrix}z & -w \\ w & z\end{pmatrix}$ | $\begin{pmatrix}z & -w \\ w^* & z^*\end{pmatrix}$ |
| 代数 | 双複素数（可換） | 四元数（非可換） |
| 行列式 | $z^2 + w^2$ | $\|z\|^2 + \|w\|^2$ |
| 零因子 | あり | なし（$Q \ne 0$なら可逆） |
| 幾何 | ユークリッドノルムを失う | ユークリッドノルムを保つ |

同じ「成分を複素数にする」という操作でも、共役を入れるかどうかで可換な双複素数と非可換な四元数に分岐します。分岐を決める指針が①でした。複素数において長さと角度を取り出す共役積を、行列の下段に書き込むこと。それが四元数の行列表現を導くヒューリスティックです。

素朴な複素化との違いは下段の共役だけで、上段の$(z\ -w)$は複素数の行列から一度も崩していません。

&&&rem ノルムの乗法性
四元数では$|pq| = |p||q|$が成り立ち、行列側では$\det(AB) = \det A\,\det B$がこれに対応します。$\det Q = |q|^2$なので、両者は同じ乗法的ノルムを見ています。双複素数の$z^2 + w^2$にはこの正定値な乗法性がありません。
&&&

# まとめ

複素数から四元数の行列表現へ至るヒューリスティックな道筋は、次のようにまとめられます。

1. **複素数の行列表現** $\begin{pmatrix}a & -b \\ b & a\end{pmatrix}$ が出発点である
2. **① 共役積** $\alpha^*\beta$ は内積と外積を同時に与え、特に$|\alpha|^2 = \alpha^*\alpha$としてユークリッド的な長さを取り出す
3. **② 素朴な複素化**は双複素数$\begin{pmatrix}z & -w \\ w & z\end{pmatrix}$を与えるが、可換性が残り、行列式$z^2 + w^2$はユークリッドノルムにならない
4. **①と②の組み合わせ**として、上段は保ったまま下段に共役を入れると、四元数の行列表現$\begin{pmatrix}z & -w \\ w^* & z^*\end{pmatrix}$が得られる
5. この表現は積を保ち、$i^2 = j^2 = k^2 = -1$、$ij = k$などの四元数の規則を再現する
6. リー群の記事の形$\begin{pmatrix}z & w \\ -w^* & z^*\end{pmatrix}$は、$w$の符号を反転した同等な規約である

複素数の幾何（共役積）と複素化という2つの操作を組み合わせた帰結として、四元数の行列表現が得られます。単位元の近傍では、この表現（およびそれと同等な規約）が$\operatorname{SU}(2)$を与え、3次元回転$\operatorname{SO}(3)$との関係へとつながっていきます[[7shi-lie2]]。
