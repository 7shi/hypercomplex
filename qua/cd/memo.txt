https://x.com/7shi/status/897002382118100993

1. 複素数の積で片方を共役にすると内積と外積が得られる。
   $$(a_1+a_2i)^*(b_1+b_2i)=(a_1b_1+a_2b_2)+(a_1b_2-a_2b_1)i$$
2. 複素数の行列表現の成分を複素数にすると双複素数が得られる。

1と2を組み合わせると、行列表現から四元数が得られることに気付いた！

四元数が複素数の自然な拡張だと言われるけど、ノルムで $|AB|=|A||B|$ と説明されるのはピンと来なかった。さっき気付いた行列表現からの組み立て方だと私にはしっくり来た。これなら四元数の積から内積と外積が得られることも、そうやって作られているからと自然に解釈できる。

ここでつぶやいていた複素数・双複素数・四元数の行列表現について数式を使ってまとめました。厳密な導出ではなく発見法的です。

【補足】

複素数a,bは、aを共役にすると内積と外積が得られます。

$$
a^*b = (a_1 - a_2 i)(b_1 + b_2 i) = (a_1 b_1 + a_2 b_2) + (a_1 b_2 - a_2 b_1) i
$$

行列表現で計算します。

$$
\begin{pmatrix} a_1 & -a_2 \\ a_2 & a_1 \end{pmatrix}^\mathsf{T}
\begin{pmatrix} b_1 \\ b_2 \end{pmatrix}
=
\begin{pmatrix} a_1 & a_2 \\ -a_2 & a_1 \end{pmatrix}
\begin{pmatrix} b_1 \\ b_2 \end{pmatrix}
=
\begin{pmatrix} a_1 b_1 + a_2 b_2 \\ -a_2 b_1 + a_1 b_2 \end{pmatrix}
$$

複素数の行列表現と同じ形の行列に入れます。

$$
\begin{pmatrix} a & -b \\ b & a \end{pmatrix}
$$

これは双複素数$a+bj$です。エルミート共役$a^*-b^*j$との積を計算します。

$$
\begin{pmatrix} a & -b \\ b & a \end{pmatrix}^\dagger
\begin{pmatrix} c \\ d \end{pmatrix}
=
\begin{pmatrix} a^* & b^* \\ -b^* & a^* \end{pmatrix}
\begin{pmatrix} c \\ d \end{pmatrix}
=
\begin{pmatrix} a^* c + b^* d \\ -(b c^* - a d^*)^* \end{pmatrix}
$$

エルミート共役ではなく、四元数と同様に$(a+bj)^*=a^*-bj$となるように調整します。

$$
\begin{pmatrix} a^* & b \\ -b & a^* \end{pmatrix}
\begin{pmatrix} c \\ d \end{pmatrix}
=
\begin{pmatrix} a^* c + b d \\ -(b c - a^* d)^* \end{pmatrix}
$$

成分の共役の位置を動かして、ベクトルは1列目と同じフォーマットにすれば、四元数になります。

$$
\begin{pmatrix} a & b \\ -b^* & a^* \end{pmatrix}
\begin{pmatrix} c \\ -d^* \end{pmatrix}
=
\begin{pmatrix} a c - b d^* \\ -(b c^* + a d)^* \end{pmatrix}
$$




---

複素数における共役との積

$$
a, b \in \mathbb{C} \\
a_1, a_2, b_1, b_2 \in \mathbb{R} \\
a = a_1 + a_2 i \\
a^* = a_1 - a_2 i \\
b = b_1 + b_2 i \\
\begin{aligned}
a^* b &= (a_1 + a_2 i)^* (b_1 + b_2 i) \\
&= (a_1 - a_2 i)(b_1 + b_2 i) \\
&= \underset{内積}{\underbrace{(a_1 b_1 + a_2 b_2)}} + \underset{外積}{\underbrace{(a_1 b_2 - a_2 b_1)}} i
\end{aligned}
$$

---

複素数 $a$ の行列表現（作用素）

$$
a \mapsto \begin{pmatrix} a_1 & -a_2 \\ a_2 & a_1 \end{pmatrix}
$$

複素数 $b$ のベクトル表現

$$
b \mapsto \begin{pmatrix} b_1 \\ b_2 \end{pmatrix}
$$

行列・ベクトル表現による共役との積

$$
\begin{aligned}
a^* b \mapsto & \begin{pmatrix} a_1 & a_2 \\ -a_2 & a_1 \end{pmatrix} \begin{pmatrix} b_1 \\ b_2 \end{pmatrix} \\
= & \begin{pmatrix} a_1 b_1 + a_2 b_2 \\ -a_2 b_1 + a_1 b_2 \end{pmatrix} \\
\mapsto & (a_1 b_1 + a_2 b_2) + (a_1 b_2 - a_2 b_1) i
\end{aligned}
$$

---

表現行列の成分を複素数にすると双複素数 $\mathbb{C}^2$ が得られる

$$
\begin{aligned}
& \begin{pmatrix} a_1 + a_2 i & -(a_3 + a_4 i) \\ a_3 + a_4 i & a_1 + a_2 i \end{pmatrix} \\
& \mapsto (a_1 + a_2 i) + (a_3 + a_4 i) j \\
&= a_1 + a_2 i + a_3 j + a_4 i j
\end{aligned}
$$

$i j = k$ とおけば

$$
a_1 + a_2 i + a_3 j + a_4 k
$$

双複素数は可換

$$
\begin{aligned}
i^2 = j^2 &= -1 \\
k^2 = (ij)^2 = i^2 j^2 &= 1 \\
ij = ji &= k \\
jk = kj &= -i \\
ki = ik &= -j
\end{aligned}
$$

参考: [Bicomplex_number](https://en.wikipedia.org/wiki/Bicomplex_number)

---

双複素数の行列表現の成分の複素数をネストした行列表現にして転置することで共役が得られる

$$
\begin{aligned}
& (a_1 + a_2 i + a_3 j + a_4 k)^* \\
\mapsto & \begin{pmatrix} a_1 & -a_2 & -a_3 & a_4 \\ a_2 & a_1 & -a_4 & -a_3 \\ a_3 & -a_4 & a_1 & -a_2 \\ a_4 & a_3 & a_2 & a_1 \end{pmatrix}^\top \\
= & \begin{pmatrix} a_1 & a_2 & a_3 & a_4 \\ -a_2 & a_1 & -a_4 & a_3 \\ -a_3 & -a_4 & a_1 & a_2 \\ a_4 & -a_3 & -a_2 & a_1 \end{pmatrix} \\
\mapsto & (a_1 - a_2 i) - (a_3 - a_4 i) j \\
= & a_1 - a_2 i - a_3 j + a_4 k
\end{aligned}
$$

---

双複素数での共役との積は複雑

$$
\begin{aligned}
a^* b
\mapsto & \begin{pmatrix} (a_1 + a_2 i)^* & (a_3 + a_4 i)^* \\ -(a_3 + a_4 i)^* & (a_1 + a_2 i)^* \end{pmatrix} \begin{pmatrix} b_1 + b_2 i \\ b_3 + b_4 i \end{pmatrix} \\
= & \begin{pmatrix} (a_1 + a_2 i)^* (b_1 + b_2 i) + (a_3 + a_4 i)^* (b_3 + b_4 i) \\ -(a_3 + a_4 i)^* (b_1 + b_2 i) + (a_1 + a_2 i)^* (b_3 + b_4 i) \end{pmatrix} \\
\mapsto & (a_1 b_1 + a_2 b_2 + a_3 b_3 + a_4 b_4) \\
& + (a_1 b_2 - a_2 b_1 + a_3 b_4 - a_4 b_3) i \\
& + (a_1 b_3 - a_3 b_1 - a_4 b_2 + a_2 b_4) j \\
& + (a_1 b_4 + a_4 b_1 - a_2 b_3 - a_3 b_2) k
\end{aligned}
$$
$$
a^* b \mapsto
\begin{pmatrix} p^* & q^* \\ -q^* & p^* \end{pmatrix}
\begin{pmatrix} r \\ s \end{pmatrix}
= \begin{pmatrix} p^* r + q^* s \\ -(q r^* - p s^*)^* \end{pmatrix}
$$

【参考】四元数の積

$$
\begin{aligned}
ab
\mapsto & \begin{pmatrix} a_1 + a_2 i & a_3 + a_4 i \\ -(a_3 + a_4 i)^* & (a_1 + a_2 i)^* \end{pmatrix} \begin{pmatrix} b_1 + b_2 i \\ -(b_3 + b_4 i)^* \end{pmatrix} \\
= & \begin{pmatrix} (a_1 + a_2 i)(b_1 + b_2 i) - (a_3 + a_4 i)(b_3 + b_4 i)^* \\ -(a_3 + a_4 i)^*(b_1 + b_2 i) - (a_1 + a_2 i)^*(b_3 + b_4 i)^* \end{pmatrix} \\
\mapsto & (a_1 b_1 - a_2 b_2 - a_3 b_3 - a_4 b_4) \\
& + (a_1 b_2 + a_2 b_1 + a_3 b_4 - a_4 b_3) i \\
& + (a_1 b_3 + a_3 b_1 + a_4 b_2 - a_2 b_4) j \\
& + (a_1 b_4 + a_4 b_1 + a_2 b_3 - a_3 b_2) k
\end{aligned}
$$
$$
ab \mapsto
\begin{pmatrix} p & q \\ -q^* & p^* \end{pmatrix}
\begin{pmatrix} r \\ -s^* \end{pmatrix}
= \begin{pmatrix} pr - qs^* \\ -(qr^* + ps)^* \end{pmatrix}
$$

---

複素数では共役により反対称性を模倣して外積を求めている。

$$
\begin{aligned}
& (a_1 + a_2 i)(b_1 + b_2 i) \\
&= \underset{ミンコフスキー計量}{\underbrace{(a_1 b_1 - a_2 b_2)}} + \underset{対称}{\underbrace{(a_1 b_2 + a_2 b_1)}} i \\
& \overset{共役}{\overbrace{(a_1 + a_2 i)^*}} (b_1 + b_2 i) \\
&= \underset{ユークリッド計量}{\underbrace{(a_1 b_1 + a_2 b_2)}} + \underset{反対称}{\underbrace{(a_1 b_2 - a_2 b_1)}} i
\end{aligned}
$$

内積には計量が現れているようだが、ここでは追求しない。

---

複素数での共役による虚部の反対称性は、順番を入れ替えることで虚部の符号が変わることで実現されている。

$$
\begin{aligned}
a^* b &= (a_1 - a_2 i)(b_1 + b_2 i) \\
b^* a &= (a_1 + a_2 i)(b_1 - b_2 i)
\end{aligned}
$$

双複素数では $i j = k$ により $k$ の符号反転が相殺されてしまうため、完全な反対称性は持たない。

$$
\begin{aligned}
a^* b &= (a_1 - a_2 i - a_3 j + a_4 k)(b_1 + b_2 i + b_3 j + b_4 k) \\
b^* a &= (a_1 + a_2 i + a_3 j + a_4 k)(b_1 - b_2 i - b_3 j + b_4 k)
\end{aligned}
$$

双複素数はそういうものだと見切りをつけて、反対称を実現する別の方法を考える。

---

既存の行列表現から共役を使わない方法を見付ける。

複素数とその共役

$$
a \mapsto \begin{pmatrix} a_1 & -a_2 \\ a_2 & a_1 \end{pmatrix} \\
a^* \mapsto \begin{pmatrix} a_1 & a_2 \\ -a_2 & a_1 \end{pmatrix}
$$

双複素数とその共役

$$
a \mapsto \begin{pmatrix} a_1 + a_2 i & -(a_3 + a_4 i) \\ a_3 + a_4 i & a_1 + a_2 i \end{pmatrix} \\
a^* \mapsto \begin{pmatrix} (a_1 + a_2 i)^* & (a_3 + a_4 i)^* \\ -(a_3 + a_4 i)^* & (a_1 + a_2 i)^* \end{pmatrix}
$$

双複素数とその共役を折衷すると、四元数の行列表現の一種が得られる。

$$
\begin{pmatrix} (a_1 + a_2 i)^* & -(a_3 + a_4 i) \\ (a_3 + a_4 i)^* & a_1 + a_2 i \end{pmatrix}
$$

---

積を確認。共役を使わなくても虚部に反対称成分が現れる。

$$
\begin{aligned}
& \begin{pmatrix} (a_1 + a_2 i)^* & -(a_3 + a_4 i) \\ (a_3 + a_4 i)^* & a_1 + a_2 i \end{pmatrix} \begin{pmatrix} (b_1 + b_2 i)^* \\ (b_3 + b_4 i)^* \end{pmatrix} \\
& = \begin{pmatrix} (a_1 b_1 - a_2 b_2 - a_3 b_3 - a_4 b_4) + \left(-a_1 b_2 - a_2 b_1 + \underset{反対称}{\underbrace{a_3 b_4 - a_4 b_3}}\right) i \\ (a_1 b_3 + a_3 b_1 + \underset{反対称}{\underbrace{a_2 b_4 - a_4 b_2}}) + \left(-a_1 b_4 - a_4 b_1 + \underset{反対称}{\underbrace{a_2 b_3 - a_3 b_2}}\right) i \end{pmatrix}
\end{aligned}
$$

四元数の内部で複素数の共役が働いている。

---

発見法的に四元数の行列表現に結びつけたが、やや複雑な対応になっている。

$$
\begin{aligned}
a + jb = & (a_1 + a_2 k) + j(b_1 + b_2 k) \\
= & a_1 + b_2 i + b_1 j + a_2 k \\
\mapsto & \begin{pmatrix} (a_1 + a_2 k)^* & -(b_1 + b_2 k) \\ (b_1 + b_2 k)^* & a_1 + a_2 k \end{pmatrix} \\
= & \begin{pmatrix} a^* & -b \\ b^* & a \end{pmatrix}
\end{aligned}
$$

四元数を双四元数に拡張すれば、この表現は素直にパウリ行列に対応する。

$$
\begin{aligned}
h^2 &= -1 \\
hi &\mapsto \sigma_1 \\
hj &\mapsto \sigma_2 \\
hk &\mapsto \sigma_3 \\
i = -h(hi) &\mapsto -i\sigma_1 \\
j = -h(hj) &\mapsto -i\sigma_2 \\
k = -h(hk) &\mapsto -i\sigma_3
\end{aligned}
$$

参考: [arxiv.org/abs/0709.2238](https://arxiv.org/abs/0709.2238)

---

四元数だけで考えると、もう少し簡単な表現がある。

$$
\begin{aligned}
a + bj = & (a_1 + a_2 i) + (b_1 + b_2 i) j \\
= & a_1 + a_2 i + b_1 j + b_2 k \\
\mapsto & \begin{pmatrix} a_1 + a_2 i & b_1 + b_2 i \\ -(b_1 + b_2 i)^* & (a_1 + a_2 i)^* \end{pmatrix} \\
= & \begin{pmatrix} a & b \\ -b^* & a^* \end{pmatrix}
\end{aligned}
$$

ケイリー=ディクソン構成に対応（積の順序は無視）

$$
\begin{pmatrix} a & b \\ -b^* & a^* \end{pmatrix} \begin{pmatrix} c \\ -d^* \end{pmatrix}
= \begin{pmatrix} ac - bd^* \\ -(bc^* + ad)^* \end{pmatrix} \\
(a, b)(c, d) = (ac - bd^*, ad + bc^*)
$$

---

一連の発言は、四元数の積から内積と外積が得られる仕組みを確認するため、行列表現の再構成を試みたものです。

$$
a, b \in \mathbb{H} \\
a_0, a_1, a_2, a_3, b_0, b_1, b_2, b_3 \in \mathbb{R} \\
\begin{aligned}
ab &= (a_0 + a_1 i + a_2 j + a_3 k)(b_0 + b_1 i + b_2 j + b_3 k) \\
&= \underset{実部の積}{\underbrace{a_0 b_0}} \\
&\quad + \underset{a の虚部の定数倍}{\underbrace{b_0 (a_1 i + a_2 j + a_3 k)}} \\
&\quad + \underset{b の虚部の定数倍}{\underbrace{a_0 (b_1 i + b_2 j + b_3 k)}} \\
&\quad - \underset{内積}{\underbrace{(a_1 b_1 + a_2 b_2 + a_3 b_3)}} \\
&\quad + \underset{外積}{\underbrace{(a_2 b_3 - a_3 b_2) i + (a_3 b_1 - a_1 b_3) j + (a_1 b_2 - a_2 b_1) k}}
\end{aligned}
$$

虚部の反対称性より、外積は虚部で構成されます。内積の符号は2ベクトルで説明されます。

---

ケイリー=ディクソン構成は、四元数を2つの複素数に分割して扱う場合に共役が付くことを表している。

$$
\begin{aligned}
& a_{01}, a_{23}, b_{01}, b_{23} \in \mathbb{C} \\
& (a_{01} + a_{23} j)(b_{01} + b_{23} j) \\
\text{双複素数}\quad & = (a_{01} b_{01} - a_{23} b_{23}) + (a_{01} b_{23} + a_{23} b_{01}) j \\
\text{四元数}\quad & = (a_{01} b_{01} - a_{23} b_{23}^*) + (a_{01} b_{23} + a_{23} b_{01}^*) j
\end{aligned}
$$

【補足】

双複素数では $i$ と $j$ が可換なため、展開しても共役は現れない。四元数では $ji = -ij$ より、複素数 $z = x + yi$ が $j$ を追い越すと共役が付く。

$$
jz = xj + y(ji) = xj - y(ij) = (x - yi)j = z^* j
$$

これを使って展開すると、$j$ の右側にあった $b_{01}, b_{23}$ にだけ共役が付く。つまりケイリー=ディクソン構成の共役は、$b$ 側の因子が $j$ を追い越したことの痕跡である。

$$
\begin{aligned}
(a_{01} + a_{23} j)(b_{01} + b_{23} j)
&= a_{01} b_{01} + a_{01} b_{23} j + a_{23} (j b_{01}) + a_{23} (j b_{23}) j \\
&= a_{01} b_{01} + a_{01} b_{23} j + a_{23} b_{01}^* j + a_{23} b_{23}^* j^2 \\
&= (a_{01} b_{01} - a_{23} b_{23}^*) + (a_{01} b_{23} + a_{23} b_{01}^*) j
\end{aligned}
$$

---

複素数の中身を展開して動きを追うと、反対称性は意外と手の込んだ組み合わせで実現されていることが分かる。

$$
a_0, a_1, a_2, a_3, b_0, b_1, b_2, b_3 \in \mathbb{R} \\
\begin{aligned}
a_{01} b_{01} &= (a_0 b_0 - a_1 b_1) + (a_0 b_1 + a_1 b_0) i \\
-a_{23} b_{23}^* &= (-a_2 b_2 - a_3 b_3) + (a_2 b_3 - a_3 b_2) i \\
a_{01} b_{23} j &= (a_0 b_2 - a_1 b_3) j + (a_0 b_3 + a_1 b_2) k \\
a_{23} b_{01}^* j &= (a_2 b_0 + a_3 b_1) j + (a_3 b_0 - a_2 b_1) k
\end{aligned}
$$

$b_{23}$ は $j k = -k j = i$  
$b_{01}^*$ は $k i = -i k = j, \ i j = -j i = k$ の反交換性を担っている。

【補足】

ここでの「反対称性」とは、外積にあたる項 $(a_2 b_3 - a_3 b_2)i,\ (a_3 b_1 - a_1 b_3)j,\ (a_1 b_2 - a_2 b_1)k$ を指す。上の展開を基底の積に戻すと、どの共役がどの反対称項を作っているかが分かる。

$i$ の反対称項は $-a_{23} b_{23}^*$ の中にまとまって現れる。これは $j, k$ 成分同士の積であり、共役 $b_{23}^*$ が $jk = -kj = i$ の符号差を作っている。

$$
(a_2 j)(b_3 k) = a_2 b_3\, jk = a_2 b_3\, i, \quad
(a_3 k)(b_2 j) = a_3 b_2\, kj = -a_3 b_2\, i
$$

もし共役がなければ双複素数の $jk = kj = -i$ となり、対称項 $-(a_2 b_3 + a_3 b_2)i$ になってしまう。

一方、$j, k$ の反対称項は1つの積に収まらず、$a_{01} b_{23} j$ と $a_{23} b_{01}^* j$ の2つに分かれて現れる。ここでは共役 $b_{01}^*$ が $ki = -ik = j$ と $ij = -ji = k$ の符号差を作っている。

$$
\begin{alignedat}{2}
(a_3 k)(b_1 i) &= a_3 b_1\, ki = a_3 b_1\, j, &\quad
(a_1 i)(b_3 k) &= a_1 b_3\, ik = -a_1 b_3\, j \\
(a_1 i)(b_2 j) &= a_1 b_2\, ij = a_1 b_2\, k, &\quad
(a_2 j)(b_1 i) &= a_2 b_1\, ji = -a_2 b_1\, k
\end{alignedat}
$$

---

四元数で2通りの行列表現を示したが、エルミート共役で対応している。

$$
\begin{pmatrix} a^* & -b \\ b^* & a \end{pmatrix}^\dagger = \begin{pmatrix} a & b \\ -b^* & a^* \end{pmatrix}
$$
