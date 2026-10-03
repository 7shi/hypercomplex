複素数の虚数単位を増やして拡張する試みから、双複素数の代数構造と零因子の性質を整理し、テンソル積の自然な導入を解説します。

# 概要

複素数に可換な第2の虚数単位を導入すると、4次元の双複素数（bicomplex number）の体系が得られます。双複素数は零因子を含むため多元体にはなりませんが、複素数同士のテンソル積$\mathbb{C} \otimes \mathbb{C}$と自然に同型対応します。本記事では、虚数単位の導入から双複素数の乗算規則と零因子の発生機構を確認し、成分ごとの積としてテンソル積代数が構成される過程を具体例を通じて明らかにします。本記事のテンソル積は、すべて実数体$\mathbb{R}$上のものです。

# 虚数単位を増やす

複素数の虚数単位$i$について$i^2=-1$が成り立ちます。ここに、同じ性質$j^2=-1$を満たす新しい虚数単位$j$を導入します。ただし、$i \ne \pm j$で、$i$と$j$は可換$ij=ji$とします。

積を計算します。（$a,b,c,d,e,f$は実数）

$$
\begin{aligned}
&(a+bi+cj)(d+ei+fj) \\
&= a(d+ei+fj) + bi(d+ei+fj) + cj(d+ei+fj) \\
&= ad + aei + afj + bdi + bei^2 + bfij + cdj + ceij + cfj^2 \\
&= ad + aei + afj + bdi - be + bfij + cdj + ceij - cf \\
&= (ad-be-cf) + (ae+bd)i + (af+cd)j + (bf+ce)ij
\end{aligned}
$$

この計算結果から、$i,j$によって生成される数は

$$
a+bi+cj+dij \quad (a,b,c,d\text{ は実数})
$$

の形で表せることが分かります。このような数の体系を**双複素数** (bicomplex number) と呼びます。[[wiki-bc]][[wiki-bc-en]]

&&&rem 三元数・四元数との違い

- $ij$は$a+bi+cj$の形では表せないことから、三元数ではありません。
- 双複素数は$ij=ji$であることから、$ij=-ji$となる四元数とは異なります。[[wiki-q]]

&&&

# 零因子の存在

$0$でない$2$つの双複素数の積が$0$になる場合があります。

$$
\begin{aligned}
&(1+ij)(1-ij) \\
&=1(1-ij)+ij(1-ij) \\
&=1-ij+ij-ijij \\
&=1-i^2j^2 \\
&=1-(-1)(-1) \\
&=1-1 \\
&=0
\end{aligned}
$$

$1+ij,\ 1-ij$のように、$0$ではないのに積が$0$になる数を**零因子**と呼びます。[[wiki-0d]]

零因子は逆数を持ちません。なぜなら、$1+ij$に何かを掛けて$1$になると仮定すれば、その両辺に$1-ij$を掛けることで矛盾が生じるためです。

$$
\begin{aligned}
(1+ij)x&=1 \\
(1-ij)(1+ij)x&=1-ij \\
0&=1-ij
\end{aligned}
$$

これは$1-ij$が$0$でないという前提に矛盾します。

&&&rem 前提条件との矛盾
もう少し式変形を進めると$ij=1,\ j=-i$となって、前提$i \ne \pm j$に矛盾します。
&&&

# テンソル積による構成

因子の順序を区別する積の演算子$\otimes$を導入して、複素数の積を計算します。（$a,b,c,d$は実数）

$$
\begin{aligned}
&(a+bi)\otimes(c+di) \\
&= a \otimes (c+di) + bi \otimes (c+di) \\
&= a \otimes c + a \otimes di + bi \otimes c + bi \otimes di
\end{aligned}
$$

これを双複素数と比較します。

$$
\begin{aligned}
&(a+bi)(c+dj) \\
&= a(c+dj) + bi(c+dj) \\
&= ac + adj + bci + bdij \\
\end{aligned}
$$

実部の基底として$1$を明示し、係数と基底 $\{1,i,j,ij\}$ を分離します。

$$
= ac(1) + ad(j) + bc(i) + bd(ij)
$$

$\otimes$でも係数と基底を分離して、係数が括り出せるという計算規則を追加します。

$$
\begin{aligned}
&= a(1) \otimes c(1) + a(1) \otimes d(i) + b(i) \otimes c(1) + b(i) \otimes d(i) \\
&= ac(1 \otimes 1) + ad(1 \otimes i) + bc(i \otimes 1) + bd(i \otimes i)
\end{aligned}
$$

以下の対応関係を認めれば、双複素数は$\otimes$による計算と一致します。

$$
1 \cong 1 \otimes 1, \quad
i \cong i \otimes 1, \quad
j \cong 1 \otimes i, \quad
ij \cong i \otimes i
$$

このように

- 因子に含まれる基底は、左右の位置を入れ替えると別のものになる
- 因子に含まれる実成分（係数）は可換で、括り出せる

という規則を持った積の演算子$\otimes$を導入することで、双複素数を構成することができます。このような$\otimes$による積を**テンソル積**と呼びます。[[wiki-tp]]

&&&rem 双線形性
このような性質を双線形性と呼びます。テンソル積と、2つの引数を取る双線形関数$B$を比較します。
$$
\begin{alignedat}{2}
ax⊗y\ &=& x⊗ay\ &= a(x⊗y) \\
B(ax,y) &=& B(x,ay) &= aB(x,y)
\end{alignedat}
$$
分配法則も双線形性に含まれます。
$$
\begin{alignedat}{2}
(ax+by)⊗z &=& ax⊗z\ +\ by⊗z\ &= a(x⊗z)+b(y⊗z) \\
B(ax+by,z) &=& B(ax,z)+B(by,z) &= aB(x,z)+bB(y,z)
\end{alignedat}
$$
&&&

テンソル積による構成によって「複素数$\otimes$複素数」という構造が明確となり、これが双複素数という名前の由来となっています。

&&&rem 双複素数とテンソル積の対応関係

- 双複素数での$i \neq j$は$(i \otimes 1) \neq (1 \otimes i)$に対応します。言い換えると、双複素数での基底の区別が、テンソル積の因子の位置の違いとして現れます。
- 双複素数では$ij=ji$ですが、テンソル積では両辺とも$i \otimes i$となるため区別がありません。言い換えると、双複素数での可換性が、テンソル積では同一な表現として現れます。

&&&

## テンソル積の積

双複素数の積を、テンソル積で書き直します。

$$
\begin{aligned}
ii&=-1 &(i \otimes 1)(i \otimes 1) &= -(1 \otimes 1) \\
jj&=-1 &(1 \otimes i)(1 \otimes i) &= -(1 \otimes 1) \\
(i)(j)&=ij &(i \otimes 1)(1 \otimes i) &= i \otimes i \\
(ij)(ij)&=1 &(i \otimes i)(i \otimes i) &= 1 \otimes 1
\end{aligned}
$$

対応関係の観察から、左因子は左因子と、右因子は右因子と積を計算すると定義します。

&&&def テンソル積の積
$$
(\alpha \otimes \beta)(\gamma \otimes \delta) = \alpha\gamma \otimes \beta\delta \quad(\alpha,\beta,\gamma,\delta\text{ は複素数})
$$
&&&

## 計算例

冒頭に挙げた双複素数での計算例をテンソル積で書き直します。

$$
\begin{aligned}
&(a+bi+cj)(d+ei+fj) \\
&\cong (a \otimes 1 + bi \otimes 1 + c \otimes i)(d \otimes 1 + ei \otimes 1 + f \otimes i) \\
&=  (a  \otimes 1)(d \otimes 1 + ei \otimes 1 + f \otimes i) \\
&\ +(bi \otimes 1)(d \otimes 1 + ei \otimes 1 + f \otimes i) \\
&\ +(c  \otimes i)(d \otimes 1 + ei \otimes 1 + f \otimes i) \\
&=   ad  \otimes 1 + aei   \otimes 1 + af  \otimes i \\
&\ + bdi \otimes 1 + bei^2 \otimes 1 + bfi \otimes i \\
&\ + cd  \otimes i + cei   \otimes i + cf  \otimes i^2 \\
&=   ad(1 \otimes 1) + ae(i \otimes 1) + af(1 \otimes i) \\
&\ + bd(i \otimes 1) - be(1 \otimes 1) + bf(i \otimes i) \\
&\ + cd(1 \otimes i) + ce(i \otimes i) - cf(1 \otimes 1) \\
&=(ad-be-cf)(1 \otimes 1) + (ae+bd)(i \otimes 1) + (af+cd)(1 \otimes i) + (bf+ce)(i \otimes i) \\
&\cong (ad-be-cf) + (ae+bd)i + (af+cd)j + (bf+ce)ij
\end{aligned}
$$

双複素数の計算結果と一致しました。

# まとめ

本記事では、複素数に可換な虚数単位を追加して得られる双複素数を通じて、零因子の性質およびテンソル積代数の自然な構成法を解説しました。

&&& 双複素数とテンソル積の同型
複素数体$\mathbb{C}$に$i$と可換な虚数単位$j\ (j^2=-1)$を付加した双複素数は、テンソル積$\mathbb{C} \otimes \mathbb{C}$と同型です。
$$
a + bi + cj + dij \cong a(1 \otimes 1) + b(i \otimes 1) + c(1 \otimes i) + d(i \otimes i)
$$
積の演算規則は以下の通りです。
$$
(\alpha \otimes \beta)(\gamma \otimes \delta) = \alpha\gamma \otimes \beta\delta
$$
&&&

&&& 零因子
双複素数には$0$でない元の積が$0$となる零因子が存在し、体の構造を持ちません。
$$
(1 + ij)(1 - ij) = 1 - (ij)^2 = 1 - (-1)(-1) = 0
$$
&&&

テンソル積によって代数系をさらに拡張することができます。

- 任意個の複素数のテンソル積 → セグレの多重複素数[[wiki-smc]][[7shi-mc]]
- 複素数と四元数のテンソル積 → 双四元数（パウリ行列が生成する代数と同型）[[wiki-q]][[wiki-bq]][[7shi-bq]]
- 複素数と八元数のテンソル積 → 双八元数[[wiki-o]][[wiki-bo]]

テンソル積の重要な応用として、多量子ビット状態を記述する量子コンピューターがあります。[[7shi-qc1]][[7shi-qc2]]
