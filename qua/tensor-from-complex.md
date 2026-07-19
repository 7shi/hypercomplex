複素数の虚数単位を増やして拡張することから、テンソル積を導入します。

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

の形で表せることが分かります。このような数の体系を**双複素数**[[wiki-bc]] (bicomplex number[[wiki-bc-en]]) と呼びます。

&&&rem

- $ij$は$a+bi+cj$の形では表せないことから、三元数ではありません。
- 双複素数は$ij=ji$であることから、$ij=-ji$となる四元数[[wiki-q]]とは異なります。

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

$1+ij,\ 1-ij$のように、$0$ではないのに積が$0$になる数を**零因子**[[wiki-0d]]と呼びます。

零因子は逆数を持ちません。なぜなら、$1+ij$に何かを掛けて$1$になると仮定すれば、その両辺に$1-ij$を掛けることで矛盾が生じるためです。

$$
\begin{aligned}
(1+ij)x&=1 \\
(1-ij)(1+ij)x&=1-ij \\
0&=1-ij
\end{aligned}
$$

これは$1-ij$が$0$でないという前提に矛盾します。

&&&rem
もう少し式変形を進めると$ij=1,\ j=-i$となって、前提$i \ne \pm j$に矛盾します。
&&&

# テンソル積による構成

非可換な積の演算子$\otimes$を導入して、複素数の積を計算します。（$a,b,c,d$は実数）

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

- 因子に含まれる基底は非可換
- 因子に含まれる実成分（係数）は可換で、括り出せる

という規則を持った積の演算子$\otimes$を導入することで、双複素数を構成することができます。このような$\otimes$による積を**テンソル積**[[wiki-tp]]と呼びます。

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

&&&rem

- 双複素数での$i \neq j$は$(i \otimes 1) \neq (1 \otimes i)$に対応します。言い換えると、双複素数での基底の区別が、テンソル積の因子の非可換性として現れます。
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

テンソル積によって代数系を拡張することができます。いくつか例を挙げます。

- 任意個の複素数のテンソル積 → セグレの多重複素数[[wiki-smc]][[7shi-mc]]
- 複素数と四元数[[wiki-q]]のテンソル積 → 双四元数[[wiki-bq]]（パウリ行列が生成する代数と同型[[7shi-bq]]）
- 複素数と八元数[[wiki-o]]のテンソル積 → 双八元数[[wiki-bo]]

テンソル積の応用として量子コンピューターがあります。[[7shi-qc1]][[7shi-qc2]]
