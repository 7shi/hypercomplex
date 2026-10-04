複素数の虚数単位を増やして拡張する試みから、双複素数の代数構造と零因子の性質を整理し、テンソル積の自然な導入を解説します。

&&& 改訂履歴
- 2026.10.04 テンソル積を左右の複素数を区別して記録する記号として導入し直し、実数上のテンソル積$\mathbb C\otimes_{\mathbb R}\mathbb C$の定義と双複素数との同型の命題と証明を追加
&&&

# 概要

複素数にもう1つ虚数単位を加えると、どのような数の体系になるでしょうか。新しい虚数単位を元の虚数単位と可換にすると、4次元の双複素数（bicomplex number）の体系が得られます。双複素数には零因子があるため、$0$でないすべての元が逆数を持つ体系にはなりません。一方で、双複素数は2つの複素数を左右に並べて組み合わせるテンソル積$\mathbb{C} \otimes \mathbb{C}$と同型になります。

前提とするのは複素数の計算だけです。本記事では、虚数単位の導入から双複素数の乗算規則と零因子の発生機構を確認したうえで、左右の複素数を区別して記録する記号$\otimes$を導入し、左右それぞれの因子で積を取るという規則によってテンソル積代数が構成される過程を具体例で確認します。本記事のテンソル積は、すべて実数体$\mathbb{R}$上のものです。テンソル積の一般的な定義には立ち入らず、$\mathbb C\otimes_{\mathbb R}\mathbb C$を基底によって具体的に構成します。

# 虚数単位を増やす

複素数の虚数単位$i$について$i^2=-1$が成り立ちます。ここに、同じ性質$j^2=-1$を満たす新しい虚数単位$j$を導入します。ただし、$i \ne \pm j$で、$i$と$j$は可換$ij=ji$とします。また、通常の数と同じく結合法則と分配法則が成り立ち、実数はすべての元と可換であるとします。

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

$ij$を含む元どうしの積でも、$ij=ji$によって因子を並べ替え、$i^2=j^2=-1$を使えば、$1,i,j,ij$の実数係数の和に整理できます。したがって、$i,j$から和と積によって得られる元は、次の形で表せます。

$$
a+bi+cj+dij \quad (a,b,c,d\text{ は実数})
$$

この表し方は一意です。

&&&prop 表し方の一意性 [prop-unique]
実数$a,b,c,d$について$a+bi+cj+dij=0$ならば、$a=b=c=d=0$です。
&&&

&&&prf
$a+bi+cj+dij=(a+bi)+(c+di)j=0$とする。$c+di\ne0$なら、複素数$c+di$の逆数を使って$j$について解くと、次のようになる。
$$
j=-\frac{a+bi}{c+di}
$$
右辺は複素数なので、$j$は複素数であることになる。複素数の中で2乗が$-1$になるのは$\pm i$だけなので$j=\pm i$となり、前提$i\ne\pm j$に反する。したがって$c=d=0$であり、$a+bi=0$から$a=b=0$である。
&&&

&&&def 双複素数
$i^2=j^2=-1$、$ij=ji$、$i\ne\pm j$を満たす$i,j$を用いて、実数$a,b,c,d$により$a+bi+cj+dij$と表される数の体系を**双複素数** (bicomplex number) と呼びます。[[wiki-bc]][[wiki-bc-en]]
&&&

&&&rem 三元数・四元数との違い

- [[prop-unique]]より、$ij$は$a+bi+cj$の形では表せません。したがって双複素数は3次元の三元数ではなく、実4次元の体系です。
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

ここで$1+ij$と$1-ij$は、どちらも$0$ではありません。$1+ij=0$なら$ij=-1$で、両辺に左から$i$を掛けると$-j=-i$、つまり$j=i$となります。同様に$1-ij=0$なら$j=-i$となります。どちらも前提$i\ne\pm j$に反します。

&&&def 零因子
$x\ne0$に対して、ある$y\ne0$が存在して$xy=0$となるとき、$x$を**零因子**と呼びます。[[wiki-0d]]
&&&

双複素数は可換なので、左右の区別は要りません。$1+ij$と$1-ij$は、どちらも零因子です。

零因子は逆数を持ちません。なぜなら、$1+ij$に何かを掛けて$1$になると仮定すれば、その両辺に$1-ij$を掛けることで矛盾が生じるためです。

$$
\begin{aligned}
(1+ij)x&=1 \\
(1-ij)(1+ij)x&=1-ij \\
0&=1-ij
\end{aligned}
$$

これは$1-ij$が$0$でないことに反します。

# テンソル積による構成

左右の複素数を区別して記録する記号$\otimes$を導入します。これは、2つの複素数を通常どおり掛け合わせる記号ではありません。各因子について分配法則が成り立つものとして、次の式を展開します。（$a,b,c,d$は実数）

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
= ac\cdot1 + ad\,j + bc\,i + bd\,ij
$$

$\otimes$でも、実数の係数は括り出せるという計算規則を追加して、係数と基底を分離します。

$$
\begin{aligned}
&= (a\cdot1) \otimes (c\cdot1) + (a\cdot1) \otimes (d\,i) + (b\,i) \otimes (c\cdot1) + (b\,i) \otimes (d\,i) \\
&= ac(1 \otimes 1) + ad(1 \otimes i) + bc(i \otimes 1) + bd(i \otimes i)
\end{aligned}
$$

2つの展開を比べると、双複素数の4つの基底$1,j,i,ij$と、$\otimes$による4つの元$1\otimes1,\ 1\otimes i,\ i\otimes1,\ i\otimes i$が同じ係数で現れています。そこで、この4つの元を基底とする空間を考えます。

&&&def 実数上のテンソル積
4つの元$1\otimes1,\ i\otimes1,\ 1\otimes i,\ i\otimes i$を独立な基底とする実ベクトル空間を$\mathbb C\otimes_{\mathbb R}\mathbb C$と書きます。複素数$a+bi,\ c+di$（$a,b,c,d$は実数）の**テンソル積**を、その元として次のように定めます。
$$
(a+bi)\otimes(c+di)=ac(1\otimes1)+ad(1\otimes i)+bc(i\otimes1)+bd(i\otimes i)
$$
&&&

基底が独立であるとは、この4つを使った実数係数の表し方が一意になることです。一般の元は4つの基底の実数係数の和であり、必ずしも1つの$\alpha\otimes\beta$の形にはなりません。以降、$\otimes_{\mathbb R}$を単に$\otimes$と書きます。[[wiki-tp]]

この定義から、$\otimes$は次の規則を持ちます。

- 因子に含まれる基底は、左右の位置を入れ替えると別のものになる（$i\otimes1\ne1\otimes i$）
- 因子に含まれる実数の係数は可換で、括り出せる

外に出せる係数は実数だけです。$i$は係数ではないので、$i\otimes1$の$i$を右の因子に移すことはできません。

各因子について分配法則が成り立ち、実数の係数を外に出せる性質を**双線形性**と呼びます。

&&&fml 双線形性
複素数$x,y,z$と実数$a,b$について、次が成り立ちます。
$$
\begin{aligned}
(ax+by)\otimes z &= a(x\otimes z)+b(y\otimes z) \\
x\otimes(ay+bz) &= a(x\otimes y)+b(x\otimes z)
\end{aligned}
$$
特に、$ax\otimes y=x\otimes ay=a(x\otimes y)$です。
&&&

&&&rem 双線形関数との比較
2つの引数を取る関数$B$が、$B(ax+by,z)=aB(x,z)+bB(y,z)$と$B(x,ay+bz)=aB(x,y)+bB(x,z)$を満たすとき、$B$を双線形関数と呼びます。テンソル積は、この性質を記号$\otimes$の計算規則として持たせたものです。

ただし、双線形性だけではテンソル積は定まりません。たとえば複素数の通常の掛け算も実数について双線形ですが、$i\cdot1=1\cdot i$となり左右を区別しません。テンソル積を定めているのは、4つの元を独立な基底とした点です。
&&&

双複素数の基底を次のように対応させると、$(a+bi)(c+dj)$の展開と$(a+bi)\otimes(c+di)$の展開は同じ係数を持ちます。

$$
1 \mapsto 1 \otimes 1, \quad
i \mapsto i \otimes 1, \quad
j \mapsto 1 \otimes i, \quad
ij \mapsto i \otimes i
$$

双複素数での基底の区別$i \neq j$は、テンソル積では因子の位置の違い$i \otimes 1 \neq 1 \otimes i$として現れます。この段階で対応しているのは、実ベクトル空間としての基底だけです。積の対応は次の節で確かめます。

テンソル積による構成から、双複素数を2つの複素数の構造を組み合わせたもの（「複素数$\otimes$複素数」）として捉えられます。

## テンソル積の積

双複素数の積が対応先でも保たれるためには、テンソル積の側で次が成り立つ必要があります。

$$
\begin{aligned}
ii&=-1 &(i \otimes 1)(i \otimes 1) &= -(1 \otimes 1) \\
jj&=-1 &(1 \otimes i)(1 \otimes i) &= -(1 \otimes 1) \\
(i)(j)&=ij &(i \otimes 1)(1 \otimes i) &= i \otimes i \\
(ij)(ij)&=1 &(i \otimes i)(i \otimes i) &= 1 \otimes 1
\end{aligned}
$$

左因子は左因子と、右因子は右因子と積を取れば、これらはすべて満たされます。そこで次のように定義します。

&&&def テンソル積の積
$\alpha\otimes\beta$の形の元どうしの積を次の式で定め、一般の元どうしの積は分配法則と実数係数に関する双線形性によって定めます。
$$
(\alpha \otimes \beta)(\gamma \otimes \delta) = \alpha\gamma \otimes \beta\delta \quad(\alpha,\beta,\gamma,\delta\text{ は複素数})
$$
&&&

&&&rem 一般の代数のテンソル積
同じ規則は、複素数に限らず、実数上の結合的な代数$A,B$のテンソル積$A\otimes B$の積にも用います。その場合、各因子の中では掛ける順序を保ちます。
&&&

この規則のもとで、$i\otimes1$と$1\otimes i$の積は順序によりません。

$$
(i\otimes1)(1\otimes i)=i\otimes i=(1\otimes i)(i\otimes1)
$$

これは$i\otimes1$と$1\otimes i$が同じ元であるという意味ではありません。異なる元でも、積の順序を交換できるという意味です。双複素数での可換性$ij=ji$は、テンソル積では左右の因子がそれぞれ別に掛け合わされることとして現れます。一般に、異なる因子から来た元$\alpha\otimes1$と$1\otimes\beta$は可換です。

以上で、双複素数とテンソル積の対応が積も保つことを示せます。

&&&prop 双複素数とテンソル積の同型 [prop-iso]
双複素数$a+bi+cj+dij$を$a(1\otimes1)+b(i\otimes1)+c(1\otimes i)+d(i\otimes i)$に移す写像$\Phi$は、実代数としての同型です。
&&&

&&&prf
$\Phi$は基底を基底に移す実線形写像なので、全単射である。積を保つことを示す。双線形性により、基底どうしの積について確かめれば十分である。

基底は$m,n\in\{0,1\}$を用いて$i^mj^n$と書け、$\Phi(i^mj^n)=i^m\otimes i^n$である。この式は任意の$m,n\ge0$でも成り立つ。実際、$i^2=-1$より$i^m=s\,i^{m'}$（$m'\in\{0,1\}$、$s=\pm1$）と書け、$j^2=-1$より同じ符号で$j^m=s\,j^{m'}$となる。$n$についても同様に$i^n=t\,i^{n'}$、$j^n=t\,j^{n'}$とすれば、双線形性により次のようになる。
$$
\Phi(i^mj^n)=st\,\Phi(i^{m'}j^{n'})=st\,(i^{m'}\otimes i^{n'})=(s\,i^{m'})\otimes(t\,i^{n'})=i^m\otimes i^n
$$

したがって、$ij=ji$を用いて次が成り立つ。
$$
\Phi(i^mj^n)\,\Phi(i^pj^q)=(i^m\otimes i^n)(i^p\otimes i^q)=i^{m+p}\otimes i^{n+q}=\Phi(i^{m+p}j^{n+q})=\Phi(i^mj^n\,i^pj^q)
$$
よって$\Phi$は積を保ち、実代数としての同型である。
&&&

以降、同型で対応する元を$\cong$で結んで書きます。

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

本記事では、複素数に可換な虚数単位を追加して得られる双複素数を通じて、零因子の性質と、実数上のテンソル積による双複素数の構成を確認しました。

&&& 双複素数とテンソル積の同型
複素数体$\mathbb{C}$に$i$と可換な虚数単位$j\ (j^2=-1)$を付加した双複素数は、実数上のテンソル積$\mathbb{C} \otimes_{\mathbb R} \mathbb{C}$と実代数として同型です。
$$
a + bi + cj + dij \cong a(1 \otimes 1) + b(i \otimes 1) + c(1 \otimes i) + d(i \otimes i)
$$
積は左右の因子ごとに取ります。
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

以下は本記事では扱わない関連話題です。テンソル積によって代数系をさらに拡張することができます。

- 任意個の複素数のテンソル積 → セグレの多重複素数[[wiki-smc]][[7shi-mc]]
- 複素数と四元数のテンソル積 → 双四元数（パウリ行列が実代数として生成する代数と同型）[[wiki-q]][[wiki-bq]][[7shi-bq]]
- 複素数と八元数のテンソル積 → 双八元数[[wiki-o]][[wiki-bo]]

テンソル積の重要な応用として、多量子ビット状態を記述する量子コンピューターがあります。量子状態の記述では、通常は複素数体上のテンソル積を用います。[[7shi-qc1]][[7shi-qc2]]
