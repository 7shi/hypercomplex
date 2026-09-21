ホップファイブレーションは、複素数から四元数、八元数へと拡張することができます。複素数による$S^3\to S^2$の写像を基礎として、四元数による$S^7\to S^4$、および八元数による$S^{15}\to S^8$のホップ写像を構成します。

複素数・実数・四元数・八元数の各代数について、回転による構成の可否、ホップ写像の直接的な拡張、およびファイバーの構造を確認します。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

# 複素数

初回の記事において、単位四元数の分解から得られる複素数のペアによって、ホップファイブレーションの表式を導きました。[[7shi-h]] ここでは改めて複素数のペアを$\alpha,\beta\in\mathbb{C}$と書きます。

&&&def 複素数ホップ写像
単位球面上の点$|\alpha|^2 + |\beta|^2 = 1$（$\alpha,\beta\in\mathbb{C}$）は、以下の写像によって2次元球面上の点に射影されます。
$$
H(\alpha, \beta) = \left(2\alpha \beta^*, |\alpha|^2 - |\beta|^2\right)
$$
$$
S^1\hookrightarrow S^3\xrightarrow{H}S^2
$$
&&&

ここでの$\alpha \beta^*$は複素数（2つの実成分）であり、$|\alpha|^2 - |\beta|^2$は実数（1つの実成分）であるため、写像の像は3次元空間内のベクトルとなります。このベクトルの長さ（ノルム）が$1$になることを確認します。

&&&prop 球面への射影
$|\alpha|^2 + |\beta|^2 = 1$のとき、$|H(\alpha,\beta)| = 1$が成り立ちます。すなわち像は単位2次元球面$S^2$上にあります。
&&&

&&&prf
ノルムの乗法性$|\alpha\beta^*| = |\alpha||\beta^*| = |\alpha||\beta|$を用いると、以下のように計算できる。

$$
\begin{aligned}
|H(\alpha, \beta)|^2 &= |2\alpha \beta^*|^2 + (|\alpha|^2 - |\beta|^2)^2 \\
&= 4|\alpha|^2|\beta|^2 + |\alpha|^4 - 2|\alpha|^2|\beta|^2 + |\beta|^4 \\
&= |\alpha|^4 + 2|\alpha|^2|\beta|^2 + |\beta|^4 \\
&= (|\alpha|^2 + |\beta|^2)^2 \\
&= 1
\end{aligned}
$$
&&&

この計算で使ったのは、ノルムの乗法性$|xy| = |x||y|$と、その後の実数の計算だけです。乗法の結合法則は使っていません。

ノルムの乗法性は、複素数だけでなく実数$\mathbb{R}$、四元数$\mathbb{H}$、八元数$\mathbb{O}$においても成立します。

&&&rem フルヴィッツの定理
実数・複素数・四元数・八元数はいずれもノルムの乗法性を持ちます。逆にフルヴィッツの定理によれば、正定値の乗法的ノルムを持つ有限次元の実単位的代数は、この4種類に限られます。
&&&

したがって、写像$H(\alpha, \beta) = (2\alpha \beta^*, |\alpha|^2 - |\beta|^2)$は、成分を実数・四元数・八元数に置き換えた$H_{\mathbb{R}}, H_{\mathbb{H}}, H_{\mathbb{O}}$として、同じ式のまま球面から球面への写像を定義します。以下、代数ごとにそのファイバーを調べ、ホップファイブレーションが得られることを確認します。

# 実数

&&&def 実数ホップ写像
単位円周上の点$\alpha^2 + \beta^2 = 1$は、以下の写像によって1次元球面（円周）上の点に射影されます。
$$
H_{\mathbb{R}}(\alpha,\beta)=(2\alpha\beta,\ \alpha^2-\beta^2)
$$
$$
S^0\hookrightarrow S^1\xrightarrow{H_{\mathbb{R}}}S^1
$$
&&&

初回の記事では[[7shi-h]]、複素数ペア$(\alpha,\beta)$を1つ上の代数である四元数$\omega$に埋め込み、単位四元数による共役作用$\omega\mathbf{k}\omega^*$で純虚四元数を回転させることで$S^3\to S^2$を構成しました。同じ「1つ上の代数を使う」という発想が、実数の場合にどこまで通用するかを確認します。

まず$H_{\mathbb{R}}$の定義を確認します。実数には虚部がなく、共役は恒等写像となります（$x^*=x$）。そのため$H$の定義に現れる$\beta^*$は$\beta$のままとなり、$H_{\mathbb{R}}(\alpha,\beta)=(2\alpha\beta,\ \alpha^2-\beta^2)$が得られます。

次に、実数ペア$(\alpha,\beta)$を1つ上の代数である複素数$\zeta=\alpha+\beta i$とみなして、その2乗を計算します。

$$
\zeta^2=(\alpha+\beta i)^2=(\alpha^2-\beta^2)+2\alpha\beta i
$$

$H_{\mathbb{R}}$の像のペアは、左側$2\alpha\beta$を虚部の係数、右側$\alpha^2-\beta^2$を実部と見なせば、$\zeta^2$そのものです。つまり実数ホップ写像は、1つ上の代数$\mathbb{C}$における2乗として実現されます。

$\zeta=\cos\theta+i\sin\theta$とおけば$\zeta^2=\cos2\theta+i\sin2\theta$であり、角度を2倍にする写像、すなわち二重被覆$\operatorname{Spin}(2)\to\operatorname{SO}(2)$そのものとして実数ホップファイブレーションを構成できます。

なお、四元数の場合の共役作用に対応する$\zeta p\zeta^*$は、$\mathbb{C}$が可換であるため任意の$p$を$p$自身に写し、回転として機能しません（これは$\operatorname{SO}(1)$が自明群であることに対応します）。実数の場合に共役作用の役割を果たしているのは、乗算そのもの、すなわち2乗$\zeta^2$です。

## ホップ写像とファイバー

写像$H_{\mathbb{R}}(\alpha, \beta)$の像について成分を確認します。
- $2\alpha \beta$は実数（1つの実成分）
- $\alpha^2 - \beta^2$は実数（1つの実成分）

合計で2つの成分を持ち、ノルムは1となるため、像は1次元球面（円周）$S^1$上の点となります。これが$S^1 \to S^1$のホップファイブレーションであり、ファイバーは0次元球面（2点）$S^0$となります。

先に見たとおり$H_{\mathbb{R}}$は角度を2倍にする写像$z\mapsto z^2$、すなわち二重被覆$\operatorname{Spin}(2)\to\operatorname{SO}(2)$ですが、ファイバーが0次元（離散）であるため、この束は被覆写像そのものでもあります。各ファイバーは核$\{\pm1\}$の剰余類であり、いずれも2点からなる$S^0$です。

&&&rem 回転群と1つ上の代数
底空間の円$S^1$の回転群は$\operatorname{SO}(2)\cong\operatorname{Spin}(2)\cong\operatorname{U}(1)$（単位複素数群）であり、これは$\mathbb{R}$の1つ上の代数である$\mathbb{C}$のノルム$1$の元の集合そのものです。

これらはリー群として互いに同型ですが、標準的な被覆準同型$\operatorname{Spin}(2)\to\operatorname{SO}(2)$は、その同型写像ではありません。前者は群同士が同型であるという主張、後者は特定の準同型の性質です。両群を$\operatorname{U}(1)$と同一視すると、被覆準同型は$z\mapsto z^2$になります。
&&&

&&&rem 複素数以降との違い
ファイバーと核の対応は実数の場合に限られます。複素数以降ではファイバーが正の次元を持つため、$S^1\hookrightarrow S^3\to S^2$は被覆ではありません（$S^2$は単連結なので、その連結な被覆は自明なものに限られます）。複素数における二重被覆は$\operatorname{Spin}(3)\cong\operatorname{SU}(2)\to\operatorname{SO}(3)$という別の写像が担い[[7shi-lie2]]、ファイバー$S^1\cong\operatorname{U}(1)$のうち部分群$\{\pm1\}$だけがその核に対応します。代数を実数に替えると、ノルム$1$の元の集合は$\operatorname{U}(1)$ではなく$\{\pm1\}$になるため、両者が重なります。
&&&

# 四元数

&&&def 四元数ホップ写像
単位球面上の点$|\alpha|^2 + |\beta|^2 = 1$（$\alpha,\beta\in\mathbb{H}$）は、以下の写像によって4次元球面上の点に射影されます。
$$
H_{\mathbb{H}}(\alpha, \beta) = \left(2\alpha \beta^*, |\alpha|^2 - |\beta|^2\right)
$$
$$
S^3\hookrightarrow S^7\xrightarrow{H_{\mathbb{H}}}S^4
$$
&&&

複素数から四元数への構成が成立したのは、単位四元数による共役作用$\omega p \omega^*$が純虚四元数の空間（3次元）上でちょうど$\operatorname{SO}(3)$の回転として作用するという、$\operatorname{Spin}(3) \cong \operatorname{SU}(2) \cong \operatorname{Sp}(1)$という特別な一致によるものです。ここで単位四元数の群$\operatorname{Sp}(1)$は、回転群$\operatorname{SO}(3)$の二重被覆$\operatorname{Spin}(3)$と同一視されています。

同じ発想を1段階ずらして、四元数ペア$(\alpha,\beta)\in\mathbb{H}^2$を1つ上の八元数$\omega$に埋め込み、単位八元数の共役作用で純虚八元数を回転させることを考えます。しかし純虚八元数の空間は7次元であるため像は$S^6$となり、目的の$S^4$と次元が一致しません。

そこで本記事では、1つ上の代数の単一の元による共役作用という形の構成ではなく、初回で得られた成分表示の結果を直接拡張する方針を取ります。

## ホップ写像

写像$H_{\mathbb{H}}(\alpha, \beta)$の像について成分を確認します。
- $2\alpha \beta^*$は四元数（4つの実成分）
- $|\alpha|^2 - |\beta|^2$は実数（1つの実成分）

合計で5つの成分を持ち、前述の通りノルムは1となるため、像は4次元球面$S^4$上の点となります。

&&&rem 回転群との関係
底空間の球面$S^4$の回転群は$\operatorname{SO}(5)\cong\operatorname{Spin}(5)/\mathbb{Z}_2$であり、$\operatorname{Spin}(5)\cong \operatorname{Sp}(2)$（四元数成分の$2\times2$ユニタリ行列群）です。これは単一の四元数や八元数による共役作用ではなく、$\mathbb{H}^2$に作用する群です。

$\operatorname{Sp}(2)$は$S^7$に推移的に作用し、その点の安定化群は$\operatorname{Sp}(1)$です。同次空間としては
$$
S^7\cong\operatorname{Sp}(2)/\operatorname{Sp}(1),\quad
S^4\cong\operatorname{Sp}(2)/\bigl(\operatorname{Sp}(1)\times\operatorname{Sp}(1)\bigr)
$$
であり、その間の射影のファイバーは$\bigl(\operatorname{Sp}(1)\times\operatorname{Sp}(1)\bigr)/\operatorname{Sp}(1)\cong S^3$と表されます。$S^7\to S^4$と回転群との関係は、この$\operatorname{Spin}(5)\cong \operatorname{Sp}(2)$の作用を用いて記述できます。

なお、ここに現れる安定化群$\operatorname{Sp}(1)$（全空間の点を固定する群）と、次節でファイバーを動かす右からの共通乗算とは、別のものです。後者はこの束を主$\operatorname{Sp}(1)$束として記述しています。
&&&

## ファイバー（結合的）

ホップファイブレーションの特徴は、写像の逆像（ファイバー）自体が球面になることです。写像が同じ値を取る（同じ点に射影される）条件として、$(\alpha, \beta)$の右からノルム1の元$q$を掛けた点$(\alpha q, \beta q)$を考えます。

四元数では乗法の結合法則が成り立ちます。右から単位四元数$q \in S^3$を掛けた場合の像を計算します。

$$
\begin{aligned}
H_{\mathbb{H}}(\alpha q, \beta q) &= \left(2(\alpha q)(\beta q)^*, |\alpha q|^2 - |\beta q|^2\right) \\
&= \left(2(\alpha q)(q^* \beta^*), |\alpha|^2 - |\beta|^2\right) \\
&= \left(2\alpha (qq^*) \beta^*, |\alpha|^2 - |\beta|^2\right) \quad (\text{結合法則}) \\
&= \left(2\alpha \beta^*, |\alpha|^2 - |\beta|^2\right) \\
&= H_{\mathbb{H}}(\alpha, \beta)
\end{aligned}
$$

$(\alpha, \beta)$と$(\alpha q, \beta q)$は同じ点に射影されることがわかります。

逆に、同じ点に射影される点が必ずこの形に書けることも確認できます。$\alpha\neq0$として、$H_{\mathbb{H}}(\alpha',\beta')=H_{\mathbb{H}}(\alpha,\beta)$を満たす点に対して$q=\alpha^{-1}\alpha'$とおくと、第2成分の一致から$|\alpha'|=|\alpha|$すなわち$|q|=1$であり、第1成分の一致$\alpha'\beta'^*=\alpha\beta^*$から$\beta'=\beta q$が従います。$\alpha=0$の場合は$\beta$を使って同様に議論できます。

したがって一つの点の逆像は$(\alpha q,\beta q)$の全体であり、単位四元数$q$の集合が$S^3$を成すため、ファイバーは$S^3$となります。

# 八元数

&&&def 八元数ホップ写像
単位球面上の点$|\alpha|^2 + |\beta|^2 = 1$（$\alpha,\beta\in\mathbb{O}$）は、以下の写像によって8次元球面上の点に射影されます。
$$
H_{\mathbb{O}}(\alpha, \beta) = \left(2\alpha \beta^*, |\alpha|^2 - |\beta|^2\right)
$$
$$
S^7\hookrightarrow S^{15}\xrightarrow{H_{\mathbb{O}}}S^8
$$
&&&

八元数では乗法の結合法則が成り立たないため、前節で用いた共通右乗算によるファイバー保存の計算を、そのまま適用することはできません。そこで、直接拡張した写像とそのファイバーを確認します。

&&&rem 部分代数の結合性
八元数でも、2つの元が生成する部分代数は結合的です。そのため、初回のように単位純虚元$p$に対して$q=u+vp$と置く場合には$(qp)q^*=p$が成り立ちます。問題になるのは、任意の$\alpha,\beta,q$に対して$(\alpha q)(q^*\beta^*)$の括弧を$\alpha(qq^*)\beta^*$へ付け替える操作の方です。
&&&

## ホップ写像

写像$H_{\mathbb{O}}(\alpha, \beta)$の像について成分を確認します。
- $2\alpha \beta^*$は八元数（8つの実成分）
- $|\alpha|^2 - |\beta|^2$は実数（1つの実成分）

合計で9つの成分を持ち、ノルムはやはり1となるため、像は8次元球面$S^8$上の点となります。ノルムの計算にはノルムの乗法性しか使っていないため、非結合性は写像の定義には影響しません。

&&&rem 回転群との関係
底空間の球面$S^8$の回転群は$\operatorname{SO}(9)\cong\operatorname{Spin}(9)/\mathbb{Z}_2$です。$\operatorname{Spin}(9)$はスピン表現を通じて$S^{15}$に推移的に作用し、その点の安定化群は$\operatorname{Spin}(7)$（21次元）です。一方、底空間$S^8$の点の安定化群は$\operatorname{Spin}(8)$であり、同次空間としては
$$
S^{15}\cong\operatorname{Spin}(9)/\operatorname{Spin}(7),\quad
S^8\cong\operatorname{Spin}(9)/\operatorname{Spin}(8)
$$
です。ファイバーは$\operatorname{Spin}(8)/\operatorname{Spin}(7)\cong S^7$と表され、この表示から構造群を$\operatorname{Spin}(8)$とする球面束が得られます。

八元数の非結合性により単位八元数の集合$S^7$自体は乗法に関して群をなさないため、構造群はファイバー自身の乗法ではなく、$\operatorname{Spin}(8)$という別の群になります。
&&&

## ファイバー（非結合的）

四元数の場合と同様に、$(\alpha, \beta)$の右からノルム1の元$q$を掛けた点$(\alpha q, \beta q)$を考えます。八元数では乗法の結合法則が成り立たないため、右から単位八元数$q \in S^7$を掛けた場合の像は次のようになります。

$$
\begin{aligned}
H_{\mathbb{O}}(\alpha q, \beta q) &= \left(2(\alpha q)(\beta q)^*, |\alpha q|^2 - |\beta q|^2\right) \\
&= \left(2(\alpha q)(q^* \beta^*), |\alpha|^2 - |\beta|^2\right)
\end{aligned}
$$

ここで、八元数の非結合性により、一般には以下のような非等式が成り立ちます。

$$
(\alpha q)(q^* \beta^*) \neq \alpha (qq^*) \beta^* = \alpha \beta^*
$$

したがって、一般には$H_{\mathbb{O}}(\alpha q, \beta q) = H_{\mathbb{O}}(\alpha, \beta)$は成り立ちません。

&&&ex
三つ組$123,145,246,365$[[7shi-oct1]]を用いて、等式が成り立たないことを確認します。

$$
\alpha=e_1/\sqrt2,\quad \beta=e_4/\sqrt2,\quad q=e_2\quad(|q|=1)
$$
とおいて計算します。

$$
\begin{aligned}
H_{\mathbb{O}}(\alpha,\beta)
&= \left(2\alpha\beta^*, 0\right)
 = \left(-e_1e_4, 0\right)
 = (-e_5, 0)
\\
H_{\mathbb{O}}(\alpha q,\beta q)
&=H_{\mathbb{O}}\left(\frac{e_3}{\sqrt2},\frac{-e_6}{\sqrt2}\right)
 =\left(2\cdot\frac{e_3}{\sqrt2}\cdot\frac{e_6}{\sqrt2},0\right)
 = (e_5, 0)
\end{aligned}
$$

よって、$H_{\mathbb{O}}(\alpha,\beta) \neq H_{\mathbb{O}}(\alpha q, \beta q)$が確認できました。
&&&

これは、八元数ホップファイブレーションにおいて、単純な右からの乗算では同じファイバーをなぞることができないことを意味しています。写像自体は$S^8$への射影として問題なく定義でき、束としても成立しますが、ファイバー$S^7$自身を構造群とする主束としては記述できません。

## 恒等式によるファイバーの具体的構成

前述の通り、単位八元数の集合$S^7$は非結合性のため乗法に関して群をなさず、構造群にはなれませんが、代わりに$\operatorname{Spin}(8)$を構造群とする球面束として$S^7\hookrightarrow S^{15}\to S^8$は成立しています。ここでは群の言葉を使わず、底空間の点を固定したときのファイバーを、$\beta$を$\alpha$の**左からの乗算**で決める形で具体的に表示します。

用いるのは、結合則を使わず交代則から従う合成代数の基本公式です。

&&&fml 合成代数の恒等式
$$
x(x^*y) = |x|^2 y,\qquad x^*(xy) = |x|^2 y
$$
2番目の式は、1番目の式で$x$を$x^*$に置き換えると得られます。
&&&

これを使って、赤道上の点$(p,\ 0)\in S^8$（$p$は単位八元数、$|p|=1$）のファイバーを求めます。

&&&prop 赤道上のファイバー
単位八元数$p$を固定すると、点$(p,\ 0)\in S^8$の逆像は
$$
H_{\mathbb{O}}^{-1}(p,\ 0)=\left\{(\alpha,\ p^*\alpha) : \alpha\in\mathbb{O},\ |\alpha|=\frac{1}{\sqrt2}\right\}
$$
であり、半径$1/\sqrt2$の球面、すなわち$S^7$と同一視されます。
&&&

&&&prf
まず、この形の点が$(p,0)$に写ることを示す。$\beta^*=\alpha^*p$、すなわち$\beta=(\alpha^*p)^*=p^*\alpha$と定める。これは$\alpha$から$\beta$を決める固定の線形写像（左乗算$L_{p^*}$）であり、$(\alpha,\beta)\mapsto(\alpha q,\beta q)$のような共通右乗算とは異なる作り方である。

$|\alpha|=1/\sqrt2$、$|p|=1$のとき、上記の恒等式により

$$
\begin{aligned}
H_{\mathbb{O}}(\alpha,p^*\alpha)
&= (2\alpha(\alpha^*p),|\alpha|^2 - |p|^2|\alpha|^2) \\
&= (2|\alpha|^2 p,|\alpha|^2(1-|p|^2)) \\
&= (p,\ 0)
\end{aligned}
$$

となる。逆に、$H_{\mathbb{O}}(\alpha,\beta)=(p,0)$を満たす点を考える。第2成分から$|\alpha|=|\beta|=1/\sqrt2$であり、第1成分$\alpha\beta^*=p/2$の両辺に左から$\alpha^*$を掛けると、上記の恒等式の2番目の式により$|\alpha|^2\beta^*=\alpha^*p/2$、すなわち$\beta^*=\alpha^*p$が一意に定まる。よって逆像は上記の集合に一致する。
&&&

このように、ファイバー$S^7$は「$\alpha$をノルム$1/\sqrt2$の球面上で自由に選び、$\beta$を固定の左乗算$L_{p^*}$で決める」という形で具体的に構成できます。これは構造群の作用を直接記述するものではなく、底空間の点を固定したときのファイバーを、左乗算のグラフとして表示したものです。八元数では共通右乗算による方法を一般には適用できないため、こうした表示がファイバーを具体的に書き下す手段になります。

&&&prop 一般の点のファイバー
一般の点$(c,\ r)\in S^8$についても、$r\neq-1$であれば
$$
|\alpha|^2=\frac{1+r}{2},\quad \beta^*=\frac{\alpha^* c}{1+r}
$$
と定めた点の全体が逆像となります。$r=-1$（南極）の逆像は$\{(0,\beta):|\beta|=1\}$です。いずれも$S^7$と同一視されます。
&&&

&&&prf
$r\neq-1$とする。$|c|^2=1-r^2$であるから、上記の$\beta$のノルムは

$$
|\beta|^2=\frac{|\alpha|^2|c|^2}{(1+r)^2}=\frac{1+r}{2}\cdot\frac{1-r^2}{(1+r)^2}=\frac{1-r}{2}
$$

となり、$|\alpha|^2+|\beta|^2=1$を満たす。第2成分は$|\alpha|^2-|\beta|^2=r$であり、第1成分は赤道の場合と同じ上記の恒等式により

$$
2\alpha\beta^*=\frac{2\alpha(\alpha^*c)}{1+r}=\frac{2|\alpha|^2c}{1+r}=c
$$

となる。よって$H_{\mathbb{O}}(\alpha,\beta)=(c,\ r)$である。逆向きの一意性も赤道の場合と同様である。$r=-1$のときは$|\alpha|^2-|\beta|^2=-1$と$|\alpha|^2+|\beta|^2=1$から$\alpha=0$、$|\beta|=1$が定まり、このとき第1成分は自動的に$0$となる。
&&&

&&&rem 四元数との関係
この構成は恒等式$x(x^*y)=|x|^2y$のみに基づくため、結合則を使わず、四元数でも同様に成り立ちます。ただし四元数では、結合性により$p^*(\alpha_0 q)=(p^*\alpha_0)q=\beta_0 q$が成り立つため、左乗算による構成は右からの共通乗算$(\alpha_0 q,\beta_0 q)$（構造群＝ファイバー自身の乗法）と式レベルで一致し、別の構成ではなく同じものになります。八元数では非結合性により、一般には$p^*(\alpha_0 q)=(p^*\alpha_0)q$が成り立たないため、この一致は崩れます。
&&&

# まとめ

初回で用いた「1つ上の代数を使って回転させる」構成は、実数の場合は共役作用ではなく乗算（2乗$\zeta^2$、二重被覆$\operatorname{Spin}(2)\to\operatorname{SO}(2)$）によって成立する一方、$\operatorname{Spin}(3)\cong \operatorname{SU}(2)\cong \operatorname{Sp}(1)$という特別な一致に依る共役作用$\omega p \omega^*$の構成は、次元の不一致や非結合性のために$S^7\to S^4$、$S^{15}\to S^8$へはそのまま拡張できません。

そこで、ノルムの乗法性に依存するホップ写像$H(\alpha, \beta) = (2\alpha \beta^*, |\alpha|^2 - |\beta|^2)$を、実数版$H_{\mathbb{R}}$・四元数版$H_{\mathbb{H}}$・八元数版$H_{\mathbb{O}}$として直接拡張しました。実数・四元数では結合性により右からの共通乗算$(\alpha q, \beta q)$がそのままファイバーを保存し、それぞれ$S^0$・$S^3$を与えます。

一方、$H_{\mathbb{O}}$では非結合性の影響が現れ、右乗算によるファイバーの保存則が成り立たなくなることを具体例で確認しました。これは八元数ホップファイブレーションが、ファイバー自身の乗法を構造群とする主束ではなく、$\operatorname{Spin}(8)$を構造群とする球面束であることに対応しています。それでも各ファイバーは、恒等式$x(x^*y)=|x|^2y$に基づく固定の左乗算$\beta=p^*\alpha$のグラフとして具体的に書き下すことができ、$S^7$全体を得ることができました。
