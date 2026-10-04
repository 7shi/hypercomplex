テレンス・タオ教授の解説記事をもとに、四元数の代数構造と球面上の回転作用を通じて球面三角法の諸公式を導出します。

&&& 改訂履歴
- 2026.10.04 マスター方程式を定理とし、右辺が$-1$になることの証明を追加、単位四元数の呼称を「回転子」に変更
&&&

&&&rem
元記事から、積の並び順・三角形の向き・右辺の符号を変更しています（マスター方程式の後のremで説明します）。結論には影響しません。
&&&

# 概要

ハミルトンの四元数は複素数の非可換な拡張であり、単位四元数を使うと3次元空間の回転を積で表せます。これにより、球面上の幾何学法則を代数的に導けます。本記事では、四元数の基本性質や随伴性から出発し、単位四元数による3次元回転作用の合成を通じて「球面三角法のマスター方程式」を導出します。さらに、この方程式を変形し、内積と共役作用の成分を比較して、球面余弦定理・正弦定理・5要素の規則を導きます。その応用として、天文学における日の出方程式も導きます。[[tao-h]]

複素数、ベクトルの内積、三角関数を前提とし、四元数の基本事項は本文で説明します。単位四元数の群と$\mathrm{SU}(2)$の同型と、共役作用が位相的な意味での被覆写像であることは既知の結果として扱い、証明には立ち入りません。リー代数との関係は補足として触れるに留めます。

# 四元数の基礎

四元数は$t + xi + yj + zk$の形の数です。ここで$t,x,y,z$は実数であり、$i,j,k$は$-1$の平方根で、以下の関係を満たします。

$$
i^2 = j^2 = k^2 = -1, \quad ij=k, \quad jk=i, \quad ki=j
$$

これらは非可換（$ij \neq ji$）ですが、結合法則は成り立ちます。四元数全体は斜体（可除環）を形成し、すべての非ゼロ四元数は逆元を持ちます。

## 共役と分解

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

- 複素数 $\mathrm{Im}(1+2i) = 2$
- 四元数 $\mathrm{Im}(1+2i+3j+4k) = 2i+3j+4k$
&&&

## 内積とノルム

内積は次のように定義されます。

$$
\langle q, r \rangle := \mathrm{Re} (q \overline{r})
$$

これは実双線形・対称・正定値であり、$1,i,j,k$は正規直交基底を形成します。ノルムは次のように与えられます。

$$
|q| = \sqrt{\langle q,q \rangle} = \sqrt{t^2 + x^2 + y^2 + z^2}
$$

ノルムは乗法的な性質$|qr| = |q| |r|$を持ちます。

単位四元数の集合$\mathrm{U}(1,\mathbb{H}) = \{ q \in \mathbb{H}: |q|=1\}$は群をなし、$\mathrm{SU}(2)$と同型であることが知られています。単位四元数$q$については$q\overline{q}=|q|^2=1$より$\overline{q}=q^{-1}$です。

## 随伴性

内積には以下の随伴性が成り立ちます。これは線形代数における随伴作用素の概念に対応します。

$$
\langle qr, s \rangle = \langle q, s\overline{r} \rangle, \quad
\langle rq, s \rangle = \langle q, \overline{r}s \rangle
$$

証明では、$\mathrm{Re}(ab)=\mathrm{Re}(ba)$を使います。積$ab$と$ba$の違いは虚部の外積に当たる項の符号だけで、実部は成分表示から一致することが確かめられます。

&&&prf 随伴性の証明
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

## オイラーの公式

$i,j,k$は$-1$の平方根であるため、実数$\theta$に対して以下のオイラーの公式が成り立ちます。

$$
e^{i\theta} = \cos \theta + i \sin \theta, \quad e^{j\theta} = \cos \theta + j \sin \theta, \quad e^{k\theta} = \cos \theta + k \sin \theta
$$

# 球面上の回転作用

単位四元数$q \in \mathrm{U}(1,\mathbb{H})$は、純虚四元数（$\mathbb{R}^3$と同一視）に対して共役によって作用します。

$$
v \mapsto q v \overline{q}
$$

この作用が$\mathbb{R}^3$の内積と向きを保ち、回転を与えることを確かめます。

&&&prop 共役作用による回転
単位四元数$q$に対して、写像$v \mapsto qv\overline{q}$は純虚四元数の3次元空間の回転（向きを保つ直交変換）です。
&&&

&&&prf
$v$を純虚四元数とすると$\overline{v}=-v$であり、共役の反準同型性より次が成り立つ。

$$
\overline{qv\overline{q}} = q\overline{v}\,\overline{q} = -qv\overline{q}
$$

よって$qv\overline{q}$も純虚である。写像は実線形で、ノルムの乗法性より$|qv\overline{q}|=|q||v||\overline{q}|=|v|$となり、ノルムを保つ。内積は$\langle v,w\rangle = \frac12(|v+w|^2-|v|^2-|w|^2)$とノルムで表せるから内積も保たれ、この写像は直交変換である。その行列式は$\pm1$のいずれかで、$q$について連続である。単位四元数全体は4次元空間の単位球面であり連結なので、行列式は一定で、$q=1$での値$1$に等しい。よって向きも保たれる。
&&&

&&&rem 共役作用の用語
群論における**共役**（conjugation）は通常$gxg^{-1}$の形を指します。単位四元数においては$\overline{q} = q^{-1}$が成り立つため、$qv\overline{q}$は$qvq^{-1}$と等価であり、共役作用と呼ばれます。なお、四元数$\overline{q}$自体も「共役（四元数）」と呼ばれるため、文脈による区別が必要です。
&&&

## 二重被覆と半角

この作用において、$q$と$-q$は同じ回転を引き起こします。

$$
(-q)v\overline{(-q)} = qv\overline{q}
$$

&&&rem 符号の相殺
直感的には、$q$とその共役$\overline{q}$とで両側から挟むことで、符号の違いがキャンセルされるためだと解釈できます。
&&&

異なる単位四元数が同じ回転を与えることを、群論の言葉で整理します。

&&&def 忠実な作用と核
群$G$の集合$X$への作用が**忠実** (faithful) であるとは、異なる群元が異なる変換を引き起こすこと、つまり作用が単射であることをいいます。また、恒等変換として作用する群元の集合を作用の**核** (kernel) と呼びます。
&&&

作用が忠実であることは、核が単位元のみ（$\{1\}$）であることと同値です。

$q=1$と$q=-1$による共役作用は、どちらも恒等変換（何もしない変換）となります。逆に、共役作用が恒等変換になるのは、$qv=vq$がすべての純虚四元数$v$について成り立つときです。$i,j,k$のすべてと可換な四元数は実数に限られるため、単位長なら$q=\pm1$です。したがって、この作用の核は非自明な元$-1$を含む$\{1, -1\}$であり、作用は忠実ではありません。$q_1,q_2$が同じ回転を与えるなら$q_2^{-1}q_1$が核に入るため、同じ回転を与えるのは常に$\pm q$の組に限られます。また、後で見るように任意の回転は$e^{\theta n/2}$の形で実現されるため、対応は全射です。以上から、$\mathrm{SU}(2)$との同型を通じて、この対応は$\mathrm{SU}(2)$から$\mathrm{SO}(3)$への2対1の全射となります。これが位相的な意味での被覆写像（**二重被覆**）であることは、既知の結果として認めます。

例えば、$i$軸周りの$\theta$回転は、$e^{i\theta/2}$による共役で表されます。

$$
\begin{aligned}
e^{i\theta/2} i e^{-i\theta/2} &= i \\
e^{i\theta/2} j e^{-i\theta/2} &= (\cos\theta) j + (\sin\theta) k \\
e^{i\theta/2} k e^{-i\theta/2} &= (\cos\theta) k - (\sin\theta) j
\end{aligned}
$$

&&&prf $i$軸周りの回転の証明
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

上の計算が示す半角の関係は、リー代数の生成子の正規化にも現れます。$e^{ti}ve^{-ti}$を$t=0$で微分すると、次のようになります。

$$
\left.\frac{d}{dt}\right|_{t=0} e^{ti} v e^{-ti} = iv - vi = [i, v]
$$

$[i, j] = 2k$、$[i, k] = -2j$なので、$i$は$i$軸周りの角速度$2$の回転を生成し、角速度$1$の回転を生成するのは$i/2$です。3次元回転群$\mathrm{SO}(3)$のリー代数$\mathfrak{so}(3)$の標準的な基底$L_x, L_y, L_z$は交換関係$[L_x, L_y] = L_z$（およびその巡回置換）を満たし、四元数の側では$[i/2, j/2] = k/2$が対応します。この正規化された基底$i/2$を使えば、角度$\theta$の回転は$e^{\theta(i/2)} = e^{i\theta/2}$と書けます。

&&&rem リー代数と半角の意味
$\theta/2$という「半分の角度」は、回転の生成子が$i$ではなく$i/2$である（基底の方が半分である）と見ることもできます。交換関係における係数$2$のずれが、指数写像における$1/2$の因子として現れています。

直感的には、共役作用$q v \overline{q}$が左右両側から半分ずつ作用すると解釈することもできます。ただし、左右の乗算はそれぞれ単独では純虚四元数を純虚四元数に移さないため、3次元の回転を2回合成するという意味ではありません。$q=e^{i\theta/2}$の場合、軸に垂直な$j,k$の平面では左右の効果が加わって角度$\theta$の回転となり、軸方向の$i$では打ち消し合って不変となります。
&&&

## 回転子と4π周期性

リー群とリー代数の対応において、生成子$X$に対応する回転変換は、指数写像$e^{\theta X}$で表されます。四元数の場合、これは$\mathrm{SU}(2)$の元に当たる単位四元数であり、回転を表すものとして**回転子**（rotor）と呼びます。回転子は、回転角を連続的に増やしていくと、1回転（$360^\circ$）では元に戻らず、2回転（$720^\circ$）して初めて元に戻るという性質を持っています。（後で具体例を見ます）

生成子$X$は回転軸の方向を表す純虚四元数です。例えば、$X=i/2$は$i$軸周りの回転を表します。軸の方向は任意に選べますが、$\theta$を回転角とするため、$X$の大きさは$1/2$に定めます。

回転軸$n$をノルム$1$の単位純虚四元数として表す場合、回転角度を$\theta$として$X=n/2$であり、回転子は$e^{\theta n/2}$となります。回転子は軸の3成分と回転角で表せますが、軸には単位長の制約があるため、独立な自由度は3です。$n^2=-1$より、指数関数において$n$は虚数単位$i,j,k$と同様に振る舞います。

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

回転子は$e^{\theta n/2} = \cos(\theta/2) + \sin(\theta/2) n$と表されることから、三角関数の引数が回転角の半分になります。そのため、回転角を連続的に$\theta=2\pi$まで増やすと、回転子は$-1$に達します。

$$
e^{2\pi n/2} = e^{\pi n} = -1
$$

2回転$\theta=4\pi$によって$1$に戻ります。

$$
e^{4\pi n/2} = e^{2\pi n} = 1
$$

これらは、既に見たように$1$と$-1$による共役作用が恒等変換であることを示しています。

&&&rem 物理的アナロジー
2回転で元に戻る回転は、数学的な概念に留まらず、腕やベルトなどで物理的に再現できます。[[wiki-trick]]
&&&

## 回転の合成と順序

回転子の積は、ベクトルに対する回転（共役作用）の合成に対応しています。あるベクトル$v$に対して、まず回転$q_1$を適用し、次に回転$q_2$を適用する場合を考えます。

1. $v \mapsto q_1 v \overline{q_1}$
2. $(q_1 v \overline{q_1}) \mapsto q_2 (q_1 v \overline{q_1}) \overline{q_2} = (q_2 q_1) v \overline{(q_2 q_1)}$

このように、回転の合成は回転子の積$q_2 q_1$として表されます。$v$を包むように操作が重なるため、$q_2 q_1$のように右から左に並べる必要があります。

# 球面三角形とマスター方程式

単位球面$S^2 \subset \mathbb{R}^3$上の球面三角形を考えます。以下では、ある開半球内にあり、各辺を短い大円弧とする非退化な球面三角形を扱います。辺の長さと内角は$0<a,b,c<\pi$、$0<\alpha,\beta,\gamma<\pi$を満たします。

空間に固定された直交座標系$(i,j,k)$と、原点を中心に自由に回転できる球体を考えます。球体とは独立して固定された針が設置されており、球体上の一点$(1,0,0)$（$i$軸上の点）を指します。

&&&rem 回転スタンドのモデル
針は球体の回転には連動しません。球体が回転スタンドに置かれており、針は回転スタンドに取り付けられているイメージです。
&&&

回転子による回転操作は以下のように対応します。回転の向きは、軸の正の側から原点を見たときの向きで表します。

* $e^{i\theta/2}$：球体を$i$軸（針の軸）周りに反時計回りで$\theta$回転させます。針が指す点は変わらず、球面の向きだけが変わります。（方向転換に相当）
* $e^{k\theta/2}$：球体を$k$軸周りに反時計回りで$\theta$回転させます。球面上の点を横方向に移動させるのに使用します。（自転に相当）

## 三角形のトレース

単位球面上には、頂点$A, B, C$を持ち、辺の長さ（中心角）が$c, a, b$、頂点の角が$\alpha, \beta, \gamma$である球面三角形が描かれています。操作の符号を揃えるため、球面を外側から見て頂点が時計回りに並ぶ向きを採用します。辺の長さと内角は鏡映で変わらないため、最終的に得られる関係式はこの向きの選択によりません。

&&&rem 単位球面の辺長と中心角
単位球面では辺の長さと中心角（ラジアン）が一致します。また、頂点が反時計回りに配置されている場合は、$i$軸周りの回転を逆向きにすれば同じ操作ができます。
&&&

球体を回転させて、各頂点を順番に針の位置まで持ってくる操作を考えます。

1. 針が点$A$を指し、点$B$が$(\cos c)i-(\sin c)j$にある状態からスタートします。点$B$は針から見て$-j$の方向にあります。
2. $e^{kc/2}$：球体を$k$軸周りに回転させ、距離$c$だけ離れた点$B$を針の位置まで持ってきます。
3. $e^{i(\pi-\beta)/2}$：球体を$i$軸周りに$\pi-\beta$回転させ、次の点$C$が針から見て$-j$の方向に来るようにします。
4. $e^{ka/2}$：球体を$k$軸周りに回転させ、距離$a$だけ離れた点$C$を針の位置まで持ってきます。
5. $e^{i(\pi-\gamma)/2}$：球体を$i$軸周りに$\pi-\gamma$回転させ、次の点$A$が針から見て$-j$の方向に来るようにします。
6. $e^{kb/2}$：球体を$k$軸周りに回転させ、距離$b$だけ離れた点$A$を針の位置まで持ってきます。
7. $e^{i(\pi-\alpha)/2}$：球体を$i$軸周りに$\pi-\alpha$回転させ、最初の向き（点$B$が$(\cos c)i-(\sin c)j$にある状態）に戻します。

三角形の内角を$\alpha, \beta, \gamma$としたとき、頂点で方向転換する角度は$\pi - \alpha$、$\pi - \beta$、$\pi - \gamma$となります。方向転換は三角形の外側で行われるため、回転角度は内角の**補角**（外角）となります。

## 球面三角法のマスター方程式

三角形を一周して元の位置と向きに戻る操作を、「球体を回転させる操作の積」として回転子で表現します。操作は右から左へ並べるため、トレースの7つの操作の積について次が成り立ちます。

&&&thm 球面三角法のマスター方程式
$$
e^{i(\pi-\alpha)/2} e^{kb/2} e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2} e^{kc/2} = -1 \tag{4}
$$
&&&

これは球面三角形の辺$a,b,c$と角$\alpha,\beta,\gamma$が満たすべき制約を課しています。この式から各種公式が導出できるため、本記事では球面三角法の「マスター方程式」と呼びます。

証明では、積が$\pm1$のどちらかに限られることを示し、三角形を縮める連続変形で符号を決めます。

&&&prf
左辺を$Q$とおく。トレースの最後に、点$A$は再び針の位置$i$に戻り、点$B$も最初の位置$(\cos c)i-(\sin c)j$に戻る。$0<c<\pi$より、この2つのベクトルは一次独立である。$Q$による共役作用は、一次独立な2つのベクトルを固定する回転なので、それらの外積の方向も固定し、恒等変換である。作用の核は$\{1,-1\}$なので、$Q=\pm1$である。

次に、三角形を連続的に変形する。三角形を含む開半球の中心で球面に接する平面へ、球の中心から射影する（心射図法）。大円は直線に写るため、球面三角形は平面三角形に写る。平面上で三角形を内部の一点に向けて縮小し、球面へ戻せば、上の範囲の条件と向きを保ったまま、三角形を連続的に縮められる。この間、辺と角は連続に変化するから、$Q$も連続に変化する。$Q$は離散集合$\{1,-1\}$に値を取るから、変形の間一定である。

最後に、縮小の極限を考える。辺の長さは$0$に近づくから、$e^{ka/2},e^{kb/2},e^{kc/2}$は$1$に近づく。三角形は縮小の中心の近くに集まり、内角はその点での接平面上の平面三角形の内角に近づくから、内角の和は$\pi$に近づく。したがって次のようになる。

$$
Q \to e^{i(\pi-\alpha)/2} e^{i(\pi-\gamma)/2} e^{i(\pi-\beta)/2}
= e^{i(3\pi-(\alpha+\beta+\gamma))/2} \to e^{i\pi} = -1
$$

$Q$は一定であったから、$Q=-1$である。
&&&

&&&rem 右辺の符号の解釈
三角形を一周すると、球体の向きは元に戻ります。回転子は回転角を連続的に$2\pi$まで増やすと$-1$に達するため、右辺の$-1$は、方向転換の外角の和として$2\pi$の回転が蓄積されたものと見ることができます。ただし、軸の異なる回転の角は単純に足せないため、これは平面に近い微小な三角形での見方であり、証明は上の連続性の議論によります。
&&&

&&&rem 式の並び順と元記事との相違点
前述のように回転子は操作順に右から左へと並べます（例：$q_1$の後に$q_2$を適用する場合は$q_2 q_1$）。元記事では、点$A$に置いた接ベクトルを進めたり回したりする操作として因子を並べており、並びが本記事と逆で、三角形も反時計回りに配置しています。本記事では、ベクトルと同じように共役作用によって球体を回転させるモデルに変更しました。また、元記事では右辺を$1$としていますが、上の証明のとおり縮小の極限では左辺が$-1$に近づくため、本記事では$-1$としています。以下の公式の導出では、式(5)と(6)の内積を取るときや$\overline{F}iF$を作るときに右辺の符号が打ち消し合うため、この符号は結果に影響しません。
&&&

証明で用いた縮小の極限を、1次の項まで計算してみます。

&&&ex 無限小三角形（ユークリッド極限）
辺の長さ$a,b,c$を$\varepsilon a, \varepsilon b, \varepsilon c$に置き換え、$\varepsilon$について1次の項までテイラー展開します。$e^{k\varepsilon x/2} \approx 1 + k\varepsilon x/2$を用います。

$$
\begin{aligned}
&e^{i(\pi-\alpha)/2} \left(1+\frac{k\varepsilon b}{2}\right) e^{i(\pi-\gamma)/2} \left(1+\frac{k\varepsilon a}{2}\right) e^{i(\pi-\beta)/2} \left(1+\frac{k\varepsilon c}{2}\right) \\
= &e^{i(3\pi - (\alpha+\beta+\gamma))/2}
   \left(1 + \frac{k\varepsilon b}{2} e^{i(\pi-\gamma)} e^{i(\pi-\beta)} + \frac{k\varepsilon a}{2} e^{i(\pi-\beta)} + \frac{k\varepsilon c}{2} \right) + O(\varepsilon^2) \\
\end{aligned}
$$

縮小した球面三角形の内角は、平面三角形の内角に収束します。ここでは$\alpha,\beta,\gamma$として平面三角形の内角を使い、1次の項の様子を見ます。平面三角形では内角の和が$\alpha+\beta+\gamma = \pi$となることを利用します。$\varepsilon$について1次の項まで残します。

$$
\begin{aligned}
&(-1) \left(1 + \frac{k\varepsilon b}{2} e^{i(\pi-\gamma)} e^{i(\pi-\beta)} + \frac{k\varepsilon a}{2} e^{i(\pi-\beta)} + \frac{k\varepsilon c}{2} \right) \\
= &-1 - \frac{k\varepsilon}{2}(b e^{i(\pi-\gamma)} e^{i(\pi-\beta)} + a e^{i(\pi-\beta)} + c)
\end{aligned}
$$

第2項の括弧の中は、複素平面で表現された3辺のベクトルの和であるため、$0$となります。（これは無限小の領域において、回転が並進として扱えることを示唆します）

よって、三角形を縮めると式(4)の左辺は$-1$に収束し、1次の項も打ち消し合います。
&&&

&&&ex 三直角の正三角形（球面の8分の1）
球面の8分の1を占める球面三角形を考えます。

$$
a = b = c = \frac{\pi}{2}, \quad \alpha = \beta = \gamma = \frac{\pi}{2}
$$

このとき各項は以下のようになります。

$$
e^{k\pi/4} = \frac{1}{\sqrt{2}}(1+k), \quad e^{i(\pi-\pi/2)/2} = e^{i\pi/4} = \frac{1}{\sqrt{2}}(1+i)
$$

これらを掛け合わせます。

$$
e^{i\pi/4} e^{k\pi/4} = \frac{1}{2}(1+i)(1+k) = \frac{1}{2}(1+i-j+k)
$$

これは正規化された純虚四元数$u = \frac{1}{\sqrt{3}}(i-j+k)$を用いて$e^{u\pi/3}$と書けます。

$$
e^{u\pi/3}
= \cos \frac{\pi}{3} + \sin \frac{\pi}{3}u
= \frac{1}{2} + \frac{1}{2}(i-j+k)
$$

式 (4) の左辺はこれの3乗となります。

$$
(e^{u\pi/3})^3 = e^{u\pi} = -1
$$

式 (4) の右辺どおり$-1$となることが確認できます。
&&&

# 球面余弦定理の導出

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

**左辺の計算**は以下の通りです。
随伴性$\langle qr, s \rangle = \langle q, s\overline{r} \rangle,\ \langle rq, s \rangle = \langle q, \overline{r}s \rangle$を用いて両端の因子を消去します。

$$
\begin{aligned}
&\langle e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2}, e^{i(\pi-\gamma)/2} e^{-ka/2} e^{i(\pi-\beta)/2} \rangle \\
= &\langle e^{ka/2}, e^{-ka/2} \rangle
= \mathrm{Re}(e^{ka/2} e^{ka/2})
= \mathrm{Re}(e^{ka}) = \cos a
\end{aligned}
$$

**右辺の計算**は以下の通りです。
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

# 球面正弦定理の導出

次に、式 (5) の左辺を$F$、右辺を$G$とおき、$F=G$から得られる関係$\overline{F}iF = \overline{G}iG$を利用します。この両辺と$k$との内積をとります。

$$
\langle \overline{F}iF, k \rangle = \langle \overline{G}iG, k \rangle
$$

これを計算することで、球面正弦定理が得られます。

&&&thm 球面正弦定理
$$
\frac{\sin \alpha}{\sin a} = \frac{\sin \beta}{\sin b} = \frac{\sin \gamma}{\sin c}
$$
&&&

&&&prf 球面正弦定理の証明
$$
\begin{aligned}
\langle \overline{F}iF, k \rangle
&= \langle e^{-i(\pi-\beta)/2} e^{-ka/2} e^{-i(\pi-\gamma)/2} i e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2}, k \rangle \\
&= \langle e^{-ka/2} i e^{ka/2}, e^{i(\pi-\beta)/2} k e^{-i(\pi-\beta)/2} \rangle \\
&= \langle (\cos a) i - (\sin a) j, -(\cos \beta) k - (\sin \beta) j \rangle \\
&= \sin a \sin \beta \\

\langle \overline{G}iG, k \rangle
&= \langle e^{kc/2} e^{i(\pi-\alpha)/2} e^{kb/2} i e^{-kb/2} e^{-i(\pi-\alpha)/2} e^{-kc/2}, k \rangle \\
&= \langle e^{kb/2} i e^{-kb/2}, e^{-i(\pi-\alpha)/2} k e^{i(\pi-\alpha)/2} \rangle \\
&= \langle (\cos b) i + (\sin b) j, -(\cos \alpha) k + (\sin \alpha) j \rangle \\
&= \sin b \sin \alpha
\end{aligned}
$$
$$
\begin{aligned}
\langle \overline{F}iF, k \rangle &= \langle \overline{G}iG, k \rangle \\
\sin a \sin \beta &= \sin b \sin \alpha \\
\frac{\sin \beta}{\sin b} &= \frac{\sin \alpha}{\sin a}
\end{aligned}
$$
&&&

証明では最初の等号を示しました。頂点の名前を巡回的に付け替えて同じ計算を行えば、$\dfrac{\sin \beta}{\sin b} = \dfrac{\sin \gamma}{\sin c}$も得られます。

こちらも辺$a,b,c$の無限小極限によって平面上の正弦定理に収束します。

$$
\frac{\sin \alpha}{a} = \frac{\sin \beta}{b} = \frac{\sin \gamma}{c}
$$

# 5要素の規則

同様に、$\overline{F}iF = \overline{G}iG$の両辺と$j$との内積をとります。

$$
\begin{aligned}
\langle \overline{F}iF, j \rangle &= \langle \overline{G}iG, j \rangle \\
\end{aligned}
$$

これを計算することで、5要素の規則が得られます。

&&&fml 5要素の規則 (five-part rules)
$$
\sin a \cos \beta = \cos b \sin c - \sin b \cos c \cos \alpha
$$
&&&

&&&prf 5要素の規則の証明
$$
\begin{aligned}
\langle \overline{F}iF, j \rangle
&= \langle e^{-i(\pi-\beta)/2} e^{-ka/2} e^{-i(\pi-\gamma)/2} i e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2}, j \rangle \\
&= \langle e^{-ka/2} i e^{ka/2}, e^{i(\pi-\beta)/2} j e^{-i(\pi-\beta)/2} \rangle \\
&= \langle (\cos a) i - (\sin a) j, -(\cos \beta) j + (\sin \beta) k \rangle \\
&= \sin a \cos \beta \\

\langle \overline{G}iG, j \rangle
&= \langle e^{kc/2} e^{i(\pi-\alpha)/2} e^{kb/2} i e^{-kb/2} e^{-i(\pi-\alpha)/2} e^{-kc/2}, j \rangle \\
&= \langle e^{i(\pi-\alpha)/2} e^{kb/2} i e^{-kb/2} e^{-i(\pi-\alpha)/2}, e^{-kc/2} j e^{kc/2} \rangle \\
&= \langle (\cos b) i - (\sin b)(\cos \alpha) j + (\sin b)(\sin \alpha) k, (\cos c) j + (\sin c) i \rangle \\
&= \cos b \sin c - \sin b \cos \alpha \cos c
\end{aligned}
$$
$$
\begin{aligned}
\langle \overline{F}iF, j \rangle &= \langle \overline{G}iG, j \rangle \\
\sin a \cos \beta &= \cos b \sin c - \sin b \cos \alpha \cos c
\end{aligned}
$$
&&&

# 応用：日の出方程式

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

これを用いて、地球上の緯度$\phi$における日の出時刻を求めることができます。まず、用いる天文学の用語を定義します。

&&&def 赤緯と時角
- **赤緯**：天の赤道（地球の赤道を天球に投影した大円）から天体がどれだけ北（または南）に離れているかを示す角度です。
- **時角**：天体が南中（子午線の上方での通過）してから、天の北極を軸に西向きにどれだけ回転したかを示す角度です。南中時に$0^\circ$となり、西（午後）へ行くとプラス、東（午前）へ行くとマイナスになります。
&&&

太陽の赤緯は季節によって変化し、夏至には最大$+23.5^\circ$、冬至には最小$-23.5^\circ$となり、春分・秋分には$0^\circ$となります。

&&&rem 白夜と極夜
北半球において、太陽の赤緯$\delta$が正のときは、自転軸の北側が太陽の方向に傾いているため、北極点$(\phi=\pi/2)$から$\delta$の範囲内$(\phi > \pi/2-\delta)$で太陽は沈みません（**白夜**）。$\delta$が負のときは、北極点から$-\delta$の範囲内$(\phi > \pi/2+\delta)$で太陽は昇りません（**極夜**）。
&&&

## 設定

北半球の$0<\phi<\pi/2$を考えます。大気による屈折、太陽の視半径、地形は無視し、太陽の中心が地平線（高度$0$）を通過する時刻を日の出・日の入りとします。天球上で以下の直角三角形$\triangle PSN$を考えます。

* 頂点$P$：天の北極
* 頂点$S$：日の出時刻の地平線上の太陽
* 頂点$N$：地平線上の真北の点
* 直角$\angle N = \pi/2$（子午線と地平線の交角）

各辺と角は以下のようになります。

* 斜辺$b$：天の北極から太陽までの角度。太陽の赤緯を$\delta$とすると、$b = PS = \pi/2 - \delta$。
* 垂直辺$c$：地平線から天の北極までの角度（天の北極の高度）。観測者の緯度$\phi$に等しい。$c = PN = \phi$。
* 角$\alpha$：天の北極における角度。日の出の時角の大きさを$\omega$とすると、真北基準の角は$\alpha = \angle P = \pi - \omega$。

&&&rem 時角の符号と大きさ
日の出の時角は負の値$-\omega$ですが、ここでは三角形の幾何学的性質を扱うため、その大きさ$\omega$を用います。
&&&

## 導出

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

$\delta=0$では$b=\pi/2$となり$\tan b$が定義されませんが、(7)に変形する前の式$0 = \cos b \sin c - \sin b \cos c \cos \alpha$に代入すれば$\cos\omega=0$となり、同じ結果が得られます。

$\omega$について解くと、次のようになります。

$$
\omega = \arccos(-\tan \phi \tan \delta)
$$

$\arccos$の定義域から、この式が使えるのは$|\tan \phi \tan \delta| \le 1$の場合です。$\tan \phi \tan \delta > 1$では白夜、$\tan \phi \tan \delta < -1$では極夜となり、日の出・日の入りはありません。$\tan \phi \tan \delta = \pm 1$は、太陽が地平線に接する境界です。

$\omega$は日の出から南中までの時角の変化量を表すため、日の出から日の入りまでの時角の変化量は$2\omega$です。時角は24時間で$2\pi$進むため、昼の長さは$24 \times \frac{2\omega}{2\pi} = \frac{24}{\pi} \omega$時間となります。

$|\tan \phi \tan \delta| \le 1$の範囲で$\arccos(-x) = \pi - \arccos x$の関係を用いると、昼の長さは次のように表されます。

$$
\begin{aligned}
\frac{24}{\pi} \omega
&= \frac{24}{\pi}\arccos(-\tan \phi \tan \delta) \\
&= \frac{24}{\pi}(\pi - \arccos(\tan \phi \tan \delta)) \\
&= 24 - \frac{24}{\pi}\arccos(\tan \phi \tan \delta) \\
\end{aligned}
$$

&&&ex 季節による昼の長さの変化
北半球$(\phi > 0)$において、以下の関係が成り立ちます。

* 春分・秋分$(\delta = 0)$：$\cos \omega = 0$より$\omega = \pi/2$（12時間）となり、昼夜の長さが等しくなります。
* 夏$(\delta > 0)$：$\cos \omega < 0$より$\omega > \pi/2$となり、昼が長くなります。
* 冬$(\delta < 0)$：$\cos \omega > 0$より$\omega < \pi/2$となり、昼が短くなります。

夏$(\delta > 0)$に白夜の境界$(\phi = \pi/2 - \delta)$では、次のようになります。

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

したがって、昼の長さは$24$時間となります。
&&&

# まとめ

本記事では、単位四元数による3次元回転作用を用いて、球面三角形を一周して元に戻る操作を回転子の積として表し、球面三角法のマスター方程式を導きました。一周後の共役作用は恒等変換なので積は$\pm1$に限られ、三角形を縮める連続変形により、その値は$1$ではなく$-1$と定まります。

マスター方程式を変形することで球面余弦定理が得られ、変形した式の両辺$F=G$から作った$\overline{F}iF=\overline{G}iG$と$k,j$との内積から、球面正弦定理と5要素の規則が得られました。さらに5要素の規則を直角三角形に適用し、日の出の時角を与える日の出方程式を導きました。

&&& 球面三角法のマスター方程式
頂点$A,B,C$の角を$\alpha,\beta,\gamma$、それぞれの対辺の長さを$a,b,c$とすると、開半球内の非退化な球面三角形について次が成り立ちます。
$$
e^{i(\pi-\alpha)/2} e^{kb/2} e^{i(\pi-\gamma)/2} e^{ka/2} e^{i(\pi-\beta)/2} e^{kc/2} = -1
$$
&&&

&&& 球面余弦定理
$$
\cos a = \cos b \cos c + \sin b \sin c \cos \alpha
$$
&&&

&&& 球面正弦定理
$$
\frac{\sin \alpha}{\sin a} = \frac{\sin \beta}{\sin b} = \frac{\sin \gamma}{\sin c}
$$
&&&

&&& 5要素の規則
$$
\sin a \cos \beta = \cos b \sin c - \sin b \cos c \cos \alpha
$$
&&&

&&& 日の出方程式
北半球の緯度$\phi$、赤緯$\delta$における日の出の時角の大きさ$\omega$は、次の式で与えられます。
$$
\cos \omega = -\tan \phi \tan \delta
$$
&&&
