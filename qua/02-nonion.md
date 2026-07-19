シルベスターの「九元数 (nonion)」の論考について、原文の流れを保ちながら、記法と説明を現代化して再構成します。[[nonion]][[sylvester]][[kindai1]]

## $2$次行列から四元数が現れる

シルベスターは、まず$2$次の行列$u, v$をとり、

$$
\det(zI + yv + xu)
$$

を$x, y, z$の$2$次式として見るところから始めます。その形を

$$
z^2 + 2bxz + 2cyz + dx^2 + 2exy + fy^2
$$

と書きます。

&&& 補足：2×2行列の成分計算
$u, v$を成分で

$$
u = \begin{pmatrix} u_1 & u_2 \\ u_3 & u_4 \end{pmatrix}, \quad
v = \begin{pmatrix} v_1 & v_2 \\ v_3 & v_4 \end{pmatrix}
$$

と書けば、

$$
M = zI + yv + xu = \begin{pmatrix} z & 0 \\ 0 & z \end{pmatrix} + y \begin{pmatrix} v_1 & v_2 \\ v_3 & v_4 \end{pmatrix} + x \begin{pmatrix} u_1 & u_2 \\ u_3 & u_4 \end{pmatrix}
$$

であり、まとめると

$$
M = \begin{pmatrix} z + v_1 y + u_1 x & v_2 y + u_2 x \\ v_3 y + u_3 x & z + v_4 y + u_4 x \end{pmatrix}.
$$

この行列式を展開すると

$$
\begin{aligned}
&z^2 + (u_1+u_4)xz + (v_1+v_4)yz + (u_1u_4-u_2u_3)x^2 \\
&+ (u_1v_4+v_1u_4-u_2v_3-v_2u_3)xy + (v_1v_4-v_2v_3)y^2
\end{aligned}
$$

が得られ、これを整理したものがシルベスターの

$$
z^2 + 2bxz + 2cyz + dx^2 + 2exy + fy^2
$$

です。
&&&

&&&rem なぜ交差項の係数に$2$が付くのか
シルベスターが$2$次形式の形で整理することを意図していたためだと考えられます。

$$
\begin{aligned}
&\begin{pmatrix} z & x & y \end{pmatrix}
\begin{pmatrix} 1 & b & c \\ b & d & e \\ c & e & f \end{pmatrix}
\begin{pmatrix} z \\ x \\ y \end{pmatrix} \\
&= z^2 + 2bzx + 2czy + dx^2 + 2exy + fy^2
\end{aligned}
$$

対称行列の非対角成分は展開時に$2$回ずつ現れるので、$xz, yz, xy$の係数がそれぞれ$2b, 2c, 2e$になります。
&&&

$$
vu + uv = 0
$$

が成り立つための必要十分条件は

$$
b = 0, \quad c = 0, \quad e = 0
$$

です。つまり、$\det(zI + yv + xu)$の交差項が消えるための条件が、$u, v$の反交換関係を導くということです。

さらに$d = 1,\ f = 1$を課し、$uv = w$とおくと、

$$
u^2 = -1, \quad v^2 = -1, \quad w^2 = -1
$$
$$
uv = -vu = w, \quad vw = -wv = u, \quad wu = -uw = v
$$

となります。したがって$1, u, v, w$は四元数系をなします。つまり

$$
\det(zI + yv + xu) = z^2 + y^2 + x^2
$$

という最も対称的な形に整えようとすると、四元数の関係式が自動的に現れるわけです。

&&&prf $uv=-vu,\ u^2=v^2=-I$
$$
\det(zI + yv + xu) = z^2 + 2bxz + 2cyz + dx^2 + 2exy + fy^2
$$

が$z^2 + y^2 + x^2$となる条件は

$$
b=c=e=0,\ d=f=1.
$$

$b=0$は$u_1+u_4=0$、$c=0$は$v_1+v_4=0$なので

$$
u_4=-u_1, \quad v_4=-v_1.
$$

$e=0$は$u_1v_4+v_1u_4-u_2v_3-v_2u_3=0$なので

$$
2u_1v_1+u_2v_3+u_3v_2=0.
$$

このとき

$$
\begin{aligned}
uv
&=\begin{pmatrix} u_1 & u_2 \\ u_3 & -u_1 \end{pmatrix}
  \begin{pmatrix} v_1 & v_2 \\ v_3 & -v_1 \end{pmatrix}
=\begin{pmatrix}
u_1v_1+u_2v_3 & u_1v_2-u_2v_1 \\
u_3v_1-u_1v_3 & u_3v_2+u_1v_1
\end{pmatrix} \\
vu
&=\begin{pmatrix} v_1 & v_2 \\ v_3 & -v_1 \end{pmatrix}
  \begin{pmatrix} u_1 & u_2 \\ u_3 & -u_1 \end{pmatrix}
=\begin{pmatrix}
v_1u_1+v_2u_3 & v_1u_2-v_2u_1 \\
v_3u_1-v_1u_3 & v_3u_2+v_1u_1
\end{pmatrix}
\end{aligned}
$$

$$
uv+vu
=\begin{pmatrix}
2u_1v_1+u_2v_3+u_3v_2 & 0 \\
0 & 2u_1v_1+u_2v_3+u_3v_2
\end{pmatrix}
=0
$$

$$
\therefore uv=-vu.
$$

次に、$d=1$は$u_1u_4-u_2u_3=1$なので、$u_4=-u_1$より

$$
-u_1^2-u_2u_3=1 \quad \Rightarrow \quad
u_1^2+u_2u_3=-1
$$

$$
u^2
=\begin{pmatrix} u_1 & u_2 \\ u_3 & -u_1 \end{pmatrix}^2
=\begin{pmatrix} u_1^2+u_2u_3 & u_1u_2 - u_2u_1 \\ u_3u_1 - u_1u_3 & u_1^2+u_2u_3 \end{pmatrix}
=\begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}
$$

$$
\therefore u^2=-I.
$$

同様に、$f=1$は$v_1v_4-v_2v_3=1$なので、$v_4=-v_1$より

$$
-v_1^2-v_2v_3=1 \quad \Rightarrow \quad
v_1^2+v_2v_3=-1
$$

$$
v^2
=\begin{pmatrix} v_1 & v_2 \\ v_3 & -v_1 \end{pmatrix}^2
=\begin{pmatrix} v_1^2+v_2v_3 & v_1v_2 - v_2v_1 \\ v_3v_1 - v_1v_3 & v_1^2+v_2v_3 \end{pmatrix}
=\begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}
$$

$$
\therefore v^2=-I.
$$
&&&

## 具体例

この条件は、たとえば

$$
v = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}, \quad
u = \begin{pmatrix} 0 & \theta \\ \theta & 0 \end{pmatrix} \quad
(\theta = \sqrt{-1})
$$

とすれば満たされます。このとき

$$
z + yv + xu = \begin{pmatrix} z & y + x\theta \\ -y + x\theta & z \end{pmatrix}
$$

したがって、次の$4$つの行列

$$
\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix},
\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix},
\begin{pmatrix} 0 & \theta \\ \theta & 0 \end{pmatrix},
\begin{pmatrix} -\theta & 0 \\ 0 & \theta \end{pmatrix}$$

は、実数係数のもとで通常の四元数系$1, i, j, k$の行列表現になっています。

&&&rem この具体例をどう選ぶか

本文ではシルベスターの与えた例をそのまま示しましたが、どうやって選んだかを考察します。

$b=c=0$に対応する「トレースが$0$」という条件は、対角成分を$0$にすれば満たせます。

$$
v = \begin{pmatrix} 0 & v_2 \\ v_3 & 0 \end{pmatrix}, \quad
u = \begin{pmatrix} 0 & u_2 \\ u_3 & 0 \end{pmatrix}
$$

行列式を$1$にするには、非対角成分の積が$-1$になる必要があります。もっとも単純な組み合わせとして$1(-1)=\theta^2=-1$を選びます。

$$
v = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}, \quad
u = \begin{pmatrix} 0 & \theta \\ \theta & 0 \end{pmatrix}
$$

なお、$vu$を計算するとシルベスターの示した行列とは符号が逆になります。

$$
vu = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}
     \begin{pmatrix} 0 & \theta \\ \theta & 0 \end{pmatrix}
   = \begin{pmatrix} \theta & 0 \\ 0 & -\theta \end{pmatrix}
$$
&&&

&&&rem なぜ$\theta$を入れるのか

ここで$\theta$を導入する理由は、$z^2 + y^2 + x^2$という正定値の対称な形を保ちたいからです。符号が混ざる二次形式ではなく、すべてがプラスの「球面的なノルム」を行列式として実現しようとすると、虚数を組み込んだこの表現が自然に現れます。
&&&

&&& 補足：パウリ行列との比較
ここで現れる2×2行列表示はパウリ行列を連想させます。実際、四元数は2×2複素行列で表現できますが、標準的なパウリ行列はそのままだと$2$乗して$+I$になります。四元数のように$2$乗して$-I$にしたいときは、虚数単位を掛けた形にする必要があります。
&&&

## 行列単位による記述

シルベスターは、現代の線形代数では$E_{00}, E_{01}, E_{10}, E_{11}$と表記される2×2行列の行列単位を$\lambda, \mu, \nu, \tau$と書き、その積を

$$
\begin{array}{c|cccc}
 & \lambda & \mu & \nu & \tau \\
\hline
\lambda & \lambda & \mu & 0 & 0 \\
\mu & 0 & 0 & \lambda & \mu \\
\nu & \nu & \tau & 0 & 0 \\
\tau & 0 & 0 & \nu & \tau
\end{array}
$$

と表します。そして

$$
\lambda + \tau = 1, \quad -\mu + \nu = i, \quad
-\theta \lambda + \theta \tau = k, \quad \theta \mu + \theta \nu = j
$$

とおくことで、四元数基底が組み立てられることを示します。原文ではこのあと、$u, v$の$8$個の成分に$5$条件がかかるので自由度は$3$つ残る、と数え上げています。さらに四元数基底$(i, j, k)$を直交変換した$(i', j', k')$に置き換える自由度もちょうど$3$つあるため、上の具体例から一般形へ移れる、という見通しを述べています。

## 四元数から九元数へ

ここからが本題です。$2$次の場合の

$$
vu = -uv, \quad u^2 = -1, \quad v^2 = -1
$$

が

$$
\det(zI + yu + xv) = z^2 + y^2 + x^2
$$

と対応したのとまったく同じ発想で、$3$次行列$u, v$と$1$の立方根$\rho=e^{2\pi i/3}$に対して

$$
vu = \rho uv, \quad u^3 = 1, \quad v^3 = 1
$$

を満たすための必要十分条件は

$$
\det(zI + yu + xv) = z^3 + y^3 + x^3
$$

である、とシルベスターは主張します。現代的には、$2$次のノルム形式から四元数が出るのと同様に、$3$次の対称なノルム形式から「九元数 (nonion)」と呼ぶべき$9$次元の非可換系が現れる、という見方ができます。

&&&rem なぜ$2$次では$-1$で、$3$次では$1$なのか
ケイリー・ハミルトンの定理を念頭に置くと見通しが立ちます。

$2$次の場合、$x=1,\ y=0$とすると

$$
\det(zI+u)=z^2+1
$$

であり、$\det$内が$0$になる条件より$z=-u$を右辺に**形式的に**代入すれば

$$
(-u)^2+1=0 \quad \Rightarrow \quad u^2=-1
$$

が得られます。いっぽう$3$次では

$$
\det(zI+u)=z^3+1
$$

なので、

$$
(-u)^3+1=0 \quad \Rightarrow \quad u^3=1.
$$

偶数次では負号が打ち消され、奇数次では負号が残るため、$2$次では$-1$、$3$次では$1$が現れるわけです。
&&&

## 九元数の具体的な生成元

シルベスターは具体例として

$$
u = \begin{pmatrix} 0 & 0 & 1 \\ \rho & 0 & 0 \\ 0 & \rho^2 & 0 \end{pmatrix}, \quad
v = \begin{pmatrix} 0 & 0 & 1 \\ \rho^2 & 0 & 0 \\ 0 & \rho & 0 \end{pmatrix}
$$

を与えます。このとき

$$
z + yu + xv
= \begin{pmatrix} z & 0 & y + x \\ \rho y + \rho^2 x & z & 0 \\ 0 & \rho^2 y + \rho x & z \end{pmatrix}
$$

であり、その行列式は

$$
z^3 + (y + x)(\rho y + \rho^2 x)(\rho^2 y + \rho x) = z^3 + y^3 + x^3.
$$

この最後の等式は、$\rho^3 = 1$と$1 + \rho + \rho^2 = 0$を使えば確認できます。

## 九元数の$9$個の基底

したがって、四元数で$1, u, v, uv$が基底になるのと同じく、ここでは$9$個の基底からなる系が得られます。

$$
\begin{alignedat}{6}
I &= \begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{pmatrix}&&,&
v &= \begin{pmatrix}
0 & 0 & 1 \\
\rho^2 & 0 & 0 \\
0 & \rho & 0
\end{pmatrix}&&,&
v^2 &= \begin{pmatrix}
0 & \rho & 0 \\
0 & 0 & \rho^2 \\
1 & 0 & 0
\end{pmatrix}&&, \\
u &= \begin{pmatrix}
0 & 0 & 1 \\
\rho & 0 & 0 \\
0 & \rho^2 & 0
\end{pmatrix}&&,&
uv &= \begin{pmatrix}
0 & \rho & 0 \\
0 & 0 & \rho \\
\rho & 0 & 0
\end{pmatrix}&&,&
uv^2 &= \begin{pmatrix}
1 & 0 & 0 \\
0 & \rho^2 & 0 \\
0 & 0 & \rho
\end{pmatrix}&&, \\
u^2 &= \begin{pmatrix}
0 & \rho^2 & 0 \\
0 & 0 & \rho \\
1 & 0 & 0
\end{pmatrix}&&,&\quad
u^2v &= \begin{pmatrix}
\rho & 0 & 0 \\
0 & \rho^2 & 0 \\
0 & 0 & 1
\end{pmatrix}&&,&\quad
u^2v^2 &= \begin{pmatrix}
0 & 0 & \rho \\
\rho & 0 & 0 \\
0 & \rho & 0
\end{pmatrix}&&.
\end{alignedat}
$$

この$9$個にさらに$1, \rho, \rho^2$を掛けた$27$個の元は閉じた群を構成します。四元数で$\pm 1, \pm u, \pm v, \pm uv$が閉じることに対応します。

&&&rem $\rho^2uv$とシフト行列
現代の言葉では、これは$3$次の一般化パウリ行列（クロック行列・シフト行列）と実質的に同じものです。[[wiki-gpm]]

たとえば

$$
uv=\begin{pmatrix}0&\rho&0\\0&0&\rho\\\rho&0&0\end{pmatrix}
$$

に$\rho^2$を掛けると$3$次のシフト行列が得られます。

$$
\rho^2uv
= \rho^2\begin{pmatrix}0&\rho&0\\0&0&\rho\\\rho&0&0\end{pmatrix}
= \begin{pmatrix}0&1&0\\0&0&1\\1&0&0\end{pmatrix}
$$
&&&

&&& 補足：現代数学・物理とのつながり
この種の関係式

$$
u^3=1,\quad v^3=1,\quad vu=\rho uv
$$

は、量子情報理論で用いられる$3$準位系 (qutrit) の演算子や、一般化クリフォード代数の基本例として現れます。$2$次ではパウリ行列、$3$次以上ではその一般化という形で、シルベスターの洞察は現代の線形代数・数理物理にそのままつながっています。[[wiki-gca]]
&&&

## 一般形の自由度

シルベスターはさらに、$3$次行列$u, v$には$18$個の定数があり、

$$\det(zI + uy + vx) = z^3 + y^3 + x^3$$

という条件は$9$個の制約を与えるので、九元数系の一般形には$18 - 9 = 9$個の任意定数が含まれるはずだ、と述べます。ただし、その$9$自由度を上の具体例からどう明示的に取り出すかは「今後の検討に委ねる」としています。

## パースへの言及と訂正

原文には

> These forms can be derived from an algebra given by Mr Charles S. Peirce (Logic of Relatives, 1870).

という一文がありますが、後注でシルベスター自身がこれを訂正しています。要点は、「自分がその事実を直接知っているわけではなく、C. S. Peirce 本人からそう聞いたのでそう書いた」ということです。つまり、パースの Logic of Relatives から本当にこの形が導けるかどうかは、シルベスター自身が確認したわけではない、という訂正です。

同じ注の中で彼は有名な一言も残しています。

> I have also a great repugnance to being made to speak of Algebras in the plural; I would as lief acknowledge a plurality of Gods as of Algebras.

「代数を複数形で語るのはどうにも抵抗がある。神々の複数性を認めるよりも、代数の複数性を認める方がまだ難しい」といった調子です。代数体系が次々に増殖していた時代への、いかにもシルベスターらしい反応です。

## 論理と代数

最後にシルベスターは「代数を論理に応用するのはすでに知られた話だが、その逆、つまり論理を代数に応用することにはもっと高い知的可能性がある」と述べます。ここで彼はパースの仕事に期待を寄せつつ、ボールが論理計算を用いて不等式の定理を見抜いた例にも触れています。

現代の観点から見ると、四元数・九元数の議論そのものが「式の対称性」から「代数構造」を読み取るという発想で進んでいます。これはまさに論理的な構成原理が代数の形を決める、という見方につながっています。また、一般化パウリ行列や一般化クリフォード代数へ向かう直感を既に含んでいます。

## まとめ

対称な多項式$z^2 + y^2 + x^2,\ z^3 + y^3 + x^3$を行列式として実現しようとすると、四元数や九元数に対応する非可換な積の法則が自然に現れます。

$2$次では反交換関係$vu = -uv$、$3$次ではねじれた交換関係$vu = \rho uv$となります。
