四元数は実部も含めた4成分で4次元空間の点を表せ、その回転は左右2つの単位四元数で挟む$r_L\,q\,r_R$という形で書けます。一方、クリフォード代数$C\ell_{4,0}(\mathbb R)$では回転子で挟んで回転を表します。本記事では、この2つの「挟み込み」の対応を、二重回転の2つの回転角を和と差に組み替える操作として与える、中嶋慧氏による方法をまとめます。[[nakajima]]

四元数とクリフォード代数の基本的な演算は既知とします。

# 左右作用と$M_4(\mathbb R)$

四元数$q$を基底$(1,i,j,k)$に関する実4次元の列ベクトルとみなせば、虚数単位$u=i,j,k$を左右から掛ける操作は$4\times4$実行列$L_u,R_u$で表されます。

$$
L_u\,q=uq,\quad R_u\,q=qu
$$

結合律から積の規則が決まり、左作用と右作用は互いに可換です。

$$
L_uL_v=L_{uv},\quad R_uR_v=R_{vu},\quad L_uR_v=R_vL_u
$$

左右の積$L_uR_v$（$u,v\in\{1,i,j,k\}$）は16個あり、一次独立で$4\times4$実行列全体$M_4(\mathbb R)$を張ります。これはテンソル積の同型$\mathbb H\otimes\mathbb H\cong M_4(\mathbb R)$の具体的な実現であり、2つのテンソル因子が左作用と右作用として現れています。左右から挟む回転$q'=r_L\,q\,r_R=L_{r_L}R_{r_R}\,q$は、この行列環の中に住んでいます。

&&&rem テンソル積の向き
$R_uR_v=R_{vu}$と積の順序が反転するため、正確にはテンソル積の第2因子は共役を通して右作用に対応します。
$$
u\otimes v\cong L_uR_{v^*}
$$
&&&

&&&rem $C\ell_{3,1}(\mathbb R)$としての$M_4(\mathbb R)$
$M_4(\mathbb R)$は、それ自身が1つのクリフォード代数でもあります。生成元を$e_1\cong L_iR_i,\ e_2\cong L_jR_i,\ e_3\cong L_kR_i,\ e_4\cong R_j$（2乗$+1,+1,+1,-1$）に取れば、ミンコフスキー計量のクリフォード代数$C\ell_{3,1}(\mathbb R)\cong M_4(\mathbb R)$と同一視できます。ただしこれは、本記事で後に使う$C\ell_{4,0}(\mathbb R)$とは異なる代数が、同じ行列環にたまたま同居しているだけで、作用のさせ方も異なります。

ここでの$\mathrm{SO}(4)$の構成は行列を四元数（ベクトル）に直接作用させる方式ですが、$C\ell_{3,1}(\mathbb R)$本来の構成は回転子でグレード1のベクトルを挟む方式で、後者を取ると$\operatorname{SO}(4)$ではなくローレンツ変換が得られます。
&&&

## $\mathrm{SO}(4)$の直接構成

ここで重要なのは、$M_4(\mathbb R)$の元が4次元ベクトル（四元数）に直接作用する行列だという点です。単位四元数による左右からの挟み込みは、1つの行列にまとまります。

$$
q'=r_L\,q\,r_R=L_{r_L}R_{r_R}\,q
$$

ノルムの乗法性より$|r_L\,q\,r_R|=|q|$なので、$L_{r_L}R_{r_R}$は長さを保つ直交行列、実際には行列式$1$の$\mathrm{SO}(4)$の回転行列そのものです（単位四元数は$1$と連続につながるため行列式は$+1$）。逆に$\mathrm{SO}(4)$のすべての元がこの形で書けることが知られています。つまり、左右からの作用は$\mathrm{SO}(4)$を**直接**構成します。

この構成を群の言葉で整理します。単位四元数のなす群は$\mathrm{SU}(2)$と同型で、$(r_L,r_R)\mapsto L_{r_L}R_{r_R}$は直積群$\mathrm{SU}(2)\times\mathrm{SU}(2)$から$\mathrm{SO}(4)$への全射準同型です。ただし、対$(r_L,r_R)$と$(-r_L,-r_R)$は積を取ると同じ行列になるため、この写像は2対1で、二重被覆の情報はここで失われます。

$$
\mathrm{SO}(4)\cong(\mathrm{SU}(2)\times\mathrm{SU}(2))/\{\pm(1,1)\}
=\mathrm{SU}(2)\otimes\mathrm{SU}(2)
$$

最後の$\otimes$は、この商が2つの$\mathrm{SU}(2)$の**テンソル積**として実現されることを表す記法です。テンソル積では$(-A)\otimes(-B)=A\otimes B$と符号が自動的に相殺するため、代数のテンソル積$\mathbb H\otimes\mathbb H\cong M_4(\mathbb R)$を単位四元数に制限すると、直積ではなくテンソル積$\mathrm{SU}(2)\otimes\mathrm{SU}(2)\cong\mathrm{SO}(4)$が現れるのです。一方、符号を区別する直積$\mathrm{SU}(2)\times\mathrm{SU}(2)$そのものが二重被覆$\operatorname{Spin}(4)$です。$M_4(\mathbb R)$の中に住んでいるのは$\operatorname{Spin}(4)$ではなく、商である$\mathrm{SO}(4)$です。

&&&prf
写像$(r_L,r_R)\mapsto L_{r_L}R_{r_R}$の核が$\{\pm(1,1)\}$であることを示す。核の元は任意の四元数$q$に対して$r_L\,q\,r_R=q$を満たす。$q=1$と取れば$r_Lr_R=1$、すなわち$r_R=r_L^{-1}$である。これを戻すと任意の$q$で$r_L\,q\,r_L^{-1}=q$、つまり$r_L$はすべての四元数と可換である。$r_L=a+bi+cj+dk$とおくと、$i$との可換性から$c=d=0$、$j$との可換性から$b=0$が従い、$r_L$は実数となる。単位四元数であることから$r_L=\pm1$、対応して$(r_L,r_R)=\pm(1,1)$である。逆に$\pm(1,1)$が恒等変換を与えることは明らかであるため、核はこの2元のみであり、準同型定理より商$(\mathrm{SU}(2)\times\mathrm{SU}(2))/\{\pm(1,1)\}$が像$\mathrm{SO}(4)$と同型になる。
&&&

符号の情報がこの商で失われている以上、$q$への作用は標準的な意味での**スピノル**でもありません。スピノルは、片側だけの作用$\psi\mapsto r\psi$が半角で効くために$2\pi$回転（$r\to-r$）で符号が反転するという性質で特徴づけられますが、$(r_L,r_R)\to(-r_L,-r_R)$としても$q$への作用が変わらない以上、この性質を満たしません。

対$(r_L,r_R)$を符号ごと1つの代数の元として保持するには、$4\times4$実行列では足りません。$\operatorname{Spin}(4)\cong\mathrm{SU}(2)\times\mathrm{SU}(2)$の忠実な実表現は最低8次元であることが知られており、その置き場所が後で見る偶部分代数$\mathbb H\oplus\mathbb H$です。

## 二重回転と2つの回転角

$r_L,r_R$の$3+3$自由度が$\mathrm{SO}(4)$の6自由度をどう覆うのかは、クリフォード代数に進む前に、四元数の計算だけで確認できます。鍵となる事実は、4次元の任意の回転が、直交基底を適切に取れば、互いに直交する2つの平面がそれぞれ独立な角度で回る**二重回転**に正規化できることです。

まず、片側だけの作用を調べます。虚数単位$i$を左から掛けると、基底は
$$
(1,i,j,k) \mapsto (i,-1,k,-j)
$$
と動きます。つまり左作用$q\mapsto e^{\alpha i}q$は、$(1,i)$平面と$(j,k)$平面を共に同じ角度$\alpha$だけ回す回転です。このように直交2平面が同じ角度で回る回転を**等傾回転**と呼びます。虚数単位$i$を右から掛けると、基底は
$$
(1,i,j,k) \mapsto (i,-1,-k,j)
$$
と動きます。右作用$q\mapsto q\,e^{\beta i}$も等傾回転ですが、$(1,i)$平面は$\beta$、$(j,k)$平面は$-\beta$と、第2平面の回る向きが逆になります。

左右を合成すると、この向きの違いによって2平面の回転角が分離します。$q'=e^{\alpha i}\,q\,e^{\beta i}$は$(1,i)$平面を$\alpha+\beta$、$(j,k)$平面を$\alpha-\beta$だけ回し、和と差は独立に動かせます。$(1,i)$平面の回転角を$\theta_1$、第2平面の回転角を$k\to j$の向きに測って$\theta_2$とすれば、
$$
\alpha=\dfrac{\theta_1-\theta_2}2,\quad
\beta=\dfrac{\theta_1+\theta_2}2
$$
となり、二重回転の四元数表現が得られます。

$$
q'=\exp\left(\frac{\theta_1-\theta_2}2i\right)q\exp\left(\frac{\theta_1+\theta_2}2i\right)
$$

&&&rem
この向きの取り方は、後でクリフォード代数の回転面$e_4e_3$と対応させるためです。
&&&

つまり、二重回転の2つの回転角$(\theta_1,\theta_2)$は、和と差$\dfrac{\theta_1\pm\theta_2}2$に組み替えられて、右と左の回転角になります。

特別な場合を確認します。

| $(\theta_1,\theta_2)$ | 四元数による表現                        | 回転の種類                     |
| --------------------- | ------------------------------- | ------------------------- |
| $(\theta,\theta)$     | $q\,e^{\theta i}$               | 右作用の等傾回転                  |
| $(\theta,-\theta)$    | $e^{\theta i}\,q$               | 左作用の等傾回転                  |
| $(\theta,0)$          | $e^{\theta i/2}\,q\,e^{\theta i/2}$ | $(1,i)$平面の単純回転（両側から半角） |

自由度を勘定します。二重回転は、回転面の取り方も含めれば$\mathrm{SO}(4)$の6自由度をすべて使います。一方、左右の単位四元数$r_L,r_R$はそれぞれ3自由度（回転軸の向き2＋回転角1）で、合わせて$3+3=6$と勘定が合います。ここでは回転面を$(1,i)$と$(j,k)$に固定して回転角だけを見ましたが、一般の回転面でも左右の回転子が具体的に書き下せること、つまり和と差への組み替えが回転角だけでなく6成分全体で成り立つことを、これからクリフォード代数を使って示します。

&&&rem
ここまでは$(1,i)$平面を固定するために虚数単位$i$を使いましたが、$i$である必然性はありません。$u^2=-1$を満たす単位純虚四元数$u=xi+yj+zk$（$x^2+y^2+z^2=1$）であれば、指数関数の展開は$i$のときと全く同じ形で成り立ちます。
$$
e^{\theta u}=\cos\theta+\sin\theta\,u
$$
したがって$i$を$u$に置き換えるだけで、$(1,u)$平面とそれに直交する平面の二重回転が同様に得られ、回転面の向き（3次元的な回転軸の向き）を自由に選べます。
&&&

$\theta_1$を$\theta_1+2\pi$に取り替えると（片方の平面だけ1回転追加）、回転としては同じですが、指数の肩が共に$\pi i$だけずれるため、$r_L,r_R$は共に$-1$倍されます。行列$L_{r_L}R_{r_R}$では相殺して見えなくなるこの符号を、対$(r_L,r_R)$は保持しています。1つの回転に対して符号違いの対が2通り対応するこの性質が**二重被覆**です。

## $\mathrm{SO}(4)$と$C\ell_{4,0}(\mathbb R)$

$\mathrm{SO}(4)$を回転子の共役作用として扱う枠組み、すなわち$\operatorname{Spin}(4)$の舞台は、4つの生成元がすべて$+1$に2乗する$C\ell_{4,0}(\mathbb R)$です。これもテンソル積で表せますが、第2因子は分解型四元数$\mathbb H'$となります。

$$
C\ell_{4,0}(\mathbb R)\cong M_2(\mathbb H)\cong\mathbb H\otimes\mathbb H'
$$

$C\ell_{4,0}(\mathbb R)$における回転は、回転子$r$による挟み込みです。

$$
v'=r^{-1}vr,\quad
r=\cos\frac\theta2+\sin\frac\theta2\,p
$$

ここでベクトル$v=q_0e_1+q_1e_2+q_2e_3+q_3e_4$は四元数$q=q_0+q_1i+q_2j+q_3k$と同一視します。回転面$p$は6成分の2ベクトルです。

$$
p=p_1e_1e_2+p_2e_1e_3+p_3e_1e_4+p_4e_4e_3+p_5e_2e_4+p_6e_3e_2
$$

&&&rem
基底は後で比較しやすいように並べています。
&&&

同じ回転が四元数側では$q'=r_L\,q\,r_R$と表されます。6成分の回転子が$3+3$成分の$(r_L,r_R)$にどう対応するのかが問題であり、二重回転の節で見た和と差への組み替えを一般の回転面へ広げる鍵が偶部分代数です。

# 偶部分代数と次元の降下

回転子はスカラー（グレード0）と2ベクトル（グレード2）の線形結合で、その積からはグレード4の$\omega=e_1e_2e_3e_4$も現れます。偶数グレードの元は積について閉じており、**偶部分代数**をなします。偶部分代数は、特定の生成元を含む2ベクトルを新たな生成元と見なすことで、生成元が1つ少ない（次元が1つ下がった）クリフォード代数と同型になります。新たな生成元の2乗は$(e_me_n)^2=-e_m^2e_n^2$で決まります。

回転子とスピン群は偶部分代数に住むため、この降下を3次元と4次元で並べると、置き場所の階梯が見えてきます。

- **3次元回転**は$C\ell_{3,0}(\mathbb R)$または$C\ell_{0,3}(\mathbb R)$の中で記述されます。偶部分代数の生成元を$e_1e_3,\ e_2e_3$と取れば、2乗は計量の符号によらず$(e_ae_3)^2=-e_a^2e_3^2=-1$となり、どちらの偶部分代数も$C\ell_{0,2}(\mathbb R)\cong\mathbb H$で、四元数そのものです。そのノルム$1$の元（単位四元数）がなす群が$\operatorname{Spin}(3)\cong\mathrm{SU}(2)$で、単位四元数による挟み込み$rqr^*$はこの構造の現れです。
- **4次元回転**は$C\ell_{4,0}(\mathbb R)$または$C\ell_{0,4}(\mathbb R)$の中で記述されます。偶部分代数の生成元を$f_a=e_1e_{a+1}$（$a=1,2,3$）と取れば、やはり$f_a^2=-1$となり、どちらの偶部分代数も$C\ell_{0,3}(\mathbb R)$と同型です。偶数グレード$\{0,2,4\}$の$1+6+1=8$次元が、$C\ell_{0,3}(\mathbb R)$の全グレードの$1+3+3+1=8$次元にちょうど対応します。

| 回転 | 代数 | 偶部分代数 | スピン群 |
| --- | --- | --- | --- |
| 3次元 | $C\ell_{3,0}(\mathbb R)$，$C\ell_{0,3}(\mathbb R)$ | $C\ell_{0,2}(\mathbb R)\cong\mathbb H$   | $\operatorname{Spin}(3)\cong\mathrm{SU}(2)$                      |
| 4次元 | $C\ell_{4,0}(\mathbb R)$，$C\ell_{0,4}(\mathbb R)$ | $C\ell_{0,3}(\mathbb R)\cong\mathbb H\oplus\mathbb H$  | $\operatorname{Spin}(4)\cong\mathrm{SU}(2)\times\mathrm{SU}(2)$  |

3次元回転が記述される$C\ell_{0,3}(\mathbb R)$が、1つ上の4次元では偶部分代数として再登場しています。$C\ell_{0,3}(\mathbb R)$は「四元数と、それに可換な$\omega$」という形をしており、4次元では$\omega^2=e_1^2e_2^2e_3^2e_4^2$となるため、どちらの計量でも$\omega^2=+1$です。これが**分解型双四元数**$\mathbb H\otimes\mathbb C'$であり、分解型複素数と同様に冪等元$\dfrac{1\pm\omega}2$によって$\mathbb H\oplus\mathbb H$へ直和分解できます。後で見るように、この2つの四元数成分がちょうど左右の回転子$(r_L,r_R)$を受け止めます。対は符号ごと保持されるため（$(-r_L,-r_R)\ne(r_L,r_R)$）、単位四元数の対がなすのは商$\mathrm{SU}(2)\otimes\mathrm{SU}(2)\cong\mathrm{SO}(4)$ではなく直積$\mathrm{SU}(2)\times\mathrm{SU}(2)\cong\operatorname{Spin}(4)$です。テンソル積として回転行列そのものを構成する$M_4(\mathbb R)$と、直積として対を符号ごと保持する$\mathbb H\oplus\mathbb H$という、2つの表現がここで分岐します。

# グレードの振り分け

前節の生成元$f_1=e_1e_2,\ f_2=e_1e_3,\ f_3=e_1e_4$を用いて、$C\ell_{4,0}(\mathbb R)$の偶部分代数を$C\ell_{0,3}(\mathbb R)$として書き直します。生成元の積は

$$
f_2f_3=e_1e_3e_1e_4=-e_3e_4=e_4e_3,\quad
f_3f_1=e_2e_4,\quad
f_1f_2=e_3e_2
$$

となります。8個の基底の振り分けをまとめます。

| $C\ell_{0,3}$のグレード | 基底                          | $C\ell_{4,0}$での正体          |
| ------------------------------ | --------------------------- | ------------------------------------- |
| 0                              | $1$                         | $1$（グレード0）                            |
| 1                              | $f_1,\ f_2,\ f_3$           | $e_1e_2,\ e_1e_3,\ e_1e_4$（グレード2）     |
| 2                              | $f_2f_3,\ f_3f_1,\ f_1f_2$  | $e_4e_3,\ e_2e_4,\ e_3e_2$（グレード2）     |
| 3                              | $f_1f_2f_3$                 | $-e_1e_2e_3e_4$（グレード4）                |

つまり、$C\ell_{4,0}$のグレード2にある回転面の6成分は、$e_1$を含むかどうかで$C\ell_{0,3}$のグレード1と2に3成分ずつ振り分けられ、グレード4の$\omega$はグレード3（擬スカラー）へ降ります。回転子が動くグレード$\{0,2,4\}$の8成分は、$C\ell_{0,3}$の全グレードをちょうど使い切ります。

# 分解型双四元数とのマッピング

$C\ell_{0,3}$のグレード2の3つの元は、2乗が$-1$で循環的な積の関係$(f_2f_3)(f_3f_1)=f_1f_2$を満たすため、グレード$\{0,2\}$（偶部分）は再び四元数$\mathbb H$をなします。グレード$\{1,3\}$はその$\omega$倍です。$C\ell_{4,0}(\mathbb R)$全体（16次元）は分解型双四元数（8次元）と同型ではないため、変換は偶部分代数からの同型写像$\varphi$として明示し、暗黙の同一視は避けます。$\varphi$はグレード2の基底を四元数の基底$i,j,k$へ、擬スカラー$\omega=e_1e_2e_3e_4$を分解型双四元数の$\omega$へ移します。

&&&def 偶部分代数から分解型双四元数への同型写像
$$
\varphi:C\ell_{4,0}^0(\mathbb R)\to\mathbb H\otimes\mathbb C',\quad
e_4e_3\mapsto i,\quad
e_2e_4\mapsto j,\quad
e_3e_2\mapsto k,\quad
e_1e_2e_3e_4\mapsto\omega
$$
&&&

行き先の$\{i,j,k,\omega\}$が生成する代数が前節の分解型双四元数で、任意の元は$X'+\omega X''$（$X',X''\in\mathbb H$）と書けます。$\omega$を掛けるとグレード1と2が入れ替わるため、$e_1$を含む2ベクトルは$\omega i,\omega j,\omega k$へ移ります。

$$
\varphi(e_1e_2)=\omega i,\quad
\varphi(e_1e_3)=\omega j,\quad
\varphi(e_1e_4)=\omega k
$$

# 射影としての準同型写像

以上の準備のもとで、回転子から$r_L,r_R$を取り出す中嶋慧氏の方法を説明します。目標は「二重回転と2つの回転角」で見た構図、すなわち左右の回転角が二重回転の2つの回転角の**和と差**になるという構図の再現です。そこで、$C\ell_{4,0}(\mathbb R)$の偶部分代数を同型写像$\varphi$で分解型双四元数へ移したうえで、回転子の情報から「和」の成分と「差」の成分を取り出す操作を作ります。その正体が、直和分解$\mathbb H\oplus\mathbb H$の各成分への射影です。詳細はこちらの資料の(3.10)式を参照してください。[[nakajima]]

## 直和の行列表現と射影行列

分解型双四元数の元$X=X'+\omega X''$（$X',X''\in\mathbb H$）を、和$X'+X''$と差$X'-X''$を対角に並べた、四元数係数の対角行列で表します。

$$
X'+\omega X''\mapsto
\begin{pmatrix}X'+X''&0\\0&X'-X''\end{pmatrix},\quad
\omega\mapsto
\begin{pmatrix}1&0\\0&-1\end{pmatrix}
$$

$\omega$の行き先は$\omega^2=1$を満たす対角行列で、$\omega$を固有値$\pm1$で対角化した形です。この対応は積を保ちます。すなわち、和成分と差成分はそれぞれ独立に積をなし、対角行列の積の規則$\operatorname{diag}(A,B)\operatorname{diag}(C,D)=\operatorname{diag}(AC,BD)$がそのまま成り立ちます。直和分解$\mathbb H\oplus\mathbb H$の2つの成分とは、この2つの対角成分$(X'+X'',\ X'-X'')$のことです。

&&&prf
$\omega$が中心的（$\mathbb H$の元と可換）で$\omega^2=1$を満たすことを用いて、$X=X'+\omega X''$と$Y=Y'+\omega Y''$の積を展開する。

$$
XY=(X'+\omega X'')(Y'+\omega Y'')
=(X'Y'+X''Y'')+\omega(X'Y''+X''Y')
$$

積$XY$の和成分と差成分を計算すると

$$
\begin{aligned}
(X'Y'+X''Y'')+(X'Y''+X''Y')&=(X'+X'')(Y'+Y'') \\
(X'Y'+X''Y'')-(X'Y''+X''Y')&=(X'-X'')(Y'-Y'')
\end{aligned}
$$

と、それぞれ$X,Y$の和成分同士・差成分同士の積に因数分解される。よって上の対応は積を保ち、積は対角成分ごとに閉じている。
&&&

この表現では、冪等元$\varepsilon=\dfrac{1+\omega}2$は文字通りの**射影行列**になります。

$$
\varepsilon=\frac{1+\omega}2\mapsto
\begin{pmatrix}1&0\\0&0\end{pmatrix}
$$

冪等性と「$\varepsilon$の上では$\omega$が$1$になる」ことは、行列の形からも、$\omega^2=1$を使った代数の計算からも確かめられます。

$$
\varepsilon^2=\frac{1+2\omega+\omega^2}4=\frac{1+\omega}2=\varepsilon,\quad
\omega\varepsilon=\frac{\omega+\omega^2}2=\varepsilon
$$

もう1つの冪等元$\dfrac{1-\omega}2\mapsto\operatorname{diag}(0,1)$は差成分への射影で、2つの射影は互いに直交し、足すと$1$になります。$X$に$\varepsilon$を掛けると

$$
\varepsilon X=\varepsilon X'+\varepsilon\omega X''=\varepsilon(X'+X'')
$$

と、和成分$X'+X''$を係数とする$\varepsilon\,\mathbb H$の元だけが生き残ります。

&&&rem $M_2(\mathbb H)$の対角部分
四元数係数の対角行列が現れるのは偶然ではありません。$C\ell_{4,0}(\mathbb R)\cong M_2(\mathbb H)$の実現では、偶部分代数がちょうど対角行列全体$\operatorname{diag}(A,B)\cong\mathbb H\oplus\mathbb H$として現れます。上の表現は、その対角部分を分解型双四元数の側から書いたものです。
&&&

## 準同型写像$T$

射影$\varepsilon X=\varepsilon(X'+X'')$で生き残るのは係数の四元数$X'+X''$だけなので、直和成分$\varepsilon\,\mathbb H$を対応$\varepsilon A\leftrightarrow A$（$A\in\mathbb H$）によって$\mathbb H$と同一視し、係数を直接読み出す写像を$T$とします。

&&&def 分解型双四元数から四元数への準同型写像
$T(X'+\omega X'')=X'+X''\quad(X',X'' \in \mathbb H)$
&&&

$T$は$\omega=1$と読み替えることに相当します。行列表現では左上の対角成分を読み出す対応$\operatorname{diag}(A,B)\mapsto A$、直和の言葉では第1成分への射影$\mathbb H\oplus\mathbb H\to\mathbb H$です。前節で見たとおり積は対角成分ごとに閉じているため、成分の読み出しは自動的に積を保ちます。

$$
T(XY)=T(X)T(Y)
$$

つまり$T$の準同型性は、「射影行列で片方の対角成分だけを見る」ことの帰結です。

もう一方の成分（差$X'-X''$）も、同じ$T$で読み出せます。$\omega$共役${}^{\dagger}$（$\omega\mapsto-\omega$、$X'+\omega X''\mapsto X'-\omega X''$）は、行列表現では2つの対角成分の入れ替え$\operatorname{diag}(A,B)\mapsto\operatorname{diag}(B,A)$なので、差成分は$T(X^{\dagger})=X'-X''$として得られます。

## 回転軸への変換

回転面$p$を$\varphi$で分解型双四元数へ移してから$T$を適用すると、回転軸$q_R$が得られます。

$$
\begin{aligned}
p&=p_1e_1e_2+p_2e_1e_3+p_3e_1e_4+p_4e_4e_3+p_5e_2e_4+p_6e_3e_2 \\
\varphi(p)&=p_1\omega i+p_2\omega j+p_3\omega k+p_4i+p_5j+p_6k \\
 &=(p_1\omega+p_4)i+(p_2\omega+p_5)j+(p_3\omega+p_6)k \\
q_R
&=T(\varphi(p)) \\
&=(p_1+p_4)i+(p_2+p_5)j+(p_3+p_6)k
\end{aligned}
$$

$\omega$共役${}^{\dagger}$で差成分に切り替えれば、回転軸$q_L$にも変換できます。

$$
\begin{aligned}
q_L
&=-T(\varphi(p)^{\dagger}) \\
&=-T(-\omega(p_1i+p_2j+p_3k)+p_4i+p_5j+p_6k) \\
&=p_1i+p_2j+p_3k-p_4i-p_5j-p_6k \\
&=(p_1-p_4)i+(p_2-p_5)j+(p_3-p_6)k
\end{aligned}
$$

これによって回転子を表します。

$$
r_L=\exp\left(-\frac\theta2T(\varphi(p)^{\dagger})\right),\quad r_R=\exp\left(\frac\theta2T(\varphi(p))\right)
$$

まとめると、$C\ell_{4,0}(\mathbb R)$の回転$v'=r^{-1}vr$と同じ回転が、四元数では次のように表されます。

$$
\begin{aligned}
q'&=r_L\,q\,r_R \\
r_L&=\cos\frac\theta2+\sin\frac\theta2\,q_L,&q_L&=(p_1-p_4)i+(p_2-p_5)j+(p_3-p_6)k \\
r_R&=\cos\frac\theta2+\sin\frac\theta2\,q_R,&q_R&=(p_1+p_4)i+(p_2+p_5)j+(p_3+p_6)k
\end{aligned}
$$

成分を見ると、対をなす成分の和が$q_R$、差が$q_L$——二重回転の節で回転角に見た「和と差への組み替え」が、回転面の6成分全体に広がった形です。

&&&prf
ベクトルと四元数の対応から確認する。$v$自体はグレード1なので$\varphi$の定義域に入らないが、奇数グレード同士の積は偶数グレードになるため、$e_1$を掛けた$e_1v$は偶部分代数の元となり、$\varphi$で移せる。$e_1^2=1$に注意すると

$$
e_1v=q_0+q_1e_1e_2+q_2e_1e_3+q_3e_1e_4,\quad
\varphi(e_1v)=q_0+\omega(q_1i+q_2j+q_3k)
$$

となり、$T(\varphi(e_1v))=q_0+q_1i+q_2j+q_3k=q$が得られる。$\varphi(e_1v)$は「$X'$が実部、$X''$が虚部」という形をしており、この形の元は$T$の像から一意に復元できるため、ベクトル$v$と四元数$q=T(\varphi(e_1v))$は1対1に対応する。

次に、偶部分代数における$e_1$による共役が、$\varphi$を通して$\omega$共役に対応することを示す。

$$
\varphi(e_1xe_1)=\varphi(x)^{\dagger}
$$

基底ごとに確認する。スカラー$1$は不変で、$\dagger$でも不変。$e_1$を含まない2ベクトル$e_me_n$（$m,n\ne1$）は、$e_1$と2回の反交換を経て可換となるため不変で、行き先$i,j,k$も$\dagger$で不変。$e_1$を含む2ベクトルは$e_1(e_1e_a)e_1=e_ae_1=-e_1e_a$と符号が反転し、行き先$\omega i,\omega j,\omega k$も$\dagger$で反転。$\omega=e_1e_2e_3e_4$は$e_1$を1つ含むため$e_1$と反交換し、$e_1\omega e_1=-\omega$と反転する。すなわち両辺はどの基底でも一致し、いずれも$X'+\omega X''\mapsto X'-\omega X''$として作用する。

以上を組み合わせる。回転子$r$と逆元$r^{-1}$は偶部分代数の元だから、回転$v'=r^{-1}vr$の両辺に$e_1$を掛けて$e_1^2=1$を挟むと

$$
e_1v'=e_1r^{-1}vr=(e_1r^{-1}e_1)(e_1v)r
$$

となる。右辺の3つの因子はすべて偶元なので$\varphi$で移せて、先の補題より

$$
\varphi(e_1v')=\varphi(r^{-1})^{\dagger}\,\varphi(e_1v)\,\varphi(r)
=(\varphi(r)^{\dagger})^{-1}\,\varphi(e_1v)\,\varphi(r)
$$

となる。ここで$\dagger$は$\omega\mapsto-\omega$による代数の自己同型なので、逆元を取る操作と可換であることを用いた。準同型$T$を適用すれば

$$
q'=T(\varphi(e_1v'))=T(\varphi(r)^{\dagger})^{-1}\,q\,T(\varphi(r))
$$

となり、クリフォード代数の回転は四元数の挟み込み$q'=r_L\,q\,r_R$として実現される。

$$
r_L=T(\varphi(r)^{\dagger})^{-1},\quad r_R=T(\varphi(r))
$$

最後に$r=\exp\left(\frac\theta2p\right)$を代入する。$\varphi,T,\dagger$はいずれも準同型で指数関数と可換だから

$$
r_R=\exp\left(\frac\theta2T(\varphi(p))\right),\quad
r_L=\exp\left(\frac\theta2T(\varphi(p)^{\dagger})\right)^{-1}
=\exp\left(-\frac\theta2T(\varphi(p)^{\dagger})\right)
$$

となり、$q_R=T(\varphi(p))$と$q_L=-T(\varphi(p)^{\dagger})$が本文の式を与える。回転子が$\cos\frac\theta2+\sin\frac\theta2\,p$と書ける場合、すなわち$p^2=-1$のときは、準同型性より$q_R^2=T(\varphi(p^2))=-1$、同様に$q_L^2=-1$となるため、$q_L,q_R$は単位虚四元数で、指数関数は$r_L=\cos\frac\theta2+\sin\frac\theta2\,q_L$、$r_R=\cos\frac\theta2+\sin\frac\theta2\,q_R$を与える。
&&&

&&&rem 左右の回転子と直和分解
証明で得た$r_R=T(\varphi(r))$と$r_L=T(\varphi(r)^{\dagger})^{-1}=T(\varphi(r)^{\dagger})^*$が示すように、偶部分代数$\mathbb H\oplus\mathbb H$の2つの直和成分は、共役を除いて左右の回転子そのものです。
&&&

&&&rem 3次元回転
4次元回転の一部としての3次元回転は、$p_1=p_2=p_3=0$（回転面が$e_1$を含まない）より$q_L=-q_R$として得られます。このとき$r_L=r_R^*$となり、挟み込み$rqr^*$の形に帰着します。
&&&

## 二重回転との整合

「二重回転と2つの回転角」で四元数だけから導いた式が、$\varphi$と$T$から再現されることを確認します。回転面を$e_1e_2$と$e_4e_3$に取れば、回転子は

$$
r=\exp\left(\frac{\theta_1}2e_1e_2+\frac{\theta_2}2e_4e_3\right)
$$

となります。$e_1e_2$は$(1,i)$平面、$e_4e_3$は$k\to j$の向きに取った第2平面の回転面に対応しており、これは前半の二重回転の設定そのものです。$\varphi$で移すと、$e_1e_2\mapsto\omega i$、$e_4e_3\mapsto i$より、指数の肩は1つにまとまります。

$$
\varphi\left(\frac{\theta_1}2e_1e_2+\frac{\theta_2}2e_4e_3\right)
=\frac{\theta_1\omega+\theta_2}2\,i
$$

証明で見たように、$r_R=T(\varphi(r))$は$\omega=1$への読み替え、$r_L=T(\varphi(r)^{\dagger})^{-1}$は$\omega=-1$への読み替えと逆元（指数の符号反転）でした。指数の肩はそれぞれ$\dfrac{\theta_1+\theta_2}2i$と$\dfrac{\theta_1-\theta_2}2i$となり、

$$
q'=\exp\left(\frac{\theta_1-\theta_2}2i\right)q\exp\left(\frac{\theta_1+\theta_2}2i\right)
$$

と、二重回転の節の式が確かに再現されます。この証明では$\omega=\pm1$成分への射影という代数の分解によって、等傾回転の合成と同じ和と差への組み替えが再現されます。

自由度の勘定も成分のレベルで裏付けられます。和$(p_1+p_4,\ p_2+p_5,\ p_3+p_6)$と差$(p_1-p_4,\ p_2-p_5,\ p_3-p_6)$から元の6成分は復元できるため、この組み替えで情報は失われません。回転角だけでなく回転面の向きも独立に選べるため、和と差の各3成分（回転面の向き2＋回転角1）が、そのまま単位四元数$r_R,r_L$の3自由度に対応します。$3+3=6$が$\mathrm{SO}(4)$の次元と一致するのは、この可逆な組み替えの反映です。

&&&rem 自己双対分解
和と差への組み替えは、2ベクトルの言葉では対をなす基底の和と差$e_1e_2\pm e_4e_3$などへの組み替えです。$\omega$を掛けると対をなす2ベクトルが入れ替わるため（$\omega\,e_1e_2=e_4e_3$）、和は不変で、差は符号が反転します。2ベクトルに$\omega$を掛ける操作はホッジ双対に相当し、不変な2ベクトルは**自己双対**、反転するものは**反自己双対**と呼ばれます。この分解はリー代数の分解$\mathfrak{so}(4)\cong\mathfrak{su}(2)\oplus\mathfrak{su}(2)$（$\operatorname{Spin}(4)\cong\mathrm{SU}(2)\times\mathrm{SU}(2)$のリー代数版）そのものです。
&&&

# まとめ

&&&
左右作用による$\mathrm{SO}(4)$の直接構成：
$$
\mathbb H\otimes\mathbb H\cong M_4(\mathbb R),\quad
\mathrm{SU}(2)\otimes\mathrm{SU}(2)\cong\mathrm{SO}(4)
$$
二重回転の分解（$e_1e_2$面で$\theta_1$、$e_4e_3$面で$\theta_2$の回転）：
$$
q'=\exp\left(\frac{\theta_1-\theta_2}2i\right)q\exp\left(\frac{\theta_1+\theta_2}2i\right)
$$
偶部分代数の階梯とスピン群：
$$
\begin{alignedat}{2}
C\ell_{3,0}^0(\mathbb R)&\cong C\ell_{0,3}^0(\mathbb R)\cong C\ell_{0,2}(\mathbb R)\cong\mathbb H,&\quad
\operatorname{Spin}(3)&\cong\mathrm{SU}(2) \\
C\ell_{4,0}^0(\mathbb R)&\cong C\ell_{0,4}^0(\mathbb R)\cong C\ell_{0,3}(\mathbb R)\cong\mathbb H\oplus\mathbb H,&\quad
\operatorname{Spin}(4)&\cong\mathrm{SU}(2)\times\mathrm{SU}(2)
\end{alignedat}
$$
偶部分代数から分解型双四元数への同型写像：
$$
\varphi:C\ell_{4,0}^0(\mathbb R)\to\mathbb H\otimes\mathbb C',\quad
e_4e_3\mapsto i,\quad
e_2e_4\mapsto j,\quad
e_3e_2\mapsto k,\quad
e_1e_2e_3e_4\mapsto\omega
$$
分解型双四元数の行列表現・射影行列と準同型写像$T$：
$$
\omega \mapsto \begin{pmatrix}1&0\\0&-1\end{pmatrix} \quad (\omega^2 \mapsto I), \quad
\varepsilon=\frac{1+\omega}2\mapsto\begin{pmatrix}1&0\\0&0\end{pmatrix} \quad (\varepsilon^2=\varepsilon)
$$
$$
X'+\omega X''\mapsto\begin{pmatrix}X'+X''&0\\0&X'-X''\end{pmatrix},\quad
T(X'+\omega X'')=X'+X''
$$
回転子の対応：
$$
r_L=\exp\left(-\frac\theta2T(\varphi(p)^{\dagger})\right),\quad
r_R=\exp\left(\frac\theta2T(\varphi(p))\right)
$$
&&&

四元数の左右から挟む作用は$\mathbb H\otimes\mathbb H\cong M_4(\mathbb R)$を張り、$\mathrm{SO}(4)$の回転行列そのものを直接構成します。

テンソル積では対の符号が相殺するため、そこに現れるのは直積$\operatorname{Spin}(4)\cong\mathrm{SU}(2)\times\mathrm{SU}(2)$ではなく商$\mathrm{SO}(4)\cong\mathrm{SU}(2)\otimes\mathrm{SU}(2)$であり、符号の情報がこの商で失われている以上、標準的な意味でのスピノルではなく、$M_4(\mathbb R)$がそのまま作用する表現空間という位置づけです。符号を区別する$\operatorname{Spin}(4)$の置き場所は、偶部分代数の階梯の上にあります。3次元回転では偶部分代数$C\ell_{0,2}(\mathbb R)\cong\mathbb H$に$\operatorname{Spin}(3)\cong\mathrm{SU}(2)$が住むように、4次元回転では偶部分代数$C\ell_{0,3}(\mathbb R)\cong\mathbb H\oplus\mathbb H$（分解型双四元数）に$\operatorname{Spin}(4)$が対$(r_L,r_R)$として住みます。

4次元の任意の回転は直交2平面の二重回転に正規化でき、その独立な2つの回転角が和と差に組み替えられて右と左の回転角になることは、等傾回転の合成として四元数だけで確認できました。準同型写像$T\circ\varphi$は、この和と差への組み替えを一般の回転面の6成分全体へ広げる装置であり、その正体は2つの冪等元（行列表現では射影行列$\operatorname{diag}(1,0),\operatorname{diag}(0,1)$）による直和成分への射影、2ベクトルの言葉では自己双対・反自己双対分解です。これにより、6成分の回転子が$3+3$成分の$(r_L,r_R)$に分離されます。
