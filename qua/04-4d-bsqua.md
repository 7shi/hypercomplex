四元数は実部も含めた4成分で4次元空間の点を表せ、その回転は左右2つの単位四元数で挟む$r_L\,q\,r_R$という形で書けます。一方、クリフォード代数$C\ell_{4,0}(\mathbb R)$では回転子で挟んで回転を表します。本記事では、この2つの「挟み込み」の対応を、二重回転の2つの回転角を和と差に組み替える操作として与える、中嶋慧氏による方法をまとめます。[[nakajima-4d-rot]]

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
p=p_{12}e_1e_2+p_{13}e_1e_3+p_{14}e_1e_4+p_{43}e_4e_3+p_{24}e_2e_4+p_{32}e_3e_2
$$

&&&rem
基底は後で比較しやすいように並べています。
&&&

同じ回転が四元数側では$q'=r_L\,q\,r_R$と表されます。6成分の回転子が$3+3$成分の$(r_L,r_R)$にどう対応するのかが問題であり、二重回転の節で見た和と差への組み替えを一般の回転面へ広げる鍵が偶部分代数です。

## 等傾回転の合成

二重回転の節で見た等傾回転を回転子で書き直し、その合成を確認します。軸は一般の単位純虚四元数$u=xi+yj+zk$（$x^2+y^2+z^2=1$）に取り、対応するベクトルを$n=xe_2+ye_3+ze_4$とします。$(1,u)$平面に対応する回転面は

$$
b=e_1n=xe_1e_2+ye_1e_3+ze_1e_4
$$

で、$e_1$と$n$が直交する単位ベクトルであることから$b^2=e_1ne_1n=-e_1^2n^2=-1$です。もう一方の回転面（直交補平面）は、擬スカラー$\omega=e_1e_2e_3e_4$を掛けると得られます。基底ごとの計算

$$
\omega\,e_1e_2=e_4e_3,\quad
\omega\,e_1e_3=e_2e_4,\quad
\omega\,e_1e_4=e_3e_2
$$

より

$$
\omega b=xe_4e_3+ye_2e_4+ze_3e_2
$$

です。$u=i$のときは$b=e_1e_2$と$\omega b=e_4e_3$となり、二重回転の節で固定した$(1,i)$平面と、$k\to j$の向きに取った第2平面が再現されます。

$\omega$は$\omega^2=1$を満たし、各生成元$e_m$と反交換するため、2ベクトルとは可換です。したがって

$$
(\omega b)^2=\omega^2b^2=-1,\quad
b\,(\omega b)=(\omega b)\,b=\omega b^2=-\omega
$$

となり、$b$と$\omega b$は可換です。積が純グレード4になる（グレード0と2の成分を持たない）ことは、2平面が互いに直交していることの反映です。可換性から指数関数は因数分解でき

$$
r=\exp\left(\frac{\theta_1}2b+\frac{\theta_2}2\,\omega b\right)
=\exp\left(\frac{\theta_1}2b\right)\exp\left(\frac{\theta_2}2\,\omega b\right)
$$

は、挟み込み$v'=r^{-1}vr$によって$b$平面を$e_1\to n$の向きに$\theta_1$、$\omega b$平面を$\theta_2$だけ回す二重回転の回転子です。各因子は他方の平面内のベクトルと可換なので、自分の平面だけを回します。

二重回転の節の結果を思い出すと、左作用$q\mapsto e^{\alpha u}q$は$(\theta_1,\theta_2)=(\alpha,-\alpha)$、右作用$q\mapsto q\,e^{\beta u}$は$(\theta_1,\theta_2)=(\beta,\beta)$の等傾回転でした。代入すれば、対応する回転子が得られます。

$$
\begin{alignedat}{2}
q&\mapsto e^{\alpha u}q,&\quad
r_-&=\exp\left(\frac\alpha2(b-\omega b)\right)=\exp\left(\frac\alpha2(1-\omega)b\right) \\
q&\mapsto q\,e^{\beta u},&\quad
r_+&=\exp\left(\frac\beta2(b+\omega b)\right)=\exp\left(\frac\beta2(1+\omega)b\right)
\end{alignedat}
$$

回転面$b$（軸$u$）は共通で、左作用か右作用かの違いは、因子$1\mp\omega$における$\omega$の符号だけに集約されます。

&&&prf 一般の軸での第2平面の向き
$u=i$の場合は二重回転の節で確認済みなので、一般の$u$について、第2平面とその向きが四元数側とクリフォード代数側で一致することを確かめる。$u$と直交する単位純虚四元数$w$を1つ取れば、$\{1,u,w,uw\}$は$\mathbb H$の正規直交基底である。$w,uw$に対応するベクトルを$\hat w,\hat m$とすれば、$(e_1,n,\hat w,\hat m)$は$\mathbb R^4$の正規直交基底で、$(u,w)=(i,j)$のとき$(e_1,e_2,e_3,e_4)$に一致する。直交対$(u,w)$の全体は連結で、正規直交基底の向き（行列式$\pm1$）は連続的に変わらないため、この基底は常に正の向きであり、積は擬スカラーに等しい。

$$
e_1n\,\hat w\hat m=\omega
$$

右から$b=e_1n$を掛けると、4つのベクトルが互いに直交することから$\hat w\hat m$と$e_1n$は可換なので

$$
\omega b=e_1n\,\hat w\hat m\,e_1n=(e_1n)^2\,\hat w\hat m=-\hat w\hat m=\hat m\hat w
$$

となる。すなわち第2の回転面$\omega b$は、$(w,uw)$平面を$uw\to w$の向きに向き付けたものである。

四元数側の作用をこの平面で確認する。$uw=-wu$と$u^2=-1$より$uwu=w$に注意すれば

$$
e^{\alpha u}\,w=\cos\alpha\,w+\sin\alpha\,uw,\quad
(uw)\,e^{\beta u}=\cos\beta\,uw+\sin\beta\,w
$$

と、左作用は$w\to uw$の向きに$\alpha$、すなわち$\omega b$の向き（$uw\to w$）では$-\alpha$回し、右作用は$\omega b$の向きに$\beta$回す。$(1,u)$平面は左右とも$1\to u$（$e_1\to n$）の向きに回るから、$\theta_2$を$\omega b$の向きで測れば、左作用は$(\theta_1,\theta_2)=(\alpha,-\alpha)$、右作用は$(\beta,\beta)$である。
&&&

合成を確認します。$b$と$\omega b$の可換性から$r_-$と$r_+$も可換で、積は指数の肩の和になります。

$$
r_-r_+=r_+r_-
=\exp\left(\frac{\alpha+\beta}2\,b+\frac{\beta-\alpha}2\,\omega b\right)
$$

これは$\theta_1=\alpha+\beta$、$\theta_2=\beta-\alpha$の二重回転の回転子であり、逆に解けば

$$
\alpha=\frac{\theta_1-\theta_2}2,\quad
\beta=\frac{\theta_1+\theta_2}2
$$

と、二重回転の節で四元数の計算から得た「和と差への組み替え」が回転子の側でも再現されます。逆向きに読めば、回転面の対$(b,\omega b)$を持つ任意の二重回転の回転子は、可換な左右の等傾回転子の積に分解されます。

$$
\exp\left(\frac{\theta_1}2b+\frac{\theta_2}2\,\omega b\right)
=\exp\left(\frac{\theta_1-\theta_2}4(1-\omega)b\right)
\exp\left(\frac{\theta_1+\theta_2}4(1+\omega)b\right)
$$

この分離を担っているのが因子$1\pm\omega$です。$\omega^2=1$より

$$
(1+\omega)(1-\omega)=1-\omega^2=0
$$

なので、$(1+\omega)$型の2ベクトルと$(1-\omega)$型の2ベクトルは、どちらの順で掛けても積が消え、互いに干渉しません。左右の情報が$\omega$の符号$\pm$で分離されるこの構造を、$\dfrac{1\pm\omega}2$を冪等元とする直和分解として整備し、一般の回転子から$(r_L,r_R)$を取り出す装置に仕上げることが、以降の節の内容です。

&&&rem 軸が異なっても可換
$\omega$が2ベクトルと可換であることから、任意の2ベクトル$B_1,B_2$に対して
$$
\{(1-\omega)B_1\}\{(1+\omega)B_2\}=(1-\omega)(1+\omega)B_1B_2=0
$$
と、どちらの順の積も$0$になります。したがって指数の肩が交換し、左作用の等傾回転子と右作用の等傾回転子は、軸が異なっていても常に可換です。これは冒頭で見た左右の作用の可換性$L_uR_v=R_vL_u$の回転子版です。
&&&

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

## グレードの振り分け

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

## 分解型双四元数とのマッピング

$C\ell_{0,3}$のグレード2の3つの元は、2乗が$-1$で循環的な積の関係$(f_2f_3)(f_3f_1)=f_1f_2$を満たすため、グレード$\{0,2\}$（偶部分）は再び四元数$\mathbb H$をなします。グレード$\{1,3\}$はその$\omega$倍です。$C\ell_{4,0}(\mathbb R)$全体（16次元）は分解型双四元数（8次元）と同型ではないため、変換は偶部分代数からの同型写像$\varphi$として明示し、暗黙の同一視は避けます。$\varphi$はグレード2の基底を四元数の基底$i,j,k$へ、擬スカラー$e_1e_2e_3e_4$を分解型双四元数の$\omega$へ移します。

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

以上の準備のもとで、回転子から$r_L,r_R$を取り出す中嶋慧氏の方法を説明します。目標は「二重回転と2つの回転角」で見た構図、すなわち左右の回転角が二重回転の2つの回転角の**和と差**になるという構図の再現です。そこで、$C\ell_{4,0}(\mathbb R)$の偶部分代数を同型写像$\varphi$で分解型双四元数へ移したうえで、回転子の情報から「和」の成分と「差」の成分を取り出す操作を作ります。その正体が、直和分解$\mathbb H\oplus\mathbb H$の各成分への射影です。

&&&ref
詳細は参考資料の(3.10)式を参照してください。[[nakajima-4d-rot]]
&&&

## 直和の行列表現と射影行列

まず冪等元を作ります。$\omega^2=1$なので、$\varepsilon=\dfrac{1+\omega}2$は冪等元です。冪等性に加えて、$\varepsilon$の上では$\omega$が$1$として振る舞うことも確かめられます。

$$
\varepsilon^2=\frac{1+2\omega+\omega^2}4=\frac{1+\omega}2=\varepsilon,\quad
\omega\varepsilon=\frac{\omega+\omega^2}2=\varepsilon
$$

もう1つの冪等元$1-\varepsilon=\dfrac{1-\omega}2$も同様に冪等で、こちらの上では逆に$\omega$が$-1$として振る舞います。

$$
(1-\varepsilon)^2=1-\varepsilon, \quad
\omega(1-\varepsilon)=-(1-\varepsilon)
$$

2つは直交し（$\varepsilon(1-\varepsilon)=0$）、足すと$1$になります。$\omega$が2つの冪等元の上でそれぞれ$\pm1$になることが、後で見る対角行列$\omega\mapsto\operatorname{diag}(1,-1)$の根拠です。

分解型双四元数の元$X=X'+\omega X''$（$X',X''\in\mathbb H$）は、この2つの冪等元の線形結合として書き直せます。

$$
X=(X'+X'')\varepsilon+(X'-X'')(1-\varepsilon)
$$

&&&prf
$2\varepsilon-1=\omega$に注意して右辺を展開する。

$$
(X'+X'')\varepsilon+(X'-X'')(1-\varepsilon)
=X'+X''(2\varepsilon-1)=X'+\omega X''
$$
&&&

つまり$X$は、和成分$X'+X''$を$\varepsilon$に、差成分$X'-X''$を$1-\varepsilon$に掛けた形に分解されます。$\varepsilon,\ 1-\varepsilon$を対角行列$\operatorname{diag}(1,0),\ \operatorname{diag}(0,1)$へ対応させる、四元数係数の対角行列でこれを表せば、$\varepsilon,\ 1-\varepsilon$は文字通りの**射影行列**になります。

$$
\varepsilon\mapsto\begin{pmatrix}1&0\\0&0\end{pmatrix},\quad
1-\varepsilon\mapsto\begin{pmatrix}0&0\\0&1\end{pmatrix},\quad
\omega=2\varepsilon-1\mapsto\begin{pmatrix}1&0\\0&-1\end{pmatrix}
$$

このとき$X=(X'+X'')\varepsilon+(X'-X'')(1-\varepsilon)$の行き先は、和$X'+X''$と差$X'-X''$を対角に並べた形になります。

$$
X'+\omega X''\mapsto
\begin{pmatrix}X'+X''&0\\0&X'-X''\end{pmatrix}
$$

この対応は積を保ちます。すなわち、和成分と差成分はそれぞれ独立に積をなし、対角行列の積の規則
$$
\begin{pmatrix}A&0\\0&B\end{pmatrix}\begin{pmatrix}C&0\\0&D\end{pmatrix}=\begin{pmatrix}AC&0\\0&BD\end{pmatrix}
$$
がそのまま成り立ちます。直和分解$\mathbb H\oplus\mathbb H$の2つの成分とは、この2つの対角成分$(X'+X'',\ X'-X'')$のことです。

$X$に$\varepsilon$を掛けると

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

もう一方の成分（差$X'-X''$）も、同じ$T$で読み出せます。$\omega$共役を$\omega^\dagger=-\omega$と表記すれば、行列表現では2つの対角成分の入れ替え$\operatorname{diag}(A,B)\mapsto\operatorname{diag}(B,A)$なので、差成分は$T(X^{\dagger})=X'-X''$として得られます。

## 回転軸への変換

回転面$p$を$\varphi$で分解型双四元数へ移してから$T$を適用すると、回転軸$q_R$が得られます。

$$
\begin{aligned}
p&=p_{12}e_1e_2+p_{13}e_1e_3+p_{14}e_1e_4+p_{43}e_4e_3+p_{24}e_2e_4+p_{32}e_3e_2 \\
\varphi(p)&=p_{12}\omega i+p_{13}\omega j+p_{14}\omega k+p_{43}i+p_{24}j+p_{32}k \\
 &=(p_{12}\omega+p_{43})i+(p_{13}\omega+p_{24})j+(p_{14}\omega+p_{32})k \\
q_R
&=T(\varphi(p)) \\
&=(p_{12}+p_{43})i+(p_{13}+p_{24})j+(p_{14}+p_{32})k
\end{aligned}
$$

$\omega$共役${}^{\dagger}$で差成分に切り替えれば、回転軸$q_L$にも変換できます。

$$
\begin{aligned}
q_L
&=-T(\varphi(p)^{\dagger}) \\
&=-T(-\omega(p_{12}i+p_{13}j+p_{14}k)+p_{43}i+p_{24}j+p_{32}k) \\
&=p_{12}i+p_{13}j+p_{14}k-p_{43}i-p_{24}j-p_{32}k \\
&=(p_{12}-p_{43})i+(p_{13}-p_{24})j+(p_{14}-p_{32})k
\end{aligned}
$$

これによって回転子を表します。

$$
r_L=\exp\left(-\frac\theta2T(\varphi(p)^{\dagger})\right),\quad r_R=\exp\left(\frac\theta2T(\varphi(p))\right)
$$

回転させるベクトル$v$はグレード1で、そのままでは偶部分代数に入らないため、$\varphi$で移せません。そこで$e_1$を掛けます。$e_1e_1$は相殺してグレード0（スカラー）になり、それ以外はグレード2（2ベクトル）になるため、偶部分代数に収まります。これによってベクトルも分解型双四元数経由で四元数に移せます。

$$
\begin{aligned}
T(\varphi(e_1v))
&=T(\varphi(e_1(v_0e_1+v_1e_2+v_2e_3+v_3e_4))) \\
&=T(\varphi(v_0+v_1e_1e_2+v_2e_1e_3+v_3e_1e_4)) \\
&=T(v_0+v_1\omega i+v_2 \omega j+v_3\omega k) \\
&=v_0+v_1i+v_2j+v_3k
\end{aligned}
$$

よって、回転は次のように分解型双四元数に移せます。

$$
\begin{aligned}
r^{-1}vr
&=\exp\left(-\frac\theta2p\right)v\exp\left(\frac\theta2p\right) \\
&\mapsto
\exp\left(-\frac\theta2T(\varphi(p)^{\dagger})\right)
T(\varphi(e_1v))
\exp\left(\frac\theta2T(\varphi(p))\right)
\end{aligned}
$$

まとめると、$C\ell_{4,0}(\mathbb R)$の回転$v'=r^{-1}vr$と同じ回転が、四元数では次のように表されます。

$$
\begin{aligned}
q'&=r_L\,q\,r_R \\
r_L&=\exp\left(\frac\theta2\,q_L\right),&q_L&=(p_{12}-p_{43})i+(p_{13}-p_{24})j+(p_{14}-p_{32})k \\
r_R&=\exp\left(\frac\theta2\,q_R\right),&q_R&=(p_{12}+p_{43})i+(p_{13}+p_{24})j+(p_{14}+p_{32})k
\end{aligned}
$$

成分を見ると、対をなす成分の和が$q_R$、差が$q_L$となり、二重回転の節で回転角に見た「和と差への組み替え」が、回転面の6成分全体に広がった形です。

自由度の勘定も成分のレベルで裏付けられます。和$(p_{12}+p_{43},\ p_{13}+p_{24},\ p_{14}+p_{32})$と差$(p_{12}-p_{43},\ p_{13}-p_{24},\ p_{14}-p_{32})$から元の6成分は復元できるため、この組み替えで情報は失われません。回転角だけでなく回転面の向きも独立に選べるため、和と差の各3成分（回転面の向き2＋回転角1）が、そのまま単位四元数$r_R,r_L$の3自由度に対応します。$3+3=6$が$\mathrm{SO}(4)$の次元と一致するのは、この可逆な組み替えの反映です。

&&&prf $q'=r_L\,q\,r_R$が同じ回転を表すこと
本文で述べたのは、回転を構成する部品（回転面$p$とベクトル$v$）がそれぞれ$T\circ\varphi$で四元数に移せるということであり、部品ごとの対応だけでは、移した部品を組み合わせた式が元と同じ回転を表すことはまだ保証されない。示すべきことは、ベクトルと四元数の同一視$v\leftrightarrow q$（本文の計算で確認したとおり$q=T(\varphi(e_1v))$）のもとで、回転後のベクトル$v'=r^{-1}vr$に対応する四元数が$r_L\,q\,r_R$に一致することである。これは「等傾回転の合成」の節の結果に帰着できる。

回転子の指数の肩を$1\pm\omega$で振り分ける。

$$
\frac\theta2p=\frac\theta4(1-\omega)p+\frac\theta4(1+\omega)p
$$

$\omega$を掛けると対をなす基底は入れ替わるから（$\omega e_4e_3=e_1e_2$など）、$(1\mp\omega)e_4e_3=\mp(1\mp\omega)e_1e_2$のように、対の片割れは$e_1$を含む側へ符号付きでまとめられる。

$$
\begin{alignedat}{2}
(1-\omega)p&=(1-\omega)b_L,&\quad
b_L&=(p_{12}-p_{43})e_1e_2+(p_{13}-p_{24})e_1e_3+(p_{14}-p_{32})e_1e_4 \\
(1+\omega)p&=(1+\omega)b_R,&\quad
b_R&=(p_{12}+p_{43})e_1e_2+(p_{13}+p_{24})e_1e_3+(p_{14}+p_{32})e_1e_4
\end{alignedat}
$$

$b_L,b_R$は$e_1n$型の回転面（軸を持つ回転面）で、その成分は本文の$q_L,q_R$の成分そのものである。「等傾回転の合成」の節で見たように$(1-\omega)$型と$(1+\omega)$型の2ベクトルはどちらの順でも積が$0$、特に可換なので、指数関数は因数分解できる。

$$
r=\exp\left(\frac\theta2p\right)
=\exp\left(\frac\theta4(1-\omega)b_L\right)\exp\left(\frac\theta4(1+\omega)b_R\right)
$$

第1因子は、$b_L=|q_L|\,\hat b$（$\hat b^2=-1$）と規格化すれば、軸$\hat b$・角$\dfrac\theta2|q_L|$の左作用の等傾回転子$\exp\left(\dfrac{\theta|q_L|}4(1-\omega)\hat b\right)$である（$q_L=0$なら恒等）。「等傾回転の合成」の節より、その作用は四元数側では左作用となる。

$$
q\mapsto\exp\left(\frac\theta2|q_L|\,\hat u\right)q=\exp\left(\frac\theta2q_L\right)q=r_L\,q
$$

ここで$\hat u$は$\hat b$に対応する軸で、$|q_L|\hat u=q_L$である。同様に第2因子は右作用$q\mapsto q\exp\left(\frac\theta2q_R\right)=q\,r_R$として作用する。回転子$r$は2つの因子の積だから、回転も2つの作用の合成となる。

$$
v'=r^{-1}vr\quad\leftrightarrow\quad q'=r_L\,q\,r_R
$$

これで本文に掲げた式が示された。$T\circ\varphi$による表示とは次のようにつながる。本文の計算のとおり$q_R=T(\varphi(p))$、$q_L=-T(\varphi(p)^{\dagger})$であり、$\varphi,T,\dagger$はいずれも準同型で指数関数と可換だから

$$
r_R=\exp\left(\frac\theta2T(\varphi(p))\right)=T(\varphi(r)),\quad
r_L=\exp\left(-\frac\theta2T(\varphi(p)^{\dagger})\right)=T(\varphi(r)^{\dagger})^{-1}
$$

が成り立つ。回転子が$\cos\frac\theta2+\sin\frac\theta2\,p$と書ける場合、すなわち$p^2=-1$のときは、準同型性より$q_R^2=T(\varphi(p^2))=-1$、同様に$q_L^2=-1$となるため、$q_L,q_R$は単位虚四元数で、指数関数は$r_L=\cos\frac\theta2+\sin\frac\theta2\,q_L$、$r_R=\cos\frac\theta2+\sin\frac\theta2\,q_R$を与える。
&&&

&&&rem 左右の回転子と直和分解
証明で得た$r_R=T(\varphi(r))$と$r_L=T(\varphi(r)^{\dagger})^{-1}=T(\varphi(r)^{\dagger})^*$が示すように、偶部分代数$\mathbb H\oplus\mathbb H$の2つの直和成分は、共役を除いて左右の回転子そのものです。
&&&

&&&rem 3次元回転
4次元回転の一部としての3次元回転は、$p_{12}=p_{13}=p_{14}=0$（回転面が$e_1$を含まない）より$q_L=-q_R$として得られます。このとき$r_L=r_R^*$となり、挟み込み$rqr^*$の形に帰着します。
&&&

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
等傾回転の回転子（$b=e_1n$は$(1,u)$平面、$\omega b$は直交補平面）：
$$
q\mapsto e^{\alpha u}q\ \leftrightarrow\ \exp\left(\frac\alpha2(1-\omega)b\right),\quad
q\mapsto q\,e^{\beta u}\ \leftrightarrow\ \exp\left(\frac\beta2(1+\omega)b\right)
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

4次元の任意の回転は直交2平面の二重回転に正規化でき、その独立な2つの回転角が和と差に組み替えられて右と左の回転角になることは、等傾回転の合成として四元数だけで確認できました。回転子の側では、左作用か右作用かの違いは因子$1\pm\omega$に集約され、一般の回転子は可換な左右の等傾回転子の積に因数分解されます。準同型写像$T\circ\varphi$は、この$1\pm\omega$の成分を四元数として読み出し、和と差への組み替えを一般の回転面の6成分全体へ広げる装置であり、その正体は2つの冪等元（行列表現では射影行列$\operatorname{diag}(1,0),\operatorname{diag}(0,1)$）による直和成分への射影、2ベクトルの言葉では自己双対・反自己双対分解です。これにより、6成分の回転子が$3+3$成分の$(r_L,r_R)$に分離されます。
