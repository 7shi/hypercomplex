[[7shi-em4]]では、時空代数$\operatorname{Cl}_{1,3}(\mathbb R)$でマクスウェル方程式を$DF=\mu_0cJ$と書き、電場と磁場が1つの2ベクトル$F$の成分であり、その分け方は時間方向$\gamma_0$の選び方によることを見ました。本記事では、時空の正規直交基底を取り替えるローレンツ変換を、回転子$x\mapsto Rx\tilde R$として扱います。特に、時間方向を変えるブーストによって、観測者による電場と磁場の分け方を調べます。2ベクトルの指数関数で回転子を作ると、2乗が$-1$の空間の2ベクトルからは三角関数による空間の回転が、2乗が$+1$の2ベクトル$\sigma_k=\gamma_k\gamma_0$からは双曲線関数によるブーストが得られます。係数の2乗の符号がここでも2つの型を分けます。回転子で$F$を変換すると電場と磁場が混ざり、$F^2$のスカラー部と擬スカラー部は不変量になります。最後に、荷電粒子の運動方程式（ローレンツ力）を$F$とベクトルの内積として書きます。

# 回転子

## 反転と回転子

[[7shi-em4]]の記号をそのまま使います。$\operatorname{Cl}_{1,3}(\mathbb R)$の生成元を$\gamma_0,\gamma_1,\gamma_2,\gamma_3$（$\gamma_0^2=1$、$\gamma_k^2=-1$）、擬スカラーを$I=\gamma_0\gamma_1\gamma_2\gamma_3$、相対ベクトルを$\sigma_k=\gamma_k\gamma_0$とします。時空の元$X$の反転$\tilde X$は、基底の積の順序を逆にする操作です。ベクトルは変えず、2ベクトルの符号を変え、$\widetilde{XY}=\tilde Y\tilde X$を満たします。

&&&def 回転子
偶部分代数の元$R$で$R\tilde R=1$を満たすものを**回転子**と呼びます。
&&&

回転子$R$と時空のベクトル$x$について、$Rx\tilde R$の2乗は

$$
(Rx\tilde R)^2=Rx\tilde RRx\tilde R=Rx^2\tilde R=x^2R\tilde R=x^2
$$

です。$x^2$はスカラーなので$R$と可換です。さらに$Rx\tilde R$はベクトルになります。

&&&prop 回転子の作用
回転子$R$と時空のベクトル$x$について、$Rx\tilde R$はベクトルです。したがって写像$x\mapsto Rx\tilde R$は線形でミンコフスキー計量を保ち、ローレンツ変換になります。
&&&

&&&prf
$Y=Rx\tilde R$は偶・奇・偶の元の積なので奇数グレードの成分だけを持ち、4次元ではグレード1と3である。$\tilde Y=\widetilde{\tilde R}\,\tilde x\,\tilde R=Rx\tilde R=Y$であり、反転はグレード1を変えずグレード3の符号を変えるので、$Y$のグレード3の成分は$0$である。
&&&

2ベクトル$B$は$\tilde B=-B$を満たすので、$\widetilde{e^B}=e^{-B}$であり、$B$は自身と可換だから$e^B\widetilde{e^B}=e^Be^{-B}=1$です。したがって2ベクトルの指数関数は回転子です。以下では、2乗が$-1$と$+1$の2ベクトルについて、作用を具体的に計算します。

## 2乗が$-1$の2ベクトル：回転

空間の面を表す2ベクトル、たとえば$I\sigma_3=\gamma_2\gamma_1$は$(I\sigma_3)^2=-1$を満たします。指数関数は複素数と同じ級数の計算で三角関数になります。

$$
e^{-I\sigma_3\theta/2}=\cos\frac\theta2-I\sigma_3\sin\frac\theta2
$$

&&&fml 空間の回転
$R=e^{-I\sigma_3\theta/2}$について

$$
R\gamma_1\tilde R=\gamma_1\cos\theta+\gamma_2\sin\theta,\qquad
R\gamma_2\tilde R=-\gamma_1\sin\theta+\gamma_2\cos\theta,
$$

$$
R\gamma_0\tilde R=\gamma_0,\qquad R\gamma_3\tilde R=\gamma_3
$$
&&&

&&&prf
$\gamma_0$と$\gamma_3$は$\gamma_2\gamma_1$と可換なので$R$とも可換であり、$R\gamma_0\tilde R=\gamma_0R\tilde R=\gamma_0$となる。$\gamma_1,\gamma_2$は$\gamma_2\gamma_1$と反可換なので$\gamma_1\tilde R=R\gamma_1$であり、$R\gamma_1\tilde R=R^2\gamma_1=e^{-I\sigma_3\theta}\gamma_1=\gamma_1\cos\theta-\gamma_2\gamma_1\gamma_1\sin\theta$である。$\gamma_2\gamma_1\gamma_1=-\gamma_2$より第1式を得る。第2式も同様である。
&&&

$(x_1,x_2)$平面の角$\theta$の回転です。$R$は$\gamma_0$と可換なので、相対ベクトル$\sigma_k=\gamma_k\gamma_0$にも$R\sigma_k\tilde R=(R\gamma_k\tilde R)\gamma_0$と同じ回転として作用します。これは$\operatorname{Cl}_{3,0}(\mathbb R)$の回転子による空間の回転（[[7shi-lie3]]）と同じものです。[[7shi-lie3]]では$Q=e^{I\sigma_3\theta/2}$による$Q^{-1}xQ$を使いました。本記事ではその逆元を$R$と書くため、指数の符号と挟む順序がともに逆になっていますが、表す回転は同じです。

## 2乗が$+1$の2ベクトル：ブースト

$\sigma_1=\gamma_1\gamma_0$は$\sigma_1^2=+1$を満たします（[[7shi-em4]]）。級数の偶数次と奇数次を分けると、三角関数の代わりに双曲線関数が現れます。

$$
e^{\sigma_1a}=\sum_{n=0}^\infty\frac{a^n\sigma_1^n}{n!}=\cosh a+\sigma_1\sinh a
$$

1と2乗が$+1$の元が張る代数は分解型複素数（[[7shi-clif2]]）であり、その指数関数は双曲線関数で書けます。

&&&fml ブースト [fml-boost]
$R=e^{\sigma_1\eta/2}$について

$$
R\gamma_0\tilde R=\gamma_0\cosh\eta+\gamma_1\sinh\eta,\qquad
R\gamma_1\tilde R=\gamma_0\sinh\eta+\gamma_1\cosh\eta,
$$

$$
R\gamma_2\tilde R=\gamma_2,\qquad R\gamma_3\tilde R=\gamma_3
$$
&&&

&&&prf
$\tilde\sigma_1=\gamma_0\gamma_1=-\sigma_1$より$\tilde R=e^{-\sigma_1\eta/2}$、$R\tilde R=1$である。$\gamma_0,\gamma_1$は$\sigma_1$と反可換なので$\gamma_0\tilde R=R\gamma_0$であり、$R\gamma_0\tilde R=R^2\gamma_0=(\cosh\eta+\gamma_1\gamma_0\sinh\eta)\gamma_0=\gamma_0\cosh\eta+\gamma_1\sinh\eta$となる。$\gamma_1$についても$\gamma_1\gamma_0\gamma_1=\gamma_0$から同様である。$\gamma_2,\gamma_3$は$\sigma_1$と可換なので変わらない。
&&&

$R\gamma_0\tilde R$を新しい時間方向$\gamma_0'$と見ると、$\gamma_0'=\cosh\eta\,(\gamma_0+\gamma_1\tanh\eta)$は、$\sigma_1$方向に速さ$v=c\tanh\eta$で動く観測者の世界線の向きです。係数$\cosh\eta=1/\sqrt{1-v^2/c^2}$はローレンツ因子で（以下、添字のない$\gamma$はローレンツ因子を表します）、パラメーター$\eta$はラピディティと呼ばれます。

$\gamma_0'$が動く観測者の世界線の向きであることは、次のように確かめられます。$\gamma_0'$の向きに進むと、$x_0=ct$が$\cosh\eta$増えるごとに$x_1$が$\sinh\eta$増えるので、$dx_1/dt=c\tanh\eta=v$です。つまりブーストは、静止した観測者を、$\sigma_1$方向に一定の速さ$v$で動く観測者に取り替える変換です。

同じ式から**時間の遅れ**が読み取れます。動く観測者が自分の時計で時間$\tau$を測る間に、その世界線は$c\tau\gamma_0'$だけ進みます。その$\gamma_0$成分は$c\tau\cosh\eta$なので、静止した観測者の時計では時間$t=\gamma\tau$が経過しています。$\gamma>1$なので、動く時計は静止した観測者から見て遅れます。ミンコフスキー計量で測った長さ$|c\tau\gamma_0'|=c\tau$が、動く観測者自身の時計の示す時間です。

同じ向きのブーストを続けると、$e^{\sigma_1a/2}e^{\sigma_1b/2}=e^{\sigma_1(a+b)/2}$よりラピディティが足し合わされます。速度で書くと$\tanh(a+b)=(\tanh a+\tanh b)/(1+\tanh a\tanh b)$であり、相対論的な速度の合成則です。速度$v_1,v_2$の合成は$(v_1+v_2)/(1+v_1v_2/c^2)$で、$v_1,v_2$が$c$より十分小さければ日常的な足し算$v_1+v_2$に近く、どちらかが$c$なら結果も$c$です。$\tanh$の値は$1$を超えないので、$c$より小さい速さをいくら合成しても$c$には届きません。光速度不変の原理は、この合成則に組み込まれています。

| 2ベクトル | 2乗 | 指数関数 | 変換 |
|---|---|---|---|
| $I\sigma_k$（空間の面） | $-1$ | $\cos$、$\sin$ | 空間の回転 |
| $\sigma_k$（時間を含む面） | $+1$ | $\cosh$、$\sinh$ | ブースト |

[[7shi-em2]]と[[7shi-em4]]では、偶部分の生成元の2乗の符号が楕円型と双曲型を分けました。同じ符号が、指数関数を三角関数と双曲線関数に分けています。空間の向きと時間の向きを保つローレンツ変換が回転とブーストの合成で得られること、および$R$と$-R$が同じ変換を与えること（二重被覆）の一般論は、本記事では扱いません。

# 場の変換

## 観測者の見る電場と磁場

電磁場$F$そのものは時空の2ベクトルで、観測者によりません。観測者が測る電場と磁場は、その観測者の時間方向を基準にした成分です。$\gamma_0'=R\gamma_0\tilde R$を時間方向とする観測者は、$\gamma_k'=R\gamma_k\tilde R$から相対ベクトル$\sigma_k'=\gamma_k'\gamma_0'=R\sigma_k\tilde R$を作り、

$$
F=\sum_kE_k'\sigma_k'+Ic\sum_kB_k'\sigma_k'
$$

と分けます。$I$は偶部分の$R$と可換なので、両辺に左から$\tilde R$、右から$R$を掛けると次のようになります。

&&&fml 観測者の見る場
$$
\tilde RFR=\sum_kE_k'\sigma_k+Ic\sum_kB_k'\sigma_k
$$
&&&

$F'=\tilde RFR$を元の$\sigma_k$で分ければ、新しい観測者の電場と磁場が読み取れます。

## ブーストによる混合

[ブースト](#fml-boost)$R=e^{\sigma_1\eta/2}$の場合を計算します。$F$の成分のうち$\sigma_1$と$I\sigma_1$は$\sigma_1$と可換で、$\sigma_2,\sigma_3,I\sigma_2,I\sigma_3$は反可換です。可換な部分$F_\parallel$は$\tilde RF_\parallel R=F_\parallel$、反可換な部分$F_\perp$は$\tilde RF_\perp R=F_\perp R^2=F_\perp e^{\sigma_1\eta}$です。

&&&fml 電場と磁場の変換
$v=c\tanh\eta$、$\gamma=\cosh\eta$、$\boldsymbol v=v\sigma_1$とすると

$$
\boldsymbol E'_\parallel=\boldsymbol E_\parallel,\qquad
\boldsymbol B'_\parallel=\boldsymbol B_\parallel,\qquad
\boldsymbol E'_\perp=\gamma(\boldsymbol E+\boldsymbol v\times\boldsymbol B)_\perp,\qquad
\boldsymbol B'_\perp=\gamma\Bigl(\boldsymbol B-\frac1{c^2}\boldsymbol v\times\boldsymbol E\Bigr)_\perp
$$

ここで$\parallel$は$\sigma_1$方向の成分、$\perp$はそれに垂直な成分です。
&&&

&&&prf
$F_\perp=E_2\sigma_2+E_3\sigma_3+Ic(B_2\sigma_2+B_3\sigma_3)$に$e^{\sigma_1\eta}=\cosh\eta+\sigma_1\sinh\eta$を右から掛ける。$\sigma_2\sigma_1=-I\sigma_3$、$\sigma_3\sigma_1=I\sigma_2$より

$$
F_\perp\sigma_1=-E_2I\sigma_3+E_3I\sigma_2+c(B_2\sigma_3-B_3\sigma_2)
$$

である（$I^2=-1$を使った）。したがって$\tilde RF_\perp R$の$\sigma_2$の係数は$E_2\cosh\eta-cB_3\sinh\eta=\gamma(E_2-vB_3)$、$\sigma_3$の係数は$\gamma(E_3+vB_2)$、$Ic\sigma_2$の係数は$\gamma(B_2+vE_3/c^2)$、$Ic\sigma_3$の係数は$\gamma(B_3-vE_2/c^2)$となる。$\boldsymbol v\times\boldsymbol B=v(0,-B_3,B_2)$、$\boldsymbol v\times\boldsymbol E=v(0,-E_3,E_2)$と比べればよい。
&&&

静止した電荷のまわりにある純粋な電場も、動く観測者からは$\boldsymbol B'_\perp=-\gamma\boldsymbol v\times\boldsymbol E/c^2$の磁場を伴って見えます。動く観測者から見ればその電荷は$-\boldsymbol v$の向きに動いており、[[7shi-em4]]で見たとおり電流を伴います。この磁場は、その電流がビオ＝サバールの法則に従って作る磁場です。

逆向きの読み方もできます。磁場$\boldsymbol B$の中を速度$\boldsymbol v$で動く電荷には、[[7shi-em1]]のローレンツ力$q\boldsymbol v\times\boldsymbol B$が働きます。電荷と一緒に動く観測者から見ると電荷は静止しているので、磁場による力は働かないはずです。実際、その観測者の電場は$\boldsymbol E'_\perp=\gamma(\boldsymbol E+\boldsymbol v\times\boldsymbol B)_\perp$であり、$\boldsymbol v\times\boldsymbol B$の力は電場による力として現れます。電気的な力と磁気的な力の区別は、観測者の選び方による相対的なものです。ベクトル解析では6成分の変換公式として覚える式ですが、ここでは$F$を回転子で挟むという1つの操作で、$\sigma_1$と反可換な成分に$e^{\sigma_1\eta}$が掛かるだけです。

## 不変量

[[7shi-em3]]で見たとおり、$F^2$はスカラー部と擬スカラー部だけを持ちます。

$$
F^2=\bigl(|\boldsymbol E|^2-c^2|\boldsymbol B|^2\bigr)+2Ic\,\boldsymbol E\cdot\boldsymbol B
$$

スカラーと擬スカラーはどちらも偶部分代数の元と可換なので

$$
F'^2=\tilde RFR\tilde RFR=\tilde RF^2R=F^2\tilde RR=F^2
$$

です。$|\boldsymbol E|^2-c^2|\boldsymbol B|^2$と$\boldsymbol E\cdot\boldsymbol B$は、回転子で結ばれた観測者、すなわち空間と時間の向きの規約を共有する観測者に共通の量です（$\boldsymbol E\cdot\boldsymbol B$は空間反転で符号を変えます）。とくに[[7shi-em3]]のヌル場（$F^2=0$）は、どの観測者から見てもヌル場です。

&&&ex 平面波のドップラー効果
$\sigma_1$方向に進む平面波$F=(1+\sigma_1)\boldsymbol E(\xi)$（$\boldsymbol E\perp\sigma_1$、$\xi=x_0-x_1$）を、同じ向きに動く観測者から見ます。$(1+\sigma_1)\sigma_1=1+\sigma_1$と、$\boldsymbol E$が$\sigma_1$と反可換であることから

$$
\tilde RFR=(1+\sigma_1)\boldsymbol E\,e^{\sigma_1\eta}=(1+\sigma_1)e^{-\sigma_1\eta}\boldsymbol E=e^{-\eta}(1+\sigma_1)\boldsymbol E
$$

となります。観測者の座標$x_0'=x\cdot\gamma_0'$、$x_1'=-x\cdot\gamma_1'$で位相を書くと、[ブースト](#fml-boost)から$x_0-x_1=e^{-\eta}(x_0'-x_1')$です。したがって観測者の見る場は

$$
F'(x_0',x_1')=e^{-\eta}(1+\sigma_1)\boldsymbol E\bigl(e^{-\eta}(x_0'-x_1')\bigr)
$$

であり、外側の因子が振幅を、引数の因子が振動数を、ともに$e^{-\eta}=\sqrt{(1-v/c)/(1+v/c)}$倍にします。これは光の縦方向のドップラー効果の因子です。
&&&

# ローレンツ力

## 固有速度

粒子と一緒に動く時計が刻む時間を**固有時**と呼びます。前節の時間の遅れで見たとおり、世界線に沿った微小な変位$dx$について$dx^2=c^2d\tau^2$であり、固有時は世界線のミンコフスキー計量による長さを$c$で割ったものです。粒子の速度が変わっても、各瞬間にこの関係が成り立ちます。

質量$m$、電荷$q$の粒子の世界線を固有時$\tau$で$x(\tau)$と表し、固有速度を$U=dx/d\tau$とします。$dx^2=c^2d\tau^2$より$U^2=c^2$です。粒子の3次元の速度を$\boldsymbol u$、ローレンツ因子を$\gamma=1/\sqrt{1-|\boldsymbol u|^2/c^2}$とすると

$$
U=\gamma\Bigl(c\gamma_0+\sum_ku_k\gamma_k\Bigr),\qquad U\gamma_0=\gamma(c+\boldsymbol u)
$$

であり、$\gamma_0$を掛けるとパラベクトルになります。

## 運動方程式

[[7shi-cla2]]の定義「ベクトルと2ベクトルの内積」にならって、2ベクトル$F$とベクトル$U$の内積を$F\cdot U=\frac12(FU-UF)$とします。これはベクトルです。この規約では$U\cdot F=-F\cdot U$なので、順序に注意します。

&&&thm ローレンツ力の成分表示
$$
m\frac{dU}{d\tau}=\frac qcF\cdot U
$$

の$\gamma_0$成分と$\gamma_k$成分は、それぞれ

$$
\frac d{d\tau}(\gamma mc)=\gamma\frac qc\boldsymbol E\cdot\boldsymbol u,\qquad
\frac d{d\tau}(\gamma m\boldsymbol u)=\gamma q(\boldsymbol E+\boldsymbol u\times\boldsymbol B)
$$

です。
&&&

&&&prf
$U$の成分を代入して$F\cdot U$を計算する。$\sigma_k\cdot\gamma_0=\frac12(\gamma_k\gamma_0\gamma_0-\gamma_0\gamma_k\gamma_0)=\gamma_k$、$\sigma_k\cdot\gamma_l=\frac12(\gamma_k\gamma_0\gamma_l-\gamma_l\gamma_k\gamma_0)$は$k=l$のとき$\gamma_0$、$k\ne l$のとき$0$である。$I\sigma_k$は$\gamma_0$と可換なので$(I\sigma_k)\cdot\gamma_0=0$であり、$(I\sigma_k)\cdot\gamma_l=-\sum_m\epsilon_{klm}\gamma_m$である（たとえば$(I\sigma_1)\cdot\gamma_2=(\gamma_3\gamma_2)\cdot\gamma_2=-\gamma_3$）。したがって磁場の寄与は

$$
-\gamma c\sum_{k,l,m}B_ku_l\epsilon_{klm}\gamma_m=\gamma c\sum_m(\boldsymbol u\times\boldsymbol B)_m\gamma_m
$$

となる。これらを集めると

$$
F\cdot U=\gamma\bigl(\boldsymbol E\cdot\boldsymbol u\bigr)\gamma_0+\gamma c\sum_k\bigl(\boldsymbol E+\boldsymbol u\times\boldsymbol B\bigr)_k\gamma_k
$$

となる。$m\,dU/d\tau$の$\gamma_0$成分は$d(\gamma mc)/d\tau$、$\gamma_k$成分は$d(\gamma mu_k)/d\tau$である。
&&&

$d\tau=dt/\gamma$で書き直すと、空間成分は$\frac d{dt}(\gamma m\boldsymbol u)=q(\boldsymbol E+\boldsymbol u\times\boldsymbol B)$、時間成分に$c$を掛けたものは$\frac d{dt}(\gamma mc^2)=q\boldsymbol E\cdot\boldsymbol u$であり、相対論的な運動量とエネルギーの式になります。電場による力と磁場による力（ローレンツ力）、および電場のする仕事が、$F$と$U$の内積という1つの式にまとまっています。

粒子の速さが$c$より十分小さければ$\gamma\approx1$で、空間成分は

$$
m\frac{d\boldsymbol u}{dt}=q(\boldsymbol E+\boldsymbol u\times\boldsymbol B)
$$

となります。これは[[7shi-em1]]のローレンツ力をニュートンの運動方程式に入れたものです。時間成分の$\gamma mc^2$は粒子のエネルギーで、$|\boldsymbol u|\ll c$では

$$
\gamma mc^2=mc^2+\frac12m|\boldsymbol u|^2+\cdots
$$

と展開できます。第2項は運動エネルギーで、時間成分の式は「運動エネルギーの変化率が電場のする仕事率$q\boldsymbol E\cdot\boldsymbol u$に等しい」ことを表します。第1項の$mc^2$は静止している粒子も持つエネルギーで、静止エネルギーと呼ばれます。磁場の力は$\boldsymbol u$に垂直なので、時間成分には現れません。

$U\cdot(F\cdot U)=0$が成り立つので、$\frac d{d\tau}U^2=2U\cdot\frac{dU}{d\tau}=0$であり、$U^2=c^2$は運動の間保たれます。また、ここでは観測者の基底を取り替えるのではなく、場と粒子の運動をともに能動的に変換します。一定の回転子$R$で$F\mapsto RF\tilde R$、$U\mapsto RU\tilde R$と変換しても$F\cdot U\mapsto R(F\cdot U)\tilde R$なので、方程式の形は変わりません。

# まとめ

ローレンツ変換を回転子$x\mapsto Rx\tilde R$として扱いました。

- **回転子**：$R\tilde R=1$を満たす偶部分の元はベクトルをベクトルに移し、ミンコフスキー計量を保ちます。2ベクトルの指数関数で作ると、2乗が$-1$の空間の2ベクトルからは三角関数による回転、2乗が$+1$の$\sigma_k$からは双曲線関数によるブーストが得られます。ラピディティは加法的で、速度の合成則を与えます。
- **場の変換**：観測者の見る場は$\tilde RFR$の成分であり、ブーストの方向と反可換な成分に$e^{\sigma_1\eta}$が掛かることから、電場と磁場の混合の公式が出ます。
- **不変量**：$F^2$のスカラー部$|\boldsymbol E|^2-c^2|\boldsymbol B|^2$と擬スカラー部$2Ic\,\boldsymbol E\cdot\boldsymbol B$は、回転子で結ばれた観測者によりません。平面波の振幅と振動数はドップラー因子$e^{-\eta}$倍になります。
- **ローレンツ力**：$m\,dU/d\tau=(q/c)F\cdot U$を座標時$t$による微分に直すと、時間成分からエネルギーの変化率$q\boldsymbol E\cdot\boldsymbol u$、空間成分から相対論的な運動量の変化率$q(\boldsymbol E+\boldsymbol u\times\boldsymbol B)$が得られます。
