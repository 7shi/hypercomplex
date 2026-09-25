[[7shi-em2]]では、マクスウェル方程式を$\mathcal DF=(\rho-\boldsymbol J/c)/\varepsilon_0$という1本の式にまとめ、共役$\bar{\mathcal D}$を掛けて波動方程式を導きました。本記事では、電荷も電流もない真空で$\mathcal DF=0$を直接解きます。進行方向の単位ベクトル$\hat{\boldsymbol k}$について$(1\pm\hat{\boldsymbol k})/2$が冪等元になることから、一定の背景場を除けば平面波は$F=(1+\hat{\boldsymbol k})\boldsymbol E$の形に限られ、$F^2=0$を満たします。正弦波では複素指数関数の$i$を擬スカラー$I$に置き換えられ、そのまま円偏光の回転が得られます。$I$と可換な$\mathcal D$の性質から、$I$を掛ける操作は電場と磁場を入れ替える双対性として働きます。最後に、$F$と反転$F^\dagger$の積がエネルギー密度とポインティングベクトルをまとめたパラベクトルになり、そのノルムが$F^2$の絶対値の2乗に等しいことからエネルギーの流れが光速を超えないこと、スカラー部の保存則がポインティングの定理になることを示します。

# 設定

[[7shi-em2]]の記号をそのまま使います。$\operatorname{Cl}_{3,0}(\mathbb R)$の生成元を$e_1,e_2,e_3$、擬スカラーを$I=e_1e_2e_3$とし、$x_0=ct$、$\mathcal D=\partial_0+D$、$F=\boldsymbol E+Ic\boldsymbol B$とします。真空のマクスウェル方程式は

$$
\mathcal DF=0
$$

です。$I$は$I^2=-1$を満たし、すべての元と可換です。

# 平面波

## 進行方向の冪等元

単位ベクトル$\hat{\boldsymbol k}$の方向に進む平面波を考えます。$F$が$\xi=x_0-\hat{\boldsymbol k}\cdot\boldsymbol x$だけの関数$F(\xi)$であるとします。

この形の意味を確かめておきます。ある時刻に$\xi$が一定の点の集まりは、$\hat{\boldsymbol k}\cdot\boldsymbol x=\text{一定}$という$\hat{\boldsymbol k}$に垂直な平面で、その上で場は一様です。この平面を**波面**と呼びます。時間が$\Delta t$だけ進むと$x_0=ct$は$c\Delta t$増えるので、同じ$\xi$の値を持つ平面は$\hat{\boldsymbol k}$の向きに$c\Delta t$だけ移動します。したがって$F(\xi)$は、波形を変えずに$\hat{\boldsymbol k}$の向きへ速さ$c$で進む波です。電荷から十分遠くでは、狭い範囲の電磁波はこの平面波で近似できます。

$F=F(\xi)$について$\partial_0F=F'$、$\partial_kF=-\hat k_kF'$より

$$
\mathcal DF=(1-\hat{\boldsymbol k})F'
$$

です。$\mathcal DF=0$は、$F'$が左から掛けた$1-\hat{\boldsymbol k}$で消えることを意味します。

$\hat{\boldsymbol k}^2=1$なので、$1$と$\hat{\boldsymbol k}$が張る部分代数は分解型複素数と同型で、[[7shi-clif2]]の冪等元が現れます。

&&&fml 進行方向の冪等元
$$
P_\pm=\frac{1\pm\hat{\boldsymbol k}}2,\qquad
P_\pm^2=P_\pm,\qquad P_+P_-=P_-P_+=0,\qquad P_++P_-=1
$$
&&&

$(1-\hat{\boldsymbol k})X=0$は$P_-X=0$、すなわち$X=P_+X$と同値です。したがって方程式からは$(P_-F)'=0$、つまり$P_-F$が$\xi$によらないことが分かります。この一定の背景場を差し引けば、残りの場は$P_-F=0$、すなわち$F=P_+F$を満たします。$\xi\to-\infty$で$F\to0$という条件や、周期的な波で平均が$0$という条件も、$P_-F=0$を保証します。以下では、背景場を除いた波の部分を扱います。

&&&thm 平面波 [thm-plane]
$F=\boldsymbol E+Ic\boldsymbol B$が$\xi=x_0-\hat{\boldsymbol k}\cdot\boldsymbol x$の関数で$F=P_+F$を満たすことは、

$$
\hat{\boldsymbol k}\cdot\boldsymbol E=0,\qquad c\boldsymbol B=\hat{\boldsymbol k}\times\boldsymbol E
$$

と同値で、このとき$F=(1+\hat{\boldsymbol k})\boldsymbol E$です。
&&&

&&&prf
$F=P_+F$は$(1-\hat{\boldsymbol k})F=0$と同値である。$\hat{\boldsymbol k}\boldsymbol E=\hat{\boldsymbol k}\cdot\boldsymbol E+I(\hat{\boldsymbol k}\times\boldsymbol E)$、$\hat{\boldsymbol k}(Ic\boldsymbol B)=Ic\,\hat{\boldsymbol k}\cdot\boldsymbol B-c\,\hat{\boldsymbol k}\times\boldsymbol B$より

$$
(1-\hat{\boldsymbol k})F=-\hat{\boldsymbol k}\cdot\boldsymbol E+\bigl(\boldsymbol E+c\,\hat{\boldsymbol k}\times\boldsymbol B\bigr)+I\bigl(c\boldsymbol B-\hat{\boldsymbol k}\times\boldsymbol E\bigr)-Ic\,\hat{\boldsymbol k}\cdot\boldsymbol B
$$

である。2ベクトル部が$0$であることは$c\boldsymbol B=\hat{\boldsymbol k}\times\boldsymbol E$であり、このとき$\hat{\boldsymbol k}\cdot\boldsymbol B=0$で擬スカラー部も$0$になる。さらにベクトル部は$\boldsymbol E+\hat{\boldsymbol k}\times(\hat{\boldsymbol k}\times\boldsymbol E)=\hat{\boldsymbol k}(\hat{\boldsymbol k}\cdot\boldsymbol E)$となり、スカラー部とともに$\hat{\boldsymbol k}\cdot\boldsymbol E=0$のときに限り$0$になる。したがって$(1-\hat{\boldsymbol k})F=0$は$\hat{\boldsymbol k}\cdot\boldsymbol E=0$と$c\boldsymbol B=\hat{\boldsymbol k}\times\boldsymbol E$の2条件と同値である。このとき$\hat{\boldsymbol k}\boldsymbol E=\hat{\boldsymbol k}\wedge\boldsymbol E=I(\hat{\boldsymbol k}\times\boldsymbol E)=Ic\boldsymbol B$なので、$F=\boldsymbol E+\hat{\boldsymbol k}\boldsymbol E=(1+\hat{\boldsymbol k})\boldsymbol E$である。
&&&

電場と磁場はともに進行方向に垂直で互いにも垂直、大きさは$|\boldsymbol E|=c|\boldsymbol B|$です。振動する量が進行方向に垂直な波を**横波**と呼び、電磁波は横波です（空気の疎密が進行方向に沿って伝わる音波は縦波です）。$c\boldsymbol B=\hat{\boldsymbol k}\times\boldsymbol E$より、$\boldsymbol E$、$\boldsymbol B$、$\hat{\boldsymbol k}$はこの順に右手系をなし、進行方向は$\boldsymbol E\times\boldsymbol B$の向きです。ベクトル解析では、これらの性質を平面波の仮定からマクスウェル方程式の4本の式に代入して1つずつ導きます。ここでは、$F$が冪等元$P_+$の像に入るという1つの条件にまとまっています。逆に、$\hat{\boldsymbol k}$に垂直な値を取る任意の$C^1$級の関数$\boldsymbol E(\xi)$について、$F=(1+\hat{\boldsymbol k})\boldsymbol E(\xi)$は$\mathcal DF=(1-\hat{\boldsymbol k})(1+\hat{\boldsymbol k})\boldsymbol E'=0$を満たす解です。波形は正弦波に限りません。

&&&rem 双曲型の正則関数
複素解析では、$z=x+iy$の任意の正則関数$f(z)$がコーシー＝リーマンの方程式を満たしました。$\mathcal DF=0$はその双曲型版で、$z$の役割を光的な変数$\xi=x_0-\hat{\boldsymbol k}\cdot\boldsymbol x$が担います。ただし$\xi$の任意の関数がそのまま解になるのではなく、冪等元$P_+$の像に入ることが必要です。空間1次元（$\mathcal D=\partial_0+e_1\partial_1$）では、$P_\pm=(1\pm e_1)/2$として

$$
F=P_+f(x_0-x_1)+P_-g(x_0+x_1)
$$

が一般解で、右向きと左向きに進む波が2つの冪等元の成分に分かれます（ダランベールの解）。分解型複素数（[[7shi-clif2]]）の冪等元による直和分解が、そのまま進行方向による波の分解になっています。楕円型の正則関数と違い、$f,g$は解析的である必要がなく、微分可能でさえあれば任意です。
&&&

## $F^2=0$

$F$の2乗を計算します。$\boldsymbol E\boldsymbol B+\boldsymbol B\boldsymbol E=2\,\boldsymbol E\cdot\boldsymbol B$より、一般の$F$について

&&&fml 場の2乗 [fml-square]
$$
F^2=\bigl(|\boldsymbol E|^2-c^2|\boldsymbol B|^2\bigr)+2Ic\,\boldsymbol E\cdot\boldsymbol B
$$
&&&

です。$F^2$はスカラー部と擬スカラー部だけを持ちます。平面波では$\boldsymbol E\hat{\boldsymbol k}=-\hat{\boldsymbol k}\boldsymbol E$より

$$
F^2=(1+\hat{\boldsymbol k})\boldsymbol E(1+\hat{\boldsymbol k})\boldsymbol E=(1+\hat{\boldsymbol k})(1-\hat{\boldsymbol k})\boldsymbol E^2=0
$$

です。スカラー部が$0$であることが$|\boldsymbol E|=c|\boldsymbol B|$に、擬スカラー部が$0$であることが$\boldsymbol E\perp\boldsymbol B$にあたります。$F^2=0$を満たす$0$でない$F$を**ヌル場**（null field）と呼びます。

&&&prop ヌル場の点ごとの形
$F=\boldsymbol E+Ic\boldsymbol B\ne0$が$F^2=0$を満たすなら、$\hat{\boldsymbol k}=\boldsymbol E\times c\boldsymbol B/|\boldsymbol E|^2$は単位ベクトルで、$F=(1+\hat{\boldsymbol k})\boldsymbol E$です。
&&&

&&&prf
[場の2乗](#fml-square)より$|\boldsymbol E|=c|\boldsymbol B|$かつ$\boldsymbol E\perp\boldsymbol B$である。$F\ne0$より$\boldsymbol E\ne0$で、$|\boldsymbol E\times c\boldsymbol B|=|\boldsymbol E|\,c|\boldsymbol B|=|\boldsymbol E|^2$だから$\hat{\boldsymbol k}$は単位ベクトルであり、$\boldsymbol E$と$\boldsymbol B$に垂直である。$\hat{\boldsymbol k}\times\boldsymbol E=\bigl((\boldsymbol E\cdot\boldsymbol E)c\boldsymbol B-(\boldsymbol E\cdot c\boldsymbol B)\boldsymbol E\bigr)/|\boldsymbol E|^2=c\boldsymbol B$となるので、[平面波](#thm-plane)で用いた代数的条件$\hat{\boldsymbol k}\cdot\boldsymbol E=0$、$c\boldsymbol B=\hat{\boldsymbol k}\times\boldsymbol E$が満たされ、$F=(1+\hat{\boldsymbol k})\boldsymbol E$である。
&&&

これは各点での代数的な特徴づけです。ヌル場が$\mathcal DF=0$を満たすことや、$\hat{\boldsymbol k}$が点によらず一定であること、場が$\xi$だけに依存することまでは意味しません。ヌル場がマクスウェル方程式を満たしても、一般には平面波とは限りません。

# 擬スカラーによる指数関数

## 円偏光

波形が三角関数で表される平面波を正弦波と呼びます。正弦波の位相を$\theta=\omega t-\boldsymbol k\cdot\boldsymbol x$と書くとき、$\omega$を角振動数、$\boldsymbol k$を波数ベクトルと呼びます。1秒あたりの振動の回数（振動数）は$f=\omega/2\pi$、波の繰り返しの長さ（波長）は$\lambda=2\pi/|\boldsymbol k|$です。前節の$\xi$の関数であるためには$\omega=c|\boldsymbol k|$、すなわち$f\lambda=c$でなければなりません。

電場の振動する向きを**偏光**と呼びます。電場が1本の直線に沿って往復するものを直線偏光、大きさを保って進行方向のまわりを回転するものを円偏光と呼びます。自然光は、さまざまな偏光が不規則に混ざったものです。

正弦波を扱うとき、ベクトル解析では複素振幅$\boldsymbol E_0e^{i(\omega t-\boldsymbol k\cdot\boldsymbol x)}$を使い、最後に実部を取ります。$\operatorname{Cl}_{3,0}(\mathbb R)$には、2乗が$-1$ですべての元と可換な$I$があるので、$i$の代わりに$I$を使った指数関数を作れます。ただしこれは実部を取る処方の書き換えではなく、電磁場全体に$e^{I\theta}$を掛けて電場と磁場の成分を混ぜる操作です。

$\hat{\boldsymbol k}$に垂直な定ベクトル$\boldsymbol E_0$と$\kappa=\omega/c>0$を取り、$\theta=\kappa\xi=\omega t-\boldsymbol k\cdot\boldsymbol x$（$\boldsymbol k=\kappa\hat{\boldsymbol k}$）とします。

&&&ex 擬スカラーによる正弦波
$$
F=(1+\hat{\boldsymbol k})\boldsymbol E_0\,e^{I\theta},\qquad e^{I\theta}=\cos\theta+I\sin\theta
$$
&&&

$e^{I\theta}$はすべての元と可換で$\partial_ae^{I\theta}=I(\partial_a\theta)e^{I\theta}$なので、$\mathcal DF=(1-\hat{\boldsymbol k})(1+\hat{\boldsymbol k})\boldsymbol E_0\,I\kappa\,e^{I\theta}=0$です。この$F$は実数の$\operatorname{Cl}_{3,0}(\mathbb R)$の元で、実部を取る操作は要りません。成分に分けると次のようになります。

&&&fml 円偏光の成分
$$
\boldsymbol E=\boldsymbol E_0\cos\theta-(\hat{\boldsymbol k}\times\boldsymbol E_0)\sin\theta,\qquad
c\boldsymbol B=\hat{\boldsymbol k}\times\boldsymbol E
$$
&&&

&&&prf
$(1+\hat{\boldsymbol k})\boldsymbol E_0=\boldsymbol E_0+I(\hat{\boldsymbol k}\times\boldsymbol E_0)$に$\cos\theta+I\sin\theta$を掛けると、$I^2=-1$より

$$
F=\bigl(\boldsymbol E_0\cos\theta-(\hat{\boldsymbol k}\times\boldsymbol E_0)\sin\theta\bigr)+I\bigl((\hat{\boldsymbol k}\times\boldsymbol E_0)\cos\theta+\boldsymbol E_0\sin\theta\bigr)
$$

である。ベクトル部が$\boldsymbol E$であり、2ベクトル部の$I$の係数が$\hat{\boldsymbol k}\times\boldsymbol E$に等しいことは、$\hat{\boldsymbol k}\times(\hat{\boldsymbol k}\times\boldsymbol E_0)=-\boldsymbol E_0$から確かめられる。
&&&

$\boldsymbol E_0$と$\hat{\boldsymbol k}\times\boldsymbol E_0$は、$\hat{\boldsymbol k}$に垂直な平面の直交する2つのベクトルで、長さは等しいです。したがって$\boldsymbol E$は、長さ$|\boldsymbol E_0|$を保ったまま、この平面内で回転します。$\hat{\boldsymbol k}$を軸とする右手系の正の回転を基準にすると、回転角は$-\theta$です。擬スカラーの指数関数は、そのまま円偏光を表します。円偏光を「右回り」「左回り」と呼び分けるときは、波を迎える側から見るか送る側から見るかで分野により流儀が逆になるので、本記事では回転の向きを式で指定するだけにします。

## 偏光面の回転としての$I$

$I$がこの平面の回転として働くことは、式の形からも読み取れます。$(1+\hat{\boldsymbol k})\hat{\boldsymbol k}=1+\hat{\boldsymbol k}$なので$(1+\hat{\boldsymbol k})I=(1+\hat{\boldsymbol k})\hat{\boldsymbol k}I=(1+\hat{\boldsymbol k})(I\hat{\boldsymbol k})$であり、$P_+$の像の上では、左から掛ける作用として$I$と2ベクトル$I\hat{\boldsymbol k}$が一致します（右から掛ける作用まで同じという意味ではありません）。$I\hat{\boldsymbol k}$は$\hat{\boldsymbol k}$に垂直な平面を表す2ベクトルで、$(I\hat{\boldsymbol k})^2=-1$です。$e^{I\theta}$は可換なので左に移せ、$(1+\hat{\boldsymbol k})$の右では$I$を$I\hat{\boldsymbol k}$に替えられます。さらに$I\hat{\boldsymbol k}$は平面内のベクトル$\boldsymbol E_0$と反可換なので

$$
(1+\hat{\boldsymbol k})\boldsymbol E_0\,e^{I\theta}=(1+\hat{\boldsymbol k})\,e^{I\theta}\boldsymbol E_0=(1+\hat{\boldsymbol k})\,e^{I\hat{\boldsymbol k}\theta}\boldsymbol E_0=(1+\hat{\boldsymbol k})\boldsymbol E_0\,e^{-I\hat{\boldsymbol k}\theta}
$$

となります。右端の$\boldsymbol E_0e^{-I\hat{\boldsymbol k}\theta}=\boldsymbol E_0\cos\theta-(\hat{\boldsymbol k}\times\boldsymbol E_0)\sin\theta$は、[[7shi-vge]]と同じく、平面の2ベクトルの指数関数による平面内の回転です。[[7shi-cla3]]で、複素数の$i$が平面の擬スカラーとして現れたのと同じ関係が、偏光面の上で成り立っています。

## 直線偏光

$e^{-I\theta}$を使えば逆向きに回る円偏光が得られます。2つの平均を取ると

$$
\frac12(1+\hat{\boldsymbol k})\boldsymbol E_0\bigl(e^{I\theta}+e^{-I\theta}\bigr)=(1+\hat{\boldsymbol k})\boldsymbol E_0\cos\theta
$$

であり、$\boldsymbol E=\boldsymbol E_0\cos\theta$の直線偏光になります。通常の複素表示で実ベクトル振幅$\boldsymbol E_0$を選んで実部を取ると、$\operatorname{Re}\bigl(\boldsymbol E_0e^{i\theta}\bigr)$は直線偏光です（複素ベクトル振幅を選べば、複素表示でも円偏光を直接表せます）。擬スカラーによる方法では、実振幅$\boldsymbol E_0$のままで指数関数そのものが円偏光になり、直線偏光は逆回りの2つの円偏光の重ね合わせです。回転の向きは指数の$I$の符号で決まります。

## 電気と磁気の双対性

円偏光で使った$e^{I\theta}$は、平面波に限らず真空の場全体に働きます。$I$は$\operatorname{Cl}_{3,0}(\mathbb R)$の中心にあり$\mathcal D$と可換なので、$\mathcal DF=0$なら任意の実数$\alpha$について

$$
\mathcal D\bigl(Fe^{I\alpha}\bigr)=(\mathcal DF)e^{I\alpha}=0
$$

です。$\alpha=\pi/2$とすると$FI=(\boldsymbol E+Ic\boldsymbol B)I=-c\boldsymbol B+I\boldsymbol E$なので

$$
\boldsymbol E\mapsto-c\boldsymbol B,\qquad c\boldsymbol B\mapsto\boldsymbol E
$$

となり、電場と磁場の役割が入れ替わります。真空のマクスウェル方程式が電場と磁場の入れ替えに対して形を保つこと（電気と磁気の双対性）は、擬スカラーを掛けるという1つの操作として表れます。一般の$\alpha$では電場と磁場が角$\alpha$で混ざり、円偏光の$e^{I\theta}$はこの双対回転を位相として使ったものと読めます。源がある場合は、源も同じ角だけ回す必要があり、磁荷（[[7shi-em1]]）を持ち込まない限りこの対称性は成り立ちません。

双対回転で$F^2$は$F^2e^{2I\alpha}$となり、スカラー部$|\boldsymbol E|^2-c^2|\boldsymbol B|^2$と擬スカラー部$2Ic\,\boldsymbol E\cdot\boldsymbol B$は複素数の実部と虚部のように混ざります。一方、次節の$FF^\dagger$は$(Fe^{I\alpha})(Fe^{I\alpha})^\dagger=Fe^{I\alpha}e^{-I\alpha}F^\dagger=FF^\dagger$で変わりません。

# エネルギーの流れ

## 場とその反転の積

$\operatorname{Cl}_{3,0}(\mathbb R)$の反転$X\mapsto X^\dagger$は、基底の積の順序を逆にする操作で、ベクトルを変えず、2ベクトルの符号を変えます。$(XY)^\dagger=Y^\dagger X^\dagger$、$I^\dagger=-I$であり、$F^\dagger=\boldsymbol E-Ic\boldsymbol B$です。

&&&fml エネルギー密度とポインティングベクトル [fml-energy]
$$
\frac{\varepsilon_0}2FF^\dagger=u+\frac{\boldsymbol S}c,\qquad
u=\frac{\varepsilon_0}2\bigl(|\boldsymbol E|^2+c^2|\boldsymbol B|^2\bigr),\qquad
\boldsymbol S=\frac1{\mu_0}\boldsymbol E\times\boldsymbol B
$$
&&&

&&&prf
$FF^\dagger=(\boldsymbol E+Ic\boldsymbol B)(\boldsymbol E-Ic\boldsymbol B)=|\boldsymbol E|^2+c^2|\boldsymbol B|^2-Ic(\boldsymbol E\boldsymbol B-\boldsymbol B\boldsymbol E)$である。$\boldsymbol E\boldsymbol B-\boldsymbol B\boldsymbol E=2\boldsymbol E\wedge\boldsymbol B=2I(\boldsymbol E\times\boldsymbol B)$より、最後の項は$2c\,\boldsymbol E\times\boldsymbol B$となる。$\varepsilon_0c=1/\mu_0c$を使えばよい。
&&&

$u$は電磁場のエネルギー密度で、電場と磁場が存在する空間には単位体積あたり$u$のエネルギーが蓄えられていると考えます。$c^2\varepsilon_0=1/\mu_0$より$u=\frac{\varepsilon_0}2|\boldsymbol E|^2+\frac1{2\mu_0}|\boldsymbol B|^2$で、第1項が電場の、第2項が磁場のエネルギーです。単位は$\mathrm{J/m^3}$です。

&&&ex 平行板コンデンサーのエネルギー
面積$S$の2枚の金属板を間隔$d$で平行に置き、電荷$+Q$と$-Q$を与えます。板の端の効果を無視すると、[[7shi-em1]]のガウスの法則から、板の間の電場は一様で大きさ$E=Q/\varepsilon_0S$、板の外では$0$です。2枚の板の電位差は$V=Ed=Qd/\varepsilon_0S$です。電荷が$q$まで充電された状態から、さらに$dq$を負の板から正の板へ運ぶには、電場に逆らって仕事$\frac{qd}{\varepsilon_0S}dq$が要ります。$0$から$Q$まで積分すると、充電に要する仕事は
$$
W=\frac{Q^2d}{2\varepsilon_0S}=\frac{\varepsilon_0}2E^2\cdot Sd
$$
です。$Sd$は電場のある領域の体積なので、この仕事はエネルギー密度$\frac{\varepsilon_0}2E^2$の電場として板の間に蓄えられていると読めます。
&&&

$\boldsymbol S$はポインティングベクトルで、エネルギーの流れの密度を表します。$\boldsymbol S$に垂直な単位面積を単位時間に通過するエネルギーが$|\boldsymbol S|$で、単位は$\mathrm{W/m^2}$です。$FF^\dagger$はスカラーとベクトルの和、つまりパラベクトルです。[場の2乗](#fml-square)$F^2$はスカラー部が正とは限らないのに対し、$FF^\dagger$のスカラー部は$\boldsymbol E$と$c\boldsymbol B$の長さの2乗の和で、つねに$0$以上です。

平面波$F=(1+\hat{\boldsymbol k})\boldsymbol E$では$(1+\hat{\boldsymbol k})^\dagger=1+\hat{\boldsymbol k}$、$(1+\hat{\boldsymbol k})^2=2(1+\hat{\boldsymbol k})$より

$$
FF^\dagger=(1+\hat{\boldsymbol k})\boldsymbol E\boldsymbol E(1+\hat{\boldsymbol k})=2|\boldsymbol E|^2(1+\hat{\boldsymbol k})
$$

なので、$u=\varepsilon_0|\boldsymbol E|^2$、$\boldsymbol S=cu\,\hat{\boldsymbol k}$です。平面波では電場と磁場のエネルギーが等しく、エネルギーは速さ$c$で進行方向に運ばれます。

一般の場でも、エネルギーが流れる速さは$c$を超えません。

&&&prop エネルギーの流れの速さ
$$
u^2-\frac{|\boldsymbol S|^2}{c^2}=\Bigl(\frac{\varepsilon_0}2\Bigr)^2\Bigl(\bigl(|\boldsymbol E|^2-c^2|\boldsymbol B|^2\bigr)^2+4c^2(\boldsymbol E\cdot\boldsymbol B)^2\Bigr)\ge0
$$

です。したがって$|\boldsymbol S|\le cu$であり、等号はヌル場（$F^2=0$）のときに限ります。
&&&

&&&prf
$\boldsymbol S/c=\varepsilon_0c\,\boldsymbol E\times\boldsymbol B$と$|\boldsymbol E\times\boldsymbol B|^2=|\boldsymbol E|^2|\boldsymbol B|^2-(\boldsymbol E\cdot\boldsymbol B)^2$より、$(2/\varepsilon_0)^2(u^2-|\boldsymbol S|^2/c^2)=(|\boldsymbol E|^2+c^2|\boldsymbol B|^2)^2-4c^2|\boldsymbol E|^2|\boldsymbol B|^2+4c^2(\boldsymbol E\cdot\boldsymbol B)^2$であり、第1項と第2項をまとめれば右辺を得る。
&&&

右辺の括弧は、[場の2乗](#fml-square)$F^2=\alpha+I\beta$の実部と虚部の2乗和$\alpha^2+\beta^2$です。パラベクトル$u+\boldsymbol S/c$のノルム$u^2-|\boldsymbol S|^2/c^2$が$F^2$の「絶対値の2乗」に等しいので、エネルギーの流れの速さ$|\boldsymbol S|/u$が$c$を超えないことは、代数の恒等式から出ます。$c$に達するのは、平面波と同じ代数的な形をしたヌル場だけです。

&&&ex 太陽光の電場
地球の位置で太陽光に垂直な面が受けるエネルギーは、およそ$1.4\times10^3\ \mathrm{W/m^2}$です。これを振幅$E_0$の直線偏光の正弦波と見なすと、$|\boldsymbol E|^2=E_0^2\cos^2\theta$の時間平均は$E_0^2/2$なので、$|\boldsymbol S|$の時間平均は$c\varepsilon_0E_0^2/2$です。これを$1.4\times10^3\ \mathrm{W/m^2}$と等しいと置くと$E_0\approx1.0\times10^3\ \mathrm{V/m}$、磁場の振幅は$E_0/c\approx3.4\times10^{-6}\ \mathrm T$となります。実際の太陽光はさまざまな振動数と偏光の混ざったものなので、これは大きさの目安です。
&&&

## ポインティングの定理

$\frac{\varepsilon_0}2FF^\dagger$に$\mathcal D$を作用させてスカラー部を取ると、$\langle\mathcal D(u+\boldsymbol S/c)\rangle_0=\partial_0u+\frac1c\nabla\cdot\boldsymbol S$です。ここで$\langle X\rangle_0$は$X$のスカラー部を表します。一方、この量はマクスウェル方程式で書き換えられます。ここでは真空に限らず、源のある場合も含めて示します。

&&&thm ポインティングの定理
マクスウェル方程式$\mathcal DF=(\rho-\boldsymbol J/c)/\varepsilon_0$の解について

$$
\partial_tu+\nabla\cdot\boldsymbol S=-\boldsymbol J\cdot\boldsymbol E
$$

が成り立ちます。
&&&

&&&prf
スカラー部には$\langle XY\rangle_0=\langle YX\rangle_0$と$\langle X^\dagger\rangle_0=\langle X\rangle_0$という性質がある。$\mathcal D=\sum_{a=0}^3e_a\partial_a$（$e_0=1$）と書くと、積の微分から

$$
\langle\mathcal D(FF^\dagger)\rangle_0=\sum_a\langle e_a(\partial_aF)F^\dagger\rangle_0+\sum_a\langle e_aF(\partial_aF^\dagger)\rangle_0
$$

である。第2項の各項は、反転を取ってから積の順序を巡回させると$\langle e_aF\,\partial_aF^\dagger\rangle_0=\langle(\partial_aF)F^\dagger e_a\rangle_0=\langle e_a(\partial_aF)F^\dagger\rangle_0$となり、第1項に等しい。したがって

$$
\langle\mathcal D(FF^\dagger)\rangle_0=2\langle(\mathcal DF)F^\dagger\rangle_0=\frac2{\varepsilon_0}\Bigl\langle\Bigl(\rho-\frac{\boldsymbol J}c\Bigr)(\boldsymbol E-Ic\boldsymbol B)\Bigr\rangle_0=-\frac2{\varepsilon_0c}\boldsymbol J\cdot\boldsymbol E
$$

である。最後の等号では、$\rho\boldsymbol E$と$\rho Ic\boldsymbol B$と$\boldsymbol J(Ic\boldsymbol B)$がスカラー部を持たないことを使った。両辺に$\varepsilon_0/2$を掛け、[エネルギー密度とポインティングベクトル](#fml-energy)を使うと$\partial_0u+\frac1c\nabla\cdot\boldsymbol S=-\frac1c\boldsymbol J\cdot\boldsymbol E$となる。$\partial_0=\frac1c\partial_t$として両辺に$c$を掛ければよい。
&&&

右辺の$\boldsymbol J\cdot\boldsymbol E$は、[[7shi-em1]]のローレンツ力で見た、電場が電流にする単位体積あたりの仕事率です。磁場による力は仕事をしないので、$\boldsymbol B$は現れません。領域$V$で積分して発散定理を使うと

$$
\frac d{dt}\int_Vu\,dV=-\oint_{\partial V}\boldsymbol S\cdot\boldsymbol n\,dS-\int_V\boldsymbol J\cdot\boldsymbol E\,dV
$$

であり、領域の中の場のエネルギーの減少は、境界から流れ出るエネルギーと電荷への仕事の和に等しくなります。電荷に渡されたエネルギーは、電荷の運動エネルギーになるか、抵抗のある導線では熱（ジュール熱）になります。ベクトル解析では、$\boldsymbol E$とアンペール＝マクスウェルの法則の内積、$\boldsymbol B$とファラデーの法則の内積を取って組み合わせ、公式$\nabla\cdot(\boldsymbol E\times\boldsymbol B)=\boldsymbol B\cdot(\nabla\times\boldsymbol E)-\boldsymbol E\cdot(\nabla\times\boldsymbol B)$で整理します。ここでは、スカラー部の対称性で2つの項が揃い、$\mathcal DF$がそのまま現れます。

運動量の保存則には、エネルギー密度とポインティングベクトルに加えて、運動量の流れを表す応力の成分が必要です。本記事ではエネルギーの保存だけを扱います。

# まとめ

真空のマクスウェル方程式$\mathcal DF=0$を、進行方向の冪等元と擬スカラーを使って調べました。

- **平面波**：$\xi=x_0-\hat{\boldsymbol k}\cdot\boldsymbol x$の関数では$\mathcal DF=(1-\hat{\boldsymbol k})F'$であり、一定の背景場を除けば$F$は冪等元$P_+=(1+\hat{\boldsymbol k})/2$の像に入ります。その形は$F=(1+\hat{\boldsymbol k})\boldsymbol E$に限られ、$\boldsymbol E\perp\hat{\boldsymbol k}$、$c\boldsymbol B=\hat{\boldsymbol k}\times\boldsymbol E$が同時に出ます。
- **$F^2=0$**：$F^2$のスカラー部$|\boldsymbol E|^2-c^2|\boldsymbol B|^2$と擬スカラー部$2Ic\,\boldsymbol E\cdot\boldsymbol B$がともに$0$になります。逆に$F^2=0$の場は、各点で平面波と同じ代数的な形をしています。
- **擬スカラーの指数関数**：$e^{I\theta}$を掛けた正弦波は実部を取らずにそのまま円偏光を表します。偏光面の上で$I$は左からの作用として2ベクトル$I\hat{\boldsymbol k}$と同じに働き、平面内の回転になります。
- **双曲型の正則関数**：真空の解は光的な変数の関数で、冪等元への分解が進行方向による波の分解になります。解析性は要りません。
- **電気と磁気の双対性**：$I$は$\mathcal D$と可換なので、$Fe^{I\alpha}$も真空の解です。$\alpha=\pi/2$で電場と磁場が入れ替わり、$FF^\dagger$は変わりません。
- **エネルギーの流れ**：$\frac{\varepsilon_0}2FF^\dagger=u+\boldsymbol S/c$はエネルギー密度とポインティングベクトルのパラベクトルです。そのノルムは$F^2$の絶対値の2乗に等しく、エネルギーの流れの速さは$c$を超えません。スカラー部についての式$\langle\mathcal D(FF^\dagger)\rangle_0=2\langle(\mathcal DF)F^\dagger\rangle_0$からポインティングの定理が出ます。
