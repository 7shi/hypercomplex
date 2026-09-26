[[7shi-em1]]では、静的な場を$F=-DP$、$P=\varphi-c\boldsymbol A$と1つのポテンシャルの微分で書き、$F$のスカラー部が$c\,\nabla\cdot\boldsymbol A$になることを見ました。本記事では、時空代数$\operatorname{Cl}_{1,3}(\mathbb R)$（[[7shi-em4]]）でこれを時間を含む場に広げます。ポテンシャルは時空のベクトル$A$になり、$DA$の2ベクトル部$D\wedge A$が電磁場$F$、スカラー部$D\cdot A$はゲージ変換によって変わる量です。$D\cdot A=0$と選ぶ（ローレンスゲージ、Lorenz gauge）と$F=DA$となり、マクスウェル方程式は$D^2A=\mu_0cJ$、すなわち成分ごとの波動方程式になります。後半では、[[7shi-em3]]の$\frac{\varepsilon_0}2FF^\dagger$を時空版のエネルギー運動量$T(a)=-\frac{\varepsilon_0}2FaF$として回収し、その発散がローレンツ力の密度になることを示します。

# 設定

[[7shi-em4]]の記号をそのまま使います。$\operatorname{Cl}_{1,3}(\mathbb R)$の生成元を$\gamma_\mu$（$\gamma_0^2=1$、$\gamma_k^2=-1$）、逆元を$\gamma^\mu$、ディラック作用素を$D=\sum_\mu\gamma^\mu\partial_\mu$とし、相対ベクトルを$\sigma_k=\gamma_k\gamma_0$、擬スカラーを$i$とします。電磁場$F=\sum_kE_k\sigma_k+ic\sum_kB_k\sigma_k$は時空の2ベクトル、電流$J=c\rho\gamma_0+\sum_kJ_k\gamma_k$は時空のベクトルで、マクスウェル方程式は

$$
DF=\mu_0cJ,\qquad\text{すなわち}\qquad D\cdot F=\mu_0cJ,\quad D\wedge F=0
$$

です。$D^2=\partial_0^2-\Delta$を$\square$と書きます。

# ベクトルポテンシャル

## 時空のポテンシャル

ポテンシャルを使うと、源を持たない2本の式を自動的に満たせます。残る式も、ゲージを適切に選べば成分ごとの波動方程式になります。

&&&def 時空のポテンシャル
スカラーポテンシャル$\varphi$とベクトルポテンシャル$\boldsymbol A=\sum_kA_k\sigma_k$から、時空のベクトル

$$
A=\varphi\,\gamma_0+c\sum_kA_k\gamma_k
$$

を作ります。
&&&

$\gamma_0$を掛けると$\gamma_0A=\varphi+c\sum_kA_k\gamma_0\gamma_k=\varphi-c\boldsymbol A$であり、[[7shi-em1]]の$P$が現れます。[[7shi-em4]]の電流$J$と$\gamma_0J=c(\rho-\boldsymbol J/c)$の関係と同じ形です。

$A$はベクトルなので、$DA$はスカラー部$D\cdot A$と2ベクトル部$D\wedge A$だけを持ちます。

&&&fml ポテンシャルの微分
$$
D\wedge A=\boldsymbol E+ic\boldsymbol B,\qquad
\boldsymbol E=-\nabla\varphi-\partial_t\boldsymbol A,\qquad
\boldsymbol B=\nabla\times\boldsymbol A,
$$

$$
D\cdot A=\frac1c\partial_t\varphi+c\,\nabla\cdot\boldsymbol A
$$
&&&

&&&prf
[[7shi-em4]]の作用素の分割$D(\gamma_0H)=\bar{\mathcal D}H$を$H=\gamma_0A=\varphi-c\boldsymbol A$に使うと、$DA=\bar{\mathcal D}(\varphi-c\boldsymbol A)=(\partial_0-D_3)(\varphi-c\boldsymbol A)$である。ここで$D_3=\sum_k\sigma_k\partial_k$は空間のディラック作用素である。展開すると

$$
DA=\partial_0\varphi-c\,\partial_0\boldsymbol A-\nabla\varphi+c\,D_3\boldsymbol A
=\bigl(\partial_0\varphi+c\,\nabla\cdot\boldsymbol A\bigr)+\bigl(-\nabla\varphi-c\,\partial_0\boldsymbol A\bigr)+ic\,\nabla\times\boldsymbol A
$$

となる。$c\,\partial_0=\partial_t$より、偶部分を$\operatorname{Cl}_{3,0}(\mathbb R)$と見たときのベクトル部が$\boldsymbol E$、$i$の項が$ic\boldsymbol B$、スカラー部が$D\cdot A$である。
&&&

$\boldsymbol E=-\nabla\varphi-\partial_t\boldsymbol A$と$\boldsymbol B=\nabla\times\boldsymbol A$は、ベクトル解析でポテンシャルから場を作る式そのものです。静的な場合と違い、電場には$-\partial_t\boldsymbol A$の項が加わります。磁場が時間変化すると電場は渦を持ち（[[7shi-em2]]のファラデーの法則）、勾配だけでは表せないためです。時間変化する磁束を囲む経路では電場の周回積分が$0$とは限らないので、固定した時刻の電場を2点間で積分した値は経路によって異なり得ます。固定した経路$C:\boldsymbol a\to\boldsymbol b$について

$$
\int_C\boldsymbol E\cdot d\boldsymbol l=\varphi(\boldsymbol a,t)-\varphi(\boldsymbol b,t)-\partial_t\int_C\boldsymbol A\cdot d\boldsymbol l
$$

であり、左辺は物理的な電場で決まりますが、右辺の$\varphi$と$\boldsymbol A$への分担はゲージによって変わります。一般のゲージ変換では$\varphi$の差そのものも変わるので、[[7shi-em1]]のように電位の差を測定される電圧と同一視することはできません。

## 同次方程式

偏微分は交換するので、$D\wedge D=\sum_{\mu,\nu}\gamma^\mu\wedge\gamma^\nu\partial_\mu\partial_\nu=0$です（$\gamma^\mu\wedge\gamma^\nu$は$\mu,\nu$について反対称）。偏微分は順序を入れ替えても同じですが、外積は順序を入れ替えると符号が変わるので、2回の外微分は打ち消し合います。ベクトル解析の恒等式$\nabla\times\nabla\varphi=0$と$\nabla\cdot(\nabla\times\boldsymbol A)=0$は、この1つの式の成分にあたります（[[7shi-em2]]）。したがって$F=D\wedge A$と置けば

$$
D\wedge F=D\wedge(D\wedge A)=0
$$

が自動的に成り立ちます。マクスウェル方程式のうち源を持たない2本（ファラデーの法則と$\nabla\cdot\boldsymbol B=0$）は、ポテンシャルを使うと恒等式になります。逆に$D\wedge F=0$なら、局所的には$F=D\wedge A$となる$A$が存在します（ポアンカレの補題）。以下では、そのようなポテンシャルを選べる領域で考え、大域的な条件には立ち入りません。

# ゲージ

## ゲージ変換

同じ$F$を与えるポテンシャルは1つではありません。スカラー関数$\chi$の勾配$D\chi$を加えても、$D\wedge D\chi=0$なので$D\wedge A$は変わりません。

&&&def ゲージ変換
$$
A\mapsto A+D\chi
$$
&&&

$D\chi=\gamma_0\partial_0\chi-\sum_k\gamma_k\partial_k\chi$なので、成分では$\varphi\mapsto\varphi+\partial_0\chi$、$c\boldsymbol A\mapsto c\boldsymbol A-\nabla\chi$です。一方、スカラー部は

$$
D\cdot(A+D\chi)=D\cdot A+\square\chi
$$

と変わります。$DA$のうち、2ベクトル部$D\wedge A$は物理的な場で、スカラー部$D\cdot A$はゲージ変換によって変わる量です。$A$の4成分のうち、$F$を決めるのは$D\wedge A$だけで、$D\cdot A$は$F$に現れません。

ゲージ変換で加える$D\chi$はスカラー関数の勾配で、$D\wedge D\chi=0$なので$F=D\wedge A$には現れません。3次元で$\boldsymbol B=\nabla\times\boldsymbol A$が$\boldsymbol A$に加えた勾配を見ないのと同じで、静的な場合にもベクトルポテンシャルに勾配を足す自由度はすでにありました。ゲージ変換で$DA$が変わるのはスカラー部$D\cdot A$だけなので、この自由度を使って$D\cdot A$を扱いやすい値、たとえば$0$に選べます。ただし、この条件だけでは、場に現れない自由度をすべて取り除けるわけではありません（下のローレンスゲージの節）。電位に定数を足しても電場が変わらないことの時空版で、足せるものが定数から勾配$D\chi$の形の関数全体に広がっています（$\chi=kx_0$とすれば$D\chi=k\gamma_0$で、電位に定数$k$を足す変換です）。

電荷に働く力は[[7shi-em1]]のローレンツ力で決まり、そこに現れるのは$\boldsymbol E$と$\boldsymbol B$だけです。したがって、ゲージ変換で結ばれた2つのポテンシャルは、同じ物理的状況を表します。ここではポテンシャルを、同じ場を複数の仕方で表せる計算上の変数として使い、ゲージは計算の都合に合わせて選びます。

## ローレンスゲージ

波動方程式$\square\chi=-D\cdot A$の解$\chi$を取ってゲージ変換すると、$D\cdot A=0$にできます。

&&&def ローレンスゲージ
$$
D\cdot A=\frac1c\partial_t\varphi+c\,\nabla\cdot\boldsymbol A=0
$$
&&&

名前は、この条件を用いたデンマークの物理学者ローレンス（L. Lorenz）に由来し、ローレンツ変換のローレンツ（H. A. Lorentz）とは別人です。日本語の文献ではどちらも「ローレンツ」と表記されることが多いですが、本シリーズでは両者を区別して、Lorenzをローレンス、Lorentzをローレンツと書きます。

このとき$DA$は2ベクトル部だけになり、$F=DA$です。ただし、$\square\chi=0$を満たす$\chi$によるゲージ変換はこの条件を保つので、ローレンスゲージでもポテンシャルが一意に決まるわけではありません。

&&&thm ポテンシャルの波動方程式
ローレンスゲージの下で、マクスウェル方程式$DF=\mu_0cJ$は

$$
\square A=\mu_0cJ
$$

と同値です。成分では

$$
\Bigl(\frac1{c^2}\partial_t^2-\Delta\Bigr)\varphi=\frac\rho{\varepsilon_0},\qquad
\Bigl(\frac1{c^2}\partial_t^2-\Delta\Bigr)\boldsymbol A=\mu_0\boldsymbol J
$$

です。
&&&

&&&prf
$F=DA$より$DF=D^2A=\square A$である。$\gamma_0$成分は$\square\varphi=\mu_0c\cdot c\rho=\rho/\varepsilon_0$、$\gamma_k$成分は$\square(cA_k)=\mu_0cJ_k$である。
&&&

ゲージを固定しない場合は、$DF=D(DA-D\cdot A)=\square A-D(D\cdot A)$です。ローレンスゲージを選ぶと、この第2項が消え、各成分についての波動方程式になります。

## 静的な場合

場と源が時間によらず、ポテンシャルも時間によらないものを選ぶと、ローレンス条件は$\nabla\cdot\boldsymbol A=0$です。[[7shi-em1]]では、ベクトルポテンシャルを積分で定めると定常電流に対して$\nabla\cdot\boldsymbol A=0$が自動的に成り立ち、それが$F$のスカラー部の消滅と同じことでした。これはローレンスゲージの静的な場合にあたります。このとき$F=DA=\bar{\mathcal D}P=-D_3P$となり、[[7shi-em1]]の$F=-DP$（$D$は空間のディラック作用素）に戻ります。上の波動方程式も、$\partial_t=0$として[[7shi-em1]]のポアソン方程式$\Delta\varphi=-\rho/\varepsilon_0$、$\Delta\boldsymbol A=-\mu_0\boldsymbol J$になります。

# エネルギー運動量

## 定義

ポテンシャルにはゲージの自由度がありますが、場のエネルギーや運動量はその選び方によらず決まらなければなりません。そこで、ポテンシャルではなく$F$そのものから作ります。

[[7shi-em3]]では、$\operatorname{Cl}_{3,0}(\mathbb R)$の反転$F^\dagger=\boldsymbol E-ic\boldsymbol B$を使って$\frac{\varepsilon_0}2FF^\dagger=u+\boldsymbol S/c$を作りました。時空代数では反転$\tilde F=-F$（2ベクトルの符号の反転）が自然な操作で、$F^\dagger$はそれを$\gamma_0$で挟んだ$F^\dagger=\gamma_0\tilde F\gamma_0=-\gamma_0F\gamma_0$です。$\gamma_0$を任意のベクトル$a$に替えて、次の線形写像を定めます。

&&&def エネルギー運動量
$$
T(a)=-\frac{\varepsilon_0}2FaF
$$
&&&

&&&prop エネルギー運動量の性質
(1) ベクトル$a$に対して$T(a)$はベクトルで、$a\cdot T(b)=b\cdot T(a)$が成り立ちます。

(2) $T(\gamma_0)\gamma_0=\frac{\varepsilon_0}2FF^\dagger=u+\boldsymbol S/c$、すなわち$T(\gamma_0)=u\gamma_0+\sum_k\frac{S_k}c\gamma_k$です。
&&&

&&&prf
(1) $FaF$は奇数個のベクトルの積の和なので、グレード1と3の成分だけを持つ。反転はグレード$k$の成分に$(-1)^{k(k-1)/2}$を掛け、グレード1は変えずグレード3の符号を変える。一方$\widetilde{FaF}=\tilde Fa\tilde F=(-F)a(-F)=FaF$なので、グレード3の成分は$0$である。対称性は、スカラー部の性質$\langle XY\rangle_0=\langle YX\rangle_0$を$X=aF$、$Y=bF$に使った$\langle aFbF\rangle_0=\langle bFaF\rangle_0$から従う。

(2) $\sigma_k$は$\gamma_0$と反可換、$i\sigma_k$は可換なので、$\gamma_0F\gamma_0=-\boldsymbol E+ic\boldsymbol B=-F^\dagger$である。したがって$T(\gamma_0)\gamma_0=-\frac{\varepsilon_0}2F(\gamma_0F\gamma_0)=\frac{\varepsilon_0}2FF^\dagger$となり、[[7shi-em3]]の結果を使えばよい。右から$\gamma_0$を掛けると、$\sigma_k\gamma_0=\gamma_k$より第2式を得る。
&&&

$T(\gamma_0)$は、時間方向$\gamma_0$を向いた観測者から見たエネルギー密度$u$と、運動量密度$\boldsymbol S/c^2$の$c$倍を成分として持つ時空のベクトルです。[[7shi-em3]]のパラベクトル$u+\boldsymbol S/c$は、このベクトルに$\gamma_0$を掛けて$\operatorname{Cl}_{3,0}(\mathbb R)$に移したものです。[[7shi-em5]]では、粒子の4元運動量を観測者の時間軸へ射影するとエネルギーが得られました。場では、$T$に観測者の時間方向$n$を入れたベクトル$T(n)$を同じ$n$へ射影した$n\cdot T(n)$が、その観測者の見るエネルギー密度です（$\gamma_0\cdot T(\gamma_0)=u$）。粒子と違い、$T(n)$を作る段階にも観測者が入ります。

$T$全体は通常のエネルギー運動量テンソル$T^{\mu\nu}=\gamma^\mu\cdot T(\gamma^\nu)$に対応します。$T(a)$の$a$は取り出す方向で、時間方向を入れると密度が、空間方向を入れると流れが出ます。エネルギーの流れは向きだけを指定すればよいので$\boldsymbol S$はベクトルでしたが、運動量はそれ自身が向きを持つので、どの向きの運動量がどの向きの面を通るかの2つを指定します。$T^{kl}$は、法線が$l$方向の面を通る$k$方向の運動量の流束で、マクスウェルの応力$\tau_{kl}$とは$T^{kl}=-\tau_{kl}$の関係にあります。

&&&fml エネルギー運動量テンソルの空間成分 [fml-stress]
$$
T^{kl}=-\tau_{kl},\qquad
\tau_{kl}=\varepsilon_0\Bigl(E_kE_l+c^2B_kB_l-\frac12\delta_{kl}\bigl(|\boldsymbol E|^2+c^2|\boldsymbol B|^2\bigr)\Bigr)
$$

ここで$\tau_{kl}$はマクスウェルの応力テンソルです。
&&&

## 保存則

&&&thm エネルギー運動量の保存
マクスウェル方程式$DF=\mu_0cJ$の解について

$$
\sum_\mu\partial_\mu T(\gamma^\mu)=-\frac1cF\cdot J
$$

が成り立ちます。ここで$F\cdot J=\frac12(FJ-JF)$で、

$$
\frac1cF\cdot J=\frac1c(\boldsymbol J\cdot\boldsymbol E)\,\gamma_0+\sum_k\bigl(\rho\boldsymbol E+\boldsymbol J\times\boldsymbol B\bigr)_k\gamma_k
$$

です。
&&&

&&&prf
ここでは非同次方程式だけでなく同次方程式も用い、$DF=D\cdot F+D\wedge F=\mu_0cJ$とする。右からの作用$\sum_\mu(\partial_\mu F)\gamma^\mu$を考える。$(\gamma^\mu\partial_\mu F)^\sim=(\partial_\mu\tilde F)\gamma^\mu=-(\partial_\mu F)\gamma^\mu$なので$\sum_\mu(\partial_\mu F)\gamma^\mu=-\widetilde{DF}=-\mu_0c\tilde J=-\mu_0cJ$である。積の微分から

$$
\sum_\mu\partial_\mu T(\gamma^\mu)=-\frac{\varepsilon_0}2\sum_\mu\bigl((\partial_\mu F)\gamma^\mu F+F\gamma^\mu\partial_\mu F\bigr)
=-\frac{\varepsilon_0\mu_0c}2(-JF+FJ)
$$

となり、$\varepsilon_0\mu_0c=1/c$から第1式を得る。第2式は、$F\cdot J$の成分を[[7shi-cla2]]の内積の規約に従って計算すれば得られる（[[7shi-em5]]のローレンツ力の計算で$U$を$J$に替えたものである）。
&&&

$\gamma_0$成分は$\partial_0u+\frac1c\nabla\cdot\boldsymbol S=-\frac1c\boldsymbol J\cdot\boldsymbol E$であり、[[7shi-em3]]のポインティングの定理です。$\gamma_k$成分は運動量の保存則で、[エネルギー運動量テンソルの空間成分](#fml-stress)を使うと

$$
\partial_t\frac{S_k}{c^2}-\sum_l\partial_l\tau_{kl}=-\bigl(\rho\boldsymbol E+\boldsymbol J\times\boldsymbol B\bigr)_k
$$

となります。場の運動量密度$\boldsymbol S/c^2$の時間変化率と運動量流束$-\tau_{ij}$の発散の和が、電荷と電流が受ける力の密度の符号を変えたものに等しいことを表します。[[7shi-em3]]で扱わなかった運動量の保存は、エネルギーの保存と同じ1本の式の別の成分として得られます。右辺の$\frac1cF\cdot J$は、[[7shi-em5]]のローレンツ力$\frac qcF\cdot U$で、点電荷$qU$を電流密度$J$に置き換えたものです。空間成分の$\rho\boldsymbol E+\boldsymbol J\times\boldsymbol B$は、[[7shi-em1]]で見たローレンツ力の密度です。[[7shi-em1]]で磁場の面が速度を面内で曲げた反対称積$\frac12\bigl((i\boldsymbol B)\boldsymbol v-\boldsymbol v(i\boldsymbol B)\bigr)$が、ここでは電流密度に働く磁気力$\boldsymbol J\times\boldsymbol B=\frac12\bigl((i\boldsymbol B)\boldsymbol J-\boldsymbol J(i\boldsymbol B)\bigr)$として$F\cdot J$に含まれています。

電荷が場から受け取る運動量は右辺に現れ、領域内の場の運動量は境界を通る運動量の流れによっても変わります。外部との受け渡しがない系では、物質と場を合わせた全運動量が保存され、電磁場は物体と同じく運動量を持って運ぶ存在です。

&&&ex 放射圧
静止した壁に$\hat{\boldsymbol k}$の向きの平面波が垂直に入射するとし、入射波のエネルギー密度を$u$とします。[[7shi-em3]]で見たとおり平面波では$\boldsymbol S=cu\,\hat{\boldsymbol k}$なので、運動量密度は$\boldsymbol S/c^2=(u/c)\hat{\boldsymbol k}$で、それが速さ$c$で壁に流れ込みます。進行方向の運動量が進行方向へ運ばれることは、平面波で$T^{kl}=u\hat k_k\hat k_l$となることに表れています。波が完全に吸収されるとき、壁が単位面積・単位時間あたりに受け取る、壁に垂直な運動量、すなわち壁が受ける圧力は$c\cdot u/c=u$です。完全に反射される場合は運動量の向きが反転するので、圧力は入射波のエネルギー密度の2倍の$2u$になります。[[7shi-em3]]の太陽光の例では、完全に吸収される場合の時間平均の圧力は$1.4\times10^3/(3.0\times10^8)\approx4.7\times10^{-6}\ \mathrm{Pa}$です。日常では感じられないほど小さい力ですが、宇宙空間の探査機の軌道には影響し、太陽光を大きな帆で受けて進む実験も行われています。
&&&

# まとめ

時空代数で、ポテンシャルとエネルギー運動量を扱いました。

- **ポテンシャル**：$A=\varphi\gamma_0+c\sum A_k\gamma_k$は$\gamma_0A=\varphi-c\boldsymbol A$を満たし、$DA$の2ベクトル部$D\wedge A$が$F$です。$D\wedge D=0$から、源を持たない2本の式は恒等式になります。
- **ゲージ**：$A\mapsto A+D\chi$は、$F$に見えない勾配の成分を足すだけの変換です。$D\wedge A$を変えず、スカラー部$D\cdot A$だけを$\square\chi$だけ変えます。$D\cdot A=0$（ローレンスゲージ）では$F=DA$、$\square A=\mu_0cJ$です。静的な場合は[[7shi-em1]]の$F=-DP$と$\nabla\cdot\boldsymbol A=0$に戻ります。
- **エネルギー運動量**：$T(a)=-\frac{\varepsilon_0}2FaF$はベクトルからベクトルへの対称な線形写像で、$T(\gamma_0)\gamma_0$が[[7shi-em3]]の$\frac{\varepsilon_0}2FF^\dagger$です。その発散$\sum\partial_\mu T(\gamma^\mu)=-\frac1cF\cdot J$の時間成分がポインティングの定理、空間成分が運動量の保存則です。
