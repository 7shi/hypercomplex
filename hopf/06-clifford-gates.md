前回まで、1・2・3量子ビットの状態空間の幾何をホップ写像で読み取ってきました。[[7shi-ent]] 今回は状態を動かす側、すなわち量子ゲートを扱います。全体位相が観測に現れないことから、1量子ビットのゲートは実質的にブロッホ球の回転であり、シリーズを通じて使ってきた回転子の共役作用にほかなりません。この見方に立つと、量子計算で基本的な役割を果たす**クリフォードゲート**が「パウリ行列の3軸を3軸に写す回転」として幾何的に特徴づけられ、位相を除けば正八面体の回転群とちょうど一致することを示します。さらに2量子ビットのCNOTゲートにも軸の行き先を追う見方がそのまま通用し、前回扱ったもつれの生成が軸の言葉で読めることを確認します。ブラケット記法は前回導入したものを引き続き使います。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

# 量子ゲートとブロッホ球の回転

1量子ビットの状態は、単位ベクトル$\Psi=(\alpha,\beta)^T$（$|\alpha|^2+|\beta|^2=1$）で表されるのでした。[[7shi-s]][[7shi-bd]] この状態を別の状態に変える操作が量子ゲートです。ノルムを保つ線形変換なので、ユニタリ行列で表されます。

&&&def 量子ゲート
状態ベクトルに左から掛けるユニタリ行列$U$（$U^\dagger U=\mathrm I$）を**量子ゲート**と呼びます。$n$量子ビットでは$2^n$次のユニタリ行列です。ゲートを順に施す列を**量子回路**と呼び、行列の積で表されます。
&&&

&&&rem 作用の順序
状態に$U_1,U_2$の順で施すと$U_2U_1\Psi$となり、式の上では右にあるゲートほど先に作用します。
&&&

観測量の期待値$\Psi^\dagger A\Psi$は全体位相$e^{i\varphi}\Psi$で変わらないため[[7shi-bd]]、ゲートも$U$と$e^{i\varphi}U$は区別できません。観測に効くのは、密度行列$\rho=\Psi\Psi^\dagger$の変換
$$
\rho\ \mapsto\ (U\Psi)(U\Psi)^\dagger=U\rho U^\dagger
$$
すなわち共役作用だけです。状態への片側作用が外積への共役作用を誘導するこの構図は、リー群シリーズでスピノルとベクトルの区別として確認したものです。[[7shi-spin]] $\rho=\frac12(\mathrm I+\boldsymbol r\cdot\boldsymbol\sigma)$[[7shi-bd]]より、共役作用はブロッホベクトル$\boldsymbol r$の変換として読めます。

その変換の正体が回転です。ユニタリ行列は位相を除けば行列式1にでき、$\mathrm{SU}(2)$の元はトレース0の反エルミート行列の指数関数、すなわちパウリ行列を軸とする回転子として書けるのでした。[[7shi-su2]]

&&&fml 量子ゲートの回転子分解
任意の1量子ビットゲートは、単位ベクトル$\boldsymbol n$と角度$\theta$、位相$\varphi$により
$$
U=e^{i\varphi}R_{\boldsymbol n}(\theta),\quad
R_{\boldsymbol n}(\theta)=\exp\Bigl(-\frac{i\theta}2\,\boldsymbol n\cdot\boldsymbol\sigma\Bigr)
=\cos\frac\theta2\,\mathrm I-i\sin\frac\theta2\,\boldsymbol n\cdot\boldsymbol\sigma
$$
と分解できます。共役作用$U(\boldsymbol r\cdot\boldsymbol\sigma)U^\dagger$は、ブロッホベクトルの軸$\boldsymbol n$周り角度$\theta$の回転（右ねじの向き）です。
&&&

共役作用が2倍角ではなく角度$\theta$そのものの回転になるのは、両側から半分ずつ作用するためです。鏡映2回の合成としての回転子の導入と指数関数表現は[[7shi-rot]]、パウリ行列による定式化は[[pauli-qua]]、$\mathrm{SU}(2)$と$\mathrm{SO}(3)$の対応としての導出は[[7shi-su2]]を参照してください。

&&&rem 回転子との記法の対応
パウリ行列による回転の記事[[pauli-qua]]では、四元数の回転子$r=\exp(\frac\theta2 n)$に対応する行列を$R=\exp(\frac\theta2 N)$（$N=i\,\boldsymbol n\cdot\boldsymbol\sigma\cong-n$）とおき、回転を$R^\dagger VR$と書きました。四元数の虚数単位との対応$\mathbf i\cong-i\sigma_x,\ \mathbf j\cong-i\sigma_y,\ \mathbf k\cong-i\sigma_z$[[7shi-bd]]をたどると$r\cong\exp(-\frac{i\theta}2\,\boldsymbol n\cdot\boldsymbol\sigma)$なので、本記事の$R_{\boldsymbol n}(\theta)=R^\dagger$は回転子$r$の行列表現そのものであり、共役作用$R_{\boldsymbol n}VR_{\boldsymbol n}^\dagger=R^\dagger VR$は$rvr^*$にあたります。
&&&

向きを1つの例で確認します。$z$軸周りの回転子は対角行列
$$
R_z(\theta)=\begin{pmatrix}e^{-i\theta/2}&0\\0&e^{i\theta/2}\end{pmatrix}
$$
で、$\sigma_x$を挟むと
$$
R_z(\theta)\,\sigma_x\,R_z(\theta)^\dagger
=\begin{pmatrix}0&e^{-i\theta}\\e^{i\theta}&0\end{pmatrix}
=\cos\theta\,\sigma_x+\sin\theta\,\sigma_y
$$
となります。$x$軸が$y$軸の方へ倒れる、右ねじの$+\theta$回転です。

&&&rem 二重被覆とゲート
$\theta=2\pi$では$R_{\boldsymbol n}(2\pi)=-\mathrm I$となり、状態ベクトルは符号を変えます（$4\pi$で元に戻る二重被覆[[7shi-su2]]）。ただしゲートとしては$-\mathrm I$は全体位相なので、$2\pi$回転は恒等ゲートと区別できません。つまり1量子ビットゲートの全体は、位相を除けば回転群$\mathrm{SO}(3)$です。
&&&

# 代表的な1量子ビットゲート

量子計算に登場する代表的なゲートを導入します。

&&&def 代表的な1量子ビットゲート
$$
X=\sigma_x,\quad Y=\sigma_y,\quad Z=\sigma_z
$$
$$
H=\frac1{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}=\frac{\sigma_x+\sigma_z}{\sqrt2},\quad
S=\begin{pmatrix}1&0\\0&i\end{pmatrix},\quad
T=\begin{pmatrix}1&0\\0&e^{i\pi/4}\end{pmatrix}
$$
$X,Y,Z$を**パウリゲート**、$H$を**アダマールゲート**、$S$を**位相ゲート**、$T$を**Tゲート**と呼びます。
&&&

基底に対する作用を見ると、それぞれの役割が分かります。$X$はビットを反転し（$X|0\rangle=|1\rangle,\ X|1\rangle=|0\rangle$）、$Z$は$|1\rangle$の符号だけを変えます。$H$は基底から重ね合わせを作ります。
$$
H|0\rangle=\frac{|0\rangle+|1\rangle}{\sqrt2},\quad
H|1\rangle=\frac{|0\rangle-|1\rangle}{\sqrt2}
$$
$S,T$は$|1\rangle$に位相$i,\ e^{i\pi/4}$を掛けます。$S=T^2,\ Z=S^2$という入れ子の関係があります。

これらを回転子として分解すると、次の表になります。

&&&fml 回転子分解の表
$$
\begin{aligned}
X&=i\,R_x(\pi) & Y&=i\,R_y(\pi) & Z&=i\,R_z(\pi)\\
H&=i\,R_{\boldsymbol n_H}(\pi)\ \ \Bigl(\boldsymbol n_H=\tfrac{(1,0,1)}{\sqrt2}\Bigr) &
S&=e^{i\pi/4}R_z\Bigl(\frac\pi2\Bigr) &
T&=e^{i\pi/8}R_z\Bigl(\frac\pi4\Bigr)
\end{aligned}
$$
ここで$R_x,R_y,R_z$は座標軸周りの回転子$R_{\boldsymbol n}(\theta)$です。
&&&

&&&ex 表の検算
$\theta=\pi$では$R_{\boldsymbol n}(\pi)=-i\,\boldsymbol n\cdot\boldsymbol\sigma$なので、$\boldsymbol n\cdot\boldsymbol\sigma=i\,R_{\boldsymbol n}(\pi)$です。$X,Y,Z$は座標軸方向、$H=\frac{\sigma_x+\sigma_z}{\sqrt2}$は$\boldsymbol n_H$方向の$\boldsymbol n\cdot\boldsymbol\sigma$なので、4つとも$\pi$回転の行です。$S$は対角行列なので直接計算できます。
$$
e^{i\pi/4}R_z\Bigl(\frac\pi2\Bigr)
=e^{i\pi/4}\begin{pmatrix}e^{-i\pi/4}&0\\0&e^{i\pi/4}\end{pmatrix}
=\begin{pmatrix}1&0\\0&i\end{pmatrix}=S
$$
$T$も同様です。
&&&

パウリゲートは「$x,y,z$軸周りの$180°$回転」、$H$は「$x$軸と$z$軸のちょうど中間の軸周りの$180°$回転」、$S,T$は「$z$軸周りの$90°,45°$回転」です。ビット反転や重ね合わせの生成といった作用が、すべてブロッホ球の回転として統一的に読めました。

## 観測量・ゲート・回転子の3役

表の上段の4つ（$X,Y,Z,H$）はどれも$\boldsymbol n\cdot\boldsymbol\sigma$の形をしています。この形の行列は、エルミートかつユニタリという二重の性格を持ちます。

&&&fml エルミートかつユニタリな行列
$2$次の行列$A$がエルミートかつユニタリで$A\neq\pm\mathrm I$であることと、単位ベクトル$\boldsymbol n$により$A=\boldsymbol n\cdot\boldsymbol\sigma$と書けることは同値です。このとき共役作用は
$$
(\boldsymbol n\cdot\boldsymbol\sigma)(\boldsymbol r\cdot\boldsymbol\sigma)(\boldsymbol n\cdot\boldsymbol\sigma)
=\{2(\boldsymbol n\cdot\boldsymbol r)\boldsymbol n-\boldsymbol r\}\cdot\boldsymbol\sigma
$$
すなわち$\boldsymbol n$周りの$\pi$回転です。
&&&

&&&prf
エルミート行列は$A=w\mathrm I+\boldsymbol r\cdot\boldsymbol\sigma$（$w,\boldsymbol r$は実）と展開できます。[[7shi-bd]] ユニタリ性は$A^2=\mathrm I$と同値で、パウリ行列の反交換関係$\sigma_a\sigma_b+\sigma_b\sigma_a=2\delta_{ab}\mathrm I$より
$$
A^2=(w^2+|\boldsymbol r|^2)\mathrm I+2w\,\boldsymbol r\cdot\boldsymbol\sigma=\mathrm I
$$
なので$w\boldsymbol r=0$です。$\boldsymbol r=0$なら$A=\pm\mathrm I$、そうでなければ$w=0,\ |\boldsymbol r|=1$となり$A=\boldsymbol n\cdot\boldsymbol\sigma$です。

挟み込みの公式は、$\boldsymbol r\cdot\boldsymbol\sigma$を$\boldsymbol n$方向と直交成分に分けます。
$$
\boldsymbol r\cdot\boldsymbol\sigma=(\boldsymbol n\cdot\boldsymbol r)(\boldsymbol n\cdot\boldsymbol\sigma)+\boldsymbol m\cdot\boldsymbol\sigma,\quad
\boldsymbol m=\boldsymbol r-(\boldsymbol n\cdot\boldsymbol r)\boldsymbol n\perp\boldsymbol n
$$
反交換関係を線形に組み合わせると$(\boldsymbol n\cdot\boldsymbol\sigma)(\boldsymbol m\cdot\boldsymbol\sigma)+(\boldsymbol m\cdot\boldsymbol\sigma)(\boldsymbol n\cdot\boldsymbol\sigma)=2(\boldsymbol n\cdot\boldsymbol m)\mathrm I$なので、直交する$\boldsymbol n\cdot\boldsymbol\sigma$と$\boldsymbol m\cdot\boldsymbol\sigma$は反交換し、また$(\boldsymbol n\cdot\boldsymbol\sigma)^2=|\boldsymbol n|^2\mathrm I=\mathrm I$です。よって
$$
(\boldsymbol n\cdot\boldsymbol\sigma)\{(\boldsymbol n\cdot\boldsymbol r)(\boldsymbol n\cdot\boldsymbol\sigma)+\boldsymbol m\cdot\boldsymbol\sigma\}(\boldsymbol n\cdot\boldsymbol\sigma)
=(\boldsymbol n\cdot\boldsymbol r)(\boldsymbol n\cdot\boldsymbol\sigma)-\boldsymbol m\cdot\boldsymbol\sigma
$$
となり、$\boldsymbol n$成分を保ち直交成分を反転する変換、すなわち$\boldsymbol n$周りの$\pi$回転です。
&&&

つまり同じ行列$\boldsymbol n\cdot\boldsymbol\sigma$が3つの役を演じます。

1. **観測量** — 期待値$\Psi^\dagger(\boldsymbol n\cdot\boldsymbol\sigma)\Psi$はブロッホベクトルの$\boldsymbol n$方向成分（$\boldsymbol n$方向の測定[[7shi-bd]]）
2. **ゲート** — ユニタリなので状態に掛けられる
3. **回転子** — 挟めば$\boldsymbol n$周りの$\pi$回転

$X,Y,Z,H$が自分自身の逆元であること（$A^2=\mathrm I$）は、幾何的には「$180°$回転は2回で元に戻る」という当たり前の事実です。

&&&rem 鏡映との関係
クリフォード代数では、単位ベクトル$n$による挟み込み$-nvn$が超平面鏡映でした。[[pauli-qua]] 符号を外した$+nvn=2(\boldsymbol n\cdot\boldsymbol v)\boldsymbol n-\boldsymbol v$は「$\boldsymbol n$を通る直線に関する鏡映」であり、3次元では$\boldsymbol n$周りの$\pi$回転と同じ変換です。ゲートの共役作用に現れるのは符号を外した後者で、行列式$+1$の回転として扱えます。
&&&

# パウリ群と四元数群

パウリゲートを掛け合わせると、$\sigma_x\sigma_y=i\sigma_z$のように$\pm1,\pm i$の係数が現れます。[[7shi-bd]] 積で閉じた集合にするには、これらの位相ごと含める必要があります。

&&&def パウリ群
$$
P_1=\{i^k\sigma\mid k=0,1,2,3,\ \sigma\in\{\mathrm I,X,Y,Z\}\}
$$
は行列の積で閉じ、位数16の群をなします。これを（1量子ビットの）**パウリ群**と呼びます。
&&&

共役作用で見ると位相$i^k$は打ち消されるため、$P_1$のブロッホ球への作用は4種類だけです。恒等と、3つの座標軸周りの$\pi$回転
$$
\begin{pmatrix}1&&\\&-1&\\&&-1\end{pmatrix},\quad
\begin{pmatrix}-1&&\\&1&\\&&-1\end{pmatrix},\quad
\begin{pmatrix}-1&&\\&-1&\\&&1\end{pmatrix}
$$
であり、この4つはクラインの四元群をなします。どの2つを合成しても残りの1つになる、最小の非自明な回転の群です。

&&&rem 四元数群
$P_1$のうち行列式1のもの、すなわち$\mathrm{SU}(2)$に属する元は
$$
P_1\cap\mathrm{SU}(2)=\{\pm\mathrm I,\ \pm i\sigma_x,\ \pm i\sigma_y,\ \pm i\sigma_z\}
$$
の8つです。対応$\mathbf i\cong-i\sigma_x,\ \mathbf j\cong-i\sigma_y,\ \mathbf k\cong-i\sigma_z$[[7shi-bd]]により、これは四元数の単位元$\{\pm1,\pm\mathbf i,\pm\mathbf j,\pm\mathbf k\}$のなす**四元数群**$Q_8$の行列表現です（対応で$i\sigma_x\cong-\mathbf i$と符号は入れ替わりますが、集合としては一致します）。パウリ群の芯にはやはり四元数がいました。
&&&

# クリフォードゲート

回転子の見方では、パウリ行列$\sigma_x,\sigma_y,\sigma_z$はブロッホ球の3つの直交軸です。ゲートの共役作用は軸を回すので、「パウリの軸をパウリの軸に写す」ゲートは特別な位置を占めます。これがクリフォードゲートです。

&&&def クリフォードゲート
共役作用でパウリ群を保つユニタリ行列
$$
UP_1U^\dagger=P_1
$$
を**クリフォードゲート**と呼びます。クリフォードゲートの全体は群をなし（パウリ群の**正規化群**）、**クリフォード群**と呼ばれます。
&&&

&&&rem クリフォード代数との関係
このクリフォード群は、本シリーズや関連記事で扱ってきたクリフォード代数と同じくW. K. クリフォードにちなみますが、指す対象は別物です。名称は、クリフォード代数の生成元がなす有限群とその正規化群を研究したBolt–Room–Wallの一連の仕事[[brw]]で「クリフォード群」の名が使われたことに遡るとされます。パウリ行列がクリフォード代数の生成元であること[[pauli-qua]]を思い出すと、パウリ群の正規化群にこの名前が付くのは自然な流れです。なおクリフォード代数の理論には、可逆元のなす「クリフォード群（リプシッツ群）」という別の群もあり、そちらとも異なる概念なので注意してください。
&&&

クリフォードゲートは、回転として次のように特徴づけられます。

&&&thm 1量子ビットのクリフォードゲート
1量子ビットゲート$U$について、次の3条件は同値です。
1. $U$はクリフォードゲート
2. $U$の共役作用（ブロッホ球の回転）は、6つの符号付き座標軸$\pm x,\pm y,\pm z$の置換
3. その回転は、6つの符号付き座標軸を頂点とする正八面体の回転群の元

位相$e^{i\varphi}$を除くと、クリフォードゲートはちょうど24個です。
&&&

&&&prf
1⇒2: $U\sigma_aU^\dagger\in P_1$は、エルミート性（$(UAU^\dagger)^\dagger=UA^\dagger U^\dagger$）とトレース（相似変換で不変）を保ちます。$P_1$の元$i^k\sigma$のうちエルミートなのは$k$が偶数の$\pm\sigma$に限られ、トレース0まで要求すると$\pm\mathrm I$も除かれて$\pm\sigma_b$だけが残ります。よって各座標軸の像は符号付きの座標軸です。

2⇔3: 6つの符号付き座標軸（単位球面上の6点）を頂点とする正多面体が正八面体であり、頂点集合を保つ回転の全体がその回転群です。

3⇒1: 軸の置換となる回転$M\in\mathrm{SO}(3)$に対し、それを共役作用として実現する$U$が位相を除いて一意に存在します（$\mathrm{SU}(2)\to\mathrm{SO}(3)$の全射性[[7shi-su2]]）。$U\sigma_aU^\dagger=\pm\sigma_b\in P_1$であり、$P_1$の一般の元$i^k\sigma$の像も$i^k(\pm\sigma_b)\in P_1$なので$UP_1U^\dagger=P_1$です。

個数: 正八面体の回転は「$x$軸の頂点の行き先が6通り、それを固定して残りを回す方法が4通り」で$6\times4=24$個です。
&&&

## HとSがすべてを生成する

$H$と$S$の共役作用を具体的に見ます。

&&&fml HとSの共役作用
$$
\begin{aligned}
SXS^\dagger&=Y,& SYS^\dagger&=-X,& SZS^\dagger&=Z\\
HXH&=Z,& HYH&=-Y,& HZH&=X
\end{aligned}
$$
ブロッホベクトルの成分では、$S$は$(x,y,z)\mapsto(-y,x,z)$（$z$軸周り$90°$）、$H$は$(x,y,z)\mapsto(z,-y,x)$（$\boldsymbol n_H$周り$180°$、$x$軸と$z$軸の交換）です。
&&&

どちらも軸を軸に写すので、$H,S$はクリフォードゲートです。しかもこの2つだけで、位相を除いた24個すべてが生成されます。実際、$S$は$z$軸の4回対称を与え、$HSH$の共役作用は$x$軸周りの$90°$回転（$y\mapsto z,\ z\mapsto-y$）で$x$軸の4回対称を与えます。直交する2本の4回軸があれば正八面体のすべての回転が合成できます。例えば$HS$の共役作用は$x\mapsto-y\mapsto z\mapsto x$と3軸を巡回させる位数3の回転（頂点ではなく面の中心を貫く3回軸）です。24個すべての列挙は検証コードで確認します。

&&&ex 正八面体の頂点と6つの状態
6つの頂点に対応する状態を挙げます。ブロッホベクトルが$\pm z$の状態は$|0\rangle,|1\rangle$、$\pm x$は$|\pm\rangle=\frac{|0\rangle\pm|1\rangle}{\sqrt2}$、$\pm y$は$|{\pm i}\rangle=\frac{|0\rangle\pm i|1\rangle}{\sqrt2}$です。[[7shi-s]] クリフォードゲートはこの6状態を（位相を除いて）互いに移し、$|0\rangle$に$H,S$を繰り返し施して得られる状態はちょうどこの6つです。例えば$H|0\rangle=|+\rangle$、$S|+\rangle=|{+i}\rangle$です。
&&&

&&&rem 対称群との同型
正八面体の回転群は4次対称群$S_4$と同型です。正八面体と双対な立方体で見ると、回転は4本の対角線（$(1,1,1)$方向など）を置換し、逆にこの置換から回転が1つに決まります。$HS$の位数3の回転は、1本の対角線を固定して残り3本を巡回させる置換にあたります。
&&&

# Tゲートと八面体対称性の破れ

$T$は$z$軸周りの$45°$回転でした。共役作用を計算すると
$$
TXT^\dagger=\frac{\sigma_x+\sigma_y}{\sqrt2}
$$
となり、$x$軸は$x$軸と$y$軸の中間方向へ写ります。これはパウリの軸ではないので、$T$は**クリフォードゲートではありません**。幾何的には、$45°$回転は正八面体の頂点を辺の中点の方向、つまり頂点でない点へ写します。正八面体の対称性からこぼれる最初のゲートです。

&&&rem 万能性
クリフォードゲートだけをいくら合成しても、位相を除けば24個の回転しか作れません。ところが$T$を加えると、$H$と$T$の生成する回転の集合は$\mathrm{SO}(3)$の中で稠密になり、任意の1量子ビットゲートを任意の精度で近似できます。さらに次節のCNOTを加えた$\{H,T,\mathrm{CNOT}\}$は、多量子ビットの任意のゲートを近似できる**万能ゲート集合**となります。[[nielsen-chuang]] 有限の対称性を破る1つのゲートが、有限群と連続群の隔たりを埋めます。
&&&

# CNOTと2量子ビットのパウリ群

2量子ビットの状態は$\Psi=(\alpha,\beta,\gamma,\delta)^T$、ゲートは4次のユニタリ行列です。[[7shi-ent]] 1量子ビットゲートのペア$U\otimes V$はテンソル積（クロネッカー積[[7shi-tensor]]）として作用しますが、それだけでは積状態が積状態にしか写らず、もつれを作れません。もつれを作る代表的なゲートがCNOTです。

&&&def CNOTゲート
$$
\mathrm{CNOT}=\begin{pmatrix}1&0&0&0\\0&1&0&0\\0&0&0&1\\0&0&1&0\end{pmatrix}
=\Pi_0\otimes\mathrm I+\Pi_1\otimes X,\quad
\Pi_0=\begin{pmatrix}1&0\\0&0\end{pmatrix},\ 
\Pi_1=\begin{pmatrix}0&0\\0&1\end{pmatrix}
$$
基底では$|00\rangle,|01\rangle$を固定し、$|10\rangle\leftrightarrow|11\rangle$を交換します。すなわち量子ビット$A$（制御）が$|1\rangle$のときだけ$B$（標的）に$X$を施す、**制御NOT**です。$\mathrm{CNOT}^2=\mathrm I$で自分自身の逆元です。
&&&

2量子ビットのパウリ群も、テンソル積で同様に定義されます。

&&&def 2量子ビットのパウリ群
$$
P_2=\{i^k\,\sigma_a\otimes\sigma_b\mid k=0,1,2,3,\ \sigma_a,\sigma_b\in\{\mathrm I,X,Y,Z\}\}
$$
は位数64の群です。位相を除いた非自明な元は$4\times4-1=15$個で、これが2量子ビットの「軸」の本数です。
&&&

軸が15本という数え方は、密度行列の展開から来ています。エルミート行列のパウリ基底展開[[7shi-bd]]はテンソル積でもそのまま通用し（$\operatorname{tr}\{(\sigma_a\otimes\sigma_b)(\sigma_c\otimes\sigma_d)\}=\operatorname{tr}(\sigma_a\sigma_c)\operatorname{tr}(\sigma_b\sigma_d)$）、2量子ビットの密度行列は
$$
\rho=\frac14\sum_{a,b}c_{ab}\,\sigma_a\otimes\sigma_b,\quad
c_{ab}=\operatorname{tr}\{\rho(\sigma_a\otimes\sigma_b)\}
$$
と展開できます（添字は$\mathrm I,x,y,z$を走り、$c_{00}=\operatorname{tr}\rho=1$）。残る15個の実係数$c_{ab}$がブロッホベクトルの2量子ビット版であり、ゲートの共役作用はこの15本の軸を回します。クリフォードゲートの定義も1量子ビットと同文で、$UP_2U^\dagger=P_2$です。

CNOTの共役作用を求めます。生成元の行き先だけ押さえれば、あとは積で決まります。

&&&fml CNOTの共役作用
$$
\begin{aligned}
X\otimes\mathrm I&\mapsto X\otimes X,&
Y\otimes\mathrm I&\mapsto Y\otimes X,&
Z\otimes\mathrm I&\mapsto Z\otimes\mathrm I\\
\mathrm I\otimes X&\mapsto \mathrm I\otimes X,&
\mathrm I\otimes Y&\mapsto Z\otimes Y,&
\mathrm I\otimes Z&\mapsto Z\otimes Z
\end{aligned}
$$
&&&

1行目の最初だけ計算を示します。$X\Pi_0=\Pi_1X,\ X\Pi_1=\Pi_0X$（$X$は$|0\rangle,|1\rangle$を入れ替える）と$\Pi_0\Pi_1=0,\ \Pi_a^2=\Pi_a$を使うと
$$
\begin{aligned}
\mathrm{CNOT}\,(X\otimes\mathrm I)\,\mathrm{CNOT}
&=(\Pi_0\otimes\mathrm I+\Pi_1\otimes X)(X\otimes\mathrm I)(\Pi_0\otimes\mathrm I+\Pi_1\otimes X)\\
&=(X\Pi_1\otimes\mathrm I+X\Pi_0\otimes X)(\Pi_0\otimes\mathrm I+\Pi_1\otimes X)\\
&=X\Pi_1\otimes X+X\Pi_0\otimes X=X\otimes X
\end{aligned}
$$
です。他の行も同様に計算でき（検証コードで全数確認します）、15本の軸はすべて符号付きの軸に写ります。**CNOTはクリフォードゲートです**。

## もつれの生成

表には、局所的なゲートには決してできないことが現れています。$U\otimes V$の共役作用は$\sigma_a\otimes\mathrm I\mapsto(U\sigma_aU^\dagger)\otimes\mathrm I$のように$A$側と$B$側の軸を別々に回すだけですが、CNOTは$A$側だけの軸$X\otimes\mathrm I$を両側にまたがる軸$X\otimes X$に写しています。

このことは、もつれの生成として確かめられます。前回の判定[[7shi-ent]]で使った行列$M=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}$は、局所ゲート$U\otimes V$のもとで$M\mapsto UMV^T$と変換するため、行列式は$\det U\det V$倍、その絶対値すなわちもつれの深さ$C=2|\alpha\delta-\beta\gamma|$は不変です。局所ゲートはもつれを作りも壊しもしません。一方CNOTは$C$を変えます。

&&&ex ベル状態の生成
$|00\rangle$に$H\otimes\mathrm I$、続いてCNOTを施します。
$$
|00\rangle
\ \xrightarrow{H\otimes\mathrm I}\ \frac{|00\rangle+|10\rangle}{\sqrt2}
\ \xrightarrow{\mathrm{CNOT}}\ \frac{|00\rangle+|11\rangle}{\sqrt2}=\Phi_+
$$
$\alpha\delta-\beta\gamma=\frac12\neq0$で、$C=1$の最大もつれです。前回の四元数ホップ像では、分離可能な状態の像がなす$S^2$の北極$|00\rangle$から、最大もつれの円上の$-\mathbf j$へ移ります。[[7shi-ent]] たった2つのゲートで、分離可能な世界からもつれの世界へ渡れました。
&&&

## 軸の行き先でベル状態を読む

同じ回路を、状態ではなく軸を追って読み直します。出発点の$|00\rangle$は、2つの観測量$Z\otimes\mathrm I$と$\mathrm I\otimes Z$の固有値$+1$の同時固有状態であり、この2条件で（位相を除いて）一意に決まります。一般に$A\Psi=\Psi$なら
$$
(UAU^\dagger)(U\Psi)=UA\Psi=U\Psi
$$
なので、回路$U$を通した状態は、軸の像$UAU^\dagger$の固有値$+1$の固有状態です。回路$U=\mathrm{CNOT}(H\otimes\mathrm I)$で2本の軸を追うと
$$
Z\otimes\mathrm I
\ \xrightarrow{H\otimes\mathrm I}\ X\otimes\mathrm I
\ \xrightarrow{\mathrm{CNOT}}\ X\otimes X,\qquad
\mathrm I\otimes Z
\ \xrightarrow{H\otimes\mathrm I}\ \mathrm I\otimes Z
\ \xrightarrow{\mathrm{CNOT}}\ Z\otimes Z
$$
となります。つまり$\Phi_+$は$X\otimes X$と$Z\otimes Z$の同時$+1$固有状態です。期待値で言えば$\langle X\otimes X\rangle=\langle Z\otimes Z\rangle=1$、一方で$A$単独の軸$\sigma\otimes\mathrm I$の期待値はすべて0（$\boldsymbol r_A=0$、最大混合[[7shi-ent]]）です。**個別の量子ビットの軸には何も見えないのに、相関の軸には完全な相関が見える**——前回、ホップ像とブロッホベクトルで見たもつれの姿が、軸の行き先を2本追うだけで再現されました。

&&&rem ゴッツマン＝クニルの定理
この「状態ではなく軸を追う」方法は$n$量子ビットに一般化できます。$|0\cdots0\rangle$は$n$本の軸$Z_i$（$i$番目だけ$Z$で他は$\mathrm I$のテンソル積）の同時$+1$固有状態で、クリフォード回路を通した後も、追跡すべきは$n$本の軸の行き先（それぞれ符号付きのパウリのテンソル積）だけです。記録量も更新の手間も$n$の多項式で済むため、**クリフォードゲートと計算基底での測定だけからなる量子回路は、古典計算機で効率よくシミュレートできます**（ゴッツマン＝クニルの定理[[gottesman]]）。もつれを作れるCNOTを含んでいても、軸の置換しか起こさない回路は古典で追えるのです。量子計算が古典計算を超えるには、$T$のようにクリフォード群からはみ出すゲートが欠かせません。
&&&

&&&rem quditへの一般化
パウリ行列$X=\sigma_x,\ Z=\sigma_z$は、クロック行列とシフト行列[[gpm]]の2次の場合でした。一般の次数$d$でも、ねじれた交換関係$ZX=\omega XZ$を満たす一般化パウリ行列$X^aZ^b$は位相を込めて群をなし、量子情報では$d$準位系（**qudit**）のパウリ群（ワイル＝ハイゼンベルク群）と呼ばれます。その正規化群としてquditのクリフォード群が定義され、本記事はちょうど$d=2$の場合です。クロック行列とシフト行列を「位置に応じた位相掛け」と「平行移動」として読む見方[[weyl-alg]]も、quditの言葉にそのまま引き継がれます。
&&&

# まとめ

これまで状態空間の幾何を追ってきたホップシリーズの延長として、状態を動かす量子ゲートを回転子として読み直しました。量子計算の基本語彙が、ブロッホ球の回転の幾何としてそのまま理解できます。

- 1量子ビットゲートは位相×回転子$e^{i\varphi}R_{\boldsymbol n}(\theta)$に分解され、共役作用はブロッホベクトルの軸$\boldsymbol n$周り角$\theta$の回転。位相を除けばゲートの全体は$\mathrm{SO}(3)$
- エルミートかつユニタリな$\boldsymbol n\cdot\boldsymbol\sigma$（$X,Y,Z,H$）は、観測量・ゲート・$\pi$回転の3役を兼ねる
- パウリ群は位相を除けばクラインの四元群、$\mathrm{SU}(2)$に属する代表は四元数群$Q_8$
- クリフォードゲート（パウリ群の正規化群）は、位相を除けば正八面体の回転群（位数24、$S_4$と同型）と一致し、$H$と$S$が生成する。頂点は6つの基本状態$|0\rangle,|1\rangle,|\pm\rangle,|{\pm i}\rangle$
- $T$（$45°$回転）は頂点を頂点でない点に写す非クリフォードゲートで、加えると万能になる
- CNOTは2量子ビットの15本の軸の符号付き置換だが局所ゲートではなく、軸の行き先を2本追うだけでベル状態の生成が読める
- 軸の追跡だけで済むことが、クリフォード回路の古典シミュレート可能性（ゴッツマン＝クニルの定理）の内容

&&&rem 参考文献
ゴッツマン＝クニルの定理と軸（スタビライザー）による記述はGottesman[[gottesman]]によります。「クリフォード群」の名称はBolt–Room–Wall[[brw]]の研究に遡るとされます。量子ゲートの標準的な事項はNielsenとChuangの教科書[[nielsen-chuang]]を参照してください。本記事の導出は本シリーズの範囲で自己完結するよう構成したものです。
&&&
