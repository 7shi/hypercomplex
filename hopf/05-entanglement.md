前回、ホップ写像を四元数による$S^7\to S^4$、八元数による$S^{15}\to S^8$に拡張しました。[[7shi-hopfext]] 今回はこの拡張を量子情報の側から読み直します。2量子ビットの状態空間が$S^7$、3量子ビットの状態空間が$S^{15}$であり、もつれの有無がホップ像の成分の消失として読み取れることを示します。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

# 2量子ビットの状態空間

これまで扱ってきた状態ベクトル$\Psi=(\alpha,\beta)^T$（$|\alpha|^2+|\beta|^2=1$）は、1つの量子ビットの状態を表します。[[7shi-s]][[7shi-bloch]] 2つの量子ビット$A,B$からなる合成系の状態空間は、それぞれの状態空間のテンソル積で与えられます。列ベクトルに対するテンソル積は、成分の総当たりの積（クロネッカー積）です。[[7shi-tp]]

&&&def 2量子ビットの状態ベクトル
$\Phi_A=(a,b)^T$と$\Phi_B=(c,d)^T$のテンソル積を
$$
\Phi_A\otimes\Phi_B=(ac,\ ad,\ bc,\ bd)^T
$$
と定めます。一般の2量子ビットの状態ベクトルは、この形とは限らない$\mathbb{C}^4$の単位ベクトル
$$
\Psi=(\alpha,\beta,\gamma,\delta)^T,\quad |\alpha|^2+|\beta|^2+|\gamma|^2+|\delta|^2=1
$$
です。実8成分と見なせば、$\Psi$は$S^7$上の点です。
&&&

&&&rem ブラケット記法
量子情報では基底ベクトルを$|0\rangle=(1,0)^T,\ |1\rangle=(0,1)^T$と書き、テンソル積の基底を並べて
$$
|00\rangle=(1,0,0,0)^T,\quad
|01\rangle=(0,1,0,0)^T,\quad
|10\rangle=(0,0,1,0)^T,\quad
|11\rangle=(0,0,0,1)^T
$$
と略記します（数字は$A,B$の順）。この記法では$\Psi=\alpha|00\rangle+\beta|01\rangle+\gamma|10\rangle+\delta|11\rangle$です。以下、具体的な状態を書くときに使います。
&&&

&&&def 分離可能ともつれ
$\Psi=\Phi_A\otimes\Phi_B$と書ける状態を**分離可能**（積状態）、書けない状態を**もつれた状態**と呼びます。
&&&

&&&rem 分解の非一意性
分離可能な状態の分解は位相の分だけ非一意です。$(e^{i\theta}\Phi_A)\otimes(e^{-i\theta}\Phi_B)=\Phi_A\otimes\Phi_B$となるため、$\Phi_A,\Phi_B$はそれぞれ位相を除いて定まります。
&&&

# 分離可能性の行列式判定

分離可能かどうかは、成分を$2\times2$行列に並べると判定できます。$\Psi$の成分を、行を$A$の基底、列を$B$の基底として
$$
M=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}
$$
と並べます。積状態$\Phi_A\otimes\Phi_B$の成分$(ac,ad,bc,bd)$を並べると
$$
M=\begin{pmatrix}ac&ad\\bc&bd\end{pmatrix}=\Phi_A\Phi_B^T
$$
となり、$M$の階数は1です。逆に階数1なら列ベクトル2本の積に分解できるため、分離可能性は$M$の階数、すなわち行列式で判定できます。

&&&thm 分離可能性の判定
$\Psi=(\alpha,\beta,\gamma,\delta)^T$が分離可能であることと、$\alpha\delta-\beta\gamma=0$は同値です。
&&&

&&&prf
分離可能なら$\alpha\delta-\beta\gamma=(ac)(bd)-(ad)(bc)=0$です。

逆に$\alpha\delta-\beta\gamma=0$とします。$\Psi\neq0$なので$M$のどちらかの行は非零です。第1行$(\alpha,\beta)$が非零なら、行列式の消失は第2行が第1行の複素数倍$(\gamma,\delta)=\lambda(\alpha,\beta)$であることを意味します（$\alpha\neq0$なら$\lambda=\gamma/\alpha$、$\beta\neq0$なら$\lambda=\delta/\beta$で、両方非零なら一致します）。よって
$$
\Psi=(1,\lambda)^T\otimes(\alpha,\beta)^T
$$
と分解でき、規格化すれば積状態です。第2行が非零の場合も同様です。
&&&

&&&ex ベル状態
もつれた状態の代表例が**ベル状態**です。
$$
\Phi_\pm=\frac{|00\rangle\pm|11\rangle}{\sqrt2},\quad
\Psi_\pm=\frac{|01\rangle\pm|10\rangle}{\sqrt2}
$$
例えば$\Phi_+$は$\alpha=\delta=\frac1{\sqrt2}$より$\alpha\delta-\beta\gamma=\frac12\neq0$で、もつれています。一方
$$
\frac{|00\rangle+|01\rangle+|10\rangle+|11\rangle}{2}=\frac{(1,1)^T}{\sqrt2}\otimes\frac{(1,1)^T}{\sqrt2}
$$
は$\alpha\delta-\beta\gamma=\frac14-\frac14=0$で、分離可能です。
&&&

# 四元数への埋め込み

$\Psi$の実8成分が$S^7$をなすことは、前回の四元数ホップ写像$H_{\mathbb{H}}$の定義域と同じです。そこで$\mathbb{C}^4$を四元数のペア$\mathbb{H}^2$と同一視します。

四元数は複素数2つに分けられます。$q=w+x\mathbf i+y\mathbf j+z\mathbf k$は、$\mathbf k=\mathbf{ij}$より
$$
q=(w+x\mathbf i)+(y+z\mathbf i)\mathbf j
$$
と書けるため、$\mathbf i$を虚数単位$i$と見なせば$\mathbb{H}=\mathbb{C}+\mathbb{C}\mathbf j$です。このとき$\mathbf{ji}=-\mathbf{ij}$から、複素数$z$と$\mathbf j$の交換則が従います。
$$
\mathbf jz=z^*\mathbf j,\quad (z\mathbf j)^*=-z\mathbf j
$$

&&&def 四元数への詰め込み
2量子ビットの状態ベクトル$\Psi=(\alpha,\beta,\gamma,\delta)^T$に対し、成分を2つずつ四元数に詰めます。
$$
q_1=\alpha+\beta\mathbf j,\quad q_2=\gamma+\delta\mathbf j
$$
$|q_1|^2+|q_2|^2=|\alpha|^2+|\beta|^2+|\gamma|^2+|\delta|^2=1$より、$(q_1,q_2)$は四元数ホップ写像$H_{\mathbb{H}}$の定義域$S^7$の点です。
&&&

&&&rem 前回の記法との対応
前回[[7shi-hopfext]]は$H_{\mathbb{H}}$の引数を$(\alpha,\beta)\in\mathbb{H}^2$と書きましたが、本記事では複素成分に$\alpha,\beta,\gamma,\delta$を使うため、四元数のペアを$(q_1,q_2)$と書きます。
&&&

&&&rem 係数を置く向き
$\mathbf j$は複素数と可換ではない（$\mathbf j\beta=\beta^*\mathbf j$）ため、係数を$\mathbf j$の左に置くか右に置くかで結果が変わります。逆向きに$q=\alpha+\mathbf j\beta$と詰めると、後で計算する$q_1q_2^*$の$\mathbf j$係数は$\beta^*\gamma-\alpha\delta^*$となり、分離可能でも消えません。実際、分離可能な状態
$$
\frac{(1,i)^T}{\sqrt2}\otimes\frac{(1,1)^T}{\sqrt2}=\frac{(1,\ 1,\ i,\ i)^T}{2}
$$
では$\beta^*\gamma-\alpha\delta^*=\frac i4+\frac i4=\frac i2\neq0$です。以下の結果は、係数を左に置く詰め方で成り立ちます。
&&&

&&&rem 初回の埋め込みとの違い
初回[[7shi-h]]では、1量子ビットの$(\alpha,\beta)$を1つの四元数$\omega$に埋め込みました（複素数の虚数単位に$\mathbf k$を使用）。今回は2量子ビットの4つの複素成分を四元数2つに詰めるという、別の埋め込みです。
&&&

# 四元数ホップ像ともつれ

$(q_1,q_2)$のホップ像を成分で計算します。交換則$\mathbf j\gamma^*=\gamma\mathbf j$と$\mathbf j\delta\mathbf j=\delta^*\mathbf j^2=-\delta^*$を使います。
$$
\begin{aligned}
q_1q_2^*&=(\alpha+\beta\mathbf j)(\gamma^*-\delta\mathbf j)\\
&=\alpha\gamma^*-\alpha\delta\mathbf j+\beta\mathbf j\gamma^*-\beta\mathbf j\delta\mathbf j\\
&=\alpha\gamma^*-\alpha\delta\mathbf j+\beta\gamma\mathbf j+\beta\delta^*\\
&=(\alpha\gamma^*+\beta\delta^*)+(\beta\gamma-\alpha\delta)\mathbf j
\end{aligned}
$$
したがってホップ像は
$$
H_{\mathbb{H}}(q_1,q_2)
=\bigl(\underbrace{2(\alpha\gamma^*+\beta\delta^*)}_{1,\mathbf i\text{ 成分}}+\underbrace{2(\beta\gamma-\alpha\delta)\mathbf j}_{\mathbf j,\mathbf k\text{ 成分}},\ \underbrace{|\alpha|^2+|\beta|^2-|\gamma|^2-|\delta|^2}_{\text{実成分}}\bigr)\in S^4
$$
となります。$\mathbf j,\mathbf k$成分の係数として、分離可能性を判定する行列式$\alpha\delta-\beta\gamma$がそのまま現れました。

&&&thm もつれとホップ像
2量子ビットの状態$\Psi$について、次の3条件は同値です。
1. $\Psi$は分離可能
2. $\alpha\delta-\beta\gamma=0$
3. ホップ像の$\mathbf j,\mathbf k$成分が消える
&&&

1と2の同値は前節の定理、2と3の同値は上の成分表示から明らかです。条件3は、$S^4$のうち$(1,\mathbf i,\text{実})$の張る3次元部分空間との交わり$S^2$に像が載ることを意味します。

分離可能な場合に像がどこに写るかを確認します。$\Psi=\Phi_A\otimes\Phi_B$（$\Phi_A=(a,b)^T,\ \Phi_B=(c,d)^T$）のとき、$B$の状態を四元数1つに詰めた$\psi=c+d\mathbf j$（$|\psi|=1$）を使うと
$$
q_1=ac+ad\mathbf j=a\psi,\quad q_2=b\psi
$$
となり、結合法則から
$$
q_1q_2^*=(a\psi)(\psi^*b^*)=a(\psi\psi^*)b^*=ab^*
$$
です。よって像は
$$
H_{\mathbb{H}}(q_1,q_2)=(2ab^*,\ |a|^2-|b|^2)=H(a,b)
$$
となります。これは初回[[7shi-h]]の複素数ホップ写像による、量子ビット$A$のブロッホ球そのものです。つまり分離可能な状態の像がなす$S^2\subset S^4$は、$S^3\to S^2$のホップファイブレーションの底空間が入れ子で再現されたものであり、$B$の状態$\psi$は像から消えます。

# 観測量とブロッホベクトル

像の成分が観測量とどう対応するかを確認します。量子ビット$A$だけを測定するときの観測量は、パウリ行列[[7shi-bloch]]を$A$側に置いたテンソル積$\sigma\otimes\mathrm I$です。例えば
$$
\sigma_x\otimes\mathrm I=\begin{pmatrix}0&0&1&0\\0&0&0&1\\1&0&0&0\\0&1&0&0\end{pmatrix}
$$
で、$\Psi=(\alpha,\beta,\gamma,\delta)^T$の前半と後半を入れ替えます。期待値を計算すると
$$
\begin{aligned}
x_A&=\Psi^\dagger(\sigma_x\otimes\mathrm I)\Psi=\alpha^*\gamma+\beta^*\delta+\gamma^*\alpha+\delta^*\beta\\
y_A&=\Psi^\dagger(\sigma_y\otimes\mathrm I)\Psi=-i(\alpha^*\gamma+\beta^*\delta-\gamma^*\alpha-\delta^*\beta)\\
z_A&=\Psi^\dagger(\sigma_z\otimes\mathrm I)\Psi=|\alpha|^2+|\beta|^2-|\gamma|^2-|\delta|^2
\end{aligned}
$$
となり、ホップ像の成分と比べると
$$
2(\alpha\gamma^*+\beta\delta^*)=x_A-iy_A,\quad
|q_1|^2-|q_2|^2=z_A
$$
です。初回の1量子ビットでも像の複素部分は$2\alpha\beta^*=x-iy$でした。[[7shi-h]] つまりホップ像の複素部分と実成分は、$A$のブロッホベクトル$\boldsymbol r_A=(x_A,y_A,z_A)$です。

残る$\mathbf j,\mathbf k$成分の大きさを
$$
C=2|\alpha\delta-\beta\gamma|
$$
とおきます。像のノルムが1であること（ノルムの乗法性[[7shi-hopfext]]）は、次の関係を与えます。

&&&fml ブロッホベクトルともつれの内訳
$$
x_A^2+y_A^2+z_A^2+C^2=1
$$
&&&

分離可能なら$C=0$で$|\boldsymbol r_A|=1$、すなわち$A$のブロッホベクトルは球面上にあります。もつれる（$C>0$）と$|\boldsymbol r_A|<1$となり、ブロッホベクトルは球の内部に落ちます。前々回[[7shi-bloch]]見たように、球の内部は混合状態です。$A$だけの測定期待値はすべて$\boldsymbol r_A$で決まるため、$A$単独の状態は密度行列$\rho_A=\frac12(\mathrm I+\boldsymbol r_A\cdot\boldsymbol\sigma)$で記述され、その純粋度は
$$
\operatorname{tr}\rho_A^2=\frac{1+|\boldsymbol r_A|^2}2=1-\frac{C^2}2
$$
です。合成系$\Psi$は純粋状態なのに、部分だけを見ると混合状態になる——これがもつれの測定における現れです。

&&&rem 量子情報の用語
$C$は量子情報で**コンカレンス**と呼ばれるもつれの尺度です。また$\rho_A$を得る操作は**部分トレース**、もつれを標準形に直す分解は**シュミット分解**と呼ばれます。本記事ではこれらの一般論には立ち入らず、ホップ像から直接読み取れる範囲を扱います。
&&&

&&&rem 全体位相
$e^{i\varphi}\Psi$は$q_1,q_2$への複素数$e^{i\varphi}$の左乗算です。$q_1q_2^*\mapsto e^{i\varphi}(q_1q_2^*)e^{-i\varphi}$となり、複素部分は不変ですが、$\mathbf j,\mathbf k$成分$w\mathbf j$は$e^{i\varphi}w\mathbf je^{-i\varphi}=e^{2i\varphi}w\mathbf j$と角度$2\varphi$だけ回転します。1量子ビットではホップ像自体が全体位相で不変でしたが、2量子ビットでは像が$(\mathbf j,\mathbf k)$平面内で動きます。ただし測定に現れる$\boldsymbol r_A$と$C$は不変です。
&&&

# 具体例

いくつかの状態のホップ像を挙げます。像の座標は$(1,\mathbf i,\mathbf j,\mathbf k;\ \text{実成分})$の順です。
$$
\begin{aligned}
\frac{|00\rangle+|01\rangle+|10\rangle+|11\rangle}2\ &\mapsto\ (1,0,0,0;\ 0) && C=0\\
\Phi_+=\frac{|00\rangle+|11\rangle}{\sqrt2}\ &\mapsto\ -\mathbf j && C=1\\
\Phi_-=\frac{|00\rangle-|11\rangle}{\sqrt2}\ &\mapsto\ +\mathbf j && C=1\\
\Psi_+=\frac{|01\rangle+|10\rangle}{\sqrt2}\ &\mapsto\ +\mathbf j && C=1\\
\Psi_-=\frac{|01\rangle-|10\rangle}{\sqrt2}\ &\mapsto\ -\mathbf j && C=1\\
\frac{|00\rangle+i|11\rangle}{\sqrt2}\ &\mapsto\ -\mathbf k && C=1
\end{aligned}
$$
1つ目は分離可能な状態で、像は$A$のブロッホベクトル（$x$方向）です。ベル状態はいずれも$\boldsymbol r_A=0$（最大混合）で、像は$(\mathbf j,\mathbf k)$平面の単位円上にあります。この円は、$C=1$の状態（最大もつれ）の像の全体です。

両者をつなぐのが子午線です。
$$
\cos t\,|00\rangle+\sin t\,|11\rangle\ \mapsto\ (0,0,-\sin2t,0;\ \cos2t),\quad C=\sin2t
$$
$t=0$では北極（分離可能な$|00\rangle$）、$t=\frac\pi4$では$-\mathbf j$（最大もつれの$\Phi_+$）に至ります。もつれの深さ$C$は、像が分離可能状態のなす$S^2$からどれだけ離れているかを表す量です。

# ファイバーの解釈

初回[[7shi-h]]のファイバー$S^1$は全体位相でした。四元数ホップ写像のファイバーは、結合法則により右からの共通乗算$(q_1q,\ q_2q)$（$q\in S^3$）で保たれます。[[7shi-hopfext]] この$S^3$が量子ビットの言葉で何にあたるかを確認します。

単位四元数$q$は$\mathbb{H}=\mathbb{C}+\mathbb{C}\mathbf j$の分解で複素数のペアと同一視でき、これは1量子ビットの状態ベクトルと同じ形です。分離可能な状態$q_1=a\psi,\ q_2=b\psi$に右から$q$を掛けると、結合法則により
$$
q_1q=a(\psi q),\quad q_2q=b(\psi q)
$$
となります。$\psi q$も単位四元数なので、これは$B$の状態を$\psi$から$\psi q$に取り替えた分離可能状態です。$A$の状態$(a,b)$は変わりません。

つまり、**ファイバーに沿う移動は$B$の状態の取り替え**です。1量子ビットの「ファイバー$S^1$＝全体位相」は、2量子ビットでは「ファイバー$S^3$＝$B$の状態空間（位相込み）」に一般化されます。分離可能な像点$S^2$上のファイバーは「$A$のブロッホベクトルを固定した分離可能状態の全体」であり、像から消えた$B$の状態がファイバー方向に格納されています。

もつれた像点でもファイバーは$S^3$で、$C$は像だけで決まるためファイバー上で一定です。例えば$\Phi_+$と$\Psi_-$はどちらも$-\mathbf j$に写り、実際$q=\mathbf j$の右乗算で移り合います。

# 3量子ビットと八元数

同じ構成を1段階上げます。3量子ビット$A,B,C$の状態空間は$\mathbb{C}^2\otimes\mathbb{C}^2\otimes\mathbb{C}^2\cong\mathbb{C}^8$で、状態ベクトルは実16成分と見なせば$S^{15}$上の点、すなわち八元数ホップ写像$H_{\mathbb{O}}$の定義域です。

$\Psi\in\mathbb{C}^8$を$A$の基底でブロックに分けます。
$$
\Psi=|0\rangle\otimes\Psi_0+|1\rangle\otimes\Psi_1,\quad \Psi_0,\Psi_1\in\mathbb{C}^4
$$
$\Psi_0,\Psi_1$は$BC$の2量子ビット分の成分です。

&&&def A|BC分離可能
$\Psi=\Phi_A\otimes X$（$X\in\mathbb{C}^4$は$BC$の合成状態）と書ける状態を**A|BC分離可能**と呼びます。これはブロックが複素数倍$\Psi_1=\lambda\Psi_0$（または$\Psi_0=0$）であることと同値です（2量子ビットの行列式判定と同じ論法を$2\times4$行列に適用）。$BC$の内部にもつれがあってもかまいません。
&&&

次に$\mathbb{C}^4\cong\mathbb{O}$の同一視を作ります。四元数のときは$\mathbb{H}=\mathbb{C}+\mathbb{C}\mathbf j$という分解を使いましたが、八元数では虚数単位$e_1$の左乗算を複素数の$i$と見なし、$\mathbb{C}$上の基底$\{1,e_2,e_4,e_6\}$への係数として定義します。

&&&def 八元数への詰め込み
$X=(z_1,z_2,z_3,z_4)^T\in\mathbb{C}^4$を八元数
$$
x=z_1\cdot1+z_2e_2+z_3e_4+z_4e_6
$$
に対応させます。ここで複素係数$z=\operatorname{Re}z+\operatorname{Im}z\,e_1$は左から掛けます。三つ組$123,145,176$[[7shi-oct1]]より$e_1e_2=e_3,\ e_1e_4=e_5,\ e_1e_6=-e_7$なので、成分で書くと
$$
x=\operatorname{Re}z_1+\operatorname{Im}z_1\,e_1
+\operatorname{Re}z_2\,e_2+\operatorname{Im}z_2\,e_3
+\operatorname{Re}z_3\,e_4+\operatorname{Im}z_3\,e_5
+\operatorname{Re}z_4\,e_6-\operatorname{Im}z_4\,e_7
$$
です。8つの実成分がちょうど1回ずつ現れるため、これは実線形同型です。
&&&

&&&rem 符号の注意
三つ組$176$は$e_1e_7=e_6$を意味するため、$e_1e_6=-e_7$と符号が反転します。$z_4$の虚部だけ$-e_7$に入るのはこのためです。
&&&

この対応が複素数の左乗算と両立すること（$\mathbb{C}$-線形性）を確認します。各基底$u\in\{1,e_2,e_4,e_6\}$について、$\{1,e_1,u,e_1u\}$は1つの三つ組で閉じる四元数部分代数の基底なので、その内部では結合法則が成り立ち
$$
a(zu)=(az)u\quad(a,z\in\mathbb{C})
$$
です。したがって$X\mapsto x$は$\mathbb{C}$-線形同型$\mathbb{C}^4\cong\mathbb{O}$であり、ブロックの複素数倍$\Psi_1=\lambda\Psi_0$は八元数の複素数倍$o_2=\lambda o_1$に移ります。

&&&rem 入れ子で定義してはいけない
四元数の$\mathbb{H}=\mathbb{C}+\mathbb{C}\mathbf j$の類推で、八元数を$(z_1+z_2e_2)+(z_3+z_4e_2)e_4$と入れ子（ケイリー＝ディクソン構成の形）で定義すると、結果が上と食い違います。$z=x+ye_1$に対して
$$
(ze_2)e_4=(xe_2+ye_3)e_4=xe_6+ye_7=z^*e_6
$$
となり（三つ組$347$より$e_3e_4=e_7$）、非結合性のために共役が混入して$\mathbb{C}$-線形になりません。基底への係数として直接定義するのはこれを避けるためです。
&&&

ブロックを八元数に詰めた$(o_1,o_2)$は$|o_1|^2+|o_2|^2=1$を満たし、$H_{\mathbb{O}}$の定義域$S^{15}$の点となります。

# 八元数ホップ像とA|BCもつれ

2量子ビットの定理を八元数で再現します。今回は成分計算の代わりに、前回[[7shi-hopfext]]使った合成代数の恒等式$x(x^*y)=|x|^2y,\ x^*(xy)=|x|^2y$を使います。両辺の共役を取って$y^*$を置き直すと、鏡像の形が得られます。
$$
(zx)x^*=|x|^2z,\quad (zx^*)x=|x|^2z
$$
これらは結合法則を使わずに成り立つ恒等式です。

&&&thm A|BCもつれとホップ像
3量子ビットの状態$\Psi$について、次の3条件は同値です。
1. $\Psi$はA|BC分離可能
2. $o_1o_2^*\in\mathbb{C}$（$1,e_1$成分のみ）
3. ホップ像$H_{\mathbb{O}}(o_1,o_2)$の$e_2,\dots,e_7$成分が消える
&&&

&&&prf
2と3の同値は、像の八元数部分が$2o_1o_2^*$そのものであることから明らかです。

1$\Rightarrow$2: $\Psi=\Phi_A\otimes X$（$\Phi_A=(a,b)^T$）なら、$X$を詰めた八元数$x$（$|x|=1$）を使って$o_1=ax,\ o_2=bx$です。$a\neq0$として$\lambda=b/a$とおくと$o_2=\lambda o_1$であり、恒等式$x(x^*y)=|x|^2y$より
$$
o_1o_2^*=o_1(\lambda o_1)^*=o_1(o_1^*\lambda^*)=|o_1|^2\lambda^*=ab^*\in\mathbb{C}
$$
となります（$|o_1|^2=|a|^2$）。$a=0$なら$o_1=0$で自明です。

2$\Rightarrow$1: $o_1o_2^*=c\in\mathbb{C}$とします。$o_2\neq0$なら、右から$o_2$を掛けて恒等式$(zx^*)x=|x|^2z$を$z=o_1,\ x=o_2$に使うと
$$
co_2=(o_1o_2^*)o_2=|o_2|^2o_1
$$
すなわち$o_1=\dfrac{c}{|o_2|^2}\,o_2$で、ブロックが複素数倍となりA|BC分離可能です。$o_2=0$なら$\Psi=|0\rangle\otimes\Psi_0$で自明です。
&&&

この証明は結合法則を一切使っていません。2量子ビットの成分計算では$\mathbf j,\mathbf k$の2実成分（1複素条件$\alpha\delta-\beta\gamma=0$）が消えましたが、八元数では$e_2,\dots,e_7$の6実成分（3複素条件）の消失がブロックの比例と同値になります。

分離可能な場合の像は、証明中の計算から
$$
H_{\mathbb{O}}(o_1,o_2)=(2ab^*,\ |a|^2-|b|^2)
$$
で、またも量子ビット$A$のブロッホ球が$S^8$の中の$S^2$として入れ子で現れます。ブロッホベクトルとの対応も2量子ビットとまったく同様で、$\sigma\otimes\mathrm I_4$の期待値$\boldsymbol r_A$と$e_2,\dots,e_7$成分の大きさ$C_{A|BC}$について
$$
x_A^2+y_A^2+z_A^2+C_{A|BC}^2=1
$$
が成り立ちます。

# 3量子ビットの具体例

3量子ビットの代表的な状態のホップ像を挙げます。像は$(2o_1o_2^*;\ |o_1|^2-|o_2|^2)$で、八元数部分と実成分の順に書きます。
$$
\begin{aligned}
|000\rangle\ &\mapsto\ (0;\ 1) && C_{A|BC}=0\\
|0\rangle\otimes\frac{|00\rangle+|11\rangle}{\sqrt2}\ &\mapsto\ (0;\ 1) && C_{A|BC}=0\\
\mathrm{GHZ}=\frac{|000\rangle+|111\rangle}{\sqrt2}\ &\mapsto\ (-e_6;\ 0) && C_{A|BC}=1\\
\mathrm{W}=\frac{|001\rangle+|010\rangle+|100\rangle}{\sqrt3}\ &\mapsto\ \left(\frac{2(e_2+e_4)}3;\ \frac13\right) && C_{A|BC}=\frac{2\sqrt2}3\\
\frac{|000\rangle+|101\rangle}{\sqrt2}\ &\mapsto\ (-e_2;\ 0) && C_{A|BC}=1\\
\frac{|000\rangle+|110\rangle}{\sqrt2}\ &\mapsto\ (-e_4;\ 0) && C_{A|BC}=1
\end{aligned}
$$
2つ目は$BC$の内部が最大にもつれた状態（$B$と$C$のベル状態）ですが、A|BCの分割では分離可能なので、$|000\rangle$と同じ北極に写ります（同じファイバー上の2点です）。

下の2つは、それぞれ$A$と$C$、$A$と$B$がベル状態を組み、残り1つが$|0\rangle$の状態です。GHZ状態と見比べると、$A$の相手が$C$だけなら$e_2$軸、$B$だけなら$e_4$軸、$BC$全体（GHZ）なら$e_6$軸と、**誰ともつれているかに応じて像の軸が分かれます**。W状態は$e_2$と$e_4$の中間方向を向き、もつれの深さも$C_{A|BC}=\frac{2\sqrt2}3<1$と最大値に届きません。

&&&rem 他の分割と3量子ビットのもつれの分類
ここで測っているのはあくまでA|BC分割のもつれです。B|AC・C|ABの分割は、テンソル積の並べ替えで同様に測れます。3つの分割すべてで$C=0$となることが、完全な分離可能性$\Psi=\Phi_A\otimes\Phi_B\otimes\Phi_C$と同値です。GHZ状態とW状態は3量子ビットのもつれの代表的な2つの型ですが、その分類（SLOCC分類）には立ち入りません。
&&&

# 複素構造の選択と$\operatorname{SU}(3)$

八元数への詰め込みで唯一恣意的だったのは、複素数の$i$の役割を担う虚数単位として$e_1$を選んだことです。この選択で虚部7次元は、$e_1$自身と、残り6次元$\{e_2,e_3,e_4,e_5,e_6,e_7\}$すなわち$\mathbb{C}$上の基底$\{e_2,e_4,e_6\}$が張る$\mathbb{C}^3$に分かれます。「虚数単位を1つ固定すると残りの6次元が$\mathbb{C}^3$になる」——これはリー群の記事[[7shi-lie5]]の結びで、$\operatorname{SU}(3)$が八元数の対称性の中に再登場する構図として予告したものです。ここでその再登場を確認します。

八元数の**自己同型**、すなわち乗法を保つ実線形写像$\varphi(xy)=\varphi(x)\varphi(y)$のうち、$e_1$を固定するものを考えます。$\varphi(e_1x)=e_1\varphi(x)$となるため、$\varphi$は複素構造と両立し、$\mathbb{C}^3$に$\mathbb{C}$-線形に作用します。この作用がどんな行列になるかを、対角的な場合で計算してみます。

&&&ex 対角位相と行列式
基底に位相を掛ける写像$\varphi\colon e_2\mapsto e^{i\theta_1}e_2,\ e_4\mapsto e^{i\theta_2}e_4,\ e_6\mapsto e^{i\theta_3}e_6$が自己同型になる条件を求めます。$e_2e_4=e_6$（三つ組$246$）を保つ必要があるので、左辺の像を計算すると
$$
\begin{aligned}
(e^{i\theta_1}e_2)(e^{i\theta_2}e_4)
&=(\cos\theta_1\,e_2+\sin\theta_1\,e_3)(\cos\theta_2\,e_4+\sin\theta_2\,e_5)\\
&=\cos\theta_1\cos\theta_2\,e_6+\cos\theta_1\sin\theta_2\,e_7+\sin\theta_1\cos\theta_2\,e_7-\sin\theta_1\sin\theta_2\,e_6\\
&=\cos(\theta_1+\theta_2)\,e_6+\sin(\theta_1+\theta_2)\,e_7\\
&=e^{-i(\theta_1+\theta_2)}e_6
\end{aligned}
$$
となります（三つ組$257,347,365$より$e_2e_5=e_7,\ e_3e_4=e_7,\ e_3e_5=-e_6$。最後の等号は$e_1e_6=-e_7$による）。これが$\varphi(e_6)=e^{i\theta_3}e_6$と一致する条件は
$$
\theta_1+\theta_2+\theta_3\equiv0\pmod{2\pi}
$$
すなわち、対角行列の行列式$e^{i(\theta_1+\theta_2+\theta_3)}$が1であることです。
&&&

位相を勝手に掛けるだけでは乗法が壊れ、行列式1という制約が乗法の保存から現れました。これは一般に成り立つ事実の対角的な断面です。$e_1$を固定する自己同型の全体は、$\mathbb{C}^3$へのユニタリかつ行列式1の作用、すなわち$\operatorname{SU}(3)$とちょうど一致します。数体系の単位球面からは決して得られなかった$\operatorname{SU}(3)$[[7shi-lie5]]が、八元数の対称性として現れる場面です。なお八元数の自己同型全体は例外型リー群$G_2$と呼ばれる14次元の群で、$\operatorname{SU}(3)$はその中で$e_1$を固定する部分群です。$G_2$は今後のリー群シリーズで扱います。

この$\operatorname{SU}(3)$がもつれの描像で果たす役割を見ます。自己同型は1を固定して虚部を虚部に写すため共役と可換で、ホップ写像とも両立します。
$$
H_{\mathbb{O}}(\varphi(o_1),\varphi(o_2))=(\varphi(2o_1o_2^*),\ |o_1|^2-|o_2|^2)
$$
$\varphi$は像の複素部分（$A$のブロッホ球$S^2$）を固定し、$e_2,\dots,e_7$成分のなす$\mathbb{C}^3$を回します。つまり$\boldsymbol r_A$と$C_{A|BC}$を保ったまま、**もつれの「型」だけを回す**操作です。最大もつれ（$C_{A|BC}=1$）の像の全体は$\mathbb{C}^3$内の単位球面$S^5$で、$\operatorname{SU}(3)$はその上に推移的に作用します（[[7shi-lie5]]で$\operatorname{SU}(3)$を$S^5$上のファイバー束と見た構図がここに現れます）。2量子ビットで最大もつれの像が$(\mathbf j,\mathbf k)$平面の$S^1$だったことの一般化です。

具体例で見た3つの軸も、$\operatorname{SU}(3)$の1つの軌道上にあります。巡回置換
$$
(z_2,z_3,z_4)\mapsto(z_4,z_2,z_3)
$$
は置換行列として行列式1の$\operatorname{SU}(3)$の元であり、像の軸を$-e_6\mapsto-e_2\mapsto-e_4\mapsto-e_6$と回します。状態の側では、$\varphi$は両ブロック$\Psi_0,\Psi_1$に同じユニタリ$U_{BC}$として作用するため、3量子ビット全体には$\mathrm I_2\otimes U_{BC}$、すなわち$BC$側だけの変換として働きます（$A$に触れないので$\boldsymbol r_A$が不変なのは当然です）。実際にこの巡回置換を状態に施すと
$$
\mathrm{GHZ}\ \mapsto\ \frac{|000\rangle+|101\rangle}{\sqrt2}\ \mapsto\ \frac{|000\rangle+|110\rangle}{\sqrt2}
$$
と、GHZ状態と2つのベル状態型がちょうど移り合います。

# 非結合性の現れとファイバー$S^7$

最後に、前回確認した非結合性の帰結[[7shi-hopfext]]が量子側でどう見えるかを整理します。

2量子ビットでは、ファイバーに沿う移動＝$B$の状態の取り替えが、単位四元数の右乗算という群作用で実現されていました。八元数では右からの共通乗算$(o_1q,\ o_2q)$がホップ像を保ちません。つまり「$BC$の状態を右乗算で取り替える」という操作は、非結合性のために写像と両立しなくなります。

代わりに前回構成した、恒等式に基づくファイバーのパラメータ表示が使えます。分離可能な像点$(c,\ r)$（$c$は複素数、$r\neq-1$）に対して
$$
|o_1|^2=\frac{1+r}2,\quad o_2^*=\frac{o_1^*c}{1+r}
$$
と定めると$H_{\mathbb{O}}(o_1,o_2)=(c,\ r)$でした。共役を取ると$o_2=\dfrac{c^*}{1+r}\,o_1$、すなわち$o_2$は$o_1$の複素数倍なので、**このファイバー$S^7$上の点はすべてA|BC分離可能**です。$o_1$を球面上で自由に動かすことは、$BC$の合成状態$X\in S^7$を自由に選ぶことにあたります。

こうして、ファイバーの系列に量子的な意味が通りました。
$$
\underbrace{S^1}_{\text{全体位相}}\ \to\ \underbrace{S^3}_{B\text{の状態}}\ \to\ \underbrace{S^7}_{BC\text{の状態}}
$$
いずれも「$A$のブロッホベクトルを固定したとき残りに許される自由度」がファイバーです。ただし$S^1,S^3$が群（全体位相の$U(1)$、単位四元数の$\operatorname{Sp}(1)$）としてファイバーを兼ねたのに対し、$S^7$は群になれず[[7shi-lie5]]、ファイバーを動かす対称性は$\operatorname{Spin}(7)$という外部の群に委ねられます。前回の左乗算による構成は、その具体的な代替手段でした。

&&&rem ファイバー内のもつれ
$BC$の内部のもつれはファイバー内の自由度であり、像には現れません。具体例で見たように、$B$と$C$がベル状態を組む$|0\rangle\otimes\frac{|00\rangle+|11\rangle}{\sqrt2}$は$|000\rangle$と同じ北極に写ります。ホップ像が検出するのはあくまでA|BC分割のもつれです。
&&&

# まとめ

複素数・四元数・八元数のホップ写像の定義域$S^3,S^7,S^{15}$は、それぞれ1・2・3量子ビットの状態空間です。状態を四元数・八元数のペアに詰めると、ホップ像の成分がもつれの情報を与えます。

- 分離可能$\iff\alpha\delta-\beta\gamma=0\iff$像の$\mathbf j,\mathbf k$成分の消失。3量子ビットではA|BC分離可能$\iff$像の$e_2,\dots,e_7$成分の消失
- 分離可能な状態の像は、量子ビット$A$のブロッホ球$S^2$の入れ子
- 像のノルム1の内訳が$x_A^2+y_A^2+z_A^2+C^2=1$を与え、もつれる（$C>0$）と$A$単独のブロッホベクトルは球の内部（混合状態）に落ちる
- ファイバーは「残りの系の状態空間」であり、$S^1$（全体位相）$\to S^3$（$B$の状態）$\to S^7$（$BC$の状態）と一般化される
- 複素構造（$e_1$）の固定は自己同型群$\operatorname{SU}(3)$を呼び込み、$\boldsymbol r_A$と$C_{A|BC}$を保ったまま、もつれの型（$e_2,e_4,e_6$の軸）を回す
- 非結合性は「$BC$の状態の右乗算による取り替え」の群構造喪失として現れ、ファイバー$S^7$は外部の対称性$\operatorname{Spin}(7)$と左乗算構成に委ねられる

&&&rem 参考文献
2量子ビットと四元数ホップファイブレーションの対応はMosseriとDandoloff[[mosseri]]、3量子ビットと八元数への拡張はBernevigとChen[[bernevig]]によります。本記事の導出は本シリーズの範囲で自己完結するよう構成し直したものです。
&&&
