複素解析をディラック作用素の核として組み直した枠組み[[7shi-cla1]][[7shi-cla2]][[7shi-cla3]]を、4次元に持ち上げます。[[7shi-cla1]]で位置ベクトルに$e_1$を掛けて複素数$z=e_1\boldsymbol x$を作った操作を$\operatorname{Cl}_{4,0}(\mathbb R)$で繰り返すと、四元数の変数$q=e_0\boldsymbol x$とフューター作用素が得られます。フューターの意味で正則な四元数値関数では、恒等写像も逆数も正則でなく、正則関数の積や合成も正則になりません。一方で、成分の調和性、コーシーの積分公式、平均値の性質は、核を取り替えるだけで生き残ります。恒等写像と逆元の正則性が失われる理由は次元の数え上げに、積の正則性が保たれない理由は非可換性に現れます。クリフォード代数を経由するのは、四元数の関係式を導き直すためではなく、作用素と積分公式の出所を2次元の場合と共通にするためです。

# 四元数変数

## ベクトル変数から四元数へ

[[7shi-cla2]]の一般次元の設定を4次元で使います。$\operatorname{Cl}_{4,0}(\mathbb R)$の生成元に$0$から添字を振って$e_0,e_1,e_2,e_3$（$e_a^2=1$、$a\ne b$なら$e_ae_b=-e_be_a$）とし、位置ベクトルとディラック作用素を

$$
\boldsymbol x=\sum_{a=0}^3x_ae_a,\qquad D=\sum_{a=0}^3e_a\partial_a
$$

とします。$D^2=\Delta$（4次元のラプラシアン）は[[7shi-cla1]]と同じ計算で成り立ちます。

[[7shi-cla1]]では、位置ベクトルに$e_1$を左から掛けて偶部分の元$z=e_1\boldsymbol x=x+Iy$を作り、作用素にも$e_1$を掛けてヴィルティンガー微分$e_1D=\partial_x+I\partial_y$を得ました。同じことを$e_0$で行います。

&&&def 四元数に対応する単位と変数
$$
i=e_0e_1,\qquad j=e_0e_2,\qquad k=e_0e_3
$$

とし、位置ベクトルに$e_0$を左から掛けたものを

$$
q=e_0\boldsymbol x=x_0+x_1i+x_2j+x_3k
$$

と書きます。
&&&

$i,j,k$はどれも2ベクトルで、偶部分代数$\operatorname{Cl}_{4,0}^0(\mathbb R)$に属します。$i^2=e_0e_1e_0e_1=-e_0^2e_1^2=-1$であり、$j,k$も同様です。また$ij=e_0e_1e_0e_2=-e_1e_2$、$ji=-e_2e_1=e_1e_2$より$ij=-ji$で、他の組も同様に反交換します。これは[[7shi-clif1]]で偶部分代数の同型$\operatorname{Cl}_{p,q}^0(\mathbb R)\cong\operatorname{Cl}_{q,p-1}(\mathbb R)$を示したときの、2乗が$+1$の生成元$e_0$との積で新しい生成元を作る構成と（掛ける順序による符号を除いて）同じもので、$i,j,k$は$\operatorname{Cl}_{0,3}(\mathbb R)$の生成元として振る舞います。

2次元では、$e_1$を掛けた$e_1e_2=I$が1つだけ現れ、偶部分は$\mathbb C$でした。4次元では、$e_0$を掛けた単位が3つ現れます。

## 四元数の成分への制限

ハミルトンの関係式$ijk=-1$は、まだ成り立っていません。擬スカラー$I_4=e_0e_1e_2e_3$を使うと

$$
ijk=e_0e_1e_0e_2e_0e_3=-e_1e_2e_0e_3=-I_4
$$

です。$I_4^2=1$であり、4次元の擬スカラーは偶数グレードの元と可換です。そこで[[7shi-clif1]]と同じく冪等元$P_\pm=(1\pm I_4)/2$を取ると、偶部分代数は2つの成分に分かれます。

&&&fml 偶部分代数の分解
$$
\operatorname{Cl}_{4,0}^0(\mathbb R)=\operatorname{Cl}_{4,0}^0(\mathbb R)P_+\oplus\operatorname{Cl}_{4,0}^0(\mathbb R)P_-\cong\mathbb H\oplus\mathbb H
$$
&&&

$P_+$の成分では$I_4=1$として働くので$ijk=-1$となり、$i,j,k$はハミルトンの四元数の単位そのものです。$P_-$の成分では$ijk=1$ですが、$k$を$-k$に取り替えれば同じく四元数になります。

作用素$e_0D$の係数$1,i,j,k$は偶数グレードなので$I_4$と可換であり、$e_0D$は2つの成分を混ぜません。したがって偶部分に値を取る関数は、$P_+$の成分と$P_-$の成分に分けて扱えます。$P_-$の成分を$i,j,-k$によって四元数と同一視すると、作用素の$\partial_3$の係数だけが反転します。したがって座標$x_3$も反転すれば、$P_+$の成分と同じ方程式になります。以下では$P_+$の成分だけを考え、関数の値を四元数$\mathbb H$に取ります。

厳密には、四元数の$1,i,j,k$と同一視するのは射影した$P_+,\ iP_+,\ jP_+,\ kP_+$であり、この成分の単位元は$P_+$です。四元数の変数$q$も$(e_0\boldsymbol x)P_+$の像です。以下の$\bar q$、$q^{-1}$、核$E$などの式は、射影前の偶部分で計算してから$P_+$の成分に射影したものと読みます。$q\bar q$がスカラーになるため、射影前の逆元を射影すれば成分の中での逆元になります。

関数を$\mathbb H$値に限ってよいのは、このように作用素が偶部分代数の直和分解を保つからです。

## フューター作用素

作用素に$e_0$を左から掛けたものを、四元数の記法で書きます。

&&&def フューター作用素
$$
\mathcal D=e_0D=\partial_0+i\partial_1+j\partial_2+k\partial_3,\qquad
\bar{\mathcal D}=\partial_0-i\partial_1-j\partial_2-k\partial_3
$$

$\mathbb H$値の関数$F$に対して、右からの作用を$G\mathcal D=\partial_0G+(\partial_1G)i+(\partial_2G)j+(\partial_3G)k$とします。$\mathcal DF=0$を満たす$F$を**左正則**、$G\mathcal D=0$を満たす$G$を**右正則**と呼びます。単に正則と言えば左正則を指します。
&&&

$e_0$は$e_0^2=1$を満たして可逆なので、$\mathcal DF=0$は$DF=0$と同値です。正則性は、$\operatorname{Cl}_{4,0}(\mathbb R)$のディラック作用素についてのモノジェニック性を、偶部分の関数に対して書き直したものにあたります。[[7shi-cla1]]で、偶部分の関数について左モノジェニック性がコーシー＝リーマンの方程式と同値になったのと同じ関係です。

右からの作用も$D$から得られます。$e_0e_0=1$、$e_0e_1=i$、$e_0e_2=j$、$e_0e_3=k$なので、$G$に$e_0$を右から掛けて$D$を右から作用させると

$$
(Ge_0)D=\sum_a(\partial_aG)e_0e_a=G\mathcal D
$$

です。

共役を取った作用素$\bar{\mathcal D}$は、$e_0$を$D$の右に置いたものです。$\ell=1,2,3$では$e_\ell e_0=-e_0e_\ell$なので、$D(e_0F)=\partial_0F+\sum_{\ell=1}^3e_\ell e_0\,\partial_\ell F=\bar{\mathcal D}F$となります。ここから$D^2=\Delta$がそのまま使えます。

&&&fml 共役との積
$$
\bar{\mathcal D}\mathcal D=\mathcal D\bar{\mathcal D}=\Delta
$$
&&&

&&&prf
$\bar{\mathcal D}F=D(e_0F)$であるから、$\bar{\mathcal D}\mathcal DF=D(e_0e_0DF)=D^2F=\Delta F$である。$\mathcal D\bar{\mathcal D}F=e_0D\,D(e_0F)=e_0\Delta(e_0F)=\Delta F$も同様である。
&&&

[[7shi-cla1]]の$D^2=\Delta$とは違って、$\mathcal D$そのものの2乗ではなく共役との積になります。これは、$e_0$を作用素の片側にだけ掛けたためです。四元数の計算では$i,j,k$の2乗が$-1$であることと反交換性から直接確かめられますが、ここでは$D^2=\Delta$の言い換えとして得られます。

&&&thm 成分の調和性
$C^2$級の正則関数の各成分は、4次元の調和関数です。
&&&

$\mathcal DF=0$なら$\Delta F=\bar{\mathcal D}\mathcal DF=0$です。右正則な関数でも、右からの作用について同じ計算をすれば同じことが言えます。

# 恒等写像と逆元

## 次元の数え上げ

[[7shi-cla1]]では$(\partial_x+I\partial_y)z=1+I^2=0$であり、恒等写像$z$は正則でした。四元数の変数では

$$
\mathcal Dq=1+i^2+j^2+k^2=1-3=-2
$$

となり、恒等写像は正則ではありません。この$-2$の出所は、ベクトル変数に戻すとはっきりします。

&&&fml 生成元の挟み込み
$\operatorname{Cl}_{n,0}(\mathbb R)$の生成元$e_1,\dots,e_n$について

$$
\sum_{a=1}^ne_ae_1e_a=(2-n)e_1
$$
&&&

&&&prf
$a=1$の項は$e_1^3=e_1$である。$a\ne1$の項は$e_ae_1e_a=-e_1e_a^2=-e_1$であり、これが$n-1$個ある。和は$e_1-(n-1)e_1=(2-n)e_1$である。
&&&

$D(e_1\boldsymbol x)=\sum_ae_ae_1\,\partial_a\boldsymbol x=\sum_ae_ae_1e_a$なので、$e_1\boldsymbol x$が$D$で消えるのは$n=2$のときに限ります。$n=4$で添字を$0$から振れば$D(e_0\boldsymbol x)=-2e_0$であり、$e_0$を左から掛けると$\mathcal Dq=-2$です。2次元で$z$が正則になるのは、生成元の反交換で出る符号$-1$の個数と、自分自身との積で出る$+1$の個数がちょうど釣り合うからです。

逆元も同じ理由で正則になりません。$\boldsymbol x$の逆元は$\boldsymbol x^{-1}=\boldsymbol x/|\boldsymbol x|^2$であり、[[7shi-cla1]]の計算を$n$次元で繰り返すと$D\boldsymbol x=n$、$D|\boldsymbol x|^{-2}=-2\boldsymbol x/|\boldsymbol x|^4$より

$$
D\boldsymbol x^{-1}=-\frac{2\boldsymbol x^2}{|\boldsymbol x|^4}+\frac n{|\boldsymbol x|^2}=\frac{n-2}{|\boldsymbol x|^2}
$$

です。四元数の共役を$\bar q=x_0-x_1i-x_2j-x_3k$とすると、$e_ke_0=-e_0e_k$より$\bar q=\boldsymbol xe_0$であり、$q\bar q=e_0\boldsymbol x^2e_0=|\boldsymbol x|^2$です。したがって$q^{-1}=\bar q/|q|^2=\boldsymbol x^{-1}e_0$であり

$$
\mathcal Dq^{-1}=e_0D(\boldsymbol x^{-1})e_0=\frac2{|q|^2}\ne0
$$

となります。[[7shi-cla3]]で積分公式の核だった$\boldsymbol x^{-1}$は、4次元ではモノジェニックではありません。恒等写像と逆元がともに正則でなくなるのは、非可換性ではなく、次元が$2$でないことによります。

&&&ex 冪
$q^2=x_0^2-x_1^2-x_2^2-x_3^2+2x_0(x_1i+x_2j+x_3k)$に$\mathcal D$を作用させると

$$
\mathcal Dq^2=-4x_0
$$

です。複素解析では$z^n$がすべて正則で、冪級数が正則関数を与えましたが、四元数では冪$q^n$が正則でないため、この構成は使えません。
&&&

## 核の取り直し

逆元の計算で$|\boldsymbol x|^{-2}$を$|\boldsymbol x|^{-n}$に替えると、$D|\boldsymbol x|^{-n}=-n\boldsymbol x/|\boldsymbol x|^{n+2}$より$D(\boldsymbol x/|\boldsymbol x|^n)=-n/|\boldsymbol x|^n+n/|\boldsymbol x|^n=0$となり、打ち消し合います。右からの作用も同様です。

&&&prop コーシー＝フューター核
$$
E(q)=\frac{\bar q}{|q|^4}=\frac{\boldsymbol x}{|\boldsymbol x|^4}\,e_0
$$

は原点を除いて左右ともに正則です。
&&&

&&&prf
$\bar q=\boldsymbol xe_0$と$|q|=|\boldsymbol x|$から等式を得る。$\boldsymbol x/|\boldsymbol x|^4$は左右ともに$D$で消えるから、$\mathcal DE=e_0D(\boldsymbol x/|\boldsymbol x|^4)e_0=0$、$E\mathcal D=(Ee_0)D=(\boldsymbol x/|\boldsymbol x|^4)D=0$である。
&&&

2次元では$\boldsymbol x/|\boldsymbol x|^2=\boldsymbol x^{-1}$だったため、核と逆元が一致していました。4次元では両者が分かれ、核は$-3$次の同次関数になります。[[7shi-cla3]]では、核と法線の積を半径$\varepsilon$の円の上で積分した値が半径によらないことを使いました。$S^3$の上では面積要素が$\varepsilon^3$に比例するので、同次な核で、球面上の積分が半径によらない$0$でない値になることを要求すると、次数は$-3$でなければなりません。$-1$次の$q^{-1}$はこの条件を満たさず、$-3$次の$E$が満たします。

# 積と合成

## 積の微分

$h_0=1,\ h_1=i,\ h_2=j,\ h_3=k$と書くと$\mathcal D=\sum_ah_a\partial_a$です。$\mathbb H$値の関数$f,g$の積に作用させると、次の式が成り立ちます。

&&&fml 積の微分
$$
\mathcal D(fg)=(\mathcal Df)g+f(\mathcal Dg)+\sum_{a=1}^3[h_a,f]\,\partial_ag
$$

ここで$[h_a,f]=h_af-fh_a$は交換子です。
&&&

&&&prf
$\partial_a(fg)=(\partial_af)g+f\partial_ag$より$\mathcal D(fg)=(\mathcal Df)g+\sum_ah_af\,\partial_ag$である。$h_af=fh_a+[h_a,f]$を代入すれば$\sum_afh_a\partial_ag=f(\mathcal Dg)$となり、残りが交換子の項である。$h_0=1$の交換子は$0$である。
&&&

複素解析では値が可換なので交換子の項が消え、正則関数の積は正則です。四元数では交換子の項が残り、$f,g$がともに正則でも$fg$は一般に正則ではありません。

## フューター変数

恒等写像の代わりに正則な1次式を探すと、次のものがあります。

&&&def フューター変数
$$
\zeta_1=x_1-x_0i,\qquad\zeta_2=x_2-x_0j,\qquad\zeta_3=x_3-x_0k
$$
&&&

$\mathcal D\zeta_1=\partial_0\zeta_1+i\partial_1\zeta_1=-i+i=0$であり、右からの作用でも$-i+i=0$です。$\zeta_2,\zeta_3$も同様で、フューター変数は左右ともに正則です。$\zeta_1=-i(x_0+x_1i)$と書けるので、$\zeta_1$は$(x_0,x_1)$平面の複素変数を$-i$倍したものです。$\zeta_1^2$は、積の微分の交換子の項が$[i,\zeta_1]\partial_1\zeta_1=0$（$\zeta_1$は$i$と可換）だけなので正則です。

&&&ex 異なるフューター変数の積
$\zeta_1\zeta_2$では、積の微分の交換子の項のうち$\partial_2\zeta_2=1$に掛かる$[j,\zeta_1]$だけが残ります。

$$
\mathcal D(\zeta_1\zeta_2)=[j,\zeta_1]=-x_0(ji-ij)=2x_0k
$$

となり、$\zeta_1\zeta_2$は正則ではありません。$\zeta_2\zeta_1$では同様に$[i,\zeta_2]=-2x_0k$が残るので、対称化した積$\zeta_1\zeta_2+\zeta_2\zeta_1$は正則です。
&&&

複素数値の正則多項式は、$z$の冪の1次結合で表されます。四元数では、フューター変数の対称化した積（各変数を指定の個数ずつ、すべての順序で掛けた積の和）が、右から四元数の係数を掛けた1次結合として、斉次な正則多項式の全体を張ります（フューター多項式）。本記事では一般の場合の証明は扱いません。

## 左正則と右正則

[[7shi-cla1]]では、偶部分の関数について左モノジェニックが正則、右モノジェニック（$FD=0$）が反正則にあたりました。本記事の右正則性は$F\mathcal D=(Fe_0)D=0$であり、元の$D$による右モノジェニック性そのものではありません。2次元でこれにあたる条件$F(\partial_x+I\partial_y)=0$は、偶部分が可換なので左からの条件$(\partial_x+I\partial_y)F=0$と一致します。四元数では、非可換性のために左正則と右正則が別の条件に分かれます。

&&&ex 左正則だが右正則でない関数
定数を右から掛けても左正則性は保たれるので、$\zeta_1j$は左正則です。一方、右からの作用は

$$
(\zeta_1j)\mathcal D=\partial_0(\zeta_1j)+\partial_1(\zeta_1j)\,i=-ij+ji=-2k
$$

となり、右正則ではありません。
&&&

左正則な関数は右から定数を掛けても左正則ですが、左から定数を掛けると一般に左正則ではなくなります。正則関数の全体は、右からの定数倍についてのみ閉じています。

## 合成

複素解析では、正則関数の合成は正則です。四元数では、合成も正則性を保ちません。

&&&ex フューター変数の合成
$\zeta_1(q)=x_1-x_0i$の変数$q$に$\zeta_2=x_2-x_0j$を代入します。$\zeta_2$の実部は$x_2$、$i$成分は$0$なので

$$
\zeta_1(\zeta_2)=0-x_2i=-x_2i,\qquad
\mathcal D(-x_2i)=-ji=k\ne0
$$

であり、正則関数どうしの合成が正則になりません。
&&&

複素解析で合成が正則性を保つのは、正則性が複素変数$z$についての微分可能性として言い換えられ、連鎖律が使えるためです。四元数値の関数でも実変数としての連鎖律は使えますが、そこから正則性が合成で保たれるとは言えません。恒等写像$q$自体が正則でないことも、この正則性が変数$q$についての微分可能性ではないことを示しています。

# フューターの定理

冪$q^n$が正則でないため、正則関数を作るには別の手段が要ります。複素解析の正則関数から四元数の正則関数を作る方法として、フューターの定理があります。

## 軸対称な関数

四元数を$q=x_0+\boldsymbol q$（$\boldsymbol q=x_1i+x_2j+x_3k$）と実部と虚部に分け、$r=|\boldsymbol q|$、$\omega=\boldsymbol q/r$とします。$\omega^2=-1$なので、$\omega$を固定すると$x_0+\omega y$（$y\in\mathbb R$）は複素平面の写しになり、$q=x_0+\omega r$はその中の点です。

複素共役で不変な開集合$U\subset\mathbb C$の上の正則関数$f(z)=u(x,y)+iv(x,y)$で、$f(\bar z)=\overline{f(z)}$を満たすもの（$u$は$y$について偶、$v$は奇）を取り、各平面の写しに持ち込みます。定義域は$\Omega_U=\{x_0+\boldsymbol q:x_0+i|\boldsymbol q|\in U\}$です。たとえば$f(z)=1/z$なら$U=\mathbb C\setminus\{0\}$、$\Omega_U$は原点を除いた$\mathbb R^4$です。

$$
\tilde f(q)=u(x_0,r)+\omega\,v(x_0,r)
$$

$\omega$は$r>0$でしか定義されないので、実軸の上（$r=0$）では$\tilde f$を連続的に延長して定めます（下の定理の証明で見るとおり、延長は滑らかです）。

$f(z)=z^n$なら、$x_0+\omega r$の冪を$x_0+iy$の冪と同じく計算できるので$\tilde f=q^n$、$f(z)=1/z$なら$\tilde f=q^{-1}$です。このように$x_0$と$r$の関数$A,B$で$A+\omega B$と書ける関数に、$\mathcal D$を作用させます。

&&&fml 軸対称な関数の微分 [fml-axial]
$A,B$を$x_0,r$の実数値$C^1$級関数とすると、$r>0$で

$$
\mathcal D(A+\omega B)=(A+\omega B)\mathcal D
=\partial_0A-\partial_rB-\frac2rB+\omega\,(\partial_rA+\partial_0B)
$$
&&&

&&&prf
$h_1,h_2,h_3=i,j,k$とする。$\partial_\ell r=x_\ell/r$より$\sum_\ell h_\ell\partial_\ell A=\omega\,\partial_rA$である。$\varphi=B/r$とおくと$\omega B=\boldsymbol q\varphi$であり

$$
\sum_{\ell=1}^3h_\ell\partial_\ell(\boldsymbol q\varphi)
=\sum_\ell h_\ell^2\,\varphi+\sum_\ell h_\ell\boldsymbol q\,\frac{x_\ell}r\,\partial_r\varphi
=-3\varphi+\frac{\boldsymbol q^2}r\,\partial_r\varphi
=-3\varphi-r\,\partial_r\varphi
$$

となる。$r\,\partial_r\varphi=\partial_rB-B/r$を代入すると$-\partial_rB-2B/r$である。$\partial_0$の項$\partial_0A+\omega\,\partial_0B$と合わせて左からの作用を得る。右からの作用では$h_\ell\boldsymbol q$が$\boldsymbol qh_\ell$に替わるが、$\sum_\ell\boldsymbol qh_\ell x_\ell=\boldsymbol q^2$なので同じ結果になる。$A,B$は実数値なので、$\omega$や$h_\ell$と可換であることを使った。
&&&

$\tilde f$に当てはめます。コーシー＝リーマンの方程式$u_x=v_y$、$u_y=-v_x$より$\omega$の係数と$\partial_0A-\partial_rB$が消え

$$
\mathcal D\tilde f=-\frac{2v}r
$$

が残ります。$f(z)=z$では$v=r$なので$\mathcal Dq=-2$であり、前節の結果に一致します。$-2B/r$の係数$2$は、$\omega$に直交する虚数単位の個数$3-1$です。複素数では虚数単位が1つだけなのでこの項は現れず、$\tilde f=f$のまま正則でした。

## 定理

$\mathcal D\tilde f$は$0$になりませんが、残ったものは調和関数です。

&&&thm フューターの定理
$f$が上の条件を満たす正則関数なら、$\Delta\tilde f$は左右ともに正則です。
&&&

&&&prf
$\mathcal D$は定数係数なので$\Delta$と可換であり、$\mathcal D\Delta\tilde f=\Delta\mathcal D\tilde f=-2\Delta(v/r)$である。$x_0$と$r$だけの関数$g$に対して、4次元のラプラシアンは$\Delta g=\partial_0^2g+\partial_r^2g+(2/r)\partial_rg$と書ける（3次元の動径方向のラプラシアンに$\partial_0^2$を加えたもの）。$g=v/r$を代入すると

$$
\partial_r^2\frac vr+\frac2r\,\partial_r\frac vr=\frac{\partial_r^2v}r
$$

より$\Delta(v/r)=(\partial_0^2v+\partial_r^2v)/r$となる。$v$は2変数の調和関数なので、これは$r>0$で$0$である。$f$は正則なので実解析的であり、$v$は第2変数について奇、$u$は偶である。したがって$v(x_0,r)/r$と$u(x_0,r)$は$r$の偶数冪で局所的に展開され、$r^2=x_1^2+x_2^2+x_3^2$の滑らかな関数として実軸まで延長される。$\tilde f=u+\boldsymbol q\,(v/r)$も同様に滑らかに延長され、$r>0$で得た等式は連続性により実軸の上でも成り立つ。右からの作用も[軸対称な関数の微分](#fml-axial)と同じ式になるので同様である。
&&&

$\mathcal D\tilde f=-2v/r$は、2次元から4次元に持ち込んだことで生じたずれです。このずれ自体が4次元の調和関数なので、ラプラシアンを1回作用させると消えます。

&&&ex 冪と逆数
$f(z)=z^n$から、$\Delta q^n$が正則関数になります。

$$
\Delta q=0,\qquad\Delta q^2=-4,\qquad\Delta q^3=-4(2q+\bar q)
$$

$\mathcal D(2q+\bar q)=2\cdot(-2)+4=0$であり、確かに正則です（$\mathcal D\bar q=1-i^2-j^2-k^2=4$）。

$f(z)=1/z$では$\tilde f=q^{-1}$であり

$$
\Delta q^{-1}=-4\,\frac{\bar q}{|q|^4}=-4E
$$

となります。$|q|^{-2}$が4次元の調和関数であることから、$\Delta(x_0/|q|^2)=2\partial_0|q|^{-2}=-4x_0/|q|^4$となり、虚部も同様です。
&&&

$1/z$から出発すると、コーシー＝フューター核が得られます。2次元では$1/z$そのものが核でしたが、4次元では逆数$q^{-1}$は$-1$次であり、核に必要な$-3$次に足りません。ラプラシアンは次数を$2$下げるので、1回作用させると$-3$次の核になります。

他の次元では、ずれの係数や、それを消すための作用素が変わります。本記事では扱いません。

# 積分公式

## 四元数による基本定理

[[7shi-cla2]]の$n$次元の両側形式を$\operatorname{Cl}_{4,0}(\mathbb R)$で使います。$M\subset\mathbb R^4$は区分的に滑らかな境界を持つ有界な領域とし、関数$F,G$は$M$とその境界を含む開集合上で定義された$\mathbb H$値の$C^1$級関数とします。核を使う場合は、その特異点を除いた開集合上でこの条件を課します。$M$の外向きの単位法線を$\boldsymbol n=\sum_an_ae_a$とし、法線にも$e_0$を左から掛けて四元数で書きます。

$$
n=e_0\boldsymbol n=n_0+n_1i+n_2j+n_3k
$$

両側形式$\int\bigl((G'D)F+G'(DF)\bigr)dV=\oint G'\boldsymbol nF\,dS$で$G'=Ge_0$とすると、$(Ge_0)D=G\mathcal D$、$Ge_0(DF)=G(\mathcal DF)$、$Ge_0\boldsymbol n=Gn$より次の形になります。

&&&fml 四元数による両側形式
$$
\int_M\bigl((G\mathcal D)F+G(\mathcal DF)\bigr)dV=\oint_{\partial M}G\,nF\,dS
$$
&&&

[[7shi-cla2]]の直方体による証明は、係数が定数であることしか使っていないので、係数を$1,i,j,k$に替えても同じ証明が通ります。一般の領域については、[[7shi-cla2]]と同じく概略です。

&&&thm コーシー＝フューターの積分定理
$M$上で$G\mathcal D=0$かつ$\mathcal DF=0$なら

$$
\oint_{\partial M}G\,nF\,dS=0
$$
&&&

## 小球面上の核

点$a$を中心とする半径$\varepsilon$の球面の上では、$q-a=\varepsilon n$（$n$は外向きの単位法線、$|n|=1$）です。$\bar nn=n\bar n=1$なので

&&&fml 球面の上の核と法線
$$
E(q-a)\,n=n\,E(q-a)=\frac{\varepsilon\bar n\,n}{\varepsilon^4}=\frac1{\varepsilon^3}
$$
&&&

となり、[[7shi-cla3]]と同じく核と法線の積はスカラーになります。単位球面$S^3$の面積（3次元の体積）は$2\pi^2$であり、面積要素は$dS=\varepsilon^3d\Omega$（$d\Omega$は$S^3$の面積要素）なので

$$
\oint_{|q-a|=\varepsilon}E(q-a)\,n\,dS=\frac1{\varepsilon^3}\cdot2\pi^2\varepsilon^3=2\pi^2
$$

であり、半径によりません。核が距離の3乗に反比例して減衰する分と、球面の面積が半径の3乗に比例して伸びる分が打ち消し合い、$S^3$の面積$2\pi^2$だけが残ります。

## コーシー＝ポンペイウの公式

[[7shi-cla3]]と同じ手順で、内部の値を境界の値と$\mathcal DF$の積分で表します。

&&&thm コーシー＝ポンペイウの公式 [thm-pompeiu]
$a$が$M$の内部の点なら

$$
F(a)=\frac1{2\pi^2}\left(\oint_{\partial M}E(q-a)\,nF\,dS-\int_ME(q-a)\,\mathcal DF\,dV\right)
$$
&&&

&&&prf
$G=E(q-a)$とし、$M$から中心$a$、半径$\varepsilon$の閉球を除いた領域$M_\varepsilon$に両側形式を当てはめる。$M_\varepsilon$上では$G\mathcal D=0$なので

$$
\int_{M_\varepsilon}G(\mathcal DF)\,dV=\oint_{\partial M}GnF\,dS-\oint_{|q-a|=\varepsilon}GnF\,dS
$$

である。小球面の項では、$M_\varepsilon$から見た外向きの法線が$-n$であることを符号に取り込んだ。小球面の上では$Gn=1/\varepsilon^3$、$dS=\varepsilon^3d\Omega$なので

$$
\oint_{|q-a|=\varepsilon}GnF\,dS=\oint_{S^3}F(a+\varepsilon n)\,d\Omega
\ \longrightarrow\ 2\pi^2F(a)\qquad(\varepsilon\to0)
$$

となる。左辺では$|G|=1/\rho^3$（$\rho=|q-a|$）、$dV=\rho^3d\rho\,d\Omega$より、被積分関数$G(\mathcal DF)\rho^3$は$a$の近くで有界であり、$\varepsilon\to0$で$M$上の積分に収束する。以上を整理すれば公式を得る。

$M$が球であっても、$M_\varepsilon$は球殻であり、[[7shi-cla3]]の円環のような1次元の基本定理だけによる証明は与えていない。両側形式の適用は[[7shi-cla2]]と同じく概略である。
&&&

&&&cor コーシー＝フューターの積分公式
$M$上で$\mathcal DF=0$なら、$M$の内部の点$a$で

$$
F(a)=\frac1{2\pi^2}\oint_{\partial M}E(q-a)\,nF\,dS
$$
&&&

[[7shi-cla3]]のコーシーの積分公式と比べると、式の形は同じで、核が$\boldsymbol x^{-1}$から$E=\bar q/|q|^4$に、$2\pi$が$2\pi^2$に替わっただけです。どちらも、核と法線の積を小さな球面の上で積分した値であり、単位球面の面積にあたります。

&&&rem 右正則な関数と基本解
両側形式で$G$の位置に$F$を、$F$の位置に$E(q-a)$を置くと、右からの作用に対する公式

$$
F(a)=\frac1{2\pi^2}\left(\oint_{\partial M}FnE(q-a)\,dS-\int_M(F\mathcal D)E(q-a)\,dV\right)
$$

が得られます。核が左右ともに正則であるため、どちらの側にも同じ核が使えます。また[[7shi-cla3]]と同じく、台がコンパクトな関数では境界の項が消え、$E/2\pi^2$が$\mathcal D$の基本解であることを表します。
&&&

## 平均値の性質

$M$を$a$を中心とする半径$R$の球とすると、境界の上で$E(q-a)n=1/R^3$が定数になります。

&&&cor 平均値の性質
$a$を中心とする半径$R$の閉球を含む開集合上で$\mathcal DF=0$なら

$$
F(a)=\frac1{2\pi^2R^3}\oint_{|q-a|=R}F\,dS=\frac1{2\pi^2}\oint_{S^3}F(a+Rn)\,d\Omega
$$
&&&

中心での値は、球面上の値の平均です。成分の調和性から、調和関数の平均値の性質としても得られます。

# まとめ

[[7shi-cla1]]の$z=e_1\boldsymbol x$、$e_1D=\partial_x+I\partial_y$を$\operatorname{Cl}_{4,0}(\mathbb R)$で繰り返すと、四元数の変数$q=e_0\boldsymbol x$とフューター作用素$\mathcal D=e_0D$が得られます。$i=e_0e_1$、$j=e_0e_2$、$k=e_0e_3$は偶部分代数$\operatorname{Cl}_{4,0}^0(\mathbb R)\cong\mathbb H\oplus\mathbb H$に属し、作用素がこの分解を保つので、関数の値を$\mathbb H$に取れます。$\bar{\mathcal D}\mathcal D=\Delta$は$D^2=\Delta$の言い換えです。

複素解析で成り立っていた性質は、次のように分かれます。

| 性質 | 2次元 | 4次元 | 原因 |
|---|---|---|---|
| 恒等写像の正則性 | $\bar\partial z=0$ | $\mathcal Dq=-2$ | 次元（$\sum_ae_ae_0e_a=(2-n)e_0$） |
| 逆元の正則性 | $\boldsymbol x^{-1}$は核 | $\mathcal Dq^{-1}=2/\lvert q\rvert^2$ | 次元（$D\boldsymbol x^{-1}=(n-2)/\lvert\boldsymbol x\rvert^2$） |
| 積の正則性 | 成り立つ | 交換子の項が残る | 非可換性 |
| 成分の調和性 | $D^2=\Delta$ | $\bar{\mathcal D}\mathcal D=\Delta$ | 壊れない |
| 積分公式 | 核$\boldsymbol x^{-1}$、$2\pi$ | 核$\bar q/\lvert q\rvert^4$、$2\pi^2$ | 核と定数だけが替わる |
| 平均値の性質 | 円周上の平均 | 球面上の平均 | 壊れない |

冪と合成も正則性を保たず、冪級数の代わりにフューター変数$\zeta_l=x_l-x_0h_l$の対称化した積が正則多項式を与えます。複素解析の正則関数$f=u+iv$を軸対称に持ち込んだ$\tilde f$は$\mathcal D\tilde f=-2v/r$のずれを持ちますが、ずれが調和関数であるため$\Delta\tilde f$は正則になります（フューターの定理）。$1/z$に当てはめると、$-1$次の$q^{-1}$から$-3$次のコーシー＝フューター核$-4E$が得られます。
