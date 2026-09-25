[[7shi-em1]]では、静的な電場と磁場を$F=\boldsymbol E+Ic\boldsymbol B$にまとめ、ディラック作用素$D$について$DF=(\rho-\boldsymbol J/c)/\varepsilon_0$が成り立つことを見ました。本記事では時間変化する場に進みます。$x_0=ct$を座標に加えた作用素$\mathcal D=\partial_0+D$を使うと、マクスウェル方程式は$\mathcal DF=(\rho-\boldsymbol J/c)/\varepsilon_0$という1本の式になり、ベクトル解析の4本の式はそのスカラー・ベクトル・2ベクトル・擬スカラー部として回収されます。この作用素は[[7shi-cla4]]の四元数解析の作用素$\mathcal D=\partial_0+i\partial_1+j\partial_2+k\partial_3$と同じ形をしていますが、係数の2乗が$-1$から$+1$に替わっています。その結果、共役との積はラプラシアンではなく波動作用素$\partial_0^2-\Delta$になり、正則関数の成分が調和関数になったことに対応して、電磁場の成分が波動方程式を満たすことが得られます。

# 設定

[[7shi-em1]]と同じく、$\operatorname{Cl}_{3,0}(\mathbb R)$の生成元を$e_1,e_2,e_3$、擬スカラーを$I=e_1e_2e_3$とし、空間のディラック作用素を$D=\sum_{k=1}^3e_k\partial_k$とします。時間$t$は光速$c$を掛けた$x_0=ct$で表し、$\partial_0=\partial/\partial x_0=\frac1c\partial_t$と書きます。

&&&def 時空の作用素
$$
\mathcal D=\partial_0+D=\partial_0+\sum_{k=1}^3e_k\partial_k,\qquad
\bar{\mathcal D}=\partial_0-D
$$
&&&

$\mathcal D$は、スカラー$1$と空間ベクトル$e_k$を係数とするパラベクトル型の作用素です。[[7shi-cla4]]のフューター作用素も、$1$と3つの四元数単位を係数とする同じ形式を持ちます。

電場$\boldsymbol E$と磁場$\boldsymbol B$は時刻と位置の関数とし、[[7shi-em1]]と同じく

$$
F=\boldsymbol E+Ic\boldsymbol B
$$

とまとめます。電荷密度$\rho$と電流密度$\boldsymbol J$も時間に依存してかまいません。関数はすべて必要な回数だけ連続微分可能とします。

# 静的な式を越えて

[[7shi-em1]]の4本の式は、電荷と電流が時間変化しない場合の法則でした。時間変化する場合には、そのうち回転の2本が修正を受けます。どちらも19世紀の電磁気学の中心的な発見です。

## 電荷の保存

[[7shi-em1]]で述べたとおり、電荷は生成も消滅もしません。領域$V$の中の電荷$Q=\int_V\rho\,dV$が変化するのは、境界$\partial V$を通って電荷が出入りするときだけです。境界を外向きに単位時間に通過する電荷は$\oint_{\partial V}\boldsymbol J\cdot\boldsymbol n\,dS$なので

$$
\frac{d}{dt}\int_V\rho\,dV=-\oint_{\partial V}\boldsymbol J\cdot\boldsymbol n\,dS
$$

です。右辺に発散定理を使い、これが任意の$V$で成り立つことから、微分形の**連続の式**

$$
\partial_t\rho+\nabla\cdot\boldsymbol J=0
$$

を得ます。電荷密度が時間変化する場所では$\nabla\cdot\boldsymbol J\ne0$です。

## 変位電流

静的な式のベクトル部はアンペールの法則$\nabla\times\boldsymbol B=\mu_0\boldsymbol J$でした。回転の発散は$0$なので、この式からは$\nabla\cdot\boldsymbol J=0$が従い、連続の式と両立しません。

具体的な場面として、コンデンサー（向かい合った2枚の金属板）を充電する回路を考えます。導線には電流が流れ、電荷は一方の板にたまっていきますが、板の間の隙間には電荷の流れがありません。導線を囲む閉曲線に[[7shi-em1]]のアンペールの法則の積分形を当てはめると、その閉曲線を縁とする曲面として導線を貫く曲面を選ぶか、板の間の隙間を通る曲面を選ぶかで、貫く電流が$I$と$0$に分かれてしまいます。隙間で起きているのは、板にたまる電荷が増え、板の間の電場が強くなることです。

マクスウェルは、アンペールの法則の電流に**変位電流**$\varepsilon_0\partial_t\boldsymbol E$を加えてこの矛盾を解消しました。ガウスの法則$\nabla\cdot\boldsymbol E=\rho/\varepsilon_0$を時間で微分すれば$\nabla\cdot(\varepsilon_0\partial_t\boldsymbol E)=\partial_t\rho=-\nabla\cdot\boldsymbol J$なので、$\boldsymbol J+\varepsilon_0\partial_t\boldsymbol E$の発散はつねに$0$になります。コンデンサーの隙間では、導線の電流$I$の代わりに変位電流が同じ大きさで流れ、どちらの曲面を選んでも結果が一致します。修正された式

$$
\nabla\times\boldsymbol B=\mu_0\boldsymbol J+\mu_0\varepsilon_0\partial_t\boldsymbol E
$$

を**アンペール＝マクスウェルの法則**と呼びます。電流がなくても、電場が時間変化すれば磁場の渦が生じます。

## 電磁誘導

コイル（巻いた導線）に磁石を近づけたり遠ざけたりすると、電池をつながなくてもコイルに電流が流れます。ファラデーが発見したこの現象を**電磁誘導**と呼びます。コイルを縁とする曲面を貫く磁束$\Phi=\int_S\boldsymbol B\cdot\boldsymbol n\,dA$が時間変化すると、コイルに沿って電荷を一周させる力の仕事（起電力）$\oint_{\partial S}\boldsymbol E\cdot d\boldsymbol x$が生じ、その大きさは$-d\Phi/dt$に等しくなります。負号は、誘導された電流の作る磁場が磁束の変化を打ち消す向きになることを表します（レンツの法則）。

静電場では$\oint\boldsymbol E\cdot d\boldsymbol x=0$でしたが、磁場が時間変化すると電場は渦を持ちます。微分形では、静電場の渦なしの式$\nabla\times\boldsymbol E=0$が

$$
\nabla\times\boldsymbol E=-\partial_t\boldsymbol B
$$

に置き換わります。これを**ファラデーの法則**と呼びます。

## マクスウェル方程式

ガウスの法則と$\nabla\cdot\boldsymbol B=0$はそのまま保たれます。$\mu_0\varepsilon_0=1/c^2$を使って4本をまとめたものが、SI単位系のマクスウェル方程式です。

&&&fml マクスウェル方程式 [fml-maxwell]
$$
\nabla\cdot\boldsymbol E=\frac\rho{\varepsilon_0},\qquad
\nabla\times\boldsymbol B-\frac1{c^2}\partial_t\boldsymbol E=\mu_0\boldsymbol J,\qquad
\nabla\times\boldsymbol E+\partial_t\boldsymbol B=0,\qquad
\nabla\cdot\boldsymbol B=0
$$
&&&

時間微分が加わったのは、2本の回転の式だけです。電場の時間微分は磁場の回転の式に、磁場の時間微分は電場の回転の式に入り、互いに相手の回転と組になっています。

# 1本の式

$\partial_0$はスカラーの作用素なので、$F$に作用させてもグレードを変えません。$\partial_0F=\partial_0\boldsymbol E+Ic\,\partial_0\boldsymbol B$のベクトル部と2ベクトル部が、$DF$の同じグレードの項に加わります。$DF$の成分は[[7shi-em1]]で見たとおりなので、$\mathcal DF$の成分は次のようになります。

&&&fml $\mathcal DF$のグレード成分 [fml-grades]
$$
\mathcal D(\boldsymbol E+Ic\boldsymbol B)
=\nabla\cdot\boldsymbol E
+\bigl(\partial_0\boldsymbol E-c\,\nabla\times\boldsymbol B\bigr)
+I\bigl(\nabla\times\boldsymbol E+c\,\partial_0\boldsymbol B\bigr)
+Ic\,\nabla\cdot\boldsymbol B
$$
&&&

時間微分は、静的な場合のアンペールの法則の位置（ベクトル部）に$\partial_0\boldsymbol E$を、渦なしの式の位置（2ベクトル部）に$c\,\partial_0\boldsymbol B$を加えます。これはマクスウェルが加えた項と、ファラデーの法則が加えた項そのものです。

&&&thm マクスウェル方程式の1本の式 [thm-maxwell]
$$
\mathcal DF=\frac1{\varepsilon_0}\Bigl(\rho-\frac{\boldsymbol J}c\Bigr)
$$

は[マクスウェル方程式](#fml-maxwell)の4本の式と同値です。
&&&

&&&prf
[$\mathcal DF$のグレード成分](#fml-grades)と右辺をグレードごとに比べる。スカラー部は$\nabla\cdot\boldsymbol E=\rho/\varepsilon_0$である。ベクトル部は$\partial_0\boldsymbol E-c\,\nabla\times\boldsymbol B=-\boldsymbol J/\varepsilon_0c$であり、$\partial_0=\frac1c\partial_t$として両辺を$-c$で割り、$1/\varepsilon_0c^2=\mu_0$を使うと$\nabla\times\boldsymbol B-\frac1{c^2}\partial_t\boldsymbol E=\mu_0\boldsymbol J$となる。右辺には2ベクトル部と擬スカラー部がないので、$\nabla\times\boldsymbol E+\partial_t\boldsymbol B=0$と$\nabla\cdot\boldsymbol B=0$を得る。逆に4本の式が成り立てば、各グレードが一致するので1本の式が成り立つ。
&&&

| グレード | 式 | 静的な場合（[[7shi-em1]]） |
|---|---|---|
| スカラー | $\nabla\cdot\boldsymbol E=\rho/\varepsilon_0$ | 同じ |
| ベクトル | $\nabla\times\boldsymbol B-\frac1{c^2}\partial_t\boldsymbol E=\mu_0\boldsymbol J$ | $\nabla\times\boldsymbol B=\mu_0\boldsymbol J$ |
| 2ベクトル | $\nabla\times\boldsymbol E+\partial_t\boldsymbol B=0$ | $\nabla\times\boldsymbol E=0$ |
| 擬スカラー | $\nabla\cdot\boldsymbol B=0$ | 同じ |

場が時間によらなければ$\mathcal DF=DF$であり、[[7shi-em1]]の静的な式に戻ります。

積分形も[[7shi-em1]]と同じく基本定理のグレード成分として得られます。時間によらない曲面$S$について曲面の基本定理を$F$に使うと、スカラー部と擬スカラー部の右辺は$\boldsymbol n\cdot(\nabla\times\boldsymbol E)$と$Ic\,\boldsymbol n\cdot(\nabla\times\boldsymbol B)$の面積分で、そこに2ベクトル部とベクトル部の式を代入します。

&&&fml 回転の2本の積分形
$$
\oint_{\partial S}\boldsymbol E\cdot d\boldsymbol x=-\frac d{dt}\int_S\boldsymbol B\cdot\boldsymbol n\,dA,\qquad
\oint_{\partial S}\boldsymbol B\cdot d\boldsymbol x=\mu_0\int_S\boldsymbol J\cdot\boldsymbol n\,dA+\frac1{c^2}\frac d{dt}\int_S\boldsymbol E\cdot\boldsymbol n\,dA
$$
&&&

左の式が電磁誘導の法則そのもの、右の式が変位電流を含むアンペールの法則の積分形です。領域の基本定理からは、ガウスの法則と閉曲面を通る磁束が$0$であることが[[7shi-em1]]と同じ形で得られます。

&&&rem 四元数からベクトル解析へ、そして再び1本の式へ
[[7shi-hist]]で見たとおり、マクスウェルは1873年の著書で四元数の記法を併記していましたが、ギブスとヘヴィサイドは四元数の積を内積とベクトル積に分け、電磁気学の法則を発散と回転の式として書き直しました。上の4本の式はその形です。幾何代数では、ベクトル値の関数に$D$を作用させると、グレードを下げる部分に発散が、上げる部分に回転（の双対）が同時に現れます。2ベクトル値の関数についても、同じ幾何積の微分によって発散と回転がまとめて得られます。分けられた発散と回転を1つの作用素に戻すと、4本の式は1本の式のグレード成分になります。通常の四元数の記法では値を実部と3成分の虚部に分けるだけですが、幾何代数ではベクトルと2ベクトルが別のグレードとして区別され、$\boldsymbol E$と$\boldsymbol B$を1つの元の別の成分として持てます。
&&&

# 係数の2乗の符号

[[7shi-cla4]]では、$\mathcal D$と共役$\bar{\mathcal D}$の積が4次元のラプラシアンになりました。同じ計算を$\operatorname{Cl}_{3,0}(\mathbb R)$で行います。

&&&fml 共役との積
係数$a_1,a_2,a_3$が座標によらない定数で、互いに反可換（$l\ne m$なら$a_la_m=-a_ma_l$）で、$a_l^2$がスカラーのとき

$$
\Bigl(\partial_0-\sum_la_l\partial_l\Bigr)\Bigl(\partial_0+\sum_la_l\partial_l\Bigr)=\partial_0^2-\sum_la_l^2\,\partial_l^2
$$

です。
&&&

&&&prf
展開すると、$\partial_0\sum_la_l\partial_l$の項は打ち消し合い、$-\sum_{l,m}a_la_m\partial_l\partial_m$が残る。$l\ne m$の項は$a_la_m+a_ma_l=0$により対ごとに消え、$l=m$の項$-a_l^2\partial_l^2$だけが残る。積の順序を逆にしても、時間微分と空間微分の交差項が同様に打ち消し合い、同じ結果になる。
&&&

[[7shi-cla4]]の係数$h_l=e_0e_l$（四元数の$i,j,k$）は$h_l^2=-1$を満たし、右辺は$\partial_0^2+\partial_1^2+\partial_2^2+\partial_3^2$、4次元のラプラシアンです。本記事の係数$e_k$は$e_k^2=+1$を満たし、右辺は

$$
\bar{\mathcal D}\mathcal D=\mathcal D\bar{\mathcal D}=\partial_0^2-\Delta
$$

となります。$\Delta$は3次元のラプラシアンです。$\partial_0^2-\Delta=\frac1{c^2}\partial_t^2-\Delta$は速さ$c$の波動作用素です。

| | 係数 | 係数の2乗 | 共役との積 | 型 |
|---|---|---|---|---|
| 四元数解析（[[7shi-cla4]]） | $h_l=e_0e_l$ | $-1$ | $\partial_0^2+\Delta$ | 楕円型 |
| 本記事 | $e_k$ | $+1$ | $\partial_0^2-\Delta$ | 双曲型 |

作用素の形は同じで、この積の計算に現れる違いは空間方向の係数の2乗の符号だけです。係数の2乗が$-1$なら共役との積がラプラシアンになり、$+1$なら波動作用素になります。同じ形の1階作用素が、前者では正則関数の理論に、後者では電磁場の波動方程式に結びつきます。

# 連続の式と波動方程式

[1本の式](#thm-maxwell)に左から$\bar{\mathcal D}$を掛けます。左辺は$(\partial_0^2-\Delta)F$です。$\partial_0^2-\Delta$はスカラーの作用素なので、ベクトル部と2ベクトル部しか持ちません。右辺を計算すると次のようになります。

$$
\bar{\mathcal D}\Bigl(\rho-\frac{\boldsymbol J}c\Bigr)
=\Bigl(\partial_0\rho+\frac1c\nabla\cdot\boldsymbol J\Bigr)
-\Bigl(\nabla\rho+\frac1c\partial_0\boldsymbol J\Bigr)
+\frac Ic\nabla\times\boldsymbol J
$$

ここで$D\boldsymbol J=\nabla\cdot\boldsymbol J+I\,\nabla\times\boldsymbol J$を使いました。

## 連続の式

$F$はベクトル部と2ベクトル部だけを持つので、左辺にスカラー部はありません。したがって右辺のスカラー部も$0$でなければなりません。$\partial_0=\frac1c\partial_t$より

&&&fml 連続の式
$$
\partial_t\rho+\nabla\cdot\boldsymbol J=0
$$
&&&

が得られます。[1本の式](#thm-maxwell)の解については、電荷の保存を独立に仮定する必要はありません。逆に、源$\rho,\boldsymbol J$を与えて場を求める立場では、連続の式は解が存在するための必要条件です。ベクトル解析では、ガウスの法則の時間微分とアンペール＝マクスウェルの法則の発散を組み合わせて導く式です。ここでは$\bar{\mathcal D}\mathcal D$がスカラーの作用素であることから、グレードの数え上げで出ます。

## 波動方程式

ベクトル部と2ベクトル部からは、$\boldsymbol E$と$\boldsymbol B$が個別に満たす方程式が得られます。$\partial_0=\frac1c\partial_t$と$1/\varepsilon_0c^2=\mu_0$を使って整理します。

&&&fml 電場と磁場の波動方程式
$$
\Bigl(\frac1{c^2}\partial_t^2-\Delta\Bigr)\boldsymbol E=-\frac1{\varepsilon_0}\nabla\rho-\mu_0\partial_t\boldsymbol J,\qquad
\Bigl(\frac1{c^2}\partial_t^2-\Delta\Bigr)\boldsymbol B=\mu_0\nabla\times\boldsymbol J
$$
&&&

電荷も電流もない真空では右辺が$0$になり、$\boldsymbol E$と$\boldsymbol B$の各成分は速さ$c$で伝わる波動方程式を満たします。電磁場の各成分が速さ$c$の波動方程式を満たすことが、1階の方程式$\mathcal DF=0$に共役$\bar{\mathcal D}$を掛けるという1回の操作から分かります。逆に、各成分が波動方程式を満たすだけでは、元のマクスウェル方程式を満たすとは限りません。

&&&rem 光は電磁波
波動方程式は、ある場所の乱れが一定の速さで周囲に伝わっていくことを表します。真空のマクスウェル方程式では、磁場の時間変化が電場の渦を生み（ファラデーの法則）、電場の時間変化が磁場の渦を生みます（変位電流）。この連鎖によって、電荷や電流から離れた場所でも電場と磁場が互いを支えながら空間を伝わっていきます。これが**電磁波**です。

伝わる速さ$c$は、[[7shi-em1]]の関係$\varepsilon_0\mu_0c^2=1$で定めた$c=1/\sqrt{\varepsilon_0\mu_0}$です。$\varepsilon_0$と$\mu_0$は静止した電荷の間の力と電流の間の力という、光とは無関係に見える実験から決まる定数ですが、そこから計算した電磁波の速さは、測定されていた光の速さと一致しました。マクスウェルはこれをもとに、光は電磁波であると結論しました（1864年）。電磁波の存在は、1888年にヘルツが電気回路からの電磁波の発生と検出によって実験的に確かめています。現在では、電波・赤外線・可視光・紫外線・X線などは、振動数（波長）の違う電磁波として統一的に理解されています。

[[7shi-em1]]では、$\boldsymbol E$と$c\boldsymbol B$の単位が揃うことを、電場と磁場をまとめる理由として使いました。そこで現れた$c$の正体が、電磁波の速さです。
&&&

[[7shi-cla4]]では、同じ操作で正則関数の各成分が調和関数になることを導きました。操作は同じで、係数の2乗の符号が違うために、結論が「調和関数」から「波動方程式の解」に替わります。

## 静的な場合

場と源が時間によらないとき、$\bar{\mathcal D}\mathcal D=-\Delta$であり、上の式は$\Delta\boldsymbol E=\nabla\rho/\varepsilon_0$、$\Delta\boldsymbol B=-\mu_0\nabla\times\boldsymbol J$です。連続の式は定常電流の条件$\nabla\cdot\boldsymbol J=0$になり、[[7shi-em1]]でスカラー部が消えるために必要だった条件と一致します。

[[7shi-em1]]では、静的な場を$D$の基本解による積分公式で表しました。時間を含む場合に$\mathcal D$の基本解を求めて同じ形の公式を作ることは、本記事では扱いません。

# まとめ

時間を$x_0=ct$として作用素$\mathcal D=\partial_0+D$を作ると、マクスウェル方程式は$\mathcal DF=(\rho-\boldsymbol J/c)/\varepsilon_0$の1本の式になりました。

- **グレード成分**：スカラー部がガウスの法則、ベクトル部がアンペール＝マクスウェルの法則、2ベクトル部がファラデーの法則、擬スカラー部が$\nabla\cdot\boldsymbol B=0$です。時間微分はスカラーの作用素なのでグレードを変えず、静的な式の回転の2本にだけ項を加えます。
- **係数の2乗の符号**：$\mathcal D$は[[7shi-cla4]]の作用素と同じ形で、係数の2乗が$-1$から$+1$に替わっています。共役との積は$\partial_0^2+\Delta$（楕円型）から$\partial_0^2-\Delta$（双曲型）に替わります。
- **連続の式と波動方程式**：$\bar{\mathcal D}$を掛けると、左辺にスカラー部がないことから連続の式が、ベクトル部と2ベクトル部から$\boldsymbol E$と$\boldsymbol B$の波動方程式が出ます。
