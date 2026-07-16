レイリー・ジーンズの公式におけるエネルギーを量子化すればプランクの公式が得られます。また、プランクの公式をすべての周波数に渡って積分すればゼータ関数が現れます。

# 黒体放射を表す3つの公式

黒体放射を数式で表す試みは、当初、低周波数と高周波数の領域が別の式で表されていました。

&&&fml レイリー・ジーンズの公式（低周波数）
$$
u(\nu, T) = \frac{8\pi\nu^2}{c^3} kT
$$
&&&

&&&fml ウィーンの公式（高周波数）
$$
u(\nu, T) = \frac{8\pi\nu^3}{c^3} \frac{1}{e^{h\nu/kT}}
$$
&&&

これらの公式を1つにまとめたのがプランクの公式です。

&&&fml プランクの公式（全周波数）
$$
u(\nu, T) = \frac{8\pi h\nu^3}{c^3} \frac{1}{e^{h\nu/kT} - 1}
$$
&&&

ここで$u(\nu, T)$は、周波数$\nu$、絶対温度$T$（ケルビン）における分光エネルギー密度を表します。また、$k$はボルツマン定数、$c$は光速度、$h$はプランク定数です。

## 3つの公式の関係

低周波極限と高周波極限において、プランクの公式から他の公式が得られます。

&&&prf レイリー・ジーンズの公式
低周波極限$h\nu \ll kT$では、$e^{h\nu/kT} \approx 1 + h\nu/kT$より

$$
u(\nu, T)
= \frac{8\pi h\nu^3}{c^3} \frac{1}{e^{h\nu/kT} - 1}
\approx \frac{8\pi h\nu^3}{c^3} \frac{1}{h\nu/kT}
= \frac{8\pi\nu^2}{c^3} kT
$$
&&&

&&&prf ウィーンの公式
高周波極限$h\nu \gg kT$では、$e^{h\nu/kT} - 1 \approx e^{h\nu/kT}$より

$$
u(\nu, T)
= \frac{8\pi h\nu^3}{c^3} \frac{1}{e^{h\nu/kT} - 1}
\approx \frac{8\pi h\nu^3}{c^3} \frac{1}{e^{h\nu/kT}}
$$
&&&

# 3つの公式の一般構造

これらの公式は共通の構造を持っています。

$$
u(\nu,T) = \frac{8\pi\nu^2}{c^3} \langle E \rangle
$$

この式は2つの部分からなります。

1. **モード密度**：$\dfrac{8\pi\nu^2}{c^3}$
周波数$\nu$における単位体積当たりのモード（振動パターン）の数を表します。

2. **平均エネルギー**：$\langle E \rangle$
周波数$\nu$におけるモードが持つ平均エネルギーです。

## モード密度の導出

モード密度は、3次元空間内に閉じ込められた電磁波の定常波の数を数えることで導出されます。一辺$L$の立方体容器内では波数ベクトルが$\pi/L$刻みで離散化されます。境界条件によって波の形が制限されるため、物理的に意味があるのは波数ベクトル成分が正（$k_x, k_y, k_z \geq 0$）となる第1象限のみであるため、全空間の$1/8$となります。

波数空間で波数の大きさが$k$と$k+dk$の間にある状態数は、球殻の体積要素

$$
d\left(\frac{4}{3}\pi k^3\right) = 4\pi k^2 dk
$$

を状態1つあたりの体積$(\pi/L)^3$と第1象限への制限$8$で割った値になります。電磁波には2つの独立した偏光自由度があるため、モード数は

$$
4\pi k^2 dk \cdot \frac{1}{(\pi/L)^3} \cdot \frac{1}{8} \cdot 2
=\frac{L^3 k^2}{\pi^2} dk
$$

となります。これを体積$L^3$で割って単位体積あたりに直し、関係式$k=2\pi\nu/c$を用いて周波数に変換すれば、$d\nu$の係数としてモード密度が得られます。

$$
\frac{k^2}{\pi^2} dk
=\left(\frac{2\pi\nu}{c}\right)^2 \frac{1}{\pi^2} \left(\frac{2\pi}{c}d\nu\right)
=\frac{8\pi\nu^2}{c^3} d\nu
$$

## 平均エネルギーの違い

各公式の違いは、平均エネルギー$\langle E \rangle$にあります。

- **レイリー・ジーンズの公式**：$\langle E \rangle = kT$

- **ウィーンの公式**：$\langle E \rangle = \dfrac{h\nu}{e^{h\nu/kT}}$

- **プランクの公式**：$\langle E \rangle = \dfrac{h\nu}{e^{h\nu/kT}-1}$

前述のように、プランクの公式の高周波極限がウィーンの公式です。

レイリー・ジーンズの公式とプランクの公式の違いは、単に低周波極限というだけに留まらず、取り得るエネルギーの違いに由来します。

レイリー・ジーンズの公式の平均エネルギーは、古典統計力学による振動子の平均エネルギーから求められます。

$$
\langle E \rangle
= \frac{\int_0^{\infty} E \cdot e^{-E/kT} dE}{\int_0^{\infty} e^{-E/kT} dE}
= kT
$$

プランクの公式の平均エネルギーは、エネルギーを離散的な値$E_n=nh\nu$（$n = 0,1,2,\ldots$）に制限することにより（**エネルギーの量子化**）、計算を積分から総和（シグマ）に変えることで求められます。

$$
\langle E \rangle
= \frac{\sum_{n=0}^{\infty} E_n \cdot e^{-E_n/kT}}{\sum_{n=0}^{\infty} e^{-E_n/kT}}
= \dfrac{h\nu}{e^{h\nu/kT}-1}
$$

歴史的には、プランクはまず実験データを説明できる式を見つけ、その後、この式を理論的に導出するために量子化という概念を導入しました。プランクは当初、$h$を実験データへの当てはめパラメーターと見なしており、その物理的意味はアインシュタインが光量子仮説を提案するまで完全には認識されませんでした。

## 連続的なエネルギー

レイリー・ジーンズの公式の平均エネルギーは、統計力学におけるボルツマン分布によって求められます。

統計力学では、熱平衡状態にある系においてエネルギー$E$を持つ状態が現れる確率は、ボルツマン因子$e^{-E/kT}$に比例します。これがボルツマン分布の基本的な考え方です。

$$
P(E) \propto e^{-E/kT}
$$

この表現を確率分布として使うため、全確率が$1$になるように規格化します。

$$
P(E) = \frac{e^{-E/kT}}{\int_0^{\infty} e^{-E/kT} dE}
$$

物理量の平均値は、その量に確率を掛けて全状態にわたって積分することで求められます。これにより平均エネルギー$\langle E \rangle$が求まります。

$$
\langle E \rangle
= \int_0^{\infty} E \cdot P(E) dE
= \frac{\int_0^{\infty} E \cdot e^{-E/kT} dE}{\int_0^{\infty} e^{-E/kT} dE}
$$

変数変換$x = \dfrac{E}{kT}$より、分母の積分を計算します。

$$
\int_0^{\infty} e^{-E/kT} dE
= kT \int_0^{\infty} e^{-x} dx
= kT [-e^{-x}]_0^{\infty}
= kT
$$

同じ変数変換で分子の積分も計算します。

$$
\int_0^{\infty} E \cdot e^{-E/kT} dE
= \int_0^{\infty} xkT \cdot e^{-x} kT \, dx
= (kT)^2 \int_0^{\infty} x \cdot e^{-x} dx
$$

積分の部分は、部分積分で計算できます。表面項（中辺第1項）が消える（$0$になる）のに注意してください。

$$
\int_0^{\infty} x \cdot e^{-x} du
= [x \cdot (-e^{-x})]_0^{\infty} + \int_0^{\infty} e^{-x} dx
= 0+[-e^{-x}]_0^{\infty}
= 1
$$

ここまでの結果を$\langle E \rangle$の式に代入します。

$$
\langle E \rangle
= \frac{\int_0^{\infty} E \cdot e^{-E/kT} dE}{\int_0^{\infty} e^{-E/kT} dE}
= \frac{(kT)^2}{kT}=kT
$$

## 離散的なエネルギー

プランクの公式の平均エネルギーを計算します。

$$
\langle E \rangle
= \frac{\sum_{n=0}^{\infty} E_n \cdot e^{-E_n/kT}}{\sum_{n=0}^{\infty} e^{-E_n/kT}}
= \frac{\sum_{n=0}^{\infty} nh\nu \cdot e^{-nh\nu/kT}}{\sum_{n=0}^{\infty} e^{-nh\nu/kT}}
$$

$x = e^{-h\nu/kT}\ (0<x<1)$とおいて、分母について等比数列の和の公式を適用します。

$$
\sum_{n=0}^{\infty} e^{-nh\nu/kT}
= \sum_{n=0}^{\infty} x^n
= \frac{1}{1-x}
$$

分子を整理します。

$$
\sum_{n=0}^{\infty} nh\nu \cdot e^{-nh\nu/kT}
= h\nu \sum_{n=0}^{\infty} n x^n
$$

総和の部分は、$x\dfrac{d}{dx}x^n=nx^n$より計算できます。

$$
\sum_{n=0}^{\infty} n x^n
= x \frac{d}{dx}\sum_{n=0}^{\infty} x^n
= x \frac{d}{dx}\frac{1}{1-x}
= \frac{x}{(1-x)^2}
$$

ここまでの結果を$\langle E \rangle$の式に代入します。

$$
\langle E \rangle
= h\nu \left( \frac{x}{1-x} \right)
= h\nu \left( \frac{1}{\frac{1}{x}-1} \right)
= \frac{h\nu}{e^{h\nu/kT}-1}
$$

# 全エネルギー密度とゼータ関数

プランクの公式は、周波数$\nu$における単位体積あたりのエネルギー密度を表します。

$$
u(\nu, T) = \frac{8\pi h\nu^3}{c^3} \frac{1}{e^{h\nu/kT} - 1}
$$

これをすべての周波数に渡って積分することで全エネルギー密度$u(T)$が得られ、ゼータ関数が現れます。

$$
u(T)
= \int_0^\infty u(\nu, T) d\nu
= \frac{48\pi k^4 T^4}{c^3 h^3} \zeta(4)
$$

## 全エネルギー密度の計算

この積分を計算するための変数変換として$x = \dfrac{h\nu}{kT}$とおけば、$\nu = \dfrac{kT}{h}x$となります。

$$
u(T)
= \int_0^\infty \frac{8\pi h}{c^3} \left(\dfrac{kT}{h}x\right)^3 \frac{1}{e^x - 1} d\left(\frac{kT}{h}x\right)
= \frac{8\pi k^4 T^4}{c^3 h^3} \int_0^\infty \frac{x^3}{e^x - 1} dx
$$

積分の中を計算します。

$$
\frac{1}{e^x - 1}
= \frac{1}{e^x(1 - e^{-x})}
= \frac{e^{-x}}{1 - e^{-x}}
= e^{-x} \sum_{n=0}^{\infty} e^{-nx}
= \sum_{n=1}^{\infty} e^{-nx}
$$

これを積分に戻して計算を進めます。途中、積分と総和の順序が交換できることを利用します。

$$
\int_0^\infty \frac{x^3}{e^x - 1} dx
= \int_0^\infty x^3 \sum_{n=1}^{\infty} e^{-nx} dx
= \sum_{n=1}^{\infty} \int_0^\infty x^3 e^{-nx} dx
$$

積分の部分に、部分積分を3回適用します。表面項はすべて消えます。

$$
\begin{alignedat}{2}
\int_0^\infty x^3 e^{-nx} dx
 &= 0 - \int_0^\infty 3x^2 \left(-\frac{1}{n}e^{-nx}\right) dx
&&= \frac{3}{n}\int_0^\infty x^2 e^{-nx} dx
\\
\int_0^\infty x^2 e^{-nx} dx
 &= 0 - \int_0^\infty 2x \left(-\frac{1}{n}e^{-nx}\right) dx
&&= \frac{2}{n}\int_0^\infty x e^{-nx} dx
\\
\int_0^\infty x e^{-nx} dx
 &= 0 - \int_0^\infty \left(-\frac{1}{n}e^{-nx}\right) dx
&&= \frac{1}{n}\int_0^\infty e^{-nx} dx
\\
\end{alignedat}
$$

最後の積分は直接計算できます。

$$
\int_0^\infty e^{-nx} dx
= \left[ -\frac{1}{n}e^{-nx} \right]_0^\infty
= \frac{1}{n}
$$

この結果から逆にたどります。

$$
\begin{alignedat}{2}
\int_0^\infty x e^{-nx} dx &= \frac{1}{n} \cdot \frac{1}{n} &&= \frac{1}{n^2} \\
\int_0^\infty x^2 e^{-nx} dx &= \frac{2}{n} \cdot \frac{1}{n^2} &&= \frac{2}{n^3} \\
\int_0^\infty x^3 e^{-nx} dx &= \frac{3}{n} \cdot \frac{2}{n^3} &&= \frac{6}{n^4} \\
\int_0^\infty \frac{x^3}{e^x - 1} dx &= \sum_{n=1}^{\infty} \frac{6}{n^4} &&= 6\zeta(4)
\end{alignedat}
$$

ここでリーマンゼータ関数が現れました。

&&&def リーマンゼータ関数
$$
\zeta(s) = \sum_{n=1}^{\infty} \frac{1}{n^s}
$$
&&&

$\zeta(4) = \dfrac{\pi^4}{90}$であることが知られているため（後述）、$u(T)$の式に代入します。

$$
u(T)
= \frac{8\pi k^4 T^4}{c^3 h^3} \int_0^\infty \frac{x^3}{e^x - 1} dx
= \frac{8\pi k^4 T^4}{c^3 h^3} \cdot 6 \cdot \frac{\pi^4}{90}
= \frac{8\pi^5 k^4 T^4}{15 c^3 h^3}
$$

これはシュテファン=ボルツマンの法則に対応します。

&&&fml シュテファン=ボルツマンの法則
$$
u(T) = \dfrac{4\sigma}{c}T^4
$$
&&&

&&&def シュテファン=ボルツマン定数
$$
\sigma=\dfrac{2\pi^5 k^4}{15 c^2 h^3}
$$
&&&

## ゼータ関数の計算

$\zeta(4) = \dfrac{\pi^4}{90}$を示すには、フーリエ級数展開によるパーセヴァルの等式を用いる方法が一般的です。

区間$[- \pi, \pi]$で定義された関数$x^2$を考えます。これは偶関数なので、フーリエ級数は余弦項のみで表されます。

$$
\begin{aligned}
x^2
&= \frac{1}{2\pi} \int_{-\pi}^\pi x^2 \, dx + \sum_{n=1}^\infty \frac{1}{\pi} \int_{-\pi}^\pi x^2 \cos(nx) \, dx \\
&= \frac{\pi^2}{3} + \sum_{n=1}^\infty \frac{4 (-1)^n}{n^2} \cos(nx)
\end{aligned}
$$

パーセヴァルの等式を用いて、フーリエ係数から$L^2$ノルムを結びつけます。

$$
\begin{aligned}
\frac{1}{2\pi} \int_{-\pi}^\pi |x^2|^2 \, dx
&= \left|\frac{\pi^2}{3}\right|^2 + \frac{1}{2} \sum_{n=1}^\infty \left|\frac{4 (-1)^n}{n^2}\right|^2
\\
\frac{1}{2\pi} \int_{-\pi}^\pi x^4 \, dx
&= \frac{\pi^4}{9} + \frac{1}{2} \sum_{n=1}^\infty \frac{16}{n^4}
\\
\frac{\pi^4}{5} &= \frac{\pi^4}{9} + 8\zeta(4) \\
8\zeta(4) &= \frac{\pi^4}{5} - \frac{\pi^4}{9} \\
\zeta(4) &= \frac{\pi^4}{8}\left(\frac{1}{5}-\frac{1}{9}\right)=\frac{\pi^4}{90}
\end{aligned}
$$
