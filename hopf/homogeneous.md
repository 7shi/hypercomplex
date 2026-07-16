実四次元座標から複素数のペアの比（同次座標）を求めることで、ホップファイブレーション[[hopf1931]]を計算します。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

# 概要

三次元球面$S^3$上の点$(x_1, x_2, x_3, x_4)\in\mathbb{R}^4$は条件$x_1^2 + x_2^2 + x_3^2 + x_4^2=1$を満たします。

この成分で表される複素数のペアの比$z = \dfrac{x_1 + i x_2}{x_3 + i x_4}\in\mathbb{C}$から二次元球面$S^2$上へ射影した座標$(\xi_1, \xi_2, \xi_3)\in\mathbb{R}^3$を求めます。

# 複素平面への射影

実四次元上の座標を複素二次元に写像します。

$$
\mathbb{R}^4\ni(x_1, x_2, x_3, x_4)\mapsto(x_1+ix_2,x_3+ix_4)\in\mathbb{C}^2
$$

これを同次座標（比として表した座標）に射影します。[[7shi-mobius]]

$$
(x_1+ix_2,x_3+ix_4)\mapsto[x_1+ix_2:x_3+ix_4]
$$

$x_3+ix_4\ne0$の場合、同次座標の代表点は成分の商$z$として表せます。

$$
z=\frac{x_1+ix_2}{x_3+ix_4}
$$

$z$の分母を実数化して、実部と虚部に分離します。

$$
\begin{aligned}
z
&= \frac{(x_1 + i x_2)(x_3 - i x_4)}{(x_3 + i x_4)(x_3 - i x_4)} \\
&= \frac{(x_1 x_3 + x_2 x_4) + i (x_2 x_3 - x_1 x_4)}{x_3^2 + x_4^2} \\
&= \frac{x_1 x_3 + x_2 x_4}{x_3^2 + x_4^2} + i\cdot\frac{x_2 x_3 - x_1 x_4}{x_3^2 + x_4^2} \\
\end{aligned}
$$

実部と虚部の係数を $u,v \in \mathbb{R}$ とおきます。

$$
u = \frac{x_1 x_3 + x_2 x_4}{x_3^2 + x_4^2}, \quad v = \frac{x_2 x_3 - x_1 x_4}{x_3^2 + x_4^2}
$$
$$
z = u + i v
$$

# 二次元球面への射影

$z$を実三次元空間の$xy$平面に埋め込みます。

$$
\mathbb{C} \ni z=u+iv \hookrightarrow (u,v,0) \in \mathbb{R}^3
$$

$(u,v,0)$を二次元球面$S^2$上に射影するため、$S^2$の北極点$(0,0,1)$と$(u,v,0)$を結ぶ直線を媒介変数$t$で表します。

$$
(x,y,z)=(ut,vt,1-t)
\qquad\left\{
\begin{aligned}
&(0,0,1)&(t&=0) \\
&(u,v,0)&(t&=1) \\
\end{aligned}
\right.
$$
これを単位球面 $x^2 + y^2 + z^2 = 1$ の式に代入して交点を求めます。
$$
\begin{aligned}
(u t)^2 + (v t)^2 + (1 - t)^2 &= 1 \\
(u^2 + v^2 + 1) t^2 - 2 t &= 0 \\
t\{(u^2 + v^2 + 1)t-2\} &= 0
\end{aligned}
$$
$$
t = 0,\ \frac{2}{u^2 + v^2 + 1}
$$

$t=\dfrac{2}{u^2 + v^2 + 1}$ のときの座標を$(\xi_1, \xi_2, \xi_3)\in\mathbb{R}^3$とします。

$$
(\xi_1, \xi_2, \xi_3)
=\left(
  \frac{2 u}{u^2 + v^2 + 1},
  \frac{2 v}{u^2 + v^2 + 1},
  1 - \frac{2}{u^2 + v^2 + 1}
\right)
$$

分母に現れる$u^2+v^2+1$を計算します。$z$の複素共役を$z^*$とします。

$$
\begin{aligned}
u^2+v^2+1
&=(u+iv)(u-iv)+1 \\
&=zz^*+1 \\
&=\left(\frac{x_1+ix_2}{x_3+ix_4}\right)\left(\frac{x_1-ix_2}{x_3-ix_4}\right)+1 \\
&=\frac{x_1^2 + x_2^2}{x_3^2 + x_4^2}+1 \\
&=\frac{x_1^2 + x_2^2 + x_3^2 + x_4^2}{x_3^2 + x_4^2} \\
&=\frac{1}{x_3^2 + x_4^2}
\end{aligned}
$$

$\xi_1, \xi_2, \xi_3$を$x_1, x_2, x_3, x_4$で表現します。

$$
\begin{aligned}
\xi_1 &= \frac{2 u}{u^2 + v^2 + 1} = 2 \left( \frac{x_1 x_3 + x_2 x_4}{x_3^2 + x_4^2} \right) (x_3^2 + x_4^2) = 2 (x_1 x_3 + x_2 x_4) \\
\xi_2 &= \frac{2 v}{u^2 + v^2 + 1} = 2 \left( \frac{x_2 x_3 - x_1 x_4}{x_3^2 + x_4^2} \right) (x_3^2 + x_4^2) = 2 (x_2 x_3 - x_1 x_4) \\
\xi_3 &= 1 - \frac{2}{u^2 + v^2 + 1} = 1 - 2(x_3^2 + x_4^2) = x_1^2 + x_2^2 - x_3^2 - x_4^2
\end{aligned}
$$

この結果はHopfの原論文§5.(1)と一致します。[[hopf1931]]

$x_3+ix_4 \rightarrow 0$の極限で$(\xi_1,\xi_2,\xi_3)\rightarrow(0,0,1)$となることから、$x_3+ix_4=0$は北極点に対応付けます。

# まとめ

$$
\mathbb{R}^4 \supset S^3 \ni
\begin{pmatrix}x_1 \\ x_2 \\ x_3 \\ x_4\end{pmatrix}
\mapsto
\begin{pmatrix}\xi_1 \\ \xi_2 \\ \xi_3\end{pmatrix}
=\begin{pmatrix}
  2 (x_1 x_3 + x_2 x_4) \\
  2 (x_2 x_3 - x_1 x_4) \\
  x_1^2 + x_2^2 - x_3^2 - x_4^2
 \end{pmatrix}
\in S^2 \subset \mathbb{R}^3
$$

# 複素数のペアでの計算

$(x_1, x_2, x_3, x_4)$から構成される複素数のペアを$\alpha=x_1+ix_2,\ \beta=x_3+ix_4$ として計算します。

$$
z=\frac{\alpha}{\beta}=\frac{\alpha\beta^*}{\beta\beta^*},\quad
u=\mathrm{Re}(z),\quad
v=\mathrm{Im}(z)
$$
$$
\begin{pmatrix}\xi_1 \\ \xi_2 \\ \xi_3\end{pmatrix}
=\begin{pmatrix}
  2 \mathrm{Re}(\alpha \beta^*) \\
  2 \mathrm{Im}(\alpha \beta^*) \\
  \alpha \alpha^* - \beta \beta^*
 \end{pmatrix}
=\begin{pmatrix}\begin{aligned}
  \alpha \beta^* &+ \beta \alpha^* \\
  -i(\alpha \beta^* &- \beta \alpha^*) \\
  \alpha \alpha^* &- \beta \beta^*
 \end{aligned}\end{pmatrix}
$$

&&&rem 四元数
この結果は四元数による導出とほぼ同じですが、第2成分の符号が異なります。[[7shi-qhopf1]]

$$
ω\mathbf{k}ω^*
=(α^*β+β^*α)\mathbf{i}-i(α^*β-β^*α)\mathbf{j}+(α^*α-β^*β)\mathbf{k}
$$

四元数と同じ結果を得るには、$xy$平面への埋め込み方を変更します。
$$
\mathbb{C} \ni z=u+iv \hookrightarrow (u,\textcolor{red}{-v},0) \in \mathbb{R}^3
$$
&&&

複素ベクトルの変換として表記すれば、実部と虚部に分離する必要がなくなります。

$$
\begin{pmatrix}\alpha \\ \beta\end{pmatrix}
\mapsto
\begin{pmatrix}\xi_1 + i\xi_2 \\ \xi_3\end{pmatrix}
=\begin{pmatrix}
  2 \alpha \beta^* \\
  \alpha \alpha^* - \beta \beta^*
 \end{pmatrix}
$$

&&&rem
$\alpha\beta^*$の形は、分母を実数化した$z=\dfrac{\alpha\beta^*}{\beta\beta^*}$の分子に表れています。
&&&
