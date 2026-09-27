実4次元座標から複素数のペアの比（同次座標）を構成し、立体射影を経由してホップファイブレーションを導出します。

シリーズ：[ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

# 概要

3次元球面$S^3$から2次元球面$S^2$への非自明な写像である**ホップファイブレーション**は、Hopfの1931年の原論文[[hopf1931]]において、同次座標（比として表した座標）と立体射影を用いて初めて構成されました。現代では四元数による定式化[[7shi-h]]も知られていますが、原論文の初等的な構成法を追うことで、複素2次元ベクトルがどのように2次元球面の点へ結びつくのか、その幾何学的な道筋が直接見通せます。

本記事では、3次元球面$S^3 \subset \mathbb R^4$の点$(x_1, x_2, x_3, x_4)$を出発点とします。まず実4次元の座標を複素数のペア$(x_1+ix_2, x_3+ix_4)$とみなし、その比$z \in \mathbb C$を取ることで同次座標の代表点を求めます。次に、複素平面$\mathbb C$を実3次元空間の$xy$平面に埋め込み、北極点からの逆立体射影によって2次元球面$S^2 \subset \mathbb R^3$上の点$(\xi_1, \xi_2, \xi_3)$を具体的に計算します。さらに、複素数のペア$(\alpha, \beta)$を用いた代数的表現に整理し、四元数による表式との符号関係を比較します。

高校数学レベルの複素数と初等幾何（直線の媒介変数表示と球面の交点計算）を前提とします。本記事では同次座標を用いたファイブレーションの初等的な座標計算に集中し、ファイバー（大円$S^1$）の幾何学的性質や結び目理論への発展は扱いません。

# 複素平面への射影

実4次元上の座標を複素2次元に写像します。

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

# 2次元球面への射影

$z$を実3次元空間の$xy$平面に埋め込みます。

$$
\mathbb{C} \ni z=u+iv \hookrightarrow (u,v,0) \in \mathbb{R}^3
$$

$(u,v,0)$を2次元球面$S^2$上に射影するため、$S^2$の北極点$(0,0,1)$と$(u,v,0)$を結ぶ直線を媒介変数$t$で表します。

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

# 複素数のペアによる表現

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
この結果は四元数による導出とほぼ同じですが、第2成分の符号が異なります。[[7shi-h]]

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

# まとめ

本記事では、3次元球面$S^3$上の点から複素数の比（同次座標）を取り出し、立体射影を組み合わせることで、Hopfの原論文の手法に基づきホップファイブレーションを計算しました。

実4次元座標による写像と、複素数のペアによる写像は次のようにまとめられます。

&&&
実座標表示：
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
複素ペア表示（$\alpha = x_1+ix_2,\ \beta = x_3+ix_4$）：
$$
\begin{pmatrix}\alpha \\ \beta\end{pmatrix}
\mapsto
\begin{pmatrix}\xi_1 + i\xi_2 \\ \xi_3\end{pmatrix}
=\begin{pmatrix}
  2 \alpha \beta^* \\
  \alpha \alpha^* - \beta \beta^*
 \end{pmatrix}
$$
&&&
