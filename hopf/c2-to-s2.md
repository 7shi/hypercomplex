複素数のペア$(z_0, z_1)$を2次元球面上の座標に対応させる方法を、なるべく初等的に説明します。

オイラーの公式を前提とします。

&&&fml オイラーの公式
$$
e^{i\theta} = \cos\theta + i\sin\theta
$$
&&&

## 複素数のペアと球面

複素数のペア$(z_0, z_1)$があり、その大きさの和が 1 になるように正規化されているとします。

$$
|z_0|^2 + |z_1|^2 = 1
$$

これは4次元空間（複素数は実数2つ分なので$2+2=4$次元）の中にある**3次元球面**の方程式です。

複素数$z_0, z_1$を極形式（絶対値と偏角）で表します。$|z_0|^2 + |z_1|^2 = 1$なので、絶対値を$\cos\theta,\ \sin\theta$とおくことができます（$0 \le \theta \le \pi/2$）。

$$
\begin{pmatrix}z_0 \\ z_1\end{pmatrix}
= \begin{pmatrix}e^{i\alpha} \cos\theta \\ e^{i\beta} \sin\theta\end{pmatrix}
$$

&&&rem 任意性
$\cos\theta,\ \sin\theta$の割り当てには任意性がありますが、本記事ではパウリ行列によって構成される標準的なブロッホ球の定義に準拠しています。[[7shi-h]]
&&&

## グローバル位相の無視（自由度の減少）

複素数のペアは実数で4自由度を持つため、3次元空間内の**2次元球面**に射影するには自由度を1つ捨てる必要があります。

そのために、この複素数のペアに対して「全体を回転させても同じとみなす」というルールを導入します。複素数の掛け算における回転とは、絶対値が1の複素数$e^{i\omega}$を掛けることです。

$$
\begin{pmatrix} z_0 \\ z_1 \end{pmatrix}
\sim e^{i\omega} \begin{pmatrix} z_0 \\ z_1 \end{pmatrix}
= \begin{pmatrix} e^{i\omega}z_0 \\ e^{i\omega}z_1 \end{pmatrix}
$$

これは、2つの数の比$z_0 : z_1$だけに注目し、全体にかかる余分な係数を無視することに相当します。

&&&rem 同値関係
$A \sim B$は、特定の条件を課すことで同一視される関係（**同値関係**）を表します。
&&&

## 代表元の選び方（第1成分の絶対値化）

同値と見なされる無数の組み合わせの中から、代表となる1つを選び出して形を固定します。

第1成分$z_0 = e^{i\alpha} \cos\theta$の位相を打ち消すため、全体に$e^{-i\alpha}$を掛けます。

$$
\begin{pmatrix} z_0 \\ z_1 \end{pmatrix}
\sim e^{-i\alpha} \begin{pmatrix} z_0 \\ z_1 \end{pmatrix}
= e^{-i\alpha} \begin{pmatrix} e^{i\alpha} \cos\theta \\ e^{i\beta} \sin\theta \end{pmatrix}
= \begin{pmatrix} \cos\theta \\ e^{i(\beta-\alpha)} \sin\theta \end{pmatrix}
$$

$-1 \le \cos\theta \le 1$ですが、符号の違いはオイラーの等式$e^{i\pi}=-1$より同値と見なされます。

$$
-\begin{pmatrix} z_0 \\ z_1 \end{pmatrix}
= e^{i\pi} \begin{pmatrix} z_0 \\ z_1 \end{pmatrix}
\sim \begin{pmatrix} z_0 \\ z_1 \end{pmatrix}
$$

符号の違いを吸収するため、$\theta$の定義域により正の値域に制限します。これは第1成分の絶対値に相当します。

$$
\cos\theta = |z_0| \ge 0 \quad \left( 0 \le \theta \le \frac{\pi}{2} \right)
$$

ここで$\phi = \beta-\alpha$とおきます。標準的なブロッホ球の定義では、第1成分を北極と南極を貫く軸として扱うため、$Z$とおきます。第2成分は複素数なので、実部$X$と虚部$Y$に分けます。

$$
\begin{pmatrix} Z \\ X+iY \end{pmatrix}
= \begin{pmatrix} \cos\theta \\ e^{i\phi} \sin\theta \end{pmatrix}
$$

座標$(X, Y, Z) \in \mathbb R^3$は以下のように求まります。これは標準的な3次元の極座標表示です。

$$
\begin{pmatrix}X \\ Y \\ Z\end{pmatrix}
= \begin{pmatrix}
\text{Re}\left( e^{i\phi} \sin\theta \right) \\
\text{Im}\left( e^{i\phi} \sin\theta \right) \\
\cos\theta
\end{pmatrix}
= \begin{pmatrix}
\cos\phi \sin\theta \\
\sin\phi \sin\theta \\
\cos\theta
\end{pmatrix}
$$

&&&rem
$z_0=0$のときは位相$\alpha$が不定ですが、その場合は$Z=0$となり、後で扱う「赤道の縮退」のケースとなります。
&&&

## 3次元座標と半球

この$(X, Y, Z)$が2次元球面$S^2$上の一点を指すことを確認します。

$$
\begin{aligned}
X^2 + Y^2 + Z^2
&= \cos^2\phi \sin^2\theta + \sin^2\phi \sin^2\theta + \cos^2\theta \\
&= (\cos^2\phi + \sin^2\phi) \sin^2\theta + \cos^2\theta \\
&= \sin^2\theta + \cos^2\theta \\
&= 1
\end{aligned}
$$

ただし、$Z \ge 0$という制限により、この段階で表される範囲は、球面の北半球になります。

## 赤道の縮退（巾着絞り）

ここで、半球の縁である赤道に注目します。赤道は$Z=0$であり、元の第1成分も$z_0=0$です。

$$
\begin{pmatrix} 0 \\ z_1 \end{pmatrix}
$$

このとき、第2成分$z_1$は（正規化により大きさは1ですが）どのような位相$e^{i\phi}$を持っていても構いません。しかし、ルール「全体を回転させても同じ」より

$$
\begin{pmatrix} 0 \\ e^{i\phi} \end{pmatrix} = e^{i\phi} \begin{pmatrix} 0 \\ 1 \end{pmatrix} \sim \begin{pmatrix} 0 \\ 1 \end{pmatrix}
$$

となり、赤道上の点 $X^2+Y^2=1,\ Z=0$ は、位相が違ってもすべて同じ点とみなされます。

つまり、座標上では円（赤道）に見えているものが、同値関係により「たった1つの点」に潰れているのです。

## 閉じた球面の完成

この状況を幾何学的にイメージしてみましょう。

1. 北半球があります。
2. その縁である赤道を、ぎゅっと絞って1点にします。

半球の口を閉じれば、それはもはや半球ではなく、閉じた完全な球面になります。赤道だった部分は、閉じた後の球面の南極になります。

## 角度の倍加（仕上げ）

この「口を閉じた球面」が単位球面に一致するように座標を補正します。

元の座標では北極$\theta=0,\ Z=1$から赤道$\theta=\pi/2,\ Z=0$までしかありませんでした。

$$
\begin{alignedat}{2}
\left. \begin{pmatrix}\cos \theta \\ e^{i\phi}\sin \theta\end{pmatrix} \right|_{\theta=0}
&= \begin{pmatrix}1 \\ 0\end{pmatrix}
&&\mapsto \left. \begin{pmatrix}
  \cos\phi \sin \theta \\
  \sin\phi \sin \theta \\
  \cos \theta
\end{pmatrix} \right|_{\theta=0}
= \begin{pmatrix}0 \\ 0 \\ 1\end{pmatrix}
\\
\left. \begin{pmatrix}\cos \theta \\ e^{i\phi}\sin \theta\end{pmatrix} \right|_{\theta=\frac{\pi}{2}}
&= \begin{pmatrix}0 \\ e^{i\phi}\end{pmatrix}
&&\mapsto \left. \begin{pmatrix}
  \cos\phi \sin \theta \\
  \sin\phi \sin \theta \\
  \cos \theta
\end{pmatrix} \right|_{\theta=\frac{\pi}{2}}
= \begin{pmatrix}\cos\phi \\ \sin\phi \\ 0\end{pmatrix}
\end{alignedat}
$$

しかし、閉じてしまった赤道を南極として扱うため、新しい球面では北極から南極まで$(0 \sim \pi)$の広がりを持つべきです。したがって、元の角度$\theta\ (0 \sim \pi/2)$を2倍に引き伸ばして、新しい角度$2\theta\ (0 \sim \pi)$に対応させる必要があります。

$$
\begin{alignedat}{2}
\left. \begin{pmatrix}\cos \theta \\ e^{i\phi}\sin \theta\end{pmatrix} \right|_{\theta=0}
&= \begin{pmatrix}1 \\ 0\end{pmatrix}
&&\mapsto \left. \begin{pmatrix}
  \cos\phi \sin 2\theta \\
  \sin\phi \sin 2\theta \\
  \cos 2\theta
\end{pmatrix} \right|_{2\theta=0}
= \begin{pmatrix}0 \\ 0 \\ 1\end{pmatrix}
\\
\left. \begin{pmatrix}\cos \theta \\ e^{i\phi}\sin \theta\end{pmatrix} \right|_{\theta=\frac{\pi}{2}}
&= \begin{pmatrix}0 \\ e^{i\phi}\end{pmatrix}
\sim \begin{pmatrix}0 \\ 1\end{pmatrix}
&&\mapsto \left. \begin{pmatrix}
  \cos\phi \sin 2\theta \\
  \sin\phi \sin 2\theta \\
  \cos 2\theta
\end{pmatrix} \right|_{2\theta=\pi}
= \begin{pmatrix}0 \\ 0 \\ -1\end{pmatrix}
\end{alignedat}
$$

以上をまとめると、$\theta \to 2\theta$によって補正した3次元座標$(x,y,z) \in \mathbb R^3$は、以下のようになります。

$$
\begin{pmatrix}X \\ Y \\ Z\end{pmatrix}
=\begin{pmatrix}
    \cos\phi \sin \theta \\
    \sin\phi \sin \theta \\
    \cos \theta
  \end{pmatrix}
\xrightarrow{\theta \to 2\theta}
  \begin{pmatrix}
    \cos\phi \sin 2\theta \\
    \sin\phi \sin 2\theta \\
    \cos 2\theta
  \end{pmatrix}
=: \begin{pmatrix}x \\ y \\ z\end{pmatrix}
$$

この補正によって、射影された球面は単位球面に一致します。

このように、複素数のペアを単位球面に射影しようとすると、必然的に俯角$\theta$が2倍になります。

&&&rem
このことをあらかじめ織り込んで、複素ベクトルの段階で角度を半分にしておくことが一般的です。
$$
\begin{pmatrix}z_0 \\ z_1\end{pmatrix}
= \begin{pmatrix}e^{i\alpha} \cos\frac{\theta}{2} \\ e^{i\beta} \sin\frac{\theta}{2}\end{pmatrix}
\sim \begin{pmatrix}\cos\frac{\theta}{2} \\ e^{i\phi} \sin\frac{\theta}{2}\end{pmatrix}
$$
&&&

## 複素数のまま計算する

最後に、この結果$(x,y,z)$を、角度$\theta, \phi$を使わずに元の複素数$z_0, z_1$だけで表してみましょう。

まず、$z$成分（高さ）について、2倍角の公式を用いて計算します。

$$
z = \cos 2\theta = \cos^2\theta - \sin^2\theta
$$

ここで、$|z_0| = \cos\theta,\ |z_1| = \sin\theta$より、

$$
z = |z_0|^2 - |z_1|^2
$$

となります。

次に、$x, y$成分（水平方向）について計算します。

$$
x + iy = (\cos\phi + i\sin\phi) \sin 2\theta = e^{i\phi} \sin 2\theta
$$

ここで、$\phi = \beta - \alpha$ と2倍角の公式より、

$$
\begin{aligned}
x + iy &= e^{i(\beta-\alpha)} \cdot 2\sin\theta\cos\theta \\
&= 2 (e^{-i\alpha}\cos\theta) (e^{i\beta}\sin\theta) \\
&= 2 z_0^* z_1
\end{aligned}
$$

となります。（$z_0^*$は$z_0$の複素共役）

まとめると、正規化された複素数のペア$(z_0, z_1)$から2次元球面上の点$(x,y,z)$への変換式は以下のようになります。

&&&fml ホップ写像
$$
\begin{pmatrix}x \\ y \\ z\end{pmatrix}
= \begin{pmatrix}
    2\text{Re}(z_0^* z_1) \\
    2\text{Im}(z_0^* z_1) \\
    |z_0|^2 - |z_1|^2
  \end{pmatrix}
= \begin{pmatrix}
    z_0^* z_1 + z_1^* z_0 \\
    -i(z_0^* z_1 - z_1^* z_0) \\
    z_0^* z_0 - z_1^* z_1
  \end{pmatrix}
$$
&&&

「グローバル位相の無視による代表元の固定」と「閉じた球面の形成」という幾何学的な要請から、この数式が自然に導かれました。

&&&rem
$z_0, z_1$の順序や符号は文献により多少異なりますが、本質的な構造は同じです。

また、$z_0^* z_1$はクリフォード代数の幾何積に相当する操作で、内積と外積（ウェッジ積）を同時に計算することに相当するため、結果には位相差のみが反映されます。
&&&
