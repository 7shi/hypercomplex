変数間に依存関係がある場合（例: $y=x^2$）、「依存関係の代入のタイミングによって、結果が変わるのではないか？」という疑問が生じることがあります。

本記事では、代入と微分の順序を入れ替えても、同じ結果が得られることを確認します。

# 操作の順序

関数$f(x,y)$について、変数間に$y=y(x)$という依存関係がある場合の微分を考えます。以下の2つの方法があります。

&&& 先に依存関係を代入
依存関係$y=y(x)$を$f(x,y)$に代入して、一変数関数を得ます。

$$
f(x,y(x))
$$

これを通常の一変数関数として微分します。

$$
\frac{df}{dx} = \frac{\partial f}{\partial x} + \frac{\partial f}{\partial y}\frac{dy}{dx}
$$

右辺は多変数関数の微分と合成関数の微分の組み合わせです。実際の計算では先に代入しているため、必ずしもこの形を経由するわけではありません。
&&&

&&& 後で依存関係を代入
$f(x,y)$を多変数関数として全微分します。

$$
df = \frac{\partial f}{\partial x}dx + \frac{\partial f}{\partial y}dy
$$

$dy = \frac{dy}{dx}dx$を代入します。これは合成関数の微分を明示的に行うことに相当します。

$$
df = \frac{\partial f}{\partial x}dx + \frac{\partial f}{\partial y}\frac{dy}{dx}dx
$$

両辺を$dx$で割れば、先に依存関係を代入したのと同じ形になります。

$$
\frac{df}{dx}=\frac{\partial f}{\partial x} + \frac{\partial f}{\partial y}\frac{dy}{dx}
$$

偏微分は変数の依存関係を考慮せず、形式的に指定された変数だけを対象にしますが、全微分を通じて依存関係を組み込むことで、最終的な結果に正しく反映されます。
&&&

以上の結果を、いくつかの具体例で確認してみましょう。

## 例 1

$$
f(x,y) = xy,\quad y=x^2
$$

&&& 先に依存関係を代入
\begin{aligned}
f(x,x^2) &= xx^2 = x^3 \\
\frac{df}{dx} &= 3x^2
\end{aligned}
&&&

&&& 後で依存関係を代入
$$
df = \frac{\partial f}{\partial x}dx + \frac{\partial f}{\partial y}dy = y\,dx + x\,dy
$$

$y=x^2,\ dy = 2x\,dx$ より

$$
\begin{aligned}
df &= x^2dx + x(2x\,dx) = 3x^2dx \\
\frac{df}{dx} &= 3x^2
\end{aligned}
$$
&&&

&&&rem ライプニッツ則
全微分によって同じ項$xy$を別々に偏微分した結果、自動的にライプニッツ則が得られることに注目してください。

$$
d(xy) = \frac{\partial(xy)}{\partial x}dx + \frac{\partial(xy)}{\partial y}dy = y\,dx + x\,dy
$$
&&&

## 例 2

$$
f(x,y) = \sin(xy),\quad y = e^x
$$

&&& 先に依存関係を代入
\begin{aligned}
f(x,e^x) &= \sin(xe^x) \\
\frac{df}{dx} &= \cos(xe^x)(e^x + xe^x)
\end{aligned}
&&&

&&& 後で依存関係を代入
$$
df = \frac{\partial f}{\partial x}dx + \frac{\partial f}{\partial y}dy = y\cos(xy)\,dx + x\cos(xy)\,dy
$$

$y=e^x,\ dy = e^xdx$ より

$$
\begin{aligned}
df &= e^x\cos(xe^x)dx + x\cos(xe^x)e^xdx = \cos(xe^x)(e^x + xe^x)dx \\
\frac{df}{dx} &= \cos(xe^x)(e^x + xe^x)
\end{aligned}
$$
&&&

## 例 3

$$
f(x,y) = x^2 + y^2,\quad y = \sin x
$$

&&& 先に依存関係を代入
\begin{aligned}
f(x,\sin(x)) &= x^2 + \sin^2 x \\
\frac{df}{dx} &= 2x + 2\sin x\cos x
\end{aligned}
&&&

&&& 後で依存関係を代入
$$
df = \frac{\partial f}{\partial x}dx + \frac{\partial f}{\partial y}dy = 2x\,dx + 2y\,dy
$$

$y=\sin x,\ dy = \cos x\,dx$ より

$$
\begin{aligned}
df &= 2x\,dx + 2\sin x\cos x\,dx = (2x + 2\sin x\cos x)dx \\
\frac{df}{dx} &= 2x + 2\sin x\cos x
\end{aligned}
$$
&&&

# 補足：オイラー＝ラグランジュ方程式

これまでの議論は、オイラー＝ラグランジュ方程式における偏微分の扱いを理解する上でも役立ちます。オイラー＝ラグランジュ方程式は、ラグランジアン$L(q, \dot{q}, t)$を用いて、運動方程式を記述します。ここで、$q$は一般化座標、$\dot{q}$は一般化速度、そして$t$は時間です。[[7shi-eleq]]

&&&fml オイラー＝ラグランジュ方程式
$$
\frac{\partial L}{\partial q} - \frac{d}{dt}\left(\frac{\partial L}{\partial \dot{q}}\right) = 0
$$
&&&

この方程式で一見奇妙に思えるのは、$\dot{q}$ が$q$の時間微分であるにもかかわらず、$L$ を$q$と$\dot{q}$で独立に偏微分している点です。しかし、ここまでの議論で見てきたように、依存関係を考慮した微分の結果は、先に依存関係を代入して微分した場合と、後で代入した場合で一致します。

つまり、オイラー＝ラグランジュ方程式における偏微分は、形式的に$q$と$\dot{q}$を独立変数として扱っているように見えますが、最終的に$\dot{q} = \frac{dq}{dt}$の関係を代入することで、時間微分を考慮した正しい運動方程式が得られます。

このように、依存関係を持つ変数に対する微分の計算手順を理解しておくことは、オイラー＝ラグランジュ方程式の理解を深める上でも役立ちます。
