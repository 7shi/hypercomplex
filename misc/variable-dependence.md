変数間に依存関係がある多変数関数において、依存関係の代入と微分の順序を交換しても同一の結果が得られる構造を全微分を通して確認し、オイラー＝ラグランジュ方程式における形式的偏微分の正当性を整理します。

# 概要

多変数関数において$y=y(x)$のように変数間に従属関係があるとき、「先に代入して1変数関数として微分する」計算と、「多変数関数として偏微分した後に依存関係（微分形式）を代入する」計算で結果が異なるのではないかという疑問が生じることがあります。また、変分法や解析力学のオイラー＝ラグランジュ方程式において、位置$q$とその時間微分$\dot{q}$を独立な変数として偏微分することに対して直観的な違和感を抱く読者も少なくありません。

本記事では、合成関数の微分法則（全微分）を用いることで、代入と微分の順序が数学的に可換であり、どちらの手順を踏んでも同一の結果が得られることを一般論として示します。次に、積の微分、合成関数、和の自乗という3つの具体例を通じて、全微分がライプニッツ則や連鎖律を自然に内包していることを計算で確かめます。最後に、オイラー＝ラグランジュ方程式における$q$と$\dot{q}$の偏微分がこの可換性の直接的な応用であることを解説します。

1変数および多変数関数の微分法、合成関数の微分則、全微分の基礎知識を前提とします。本記事では2変数関数と1変数の依存関係における微分計算に焦点を当て、陰関数定理の厳密な証明や多次元多様体上の外微分代数は扱いません。

# 操作の順序

関数$f(x,y)$について、変数間に$y=y(x)$という依存関係がある場合の微分を考えます。以下の2つの方法があります。

&&& 先に依存関係を代入する手順
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

&&& 全微分後に依存関係を代入する手順
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

## 積の微分とライプニッツ則

$$
f(x,y) = xy,\quad y=x^2
$$

&&& 先に代入する場合
$$
\begin{aligned}
f(x,x^2) &= xx^2 = x^3 \\
\frac{df}{dx} &= 3x^2
\end{aligned}
$$
&&&

&&& 全微分後に代入する場合
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

## 指数関数と合成関数の微分

$$
f(x,y) = \sin(xy),\quad y = e^x
$$

&&& 先に代入する場合
$$
\begin{aligned}
f(x,e^x) &= \sin(xe^x) \\
\frac{df}{dx} &= \cos(xe^x)(e^x + xe^x)
\end{aligned}
$$
&&&

&&& 全微分後に代入する場合
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

## 和の自乗と三角関数の微分

$$
f(x,y) = x^2 + y^2,\quad y = \sin x
$$

&&& 先に代入する場合
$$
\begin{aligned}
f(x,\sin(x)) &= x^2 + \sin^2 x \\
\frac{df}{dx} &= 2x + 2\sin x\cos x
\end{aligned}
$$
&&&

&&& 全微分後に代入する場合
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

# オイラー＝ラグランジュ方程式との関連

これまでの議論は、オイラー＝ラグランジュ方程式における偏微分の扱いを理解する上でも役立ちます。オイラー＝ラグランジュ方程式は、ラグランジアン$L(q, \dot{q}, t)$を用いて、運動方程式を記述します。ここで、$q$は一般化座標、$\dot{q}$は一般化速度、そして$t$は時間です。[[7shi-eleq]]

&&&fml オイラー＝ラグランジュ方程式
$$
\frac{\partial L}{\partial q} - \frac{d}{dt}\left(\frac{\partial L}{\partial \dot{q}}\right) = 0
$$
&&&

この方程式で一見奇妙に思えるのは、$\dot{q}$ が$q$の時間微分であるにもかかわらず、$L$ を$q$と$\dot{q}$で独立に偏微分している点です。しかし、ここまでの議論で見てきたように、依存関係を考慮した微分の結果は、先に依存関係を代入して微分した場合と、後で代入した場合で一致します。

つまり、オイラー＝ラグランジュ方程式における偏微分は、形式的に$q$と$\dot{q}$を独立変数として扱っているように見えますが、最終的に$\dot{q} = \frac{dq}{dt}$の関係を代入することで、時間微分を考慮した正しい運動方程式が得られます。

このように、依存関係を持つ変数に対する微分の計算手順を理解しておくことは、オイラー＝ラグランジュ方程式の理解を深める上でも役立ちます。

# まとめ

本記事では、変数間に依存関係がある関数の微分において、代入と微分の順序を入れ替えても全微分を通じて完全に同一の結果が得られることを確認しました。

&&& 微分と代入の可換性
多変数関数$f(x,y)$と従属関係$y=y(x)$に対し、先に代入して1変数関数として微分する操作と、全微分$df=\frac{\partial f}{\partial x}dx+\frac{\partial f}{\partial y}dy$を求めてから$dy=\frac{dy}{dx}dx$を代入する操作は、合成関数の微分公式を通じて一致します。
$$
\frac{df}{dx} = \frac{\partial f}{\partial x} + \frac{\partial f}{\partial y}\frac{dy}{dx}
$$
&&&

&&& オイラー＝ラグランジュ方程式における形式的偏微分
オイラー＝ラグランジュ方程式においてラグランジアン$L(q,\dot{q},t)$を$q$と$\dot{q}$で独立に偏微分できるのは、この可換性に基づき、形式的に独立変数として偏微分した後に$\dot{q}=\frac{dq}{dt}$という関係を結びつける操作が正当化されているためです。
$$
\frac{\partial L}{\partial q} - \frac{d}{dt}\left(\frac{\partial L}{\partial \dot{q}}\right) = 0
$$
&&&
