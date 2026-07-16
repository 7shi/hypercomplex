部分積分は積分計算の重要なテクニックですが、関数の微分を用いた表記を利用すれば、段階を踏みやすくなります。具体的な例を交えながら説明します。

# 基本的な考え方

まず、二つの関数$f=f(x)$, $g=g(x)$について考えます。ライプニッツ則から始めることで、部分積分の本質が見えてきます。

$$
\begin{aligned}
d(fg) &= (df)g+f(dg) \\
f\,dg &= d(fg)-g\,df
\end{aligned}
$$

この式の両辺を積分することで、部分積分の公式が導かれます。

&&&fml 部分積分の公式（関数の微分）
$$
\int f\,dg = fg - \int g\,df
$$
&&&

ここで、$df=\frac{df}{dx}dx=f'dx$,$dg=\frac{dg}{dx}dx=g'dx$という関数の微分[[wiki-df]]の関係を用いれば、通常の形式が得られます。

&&&fml 部分積分の公式（導関数）
$$
\int fg'\,dx = fg - \int f'g\,dx
$$
&&&

## 例 1

基本的な例として、$f(x)=x$,$g(x)=\sin x$の場合を見てみましょう。

$$
\begin{aligned}
\int x\,d(\sin x) &= x\sin x - \int \sin x\,dx \\
\int x\cos x\,dx &= x\sin x + \cos x + C
\end{aligned}
$$

関数の微分を用いたおかげで、右辺に$\sin x$を含むことが自然に理解できます。

通常の計算では、2行目の左辺から始めることになります。

$$
\begin{aligned}
\int x\cos x\,dx
&= \int x\,d(\sin x) \\
&= x\sin x - \int \sin x\,dx \\
&= x\sin x + \cos x + C
\end{aligned}
$$

&&&rem
$d$の中に入れる際、暗算でできるような簡単な積分を行います。
$$
\cos x\,dx=(\sin x)'dx=d(\sin x)
$$
この書き方は$\frac{d}{dx}$に$dx$を掛けて分母を取り払ったものだと解釈できます。
$$
d(\sin x)=\frac{d(\sin x)}{dx}\,dx=\left(\frac{d}{dx}\sin x\right)dx=(\sin x)'dx
$$
&&&

### 補足

関数の微分を使わない場合、次のような書き方が考えられます。

$$
\begin{aligned}
\int x\cos x\,dx
&= \int x(\sin x)'dx \\
&= x\sin x - \int (x)'\sin x\,dx \\
&= x\sin x - \int \sin x\,dx \\
&= x\sin x + \cos x + C
\end{aligned}
$$

## 例 2

次に、$f(x)=x^2$,$g(x)=\cos x$という、やや複雑な例を見てみましょう。

$$
\begin{aligned}
\int x^2\,d(\cos x) &= x^2\cos x - \int \cos x\,d(x^2) \\
\int x^2(-\sin x)\,dx &= x^2\cos x - \int 2x\cos x\,dx \\
\int x^2\sin x\,dx &= -x^2\cos x + \int 2x\cos x\,dx \\
&= -x^2\cos x + \int 2x\,d(\sin x) \\
&= -x^2\cos x + 2x\sin x - \int 2\sin x\,dx \\
&= -x^2\cos x + 2x\sin x + 2\cos x + C
\end{aligned}
$$

通常の計算では、3行目の左辺から始めることになります。

$$
\begin{aligned}
\int x^2\sin x\,dx
&= \int x^2\,d(-\cos x) \\
&= x^2(-\cos x) - \int (-\cos x)\,d(x^2) \\
&= -x^2\cos x + \int 2x\cos x\,dx \\
&= -x^2\cos x + \int 2x\,d(\sin x) \\
&= -x^2\cos x + 2x\sin x - \int 2\sin x\,dx \\
&= -x^2\cos x + 2x\sin x + 2\cos x + C
\end{aligned}
$$

この例からわかるように、関数の微分を用いることで、各ステップでどのような変形が行われているのかが明確になります。特に、$\cos x$と$\sin x$の関係、そして$x^2$の微分が計算の中でどのように現れるかが見通しやすくなっています。

## 注意点

関数の微分を使う際、次数が上がるような変形は避けるべきです。

$$
\begin{aligned}
\int x^2\sin x\,dx
&= \int \sin x\,d\left(\frac{x^3}{3}\right) \\
&= \sin x\frac{x^3}{3} - \int \frac{x^3}{3}\,d(\sin x) \\
&= \sin x\frac{x^3}{3} - \frac{1}{3}\int x^3\cos x\,dx
\end{aligned}
$$

この変形では$x$の次数が上がってしまい、積分を外すことができません。次数を下げる方向に計算を進めることが重要です。

# まとめ

関数の微分を用いた表記を使うことで、部分積分の計算において結果に現れる項の形が予測しやすくなります。
