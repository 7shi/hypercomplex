微分のライプニッツ則に基づく関数の微分（微分形式）表記を用いることで、部分積分の各変形ステップと項の対応関係を明快に整理します。

# 概要

部分積分の公式 $\int fg'\,dx = fg - \int f'g\,dx$ は積の微分法の逆演算として広く知られていますが、どの関数を積分しどの関数を微分するか、なぜ符号が反転するのかの操作手順が機械的になりがちで、計算の見通しを失いやすい論点です。

本記事では、1階微分形式としての「関数の微分」$df = f'(x)dx$ を導入し、ライプニッツ則（積の微分法則）$d(fg) = f\,dg + g\,df$ と部分積分の関係を直接結びつけます。これにより、暗記に頼らず自然に変形を進められるようにします。続いて、1次式および2次式と三角関数の積の具体例を通して、関数の微分を用いた記法がどのように計算の段階を明快にするかを解説します。最後に、被積分関数の選択において次数を上げる方向を選んではならない注意点を確認します。

1変数関数の微分法、積の微分則、および不定積分の基礎知識を前提とします。本記事では初等関数の不定積分における部分積分法に焦点を当て、定積分の境界項処理や微分形式の外微分代数・ストークスの定理などの一般化は扱いません。

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

ここで、$df=\frac{df}{dx}dx=f'dx$, $dg=\frac{dg}{dx}dx=g'dx$という関数の微分の関係を用いれば、通常の形式が得られます。[[wiki-df]]

&&&fml 部分積分の公式（導関数）
$$
\int fg'\,dx = fg - \int f'g\,dx
$$
&&&

## 1次式と三角関数の積

基本的な例として、$f(x)=x$, $g(x)=\sin x$の場合を見てみましょう。

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

&&&rem 微分記号と関数の微分
$d$の中に入れる際、暗算でできるような簡単な積分を行います。
$$
\cos x\,dx=(\sin x)'dx=d(\sin x)
$$
この書き方は$\frac{d}{dx}$に$dx$を掛けて分母を取り払ったものだと解釈できます。
$$
d(\sin x)=\frac{d(\sin x)}{dx}\,dx=\left(\frac{d}{dx}\sin x\right)dx=(\sin x)'dx
$$
&&&

### 導関数表記との比較

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

## 2次式と三角関数の積

次に、$f(x)=x^2$, $g(x)=\cos x$という、やや複雑な例を見てみましょう。

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

## 被積分関数の選択と次数の変化

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

本記事では、関数の微分$dg = g'dx$を用いた表記で、ライプニッツ則から部分積分の公式を導き、1次式・2次式と三角関数の積の計算に適用しました。関数の微分を用いた表記を使うことで、部分積分の計算において結果に現れる項の形が予測しやすくなります。また、次数が上がる変形を避け、次数を下げる方向に計算を進めることが重要です。

&&& ライプニッツ則と部分積分
$$
d(fg) = f\,dg + g\,df, \quad
\int f\,dg = fg - \int g\,df
$$
&&&

&&& 導関数表記との対応
$$
\int fg'\,dx = \int f\,dg = fg - \int g\,df = fg - \int f'g\,dx
$$
&&&
