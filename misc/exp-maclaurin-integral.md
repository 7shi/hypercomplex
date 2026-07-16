ピカールの逐次近似法の説明で指数関数を求める例がよく使われます[[rn-picard]][[wiki-picard]]。普通に計算すれば再帰的な積分からマクローリン展開と同じ結果が現れます。

# マクローリン展開

指数関数$e^x$をマクローリン展開します。[[masuo-maclaurin]]

$$
e^x = 1 + x + \frac{x^2}{2} + \frac{x^3}{3!} + \frac{x^4}{4!} + \dots = \sum_{n=0}^{\infty}\frac{x^n}{n!}
$$

これと同じ形が積分で得られることを確認します。

# 定積分

$e^x$を$0$から$x$まで積分します。

$$
\int_0^x e^{x'} dx' = [e^{x'}]_0^x = e^x - 1
$$

これを整理します。

$$
e^x = 1 + \int_0^x e^{x'} dx'
$$

右辺の積分内の$e^{x'}$に、変数を調整して右辺全体を再帰的に代入します。

$$
\begin{aligned}
e^x
&= 1 + \int_0^x \left(1 + \int_0^{x'} e^{x''} dx''\right) dx' \\
&= 1 + \int_0^x 1\,dx' + \int_0^x \int_0^{x'} e^{x''} dx'' dx'
\end{aligned}
$$

同様に続けます。

$$
\begin{aligned}
e^x
&= 1 + \int_0^x 1\,dx' + \int_0^x \int_0^{x'} \left(1 + \int_0^{x''} e^{x'''} dx'''\right) dx'' dx' \\
&= 1 + \int_0^x 1\,dx' + \int_0^x \int_0^{x'} 1\,dx'' dx' + \int_0^x \int_0^{x'} \int_0^{x''} e^{x'''} dx''' dx'' dx'
\end{aligned}
$$

これを繰り返せば、一般には次のようになります。

$$
e^x = 1 + \int_0^x 1\,dx' + \int_0^x \int_0^{x'} 1\,dx'' dx' + \int_0^x \int_0^{x'} \int_0^{x''} 1\,dx''' dx'' dx' + \cdots
$$

$1$の多重積分を計算します。

$$
e^x = 1 + x + \frac{x^2}{2} + \frac{x^3}{3!} + \frac{x^4}{4!} + \dots
$$

これはマクローリン展開に一致します。

積分からマクローリン展開と同じ結果が現れるのは興味深いですが、計算自体はピカールの逐次近似法を使った方が簡単です。
