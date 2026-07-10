前回[[7shi-h]]、四元数の回転でホップファイブレーションを導出しました。その続きとして、具体的に成分を計算します。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

# パラメーターの調整

単位四元数による回転子$ω$と$q$を定義します。

&&&def 単位四元数による回転子
$$\begin{aligned}
ω&=(α_0-α_1\mathbf{k})+(β_0\mathbf{j}-β_1\mathbf{i})&&(α_0^2+α_1^2+β_0^2+β_1^2=1) \\
 &=(α_0-α_1\mathbf{k})+\mathbf{j}(β_0-β_1\mathbf{k}) \\
q&=u+v\mathbf{k}&&(u^2+v^2=1)
\end{aligned}$$
$q\mathbf kq^*=\mathbf k$より$(ωq)\mathbf k(ωq)^*=ω\mathbf kω^*$
&&&

$ωq$を計算します。

$$\begin{aligned}
ωq&=\{(α_0-α_1\mathbf{k})+\mathbf{j}(β_0-β_1\mathbf{k})\}(u+v\mathbf{k}) \\
&=(α_0-α_1\mathbf{k})(u+v\mathbf{k})+\mathbf{j}(β_0-β_1\mathbf{k})(u+v\mathbf{k})
\end{aligned}$$

$u,v$を調整すれば、$\mathbf k$を消すことができます。

$$\begin{aligned}
α_0-α_1\mathbf{k}&=γ(α_0'-α_1'\mathbf{k})&&(α_0'^2+α_1'^2=1) \\
u+v\mathbf{k}&=α_0'+α_1'\mathbf{k}
\end{aligned}$$
$$\begin{aligned}
ωq&=γ(α_0'-α_1'\mathbf{k})(α_0'+α_1'\mathbf{k})+\mathbf{j}(β_0-β_1\mathbf{k})(α_0'+α_1'\mathbf{k}) \\
  &=γ(α_0'^2+α_1'^2)+\mathbf{j}\{(β_0α_0'+β_1α_1')+(β_0α_1'-β_1α_0')\mathbf{k}\} \\
  &=γ+(β_0α_0'+β_1α_1')\mathbf{j}+(β_0α_1'-β_1α_0')\mathbf{i} \\
\end{aligned}$$

つまり$ω$と同じ回転を表す$ωq$のうち、$\mathbf k$の項がないものが存在します。

## 三角関数表示

成分計算だけではイメージが湧きにくいため、三角関数でパラメーター表示します。[[7shi-c]]

&&&def 回転子の三角関数表示
\begin{aligned}
ω&=\cosθ\,(\cos a-\sin a\,\mathbf{k})+\sinθ\,\mathbf{j}(\cos b-\sin b\,\mathbf{k}) \\
 &=\cosθ\,e^{-a\,\mathbf{k}}+\sinθ\,\mathbf{j}\,e^{-b\,\mathbf{k}} \\
q&=\cos c+\sin c\,\mathbf{k} \\
 &=e^{c\,\mathbf{k}}
\end{aligned}
&&&

$ωq$を計算します。

$$\begin{aligned}
ωq&=(\cosθ\,e^{-a\,\mathbf{k}}+\sinθ\,\mathbf{j}\,e^{-b\,\mathbf{k}})e^{c\,\mathbf{k}} \\
  &=\cosθ\,e^{(c-a)\,\mathbf{k}}+\sinθ\,\mathbf{j}\,e^{(c-b)\,\mathbf{k}}
\end{aligned}$$

$a=c$のとき

$$\begin{aligned}
ωq&=\cosθ+\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}} \\
  &=\cosθ+\sinθ\,\mathbf{j}\left\{\cos(a-b)-\sin(a-b)\,\mathbf{k}\right\} \\
  &=\cosθ+\sinθ\left\{\cos(a-b)\,\mathbf{j}-\sin(a-b)\,\mathbf{i}\right\}
\end{aligned}$$

よって、$q$が$e^{-a\,\mathbf{k}}$を打ち消すとき、$\mathbf k$の項が消えます。これを$ω'$とします。

&&&def k のない回転子
\begin{aligned}
ω'&=\cosθ+\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}} \\
\end{aligned}
&&&

## 回転

$ω'$で$\mathbf k$を回転させます。

$$\begin{aligned}
ω'\mathbf kω'^*
&=(\cosθ+\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}})\mathbf{k}(\cosθ+\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}})^* \\
&=(\cosθ+\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}})\mathbf{k}(\cosθ-\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}}) \\
&=(\cosθ+\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}})^2\,\mathbf{k} \\
&=(\cos^2θ+2\cosθ\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}}
  +\sin^2θ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}}\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}})\mathbf{k} \\
&=(\cos^2θ+2\cosθ\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}}
  +\sin^2θ\,\mathbf{j}^2\,e^{-(a-b)\,\mathbf{k}}e^{(a-b)\,\mathbf{k}})\mathbf{k} \\
&=2\cosθ\sinθ\,\mathbf{i}\left\{\cos(a-b)+\sin(a-b)\mathbf{k}\right\}
  +(\cos^2θ-\sin^2θ)\mathbf{k} \\
&=\sin2θ\left\{\cos(a-b)\mathbf{i}-\sin(a-b)\mathbf{j}\right\}+\cos2θ\,\mathbf{k} \\
&=\sin2θ\left\{\cos(b-a)\mathbf{i}+\sin(b-a)\mathbf{j}\right\}+\cos2θ\,\mathbf{k} \\
\end{aligned}$$

左右から回転子で挟むため、回転角が2倍の$2θ$になります。

&&&rem
一般的には回転子の側で角度を半分にして、回転角が$θ$となるように調整します。
$$ω''=\cos\fracθ2+\sin\fracθ2\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}}$$
ここでは敢えて、理由を説明する前に天下りで導入するのを避けました。
&&&

三角関数の引数を$b-a$とした理由は、次で説明します。

## 内積と外積

複素数$e^{i(b-a)}$を計算します。

$$\begin{aligned}
e^{i(b-a)}
&=e^{ib}e^{-ia} \\
&=(e^{ia})^*e^{ib} \\
&=(\cos a-i\sin a)(\cos b+i\sin b) \\
&=(\cos a\cos b+\sin a\sin b)+i(\cos a\sin b-\sin a\cos b)
\end{aligned}$$

これは加法定理ですが、様子を見るため実部を$x$座標、虚部を$y$座標に割り当てます。

$$\begin{alignedat}{2}
e^{ia}&=\cos a+i\sin a&&=a_x+ia_y \\
e^{ib}&=\cos b+i\sin b&&=b_x+ib_y
\end{alignedat}$$
$$e^{i(b-a)}=(a_xb_x+a_yb_y)+i(a_xb_y-a_yb_x)$$

位相差$b-a$を計算することで、実部は内積、虚部は（ウェッジ積の意味での）外積が得られると解釈できます。

&&&def 複素数による内積と外積の計算
$$a^*b=a\cdot b+i\norm{a\wedge b}$$
&&&

外積は演算の順序で符号が変わります。偏角を昇順で並べたものを順方向とします。偏角$a,b$を昇順とすれば、位相差は$b-a$として計算するのが自然です。

&&&ex 位相差
$20°, 50°$の順に並んでいれば、位相差は$50°-20°=30°$
&&&

# 状態ベクトルとパウリ行列

同じことを状態ベクトルとパウリ行列で計算します。

$$\begin{aligned}
ω&=\pmatrix{α\\β}
  =\pmatrix{α_0+iα_1\\β_0+iβ_1}
  =\pmatrix{\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}}
  =e^{ia}\pmatrix{\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}} \\
q&=e^{-ia} \\
ωq&=\pmatrix{\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}}=:ω'
\end{aligned}$$

$ωq$は括り出した係数の消去に相当することが、形式的な操作から浮き彫りになります。

パウリ行列でブロッホ球上の座標に変換します。

$$\begin{aligned}
ω'^\dagger σ_x ω'
&=\pmatrix{\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}}\pmatrix{0&1\\1&0}\pmatrix{\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}} \\
&=\pmatrix{\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}}\pmatrix{\sin\fracθ2\,e^{i(b-a)} \\ \cos\fracθ2} \\
&=\cos\fracθ2\sin\fracθ2\,e^{i(b-a)} + \sin\fracθ2\,e^{-i(b-a)}\cos\fracθ2 \\
&=\cos\fracθ2\sin\fracθ2\,(e^{i(b-a)}+e^{-i(b-a)}) \\
&=\cos\fracθ2\sin\fracθ2\left\{\cos(b-a)+i\sin(b-a)+\cos(b-a)-i\sin(b-a)\right\} \\
&=2\cos\fracθ2\sin\fracθ2\cos(b-a) \\
&=\sinθ\cos(b-a) \\

ω'^\dagger σ_y ω'
&=\pmatrix{\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}}\pmatrix{0&-i\\i&0}\pmatrix{\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}} \\
&=\pmatrix{\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}}\pmatrix{-i\sin\fracθ2\,e^{i(b-a)} \\ i\cos\fracθ2} \\
&=-i\cos\fracθ2\sin\fracθ2\,e^{i(b-a)} + i\sin\fracθ2\,e^{-i(b-a)}\cos\fracθ2 \\
&=-i\cos\fracθ2\sin\fracθ2\,(e^{i(b-a)}-e^{-i(b-a)}) \\
&=-i\cos\fracθ2\sin\fracθ2\left\{\cos(b-a)+i\sin(b-a)-\cos(b-a)+i\sin(b-a)\right\} \\
&=2\cos\fracθ2\sin\fracθ2\sin(b-a) \\
&=\sinθ\sin(b-a) \\

ω'^\dagger σ_z ω'
&=\pmatrix{\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}}\pmatrix{1&0\\0&-1}\pmatrix{\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}} \\
&=\pmatrix{\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}}\pmatrix{\cos\fracθ2 \\ -\sin\fracθ2\,e^{i(b-a)}} \\
&=\cos^2\fracθ2-\sin^2\fracθ2\,e^{-i(b-a)}e^{i(b-a)} \\
&=\cosθ
\end{aligned}$$

四元数は全成分を一度に計算するためごちゃごちゃしますが、こちらは成分ごとに計算するため冗長になります。どちらが簡単だとは一概には言えませんが、結果だけまとめるのにベクトルはシンプルです。

$$\begin{aligned}
ω'&=\pmatrix{\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}}
\mapsto\pmatrix{\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ}
\end{aligned}$$

$θ$が2倍になって、指数関数が表す$\cos$と$\sin$が別成分に分離されました。これを**ブロッホベクトル**と呼びます。

$ω$を計算しても同じ結果になります。

$$\begin{aligned}
ω&=\pmatrix{\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}}
\mapsto\pmatrix{\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ}
\end{aligned}$$

このように、ブロッホベクトルには$a,b$個々の位相ではなく、位相差$b-a$だけが現れます。量子情報では位相差だけが観測されるという重要な事実を反映しています。

&&&rem
$ωq$で見たように、状態ベクトル全体に位相を掛けてもブロッホベクトルには反映されません。
&&&

## まとめて計算

パウリ行列を使った方法でも、全成分をまとめて計算することはできます。

$A,B$を複素ベクトルとして、$A^\dagger B$は内積となります。エルミート共役を逆にした$AB^\dagger$は外積と呼ばれます。

&&&rem 外積[[7shi-p]]
既に出て来たウェッジ積や、ベクトル解析で一般的なベクトル積（クロス積）とは別の外積で、テンソル積（クロネッカー積）の特殊な場合です。ややこしいですが、英語では区別があります。

英語|日本語|備考
----|----|----
outer product<br>tensor direct product|外積<br>直積[[wikipedia-dp]][[dp]]|テンソル積（クロネッカー積）の特殊な場合<br>直積も他の用法と紛らわしくマイナー
exterior product<br>wedge product|外積、外部積<br>ウェッジ積、楔積|外積代数で使用、次元に依存しない<br>3次元ではベクトル積のホッジ双対
<br>vector product<br>cross product|外積<br>ベクトル積<br>クロス積|ベクトル解析で使用<br>次元に依存する（3次元と7次元）<br>四元数・八元数の積の虚部から得られる
&&&

（テンソル積の意味での）外積はベクトルから行列を作る計算で、自身との外積は複素共役との積を列挙したものとなります。

$$\begin{aligned}
ω&=\pmatrix{α\\β}=\pmatrix{\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}} \\
ωω^\dagger
&=\pmatrix{α\\β}\pmatrix{α^*&β^*} \\
&=\pmatrix{αα^*&αβ^*\\βα^*&ββ^*} \\
&=\pmatrix{
  \cos\fracθ2\,e^{ia}\cos\fracθ2\,e^{-ia} & \cos\fracθ2\,e^{ia}\sin\fracθ2\,e^{-ib} \\
  \sin\fracθ2\,e^{ib}\cos\fracθ2\,e^{-ia} & \sin\fracθ2\,e^{ib}\sin\fracθ2\,e^{-ib}} \\
&=\pmatrix{
  \cos^2\fracθ2 & \cos\fracθ2\sin\fracθ2\,e^{-i(b-a)} \\
  \cos\fracθ2\sin\fracθ2\,e^{i(b-a)} & \sin^2\fracθ2}
\end{aligned}$$

成分に見覚えのある形が現れています。半角公式[[half]]などを使って計算を進めます。

$$\begin{aligned}
ωω^\dagger
&=\frac12\pmatrix{
  1+\cosθ & \sinθ\left\{\cos(b-a)-i\sin(b-a)\right\} \\
  \sinθ\left\{\cos(b-a)+i\sin(b-a)\right\} & 1-\cosθ}
\end{aligned}$$

対角成分と非対角成分で、それぞれ符号が異なる項が現れています。これらを分離すれば、パウリ行列の線形結合が現れます。

$$\begin{aligned}
ωω^\dagger
&=\frac12\left\{
  \pmatrix{1&0 \\ 0&1} \\
  +\sinθ\cos(b-a)\pmatrix{0&1 \\ 1&0} \\
  +\sinθ\sin(b-a)\pmatrix{0&-i \\ i&0} \\
  +\cosθ\pmatrix{1&0 \\ 0&-1}\right\} \\
&=\frac12\left\{ \mathrm{I} + \sinθ\cos(b-a)σ_x + \sinθ\sin(b-a)σ_y + \cosθ\,σ_z \right\}
\end{aligned}$$

これは**密度行列**と呼ばれます。スケールを$2$倍してパウリ行列の係数を抜き出せば、ブロッホベクトルと一致します。

$$\begin{aligned}
ω&=\pmatrix{\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}}
\mapsto\pmatrix{\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ}
\end{aligned}$$

&&&rem 始点
四元数による回転の始点を$k$としたのは、外積から計算される成分と一致させるためです。
&&&

&&&rem 密度行列
ここでの密度行列は1つのブロッホベクトルを変換したもので、密度というよりベクトルの別表現に見えます。このように1つのブロッホベクトルしか含まない状態を**純粋状態**と呼びます。それに対して、複数のブロッホベクトルに対応する密度行列を重み付けして足し合わせることがあり、**混合状態**と呼びます。

純粋状態の密度行列に対応するブロッホベクトルはブロッホ球の表面を指しますが、混合状態ではブロッホ球の内部を指します。[[quantumuniverse]]
&&&

# まとめ

四元数での回転によるモデルを導入して、線形代数でその結果を利用する流れを示しました。量子情報への導入とすることを意識しています。

他にも射影平面による定式化があり、こちらの方がより一般的です。[[dim7-8]]
