前回[[7shi-h]]、四元数の回転でホップファイブレーションを導出しました。その続きとして、具体的に成分を計算します。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

改訂履歴：

- 2026.07.14 回転子の定義を調整、回転と内積・外積を削除、基準依存性を追記

# パラメーターの調整

単位四元数による回転子$ω$と$q$を定義します。

&&&def 単位四元数による回転子
$$\begin{aligned}
ω&=(α_0\mathbf{k}+α_1)+\mathbf{j}(β_0\mathbf{k}+β_1)&&(α_0^2+α_1^2+β_0^2+β_1^2=1) \\
q&=u+v\mathbf{k}&&(u^2+v^2=1)
\end{aligned}$$
&&&

$q\mathbf kq^*=\mathbf k$より、$q$が$\mathbf k$を固定するファイバーを構成するため、

$$
(ωq)\mathbf k(ωq)^*=ω\mathbf kω^*
$$

となります。この自由度を使って、$ω$の$\mathbf k$の項を消去した標準形を作ります。

$ωq$を計算します。

$$\begin{aligned}
ωq&=\{(α_0\mathbf{k}+α_1)+\mathbf{j}(β_0\mathbf{k}+β_1)\}(u+v\mathbf{k}) \\
&=(α_0\mathbf{k}+α_1)(u+v\mathbf{k})+\mathbf{j}(β_0\mathbf{k}+β_1)(u+v\mathbf{k})
\end{aligned}$$

$u,v$を調整すれば、$\mathbf k$を消すことができます。

$$\begin{aligned}
α_0\mathbf{k}+α_1&=γ(α_0'\mathbf{k}+α_1')&&(α_0'^2+α_1'^2=1) \\
u+v\mathbf{k}&=-α_0'\mathbf{k}+α_1'
\end{aligned}$$
$$\begin{aligned}
ωq&=γ(α_0'\mathbf{k}+α_1')(-α_0'\mathbf{k}+α_1')+\mathbf{j}(β_0\mathbf{k}+β_1)(-α_0'\mathbf{k}+α_1') \\
  &=γ(α_0'^2+α_1'^2)+\mathbf{j}\{(β_0α_0'+β_1α_1')+(β_0α_1'-β_1α_0')\mathbf{k}\} \\
  &=γ+(β_0α_0'+β_1α_1')\mathbf{j}+(β_0α_1'-β_1α_0')\mathbf{i} \\
\end{aligned}$$

つまり$ω$と同じ回転を表す$ωq$のうち、$\mathbf k$の項がないものが存在します。

## 三角関数表示

成分計算だけではイメージが湧きにくいため、三角関数でパラメーター表示します。[[7shi-c]]

&&&def 回転子の三角関数表示
$$
\begin{aligned}
ω&=\cosθ\,(\sin a\,\mathbf{k}+\cos a)+\sinθ\,\mathbf{j}(\sin b\,\mathbf{k}+\cos b) \\
 &=\cosθ\,e^{a\,\mathbf{k}}+\sinθ\,\mathbf{j}\,e^{b\,\mathbf{k}} \\
q&=\cos c+\sin c\,\mathbf{k} \\
 &=e^{c\,\mathbf{k}}
\end{aligned}
$$
&&&

$ωq$を計算します。

$$
\begin{aligned}
ωq&=(\cosθ\,e^{a\,\mathbf{k}}+\sinθ\,\mathbf{j}\,e^{b\,\mathbf{k}})e^{c\,\mathbf{k}} \\
  &=\cosθ\,e^{(a+c)\,\mathbf{k}}+\sinθ\,\mathbf{j}\,e^{(b+c)\,\mathbf{k}}
\end{aligned}
$$

$c=-a$のとき

$$\begin{aligned}
ωq&=\cosθ+\sinθ\,\mathbf{j}\,e^{(b-a)\,\mathbf{k}} \\
  &=\cosθ+\sinθ\,\mathbf{j}\left\{\cos(b-a)+\sin(b-a)\,\mathbf{k}\right\} \\
  &=\cosθ+\sinθ\left\{\cos(b-a)\,\mathbf{j}+\sin(b-a)\,\mathbf{i}\right\}
\end{aligned}$$

よって、$q$が$e^{a\,\mathbf{k}}$を打ち消すとき、$\mathbf k$の項が消えます。これを$ω'$とします。

&&&def k のない回転子
$$
\begin{aligned}
ω'&=\cosθ+\sinθ\,\mathbf{j}\,e^{(b-a)\,\mathbf{k}} \\
\end{aligned}
$$
&&&

# 状態ベクトルとパウリ行列

同じことを状態ベクトルとパウリ行列で計算します。回転子を表す2成分の複素ベクトルは**スピノル**とも呼ばれます。

$$\begin{aligned}
ω&=\begin{pmatrix}α\\β\end{pmatrix}
  =\begin{pmatrix}α_0+iα_1\\β_0+iβ_1\end{pmatrix}
  =\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix}
  =e^{ia}\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix} \\
q&=e^{-ia} \\
ωq&=\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix}=:ω'
\end{aligned}$$

$ωq$は括り出した係数の消去に相当することが、形式的な操作から浮き彫りになります。

パウリ行列でブロッホ球上の座標に変換します。

$$\begin{aligned}
ω'^\dagger \sigma_x ω'
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}0&1\\1&0\end{pmatrix}\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix} \\
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}\sin\fracθ2\,e^{i(b-a)} \\ \cos\fracθ2\end{pmatrix} \\
&=\cos\fracθ2\sin\fracθ2\,e^{i(b-a)} + \sin\fracθ2\,e^{-i(b-a)}\cos\fracθ2 \\
&=\cos\fracθ2\sin\fracθ2\,(e^{i(b-a)}+e^{-i(b-a)}) \\
&=\cos\fracθ2\sin\fracθ2\left\{\cos(b-a)+i\sin(b-a)+\cos(b-a)-i\sin(b-a)\right\} \\
&=2\cos\fracθ2\sin\fracθ2\cos(b-a) \\
&=\sinθ\cos(b-a) \\

ω'^\dagger \sigma_y ω'
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}0&-i\\i&0\end{pmatrix}\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix} \\
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}-i\sin\fracθ2\,e^{i(b-a)} \\ i\cos\fracθ2\end{pmatrix} \\
&=-i\cos\fracθ2\sin\fracθ2\,e^{i(b-a)} + i\sin\fracθ2\,e^{-i(b-a)}\cos\fracθ2 \\
&=-i\cos\fracθ2\sin\fracθ2\,(e^{i(b-a)}-e^{-i(b-a)}) \\
&=-i\cos\fracθ2\sin\fracθ2\left\{\cos(b-a)+i\sin(b-a)-\cos(b-a)+i\sin(b-a)\right\} \\
&=2\cos\fracθ2\sin\fracθ2\sin(b-a) \\
&=\sinθ\sin(b-a) \\

ω'^\dagger \sigma_z ω'
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}1&0\\0&-1\end{pmatrix}\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix} \\
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}\cos\fracθ2 \\ -\sin\fracθ2\,e^{i(b-a)}\end{pmatrix} \\
&=\cos^2\fracθ2-\sin^2\fracθ2\,e^{-i(b-a)}e^{i(b-a)} \\
&=\cosθ
\end{aligned}$$

四元数は全成分を一度に計算するためごちゃごちゃしますが、こちらは成分ごとに計算するため冗長になります。どちらが簡単だとは一概には言えませんが、結果だけまとめるのにベクトルはシンプルです。

$$\begin{aligned}
ω'&=\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix}
\mapsto\begin{pmatrix}\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ\end{pmatrix}
\end{aligned}$$

$θ$が2倍になって、指数関数が表す$\cos$と$\sin$が別成分に分離されました。これを**ブロッホベクトル**と呼びます。

$ω$を計算しても同じ結果になります。

$$\begin{aligned}
ω&=\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix}
\mapsto\begin{pmatrix}\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ\end{pmatrix}
\end{aligned}$$

このように、ブロッホベクトルには$a,b$個々の位相ではなく、位相差$b-a$だけが現れます。量子情報では位相差だけが観測されるという重要な事実を反映しています。

## 全体位相とファイバー

状態ベクトル全体に位相$e^{iφ}$を掛けると、

$$
e^{iφ}ω=\begin{pmatrix}\cos\fracθ2\,e^{i(a+φ)} \\ \sin\fracθ2\,e^{i(b+φ)}\end{pmatrix}
$$

となり、$a,b$がともに$φ$だけ動きます。位相差$b-a$は変わらないため、ブロッホベクトルも変わりません。

これは$ωq$におけるファイバー自由度$q$に対応します。つまり状態ベクトルに全体位相を掛ける操作は、ホップファイブレーションのファイバーに沿って点を動かす操作に対応します。ブロッホベクトルがこの操作で変わらないのは、ファイバー上の点$ωq$（$q$を動かして得られる各点）がすべて同じ$S^2$上の点に写ることの表れです。

&&&rem
「ファイバー」という用語は、ホップファイブレーションの構造として抽象的に円周$S^1$を指す場合と、$S^3$内に埋め込まれた具体的な部分集合$\{ωq:q\in S^1\}$を指す場合の両方に使われます。$q$は$S^1$上を動くパラメーターで、$ωq$は$S^3$上の点です。
&&&

## まとめて計算

パウリ行列を使った方法でも、全成分をまとめて計算することはできます。

$A,B$を複素ベクトルとして、$A^\dagger B$は内積となります。エルミート共役を逆にした$AB^\dagger$は外積と呼ばれます。

&&&rem 外積[[7shi-p]]
ベクトル解析で一般的なベクトル積（クロス積）や、ウェッジ積とは別の外積で、テンソル積（クロネッカー積）の特殊な場合です。ややこしいですが、英語では区別があります。

英語|日本語|備考
----|----|----
outer product<br>tensor direct product|外積<br>直積[[wikipedia-dp]][[dp]]|テンソル積（クロネッカー積）の特殊な場合<br>直積も他の用法と紛らわしくマイナー
exterior product<br>wedge product|外積、外部積<br>ウェッジ積、楔積|外積代数で使用、次元に依存しない<br>3次元ではベクトル積のホッジ双対
<br>vector product<br>cross product|外積<br>ベクトル積<br>クロス積|ベクトル解析で使用<br>次元に依存する（3次元と7次元）<br>四元数・八元数の積の虚部から得られる
&&&

（テンソル積の意味での）外積はベクトルから行列を作る計算で、自身との外積は複素共役との積を列挙したものとなります。

$$\begin{aligned}
ω&=\begin{pmatrix}α\\β\end{pmatrix}=\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix} \\
ωω^\dagger
&=\begin{pmatrix}α\\β\end{pmatrix}\begin{pmatrix}α^*&β^*\end{pmatrix} \\
&=\begin{pmatrix}αα^*&αβ^*\\βα^*&ββ^*\end{pmatrix} \\
&=\begin{pmatrix}
  \cos\fracθ2\,e^{ia}\cos\fracθ2\,e^{-ia} & \cos\fracθ2\,e^{ia}\sin\fracθ2\,e^{-ib} \\
  \sin\fracθ2\,e^{ib}\cos\fracθ2\,e^{-ia} & \sin\fracθ2\,e^{ib}\sin\fracθ2\,e^{-ib}\end{pmatrix} \\
&=\begin{pmatrix}
  \cos^2\fracθ2 & \cos\fracθ2\sin\fracθ2\,e^{-i(b-a)} \\
  \cos\fracθ2\sin\fracθ2\,e^{i(b-a)} & \sin^2\fracθ2\end{pmatrix}
\end{aligned}$$

成分に見覚えのある形が現れています。半角公式[[half]]などを使って計算を進めます。

$$\begin{aligned}
ωω^\dagger
&=\frac12\begin{pmatrix}
  1+\cosθ & \sinθ\left\{\cos(b-a)-i\sin(b-a)\right\} \\
  \sinθ\left\{\cos(b-a)+i\sin(b-a)\right\} & 1-\cosθ\end{pmatrix}
\end{aligned}$$

対角成分と非対角成分で、それぞれ符号が異なる項が現れています。これらを分離すれば、パウリ行列の線形結合が現れます。

$$\begin{aligned}
ωω^\dagger
&=\frac12\left\{
  \begin{pmatrix}1&0 \\ 0&1\end{pmatrix} \\
  +\sinθ\cos(b-a)\begin{pmatrix}0&1 \\ 1&0\end{pmatrix} \\
  +\sinθ\sin(b-a)\begin{pmatrix}0&-i \\ i&0\end{pmatrix} \\
  +\cosθ\begin{pmatrix}1&0 \\ 0&-1\end{pmatrix}\right\} \\
&=\frac12\left\{ \mathrm{I} + \sinθ\cos(b-a)\sigma_x + \sinθ\sin(b-a)\sigma_y + \cosθ\,\sigma_z \right\}
\end{aligned}$$

これは**密度行列**と呼ばれます。スケールを$2$倍してパウリ行列の係数を抜き出せば、ブロッホベクトルと一致します。

$$\begin{aligned}
ω&=\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix}
\mapsto\begin{pmatrix}\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ\end{pmatrix}
\end{aligned}$$

&&&rem 始点
四元数による回転の始点を$\mathbf k$としたのは、外積から計算される成分と一致させるためです。
&&&

&&&rem 密度行列
ここでの密度行列は1つのブロッホベクトルを変換したもので、密度というよりベクトルの別表現に見えます。このように1つのブロッホベクトルしか含まない状態を**純粋状態**と呼びます。それに対して、複数のブロッホベクトルに対応する密度行列を重み付けして足し合わせることがあり、**混合状態**と呼びます。

純粋状態の密度行列に対応するブロッホベクトルはブロッホ球の表面を指しますが、混合状態ではブロッホ球の内部を指します。[[quantumuniverse]]
&&&

# 全体位相の基準依存性

ホップファイブレーションが定める束$S^1\hookrightarrow S^3\xrightarrow{\text{Hopf}}S^2$は、単純な直積$S^2\times S^1$ではありません。任意の「$S^3$上の点」を「$S^2$上の点」と「$S^1$上の位相」の組に一意に分解できないためです。具体的に確認します。

先ほど見た通り、状態ベクトルから位相を括り出すことで、片方の成分を実数にできます。

$$
ω
=\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix}
=e^{ia}\begin{pmatrix}\cos\fracθ2\\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix}
=e^{ib}\begin{pmatrix}\cos\fracθ2\,e^{-i(b-a)}\\ \sin\fracθ2\end{pmatrix}
$$

ただし、$ω$の第1成分が$0$になるとき、$e^{ia}$は不明のため括り出すことはできません。同様に、$ω$の第2成分が$0$になるとき、$e^{ib}$を括り出すことはできません。これらを踏まえて定義域から除外した上で、括り出したベクトルに名前を付けます。

&&&def 北側・南側の基準ベクトル
$0 \le θ/2 < 2π$より$0 \le θ < 4π$とします。
$$
\begin{aligned}
χ_N&=\begin{pmatrix}\cos\fracθ2\\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix} && (θ \ne π,3π) \\
χ_S&=\begin{pmatrix}\cos\fracθ2\,e^{-i(b-a)}\\ \sin\fracθ2\end{pmatrix} && (θ \ne 0,2π)
\end{aligned}
$$
&&&

&&&rem
全体位相の括り出しに、基準依存性があるということです。また、$a$だけでは南極（$θ=π,3π$）が、$b$だけでは北極（$θ=0,2π$）が定義域から抜けるため、単一の基準では$S^2$全体を覆えません。
&&&

両方の定義域が重なる領域（$θ\neq0,π,2π,3π$）では、$χ_N,χ_S$の定義を見比べると

$$
χ_N=e^{i(b-a)}χ_S
$$

という関係があります。つまり$b-a$は、南側の位相$b$から北側の位相$a$へ乗り換えるときに掛かる補正（位相のずれ）です。この補正が何を表すかは、ブロッホベクトルの成分から分かります。

$$
\begin{pmatrix}\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ\end{pmatrix}
$$

これは天頂角$θ$、方位角$b-a$による極座標表示なので、補正$b-a$はブロッホベクトルの方位角と一致します。

方位角は点によって値が変わります。たとえば方位角が$0$の点では補正は$0$ですが、方位角が$π$の点では補正は$π$になり、南北の基準がちょうど逆転します。このように、南北の位相の基準を一定の補正だけで揃えることはできず、位置ごとに異なる補正が必要になります。

単一の基準では特異点のため$S^2$全体を覆えず、複数の基準を組み合わせても位相を一意に揃えることはできません。これが「$S^3$を$S^2\times S^1$という単純な直積とみなせない」ことの具体的な中身です。

# まとめ

四元数での回転によるモデルを導入して、線形代数でその結果を利用する流れを示しました。量子情報への導入とすることを意図しています。

&&&rem
ホップファイブレーションの定式化はいくつかあり、射影直線によるものが一般的です。[[dim7-8]]
&&&
