四元数のファイバー自由度を利用して回転子を標準形に変形し、状態ベクトル（スピノル）とそのエルミート共役のテンソル積からブロッホベクトルと密度行列を導出します。

シリーズ：[ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

&&& 改訂履歴
- 2026.10.04 三角関数表示の位相の符号を前回の複素数への対応に合わせ、回転子と状態ベクトルが全体位相$i$を除いて対応することを明記、非自明性を大域的・連続的な同一視の不可能性として述べ直し、その証明は扱わない旨を明記、複素数側の位相を$ζ$で区別、密度行列を「エルミート共役とのテンソル積（外積）」とし純粋状態・混合状態の説明を修正、標準形の$γ$の定義と回転角の説明を補足
- 2026.09.22 北側・南側の基準ベクトルの定義域を天頂角$0\leθ\leπ$に修正
- 2026.07.14 回転子の定義を調整、内積・外積を削除、スピノルを補足説明、基準依存性を追記
&&&

# 概要

前回は、単位四元数による回転を通してホップファイブレーションを構成し、複素数表示を経てパウリ行列を導出しました。[[7shi-h]]

しかし、四元数の回転子$ω$には始点$\mathbf{k}$を固定する回転（円周$S^1$）の自由度が残されており、このファイバー方向の自由度が状態ベクトルの全体位相とどう結びついているのか、またベクトルとスピノルの関係が代数的にどう捉えられるかという問いが残されていました。

本記事では、まずファイバー自由度$q$を用いて回転子$ω$の$\mathbf{k}$成分を消去した標準形を求め、三角関数表示を通して、これが状態ベクトル（スピノル）の全体位相の消去に対応することを示します。次に、パウリ行列を用いてブロッホ球上の座標（ブロッホベクトル）を成分ごとに計算し、さらにエルミート共役とのテンソル積（外積）$\Psi\Psi^\dagger$として一括計算することで、密度行列が自然に現れることを見ます。これにより、ブロッホベクトルがスピノルのテンソル積から得られる派生的な対象であることを明らかにします。最後に、状態ベクトルから全体位相を括り出す操作には基準依存性があり、単一の基準では球面全体を覆えないことを確かめます。これは、束が自明な直積$S^2\times S^1$にならないことの具体的な現れです。

前回の四元数によるホップファイブレーションの構成を前提とします。本記事では1量子ビットに相当する状態ベクトルとブロッホ球、密度行列の幾何学的関係の解明に集中し、複数量子ビットのもつれや一般の量子測定理論の詳細は扱いません。

# パラメーターの調整

単位四元数による回転子$ω$と$q$を定義します。

&&&def 単位四元数による回転子
$$\begin{aligned}
ω&=(α_0\mathbf{k}+α_1)+\mathbf{j}(β_0\mathbf{k}+β_1)&&(α_0^2+α_1^2+β_0^2+β_1^2=1) \\
q&=u+v\mathbf{k}&&(u^2+v^2=1)
\end{aligned}$$
&&&

$q\mathbf kq^*=\mathbf k$より、$q$が$\mathbf k$を固定するファイバーを構成するため、次が成り立ちます。

$$
(ωq)\mathbf k(ωq)^*=ω\mathbf kω^*
$$

この自由度を使って、$ω$の$\mathbf k$の項を消去した標準形を作ります。

$ωq$を計算します。

$$\begin{aligned}
ωq&=\{(α_0\mathbf{k}+α_1)+\mathbf{j}(β_0\mathbf{k}+β_1)\}(u+v\mathbf{k}) \\
&=(α_0\mathbf{k}+α_1)(u+v\mathbf{k})+\mathbf{j}(β_0\mathbf{k}+β_1)(u+v\mathbf{k})
\end{aligned}$$

$u,v$を調整すれば、$\mathbf k$を消すことができます。$γ=\sqrt{α_0^2+α_1^2}$とします。$γ=0$のときは、もとの$ω$にすでに$\mathbf k$の項がないため、$q=1$とすれば足ります。$γ>0$のときは$α_0'=α_0/γ,\ α_1'=α_1/γ$とおき、次のように$q$を選びます。

$$\begin{aligned}
α_0\mathbf{k}+α_1&=γ(α_0'\mathbf{k}+α_1')&&(α_0'^2+α_1'^2=1) \\
u+v\mathbf{k}&=-α_0'\mathbf{k}+α_1'
\end{aligned}$$
$$\begin{aligned}
ωq&=γ(α_0'\mathbf{k}+α_1')(-α_0'\mathbf{k}+α_1')+\mathbf{j}(β_0\mathbf{k}+β_1)(-α_0'\mathbf{k}+α_1') \\
  &=γ(α_0'^2+α_1'^2)+\mathbf{j}\{(β_0α_0'+β_1α_1')+(β_0α_1'-β_1α_0')\mathbf{k}\} \\
  &=γ+(β_0α_0'+β_1α_1')\mathbf{j}+(β_0α_1'-β_1α_0')\mathbf{i} \\
\end{aligned}$$

つまり、$\mathbf k$を$ω$と同じ点へ移す$ωq$のうち、$\mathbf k$の項がないものが存在します。

## 三角関数表示

成分計算だけではイメージが湧きにくいため、三角関数でパラメーター表示します。[[7shi-coord]]

&&&def 回転子の三角関数表示
$$
\begin{aligned}
ω&=\cosθ\,(\cos a-\sin a\,\mathbf{k})+\sinθ\,\mathbf{j}(\cos b-\sin b\,\mathbf{k}) \\
 &=\cosθ\,e^{-a\,\mathbf{k}}+\sinθ\,\mathbf{j}\,e^{-b\,\mathbf{k}} \\
q&=\cos c+\sin c\,\mathbf{k} \\
 &=e^{c\,\mathbf{k}}
\end{aligned}
$$
&&&

指数の位相に負号を付けているのは、後で状態ベクトルの位相$e^{ia},e^{ib}$と対応させるためです。

$ωq$を計算します。

$$
\begin{aligned}
ωq&=(\cosθ\,e^{-a\,\mathbf{k}}+\sinθ\,\mathbf{j}\,e^{-b\,\mathbf{k}})e^{c\,\mathbf{k}} \\
  &=\cosθ\,e^{(c-a)\,\mathbf{k}}+\sinθ\,\mathbf{j}\,e^{(c-b)\,\mathbf{k}}
\end{aligned}
$$

$c=a$とすると、次のようになります。

$$\begin{aligned}
ωq&=\cosθ+\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}} \\
  &=\cosθ+\sinθ\,\mathbf{j}\left\{\cos(a-b)+\sin(a-b)\,\mathbf{k}\right\} \\
  &=\cosθ+\sinθ\left\{\cos(a-b)\,\mathbf{j}+\sin(a-b)\,\mathbf{i}\right\}
\end{aligned}$$

よって、$q$が$e^{-a\,\mathbf{k}}$を打ち消すとき、$\mathbf k$の項が消えます。これを$ω'$とします。

&&&def k のない回転子
$$
\begin{aligned}
ω'&=\cosθ+\sinθ\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}} \\
\end{aligned}
$$
&&&

## 回転角の調整

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
&=\sin2θ\left\{\cos(b-a)\mathbf{i}+\sin(b-a)\mathbf{j}\right\}+\cos2θ\,\mathbf{k}
\end{aligned}$$

左右から回転子で挟むため、回転子のパラメーター$θ$に対して、実際の回転角は$2θ$になります。

最終的に得られる回転角が$θ$になるように、以降では回転子の側で角度を半分にして調整します。

$$
ω''=\cos\fracθ2+\sin\fracθ2\,\mathbf{j}\,e^{(a-b)\,\mathbf{k}}
$$

&&&rem 回転角の半角調整
最初から半角を用いるのではなく、回転子のパラメーターの2倍が実際の回転角になることを計算で確かめてから、$θ/2$に置き換えました。
&&&

# 状態ベクトルとパウリ行列

同じことを状態ベクトルとパウリ行列で計算します。この複素ベクトルは**スピノル**とも呼ばれます。

半角に調整した回転子$\cos\fracθ2\,e^{-a\mathbf k}+\sin\fracθ2\,\mathbf j\,e^{-b\mathbf k}$から、前回の対応$α=α_0+iα_1,\ β=β_0+iβ_1$で複素数を作ると、次のようになります。

$$
\begin{pmatrix}α\\β\end{pmatrix}
=\begin{pmatrix}α_0+iα_1\\β_0+iβ_1\end{pmatrix}
=i\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix}
$$

係数$i$は両成分に共通に掛かる定数の位相で、$(e^{iδ}\Psi)^\dagger σ(e^{iδ}\Psi)=\Psi^\dagger σ\Psi$よりパウリ行列による座標には影響しません。そこで、これを除いたものを状態ベクトル$\Psi$とします。

$$\begin{aligned}
\Psi&=\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix}
  =e^{ia}\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix} \\
ζ&=e^{-ia} \\
\Psi ζ&=\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix}=:\Psi'
\end{aligned}$$

四元数の$e^{c\mathbf k}$を右から掛けることは、前回の対応では複素数の組に$e^{-ic}$を掛けることに当たります。したがって$\mathbf k$の項を消す四元数の$q=e^{a\mathbf k}$は複素数の$ζ=e^{-ia}$に対応し、$\Psi ζ$は括り出した係数の消去に相当します。また、$ω''$から前回の対応で作った複素数の組は$i\Psi'$となり、$\Psi'$と全体位相$i$だけ異なります。

パウリ行列でブロッホ球上の座標に変換します。

$$\begin{aligned}
\Psi'^\dagger \sigma_x \Psi'
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}0&1\\1&0\end{pmatrix}\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix} \\
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}\sin\fracθ2\,e^{i(b-a)} \\ \cos\fracθ2\end{pmatrix} \\
&=\cos\fracθ2\sin\fracθ2\,e^{i(b-a)} + \sin\fracθ2\,e^{-i(b-a)}\cos\fracθ2 \\
&=\cos\fracθ2\sin\fracθ2\,(e^{i(b-a)}+e^{-i(b-a)}) \\
&=\cos\fracθ2\sin\fracθ2\left\{\cos(b-a)+i\sin(b-a)+\cos(b-a)-i\sin(b-a)\right\} \\
&=2\cos\fracθ2\sin\fracθ2\cos(b-a) \\
&=\sinθ\cos(b-a) \\

\Psi'^\dagger \sigma_y \Psi'
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}0&-i\\i&0\end{pmatrix}\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix} \\
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}-i\sin\fracθ2\,e^{i(b-a)} \\ i\cos\fracθ2\end{pmatrix} \\
&=-i\cos\fracθ2\sin\fracθ2\,e^{i(b-a)} + i\sin\fracθ2\,e^{-i(b-a)}\cos\fracθ2 \\
&=-i\cos\fracθ2\sin\fracθ2\,(e^{i(b-a)}-e^{-i(b-a)}) \\
&=-i\cos\fracθ2\sin\fracθ2\left\{\cos(b-a)+i\sin(b-a)-\cos(b-a)+i\sin(b-a)\right\} \\
&=2\cos\fracθ2\sin\fracθ2\sin(b-a) \\
&=\sinθ\sin(b-a) \\

\Psi'^\dagger \sigma_z \Psi'
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}1&0\\0&-1\end{pmatrix}\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix} \\
&=\begin{pmatrix}\cos\fracθ2 & \sin\fracθ2\,e^{-i(b-a)}\end{pmatrix}\begin{pmatrix}\cos\fracθ2 \\ -\sin\fracθ2\,e^{i(b-a)}\end{pmatrix} \\
&=\cos^2\fracθ2-\sin^2\fracθ2\,e^{-i(b-a)}e^{i(b-a)} \\
&=\cosθ
\end{aligned}$$

四元数は全成分を一度に計算するため式が込み入りますが、こちらは成分ごとに計算するため冗長になります。どちらが簡単だとは一概には言えませんが、結果だけまとめるのにベクトルはシンプルです。

$$\begin{aligned}
\Psi'&=\begin{pmatrix}\cos\fracθ2 \\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix}
\mapsto\begin{pmatrix}\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ\end{pmatrix}
\end{aligned}$$

スピノルに現れる半角$θ/2$から天頂角$θ$が得られ、指数関数が表す$\cos$と$\sin$が別成分に分離されました。これを**ブロッホベクトル**と呼びます。

$\Psi$を計算しても同じ結果になります。

$$\begin{aligned}
\Psi&=\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix}
\mapsto\begin{pmatrix}\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ\end{pmatrix}
\end{aligned}$$

このように、ブロッホベクトルには$a,b$個々の位相ではなく、相対位相$b-a$だけが現れます。全体位相はブロッホベクトルに影響しません。

## 全体位相とファイバー

状態ベクトル全体に位相$e^{iφ}$を掛けると、次のようになります。

$$
e^{iφ}\Psi=\begin{pmatrix}\cos\fracθ2\,e^{i(a+φ)} \\ \sin\fracθ2\,e^{i(b+φ)}\end{pmatrix}
$$

$a,b$がともに$φ$だけ動きます。位相差$b-a$は変わらないため、ブロッホベクトルも変わりません。

これは$\Psi ζ$における位相$ζ$、すなわち四元数のファイバー自由度$q$に対応します。つまり状態ベクトルに全体位相を掛ける操作は、ホップファイブレーションのファイバーに沿って点を動かす操作に対応します。ブロッホベクトルがこの操作で変わらないのは、ファイバー上の点$\Psi ζ$（$ζ$を動かして得られる各点）がすべて同じ$S^2$上の点に写ることの表れです。

&&&rem ファイバーという用語の二面性
「ファイバー」という用語は、ホップファイブレーションの構造として抽象的に円周$S^1$を指す場合と、$S^3$内に埋め込まれた具体的な部分集合$\{\Psi ζ:|ζ|=1\}$を指す場合の両方に使われます。$ζ$は$S^1$上を動くパラメーターで、$\Psi ζ$は$S^3$上の点です。
&&&

## まとめて計算

ここまでは成分ごとに座標を求めましたが、パウリ行列を使えば全成分をまとめて計算することもできます。この計算は、スピノル$\Psi$とそのエルミート共役のテンソル積（外積）によってブロッホベクトルを作り出す操作になっています。

$A,B$を複素ベクトルとして、$A^\dagger B$は内積となります。エルミート共役を逆にした$AB^\dagger$は外積と呼ばれます。

&&&rem 外積の用語の区別
ベクトル解析で一般的なベクトル積（クロス積）や、ウェッジ積とは別の外積で、テンソル積（クロネッカー積）の特殊な場合です。日本語では同じ訳語が用いられますが、英語では区別があります。

英語|日本語|備考
----|----|----
outer product<br>tensor direct product|外積<br>直積|テンソル積（クロネッカー積）の特殊な場合<br>直積も他の用法と紛らわしくマイナー
exterior product<br>wedge product|外積、外部積<br>ウェッジ積、楔積|外積代数で使用、次元に依存しない<br>3次元ではベクトル積のホッジ双対
<br>vector product<br>cross product|外積<br>ベクトル積<br>クロス積|ベクトル解析で使用<br>次元に依存する（3次元と7次元）<br>純虚四元数・純虚八元数どうしの積の虚部から得られる

外積の多様な用法については以前の記事を参照してください。[[7shi-p]]

直積の用語については以下を参照してください。[[wiki-dp]][[hyodo-dp]]
&&&

（テンソル積の意味での）外積はベクトルから行列を作る計算で、エルミート共役との外積は成分と複素共役との積を列挙したものとなります。以下では記号を簡潔にするため、改めて$\Psi$の成分を$α,β$と書きます。これらは、もとの四元数から前回の対応で直接得た成分の$-i$倍です。

$$\begin{aligned}
\Psi&=\begin{pmatrix}α\\β\end{pmatrix}=\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix} \\
\Psi\Psi^\dagger
&=\begin{pmatrix}α\\β\end{pmatrix}\begin{pmatrix}α^*&β^*\end{pmatrix} \\
&=\begin{pmatrix}αα^*&αβ^*\\βα^*&ββ^*\end{pmatrix} \\
&=\begin{pmatrix}
  \cos\fracθ2\,e^{ia}\cos\fracθ2\,e^{-ia} & \cos\fracθ2\,e^{ia}\sin\fracθ2\,e^{-ib} \\
  \sin\fracθ2\,e^{ib}\cos\fracθ2\,e^{-ia} & \sin\fracθ2\,e^{ib}\sin\fracθ2\,e^{-ib}\end{pmatrix} \\
&=\begin{pmatrix}
  \cos^2\fracθ2 & \cos\fracθ2\sin\fracθ2\,e^{-i(b-a)} \\
  \cos\fracθ2\sin\fracθ2\,e^{i(b-a)} & \sin^2\fracθ2\end{pmatrix}
\end{aligned}$$

成分に見覚えのある形が現れています。半角公式などを使って計算を進めます。[[half]]

$$\begin{aligned}
\Psi\Psi^\dagger
&=\frac12\begin{pmatrix}
  1+\cosθ & \sinθ\left\{\cos(b-a)-i\sin(b-a)\right\} \\
  \sinθ\left\{\cos(b-a)+i\sin(b-a)\right\} & 1-\cosθ\end{pmatrix}
\end{aligned}$$

対角成分と非対角成分で、それぞれ符号が異なる項が現れています。これらを分離すれば、パウリ行列の線形結合が現れます。

&&&fml 密度行列のパウリ分解
$$\begin{aligned}
\Psi\Psi^\dagger
&=\frac12\left\{
  \begin{pmatrix}1&0 \\ 0&1\end{pmatrix} \\
  +\sinθ\cos(b-a)\begin{pmatrix}0&1 \\ 1&0\end{pmatrix} \\
  +\sinθ\sin(b-a)\begin{pmatrix}0&-i \\ i&0\end{pmatrix} \\
  +\cosθ\begin{pmatrix}1&0 \\ 0&-1\end{pmatrix}\right\} \\
&=\frac12\left\{ \mathrm{I} + \sinθ\cos(b-a)\sigma_x + \sinθ\sin(b-a)\sigma_y + \cosθ\,\sigma_z \right\}
\end{aligned}$$
&&&

これは**密度行列**と呼ばれます。スケールを$2$倍してパウリ行列の係数を抜き出せば、ブロッホベクトルと一致します。

$$\begin{aligned}
\Psi&=\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix}
\mapsto\begin{pmatrix}\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ\end{pmatrix}
\end{aligned}$$

&&&rem 始点
四元数による回転の始点を$\mathbf k$としたのは、外積から計算される成分と一致させるためです。
&&&

$\Psi\Psi^\dagger$は密度行列という形を取っていますが、表している情報は、単位行列の部分とスケールの違いを除けばブロッホベクトルと同じです。つまりテンソル積$\Psi\Psi^\dagger$は、ブロッホベクトルの別の表現に過ぎません。

この見方に立つと、ブロッホベクトルは**スピノルのテンソル積**から得られる派生的な対象であり、スピノル$\Psi$の方がより基本的な対象だと言えます。実際、$\Psi$に全体位相$e^{iφ}$を掛けても$\Psi\Psi^\dagger$は変わらないため、$\Psi$はブロッホベクトルよりも多くの情報（全体位相）を持っています。逆にブロッホベクトルからスピノルを復元しようとしても、全体位相の分だけ情報が不足していて一意に定まりません。

&&&rem 密度行列の状態
1つの状態ベクトルから作られる密度行列$\Psi\Psi^\dagger$を**純粋状態**と呼びます。それに対して、異なる純粋状態の密度行列を確率で重み付けして足し合わせたものを**混合状態**と呼びます。

純粋状態の密度行列に対応するブロッホベクトルはブロッホ球の表面を指しますが、混合状態ではブロッホ球の内部を指します。[[quantumuniverse]]
&&&

# 全体位相の基準依存性

ホップファイブレーションが定める束$S^1\hookrightarrow S^3\xrightarrow{H}S^2$は自明束ではありません。すなわち、射影と整合する形で、$S^3$を直積$S^2\times S^1$と大域的かつ連続的に同一視することはできません。この証明は本記事では扱いませんが、全体位相の括り出しを通して、この性質がどのように現れるかを具体的に確認します。

先ほど見た通り、状態ベクトルから位相を括り出すことで、片方の成分を実数にできます。

$$
\Psi
=\begin{pmatrix}\cos\fracθ2\,e^{ia} \\ \sin\fracθ2\,e^{ib}\end{pmatrix}
=e^{ia}\begin{pmatrix}\cos\fracθ2\\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix}
=e^{ib}\begin{pmatrix}\cos\fracθ2\,e^{-i(b-a)}\\ \sin\fracθ2\end{pmatrix}
$$

ただし、$\Psi$の第1成分が$0$になるとき、その成分から位相$e^{ia}$を決められず、この成分を基準にした代表の選び方が定まりません。同様に、$\Psi$の第2成分が$0$になるとき、第2成分を基準にした代表の選び方は定まりません。これらを踏まえて定義域から除外した上で、括り出したベクトルに名前を付けます。

&&&def 北側・南側の基準ベクトル
$θ$はブロッホ球の天頂角として$0\leθ\leπ$とします。
$$
\begin{aligned}
χ_N&=\begin{pmatrix}\cos\fracθ2\\ \sin\fracθ2\,e^{i(b-a)}\end{pmatrix} && (0\leθ<π) \\
χ_S&=\begin{pmatrix}\cos\fracθ2\,e^{-i(b-a)}\\ \sin\fracθ2\end{pmatrix} && (0<θ\leπ)
\end{aligned}
$$
&&&

&&&rem 全体位相の基準依存性
全体位相の括り出しに、基準依存性があるということです。また、$a$だけでは南極（$θ=π$）が、$b$だけでは北極（$θ=0$）が定義域から抜けるため、ここで用いた各基準だけでは$S^2$全体を覆えません。
&&&

両方の定義域が重なる領域（$0<θ<π$）では、$χ_N,χ_S$の定義を見比べると次の関係があります。

$$
χ_N=e^{i(b-a)}χ_S
$$

つまり$b-a$は、南側の位相$b$から北側の位相$a$へ乗り換えるときに掛かる補正（位相のずれ）です。この補正が何を表すかは、ブロッホベクトルの成分から分かります。

$$
\begin{pmatrix}\sinθ\cos(b-a) \\ \sinθ\sin(b-a) \\ \cosθ\end{pmatrix}
$$

これは天頂角$θ$、方位角$b-a$による極座標表示なので、補正$b-a$はブロッホベクトルの方位角と一致します。

方位角は点によって値が変わります。たとえば方位角が$0$の点では補正は$0$ですが、方位角が$π$の点では補正は$π$になり、南北の基準がちょうど逆転します。このように、南北の位相の基準を一定の補正だけで揃えることはできず、位置ごとに異なる補正が必要になります。

ここで用いた各基準は特異点のため単独では$S^2$全体を覆えず、南北の基準の間には位置ごとに異なる補正が必要になります。これは「$S^3$を$S^2\times S^1$と大域的かつ連続的に同一視できない」ことの具体的な現れです。ただし、ここで確かめたのは南北の2つの基準についてだけです。非自明性の証明には、赤道上で補正$e^{i(b-a)}$が1周することと、この周回数が南北それぞれの領域全体での基準の連続的な変更では変わらないことを用います。本記事ではその証明には立ち入りません。

# まとめ

四元数の回転子$ω$はファイバー方向の自由度$q$を持ち、この自由度を使って$\mathbf k$の項を持たない標準形$ω'$を示しました。$q$による$ωq$の調整は、状態ベクトルから全体位相$e^{iφ}$を括り出す操作と同じ構造であり、回転角を半角に調整した標準形$ω''$は、全体位相$i$を除いて状態ベクトル$(\cos\fracθ2,\ \sin\fracθ2\,e^{i(b-a)})$に対応します。パウリ行列によってブロッホベクトルへ変換すると、四元数による回転と状態ベクトルによる記述が同じ数学的構造を共有していることが分かります。

成分ごとの計算と、エルミート共役とのテンソル積$\Psi\Psi^\dagger$による一括計算の両方を行い、後者が密度行列の形を与え、単位行列の部分とスケールを除けばブロッホベクトルと一致することを確認しました。ここから、ブロッホベクトルはスピノル$\Psi$のテンソル積から得られる派生的な対象であり、全体位相の情報を保持するスピノルの方がより基本的な対象であることが分かりました。

さらに全体位相の括り出しには基準依存性があり、北側・南側いずれの基準も球面上の一点（それぞれ南極・北極）では定義できないこと、両者を結ぶ補正がブロッホベクトルの方位角そのものであることを示しました。これは、ホップファイブレーションが与える束$S^1\hookrightarrow S^3\to S^2$が単純な直積$S^2\times S^1$にはならないという事実の、具体的な現れです。

このように、四元数での回転によるモデルを導入して、線形代数でその結果を利用する流れを示しました。量子情報への導入とすることを意図しています。
