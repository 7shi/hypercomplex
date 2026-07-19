一般次元における球座標系のバリエーションを、成分のパターン分けによって考えます。例として量子情報で使われる2準位系と3準位系の状態ベクトルを紹介します。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

# 概要

実ユークリッド空間のベクトルは長さ$r$と方向を表す単位ベクトル$\vec u$の積として表せます。

&&&def 実ユークリッド空間のベクトル
$$r\vec u\quad(r\in\mathbb R,\ n\in\mathbb N,\ \vec u\in\mathbb R^n,\ \norm{\vec u}=1)$$
&&&

一般次元で$\vec{u}$は単位超球面上の一点を指します。成分は三角関数で表現できます。[[wiki-nsphere]]

&&&def 球座標系から直交座標系への変換
$$\begin{aligned}
x_{1}&=r\cos(\phi _{1})\\x_{2}&=r\sin(\phi _{1})\cos(\phi _{2})\\x_{3}&=r\sin(\phi _{1})\sin(\phi _{2})\cos(\phi _{3})\\&\vdots \\x_{n-1}&=r\sin(\phi _{1})\cdots \sin(\phi _{n-2})\cos(\phi _{n-1})\\x_{n}&=r\sin(\phi _{1})\cdots \sin(\phi _{n-2})\sin(\phi _{n-1})\
\end{aligned}$$
&&&

3次元であれば角度は天頂角や方位角として図形的に理解できますが、4次元以上では図形的な理解は困難です。また、角度の取り方には任意性があり、別の表し方も可能です。

図形的な理解は度外視した上で、三角関数の性質から成分をパターン分けして、バリエーションを考えます。長さ$r$は外して、単位ベクトルに絞ります。

# 2次元

まず2次元について考えます。単位ベクトルは直交座標系の成分の2乗和が1になります。

&&&def 2次元の単位ベクトル
$$\vec{u}_2=\pmatrix{x_1\\x_2}\quad(x_1^2+x_2^2=1)$$
&&&

2乗和が1になることから、三角関数が使えます。

&&&thm 三角関数の2乗和
$$\sin^2ϕ+\cos^2ϕ=1$$
&&&

$x_1,x_2$のどちらを$\sinϕ,\cosϕ$に割り当てるかは任意性があります。また、平方根となるため符号にも任意性があります。円を三角関数でパラメーター表示するときの慣習に従って定義します。

&&&def 単位円上の一点
$$\pmatrix{x_1\\x_2}=\pmatrix{\cosϕ\\\sinϕ}$$
&&&

このように座標を取れば、複素平面で$ϕ=0$のとき単位元$1$になります。また、オイラーの公式によって指数関数で表せます。

&&&def 単位複素数
$$x_1+ix_2=\cosϕ+i\sinϕ=e^{iϕ}$$
&&&

# 3次元

単位ベクトルは直交座標系の成分の2乗和が1になります。次元の数だけ成分があるため、2次元より増えます。

&&&def 3次元の単位ベクトル
$$\vec{u}_3=\pmatrix{x_1\\x_2\\x_3}\quad(x_1^2+x_2^2+x_3^2=1)$$
&&&

これを三角関数の2乗和に当てはめると、3成分を2つに分けることになります。cosに対応する成分の数でパターン分けします。

&&& 成分のパターン分け
$$
(1)\left\{\begin{aligned}x_1^2&=\cos^2ϕ\\x_2^2+x_3^2&=\sin^2ϕ\end{aligned}\right.
\qquad
(2)\left\{\begin{aligned}x_1^2+x_2^2&=\cos^2ϕ\\x_3^2&=\sin^2ϕ\end{aligned}\right.
$$
&&&

通常は(1)を使うことがほとんどです。sinで括った中に2乗和の関係が現れるため、そこにネストして三角関数の関係が入ります。

&&& 三角関数のネスト
$$\begin{aligned}
\pmatrix{x_1\\x_2\\x_3}
&=\pmatrix{\cosϕ \\ \sinϕ\pmatrix{x_2'\\x_3'}}&&(x_2'^2+x_3'^2=1) \\
&=\pmatrix{\cosϕ \\ \sinϕ\pmatrix{\cosϕ' \\ \sinϕ'}} \\
&=\pmatrix{\cosϕ \\ \sinϕ\cosϕ' \\ \sinϕ\sinϕ'}
\end{aligned}$$
&&&

このようにcosに対応する成分を1つ取って、残りをsinで括ってネストさせることを一般次元に拡張したのが定義2です。

# 4次元

2通りのパターンを示します。これらは三次元球面の座標系として扱われます。[[wiki-s3]]

## 1:(1:2)

定義2の方式です。2段階にネストします。

&&&def 4次元の単位ベクトル (1:(1:2))
$$\begin{aligned}
\vec{u}_4
&=\pmatrix{x_1\\x_2\\x_3\\x_4}&&(x_1^2+x_2^2+x_3^2+x_4^2=1) \\
&=\pmatrix{\cosϕ_1 \\ \sinϕ_1\pmatrix{x_2'\\x_3'\\x_4'}}&&(x_2'^2+x_3'^2+x_4'^2=1) \\
&=\pmatrix{\cosϕ_1 \\ \sinϕ_1\pmatrix{\cosϕ_2 \\ \sinϕ_2 \pmatrix{x_3''\\x_4''}}}&&(x_3''^2+x_4''^2=1) \\
&=\pmatrix{\cosϕ_1 \\ \sinϕ_1\pmatrix{\cosϕ_2 \\ \sinϕ_2 \pmatrix{\cosϕ_3 \\ \sinϕ_4}}} \\
&=\pmatrix{\cosϕ_1 \\ \sinϕ_1\cosϕ_2 \\ \sinϕ_1\sinϕ_2\cosϕ_3 \\ \sinϕ_1\sinϕ_2\sinϕ_4}
\end{aligned}$$
&&&

2行目は、単位四元数によって回転角度と回転軸を指定する方式に相当します。[[7shi-h]]

&&&def 単位四元数による回転子
$$\begin{aligned}
q&=\cos\frac\theta2+\sin\frac\theta2(x\,\mathbf{i}+y\,\mathbf{j}+z\,\mathbf{k})&&(x^2+y^2+z^2=1)
\end{aligned}$$
&&&

## 2:2

2成分ごとに分割するパターンを考えます。

&&&def 4次元の単位ベクトル (2:2)
$$\begin{aligned}
\vec{u}_4
&=\pmatrix{x_1\\x_2\\x_3\\x_4}&&(x_1^2+x_2^2+x_3^2+x_4^2=1) \\
&=\pmatrix{\cosϕ_1\pmatrix{x_1'\\x_2'} \\ \sinϕ_1\pmatrix{x_3'\\x_4'}}&&(x_1'^2+x_2'^2=x_3'^2+x_4'^2=1) \\
&=\pmatrix{\cosϕ_1\pmatrix{\cosϕ_2 \\ \sinϕ_2} \\ \sinϕ_1\pmatrix{\cosϕ_3 \\ \sinϕ_3}} \\
&=\pmatrix{\cosϕ_1\cosϕ_2 \\ \cosϕ_1\sinϕ_2 \\ \sinϕ_1\cosϕ_3 \\ \sinϕ_1\sinϕ_3}
\end{aligned}$$
&&&

これを複素数2成分のベクトルとして表したのが、量子情報などで使われる2準位系の状態ベクトルです。

&&&def 2準位系の状態ベクトル
$$
\pmatrix{\cosϕ_1\pmatrix{\cosϕ_2 \\ \sinϕ_2} \\ \sinϕ_1\pmatrix{\cosϕ_3 \\ \sinϕ_3}}
\mapsto\pmatrix{\cosϕ_1(\cosϕ_2+i\sinϕ_2) \\ \sinϕ_1(\cosϕ_3+i\sinϕ_3)} \\
=\pmatrix{\cosϕ_1\,e^{iϕ_2} \\ \sinϕ_1\,e^{iϕ_3}}
$$
&&&

状態ベクトル全体に単位複素数を掛けても測定結果に影響しないため、第1成分を実数化することが行われます。これは次元を下げる操作で、ホップファイブレーションに関係します。（$ϕ_1$のスケールが異なるため、そのものではありません）[[7shi-s]]

&&&def 第1成分の実数化
$$\begin{aligned}
\pmatrix{\cosϕ_1\,e^{iϕ_2} \\ \sinϕ_1\,e^{iϕ_3}}
&=e^{iϕ_2}\pmatrix{\cosϕ_1 \\ \sinϕ_1\,e^{i(ϕ_3-ϕ_2)}}
\mapsto\pmatrix{\cosϕ_1 \\ \sinϕ_1\,e^{iα}}
\end{aligned}$$
&&&

アダマール積（成分ごとの積）によって、2次元の球座標系と各成分の位相に分解できます。

&&&def アダマール積による分解
$$\begin{aligned}
\pmatrix{\cosϕ_1\,e^{iϕ_2} \\ \sinϕ_1\,e^{iϕ_3}}
&=\pmatrix{\cosϕ_1 \\ \sinϕ_1}∘\pmatrix{e^{iϕ_2} \\ e^{iϕ_3}}
\end{aligned}$$
&&&

# 6次元

5次元は飛ばして、6次元の量子情報に関連する例を挙げます。（中嶋慧氏よりご教示）

&&&def 6次元の単位ベクトル (2:(2:2))
$$\begin{aligned}
\vec{u}_6
&=\pmatrix{x_1\\x_2\\x_3\\x_4\\x_5\\x_6}&&(x_1^2+x_2^2+x_3^2+x_4^2+x_5^2+x_6^2=1) \\
&=\pmatrix{\cosϕ_1\pmatrix{x_1'\\x_2'} \\ \sinϕ_1\pmatrix{x_3'\\x_4'\\x_5'\\x_6'}}&&(x_1'^2+x_2'^2=x_3'^2+x_4'^2+x_5'^2+x_6'^2=1) \\
&=\pmatrix{\cosϕ_1\pmatrix{\cosϕ_2 \\ \sinϕ_2} \\ \sinϕ_1\pmatrix{\cosϕ_3\pmatrix{x_3''\\x_4''} \\ \sinϕ_3\pmatrix{x_5''\\x_6''}}}&&(x_3''^2+x_4''^2=x_5''^2+x_6''^2=1) \\
&=\pmatrix{\cosϕ_1\pmatrix{\cosϕ_2 \\ \sinϕ_2} \\ \sinϕ_1\pmatrix{\cosϕ_3\pmatrix{\cosϕ_4 \\ \sinϕ_4} \\ \sinϕ_3\pmatrix{\cosϕ_5 \\ \sinϕ_5}}} \\
&=\pmatrix{\cosϕ_1\cosϕ_2 \\ \cosϕ_1\sinϕ_2 \\ \sinϕ_1\cosϕ_3\cosϕ_4 \\ \sinϕ_1\cosϕ_3\sinϕ_4 \\ \sinϕ_1\sinϕ_3\cosϕ_5 \\ \sinϕ_1\sinϕ_3\cosϕ_5}
\end{aligned}$$
&&&

これを複素化すれば3準位系の状態ベクトルになります。

&&& 3準位系の状態ベクトル
$$\begin{alignedat}{2}
&\pmatrix{\cosϕ_1\pmatrix{\cosϕ_2 \\ \sinϕ_2} \\ \sinϕ_1\pmatrix{\cosϕ_3\pmatrix{\cosϕ_4 \\ \sinϕ_4} \\ \sinϕ_3\pmatrix{\cosϕ_5 \\ \sinϕ_5}}}
&\mapsto&\pmatrix{\cosϕ_1(\cosϕ_2+i\sinϕ_2) \\ \sinϕ_1\cosϕ_3(\cosϕ_4+i\sinϕ_4) \\ \sinϕ_1\sinϕ_3(\cosϕ_5+i\sinϕ_5)} \\
&&=&\pmatrix{\cosϕ_1\,e^{iϕ_2} \\ \sinϕ_1\cosϕ_3\,e^{iϕ_4} \\ \sinϕ_1\sinϕ_3\,e^{iϕ_5}}
\end{alignedat}$$
&&&

2準位系と同様に次元を1つ下げることができます。

&&& 第1成分の実数化
$$\begin{aligned}
&\pmatrix{\cosϕ_1\,e^{iϕ_2} \\ \sinϕ_1\cosϕ_3\,e^{iϕ_4} \\ \sinϕ_1\sinϕ_3\,e^{iϕ_5}}
=e^{iϕ_2}\pmatrix{\cosϕ_1 \\ \sinϕ_1\cosϕ_3\,e^{i(ϕ_4-ϕ_2)} \\ \sinϕ_1\sinϕ_3\,e^{i(ϕ_5-ϕ_2)}}
\mapsto\pmatrix{\cosϕ_1 \\ \sinϕ_1\cosϕ_3\,e^{iα} \\ \sinϕ_1\sinϕ_3\,e^{iβ}} \\
\end{aligned}$$
&&&

アダマール積によって、3次元の球座標系と各成分の位相に分解できます。

&&& アダマール積による分解
$$\begin{aligned}
&\pmatrix{\cosϕ_1\,e^{iϕ_2} \\ \sinϕ_1\cosϕ_3\,e^{iϕ_4} \\ \sinϕ_1\sinϕ_3\,e^{iϕ_5}}
=\pmatrix{\cosϕ_1 \\ \sinϕ_1\cosϕ_3 \\ \sinϕ_1\sinϕ_3}
∘\pmatrix{e^{iϕ_2} \\ e^{iϕ_4} \\ e^{iϕ_5}}
\end{aligned}$$
&&&

# まとめ

高次元の球座標系を考えるときは、2次元に行き着くまで分割します。
