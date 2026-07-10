実数体・複素数体・四元数体と分解型の代数を基本構造として、テンソル積による拡張公式を繰り返し適用することで、クリフォード代数$\operatorname{Cl}_{p,q}(\mathbb R)$を行列環として分類します。拡張公式は前回の記事で導出したものを使用します。[[qua-tensor]]

# 基本的な代数構造

クリフォード代数には、5つの基本的な代数構造があります。

$$
\begin{alignedat}{3}
&\mathbb{R} &&\cong \operatorname{Cl}_{0,0}(\mathbb{R}) &&\quad\text{（実数体）} \\
&\mathbb{C} &&\cong \operatorname{Cl}_{0,1}(\mathbb{R}) &&\quad\text{（複素数体）} \\
&\mathbb{C}' &&\cong \operatorname{Cl}_{1,0}(\mathbb{R}) &&\quad\text{（分解型複素数）} \\
&\mathbb{H} &&\cong \operatorname{Cl}_{0,2}(\mathbb{R}) &&\quad\text{（四元数体）} \\
&\mathbb{H}' &&\cong \operatorname{Cl}_{2,0}(\mathbb{R}) \cong \operatorname{Cl}_{1,1}(\mathbb{R}) &&\quad\text{（分解型四元数）}
\end{alignedat}
$$

&&&rem
$\mathbb H,\mathbb H'$とクリフォード代数の対応は、前回の記事で確認しました。[[qua-tensor]]

$\mathbb C,\mathbb C'$は生成元が1個のクリフォード代数に対応します。$i^2=-1$となる虚数単位$i$を生成元とすれば$\operatorname{Cl}_{0,1}(\mathbb R)$、$j^2=1$となる虚数単位$j$を生成元とすれば$\operatorname{Cl}_{1,0}(\mathbb R)$が得られます。
&&&

これらの基本構造は、より高次のクリフォード代数を構成する際の基礎となります。分解型の2つの代数は実数の直和や実行列環として表せるため、以下で確認します。

## 分解型複素数の直和分解

分解型複素数 $\mathbb{C}'$ は、通常の複素数とは異なる性質を持ちます。その最も基本的な特徴は、虚数単位 $j$ が $j^2 = +1$ を満たすことです（通常の虚数単位 $i$ は $i^2 = -1$ を満たす）。

この代数は、2つの冪等元
$$
e = \dfrac{1-j}{2},\ e^* = \dfrac{1+j}{2}
$$
を用いて直和分解することができます。

この冪等元は以下の重要な性質を持ちます：

- 冪等性：$e^2 = e$, $(e^*)^2 = e^*$
- 直交性：$ee^* = e^*e = 0$

これにより、任意の分解型複素数 $z = x + jy$ は以下のように一意的に分解できます：

$$z = (x-y)e + (x+y)e^*$$

この直和分解を用いると、2つの分解型複素数の積は非常に簡単に計算できます：

$$
(a_1e + b_1e^*)(a_2e + b_2e^*) = a_1a_2e + b_1b_2e^*
$$

この計算規則は成分ごとの積に対応しており、分解型複素数が実数の直和と同型であることを示しています：

$$\mathbb{C}' \cong \mathbb{R} \oplus \mathbb{R} \cong 2\mathbb{R}$$

&&&rem
以下、代数$A$の直和$A \oplus A$を$2A$と略記します。
&&&

## 分解型四元数の行列表現

分解型四元数 $\mathbb{H}'$ は2次の実行列環と同型です：

$$\mathbb{H}' \cong \mathbb{R}(2)$$

&&&rem
以下、$\mathbb F$を成分とする$n$次の全行列環を$\mathbb F(n)$と表記します。
&&&

その生成元は以下の行列で表されます：

$$
i \cong \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \quad
j \cong \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
k = ij \cong \begin{pmatrix} -1 & 0 \\ 0 & 1 \end{pmatrix}
$$

この行列表現により、任意の分解型四元数 $q = a + bi + cj + dk$ は以下の2×2実行列として具体的に表現できます：

$$
q \cong \begin{pmatrix} a-d & -b+c \\ b+c & a+d \end{pmatrix}
$$

この表現により、分解型四元数の演算は2×2実行列の演算として具体的に計算することができます。

# テンソル積による拡張公式

クリフォード代数は、$\mathbb H,\mathbb H'$とのテンソル積によって拡張できます。前回導出した公式を再掲します。[[qua-tensor]]

&&&fml クリフォード代数の$⊗\mathbb H$による拡張
$$
\operatorname{Cl}_{p,q}(\mathbb R) ⊗ \mathbb H
\cong \operatorname{Cl}_{q,p+2}(\mathbb R)
$$
&&&

&&&fml クリフォード代数の$⊗\mathbb H'$による拡張
$$
\operatorname{Cl}_{p,q}(\mathbb R) ⊗ \mathbb H'
\cong \operatorname{Cl}_{q+2,p}(\mathbb R)
\cong \operatorname{Cl}_{p+1,q+1}(\mathbb R)
$$
&&&

&&&rem
テンソル積によって、クリフォード代数としての生成元は2個増えます。右辺で$p,q$の位置が入れ替わるのは、拡張時に元の生成元の計量が反転することに由来します。詳細は前回の記事を参照してください。[[qua-tensor]]
&&&

これらの公式を基本構造に繰り返し適用すれば、より高次のクリフォード代数を系統的に構成できます。

## テンソル積の計算規則

公式を適用した結果を行列環として表すために、以下の同型を使用します。

1. 可換性：

$$A \otimes B \cong B \otimes A$$

2. 行列環とのテンソル積：

$$\mathbb{F}(m) \otimes \mathbb{G}(n) \cong (\mathbb{F} \otimes \mathbb{G})(mn)$$

特に$\mathbb F \otimes \mathbb R \cong \mathbb F$より$\mathbb{F}(m) \otimes \mathbb{R}(n) \cong \mathbb{F}(mn)$となるため、$\mathbb H' \cong \mathbb R(2)$とのテンソル積は行列表現の次数を2倍にします。

3. 直和との分配性：

$$2A \otimes B \cong 2(A \otimes B)$$

4. 四元数とのテンソル積：

$$
\begin{aligned}
\mathbb{H} \otimes \mathbb{H} &\cong \mathbb{H}' \otimes \mathbb{H}' \cong \mathbb{R}(4) \\
\mathbb{H}' \otimes \mathbb{H} &\cong \mathbb{R}(2) \otimes \mathbb{H} \cong \mathbb{H}(2) \\
\mathbb{C} \otimes \mathbb{H} &\cong \mathbb{C} \otimes \mathbb{H}' \cong \mathbb{C}(2)
\end{aligned}
$$

&&&rem
1行目：$\mathbb{H} \otimes \mathbb{H} \cong \mathbb{H}' \otimes \mathbb{H}'$は前回示した同型です[[qua-tensor]]。$\mathbb H' \cong \mathbb R(2)$より、$\mathbb{H}' \otimes \mathbb{H}' \cong \mathbb{R}(2) \otimes \mathbb{R}(2) \cong \mathbb{R}(4)$と計算できます。

3行目：$\mathbb{C} \otimes \mathbb{H}$は双四元数であり、パウリ行列との対応によって$\mathbb{C}(2)$と同型になることを以前の記事で確認しました[[pauli-qua]]。一方、$\mathbb{C} \otimes \mathbb{H}' \cong \mathbb{C} \otimes \mathbb{R}(2) \cong \mathbb{C}(2)$となるため、両者は同型です。
&&&

これらの規則は組み合わせて使用できます。例えば：

$$\mathbb{H}(2) \otimes \mathbb{H} \cong (\mathbb{H} \otimes \mathbb{H})(2) \cong \mathbb{R}(4)(2) \cong \mathbb{R}(8)$$

# 系列による代数構造の導出

基本構造に2つの公式を交互に適用して、分類表の第1行（$p=0$）と第1列（$q=0$）を導出します。符号数$(p,0)$に公式1を適用すると$(0,p+2)$へ、$(0,q)$に公式2を適用すると$(q+2,0)$へ移るため、交互に適用することで$p$軸と$q$軸を往復しながら次数が上がっていきます。

## 複素数体からの系列

複素数体 $\mathbb{C} \cong \operatorname{Cl}_{0,1}(\mathbb{R})$ を出発点として、テンソル積による拡張を繰り返すと、以下の系列が得られます：

$$
\begin{aligned}
\operatorname{Cl}_{0,1}(\mathbb{R}) &\cong \mathbb{C} \\
\operatorname{Cl}_{3,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,1}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{C} \otimes \mathbb{R}(2) \cong \mathbb{C}(2) \\
\operatorname{Cl}_{0,5}(\mathbb{R}) &\cong \operatorname{Cl}_{3,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{C}(2) \otimes \mathbb{H} \cong \mathbb{C}(4) \\
\operatorname{Cl}_{7,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,5}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{C}(4) \otimes \mathbb{R}(2) \cong \mathbb{C}(8)
\end{aligned}
$$

## 分解型複素数からの系列

分解型複素数 $\mathbb{C}' \cong \operatorname{Cl}_{1,0}(\mathbb{R})$ からは、以下の系列が導かれます：

$$
\begin{aligned}
\operatorname{Cl}_{1,0}(\mathbb{R}) &\cong \mathbb{C}' \cong 2\mathbb{R} \\
\operatorname{Cl}_{0,3}(\mathbb{R}) &\cong \operatorname{Cl}_{1,0}(\mathbb{R}) \otimes \mathbb{H} \cong 2\mathbb{R} \otimes \mathbb{H} \cong 2\mathbb{H} \\
\operatorname{Cl}_{5,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,3}(\mathbb{R}) \otimes \mathbb{H}' \cong 2\mathbb{H} \otimes \mathbb{R}(2) \cong 2\mathbb{H}(2) \\
\operatorname{Cl}_{0,7}(\mathbb{R}) &\cong \operatorname{Cl}_{5,0}(\mathbb{R}) \otimes \mathbb{H} \cong 2\mathbb{H}(2) \otimes \mathbb{H} \cong 2\mathbb{R}(8)
\end{aligned}
$$

## 四元数体からの系列

四元数体 $\mathbb{H} \cong \operatorname{Cl}_{0,2}(\mathbb{R})$ からは、以下の系列が生成されます：

$$
\begin{aligned}
\operatorname{Cl}_{0,2}(\mathbb{R}) &\cong \mathbb{H} \\
\operatorname{Cl}_{4,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,2}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{H} \otimes \mathbb{R}(2) \cong \mathbb{H}(2) \\
\operatorname{Cl}_{0,6}(\mathbb{R}) &\cong \operatorname{Cl}_{4,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{H}(2) \otimes \mathbb{H} \cong \mathbb{R}(8) \\
\operatorname{Cl}_{8,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,6}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{R}(8) \otimes \mathbb{R}(2) \cong \mathbb{R}(16)
\end{aligned}
$$

## 分解型四元数からの系列

分解型四元数 $\mathbb{H}' \cong \operatorname{Cl}_{2,0}(\mathbb{R})$ からは、以下の系列が得られます：

$$
\begin{aligned}
\operatorname{Cl}_{2,0}(\mathbb{R}) &\cong \mathbb{H}' \cong \mathbb{R}(2) \\
\operatorname{Cl}_{0,4}(\mathbb{R}) &\cong \operatorname{Cl}_{2,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{R}(2) \otimes \mathbb{H} \cong \mathbb{H}(2) \\
\operatorname{Cl}_{6,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,4}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{H}(2) \otimes \mathbb{R}(2) \cong \mathbb{H}(4) \\
\operatorname{Cl}_{0,8}(\mathbb{R}) &\cong \operatorname{Cl}_{6,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{H}(4) \otimes \mathbb{H} \cong \mathbb{R}(16)
\end{aligned}
$$

## 系列の整理

これらの系列をまとめると、以下の表が得られます：

$$
\begin{array}{|c|ccccccccc|}
\hline
p \backslash q & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 \\
\hline
0 & \mathbb{R} & \mathbb{C} & \mathbb{H} & 2\mathbb{H} & \mathbb{H}(2) & \mathbb{C}(4) & \mathbb{R}(8) & 2\mathbb{R}(8) & \mathbb{R}(16) \\
1 & 2\mathbb{R} & & & & & & & & \\
2 & \mathbb{R}(2) & & & & & & & & \\
3 & \mathbb{C}(2) & & & & & & & & \\
4 & \mathbb{H}(2) & & & & & & & & \\
5 & 2\mathbb{H}(2) & & & & & & & & \\
6 & \mathbb{H}(4) & & & & & & & & \\
7 & \mathbb{C}(8) & & & & & & & & \\
8 & \mathbb{R}(16) & & & & & & & & \\
\hline
\end{array}
$$

この表の空白部分は、公式2の$\operatorname{Cl}_{p+1,q+1}(\mathbb{R})$の形を用いて埋めることができます。

# 分類表の完成

## 対角線方向の充填

公式2は、$\mathbb H'$とのテンソル積によって符号数が$(p,q)$から$(p+1,q+1)$へ拡張できることを示しています。$\mathbb H' \cong \mathbb R(2)$より：

$$\operatorname{Cl}_{p+1,q+1}(\mathbb{R}) \cong \operatorname{Cl}_{p,q}(\mathbb{R}) \otimes \mathbb{R}(2)$$

つまり、表を対角線（右下）方向に1つ進むごとに行列の次数が2倍になり、型（$\mathbb{R},\mathbb{C},\mathbb{H}$の種別と直和の有無）は変わりません。系列によって得られた第1行と第1列を起点として、対角線方向に空白部分を埋めることができます。例えば：

1. $p=1, q=1$ の場合：
   $\operatorname{Cl}_{1,1}(\mathbb{R}) \cong \operatorname{Cl}_{0,0}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{R} \otimes \mathbb{R}(2) \cong \mathbb{R}(2)$

2. $p=1, q=2$ の場合：
   $\operatorname{Cl}_{1,2}(\mathbb{R}) \cong \operatorname{Cl}_{0,1}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{C} \otimes \mathbb{R}(2) \cong \mathbb{C}(2)$

3. $p=2, q=1$ の場合：
   $\operatorname{Cl}_{2,1}(\mathbb{R}) \cong \operatorname{Cl}_{1,0}(\mathbb{R}) \otimes \mathbb{H}' \cong 2\mathbb{R} \otimes \mathbb{R}(2) \cong 2\mathbb{R}(2)$

## 完全な分類表

このようにして得られる完全な分類表は以下のとおりです。[[wiki-clif]]

$$
\begin{array}{|c|ccccccccc|}
\hline
p \backslash q & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 \\
\hline
0 & \mathbb{R} & \mathbb{C} & \mathbb{H} & 2\mathbb{H} & \mathbb{H}(2) & \mathbb{C}(4) & \mathbb{R}(8) & 2\mathbb{R}(8) & \mathbb{R}(16) \\
1 & 2\mathbb{R} & \mathbb{R}(2) & \mathbb{C}(2) & \mathbb{H}(2) & 2\mathbb{H}(2) & \mathbb{H}(4) & \mathbb{C}(8) & \mathbb{R}(16) & 2\mathbb{R}(16) \\
2 & \mathbb{R}(2) & 2\mathbb{R}(2) & \mathbb{R}(4) & \mathbb{C}(4) & \mathbb{H}(4) & 2\mathbb{H}(4) & \mathbb{H}(8) & \mathbb{C}(16) & \mathbb{R}(32) \\
3 & \mathbb{C}(2) & \mathbb{R}(4) & 2\mathbb{R}(4) & \mathbb{R}(8) & \mathbb{C}(8) & \mathbb{H}(8) & 2\mathbb{H}(8) & \mathbb{H}(16) & \mathbb{C}(32) \\
4 & \mathbb{H}(2) & \mathbb{C}(4) & \mathbb{R}(8) & 2\mathbb{R}(8) & \mathbb{R}(16) & \mathbb{C}(16) & \mathbb{H}(16) & 2\mathbb{H}(16) & \mathbb{H}(32) \\
5 & 2\mathbb{H}(2) & \mathbb{H}(4) & \mathbb{C}(8) & \mathbb{R}(16) & 2\mathbb{R}(16) & \mathbb{R}(32) & \mathbb{C}(32) & \mathbb{H}(32) & 2\mathbb{H}(32) \\
6 & \mathbb{H}(4) & 2\mathbb{H}(4) & \mathbb{H}(8) & \mathbb{C}(16) & \mathbb{R}(32) & 2\mathbb{R}(32) & \mathbb{R}(64) & \mathbb{C}(64) & \mathbb{H}(64) \\
7 & \mathbb{C}(8) & \mathbb{H}(8) & 2\mathbb{H}(8) & \mathbb{H}(16) & \mathbb{C}(32) & \mathbb{R}(64) & 2\mathbb{R}(64) & \mathbb{R}(128) & \mathbb{C}(128) \\
8 & \mathbb{R}(16) & \mathbb{C}(16) & \mathbb{H}(16) & 2\mathbb{H}(16) & \mathbb{H}(32) & \mathbb{C}(64) & \mathbb{R}(128) & 2\mathbb{R}(128) & \mathbb{R}(256) \\
\hline
\end{array}
$$

## 8周期性

対角線方向に型が変わらないことから、型は$p-q$だけで決まります。さらに、公式1と公式2を続けて適用すると：

$$
\begin{alignedat}{3}
\operatorname{Cl}_{p,q}(\mathbb R) ⊗ \mathbb H ⊗ \mathbb H'
&\cong \operatorname{Cl}_{q,p+2}(\mathbb R) ⊗ \mathbb H'
&&\cong \operatorname{Cl}_{p+4,q}(\mathbb R) \\
\operatorname{Cl}_{p,q}(\mathbb R) ⊗ \mathbb H' ⊗ \mathbb H
&\cong \operatorname{Cl}_{q+2,p}(\mathbb R) ⊗ \mathbb H
&&\cong \operatorname{Cl}_{p,q+4}(\mathbb R)
\end{alignedat}
$$

$\mathbb H ⊗ \mathbb H' \cong \mathbb H(2)$より、$p$または$q$を$4$増やすことは$\mathbb H(2)$とのテンソル積に対応します：

$$
\operatorname{Cl}_{p+4,q}(\mathbb{R})
\cong \operatorname{Cl}_{p,q+4}(\mathbb{R})
\cong \operatorname{Cl}_{p,q}(\mathbb{R}) \otimes \mathbb{H}(2)
$$

これをもう一度適用すれば、$\mathbb H(2) ⊗ \mathbb H(2) \cong (\mathbb H ⊗ \mathbb H)(4) \cong \mathbb R(16)$より、周期$8$で同じ型に戻ります。

&&&thm 8周期性
$$
\operatorname{Cl}_{p+8,q}(\mathbb{R})
\cong \operatorname{Cl}_{p,q+8}(\mathbb{R})
\cong \operatorname{Cl}_{p,q}(\mathbb{R}) \otimes \mathbb{R}(16)
$$
&&&

したがって、型は$p-q \bmod 8$のみで決まり、行列の次数は全体の次元$2^{p+q}$から定まります。分類表から型を読み取ると、以下のようになります：

$$
\begin{array}{c|cccccccc}
p-q \bmod 8 & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 \\
\hline
\text{型} & \mathbb{R} & 2\mathbb{R} & \mathbb{R} & \mathbb{C} & \mathbb{H} & 2\mathbb{H} & \mathbb{H} & \mathbb{C}
\end{array}
$$

&&&rem
分類表に見られる規則性、すなわち対角線に沿った同一パターンの繰り返し、右下に向かう次数の指数的な増大、$2\mathbb{R},2\mathbb{H}$のような直和の規則的な出現は、すべて対角線方向の充填とこの8周期性から説明できます。

この周期性は実クリフォード代数の**ボット周期性**と呼ばれ、K理論におけるボット周期性定理と対応することが知られています。[[wiki-clif]]
&&&

# まとめ

本記事では、5つの基本構造から出発し、テンソル積による2つの拡張公式を繰り返し適用することで、クリフォード代数$\operatorname{Cl}_{p,q}(\mathbb R)$を行列環として分類しました。

&&&
基本構造：
$$
\begin{alignedat}{2}
&\mathbb{C}' &&\cong 2\mathbb{R} \\
&\mathbb{H}' &&\cong \mathbb{R}(2)
\end{alignedat}
$$
拡張公式：
$$
\begin{alignedat}{2}
\operatorname{Cl}_{p,q}(\mathbb R) &⊗ \mathbb H  &&\cong \operatorname{Cl}_{q,p+2}(\mathbb R) \\
\operatorname{Cl}_{p,q}(\mathbb R) &⊗ \mathbb H' &&\cong \operatorname{Cl}_{q+2,p}(\mathbb R) \cong \operatorname{Cl}_{p+1,q+1}(\mathbb R)
\end{alignedat}
$$
8周期性：
$$
\operatorname{Cl}_{p+8,q}(\mathbb{R})
\cong \operatorname{Cl}_{p,q+8}(\mathbb{R})
\cong \operatorname{Cl}_{p,q}(\mathbb{R}) ⊗ \mathbb{R}(16)
$$
&&&

完成した分類表では、型（$\mathbb{R},\mathbb{C},\mathbb{H}$の種別と直和の有無）は$p-q \bmod 8$のみで決まり、行列の次数は全体の次元$2^{p+q}$から定まります。
