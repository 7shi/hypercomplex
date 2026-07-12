実数体・複素数体・四元数体と分解型の代数を基本構造として、テンソル積による拡張公式を繰り返し適用することで、クリフォード代数$\operatorname{Cl}_{p,q}(\mathbb R)$を行列環として分類します。拡張公式は前回の記事で導出したものを使用します。[[qua-tensor]] 完成した分類表からは、擬スカラーによる型の判別、偶部分代数、複素化と2周期性、ピノル表現といった構造を読み取ります。

&&& 凡例
以下、代数$A$の直和$A \oplus A$を$2A$と略記し（例：$2\mathbb H \cong \mathbb H \oplus \mathbb H$）、$\mathbb F$を成分とする$n$次の全行列環$M_n(\mathbb F)$を$\mathbb F(n)$と表記します（例：$\mathbb R(2) \cong M_2(\mathbb R)$）。
&&&

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

これにより、任意の分解型複素数 $z = x + jy$ は以下のように一意的に分解できます。

$$z = (x-y)e + (x+y)e^*$$

この直和分解を用いると、2つの分解型複素数の積は非常に簡単に計算できます。

$$
(a_1e + b_1e^*)(a_2e + b_2e^*) = a_1a_2e + b_1b_2e^*
$$

この計算規則は成分ごとの積に対応しており、分解型複素数が実数の直和と同型であることを示しています。

$$\mathbb{C}' \cong \mathbb{R} \oplus \mathbb{R} \cong 2\mathbb{R}$$

&&&rem
代数$A$の直和$A \oplus A$を$2A$と略記します。
&&&

## 分解型四元数の行列表現

分解型四元数 $\mathbb{H}'$ は2次の実行列環と同型です。

$$\mathbb{H}' \cong \mathbb{R}(2)$$

&&&rem
$\mathbb F$を成分とする$n$次の全行列環$M_n(\mathbb F)$を$\mathbb F(n)$と表記します。
&&&

その生成元は以下の行列で表されます。

$$
i \cong \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \quad
j \cong \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
k = ij \cong \begin{pmatrix} -1 & 0 \\ 0 & 1 \end{pmatrix}
$$

この行列表現により、任意の分解型四元数 $q = a + bi + cj + dk$ は以下の2×2実行列として具体的に表現できます。

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

これらの規則は組み合わせて使用できます。例を示します。

$$\mathbb{H}(2) \otimes \mathbb{H} \cong (\mathbb{H} \otimes \mathbb{H})(2) \cong \mathbb{R}(4)(2) \cong \mathbb{R}(8)$$

# 系列による代数構造の導出

基本構造に2つの公式を交互に適用して、分類表の第1行（$p=0$）と第1列（$q=0$）を導出します。符号数$(p,0)$に公式1を適用すると$(0,p+2)$へ、$(0,q)$に公式2を適用すると$(q+2,0)$へ移るため、交互に適用することで$p$軸と$q$軸を往復しながら次数が上がっていきます。

## 複素数体からの系列

複素数体 $\mathbb{C} \cong \operatorname{Cl}_{0,1}(\mathbb{R})$ を出発点として、テンソル積による拡張を繰り返すと、以下の系列が得られます。

$$
\begin{aligned}
\operatorname{Cl}_{0,1}(\mathbb{R}) &\cong \mathbb{C} \\
\operatorname{Cl}_{3,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,1}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{C} \otimes \mathbb{R}(2) \cong \mathbb{C}(2) \\
\operatorname{Cl}_{0,5}(\mathbb{R}) &\cong \operatorname{Cl}_{3,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{C}(2) \otimes \mathbb{H} \cong \mathbb{C}(4) \\
\operatorname{Cl}_{7,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,5}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{C}(4) \otimes \mathbb{R}(2) \cong \mathbb{C}(8)
\end{aligned}
$$

## 分解型複素数からの系列

分解型複素数 $\mathbb{C}' \cong \operatorname{Cl}_{1,0}(\mathbb{R})$ からは、以下の系列が導かれます。

$$
\begin{aligned}
\operatorname{Cl}_{1,0}(\mathbb{R}) &\cong \mathbb{C}' \cong 2\mathbb{R} \\
\operatorname{Cl}_{0,3}(\mathbb{R}) &\cong \operatorname{Cl}_{1,0}(\mathbb{R}) \otimes \mathbb{H} \cong 2\mathbb{R} \otimes \mathbb{H} \cong 2\mathbb{H} \\
\operatorname{Cl}_{5,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,3}(\mathbb{R}) \otimes \mathbb{H}' \cong 2\mathbb{H} \otimes \mathbb{R}(2) \cong 2\mathbb{H}(2) \\
\operatorname{Cl}_{0,7}(\mathbb{R}) &\cong \operatorname{Cl}_{5,0}(\mathbb{R}) \otimes \mathbb{H} \cong 2\mathbb{H}(2) \otimes \mathbb{H} \cong 2\mathbb{R}(8)
\end{aligned}
$$

## 四元数体からの系列

四元数体 $\mathbb{H} \cong \operatorname{Cl}_{0,2}(\mathbb{R})$ からは、以下の系列が生成されます。

$$
\begin{aligned}
\operatorname{Cl}_{0,2}(\mathbb{R}) &\cong \mathbb{H} \\
\operatorname{Cl}_{4,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,2}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{H} \otimes \mathbb{R}(2) \cong \mathbb{H}(2) \\
\operatorname{Cl}_{0,6}(\mathbb{R}) &\cong \operatorname{Cl}_{4,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{H}(2) \otimes \mathbb{H} \cong \mathbb{R}(8) \\
\operatorname{Cl}_{8,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,6}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{R}(8) \otimes \mathbb{R}(2) \cong \mathbb{R}(16)
\end{aligned}
$$

## 分解型四元数からの系列

分解型四元数 $\mathbb{H}' \cong \operatorname{Cl}_{2,0}(\mathbb{R})$ からは、以下の系列が得られます。

$$
\begin{aligned}
\operatorname{Cl}_{2,0}(\mathbb{R}) &\cong \mathbb{H}' \cong \mathbb{R}(2) \\
\operatorname{Cl}_{0,4}(\mathbb{R}) &\cong \operatorname{Cl}_{2,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{R}(2) \otimes \mathbb{H} \cong \mathbb{H}(2) \\
\operatorname{Cl}_{6,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,4}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{H}(2) \otimes \mathbb{R}(2) \cong \mathbb{H}(4) \\
\operatorname{Cl}_{0,8}(\mathbb{R}) &\cong \operatorname{Cl}_{6,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{H}(4) \otimes \mathbb{H} \cong \mathbb{R}(16)
\end{aligned}
$$

## 系列の整理

これらの系列をまとめると、以下の表が得られます。

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

公式2は、$\mathbb H'$とのテンソル積によって符号数が$(p,q)$から$(p+1,q+1)$へ拡張できることを示しています。$\mathbb H' \cong \mathbb R(2)$なので、次の式が成り立ちます。

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

対角線方向に型が変わらないことから、型は$p-q$だけで決まります。さらに、公式1と公式2を続けて適用します。

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

$\mathbb H ⊗ \mathbb H' \cong \mathbb H(2)$より、$p$または$q$を$4$増やすことは$\mathbb H(2)$とのテンソル積に対応します。

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

したがって、型は$p-q \bmod 8$のみで決まり、行列の次数は全体の次元$2^{p+q}$から定まります。分類表から型を読み取ると、以下のようになります。

$$
\begin{array}{c|cccccccc}
p-q \bmod 8 & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 \\
\hline
\text{型} & \mathbb{R} & 2\mathbb{R} & \mathbb{R} & \mathbb{C} & \mathbb{H} & 2\mathbb{H} & \mathbb{H} & \mathbb{C}
\end{array}
$$

&&&rem
分類表に見られる規則性、すなわち対角線に沿った同一パターンの繰り返し、右下に向かう次数の指数的な増大、$2\mathbb{R},2\mathbb{H}$のような直和の規則的な出現は、すべて対角線方向の充填とこの8周期性から説明できます。直和や$\mathbb C$型がこの位置に現れる構造的な理由は、次節で擬スカラーを用いて説明します。

この周期性は実クリフォード代数の**ボット周期性**と呼ばれ、K理論におけるボット周期性定理と対応することが知られています。[[wiki-clif]]
&&&

# 擬スカラーと型の判別

分類表の型の並びには、$2\mathbb R,2\mathbb H$（直和型）が$p-q\equiv1\pmod4$に、$\mathbb C$型が$p-q\equiv3\pmod4$に現れるという規則性があります。これは、最高グレードの基底である**擬スカラー**

$$\omega = e_1e_2\cdots e_n \quad (n=p+q)$$

の性質、すなわち中心性（すべての元と可換かどうか）と$\omega^2$の符号から説明できます。

生成元$e_i$を$\omega$の中を通り抜けさせると、自分自身とは可換で、残りの$n-1$個の生成元とは反交換するため

$$e_i\,\omega=(-1)^{n-1}\omega\,e_i$$

となります。つまり$n$が奇数のとき、$\omega$はすべての生成元と可換、すなわち中心的です。一方、$\omega^2$は生成元を並べ替えて計算できます。並べ替えに必要な転倒数は$n(n-1)/2$なので

$$\omega^2=(-1)^{n(n-1)/2}\,e_1^2e_2^2\cdots e_n^2=(-1)^{n(n-1)/2+q}$$

となり、この符号は$p-q \bmod 4$だけで決まります。

$$
\omega^2=
\begin{cases}
+1 & (p-q\equiv0,1 \pmod 4)\\
-1 & (p-q\equiv2,3 \pmod 4)
\end{cases}
$$

&&&prf
符号数を$(p,q)\to(p+1,q+1)$と変えると、$n\to n+2$より指数の変化は$\{(n+2)(n+1)-n(n-1)\}/2+1=2n+2$で偶数となり、符号は変わらない。また$(p,q)\to(p+4,q)$では指数の変化は$\{(n+4)(n+3)-n(n-1)\}/2=4n+6$で偶数となり、やはり符号は変わらない。したがって符号は$p-q \bmod 4$のみに依存し、代表例$\operatorname{Cl}_{0,0},\operatorname{Cl}_{1,0},\operatorname{Cl}_{2,0},\operatorname{Cl}_{3,0}$での値$+1,+1,-1,-1$から上式が得られる。
&&&

$n$の偶奇と$\omega^2$の符号によって、中心（すべての元と可換な元のなす部分代数）が決まります：

- **$n$が偶数**：$\omega$は生成元と反交換するため中心に入らず、中心は$\mathbb R$のみです。型は$\mathbb R$か$\mathbb H$ですが、どちらになるか（$p-q \bmod 8$の区別）は中心だけでは決まりません。
- **$n$が奇数で$\omega^2=-1$**（$p-q\equiv3\pmod4$）：中心は$\mathbb R+\mathbb R\omega\cong\mathbb C$です。代数全体が複素行列環となり、$\mathbb C$型が現れます。
- **$n$が奇数で$\omega^2=+1$**（$p-q\equiv1\pmod4$）：中心は$\mathbb R+\mathbb R\omega\cong\mathbb C'\cong2\mathbb R$です。分解型複素数と同じ冪等元$\dfrac{1\pm\omega}2$によって代数全体が直和分解され、$2\mathbb R$型・$2\mathbb H$型が現れます。

型の表に中心を並べると、直和型と$\mathbb C$型の出現位置が中心の構造とちょうど一致していることが確認できます。

$$
\begin{array}{c|cccccccc}
p-q \bmod 8 & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 \\
\hline
\text{型} & \mathbb{R} & 2\mathbb{R} & \mathbb{R} & \mathbb{C} & \mathbb{H} & 2\mathbb{H} & \mathbb{H} & \mathbb{C} \\
\text{中心} & \mathbb R & 2\mathbb R & \mathbb R & \mathbb C & \mathbb R & 2\mathbb R & \mathbb R & \mathbb C
\end{array}
$$

&&&ex 冪等元による直和分解
最小の例は$\mathbb C'\cong\operatorname{Cl}_{1,0}(\mathbb R)$です。擬スカラー$\omega=j$は$\omega^2=+1$を満たし、冒頭で確認した冪等元$e=\dfrac{1-j}2,\ e^*=\dfrac{1+j}2$による直和分解$\mathbb C'\cong2\mathbb R$は、この仕組みの現れです。

次の例は$\operatorname{Cl}_{0,3}(\mathbb R)\cong2\mathbb H$（$p-q\equiv5\pmod8$）です。擬スカラー$\omega=e_1e_2e_3$は中心的で$\omega^2=+1$となり、冪等元$\dfrac{1\pm\omega}2$によって$\mathbb H\oplus\mathbb H$へ直和分解されます。
&&&

# 偶部分代数

偶数グレードの元は積について閉じており、**偶部分代数**$\operatorname{Cl}_{p,q}^0(\mathbb R)$をなします。次元は全体の半分の$2^{n-1}$です。偶部分代数は、生成元が1つ少ないクリフォード代数と同型になります。

&&&fml 偶部分代数の同型
$$
\operatorname{Cl}_{p,q}^0(\mathbb R)
\cong\operatorname{Cl}_{p,q-1}(\mathbb R)\ (q\ge1),\quad
\operatorname{Cl}_{p,q}^0(\mathbb R)
\cong\operatorname{Cl}_{q,p-1}(\mathbb R)\ (p\ge1)
$$
&&&

&&&prf
$e_n^2=-1$となる生成元$e_n$を1つ選び（$q\ge1$より存在する）、残りの生成元との積$f_i=e_ie_n$（$i=1,\dots,n-1$）を新たな生成元とする。$f_i$はグレード2なので偶部分代数に属し、2乗と反交換関係は

$$
f_i^2=e_ie_ne_ie_n=-e_i^2e_n^2=e_i^2,\quad
f_if_j=e_ie_ne_je_n=-e_ie_je_n^2=e_ie_j=-f_jf_i\ (i\ne j)
$$

となる。すなわち$f_i$は元の$e_i$と同じ2乗を持ち互いに反交換するため、符号数$(p,q-1)$のクリフォード代数を生成する。$f_i$の偶数個の積は$e_n$を含まない偶グレードの基底、奇数個の積は$e_n$を含む偶グレードの基底を（符号を除いて）与えるため、生成される代数は偶部分代数全体と一致し、次元$2^{n-1}$も一致する。

$e_n^2=+1$となる生成元を選んだ場合（$p\ge1$）は、$f_i^2=-e_i^2$とすべての2乗の符号が反転するため、符号数$(q,p-1)$のクリフォード代数が得られる。
&&&

2つの公式を組み合わせると、分類表の対称性が得られます。

$$\operatorname{Cl}_{p,q}(\mathbb R)\cong\operatorname{Cl}_{q+1,p-1}(\mathbb R)\quad(p\ge1)$$

例えば$\operatorname{Cl}_{4,0}(\mathbb R)\cong\operatorname{Cl}_{1,3}(\mathbb R)\cong\mathbb H(2)$です。また、符号数を反転した$\operatorname{Cl}_{p,q}(\mathbb R)$と$\operatorname{Cl}_{q,p}(\mathbb R)$は一般には同型ではありませんが（例：$\mathbb H'\not\cong\mathbb H$）、偶部分代数は共通です。

$$\operatorname{Cl}_{p,q}^0(\mathbb R)\cong\operatorname{Cl}_{q,p}^0(\mathbb R)$$

偶部分代数の中身は分類表から読み取ることができます。例を示します。

$$
\begin{alignedat}{2}
\operatorname{Cl}_{3,0}^0(\mathbb R)&\cong\operatorname{Cl}_{0,2}(\mathbb R)&&\cong\mathbb H \\
\operatorname{Cl}_{4,0}^0(\mathbb R)&\cong\operatorname{Cl}_{0,3}(\mathbb R)&&\cong2\mathbb H
\end{alignedat}
$$

&&&rem
回転群の二重被覆であるスピン群は偶部分代数の中に住んでおり、この同型は$\operatorname{Spin}(3)\cong\mathrm{SU}(2)$（単位四元数）や$\operatorname{Spin}(4)\cong\mathrm{SU}(2)\times\mathrm{SU}(2)$（単位四元数の対）の置き場所を与えます。
&&&

# 複素化と2周期性

係数を複素数体に拡張すると、分類は大幅に単純になります。複素数体上では生成元$e$を$ie$に置き換えることで$e^2$の符号を反転できるため、計量の符号の区別が消え、複素クリフォード代数$\operatorname{Cl}_n(\mathbb C)$は生成元の個数$n$だけで決まります。

$$\operatorname{Cl}_{p,q}(\mathbb R)\otimes\mathbb C\cong\operatorname{Cl}_{p+q}(\mathbb C)$$

複素化の計算規則は以下のとおりです。

$$
\mathbb R\otimes\mathbb C\cong\mathbb C,\quad
\mathbb C\otimes\mathbb C\cong2\mathbb C,\quad
\mathbb H\otimes\mathbb C\cong\mathbb C(2)
$$

&&&rem
3つ目は双四元数として確認済みです。[[pauli-qua]] 2つ目（双複素数）では、$u=i\otimes i$が$u^2=1$を満たすため、分解型複素数と同じ冪等元$\dfrac{1\pm u}2$による直和分解が起こります。
&&&

拡張公式$\operatorname{Cl}_{p,q}(\mathbb R)\otimes\mathbb H\cong\operatorname{Cl}_{q,p+2}(\mathbb R)$の両辺を複素化すると、$\mathbb H\otimes\mathbb C\cong\mathbb C(2)$より次の周期性が得られます。

&&&thm 2周期性
$$
\operatorname{Cl}_{n+2}(\mathbb C)
\cong\operatorname{Cl}_n(\mathbb C)\otimes\mathbb C(2)
$$
&&&

$\operatorname{Cl}_0(\mathbb C)\cong\mathbb C$と$\operatorname{Cl}_1(\mathbb C)\cong\mathbb C\otimes\mathbb C\cong2\mathbb C$を出発点として、複素の分類表が完成します。

$$
\begin{array}{c|ccccccc}
n & 0 & 1 & 2 & 3 & 4 & 5 & \cdots \\
\hline
\operatorname{Cl}_n(\mathbb C) & \mathbb C & 2\mathbb C & \mathbb C(2) & 2\mathbb C(2) & \mathbb C(4) & 2\mathbb C(4) & \cdots
\end{array}
$$

実の場合の8周期性に対して、複素では周期2で同じ型に戻り、型は$n \bmod 2$のみで決まります。

&&&rem
擬スカラーによる説明もそのまま通用します。複素数体上では$\omega^2=-1$であっても$i\omega$と取り直せば2乗が$+1$になるため、$n$が奇数なら必ず中心が$2\mathbb C$型となり、直和分解が起こります。実の8周期性が実K理論、複素の2周期性が複素K理論のボット周期性に対応します。[[wiki-clif]]
&&&

# ピノルとスピノル

分類表の行列環$\mathbb F(n)$は、$n$次元の列ベクトル空間$\mathbb F^n$に作用します。この作用に自明でない不変部分空間はなく、クリフォード代数の既約表現を与えます。

- 単純型$\mathbb F(n)$：既約表現は$\mathbb F^n$のただ1つです。
- 直和型$2\mathbb F(n)$：2つの直和成分がそれぞれ$\mathbb F^n$に作用するため、既約表現は2つあります。両者は中心的な擬スカラーが$\omega=+1,-1$のどちらとして作用するか（冪等元$\dfrac{1\pm\omega}2$のどちらの成分か）で区別されます。

&&&rem ピノルとスピノル
クリフォード代数全体の既約表現の空間の元を**ピノル**、偶部分代数の既約表現の空間の元を**スピノル**と呼び分けます。名称は、鏡映を含む直交群$\mathrm O(p,q)$の二重被覆である$\operatorname{Pin}(p,q)$群と、回転群$\mathrm{SO}(p,q)$の二重被覆である$\operatorname{Spin}(p,q)$群に由来します。
&&&

&&&rem 用語の慣習
物理の文献では、この区別を緩めて、クリフォード代数全体の既約表現の元も慣例的に**スピノル**と呼ぶことがよくあります。ピノル空間を偶部分代数に制限すると、スピノル表現そのもの、または2つのスピノル表現の直和が得られるため、両者は同じ空間をどちらの立場から見るかの違いです。
&&&

例えば$\operatorname{Cl}_{3,0}(\mathbb R)\cong\mathbb C(2)$のピノル空間は$\mathbb C^2$です。これはパウリ行列が作用する空間であり、偶部分代数$\operatorname{Cl}_{3,0}^0(\mathbb R)\cong\mathbb H$に制限すれば、$\mathrm{SU}(2)$（単位四元数）が作用する2成分スピノル（パウリスピノル）となります。[[pauli-qua]]

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
偶部分代数：
$$
\operatorname{Cl}_{p,q}^0(\mathbb R)
\cong\operatorname{Cl}_{p,q-1}(\mathbb R)
\cong\operatorname{Cl}_{q,p-1}(\mathbb R)
\quad(p,q\ge1)
$$
複素化と2周期性：
$$
\operatorname{Cl}_{p,q}(\mathbb R)\otimes\mathbb C\cong\operatorname{Cl}_{p+q}(\mathbb C),\quad
\operatorname{Cl}_{n+2}(\mathbb C)\cong\operatorname{Cl}_n(\mathbb C)\otimes\mathbb C(2)
$$
&&&

完成した分類表では、型（$\mathbb{R},\mathbb{C},\mathbb{H}$の種別と直和の有無）は$p-q \bmod 8$のみで決まり、行列の次数は全体の次元$2^{p+q}$から定まります。直和型と$\mathbb C$型の出現位置は擬スカラー$\omega$の中心性と$\omega^2$の符号に対応します。さらに、行列環が作用する列ベクトル空間はピノル・スピノルとして、スピン群を通じて回転群の表現論につながります。
