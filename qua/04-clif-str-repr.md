クリフォード代数をテンソル積で拡張して、行列表現の表を作成します。

# 基本的な代数構造

クリフォード代数には、5つの基本的な代数構造があります。

$$
\begin{alignedat}{3}
&\mathbb{R} &&\cong \operatorname{Cl}_{0,0}(\mathbb{R}) &&\quad\text{（実数体）} \\
&\mathbb{C} &&\cong \operatorname{Cl}_{0,1}(\mathbb{R}) &&\quad\text{（複素数体）} \\
&\mathbb{C}' &&\cong \operatorname{Cl}_{1,0}(\mathbb{R}) &&\quad\text{（分解型複素数体）} \\
&\mathbb{H} &&\cong \operatorname{Cl}_{0,2}(\mathbb{R}) &&\quad\text{（四元数体）} \\
&\mathbb{H}' &&\cong \operatorname{Cl}_{2,0}(\mathbb{R}) \cong \operatorname{Cl}_{1,1}(\mathbb{R}) &&\quad\text{（分解型四元数体）}
\end{alignedat}
$$

これらの基本構造は、より高次のクリフォード代数を構成する際の基礎となります。



## 1.2 分解型代数の特徴

### 分解型複素数の構造

分解型複素数 $\mathbb{C}'$ は、通常の複素数体とは異なる性質を持ちます。その最も基本的な特徴は、虚数単位 $j$ が $j^2 = +1$ を満たすことです（通常の虚数単位 $i$ は $i^2 = -1$ を満たす）。

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

$\mathbb{C}' \cong \mathbb{R} \oplus \mathbb{R} \cong 2\mathbb{R}$

### 分解型四元数の表現

分解型四元数 $\mathbb{H}'$ は2×2実行列環と同型です：

$$\mathbb{H}' \cong \mathbb{R}(2) \cong \operatorname{Cl}_{2,0}(\mathbb{R}) \cong \operatorname{Cl}_{1,1}(\mathbb{R})$$

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

# 2. テンソル積による拡張

## 2.1 基本的な拡張公式

クリフォード代数は、テンソル積を用いて以下の公式により拡張することができます：

$$
\begin{alignedat}{2}
\operatorname{Cl}_{p,q}(\mathbb{R}) &\otimes \mathbb{H}  &&\cong \operatorname{Cl}_{q,p+2}(\mathbb{R}) \\
\operatorname{Cl}_{p,q}(\mathbb{R}) &\otimes \mathbb{H}' &&\cong \operatorname{Cl}_{q+2,p}(\mathbb{R}) \cong \operatorname{Cl}_{p+1,q+1}(\mathbb{R})
\end{alignedat}
$$

これらの公式は、基本的な代数構造から出発して、より高次の代数を系統的に構成する方法を与えます。

## 2.2 テンソル積の性質

テンソル積による拡張では、以下の性質が重要な役割を果たします：

1. 行列環との関係：

$$\mathbb{F}(n) \otimes \mathbb{R}(2) \cong \mathbb{F}(2n)\quad(\mathbb{F}\in\{\mathbb{R},\mathbb{C},\mathbb{H}\})$$

この関係は、行列表現の次元が拡張によって2倍になることを示しています。

2. 四元数との関係：

$$
\begin{aligned}
\mathbb{H}' \otimes \mathbb{H} &\cong \mathbb{R}(2) \otimes \mathbb{H} \cong \mathbb{H}(2) \\
\mathbb{C} \otimes \mathbb{H} &\cong \mathbb{C} \otimes \mathbb{H}' \cong \mathbb{C}(2)
\end{aligned}
$$

## 2.3 同型対応の法則

テンソル積による拡張において、以下の3つの重要な同型対応が成り立ちます：

$$
\begin{aligned}
A \otimes B &\cong B \otimes A \\
\mathbb{H} \otimes \mathbb{H} &\cong \mathbb{H}' \otimes \mathbb{H}' \\
\mathbb{C} \otimes \mathbb{H} &\cong \mathbb{C} \otimes \mathbb{H}'
\end{aligned}
$$

これらの同型関係は、テンソル積の表現を得るために使用します。

# 3. 系列による代数構造の導出

## 3.1 複素数体からの系列

複素数体 $\mathbb{C} \cong \operatorname{Cl}_{0,1}(\mathbb{R})$ を出発点として、テンソル積による拡張を繰り返すと、以下の系列が得られます：

$$
\begin{aligned}
\operatorname{Cl}_{0,1}(\mathbb{R}) &\cong \mathbb{C} \\
\operatorname{Cl}_{3,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,1}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{C} \otimes \mathbb{R}(2) \cong \mathbb{C}(2) \\
\operatorname{Cl}_{0,5}(\mathbb{R}) &\cong \operatorname{Cl}_{3,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{C}(2) \otimes \mathbb{H} \cong \mathbb{C}(4) \\
\operatorname{Cl}_{7,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,5}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{C}(4) \otimes \mathbb{R}(2) \cong \mathbb{C}(8)
\end{aligned}
$$

## 3.2 分解型複素数体からの系列

分解型複素数体 $\mathbb{C}' \cong \operatorname{Cl}_{1,0}(\mathbb{R})$ からは、以下の系列が導かれます：

$$
\begin{aligned}
\operatorname{Cl}_{1,0}(\mathbb{R}) &\cong \mathbb{C}' \cong 2\mathbb{R} \\
\operatorname{Cl}_{0,3}(\mathbb{R}) &\cong \operatorname{Cl}_{1,0}(\mathbb{R}) \otimes \mathbb{H} \cong 2\mathbb{R} \otimes \mathbb{H} \cong 2\mathbb{H} \\
\operatorname{Cl}_{5,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,3}(\mathbb{R}) \otimes \mathbb{H}' \cong 2\mathbb{H} \otimes \mathbb{R}(2) \cong 2\mathbb{H}(2) \\
\operatorname{Cl}_{0,7}(\mathbb{R}) &\cong \operatorname{Cl}_{5,0}(\mathbb{R}) \otimes \mathbb{H} \cong 2\mathbb{H}(2) \otimes \mathbb{H} \cong 2\mathbb{R}(8)
\end{aligned}
$$

## 3.3 四元数体からの系列

四元数体 $\mathbb{H} \cong \operatorname{Cl}_{0,2}(\mathbb{R})$ からは、以下の系列が生成されます：

$$
\begin{aligned}
\operatorname{Cl}_{0,2}(\mathbb{R}) &\cong \mathbb{H} \\
\operatorname{Cl}_{4,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,2}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{H} \otimes \mathbb{R}(2) \cong \mathbb{H}(2) \\
\operatorname{Cl}_{0,6}(\mathbb{R}) &\cong \operatorname{Cl}_{4,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{H}(2) \otimes \mathbb{H} \cong \mathbb{R}(8) \\
\operatorname{Cl}_{8,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,6}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{R}(8) \otimes \mathbb{R}(2) \cong \mathbb{R}(16)
\end{aligned}
$$

## 3.4 分解型四元数体からの系列

分解型四元数体 $\mathbb{H}' \cong \operatorname{Cl}_{2,0}(\mathbb{R})$ からは、以下の系列が得られます：

$$
\begin{aligned}
\operatorname{Cl}_{2,0}(\mathbb{R}) &\cong \mathbb{H}' \cong \mathbb{R}(2) \\
\operatorname{Cl}_{0,4}(\mathbb{R}) &\cong \operatorname{Cl}_{2,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{R}(2) \otimes \mathbb{H} \cong \mathbb{H}(2) \\
\operatorname{Cl}_{6,0}(\mathbb{R}) &\cong \operatorname{Cl}_{0,4}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{H}(2) \otimes \mathbb{R}(2) \cong \mathbb{H}(4) \\
\operatorname{Cl}_{0,8}(\mathbb{R}) &\cong \operatorname{Cl}_{6,0}(\mathbb{R}) \otimes \mathbb{H} \cong \mathbb{H}(4) \otimes \mathbb{H} \cong \mathbb{R}(16)
\end{aligned}
$$

## 3.5 系列の整理と表による分類

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

この表の空白部分は、次章で導入する$(p+1,q+1)$の公式を用いて埋めることができます。

# 4. テンソル積による完全分類

## 4.1 (p+1,q+1)の拡張公式

クリフォード代数は以下の公式によって拡張することができます：

$$\operatorname{Cl}_{p,q}(\mathbb{R}) \otimes \mathbb{H}' \cong \operatorname{Cl}_{p+1,q+1}(\mathbb{R})$$

この公式は、各クリフォード代数を分解型四元数とのテンソル積によって拡張できることを示しています。

## 4.2 空白部分の導出

前章で得られた表の空白部分は、この拡張公式を繰り返し適用することで埋めることができます。例えば：

1. $p=1, q=1$ の場合：
   $\operatorname{Cl}_{1,1}(\mathbb{R}) \cong \operatorname{Cl}_{0,0}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{R} \otimes \mathbb{R}(2) \cong \mathbb{R}(2)$

2. $p=1, q=2$ の場合：
   $\operatorname{Cl}_{1,2}(\mathbb{R}) \cong \operatorname{Cl}_{0,1}(\mathbb{R}) \otimes \mathbb{H}' \cong \mathbb{C} \otimes \mathbb{R}(2) \cong \mathbb{C}(2)$

3. $p=2, q=1$ の場合：
   $\operatorname{Cl}_{2,1}(\mathbb{R}) \cong \operatorname{Cl}_{1,0}(\mathbb{R}) \otimes \mathbb{H}' \cong 2\mathbb{R} \otimes \mathbb{R}(2) \cong 2\mathbb{R}(2)$

## 4.3 完全な分類表

このようにして得られる完全な分類表は以下のとおりです：

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

## 4.4 周期性と規則性の解析

この完全な分類表からは、以下の重要な規則性が読み取れます：

1. 対角線に沿って同様のパターンが現れます
2. 行列環の次元は右下に向かうにつれて指数的に大きくなります
3. $\mathbb{R}, \mathbb{C}, \mathbb{H}$ の順序パターンが繰り返し現れます
4. 2倍化（例：$2\mathbb{H}, 2\mathbb{R}(8)$）が規則的に現れます

これらの規則性は、クリフォード代数の構造を理解する上で重要な指針となります。
