シリーズ：[四元数の行列表現](https://mathlog.info/series/PXPuUuQLYZk6HHho9eP8)

四元数の実行列表現は、演算規則から直接構成することができます。作用を受けるベクトル表現の要素の並び順を調整することで、複素行列表現に変換できることを見ていきます。

# はじめに

四元数$q = a + bi + cj + dk$は、実数$a, b, c, d$と虚数単位$i, j, k$を用いて表される数体系です。これらの虚数単位は$i^2 = j^2 = k^2 = ijk = -1$という特徴的な乗算規則を持ちます。四元数の演算、特に乗算は複雑に見えますが、線形代数の行列を用いることで、より直感的に理解することができます。

本記事では、四元数を4次元の実ベクトルとして表現して、その乗算（特に虚数単位$i, j, k$の作用）を行列で表します。ベクトルの要素の並び順は任意性がありますが、そのうち特定の並び順を選ぶことで、$4×4$の実行列が$2×2$の複素行列へと変換できることを示します。その特別な場合として、得られる複素行列表現がパウリ行列と簡潔な形で対応付けられるケースを調べます。

# 四元数の4次元実ベクトル表現と行列作用

四元数を、その係数を用いて4次元の実ベクトルに対応付けます。

&&&def 四元数と実ベクトル表現
$$
a + bi + cj + dk \mapsto \begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}
$$
&&&

実ベクトル表現に対して、左から虚数単位$i, j, k$を掛ける操作は、線形変換として$4×4$の実行列で表すことができます。導出に必要な演算規則を示します。

&&&def 四元数の乗算規則
$$
i^2 = j^2 = k^2 = -1
$$
$$
ij = -ji = k, \quad
jk = -kj = i, \quad
ki = -ik = j
$$
&&&

## 左からの$i$の作用

$$
\begin{aligned}
i(a + bi + cj + dk)
&= ai + bi^2 + cij + dik \\
&= -b + ai - dj + ck \\
&\mapsto \begin{pmatrix} -b \\ a \\ -d \\ c \end{pmatrix}
\end{aligned}
$$

この作用を行列$M_i$で表せば、以下のようになります。
$$
M_i
\begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}
= \begin{pmatrix}
0 & -1 & 0 & 0 \\
1 & 0 & 0 & 0 \\
0 & 0 & 0 & -1 \\
0 & 0 & 1 & 0
\end{pmatrix}
\begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}
=
\begin{pmatrix} -b \\ a \\ -d \\ c \end{pmatrix}
$$

## 左からの$j$の作用

$$
\begin{aligned}
j(a + bi + cj + dk)
&= aj + bji + cj^2 + djk \\
&= -c + di + aj - bk \\
&\mapsto \begin{pmatrix} -c \\ d \\ a \\ -b \end{pmatrix}
\end{aligned}
$$
この作用を行列$M_j$で表せば、以下のようになります。

$$
M_j
\begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}
= \begin{pmatrix}
0 & 0 & -1 & 0 \\
0 & 0 & 0 & 1 \\
1 & 0 & 0 & 0 \\
0 & -1 & 0 & 0
\end{pmatrix}
\begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}
=
\begin{pmatrix} -c \\ d \\ a \\ -b \end{pmatrix}
$$

## 左からの$k$の作用

$$
\begin{aligned}
k(a + bi + cj + dk)
&= ak + bki + ckj + dk^2 \\
&= -d - ci + bj + ak \\
&\mapsto \begin{pmatrix} -d \\ -c \\ b \\ a \end{pmatrix}
\end{aligned}
$$
この作用を行列$M_k$で表せば、以下のようになります。

$$
M_k
\begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}
= \begin{pmatrix}
0 & 0 & 0 & -1 \\
0 & 0 & -1 & 0 \\
0 & 1 & 0 & 0 \\
1 & 0 & 0 & 0
\end{pmatrix}
\begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}
=
\begin{pmatrix} -d \\ -c \\ b \\ a \end{pmatrix}
$$

# 複素行列表現とベクトルの要素の並べ替え

ベクトルに対する作用として行列を構成したため、ベクトル表現の要素を並べ替えると、それぞれ異なる$4×4$実行列表現を与えます。そのうち$2×2$複素行列表現に変換できる場合があります。

$4×4$の実行列を$2×2$のブロックに区切ります。各ブロックが以下の複素数の行列表現のパターンに一致すれば、$4×4$の実行列は$2×2$の複素行列に変換できます。

$$
\alpha + i\beta \cong
\begin{pmatrix} \alpha & -\beta \\ \beta & \alpha \end{pmatrix}
\quad(α,β\in\mathbb{R})
$$

しかし、先ほど求めた$M_i,M_j,M_k$のうち、複素行列に変換できるのは$M_i$のみで、$M_j,M_k$は条件を満たしません。

ベクトル表現の要素の並び順は$4! = 24$通りです。総当たりで確認したところ、そのうち半分の12通りで、構成した実行列表現が複素行列表現に変換できることが分かりました。[[7shi-qrc-colab]]

- (a,b,d,c), (a,c,b,d), (a,d,c,b), (b,a,c,d), (b,c,d,a), (b,d,a,c), (c,a,d,b), (c,b,a,d), (c,d,b,a), (d,a,b,c), (d,b,c,a), (d,c,a,b)

このうち特徴的な結果となる組み合わせを1つ選んで、複素行列への変換の例を示します。

## ベクトル要素の並べ替え

ベクトル表現の要素を$(d,a,b,c)$の順に並べ変えます。

$$
\begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}
→ \begin{pmatrix} d \\ a \\ b \\ c \end{pmatrix}
$$

$(a,b,c,d)$の要素を1つずつずらして循環させた形であることに注目すれば、$i,j,k$を左から作用させたベクトルが構成できます。

$$
i\times: \begin{pmatrix} -b \\ a \\ -d \\ c \end{pmatrix}
→ \begin{pmatrix} c \\ -b \\ a \\ -d \end{pmatrix},\quad
j\times: \begin{pmatrix} -c \\ d \\ a \\ -b \end{pmatrix}
→ \begin{pmatrix} -b \\ -c \\ d \\ a \end{pmatrix},\quad
k\times: \begin{pmatrix} -d \\ -c \\ b \\ a \end{pmatrix}
→ \begin{pmatrix} a \\ -d \\ -c \\ b \end{pmatrix}
$$

この結果を再現するように、作用の行列表現$M_i',M_j',M_k'$を求めます。

$$
\begin{alignedat}{2}
M_i'
\begin{pmatrix} d \\ a \\ b \\ c \end{pmatrix}
&= \begin{pmatrix}
0 & 0 & 0 & 1 \\
0 & 0 & -1 & 0 \\
0 & 1 & 0 & 0 \\
-1 & 0 & 0 & 0
\end{pmatrix}
\begin{pmatrix} d \\ a \\ b \\ c \end{pmatrix}
&&=
\begin{pmatrix} c \\ -b \\ a \\ -d \end{pmatrix}
\\
M_j'
\begin{pmatrix} d \\ a \\ b \\ c \end{pmatrix}
&= \begin{pmatrix}
0 & 0 & -1 & 0 \\
0 & 0 & 0 & -1 \\
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0
\end{pmatrix}
\begin{pmatrix} d \\ a \\ b \\ c \end{pmatrix}
&&=
\begin{pmatrix} -b \\ -c \\ d \\ a \end{pmatrix}
\\
M_k'
\begin{pmatrix} d \\ a \\ b \\ c \end{pmatrix}
&= \begin{pmatrix}
0 & 1 & 0 & 0 \\
-1 & 0 & 0 & 0 \\
0 & 0 & 0 & -1 \\
0 & 0 & 1 & 0
\end{pmatrix}
\begin{pmatrix} d \\ a \\ b \\ c \end{pmatrix}
&&=
\begin{pmatrix} a \\ -d \\ -c \\ b \end{pmatrix}
\end{alignedat}
$$

## 複素行列への変換

これらの$4×4$実行列$M_i', M_j', M_k'$を$2×2$のブロックに分け、$\begin{pmatrix} \alpha & -\beta \\ \beta & \alpha \end{pmatrix} → \alpha + i\beta$の規則で変換して、$I_H,J_H,K_H$とします。（添え字は四元数全体の集合を表す$\mathbb{H}$に由来）

$$
\begin{alignedat}{3}
M_i' &= \left(
  \begin{array}{cc|cc}
    0 & 0 & 0  & 1 \\
    0 & 0 & -1 & 0 \\ \hline
    0 & 1 & 0  & 0 \\
    -1 & 0 & 0  & 0
  \end{array}
\right)
&&\ → \ \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix} &&=: I_H
\\
M_j' &= \left(
    \begin{array}{cc|cc}
        0 & 0 & -1 & 0 \\
        0 & 0 & 0 & -1 \\ \hline
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0
    \end{array}
\right)
&&\ → \ \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} &&=: J_H
\\
M_k' &= \left(
    \begin{array}{cc|cc}
        0 & 1 & 0 & 0 \\
        -1 & 0 & 0 & 0 \\ \hline
        0 & 0 & 0 & -1 \\
        0 & 0 & 1 & 0
    \end{array}
\right)
&&\ → \ \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix} &&=: K_H
\end{alignedat}
$$

行列の作用を受けるベクトルは2成分ずつ複素数に変換します。

$$
\begin{pmatrix} d \\ a \\ b \\ c \end{pmatrix}
→ \begin{pmatrix} d+ia \\ b+ic \end{pmatrix}
$$

この形式で、複素ベクトルとして期待される結果を確認します。

$$
i\times: \begin{pmatrix} c \\ -b \\ a \\ -d \end{pmatrix}
→ \begin{pmatrix} c-ib \\ a-id \end{pmatrix},\quad
j\times: \begin{pmatrix} -b \\ -c \\ d \\ a \end{pmatrix}
→ \begin{pmatrix} -b-ic \\ d+ia \end{pmatrix},\quad
k\times: \begin{pmatrix} a \\ -d \\ -c \\ b \end{pmatrix}
→ \begin{pmatrix} a-id \\ -c+ib \end{pmatrix}
$$

$I_H,J_H,K_H$の作用を計算すれば、期待される結果と一致します。

$$
\begin{alignedat}{3}
I_H
\begin{pmatrix} d+ia \\ b+ic \end{pmatrix}
&= &\begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}
&\begin{pmatrix} d+ia \\ b+ic \end{pmatrix}
&&=
\begin{pmatrix} c-ib \\ a-id \end{pmatrix}
\\
J_H
\begin{pmatrix} d+ia \\ b+ic \end{pmatrix}
&= &\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
&\begin{pmatrix} d+ia \\ b+ic \end{pmatrix}
&&=
\begin{pmatrix} -b-ic \\ d+ia \end{pmatrix}
\\
K_H
\begin{pmatrix} d+ia \\ b+ic \end{pmatrix}
&= &\begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
&\begin{pmatrix} d+ia \\ b+ic \end{pmatrix}
&&=
\begin{pmatrix} a-id \\ -c+ib \end{pmatrix}
\end{alignedat}
$$

## 複素行列表現の性質

得られた$2×2$複素行列は、元の四元数の虚数単位$i, j, k$の代数構造を保持しています。例えば、四元数の積$ij = k$に対応する行列の積を確認します。

$$
I_H J_H
=\begin{pmatrix}  0 & -i \\ -i & 0 \end{pmatrix}
 \begin{pmatrix}  0 & -1 \\  1 & 0 \end{pmatrix}
=\begin{pmatrix} -i &  0 \\  0 & i \end{pmatrix}
=K_H
$$

期待される結果を再現しました。同様に、他の積（例：$J_H I_H = -K_H$）や2乗（例：${I_H}^2 = {J_H}^2 = {K_H}^2 = -I$、ここで$I$は単位行列）も、四元数の規則と一致することが確認できます。

## ベクトル表現の任意性

行列の作用を受けるベクトルは、定数倍の任意性があります。（$c\in\mathbb{C},\ c\ne0$とする）

$$
M\mathbf v=\mathbf v' \iff M(c\,\mathbf v)=c\,\mathbf v'
$$

ここまで使用した2次元複素ベクトル$\begin{pmatrix} d+ia \\ b+ic \end{pmatrix}$を$-i$倍します。

$$
-i\begin{pmatrix} d+ia \\ b+ic \end{pmatrix}
=\begin{pmatrix} a-id \\ c-ib \end{pmatrix}
$$

このベクトル表現は行列表現の線形結合の第1列と一致するため、対応が付けやすいです。

&&&fml 四元数の複素行列表現の線形結合
$$
\begin{aligned}
aI+bI_H+cJ_H+dK_H
&= a\begin{pmatrix}  1 &  0 \\  0 & 1 \end{pmatrix}
+ b\begin{pmatrix}  0 & -i \\ -i & 0 \end{pmatrix}
+ c\begin{pmatrix}  0 & -1 \\  1 & 0 \end{pmatrix}
+ d\begin{pmatrix} -i &  0 \\  0 & i \end{pmatrix} \\
&=\begin{pmatrix}a-id & -c-ib \\ c-ib & a+id \end{pmatrix}
\end{aligned}
$$
&&&

## パウリ行列との対応

得られた$2×2$複素行列$I_H,J_H,K_H$に$i$を掛ければ、パウリ行列$\sigma_1, \sigma_2, \sigma_3$が得られます。[[7shi-qp]]

$$
\begin{alignedat}{3}
i I_H
&= i \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}
&&= \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
&&= \sigma_1
\\
i J_H
&= i \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
&&= \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}
&&= \sigma_2
\\
i K_H
&= i \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
&&= \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
&&= \sigma_3
\end{alignedat}
$$

このように、四元数の虚数単位$i, j, k$は、$(d,a,b,c)$の並び順から構成した複素行列表現を通じて、パウリ行列と対応付けられます。

# まとめ

四元数$q = a + bi + cj + dk$の演算は、4次元の実ベクトル空間における線形変換として捉えることができます。虚数単位$i, j, k$の左からの乗算は、特定の$4×4$実行列$M_i, M_j, M_k$に対応します。

&&&def 四元数の実行列表現 $(a,b,c,d)$
$$
M_i = \begin{pmatrix}
0 & -1 & 0 & 0 \\
1 & 0 & 0 & 0 \\
0 & 0 & 0 & -1 \\
0 & 0 & 1 & 0
\end{pmatrix},\quad
M_j = \begin{pmatrix}
0 & 0 & -1 & 0 \\
0 & 0 & 0 & 1 \\
1 & 0 & 0 & 0 \\
0 & -1 & 0 & 0
\end{pmatrix},\quad
M_k = \begin{pmatrix}
0 & 0 & 0 & -1 \\
0 & 0 & -1 & 0 \\
0 & 1 & 0 & 0 \\
1 & 0 & 0 & 0
\end{pmatrix}
$$
&&&

ベクトルの要素を$(d,a,b,c)$という並び順（四元数の係数を$k,1,i,j$の順でベクトル化）に変更して、対応する$4×4$実行列$M'_i, M'_j, M'_kを2×2$のブロックに分解します。

&&&def 四元数の実行列表現 $(d,a,b,c)$
$$
M_i' = \left(
  \begin{array}{cc|cc}
    0 & 0 & 0  & 1 \\
    0 & 0 & -1 & 0 \\ \hline
    0 & 1 & 0  & 0 \\
    -1 & 0 & 0  & 0
  \end{array}
\right),\quad
M_j' = \left(
    \begin{array}{cc|cc}
        0 & 0 & -1 & 0 \\
        0 & 0 & 0 & -1 \\ \hline
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0
    \end{array}
\right),\quad
M_k' = \left(
    \begin{array}{cc|cc}
        0 & 1 & 0 & 0 \\
        -1 & 0 & 0 & 0 \\ \hline
        0 & 0 & 0 & -1 \\
        0 & 0 & 1 & 0
    \end{array}
\right)
$$
&&&

各ブロックは複素数の行列表現に合致します。

&&&def 複素数の行列表現
$$
\alpha + i\beta \cong
  \alpha\begin{pmatrix}1 &  0 \\ 0 & 1\end{pmatrix}
+ \beta \begin{pmatrix}0 & -1 \\ 1 & 0\end{pmatrix}
= \begin{pmatrix} \alpha & -\beta \\ \beta & \alpha \end{pmatrix}
\quad(α,β\in\mathbb{R})
$$
&&&

各ブロックを複素数に変換すれば、$2×2$の複素行列$I_H,J_H,K_H$が得られます。

&&&def 四元数の複素行列表現
$$
I_H = \begin{pmatrix}  0 & -i \\ -i & 0 \end{pmatrix}, \quad
J_H = \begin{pmatrix}  0 & -1 \\  1 & 0 \end{pmatrix}, \quad
K_H = \begin{pmatrix} -i &  0 \\  0 & i \end{pmatrix}
$$
&&&

この$2×2$複素行列表現は、四元数の代数構造（乗算規則、非可換性など）を忠実に再現します。これにより、四元数の理論を複素行列の枠組みで扱うことが可能になり、コンピュータグラフィックス（3次元回転）など、様々な分野での応用につながります。

この$(d,a,b,c)$という並び順は、それによって得られた$I_H,J_H,K_H$に$i$を掛けることで、パウリ行列が構成できるように選ばれています。

&&&def 四元数の行列表現とパウリ行列の関係
$$
\sigma_1 = iI_H = \begin{pmatrix} 0 &  1 \\ 1 &  0 \end{pmatrix},\quad
\sigma_2 = iJ_H = \begin{pmatrix} 0 & -i \\ i &  0 \end{pmatrix},\quad
\sigma_3 = iK_H = \begin{pmatrix} 1 &  0 \\ 0 & -1 \end{pmatrix}
$$
両辺に$i^{-1}=-i$を掛けた逆の関係：
$$
I_H=-i\sigma_1,\quad J_H=-i\sigma_2,\quad K_H=-i\sigma_3
$$
&&&

この$-i\sigma_i\ (i=1,2,3)$という形は、四元数が量子力学におけるスピンや、ローレンツ群の表現論など、物理学の基本的な記述と密接に関連していることを示唆しています。
