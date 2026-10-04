四元数の$4×4$実行列による表現から、ベクトルの成分の並び順を調整することで2次複素行列表現およびパウリ行列を導出します。

シリーズ：[四元数の行列表現](https://mathlog.info/series/PXPuUuQLYZk6HHho9eP8)

&&& 改訂履歴
- 2026.10.04 ブロックの複素数化と一般の四元数の行列表現が積を保つことを示し、並び順の選び方と関連記事の表現との関係を補い、まとめの応用への言及を削除した
&&&

# 概要

四元数は、一般に$q = a + bi + cj + dk\ (a,b,c,d\in\mathbb R)$と表されます。虚数単位$i, j, k$は$i^2 = j^2 = k^2 = ijk = -1$という乗算規則を持ちます。四元数の左乗算は、係数ベクトルに作用する実線形変換として表せます。

本記事では、四元数を4次元の実ベクトルとして表現し、左からの虚数単位の作用を行列として表します。ベクトルの成分の並び順には任意性がありますが、特定の並び順を選ぶことで$4×4$の実行列が$2×2$の複素行列へと変換できることを示します。さらに、得られた行列を複素数の$i$倍すると、標準的なパウリ行列になることを確認します。

行列の積と複素数の基本演算を前提とします。複素数の実行列表現は本文中で確認します。回転や量子力学への応用は扱いません。

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

# 複素行列表現とベクトルの成分の並べ替え

ベクトルに対する作用として行列を構成したため、ベクトル表現の成分を並べ替えると、それぞれ異なる$4×4$実行列表現を与えます。そのうち$2×2$複素行列表現に変換できる場合があります。

複素数$α+iβ$を掛ける操作は、実部と虚部を並べた実ベクトルに対する$2×2$実行列の作用として表せます。

&&&def 複素数の実行列表現
$$
R(α+iβ)=\begin{pmatrix} α & -β \\ β & α \end{pmatrix}\quad(α,β\in\mathbb{R})
$$
&&&

ベクトル$(x,y)^\mathsf T$を複素数$x+iy$とみなすと、この行列の作用は複素数の積に対応します。

&&&fml 行列の作用と複素数の積
$$
\begin{pmatrix} α & -β \\ β & α \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix}
=\begin{pmatrix} αx-βy \\ βx+αy \end{pmatrix},\quad
(α+iβ)(x+iy)=(αx-βy)+i(βx+αy)
$$
&&&

$R$は加法と乗法を保ちます。

$$
R(z+w)=R(z)+R(w),\quad R(zw)=R(z)R(w)
$$

$4×4$の実行列を$2×2$のブロックに区切り、各ブロックがこの形になる場合を考えます。ブロック行列の積は各ブロックの積と和で計算されるため、各ブロックを$R$で複素数に置き換える操作も、行列の加法と乗法を保ちます。こうして、各ブロックがこの形になる実行列の部分代数は、$2×2$複素行列の代数と同型になります。以下、この同型による対応を$\cong$で表します。

しかし、先ほど求めた$M_i,M_j,M_k$のうち、複素行列に変換できるのは$M_i$のみで、$M_j,M_k$は条件を満たしません。

ベクトル表現の成分の並び順は$4! = 24$通りです。それぞれの並べ替えを$M_i,M_j,M_k$の3つに施し、すべての$2×2$ブロックが上の形になるかを総当たりで調べたところ、半分の12通りで複素行列表現に変換できることが分かりました。[[7shi-colab-qrc]]

- (a,b,d,c), (a,c,b,d), (a,d,c,b), (b,a,c,d), (b,c,d,a), (b,d,a,c), (c,a,d,b), (c,b,a,d), (c,d,b,a), (d,a,b,c), (d,b,c,a), (d,c,a,b)

このうち、得られる複素行列を複素数の$i$倍すると標準的なパウリ行列になる並び順$(d,a,b,c)$を選んで、複素行列への変換の例を示します。

## ベクトル成分の並べ替え

ベクトル表現の成分を$(d,a,b,c)$の順に並べ替えます。並べ替え行列$P$を用いると、次のように書けます。

$$
P\begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}
=\begin{pmatrix} d \\ a \\ b \\ c \end{pmatrix},\quad
P=\begin{pmatrix} 0&0&0&1 \\ 1&0&0&0 \\ 0&1&0&0 \\ 0&0&1&0 \end{pmatrix}
$$

並べ替えは、入力のベクトルだけでなく、左乗算した後の出力のベクトルにも同じように施します。出力の第4成分を先頭へ移し、残りの成分を順に1つずつ後ろへずらすと、$i,j,k$を左から作用させたベクトルが得られます。

$$
i\times: \begin{pmatrix} -b \\ a \\ -d \\ c \end{pmatrix}
→ \begin{pmatrix} c \\ -b \\ a \\ -d \end{pmatrix},\quad
j\times: \begin{pmatrix} -c \\ d \\ a \\ -b \end{pmatrix}
→ \begin{pmatrix} -b \\ -c \\ d \\ a \end{pmatrix},\quad
k\times: \begin{pmatrix} -d \\ -c \\ b \\ a \end{pmatrix}
→ \begin{pmatrix} a \\ -d \\ -c \\ b \end{pmatrix}
$$

この結果を再現するように、作用の行列表現$M_i',M_j',M_k'$を求めます。入力と出力の両方を$P$で並べ替えるので、$M_r'=PM_rP^{-1}\ (r=i,j,k)$となります。

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

以下、行列の成分やスカラー倍に現れる$i$は複素数の虚数単位です。四元数の虚数単位$i$の左乗算は、ここで得る$I_H$で表します。

これらの$4×4$実行列$M_i', M_j', M_k'$を$2×2$のブロックに分け、$\begin{pmatrix} α & -β \\ β & α \end{pmatrix} \mapsto α + iβ$の規則で変換して、$I_H,J_H,K_H$とします。添え字は四元数全体の集合を表す$\mathbb{H}$に由来します。

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

行列の作用を受けるベクトルには、定数倍の任意性があります。入力と出力のベクトルを同じ非零の複素数倍に変えても、作用を表す行列は変わりません。$λ\in\mathbb{C},\ λ\ne0$として、次が成り立ちます。

$$
M\mathbf v=\mathbf v' \iff M(λ\mathbf v)=λ\mathbf v'
$$

ここまで使用した2次元複素ベクトル$\begin{pmatrix} d+ia \\ b+ic \end{pmatrix}$では、四元数$1$（$a=1,\ b=c=d=0$）のベクトルが$(i,0)^\mathsf T$になります。これを$(1,0)^\mathsf T$にそろえるため、$-i$倍します。

$$
-i\begin{pmatrix} d+ia \\ b+ic \end{pmatrix}
=\begin{pmatrix} a-id \\ c-ib \end{pmatrix}
$$

このベクトル表現は、次の行列表現の線形結合の第1列と一致します。

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

一般の四元数を、この線形結合とベクトルで表します。

&&&def 四元数の複素行列表現とベクトル表現
$$
ρ(a+bi+cj+dk)=aI+bI_H+cJ_H+dK_H
$$
$$
w(a+bi+cj+dk)=\begin{pmatrix} a-id \\ c-ib \end{pmatrix}
$$
&&&

$I_H,J_H,K_H$はそれぞれ$i,j,k$の左乗算を$w$を通じて表しているので、実線形性から、任意の四元数$x,y$について次が成り立ちます。

$$
ρ(x)w(y)=w(xy)
$$

特に$w(1)=(1,0)^\mathsf T$なので$ρ(x)w(1)=w(x)$となり、$w(x)$は$ρ(x)$の第1列に一致します。

## 積の保存と単射性

行列$ρ(x)$を四元数$x$の表現と呼ぶには、一般の四元数の積が行列の積に対応している必要があります。

&&&prop 積の保存と単射性
任意の四元数$x,y$に対して、次が成り立ちます。
$$
ρ(xy)=ρ(x)ρ(y),\quad ρ(1)=I
$$
また、$ρ$は単射です。
&&&

&&&prf
任意の四元数$z$に対して、$ρ(x)w(y)=w(xy)$と四元数の結合法則より次が成り立つ。

$$
ρ(x)ρ(y)w(z)=ρ(x)w(yz)=w(x(yz))=w((xy)z)=ρ(xy)w(z)
$$

$z$の係数は任意の実数を取れるので、$w(z)$は$\mathbb C^2$全体を動く。したがって$ρ(x)ρ(y)=ρ(xy)$である。$ρ(1)=I$は定義から明らかである。また、$ρ(q)$の第1列$w(q)$から係数$a,b,c,d$が復元できるので、$ρ$は単射である。
&&&

&&&rem 関連記事の表現との関係
ベクトルの第1成分を共役にする対応$v(a+bi+cj+dk)=(a-bi,\ c+di)^\mathsf T$から複素行列表現を構成すると、虚数単位の左乗算の行列$M_I,M_J,M_K$が得られます。本記事の$w$はこれとは異なるベクトル表現ですが、行列は$I_H=-M_K,\ J_H=M_J,\ K_H=M_I$の関係にあり、その記事で最後にパウリ行列に合わせて選んだ表現と一致します。なお、本記事の小文字の添字の$M_i,M_j,M_k$は$4×4$実行列、その記事の大文字の添字の$M_I,M_J,M_K$は$2×2$複素行列です。[[7shi-qcm]]
&&&

## パウリ行列との対応

得られた$2×2$複素行列$I_H,J_H,K_H$に$i$を掛ければ、パウリ行列$σ_1, σ_2, σ_3$が得られます。[[7shi-qp]]

$$
\begin{alignedat}{3}
i I_H
&= i \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix}
&&= \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
&&= σ_1
\\
i J_H
&= i \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
&&= \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}
&&= σ_2
\\
i K_H
&= i \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
&&= \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
&&= σ_3
\end{alignedat}
$$

このように、四元数の虚数単位$i, j, k$は、$(d,a,b,c)$の並び順から構成した複素行列表現を通じて、パウリ行列と対応付けられます。

# まとめ

四元数$q = a + bi + cj + dk$の演算は、4次元の実ベクトル空間における線形変換として捉えることができます。虚数単位$i, j, k$の左からの乗算は、$4×4$実行列$M_i, M_j, M_k$に対応します。ベクトルの成分を$(d,a,b,c)$という並び順（四元数の係数を$k,1,i,j$の順でベクトル化）に変更して、対応する$4×4$実行列$M'_i, M'_j, M'_k$を$2×2$のブロックに分解します。

&&& 四元数の実行列表現 $(d,a,b,c)$
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

&&& 複素数の行列表現
$$
R(α + iβ) =
  α\begin{pmatrix}1 &  0 \\ 0 & 1\end{pmatrix}
+ β \begin{pmatrix}0 & -1 \\ 1 & 0\end{pmatrix}
= \begin{pmatrix} α & -β \\ β & α \end{pmatrix}
\quad(α,β\in\mathbb{R})
$$
&&&

各ブロックを複素数に変換すれば、$2×2$の複素行列$I_H,J_H,K_H$が得られます。

&&& 虚数単位の複素行列表現
$$
I_H = \begin{pmatrix}  0 & -i \\ -i & 0 \end{pmatrix}, \quad
J_H = \begin{pmatrix}  0 & -1 \\  1 & 0 \end{pmatrix}, \quad
K_H = \begin{pmatrix} -i &  0 \\  0 & i \end{pmatrix}
$$
&&&

一般の四元数はこれらの線形結合で表され、得られた表現は四元数の積を保つ単射です。

&&& 四元数の複素行列表現
$$
ρ(a+bi+cj+dk)=\begin{pmatrix}a-id & -c-ib \\ c-ib & a+id \end{pmatrix},\quad
ρ(xy)=ρ(x)ρ(y)
$$
&&&

この$(d,a,b,c)$という並び順は、それによって得られた$I_H,J_H,K_H$に$i$を掛けることで、パウリ行列が構成できるように選ばれています。

&&& 四元数の行列表現とパウリ行列の関係
$$
σ_1 = iI_H = \begin{pmatrix} 0 &  1 \\ 1 &  0 \end{pmatrix},\quad
σ_2 = iJ_H = \begin{pmatrix} 0 & -i \\ i &  0 \end{pmatrix},\quad
σ_3 = iK_H = \begin{pmatrix} 1 &  0 \\ 0 & -1 \end{pmatrix}
$$
両辺に$i^{-1}=-i$を掛けると、以下の逆の関係が得られます。
$$
I_H=-iσ_1,\quad J_H=-iσ_2,\quad K_H=-iσ_3
$$
&&&
