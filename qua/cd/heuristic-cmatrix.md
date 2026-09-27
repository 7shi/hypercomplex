四元数の複素行列表現について、ベクトルの第1成分の虚部に反転を加える発想から素朴かつ発見的に構成します。

シリーズ：[四元数の行列表現](https://mathlog.info/series/PXPuUuQLYZk6HHho9eP8)

# 概要

四元数を2次元複素ベクトルとして素直に表現しようとすると、虚数単位$j$を掛ける操作に複素共役が現れ、通常の複素行列（線形変換）としては表現できない壁に突き当たります。本記事では、あらかじめベクトルの第1成分に虚部の符号反転（共役）という「ひねり」を組み込むことで、すべての虚数単位の左乗算を正則な複素行列として実現する構成法を解説します。さらに、得られた表現とパウリ行列との対応関係[[7shi-qp]]を整理します。

# 四元数の複素ベクトル表現

四元数$a+bi+cj+dk\ (a,b,c,d\in\mathbb R)$を複素数の組として表現します。

$$a+bi+cj+dk = (a+bi) + (c+di)j$$

これを単純に複素ベクトルに対応させます。

$$
(a+bi) + (c+di)j \mapsto \begin{pmatrix}a+bi \\ c+di\end{pmatrix}
$$

このベクトル表現によって、四元数の基本操作がどのように表現されるか確認してみます。

## 左から$i$を掛ける操作

四元数に左から$i$を掛ける操作を考えます。計算してベクトルに対応させます。

$$
\begin{aligned}
i(a+bi+cj+dk)
&= ai-b+ck-dj \\
&= (-b+ai)+(-d+ci)j \\
&\mapsto \begin{pmatrix}-b+ai \\ -d+ci\end{pmatrix}
\end{aligned}
$$
この操作をベクトルの変換として考えます。

$$
\begin{pmatrix}a+bi \\ c+di\end{pmatrix}
\xrightarrow{i \times}
\begin{pmatrix}-b+ai \\ -d+ci\end{pmatrix}
$$

これは各成分に$i$を掛けているだけのため、行列による変換は以下のようになります。

$$
\begin{pmatrix} i & 0 \\ 0 & i \end{pmatrix}
\begin{pmatrix} a+bi \\ c+di \end{pmatrix}
= \begin{pmatrix}-b+ai \\ -d+ci\end{pmatrix}
$$

## 左から$j$を掛ける操作

次に、四元数に左から$j$を掛ける操作を考えます。

$$
\begin{aligned}
j(a+bi+cj+dk)
&= aj-bk-c+di \\
&= (-c+di)+(a-bi)j \\
&\mapsto \begin{pmatrix}-c+di \\ a-bi\end{pmatrix}
\end{aligned}
$$

この操作をベクトルの変換として考えます。

$$
\begin{pmatrix}a+bi \\ c+di\end{pmatrix}
\xrightarrow{j \times}
\begin{pmatrix}-c+di \\ a-bi\end{pmatrix}
$$

この変換を行列で表現できるでしょうか。

$$
\begin{pmatrix} ? & ? \\ ? & ? \end{pmatrix}
\begin{pmatrix} a+bi \\ c+di \end{pmatrix}
\stackrel{?}{=} \begin{pmatrix} -c+di \\ a-bi \end{pmatrix}
$$

成分の対応には複素共役操作が含まれており、複素共役は線形変換ではないため、行列を掛けるだけでは表現できません。

# ひねりを加えた表現

最初から$a-bi$が含まれるように、四元数からベクトル表現へのマッピングに「ひねり」を加えます。（第1成分の虚部の符号反転）

$$
a+bi+cj+dk \mapsto \begin{pmatrix}a-bi \\ c+di\end{pmatrix}
$$

この表現では、左から$j$を掛ける操作は次のようになります。

$$
\begin{aligned}
j(a+bi+cj+dk)
&= aj-bk-c+di \\
&= (-c+di)+(a-bi)j \\
&\mapsto \begin{pmatrix}-c-di \\ a-bi\end{pmatrix}
\end{aligned}
$$

この操作をベクトルの変換として考えます。

$$
\begin{pmatrix}a-bi \\ c+di\end{pmatrix}
\xrightarrow{j \times}
\begin{pmatrix}-c-di \\ a-bi\end{pmatrix}
$$

この変換には複素共役が含まれないため、行列で変換が可能です。

$$
\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
\begin{pmatrix} a-bi \\ c+di \end{pmatrix}
= \begin{pmatrix} -c-di \\ a-bi \end{pmatrix}
$$

このように、ベクトル表現へのマッピングに最初から複素共役の形を埋め込むことで、線形変換として扱うことが可能になります。

## 完全な行列表現

このマッピングによって、$i$および$k$による左乗算も複素行列で表現できます。

$$
\begin{aligned}
i(a+bi+cj+dk) &= (-b+ai)+(-d+ci)j \\
\begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
\begin{pmatrix} a-bi \\ c+di \end{pmatrix}
&= \begin{pmatrix} -b-ai \\ -d+ci \end{pmatrix}
\end{aligned}
$$
$$
\begin{aligned}
k(a+bi+cj+dk) &= (-d-ci)+(b+ai)j \\
\begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix}
\begin{pmatrix} a-bi \\ c+di \end{pmatrix}
&= \begin{pmatrix} -d+ci \\ b+ai \end{pmatrix}
\end{aligned}
$$

ここまでの結果をまとめます。

$$
\begin{alignedat}{3}
i\times&:\quad &&M_I &&= \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix} \\
j\times&:\quad &&M_J &&= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \\
k\times&:\quad &&M_K &&= \begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix}
\end{alignedat}
$$

これらの行列は四元数の乗算規則を満たしています。（$1$は単位行列$I$に対応）

$$
\begin{alignedat}{3}
&i^2=j^2=k^2=-1\ &&\Leftrightarrow\ &&{M_I}^2={M_J}^2={M_K}^2= -I \\
&ij=-ji=k\ &&\Leftrightarrow\ &&M_IM_J=-M_JM_I=M_K \\
&jk=-kj=i\ &&\Leftrightarrow\ &&M_JM_K=-M_KM_J=M_I \\
&ki=-ik=j\ &&\Leftrightarrow\ &&M_KM_I=-M_IM_K=M_J
\end{alignedat}
$$

行列表現を線形結合させると、第1列がベクトル表現に一致することが分かります。

$$
\begin{aligned}
a+bi+cj+dk
\ \Leftrightarrow\ &aI+bM_I+cM_J+dM_K \\
=\ &a\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}
   +b\begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
   +c\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
   +d\begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix} \\
=\ &\begin{pmatrix} a-bi & -(c-di) \\ c+di & a+bi \end{pmatrix} \\
\mapsto\ &\begin{pmatrix} a-bi \\ c+di \end{pmatrix} \\
\end{aligned}
$$

# まとめ

本記事では、複素共役をあらかじめ埋め込む発見的手法によって四元数の複素行列表現を構成しました。

&&& 複素ベクトル表現のひねり
第1成分の虚部の符号を反転させることで、複素共役を伴う非線形な操作を通常の行列の線形変換として扱えるようになります。
$$
a+bi+cj+dk \mapsto \begin{pmatrix}a-bi \\ c+di\end{pmatrix}
$$
&&&

&&& 四元数の複素行列表現
$$
a+bi+cj+dk \iff \begin{pmatrix} a-bi & -(c-di) \\ c+di & a+bi \end{pmatrix}
$$
各虚数単位の左乗算は、四元数の乗法規則を満たす$2 \times 2$複素行列として表現されます。
&&&

&&&rem パウリ行列との対応関係
四元数の行列表現には任意性があり、基底を組み替えることでパウリ行列（それぞれに$i$を掛けたもの）に対応させることができます。[[7shi-qp]]
$$
\begin{alignedat}{3}
I_H &:=-&&M_K &&= \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix} \\
J_H &:= &&M_J &&= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \\
K_H &:= &&M_I &&= \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
\end{alignedat}
$$
$I_H, J_H$を固定すれば、$K_H = I_H J_H = -M_K M_J = M_I$より自動的に定まります。
&&&
