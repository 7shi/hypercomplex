四元数の複素行列表現について、ベクトルの第1成分の虚部に反転を加える発想から素朴かつ発見的に構成します。

シリーズ：[四元数の行列表現](https://mathlog.info/series/PXPuUuQLYZk6HHho9eP8)

&&& 改訂履歴
- 2026.10.04 複素共役が複素線形でないことを写像の性質として示し、ベクトル表現と行列表現を定義して積の保存の証明を追加し、パウリ行列に合わせた表現を独立した節に移した
&&&

# 概要

四元数を2次元複素ベクトルとして素直に表現しようとすると、虚数単位$j$を掛ける操作に複素共役が現れ、複素行列による線形変換としては表せないことが分かります。本記事では、あらかじめベクトルの第1成分に虚部の符号反転（共役）という「ひねり」を組み込むことで、すべての虚数単位の左乗算を複素行列による線形変換として表します。得られた行列から四元数の複素行列表現を定め、一般の四元数の積を保つことを確かめます。さらに、得られた表現とパウリ行列との対応関係を整理します。[[7shi-qp]]

四元数の乗算規則（$i^2=j^2=k^2=ijk=-1$）と複素共役を前提とします。パウリ行列の物理的な応用は扱いません。

# 四元数の複素ベクトル表現

四元数$a+bi+cj+dk\ (a,b,c,d\in\mathbb R)$を複素数の組として表現します。

$$a+bi+cj+dk = (a+bi) + (c+di)j$$

これを単純に複素ベクトルに対応させます。

$$
(a+bi) + (c+di)j \mapsto \begin{pmatrix}a+bi \\ c+di\end{pmatrix}
$$

この対応によって、四元数の基本操作がどのように表現されるか確認してみます。

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

成分の対応には複素共役操作が含まれています。この変換を$z_1=a+bi,\ z_2=c+di$の関数として書くと、次のようになります。

$$
F\begin{pmatrix}z_1 \\ z_2\end{pmatrix}=\begin{pmatrix}-z_2^* \\ z_1^*\end{pmatrix}
$$

ベクトル$z$を$iz$に置き換えると、共役によって$i$の符号が反転します。

$$
F(iz)=-iF(z)
$$

したがって$F$は複素線形ではなく、複素行列を掛けるだけでは表現できません。

# ひねりを加えた表現

最初から$a-bi$が含まれるように、四元数とベクトルの対応に「ひねり」を加えます。第1成分の虚部の符号を反転させ、次のベクトルに対応させます。

&&&def ひねりを加えたベクトル表現
$$
v:\mathbb H\to\mathbb C^2,\quad a+bi+cj+dk \mapsto \begin{pmatrix}a-bi \\ c+di\end{pmatrix}
$$
&&&

&&&rem 第1成分を選ぶ理由
第2成分を共役にして$(a+bi,\ c-di)^\mathsf T$に対応させる方法もあります。その場合、得られる行列は本記事の行列の成分をそれぞれ共役にしたものになります。ここでは、$j$を掛けた結果の第2成分に現れる$a-bi$を入力の側に用意しておくため、第1成分を共役にします。
&&&

この表現では、左から$j$を掛ける操作は次のようになります。

$$
\begin{aligned}
j(a+bi+cj+dk)
&= aj-bk-c+di \\
&= (-c+di)+(a-bi)j \\
v(j(a+bi+cj+dk)) &= \begin{pmatrix}-c-di \\ a-bi\end{pmatrix}
\end{aligned}
$$

この操作をベクトルの変換として考えます。

$$
\begin{pmatrix}a-bi \\ c+di\end{pmatrix}
\xrightarrow{j \times}
\begin{pmatrix}-c-di \\ a-bi\end{pmatrix}
$$

この変換は成分の入れ替えと符号の反転だけで、複素共役を含みません。そのため行列で表せます。

$$
\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
\begin{pmatrix} a-bi \\ c+di \end{pmatrix}
= \begin{pmatrix} -c-di \\ a-bi \end{pmatrix}
$$

このように、四元数とベクトルの対応に最初から複素共役の形を埋め込むことで、線形変換として扱うことが可能になります。

この対応の変更は実線形ですが、複素線形ではありません。元の複素ベクトル空間で基底を取り替えるのではなく、四元数と複素ベクトルの対応の仕方そのものを変えています。

## 完全な行列表現

この対応によって、$i$および$k$による左乗算も複素行列で表現できます。

$$
\begin{aligned}
i(a+bi+cj+dk) &= (-b+ai)+(-d+ci)j \\
\begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
\begin{pmatrix} a-bi \\ c+di \end{pmatrix}
&= \begin{pmatrix} -b-ai \\ -d+ci \end{pmatrix}
\end{aligned}
$$

$i$の左乗算は、単純な対応では$\operatorname{diag}(i,i)$でしたが、ここでは第1成分が$-i$倍になります。第1成分は共役を取っているため、$(i(a+bi))^*=-i(a+bi)^*$となるからです。第2成分は共役にしていないので$i$倍のままです。

$$
\begin{aligned}
k(a+bi+cj+dk) &= (-d-ci)+(b+ai)j \\
\begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix}
\begin{pmatrix} a-bi \\ c+di \end{pmatrix}
&= \begin{pmatrix} -d+ci \\ b+ai \end{pmatrix}
\end{aligned}
$$

ここまでの結果をまとめます。

&&&fml 虚数単位の左乗算の行列
$$
\begin{alignedat}{3}
i\times&:\quad &&M_I &&= \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix} \\
j\times&:\quad &&M_J &&= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \\
k\times&:\quad &&M_K &&= \begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix}
\end{alignedat}
$$
&&&

これらの行列は四元数の乗算規則を満たしています。単位元$1$は単位行列$I$に対応させます。

$$
\begin{alignedat}{3}
&i^2=j^2=k^2=-1\ &&\Leftrightarrow\ &&{M_I}^2={M_J}^2={M_K}^2= -I \\
&ij=-ji=k\ &&\Leftrightarrow\ &&M_IM_J=-M_JM_I=M_K \\
&jk=-kj=i\ &&\Leftrightarrow\ &&M_JM_K=-M_KM_J=M_I \\
&ki=-ik=j\ &&\Leftrightarrow\ &&M_KM_I=-M_IM_K=M_J
\end{alignedat}
$$

基底に対応する行列を実線形結合すると、次の行列が得られます。

$$
\begin{aligned}
&aI+bM_I+cM_J+dM_K \\
=\ &a\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}
   +b\begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
   +c\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
   +d\begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix} \\
=\ &\begin{pmatrix} a-bi & -(c-di) \\ c+di & a+bi \end{pmatrix}
\end{aligned}
$$

この行列を四元数$a+bi+cj+dk$に対応させます。

&&&def 四元数の複素行列表現
$$
M:\mathbb H\to M_2(\mathbb C),\quad a+bi+cj+dk\mapsto aI+bM_I+cM_J+dM_K
$$
&&&

定義より$M_I=M(i),\ M_J=M(j),\ M_K=M(k)$です。$M_I,M_J,M_K$はそれぞれ$i,j,k$の左乗算を$v$を通じて表しているので、実線形性から、任意の四元数$x,y$について次が成り立ちます。

$$
M(x)v(y)=v(xy)
$$

特に$v(1)=(1,0)^\mathsf T$なので$M(x)v(1)=v(x)$となります。$M(x)$の第1列がベクトル表現$v(x)$に一致するのはこのためです。

# 積の保存

虚数単位どうしの乗算規則は確認しました。行列$M(x)$を四元数$x$の表現と呼ぶには、一般の四元数の積が行列の積に対応している必要があります。

&&&prop 積の保存
任意の四元数$x,y$に対して、次が成り立ちます。
$$
M(xy)=M(x)M(y),\quad M(1)=I
$$
&&&

&&&prf
$M(1)=I$は定義から明らかである。基底$1,i,j,k$のどの2つの積についても$M(xy)=M(x)M(y)$が成り立つことは、$M(1)=I$と上で確認した乗算規則から分かる。$M$は実線形であり、四元数の積と行列の積はどちらも実双線形なので、$x,y$を基底の実線形結合に展開すれば、任意の$x,y$についても$M(xy)=M(x)M(y)$が成り立つ。
&&&

# パウリ行列に合わせた表現

四元数の行列表現には任意性があります。パウリ行列の形に合わせるため、$j$の像を保ったまま、$i,k$の像を取り替えた別の行列表現を選びます。[[7shi-qp]]

&&&def パウリ行列に合わせた虚数単位の行列
$$
\begin{alignedat}{3}
I_H &:=-&&M_K &&= \begin{pmatrix} 0 & -i \\ -i & 0 \end{pmatrix} \\
J_H &:= &&M_J &&= \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \\
K_H &:= &&M_I &&= \begin{pmatrix} -i & 0 \\ 0 & i \end{pmatrix}
\end{alignedat}
$$
&&&

&&&rem $K_H$の決まり方
$I_H,J_H$を決めれば、残りは$K_H=I_HJ_H=-M_KM_J=M_I$より自動的に決まります。
&&&

それぞれに複素数の虚数単位$i$を掛けると、パウリ行列$\sigma_1,\sigma_2,\sigma_3$が得られます。ここでの$i$は行列の成分に使う複素数としての$i$で、行列$M_I$を掛けることではありません。

$$
iI_H=\sigma_1=\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix},\quad
iJ_H=\sigma_2=\begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix},\quad
iK_H=\sigma_3=\begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$

# まとめ

本記事では、四元数を複素ベクトルに対応させる際に第1成分を共役にすることで、四元数の左乗算を複素行列で表しました。この対応の変更自体は複素線形ではありませんが、共役を含んでいた左乗算は、変更後の座標で複素線形な行列として表せます。得られた行列表現は一般の四元数の積を保ち、$j$の像を保ったまま$i,k$の像を取り替えることで、パウリ行列に対応する表現が得られます。

&&& ベクトル表現
$$
v(a+bi+cj+dk)=\begin{pmatrix}a-bi \\ c+di\end{pmatrix}
$$
&&&

&&& 行列表現
$$
M(a+bi+cj+dk)=\begin{pmatrix} a-bi & -(c-di) \\ c+di & a+bi \end{pmatrix}
$$
&&&

&&& 積の保存
$$
M(x)v(y)=v(xy),\quad M(xy)=M(x)M(y)
$$
&&&
