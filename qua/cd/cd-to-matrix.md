複素数の組から四元数を合成するケイリー＝ディクソンの構成法を出発点として、共役を伴う複素行列表現を発見的に導出します。

シリーズ：[四元数の行列表現](https://mathlog.info/series/PXPuUuQLYZk6HHho9eP8)

&&& 改訂履歴
- 2026.10.04 目標値の調整の根拠を補い、行列とベクトルへの写像を定義して積の保存の証明を追加し、既出の行列表現との一致を独立した節に移した
&&&

# 概要

四元数は2つの複素数の順序対として表すことができ、その乗法規則はケイリー＝ディクソンの構成法によって与えられます。[[wiki-cd]]

本記事では、この積公式を目標に据え、行列とベクトルの積として書き表すことを試みます。未知の行列成分を仮に置き、共役の不一致を1つずつ修正していくことで、第1成分に複素共役が現れる四元数の複素行列表現を段階的に導出します。積公式に現れる共役が、行列とベクトルの成分にどう配置されるかを観察します。最後に、得られた行列が四元数の積を保つことを確かめ、別の発想で構成した既出の行列表現と一致することを確認します。

四元数の乗算規則（$i^2=j^2=k^2=ijk=-1$）と複素共役を前提とします。扱うのは複素数から四元数を構成する段階だけで、四元数から八元数への拡張は扱いません。

# ケイリー＝ディクソンの構成法

ケイリー＝ディクソンの構成法は、代数を合成して高次元の代数を得るための規則です。[[wiki-cd]]

&&&rem 扱う構成法の範囲
本記事では、複素数の組から四元数を合成する構成法だけを取り扱います。
&&&

四元数$a+bi+cj+dk\ (a,b,c,d\in\mathbb R)$を2つの複素数$p=a+bi,\ q=c+di$の組として表現します。複素数全体$\mathbb C$は、虚数単位$i$の表記を維持したまま$\{a+bi\mid a,b\in\mathbb R\}\subset\mathbb H$と同一視します。

$$
a+bi+cj+dk=(a+bi)+(c+di)j=p+qj
$$

$j$と複素数の順序を入れ替えると、複素数に共役が付きます。共役を$*$で表記します。

&&&fml $j$と複素数の交換規則
$$
jp=p^*j\quad(p=a+bi\in\mathbb C;\ a,b\in\mathbb R)
$$
&&&

&&&prf
$$
jp=j(a+bi)=aj+bji=aj-bij=(a-bi)j=p^*j
$$
&&&

これを利用して、複素数$r,s$で構成される四元数$r+sj$との積を計算します。

$$
\begin{aligned}
(p+qj)(r+sj)
&=pr+psj+qjr+qjsj \\
&=pr+psj+qr^*j+qs^*j^2 \\
&=(pr-qs^*)+(ps+qr^*)j
\end{aligned}
$$

逆に、複素数の順序対全体に成分ごとの加法と次の積を定めることで、四元数を構成できます。これを**ケイリー＝ディクソンの構成法**と呼びます。

&&&def ケイリー＝ディクソンの構成法（複素数 → 四元数）
$$
(p,q)(r,s)=(pr-qs^*,ps+qr^*)\quad(p,q,r,s\in\mathbb C)
$$
&&&

&&&rem 八元数への拡張時の注意
ここでは成分が複素数なので、その積の順序を入れ替えられます。四元数から八元数を構成する場合は成分が非可換になるため、積の順序を指定する必要があります。この式の成分を単に四元数へ置き換えることはできません。
&&&

# 目標設定

ケイリー＝ディクソンの構成法を出発点にして、四元数の行列表現を導出します。左から掛ける四元数$p+qj$を固定し、右の四元数$r+sj$を列ベクトルとみなして、積を行列とベクトルの積として表すことを目指します。イメージを示します。

$$
\begin{pmatrix} ? & ? \\ ? & ? \end{pmatrix} \begin{pmatrix} r \\ s \end{pmatrix}
\stackrel{?}{=}\begin{pmatrix} pr-qs^* \\ ps+qr^* \end{pmatrix}
$$

## 仮計算

未知の成分に変数を割り当てて計算します。行列は$p,q$だけで決まり、任意の$r,s$に対して働くものとします。

$$
\begin{pmatrix}w & x \\ y & z\end{pmatrix}
\begin{pmatrix}r \\ s\end{pmatrix}
=\begin{pmatrix}wr+xs \\ yr+zs\end{pmatrix}
$$

この結果を目標の成分と比較します。

$$
\begin{pmatrix}wr+xs \\ yr+zs\end{pmatrix}
\stackrel{?}{=}\begin{pmatrix} pr-qs^* \\ qr^*+ps \end{pmatrix}
$$

両辺を見比べ、共役の不一致をひとまず無視して、次の行列を候補とします。

$$
\begin{pmatrix}w & x \\ y & z\end{pmatrix}
=\begin{pmatrix}p & -q \\ q & p\end{pmatrix}
$$

改めて計算して比較します。

$$
\begin{pmatrix}p & -q \\ q & p\end{pmatrix}
\begin{pmatrix}r \\ s\end{pmatrix}
=\begin{pmatrix}pr-qs \\ qr+ps\end{pmatrix}
\stackrel{?}{=} \begin{pmatrix} pr-qs^* \\ qr^*+ps \end{pmatrix}
$$

## 第2成分の一致

まず第2成分を合わせることに焦点を当てます。行列計算の第2成分 $qr + ps$ が目標の $qr^*+ps$ と一致するには、$r$ を $r^*$ に置き換える必要があります。

そこで四元数のベクトル表現を修正します。

$$
\begin{pmatrix} p & -q \\ q & p \end{pmatrix} \begin{pmatrix} r^* \\ s \end{pmatrix}
=\begin{pmatrix} pr^* - qs \\ qr^* + ps \end{pmatrix}
\stackrel{?}{=} \begin{pmatrix} pr-qs^* \\ qr^*+ps \end{pmatrix}
$$

これにより第2成分は目標と一致しましたが、第1成分にはまだ不一致があります。

## 目標値の調整

入力の表現を$\begin{pmatrix} r^* \\ s \end{pmatrix}$に変更したので、積の結果も同じ規則で表します。積を$u+wj\ (u,w\in\mathbb C)$とすれば、そのベクトル表現は$\begin{pmatrix} u^* \\ w \end{pmatrix}$です。

したがって、目標値の第1成分の共役を取ります。

$$
\begin{pmatrix} pr-qs^* \\ qr^*+ps \end{pmatrix}
\rightarrow \begin{pmatrix} (pr-qs^*)^* \\ qr^*+ps \end{pmatrix}
=\begin{pmatrix} p^*r^*-q^*s \\ qr^*+ps \end{pmatrix}
$$

## 行列表現の調整

改めてここまでの計算式と比較します。

$$
\begin{pmatrix} p & -q \\ q & p \end{pmatrix} \begin{pmatrix} r^* \\ s \end{pmatrix}
=\begin{pmatrix} pr^* - qs \\ qr^* + ps \end{pmatrix}
\stackrel{?}{=} \begin{pmatrix} p^*r^*-q^*s \\ qr^*+ps \end{pmatrix}
$$

目標値と一致させるため、行列の第1行の共役を取ります。

$$
\begin{pmatrix} p^* & -q^* \\ q & p \end{pmatrix} \begin{pmatrix} r^* \\ s \end{pmatrix}
=\begin{pmatrix} p^*r^* - q^*s \\ qr^* + ps \end{pmatrix}
=\begin{pmatrix} p^*r^*-q^*s \\ qr^*+ps \end{pmatrix}
$$

これで、入力と出力に同じベクトル表現を用いた等式が得られました。

ここで得られた対応を、四元数から行列とベクトルへの写像として定めます。

&&&def 四元数の複素行列表現とベクトル表現
$$
M:\mathbb H\to M_2(\mathbb C),\quad p+qj\mapsto\begin{pmatrix}p^*&-q^*\\q&p\end{pmatrix}
$$
$$
v:\mathbb H\to\mathbb C^2,\quad p+qj\mapsto\begin{pmatrix}p^*\\q\end{pmatrix}
$$
&&&

ベクトル$v(x)$は行列$M(x)$の第1列です。ここまでの計算は、四元数$x=p+qj,\ y=r+sj$に対して次のように書けます。

$$
M(x)v(y)=v(xy)
$$

# 積の保存

行列とベクトルの積で四元数の積を表せることが分かりました。行列$M(x)$を四元数$x$の表現と呼ぶには、行列どうしの積が四元数の積に対応している必要があります。

&&&prop 積の保存
任意の四元数$x,y$に対して、次が成り立つ。
$$
M(x)M(y)=M(xy),\quad M(1)=I
$$
&&&

&&&prf
任意の四元数$z$に対して、$M(x)v(y)=v(xy)$と四元数の結合法則より次が成り立つ。

$$
M(x)M(y)v(z)=M(x)v(yz)=v(x(yz))=v((xy)z)=M(xy)v(z)
$$

$z=r+sj$の$r,s$は任意の複素数を取れるので、$v(z)$は$\mathbb C^2$全体を動く。したがって$M(x)M(y)=M(xy)$である。また$1=1+0j$より$M(1)=I$である。
&&&

# 既出の行列表現との一致

$p=a+bi,\ q=c+di$を代入して、行列を実数の成分で書きます。

$$
M((a+bi)+(c+di)j)
=\begin{pmatrix} (a+bi)^* & -(c+di)^* \\ c+di & a+bi \end{pmatrix}
=\begin{pmatrix} a-bi & -(c-di) \\ c+di & a+bi \end{pmatrix}
$$

これは、ベクトルの第1成分を共役にする発想から構成した行列表現と一致します。[[7shi-qcm]]

ベクトル表現の選び方によって、行列内の共役や符号の配置は変わります。

# まとめ

本記事では、ケイリー＝ディクソンの構成法を行列とベクトルの演算として再構成し、四元数の複素行列表現を導出しました。得られた行列は四元数の積を保ち、既出の行列表現と一致します。

&&& 四元数・順序対・行列・ベクトルの対応
$p,q\in\mathbb C$とする。
$$
(p,q)\mapsto p+qj
$$
$$
M(p+qj)=\begin{pmatrix} p^* & -q^* \\ q & p \end{pmatrix},\quad
v(p+qj)=\begin{pmatrix} p^* \\ q \end{pmatrix}
$$
&&&

&&&rem 共役とベクトルの関係
- 積公式に現れる共役は、$ji=-ij$に由来します。今回のベクトル表現では、その共役を行列の第1行に配置することで、左乗算を複素線形な変換として表しています。
- ベクトルは行列の第1列のみを取り出した形です。
&&&

&&& 同じ積の3つの表し方
$p,q,r,s\in\mathbb C$とする。
$$
(p+qj)(r+sj)=(pr-qs^*)+(ps+qr^*)j
$$
$$
(p,q)(r,s)=(pr-qs^*,ps+qr^*)
$$
$$
\begin{pmatrix} p^* & -q^* \\ q & p \end{pmatrix} \begin{pmatrix} r^* \\ s \end{pmatrix}
=\begin{pmatrix} (pr - qs^*)^* \\ qr^* + ps \end{pmatrix}
$$
&&&

&&& 積の保存
$$
M(x)M(y)=M(xy)
$$
&&&
