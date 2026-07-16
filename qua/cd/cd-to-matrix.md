シリーズ：[四元数の行列表現](https://mathlog.info/series/PXPuUuQLYZk6HHho9eP8)

複素数の組から四元数を合成するためのケイリー＝ディクソンの構成法を確認して、そこから発見的に行列表現を導出します。

# ケイリー＝ディクソンの構成法

ケイリー＝ディクソンの構成法は、代数を合成して高次元の代数を得るための規則です。[[wiki-cd]]

&&&rem
本稿では、複素数の組から四元数を合成する構成法だけを取り扱います。
&&&

四元数$a+bi+cj+dk$を2つの複素数$p=a+bi,\ q=c+di\ (a,b,c,d\in\mathbb R)$の組として表現します。複素数の虚数単位$i$は、表記を維持したまま四元数の虚数単位$i$として埋め込みます。($\mathbb C\hookrightarrow\mathbb H$)

$$
a+bi+cj+dk=(a+bi)+(c+di)j=p+qj
$$

$j$が（四元数の部分代数としての）複素数と交換するとき、複素数を共役に変化させます。共役を$*$で表記します。

&&&fml
$$
jp=p^*j\quad(p=a+bi\in\mathbb H;\ a,b\in\mathbb R)
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

複素数の組を順序対として表記したものが、ケイリー＝ディクソンの構成法と呼ばれます。

&&&def ケイリー＝ディクソンの構成法（複素数 → 四元数）
$$
(p,q)(r,s)=(pr-qs^*,ps+qr^*)\quad(p,q,r,s\in\mathbb C)
$$
&&&

&&&rem
これは四元数の構成に特化しています。四元数から八元数を構成する場合、このままでは適用できません。
&&&

# 目標設定

ケイリー＝ディクソンの構成法を出発点にして、四元数の行列表現を導出します。イメージを示します。

$$
\begin{pmatrix} ? & ? \\ ? & ? \end{pmatrix} \begin{pmatrix} r \\ s \end{pmatrix}
\stackrel{?}{=}\begin{pmatrix} pr-qs^* \\ ps+qr^* \end{pmatrix}
$$

## 仮計算

未知の成分に変数を割り当てて計算します。

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

両辺を見比べることで、共役の不一致はひとまず無視して対応を推定します。

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

そこで四元数のベクトル表現を修正します：

$$
\begin{pmatrix} p & -q \\ q & p \end{pmatrix} \begin{pmatrix} r^* \\ s \end{pmatrix}
=\begin{pmatrix} pr^* - qs \\ qr^* + ps \end{pmatrix}
\stackrel{?}{=} \begin{pmatrix} pr-qs^* \\ qr^*+ps \end{pmatrix}
$$

これにより第2成分は目標と一致しましたが、第1成分にはまだ不一致があります。

## 目標値の調整

$\begin{pmatrix} r^* \\ s \end{pmatrix}$の第1成分は共役になっていますが、これはベクトル表現に共通すると考えられます。

そこで目標値の第1成分の共役を取ります。

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
\stackrel{?}{=} \begin{pmatrix} p^*r^*-q^*s \\ qr^*+ps \end{pmatrix}
$$

無事に目標値と一致しました！

# まとめ

&&&def 四元数・順序対・行列・ベクトルの対応
$$
p+qj
\iff (p,q)
\iff \begin{pmatrix} p^* & -q^* \\ q & p \end{pmatrix}
\iff \begin{pmatrix} p^* \\ q \end{pmatrix}
\quad (p,q\in\mathbb C (\hookrightarrow \mathbb H))
$$
&&&

&&&rem
- 行列・ベクトル表現の第1行の共役が、四元数の非可換性を反映しています。
- ベクトルは行列の第1列のみを取り出した形です。
&&&

&&&fml 四元数の計算・ケイリー＝ディクソンの構成法・行列とベクトルの計算
$$
p,q,r,s\in\mathbb C (\hookrightarrow \mathbb H)
$$
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

$p=a+bi,q=c+di$の関係を展開して、四元数と行列を比較します。

$$
(a+bi)+(c+di)j
\iff \begin{pmatrix} (a+bi)^* & -(c+di)^* \\ c+di & a+bi \end{pmatrix}
$$

これは発見的に構成した行列表現と一致します。[[7shi-qcm]]
