四元数のテンソル積[[7shi-tp]]からクリフォード代数を構成する方法を確認します。同様に分解型四元数でも確認します。

# 四元数のテンソル積

四元数（全体の集合$\mathbb{H}$）は実数体上の4次元の除法代数であり、その基底は通常$\{1,i,j,k\}$で表されます。これらの基底は以下の関係式を満たします。

&&&def 四元数の基底の関係式
$$
i^2 = j^2 = k^2 = -1,\ k=ij=-ji
$$
&&&

四元数のテンソル積（全体の集合$\mathbb{H}⊗\mathbb{H}$）は16次元の代数であり、その基底は四元数の基底のテンソル積$\{1,i,j,k\}⊗\{1,i,j,k\}$から得られます。

&&&def $\mathbb{H}⊗\mathbb{H}$の基底
\begin{alignedat}{4}
\{
1&⊗1,\ &1&⊗i,\ &1&⊗j,\ &1&⊗k, \\
i&⊗1,&i&⊗i,&i&⊗j,&i&⊗k, \\
j&⊗1,&j&⊗i,&j&⊗j,&j&⊗k, \\
k&⊗1,&k&⊗i,&k&⊗j,&k&⊗k
\}
\end{alignedat}
&&&

# クリフォード代数の生成元

$\mathbb H$と$\mathbb H⊗\mathbb H$をクリフォード代数として扱う場合の、適切な生成元の選択とそれによって得られる代数構造を説明します。

&&&def クリフォード代数の生成元
グレード1の基底を生成元とする。$n$個の生成元$\{e_1,e_2,\cdots,e_n\}$から$2^n$個の基底が生成される。

生成元は互いに反交換関係を満たす。
$$
e_ie_j=-e_je_i\quad (i \ne j)
$$
&&&

$\mathbb H$の基底$\{1,i,j,k\}$は$4=2^2$個であることから、クリフォード代数としての生成元は2個です。虚数単位$i,j,k$は代数的な性質が同一なため、任意の2個を選択して生成元とすることができます。

$\{i,j\}$を生成元とすれば、$ij=k$より$k$はグレード2の基底に対応します。

&&&ex 四元数とクリフォード代数の基底の対応
$$
\begin{array}{c|cccc}
\mathbb H&1&i&j&k \\
\hline
\operatorname{Cl}&1&e_1&e_2&e_1e_2
\end{array}
$$
&&&

本記事ではこの組み合わせを使用します。

## $\mathbb H⊗\mathbb H$

$\mathbb{H}⊗\mathbb{H}$の基底は$4×4=16=2^4$個であることから、クリフォード代数としての生成元は4個です。

$\{1,i,j,k\}⊗\{1,i,j,k\}$から反交換関係を満たす基底の集合を総当たりで探索したところ、元の個数が3個と5個に分かれました。[[7shi-colab-cl]]

クリフォード代数として必要な生成元は4個であることから、3個では不足するため、5個の元からなる集合を選択します。

$$\{1⊗i,\ 1⊗j,\ i⊗k,\ j⊗k,\ k⊗k\} \tag{1}$$
$$\{1⊗i,\ 1⊗k,\ i⊗j,\ j⊗j,\ k⊗j\} \tag{2}$$
$$\{1⊗j,\ 1⊗k,\ i⊗i,\ j⊗i,\ k⊗i\} \tag{3}$$
$$\{i⊗1,\ j⊗1,\ k⊗i,\ k⊗j,\ k⊗k\} \tag{4}$$
$$\{i⊗1,\ j⊗i,\ j⊗j,\ j⊗k,\ k⊗1\} \tag{5}$$
$$\{i⊗i,\ i⊗j,\ i⊗k,\ j⊗1,\ k⊗1\} \tag{6}$$

これらは虚数単位$i,j,k$の置換や、左右の因子の交換によって移り合うため、本質的に同じ代数構造を与えます。

## 生成元の候補の選択

模式化しやすさの観点から、$(3)$を生成元の候補として選択します。

$$\{1⊗j,\ 1⊗k,\ i⊗i,\ j⊗i,\ k⊗i\} \tag{3}$$

&&&rem 模式化しやすさ
次のセクションの$(7)$のこと
&&&

$(3)$の性質を調べます。

&&&prop
$(3)$は$\mathbb H⊗\mathbb H$の基底を生成する。
&&&

&&&prf
$\mathbb H$の基底$\{1,i,j,k\}$は$\{i,j\}$から生成される。
$$
i^4=1,\ ij=k
$$

$\mathbb H⊗\mathbb H$は、$\{i,j\}⊗1$から$\{1,i,j,k\}⊗1$、$1⊗\{i,j\}$から$1⊗\{1,i,j,k\}$が生成され、これらの積からすべての基底が生成される。

$(3)$は$1⊗j$を含んでおり、残りは以下のように生成される。
$$
\begin{aligned}
(k⊗i)(j⊗i)&=i⊗1 \\
(i⊗i)(k⊗i)&=j⊗1 \\
(1⊗j)(1⊗k)&=1⊗i
\end{aligned}
$$

よって$(3)$から$\mathbb H⊗\mathbb H$の基底が生成される。
&&&

&&&prop
$(3)$において、任意の4元から残りの1元が生成される。
&&&

&&&prf
直接計算により示す。
$$
\begin{aligned}
(1⊗k)(1⊗j)(i⊗i)(j⊗i)&=k⊗i \\
(1⊗k)(i⊗i)(j⊗i)(k⊗i)&=1⊗j \\
(i⊗i)(j⊗i)(k⊗i)(1⊗j)&=1⊗k \\
(j⊗i)(k⊗i)(1⊗k)(1⊗j)&=i⊗i \\
(k⊗i)(1⊗k)(1⊗j)(i⊗i)&=j⊗i
\end{aligned}
$$
&&&

## 生成元の選択

互いに反交換ですべての基底を生成する最小の元の集合は、クリフォード代数の生成元となります。命題1,2より、$(3)$から任意の4元を選択すれば、クリフォード代数の生成元となります。

代数的な性質は生成元の選択に依存しませんが、どのようにすればテンソル積の構造が解釈しやすいかを検討します。

$(3)$の構造を捉えるため模式化します。

$$
\begin{array}{c|ccccc}
\text{左因子}&\color{red}{i}&\color{red}{j}&\color{red}{k}&1&1 \\
&⊗&⊗&⊗&⊗&⊗ \\
\text{右因子}&i&i&\color{red}{i}&\color{red}{j}&\color{red}{k}
\end{array} \tag{7}
$$

&&&rem
赤字部分は、左因子と右因子に$\{i,j,k\}$が現れるように並べ替えた様子を示します。
&&&

左因子が元となった$\mathbb H$で、それに右因子を付加して拡張していると解釈します。$\mathbb H$のクリフォード代数としての生成元を$\{i,j\}$とすれば、左因子に$k$を含まないように生成元を選択することで、テンソル積による拡張の様子が分かりやすくなります。また、右因子の$k$は$ij$に書き換えます。

$$
\begin{array}{c|cccc}
\text{左因子}&\color{red}{i}&\color{red}{j}&1&1 \\
&⊗&⊗&⊗&⊗ \\
\text{右因子}&i&i&j&ij \\
\hline
\operatorname{Cl}&e_1&e_2&e_3&e_4
\end{array}
$$

&&&rem
赤字部分は、拡張前の$\mathbb H$におけるクリフォード代数の生成元です。それ以外が、テンソル積によって拡張された部分です。
&&&

本記事では、この組み合わせを$\mathbb H⊗\mathbb H$のクリフォード代数の生成元として使用します。

# 計量と符号数

生成元をそれぞれ2乗します。

$$
\begin{alignedat}{7}
(i&⊗i)^2 &&= &i^2&⊗i^2    &&= &(-1)&⊗(-1) &&= &\color{red}{ 1}&(1⊗1) \\
(j&⊗i)^2 &&= &j^2&⊗i^2    &&= &(-1)&⊗(-1) &&= &\color{red}{ 1}&(1⊗1) \\
(1&⊗j)^2 &&= &1^2&⊗j^2    &&= &   1&⊗(-1) &&= &\color{red}{-1}&(1⊗1) \\
(1&⊗k)^2 &&= &1^2&⊗(ij)^2 &&= &   1&⊗(-1) &&= &\color{red}{-1}&(1⊗1)
\end{alignedat}
$$

これら2乗の係数（赤字部分）を**計量**、計量の値ごとの生成元の個数を**符号数**と呼びます。本記事では$1,-1$の順に符号数を数えます。

- 計量$1$が2個、計量$-1$が2個 → 符号数$(2,2)$

&&&rem
クリフォード代数としての性質は符号数にのみ依存するため、通常、生成元の具体的な選択ではなく、符号数のみが添え字で示されます。

- 符号数$(p,q)$ → $\operatorname{Cl}_{p,q}(\mathbb R)$

資料によっては$-1,1$の順に符号数を数えるものがあります。どちらを使用しているかは確認が必要です。
&&&

計量を含めて、テンソル積による拡張の様子を示します。

$$
\begin{array}{ccc}

\begin{array}{c|cc}
\mathbb H&\color{red}{i}&\color{red}{j} \\
\hline
\operatorname{Cl}&e_1&e_2 \\
\hline
\text{計量}&-1&-1
\end{array} &

\xrightarrow{⊗\mathbb H} &

\begin{array}{c|cccc}
\mathbb H&\color{red}{i}&\color{red}{j}&1&1 \\
&⊗&⊗&⊗&⊗ \\
\mathbb H&i&i&j&ij \\
\hline
\operatorname{Cl}&e_1&e_2&e_3&e_4 \\
\hline
\text{計量}&1&1&-1&-1
\end{array} \\

\operatorname{Cl}_{0,2}(\mathbb R) & &
\operatorname{Cl}_{2,2}(\mathbb R)
\end{array}
$$

&&&rem
拡張の際、赤字部分の$i,j$の計量が$⊗i$によって反転しています。また、拡張された$e_3,e_4$の計量は$-1$です。
&&&

## 一般化と公式

$\operatorname{Cl}_{2,2}(\mathbb R)$を更に拡張しても、同じ構造が繰り返されます。

$$
\begin{array}{ccc}

\begin{array}{c|cccc}
\operatorname{Cl}&e_1&e_2&e_3&e_4 \\
\hline
\text{計量}&1&1&-1&-1
\end{array} &

\xrightarrow{⊗\mathbb H} &

\begin{array}{c|cccccc}
\operatorname{Cl}&e_1&e_2&e_3&e_4&1&1 \\
&⊗&⊗&⊗&⊗&⊗&⊗ \\
\mathbb H&i&i&i&i&j&ij \\
\hline
\text{計量}&-1&-1&1&1&-1&-1
\end{array} \\

\operatorname{Cl}_{2,2}(\mathbb R) & &
\operatorname{Cl}_{2,4}(\mathbb R)
\end{array}
$$

&&&rem
$e_1,e_2,e_3,e_4$の計量が$⊗i$によって反転して、追加された2個の生成元$\{1⊗j,\ 1⊗ij\}$の計量は$-1$です。
&&&

この構造を一般化します。

&&& クリフォード代数の$⊗\mathbb H$による拡張の一般化
$$
\begin{array}{ccc}

\begin{array}{c|cccccc}
\operatorname{Cl}&e_1&\cdots&e_p&e_{p+1}&\cdots&e_{p+q} \\
\hline
\text{計量}&1&\cdots&1&-1&\cdots&-1
\end{array} &

\xrightarrow{⊗\mathbb H} &

\begin{array}{c|cccccccc}
\operatorname{Cl}&e_1&\cdots&e_p&e_{p+1}&\cdots&e_{p+q}&1&1 \\
&⊗&\cdots&⊗&⊗&\cdots&⊗&⊗&⊗ \\
\mathbb H&i&\cdots&i&i&\cdots&i&j&ij \\
\hline
\text{計量}&-1&\cdots&-1&1&\cdots&1&-1&-1
\end{array} \\

\operatorname{Cl}_{p,q}(\mathbb R) & &
\operatorname{Cl}_{q,p+2}(\mathbb R)
\end{array}
$$
&&&

この結果を公式の形にまとめます。[[wiki-clif]]

&&&fml クリフォード代数の$⊗\mathbb H$による拡張
$$
\operatorname{Cl}_{p,q}(\mathbb R) ⊗ \mathbb H
\cong \operatorname{Cl}_{q,p+2}(\mathbb R)
$$
&&&

&&&rem
右辺の$\operatorname{Cl}_{q,p+2}(\mathbb R)$は$p,q$の位置が入れ替わります。これは$⊗i$による計量の反転に由来します。

符号数の増分$(0,2)$は$\mathbb H \cong \operatorname{Cl}_{0,2}(\mathbb R)$に由来します。
&&&

# 分解型四元数

四元数を一部変更して、$j$を$j^2=1$となる実数ではない虚数単位としたものが分解型四元数です。[[wiki-sq]]

&&&def 分解型四元数の基底の関係式
$$
i^2 = -1,\ j^2 = 1,\ k=ij=-ji \quad (j\ne\pm1)
$$
&&&

この定義から$k^2=1$が導かれます。

&&&prf
$$
k^2=(ij)(ij)=i(ji)j=i(-ij)j=-(ii)(jj)=-(-1)(1)=1
$$
&&&

分解型四元数全体の集合を$\mathbb H'$と表記します。

&&&rem 分解型符号数
分解型四元数は、共役との積によって定義された擬ノルムの2乗により、計量が決まります。（クリフォード代数とは異なり、基底の2乗がそのまま計量とはなりません）

$$
\begin{aligned}
\lVert a+bi+cj+dk\rVert^2
&=(a+bi+cj+dk)(a+bi+cj+dk)^* \\
&=(a+bi+cj+dk)(a-bi-cj-dk) \\
&=a^2+b^2-c^2-d^2
\end{aligned}
$$

擬ノルムの2乗に現れる符号数は$(2,2)$となり正と負の個数が等しくなります。このような符号数を**分解型**と呼び、代数名の由来となっています。
&&&

## クリフォード代数

虚数単位によって2乗の値が変わるため、どれをクリフォード代数の生成元として使うかで符号数が変わります。

1. $\operatorname{Cl}_{1,1}(\mathbb R)$: $\{i,j\},\{i,k\}$
2. $\operatorname{Cl}_{2,0}(\mathbb R)$: $\{j,k\}$

クリフォード代数としての性質は符号数にのみ依存するため、$\operatorname{Cl}_{1,1}(\mathbb R)$の生成元としては$\{i,j\}$のみを使用します。

&&&rem
グレード2の基底まで含めれば代数として同型です。
$$
\mathbb H' \cong \operatorname{Cl}_{1,1}(\mathbb R) \cong \operatorname{Cl}_{2,0}(\mathbb R)
$$
&&&

## $\mathbb H'⊗\mathbb H'$

2×2=4種類の組み合わせを確認します。

$$
\begin{array}{ccc}

\begin{array}{c|cc}
\mathbb H'&\color{red}{i}&\color{red}{j} \\
\hline
\text{計量}&-1&1
\end{array} &

\xrightarrow{⊗\{i,j\}} &

\begin{array}{c|cccc}
\mathbb H'&\color{red}{i}&\color{red}{j}&1&1 \\
&⊗&⊗&⊗&⊗ \\
\mathbb H'&i&i&j&ij \\
\hline
\text{計量}&1&-1&1&1
\end{array} \\

\operatorname{Cl}_{1,1}(\mathbb R) & &
\operatorname{Cl}_{3,1}(\mathbb R)

\\ \ \\

\begin{array}{c|cc}
\mathbb H'&\color{red}{j}&\color{red}{k} \\
\hline
\text{計量}&1&1
\end{array} &

\xrightarrow{⊗\{i,j\}} &

\begin{array}{c|cccc}
\mathbb H'&\color{red}{j}&\color{red}{k}&1&1 \\
&⊗&⊗&⊗&⊗ \\
\mathbb H'&i&i&j&ij \\
\hline
\text{計量}&-1&-1&1&1
\end{array} \\

\operatorname{Cl}_{2,0}(\mathbb R) & &
\operatorname{Cl}_{2,2}(\mathbb R)

\\ \ \\

\begin{array}{c|cc}
\mathbb H'&\color{red}{i}&\color{red}{j} \\
\hline
\text{計量}&-1&1
\end{array} &

\xrightarrow{⊗\{j,k\}} &

\begin{array}{c|cccc}
\mathbb H'&\color{red}{i}&\color{red}{j}&1&1 \\
&⊗&⊗&⊗&⊗ \\
\mathbb H'&j&j&k&jk \\
\hline
\text{計量}&-1&1&1&-1
\end{array} \\

\operatorname{Cl}_{1,1}(\mathbb R) & &
\operatorname{Cl}_{2,2}(\mathbb R)

\\ \ \\

\begin{array}{c|cc}
\mathbb H'&\color{red}{j}&\color{red}{k} \\
\hline
\text{計量}&1&1
\end{array} &

\xrightarrow{⊗\{j,k\}} &

\begin{array}{c|cccc}
\mathbb H'&\color{red}{j}&\color{red}{k}&1&1 \\
&⊗&⊗&⊗&⊗ \\
\mathbb H'&j&j&k&jk \\
\hline
\text{計量}&1&1&1&-1
\end{array} \\

\operatorname{Cl}_{2,0}(\mathbb R) & &
\operatorname{Cl}_{3,1}(\mathbb R)
\end{array}
$$

結果をまとめます。

$$
\mathbb H'⊗\mathbb H'
\cong \operatorname{Cl}_{3,1}(\mathbb R)
\cong \operatorname{Cl}_{2,2}(\mathbb R)
$$

この結果から、以下の関係が分かります。

&&&thm 同型対応
$$
\mathbb H⊗\mathbb H \cong \mathbb H'⊗\mathbb H'
$$
&&&

&&&prf
$$
\mathbb H⊗\mathbb H
\cong \operatorname{Cl}_{2,2}(\mathbb R)
\cong \mathbb H'⊗\mathbb H'
$$
&&&

## 一般化と公式

この構造を一般化します。

&&& クリフォード代数の$⊗\mathbb H'$による拡張の一般化
$$
\begin{array}{ccc}

\begin{array}{c|cccccc}
\operatorname{Cl}&e_1&\cdots&e_p&e_{p+1}&\cdots&e_{p+q} \\
\hline
\text{計量}&1&\cdots&1&-1&\cdots&-1
\end{array} &

\xrightarrow{⊗\{i,j\}} &

\begin{array}{c|cccccccc}
\operatorname{Cl}&e_1&\cdots&e_p&e_{p+1}&\cdots&e_{p+q}&1&1 \\
&⊗&\cdots&⊗&⊗&\cdots&⊗&⊗&⊗ \\
\mathbb H'&i&\cdots&i&i&\cdots&i&j&ij \\
\hline
\text{計量}&-1&\cdots&-1&1&\cdots&1&1&1
\end{array} \\

\operatorname{Cl}_{p,q}(\mathbb R) & &
\operatorname{Cl}_{q+2,p}(\mathbb R)

\\ \ \\

&\xrightarrow{⊗\{j,k\}} &

\begin{array}{c|cccccccc}
\operatorname{Cl}&e_1&\cdots&e_p&e_{p+1}&\cdots&e_{p+q}&1&1 \\
&⊗&\cdots&⊗&⊗&\cdots&⊗&⊗&⊗ \\
\mathbb H'&j&\cdots&j&j&\cdots&j&k&jk \\
\hline
\text{計量}&1&\cdots&1&-1&\cdots&-1&1&-1
\end{array} \\

& &
\operatorname{Cl}_{p+1,q+1}(\mathbb R)

\end{array}
$$
&&&

この結果を公式の形にまとめます。[[wiki-clif]]

&&&fml クリフォード代数の$⊗\mathbb H'$による拡張
$$
\begin{aligned}
\operatorname{Cl}_{p,q}(\mathbb R) ⊗ \mathbb H'
&\cong \operatorname{Cl}_{q+2,p}(\mathbb R) \\
&\cong \operatorname{Cl}_{p+1,q+1}(\mathbb R)
\end{aligned}
$$
&&&

&&&rem
中辺の$\operatorname{Cl}_{q+2,p}(\mathbb R)$は$p,q$の位置が入れ替わります。これは$⊗i$による計量の反転に由来します。
右辺の$\operatorname{Cl}_{p+1,q+1}(\mathbb R)$は$p,q$の位置が維持されます。これは$⊗j$によって計量が変化しないことに由来します。

符号数の増分$(2,0),(1,1)$は$\mathbb H' \cong \operatorname{Cl}_{2,0}(\mathbb R) \cong \operatorname{Cl}_{1,1}(\mathbb R)$に由来します。
&&&

## $\mathbb H'⊗\mathbb H,\ \mathbb H⊗\mathbb H'$

公式1より

$$
\begin{aligned}
\mathbb H'⊗\mathbb H

&\cong \operatorname{Cl}_{2,0}(\mathbb R)⊗\mathbb H
\cong \operatorname{Cl}_{0,4}(\mathbb R) \\

&\cong \operatorname{Cl}_{1,1}(\mathbb R)⊗\mathbb H
\cong \operatorname{Cl}_{1,3}(\mathbb R) \\
\end{aligned}
$$

公式2より

$$
\begin{aligned}
\mathbb H⊗\mathbb H'
\cong \operatorname{Cl}_{0,2}(\mathbb R)⊗\mathbb H'
&\cong \operatorname{Cl}_{4,0}(\mathbb R) \\
&\cong \operatorname{Cl}_{1,3}(\mathbb R)
\end{aligned}
$$

どちらも$\operatorname{Cl}_{1,3}(\mathbb R)$を含みますが、これは同型対応におけるテンソル積の可換性を反映しています。

$$
\mathbb{H} ⊗ \mathbb{H}'
\cong \mathbb{H}' ⊗ \mathbb{H}
\cong \operatorname{Cl}_{1,3}(\mathbb R)
\cong \operatorname{Cl}_{0,4}(\mathbb R)
\cong \operatorname{Cl}_{4,0}(\mathbb R)
$$

# まとめ

本記事では、$\mathbb H$と$\mathbb H'$とそれらのテンソル積によって構成されるクリフォード代数の構造を分析しました。

テンソル積によって$\mathbb H,\mathbb H'$を付加することで、クリフォード代数としての生成元は2個増えます。

&&&
$$
\begin{alignedat}{2}
&\mathbb H &&\cong \operatorname{Cl}_{0,2}(\mathbb R) \\
&\mathbb H' &&\cong \operatorname{Cl}_{2,0}(\mathbb R) \cong \operatorname{Cl}_{1,1}(\mathbb R)
\end{alignedat}
$$
$$
\begin{alignedat}{2}
\operatorname{Cl}_{p,q}(\mathbb R) &⊗ \mathbb H  &&\cong \operatorname{Cl}_{q,p+2}(\mathbb R) \\
\operatorname{Cl}_{p,q}(\mathbb R) &⊗ \mathbb H' &&\cong \operatorname{Cl}_{q+2,p}(\mathbb R) \cong \operatorname{Cl}_{p+1,q+1}(\mathbb R)
\end{alignedat}
$$
$$
\begin{alignedat}{6}
\mathbb{H} &⊗ \mathbb{H} &&\cong \mathbb{H}' ⊗ \mathbb{H}' &&\cong \operatorname{Cl}_{3,1}(\mathbb R) \cong \operatorname{Cl}_{2,2}(\mathbb R) \\
\mathbb{H} &⊗ \mathbb{H}' &&\cong \mathbb{H}' ⊗ \mathbb{H} &&\cong \operatorname{Cl}_{1,3}(\mathbb R) \cong \operatorname{Cl}_{0,4}(\mathbb R) \cong \operatorname{Cl}_{4,0}(\mathbb R)
\end{alignedat}
$$

$$
\begin{array}{|c|ccccc|}
\hline
p \backslash q & 0 & 1 & 2 & 3 & 4 \\
\hline
0 & & & \mathbb H & & \mathbb H⊗\mathbb H' \\
1 & & \mathbb H' & & \mathbb H⊗\mathbb H' \\
2 & \mathbb H' & & \mathbb H⊗\mathbb H \\
3 & & \mathbb H⊗\mathbb H \\
4 & \mathbb H⊗\mathbb H' \\
\hline
\end{array}
$$
&&&
