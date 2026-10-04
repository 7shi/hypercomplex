四元数および分解型四元数のテンソル積から実クリフォード代数を構成する手順を整理し、テンソル積による生成元の拡張規則を体系化します。

&&& 改訂履歴
- 2026.10.04 同型を関係式・生成・次元の一致から導く命題を追加、「3個と5個」を極大な反交換集合と明記、拡張公式の一般的な証明を追加
&&&

# 概要

実クリフォード代数$\operatorname{Cl}_{p,q}(\mathbb R)$は、四元数$\mathbb H$や分解型四元数$\mathbb H'$とのテンソル積によって、生成元を2個ずつ増やして構成できることが知られています。本記事では、テンソル積の基礎概念を踏まえ、$\mathbb H$および$\mathbb H'$の基底から反交換関係を満たす生成元を体系的に選択し、$\mathbb H \otimes \mathbb H$や$\mathbb H \otimes \mathbb H'$がそれぞれどのクリフォード代数と同型になるかを具体的に導出します。さらに、一般の$\operatorname{Cl}_{p,q}(\mathbb R)$に対して$\otimes \mathbb H$や$\otimes \mathbb H'$を作用させた際の符号数（計量）の変化公式を導き、テンソル積による拡張構造を整理します。[[7shi-tp]]

# 四元数のテンソル積

本記事のテンソル積はすべて実数上の通常のテンソル積$\otimes_{\mathbb R}$とし、添え字を省略します。また、本記事の同型はすべて実結合代数としての同型で、グレードや生成元の空間の保存は要求しません。テンソル積の基本的な計算規則は先行記事を参照してください。[[7shi-tp]]

四元数（全体の集合$\mathbb{H}$）は実数体上の4次元の除法代数であり、その基底は通常$\{1,i,j,k\}$で表されます。これらの基底は以下の関係式を満たします。

&&&def 四元数の基底の関係式
$$
i^2 = j^2 = k^2 = -1,\ k=ij=-ji
$$
&&&

四元数のテンソル積（全体の集合$\mathbb{H}⊗\mathbb{H}$）は16次元の代数であり、その基底は四元数の基底のテンソル積$\{1,i,j,k\}⊗\{1,i,j,k\}$から得られます。

&&&ex $\mathbb{H}⊗\mathbb{H}$の基底
$$
\begin{alignedat}{4}
\{
1&⊗1,\ &1&⊗i,\ &1&⊗j,\ &1&⊗k, \\
i&⊗1,&i&⊗i,&i&⊗j,&i&⊗k, \\
j&⊗1,&j&⊗i,&j&⊗j,&j&⊗k, \\
k&⊗1,&k&⊗i,&k&⊗j,&k&⊗k
\}
\end{alignedat}
$$
&&&

# クリフォード代数の生成元

$\mathbb H$と$\mathbb H⊗\mathbb H$をクリフォード代数として扱う場合の、適切な生成元の選択とそれによって得られる代数構造を説明します。

&&&def クリフォード代数の生成元
$n=p+q$個の元$\{e_1,e_2,\cdots,e_n\}$が次の関係だけを満たすように生成される実代数を$\operatorname{Cl}_{p,q}(\mathbb R)$とする。
$$
e_i^2=\begin{cases}1&(i\le p)\\-1&(i>p)\end{cases},\qquad
e_ie_j=-e_je_i\quad (i \ne j)
$$
このとき、$1$と$e_{i_1}\cdots e_{i_r}\ (i_1<\cdots<i_r)$の$2^n$個が基底をなす。グレード1の基底$e_i$を生成元と呼ぶ。
&&&

ある代数の中で同じ関係を満たす元を見つけても、これらの単項式が一次独立になるとは限りません。たとえば、後で見る5個の反交換する元は、積が定数になるため16次元にしかなりません。そこで、代数との同型を次のように判定します。

&&&prop 生成元による同型の判定 [prop-cl-iso]
実代数$A$の$n$個の元$e_1,\cdots,e_n$が上の関係式を満たし、$A$を生成し、$\dim A=2^n$ならば、$A\cong\operatorname{Cl}_{p,q}(\mathbb R)$である。
&&&

&&&prf
$\operatorname{Cl}_{p,q}(\mathbb R)$は関係式だけで定まる代数なので、生成元を$e_i$に写す代数準同型$\operatorname{Cl}_{p,q}(\mathbb R)\to A$が存在する。像は$e_i$が生成する代数であり、$A$全体である。よってこの準同型は全射で、両辺の次元がともに$2^n$であるから同型である。
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

$\{1,i,j,k\}⊗\{1,i,j,k\}$から$1⊗1$を除く15個の基底について、互いに反交換し、包含関係で極大となる集合を総当たりで探索したところ、元の個数が3個と5個に分かれました。[[7shi-colab-cl]]

互いに反交換し2乗が$\pm1$となる3個の元が生成する代数は、単項式$2^3=8$個で張られるため、16次元の$\mathbb H⊗\mathbb H$全体を生成できません。3個では不足するため、5個の元からなる集合を選択します。

$$\{1⊗i,\ 1⊗j,\ i⊗k,\ j⊗k,\ k⊗k\} \tag{1}$$
$$\{1⊗i,\ 1⊗k,\ i⊗j,\ j⊗j,\ k⊗j\} \tag{2}$$
$$\{1⊗j,\ 1⊗k,\ i⊗i,\ j⊗i,\ k⊗i\} \tag{3}$$
$$\{i⊗1,\ j⊗1,\ k⊗i,\ k⊗j,\ k⊗k\} \tag{4}$$
$$\{i⊗1,\ j⊗i,\ j⊗j,\ j⊗k,\ k⊗1\} \tag{5}$$
$$\{i⊗i,\ i⊗j,\ i⊗k,\ j⊗1,\ k⊗1\} \tag{6}$$

これらは虚数単位$i,j,k$の置換や、左右の因子の交換によって移り合うため、本質的に同じ代数構造を与えます。

&&&rem 5個の組の形
$a⊗b$と$c⊗d$の積では、左因子どうしと右因子どうしの交換で符号が決まります。$1$は何とでも可換で、異なる虚数単位は反交換するので、2つの元が反交換するのは、左右の因子のうちちょうど一方が反交換する場合です。

探索の結果、5個の組は左右の交換を除いて、$\{a,b,c\}=\{i,j,k\}$として
$$\{1⊗a,\ 1⊗b,\ i⊗c,\ j⊗c,\ k⊗c\}$$
の形に限られました。$c$の選び方の3通りと左右の交換の2通りで、$(1)$〜$(6)$の6個になります。
&&&

## 生成元の候補の選択

次のセクションの模式図$(7)$で左右の因子の構造を整理しやすいことから、$(3)$を生成元の候補として選択します。

$$\{1⊗j,\ 1⊗k,\ i⊗i,\ j⊗i,\ k⊗i\} \tag{3}$$

$(3)$の性質を調べます。

&&&prop 基底の生成 [prop-tensor-span]
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

&&&prop 4元による残りの元の生成 [prop-tensor-four]
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

[[prop-tensor-span]]と[[prop-tensor-four]]より、$(3)$から任意の4元を選択すれば、それらは$\mathbb H⊗\mathbb H$を生成します。これらは互いに反交換し、2乗は$\pm1$であり（次節で計算します）、$\dim(\mathbb H⊗\mathbb H)=16=2^4$なので、[[prop-cl-iso]]により、対応する符号数のクリフォード代数の生成元となります。

基礎にある結合代数は生成元の選択によって変わりませんが、クリフォード代数としての符号数は選択によって変わります。どのようにすればテンソル積の構造が解釈しやすいかを検討します。

$(3)$の構造を捉えるため模式化します。

$$
\begin{array}{c|ccccc}
\text{左因子}&\color{red}{i}&\color{red}{j}&\color{red}{k}&1&1 \\
&⊗&⊗&⊗&⊗&⊗ \\
\text{右因子}&i&i&\color{red}{i}&\color{red}{j}&\color{red}{k}
\end{array} \tag{7}
$$

&&&rem 赤字因子の並び
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

&&&rem 生成元の模式的選択
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

これら2乗の係数（赤字部分）を、生成元の**計量**（正規直交基底での計量の対角成分）、計量の値ごとの生成元の個数を**符号数**と呼びます。本記事では$1,-1$の順に符号数を数えます。

- 計量$1$が2個、計量$-1$が2個 → 符号数$(2,2)$

&&&rem 別の4元の選択と符号数
$(3)$の5個のうち、$1⊗j,1⊗k$の2乗は$-1$、残りの3個の2乗は$1$です。$1⊗j$か$1⊗k$を除いた4元は符号数$(3,1)$、それ以外を除いた4元は符号数$(2,2)$になります。同じ代数$\mathbb H⊗\mathbb H$が、選び方によって異なる符号数のクリフォード代数として表されます。
&&&

&&&rem 符号数によるクリフォード代数の指定
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

&&&rem 計量の反転
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

&&&rem 反転と追加生成元の計量
$e_1,e_2,e_3,e_4$の計量が$⊗i$によって反転して、追加された2個の生成元$\{1⊗j,\ 1⊗ij\}$の計量は$-1$です。
&&&

この構造を一般化します。

&&&fml クリフォード代数の$⊗\mathbb H$による拡張の一般化
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

この構成が実際にクリフォード代数の生成元を与えることを、$\mathbb H$と$\mathbb H'$に共通する形で示します。

&&&prop 生成元の拡張 [prop-ext]
$e_1,\cdots,e_n$を$\operatorname{Cl}_{p,q}(\mathbb R)$の生成元（$n=p+q$、$e_a^2=\varepsilon_a$）とする。実代数$B$の元$u,v$が$u^2=s,\ v^2=t\ (s,t\in\{1,-1\})$、$uv=-vu$を満たし、$1,u,v,uv$が$B$の基底をなすとする。このとき、次の$n+2$個の元は互いに反交換し、2乗は
$$
\begin{aligned}
E_a&=e_a⊗u, & E_a^2&=s\varepsilon_a \\
F&=1⊗v, & F^2&=t \\
G&=1⊗uv, & G^2&=-st
\end{aligned}
$$
となり、$\operatorname{Cl}_{p,q}(\mathbb R)⊗B$を生成する。
&&&

&&&prf
$E_a,E_b\ (a\ne b)$は、$E_aE_b=s\,e_ae_b⊗1$であり、$e_ae_b=-e_be_a$なので反交換する。$E_a$と$F$、$E_a$と$G$、$F$と$G$は、それぞれ右因子の組$(u,v),(u,uv),(v,uv)$が反交換し、左因子が可換なので反交換する。

2乗は次のとおり。
$$
E_a^2=e_a^2⊗u^2=s\varepsilon_a,\quad F^2=1⊗v^2=t,\quad
G^2=1⊗(uv)^2=1⊗(-u^2v^2)=-st
$$

生成することを示す。$FG=1⊗vuv=1⊗(-uv^2)=-t(1⊗u)$より$1⊗u$が得られ、$E_a(1⊗u)=s(e_a⊗1)$より$e_a⊗1$が得られる。$e_a⊗1$の積は$\operatorname{Cl}_{p,q}(\mathbb R)$の基底を与え、これらに$1⊗1,1⊗u,1⊗v,1⊗uv$を掛ければ、テンソル積の基底がすべて得られる。

したがって、$n+2$個の元は$\operatorname{Cl}_{p,q}(\mathbb R)⊗B$を生成する。
&&&

次元は$\dim(\operatorname{Cl}_{p,q}(\mathbb R)⊗B)=2^n\cdot4=2^{n+2}$なので、[[prop-cl-iso]]により、$E_a,F,G$は$\operatorname{Cl}_{p,q}(\mathbb R)⊗B$の生成元となります。$(u,v)$と$(s,t)$の選び方は次のとおりです。

- $\mathbb H$の$(u,v)=(i,j)$：$(s,t)=(-1,-1)$
- $\mathbb H'$の$(u,v)=(i,j)$：$(s,t)=(-1,1)$
- $\mathbb H'$の$(u,v)=(j,k)$：$(s,t)=(1,1)$

この結果を公式の形にまとめます。[[wiki-clif]]

&&&fml クリフォード代数の$⊗\mathbb H$による拡張 [fml-ext-h]
$$
\operatorname{Cl}_{p,q}(\mathbb R) ⊗ \mathbb H
\cong \operatorname{Cl}_{q,p+2}(\mathbb R)
$$
&&&

&&&rem 符号数の入れ替わり
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
分解型四元数は、共役との積によって定義された二次形式$N(x)=xx^*$により、計量が決まります。（クリフォード代数とは異なり、基底の2乗がそのまま計量とはなりません）

$$
\begin{aligned}
N(a+bi+cj+dk)
&=(a+bi+cj+dk)(a+bi+cj+dk)^* \\
&=(a+bi+cj+dk)(a-bi-cj-dk) \\
&=a^2+b^2-c^2-d^2
\end{aligned}
$$

$N$は4次元実空間全体の二次形式で、その符号数は$(2,2)$となり正と負の個数が等しくなります。このような符号数を**分解型**と呼び、代数名の由来となっています。クリフォード代数の符号数は、選択した生成元の空間の二次形式の符号数であり、$N$の符号数とは別のものです。
&&&

## クリフォード代数

虚数単位によって2乗の値が変わるため、どれをクリフォード代数の生成元として使うかで符号数が変わります。

1. $\operatorname{Cl}_{1,1}(\mathbb R)$: $\{i,j\},\{i,k\}$
2. $\operatorname{Cl}_{2,0}(\mathbb R)$: $\{j,k\}$

クリフォード代数としての性質は符号数にのみ依存するため、$\operatorname{Cl}_{1,1}(\mathbb R)$の生成元としては$\{i,j\}$のみを使用します。

&&&rem 分解型四元数の同型対応
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
$\mathbb H'⊗\mathbb H'$は、生成元の選び方によって$\operatorname{Cl}_{3,1}(\mathbb R)$とも$\operatorname{Cl}_{2,2}(\mathbb R)$とも同型である。したがって、次が成り立つ。
$$
\mathbb H⊗\mathbb H
\cong \operatorname{Cl}_{2,2}(\mathbb R)
\cong \mathbb H'⊗\mathbb H'
$$
&&&

## 一般化と公式

この構造を一般化します。

&&&fml クリフォード代数の$⊗\mathbb H'$による拡張の一般化
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

この構成が生成元を与えることは、[[prop-ext]]の$(u,v)=(i,j),(j,k)$の場合です。この結果を公式の形にまとめます。[[wiki-clif]]

&&&fml クリフォード代数の$⊗\mathbb H'$による拡張 [fml-ext-hs]
$$
\begin{aligned}
\operatorname{Cl}_{p,q}(\mathbb R) ⊗ \mathbb H'
&\cong \operatorname{Cl}_{q+2,p}(\mathbb R) \\
&\cong \operatorname{Cl}_{p+1,q+1}(\mathbb R)
\end{aligned}
$$
&&&

&&&rem 計量変化と符号数の増分
中辺の$\operatorname{Cl}_{q+2,p}(\mathbb R)$は$p,q$の位置が入れ替わります。これは$⊗i$による計量の反転に由来します。
右辺の$\operatorname{Cl}_{p+1,q+1}(\mathbb R)$は$p,q$の位置が維持されます。これは$⊗j$によって計量が変化しないことに由来します。

符号数の増分$(2,0),(1,1)$は$\mathbb H' \cong \operatorname{Cl}_{2,0}(\mathbb R) \cong \operatorname{Cl}_{1,1}(\mathbb R)$に由来します。
&&&

## $\mathbb H'⊗\mathbb H,\ \mathbb H⊗\mathbb H'$

[[fml-ext-h]]より

$$
\begin{aligned}
\mathbb H'⊗\mathbb H

&\cong \operatorname{Cl}_{2,0}(\mathbb R)⊗\mathbb H
\cong \operatorname{Cl}_{0,4}(\mathbb R) \\

&\cong \operatorname{Cl}_{1,1}(\mathbb R)⊗\mathbb H
\cong \operatorname{Cl}_{1,3}(\mathbb R) \\
\end{aligned}
$$

[[fml-ext-hs]]より

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

&&& 四元数・分解型四元数と同型対応
四元数および分解型四元数は、低次元のクリフォード代数と同型です。
$$
\begin{alignedat}{2}
&\mathbb H &&\cong \operatorname{Cl}_{0,2}(\mathbb R) \\
&\mathbb H' &&\cong \operatorname{Cl}_{2,0}(\mathbb R) \cong \operatorname{Cl}_{1,1}(\mathbb R)
\end{alignedat}
$$
&&&

&&& テンソル積によるクリフォード代数の拡張公式
$$
\begin{alignedat}{2}
\operatorname{Cl}_{p,q}(\mathbb R) &⊗ \mathbb H  &&\cong \operatorname{Cl}_{q,p+2}(\mathbb R) \\
\operatorname{Cl}_{p,q}(\mathbb R) &⊗ \mathbb H' &&\cong \operatorname{Cl}_{q+2,p}(\mathbb R) \cong \operatorname{Cl}_{p+1,q+1}(\mathbb R)
\end{alignedat}
$$
$$
\begin{alignedat}{6}
\mathbb{H} &⊗ \mathbb{H} &&\cong \mathbb{H}' ⊗ \mathbb{H}' &&\cong \operatorname{Cl}_{3,1}(\mathbb R) \cong \operatorname{Cl}_{2,2}(\mathbb R) \\
\mathbb{H} &⊗ \mathbb{H}' &&\cong \mathbb{H}' ⊗ \mathbb{H} &&\cong \operatorname{Cl}_{1,3}(\mathbb R) \cong \operatorname{Cl}_{0,4}(\mathbb R) \cong \operatorname{Cl}_{4,0}(\mathbb R)
\end{alignedat}
$$
&&&

&&& テンソル積とクリフォード代数の対応表
$p+q=2,4$の場合を示します。
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
