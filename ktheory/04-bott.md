球面のK群を貼り合わせ関数の安定なホモトピー類として捉え、クリフォード加群から貼り合わせ関数を作って、加群の群を制限の像で割った群を分類表から計算し、実K群の8周期の表と比べます。

# 概要

以前の記事では、球面のK群が実では周期8、複素では周期2で繰り返すこと（ボット周期性）を主張として紹介し、その周期がクリフォード代数の分類表の周期と一致することを見ました。生成元は、$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群の上のパラベクトル$x_0+\sum_lx_le_l$を貼り合わせ関数として作れること、さらにK群そのものが「$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群の群を、$\operatorname{Cl}_{0,k}(\mathbb R)$の加群の制限で割る」ことで計算できること（アティヤ＝ボット＝シャピロの定理）も、主張として述べるに留めました。分類表を導いた記事で「K理論におけるボット周期性定理と対応することが知られています」と述べた対応の中身は、この定理にあります。[[7shi-kth1]][[7shi-clif1]]

前回の記事では、球面上のベクトル束を赤道での貼り合わせ関数で表し、$S^2$上の複素直線束が回転数で分類されることを示しました。一方、階数が高い束については、行列式の回転数が同型で変わらないことまでを示し、「同じ階数で行列式の回転数が等しければ同型」は主張に留めました。四元数・八元数のホップ束が自明束を足しても自明にならないことも扱いませんでした。[[7shi-kth3]]

本記事で扱う問いは次のとおりです。

> 球面のK群の元を貼り合わせ関数の言葉で表すと、群の演算は何になるか。どの貼り合わせ関数がK群の$0$を与え、残りはクリフォード代数の分類表からどう読めるか。

前提は次のとおりです。

- **以前の記事から**：K群と簡約K群の定義、コンパクトな空間では同値関係の$G$を自明束に取り替えられること、球面のK群の表、パラベクトルによる貼り合わせ関数、クリフォード代数の分類表の第1行から読んだ既約加群の次元$a_k$と、直和型の2つの既約加群が擬スカラーの符号で区別されること[[7shi-kth1]][[7shi-kth2]][[7shi-clif1]]
- **前回の記事から**：貼り合わせの向き（$D_-$の$(\boldsymbol x,v)$を$D_+$の$(\boldsymbol x,g(\boldsymbol x)v)$と同一視する）、ホモトピックな貼り合わせ関数が同型な束を与えること（複素では逆も成り立つこと）、回転数とその性質、$H^n=E_{z^n}$、ホップ束の貼り合わせ関数が単位数の左からの積$L_c$であること[[7shi-kth3]]

本記事は次の順に進みます。

1. 自明束を足すことを貼り合わせ関数の言葉に移し、複素の簡約K群の元を貼り合わせ関数の安定なホモトピー類として表す。平面の回転を使って、K群の和が貼り合わせ関数の積にあたることを示す
2. $S^2$の場合に、行列式の回転数がK群から整数への全射な準同型を与えることと、直線束の直和の類が回転数だけで決まることを示す
3. クリフォード加群から貼り合わせ関数を作り、1つ多い生成元の作用まで延長できる加群は自明な束を与えることを示す
4. 加群の群を延長できる加群の制限で割った群を、分類表から計算して$\mathbb Z_2,\mathbb Z_2,0,\mathbb Z,0,0,0,\mathbb Z$を得る。代数の側の8周期性から、この表が周期8で繰り返すことを見る

証明するのは、安定なホモトピー類による表示、回転による変形とその帰結、$S^2$上の直線束の直和の分類、延長できる加群が自明な束を与えること、制限の計算です。有限次元の加群が既約加群の直和に一意に分かれること、$\operatorname{SL}(n,\mathbb C)$の中の閉曲線が1点に縮められることは主張として使います。アティヤ＝ボット＝シャピロの定理（加群から作った束がK群全体を与え、それ以外に関係がないこと）とボット周期性定理の証明、四元数・八元数のホップ束が自明束を足しても自明にならないことの証明は、本記事では扱いません。

# 球面のK群と貼り合わせ関数

## 自明束を足す

以前の記事では、K群の元$[E]-[F]$を、ベクトル束の組$(E,F)$を$E\oplus F'\oplus G\cong E'\oplus F\oplus G$で同一視したものとして定義しました。球面のようなコンパクトな空間では$G$を自明束に取り替えられ、簡約K群の元は階数$n$の束$E$を使って$[E]-n$と書けます。前回の記事では、球面上のベクトル束を貼り合わせ関数で表しました。この2つを合わせると、K群の元とその間の等式を、貼り合わせ関数だけで書けるはずです。[[7shi-kth1]][[7shi-kth3]]

鍵になるのは、自明束を足す操作です。$E_g\oplus\underline{\mathbb C}^a$の貼り合わせ関数は、$g$の右下に単位行列を並べた$\operatorname{diag}(g,I_a)$です。以下、複素ベクトル束を中心に述べ、実ベクトル束との違いはその都度断ります。

&&&def 安定なホモトピー
$g:S^{k-1}\to\operatorname{GL}(n,\mathbb C)$と$g':S^{k-1}\to\operatorname{GL}(n',\mathbb C)$について、$n+a=n'+b$を満たす$a,b\ge0$で
$$
\operatorname{diag}(g,I_a)\simeq\operatorname{diag}(g',I_b)
$$
となるものがあるとき、$g$と$g'$は**安定にホモトピック**であると言います。
&&&

自明束を足して同型になることを、単位行列を並べてホモトピックになることに置き換えたものです。$\operatorname{diag}(g,I_a)\simeq\operatorname{diag}(g',I_b)$なら、両方にさらに$I_c$を並べても同じホモトピーが使えるので、$a,b$はいくらでも大きく取り直せます。

&&&prop 簡約K群と安定なホモトピー [prop-stable]
$S^k$上の複素ベクトル束$E_g,E_{g'}$（階数$n,n'$）について、$\tilde K(S^k)$の中で$[E_g]-n=[E_{g'}]-n'$となることと、$g$と$g'$が安定にホモトピックであることは同値です。
&&&

&&&prf
$[E_g]-n=[E_{g'}]-n'$は、K群の定義と$G$を自明束に取り替えられることから、ある$N$で$E_g\oplus\underline{\mathbb C}^{n'+N}\cong E_{g'}\oplus\underline{\mathbb C}^{n+N}$となることと同値である。両辺の貼り合わせ関数は$\operatorname{diag}(g,I_{n'+N})$と$\operatorname{diag}(g',I_{n+N})$で、前回の記事の定理「複素ベクトル束と貼り合わせ関数」により、同型はこの2つがホモトピックであることと同値である。したがって左辺が成り立てば安定にホモトピックである。[[7shi-kth3]]

逆に$n+a=n'+b$で$\operatorname{diag}(g,I_a)\simeq\operatorname{diag}(g',I_b)$なら、$E_g\oplus\underline{\mathbb C}^a\cong E_{g'}\oplus\underline{\mathbb C}^b$である。両辺に$\underline{\mathbb C}^{n+b}$を足すと、$a+n+b=n'+2b$より$E_g\oplus\underline{\mathbb C}^{n'+2b}\cong E_{g'}\oplus\underline{\mathbb C}^{n+2b}$となり、$N=2b$として上の条件を満たす。
&&&

球面上のベクトル束はすべて貼り合わせで得られるので（前回の記事で主張として使いました）、$\tilde K(S^k)$の元はすべて$[E_g]-n$の形をしています。[[prop-stable]]により、$\tilde K(S^k)$は貼り合わせ関数を安定なホモトピーで分類したものと同じです。和は直和なので、$([E_g]-n)+([E_h]-m)=[E_{\operatorname{diag}(g,h)}]-(n+m)$です。[[7shi-kth3]]

実ベクトル束でも、安定にホモトピックな貼り合わせ関数が$\widetilde{KO}(S^k)$の同じ元を与えることは同様に成り立ちます（前回の記事の命題「貼り合わせ関数の変形」は実でも成り立つため）。逆向きには、前回の記事で見た定数行列の違いが残ります。以下の議論で実ベクトル束に使うのは、ホモトピックなら同型という向きだけです。[[7shi-kth3]]

## 平面の回転による変形

和$\operatorname{diag}(g,h)$は、2つのブロックを並べただけで、$g$と$h$は混ざりません。ところが、平面の回転を使うと、2つのブロックを1つにまとめられます。以前の記事では、メビウスの帯2枚の貼り合わせ行列$-I_2$を回転$R(φ)$で$I_2$に変形しました。同じ回転を、ブロックごとに行います。[[7shi-kth1]]

$n$次の単位行列$I$を成分とするブロックの回転行列を、次のように置きます。

$$
R(t)=\begin{pmatrix}\cos t\,I&-\sin t\,I\\\sin t\,I&\cos t\,I\end{pmatrix}
$$

$R(0)$は単位行列で、$R(\pi/2)$は2つのブロックを入れ替えます。

&&&lem 回転による変形 [lem-rotation]
$g,h:S^{k-1}\to\operatorname{GL}(n,\mathbb F)$（$\mathbb F=\mathbb R$または$\mathbb C$）について、次のホモトピーがあります。

$$
\operatorname{diag}(g,h)\simeq\operatorname{diag}(gh,I)
$$
&&&

&&&prf
$0\le t\le\pi/2$について$G_t=\operatorname{diag}(g,I)\,R(t)\,\operatorname{diag}(I,h)\,R(t)^{-1}$と置く。各因子は可逆なので$G_t$は可逆で、$(t,\boldsymbol x)$について連続である。$t=0$では$G_0=\operatorname{diag}(g,I)\operatorname{diag}(I,h)=\operatorname{diag}(g,h)$である。$t=\pi/2$では、次の計算から$R(\pi/2)\operatorname{diag}(I,h)R(\pi/2)^{-1}=\operatorname{diag}(h,I)$となる。

$$
\begin{pmatrix}0&-I\\I&0\end{pmatrix}\begin{pmatrix}I&0\\0&h\end{pmatrix}\begin{pmatrix}0&I\\-I&0\end{pmatrix}=\begin{pmatrix}0&-h\\I&0\end{pmatrix}\begin{pmatrix}0&I\\-I&0\end{pmatrix}=\begin{pmatrix}h&0\\0&I\end{pmatrix}
$$

したがって$G_{\pi/2}=\operatorname{diag}(g,I)\operatorname{diag}(h,I)=\operatorname{diag}(gh,I)$である。
&&&

$R(t)$は実の行列で、行列式は$1$です。$h$の入ったブロックを、回転で$g$の隣へ運んで掛け合わせています。

&&&cor K群の和と貼り合わせ関数の積 [cor-product]
$g,h:S^{k-1}\to\operatorname{GL}(n,\mathbb C)$について、次が成り立ちます。

$$
([E_g]-n)+([E_h]-n)=[E_{gh}]-n,\qquad-([E_g]-n)=[E_{g^{-1}}]-n
$$

実ベクトル束の$\widetilde{KO}(S^k)$でも同じです。
&&&

&&&prf
[[lem-rotation]]より$E_g\oplus E_h\cong E_{gh}\oplus\underline{\mathbb C}^n$なので、$[E_g]+[E_h]=[E_{gh}]+n$である。両辺から$2n$を引けば1つ目の式を得る。$h=g^{-1}$とすると$gh=I$で、$([E_g]-n)+([E_{g^{-1}}]-n)=[E_I]-n=0$となり、2つ目の式を得る。実ベクトル束でも、[[lem-rotation]]と前回の記事の命題「貼り合わせ関数の変形」から同じ計算ができる。[[7shi-kth3]]
&&&

安定化してK群で見ると、ブロックを並べる和が、貼り合わせ関数の積に置き換わります。逆元は逆行列を貼り合わせ関数にした束です。

# S²のK群

## 直線束の和

$S^2$では赤道が円周で、前回の記事で複素直線束を$H^m=E_{z^m}$（$m\in\mathbb Z$）に分類しました。直線束の貼り合わせ関数は$\mathbb C^\times$に値を取るので、[[lem-rotation]]を$n=1$で使えます。[[7shi-kth3]]

&&&ex ホップ束の2つの和 [ex-hopf-sums]
$g=h=z$とすると、[[lem-rotation]]は$\operatorname{diag}(z,z)\simeq\operatorname{diag}(z^2,1)$を与えます。したがって次のようになります。

$$
H\oplus H\cong H^2\oplus\underline{\mathbb C}
$$

$g=z$、$h=z^{-1}$とすると$\operatorname{diag}(z,z^{-1})\simeq I_2$で、$H\oplus H^{-1}\cong\underline{\mathbb C}^2$です。
&&&

1つ目の同型の両辺は、どちらも行列式が$z^2$で、行列式の回転数は$2$です。前回の記事では、同じ階数で行列式の回転数が等しい束は同型であることを主張に留めました。[[ex-hopf-sums]]は、その具体例を回転で示したものにあたります。[[7shi-kth3]]

同じ議論を繰り返すと、直線束の直和は、行列式の回転数だけで決まります。

&&&prop 直線束の直和 [prop-line-sums]
$S^2$上の複素直線束$L_1,\dots,L_n$の回転数を$m_1,\dots,m_n$とし、$m=m_1+\dots+m_n$と置きます。このとき次が成り立ちます。

$$
L_1\oplus\dots\oplus L_n\cong H^m\oplus\underline{\mathbb C}^{n-1},\qquad[H^m]-1=m([H]-1)
$$
&&&

&&&prf
前回の記事の定理「複素直線束の分類」により$L_j\cong H^{m_j}$だから、貼り合わせ関数は$\operatorname{diag}(z^{m_1},\dots,z^{m_n})$としてよい。[[lem-rotation]]を右下の2つのブロックに順に使うと、$\operatorname{diag}(z^{m_{n-1}},z^{m_n})\simeq\operatorname{diag}(z^{m_{n-1}+m_n},1)$である。1になった成分は右端に並ぶので、残りの左側のブロックに同じ変形を繰り返せば$\operatorname{diag}(z^m,1,\dots,1)$に達する。[[7shi-kth3]]

後半は[[cor-product]]を$g=z^{m-1}$、$h=z$に使えば$[H^m]-1=([H^{m-1}]-1)+([H]-1)$となり、$m$についての帰納法で$m\ge0$の場合が従う。$m<0$の場合は[[cor-product]]の逆元の式から$[H^{-m}]-1=-([H^m]-1)$を使う。
&&&

## 行列式の回転数による準同型

直線束の直和では、K群の元が行列式の回転数で決まることが分かりました。行列式の回転数そのものは、どの束についても定義できます。これがK群の上の関数になるかを確かめます。

&&&thm S²のK群と行列式の回転数 [thm-s2]
$\tilde K(S^2)$の元$[E_g]-n$に$\operatorname{wind}(\det g)$を対応させる写像$φ:\tilde K(S^2)\to\mathbb Z$は、全射な準同型です。直線束の直和から作られる元$[L_1\oplus\dots\oplus L_n]-n$は$(m_1+\dots+m_n)([H]-1)$に等しく、$φ$の値で決まります。
&&&

&&&prf
$[E_g]-n=[E_{g'}]-n'$なら、[[prop-stable]]により$\operatorname{diag}(g,I_a)\simeq\operatorname{diag}(g',I_b)$となる$a,b$がある。ホモトピーの行列式を取れば$\det g\simeq\det g'$（$\mathbb C^\times$の中で）となり、前回の記事の命題「回転数の性質」の3から回転数は等しい。よって$φ$は矛盾なく定まる。$\det\operatorname{diag}(g,h)=\det g\det h$と命題「回転数の性質」の1から、$φ$は和を和に移す。$φ([H^m]-1)=m$なので全射である。後半は[[prop-line-sums]]である。[[7shi-kth3]]
&&&

$φ$が単射であること、すなわち$\tilde K(S^2)\cong\mathbb Z$を言うには、どの束も直線束の直和と安定同値であることが要ります。これには次の事実を使います。

$\operatorname{SL}(n,\mathbb C)$（行列式が$1$の複素行列の全体）の中の閉曲線は、$\operatorname{SL}(n,\mathbb C)$の中で$I$に連続に縮められることが知られています（$\operatorname{SL}(n,\mathbb C)$が単連結であること）。単連結性の証明は本記事では扱いません。これを認めると、前回の記事で主張に留めた分類が得られます。[[7shi-kth3]]

&&&prop S²上の複素ベクトル束 [prop-s2-classify]
$\operatorname{SL}(n,\mathbb C)$の単連結性を認めると、$g:S^1\to\operatorname{GL}(n,\mathbb C)$について、$m=\operatorname{wind}(\det g)$として次が成り立ちます。

$$
E_g\cong H^m\oplus\underline{\mathbb C}^{n-1}
$$

とくに、同じ階数の$S^2$上の複素ベクトル束は、行列式の回転数が等しければ同型です。また、[[thm-s2]]の$φ$は同型で、$\tilde K(S^2)\cong\mathbb Z$の生成元は$[H]-1$です。
&&&

&&&prf
$h=g\cdot\operatorname{diag}(\det g,1,\dots,1)^{-1}$は$\operatorname{SL}(n,\mathbb C)$に値を取る閉曲線だから、単連結性により$h\simeq I$である。したがって$g=h\cdot\operatorname{diag}(\det g,1,\dots,1)\simeq\operatorname{diag}(\det g,1,\dots,1)$となる。前回の記事の命題「回転数の性質」の1と4から$\det g/z^m\simeq1$、すなわち$\det g\simeq z^m$（$\mathbb C^\times$の中で）なので、$g\simeq\operatorname{diag}(z^m,1,\dots,1)$であり、$E_g\cong H^m\oplus\underline{\mathbb C}^{n-1}$である。後半の同型は、行列式の回転数が等しい2つの束がどちらも同じ$H^m\oplus\underline{\mathbb C}^{n-1}$に同型になることから従う。[[7shi-kth3]]

$\tilde K(S^2)$の元は$[E_g]-n=[H^m]-1=m([H]-1)$（[[prop-line-sums]]）と書けるので、$φ$の値$m$が$0$なら元も$0$であり、$φ$は単射である。[[thm-s2]]の全射性と合わせて$φ$は同型である。
&&&

&&&ex 接束の類
前回の記事では、複素直線束として$TS^2\cong H^{-2}$であることを見ました。[[prop-line-sums]]から、$\tilde K(S^2)$の中で$[TS^2]-1=-2([H]-1)$です。実ベクトル束としては$TS^2\oplus\underline{\mathbb R}\cong\underline{\mathbb R}^3$で、$\widetilde{KO}(S^2)$の中では$[TS^2]-2=0$でした。同じ束が、複素のK群では生成元の$-2$倍として残り、実のK群では消えます。[[7shi-kth3]]
&&&

&&&rem テンソル積と積
前回の記事では、複素直線束のテンソル積を定義し、$E_g\otimes E_h=E_{gh}$を見ました。一般のベクトル束にもテンソル積が定義でき、K群はそれを積として環になります。[[ex-hopf-sums]]の$H\oplus H\cong H\otimes H\oplus\underline{\mathbb C}$は、この環の中で$[H]^2-2[H]+1=0$、すなわち$([H]-1)^2=0$と書けます。一般のテンソル積とK群の環の構造には、本記事では立ち入りません。[[7shi-kth3]]
&&&

# クリフォード加群による生成元

## 加群が与える束

$S^2$の生成元は、ホップ束$H$から作った$[H]-1$でした。以前の記事では、一般の$S^k$でも、$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群の上のパラベクトルを貼り合わせ関数にすると生成元が得られると述べました。どの加群がどの元を与えるかを調べるために、加群ごとに束を定めます。[[7shi-kth1]]

以前の記事にならい、$\operatorname{Cl}_{0,k-1}(\mathbb R)$の生成元$e_l$（$l=1,\dots,k-1$）の加群$W$への作用を$J_l$と書きます。$J_l^2=-I$、$l\ne m$なら$J_lJ_m=-J_mJ_l$です。[[7shi-kth2]]

&&&def 加群が与える束 [def-module-bundle]
$\operatorname{Cl}_{0,k-1}(\mathbb R)$の有限次元の加群$W$に対して、赤道$S^{k-1}$の点$\boldsymbol x=(x_0,\dots,x_{k-1})$で

$$
g_W(\boldsymbol x)=x_0I+x_1J_1+\dots+x_{k-1}J_{k-1}
$$

と置き、$g_W:S^{k-1}\to\operatorname{GL}(W)$を貼り合わせ関数とする$S^k$上の実ベクトル束を$E_W$と書きます。$W$が複素ベクトル空間で、$J_l$が複素線形なら、$E_W$は複素ベクトル束です。
&&&

$g_W(\boldsymbol x)$が可逆であることは、以前の記事で見たとおり、$\boldsymbol v=\sum_lx_lJ_l$について$\boldsymbol v^2=-\sum_lx_l^2I$となることから従います。$(x_0I+\boldsymbol v)(x_0I-\boldsymbol v)=|\boldsymbol x|^2I=I$です。[[7shi-kth1]]

複素の場合の加群は、生成元が複素線形に作用するものです。これは複素クリフォード代数$\operatorname{Cl}_{k-1}(\mathbb C)\cong\operatorname{Cl}_{0,k-1}(\mathbb R)\otimes_{\mathbb R}\mathbb C$の加群と同じものです。以下、$\dim W$は、実の加群では実次元、複素の加群では複素次元を表します。[[7shi-clif1]]

&&&prop 加群の同型と直和 [prop-module-sum]
$W\cong W'$（加群として同型）なら$E_W\cong E_{W'}$です。また、$E_{W\oplus W'}=E_W\oplus E_{W'}$です。
&&&

&&&prf
加群の同型$P:W\to W'$は$PJ_l=J'_lP$を満たすので、$g_{W'}=Pg_WP^{-1}$である。前回の記事の同型の条件$g'=A_+gA_-^{-1}$で$A_\pm=P$（定数）と取れば$E_W\cong E_{W'}$となる。直和の加群では$J_l$がブロック対角$\operatorname{diag}(J_l,J'_l)$で作用するので、$g_{W\oplus W'}=\operatorname{diag}(g_W,g_{W'})$である。[[7shi-kth3]]
&&&

## 延長できる加群

$W$の上に、もう1つの生成元の作用$J_k$（$J_k^2=-I$で、$J_1,\dots,J_{k-1}$と反交換する）があるとします。このとき$W$は$\operatorname{Cl}_{0,k}(\mathbb R)$の加群で、その作用を$\operatorname{Cl}_{0,k-1}(\mathbb R)$に制限したものが元の$W$です。このような$W$を、$\operatorname{Cl}_{0,k}(\mathbb R)$の加群に**延長できる**と言います。

&&&prop 延長できる加群は自明な束を与える [prop-extend]
$W$が$\operatorname{Cl}_{0,k}(\mathbb R)$の加群に延長できるなら、$E_W$は自明束です。複素の場合も同様です。
&&&

&&&prf
$0\le t\le\pi/2$について次のように置く。

$$
g_t(\boldsymbol x)=\cos t\,g_W(\boldsymbol x)+\sin t\,J_k
$$

$a=x_0\cos t$、$\boldsymbol u=\cos t\sum_lx_lJ_l+\sin t\,J_k$とすると$g_t=aI+\boldsymbol u$である。$J_1,\dots,J_k$は互いに反交換して2乗が$-I$なので、$\boldsymbol u^2=-(\cos^2t\sum_lx_l^2+\sin^2t)I$となる。したがって

$$
(aI+\boldsymbol u)(aI-\boldsymbol u)=a^2I-\boldsymbol u^2=\Bigl(\cos^2t\sum_{i=0}^{k-1}x_i^2+\sin^2t\Bigr)I=I
$$

であり、$g_t(\boldsymbol x)$は可逆である。$g_0=g_W$、$g_{\pi/2}=J_k$（定数）だから、前回の記事の命題「貼り合わせ関数の変形」により$E_W\cong E_{J_k}$である。定数$J_k$による貼り合わせは、同型の条件で$A_+=J_k^{-1}$、$A_-=I$と取れば$E_I$、すなわち自明束に同型である。[[7shi-kth3]]
&&&

$g_t$は、赤道の点$\boldsymbol x\in S^{k-1}$を1つ大きい球面$S^k$の中で北極の方向$(0,\dots,0,1)$へ動かし、その点のパラベクトルを取ったものです。新しい生成元$J_k$の方向が空いているので、赤道全体を1点に寄せられます。

[[prop-module-sum]]と[[prop-extend]]から、K群の元$[E_W]-\dim W$は加群の直和について足し算になり、延長できる加群では$0$になります。低い次元で確かめます。

&&&ex メビウスの帯（k=1）
$\operatorname{Cl}_{0,0}(\mathbb R)=\mathbb R$で、生成元はありません。$W=\mathbb R$なら$g_W(x_0)=x_0=\pm1$で、$E_W$はメビウスの帯$M$です。$W=\mathbb R^2$に$J_1=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$を作用させると$\operatorname{Cl}_{0,1}(\mathbb R)\cong\mathbb C$の加群になるので、$E_{\mathbb R^2}=M\oplus M$は自明です。[[prop-extend]]の$g_t$は、$x_0=1$では$R(t)$、$x_0=-1$では$R(\pi-t)$で、以前の記事でメビウスの帯2枚をほどいた回転と同じものです。[[7shi-kth1]]
&&&

&&&ex 複素数のホップ束（k=2）
$\operatorname{Cl}_{0,1}(\mathbb R)\cong\mathbb C$の加群$W=\mathbb C$（$J_1=i$）では$g_W(x_0,x_1)=x_0+x_1i$で、$E_W$はホップ束$H$を実ベクトル束と見たものです。$\mathbb C^2=\mathbb H$は$\mathbf i,\mathbf j$の左からの積で$\operatorname{Cl}_{0,2}(\mathbb R)\cong\mathbb H$の加群になるので、$E_{\mathbb C^2}=H\oplus H$は実ベクトル束として自明で、$\widetilde{KO}(S^2)$の中で$2([H]-2)=0$です。

複素の加群としては、$\operatorname{Cl}_1(\mathbb C)\cong2\mathbb C$に既約加群が2つあり、$J_1=i$のもの$W$と$J_1=-i$のもの$\overline W$です。赤道の点を$z=x_0+x_1i$と書くと$g_W=z$、$g_{\overline W}=\bar z=z^{-1}$なので、$E_W=H$、$E_{\overline W}=H^{-1}$です。$W\oplus\overline W=\mathbb C^2$は$J_1=\operatorname{diag}(i,-i)$、$J_2=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$で$\operatorname{Cl}_2(\mathbb C)\cong M_2(\mathbb C)$の加群に延長できるので、$H\oplus H^{-1}$は自明です。これは[[ex-hopf-sums]]で回転から得た同型と一致します。
&&&

## 直和型の2つの既約加群

$\operatorname{Cl}_{0,3}(\mathbb R)\cong2\mathbb H$と$\operatorname{Cl}_{0,7}(\mathbb R)\cong2M_8(\mathbb R)$は直和型で、既約加群が2つあります。以前の記事では、四元数の$\mathbf i,\mathbf j,\mathbf k$を左から掛ける作用と、符号を反転した$-\mathbf i,-\mathbf j,-\mathbf k$を掛ける作用が、擬スカラーの符号で区別される2つの既約加群であることを見ました。八元数でも同じです。[[7shi-kth2]]

一般に、既約加群$W$の生成元の作用$J_l$の符号をすべて反転した$-J_l$も同じクリフォード関係を満たし、加群になります。これを$\overline W$と書きます。擬スカラー$ω=J_1\cdots J_{k-1}$は$(-1)^{k-1}ω$に変わります。$\operatorname{Cl}_{0,3}(\mathbb R)$や$\operatorname{Cl}_{0,7}(\mathbb R)$のような直和型（$k-1\equiv3,7\pmod8$）では、既約加群の上で擬スカラーが$\pm I$として作用し、その符号が反転するので、$\overline W$は$W$と異なる既約加群です。$\operatorname{Cl}_{0,3}(\mathbb R)$と$\operatorname{Cl}_{0,7}(\mathbb R)$の2つの既約加群は$W$と$\overline W$で尽くされます。

&&&ex 四元数・八元数のホップ束の共役（k=4, 8）
$W=\mathbb H$（$J_l=L_{\mathbf i},L_{\mathbf j},L_{\mathbf k}$）では$g_W(c)=L_c$で、$E_W$は前回の記事の四元数のホップ束です。$\overline W$では、$c=x_0+x_1\mathbf i+x_2\mathbf j+x_3\mathbf k$の共役$c^*$を使って$g_{\overline W}(c)=L_{c^*}$です。八元数の$W=\mathbb O$と$\overline W$でも、$g_W(c)=L_c$、$g_{\overline W}(c)=L_{c^*}$です。[[7shi-kth3]]

$W\oplus\overline W$は、次の作用で1つ多い生成元の加群に延長できます。

$$
J_l\mapsto\begin{pmatrix}J_l&0\\0&-J_l\end{pmatrix}\ (l=1,\dots,k-1),\qquad J_k=\begin{pmatrix}0&-I\\I&0\end{pmatrix}
$$

$J_k^2=-I$で、ブロックを計算すると$\operatorname{diag}(J_l,-J_l)J_k=-J_k\operatorname{diag}(J_l,-J_l)$です。$k=8$の場合は、以前の記事で$\operatorname{Cl}_{0,8}(\mathbb R)$の既約加群を作った倍加そのものです。[[prop-extend]]から$E_W\oplus E_{\overline W}$は自明で、K群の中で次が成り立ちます。[[7shi-kth2]]

$$
[E_{\overline W}]-\dim W=-([E_W]-\dim W)
$$

共役を取る貼り合わせ関数は、符号が逆の元を与えます。$S^2$の$H^{-1}$と$H$の関係と同じです。
&&&

# 分類表からの計算

## 加群の群を制限で割る

[[prop-module-sum]]と[[prop-extend]]により、$W\mapsto[E_W]-\dim W$は加群の直和を和に移し、延長できる加群を$0$に移します。そこで、加群にも形式的な差を導入して群を作り、延長できる加群を$0$と見なすと、K群への写像が得られます。

$\operatorname{Cl}_{0,j}(\mathbb R)$の有限次元の加群の同型類は、直和について可換なモノイドをなします。以前の記事で1点のK群$KO(\mathrm{pt})\cong\mathbb Z$を作ったのと同じく、形式的な差で群にしたものを$\mathfrak M_j$と書きます。有限次元の加群は既約加群の直和に分かれ、各既約加群が何個現れるかは加群によって一意に決まることが知られています。これを認めると、$\mathfrak M_j$は既約加群の類を基底とする自由アーベル群（既約加群ごとの個数を整数で記録する群）で、単純型なら$\mathbb Z$、直和型なら$\mathbb Z^2$です。[[7shi-kth1]]

$\operatorname{Cl}_{0,k}(\mathbb R)$の加群を、$J_k$を忘れて$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群と見る操作を**制限**と呼び、それが定める準同型を$r:\mathfrak M_k\to\mathfrak M_{k-1}$と書きます。$r$の像が、延長できる加群の差で書ける元の全体です。

&&&def 制限で割った群
制限の像で割った次の群を考えます。

$$
\mathcal A_k=\mathfrak M_{k-1}/r(\mathfrak M_k)
$$
複素クリフォード代数の加群から同じように作った群を$\mathcal A_k^{\mathbb C}$と書きます。
&&&

[[prop-extend]]により、$[W]\mapsto[E_W]-\dim W$は準同型$\mathcal A_k\to\widetilde{KO}(S^k)$（複素では$\mathcal A_k^{\mathbb C}\to\tilde K(S^k)$）を定めます。以前の記事で述べたアティヤ＝ボット＝シャピロの定理は、これが同型であることを主張しています。[[7shi-kth1]]

## 制限による分かれ方

$\mathcal A_k$を計算するには、$\operatorname{Cl}_{0,k}(\mathbb R)$の各既約加群を制限したとき、$\operatorname{Cl}_{0,k-1}(\mathbb R)$の既約加群がいくつ現れるかが分かればよいです。$\operatorname{Cl}_{0,k-1}(\mathbb R)$が単純型なら既約加群は1つなので、次元の比$a_k/a_{k-1}$がそのまま個数です。直和型なら、2つの既約加群$W,\overline W$の個数を区別する必要があります。

&&&lem 直和型への制限 [lem-restrict]
$\operatorname{Cl}_{0,k-1}(\mathbb R)$が直和型（$k-1\equiv3,7\pmod8$）とします。$\operatorname{Cl}_{0,k}(\mathbb R)$の加群を$\operatorname{Cl}_{0,k-1}(\mathbb R)$に制限すると、$W$と$\overline W$は同じ個数ずつ現れます。
&&&

&&&prf
$ω=J_1\cdots J_{k-1}$とする。$k-1$は奇数なので、$J_l$（$l\le k-1$）は$ω$の自分以外の$k-2$個（偶数個）の因子と反交換し、自分自身とは可換だから、$J_lω=ωJ_l$である。一方$J_k$は$k-1$個（奇数個）の因子すべてと反交換するので、$J_kω=-ωJ_k$である。

既約加群$W$と$\overline W$の上で$ω$はそれぞれ$\pm1$（符号は互いに逆）として作用するので、制限した加群の上で$ω^2=I$であり、空間は$ω=1$と$ω=-1$の固有空間$V_+,V_-$の直和に分かれる。$J_l$（$l\le k-1$）は$ω$と可換なので$V_\pm$を保ち、$V_\pm$はそれぞれ一方の既約加群だけの直和である。$J_k$は$ω$と反交換するので$V_+$を$V_-$に、$V_-$を$V_+$に移し、可逆だから$\dim V_+=\dim V_-$である。$W$と$\overline W$の次元は等しいので、個数も等しい。
&&&

$\operatorname{Cl}_{0,k}(\mathbb R)$が直和型（$k\equiv3,7\pmod8$）で制限先が単純型の場合は、2つの既約加群のどちらも、次元の比の個数だけ同じ既約加群を含みます。

## 実K群の表

以前の記事の既約加群の次元$a_k$と、分類表の型から、$\mathcal A_k$を計算します。$r$の行列の列は$\operatorname{Cl}_{0,k}(\mathbb R)$の既約加群、行は$\operatorname{Cl}_{0,k-1}(\mathbb R)$の既約加群に対応し、成分は制限に現れる個数です。[[7shi-kth2]]

$$
\begin{array}{c|cccc|c|c}
k&\operatorname{Cl}_{0,k-1}(\mathbb R)&a_{k-1}&\operatorname{Cl}_{0,k}(\mathbb R)&a_k&r&\mathcal A_k\\\hline
1&\mathbb R&1&\mathbb C&2&(2)&\mathbb Z_2\\
2&\mathbb C&2&\mathbb H&4&(2)&\mathbb Z_2\\
3&\mathbb H&4&2\mathbb H&4&(1\ 1)&0\\
4&2\mathbb H&4&M_2(\mathbb H)&8&\binom11&\mathbb Z\\
5&M_2(\mathbb H)&8&M_4(\mathbb C)&8&(1)&0\\
6&M_4(\mathbb C)&8&M_8(\mathbb R)&8&(1)&0\\
7&M_8(\mathbb R)&8&2M_8(\mathbb R)&8&(1\ 1)&0\\
8&2M_8(\mathbb R)&8&M_{16}(\mathbb R)&16&\binom11&\mathbb Z
\end{array}
$$

各行の計算は次のとおりです。

- **$k=1,2$**：制限先の既約加群は1つで、制限すると2個になります。$\mathfrak M_{k-1}=\mathbb Z$を像$2\mathbb Z$で割って$\mathbb Z_2$です。
- **$k=3,7$**：$\operatorname{Cl}_{0,k}(\mathbb R)$の2つの既約加群は、どちらも制限すると制限先の既約加群1個になります。像は$\mathbb Z$全体で、$\mathcal A_k=0$です。
- **$k=4,8$**：制限先は$W,\overline W$の2つを持つ直和型で、$\mathfrak M_{k-1}=\mathbb Z^2$です。[[lem-restrict]]により$\operatorname{Cl}_{0,k}(\mathbb R)$の既約加群は$W\oplus\overline W$に制限され、像は$(1,1)$の倍数です。$\mathbb Z^2/\mathbb Z(1,1)\cong\mathbb Z$で、$[W]$の類が生成元です（$[\overline W]=-[W]$）。
- **$k=5,6$**：次元が等しく、制限すると既約加群1個です。$\mathcal A_k=0$です。

得られた$\mathbb Z_2,\mathbb Z_2,0,\mathbb Z,0,0,0,\mathbb Z$は、以前の記事で主張として紹介した$\widetilde{KO}(S^k)$の表と一致します。$\mathcal A_k$の生成元は、$k=1$で$\mathbb R$（メビウスの帯）、$k=2$で$\mathbb C$（ホップ束）、$k=4$で$\mathbb H$（四元数のホップ束）、$k=8$で$\mathbb O$（八元数のホップ束）の類で、アティヤ＝ボット＝シャピロの定理ではK群の生成元に対応します。[[7shi-kth1]]

以前の記事では、$\mathbb Z$が現れる$k=4,8$は$\operatorname{Cl}_{0,k-1}(\mathbb R)$が直和型になる位置であることを見て、$k=1,2$の$\mathbb Z_2$の現れ方には既約加群の次元が関わると述べるに留めました。表から、その仕組みが読み取れます。$\mathbb Z$は、直和型の2つの既約加群のうち、制限で打ち消されない片方の向きとして現れます。$\mathbb Z_2$は、1つ多い生成元を加えたときに既約加群の次元が2倍になり、制限で2個ずつしか消えないことから現れます。$k=5,6,7$で次元が変わらないときは、どの既約加群も延長でき、何も残りません。[[7shi-kth1]]

&&&thm アティヤ＝ボット＝シャピロの定理
$[W]\mapsto[E_W]-\dim W$は、同型$\mathcal A_k\cong\widetilde{KO}(S^k)$と$\mathcal A_k^{\mathbb C}\cong\tilde K(S^k)$を与えます。
&&&

本記事で示したのは、この写像が矛盾なく定まることと、延長できる加群が$0$になることです。写像が同型であること、とくに$k=4,8$で四元数・八元数のホップ束の類が$0$でないこと（自明束を足しても自明にならないこと）は、本記事では扱いません。$k=1$の同型は、以前の記事で示した$\widetilde{KO}(S^1)\cong\mathbb Z_2$です。[[7shi-kth1]]

## 複素K群の表

複素でも同じ計算ができます。$\operatorname{Cl}_j(\mathbb C)$は、$j$が偶数なら$M_{2^{j/2}}(\mathbb C)$、奇数なら$2M_{2^{(j-1)/2}}(\mathbb C)$でした。既約加群の複素次元は、$j=0,1,2,3,\dots$で$1,1,2,2,\dots$です。[[7shi-clif1]]

- **$k$が偶数**：$\operatorname{Cl}_{k-1}(\mathbb C)$は直和型で、既約加群は$W,\overline W$の2つです。$\operatorname{Cl}_k(\mathbb C)$の既約加群は次元が2倍で、[[lem-restrict]]と同じ議論で$W\oplus\overline W$に制限されます。$\mathcal A_k^{\mathbb C}=\mathbb Z^2/\mathbb Z(1,1)\cong\mathbb Z$です。
- **$k$が奇数**：$\operatorname{Cl}_{k-1}(\mathbb C)$は単純型で、$\operatorname{Cl}_k(\mathbb C)$の2つの既約加群はどちらも同じ次元で、制限すると既約加群1個です。$\mathcal A_k^{\mathbb C}=0$です。

複素では、擬スカラーに虚数単位を掛けて2乗を$1$にそろえられるので、[[lem-restrict]]の議論は$k-1$が奇数なら常に使えます。結果は$\mathbb Z,0$が交互に並ぶ周期2の表で、$\tilde K(S^k)$の表と一致します。

$k=2$では、$\mathcal A_2^{\mathbb C}\cong\mathbb Z$の生成元$[W]$は$[H]-1$に移り、[[thm-s2]]の$φ$で$1$に移ります。$φ$との合成が$\mathbb Z$から$\mathbb Z$への同型なので、$\mathcal A_2^{\mathbb C}\to\tilde K(S^2)$は単射です。$S^2$の場合、この写像が同型であることは、単連結性を認めて得た$\tilde K(S^2)\cong\mathbb Z$と同じ内容です。

# ボット周期性

## 代数の側の周期

$\mathcal A_k$の計算に使ったのは、$\operatorname{Cl}_{0,k-1}(\mathbb R)$と$\operatorname{Cl}_{0,k}(\mathbb R)$の型（単純型か直和型か）と、既約加群の次元の比$a_k/a_{k-1}$だけでした。分類表を導いた記事では、次の8周期性を示しました。[[7shi-clif1]]

$$
\operatorname{Cl}_{0,k+8}(\mathbb R)\cong\operatorname{Cl}_{0,k}(\mathbb R)\otimes M_{16}(\mathbb R)
$$

$\otimes M_{16}(\mathbb R)$は型を変えず、既約加群の次元を16倍にします（$a_{k+8}=16a_k$）。したがって$k$と$k+8$で、型も次元の比も同じです。制限先が直和型の場合に2つの既約加群がどう配分されるかは、[[lem-restrict]]により$k-1\equiv3,7\pmod8$で常に同数ずつと決まります。こうして制限の行列も周期8で繰り返し、次が成り立ちます。

$$
\mathcal A_{k+8}\cong\mathcal A_k,\qquad\mathcal A_{k+2}^{\mathbb C}\cong\mathcal A_k^{\mathbb C}
$$

複素の周期2は、$\operatorname{Cl}_{n+2}(\mathbb C)\cong\operatorname{Cl}_n(\mathbb C)\otimes_{\mathbb C}M_2(\mathbb C)$から同様に従います。[[7shi-clif1]]

## 位相の側の周期

アティヤ＝ボット＝シャピロの定理を認めると、代数の側の周期がそのまま球面のK群に移ります。

&&&thm ボット周期性
$k\ge1$について、次の同型があります。

$$
\widetilde{KO}(S^{k+8})\cong\widetilde{KO}(S^k),\qquad\tilde K(S^{k+2})\cong\tilde K(S^k)
$$
&&&

ボット周期性定理は、アティヤ＝ボット＝シャピロの定理より先に位相の側で証明されたもので、証明は本記事では扱いません。本記事の計算は、分類表の周期8と球面のK群の周期8が、どちらもクリフォード代数の加群の制限の表から読めることを示しています。分類表を導いた記事で「対応することが知られています」と述べた対応は、この形で具体的になります。[[7shi-clif1]]

&&&rem アダムスの定理
以前の記事では、$S^{n-1}$の各点で線形独立な接ベクトル場の最大本数がラドン＝フルヴィッツ数$\rho(n)-1$であること、単位元つきの連続な積を持てる球面が$S^0,S^1,S^3,S^7$に限られることを、アダムスの定理として主張に留めました。どちらも、クリフォード加群から作った構成が最良であることを述べるもので、証明にはK理論の演算（アダムス演算）が用いられます。接ベクトル場の議論では実射影空間のKO群などを使い、ホップ不変量1の問題（積を持つ球面の問題）については、アティヤとアダムスによるK理論を用いた証明が知られています。これらの議論には、本記事では立ち入りません。[[7shi-kth2]][[7shi-lie7]]
&&&

# まとめ

球面の複素K群の元を、貼り合わせ関数の安定なホモトピー類として表しました。実の場合にも、ホモトピックな貼り合わせ関数は同じK群の元を与えます。平面の回転によって、K群の和は貼り合わせ関数の積にあたります。$S^2$では、行列式の回転数が$\tilde K(S^2)$から$\mathbb Z$への全射な準同型を与え、直線束の直和の類はこの値で決まります。$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群$W$から貼り合わせ関数$g_W$を作ると、1つ多い生成元まで延長できる加群は自明な束を与えます。加群の群を延長できる加群の制限で割った群を分類表から計算すると、実では$\mathbb Z_2,\mathbb Z_2,0,\mathbb Z,0,0,0,\mathbb Z$、複素では$0,\mathbb Z$の繰り返しとなり、球面のK群の表と一致します。

&&& 回転による変形
$R(t)$をブロックの回転として、次のホモトピーがあります。

$$
\operatorname{diag}(g,h)\simeq\operatorname{diag}(gh,I)
$$

K群では$([E_g]-n)+([E_h]-n)=[E_{gh}]-n$です。$S^2$では$H\oplus H\cong H^2\oplus\underline{\mathbb C}$、$[H^m]-1=m([H]-1)$です。
&&&

&&& 加群が与える貼り合わせ関数
$\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群$W$の生成元の作用$J_l$から、次の貼り合わせ関数を作ります。

$$
g_W(\boldsymbol x)=x_0I+x_1J_1+\dots+x_{k-1}J_{k-1}
$$

$J_k$まで延長できれば、$\cos t\,g_W+\sin t\,J_k$により$g_W$は定数にホモトピックです。
&&&

&&& 制限で割った群
$\mathfrak M_j$を$\operatorname{Cl}_{0,j}(\mathbb R)$の加群の群、$r$を制限として、次のように定めます。

$$
\mathcal A_k=\mathfrak M_{k-1}/r(\mathfrak M_k)
$$

$k=1,\dots,8$で$\mathbb Z_2,\mathbb Z_2,0,\mathbb Z,0,0,0,\mathbb Z$です。
&&&

K群の側と加群の側の対応を並べます。

| 球面のK群 | クリフォード加群 |
|:---|:---|
| 貼り合わせ関数$g_W$による類$[E_W]-\dim W$ | $\operatorname{Cl}_{0,k-1}(\mathbb R)$の加群$W$ |
| 束の直和 | 加群の直和 |
| $E_W$が自明（したがって$[E_W]-\dim W=0$） | $\operatorname{Cl}_{0,k}(\mathbb R)$に延長できる加群 |
| $\mathbb Z$（$k=4,8$）と符号の逆の元 | 直和型の2つの既約加群$W,\overline W$ |
| $\mathbb Z_2$（$k=1,2$） | 延長側の既約加群の次元が制限先の2倍になること |
| ボット周期性（周期8・2） | $\operatorname{Cl}_{0,k+8}\cong\operatorname{Cl}_{0,k}\otimes M_{16}(\mathbb R)$、$\operatorname{Cl}_{n+2}(\mathbb C)\cong\operatorname{Cl}_n(\mathbb C)\otimes M_2(\mathbb C)$ |
