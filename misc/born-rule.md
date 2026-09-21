2準位系の状態ベクトル$\Psi=(α,β)^T$から密度行列$ρ=\Psi\Psi^\dagger$とブロッホベクトル$\boldsymbol r$を導き、方向$\boldsymbol n$への測定の期待値がトレース$\mathrm{tr}(ρ(\boldsymbol n\cdot\boldsymbol σ))=\boldsymbol r\cdot\boldsymbol n$で与えられること、そこから確率$p_\pm=\dfrac{1\pm\boldsymbol r\cdot\boldsymbol n}2$が求まることは既に示されています[[7shi-bloch]]。ここでは、同じ確率を状態ベクトル$\Psi$の内積から直接求める方法（ボルンの規則）を示します。

シリーズ: [ホップファイブレーション](https://mathlog.info/series/sKmD4S7IQSBnq4CvOVlU)

# 状態ベクトルの成分と確率

$z$軸方向の測定を例に、確率が状態ベクトルの成分からどう求まるかを確認します。ブロッホベクトルの$z$成分と規格化条件は

$$
z=αα^*-ββ^*,\quad αα^*+ββ^*=1
$$

なので、これらを$p_+=\dfrac{1+z}2,\ p_-=\dfrac{1-z}2$に代入すれば

$$
\begin{aligned}
p_+&=\frac{1+z}2=\frac{(αα^*+ββ^*)+(αα^*-ββ^*)}2=αα^*=|α|^2 \\
p_-&=\frac{1-z}2=\frac{(αα^*+ββ^*)-(αα^*-ββ^*)}2=ββ^*=|β|^2
\end{aligned}
$$

となります。

&&&fml $z$軸方向の測定確率
規格化された状態ベクトル$\Psi=(α,β)^T$を$z$軸方向に測定したときの確率は、成分の絶対値の2乗として得られます。
$$
p_+=|α|^2,\quad p_-=|β|^2
$$
&&&

この$z$軸方向の測定では、確率$|α|^2,|β|^2$は密度行列$ρ=\Psi\Psi^\dagger$の対角成分そのものです。測定方向を変えた場合は、同じ行列表現の対角成分がそのまま確率になるわけではありません。

# 基底状態との内積

成分$α,β$は、$z$軸の基底状態$\Psi_\uparrow=(1,0)^T,\ \Psi_\downarrow=(0,1)^T$との内積として取り出せます。

$$
\Psi_\uparrow^\dagger\Psi=\begin{pmatrix}1 & 0\end{pmatrix}\begin{pmatrix}α \\ β\end{pmatrix}=α, \quad
\Psi_\downarrow^\dagger\Psi=\begin{pmatrix}0 & 1\end{pmatrix}\begin{pmatrix}α \\ β\end{pmatrix}=β
$$

これを使えば、確率は内積の絶対値の2乗という形に書き直せます。

$$
p_+=|\Psi_\uparrow^\dagger\Psi|^2
,\quad
p_-=|\Psi_\downarrow^\dagger\Psi|^2
$$

&&&rem 半角公式との対応
極角$θ$を使って$z=\cosθ$とおき、半角の公式から$p_+=\cos^2\dfracθ2$を導く方法が既に示されています[[7shi-bloch]]。そこでの$θ$はブロッホベクトル$\boldsymbol r$と測定方向$\boldsymbol n$のなす角であり、$\boldsymbol n$を$z$軸に選べば$\boldsymbol r$の極角（$0\leθ\leπ$）と一致します。

このとき
$$
|α|^2=\frac{1+\cosθ}2=\cos^2\fracθ2
$$
なので$|\Psi_\uparrow^\dagger\Psi|=|α|=\cos\dfracθ2$となり、$p_+=|\Psi_\uparrow^\dagger\Psi|^2=\cos^2\dfracθ2$が同じ式に一致します。つまり半角の公式による2乗は、内積の絶対値の2乗に対応します。一般には$α=\cos\dfracθ2\,e^{ia}$であって$α$自身が実数とは限りませんが、絶対値を取るためこの対応は位相によらず成り立ちます。
&&&

# ボルンの規則

内積で書き直した形は、測定方向を$z$軸に限定しません。ブロッホベクトルの内積$\boldsymbol r\cdot\boldsymbol n$が期待値を与えたのと対をなして、状態ベクトルの内積$\Psi_i^\dagger\Psi$の絶対値の2乗は確率を与えます。

任意の単位ベクトル$\boldsymbol n$を選ぶと、ブロッホベクトルが$\pm\boldsymbol n$を指す純粋状態が定まります。それらを表す規格化された状態ベクトル$\Psi_\pm$は、各々の全体位相を除いて定まり、互いに直交します。確率は基底状態と状態のどちらの全体位相にも依存しません。

$$
|(e^{iχ}\Psi_\pm)^\dagger(e^{iφ}\Psi)|^2=|\Psi_\pm^\dagger\Psi|^2
$$

この$\Psi_\pm$を基底状態として、状態$\Psi$から測定結果$i$を得る確率は$z$軸の場合と同じ形の式で与えられます。これを**ボルンの規則**と呼びます。なお$\Psi_\pm$は測定方向$\boldsymbol n$で決まる基底であり、状態$ρ$の固有分解に現れる固有ベクトル[[7shi-bloch]]とは一般に異なります。

&&&fml ボルンの規則
規格化された2準位系の純粋状態$\Psi$と、測定方向$\boldsymbol n$に対応する正規直交基底$\Psi_+,\Psi_-$を考えます。この基底は$\boldsymbol n\cdot\boldsymbolσ$の固有値$+1,-1$に対応する固有ベクトルからなります。結果$i$が得られる確率$p_i$は、対応する基底状態$\Psi_i$との内積の絶対値の2乗として得られます。
$$
p_i=|\Psi_i^\dagger\Psi|^2\quad(i=+,-)
$$
&&&

$z$軸で確かめた形が任意の方向で成り立つことは、密度行列とトレースを経由して確認できます。[[7shi-bloch]]

&&&prf 任意の測定方向での確認
$|\boldsymbol n|=1$とし、ブロッホベクトルが$\pm\boldsymbol n$である規格化された状態ベクトルを$\Psi_\pm$とする。対応する密度行列は
$$
P_\pm=\Psi_\pm\Psi_\pm^\dagger=\frac12(\mathrm I\pm\boldsymbol n\cdot\boldsymbolσ)
$$
である。$ρ=\Psi\Psi^\dagger$より
$$
\begin{aligned}
|\Psi_\pm^\dagger\Psi|^2
&=(\Psi_\pm^\dagger\Psi)(\Psi^\dagger\Psi_\pm) \\
&=\Psi_\pm^\dagger ρ\,\Psi_\pm \\
&=\mathrm{tr}(ρP_\pm) \\
&=\frac{\mathrm{tr}\,ρ\pm\mathrm{tr}(ρ(\boldsymbol n\cdot\boldsymbolσ))}2 \\
&=\frac{1\pm\boldsymbol r\cdot\boldsymbol n}2
\end{aligned}
$$
となり、期待値から求めた確率$p_\pm$と一致する。
&&&

これは既に採用した測定確率の式と内積表示が同値であることの確認であり、測定についての前提なしにボルンの規則を導出しているわけではありません。

&&&ex $x$軸方向の測定
$x$軸の基底状態は、$σ_x$の固有ベクトル
$$
\Psi_\pm=\frac1{\sqrt2}\begin{pmatrix}1\\\pm1\end{pmatrix}
$$
です。規格化条件$αα^*+ββ^*=1$と$x=α^*β+β^*α$より
$$
\begin{aligned}
p_\pm
&=|\Psi_\pm^\dagger\Psi|^2 \\
&=\frac12|α\pmβ|^2 \\
&=\frac12(α\pmβ)(α^*\pmβ^*) \\
&=\frac12(αα^*+ββ^*\pmα^*β\pmβ^*α) \\
&=\frac{1\pm x}2
\end{aligned}
$$
となり、$z$軸で見た$p_\pm=\dfrac{1\pm z}2$と同じ形が$x$軸でも成り立つことが確認できます。
&&&

# 位相の役割

$z$軸方向の測定確率には、成分の絶対値$|α|,|β|$だけが現れます。一方、$x$軸方向の測定確率には$α^*β+β^*α$が現れ、両成分の相対位相も関わります。

$α=|α|e^{ia},\ β=|β|e^{ib}$と書けば$α^*β=|α||β|e^{i(b-a)}$です。このように、成分の絶対値に加えて相対位相$b-a$を保持することで、$z$軸以外の方向の測定確率も求められます。全体位相を掛けても$b-a$は変わらないため、どの方向の測定確率も変わりません。

# まとめ

測定の確率は、密度行列側では期待値$\boldsymbol r\cdot\boldsymbol n$からスケール調整によって、状態ベクトル側では基底状態との内積の絶対値の2乗（ボルンの規則）として求まります。両者は同じ確率を異なる経路で与えます。状態ベクトルでは、成分の絶対値と相対位相が測定確率を決め、全体位相は確率に影響しません。
