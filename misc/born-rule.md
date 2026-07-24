2準位系の状態ベクトル$\Psi=(α,β)$から密度行列$ρ=\Psi\Psi^\dagger$とブロッホベクトル$\boldsymbol r$を導き、方向$\boldsymbol n$への測定の期待値がトレース$\mathrm{tr}(ρ(\boldsymbol n\cdot\boldsymbol σ))=\boldsymbol r\cdot\boldsymbol n$で与えられること、そこから確率$p_\pm=\dfrac{1\pm\boldsymbol r\cdot\boldsymbol n}2$が求まることは既に示されています[[7shi-bloch]]。ここでは、同じ確率を状態ベクトル$\Psi$の内積から直接求める方法（ボルンの規則）を示します。

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

となります。確率$|α|^2,|β|^2$は密度行列$ρ=\Psi\Psi^\dagger$の対角成分そのものです。

# 基底状態との内積

成分$α,β$は、$z$軸の基底状態$\Psi_\uparrow=(1,0),\ \Psi_\downarrow=(0,1)$との内積として取り出せます。

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

&&&rem
極角$θ$を使って$z=\cosθ$とおき、半角の公式から$p_+=\cos^2\dfracθ2$を導く方法が既に示されています[[7shi-bloch]]。$α=\cos\dfracθ2$とおけば$|\Psi_\uparrow^\dagger\Psi|=|α|=\cos\dfracθ2$となるので、$p_+=|\Psi_\uparrow^\dagger\Psi|^2=\cos^2\dfracθ2$となり、同じ式に一致します。つまり半角の公式による2乗は、内積の絶対値の2乗に対応します。
&&&

# ボルンの規則

内積で書き直した形は、測定方向を$z$軸に限定しません。ブロッホベクトルの内積$\boldsymbol r\cdot\boldsymbol n$が期待値を与えたのと対をなして、状態ベクトルの内積$\Psi_i^\dagger\Psi$は確率を与えます。

任意の方向$\boldsymbol n$を選べば、ブロッホベクトルが$\pm\boldsymbol n$を指す純粋状態の状態ベクトル$\Psi_i$が対応する基底状態として定まり、状態$\Psi$から測定結果$i$を得る確率は同じ形の式で与えられます。これを**ボルンの規則**と呼びます。

&&&def ボルンの規則
状態$\Psi$を測定したとき、結果$i$が得られる確率$p_i$は、対応する基底状態$\Psi_i$との内積の絶対値の2乗として得られる。
$$
p_i=|\Psi_i^\dagger\Psi|^2
$$
&&&

&&&ex $x$軸方向の測定
$x$軸の基底状態は、$σ_x$の固有ベクトル
$$
\Psi_\pm=\frac1{\sqrt2}\begin{pmatrix}1\\\pm1\end{pmatrix}
$$
です。$x=α^*β+β^*α$より
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

ここまで見たように、確率$p_+,p_-$は測定方向$\boldsymbol n$を指定して初めて取り出せる値です。$\Psi$の成分$α,β$が実数ではなく位相を持つ複素数であるのは、方向に関わらず定まる$\Psi$自身の中に、どの方向を指定しても確率を引き出せるだけの情報をあらかじめ位相として持たせるためです。

# まとめ

測定の確率は、密度行列側では期待値$\boldsymbol r\cdot\boldsymbol n$からスケール調整によって、状態ベクトル側では基底状態との内積の絶対値の2乗（ボルンの規則）として求まります。両者は同じ確率を異なる経路で与えており、状態ベクトルの位相は、この内積を通じてどの方向の確率も引き出せるようにするための情報です。
