# 検討メモ：時空代数とディラック方程式

パウリ方程式からディラック方程式への流れを、実クリフォード代数（Hestenes形式）で扱うシリーズ`dirac`の検討メモです。公開向けの構想は[README.md](README.md)にまとめ、本メモには構成の判断と、記事化の前に確かめる事項を記録します。

[PLAN.md](../PLAN.md)の「時空代数とディラック方程式」に対応します。[幾何代数による電磁気学](../em/README.md)（全7回、完結）が委ねた話題（スピノル、二重被覆、$\operatorname{SL}(2,\mathbb C)$、パウリ方程式・ディラック方程式）を引き受けます（[em/MEMO.md](../em/MEMO.md)の「確定した方針」）。

# 確定した方針

- **Hestenes形式で書く**。スピノルを列ベクトルではなく偶部分代数の元$\psi$として扱い、ディラック方程式を$\operatorname{Cl}_{1,3}(\mathbb R)$の中で書く。行列形式（ガンマ行列と$\mathbb C^4$）は照合のために翻訳するが、主役にはしない
- **符号数は$\operatorname{Cl}_{1,3}(\mathbb R)$**（$\gamma_0^2=+1$、$\gamma_k^2=-1$）。emシリーズと揃える（理由は[em/MEMO.md](../em/MEMO.md)の「符号数の選択」）。マヨラナ形式（$\operatorname{Cl}_{3,1}(\mathbb R)\cong M_4(\mathbb R)$の実スピノル）は`&&&rem`で触れるに留める
- **パウリから始める**。前半は$\operatorname{Cl}_{3,0}(\mathbb R)$でパウリスピノルとパウリ方程式を扱い、後半で$\operatorname{Cl}_{1,3}(\mathbb R)$に上げる。emシリーズと同じく、3次元の形式を時空代数の偶部分として回収する流れにする
- **emシリーズとの境界**。ローレンツ変換（ベクトル・2ベクトルに作用する回転子）、電磁場$F$、ポテンシャル$A$とゲージは[[7shi-em5]]・[[7shi-em6]]で済んでいるものとして後方参照する。本シリーズはスピノル（回転子の片側作用を受ける量）から先を扱う
- **qiとの境界**。[qi/PLAN.md](../qi/PLAN.md)（検討段階）は量子ビットを単位四元数として扱う案を持つ。本シリーズは空間・時空の関数としてのスピノル場と、それが従う方程式を扱い、有限次元の状態とゲート（量子情報）には立ち入らない。qiは方針未定なので、本シリーズからqiに依存しない
- **扱わないもの**。水素原子の厳密解、場の量子化（第二量子化）、曲がった時空のディラック作用素（[K理論と指数定理](../ktheory/README.md)のシリーズと[PLAN.md](../PLAN.md)の「ゲージ理論重力」に委ねる）、多粒子（多粒子時空代数、qiの検討に委ねる）

## 記号

- **擬スカラーは$\omega$**。[STYLE.md](../STYLE.md)の原則どおりで、emシリーズの$i$は引き継がない。量子力学ではシュレーディンガー方程式$i\hbar\partial_t\Psi$の虚数単位$i$が定着しており、[[7shi-clif5]]もこの形で導いている。複素数の虚数単位は$i$のまま使い、擬スカラー$\omega=\gamma_0\gamma_1\gamma_2\gamma_3=\sigma_1\sigma_2\sigma_3$と区別する。clif-analysisの$\omega$（[[7shi-cla1]]の$e_1e_2$、[[7shi-cla5]]の$e_0\cdots e_{n-1}$）とも一致する
- **スピノルは$\psi$、列ベクトル・左イデアルの元は$\Psi$**。Hestenes形式の慣例に従い、偶部分代数の元を$\psi$と書く。行列形式の列ベクトル（左イデアルの元）は$\Psi$とする。[[7shi-cla1]]の命題「偶部分と極小左イデアルの対応」は偶部分の$F$から左イデアルの$\psi=FP$を作っており、文字の役割が逆になるので、引くときに断る
- **行列にはハットを付ける**。代数の元$\sigma_k$・$\gamma_\mu$と、それを表現する行列$\hat\sigma_k$・$\hat\gamma_\mu$を区別する。翻訳の節で両者が並ぶため。ハットは行列表示と列スピノルに作用する演算子にだけ付け、単位ベクトルには付けない（$\boldsymbol n$、$\boldsymbol b$。単位ベクトルと断る）
- **作用素**。時空のディラック作用素は[[7shi-em4]]と同じく$D=\sum_\mu\gamma^\mu\partial_\mu$（相反基底$\gamma^0=\gamma_0$、$\gamma^k=-\gamma_k$）、$x_0=ct$、$D^2=\partial_0^2-\Delta=\square$。空間の作用素は、[[7shi-em6]]の証明で使った$D_3=\sum_k\sigma_k\partial_k$と書き、01から$D$と区別しておく（emの01〜03は空間の作用素を$D$と書いたが、本シリーズは02で時空に上がるため）
- **3次元の生成元は$\sigma_k$**。01から$\operatorname{Cl}_{3,0}(\mathbb R)$の生成元を$\sigma_k$と書き、02で[[7shi-em4]]の相対ベクトル$\sigma_k=\gamma_k\gamma_0$として回収する。em/01〜03の$e_k$は使わない
- **定数は$\hbar$と$c$を残す**。非相対論極限で$c\to\infty$を取るため、自然単位系にはしない。emシリーズのSI単位系と揃う。電荷は一般に$q$と書き、電子では$q=-e$とする。ポテンシャルは[[7shi-em6]]の$A=\varphi\gamma_0+c\sum_kA_k\gamma_k$（単位は$\varphi$と同じ）をそのまま使うので、結合項には$1/c$が付く（04の検討事項）

# シリーズ構成

全5回で完結しました（全回レビュー済み）。01は$\operatorname{Cl}_{3,0}(\mathbb R)$の回、02〜03で時空代数に上げてディラック方程式に至り、04〜05で電磁場との結合と非相対論極限を扱ってパウリ方程式に戻ります。04と05は分量が多いので分けました。

| # | ファイル | slug | 内容 | 到達点 |
|---|---|---|---|---|
| 01 | `01-pauli.md` | `7shi-dirac1` | パウリスピノルとパウリ方程式 | 列ベクトル$\mathbb C^2$は偶部分代数$\operatorname{Cl}_{3,0}^0(\mathbb R)\cong\mathbb H$の元$\psi=\sqrt\rho\,R$に移る。虚数単位$i$は右からの$\omega\sigma_3=\sigma_1\sigma_2$、スピンの向きは$\psi\sigma_3\tilde\psi$（ホップ写像）。一様磁場中のスピンの歳差 |
| 02 | `02-spinor.md` | `7shi-dirac2` | 時空スピノルとSL(2,C) | $\operatorname{Cl}_{1,3}^0(\mathbb R)\cong\operatorname{Cl}_{3,0}(\mathbb R)\cong M_2(\mathbb C)$（双四元数）。回転子の群$\operatorname{Spin}^+(1,3)\cong\operatorname{SL}(2,\mathbb C)$が$\operatorname{SO}^+(1,3)$を二重に覆う。ディラックスピノルは$\psi=\sqrt\rho\,e^{\omega\beta/2}R$ |
| 03 | `03-dirac.md` | `7shi-dirac3` | ディラック方程式 | $\square$の平方根として$D$を取り、行列形式$i\hbar\hat\gamma^\mu\partial_\mu\Psi=mc\Psi$を$\hbar D\psi\,\omega\sigma_3=mc\,\psi\gamma_0$に翻訳する。平面波解と正負のエネルギー |
| 04 | `04-coupling.md` | `7shi-dirac4` | 電磁場との結合とゲージ | 結合項$-qA\psi$、ゲージ変換は右からの$e^{\omega\sigma_3\chi}$（スピン軸まわりの回転）。流れ$J=\psi\gamma_0\tilde\psi$の保存 |
| 05 | `05-limit.md` | `7shi-dirac5` | 非相対論極限と$g=2$ | 大きい成分と小さい成分への分解からパウリ方程式を回収し、磁気モーメントの係数$g=2$を得る |

# 各回の検討事項

## 01 パウリスピノルとパウリ方程式

- **対応の作り方**。列ベクトル$(a_0+ia_3,\,-a_2+ia_1)^T$を$\psi=a_0+a_k\,\omega\sigma_k$に対応させる（Doran–Lasenbyの規約）。行列$\hat\sigma_k$の作用は$\sigma_k\psi\sigma_3$、$i$倍は$\psi\,\omega\sigma_3$になる。$\operatorname{Cl}_{3,0}^0(\mathbb R)\cong\mathbb H$なので、パウリスピノルは四元数そのもの
- **左イデアルとの関係**。[[7shi-lie3]]・[[7shi-ideal]]はスピノルを射影$P=(1+\sigma_3)/2$による左イデアル$\operatorname{Cl}_{3,0}(\mathbb R)P$の元として取った。$\psi\mapsto\Psi=\psi P$は偶部分代数から左イデアルへの全単射で、右からの$\omega\sigma_3$は$P$の上で$\omega$（中心の元、$M_2(\mathbb C)$の虚数単位）に一致する。これは[[7shi-cla1]]の命題「偶部分と極小左イデアルの対応」（$\operatorname{Cl}_{2,0}$、$P=(1+e_1)/2$）の3次元版で、同じ証明の流れで示せる。偶部分代数による表示は、左イデアルによる表示の言い直しとして導入する
- **虚数単位の意味**。$\psi\mapsto\psi e^{\omega\sigma_3\alpha}$は$\psi\sigma_3\tilde\psi$を変えない。大域位相は「スピン軸まわりの回転」で、[[7shi-h]]のホップファイバーにあたる。$i$が右から掛かる2ベクトルであることを、本シリーズを通じた最初の論点にする。[[7shi-em3]]で複素指数関数の$j$の役割を擬スカラーが担ったのに対し、ここで$i$の役割を担うのは2ベクトル$\omega\sigma_3$であり、この違いも述べる
- **観測量**。$\rho=\psi\tilde\psi$（確率密度）、$\psi\sigma_3\tilde\psi=\rho\boldsymbol s$（スピンの向き）。$\boldsymbol s$は[[7shi-h]]・[[7shi-bloch]]のブロッホベクトルと同じもの（対応は下記「既存シリーズとの関係」）。期待値は行列のエルミート共役を反転$\tilde\psi$に替えて書ける
- **パウリ方程式**。$i\hbar\partial_t\Psi=\hat H\Psi$を翻訳する。磁場の項$-\frac{q\hbar}{2m}\hat{\boldsymbol\sigma}\cdot\boldsymbol B$は$\boldsymbol B\psi\sigma_3$の形になる。例として一様磁場中のスピンの歳差を、回転子$R(t)$の時間発展として解く（係数$g$はここでは現象論的に置き、05で回収する）。回転の向きは[[7shi-em5]]の$R=e^{-\omega\sigma_3\theta/2}$（$x_1$から$x_2$へ角$\theta$）に揃える
- **前提知識**。シュレーディンガー方程式と確率解釈は[[7shi-clif5]]・[[7shi-born]]程度を前提とし、スピンは本記事で導入する。物理の用語は代数の構造に結び付けて導入する（emシリーズと同じ方針）

## 02 時空スピノルとSL(2,C)

- **偶部分代数**。$\sigma_k=\gamma_k\gamma_0$により$\operatorname{Cl}_{1,3}^0(\mathbb R)\cong\operatorname{Cl}_{3,0}(\mathbb R)$（[[7shi-em4]]で既出）。さらに$\operatorname{Cl}_{3,0}(\mathbb R)\cong M_2(\mathbb C)$は[[7shi-bq]]の双四元数そのもので、ディラックスピノル（複素4成分＝実8成分）は双四元数1つにあたる
- **二重被覆**。回転子$R\tilde R=1$の群$\operatorname{Spin}^+(1,3)$が$\operatorname{SL}(2,\mathbb C)$に同型で、$x\mapsto Rx\tilde R$が$\operatorname{SO}^+(1,3)$への2対1の写像になる。[[7shi-em5]]は変換そのものまでを扱ったので、ここで群の構造を扱う。空間の回転に制限すると[[7shi-lie3]]・[[7shi-cover]]の$\operatorname{SU}(2)\to\operatorname{SO}(3)$に戻る
- **スピノルの分解**。$\psi\tilde\psi=\rho e^{\omega\beta}$（スカラー＋擬スカラー）から$\psi=\sqrt\rho\,e^{\omega\beta/2}R$。$\beta$（Yvon–Takabayashi角）の物理的な意味は定まっていないので、深入りせず`&&&rem`に留める
- **観測量**。$\psi\gamma_\mu\tilde\psi=\rho e_\mu$で、回転子が基底$\gamma_\mu$を回した枠$e_\mu$が現れる。$e_0$が流れ、$e_3$がスピンの向き。01の$\psi\sigma_3\tilde\psi$の時空版
- **ワイルスピノル**（検討中）。擬スカラー$\omega$が行列形式の$\hat\gamma_5$の役割を担うこと、カイラリティを右からの射影で分けることに触れるか。触れるなら`&&&rem`で

## 03 ディラック方程式

- **導入の動機**。$D^2=\square$（[[7shi-em4]]）から、クライン＝ゴルドン方程式$(\square+m^2c^2/\hbar^2)\phi=0$の平方根として1階の方程式を探す。[[7shi-cla1]]がラプラシアンの平方根として$D$を組んだのと同じ発想で、係数の2乗の符号で双曲型になる（[[7shi-em7]]）
- **翻訳**。行列形式$i\hbar\hat\gamma^\mu\partial_\mu\Psi=mc\Psi$を$\hbar D\psi\,\omega\sigma_3=mc\,\psi\gamma_0$に移す。右から掛かる$\gamma_0$と$\omega\sigma_3$が、左からの作用（ローレンツ変換$\psi\mapsto R\psi$）と可換であることを要点にする。符号と右からの因子の順序は検証コードで行列形式と照合する
- **共変性**。$\psi\mapsto R\psi$で方程式の形が保たれることを、$D$の変換と合わせて示す
- **平面波解**。$\psi=\psi_0e^{-\omega\sigma_3\,p\cdot x/\hbar}$、$p^2=m^2c^2$。静止系の解から回転子（ブースト）で運動する解を作る。負のエネルギーの解は右からの$\sigma_3$などで表せるかを確かめる。解釈（空孔理論・反粒子）には深入りしない
- **質量0の場合**。$m=0$では$D\psi=0$となり、真空のマクスウェル方程式$DF=0$（[[7shi-em4]]）と同じ作用素の核になる。違いは値のグレード（$F$は2ベクトル、$\psi$は偶部分全体）だけである。静的なら$D_3\psi=0$で、[[7shi-cla5]]の$n=3$のモノジェニック関数（楕円型）に戻る。クリフォード解析からemを経て本シリーズへ、同じ作用素の核を値の空間と符号数を替えて追う流れとして`&&&rem`に置く
- **ジッターベヴェーグング**（検討中）。正負のエネルギーの重ね合わせで流れが振動する現象。`&&&rem`に留めるか省く

## 04 電磁場との結合とゲージ

- **結合**。$\hbar D\psi\,\omega\sigma_3-\frac qcA\psi=mc\,\psi\gamma_0$。[[7shi-em6]]の$A$は$c\sum_\mu A^\mu\gamma_\mu$（$A^\mu=(\varphi/c,\boldsymbol A)$）にあたるので、最小結合$p_\mu\to p_\mu-qA_\mu$は$\frac qcA$の形になる。[[7shi-em5]]のローレンツ力$m\,dU/d\tau=\frac qcF\cdot U$と同じ$q/c$である。結合項の符号は行列形式との照合で確定する
- **ゲージ変換**。[[7shi-em6]]の$A\mapsto A+D\chi$に対して$\psi\mapsto\psi e^{\omega\sigma_3\alpha}$、$\alpha=-q\chi/\hbar c$（上の符号で手計算。$\omega\sigma_3=\gamma_2\gamma_1$が$\gamma_0$と可換なので質量項は変わらない）。右からの回転なので$e_0$と$e_3$は変わらず、$e_1,e_2$がスピン軸まわりに回る。[[7shi-em6]]のゲージの自由度が、スピノルでは枠の回転として見える。01の大域位相を局所化したもの
- **保存則**。$D\cdot J=0$（$J=\psi\gamma_0\tilde\psi$）を方程式から導く。[[7shi-em2]]の連続の式と対比する
- **ローレンツ力**。流れの運動方程式に$F$が現れることを、エーレンフェストの定理（[[7shi-clif5]]）の時空版として示せるか（検討中）

## 05 非相対論極限と$g=2$

- **分解**。$\psi$を$\gamma_0$と可換な部分と反可換な部分（大きい成分・小さい成分）に分け、$c\to\infty$で小さい成分を消去する
- **到達点**。パウリ方程式の磁場の項が係数$g=2$で現れる。01で現象論的に置いた係数を回収する
- **$g=2$の出所**。消去の後に残る運動量の項（$D_3$とベクトルポテンシャル$\boldsymbol A$を組んだ作用素）の2乗の2ベクトル部が、磁場の2ベクトル$D_3\wedge\boldsymbol A=\omega\,\nabla\times\boldsymbol A=\omega\boldsymbol B$（[[7shi-em6]]の$D\wedge A$の空間部分）を与える。幾何積$D_3\boldsymbol A=D_3\cdot\boldsymbol A+D_3\wedge\boldsymbol A$の外積の部分がスピンと磁場の結合になる、という読み方を軸にする（係数は検証コードで確認）
- **04との統合**。分量が少なければ04の最後の節にする

# 記事化で確定した事項

全5回で執筆し、全回をレビューして指摘を反映した（04と05は分けた）。検討中だった項目の扱いと、執筆・レビューで確定した事項は次のとおり。

- **回転子の規約**：作用はemシリーズに合わせて$x\mapsto Rx\tilde R$、四元数との対応は$\mathbf i=-\omega\sigma_1$、$\mathbf j=-\omega\sigma_2$、$\mathbf k=-\omega\sigma_3$（[[7shi-bq]]と同じ）。[[7shi-lie3]]は$x\mapsto R^{-1}xR$で、同じ回転子が互いに逆元になる。本シリーズの内部では整合が取れているので統一は保留し、01の「記号」に断りを入れた（経緯は[lie/NOTATION.md](../lie/NOTATION.md)）。レビューが規約差を指摘しても、統一は提案せず、この断りで足りるかで判断する
- **共変成分の記号**：04で、$A$の共変成分を$a_\mu=\gamma_\mu\cdot A/c$と書く。空間成分$A_k$との衝突を避けるため
- **適用範囲の明示**：スピンは1/2粒子に限り$\hbar/2$の係数を落とさない。分解と枠は$\psi\tilde\psi\ne0$のスピノルに対するもので、非零でも$\psi\tilde\psi=0$となる例（ワイルスピノルなど）は対象外。一様磁場の歳差は状態がスピン部分と軌道部分に分離している場合（01）、近似の前提は正エネルギー側の低エネルギー状態（05）
- **左右の作用**：左乗算（ローレンツ変換）と右乗算（複素構造$\omega\sigma_3$など）は可換で、右乗算は複素構造を与える。一般の右乗算は観測量も変える（$e^{\omega\sigma_3\alpha}$だけが$e_0,e_3$を保つ）。左右の可換性は値に対する代数作用の話で、場の変換では引数も変わる
- **$J$の保存**：03で$J$の正値性までを示して保存は先送りし、04で方程式から示す

- **行列の規約**：03は標準（Bjorken–Drell）のディラック表現（上付き$\hat\gamma^k$の右上ブロックが$\hat\sigma_k$）を使う。対応は$\psi=\phi+\eta\sigma_3\mapsto(|\phi\rangle,|\eta\rangle)^T$で、$\hat\gamma_\mu\leftrightarrow\gamma_\mu\psi\gamma_0$、$\hat\gamma_5\leftrightarrow\psi\sigma_3$。Doran–Lasenbyは下付きの符号が逆の表現を使うので、引用時に注意
- **ワイルスピノル**：03の`&&&rem`（右からの射影$(1\pm\sigma_3)/2$）。02では扱わない
- **負のエネルギーと電荷の反転**：右からの$\sigma_1$。03で負のエネルギーの平面波（$\beta=\pi$、流れは同じ）、04で電荷$q\to-q$として扱う
- **ジッターベヴェーグング、ローレンツ力（エーレンフェストの定理の時空版）**：扱わない
- **$g=2$とサイクロトロン運動**：01の最後で、$g=2$ならスピンの回転子の方程式がem/05の固有速度の回転子の方程式と同じ形になることを示し、05で回収する

# 検証コード

`check/`に置き、`uv run`で実行する（一覧は[check/README.md](check/README.md)）。emシリーズで一般化した`src/common/clifford.py`（$\operatorname{Cl}_{p,q}$のビットマスク実装と$D$）を使う。記事化の前に確かめた事項：

- 01の対応規約（$\hat\sigma_k\Psi\leftrightarrow\sigma_k\psi\sigma_3$、$i\Psi\leftrightarrow\psi\,\omega\sigma_3$）と、$\psi\mapsto\psi P$の全単射
- 03の翻訳：ディラック表現の$\hat\gamma_\mu$で行列形式と$\hbar D\psi\,\omega\sigma_3=mc\,\psi\gamma_0$が同値になる対応写像
- 03の平面波解と$p^2=m^2c^2$、負のエネルギーの解の形
- 04の結合項とゲージ変換の符号・係数（emシリーズの単位系で）
- 05の非相対論極限で$g=2$が出ること

# 既存シリーズとの関係

本シリーズは、[クリフォード解析](../clif-analysis/README.md)→[幾何代数による電磁気学](../em/README.md)と続いた「ディラック作用素$D$の核」の系列の3番目にあたる。clif-analysisは$D^2=\Delta$（楕円型）の核として正則関数を、emは$D^2=\square$（双曲型）で電磁場を扱った。本シリーズは同じ$D$を、値を偶部分代数（スピノル）に取り、右から掛かる構造（$\omega\sigma_3$と$\gamma_0$）と質量項を加えて使う。作用素と符号数はemからそのまま引き継ぎ、新しく加わるのは値の空間と右からの作用である。

## 幾何代数による電磁気学（em）

前提として最も多く引く。emが境界として委ねた話題（[em/MEMO.md](../em/MEMO.md)の「確定した方針」）を引き受ける。

- **[[7shi-em4]]**：時空代数$\operatorname{Cl}_{1,3}(\mathbb R)$、相反基底と$D$、$D^2=\square$、$\sigma_k=\gamma_k\gamma_0$による$\operatorname{Cl}_{1,3}^0(\mathbb R)\cong\operatorname{Cl}_{3,0}(\mathbb R)$、$\omega\sigma_3=\gamma_2\gamma_1$などの2ベクトルの対応。符号数の`&&&rem`（複素化の中で$\operatorname{Cl}_{1,3}$の生成元を虚数単位倍すると$\operatorname{Cl}_{3,1}$の生成元になる）と例「パウリ行列からの構成」（ワイル表現のガンマ行列）も02・03のマヨラナ形式・行列形式の断りで引く。定義は再掲せず後方参照する
- **[[7shi-em5]]**：回転子の定義（偶部分で$R\tilde R=1$）、$x\mapsto Rx\tilde R$、回転$e^{-\omega\sigma_3\theta/2}$とブースト$e^{\sigma_1\eta/2}$。emは変換そのものと、観測者の取り替え（基底$\gamma_\mu'=R\gamma_\mu\tilde R$、成分は$\tilde RFR$）までを扱った。本シリーズは$\psi\mapsto R\psi$（片側の作用）を加え、同じ$R$と$-R$がベクトルには同じ変換を与えること（二重被覆）を02で扱う。能動・受動の区別はemの書き方に合わせる
- **[[7shi-em6]]**：ポテンシャル$A=\varphi\gamma_0+c\sum_kA_k\gamma_k$、$F=D\wedge A$、ゲージ変換$A\mapsto A+D\chi$、ローレンス条件。04の結合とゲージ変換はこの$A$と$\chi$をそのまま使う
- **[[7shi-em2]]・[[7shi-em3]]**：連続の式（04の$D\cdot J=0$と対比）、複素指数関数の$j$を擬スカラーが担うこと（01で$i$を2ベクトルが担うことと対比）
- **[[7shi-em7]]**：双曲型の基本解。質量のある場合の基本解（クライン＝ゴルドン）は本シリーズでも扱わない
- **記号の差**：emの擬スカラー$i$は本シリーズでは$\omega$、emの複素数の虚数単位$j$は$i$になる。emの式を引くときは$F=\boldsymbol E+\omega c\boldsymbol B$のように書き換え、01で1回断る。$D$・$\gamma_\mu$・$\gamma^\mu$・$\sigma_k$・反転$\tilde{\ }$・単位系（SI、$x_0=ct$）・電荷$q$・ローレンツ力の$\frac qcF\cdot U$はemと同じ
- **emの未解決課題**：リエナール＝ヴィーヘルトのポテンシャルと物質中の場（[PLAN.md](../PLAN.md)）は本シリーズでも扱わない

## クリフォード解析（clif-analysis）

- **[[7shi-cla1]]**：ラプラシアンの平方根としての$D$（03の導入の動機の前例）、命題「偶部分と極小左イデアルの対応」（01の$\psi\mapsto\psi P$の前例）、擬スカラー$\omega$の記号。文字$\psi$の役割が逆である点は「記号」の項のとおり
- **[[7shi-cla4]]**：$e_0$を掛けて偶部分に移す操作。[[7shi-em4]]の$\gamma_0$を掛ける操作と同じ構造で、03の質量項$\psi\gamma_0$でも$\gamma_0$が基準方向として現れる（こちらは右から掛かる）。この対比に触れるかは記事化の際に判断する
- **[[7shi-cla5]]**：一般次元のモノジェニック関数。03の質量0・静的な場合の`&&&rem`で$n=3$の場合として引く
- **[[7shi-cla6]]・[[7shi-em7]]**：楕円型と双曲型の性質の仕分け。本シリーズでは新たな仕分けはせず、必要なら後方参照する
- **未解決課題**：球面モノジェニックス、奇数次元のフューター＝ソーの定理、八元数解析（[PLAN.md](../PLAN.md)）は本シリーズでも扱わない。球面モノジェニックスは水素原子の角度部分（スピノル球面調和関数）と関係するが、水素原子を扱わないので立ち入らない

## スピノル・回転子を扱う既存記事（lie・hopf・qua・clif）

これらは記号の規約がシリーズごとに異なる。本シリーズはemの規約に合わせ、引くときに対応を断る。

- **[[7shi-lie3]]・[[7shi-ideal]]・[[7shi-ladder]]**：射影$P=(1+\sigma_3)/2$による左イデアルとしてのスピノル、昇降演算子と純粋スピノル。01で偶部分代数による表示と結ぶ。[[7shi-lie3]]はスピノルの列ベクトルを$\omega$と書き、回転を$R^{-1}xR$（emとは左右が逆）と書いている。四元数との対応も$i\leftrightarrow\sigma_1\sigma_2$、$k\leftrightarrow\sigma_2\sigma_3$で、[[7shi-bq]]とは異なる
- **[[7shi-bq]]**：$\operatorname{Cl}_{3,0}(\mathbb R)\cong M_2(\mathbb C)$と双四元数、四元数との対応$\mathbf i\cong-i\sigma_x$、$\mathbf j\cong-i\sigma_y$、$\mathbf k\cong-i\sigma_z$（本シリーズの記号では$\mathbf k\cong-\omega\sigma_3$）。01の「パウリスピノルは四元数」と02の「ディラックスピノルは双四元数」はこの対応で述べる
- **[[7shi-h]]・[[7shi-bloch]]・[[7shi-s]]・[[7shi-qgate]]**：ホップ写像$\omega\mapsto\omega\mathbf k\omega^*$（hopf/01の$\omega$は単位四元数）、ブロッホベクトルと密度行列、大域位相とファイバー、ゲートの位相×回転子分解。[[7shi-bq]]の対応では$\sigma_3=\omega\mathbf k$、反転が四元数の共役にあたるので、$\psi\sigma_3\tilde\psi=\omega\,(\psi\mathbf k\psi^*)$となり、01のスピンの向きはホップ写像の像の双対（2ベクトルとベクトルの読み替え）になる（符号は検証コードで確認）。スピンの向きをベクトル、四元数の虚部を2ベクトルとして区別できることは[PLAN.md](../PLAN.md)の「グレード」の論点そのものなので、01で明示する
- **[[7shi-cover]]・[[7shi-lie4]]・[[7shi-4drot]]**：二重被覆の位相的な直観、$\operatorname{Spin}(4)$。02の$\operatorname{Spin}^+(1,3)\to\operatorname{SO}^+(1,3)$は、コンパクトな場合の二重被覆を非コンパクトな場合に広げたものとして引く。[lie/README.md](../lie/README.md)が時空代数に委ねた非コンパクト群の構造は、[[7shi-em5]]が変換までを扱い、群の構造を本シリーズの02が扱う
- **[[7shi-clif1]]**：$\operatorname{Cl}_{1,3}(\mathbb R)\cong M_2(\mathbb H)$、$\operatorname{Cl}_{3,1}(\mathbb R)\cong M_4(\mathbb R)$、偶部分代数の同型。02とマヨラナ形式の`&&&rem`で引く
- **[[7shi-clif5]]・[[7shi-born]]**：シュレーディンガー方程式の導出、ボルンの規則。01の前提とする量子力学

## 量子情報（qi、検討段階）

[qi/PLAN.md](../qi/PLAN.md)の案C・案Dは、量子ビットを$\operatorname{Cl}_{3,0}^0(\mathbb R)\cong\mathbb H$の元（回転子）として扱う。01の「パウリスピノルは四元数」と同じ対応であり、qiが記事化されれば01を後方参照できる関係にある。qiは方針未定なので、本シリーズは01でこの対応を自前で示し、qiに依存しない。有限次元の状態・ゲート・多粒子はqiに任せ、本シリーズは空間・時空の関数としてのスピノル場を扱う。

## 全体構想の他テーマ

- **ゲージ理論重力（GTG）**：時空代数上のゲージ理論で、ディラック方程式の共変微分が出発点の1つになる。本シリーズの04（ゲージ変換を右からの回転として読む）が前提として引かれうる
- **[K理論と指数定理](../ktheory/README.md)**：曲がった空間（球面）上のディラック作用素とスピン構造。本シリーズは平坦な時空に限る
