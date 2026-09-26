# 検討メモ：幾何代数による電磁気学

ディラック作用素を軸に、電磁気学と特殊相対論を幾何代数で組み直すシリーズ`em`の検討メモです。公開向けの構想は[README.md](README.md)にまとめ、本メモには構成の判断と、記事化の前に確かめる事項を記録します。

[PLAN.md](../PLAN.md)の「幾何代数ベースの電磁気学・相対性理論」に対応します。[クリフォード解析](../clif-analysis/README.md)（全6回、完結）の物理側への展開にあたり、$\operatorname{Cl}_{n,0}(\mathbb R)$のディラック作用素$D=\sum e_a\partial_a$（$D^2=\Delta$）で得た基本定理・核・積分公式を、そのまま電磁場に使います。

擬スカラーの記号は2026-09-26に$I$から$i$に改めました。複素数の虚数単位は$j$、電流は$I$と書きます（clif-analysisの擬スカラーは$\omega$）。以下の記録は旧表記のままです。

[四元数が脇役になった歴史的経緯](../qua/history.md)は、マクスウェルが四元数を併記し、ギブスとヘヴィサイドが積を内積とベクトル積に分離した経緯を扱っています。本シリーズはその経緯への回答として、分離された発散と回転を1つの作用素に戻し、マクスウェル方程式の4本の式を1本の式のグレード成分として回収します。

# 確定した方針

- **舞台は$\operatorname{Cl}_{3,0}(\mathbb R)$から$\operatorname{Cl}_{1,3}(\mathbb R)$へ進む**。前半（01〜03）は$\operatorname{Cl}_{3,0}(\mathbb R)$のパラベクトル形式で、clif-analysisと地続きに書く。後半（04〜07）は時空代数$\operatorname{Cl}_{1,3}(\mathbb R)$（$\gamma_0^2=+1$、$\gamma_k^2=-1$）に上げ、前半の形式をその偶部分として回収する
- **符号数は$\operatorname{Cl}_{1,3}$**（$\operatorname{Cl}_{3,1}$ではない）。理由は下記「符号数の選択」
- **単位系はSI、時間座標は$x_0=ct$**。「見慣れた4本の式がグレード成分として出る」ことを示すには、読者の知る式と照合できる必要がある。電磁場は$F=\boldsymbol E+Ic\boldsymbol B$、$\partial_0=\partial/\partial(ct)$。クーロン定数$1/4\pi\varepsilon_0$の$4\pi$は$|S^2|$として読む
- **記号は$D$を使う**。clif-analysisの判断（$\nabla$は3次元のベクトル解析に固有）を引き継ぐ。$\nabla\cdot$・$\nabla\times$はベクトル解析の式と照合する場面に限る
- **物理の用語は代数の構造に結び付けて導入する**。電場・磁場・電荷・電流は既知の量として受け取ってよいが、ゲージ・4元電流・エネルギー運動量テンソルなどは、グレード・擬スカラー・基本定理といった既出の構造から定義する
- **「時空代数とディラック方程式」との境界**。本シリーズは古典的な場と、ローレンツ変換（ベクトル・2ベクトルに作用する回転子）までを扱う。スピノル、二重被覆、$\operatorname{SL}(2,\mathbb C)$、パウリ方程式・ディラック方程式は同テーマに委ねる。[lie/README.md](../lie/README.md)が時空代数に委ねた非コンパクト群のうち、変換そのものは05で扱い、群の構造（二重被覆など）は引き続き委ねる

## 符号数の選択

当初は$\operatorname{Cl}_{3,1}(\mathbb R)$（計量$+,+,+,-$）を推していた。空間の生成元の2乗が$+1$のままでclif-analysisの$\operatorname{Cl}_{n,0}$と規約が揃い、[Spin(4)とSO(4)](../lie/04-spin4.md)・[四元数の左作用・右作用とスピノルの左イデアル](../qua/spinor-ideal.md)とも同じ代数になるためである。検討の結果、次の理由で$\operatorname{Cl}_{1,3}(\mathbb R)$に改めた。

- **パラベクトルの計量が$(+,-,-,-)$**。$\operatorname{Cl}_{3,0}$のパラベクトル$p=p_0+\boldsymbol p$では$p\bar p=p_0^2-|\boldsymbol p|^2$。前半の作用素も$\bar{\mathcal D}\mathcal D=\partial_0^2-\Delta$で、この符号である。$\operatorname{Cl}_{3,1}$では$D^2=\Delta-\partial_0^2$と全体の符号が逆になる
- **パウリ行列の素直な拡張**。$\gamma(p)=\begin{pmatrix}0&p\\\bar p&0\end{pmatrix}$とすると$\gamma(p)^2=p\bar p$で、$\sigma^\mu=(1,\sigma_k)$・$\bar\sigma^\mu=(1,-\sigma_k)$によるワイル表現のガンマ行列そのものになる。$\operatorname{Cl}_{3,1}$にするには下のブロックを$-\bar p$に変える符号が1つ要る
- **時空分割がclif-analysisの操作と一致する**。$\sigma_k=\gamma_k\gamma_0$で$\sigma_k^2=-\gamma_k^2\gamma_0^2=+1$、$\operatorname{Cl}_{1,3}^0(\mathbb R)\cong\operatorname{Cl}_{3,0}(\mathbb R)$。[[7shi-cla4]]の$h_l=e_0e_l$と同じく、2乗が$+1$の生成元を掛けて偶部分に移す。違いは空間の生成元の2乗の符号だけで、$h_l^2=-1$（楕円型）と$\sigma_k^2=+1$（双曲型）が分かれる。$\operatorname{Cl}_{3,0}$の擬スカラー$\sigma_1\sigma_2\sigma_3$は$\gamma_0\gamma_1\gamma_2\gamma_3$に一致するので、前半の$I$がそのまま使える
- **標準の文献と揃う**。通常のガンマ行列（$\gamma_0^2=+1$、$\gamma_k^2=-1$）とHestenesの時空代数は$\operatorname{Cl}_{1,3}$である

代償と対処：

- 時空のベクトルとしての$\gamma_k$は2乗が$-1$で、$\boldsymbol x^2=|\boldsymbol x|^2$の規約から外れる。空間のベクトルは相対ベクトル$\sigma_k$（2乗$+1$）として保たれるので、前半の規約はそのまま生きる。04で1回説明すれば足りる
- lie/04・qua/spinor-idealの$\operatorname{Cl}_{3,1}$は$\operatorname{SO}(4)$とマヨラナ型の実表現の文脈であり、電磁気学とは直接ぶつからない。04の`&&&rem`で、$\operatorname{Cl}_{3,1}\cong M_4(\mathbb R)$と$\operatorname{Cl}_{1,3}\cong M_2(\mathbb H)$の違い（複素化すればともに$M_4(\mathbb C)$）と、$\operatorname{Cl}_{3,1}$の生成元が$i\gamma^\mu$にあたることを断る
- ディラック方程式の側も$\operatorname{Cl}_{1,3}$のHestenes形式で書けば一貫する。マヨラナ形式を取るかは同テーマの検討に委ねる

# シリーズ構成

全7回です（全回の下書きを作成済み、レビューは未実施）。01〜03は[クリフォード解析](../clif-analysis/README.md)の結果を3次元で使う回、04〜06は符号数を替える回、07は[[7shi-cla6]]と対をなす総括です。

| # | ファイル（予定） | slug | 内容 | 到達点 |
|---|---|---|---|---|
| 01 | `01-static.md` | `7shi-em1` | 静電場と静磁場 | クーロン場とビオ＝サバールの法則は、3次元の核による積分公式のベクトル部と2ベクトル部。静的な方程式はその$D$微分 |
| 02 | `02-maxwell.md` | `7shi-em2` | マクスウェル方程式 | 1本の式$\mathcal DF=(\rho-\boldsymbol J/c)/\varepsilon_0$のグレード成分が4本の式。係数の2乗の符号で$\Delta$と$\partial_0^2-\Delta$が分かれる |
| 03 | `03-waves.md` | `7shi-em3` | 電磁波とエネルギーの流れ | 平面波は$F^2=0$、複素指数関数の虚数単位$j$の役割を擬スカラー$i$が担う。$\frac{\varepsilon_0}2FF^\dagger$がエネルギー密度とポインティングベクトルのパラベクトル |
| 04 | `04-spacetime.md` | `7shi-em4` | 時空代数とマクスウェル方程式 | $DF=J$。$\gamma_0$を掛けると02の式に戻る。楕円型と双曲型は空間の生成元の2乗の符号で分かれる |
| 05 | `05-lorentz.md` | `7shi-em5` | ローレンツ変換と回転子 | 2乗が$+1$の2ベクトルの指数関数がブースト。$F\mapsto RF\tilde R$で電場と磁場が混ざる。ローレンツ力 |
| 06 | `06-potential.md` | `7shi-em6` | ポテンシャルとエネルギー運動量 | $F=D\wedge A$、$D\cdot A$がゲージの自由度、$D^2A=J$。エネルギー運動量の時空版 |
| 07 | `07-retarded.md` | `7shi-em7` | 双曲型の基本解と性質の仕分け | 遅延ポテンシャル。clif-analysisの性質のうち、符号数を替えて残るものと失われるもの |

03（電磁波）は04の前、$\operatorname{Cl}_{3,0}$の側に置く。$F^2=0$などは時空版のほうが自然だが、前半を$\operatorname{Cl}_{3,0}$で完結させたほうがclif-analysisの読者は入りやすい。04以降で時空版に言い直すときは、03への後方参照にする。

# 各回の検討事項

## 01 静電場と静磁場

- **導入の順序**。クーロンの法則とビオ＝サバールの法則（実験則）を出発点とし、両者が1つの式$F(\boldsymbol x)=\frac1{4\pi\varepsilon_0}\int\frac{\boldsymbol x-\boldsymbol y}{|\boldsymbol x-\boldsymbol y|^3}\bigl(\rho-\boldsymbol J/c\bigr)(\boldsymbol y)\,dV$にまとまることを示す。核は[[7shi-cla5]]の$n=3$の核、$4\pi=|S^2|$。$D$を掛けると基本解の性質（[[7shi-cla3]]・[[7shi-cla5]]の`&&&rem`）から$DF=(\rho-\boldsymbol J/c)/\varepsilon_0$となり、そのグレード成分が静的な4本の式になる
- **磁場は2ベクトル**。$(\boldsymbol x-\boldsymbol y)\wedge\boldsymbol J$から$Ic\boldsymbol B$が出る。[[7shi-hist]]の「ベクトル積は擬ベクトル」の扱いを、$\boldsymbol B$がもともと2ベクトルであることで解消する
- **確認済み**。積分の$(\boldsymbol x-\boldsymbol y)\cdot\boldsymbol J$（スカラー部）は、$\boldsymbol x-\boldsymbol y=|\boldsymbol x-\boldsymbol y|^3\nabla_{\boldsymbol y}|\boldsymbol x-\boldsymbol y|^{-1}$と部分積分により$\nabla\cdot\boldsymbol J=0$から消える（源は有界な台を仮定し、特異点は小球を除いて処理）。記事では`&&&prop`とし、ポテンシャルの節でこれが$\nabla\cdot\boldsymbol A=0$と同じことだと述べた。$\nabla\cdot\boldsymbol J\ne0$ならスカラー部$\Phi$が残り、アンペールの法則が$\nabla\Phi$だけずれることを`&&&rem`にした（[check/01-static.py](check/01-static.py)）
- **積分形**。ガウスの法則とアンペールの法則の積分形は、[[7shi-cla2]]の領域の基本定理と曲面版のグレード成分
- **ポテンシャル（決定）**。01の最後の節に入れた。$P=\varphi-c\boldsymbol A$で$F=-DP$、スカラー部は$c\,\nabla\cdot\boldsymbol A$。ゲージの自由度は「本記事では扱いません」とし、06で時空版として扱う
- **核の表記（決定）**。[[7shi-cla5]]の$E=\boldsymbol x/|\boldsymbol x|^n$は電場$\boldsymbol E$と紛らわしいので、01では核を分数のまま書き、係数$1/4\pi\varepsilon_0$を積分の外に置いた

## 02 マクスウェル方程式

- $\mathcal D=\partial_0+D$（$D=\sum e_k\partial_k$、$\operatorname{Cl}_{3,0}$）。[[7shi-cla4]]の$\mathcal D=\partial_0+\sum h_l\partial_l$と同じ形で、係数の2乗が$-1$から$+1$に替わる。$\bar{\mathcal D}\mathcal D$は$\Delta$ではなく$\partial_0^2-\Delta$
- $\mathcal D(\boldsymbol E+Ic\boldsymbol B)=(\rho-\boldsymbol J/c)/\varepsilon_0$の成分：スカラー部が$\nabla\cdot\boldsymbol E=\rho/\varepsilon_0$、ベクトル部が$\frac1c\partial_t\boldsymbol E-c\nabla\times\boldsymbol B=-\mu_0c\boldsymbol J$、2ベクトル部が$\nabla\times\boldsymbol E+\partial_t\boldsymbol B=0$、擬スカラー部が$\nabla\cdot\boldsymbol B=0$（[check/memo-signature.py](check/memo-signature.py)で確認済み）
- $\bar{\mathcal D}$を掛けてスカラー部を取ると連続の式、全体では波動方程式
- [[7shi-hist]]への回答をここで明示する

## 03 電磁波とエネルギーの流れ

- 平面波$F=(1+\hat{\boldsymbol k})\boldsymbol E_0\exp(\cdots)$の形。$F^2=0$の成分が$\boldsymbol E\perp\boldsymbol B$と$|\boldsymbol E|=c|\boldsymbol B|$
- 複素指数関数の虚数単位$j$を擬スカラー$i$に置き換えられ、円偏光が平面内の回転になる（[[7shi-cla3]]の「虚数単位は擬スカラー」）。$i$は$\operatorname{Cl}_{3,0}$の中心にあるので可換性の問題は起きない
- $\frac{\varepsilon_0}2FF^\dagger=\frac{\varepsilon_0}2(|\boldsymbol E|^2+c^2|\boldsymbol B|^2)+\varepsilon_0c\,\boldsymbol E\times\boldsymbol B$（確認済み）。マクスウェル方程式からの保存則の導出

## 04 時空代数とマクスウェル方程式

- $\operatorname{Cl}_{1,3}(\mathbb R)$、$D=\sum\gamma^\mu\partial_\mu$（$\gamma^0=\gamma_0$、$\gamma^k=-\gamma_k$）。$\gamma_0D=\partial_0+\sum\sigma_k\partial_k$で、02の$\mathcal D$に一致する（確認済み）
- $F$は時空の2ベクトル（$\binom42=6$成分）で、$\sigma_k=\gamma_k\gamma_0$の側が$\boldsymbol E$、$I\sigma_k$の側が$c\boldsymbol B$。02の$F$と同じ元
- $DF=J$を$D\cdot F=J$と$D\wedge F=0$の2本に分ける。02の単位に合わせると$DF=\mu_0cJ$、$J=c\rho\gamma_0+\sum J_k\gamma_k$（確認済み）。係数の置き方（$\mu_0c$を$J$に含めるか）は記事で決める
- 「$\gamma_0$を掛ける」が[[7shi-cla4]]の「$e_0$を掛ける」と同じ操作で、空間の生成元の2乗の符号だけが違うこと。この対比をシリーズの軸として本文で明示する
- 符号数の`&&&rem`（上記「符号数の選択」の代償と対処）
- **決定**。係数は$DF=\mu_0cJ$（$J=c\rho\gamma_0+\sum J_k\gamma_k$）とした。符号数の`&&&rem`を書く際、$\operatorname{Cl}_{3,1}$でも$(\gamma_k\gamma_0)^2=+1$で偶部分は$\operatorname{Cl}_{3,0}$と同型である点に注意。違いは$\gamma_0^2=-1$による$D^2$の全体の符号と、$e_0$と同じ「2乗が$+1$の生成元を掛ける」操作が保てるかどうか

## 05 ローレンツ変換と回転子

- 空間の2ベクトル（2乗$-1$）の指数関数が回転、$\sigma_k=\gamma_k\gamma_0$（2乗$+1$）の指数関数がcosh・sinhによるブースト。2乗の符号の軸がここでも効く
- $F\mapsto RF\tilde R$で電場と磁場が混ざる。不変量$F^2=(|\boldsymbol E|^2-c^2|\boldsymbol B|^2)+2cI\,\boldsymbol E\cdot\boldsymbol B$（スカラー部と擬スカラー部）
- ローレンツ力$\dot p=qF\cdot v$
- $R$と$-R$が同じ変換を与えること（二重被覆）は「本記事では扱いません」とする
- 下書きでは、平面波を進行方向にブーストすると振幅と振動数が$e^{-\varphi}$倍（ドップラー因子）になる`&&&ex`を加えた。観測者の見る場は$\tilde RFR$（$\gamma_0'=R\gamma_0\tilde R$）とし、ローレンツ力は$m\,dv/d\tau=\frac qcF\cdot v$（$F\cdot v=\frac12(Fv-vF)$）
- 関連候補：[分解型複素数の冪等元と直和分解](../clif/02-split-complex.md)（2乗が$+1$の単位）。既存記事にcosh・sinhによる回転を扱ったものはない（`qua/05-dual-qua.md`に語が出るのみ）

## 06 ポテンシャルとゲージ

- $F=D\wedge A$。$D\cdot A$がゲージの自由度で、$D\cdot A=0$（ローレンスゲージ）なら$F=DA$、$D^2A=J$
- 静的な場合は[[7shi-cla5]]の$D|\boldsymbol x|^{2-n}=(2-n)\boldsymbol x/|\boldsymbol x|^n$（$n=3$）に戻る
- エネルギー運動量は$T(a)=-\frac{\varepsilon_0}2FaF$で確認済み。$T(\gamma_0)\gamma_0=\frac{\varepsilon_0}2FF^\dagger$（$F^\dagger=\gamma_0\tilde F\gamma_0$）、$\sum\partial_\mu T(\gamma^\mu)=-\frac1cF\cdot J$。03で先送りした運動量の保存をここで回収する。記事の後半がエネルギー運動量になったので、タイトルを「ポテンシャルとエネルギー運動量」に改めた

## 07 双曲型の基本解と性質の仕分け

- $\partial_0^2-\Delta$の遅延基本解$\delta(t-r/c)/4\pi r$と遅延ポテンシャル。$D$の基本解は$D$を掛けて得る（$D^2=\square$）
- [[7shi-cla6]]の仕分けの双曲型版。平均値の性質・最大値原理は失われ、光円錐上に台を持つ基本解（ホイヘンスの原理）に替わる。台が円錐上に集まるかどうかは空間の次元に依存する（空間2次元では集まらない）ので、次元の軸も再び現れる
- **決定**。キルヒホッフの公式は主張に留め（多項式の初期値で直接確かめられることだけ述べる）、リエナール＝ヴィーヘルトのポテンシャルは「本記事では扱いません」とした。仕分けは「形を保つもの（成分の方程式・基本定理・保存則）」「替わるもの（積分公式・一意性の根拠）」「失われるもの（平均値の性質・最大値原理・リウヴィル・一致の定理）」「次元で変わるもの（ホイヘンスの原理）」の4区分

# 検証コード

`check/`に置き、`uv run`で実行する（[check/README.md](check/README.md)）。clif-analysisの`check/04`〜`06`にあった$\operatorname{Cl}_{n,0}$のビットマスク実装（`MV`/`Alg`）を、負の2乗を持つ生成元に対応させて$\operatorname{Cl}_{p,q}$に一般化し、`src/common/clifford.py`に移した。$D$は逆基底$e^a=e_a^{-1}$で組む。clif-analysisの`check/04`〜`06`は移行後の実装で移行前と同じ出力になることを確認した。構成案の手計算は[check/memo-signature.py](check/memo-signature.py)で、各回の数式は`check/01-static.py`〜`check/07-retarded.py`で確認した。

# 関連記事

- [クリフォード解析](../clif-analysis/README.md) — 基本定理（[[7shi-cla2]]）、積分公式と基本解（[[7shi-cla3]]）、パラベクトル変数（[[7shi-cla4]]）、核とポテンシャル・$|S^{n-1}|$（[[7shi-cla5]]）、性質の仕分け（[[7shi-cla6]]）
- [四元数が脇役になった歴史的経緯](../qua/history.md) — 内積とベクトル積の分離、擬ベクトル。本シリーズはその回答にあたる
- [ベクトルの幾何積と指数関数](../vec-oct/geometric-product-exp.md) — 幾何積の分解
- [パウリ行列と四元数：双四元数による橋渡し](../qua/01-pauli-qua.md)・[クリフォード代数とSpin(3)](../lie/03-spin.md) — $\operatorname{Cl}_{3,0}$とパウリ行列
- [Spin(4)とSO(4)](../lie/04-spin4.md)・[四元数の左作用・右作用とスピノルの左イデアル](../qua/spinor-ideal.md) — $\operatorname{Cl}_{3,1}$を使う既存記事（符号数の`&&&rem`で触れる）
- [クリフォード代数の構造と表現](../clif/01-representation.md) — $\operatorname{Cl}_{1,3}\cong M_2(\mathbb H)$、$\operatorname{Cl}_{3,1}\cong M_4(\mathbb R)$、偶部分代数の同型
- [射影幾何代数と共形幾何代数](../clif/pga-cga.md) — ローレンツ計量$\operatorname{Cl}_{4,1}$を使う既存記事
