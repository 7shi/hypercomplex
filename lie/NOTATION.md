# 回転子・四元数の規約の状況

`vec-oct`、`lie`、`em`、`dirac`、`qua/cd` で、回転の作用と四元数との対応の書き方が異なる。統一に向けて議論中で、現状と論点、方針の候補を記録する。

## 経緯

`dirac/01` のレビュー結果を検討する中で、参照記事の `lie/03` と本文の回転子の作用の規約が異なることが議論になった（`lie/03` は $R^{-1}xR$、`dirac` は `em` に合わせた $Rx\tilde R$）。`lie/03` に沿った書き方は、鏡映の積とクリフォード代数の指数関数の読みから自然だが、`dirac` は `em`・`qua/cd` と一貫した規約で書かれており、`dirac` の中では整合が取れている。

統一には `em`・`dirac` 全体の書き換えが必要で大掛かりなため、いったん保留とし、`dirac` は現行の規約のまま進めた。`dirac/01` の「記号」には、`lie/03` との規約の差を断る一文を入れてある。

その後、規約の出どころである `vec-oct/02` を確認し、著者の流儀を `vec-oct/notation.md`（記事「回転の向きと積の順序の規約」、slug `7shi-vnot`）にまとめた。流儀の内容と導出はそちらを前提とし、本メモでは各記事の状況と統一の論点だけを扱う。

## 現状

- `vec-oct/02`：回転子を $r=mn=\exp(\frac\theta2p)$、$p=m\wedge n/|m\wedge n|$ とし、作用は $v\mapsto r^{-1}vr$。四元数との対応は $ae_3e_2+be_1e_3+ce_2e_1\cong ai+bj+ck$（負号あり、向きを保つ）。その結果 $p\cong-q$ となり、四元数の $rvr^*$ とは表面上挟み方が逆になると説明している。鏡映から導く体裁だが、重点は指数関数の対応にある。未公開。
- `lie/02`：以前は四元数の行列表現を $i,j,k\Leftrightarrow i\sigma_3,i\sigma_2,i\sigma_1$（`lie/03` で添字が逆順と示す形）で天下りに導入していた。行列式が絶対値の2乗になる要請から形を絞り込み、標準形 $I_H,J_H,K_H=-i\sigma_1,-i\sigma_2,-i\sigma_3$ を導く形に改訂した（パウリ行列は使わず、$\beta=-i$ の符号の選び方は先送り）。ローカルでは改訂済み、Mathlog未反映。
- `lie/03`：回転子を鏡映の積 $R=uv$ から導き、作用は $x\mapsto R^{-1}xR$。以前は四元数との対応が $i,j,k\Leftrightarrow\sigma_1\sigma_2,\sigma_3\sigma_1,\sigma_2\sigma_3$（添字逆順）で、スピノルの例に $-\theta$ が出ていた。`lie/02` に合わせて $i,j,k\Leftrightarrow-i\sigma_1,-i\sigma_2,-i\sigma_3=\sigma_3\sigma_2,\sigma_1\sigma_3,\sigma_2\sigma_1$ に改訂し、$q\Leftrightarrow R^{-1}$ で $qxq^{-1}$ と $R^{-1}xR$ が一致する形にした。`lie/02` で先送りした符号の理由（向きをそろえるため）をremで回収。ローカルでは改訂済み、Mathlog未反映。
- `qua/cd`（`bq` を含む）：$I_H,J_H,K_H=-i\sigma_1,-i\sigma_2,-i\sigma_3$。`lie/02` とは別の行列表現で、添字の順を保つ。
- `em`：回転子は偶部分の元で $R\tilde R=1$ を満たすもの、作用は $x\mapsto Rx\tilde R$。正の回転の回転子は $e^{-i\sigma_3\theta/2}$ のように指数に負号が付く（`em` では擬スカラーを $i=\gamma_0\gamma_1\gamma_2\gamma_3$ と書く）。未公開。
- `dirac`：`em` に合わせて $x\mapsto Rx\tilde R$、四元数との対応は $\mathbf i=-\omega\sigma_1$ など（`qua/cd` と同じ）。未公開。`dirac/01` の「記号」に、`lie/03` との規約の差の断りを入れた。

## 回転を扱う記事の一覧

回転の作用を式で扱う記事を、公開・未公開に分けて整理する。公開日は `articles.tsv` による。型の名前は次のとおり。

- 流儀型：`vec-oct/notation` の規約。クリフォード代数（またはパウリ行列）で $\tilde RxR$（$r^{-1}vr$、$R^\dagger VR$）と書き、指数は $\exp(+\theta B/2)$。
- 文献流：$Rx\tilde R$ と書き、指数は $\exp(-\theta B/2)$。
- 四元数型：四元数・八元数だけで $qvq^{-1}$・$rvr^*$ と書く。`vec-oct/notation` の四元数の節と一致する。
- 対応のみ：作用は扱わず、四元数とクリフォード代数・パウリ行列との対応だけを扱う。

### 公開済み

| 記事 | 公開日 | 型 | 作用・対応 | 流儀との関係 |
|---|---|---|---|---|
| `hopf/01-quaternion` | 2024/05/14 | 四元数型 | $qpq^*$ | 一致 |
| `hopf/02-spinor-tensor` | 2024/05/24 | 四元数型 | $q\mathbf kq^*$ | 一致 |
| `oct/02-7d-3rot` | 2024/06/26 | 四元数型 | 八元数の $rxr^*$ | 一致 |
| `qua/01-pauli-qua` | 2024/07/11 | 流儀型 | $R^\dagger VR$、$R=\exp(\frac\theta2N)$、$N=i\boldsymbol n\cdot\boldsymbol\sigma\cong-n$。四元数は $rvr^*$、$R^\dagger\cong r$ | 一致 |
| `vec-oct/geometric-product-exp` | 2024/07/24 | 対応のみ | $ab=\|a\|\|b\|e^{\theta e_1e_2}$、$i\cong e_1e_2$、四元数は $i\cong e_3e_2$（負号あり）。2次元の外積を擬スカラー $e_1e_2$ の係数（符号付き面積）として扱う（Mathlog未反映） | 一致（`vec-oct/notation` の外積の節の出典） |
| `qua/cd/matrix-to-pauli` | 2025/05/03 | 対応のみ | $I_H=-i\sigma_1$ など（負号あり） | 一致 |
| `qua/spherical-trig` | 2025/11/29 | 四元数型 | $qvq^{-1}$ | 一致 |
| `lie/02-su2-so3` | 2026/07/10 | 四元数型 | $qxq^{-1}$、等傾回転。行列表現を標準形 $-i\sigma_k$ に改訂（Mathlog未反映） | 一致（改訂前は行列表現が `lie/03` の添字逆順の起点） |
| `qua/04-4d-bsqua` | 2026/07/20 | 流儀型 | $r^{-1}vr$、$r=\exp(B/2)$。四元数は $r_Lqr_R$。$\rho_{rs}=\rho_s\circ\rho_r$ の断りあり | 一致 |
| `lie/03-spin` | 2026/07/24 | 流儀型 | $R^{-1}xR$、$R=\exp(\varphi\sigma_1\sigma_2)$。四元数は $qxq^{-1}$、対応を $-i\sigma_k$ に改訂（Mathlog未反映） | 一致（改訂前は対応が添字逆順） |

公開済みの記事は、`lie/02`・`lie/03` の改訂（Mathlog未反映）を含めて、すべて流儀と一致している。文献流の公開済み記事はない。

### 未公開

| 記事 | 型 | 作用・対応 | 流儀との関係 |
|---|---|---|---|
| `vec-oct/notation` | 流儀型 | 流儀そのものの説明。$qvq^{-1}$ から双四元数を介して $r^{-1}xr$ を導く。対応は $i\cong e_3e_2\cong-i\sigma_1$ | 流儀の定義 |
| `vec-oct/02-rotation` | 流儀型 | $r^{-1}vr$、$r=mn$。四元数は $rvr^*$、$ae_3e_2+be_1e_3+ce_2e_1\cong ai+bj+ck$、4Dは $r_Lqr_R$ | 一致。流儀を最も忠実に記事にした基準 |
| `lie/04-spin4` | 流儀型 | $gvg^{-1}$（$g=R^{-1}$ と断り、四元数に向きをそろえる）。四元数は $qxq^{-1}$、$pxq^{-1}$ | 一致（`lie/03` の $R^{-1}vR$ を受ける） |
| `lie/05-bch-adjoint`・`lie/07-s7`・`lie/08-g2`・`lie/double-cover` | 四元数型 | $\rho_q(x)=qxq^{-1}$ | 一致 |
| `hopf/06-clifford-gates` | 流儀型 | $R^\dagger VR$（`qua/01` を引用）、$R_{\boldsymbol n}(\theta)=R^\dagger$ が $r$ の行列表現と断る | 一致 |
| `qua/05-dual-qua`・`qua/06-slerp` | 四元数型 | 二重四元数による剛体変換、$qvq^{-1}$ の補間 | 一致 |
| `qua/cd/bcmp-qua` | 対応のみ | $-i\sigma_k$（負号あり） | 一致 |
| `qua/history` | 文献流 | $R\boldsymbol v\tilde R$（歴史記事での一般的な紹介）。四元数は $q\boldsymbol vq^*$ | 作用が異なる |
| `em/05-lorentz` | 文献流 | $Rx\tilde R$、$R=e^{-i\sigma_3\theta/2}$（回転・ブースト） | 作用が異なる |
| `dirac/01`〜`04` | 文献流 | $Rx\tilde R$、$\psi\gamma_\mu\tilde\psi=\rho e_\mu$。四元数との対応は $\mathbf i=-\omega\sigma_1$。`dirac/03` は座標変換で $x'=\tilde RxR$ も使う | 作用が異なる、対応は一致 |
| `clif/pga-cga` | 文献流 | $X\mapsto MX\widetilde M$、$M=rq$（合成の順）。PGA・CGAの剛体変換 | 作用が異なる |
| `ktheory/06`〜`08` | 文献流 | $Ue_0\tilde U=\hat\theta$ など（枠の読み取り）。`dirac` と同じ向き | 作用が異なる |

未公開の記事のうち、流儀と異なるのは文献流の `qua/history`・`em/05`・`dirac/01`〜`04`・`clif/pga-cga`・`ktheory/06`〜`08` で、いずれも作用の向きだけが異なる。四元数との対応は、未公開の記事ではすべて負号付きで流儀と一致している。

## 著者の流儀

流儀の内容は `vec-oct/notation.md` にまとめた。要点だけを挙げる。

- 原則：複素数の積にならい反時計回りを正とし、反時計回りを与える元を負号なしで指数の肩に乗せる。
- 四元数：軸の単位は左から反時計回りに作用し、回転は $qvq^{-1}$（$q=e^{\theta u/2}$）。
- クリフォード代数：双四元数を介して $e_k\cong\sigma_k\cong hi,hj,hk$ と対応させると、2ベクトルは $i\cong e_3e_2$ と負号付きで四元数に移る。$q\cong r^{-1}$（$r=e^{\theta p/2}$）となり、回転は $r^{-1}xr$。右作用はこの結果として現れる。
- 外積：$m\wedge n$ は「$m$ から $n$ へ」の符号付き面積（2次元は $e_1e_2$ の係数、3次元以上は大きさと単位2ベクトル $p$ に分ける）。$r=mn$、$m\,r=n$、鏡映 $-n(-mvm)n=r^{-1}vr$、合成 $r_1r_2$ がすべて左から右へそろう。
- 文献流 $Rx\tilde R$ との違いは指数の肩の負号（$R=\tilde r$）に帰着する。負号は文献流では指数に、流儀では四元数との対応に置かれる。
- 擬スカラー $\omega$ は可換なので、グレードの違いは向きに影響しない。

`vec-oct/geometric-product-exp` は以前、2次元でも $|a||b|\sin\theta=|a\wedge b|$ と大きさで対応させ、$\sin\theta\ge0$ の制限を付けていたが、符号付き面積の読みに改めた（Mathlogへの反映は後述）。

## `lie/03` の向きの食い違い（改訂前）

`lie/02`・`lie/03` の改訂で解消した。経緯として残す。


- 負号を避けて添字を逆順に取る（`lie/02`・`lie/03`）と、四元数の $(i,j,k)$ に対応する軸は $(\sigma_3,\sigma_2,\sigma_1)$ で、$\sigma$ の右手系から見ると向きが反転する。
- その結果、$r=\exp(\theta\,\mathbf i/2)$ による四元数の共役作用 $rqr^{-1}$ が反時計回りでも、$\sigma$ の右手系では時計回りになる。`lie/03` の $R^{-1}xR$ は流儀どおりだが、対応だけが負号を避けた形のため、四元数の側と向きが食い違う。
- `lie/03` の外積の節で $-\theta$ が出るのも、この向きの反転の現れである。
- 以前は「$R^{-1}xR$ は $R=uv$ を『$u$ のあとに $v$』の積と読む約束の帰結で、右作用・左作用の問題ではない」と整理していた。流儀を踏まえると、鏡映の読みは後付けで、右作用は四元数との対応から結果として現れるものである。

## 規約の比較

| | クリフォード代数の作用 | 四元数との対応 | 指数 | 流儀との関係 |
|---|---|---|---|---|
| `vec-oct/02` | $r^{-1}vr$ | $k\cong e_2e_1$（負号あり） | $\exp(+\theta p/2)$ | 一致 |
| `qua/cd`・`bq` | ― | $K_H=-i\sigma_3$（負号あり） | ― | 対応は一致 |
| `lie/03` | $R^{-1}xR$ | $i\Leftrightarrow\sigma_3\sigma_2=-i\sigma_1$（改訂後） | $\exp(+\theta B/2)$ | 一致（改訂前は対応が添字逆順） |
| `em`・`dirac` | $Rx\tilde R$ | $\mathbf i=-\omega\sigma_1$（負号あり） | $\exp(-\theta B/2)$ | 対応は一致、作用が異なる |

四元数との対応は、`vec-oct`・`qua/cd`・`dirac` がすでに負号付きで一致している。`lie/03` の改訂後、流儀に統一する場合のずれは `em`・`dirac` の作用（$Rx\tilde R$）だけである。改訂前の `lie/03` のremで符号の違いが説明しにくかったのは、作用は流儀どおりで対応だけが負号を避けた形になっていたためである。

## 論点

- 流儀に統一する場合、`em`・`dirac` の作用を $\tilde RxR$ に改める。`em`・`dirac` は未公開なので書き換えの余地はあるが、次の波及がある。
    - $\psi=\sqrt\rho\,R$、$\psi\gamma_\mu\tilde\psi=\rho\,e_\mu$、$e_\mu=R\gamma_\mu\tilde R$ の形が崩れ、$\psi=\sqrt\rho\,\tilde R$ とするか、観測量を $\tilde\psi\gamma_\mu\psi$ 側に書き直す必要がある。
    - ローレンツ変換が「スピノルに左から $R$」、ゲージ変換が「右からの回転子」という左右の分担（`dirac/02`〜`04`）が入れ替わる。
    - `em/05` の回転子・ブーストの式と、`dirac` 全体の符号の向きが変わる。
    - $Rx\tilde R$ は、ヘステネスやドーラン＝ラゼンビーなど標準的な文献の書き方と一致する。流儀に合わせると文献との差を断る必要がある。
- 四元数との対応は `qua/cd`・`bq`（`dirac/01` も引用）の $-i\sigma_k$ をそのまま使える。
- `lie/02`・`lie/03` は標準形（$-i\sigma_k$）に改訂し、未公開の `lie/06`（成分内＝$k$、成分間＝$j$、積＝$i$、$aI-i(b\sigma_1+c\sigma_2+d\sigma_3)$）と `lie/04`（$x=u+vj$ の記号）も追従させた。`lie/04` の $gvg^{-1}$（$g=R^{-1}$）、`lie/05`・`lie/07`・`lie/08`・`lie/double-cover` の $qxq^{-1}$、`lie/10`・`ladder-spinor` のパウリ行列は四元数との対応に依存しないため変更不要。lieシリーズはこれで流儀に統一された。

## Mathlogへの反映

ローカルで改訂した公開済み記事のうち、Mathlogに未反映のものをまとめる。反映は `mathlog_fix.md` を書いて `bash src/mathlog_fix.sh` で行う（手順は `SLUG.md` の「更新手順」）。

- 参照（`[[slug]]`）を変えた記事は `--refs` で実行し、Mathlogの参考文献パネルも直したうえで、パネルをコピーして `refs/{ID}.toml` を取り込み直す。本文のみの修正は `--no-refs` でよい。オプションは実行単位で1つなので、両方がある場合は `mathlog_fix.md` を分けて2回実行するか、全体を `--refs` で実行する。
- 処理後に `mathlog_fix.md` を削除し、`make all` で `refs.toml` などを再生成する。
- `lie/02` の参照の変更により、現在 `reftools check` が `lie/02` について「未定義：`7shi-nonion`・`7shi-qp`」「未使用：`7shi-qcm`」を報告している。Mathlog側の参考文献パネルを直して取り込めば解消する。

### 未反映の記事

| 記事 | Mathlog | 内容 | 参照の変更 |
|---|---|---|---|
| `vec-oct/geometric-product-exp` | `YZmxak6ObeP6rLQnyU2V` | 2次元の外積を符号付き面積として扱い、角度の制限を3次元に限定 | なし（`--no-refs`） |
| `lie/02-su2-so3` | `Utdur1fLLzrWVHOJHifj` | 行列表現を行列式の要請から導く標準形に変更。$\mathfrak{su}(2)$ の分解、成分内・成分間の回転（$k$ と $j$）、等傾回転のremを追従。改訂履歴を追加 | `7shi-qp`・`7shi-nonion` を追加、`7shi-qcm` を削除（`--refs`） |
| `lie/03-spin` | `DRbXTeeL31pDcZyG6ml7` | 四元数との対応を $-i\sigma_k$ に変更。スピノルの例を $q=\exp(k\theta/2)$ にして全角 $+\theta$ に、$\mathfrak{su}(2)$ の一般形、クリフォード代数の対応・偶部分代数の乗積規則・$\operatorname{Spin}(3)$ の基底を追従。符号の選び方のremを追加し、共役作用との符号の違いのremを対応の説明に置き換え。改訂履歴を追加 | なし（`--no-refs`） |

### `mathlog_fix.md` の下書き

```markdown
## vec-oct/geometric-product-exp.md — https://mathlog.info/articles/YZmxak6ObeP6rLQnyU2V
- 本文：「2次元の外積」の節と「外積における角度の範囲」のremを差し替え（参照の変更なし）

## lie/02-su2-so3.md — https://mathlog.info/articles/Utdur1fLLzrWVHOJHifj
- 本文：改訂履歴、前提、「四元数の行列表現」の節、su(2)の分解、「SO(2)の複素ユニタリ化」の後半、「等傾回転の相殺」の分解とremを差し替え
- 7shi-qcm: 削除
- 7shi-qp: 追加
- 7shi-nonion: 追加

## lie/03-spin.md — https://mathlog.info/articles/DRbXTeeL31pDcZyG6ml7
- 本文：改訂履歴、概要の前提（su(2)の一般形、対角行列の作用）、スピノルの冒頭と例、「四元数との対応」の節、クリフォード代数の冒頭の対応、回転子のrem、偶部分代数の乗積規則とSpin(3)、グレード1への共役作用のremを差し替え（参照の変更なし）
```
