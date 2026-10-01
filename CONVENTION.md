# 回転子・四元数の規約の状況

`vec-oct`、`lie`、`em`、`dirac`、`qua` で、回転の作用と四元数との対応の書き方がそろっていなかった。規約を「四元数の回転子をクリフォード代数にそのまま移す」形に統一した。経緯と各記事の状況、Mathlogへの反映を記録する。

## 経緯

`dirac/01` のレビュー結果を検討する中で、参照記事の `lie/03` と本文の回転子の作用の規約が異なることが議論になった（`lie/03` は $R^{-1}xR$、`dirac` は `em` に合わせた $Rx\tilde R$）。

規約の出どころである `vec-oct/02` を確認し、著者の流儀を `vec-oct/convention.md`（記事「回転の向きと積の順序の規約」、slug `7shi-vconv`）にまとめた。当初は、反時計回りを与える2ベクトル $p$（回転面の向き）を負号なしで指数の肩に乗せて $r=\exp(\frac\theta2p)$ とし、作用を $r^{-1}xr$ とする形を流儀としていた。

その後、この $r$ は四元数の回転子 $q$ の像ではなく、$q^{-1}$ の像であることが問題になった。四元数の回転子 $q=e^{\frac\theta2k}$ を双四元数経由の対応 $k\cong e_2e_1$ でそのまま移すと $R=e^{\frac\theta2e_2e_1}=e^{-\frac\theta2e_1e_2}$ で、作用は四元数と同じ並びの $Rx\tilde R$ になる。指数の肩の見た目をそろえることにこだわった結果、回転子が逆元に取り替わっていた。そこで、四元数の回転子をそのまま移す $Rx\tilde R$ を規約とし、`vec-oct/convention` を書き直した。結果として、ヘステネスやドーラン＝ラゼンビーなど標準的な文献の書き方と一致する。

## 規約

内容と導出は `vec-oct/convention.md` を前提とする。要点だけを挙げる。

- 原則：複素数の積にならい反時計回りを正とし、反時計回りを与える元を負号なしで指数の肩に乗せる。
- 四元数：軸の単位は左から反時計回りに作用し、回転は $qvq^{-1}$（$q=e^{\theta u/2}$）。
- クリフォード代数：双四元数を介して $e_k\cong\sigma_k\cong hi,hj,hk$ と対応させると、2ベクトルは $i\cong e_3e_2\cong-i\sigma_1$ などと負号付きで四元数に移る。四元数の回転子はそのまま $R$ に移り、回転は $Rx\tilde R$。
- 位相の符号：回転面の向き $p=m\wedge n/|m\wedge n|$（$m$ から $n$ への符号付き面積）の解釈は変えない。左から反時計回りに作用するのは反転 $\tilde p=-p$ なので、$R=\exp(\frac\theta2\tilde p)=\exp(-\frac\theta2p)$ と、$p$ で書くと位相に負号が付く。
- 外積との整合：$R=nm$、$(nm)m=n$、鏡映 $-n(-mvm)n=Rv\tilde R$、合成 $R_2R_1$ は、写像の合成と同じく右から左へ読む。ウェッジ積だけが面の向きとして左から右へ読み、その食い違いが位相の負号になる。
- 2次元の複素数：$ve_1=x+y\,e_2e_1\leftrightarrow z$（$i\cong e_2e_1$）と対応させると、左からの $e^{i\theta}z$ に一致する。複素数の $i$ と四元数の $k$ が同じ $e_2e_1$ に移る。これは $k$ を対角（成分内）の単位に取る根拠にはなるが、$J_H$ が複素数 $i$ の実行列と同じ形になる理由には使えない（行列成分の $i$ は擬スカラー $\omega$ にあたる）。
- 擬スカラー $\omega$ は可換なので、グレードの違いは向きに影響しない。
- 負号を避けた $r=\tilde R=\exp(\frac\theta2p)$、$r^{-1}xr$ の書き方は、$r$ が $q^{-1}$ の像で、挟み方と合成の並びが四元数と逆になる（右作用）。`vec-oct/convention` の「別の書き方との関係」で扱う。

## 回転を扱う記事の一覧

公開日は `articles.tsv` による。型の名前は次のとおり。

- 規約型：クリフォード代数（またはパウリ行列）で $Rx\tilde R$（$RVR^\dagger$）と書き、$R$ は四元数の回転子の像。$p$ で書くと $R=\exp(-\frac\theta2p)$。文献の書き方と一致する。
- 四元数型：四元数・八元数だけで $qvq^{-1}$・$rvr^*$ と書く。
- 対応のみ：作用は扱わず、四元数とクリフォード代数・パウリ行列との対応だけを扱う。

### 公開済み

| 記事 | 公開日 | 型 | 作用・対応 | 状況 |
|---|---|---|---|---|
| `hopf/01-quaternion` | 2024/05/14 | 四元数型 | $qpq^*$ | 変更不要 |
| `hopf/02-spinor-tensor` | 2024/05/24 | 四元数型 | $q\mathbf kq^*$ | 変更不要 |
| `oct/02-7d-3rot` | 2024/06/26 | 四元数型 | 八元数の $rxr^*$ | 変更不要 |
| `qua/01-pauli-qua` | 2024/07/11 | 規約型 | $RVR^\dagger$、$R=\exp(-\frac{i\theta}2\boldsymbol n\cdot\boldsymbol\sigma)\cong r$。四元数は $rvr^*$。鏡映は $NM$ | 改訂済み（以前は $R^\dagger VR$、$R\cong r^*$）、Mathlog未反映 |
| `vec-oct/geometric-product-exp` | 2024/07/24 | 対応のみ | $ab=\|a\|\|b\|e^{\theta e_1e_2}$、2次元は $i\cong e_2e_1$（$ve_1$ 型の対応、$ba\cong a^*b$）、四元数は $i\cong e_3e_2$（負号あり） | 改訂済み（外積の符号付き面積、2次元の対応を $i\cong e_1e_2$ から変更）、Mathlog未反映 |
| `qua/cd/matrix-to-pauli` | 2025/05/03 | 対応のみ | $I_H=-i\sigma_1$ など（負号あり） | 変更不要 |
| `qua/spherical-trig` | 2025/11/29 | 四元数型 | $qvq^{-1}$ | 変更不要 |
| `lie/02-su2-so3` | 2026/07/10 | 四元数型 | $qxq^{-1}$、等傾回転。行列表現を標準形 $-i\sigma_k$ に改訂 | 改訂済み、Mathlog未反映 |
| `qua/04-4d-bsqua` | 2026/07/20 | 規約型 | $rvr^{-1}$、$r=\exp(-B/2)$。四元数は $r_Lqr_R$、$r_L=T(\varphi(r)^\dagger)$、$r_R=T(\varphi(r))^{-1}$ | 改訂済み（以前は $r^{-1}vr$、$r=\exp(B/2)$、逆元が $r_L$ 側）、Mathlog未反映 |
| `lie/03-spin` | 2026/07/24 | 規約型 | $RxR^{-1}$、$R=vu=\exp(-\varphi\sigma_1\sigma_2)$。四元数は $qxq^{-1}$、対応は $-i\sigma_k$ | 改訂済み（以前は対応が添字逆順、その後 $R^{-1}xR$）、Mathlog未反映 |

### 未公開

| 記事 | 型 | 作用・対応 | 状況 |
|---|---|---|---|
| `vec-oct/convention` | 規約型 | 規約そのものの説明。$qvq^{-1}$ の回転子をそのまま移して $Rx\tilde R$ | 書き直し済み |
| `vec-oct/02-rotation` | 規約型 | $rvr^{-1}$、$r=nm=\exp(-\frac\theta2p)$。四元数は $rvr^*$（$r=-nm$）、4Dは $r_Lqr_R$（$r_L=T(r^\dagger)$、$r_R=T(\tilde r)$） | 改訂済み |
| `lie/04-spin4` | 規約型 | $gvg^{-1}$。四元数は $qxq^{-1}$、$pxq^{-1}$ | 断り書きを簡略化 |
| `lie/05-bch-adjoint`・`lie/07-s7`・`lie/08-g2`・`lie/double-cover` | 四元数型 | $\rho_q(x)=qxq^{-1}$ | 変更不要 |
| `hopf/06-clifford-gates` | 規約型 | $R_{\boldsymbol n}(\theta)VR_{\boldsymbol n}^\dagger$（`qua/01` の $R$ そのもの） | 断り書きを改訂 |
| `qua/05-dual-qua`・`qua/06-slerp` | 四元数型 | 二重四元数による剛体変換、$qvq^{-1}$ の補間 | `qua/06` の `vec-oct/02` への断りを改訂 |
| `qua/cd/bcmp-qua` | 対応のみ | $p=a+bj$ から $i,j,k\mapsto i\sigma_3,i\sigma_2,i\sigma_1$（改訂前の `lie/02` と同じ論法）。標準形とは $i,j,k\mapsto-K_H,-J_H,-I_H$ の関係と断る | 構成上の表現として残す |
| `qua/history` | 規約型 | $R\boldsymbol v\tilde R$。四元数は $q\boldsymbol vq^*$ | 変更不要 |
| `em/05-lorentz` | 規約型 | $Rx\tilde R$、$R=e^{-i\sigma_3\theta/2}$（回転・ブースト） | 変更不要 |
| `dirac/01`〜`04` | 規約型 | $Rx\tilde R$、$\psi\gamma_\mu\tilde\psi=\rho e_\mu$。四元数との対応は $\mathbf i=-\omega\sigma_1$。`dirac/03` は座標変換で $x'=\tilde RxR$ も使う | `dirac/01` の「記号」の断りを改訂 |
| `clif/pga-cga` | 規約型 | $X\mapsto MX\widetilde M$、$M=rq$ | 変更不要 |
| `ktheory/06`〜`08` | 規約型 | $Ue_0\tilde U=\hat\theta$ など | 変更不要 |

`lie/10`・`lie/ladder-spinor` のパウリ行列は四元数との対応に依存しないため変更不要。

## 決定事項：ベクトルを数に移すときの左右

ベクトルを数（複素数・四元数）に移すとき、基準の生成元を左右どちらから掛けるかは、記事の目的に応じて次のように決めた。

- **回転・四元数と結び付ける2次元**（`geometric-product-exp`・`vec-oct/02`・`convention`）：右から掛ける $ve_1$（$i\cong e_2e_1$、$ba\cong a^*b$、回転子 $nm\cong m^*n$）。左から掛ける $e_1v$（$i\cong e_1e_2$、$ab\cong a^*b$）でも対応できることを述べたうえで、四元数の $k\cong e_2e_1$（右手系の要請）とそろえるために $e_2e_1$ を選ぶ、という論法で書く。
- **4次元の四元数への対応**（`qua/04` の $Q(v)=T(\varphi(e_1v))$）：左から掛ける $e_1v$ のまま。$\varphi$ の対応 $e_1e_2,e_1e_3,e_1e_4\mapsto\omega i,\omega j,\omega k$ に合わせると $e_1v$ が自然で、$ve_1$ にすると $Q$ が共役の四元数になり、左右の回転子も入れ替わる。
- **解析**（`clif-analysis` の $z=e_1\boldsymbol x$、$\omega=e_1e_2$、$e_1D=2\bar\partial$、`clif-analysis/04`・`ktheory`・`em/04` の $q=e_0\boldsymbol x$、$h_l=e_0e_l$）：左から掛けるまま（現状維持で決定）。左モノジェニックが正則にあたるため。$e_2e_1$ にすると右モノジェニックが正則にあたる鏡像になり、クリフォード解析の文献の標準（左モノジェニック）からも外れる。

## レビュー予定

規約の書き換えはまだ完成していない。書き換えが固まったら、次の記事を `uv run review` で見直す（`<stem>-prompt.md` の準備までは Claude が行い、実行は著者が判断する。手順は `REVIEW.md`）。既存の `<stem>.txt` は書き換え前の版に対するもの。

| 記事 | 状態 | 重点 |
|---|---|---|
| `vec-oct/convention` | 新規 | 四元数の回転子をそのまま移す筋立て、位相の負号と面の向きの説明、2次元で複素数に移す2通りの対応（$e_1v$ で $ab\cong a^*b$、$ve_1$ で $ba\cong a^*b$）と四元数に拡張するときに $i\cong e_2e_1$ を選ぶ論法、「別の書き方との関係」の節 |
| `vec-oct/02-rotation` | 全面改訂（既存 `.txt` は旧版） | 回転子 $r=nm$・$rvr^{-1}$、通常の鏡映 $-nvn$ への置き換え、2次元の複素数との2通りの対応と $e_2e_1$ を選ぶ理由、4次元の成分の対応（基底の反転）と $r_L=T(r^\dagger)$・$r_R=T(\tilde r)$ |
| `vec-oct/geometric-product-exp` | 公開済み・改訂 | $i\cong e_2e_1$、$ba\cong a^*b$、$e_1e_2$ でも対応できること（$ab\cong a^*b$）と四元数とそろえる理由の一文、rem「回転子との関係」、外積の符号付き面積 |
| `lie/02-su2-so3` | 公開済み・改訂 | 行列式の要請からの行列表現の導出、成分内・成分間の回転 |
| `lie/03-spin` | 公開済み・改訂 | 対応 $-i\sigma_k$、回転子 $R=vu$・$RxR^{-1}$、例とremの符号 |
| `qua/01-pauli-qua` | 公開済み・改訂 | $R=\exp(-\frac{i\theta}2N)\cong r$、$RVR^\dagger$、鏡映の節の $NM$ |
| `qua/04-4d-bsqua` | 公開済み・改訂（既存 `.txt` は旧版） | $r=\exp(-B/2)$・$rvr^{-1}$、$r_L=T(\varphi(r)^\dagger)$・$r_R=T(\varphi(r))^{-1}$、rem「直和成分と左右の回転子」 |

公開済みの記事は、レビューの指摘を反映してから Mathlog に反映するとよい。

## Mathlogへの反映

ローカルで改訂した公開済み記事のうち、Mathlogに未反映のものをまとめる。反映は `mathlog_fix.md` を書いて `bash src/mathlog_fix.sh` で行う（手順は `SLUG.md` の「更新手順」）。

- 参照（`[[slug]]`）を変えた記事は `--refs` で実行し、Mathlogの参考文献パネルも直したうえで、パネルをコピーして `refs/{ID}.toml` を取り込み直す。本文のみの修正は `--no-refs` でよい。オプションは実行単位で1つなので、両方がある場合は `mathlog_fix.md` を分けて2回実行するか、全体を `--refs` で実行する。
- 処理後に `mathlog_fix.md` を削除し、`make all` で `refs.toml` などを再生成する。
- `lie/02` の参照の変更により、現在 `reftools check` が `lie/02` について「未定義：`7shi-nonion`・`7shi-qp`」「未使用：`7shi-qcm`」を報告している。Mathlog側の参考文献パネルを直して取り込めば解消する。

### 反映の計画

参照の変更がある `lie/02` だけを先に `--refs` で反映し、残りの4本は本文のみなので、まとめて `--no-refs` で反映する。

1. **`lie/02`（`--refs`）**：リポジトリ直下の `mathlog_fix.md` に `lie/02` の1件だけを書いてある（本文の差し替え箇所と、参考文献パネルの `7shi-qcm` 削除・`7shi-qp` 追加・`7shi-nonion` 追加）。`bash src/mathlog_fix.sh --refs` で反映し、参考文献パネルを取り込んで `refs/Utdur1fLLzrWVHOJHifj.toml` を更新する。
2. **後始末**：`mathlog_fix.md` を削除し、`make all` で再生成する。`reftools check` の `lie/02` の警告が消えたことを確かめ、下の表から `lie/02` を外す。
3. **残りの4本（`--no-refs`）**：`qua/01`・`vec-oct/geometric-product-exp`・`qua/04`・`lie/03` を、下の下書きから `lie/02` を除いた内容で `mathlog_fix.md` に書き直し、`bash src/mathlog_fix.sh --no-refs` で反映する。レビュー予定の記事なので、レビューの指摘を反映してから行う。反映後に `mathlog_fix.md` を削除し、`make all` を実行する。

### 未反映の記事

| 記事 | Mathlog | 内容 | 参照の変更 |
|---|---|---|---|
| `qua/01-pauli-qua` | `lZ1X3t6exNS3NrNArqji` | パウリ行列による回転を、四元数の生成子をそのまま移した $RVR^\dagger$（$R=\exp(-\frac{i\theta}2N)\cong r$）に変更。鏡映の節を $NM$ に。改訂履歴を追加 | なし（`--no-refs`） |
| `vec-oct/geometric-product-exp` | `YZmxak6ObeP6rLQnyU2V` | 2次元の外積を符号付き面積として扱い、角度の制限を3次元に限定。複素数との対応を $i\cong e_2e_1$ に変更（$ve_1\leftrightarrow v_1+iv_2$、幾何積 $ba$ が $a^*b$ に対応）し、四元数の $k\cong e_2e_1$ との一致を追記。左から $e_1$ を掛けて $i\cong e_1e_2$ とする対応（$ab\cong a^*b$）もあり、四元数とそろえるために $e_2e_1$ を選ぶことを一文で追記。改訂履歴を追加 | なし（`--no-refs`） |
| `lie/02-su2-so3` | `Utdur1fLLzrWVHOJHifj` | 行列表現を行列式の要請から導く標準形に変更。$\mathfrak{su}(2)$ の分解、成分内・成分間の回転（$k$ と $j$）、等傾回転のremを追従。改訂履歴を追加 | `7shi-qp`・`7shi-nonion` を追加、`7shi-qcm` を削除（`--refs`） |
| `qua/04-4d-bsqua` | `asrMOxuKsJIdPqANfOs3` | 回転を $rvr^{-1}$、$r=\exp(-B/2)$ に変更。等傾回転の回転子、射影による証明（$r_L=T(\varphi(r)^\dagger)$、$r_R=T(\varphi(r))^{-1}$）、直和成分のrem、まとめを追従。改訂履歴を追加 | なし（`--no-refs`） |
| `lie/03-spin` | `DRbXTeeL31pDcZyG6ml7` | 四元数との対応を $-i\sigma_k$ に変更。スピノルの例を $q=\exp(k\theta/2)$ にして全角 $+\theta$ に。回転子による回転を $R=vu$、$RxR^{-1}$ に変更し、例・共役作用との対応のrem・まとめを追従。符号の選び方のremを追加。改訂履歴を追加 | なし（`--no-refs`） |

### `mathlog_fix.md` の下書き

```markdown
## qua/01-pauli-qua.md — https://mathlog.info/articles/lZ1X3t6exNS3NrNArqji
- 本文：改訂履歴、概要の5、「四元数からパウリ行列への変換」の後半（生成子の定義）、「パウリ行列による回転の表現」、「2次の元の回転」、「回転表現のまとめ」、「回転の生成子の分解」の積と回転表現を差し替え（参照の変更なし）

## vec-oct/geometric-product-exp.md — https://mathlog.info/articles/YZmxak6ObeP6rLQnyU2V
- 本文：改訂履歴、「指数関数による表現」の冒頭、「2次元の外積」の節、「複素数との対応」の節（rem「回転子との関係」を含む）、「外積における角度の範囲」のrem、「四元数との対応」の定義の直後、まとめを差し替え（参照の変更なし）

## lie/02-su2-so3.md — https://mathlog.info/articles/Utdur1fLLzrWVHOJHifj
- 本文：改訂履歴、前提、「四元数の行列表現」の節、su(2)の分解、「SO(2)の複素ユニタリ化」の後半、「等傾回転の相殺」の分解とremを差し替え
- 7shi-qcm: 削除
- 7shi-qp: 追加
- 7shi-nonion: 追加

## qua/04-4d-bsqua.md — https://mathlog.info/articles/asrMOxuKsJIdPqANfOs3
- 本文：改訂履歴、概要、「SO(4)とCl_{4,0}(R)」の回転と回転子、単純回転の例、「等傾回転の合成」の回転子、「純虚四元数への変換」のQの式、公式「回転子の対応」、remの例、証明1・2、rem「直和成分と左右の回転子」、まとめと枠「回転子の対応」を差し替え（参照の変更なし）

## lie/03-spin.md — https://mathlog.info/articles/DRbXTeeL31pDcZyG6ml7
- 本文：改訂履歴、概要の前提（su(2)の一般形、対角行列の作用）、スピノルの冒頭と例、「四元数との対応」の節、クリフォード代数の冒頭の対応、「単位ベクトルの積による構成」の回転子による回転・例・rem、偶部分代数の乗積規則とSpin(3)、グレード1への共役作用のrem、まとめを差し替え（参照の変更なし）
```
