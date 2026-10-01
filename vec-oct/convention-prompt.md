---
refs:
  - vec-oct/geometric-product-exp.md
  - qua/01-pauli-qua.md
---

添付は「ベクトルから八元数まで」シリーズの番外記事で、クリフォード代数と四元数で回転を書くときの向きと積の順序について、シリーズで採る規約とその整合性をまとめたものです。四元数の回転子$q=e^{\frac\theta2k}$を双四元数経由でそのまま$\operatorname{Cl}_{3,0}(\mathbb R)$に移すと$R=e^{\frac\theta2e_2e_1}=e^{-\frac\theta2e_1e_2}$、作用は$Rx\tilde R$になる、という筋立てです。参照記事`vec-oct/geometric-product-exp.md`は幾何積の指数関数表示$ab=|a||b|e^{\theta p}$と複素数・四元数との対応（$i\cong e_2e_1$、$k\cong e_2e_1$）、`qua/01-pauli-qua.md`は双四元数を介したパウリ行列と四元数の対応と、パウリ行列による回転$RVR^\dagger$を扱う公開済み記事です。以下の観点でレビューしてください。

- 数式・数学的主張の正しさ
    - 2ベクトルと四元数の対応$i\cong e_3e_2$、$j\cong e_1e_3$、$k\cong e_2e_1$（$\cong-i\sigma_k$）、$(e_2e_3,e_3e_1,e_1e_2)$では$ij=-k$となること、$Rx_\parallel\tilde R=R^2x_\parallel=x_\parallel e^{\theta p}$、$(nm)m=n$、$-n(-mvm)n=Rv\tilde R$の符号
    - 「2次元と複素数」の節：$e_1v$（$i\cong e_1e_2$）では$ab\leftrightarrow a^*b$、$ve_1$（$i\cong e_2e_1$）では$ba\leftrightarrow a^*b$となること、$e_1(Rv\tilde R)=\tilde R^2(e_1v)$、$(Rv\tilde R)e_1=R^2(ve_1)$、回転子$R$がそれぞれ$e^{-i\theta/2}$、$e^{i\theta/2}$に移ること
    - 「別の書き方との関係」：$r=\tilde R$が$q^{-1}$の像であること、合成・ベクトルの積・鏡映の読みがすべて逆になること
- 説明の分かりやすさ、論理の飛躍の有無
    - 位相の負号の説明（回転面の向き$p$と、左から反時計回りを与える$\tilde p$の食い違い）が、ウェッジ積$m\wedge n$を「$m$から$n$へ」と読む約束と矛盾なくつながっているか
    - 2次元では$i\cong e_1e_2$の方が素直だと認めたうえで、四元数に拡張するときに$i\cong e_2e_1$を選ぶ論法に説得力があるか
- 用語・記法の一貫性
    - `vec-oct/geometric-product-exp.md`との照合：$i\cong e_2e_1$を選ぶ理由、$ba\cong a^*b$、回転子$nm$の説明、外積の符号付き面積と角度の範囲（2次元で$-\pi<\theta\le\pi$、3次元以上で$0\le\theta\le\pi$）
    - `qua/01-pauli-qua.md`との照合：双四元数の虚数単位$h$、$e_k\cong\sigma_k\cong hi,hj,hk$、$R=\exp(-\frac{i\theta}2N)\cong r$と本記事の$R$の関係
- 日本語表現の自然さ
- 構成上の制約
    - 規約の整理を目的とする記事で、前提を基礎から説明しない方針です。前提の範囲（概要の箇条書き）と本文で実際に使っている事項が対応しているか
    - 後続記事を予告する書き方が残っていないか

指摘事項を箇条書きで挙げてください。問題がなければその旨を書いてください。
