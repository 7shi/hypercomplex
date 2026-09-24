---
refs:
  - qua/04-4d-bsqua.md
  - qua/spherical-trig.md
---

添付は「四元数」シリーズ第5回の記事です。球面三角法の記事で観察された「無限小の領域では回転が並進として扱える」を、$\varepsilon^2=0$ の形式的基底による厳密な読み替えに昇格させ、二重数・二重四元数を定義して3次元剛体変換 $\mathrm{SE}(3)$ の挟み込み表現とねじ運動を導き、最後に退化クリフォード代数 $\operatorname{Cl}_{3,0,1}(\mathbb R)$ の偶部分代数との同型を示します。読者は四元数の共役・指数関数・共役作用による回転を既知とします。以下の観点でレビューしてください。

- 数式・数学的主張の正しさ。特に次の導出を検算してください。
    - 角を $\varepsilon$ 倍した回転子 $e^{k\varepsilon a/2}$ が接平面の点 $i + \varepsilon x$ を $aj$ だけ並進させること、および $e^{j\varepsilon b/2}$ が $-bk$ 方向になること（符号）。
    - 剛体変換のサンドイッチ公式 $\sigma(1+\varepsilon x)\overline{\sigma^*} = 1 + \varepsilon(rx\overline r + t)$ の展開と、$\overline\sigma$ で挟むと並進が消えるという対比計算。
    - 合成則 (2) と、$\mathrm{SE}(3)$ の合成規則 $x \mapsto r_2(r_1x\overline{r_1}+t_1)\overline{r_2}+t_2$ の一致。
    - 同型表 $i \leftrightarrow e_3e_2$、$j \leftrightarrow e_1e_3$、$k \leftrightarrow e_2e_1$、$\varepsilon \leftrightarrow e_0e_1e_2e_3$ が、$i^2=j^2=k^2=ijk=-1$ と $\varepsilon$ の可換性・冪零性をすべて再現するか（表の各行各列の対応 $\varepsilon i \leftrightarrow e_0e_1$ などの符号を含めて）。
    - 擬スカラー $E$ が偶部分代数の中で中心的だという主張（1ベクトルと反交換、2ベクトルと可換）。
    - 二重角とねじ軸による指数表示 (3)、プリュッカー座標の $m = x_0 \times l$ と $x_0 = l \times m$ の関係。
- 説明の分かりやすさ、論理の飛躍の有無。特に、基点を $i$ から $1$ に取り替える箇所で、なぜ通常の共役 $\overline\sigma$ では駄目で第3の共役 $\overline{\sigma^*}$ が要るのかという動機づけが、読者の側から追える順序になっているか。
- 用語・記法の一貫性。参照記事との照合を含めて確認してください。
    - $\varepsilon$ の衝突。本稿の $\varepsilon$（$\varepsilon^2=0$ の形式的基底）、球面三角法の記事の $\varepsilon$（無限小の展開パラメーター）、4次元回転の記事の $\varepsilon = \frac{1+\omega}2$（冪等元）の3つを「記号について」で区別していますが、その説明が参照記事の実際の用法と合っているか。
    - 擬スカラーの三分法の表（$\operatorname{Cl}_{3,1}$・$\operatorname{Cl}_{4,0}$・$\operatorname{Cl}_{3,0,1}$ と $\omega^2$）が、4次元回転の記事の擬スカラー $\omega$、分解型双四元数、冪等元 $\frac{1\pm\omega}2$ による $\mathbb H \oplus \mathbb H$ 分解、$\mathrm{Spin}(4)$ の記述と整合しているか。符号数の記法（$(p,q)$ と $(p,q,r)$ の混在）、$\mathbb C'$ の記号も確認してください。
    - 内積 $\langle q,r \rangle := \mathrm{Re}(q\overline r)$、共役 $\overline{\phantom q}$・$*$ の使い分けが、参照記事の記法と衝突していないか。
- 日本語表現の自然さ。あわせて次のプロジェクト方針に反する箇所を指摘してください。
    - `——`（em dash）を使わない。箇条書きの見出しと説明の区切りは全角コロン「：」。
    - 「正体」「種明かし」のような砕けた比喩を避ける。
- `&&&` ブロックの type の使い分けが内容に合っているか。単なる `&&&rem` に丸めず、`&&&def`・`&&&thm`・`&&&fml`・`&&&ex`・`&&&prf` を使い分ける方針です。無題の `&&&rem` が続く箇所や、定理に見出しが付いていない箇所が適切か。

指摘事項を箇条書きで挙げてください。問題がなければその旨を書いてください。
