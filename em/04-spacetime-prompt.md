---
refs:
  - em/02-maxwell.md
  - clif-analysis/04-quaternion.md
---

添付は「幾何代数による電磁気学」シリーズ第4回の記事です。第1〜3回の$\operatorname{Cl}_{3,0}(\mathbb R)$の形式を、時空代数$\operatorname{Cl}_{1,3}(\mathbb R)$（$\gamma_0^2=1$、$\gamma_k^2=-1$）の偶部分として回収し、マクスウェル方程式を$DF=\mu_0cJ$と書きます。$\gamma_0$を掛ける操作を「クリフォード解析」第4回の$e_0$を掛ける操作と対比し、空間の生成元の2乗の符号が楕円型と双曲型を分けることをシリーズの軸として示します。参照記事として、$\mathcal D$・$\bar{\mathcal D}$と1本の式の出所である第2回（未レビューの下書き）と、$q=e_0\boldsymbol x$、$\mathcal D=e_0D$、$h_l=e_0e_l$、$D(e_0F)=\bar{\mathcal D}F$の出所である「四元数解析」を添付します。以下の観点でレビューしてください。

- 数式・数学的主張の正しさ
    - `&&&prop 偶部分代数と$\operatorname{Cl}_{3,0}(\mathbb R)$`の証明（符号の数え方、基底が「符号を除いて一致する」ことからの同型の結論）が十分か
    - `&&&prop 連続の式`の証明で、$\gamma^\mu\cdot(\gamma^\nu\cdot F)$の反対称性から$D\cdot(D\cdot F)=0$を導く部分が正確か
    - `&&&rem $\operatorname{Cl}_{1,3}(\mathbb R)$を選ぶ理由`の各項目が正確か。とくに$\operatorname{Cl}_{3,1}(\mathbb R)$でも$(\gamma_k\gamma_0)^2=+1$で偶部分は$\operatorname{Cl}_{3,0}(\mathbb R)$と同型であり、違いを$\gamma_0^2=-1$による$D^2$の全体の符号に求めている点、$\operatorname{Cl}_{1,3}\cong M_2(\mathbb H)$・$\operatorname{Cl}_{3,1}\cong M_4(\mathbb R)$と「$\operatorname{Cl}_{3,1}$の生成元は$\operatorname{Cl}_{1,3}$の生成元の$i$倍」という述べ方、ワイル表現のガンマ行列との対応
- 説明の分かりやすさ、論理の飛躍の有無
    - 逆元$\gamma^\mu$を係数に使う理由の説明と、$D\cdot F$・$D\wedge F$の2本への分解の表が分かりやすいか
- 用語・記法の一貫性
    - 第2回との照合：$\mathcal D$・$\bar{\mathcal D}$、$\mathcal DF=(\rho-\boldsymbol J/c)/\varepsilon_0$、$\bar{\mathcal D}\mathcal D=\partial_0^2-\Delta$を$D^2$の言い換えとする説明
    - 「四元数解析」との照合：$h_l=e_0e_l$と$\sigma_k=\gamma_k\gamma_0$で掛ける順序が逆である点が混乱を招かないか、楕円型・双曲型の対比の表の記載（「共役との積」の列）が正確か
    - 前半で$\operatorname{Cl}_{3,0}(\mathbb R)$のベクトルを$e_k$と書いていたものを、本記事から$\sigma_k$と書く切り替えが明示されているか
- 日本語表現の自然さ
- プロジェクト方針
    - 後続記事の内容を予告する前方参照がないか（「分け方は時間方向$\gamma_0$の選び方に依存します」がローレンツ変換の予告になっていないか）
    - `&&&`ブロックの種類の使い分けと、本文・主張はですます調、`&&&prf`内だけである調になっているか

次の計算は別途、記号計算で検証済みです：逆基底と$Dx=4$・$D^2=\partial_0^2-\Delta$、$\sigma_k$の関係式と$\sigma_1\sigma_2\sigma_3=\gamma_0\gamma_1\gamma_2\gamma_3$、偶部分代数の閉性、$x\gamma_0=x_0+\boldsymbol x$、$\gamma_0D=\mathcal D$・$D(\gamma_0H)=\bar{\mathcal D}H$、$I\sigma_k$の基底表示、$\gamma_0J=c(\rho-\boldsymbol J/c)$、$\gamma_0(D\cdot F)$・$\gamma_0(D\wedge F)$の成分、$D\cdot F=\mu_0cJ$、$D\cdot(D\cdot F)=0$、ワイル表現の符号。検算よりも、論証の構成と述べ方の正確さを優先してください。

指摘事項を箇条書きで挙げてください。問題がなければその旨を書いてください。
