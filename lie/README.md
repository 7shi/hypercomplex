# リー群とリー代数

リー群とリー代数について、数学的厳密さよりも直観性を重視して説明する記事シリーズです。

## 記事

- [01-u1-so2.md](01-u1-so2.md) — 複素数から構成されるリー群$\operatorname{U}(1)$を例に、リー群とリー代数の基本概念（指数写像、基底、括弧積）を説明します。複素数の行列表現を通じて$\operatorname{SO}(2)$との同型を示します。
- [02-su2-so3.md](02-su2-so3.md) — 四元数から構成されるリー群$\operatorname{Sp}(1) \cong \operatorname{SU}(2)$。非可換性の反映として初めて$0$でない括弧積が現れ、外積と一致することを見ます。行列表現の形が「絶対値の2乗＝行列式」の要求から導かれること（シルベスター流の導出）、リー群の条件を微分して$\mathfrak{su}(2)$（トレース$0$の反エルミート行列）を求め基底を得ること、共役作用$qxq^{-1}$による$\operatorname{SO}(3)$の2重被覆、リー代数の同型$\mathfrak{su}(2) \cong \mathfrak{so}(3)$について。
- [03-spin.md](03-spin.md) — 四元数の裏の顔であるクリフォード代数と$\operatorname{Spin}$群。$\operatorname{SO}(3)$の2重被覆としての$\operatorname{Spin}(3)$。$\mathfrak{su}(2)$（反エルミート行列）から$i$をくくり出すことによるパウリ行列の導出、3次元空間のベクトルとしてのパウリ行列（グレード1）と、その積（グレード2）が四元数になること、偶部分代数。さらに4次元への拡張として、独立した左右からの作用$pxq^{-1}$が等角回転の制限を解除して$\operatorname{SO}(4)$の2重被覆$\operatorname{Spin}(4) \cong \operatorname{Sp}(1) \times \operatorname{Sp}(1)$になることと、リー代数の直和分解$\mathfrak{so}(4) \cong \mathfrak{sp}(1) \oplus \mathfrak{sp}(1)$。八元数における非結合性の伏線となる「左作用と右作用の分離」について。
- [04-su3.md](04-su3.md) — $\operatorname{SO}(2) \to \operatorname{SO}(3)$の構造（座標平面ごとの$\mathfrak{so}(2)$の束、対角自由度なし）を分析し、それを複素化して$\operatorname{SU}(3)$を構成します。数体系に頼らない初めてのリー群です。実と複素の次元勘定の対比、ゲルマン行列と$\mathfrak{so}(3) \subset \mathfrak{su}(3)$、リー群とリー代数が同じ基底の空間に同居する$\operatorname{SU}(2)$までの縮退が解けて両者が分離すること、構造定数、対角位相によるランクの増加、平面回転が$\mathfrak{su}(2)$に太ることによる編み込み構造、$\mathbb{C}^3$への作用と共役作用、対応する数体系が存在しないこと（フルヴィッツの定理）について。

## 今後の予定

- 05 — BCH公式と随伴作用。指数法則のずれを括弧積の言葉で定式化するベイカー・キャンベル・ハウスドルフの公式と、共役作用のリー代数版である随伴作用$\operatorname{Ad}/\operatorname{ad}$。
- 06 — 八元数と$S^7$。積と逆元を持ちながら結合性が崩れるためリー群にならない球面。ムーファン・ループ、フルヴィッツの定理。
- 07 — $G_2$：八元数の自己同型群。リー代数としての微分全体、$G_2 \subset \operatorname{SO}(7)$、虚数単位を固定する部分群として$\operatorname{SU}(3)$が再登場すること。
- 08 — $\operatorname{Spin}(7)$と三対性。左作用が生成するクリフォード代数とスピノル表現、$\operatorname{Spin}(8)$の三対性。ルートの八元数記事との合流点。

## 記事構成の意図

本シリーズは、単にリー群を羅列するのではなく、以下のような全体を貫く構造的なストーリーを持っています。

- **数体系によるパラレルの終着**: 01と02では、「複素数 $\to \operatorname{SO}(2)$」と「四元数 $\to \operatorname{SO}(3)$」という美しい並行構造を描きます。しかし、この「数体系から直接回転群を作る」という単純なパターンは八元数で破綻します（結合則の喪失によりリー群にならないため）。
- **一般化へのバトンタッチ**: そこで03において、「クリフォード代数」と「$\operatorname{Spin}$群」を導入します。これが高次元の回転群を記述する枠組みへのバトンタッチとなります。また、ここで導入する四元数の「左右の独立した作用」は、のちの八元数における非結合性への最大の伏線です。
- **拡張の軸の転換**: 数体系の道が途絶えるため、04の $\operatorname{SU}(3)$ では「行列のサイズを上げる」という別の軸へ乗り換えて群を拡張します。
- **壮大な伏線回収**: 後半の06以降では八元数の対称性を扱います。結合律の崩れによって左右の作用が複雑に絡み合い、そこにクリフォード代数が生成する $\operatorname{Spin}(6)$ や $\operatorname{Spin}(7)$ といった群が立ち現れます。そして、その部分群として前半で構成した $\operatorname{SU}(3)$ や $G_2$ が再登場し、すべての物語が統合されます。

## 検証コード

プロジェクトルートで `uv sync` により環境を構築し、`uv run lie/<ファイル名>` で実行します。

- [check-02-quaternion.py](check-02-quaternion.py) — 02-su2-so3.md の数式の数値検証。四元数の乗法表と2×2複素行列表現、シルベスター条件（トレース$0$・行列式$1$・互いに反交換）、オイラーの公式の四元数版、指数法則の反例、括弧積と外積の関係$[p,q]=2(p\times q)$、$\operatorname{SU}(2)$条件と$\mathfrak{su}(2)$の一般形$i(x\sigma_1 + y\sigma_2 + z\sigma_3)$、共役作用による回転と2重被覆、$\mathfrak{su}(2)\cong\mathfrak{so}(3)$の構造定数の一致を確認します。
- [check-03-spin.py](check-03-spin.py) — 03-spin.md の数式の数値検証。$\mathfrak{su}(2)$の一般形が「エルミート行列の$i$倍」となりパウリ行列に分解されること、$i \Leftrightarrow i\sigma_3,\ j \Leftrightarrow i\sigma_2,\ k \Leftrightarrow i\sigma_1$の対応、パウリ行列の長さ・直交性（$\sigma_a^2 = I$・反交換）、グレード2の因数分解$i = \sigma_1\sigma_2,\ j = \sigma_3\sigma_1,\ k = \sigma_2\sigma_3$と$(\sigma_1\sigma_2)^2 = -I$、偶部分代数が積について閉じ奇数グレードを含まないこと、両側作用$pxq^{-1}$の平面ごとの回転角$\alpha \mp \beta$、一般の$p, q$で4×4直交行列（$\det = 1$）となること、$(-p, -q)$が同じ回転を与えること（2重被覆）、左作用と右作用の生成子が可換で合わせて6次元（$\mathfrak{so}(4) \cong \mathfrak{sp}(1) \oplus \mathfrak{sp}(1)$）となることを確認します。
- [check-04-su3.py](check-04-su3.py) — 04-su3.md の数式の数値検証。ゲルマン行列の性質（エルミート・トレース$0$・正規化$\operatorname{tr}(\lambda_a\lambda_b)=2\delta_{ab}$）、$L_x, L_y, L_z$が$\mathfrak{so}(2)$生成子のブロック配置であること、パウリ行列の左上ブロック埋め込み、$i\lambda_{2,5,7}$が$\pm L$に一致すること（$\mathfrak{so}(3) \subset \mathfrak{su}(3)$）、括弧積の閉性、$\exp$が$\operatorname{SU}(3)$・$\operatorname{SO}(3)$に入ること、$(i\lambda_1)^2 \ne -I$および$\exp(i\theta\lambda_1)$が$\{I, i\lambda_a\}$の実線形結合で書けないこと（$\operatorname{SU}(2)$の縮退との対比）、構造定数$f_{abc}$の値と完全反対称性、$[\lambda_3,\lambda_8]=0$、3つの$\mathfrak{su}(2)$の組（I・V・Uスピン）、$\operatorname{SO}(2) \subset \operatorname{SO}(3)$と$\operatorname{U}(1) \subset \operatorname{SU}(2) \subset \operatorname{SU}(3)$のブロック埋め込み、共役作用が$\mathfrak{su}(3)$を保ち8次元の回転を与えることを確認します。
