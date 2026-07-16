Wikipedia英語版のオイラー＝ラグランジュ方程式[[wiki-en-eleq]]の導出に使われている表記法が解析力学で馴染みのある形式とは異なるため、解析力学に寄せて解釈してみました。

# 作用積分

自然界で実現する運動は、作用$S$と呼ばれる物理量を極値（多くの場合は最小値）にする経路$q(t)$に沿って起こります。これは**最小作用の原理**として知られ、物理学において基本的な原理となっています。

この原理を数学的に定式化するため、時刻$t_1$と$t_2$の間の運動を考えます。境界条件$q(t_1) = q_1,\ q(t_2) = q_2$によって、出発点と到達点を固定化します。

![経路](/uploads/mathdown/F8VnX2adSOAanrv55B27.png =300)

ある時刻における系の状態は、ラグランジアン$L(q,\dot q,t)$という2回微分可能な関数で記述されます。例えば、ボールが投げ上げられた様子を想像してみてください。ある瞬間を切り取って「スナップショット」を撮れば、ボールの位置$q$と速度$\dot{q}$が写ります。$L(q,\dot q,t)$はこれらをパラメーターとして受け取って、系の状態をエネルギーとして計算します。

各時刻におけるラグランジアンの時間積分は**作用積分**と呼ばれます。実現する経路$q(t)$は、作用積分を極値化します。

&&&def 作用積分
$$
S[q] = \int_{t_1}^{t_2} L(q(t),\dot{q}(t),t)\, dt
$$
&&&

つまり、ラグランジアンで各時刻における状態を独立に計算し、それらを作用積分で繋ぎ合わせることで運動を再現するというアプローチです。

&&&rem 作用積分と力積
作用積分は、ラグランジアンで計算されるエネルギーの時間積分です。力の時間積分である力積に似ていますが、エネルギーなら何でも良いわけではなく、ラグランジアンに限定されます。
&&&

# ラグランジアンの微分

もし$q(t)$が境界条件の下で作用積分の極値を与えるならば、経路の変化における作用積分の変化率は$0$になります。

&&&rem
放物線（二次関数）において、頂点における接線の傾き（微分値）が$0$になることと同じような状況を想定しています。
&&&

変化後の経路を$q(t) + \varepsilon \eta(t)$とします。ここで、$\varepsilon$は任意の実数で、$\eta(t)$は境界条件$\eta(t_1) = \eta(t_2) = 0$を満たす任意の微分可能関数です。

![経路の変化](/uploads/mathdown/iVeKtVR1PcLMntPg3WxS.png =300)

&&&rem
$\varepsilon$は任意の実数としていますが、計算に寄与するのは$0$付近の微小値です。
&&&

経路の変化に伴うラグランジアン$L(q+\varepsilon\eta,\dot{q}+\varepsilon\dot{\eta},t)$の変化を調べるため、まず全微分を計算します。

$$
\begin{aligned}
dL &= \frac{\partial L}{\partial q}d(q+\varepsilon\eta) + \frac{\partial L}{\partial \dot{q}}d(\dot{q}+\varepsilon\dot{\eta}) + \frac{\partial L}{\partial t}dt \\
&= \frac{\partial L}{\partial q}(dq + \eta d\varepsilon + \varepsilon d\eta) + \frac{\partial L}{\partial \dot{q}}(d\dot{q} + \dot{\eta} d\varepsilon + \varepsilon d\dot{\eta}) + \frac{\partial L}{\partial t}dt
\end{aligned}
$$

次に、この式の両辺を$d\varepsilon$で割ります。$L$を$\varepsilon$で微分するということです。

$$
\frac{dL}{d\varepsilon} = \frac{\partial L}{\partial q}\left(\frac{dq}{d\varepsilon} + \eta + \varepsilon\frac{d\eta}{d\varepsilon}\right) + \frac{\partial L}{\partial \dot{q}}\left(\frac{d\dot{q}}{d\varepsilon} + \dot{\eta} + \varepsilon\frac{d\dot{\eta}}{d\varepsilon}\right) + \frac{\partial L}{\partial t}\frac{dt}{d\varepsilon}
$$

$q,\eta,\dot q,\dot\eta,t$は$\varepsilon$には依存しないことから、微分が$0$になります。

$$
\frac{dq}{d\varepsilon}
= \frac{d\eta}{d\varepsilon}
= \frac{d\dot q}{d\varepsilon}
= \frac{d\dot \eta}{d\varepsilon}
= \frac{dt}{d\varepsilon}
= 0
$$

これらを適用すれば、ラグランジアンの$\varepsilon$による微分は以下のようになります。

$$
\frac{dL}{d\varepsilon} = \frac{\partial L}{\partial q}\eta + \frac{\partial L}{\partial \dot{q}}\dot{\eta}
$$

# 作用積分の微分

経路変化後の作用積分

$$
S[q+\varepsilon\eta] = \int_{t_1}^{t_2} L(q(t)+\varepsilon\eta(t), \dot{q}(t)+\varepsilon\dot{\eta}(t), t)\, dt
$$

の$\varepsilon$についての微分は以下のように計算されます。

$$
\begin{aligned}
\frac{d}{d \varepsilon} S[q+\varepsilon\eta]
&= \frac{d}{d\varepsilon}\int_{t_1}^{t_2} L(q(t)+\varepsilon\eta(t), \dot{q}(t)+\varepsilon\dot{\eta}(t), t) \, dt \\
&= \int_{t_1}^{t_2} \frac{d}{d\varepsilon} L(q(t)+\varepsilon\eta(t), \dot{q}(t)+\varepsilon\dot{\eta}(t), t) \, dt \\
&= \int_{t_1}^{t_2} \left[
  \eta(t) \frac{\partial}{\partial {q} }L(q(t)+\varepsilon\eta(t), \dot{q}(t)+\varepsilon\dot{\eta}(t), t) \right. \\
&\qquad\quad \left. + \dot{\eta}(t) \frac{\partial}{\partial \dot{q}}L(q(t)+\varepsilon\eta(t), \dot{q}(t)+\varepsilon\dot{\eta}(t), t)
\right] dt
\end{aligned}
$$

$\varepsilon = 0$のとき、$S[q+\varepsilon\eta]$は極値を取るため

$$
\left.\frac{d}{d \varepsilon} S[q+\varepsilon\eta]\right|_{\varepsilon=0}
= \int_{t_1}^{t_2} \left[ \eta(t) \frac{\partial}{\partial q}L(q(t),\dot{q}(t),t) + \dot{\eta}(t) \frac{\partial}{\partial \dot{q}}L(q(t),\dot{q}(t),t) \,\right]\,dt = 0
$$

以降、$L(q(t),\dot{q}(t),t)$の引数は変化しないため省略します。

$$
\left.\frac{d}{d \varepsilon} S[q+\varepsilon\eta]\right|_{\varepsilon=0}
= \int_{t_1}^{t_2} \left[
  \eta(t) \frac{\partial L}{\partial q}
+ \color{red}{\dot{\eta}(t) \frac{\partial L}{\partial \dot{q}}}\color{black}{}
\right]\,dt = 0
$$

積分の第2項を部分積分します。[[7shi-ibp]]

$$
\int_{t_1}^{t_2}
\color{red}{\dot{\eta}(t) \frac{\partial L}{\partial \dot{q}}}\color{black}{}
\,dt
= \left[ \eta(t) \frac{\partial L}{\partial \dot{q}} \right]_{t_1}^{t_2} - \int_{t_1}^{t_2} \eta(t) \frac{d}{dt} \frac{\partial L}{\partial \dot{q}}\,dt
$$

境界条件$\eta(t_1) = \eta(t_2) = 0$より右辺第1項（表面項）は消えます。

$$
\int_{t_1}^{t_2} \dot{\eta}(t) \frac{\partial L}{\partial \dot{q}} \,dt
= - \int_{t_1}^{t_2} \eta(t) \frac{d}{dt} \frac{\partial L}{\partial \dot{q}}\,dt
$$

これを元の式に適用すれば

$$
\left.\frac{d}{d \varepsilon} S[q+\varepsilon\eta]\right|_{\varepsilon=0}
= \int_{t_1}^{t_2} \eta(t) \left[ \frac{\partial L}{\partial q} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}} \right] dt = 0
$$

$\eta(t)$は任意関数であることから、括弧の中が$0$になることが分かります。これにより**オイラー＝ラグランジュ方程式**が得られます。

&&&fml オイラー＝ラグランジュ方程式
$$
\frac{\partial L}{\partial q} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}} = 0
$$
&&&

この方程式は、$L = T - V$ （運動エネルギーから位置エネルギーを引いたもの）の場合、ニュートンの運動方程式$F=ma$と等価になります。

$$
\begin{aligned}
\frac{\partial L}{\partial q} &= \frac{\partial V}{\partial q} = F \\
\frac{d}{dt} \frac{\partial L}{\partial \dot{q}} &= \frac{d}{dt} \frac{\partial T}{\partial \dot{q}} = \frac{d}{dt} m \dot{q} = ma
\end{aligned}
$$

&&&rem
直感的に言えば、$L=T-V$によって表されるラグランジアンは、運動エネルギーと位置エネルギーのバランスを取るように定義されています。できるだけ位置エネルギーが低いところを通って、かつ、無駄に速く動かないような経路が選ばれます。その結果、運動方程式が再現されます。

$H=T+V$によって表される全エネルギーとは異なる概念であることに注意が必要です。なお、$H$から運動方程式を再現するのは正準方程式です。
&&&

# 変数の依存性に関する考察

$\varepsilon$パラメーターを導入することで、位置座標$q$とその時間微分$\dot{q}$の関係が適切に保たれており、変分法の導出においてしばしば混乱が生じる独立変数であるかのような扱いを回避しています。

&&&rem
依存関係のある$q$と$\dot q$で別々に偏微分するのが納得いかない場合は、こちらの記事を参照してください。[[7shi-diff]]
&&&

具体的には、変化後の経路を$q(t) + \varepsilon\eta(t)$とすることで、変化後の速度はその時間微分として自動的に決まります。

$$
\frac{d}{dt}(q(t) + \varepsilon\eta(t)) = \dot{q}(t) + \varepsilon\dot{\eta}(t)
$$

また、変分の微小性を定量的に扱うことができます。

$$
\left.\frac{d}{d\varepsilon}S[q + \varepsilon\eta]\right|_{\varepsilon=0}
$$

物理学でよく用いられる変分法の計算と比較します。

$$
\begin{aligned}
\left.\frac{d}{d\varepsilon}S[q + \varepsilon\eta]\right|_{\varepsilon=0}
&= \left. \int_{t_1}^{t_2} \frac{d}{d\varepsilon} L(q(t)+\varepsilon\eta(t), \dot{q}(t)+\varepsilon\dot{\eta}(t), t) \, dt \right|_{\varepsilon=0} \\
&= \int_{t_1}^{t_2} \left[ \frac{\partial L}{\partial {q}}\eta + \frac{\partial L}{\partial \dot{q}} \dot{\eta} \right] dt \\
&= \int_{t_1}^{t_2} \left[ \frac{\partial L}{\partial {q}} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}} \right] \eta \, dt = 0
\end{aligned} \tag{1}
$$
$$
\begin{aligned}
\delta S
&= \int_{t_1}^{t_2} \left[ L(q(t)+\delta q(t), \dot{q}(t)+\delta\dot{q}(t), t) - L(q(t), \dot{q}(t), t)
\right] dt \\
&\approx \int_{t_1}^{t_2} \left[ \frac{\partial L}{\partial q} \delta q + \frac{\partial L}{\partial \dot{q}} \delta \dot{q} \right] dt \\
&= \int_{t_1}^{t_2} \left[ \frac{\partial L}{\partial q} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}} \right] \delta q \, dt = 0
\end{aligned} \tag{2}
$$

- $(1)$では微分によって$\varepsilon\eta,\varepsilon\dot{\eta}$から$\varepsilon$が取り除かれて$\eta,\dot{\eta}$が残ります。
- 形式的には$(1)$の$\eta,\dot{\eta}$と$(2)$の$\delta q,\delta \dot q$が対応しますが、前者は任意の関数であり、微小であるという制約は課されません。既に見たように、それ自体はオイラー＝ラグランジュ方程式の導出に影響しません。
- $(1)$ではラグランジアンを$\varepsilon$によって微分します。$(2)$ではテイラー展開から高次の微小量を取り除くため、近似として$\approx$となっていますが、これは$(1)$における$\varepsilon=0$に対応します。
