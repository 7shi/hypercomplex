実数パラメーター$\varepsilon$を用いた変分の定式化を通して、位置とその時間微分の独立性や微小量の扱いを明確にしながら、オイラー＝ラグランジュ方程式を厳密に導出します。

# 概要

変分法におけるオイラー＝ラグランジュ方程式の導出では、「仮想変位$\delta q$の微小性」や「位置$q$と速度$\dot{q}$の偏微分における独立性」があいまいなまま形式的に計算されやすく、直観的な理解や数学的な正当化において疑問が生じやすい論点となっています。英語圏の文献などでは、変分を実数パラメーター$\varepsilon$でパラメーター表示する定式化がしばしば用いられます。[[wiki-en-eleq]]

本記事では、経路の変分を実数パラメーター$\varepsilon$と、端点条件を満たす任意の$C^1$級関数$\eta(t)$の積$\varepsilon\eta(t)$として導入することで、変分の微小性を実数微分の極限$\varepsilon\to 0$として厳密に定式化します。まず作用積分とラグランジアンの定義を確認し、合成関数の微分法則（連鎖律）によってラグランジアンの$\varepsilon$微分を求めます。次に、作用積分の停留条件から部分積分と変分法の基本補題を用いてオイラー＝ラグランジュ方程式を導出します。最後に、この$\varepsilon$法と物理で慣用される$\delta$変分法の対応を比較し、微小量の扱いの差異を明確にします。

多変数関数の微分、合成関数の微分法則（連鎖律）、および部分積分法を前提とします。本記事では1次元の単一粒子系におけるオイラー＝ラグランジュ方程式の基礎的な導出に焦点を当て、拘束条件付き変分（ラグランジュの未定乗数法）や多自由度系・場の理論への一般化は扱いません。

# 作用積分

固定した始点と終点を結ぶ経路のうち、自然界で実現する運動の経路$q(t)$は、作用$S$と呼ばれる物理量を停留させる経路です。これは**作用の停留原理**と呼ばれ、慣用的には**最小作用の原理**とも呼ばれます。物理学において基本的な原理となっています。

この原理を数学的に定式化するため、時刻$t_1$と$t_2$（$t_1<t_2$）の間の運動を考えます。境界条件$q(t_1) = q_1,\ q(t_2) = q_2$によって、出発点と到達点を固定します。

![経路](/uploads/mathdown/F8VnX2adSOAanrv55B27.png =300)

ある時刻における系の状態は、ラグランジアン$L(q,\dot q,t)$という関数で記述されます。例えば、ボールが投げ上げられた様子を想像してみてください。ある瞬間を切り取って「スナップショット」を撮れば、ボールの位置$q$と速度$\dot{q}$が写ります。$L(q,\dot q,t)$はこれらをパラメーターとして受け取って、系の状態をエネルギーとして計算します。

&&&rem 独立変数としての速度引数
$L(q,\dot q,t)$の引数$\dot q$は、この段階では$q$の時間微分$dq/dt$と結び付けられているわけではなく、単に「$\dot q$」という名前を持つ独立な変数として扱われています。

静止画には速度が写らないと思われるかもしれません。ここでのスナップショットは、速度を位置の時間変化から読み取るのではなく、位置と並ぶパラメーターとしてその瞬間の状態に含めて与える、という立場を表しています。同じ位置にあっても速度は経路によって異なるため、位置と速度は独立に値を取ります。
&&&

各時刻におけるラグランジアンの時間積分は**作用積分**と呼ばれます。実現する経路$q(t)$は、作用積分を停留させます。以下、$L$は考える経路の近傍で$C^2$級、経路$q$は$[t_1,t_2]$上で$C^2$級とします。

&&&def 作用積分
$$
S[q] = \int_{t_1}^{t_2} L\bigl(q(t),\dot{q}(t),t\bigr)\, dt
$$
&&&

先に述べた通り、$L(q,\dot q,t)$における$q,\dot q,t$は、それぞれ独立な名前を持つ3つの引数にすぎません。一方、$S[q,\dot q]$ではなく$S[q]$と表記されているのは、$S$の引数が$q$（経路を表す関数）ただ1つであることを表しています。被積分関数$L\bigl(q(t),\dot{q}(t),t\bigr)$に現れる$q(t)$と$\dot{q}(t)$は、この1つの引数$q$を「時刻$t$で評価する」「微分してから時刻$t$で評価する」ことで得られる量であり、$\dot{q}(t) = dq(t)/dt$という関係を通じて、いずれも同じ引数$q$から導かれています。この2つの段階の違いを区別しておくと、後の議論が見通しやすくなります。

つまり、ラグランジアンで各時刻における状態を独立に計算し、それらを作用積分で繋ぎ合わせることで運動を再現するというアプローチです。

&&&rem 作用積分と力積
作用積分は、ラグランジアンで計算されるエネルギーの時間積分です。力の時間積分である力積に似ていますが、エネルギーなら何でも良いわけではなく、ラグランジアンに限定されます。
&&&

# ラグランジアンの微分

もし$q(t)$が境界条件の下で作用積分を停留させるならば、経路の変化における作用積分の変化率は$0$になります。

&&&rem 関数の極値と微分係数のアナロジー
放物線（2次関数）において、頂点における接線の傾き（微分値）が$0$になることと同じような状況を想定しています。ただし、傾きが$0$になることは極値の必要条件であり、傾きが$0$になる点（停留点）が極値とは限りません。
&&&

変化後の経路を$q_\varepsilon(t) = q(t) + \varepsilon \eta(t)$とします。ここで、$\varepsilon$は実数で、$\eta(t)$は境界条件$\eta(t_1) = \eta(t_2) = 0$を満たす任意の$C^1$級関数です。$\varepsilon$は任意の実数としてかまいませんが、評価するのは$\varepsilon=0$における微分係数のため、実際に使うのは$0$付近の値だけです。$|\varepsilon|$は十分小さく取り、変化後の経路が$L$の定義域内にあるものとします。

![経路の変化](/uploads/mathdown/iVeKtVR1PcLMntPg3WxS.png =300)

変化後の速度は、変化後の経路の時間微分として決まります。

$$
\dot q_\varepsilon(t) = \dot q(t) + \varepsilon\dot\eta(t)
$$

経路の変化に伴うラグランジアン$L\bigl(q_\varepsilon(t),\dot q_\varepsilon(t),t\bigr)$の変化を調べます。時刻$t$を固定すると、$q_\varepsilon(t)$を$\varepsilon$で微分すると$\eta(t)$、$\dot q_\varepsilon(t)$を$\varepsilon$で微分すると$\dot\eta(t)$になり、第3引数の$t$は$\varepsilon$に依存しません。したがって、合成関数の微分法則（連鎖律）により次のように微分できます。

$$
\frac{\partial}{\partial\varepsilon}L\bigl(q_\varepsilon(t),\dot q_\varepsilon(t),t\bigr)
= \frac{\partial L}{\partial q}\bigl(q_\varepsilon(t),\dot q_\varepsilon(t),t\bigr)\,\eta(t)
+ \frac{\partial L}{\partial \dot q}\bigl(q_\varepsilon(t),\dot q_\varepsilon(t),t\bigr)\,\dot\eta(t)
$$

&&&rem 偏微分の意味
ここで用いている$\partial L/\partial q$や$\partial L/\partial \dot{q}$は、$L$を独立な3引数の関数とみなした通常の偏微分を、括弧内の点で評価したものです。$\dot q = dq/dt$という関係は偏微分の段階では使っておらず、$q$と$\dot q$を独立な引数として扱うことに矛盾は生じません。
&&&

# 作用積分の微分

経路変化後の作用積分は次のように表されます。

$$
S[q+\varepsilon\eta] = \int_{t_1}^{t_2} L\bigl(q_\varepsilon(t), \dot{q}_\varepsilon(t), t\bigr)\, dt
$$

この作用積分の$\varepsilon=0$における微分係数を、**第一変分**と呼びます。

&&&def 第一変分
$$
\delta S[q;\eta] = \left.\frac{d}{d\varepsilon} S[q+\varepsilon\eta]\right|_{\varepsilon=0}
$$
&&&

$q$が作用積分を停留させるとは、端点条件を満たす任意の$\eta$について$\delta S[q;\eta]=0$となることです。

$\varepsilon$についての微分は以下のように計算されます。滑らかさの仮定により、被積分関数とその$\varepsilon$偏導関数は連続なので、微分と積分を交換できます。

$$
\begin{aligned}
\frac{d}{d \varepsilon} S[q+\varepsilon\eta]
&= \frac{d}{d\varepsilon}\int_{t_1}^{t_2} L\bigl(q_\varepsilon(t), \dot{q}_\varepsilon(t), t\bigr) \, dt \\
&= \int_{t_1}^{t_2} \frac{\partial}{\partial\varepsilon} L\bigl(q_\varepsilon(t), \dot{q}_\varepsilon(t), t\bigr) \, dt \\
&= \int_{t_1}^{t_2} \left[
  \eta(t) \frac{\partial L}{\partial q}\bigl(q_\varepsilon(t), \dot{q}_\varepsilon(t), t\bigr)
+ \dot{\eta}(t) \frac{\partial L}{\partial \dot{q}}\bigl(q_\varepsilon(t), \dot{q}_\varepsilon(t), t\bigr)
\right] dt
\end{aligned}
$$

$\varepsilon=0$とすれば、第一変分が得られます。

$$
\delta S[q;\eta]
= \int_{t_1}^{t_2} \left[ \eta(t) \frac{\partial L}{\partial q}\bigl(q(t),\dot{q}(t),t\bigr) + \dot{\eta}(t) \frac{\partial L}{\partial \dot{q}}\bigl(q(t),\dot{q}(t),t\bigr) \,\right]\,dt
$$

以降、偏導関数はすべて$(q(t),\dot q(t),t)$で評価するものとし、評価点の表記を省きます。

$$
\delta S[q;\eta]
= \int_{t_1}^{t_2} \left[
  \eta(t) \frac{\partial L}{\partial q}
+ \color{red}{\dot{\eta}(t) \frac{\partial L}{\partial \dot{q}}}\color{black}{}
\right]\,dt
$$

積分の第2項を部分積分します。[[7shi-ibp]]

$$
\int_{t_1}^{t_2}
\color{red}{\dot{\eta}(t) \frac{\partial L}{\partial \dot{q}}}\color{black}{}
\,dt
= \left[ \eta(t) \frac{\partial L}{\partial \dot{q}} \right]_{t_1}^{t_2} - \int_{t_1}^{t_2} \eta(t) \frac{d}{dt} \frac{\partial L}{\partial \dot{q}}\,dt
$$

$L$と$q$は$C^2$級なので、$\frac{\partial L}{\partial\dot q}$は$t$の$C^1$級関数となり、部分積分を適用できます。境界条件$\eta(t_1) = \eta(t_2) = 0$より右辺第1項（境界項）は消えます。

$$
\int_{t_1}^{t_2} \dot{\eta}(t) \frac{\partial L}{\partial \dot{q}} \,dt
= - \int_{t_1}^{t_2} \eta(t) \frac{d}{dt} \frac{\partial L}{\partial \dot{q}}\,dt
$$

これを第一変分の式に代入します。

$$
\delta S[q;\eta]
= \int_{t_1}^{t_2} \eta(t) \left[ \frac{\partial L}{\partial q} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}} \right] dt
$$

$q$が停留経路であれば、この値は端点条件を満たす任意の$\eta$について$0$になります。括弧の中は$t$の連続関数なので、次の補題により括弧の中が$0$になることが分かります。

&&&lem 変分法の基本補題
$A(t)$を$[t_1,t_2]$上の連続関数とします。$\eta(t_1)=\eta(t_2)=0$を満たす任意の$C^1$級関数$\eta$について次が成り立つならば、$[t_1,t_2]$上で$A(t)=0$です。

$$
\int_{t_1}^{t_2}A(t)\eta(t)\,dt=0
$$
&&&

&&&prf
ある$t_0\in(t_1,t_2)$で$A(t_0)>0$と仮定する。$A$は連続なので、$t_0$を含む区間$[a,b]\subset(t_1,t_2)$（$a<b$）で$A(t)>0$となる。$\eta$を次のように定める。

$$
\eta(t)=
\begin{cases}
(t-a)^2(b-t)^2 & (a\le t\le b) \\
0 & (\text{それ以外})
\end{cases}
$$

この$\eta$は$C^1$級で端点で$0$になり、$(a,b)$で正である。したがって$\int_{t_1}^{t_2}A\eta\,dt=\int_a^b A\eta\,dt>0$となり、仮定に反する。$A(t_0)<0$の場合も同様である。よって$(t_1,t_2)$上で$A=0$であり、連続性から端点でも$A=0$である。
&&&

これにより**オイラー＝ラグランジュ方程式**が得られます。

&&&fml オイラー＝ラグランジュ方程式
$$
\frac{\partial L}{\partial q} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}} = 0
$$
&&&

&&&ex ニュートンの運動方程式
$L = T - V$ （運動エネルギーから位置エネルギーを引いたもの）で、$T=\frac12m\dot q^2$、$V=V(q)$の場合を考えます。各項は次のようになります。

$$
\begin{aligned}
\frac{\partial L}{\partial q} &= -\frac{\partial V}{\partial q} = F \\
\frac{d}{dt} \frac{\partial L}{\partial \dot{q}} &= \frac{d}{dt} \frac{\partial T}{\partial \dot{q}} = \frac{d}{dt} m \dot{q} = ma
\end{aligned}
$$

これらをオイラー＝ラグランジュ方程式に代入すると$F-ma=0$となり、ニュートンの運動方程式$F=ma$と等価になります。
&&&

&&&rem ラグランジアンの物理的直観
$L=T-V$によって表されるラグランジアンは、運動エネルギーと位置エネルギーの差です。作用の停留条件を求めると、慣性の項と力の項が釣り合う運動方程式が得られます。これは、運動エネルギーや位置エネルギーをそれぞれ小さくするという原理ではありません。

$H=T+V$によって表される全エネルギーとは異なる概念であることに注意が必要です。なお、$H$から運動方程式を再現するのは正準方程式です。
&&&

# 変数の依存性に関する考察

$L$の偏微分では、位置・速度・時刻を独立な引数として扱います。一方、経路を代入した後の位置と速度は、$\dot q=dq/dt$で結ばれています。変化後の経路を$q_\varepsilon=q+\varepsilon\eta$と置くことで、速度の変化も$\dot q_\varepsilon=\dot q+\varepsilon\dot\eta$と定まり、この関係を保ったまま変分できます。

&&&rem 変数の依存関係の補足
ラグランジアンの独立な引数と、経路上での位置・速度の関係を区別する説明は、こちらの記事を参照してください。[[7shi-diff]]
&&&

また、変分の微小性は、第一変分$\delta S[q;\eta]$を定める$\varepsilon\to 0$の極限として、定量的に扱うことができます。

物理学でよく用いられる変分法の計算と比較します。物理学では、時刻を固定した経路の変化を$\delta q(t)$とし、端点条件$\delta q(t_1)=\delta q(t_2)=0$を課します。速度の変化はその時間微分$\delta\dot q=\frac{d}{dt}\delta q$で、作用の差の1次の部分を$\delta S$と書きます。停留経路では、次のようになります。

$$
\begin{aligned}
\delta S[q;\eta]
&= \left. \int_{t_1}^{t_2} \frac{\partial}{\partial\varepsilon} L\bigl(q(t)+\varepsilon\eta(t), \dot{q}(t)+\varepsilon\dot{\eta}(t), t\bigr) \, dt \right|_{\varepsilon=0} \\
&= \int_{t_1}^{t_2} \left[ \frac{\partial L}{\partial {q}}\eta + \frac{\partial L}{\partial \dot{q}} \dot{\eta} \right] dt \\
&= \int_{t_1}^{t_2} \left[ \frac{\partial L}{\partial {q}} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}} \right] \eta \, dt = 0
\end{aligned} \tag{1}
$$
$$
\begin{aligned}
S[q+\delta q]-S[q]
&= \int_{t_1}^{t_2} \left[ L\bigl(q(t)+\delta q(t), \dot{q}(t)+\delta\dot{q}(t), t\bigr) - L\bigl(q(t), \dot{q}(t), t\bigr)
\right] dt \\
&= \delta S + (\text{2次以上の項}) \\
\delta S
&= \int_{t_1}^{t_2} \left[ \frac{\partial L}{\partial q} \delta q + \frac{\partial L}{\partial \dot{q}} \delta \dot{q} \right] dt \\
&= \int_{t_1}^{t_2} \left[ \frac{\partial L}{\partial q} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}} \right] \delta q \, dt = 0
\end{aligned} \tag{2}
$$

- $(1)$では微分によって$\varepsilon\eta,\varepsilon\dot{\eta}$から$\varepsilon$が取り除かれて$\eta,\dot{\eta}$が残ります。
- 形式的には$(1)$の$\eta,\dot{\eta}$と$(2)$の$\delta q,\delta \dot q$が対応しますが、前者は任意の関数であり、微小であるという制約は課されません。$\eta$自体を小さくする必要はなく、各$\eta$を固定したうえで$\varepsilon\to0$とすることで、その方向の第一変分が定まります。
- $(2)$では作用の差から2次以上の項を取り除いて$\delta S$とします。これは$(1)$における$\lim_{\varepsilon\to0}\frac{S[q+\varepsilon\eta]-S[q]}{\varepsilon}$に対応し、$\delta q=\varepsilon\eta$とすれば$\delta S=\varepsilon\,\delta S[q;\eta]$となります。

&&&rem εとηの分離による厳密化
関数$\eta$が「微小である」ことを直接定義しようとすると曖昧になりがちですが、実数$\varepsilon$が$0$に近づくことは通常の解析学で厳密に扱えます。$\delta$記法も、作用の差の1次の部分を第一変分として定義すれば厳密に扱えます。$\varepsilon$と$\eta$を分離する意味は、この第一変分を実数の言葉に還元し、通常の実数微分として明示できる点にあります。
&&&

# まとめ

本記事では、経路の変分を実数パラメーター$\varepsilon$と、端点で$0$になる任意の$C^1$級関数$\eta(t)$の積として$q(t)+\varepsilon\eta(t)$と表し、作用積分を$\varepsilon$の関数として扱うことで、オイラー＝ラグランジュ方程式を導出しました。変化後の速度は変化後の経路の時間微分として自動的に決まるため、位置と速度の関係は保たれたままです。停留条件は$\varepsilon=0$における通常の微分係数（第一変分）が$0$になる条件となり、部分積分と変分法の基本補題からオイラー＝ラグランジュ方程式が得られます。$L=T-V$で$T=\frac12m\dot q^2$、$V=V(q)$の場合、この方程式はニュートンの運動方程式$F=ma$と等価です。

物理学でよく用いられる$\delta$による変分の計算とは形式的に対応しますが、$\eta$には微小であるという制約は課されません。変分の微小性は実数$\varepsilon$が$0$に近づく極限として扱われ、第一変分を通常の微分として計算できます。

&&& 作用積分の第一変分
停留経路では、次が成り立ちます。

$$
\delta S[q;\eta]
= \int_{t_1}^{t_2} \left[ \frac{\partial L}{\partial {q}} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}} \right] \eta \, dt = 0
$$
&&&

&&& オイラー＝ラグランジュ方程式
$$
\frac{\partial L}{\partial q} - \frac{d}{dt} \frac{\partial L}{\partial \dot{q}} = 0
$$
&&&
