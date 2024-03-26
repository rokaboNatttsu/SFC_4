- [1. モデルでやりたいことリスト](#1-モデルでやりたいことリスト)
- [2. 行動方程式、上から計算](#2-行動方程式上から計算)
- [3. 定義式](#3-定義式)
  - [3.1. 名目値と実質値と価格の関係の恒等式](#31-名目値と実質値と価格の関係の恒等式)
- [4. 会計恒等式](#4-会計恒等式)
  - [4.1. 取引フロー表](#41-取引フロー表)
  - [4.2. バランスシート表](#42-バランスシート表)
  - [4.3. 完全統合表](#43-完全統合表)
  - [4.4. フローの整合性](#44-フローの整合性)
  - [4.5. ストックとフローの接続の整合性](#45-ストックとフローの接続の整合性)
  - [4.6. ストックの整合性](#46-ストックの整合性)
- [5. パラメータの値](#5-パラメータの値)
- [6. 今後発展させる方向性について](#6-今後発展させる方向性について)


将来的に、余裕があれば、このレジュメに各行動方程式とパラメータの説明を加えたい

# 1. モデルでやりたいことリスト

- 政府支出の増加率の変化の影響をシミュレーションしてまとめる
  - 基準シナリオ
  - 政府支出一定シナリオ
  - 技術進歩がネックにならないシナリオ
  - 技術進歩がネックになるシナリオ
  - 政府支出急増シナリオ
- カレツキアンモデルの賃金主導型と利潤主導型が再現されているかどうか確認
  - 再現されないなら、それはなぜかを考察
  - 賃金主導型成長だけ再現されるなら、その理由を考察
    - 賃金主導型成長が起こる条件を帰納的に探す
  - 利潤主導型成長だけ再現されるなら、その理由を考察
    - 利潤主導型成長が起こる条件を帰納的に探す
  - 現実社会が、利潤主導型と賃金主導型のどちらできれいに説明できるか、文献をあたる


# 2. 行動方程式、上から計算

- $\zeta_2 = \zeta_{2-1} (1 + \zeta_3 \cdot abs(randn()))$
- $W_f = (1.0 - \epsilon_1)W_{f-1} + \epsilon_1 \cdot \max(0.0, \epsilon_3(W_{f-1} + \Pi_{-1} - I_{-1}))$
- $W_b = (1.0 - \epsilon_1)W_{b-1} + \epsilon_1 \cdot \max(0, \epsilon_2(\Pi_{b-1} + r L_{-1}) + \epsilon_4 NW_{b-1})$
- $u^e = ((1-\lambda_e)u^e_{-1} + \lambda_e u_{-1})u_{-1}/(u_{-1}+\Delta u_{-1})$
- $p = p_{-1} \exp\{\mu_3 (u_{-1} - u^T)\}$
- $w_g = w_{g-1}(1 + \delta_1 - \delta_2 \frac{p-p_{-1}}{p_{-1}})$
- if $E_{b-1} < \mu_1(L_{-1} + H_{b-1} + E_{b-1})$
  - $p_e = p_{e-1} \exp(\mu_2 \cdot abs(randn()))$
- else
  - $p_e = p_{e-1} \exp(-\mu_2 \cdot abs(randn()))$
- end
- $g^D = g^D_{-1}(1 + \delta_1 - \delta_2 \frac{p-p_{-1}}{p_{-1}})$
- $C_w^D = (1 - \alpha_5) C_{w-1}^D + \alpha_5 \cdot \max\{0, \alpha_1 (W_{-1}-T_{iw-1}-T_{ew-1}-r L_{w-1}) + \alpha_2 (M_{w-1} + H_{w-1} - L_{w-1})\}$
- $C_i^D = (1 - \alpha_5) C_{i-1}^D + \alpha_5 \cdot \max(0, \alpha_3 (\Pi_{i-1}-T_{ii-1}-T_{ei-1}) + \alpha_4 (M_{i-1} + E_{i-1}))$
- $C_w = \frac{C_w^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{-1}, \zeta_2)})}$
- $C_i = \frac{C_i^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{-1}, \zeta_2)})}$
- $G = \frac{G^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{-1}, \zeta_2)})}$
- $T_{ef} = \gamma_1 K_{-1}$
- $T_{ei} = \gamma_2 (M_{i-1} + E_{i-1})$
- $T_{ew} = \gamma_2 (M_{w-1} + H_{w-1})$
- $i = btw(0, (u_{-1} - u^T)k_{-1} + \beta k_{-1}, \frac{M_{f-1} - L_{f-1}}{p})$
- $T_{ii} = \tau_1 \Pi_{i-1}$
- $T_{iw} = \tau_1 W_{-1}$
- $T_{ff} = \tau_2 (C+I+G-W_f-T_{ef}-r L_{f-1})$
- $\Pi_i = max(0, \theta_1 \frac{E_{i-1}}{E_{-1}} (\Pi - I) + \theta_2 (M_{f-1} - L_{f-1}))$
- $\Pi_b = max(0, \theta_1 \frac{E_{b-1}}{E_{-1}} (\Pi - I) + \theta_2 (M_{f-1} - L_{f-1}))$
- $T_{fb} = max(0, \tau_2 (\Pi_b + r L_{-1} - W_b))$
- $H_w = \iota_1 C_w$
- $L_w = \iota_2 W$
- $M_w = NL_w - \Delta H_w + \Delta L_w$
- $E_i = min(\iota_3 NW_i, E)$
- $M_i = NW_i - E_i$
- $L_f^D = (1 - \iota_4)L_{f-1}^D + \iota_4 max(0, I + W_f + T_{ef} + T_{ff} + r L_{f-1} - M_{f-1} - NL_f)$
- $L_f = btw(0, L_f^D, max(0, \iota_5 NL_f))$
- $\Delta M_f = NL_f + \Delta L_f$ 新規株式発行や自社株買いを行わないという仮定

# 3. 定義式

- $u = \frac{c + g}{\zeta_1 k_{-1}}$
- $NL_w = -C_w+W-T_{iw}-T_{ew}-r L_{w-1}$
- $NL_i = -C_i-T_{ii}-T_{ei}+\Pi_i$
- $NL_f = -I + \Pi_f$
- $NL_b = -W_b - T_{fb} + \Pi_b + r L_{-1}$
- $NL_g = -G-W_g+T_i+T_e+T_f$

## 3.1. 名目値と実質値と価格の関係の恒等式

- $G = p g$
- $C_w = p c_w$
- $C_i = p c_i$
- $K = p k$
- $I = p i$
- $E_i = p_e e_i$
- $E = p_e e$
- $E_b = p_e e_b$

# 4. 会計恒等式

## 4.1. 取引フロー表

|                      |    労働者     |      資本家       |  企業(経常)  |  企業(資本)   |       銀行        |  統合政府   | 合計 |
| :------------------- | :-----------: | :---------------: | :----------: | :-----------: | :---------------: | :---------: | :--: |
| 消費                 |    $-C_w$     |      $-C_i$       |     $+C$     |               |                   |             | $0$  |
| 投資                 |               |                   |     $+I$     |     $-I$      |                   |             | $0$  |
| 政府支出（賃金除く） |               |                   |     $+G$     |               |                   |    $-G$     | $0$  |
| 賃金                 |     $+W$      |                   |    $-W_f$    |               |      $-W_b$       |   $-W_g$    | $0$  |
| 所得税               |   $-T_{iw}$   |     $-T_{ii}$     |              |               |                   |   $+T_i$    | $0$  |
| 資産税               |   $-T_{ew}$   |     $-T_{ei}$     |  $-T_{ef}$   |               |                   |   $+T_e$    | $0$  |
| 法人税               |               |                   |  $-T_{ff}$   |               |     $-T_{fb}$     |   $+T_f$    | $0$  |
| 企業利潤             |               |     $+\Pi_i$      |    $-\Pi$    |   $+\Pi_f$    |     $+\Pi_b$      |             | $0$  |
| 借入金金利           | $-r L_{w-1}$  |                   | $-r L_{f-1}$ |               |    $+r L_{-1}$    |             | $0$  |
| [メモ:Net Lending]   |    $NL_w$     |      $NL_i$       |              |    $NL_f$     |      $NL_b$       |   $NL_g$    | $0$  |
| 預金                 | $-\Delta M_w$ |   $-\Delta M_i$   |              | $-\Delta M_f$ |    $+\Delta M$    |             | $0$  |
| 借入                 | $+\Delta L_w$ |                   |              | $+\Delta L_f$ |    $-\Delta L$    |             | $0$  |
| 株式                 |               | $-p_e \Delta e_i$ |              |               | $-p_e \Delta e_b$ |             | $0$  |
| 現金                 | $-\Delta H_w$ |                   |              |               |   $-\Delta H_b$   | $+\Delta H$ | $0$  |
| 合計                 |      $0$      |        $0$        |     $0$      |      $0$      |        $0$        |     $0$     |      |

## 4.2. バランスシート表

|        | 労働者  | 資本家  |  企業   |  銀行   | 統合政府 | 合計 |
| :----- | :-----: | :-----: | :-----: | :-----: | :------: | :--: |
| 資本   |         |         |  $+K$   |         |          | $+K$ |
| 預金   | $+M_w$  | $+M_i$  | $+M_f$  |  $-M$   |          | $0$  |
| 株式   |         | $+E_i$  |  $-E$   | $+E_b$  |          | $0$  |
| 借入   | $-L_w$  |         | $-L_f$  |  $+L$   |          | $0$  |
| 現金   | $+H_w$  |         |         | $+H_b$  |   $-H$   | $0$  |
| 純資産 | $-NW_w$ | $-NW_i$ | $-NW_f$ | $-NW_b$ | $-NW_g$  | $-K$ |
| 合計   |   $0$   |   $0$   |   $0$   |   $0$   |   $0$    | $0$ |

## 4.3. 完全統合表

|                        |    労働者     |           資本家            |            企業            |            銀行             |  統合政府   |           合計           |
| :--------------------- | :-----------: | :-------------------------: | :------------------------: | :-------------------------: | :---------: | :----------------------: |
| 期首純資産             |  $NW_{w-1}$   |         $NW_{i-1}$          |         $NW_{f-1}$         |         $NW_{b-1}$          | $NW_{g-1}$  |         $K_{-1}$         |
| 資本のキャピタルゲイン |               |                             |  $+\Delta p \cdot k_{-1}$  |                             |             | $+\Delta p \cdot k_{-1}$ |
| 資本の増減             |               |                             |    $+p \cdot \Delta k$     |                             |             |   $+p \cdot \Delta k$    |
| 預金の増減             | $+\Delta M_w$ |        $+\Delta M_i$        |       $+\Delta M_f$        |         $-\Delta M$         |             |           $0$            |
| 株式のキャピタルゲイン |               | $+\Delta p_e \cdot e_{i-1}$ | $-\Delta p_e \cdot e_{-1}$ | $+\Delta p_e \cdot e_{b-1}$ |             |           $0$            |
| 株式の増減             |               |   $+p_e \cdot \Delta e_i$   |                            |   $+p_e \cdot \Delta e_b$   |             |           $0$            |
| 借入の増減             | $-\Delta L_w$ |                             |       $-\Delta L_f$        |         $+\Delta L$         |             |           $0$            |
| 現金の増減             | $+\Delta H_w$ |                             |                            |        $+\Delta H_b$        | $-\Delta H$ |           $0$            |
| 期末純資産             |    $NW_w$     |           $NW_i$            |           $NW_f$           |           $NW_b$            |   $NW_g$    |           $K$            |

## 4.4. フローの整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $-C_w+W-T_{iw}-T_{ew}-r L_{w-1} = \Delta M_w - \Delta L_w + \Delta H_w$
- [x] $-C_i-T_{ii}-T_{ei}+P_i = \Delta M_i + p_e \Delta e_i$
- [x] $\Pi = C+I+G-W_f-T_{ef}-T_{ff}-r L_{f-1}$
- [x] $-I + \Pi_f = \Delta M_f - \Delta L_f$
- [x] $-W_b - T_{fb} + \Pi_b + r L_{-1} = -\Delta M + \Delta L + p_e \Delta e_b + \Delta H_b$
- [ ] $-G-W_g+T_i+T_e+T_f = -\Delta H$
- [x] $C = C_w + C_i$
- [x] $W = W_f + W_g + W_b$
- [x] $T_i = T_{iw} + T_{ii}$
- [x] $T_e = T_{ew} + T_{ei} + T_{ef}$
- [x] $T_f = T_{ff} + T_{fb}$
- [x] $\Pi_f = \Pi - \Pi_i - \Pi_b$
- [ ] $NL_w + NL_i + NL_f + NL_b + NL_g = 0$
- [x] $\Delta M = \Delta M_w + \Delta M_i + \Delta M_f$
- [x] $\Delta L = \Delta L_w + \Delta L_f$
- [ ] $\Delta e_i + \Delta e_b = 0$
- [x] $\Delta H = \Delta H_w + \Delta H_b$

## 4.5. ストックとフローの接続の整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $k = (1 - \beta)k_{-1} + i$
- [x] $\Delta k = k - k_{-1}$
- [x] $K = K_{-1} + \Delta K = p_{-1} k_{-1} + \Delta p \cdot k_{-1} + p \cdot \Delta k$
- [x] $\Delta M_w = M_w - M_{w-1}$
- [x] $\Delta M_i = M_i - M_{i-1}$
- [x] $M_f = M_{f-1} + \Delta M_f$
- [x] $M = M_{-1} + \Delta M$
- [x] $\Delta L_w = L_w - L_{w-1}$
- [x] $\Delta L_f = L_f - L_{f-1}$
- [ ] $L = L_{-1} + \Delta L$
- [x] $E_i = E_{i-1} + \Delta E_i = p_{e-1} e_{i-1} + \Delta p_e \cdot e_{i-1} + p_e \cdot \Delta e_i$
- [x] $E_b = E_{b-1} + \Delta E_b = p_{e-1} e_{b-1} + \Delta p_e \cdot e_{b-1} + p_e \cdot \Delta e_b$
- [x] $E = E_{-1} + \Delta E = p_{e-1} e_{-1} + \Delta p_e \cdot e_{-1} + p_e \cdot \Delta e$
- [x] $\Delta H_w = H_w - H_{w-1}$
- [x] $H_b = H_{b-1} + \Delta H_b$
- [x] $H = H_{-1} + \Delta H$
- [x] $NW_w = NW_{w-1} + \Delta M_w - \Delta L_w + \Delta H_w (=NW_{w-1} + NL_w)$
- [x] $NW_i = NW_{i-1} + \Delta M_i + \Delta p_e \cdot e_{i-1}  + p_e \cdot \Delta e_i (= NW_{i-1} + NL_i + \Delta p_e \cdot e_{i-1})$
- [x] $NW_f = NW_{f-1} + \Delta K + \Delta M_f - \Delta E - \Delta L_f (= NW_{f-1} + NL_f - \Delta E + \Delta K)$
- [x] $NW_b = NW_{b-1} - \Delta M + \Delta E_b + \Delta L + \Delta H_b (= NW_{b-1} + NL_b + \Delta p_e \cdot e_{b-1})$
- [x] $NW_g = NW_{g-1} - \Delta H (= NW_{g-1} + NL_g)$

## 4.6. ストックの整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $M = M_w + M_i + M_f$
- [ ] $E = E_i + E_b$
- [x] $e = e_i + e_b$
- [x] $L = L_w + L_f$
- [x] $H = H_w + H_f$
- [ ] $NW_w + NW_i + NW_f + NW_b + NW_g = 0$
- [ ] $NW_w = M_w - L_w + H_w$
- [ ] $NW_i = M_i + E_i$
- [ ] $NW_f = K + M_f - E - L_f$
- [ ] $NW_b = -M + E_b + L + H_b$
- [ ] $NW_g = -H$

# 5. パラメータの値

- $\alpha_1 = 0.9$
- $\alpha_2 = 0.1$
- $\alpha_3 = 0.1$
- $\alpha_4 = 0.05$
- $\alpha_5 = 0.5$
- $\beta = 0.05$
- $\gamma_1 = 0.015$
- $\gamma_2 = 0.02$
- $\delta_1 = 0.02$
- $\delta_2 = 0.2$
- $\epsilon_1 = 1.0$
- $\epsilon_2 = 0.7$
- $\epsilon_3 = 0.7$
- $\epsilon_4 = 0.02$
- $\zeta_1 = 1.0$
- $\zeta_2 = 300.0$ (シミュレーションのための初期値)
- $\zeta_3 = 0.0232$
- $\theta_1 = 0.2$
- $\theta_2 = 0.1$
- $\iota_1 = 0.1$
- $\iota_2 = 1.0$
- $\iota_3 = 0.5$
- $\iota_4 = 0.5$
- $\iota_5 = 10.0$
- $\mu_1 = 0.3$
- $\mu_2 = 0.05$
- $\mu_3 = 0.5$
- $\tau_1 = 0.3$
- $\tau_2 = 0.2$
- $mu = 0.2$
- $u^T = 0.8$
- $G_0 = 100$
- $r = 0.01$

# 6. 今後発展させる方向性について

- 期待値を導入する。期待値は、 $x^e=((1 - \lambda_e) x^e_{-1} + \lambda_e x_{-1})\frac{x_{-1}}{x_{-1}-\Delta x_{-1}}$ の形で計算する。経済が安定して拡大しているときは消費性向が増加し、経済が停滞しているときは消費性向が減少する傾向が、生み出される見込み
- パラメータのカリブレーションと推定をする
- 新規株式発行と自社株買いを追加する。売上規模(名目値)と発行部数が比例する状態が望ましいかもしれない。企業の資金調達方法が追加される
- 時価変動と配当であらわされる収益率が高いほど、ポートフォリオ配分割合が増えるような行動方程式を採用する。
- 政府部門内で資本ストックを作る。公的固定資本形成を作る
- 技術投資や研究投資（Spring-8みたいなの作ったり、学者に給料を払ったり）の増加が、技術進歩速度を上げることを記述できるように、政府と企業の投資と賃金にそれ専用の変数を導入する
- AB化する
