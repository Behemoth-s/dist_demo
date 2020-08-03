# 精馏 模拟 面向方程demo
# 两组分模型
# 相平衡使用Wilson 方程计算活度系数
# 使用Antoine 和 Wilson方程共同构成对温度的约束
# 本模型为模拟模型
# 目标为 每一块塔板上的气液相差为0
using Ipopt
using JuMP
using MathOptInterface
const MOI = MathOptInterface
x_ss = [0.97287970129754
   0.95636038934316
   0.93661040294083
   0.91321282273197
   0.88585270419340
   0.85442258268073
   0.81914368999749
   0.78066600262237
   0.74009371717436
   0.69889223479787
   0.65867789383232
   0.62095142736705
   0.58686770454651
   0.55711313643268
   0.53190656888283
   0.51108997385426
   0.49425664303546
   0.47785779539230
   0.45662374470122
   0.42987184311163
   0.39729919963226
   0.35923134730789
   0.31678635583970
   0.27183384568409
   0.22669854804022
   0.18369476978786
   0.14468004616214
   0.11079600651024
   0.08244638861628
   0.05944779752760
   0.04124685281752
   0.02712029870246]
T_ss = 100 * [3.54170894061095
    3.54384309151657
    3.54642323740009
    3.54952180753185
    3.55320509380443
    3.55751925276361
    3.56247139366757
    3.56800942498801
    3.57400825826251
    3.58027163739042
    3.58655486133578
    3.59260419440273
    3.59819966694026
    3.60318627239385
    3.60748534712935
    3.61108752771922
    3.61403480708638
    3.61693563210643
    3.62073526774064
    3.62559227397992
    3.63161163785833
    3.63879364907630
    3.64698872172953
    3.65588344068293
    3.66503752584723
    3.67396743659234
    3.68224472290915
    3.68956911381347
    3.69579402098210
    3.70090876480212
    3.70499767004852
    3.70819628750959]

nh = 32 # 塔板数 包括冷凝器和再沸器
ne = 2 # 组分数
P = 101325 # 塔压，本案例未考虑压降
# P = [ Ps0 + dp * (i-1) for i in 1:nh]
Feed =  24.0 / 60.0 # 进料量
x_Feed = 0.5 # 进料组成， 平衡态进料
D = 0.5 * Feed # 塔顶采出量
rr = 2.49 # 回流比

L = rr * D     # 精馏液相流量
V = L + D      # 精馏气相流量
FL = Feed + L  # 提馏段液相流量

# antoine系数
DIPPR = [   5.1087E1 8.7829E1;
            -5.2264E3 -6.9964E3;
            -4.2278E0 -9.8802E0;
            9.7554E-18 7.2099E-6;
            6.0000E0 2.0000E0];

# Wilson系数
L12 = 1.618147731;  # 组分1对组分2
L21 = 0.50253532;   # 组分2对组分1

nfeed = 17  # 进料板

# 创建模型
m = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 600,  "tol" => 1e-6, "print_level" => 5,  "output_file" => "ipopt.out", "file_print_level" => 6))

# 变量个数 nh * 7  = 224
@variable(m, 0 <= x[i in 1:nh] <= 1)  # 组分A的液相浓度
@variable(m, 0 <= y[i in 1:nh] <= 1)  # 组分A的气相浓度
@variable(m, T[i in 1:nh])            # 温度
@variable(m, gammaA[i in 1:nh])       # 组分A的活度系数
@variable(m, gammaB[i in 1:nh])       # 组分B的活度系数
@variable(m, PsatA[i in 1:nh])        # 组分A的饱和蒸汽压
@variable(m, PsatB[i in 1:nh])        # 组分B的饱和蒸汽压

# 方程数 nh* 4  = 192
# 其中目标函数包含 32 个 eqdt == 0
# 192 + 32 = 224

# 使用Wilson 方程计算活度系数, A表示组分A, 及x和y对应的组分
@NLconstraint(m, defgammaA[i in 1:nh], gammaA[i] == exp(-log(x[i] + L12 * (1 - x[i])) + (1 - x[i]) * (L12 / (x[i] + L12 * (1 - x[i])) - (L21 / (L21 * x[i] + (1 - x[i]))))))
@NLconstraint(m, defgammaB[i in 1:nh], gammaB[i] == exp(-log((1 - x[i]) + L21 * x[i]) + x[i] * (L21 / ((1 - x[i]) + L21 * x[i]) - (L12 / (L12 * (1 - x[i]) + x[i])))))

# 使用Antoine方程计算饱和蒸汽压, A表示组分A, 及x和y对应的组分
@NLconstraint(m, defPstaA[i in 1:nh], PsatA[i] == exp(DIPPR[1,1] + DIPPR[2,1] / T[i] + DIPPR[3,1] * log(T[i]) + DIPPR[4,1] * (T[i]^DIPPR[5,1])))
@NLconstraint(m, defPstaB[i in 1:nh], PsatB[i] == exp(DIPPR[1,2] + DIPPR[2,2] / T[i] + DIPPR[3,2] * log(T[i]) + DIPPR[4,2] * (T[i]^DIPPR[5,2])))

# 根据液相活度和Ps, 计算组分分压 x -> y
@NLconstraint(m, defy[i in 1:nh],  y[i] == x[i] * gammaA[i] * (PsatA[i] / P))

# 组分A的质量平衡, D冷凝器, dist精馏段, feed进料, strip提馏段, W 再沸器
@NLconstraint(m, mass_eq_D,y[2] == x[1])
@NLconstraint(m, mass_eq_dist[i in 2:nfeed], L * (x[i - 1] - x[i]) == V * (y[i] - y[i + 1]))
@NLconstraint(m, mass_eq_feed, Feed * x_Feed + L * x[nfeed - 1] == FL * x[nfeed] - V * (y[nfeed] - y[nfeed + 1]))
@NLconstraint(m, mass_eq_strip[i in nfeed + 1:nh - 1], FL * (x[i - 1] - x[i]) == V * (y[i] - y[i + 1]))
@NLconstraint(m, mass_eq_W, FL * x[nh - 1] - (Feed - D) * x[nh] - V * y[nh] == 0)

# Pa + Pb == P
@NLexpression(m, eqdt[i in 1:nh], ((x[i] * gammaA[i] * PsatA[i]) + ((1 - x[i]) * gammaB[i] * PsatB[i]) - P) / P)

# 目标 min sum(dp)
@NLobjective(m, Min, sum(eqdt[i]^2 for i in 1:nh))

# 设置初值
set_start_value.(x, x_ss)
set_start_value.(T, T_ss)

optimize!(m)
print(DataFrame([JuMP.value.(x),JuMP.value.(y), JuMP.value.(T)], [:x,:y,:T]))